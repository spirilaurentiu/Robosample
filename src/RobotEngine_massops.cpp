// ============================================================================
//  RobotEngine_massops.cpp - mass-matrix operators shared by HMC momentum
//  seeding, kinetic energy, and the Fixman log-det (INV-5).
//
//  One of four cohesive splits of the former RobotEngine.cpp (SPLIT-R3, pure
//  code motion); see RobotEngine_kinematics.cpp's header comment for the
//  Simbody provenance / validation notes shared by all four
//  RobotEngine_*.cpp translation units.
// ============================================================================

#include "RobotEngine.hpp"

#include <vector>

#include "RobotEngine_internal.hpp"
#include "math/hinge_linalg.hpp"
#include "robot_math.hpp"

using robo::ArticulatedInertia;
using robo::Mat33;
using robo::PhiMatrix;
using robo::Quaternion;
using robo::Real;
using robo::Rotation;
using robo::SpatialInertia;
using robo::SpatialVec;
using robo::SymMat33;
using robo::Transform;
using robo::Vec3;
using robo::Vec4;

using robo::detail::pseudoLogDet;
using robo::detail::spatialDot;
using robo::detail::symSqrt;
using robo::detail::symSqrtInv;

// ============================================================================
//  M^-1 f   (two passes, no velocity terms)   RigidBodyNodeSpec.cpp
// ============================================================================
void RobotEngine::multiplyByMInv(const RobotModel& m, RobotState& s, const Real* in, Real* out) {
    const SpatialVec* H = s.H();
    const SpatialVec* G = s.G();
    const Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    SpatialVec* Z = s.Z();
    SpatialVec* ZPlus = s.zPlus();
    Real* eps = s.eps();
    SpatialVec* A_GB = s.A_GB();

    for (int b = m.numBodies - 1; b >= 1; --b) {
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        SpatialVec z(Vec3(0), Vec3(0));
        for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
            const int c = m.bodyChildren[ci];
            z += Phi[c] * ZPlus[c];
        }
        SpatialVec zp = z;
        for (int j = 0; j < dof; ++j) {
            eps[uOff + j] = in[uOff + j] - spatialDot(H[uOff + j], z);
        }
        for (int j = 0; j < dof; ++j) {
            zp += G[uOff + j] * eps[uOff + j];
        }
        Z[b] = z;
        ZPlus[b] = zp;
    }

    A_GB[0] = SpatialVec(Vec3(0), Vec3(0));
    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const SpatialVec APlus = (~Phi[b]) * A_GB[p];
        const Real* DI = &DIpool[m.bodyUSqIndex[b]];
        SpatialVec acc = APlus;
        for (int i = 0; i < dof; ++i) {
            Real v = 0;
            for (int j = 0; j < dof; ++j) {
                v += DI[i * dof + j] * eps[uOff + j];
            }
            v -= spatialDot(G[uOff + i], APlus);
            out[uOff + i] = v;
            acc += H[uOff + i] * v;
        }
        A_GB[b] = acc;
    }
}

// ============================================================================
//  sqrt(M^-1) f   (HMC velocity seeding)   RigidBodyNodeSpec.cpp
// ============================================================================
void RobotEngine::multiplyBySqrtMInv(const RobotModel& m, RobotState& s, const Real* in, Real* out) {
    const SpatialVec* H = s.H();
    const SpatialVec* G = s.G();
    const Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    SpatialVec* V = s.V_GB(); // reuse as scratch outward velocity

    V[0] = SpatialVec(Vec3(0), Vec3(0));
    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const SpatialVec VPlus = (~Phi[b]) * V[p];
        Real sqrtDI[36];
        symSqrt(&DIpool[m.bodyUSqIndex[b]], dof, sqrtDI);
        SpatialVec acc = VPlus;
        for (int i = 0; i < dof; ++i) {
            Real v = 0;
            for (int j = 0; j < dof; ++j) {
                v += sqrtDI[i * dof + j] * in[uOff + j];
            }
            v -= spatialDot(G[uOff + i], VPlus);
            out[uOff + i] = v;
            acc += H[uOff + i] * v;
        }
        V[b] = acc;
    }
}

// ============================================================================
//  sqrt(M) f   (exact inverse of multiplyBySqrtMInv; for NMA Route B acceptance)
// ============================================================================
// multiplyBySqrtMInv maps in -> out = M^(-1/2) in via the outward recurrence
//     VPlus_b = ~Phi_b V_p ;  out_b = sqrtDI_b in_b - ~G_b VPlus_b ;
//     V_b     = VPlus_b + H_b out_b .
// Inverting body-by-body (out_b is this sweep's INPUT, in_b its OUTPUT):
//     in_b = sqrtDI_b^{-1} ( out_b + ~G_b VPlus_b ) ,
// while V_b is rebuilt from the SAME (given) out_b, so V_b is identical to the
// forward sweep's. Substituting the forward out_b shows in_b is recovered
// exactly: sqrtDI^{-1}(sqrtDI in_b - ~G VPlus + ~G VPlus) = in_b. Hence this is
// the exact algebraic inverse for any H,G,DI,Phi -- no approximation. It uses a
// LOCAL scratch velocity (NOT s.V_GB()), so it is safe to call after a trajectory
// when V_GB holds the real end velocities needed by calcKineticEnergy.
void RobotEngine::multiplyBySqrtM(const RobotModel& m, RobotState& s, const Real* in, Real* out) {
    const SpatialVec* H = s.H();
    const SpatialVec* G = s.G();
    const Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    std::vector<SpatialVec> V(m.numBodies, SpatialVec(Vec3(0), Vec3(0)));

    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const SpatialVec VPlus = (~Phi[b]) * V[p];
        Real sqrtDIinv[36];
        symSqrtInv(&DIpool[m.bodyUSqIndex[b]], dof, sqrtDIinv);
        // tmp_i = in_b[i] + ~G_b[i] . VPlus     (in[] plays the forward out_b)
        Real tmp[6];
        for (int i = 0; i < dof; ++i) {
            tmp[i] = in[uOff + i] + spatialDot(G[uOff + i], VPlus);
        }
        // out_b = sqrtDI^{-1} tmp
        for (int i = 0; i < dof; ++i) {
            Real v = 0;
            for (int j = 0; j < dof; ++j) {
                v += sqrtDIinv[i * dof + j] * tmp[j];
            }
            out[uOff + i] = v;
        }
        // V_b = VPlus + H_b in_b   (rebuilt from the given out_b == in[], matching forward)
        SpatialVec acc = VPlus;
        for (int i = 0; i < dof; ++i) {
            acc += H[uOff + i] * in[uOff + i];
        }
        V[b] = acc;
    }
}

robo::Real RobotEngine::calcLogDetM(const RobotModel& m, const RobotState& s) {
    using robo::Real;
    using robo::SpatialVec;
    const SpatialVec* H = s.H();
    const robo::ArticulatedInertia* P = s.P();

    Real logDet = Real(0);
    Real D[36]; // per-body dof x dof, dof <= 6; the global M is never formed
    for (int b = 1; b < m.numBodies; ++b) {
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        if (dof == 0) {
            continue; // Weld: contributes det = 1 -> log 0
        }
        // D_b = ~H_b * P_b * H_b  (exactly as in realizeArticulatedBodyInertias)
        for (int j = 0; j < dof; ++j) {
            const SpatialVec PHj = P[b] * H[uOff + j];
            for (int i = 0; i < dof; ++i) {
                D[i * dof + j] = spatialDot(H[uOff + i], PHj);
            }
        }
        logDet += pseudoLogDet(D, dof); // pseudo-ln-det(D_b): null directions contribute 0 (CC1)
    }
    return logDet; // = ln|M_phi|
}

// ============================================================================
//  KINETIC ENERGY = 1/2 u^T M u  via the spatial sweep (Mk_G, V_GB).
// ============================================================================
Real RobotEngine::calcKineticEnergy(const RobotModel& m, const RobotState& s) {
    const SpatialInertia* Mk = s.Mk_G();
    const SpatialVec* V = s.V_GB();
    Real ke = 0;
    for (int b = 1; b < m.numBodies; ++b) {
        const SpatialVec MV = Mk[b] * V[b];
        ke += spatialDot(V[b], MV);
    }
    return Real(0.5) * ke;
}
