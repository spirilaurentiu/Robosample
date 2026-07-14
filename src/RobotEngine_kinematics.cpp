// ============================================================================
//  RobotEngine_kinematics.cpp - SimTK-free articulated-body dynamics:
//  position/velocity kinematics and the q<->u differential maps.
//
//  This is a faithful port of the Simbody recursions (no MobilizedBody /
//  SimbodyMatterSubsystem / State / vtables). The SoA driver walks the body
//  tree (model.bodyParent[b] < b: outward = forward, inward = reverse) and
//  branches on model.bodyJoint[b] instead of virtual dispatch. The spatial
//  *value* algebra (SpatialVec, ArticulatedInertia, PhiMatrix, Transform,
//  SpatialInertia) is SimTK's own header-only, vtable-free code, so only the
//  recursion is reimplemented.
//
//  Source map (Simbody01/Simbody/src):
//    realizePosition / realizeVelocity ......... RigidBodyNodeSpec.h
//    calcBodyTransforms / H_PB_G ............... RigidBodyNodeSpec.{h,cpp}
//    Phi / gyroscopic / coriolis .............. RigidBodyNode.cpp
//    realizeArticulatedBodyInertias ........... RigidBodyNodeSpec.cpp
//    calcUDot pass1/2 ......................... RigidBodyNodeSpec.cpp
//    multiplyByMInv pass1/2 ................... RigidBodyNodeSpec.cpp
//    multiplyBySqrtMInvPassOutward ............ RigidBodyNodeSpec.cpp
//    per-joint H_FM / X_FM / qdot ............. RigidBodyNodeSpec_{Torsion,Translation,Free}.h
//
//  VALIDATION: Simbody has been REMOVED from the build (no find_package, no
//  SimTK/ Simbody01 tree, no live MobilizedBody code anywhere in src/ or tests/).
//  The historical "diff every operator against the still-present Simbody build"
//  gate therefore NO LONGER EXISTS and was never wired into the suite. The
//  oracle-of-record for these recursions is now:
//    * analytic closed forms (single/free rigid body, physical pendulum, KE),
//    * finite-difference cross-checks (H<->X_FM, V<->X_GB, A<->V, HDot<->H,
//      logDetM<->Cholesky of the forward-built dense M),
//    * operator round-trips (M*M^-1==I, sqrtMInv/sqrtM mutual inverse,
//      u^T M u == w^T w for u = multiplyBySqrtMInv(w)),
//    * statistical ensembles (equipartition, Fixman/Boltzmann marginal),
//    * frozen Simbody-golden CONSTANTS transcribed into TestSpatialAlgebra /
//      TestInertia (numbers, not a live diff), and
//    * a live OpenMM/ParmEd differential for POTENTIAL energy only.
//  Field/routine names below still track their Simbody originals (source map
//  above) purely as provenance, not as a runtime comparison.
//
//  This TU is one of four cohesive splits of the former RobotEngine.cpp
//  (SPLIT-R3, pure code motion): RobotEngine_kinematics.cpp (this file,
//  carries the provenance/validation note above for all four),
//  RobotEngine_dynamics.cpp, RobotEngine_massops.cpp, RobotEngine_reaction.cpp.
//  All four compile into the SAME RobotEngine class declared in RobotEngine.hpp.
// ============================================================================

#include "RobotEngine.hpp"

#include "JointKernels.hpp"
#include "RobotEngine_internal.hpp"
#include "RobotModel.hpp"
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

// ============================================================================
//  KINEMATICS
// ============================================================================
void RobotEngine::realizePosition(const RobotModel& m, RobotState& s) {
    Transform* X_GB = s.X_GB();
    Transform* X_FM = s.X_FM();
    Transform* X_PB = s.X_PB();
    PhiMatrix* Phi = s.Phi();
    SpatialInertia* Mk = s.Mk_G();
    Vec3* comG = s.comG();
    SpatialVec* HFM = s.H_FM();
    SpatialVec* H = s.H();
    const Real* q = s.q();

    X_GB[0] = Transform(); // Ground
    Mk[0] = SpatialInertia();

    // Kinetic-metric preconditioning: fictitious per-body inertia scale applied
    // to Mk_G ONLY (draw + KE + Fixman all read Mk_G, so they stay consistent).
    // Empty => all 1.0 (physical). Hoisted out of the hot loop.
    const Real* massScale = m.bodyMassScale.empty() ? nullptr : m.bodyMassScale.data();

    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const JointType jt = m.bodyJoint[b];
        const int qOff = m.bodyQIndex[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];

        // X_FM, then X_PB = X_PF * X_FM * X_MB, X_GB = X_GP * X_PB.
        X_FM[b] = robo::jointX_FM(jt, q, qOff);
        const Transform X_MB = ~m.X_BM[b];
        const Transform X_FB = X_FM[b] * X_MB;
        X_PB[b] = m.X_PF[b] * X_FB;
        X_GB[b] = X_GB[p] * X_PB[b];

        // H_FM (in F) and H_PB_G (in Ground): RigidBodyNodeSpec.cpp.
        jointH_FM(jt, X_FM[b], &HFM[uOff]);
        const Rotation R_GF = X_GB[p].R() * m.X_PF[b].R();
        // r_MB is the vector from Mo to Bo expressed in M, i.e. X_MB.p() (Simbody
        // getX_MB().p()), NOT X_BM.p(). X_BM.p() is Mo's position in B (the
        // vector Bo->Mo expressed in B): wrong frame AND wrong sign. Using it
        // corrupts the LINEAR rows of the hinge matrix for every body whose M
        // frame is offset from Bo (all non-root torsion bodies, offset = bond
        // length), giving wrong cross-mobilizer velocities -> non-conserving
        // dynamics. The position path never touches H, so the OpenMM energy
        // check cannot catch this. X_MB == ~X_BM is already computed above.
        const Vec3 r_MB = X_MB.p();
        const Rotation& R_FM = X_FM[b].R();
        const Vec3 r_MB_F = R_FM * r_MB;
        for (int j = 0; j < dof; ++j) {
            const SpatialVec& h = HFM[uOff + j];
            // H_MB_F: top row 0, bottom row = -r_MB_F % h.angular.
            const SpatialVec hpb(h[0], h[1] - (r_MB_F % h[0]));
            H[uOff + j] = SpatialVec(R_GF * hpb[0], R_GF * hpb[1]);
        }

        // Joint-independent kinematics (Phi, COM_G, Mk_G).  RigidBodyNode.cpp.
        const Vec3 p_PB_G = X_GB[p].R() * X_PB[b].p();
        Phi[b] = PhiMatrix(p_PB_G);
        const Rotation& R_GB = X_GB[b].R();
        const Vec3& p_GB = X_GB[b].p();
        const robo::UnitInertia G_Bo_G = m.bodyUnitInertia_B[b].reexpress(~R_GB);
        const Vec3 p_BBc_G = R_GB * m.bodyCom_B[b];
        comG[b] = p_GB + p_BBc_G;
        // bodyMass is the only scaled argument: UnitInertia is mass-normalised, so
        // mass*scale propagates the factor through the whole spatial inertia (mass,
        // first moment, and inertia tensor) -> a uniformly heavier, SPD-valid copy
        // of the same rigid body. Bias-free (see RobotModel::bodyMassScale).
        const Real massB = massScale ? (m.bodyMass[b] * massScale[b]) : m.bodyMass[b];
        Mk[b] = SpatialInertia(massB, p_BBc_G, G_Bo_G);
    }

    // Per-atom Ground positions and stations (R_GB * station_B). The transfer
    // currency + the input to the OpenMM force bridge. Each atom is independent (no
    // reduction), so this vectorizes cleanly: __restrict removes the aliasing barrier
    // between the destination arrays and the (const) source arrays, and `omp simd`
    // (SIMD-only, -fopenmp-simd) lets the compiler pack the R*station FMAs across atoms.
    // Bitwise-identical to the scalar form.
    Vec3* __restrict posG = s.atomPosG();
    Vec3* __restrict stG = s.atomStationG();
    const int* __restrict atomBody = m.atomBody.data();
    const Vec3* __restrict station = m.atomStation_B.data();
#pragma omp simd
    for (int a = 0; a < m.numAtoms; ++a) {
        const int b = atomBody[a];
        const Vec3 st = X_GB[b].R() * station[a];
        stG[a] = st;
        posG[a] = X_GB[b].p() + st;
    }
}

void RobotEngine::realizeVelocity(const RobotModel& m, RobotState& s) {
    calcQDot(m, s, s.qdot());

    const SpatialVec* H = s.H();
    const SpatialVec* HFM = s.H_FM();
    SpatialVec* V_FM = s.V_FM();
    SpatialVec* V_PB_G = s.V_PB_G();
    SpatialVec* V_GB = s.V_GB();
    SpatialVec* gyro = s.gyro();
    SpatialVec* a_tot = s.coriolisA();
    SpatialVec* a_mob = s.mobCoriolisA();
    const PhiMatrix* Phi = s.Phi();
    const SpatialInertia* Mk = s.Mk_G();
    const Transform* X_GB = s.X_GB(); // for R_GF (centripetal HDot_MB_F term)
    const Transform* X_FM = s.X_FM(); // for r_MB_F  (centripetal HDot_MB_F term)
    const Real* u = s.u();

    V_GB[0] = SpatialVec(Vec3(0), Vec3(0));
    a_tot[0] = SpatialVec(Vec3(0), Vec3(0));

    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const JointType jt = m.bodyJoint[b];

        SpatialVec vfm(Vec3(0), Vec3(0));
        SpatialVec vpb(Vec3(0), Vec3(0));
        for (int j = 0; j < dof; ++j) {
            vfm += HFM[uOff + j] * u[uOff + j];
            vpb += H[uOff + j] * u[uOff + j];
        }
        V_FM[b] = vfm;
        V_PB_G[b] = vpb;

        // V_GB = ~Phi * V_GP + V_PB_G.
        const SpatialVec V_GBb = (~Phi[b]) * V_GB[p] + vpb;
        V_GB[b] = V_GBb;
        const Vec3& w_GB = V_GBb[0];
        const Vec3& v_GB = V_GBb[1];

        // Gyroscopic force b = mass * [ w x (G w) ; w x (w x c) ]  (Ground),
        // with G the unit inertia about Bo and c the COM offset (RigidBodyNode.cpp).
        const SpatialInertia& Mkb = Mk[b];
        const Vec3 Iw = Mkb.getUnitInertia() * w_GB;
        const Vec3 c = Mkb.getMassCenter();
        gyro[b] = Mkb.getMass() * SpatialVec(w_GB % Iw, w_GB % (w_GB % c));

        // Mobilizer coriolis (bias) acceleration  a_mob = d(~Phi)/dt V_GP + Hdot_G u.
        //
        // EARLIER BUG: this dropped Hdot_G u, arguing H_FM is constant in the M
        // frame. Wrong frame. H_PB_G is H_FM expressed in GROUND, and the
        // inboard frame F is fixed in the PARENT body, so R_GF (hence H_PB_G)
        // rotates with the parent angular velocity w_GP. Thus
        //   Hdot_G u = w_GP x (H_PB_G u) = w_GP x V_PB_G   (NONZERO).
        // Omitting it drops the angular Coriolis coupling entirely and half the
        // linear Coriolis, so the integrator stops conserving energy: velocities
        // grow every step from a thermal start until the trajectory explodes.
        // Full bias (equivalent to the classic centripetal + 2*Coriolis form
        //  w_GP x (w_GP x p_PB) + 2 w_GP x v_PB_G):
        //   a_mob = [ w_GP x w_PB_G ;
        //             w_GP x (v_GB - v_GP) + w_GP x v_PB_G ]
        const Vec3& w_GP = V_GB[p][0];
        const Vec3& v_GP = V_GB[p][1];
        const Vec3& w_PB_G = V_PB_G[b][0]; // cross-mobilizer angular vel in G
        const Vec3& v_PB_G = V_PB_G[b][1]; // cross-mobilizer linear  vel in G

        // MISSING-TERM FIX: the bias above is only complete when the mobilized
        // (M) frame origin coincides with the body (B) origin. Simbody's
        // calcParentToChildVelocityJacobianInGroundDot carries an extra HDot_MB_F
        // contribution whenever r_MB = X_MB.p() != 0 (true for every non-root
        // Torsion/torsion body, where r_MB is the bond length). With the per-joint
        // H_FM constant in F (Torsion/Free/Translation/Weld all have HDot_FM = 0),
        // that contribution reduces to the classic centripetal term
        //     R_GF * ( w_FM x (w_FM x r_MB_F) ),   r_MB_F = R_FM * r_MB,
        // and belongs in the LINEAR row of A_mob. Omitting it leaves the joint
        // bias acceleration inconsistent with H/Hdot and slowly pumps energy.
        // NOTE r_MB = X_MB.p() = (~X_BM).p(), the Mo->Bo vector in M (see the
        // matching fix in realizePosition); X_BM.p() would be wrong here too.
        const Vec3& w_FM = V_FM[b][0];                      // cross-mobilizer ang. vel in F
        const Vec3 r_MB_F = X_FM[b].R() * (~m.X_BM[b]).p(); // Mo->Bo, expressed in F
        const Rotation R_GF = X_GB[p].R() * m.X_PF[b].R();  // F orientation in Ground
        const Vec3 centripetal_G = R_GF * (w_FM % (w_FM % r_MB_F));

        // For the seven joints with H_FM constant in F, dH_FM/dt = 0 and the
        // mobilizer bias is exactly the line above (this branch is skipped),
        // so those joints stay BIT-IDENTICAL to the prior engine. The three
        // q-dependent joints (BendStretch/SphericalCoords/FreeLine) add the
        // intrinsic dH_FM/dt*u term that the constant-H derivation dropped:
        //   extraAng_F = sum_j Hdot_ang(j) u_j
        //   extraLin_F = sum_j (Hdot_lin(j) - r_MB_F x Hdot_ang(j)) u_j
        // (The -rdot_MB_F x H_ang(j) part of the full d/dt(H_MB_F) is already
        //  captured by centripetal_G, so it is NOT re-added here.)
        Vec3 extraAng_F(0, 0, 0);
        Vec3 extraLin_F(0, 0, 0);

        if (!RobotModel::jointHasConstantHFM(jt)) {
            std::array<SpatialVec, 5> Hdot; // max dof = 5 (FreeLine)
            jointHDot_FM(jt, X_FM[b], V_FM[b], Hdot.data());
            for (int j = 0; j < dof; ++j) {
                const Real uj = u[uOff + j];
                extraAng_F += Hdot[j][0] * uj;
                extraLin_F += (Hdot[j][1] - (r_MB_F % Hdot[j][0])) * uj;
            }
        }

        const SpatialVec A_mob(w_GP % w_PB_G + R_GF * extraAng_F,
                               w_GP % (v_GB - v_GP) + w_GP % v_PB_G + centripetal_G + R_GF * extraLin_F);
        a_mob[b] = A_mob;

        // Total coriolis accel a = ~Phi * a_parent + A_mob.
        a_tot[b] = (~Phi[b]) * a_tot[p] + A_mob;
    }
}

void RobotEngine::calcQDot(const RobotModel& m, const RobotState& s, Real* qdotOut) {
    const Real* u = s.u();
    const Transform* X_FM = s.X_FM();
    const Real* q = s.q();
    for (int b = 1; b < m.numBodies; ++b) {
        const int qOff = m.bodyQIndex[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        robo::jointQDot(m.bodyJoint[b], q, qOff, u, uOff, dof, X_FM[b], qdotOut);
    }
}

// ---- calcQDotDot: qddot = N qddot-coupling. Scalar joints: udot; quaternion ----
//      bodies (Ball/FreeLine/Free): quaternion second derivative. ----
void RobotEngine::calcQDotDot(const RobotModel& m, RobotState& s) {
    const Real* u = s.u();
    const Real* udot = s.udot();
    const Real* q = s.q();
    const Transform* X_FM = s.X_FM();
    Real* qdd = s.qdotdot();
    for (int b = 1; b < m.numBodies; ++b) {
        const int qOff = m.bodyQIndex[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const JointType jt = m.bodyJoint[b];
        robo::jointQDotDot(jt, q, qOff, u, udot, uOff, dof, X_FM[b], qdd);
    }
}

// ============================================================================
//  TRANSFER: internal bodies -> per-atom Cartesian (accept path)
// ============================================================================
void RobotEngine::fillAtomPositionsFromBodies(const RobotModel& m, RobotState& s) {
    const Transform* X_GB = s.X_GB();
    Vec3* posG = s.atomPosG();
    // Atoms that are Cartesian-integrated inside the proposal (solvent-relaxing
    // NCMC) own their Ground positions directly in posG -- they are NOT placed by
    // any rigid body, so the body->atom fill must leave them untouched (otherwise
    // each step would snap them back onto the welded body and undo the relaxation).
    // mask == nullptr in the welded engine, so this is a no-op there.
    const char* cartMask = s.cartSolventMask();
    const int* __restrict atomBody = m.atomBody.data();
    const Vec3* __restrict station = m.atomStation_B.data();
    Vec3* __restrict posGr = posG;
    if (cartMask == nullptr) {
        // Welded engine (the common case): no Cartesian-solvent atoms to preserve, so
        // the whole loop vectorizes (see realizePosition for the same pattern).
#pragma omp simd
        for (int a = 0; a < m.numAtoms; ++a) {
            const int b = atomBody[a];
            posGr[a] = X_GB[b].p() + X_GB[b].R() * station[a];
        }
    } else {
        for (int a = 0; a < m.numAtoms; ++a) {
            if (cartMask[a]) {
                continue; // Cartesian-solvent atom owns posG directly; leave untouched
            }
            const int b = atomBody[a];
            posGr[a] = X_GB[b].p() + X_GB[b].R() * station[a];
        }
    }
}

// ============================================================================
//  Quaternion renormalization (the no-constraint Verlet "projection").
// ============================================================================
void RobotEngine::normalizeQuaternions(const RobotModel& m, RobotState& s) {
    Real* q = s.q();
    for (int qs : m.quaternionQStart) {
        const Real n2 = q[qs] * q[qs] + q[qs + 1] * q[qs + 1] + q[qs + 2] * q[qs + 2] + q[qs + 3] * q[qs + 3];
        if (n2 <= Real(1e-24)) {
            // Degenerate (near-zero) quaternion: reset to identity instead of
            // dividing by ~0. Should not happen now that q is seeded to identity,
            // but this guarantees no 1/sqrt(0) -> NaN can ever recur.
            q[qs] = Real(1);
            q[qs + 1] = q[qs + 2] = q[qs + 3] = Real(0);
            continue;
        }
        const Real inv = Real(1) / std::sqrt(n2);
        q[qs] *= inv;
        q[qs + 1] *= inv;
        q[qs + 2] *= inv;
        q[qs + 3] *= inv;
    }
}
