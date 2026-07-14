// ============================================================================
//  RobotEngine_dynamics.cpp - articulated-body inertia factorization and
//  forward dynamics (forces -> udot).
//
//  One of four cohesive splits of the former RobotEngine.cpp (SPLIT-R3, pure
//  code motion); see RobotEngine_kinematics.cpp's header comment for the
//  Simbody provenance / validation notes shared by all four
//  RobotEngine_*.cpp translation units.
// ============================================================================

#include "RobotEngine.hpp"

#include <stdexcept>
#include <string>

#include "RobotEngine_internal.hpp"
#include "math/hinge_linalg.hpp"
#include "math/robo_debug.hpp"
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

using robo::detail::eigDecompAndTol;
using robo::detail::invertDense;
using robo::detail::spatialDot;

// ============================================================================
//  ARTICULATED-BODY INERTIAS (inward)   RigidBodyNodeSpec.cpp
// ============================================================================
// Position-only articulated-inertia factorization: P, PPlus, D, DI, G (the per-body
// Jacobi eigensolve / invertDense). Pure function of q (masses + geometry) -- NO velocity
// dependence -- so the verlet corrector hoists it to once/step (see the wrapper below and
// docs/specs/gpu-cartesian-kinematics/03-aba-parallelization.md Sec.0.5). The velocity-
// dependent centrifugal seed is split into seedArticulatedCentrifugal.
void RobotEngine::factorizeArticulatedInertias(const RobotModel& m, RobotState& s) {
    ArticulatedInertia* P = s.P();
    ArticulatedInertia* PPlus = s.PPlus();
    const SpatialVec* H = s.H();
    SpatialVec* G = s.G();
    Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    const SpatialInertia* Mk = s.Mk_G();

    for (int b = m.numBodies - 1; b >= 1; --b) {
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];

        // P = Mk_G + sum_children Phi * PPlus_child * ~Phi.
        ArticulatedInertia Pb(Mk[b]);
        for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
            const int c = m.bodyChildren[ci];
            Pb += PPlus[c].shift(Phi[c].l());
        }
        P[b] = Pb;

        if (dof == 0) { // Weld: nothing felt through the (zero) mobility space
            PPlus[b] = Pb;
            continue;
        }

        // PH columns, D = ~H*P*H, DI = D^-1, G = P*H*DI.
        SpatialVec PH[6];
        for (int j = 0; j < dof; ++j) {
            PH[j] = Pb * H[uOff + j];
        }
        Real D[36];
        for (int i = 0; i < dof; ++i) {
            for (int j = 0; j < dof; ++j) {
                D[i * dof + j] = spatialDot(H[uOff + i], PH[j]);
            }
        }
        Real* DI = &DIpool[m.bodyUSqIndex[b]];
        const int numLocked = invertDense(D, dof, DI);

        // STEP 3 fail-loud gate (docs/specs/singular-dof-fixman.md, Review
        // outcome): invertDense null-locked >=1 direction of this body's
        // hinge inertia. That is EXPECTED and safe only for a body whose
        // D_b is PROVABLY a run-constant (locked at this q => locked for
        // EVERY q, so freezing it is bias-free, C3) -- never merely
        // "small at this reference config" (B2's rejected false-positive
        // risk: an angle-flexible multi-dof joint, e.g. BendStretch/
        // SphericalCoords/Cartesian/FreeLine, can be genuinely collinear
        // only at isolated q and must NOT be treated as structural).
        //
        // A LEAF (no children) Torsion (1-dof) body qualifies unconditionally:
        // D_b = ~H_b P_b H_b with H_b (Ground) = R_GF * H_FM, H_FM CONSTANT
        // for Torsion (jointHasConstantHFM) and R_GF depending only on
        // ancestors (upstream of this body's own q); P_b == Mk_G[b] exactly
        // (leaf: no PPlus_child term). Writing Mk_G[b] = R_GB Mk_B R_GB^T
        // with R_GB = R_GF R_FM(q_own) R_MB, the R_GF factors cancel by
        // orthogonality in D_b = H_FM^T [R_FM(q_own) R_MB Mk_B R_MB^T
        // R_FM(q_own)^T] H_FM, which no longer contains R_GF at all -- D_b
        // is therefore independent of every ancestor config AND of the
        // body's own q (a revolute joint's own-axis inertia is invariant to
        // rotation about that same axis, the standard rigid-body fact).
        // Confirmed empirically to ~1e-14 by the Cyclic1APQPhantomLogDetIsRunConstant
        // LEMMA test (tests/TestRoboticsOracleMolecule.cpp) on 1APQ's three
        // leaf single-atom on-axis Torsion phantoms. STEP 4 (a build-time
        // weld removing these bodies outright, making this gate unreachable
        // for them) is DEFERRED -- see the coder checkpoint; this runtime
        // recognition is the provably-safe substitute while it is deferred.
        // Any OTHER locked body (non-leaf, non-Torsion, or multi-dof) is NOT
        // this shape and throws instead of silently locking (Rule 11).
        if (numLocked > 0
            && !(m.bodyChildrenBeg[b] == m.bodyChildrenEnd[b] && m.bodyJoint[b] == JointType::Torsion)) {
            Real eig[6];
            Real eigV[36];
            Real eigTol;
            eigDecompAndTol(D, dof, eig, eigV, eigTol);
            Real minEig = std::abs(eig[0]);
            for (int k = 1; k < dof; ++k) {
                minEig = std::min(minEig, std::abs(eig[k]));
            }
            std::string msg = "realizeArticulatedBodyInertias: invertDense null-locked "
                              + std::to_string(numLocked) + " direction(s) of D_b on body "
                              + std::to_string(b)
                              + " (JointType=" + std::to_string(static_cast<int>(m.bodyJoint[b]))
                              + ", dof=" + std::to_string(dof) + ", atoms=[";
            // bodyAtomsBeg/End are only populated for a fully-built (real-molecule)
            // RobotModel; a hand-built synthetic test model (tests/RobotBuilders.hpp
            // buildForest without attachAtoms) leaves them empty -- guard the index
            // so a genuine gate failure never masks itself behind an out-of-bounds
            // read while formatting the diagnostic.
            if (static_cast<std::size_t>(b) < m.bodyAtomsBeg.size()) {
                for (int a = m.bodyAtomsBeg[b]; a < m.bodyAtomsEnd[b]; ++a) {
                    if (a != m.bodyAtomsBeg[b]) {
                        msg += ",";
                    }
                    msg += std::to_string(m.bodyAtoms[a]);
                }
            } else {
                msg += "unavailable: synthetic model with no atom map";
            }
            msg += "], min-eig(D)=" + std::to_string(minEig)
                   + ") that is not a recognized structural phantom (leaf Torsion) -- "
                     "see docs/specs/singular-dof-fixman.md STEP 3";
            throw std::runtime_error(msg);
        }

#if ROBO_DEBUG
        {
            bool bad = false;
            for (int k = 0; k < dof * dof; ++k) {
                if (!robodbg::fin(D[k]) || !robodbg::fin(DI[k])) {
                    bad = true;
                    break;
                }
            }
            // Also flag a wildly large inverse (near-singular D), which is the
            // precursor to a blow-up even before it becomes a literal NaN.
            Real maxDI = 0;
            for (int k = 0; k < dof * dof; ++k) {
                maxDI = std::max(maxDI, std::abs(DI[k]));
            }
            if ((bad || maxDI > Real(1e8)) && !robodbg::firstNanDumped) {
                std::cout << "\n--- realizeABI: suspect hinge inverse at body " << b << " joint "
                          << static_cast<int>(m.bodyJoint[b]) << " dof " << dof << " (max|DI|=" << maxDI
                          << ", eval#" << robodbg::evalCount << " step#" << robodbg::stepCount << ")\n";
                std::cout << "  Mk: mass=" << Mk[b].getMass() << " com=" << Mk[b].getMassCenter() << "\n";
                std::cout << "  D =\n";
                for (int i = 0; i < dof; ++i) {
                    std::cout << "    ";
                    for (int j = 0; j < dof; ++j) {
                        std::cout << std::setw(13) << D[i * dof + j] << ' ';
                    }
                    std::cout << "\n";
                }
                std::cout << "  DI =\n";
                for (int i = 0; i < dof; ++i) {
                    std::cout << "    ";
                    for (int j = 0; j < dof; ++j) {
                        std::cout << std::setw(13) << DI[i * dof + j] << ' ';
                    }
                    std::cout << "\n";
                }
                std::cout << "  H cols (Ground):\n";
                for (int j = 0; j < dof; ++j) {
                    std::cout << "    H[" << j << "] = " << H[uOff + j] << "\n";
                }
                std::cout << std::flush;
            }
        }
#endif

        for (int j = 0; j < dof; ++j) {
            SpatialVec gj(Vec3(0), Vec3(0));
            for (int k = 0; k < dof; ++k) {
                gj += PH[k] * DI[k * dof + j];
            }
            G[uOff + j] = gj;
        }

        // PPlus = P - G*~PH, reconstructed in ArticulatedInertia block form
        // exactly as Simbody does (sum of outer products over columns).
        Mat33 massMoment(0);
        Mat33 mass(0);
        Mat33 inertia(0);
        for (int j = 0; j < dof; ++j) {
            const Vec3& Ga = G[uOff + j][0];
            const Vec3& Gl = G[uOff + j][1];
            const Vec3& PHa = PH[j][0];
            const Vec3& PHl = PH[j][1];
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    massMoment(r, c) += Ga[r] * PHl[c];
                    mass(r, c) += Gl[r] * PHl[c];
                    inertia(r, c) += Ga[r] * PHa[c];
                }
            }
        }
        const SymMat33 symMass(mass(0, 0),
                               (mass(1, 0) + mass(0, 1)) / 2,
                               mass(1, 1),
                               (mass(2, 0) + mass(0, 2)) / 2,
                               (mass(2, 1) + mass(1, 2)) / 2,
                               mass(2, 2));
        const SymMat33 symInertia(inertia(0, 0),
                                  (inertia(1, 0) + inertia(0, 1)) / 2,
                                  inertia(1, 1),
                                  (inertia(2, 0) + inertia(0, 2)) / 2,
                                  (inertia(2, 1) + inertia(1, 2)) / 2,
                                  inertia(2, 2));
        PPlus[b] = Pb - ArticulatedInertia(symMass, massMoment, symInertia);
    }
}

// Velocity-dependent articulated centrifugal seed: abcf = P * a_mob + gyro (seed for
// calcUDot pass1). This is the ONLY part of the old realizeArticulatedBodyInertias that
// depends on u; splitting it out lets the corrector iterate just this + calcUDot while the
// (expensive, position-only) factorization above runs once/step. PRECONDITION: realizeVelocity
// current (a_mob/gyro) and factorizeArticulatedInertias current (P).
void RobotEngine::seedArticulatedCentrifugal(const RobotModel& m, RobotState& s) {
    // CRITICAL: Simbody (RigidBodyNode.cpp, realizeArticulatedBodyVelocityCache) forms this
    // from the MOBILIZER coriolis acceleration (the per-joint incremental term A), NOT the
    // TOTAL coriolis acceleration (a = ~Phi*a_parent + A). Using the total here adds a
    // spurious, velocity^2-scaled centrifugal force on every non-root body; it propagates
    // inward through the pass-1 Phi*zPlus sum, corrupts udot, and feeds back through Verlet
    // as monotonic energy injection -> blow-up.
    const ArticulatedInertia* P = s.P();
    const SpatialVec* a_mob = s.mobCoriolisA();
    const SpatialVec* gyro = s.gyro();
    SpatialVec* abcf = s.abCentrifugal();
    for (int b = 1; b < m.numBodies; ++b) {
        abcf[b] = P[b] * a_mob[b] + gyro[b];
    }
}

// Backward-compatible full pass = factorization + centrifugal seed, in the original order.
// Callers other than the verlet corrector (reinitialize, captureReactionSnapshot, ...) use
// this and are unchanged; the corrector (RobotIntegrator.hpp) calls the two halves
// separately, hoisting the factorization out of the u-only iteration.
void RobotEngine::realizeArticulatedBodyInertias(const RobotModel& m, RobotState& s) {
    factorizeArticulatedInertias(m, s);
    seedArticulatedCentrifugal(m, s);
}

// ============================================================================
//  FORWARD DYNAMICS  (forces -> udot)   RigidBodyNodeSpec.cpp
// ============================================================================
void RobotEngine::calcUDot(const RobotModel& m, RobotState& s) {
    const SpatialVec* H = s.H();
    const SpatialVec* G = s.G();
    const ArticulatedInertia* P = s.P();
    const Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    const SpatialVec* abcf = s.abCentrifugal();
    const SpatialVec* a_mob = s.mobCoriolisA();
    const SpatialVec* bodyF = s.bodyForceG();
    const Real* jointF = s.mobilityForce();

    SpatialVec* Z = s.Z();
    SpatialVec* ZPlus = s.zPlus();
    Real* eps = s.eps();
    SpatialVec* A_GB = s.A_GB();
    Real* udot = s.udot();

    // Pass 1 inward: z = (P a + b) - F + sum Phi*zPlus_child ; eps = f - ~H z ;
    // zPlus = z + G eps.
    for (int b = m.numBodies - 1; b >= 1; --b) {
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        SpatialVec z = abcf[b] - bodyF[b];
        for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
            const int c = m.bodyChildren[ci];
            z += Phi[c] * ZPlus[c];
        }
        for (int j = 0; j < dof; ++j) {
            eps[uOff + j] = jointF[uOff + j] - spatialDot(H[uOff + j], z);
        }
        SpatialVec zp = z;
        for (int j = 0; j < dof; ++j) {
            zp += G[uOff + j] * eps[uOff + j];
        }
        Z[b] = z;
        ZPlus[b] = zp;
    }

    // Pass 2 outward: APlus = ~Phi A_parent; udot = DI eps - ~G APlus;
    // A_GB = APlus + H udot + a_mob.
    A_GB[0] = SpatialVec(Vec3(0), Vec3(0));
    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const SpatialVec APlus = (~Phi[b]) * A_GB[p];
        const Real* DI = &DIpool[m.bodyUSqIndex[b]];
        SpatialVec acc = APlus;
        for (int i = 0; i < dof; ++i) {
            Real ud = 0;
            for (int j = 0; j < dof; ++j) {
                ud += DI[i * dof + j] * eps[uOff + j];
            }
            ud -= spatialDot(G[uOff + i], APlus); // -~G*APlus, row i
            udot[uOff + i] = ud;
            acc += H[uOff + i] * ud;
        }
        A_GB[b] = acc + a_mob[b];
    }
}
