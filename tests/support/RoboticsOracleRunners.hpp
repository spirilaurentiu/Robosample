#pragma once
// ============================================================================
//  RoboticsOracleRunners.hpp -- the four staged-comparison runner functions
//  (runOracleCase/runOracleMultiCase/runOracleAggregateCase/runOracleFuzzCase)
//  and their shared helpers, promoted from TestRoboticsOracle.cpp (TEST-005)
//  so every case-family split binary
//  (TestRoboticsOracleSingleState/MultiSystem/Aggregate/Fuzz.cpp) shares ONE
//  copy. Per RoboticsOracleLoader.hpp's own banner: "this file only changes
//  where the struct's data comes from, never what is compared" -- this
//  promotion is pure code motion (linkage location only, `inline` instead of
//  the origin file's anonymous namespace); no runner body, tolerance, or
//  comparison changed. docs/specs/robotics-oracle-differential.md is the
//  authoritative spec these runners implement (staged comparison, §6;
//  tolerances, §6/§8).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>

#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"

#include "../RobotBuilders.hpp"
#include "../RobotLinearAlgebra.hpp"
#include "../TestHelpers.hpp"

#include "../RoboticsOracleLoader.hpp"
#include "../fixtures/robotics_oracle/RoboticsOracleTypes.hpp"

using namespace robo;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::NearMat33;
using rtest::NearVec3;

namespace rtest {

#ifndef ROBOTICS_ORACLE_FIXTURE_DIR
#error "ROBOTICS_ORACLE_FIXTURE_DIR must be defined by the build (tests/fixtures/robotics_oracle)"
#endif
// CMake-defined (test-scoped, TestRoboticsOracle target only): the on-disk
// location of the <Case>.npz + <Case>.manifest.json pairs (§9).
const std::string kFixtureDir = ROBOTICS_ORACLE_FIXTURE_DIR;

// ---------------------------------------------------------------------------
//  §6 stage tolerances (the spec's own table, not the generic TestHelpers
//  tiers -- reusing the tiers here would either be looser than the spec
//  mandates (kLoose=1e-9 for a 1e-10 stage) or not map cleanly, so these are
//  named locally and cite the spec section they encode).
// ---------------------------------------------------------------------------
constexpr Real kStage1Tol = 1e-10; // X_GB
constexpr Real kFrameGateTol = 1e-9; // §4.2 X_FM frame-equality gate
constexpr Real kStage2Tol = 1e-9;  // V_GB, qdot
constexpr Real kStage3Tol = 1e-8;  // P, PPlus, DI, G (element-wise, gated)
constexpr Real kLockGateTol = 1e-9; // §6.1 min-eig(D) cross-engine agreement
constexpr Real kStage4Tol = 1e-8;  // Z, zPlus, eps, udot, A_GB
constexpr Real kStage5Tol = 1e-8;  // dense M, logDetM, reactions
constexpr Real kStage5EigTol = 1e-6; // eig(M)

// §2b/§6.1: the port's Jacobi null-space lock threshold (src/RobotEngine.cpp,
// tests/RobotLinearAlgebra.hpp::invertDense). Element-wise DI/PPlus/G diffs
// are only valid strictly above this.
constexpr Real kLockTol = 1e-12;

// §8.1 mechanism 2 (J1-B): the conditioning-stress case's report-body A_GB is
// computed FROM the ill-conditioned D (cond(D) ~ 1e6-1e8 by design), so the
// port-vs-Simbody hinge-inverse method difference (Jacobi vs factorization)
// amplifies by ~cond(D) into this end-product, same mechanism as the DI
// element error the spec says "rides at the 1e-8 boundary" for this case.
// Observed residual here is ~1.5e-7 (cond(D)~9e7 * eps~2e-16 * O(100) scale
// ~2e-6 theoretical ceiling) -- this tolerance is 100x looser than the
// standard kStage4Tol (still 10x tighter than the theoretical ceiling, so it
// stays discriminating against a real transcription bug) and applies ONLY to
// this one aggregate case's report-body acceleration, never to the
// well-conditioned structural cases above.
constexpr Real kConditioningAccTol = 1e-6;

// §8.3 fuzz batch: per-body udot/A_GB end-product tolerance. Looser than
// kStage4Tol (1e-8) because the fuzz batch's topologies (unlike the
// hand-tuned structural cases) are NOT authored to stay well-conditioned --
// the only generation-time guarantee is the §6.1 mandatory filter rejecting
// states at/below 10*lockTol (the discrete lock discontinuity), which does
// NOT bound the smooth-regime condition number (§8.1 mechanism 2). The
// ConditioningStress case already established that a deliberately
// ill-conditioned (cond(D)~1e8) end-product can carry ~1.5e-7 residual
// against the 1e-8 base (kConditioningAccTol=1e-6, 100x looser); the fuzz
// batch samples a WIDER swath of the depth x cond plane (§8.1), so its
// tolerance is set one further order of magnitude looser -- still 4+ orders
// tighter than a transcription-bug-sized (O(1) relative) divergence.
constexpr Real kFuzzEndProductTol = 1e-5;

// §8.3/§9 logDetM, ConditioningStress and fuzz-batch sites only: logDetM sums
// Sigma ln|D_k| over every body with no cancellation (unlike udot/A_GB), so a
// large-|logDetM| system (e.g. FuzzTopo_06's worst state, |logDetM|~28.9)
// accumulates conditioning-limited FP noise across O(numBodies) terms the
// same way kConditioningAccTol's report-body A_GB does across cond(D) --
// adjudicated benign (2.766e-7 abs = 9.57e-9 rel, i.e. it clears relative
// kStage5Tol=1e-8 by only ~4.3%), not a port bug (calcLogDetM confirmed
// term-for-term correct vs Simbody's calcDetMPass2Outward). One order looser
// than kStage5Tol gives headroom for the fuzz batch's unbounded cond x depth
// sampling while staying 6+ orders below an O(1)-relative transcription bug.
constexpr Real kLogDetMFuzzTol = 1e-7;

// ---------------------------------------------------------------------------
//  Fixture (row-major flat double[]) -> port value type conversions.
// ---------------------------------------------------------------------------
inline Transform transformFromFixture(const double rot[9], const double p[3]) {
    Rotation r(Mat33(rot[0], rot[1], rot[2], rot[3], rot[4], rot[5], rot[6], rot[7], rot[8]));
    return Transform(r, Vec3(p[0], p[1], p[2]));
}
inline UnitInertia unitInertiaFromFixture(const double packed[6]) {
    return UnitInertia(packed[0], packed[1], packed[2], packed[3], packed[4], packed[5]);
}
inline SpatialVec spatialVecFromFixture(const double ang[3], const double lin[3]) {
    return SpatialVec(Vec3(ang[0], ang[1], ang[2]), Vec3(lin[0], lin[1], lin[2]));
}

inline ::testing::AssertionResult NearScalar(Real a, Real b, Real tol, const char* what) {
    const Real d = std::abs(a - b);
    if (d <= tol) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure() << what << " differ by " << d << " (tol " << tol << "): " << a
                                         << " vs " << b;
}

// Compare a port ArticulatedInertia (angAng=J packed 6, angLin=F full 9,
// linLin=M packed 6) against the fixture's packed blocks (same packing,
// RoboticsOracleTypes.hpp banner).
inline void expectArticulatedInertiaNear(const char* label,
                                  const ArticulatedInertia& port,
                                  const double jFix[6],
                                  const double fFix[9],
                                  const double mFix[6],
                                  Real tol) {
    for (int i = 0; i < 6; ++i) {
        EXPECT_TRUE(NearScalar(port.angAng.elems[static_cast<std::size_t>(i)], jFix[i], tol, label))
            << label << " J[" << i << "]";
    }
    for (int i = 0; i < 9; ++i) {
        EXPECT_TRUE(NearScalar(port.angLin.elems[static_cast<std::size_t>(i)], fFix[i], tol, label))
            << label << " F[" << i << "]";
    }
    for (int i = 0; i < 6; ++i) {
        EXPECT_TRUE(NearScalar(port.linLin.elems[static_cast<std::size_t>(i)], mFix[i], tol, label))
            << label << " M[" << i << "]";
    }
}

// Independent forward-Jacobian dense M (port's OWN V_GB/Mk_G, but computed by
// a totally different path than the ABI recursion -- see TestMassMatrix.cpp's
// buildDenseM, duplicated here verbatim; this is the "correctness note" the
// spec §3 calls out: the oracle diff is what makes THIS independence externally
// meaningful, since Simbody's calcM never touches the port's V_GB/Mk_G at all).
inline std::vector<Real> buildDenseMPort(const RobotModel& m, RobotState& s) {
    const int n = m.nu;
    const int bodies = m.numBodies;
    const SpatialInertia* mk = s.Mk_G();

    std::vector<std::vector<SpatialVec>> col(static_cast<std::size_t>(n),
                                             std::vector<SpatialVec>(static_cast<std::size_t>(bodies)));
    std::vector<Real> uSave(s.u(), s.u() + n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            s.u()[i] = (i == j) ? Real(1) : Real(0);
        }
        RobotEngine::realizeVelocity(m, s);
        const SpatialVec* v = s.V_GB();
        for (int b = 0; b < bodies; ++b) {
            col[static_cast<std::size_t>(j)][static_cast<std::size_t>(b)] = v[b];
        }
    }
    std::copy(uSave.begin(), uSave.end(), s.u());
    RobotEngine::realizeVelocity(m, s); // restore V_GB for the real u before returning

    std::vector<Real> dense(static_cast<std::size_t>(n * n), 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Real acc = 0;
            for (int b = 1; b < bodies; ++b) {
                const SpatialVec mv = mk[b] * col[static_cast<std::size_t>(j)][static_cast<std::size_t>(b)];
                acc += dot(col[static_cast<std::size_t>(i)][static_cast<std::size_t>(b)], mv);
            }
            dense[static_cast<std::size_t>((i * n) + j)] = acc;
        }
    }
    return dense;
}

// ---------------------------------------------------------------------------
//  The staged comparison, run once per (JointType, fixture case). `body` is
//  always 1 (single-body robot: Ground=0, the one body=1).
// ---------------------------------------------------------------------------
inline void runOracleCase(JointType jt, const robotics_oracle::OracleCase& fix) {
    ASSERT_GT(fix.numStates, 0) << "empty fixture -- regenerate";

    BodySpec spec;
    spec.parent = 0;
    spec.joint = jt;
    spec.X_PF = transformFromFixture(fix.X_PF_R, fix.X_PF_p);
    spec.X_BM = transformFromFixture(fix.X_BM_R, fix.X_BM_p);
    spec.mass = fix.mass;
    spec.com_B = Vec3(fix.com_B[0], fix.com_B[1], fix.com_B[2]);
    spec.inertia_B = unitInertiaFromFixture(fix.unitInertia_B);

    RobotModel m = buildForest({spec});
    const int body = 1;

    // ---- §5 Scope-A correspondence: DOF layout must agree before any
    //      numeric comparison, or every later mismatch is uninterpretable. ----
    ASSERT_EQ(m.bodyNQ[body], fix.nq) << "nq mismatch vs Simbody getNumQ";
    ASSERT_EQ(m.bodyNU[body], fix.nu) << "nu mismatch vs Simbody getNumU";
    if (RobotModel::jointUsesQuaternion(jt)) {
        ASSERT_FALSE(m.quaternionQStart.empty());
        EXPECT_EQ(m.quaternionQStart.front(), m.bodyQIndex[body])
            << "quaternion block must be the first 4 q (§5 Scope A)";
    }

    RobotState s;
    s.allocateFull(m);

    for (int st = 0; st < fix.numStates; ++st) {
        const robotics_oracle::OracleState& in = fix.states[st];
        SCOPED_TRACE(::testing::Message() << "state=" << in.label << " (" << st << ")");
        ASSERT_EQ(in.nq, m.nq);
        ASSERT_EQ(in.nu, m.nu);

        for (int i = 0; i < in.nq; ++i) {
            s.q()[i] = in.q[i];
        }
        for (int i = 0; i < in.nu; ++i) {
            s.u()[i] = in.u[i];
        }
        s.bodyForceG()[0] = SpatialVec(Vec3(0), Vec3(0));
        s.bodyForceG()[body] = spatialVecFromFixture(in.bodyForceTorque, in.bodyForceForce);
        for (int i = 0; i < in.nu; ++i) {
            s.mobilityForce()[i] = in.mobilityForce[i];
        }

        // ---- Stage 1: kinematics ----
        RobotEngine::realizePosition(m, s);
        EXPECT_TRUE(NearMat33(s.X_GB()[body].R(),
                              Mat33(in.X_GB_R[0], in.X_GB_R[1], in.X_GB_R[2], in.X_GB_R[3], in.X_GB_R[4],
                                    in.X_GB_R[5], in.X_GB_R[6], in.X_GB_R[7], in.X_GB_R[8]),
                              kStage1Tol))
            << "X_GB.R";
        EXPECT_TRUE(NearVec3(s.X_GB()[body].p(), Vec3(in.X_GB_p[0], in.X_GB_p[1], in.X_GB_p[2]), kStage1Tol))
            << "X_GB.p";

        // §4.2 frame-equality gate: assert BEFORE trusting any frame-dependent
        // (stage 3/4) comparison that the port picked the same F/M convention
        // Simbody's native mobilizer uses.
        const bool frameOk =
            NearMat33(s.X_FM()[body].R(),
                     Mat33(in.X_FM_R[0], in.X_FM_R[1], in.X_FM_R[2], in.X_FM_R[3], in.X_FM_R[4], in.X_FM_R[5],
                           in.X_FM_R[6], in.X_FM_R[7], in.X_FM_R[8]),
                     kFrameGateTol)
            && NearVec3(s.X_FM()[body].p(), Vec3(in.X_FM_p[0], in.X_FM_p[1], in.X_FM_p[2]), kFrameGateTol);
        ASSERT_TRUE(frameOk) << "X_FM frame-equality gate failed -- port and Simbody disagree on the "
                                "joint's F/M convention; no frame-dependent comparison below is valid";

        // ---- Stage 2: velocity ----
        RobotEngine::realizeVelocity(m, s);
        EXPECT_TRUE(NearVec3(s.V_GB()[body].angular, Vec3(in.V_GB_ang[0], in.V_GB_ang[1], in.V_GB_ang[2]),
                             kStage2Tol))
            << "V_GB.angular";
        EXPECT_TRUE(
            NearVec3(s.V_GB()[body].linear, Vec3(in.V_GB_lin[0], in.V_GB_lin[1], in.V_GB_lin[2]), kStage2Tol))
            << "V_GB.linear";
        // FreeLine (§8.2 #3): the u->qdot N-map is rank-deficient (nq=7,
        // nu=5, spin about the body's own line suppressed) and WHICH
        // component absorbs the suppression is convention-dependent -- never
        // diff the raw qdot quaternion block for this joint. Every other
        // joint's qdot is externally re-derivable (kinematic N-map), so it is
        // compared normally.
        std::vector<Real> qdot(static_cast<std::size_t>(m.nq));
        RobotEngine::calcQDot(m, s, qdot.data());
        if (jt != JointType::FreeLine) {
            for (int i = 0; i < in.nq; ++i) {
                EXPECT_TRUE(NearScalar(qdot[static_cast<std::size_t>(i)], in.qdot[i], kStage2Tol, "qdot"))
                    << "qdot[" << i << "]";
            }
        }

        // ---- Stage 3: articulated body inertia (frame-dependent) ----
        RobotEngine::realizeArticulatedBodyInertias(m, s);

        // §8.2 #6: a 0-dof (Rigid/Weld) body has no D/DI/G/lock-gate to
        // compare -- P/PPlus (still populated, reflecting the transmitted
        // subtree inertia) ARE still compared. jacobiSymEig(n=0, ...) would
        // read an uninitialized eigenvalue slot, so this is a real guard, not
        // stylistic.
        const int dof = m.bodyNU[body];
        if (dof > 0) {
            // §6.1 lock gate: derive min-eig(D)_port from the port's own DI
            // (the port never stores D directly -- only DI persists in
            // RobotState). In the smooth regime (away from the lock) DI ==
            // D^-1 exactly, so eig(D) = 1/eig(DI); this is the ONLY way to
            // obtain a min-eig(D)_port without a disasm src/include change
            // (HARD CONSTRAINT).
            const Real* diPort = &s.DI()[m.bodyUSqIndex[body]];
            Real eigD[6], eigV[36];
            robo_linalg::jacobiSymEig(diPort, dof, eigD, eigV);
            Real maxEigDI = eigD[0];
            for (int k = 1; k < dof; ++k) {
                maxEigDI = std::max(maxEigDI, eigD[k]);
            }
            ASSERT_GT(maxEigDI, Real(0)) << "DI is not positive -- port hinge inertia is degenerate";
            const Real minEigDPort = Real(1) / maxEigDI;

            EXPECT_TRUE(NearScalar(minEigDPort, in.minEigD,
                                   kLockGateTol * std::max(Real(1), std::abs(in.minEigD)), "min-eig(D)"))
                << "min-eig(D) cross-engine agreement (§6.1 step 1)";
            ASSERT_GT(minEigDPort, kLockTol) << "port min-eig(D) at/below the null-space lock -- element-wise "
                                                "DI/PPlus/G diff is by-design invalid here (§6.1 step 2/§2b)";
            ASSERT_GT(in.minEigD, kLockTol) << "Simbody min-eig(D) at/below the lock";
        }

        // P/PPlus are always compared (well-defined even at dof==0);
        // DI/G element-wise only when there is a hinge block to compare.
        expectArticulatedInertiaNear("P", s.P()[body], in.P_J, in.P_F, in.P_M, kStage3Tol);
        expectArticulatedInertiaNear("PPlus", s.PPlus()[body], in.PPlus_J, in.PPlus_F, in.PPlus_M, kStage3Tol);
        if (dof > 0) {
            const Real* diPort = &s.DI()[m.bodyUSqIndex[body]];
            for (int i = 0; i < dof; ++i) {
                for (int j = 0; j < dof; ++j) {
                    EXPECT_TRUE(NearScalar(diPort[(i * dof) + j], in.DI[(i * dof) + j], kStage3Tol, "DI"))
                        << "DI[" << i << "," << j << "]";
                }
            }
            const SpatialVec* gPort = &s.G()[m.bodyUIndex[body]];
            for (int j = 0; j < dof; ++j) {
                EXPECT_TRUE(NearVec3(gPort[j].angular, Vec3(in.G_ang[j][0], in.G_ang[j][1], in.G_ang[j][2]),
                                     kStage3Tol))
                    << "G[" << j << "].angular";
                EXPECT_TRUE(
                    NearVec3(gPort[j].linear, Vec3(in.G_lin[j][0], in.G_lin[j][1], in.G_lin[j][2]), kStage3Tol))
                    << "G[" << j << "].linear";
            }
        }

        // ---- Stage 4: acceleration pass ----
        RobotEngine::calcUDot(m, s);
        EXPECT_TRUE(NearVec3(s.Z()[body].angular, Vec3(in.Z_ang[0], in.Z_ang[1], in.Z_ang[2]), kStage4Tol))
            << "Z.angular";
        EXPECT_TRUE(NearVec3(s.Z()[body].linear, Vec3(in.Z_lin[0], in.Z_lin[1], in.Z_lin[2]), kStage4Tol))
            << "Z.linear";
        EXPECT_TRUE(
            NearVec3(s.zPlus()[body].angular, Vec3(in.ZPlus_ang[0], in.ZPlus_ang[1], in.ZPlus_ang[2]),
                     kStage4Tol))
            << "zPlus.angular";
        EXPECT_TRUE(
            NearVec3(s.zPlus()[body].linear, Vec3(in.ZPlus_lin[0], in.ZPlus_lin[1], in.ZPlus_lin[2]),
                     kStage4Tol))
            << "zPlus.linear";
        const Real* epsPort = &s.eps()[m.bodyUIndex[body]];
        const Real* udotPort = &s.udot()[m.bodyUIndex[body]];
        for (int i = 0; i < dof; ++i) {
            EXPECT_TRUE(NearScalar(epsPort[i], in.eps[i], kStage4Tol, "eps")) << "eps[" << i << "]";
            EXPECT_TRUE(NearScalar(udotPort[i], in.udot[i], kStage4Tol, "udot")) << "udot[" << i << "]";
        }
        EXPECT_TRUE(NearVec3(s.A_GB()[body].angular, Vec3(in.A_GB_ang[0], in.A_GB_ang[1], in.A_GB_ang[2]),
                             kStage4Tol))
            << "A_GB.angular";
        EXPECT_TRUE(NearVec3(s.A_GB()[body].linear, Vec3(in.A_GB_lin[0], in.A_GB_lin[1], in.A_GB_lin[2]),
                             kStage4Tol))
            << "A_GB.linear";

        // ---- Stage 5: dense M / logDetM / reactions ----
        // M is system-wide (nu x nu); for these single-body robots nu == dof,
        // but use m.nu explicitly since that is the quantity buildDenseMPort
        // actually sizes (this test never has a multi-body system, so the two
        // coincide, but the loop bound should say what it means).
        const int nu = m.nu;
        const std::vector<Real> denseM = buildDenseMPort(m, s);
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nu; ++j) {
                EXPECT_TRUE(
                    NearScalar(denseM[static_cast<std::size_t>((i * nu) + j)], in.Mdense[(i * nu) + j],
                              kStage5Tol, "M"))
                    << "M[" << i << "," << j << "]";
            }
        }
        Real dPort[6], dFix[6], vPort[36], vFix[36];
        robo_linalg::jacobiSymEig(denseM.data(), nu, dPort, vPort);
        robo_linalg::jacobiSymEig(in.Mdense, nu, dFix, vFix);
        std::sort(dPort, dPort + nu);
        std::sort(dFix, dFix + nu);
        for (int k = 0; k < nu; ++k) {
            EXPECT_TRUE(NearScalar(dPort[k], dFix[k], kStage5EigTol, "eig(M)")) << "eig(M)[" << k << "]";
        }

        EXPECT_TRUE(NearScalar(RobotEngine::calcLogDetM(m, s), in.logDetM,
                               kStage5Tol * std::max(Real(1), std::abs(in.logDetM)), "logDetM"));

        std::vector<SpatialVec> reactBo(static_cast<std::size_t>(m.numBodies));
        std::vector<SpatialVec> reactMo(static_cast<std::size_t>(m.numBodies));
        RobotEngine::calcMobilizerReactionForces(m, s, reactBo.data(), reactMo.data());
        EXPECT_TRUE(NearVec3(reactBo[body].angular,
                             Vec3(in.reactionBoAng[0], in.reactionBoAng[1], in.reactionBoAng[2]), kStage5Tol))
            << "reaction@Bo.angular";
        EXPECT_TRUE(NearVec3(reactBo[body].linear,
                             Vec3(in.reactionBoLin[0], in.reactionBoLin[1], in.reactionBoLin[2]), kStage5Tol))
            << "reaction@Bo.linear";
        EXPECT_TRUE(NearVec3(reactMo[body].angular,
                             Vec3(in.reactionMoAng[0], in.reactionMoAng[1], in.reactionMoAng[2]), kStage5Tol))
            << "reaction@Mo.angular";
        EXPECT_TRUE(NearVec3(reactMo[body].linear,
                             Vec3(in.reactionMoLin[0], in.reactionMoLin[1], in.reactionMoLin[2]), kStage5Tol))
            << "reaction@Mo.linear";
    }
}

// ---------------------------------------------------------------------------
//  Phase 1b: multi-body structural cases (mixed chain, forest, wide-star hub,
//  zero-DOF Rigid mid-chain, duplicate molecules, applied-force
//  discriminators) and the two §8.1 stress cases (aggregate invariants only).
// ---------------------------------------------------------------------------

// Dynamic-size cyclic-Jacobi symmetric eigensolver, local to THIS test file
// only. tests/RobotLinearAlgebra.hpp::jacobiSymEig is hard-capped at n<=6
// (fixed Real a[36] -- correct and load-bearing there: it is the per-body
// hinge-block D solver and dof<=6 always) and is shared by 6+ other test
// TUs, so it is NOT extended here (Rule 3: surgical, don't touch widely-used
// shared infra for a need specific to this file's system-wide dense-M eig
// check, where nu can exceed 6 for a multi-body case).
inline void jacobiSymEigDyn(const std::vector<Real>& Ain, int n, std::vector<Real>& dOut) {
    std::vector<Real> a = Ain;
    for (int sweep = 0; sweep < 100; ++sweep) {
        Real off = 0;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                off += a[static_cast<std::size_t>((p * n) + q)] * a[static_cast<std::size_t>((p * n) + q)];
            }
        }
        if (off <= Real(1e-30)) {
            break;
        }
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const Real apq = a[static_cast<std::size_t>((p * n) + q)];
                if (apq == Real(0)) {
                    continue;
                }
                const Real app = a[static_cast<std::size_t>((p * n) + p)];
                const Real aqq = a[static_cast<std::size_t>((q * n) + q)];
                const Real theta = (aqq - app) / (2 * apq);
                const Real t =
                    (theta >= 0 ? Real(1) : Real(-1)) / (std::abs(theta) + std::sqrt((theta * theta) + 1));
                const Real cs = Real(1) / std::sqrt((t * t) + 1);
                const Real sn = t * cs;
                for (int i = 0; i < n; ++i) {
                    const Real aip = a[static_cast<std::size_t>((i * n) + p)];
                    const Real aiq = a[static_cast<std::size_t>((i * n) + q)];
                    a[static_cast<std::size_t>((i * n) + p)] = (cs * aip) - (sn * aiq);
                    a[static_cast<std::size_t>((i * n) + q)] = (sn * aip) + (cs * aiq);
                }
                for (int i = 0; i < n; ++i) {
                    const Real api = a[static_cast<std::size_t>((p * n) + i)];
                    const Real aqi = a[static_cast<std::size_t>((q * n) + i)];
                    a[static_cast<std::size_t>((p * n) + i)] = (cs * api) - (sn * aqi);
                    a[static_cast<std::size_t>((q * n) + i)] = (sn * api) + (cs * aqi);
                }
            }
        }
    }
    dOut.assign(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        dOut[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>((i * n) + i)];
    }
}

// Build a rtest::BodySpec list from a fixture's per-body model spec (§4.1/§5
// Scope-A correspondence: same body order, same jointType, so q/u align by
// construction).
inline std::vector<BodySpec> buildSpecsFromFixture(const robotics_oracle::BodyModelSpec* bodies, int numBodies) {
    std::vector<BodySpec> specs(static_cast<std::size_t>(numBodies));
    for (int b = 0; b < numBodies; ++b) {
        const robotics_oracle::BodyModelSpec& bm = bodies[b];
        BodySpec s;
        s.parent = bm.parent;
        s.joint = static_cast<JointType>(bm.joint);
        s.X_PF = transformFromFixture(bm.X_PF_R, bm.X_PF_p);
        s.X_BM = transformFromFixture(bm.X_BM_R, bm.X_BM_p);
        s.mass = bm.mass;
        s.com_B = Vec3(bm.com_B[0], bm.com_B[1], bm.com_B[2]);
        s.inertia_B = unitInertiaFromFixture(bm.unitInertia_B);
        specs[static_cast<std::size_t>(b)] = s;
    }
    return specs;
}

inline void runOracleMultiCase(const robotics_oracle::OracleMultiCase& fix) {
    ASSERT_GT(fix.numStates, 0) << "empty fixture -- regenerate";
    ASSERT_GT(fix.numBodies, 0);

    const std::vector<BodySpec> specs = buildSpecsFromFixture(fix.bodies, fix.numBodies);
    RobotModel m = buildForest(specs);

    // §5 Scope-A correspondence: total DOF layout must agree before any
    // numeric comparison (per-body ground truth already covered by the 10
    // single-joint tests; a wrong qIndex/uIndex OFFSET combination here would
    // still show up as a stage-1 X_GB mismatch on a later body).
    ASSERT_EQ(m.nq, fix.states[0].nq) << "total nq mismatch vs Simbody";
    ASSERT_EQ(m.nu, fix.states[0].nu) << "total nu mismatch vs Simbody";

    RobotState s;
    s.allocateFull(m);

    for (int st = 0; st < fix.numStates; ++st) {
        const robotics_oracle::OracleMultiState& in = fix.states[st];
        SCOPED_TRACE(::testing::Message() << "state=" << in.label << " (" << st << ")");
        ASSERT_EQ(in.nq, m.nq);
        ASSERT_EQ(in.nu, m.nu);

        for (int i = 0; i < in.nq; ++i) {
            s.q()[i] = in.q[i];
        }
        for (int i = 0; i < in.nu; ++i) {
            s.u()[i] = in.u[i];
        }
        s.bodyForceG()[0] = spatialVecFromFixture(in.groundForceTorque, in.groundForceForce);
        for (int b = 0; b < fix.numBodies; ++b) {
            s.bodyForceG()[b + 1] = spatialVecFromFixture(in.bodyForceTorque[b], in.bodyForceForce[b]);
        }
        for (int i = 0; i < in.nu; ++i) {
            s.mobilityForce()[i] = in.mobilityForce[i];
        }

        // ---- Stage 1 + §4.2 frame-equality gate, per body ----
        RobotEngine::realizePosition(m, s);
        for (int b = 0; b < fix.numBodies; ++b) {
            const int body = b + 1;
            const robotics_oracle::OracleBodyOutput& bo = in.body[b];
            SCOPED_TRACE(::testing::Message() << "body=" << body);
            EXPECT_TRUE(NearMat33(s.X_GB()[body].R(),
                                  Mat33(bo.X_GB_R[0], bo.X_GB_R[1], bo.X_GB_R[2], bo.X_GB_R[3], bo.X_GB_R[4],
                                        bo.X_GB_R[5], bo.X_GB_R[6], bo.X_GB_R[7], bo.X_GB_R[8]),
                                  kStage1Tol))
                << "X_GB.R";
            EXPECT_TRUE(
                NearVec3(s.X_GB()[body].p(), Vec3(bo.X_GB_p[0], bo.X_GB_p[1], bo.X_GB_p[2]), kStage1Tol))
                << "X_GB.p";
            const bool frameOk =
                NearMat33(s.X_FM()[body].R(),
                         Mat33(bo.X_FM_R[0], bo.X_FM_R[1], bo.X_FM_R[2], bo.X_FM_R[3], bo.X_FM_R[4],
                               bo.X_FM_R[5], bo.X_FM_R[6], bo.X_FM_R[7], bo.X_FM_R[8]),
                         kFrameGateTol)
                && NearVec3(s.X_FM()[body].p(), Vec3(bo.X_FM_p[0], bo.X_FM_p[1], bo.X_FM_p[2]), kFrameGateTol);
            ASSERT_TRUE(frameOk) << "X_FM frame-equality gate failed for body " << body;
        }

        // ---- Stage 2 ----
        RobotEngine::realizeVelocity(m, s);
        std::vector<Real> qdot(static_cast<std::size_t>(m.nq));
        RobotEngine::calcQDot(m, s, qdot.data());
        for (int b = 0; b < fix.numBodies; ++b) {
            const int body = b + 1;
            const robotics_oracle::OracleBodyOutput& bo = in.body[b];
            SCOPED_TRACE(::testing::Message() << "body=" << body);
            EXPECT_TRUE(NearVec3(s.V_GB()[body].angular, Vec3(bo.V_GB_ang[0], bo.V_GB_ang[1], bo.V_GB_ang[2]),
                                 kStage2Tol))
                << "V_GB.angular";
            EXPECT_TRUE(NearVec3(s.V_GB()[body].linear, Vec3(bo.V_GB_lin[0], bo.V_GB_lin[1], bo.V_GB_lin[2]),
                                 kStage2Tol))
                << "V_GB.linear";
            // §8.2 #3: skip the raw qdot compare for FreeLine.
            if (m.bodyJoint[body] != JointType::FreeLine) {
                const int qOff = m.bodyQIndex[body];
                const int nq = m.bodyNQ[body];
                for (int i = 0; i < nq; ++i) {
                    EXPECT_TRUE(NearScalar(qdot[static_cast<std::size_t>(qOff + i)], bo.qdot[i], kStage2Tol,
                                           "qdot"))
                        << "qdot[" << i << "]";
                }
            }
        }

        // ---- Stage 3 ----
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        for (int b = 0; b < fix.numBodies; ++b) {
            const int body = b + 1;
            const robotics_oracle::OracleBodyOutput& bo = in.body[b];
            SCOPED_TRACE(::testing::Message() << "body=" << body);
            const int dof = m.bodyNU[body];
            if (dof > 0) {
                const Real* diPort = &s.DI()[m.bodyUSqIndex[body]];
                Real eigD[6], eigV[36];
                robo_linalg::jacobiSymEig(diPort, dof, eigD, eigV);
                Real maxEigDI = eigD[0];
                for (int k = 1; k < dof; ++k) {
                    maxEigDI = std::max(maxEigDI, eigD[k]);
                }
                ASSERT_GT(maxEigDI, Real(0)) << "DI is not positive -- port hinge inertia is degenerate";
                const Real minEigDPort = Real(1) / maxEigDI;
                EXPECT_TRUE(NearScalar(minEigDPort, bo.minEigD,
                                       kLockGateTol * std::max(Real(1), std::abs(bo.minEigD)), "min-eig(D)"))
                    << "min-eig(D) cross-engine agreement (§6.1 step 1)";
                ASSERT_GT(minEigDPort, kLockTol) << "port min-eig(D) at/below the null-space lock";
                ASSERT_GT(bo.minEigD, kLockTol) << "Simbody min-eig(D) at/below the lock";
            }
            expectArticulatedInertiaNear("P", s.P()[body], bo.P_J, bo.P_F, bo.P_M, kStage3Tol);
            expectArticulatedInertiaNear("PPlus", s.PPlus()[body], bo.PPlus_J, bo.PPlus_F, bo.PPlus_M, kStage3Tol);
            if (dof > 0) {
                const Real* diPort = &s.DI()[m.bodyUSqIndex[body]];
                for (int i = 0; i < dof; ++i) {
                    for (int j = 0; j < dof; ++j) {
                        EXPECT_TRUE(
                            NearScalar(diPort[(i * dof) + j], bo.DI[(i * dof) + j], kStage3Tol, "DI"))
                            << "DI[" << i << "," << j << "]";
                    }
                }
                const SpatialVec* gPort = &s.G()[m.bodyUIndex[body]];
                for (int j = 0; j < dof; ++j) {
                    EXPECT_TRUE(NearVec3(gPort[j].angular, Vec3(bo.G_ang[j][0], bo.G_ang[j][1], bo.G_ang[j][2]),
                                         kStage3Tol))
                        << "G[" << j << "].angular";
                    EXPECT_TRUE(
                        NearVec3(gPort[j].linear, Vec3(bo.G_lin[j][0], bo.G_lin[j][1], bo.G_lin[j][2]),
                                 kStage3Tol))
                        << "G[" << j << "].linear";
                }
            }
        }

        // ---- Stage 4 ----
        RobotEngine::calcUDot(m, s);
        for (int b = 0; b < fix.numBodies; ++b) {
            const int body = b + 1;
            const robotics_oracle::OracleBodyOutput& bo = in.body[b];
            SCOPED_TRACE(::testing::Message() << "body=" << body);
            EXPECT_TRUE(NearVec3(s.Z()[body].angular, Vec3(bo.Z_ang[0], bo.Z_ang[1], bo.Z_ang[2]), kStage4Tol))
                << "Z.angular";
            EXPECT_TRUE(NearVec3(s.Z()[body].linear, Vec3(bo.Z_lin[0], bo.Z_lin[1], bo.Z_lin[2]), kStage4Tol))
                << "Z.linear";
            EXPECT_TRUE(NearVec3(s.zPlus()[body].angular, Vec3(bo.ZPlus_ang[0], bo.ZPlus_ang[1], bo.ZPlus_ang[2]),
                                 kStage4Tol))
                << "zPlus.angular";
            EXPECT_TRUE(NearVec3(s.zPlus()[body].linear, Vec3(bo.ZPlus_lin[0], bo.ZPlus_lin[1], bo.ZPlus_lin[2]),
                                 kStage4Tol))
                << "zPlus.linear";
            const int dof = m.bodyNU[body];
            const Real* epsPort = &s.eps()[m.bodyUIndex[body]];
            const Real* udotPort = &s.udot()[m.bodyUIndex[body]];
            for (int i = 0; i < dof; ++i) {
                EXPECT_TRUE(NearScalar(epsPort[i], bo.eps[i], kStage4Tol, "eps")) << "eps[" << i << "]";
                EXPECT_TRUE(NearScalar(udotPort[i], bo.udot[i], kStage4Tol, "udot")) << "udot[" << i << "]";
            }
            EXPECT_TRUE(
                NearVec3(s.A_GB()[body].angular, Vec3(bo.A_GB_ang[0], bo.A_GB_ang[1], bo.A_GB_ang[2]), kStage4Tol))
                << "A_GB.angular";
            EXPECT_TRUE(
                NearVec3(s.A_GB()[body].linear, Vec3(bo.A_GB_lin[0], bo.A_GB_lin[1], bo.A_GB_lin[2]), kStage4Tol))
                << "A_GB.linear";
        }

        // §8.2 #7: force on Ground must not leak into udot -- an explicit
        // invariant on top of the fixture diff above (which would only show
        // it as an incidental zero-vs-zero match).
        if (std::string(in.label) == "force-on-ground") {
            for (int i = 0; i < m.nu; ++i) {
                EXPECT_NEAR(s.udot()[i], Real(0), kStage4Tol) << "force on Ground leaked into udot[" << i << "]";
            }
        }

        // ---- Stage 5: dense M / logDetM / reactions ----
        const int nu = m.nu;
        const std::vector<Real> denseM = buildDenseMPort(m, s);
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nu; ++j) {
                EXPECT_TRUE(NearScalar(denseM[static_cast<std::size_t>((i * nu) + j)], in.Mdense[(i * nu) + j],
                                       kStage5Tol, "M"))
                    << "M[" << i << "," << j << "]";
            }
        }
        std::vector<Real> dPort, dFix;
        jacobiSymEigDyn(denseM, nu, dPort);
        jacobiSymEigDyn(std::vector<Real>(in.Mdense, in.Mdense + (nu * nu)), nu, dFix);
        std::sort(dPort.begin(), dPort.end());
        std::sort(dFix.begin(), dFix.end());
        for (int k = 0; k < nu; ++k) {
            EXPECT_TRUE(NearScalar(dPort[static_cast<std::size_t>(k)], dFix[static_cast<std::size_t>(k)],
                                   kStage5EigTol, "eig(M)"))
                << "eig(M)[" << k << "]";
        }

        EXPECT_TRUE(NearScalar(RobotEngine::calcLogDetM(m, s), in.logDetM,
                               kStage5Tol * std::max(Real(1), std::abs(in.logDetM)), "logDetM"));

        std::vector<SpatialVec> reactBo(static_cast<std::size_t>(m.numBodies));
        std::vector<SpatialVec> reactMo(static_cast<std::size_t>(m.numBodies));
        RobotEngine::calcMobilizerReactionForces(m, s, reactBo.data(), reactMo.data());
        for (int b = 0; b < fix.numBodies; ++b) {
            const int body = b + 1;
            const robotics_oracle::OracleBodyOutput& bo = in.body[b];
            SCOPED_TRACE(::testing::Message() << "body=" << body);
            EXPECT_TRUE(NearVec3(reactBo[body].angular, Vec3(bo.reactionBoAng[0], bo.reactionBoAng[1],
                                                             bo.reactionBoAng[2]), kStage5Tol))
                << "reaction@Bo.angular";
            EXPECT_TRUE(NearVec3(reactBo[body].linear, Vec3(bo.reactionBoLin[0], bo.reactionBoLin[1],
                                                            bo.reactionBoLin[2]), kStage5Tol))
                << "reaction@Bo.linear";
            EXPECT_TRUE(NearVec3(reactMo[body].angular, Vec3(bo.reactionMoAng[0], bo.reactionMoAng[1],
                                                             bo.reactionMoAng[2]), kStage5Tol))
                << "reaction@Mo.angular";
            EXPECT_TRUE(NearVec3(reactMo[body].linear, Vec3(bo.reactionMoLin[0], bo.reactionMoLin[1],
                                                            bo.reactionMoLin[2]), kStage5Tol))
                << "reaction@Mo.linear";
        }
    }
}

// ---------------------------------------------------------------------------
//  §8.1 stress cases: AGGREGATE invariants + end-products only (§6.1/§9,
//  J1-B) -- never element-wise DI/PPlus/G here, even though both engines stay
//  above the lock, because per-element DI error rides at the stage-3
//  tolerance boundary for the conditioning case (spec's own words).
// ---------------------------------------------------------------------------
inline void runOracleAggregateCase(const robotics_oracle::OracleAggregateCase& fix) {
    ASSERT_GT(fix.numStates, 0) << "empty fixture -- regenerate";
    ASSERT_GT(fix.numBodies, 0);

    const std::vector<BodySpec> specs = buildSpecsFromFixture(fix.bodies, fix.numBodies);
    RobotModel m = buildForest(specs);

    RobotState s;
    s.allocateFull(m);

    for (int st = 0; st < fix.numStates; ++st) {
        const robotics_oracle::OracleAggregateState& in = fix.states[st];
        SCOPED_TRACE(::testing::Message() << "state=" << in.label << " (" << st << ")");
        ASSERT_EQ(in.nq, m.nq) << "total nq mismatch vs Simbody";
        ASSERT_EQ(in.nu, m.nu) << "total nu mismatch vs Simbody";

        for (int i = 0; i < in.nq; ++i) {
            s.q()[i] = in.q[i];
        }
        for (int i = 0; i < in.nu; ++i) {
            s.u()[i] = in.u[i];
            s.mobilityForce()[i] = 0;
        }
        s.bodyForceG()[0] = SpatialVec(Vec3(0), Vec3(0));
        for (int b = 0; b < fix.numBodies; ++b) {
            s.bodyForceG()[b + 1] = spatialVecFromFixture(in.bodyForceTorque[b], in.bodyForceForce[b]);
        }

        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        RobotEngine::calcUDot(m, s);

        EXPECT_TRUE(NearScalar(RobotEngine::calcLogDetM(m, s), in.logDetM,
                               kLogDetMFuzzTol * std::max(Real(1), std::abs(in.logDetM)), "logDetM"));

        const Real kePort = RobotEngine::calcKineticEnergy(m, s);
        EXPECT_TRUE(
            NearScalar(kePort, in.totalKE, kStage4Tol * std::max(Real(1), std::abs(in.totalKE)), "KE"));

        Real ss = 0;
        for (int i = 0; i < m.nu; ++i) {
            ss += s.udot()[i] * s.udot()[i];
        }
        const Real normUdotPort = std::sqrt(ss);
        EXPECT_TRUE(NearScalar(normUdotPort, in.normUdot, kStage4Tol * std::max(Real(1), std::abs(in.normUdot)),
                               "||udot||"));

        // Worst (smallest) min-eig(D) over all bodies -- same DI-derivation
        // trick as the structural cases' §6.1 gate, just minimized instead of
        // read at a single body.
        Real worst = Real(1e300);
        for (int body = 1; body < m.numBodies; ++body) {
            const int dof = m.bodyNU[body];
            if (dof == 0) {
                continue;
            }
            const Real* diPort = &s.DI()[m.bodyUSqIndex[body]];
            Real eigD[6], eigV[36];
            robo_linalg::jacobiSymEig(diPort, dof, eigD, eigV);
            Real maxEigDI = eigD[0];
            for (int k = 1; k < dof; ++k) {
                maxEigDI = std::max(maxEigDI, eigD[k]);
            }
            ASSERT_GT(maxEigDI, Real(0)) << "DI is not positive at body " << body;
            worst = std::min(worst, Real(1) / maxEigDI);
        }
        EXPECT_TRUE(NearScalar(worst, in.minEigD, kLockGateTol * std::max(Real(1), std::abs(in.minEigD)),
                               "min-eig(D) (worst)"));
        // §8.1: the stress case is deliberately tuned to stay above the lock
        // on BOTH engines -- assert that design intent held, not just diffed.
        ASSERT_GT(worst, kLockTol) << "port worst min-eig(D) at/below the lock -- stress case must stay smooth";
        ASSERT_GT(in.minEigD, kLockTol) << "Simbody worst min-eig(D) at/below the lock";

        const int rb = in.reportBody;
        ASSERT_GT(rb, 0);
        ASSERT_LT(rb, m.numBodies);
        EXPECT_TRUE(NearMat33(s.X_GB()[rb].R(),
                              Mat33(in.reportX_GB_R[0], in.reportX_GB_R[1], in.reportX_GB_R[2],
                                    in.reportX_GB_R[3], in.reportX_GB_R[4], in.reportX_GB_R[5],
                                    in.reportX_GB_R[6], in.reportX_GB_R[7], in.reportX_GB_R[8]),
                              kStage1Tol))
            << "report body X_GB.R";
        EXPECT_TRUE(NearVec3(s.X_GB()[rb].p(),
                             Vec3(in.reportX_GB_p[0], in.reportX_GB_p[1], in.reportX_GB_p[2]), kStage1Tol))
            << "report body X_GB.p";
        EXPECT_TRUE(NearVec3(s.V_GB()[rb].angular,
                             Vec3(in.reportV_GB_ang[0], in.reportV_GB_ang[1], in.reportV_GB_ang[2]), kStage2Tol))
            << "report body V_GB.angular";
        EXPECT_TRUE(NearVec3(s.V_GB()[rb].linear,
                             Vec3(in.reportV_GB_lin[0], in.reportV_GB_lin[1], in.reportV_GB_lin[2]), kStage2Tol))
            << "report body V_GB.linear";
        // A_GB uses kConditioningAccTol, not kStage4Tol: it is computed FROM
        // the ill-conditioned D by design (§8.1 mechanism 2), so the
        // port-vs-Simbody hinge-inverse method difference amplifies by
        // ~cond(D) here, same mechanism as the DI element error the spec
        // says "rides at the 1e-8 boundary" for this case.
        EXPECT_TRUE(NearVec3(s.A_GB()[rb].angular,
                             Vec3(in.reportA_GB_ang[0], in.reportA_GB_ang[1], in.reportA_GB_ang[2]),
                             kConditioningAccTol))
            << "report body A_GB.angular";
        EXPECT_TRUE(NearVec3(s.A_GB()[rb].linear,
                             Vec3(in.reportA_GB_lin[0], in.reportA_GB_lin[1], in.reportA_GB_lin[2]),
                             kConditioningAccTol))
            << "report body A_GB.linear";
    }
}

// ---------------------------------------------------------------------------
//  §8.3 randomized fuzz batch: seeded, frozen, clone-generated (Architecture
//  B, spec banner) -- this test draws NO randomness itself, it only replays
//  the baked (q,u) states and diffs the SAME §9 aggregate invariants as
//  runOracleAggregateCase PLUS the per-body frame-invariant end-products
//  (udot, A_GB) the fuzz batch adds (§8.3 task 1/4). Never element-wise
//  DI/PPlus/G (same J1-B reasoning as the hand-tuned stress cases, only more
//  so: the fuzz topologies are not authored to stay well-conditioned).
// ---------------------------------------------------------------------------
inline void runOracleFuzzCase(const robotics_oracle::OracleFuzzCase& fix) {
    ASSERT_GT(fix.numStates, 0) << "empty fuzz fixture -- regenerate";
    ASSERT_GT(fix.numBodies, 0);
    ASSERT_GE(fix.rngSeed, 0) << "fuzz fixture is missing its §8.3/§9 rng_seed provenance";
    ASSERT_GE(fix.resampleCount, 0) << "fuzz fixture is missing its §8.3/§9 resample_count provenance";

    const std::vector<BodySpec> specs = buildSpecsFromFixture(fix.bodies, fix.numBodies);
    RobotModel m = buildForest(specs);

    RobotState s;
    s.allocateFull(m);

    for (int st = 0; st < fix.numStates; ++st) {
        const robotics_oracle::OracleFuzzState& in = fix.states[st];
        SCOPED_TRACE(::testing::Message() << "fuzz state=" << in.label << " (" << st << ")");
        ASSERT_EQ(in.nq, m.nq) << "total nq mismatch vs Simbody";
        ASSERT_EQ(in.nu, m.nu) << "total nu mismatch vs Simbody";

        for (int i = 0; i < in.nq; ++i) {
            s.q()[i] = in.q[i];
        }
        for (int i = 0; i < in.nu; ++i) {
            s.u()[i] = in.u[i];
            s.mobilityForce()[i] = 0;
        }
        for (int b = 0; b <= fix.numBodies; ++b) {
            s.bodyForceG()[b] = SpatialVec(Vec3(0), Vec3(0)); // §8.3: fuzz battery applies zero force
        }

        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        RobotEngine::calcUDot(m, s);

        // ---- §9 aggregate invariants (same pattern as runOracleAggregateCase) ----
        EXPECT_TRUE(NearScalar(RobotEngine::calcLogDetM(m, s), in.logDetM,
                               kLogDetMFuzzTol * std::max(Real(1), std::abs(in.logDetM)), "logDetM"));

        const Real kePort = RobotEngine::calcKineticEnergy(m, s);
        EXPECT_TRUE(
            NearScalar(kePort, in.totalKE, kStage4Tol * std::max(Real(1), std::abs(in.totalKE)), "KE"));

        Real ss = 0;
        for (int i = 0; i < m.nu; ++i) {
            ss += s.udot()[i] * s.udot()[i];
        }
        const Real normUdotPort = std::sqrt(ss);
        EXPECT_TRUE(NearScalar(normUdotPort, in.normUdot, kStage4Tol * std::max(Real(1), std::abs(in.normUdot)),
                               "||udot||"));

        Real worst = Real(1e300);
        for (int body = 1; body < m.numBodies; ++body) {
            const int dof = m.bodyNU[body];
            if (dof == 0) {
                continue;
            }
            const Real* diPort = &s.DI()[m.bodyUSqIndex[body]];
            Real eigD[6], eigV[36];
            robo_linalg::jacobiSymEig(diPort, dof, eigD, eigV);
            Real maxEigDI = eigD[0];
            for (int k = 1; k < dof; ++k) {
                maxEigDI = std::max(maxEigDI, eigD[k]);
            }
            ASSERT_GT(maxEigDI, Real(0)) << "DI is not positive at body " << body;
            worst = std::min(worst, Real(1) / maxEigDI);
        }
        EXPECT_TRUE(NearScalar(worst, in.minEigD, kLockGateTol * std::max(Real(1), std::abs(in.minEigD)),
                               "min-eig(D) (worst)"));
        // §8.3 mandatory singularity filter: the generator already rejected
        // any draw whose Simbody-side worst min-eig(D) < 10*lockTol, so a
        // port value at/below kLockTol here is a genuine port-vs-Simbody
        // divergence (the port locking where Simbody did not), never filter
        // slack -- fail loud rather than silently accept it.
        ASSERT_GT(worst, kLockTol) << "port worst min-eig(D) at/below the lock in a filtered fuzz state";
        ASSERT_GT(in.minEigD, kLockTol) << "Simbody worst min-eig(D) at/below the lock in a filtered fuzz state";

        // ---- §8.3 task 1/4: per-body frame-invariant end-products ----
        for (int b = 0; b < fix.numBodies; ++b) {
            const int body = b + 1;
            SCOPED_TRACE(::testing::Message() << "body=" << body);
            const int dof = m.bodyNU[body];
            const int uOff = m.bodyUIndex[body];
            for (int i = 0; i < dof; ++i) {
                const Real refUdot = in.udot[uOff + i];
                EXPECT_TRUE(NearScalar(s.udot()[uOff + i], refUdot,
                                       kFuzzEndProductTol * std::max(Real(1), std::abs(refUdot)), "udot"))
                    << "udot[" << i << "]";
            }
            const Vec3 refAGBAng(in.A_GB_ang[b][0], in.A_GB_ang[b][1], in.A_GB_ang[b][2]);
            const Vec3 refAGBLin(in.A_GB_lin[b][0], in.A_GB_lin[b][1], in.A_GB_lin[b][2]);
            EXPECT_TRUE(NearVec3(s.A_GB()[body].angular, refAGBAng,
                                 kFuzzEndProductTol * std::max(Real(1), refAGBAng.norm())))
                << "A_GB.angular";
            EXPECT_TRUE(NearVec3(s.A_GB()[body].linear, refAGBLin,
                                 kFuzzEndProductTol * std::max(Real(1), refAGBLin.norm())))
                << "A_GB.linear";
        }
    }
}

} // namespace rtest
