// ============================================================================
//  TestRoboticsOracleSingleState.cpp -- the live-Simbody differential oracle,
//  single-state family (docs/specs/robotics-oracle-differential.md), split
//  out of TestRoboticsOracle.cpp (TEST-005). One `TEST` per (JointType,
//  fixture case), run through the shared `runOracleCase` staged comparison
//  (tests/support/RoboticsOracleRunners.hpp) -- see that header for the
//  staged-comparison contract (§6) and tolerance table.
//
//  Also carries the three PORT-ONLY invariants that build a RobotModel
//  directly (no fixture, no Simbody comparison) rather than driving
//  runOracleCase/Multi/Aggregate/Fuzz: PortOnlySingularHinge (§8.2 #2/§6.1),
//  PortOnlySingularHingeNonTorsionThrows (the STEP 3 fail-loud gate), and
//  PortOnlyGroundOnlySystem (§8.2 #10). None of the four case-family runners
//  fit these (they never load a fixture), and their single-body-or-empty
//  construction style is closest to this family, so they live here rather
//  than in a fifth binary.
// ============================================================================
#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>

#include "support/RoboticsOracleRunners.hpp"

using rtest::kFixtureDir;
using rtest::kStage4Tol;
using rtest::runOracleCase;

TEST(RoboticsOracle, Torsion) {
    runOracleCase(JointType::Torsion, robotics_oracle_loader::loadCase(kFixtureDir, "Torsion"));
}

TEST(RoboticsOracle, Free) {
    runOracleCase(JointType::Free, robotics_oracle_loader::loadCase(kFixtureDir, "Free"));
}

TEST(RoboticsOracle, Rigid) {
    runOracleCase(JointType::Rigid, robotics_oracle_loader::loadCase(kFixtureDir, "Rigid"));
}

TEST(RoboticsOracle, Slider) {
    runOracleCase(JointType::Slider, robotics_oracle_loader::loadCase(kFixtureDir, "Slider"));
}

TEST(RoboticsOracle, Cylinder) {
    runOracleCase(JointType::Cylinder, robotics_oracle_loader::loadCase(kFixtureDir, "Cylinder"));
}

TEST(RoboticsOracle, Cartesian) {
    runOracleCase(JointType::Cartesian, robotics_oracle_loader::loadCase(kFixtureDir, "Cartesian"));
}

TEST(RoboticsOracle, Ball) {
    runOracleCase(JointType::Ball, robotics_oracle_loader::loadCase(kFixtureDir, "Ball"));
}

// §8.2 #5: extreme anisotropic inertia (needle) -- pairs with the
// conditioning-stress case but stays well-conditioned enough (cond(D)~1e4)
// to be a normal element-wise structural case, unlike the aggregate-only
// stress case.
TEST(RoboticsOracle, BallNeedle) {
    runOracleCase(JointType::Ball, robotics_oracle_loader::loadCase(kFixtureDir, "BallNeedle"));
}

// §6.1 priority target: q-dependent H_FM, native MobilizedBody::BendStretch.
TEST(RoboticsOracle, BendStretch) {
    runOracleCase(JointType::BendStretch, robotics_oracle_loader::loadCase(kFixtureDir, "BendStretch"));
}

// §6.1 priority target: q-dependent H_FM, native MobilizedBody::SphericalCoords
// (default ctor: radialAxis=Z, no negation, zero offsets).
TEST(RoboticsOracle, SphericalCoords) {
    runOracleCase(JointType::SphericalCoords, robotics_oracle_loader::loadCase(kFixtureDir, "SphericalCoords"));
}

// §6.1 priority target: native MobilizedBody::FreeLine; §8.2 #3 -- qdot raw
// comparison is skipped inside runOracleCase for this joint.
TEST(RoboticsOracle, FreeLine) {
    runOracleCase(JointType::FreeLine, robotics_oracle_loader::loadCase(kFixtureDir, "FreeLine"));
}
// ---------------------------------------------------------------------------
//  §8.2 #2 / §6.1: PORT-ONLY invariant -- a genuinely sub-1e-12 hinge
//  inertia. Un-buildable as a differential: Simbody either rejects a
//  near-singular (non-SPD) hinge inertia at realize or reports ~1/lambda
//  (~1e12+), a step-discontinuity divergence that is by-design (§2b) and
//  must never be diffed against the port's null-space lock. Documented here
//  as port-only per the spec's explicit instruction, not silently dropped.
// ---------------------------------------------------------------------------
TEST(RoboticsOracle, PortOnlySingularHinge) {
    rtest::Rng rng(0xC0FFEEULL);
    BodySpec spec;
    spec.parent = 0;
    spec.joint = JointType::Torsion;
    spec.X_PF = Transform(rng.rotation(), rng.vec3());
    spec.X_BM = Transform(rng.rotation(), rng.vec3());
    // Absolute rotational inertia about the hinge axis ~ mass*unitInertia
    // ~ 1e-13*0.5 = 5e-14, well below the 1e-12 lock threshold.
    spec.mass = Real(1e-13);
    spec.com_B = Vec3(0);
    spec.inertia_B = UnitInertia(Real(0.5), Real(0.5), Real(0.5));

    RobotModel m = buildForest({spec});
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = Real(0.3);
    s.u()[0] = Real(0.5);
    s.bodyForceG()[0] = SpatialVec(Vec3(0), Vec3(0));
    s.bodyForceG()[1] = SpatialVec(Vec3(0, 0, 0.02), Vec3(0.01, -0.01, 0.02));
    s.mobilityForce()[0] = Real(0.03);

    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    RobotEngine::calcUDot(m, s);

    // The null-space lock (src/RobotEngine.cpp invertDense, 1-dof closed
    // form) zeroes DI EXACTLY when |D| <= 1e-12 -- not a smooth clamp.
    EXPECT_EQ(s.DI()[m.bodyUSqIndex[1]], Real(0)) << "expected the null-space lock to zero DI exactly";

    // Finite everywhere: no NaN/Inf leaked out of the lock branch.
    EXPECT_TRUE(std::isfinite(s.eps()[0]));
    EXPECT_TRUE(std::isfinite(s.udot()[0]));
    EXPECT_TRUE(std::isfinite(s.Z()[1].angular.norm()));
    EXPECT_TRUE(std::isfinite(s.Z()[1].linear.norm()));
    EXPECT_TRUE(std::isfinite(s.A_GB()[1].angular.norm()));
    EXPECT_TRUE(std::isfinite(s.A_GB()[1].linear.norm()));

    // Zero on the locked direction: DI==0 => udot = eps*DI = 0 regardless of
    // the applied force/torque, i.e. the DOF is frozen at its initial value.
    EXPECT_NEAR(s.udot()[0], Real(0), kStage4Tol) << "locked (only) direction must carry zero udot";
}

// ---------------------------------------------------------------------------
//  STEP 3 fail-loud gate (docs/specs/singular-dof-fixman.md, Review outcome):
//  invertDense null-locking a direction is safe (silently) ONLY for a leaf
//  Torsion body (D_b provably a run-constant, see the code comment at the
//  gate site). A Slider is NOT that shape -- same recipe as
//  PortOnlySingularHinge (near-zero mass, so its 1-dof translational D also
//  locks below 1e-12) but a different joint type -- so realizeArticulatedBodyInertias
//  MUST throw instead of silently locking-and-continuing. This is the
//  positive counterpart to PortOnlySingularHinge (which proves the gate does
//  NOT false-positive on the recognized shape); together they pin the gate's
//  discriminating boundary.
// ---------------------------------------------------------------------------
TEST(RoboticsOracle, PortOnlySingularHingeNonTorsionThrows) {
    rtest::Rng rng(0xC0FFEEULL);
    BodySpec spec;
    spec.parent = 0;
    spec.joint = JointType::Slider;
    spec.X_PF = Transform(rng.rotation(), rng.vec3());
    spec.X_BM = Transform(rng.rotation(), rng.vec3());
    spec.mass = Real(1e-13); // same near-zero mass as PortOnlySingularHinge -> D also locks
    spec.com_B = Vec3(0);
    spec.inertia_B = UnitInertia(Real(0.5), Real(0.5), Real(0.5));

    RobotModel m = buildForest({spec});
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = Real(0.1);
    s.u()[0] = Real(0.0);

    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeVelocity(m, s);
    EXPECT_THROW(RobotEngine::realizeArticulatedBodyInertias(m, s), std::runtime_error)
        << "a locked Slider (not a leaf Torsion) must fail loud, not silently lock-and-continue";
}

// ---------------------------------------------------------------------------
//  §8.2 #10: Ground-only system (nq=nu=0 total) -- guards the comparator and
//  the corrector (den+1e-30) against empty-array UB. PORT-ONLY: both sides
//  trivially agree at 0 (no dynamics to differ on), so the value this case
//  adds is the port's own robustness to n=0, not an external pin -- a
//  dedicated Simbody empty-system build would add generator complexity for
//  no discriminating power (deliberate scope call, Rule 2).
// ---------------------------------------------------------------------------
TEST(RoboticsOracle, PortOnlyGroundOnlySystem) {
    RobotModel m = buildForest({});
    ASSERT_EQ(m.numBodies, 1) << "Ground only";
    ASSERT_EQ(m.nq, 0);
    ASSERT_EQ(m.nu, 0);

    RobotState s;
    s.allocateFull(m);

    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    RobotEngine::calcUDot(m, s);

    EXPECT_EQ(RobotEngine::calcLogDetM(m, s), Real(0));
    EXPECT_EQ(RobotEngine::calcKineticEnergy(m, s), Real(0));
}
