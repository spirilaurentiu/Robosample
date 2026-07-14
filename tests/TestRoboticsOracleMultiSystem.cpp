// ============================================================================
//  TestRoboticsOracleMultiSystem.cpp -- the live-Simbody differential oracle,
//  multi-system + body-output family (docs/specs/robotics-oracle-
//  differential.md), split out of TestRoboticsOracle.cpp (TEST-005). One
//  `TEST` per multi-body structural fixture case (mixed chain, forest,
//  wide-star hub, zero-DOF mid-chain, duplicate molecules, applied-force
//  discriminators), run through the shared `runOracleMultiCase` staged
//  comparison (tests/support/RoboticsOracleRunners.hpp). Body-output fields
//  (`expectedBodyOutputFields`, RoboticsOracleLoader.hpp) travel with this
//  family, not as a separate binary -- checkNoSilentGap's `kind=="multi"`
//  branch is where they are consumed.
// ============================================================================
#include <gtest/gtest.h>

#include "support/RoboticsOracleRunners.hpp"

using rtest::kFixtureDir;
using rtest::runOracleMultiCase;

// ---- Phase 1b: multi-body structural cases (§8) ----

TEST(RoboticsOracle, MixedChain) {
    runOracleMultiCase(robotics_oracle_loader::loadMultiCase(kFixtureDir, "MixedChain"));
}

TEST(RoboticsOracle, Forest) {
    runOracleMultiCase(robotics_oracle_loader::loadMultiCase(kFixtureDir, "Forest"));
}

// §8.2 #4: wide star hub (>=8 children) -- summation-order determinism.
TEST(RoboticsOracle, WideStarHub) {
    runOracleMultiCase(robotics_oracle_loader::loadMultiCase(kFixtureDir, "WideStarHub"));
}

// §8.2 #6: zero-DOF Rigid (Weld) body mid-chain.
TEST(RoboticsOracle, RigidMidChain) {
    runOracleMultiCase(robotics_oracle_loader::loadMultiCase(kFixtureDir, "RigidMidChain"));
}

// §8.2 #9: duplicate molecules in a forest -- qIndex/uIndex offset aliasing.
TEST(RoboticsOracle, DuplicateMolecules) {
    runOracleMultiCase(robotics_oracle_loader::loadMultiCase(kFixtureDir, "DuplicateMolecules"));
}

// §8.2 #7: applied-force discriminators (pure-torque, pure-force,
// force-on-Ground).
TEST(RoboticsOracle, ForceDiscriminators) {
    runOracleMultiCase(robotics_oracle_loader::loadMultiCase(kFixtureDir, "ForceDiscriminators"));
}
