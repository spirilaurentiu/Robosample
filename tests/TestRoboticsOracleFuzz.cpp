// ============================================================================
//  TestRoboticsOracleFuzz.cpp -- the live-Simbody differential oracle,
//  randomized fuzz family (docs/specs/robotics-oracle-differential.md §8.3),
//  split out of TestRoboticsOracle.cpp (TEST-005). Two batches of `TEST`s,
//  both run through the shared `runOracleFuzzCase` staged comparison
//  (tests/support/RoboticsOracleRunners.hpp): batch 1 (FuzzStates_*) replays
//  seeded random (q,u) states over each existing Scope-A topology; batch 2
//  (FuzzTopo_*) replays a seeded batch of random topologies sampling the
//  depth x cond plane the two hand-tuned stress cases only sample at two
//  points. Every case is seeded/frozen/clone-generated -- this file draws no
//  randomness itself.
// ============================================================================
#include <gtest/gtest.h>

#include "support/RoboticsOracleRunners.hpp"

using rtest::kFixtureDir;
using rtest::runOracleFuzzCase;

// ---------------------------------------------------------------------------
//  §8.3 randomized fuzz batch, batch 1: N in [32,64] seeded random (q,u)
//  states over each existing Scope-A topology (all 10 single joints +
//  MixedChain + Forest) -- catches state-dependent H_FM(q)/Coriolis/
//  velocity-coupled udot bugs the 3 hand-authored states (rest/random/force)
//  cannot, per §8.3's own rationale.
// ---------------------------------------------------------------------------
TEST(RoboticsOracle, FuzzStatesRigid) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_Rigid"));
}
TEST(RoboticsOracle, FuzzStatesTorsion) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_Torsion"));
}
TEST(RoboticsOracle, FuzzStatesSlider) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_Slider"));
}
TEST(RoboticsOracle, FuzzStatesCylinder) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_Cylinder"));
}
TEST(RoboticsOracle, FuzzStatesBendStretch) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_BendStretch"));
}
TEST(RoboticsOracle, FuzzStatesCartesian) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_Cartesian"));
}
TEST(RoboticsOracle, FuzzStatesBall) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_Ball"));
}
TEST(RoboticsOracle, FuzzStatesSphericalCoords) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_SphericalCoords"));
}
TEST(RoboticsOracle, FuzzStatesFreeLine) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_FreeLine"));
}
TEST(RoboticsOracle, FuzzStatesFree) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_Free"));
}
TEST(RoboticsOracle, FuzzStatesMixedChain) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_MixedChain"));
}
TEST(RoboticsOracle, FuzzStatesForest) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzStates_Forest"));
}

// ---------------------------------------------------------------------------
//  §8.3 randomized fuzz batch, batch 2: a small seeded batch of random
//  topologies (depth in [2,60], random branching/joint-type/mass-ratio,
//  random-but-fixed frames) -- samples the depth x cond plane (§8.1) the two
//  hand-tuned stress cases (DepthChainStress/ConditioningStress) only sample
//  at two points.
// ---------------------------------------------------------------------------
TEST(RoboticsOracle, FuzzTopo00) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_00"));
}
TEST(RoboticsOracle, FuzzTopo01) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_01"));
}
TEST(RoboticsOracle, FuzzTopo02) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_02"));
}
TEST(RoboticsOracle, FuzzTopo03) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_03"));
}
TEST(RoboticsOracle, FuzzTopo04) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_04"));
}
TEST(RoboticsOracle, FuzzTopo05) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_05"));
}
// §8.3 task 6 (counterexample capture -- adjudicated by hostile review, not
// a disasm src/include fix): this case's `logDetM` diverges from Simbody at
// 5 of its 10 states (worst: state index 7, residual 2.7663761414942201e-07
// at |logDetM|~28.9, i.e. 9.57e-9 RELATIVE) by ~27x the old §6 stage-5
// ABSOLUTE tolerance of 1e-8 -- the divergence is an absolute-tolerance
// artifact on a large-|logDetM| system, not an O(1)-relative gap. `min-eig(D)`
// stays ~1e-3 on both engines throughout (nowhere near the 1e-12 lock, §2b),
// so it is not the null-space-lock discontinuity either. Root cause,
// confirmed by hostile line-by-line review of `RobotEngine::calcLogDetM`
// against Simbody's `calcDetMPass2Outward` (term-for-term match): benign
// conditioning-limited FP accumulation -- logDetM sums Sigma ln|D_k| over
// all bodies with no cancellation (unlike udot/A_GB), and this regime
// (11 bodies, quaternion-joint-heavy: 3 Free + 3 FreeLine + 2 SphericalCoords,
// per-body mass spanning ~2.5e5) is exactly the "§8.1 scaling-blindness" the
// fuzz batch exists to surface. Verdict: BENIGN, not a port bug -- kept
// ENABLED as a regression guard under the relative kLogDetMFuzzTol (see its
// definition above), which gives ~10x headroom over this state's residual
// while staying 6+ orders below an O(1)-relative transcription-bug-sized
// divergence. The exact (topology, state) pair remains reproducible without
// this whole 10-state batch via the standalone one-state fixture
// "FuzzCounterexample_LogDetM_Topo06State7" (same seed, same RNG replay,
// Robosample/tools/gen_robotics_oracle.cpp::writeFuzzCounterexample).
TEST(RoboticsOracle, FuzzTopo06) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_06"));
}
TEST(RoboticsOracle, FuzzTopo07) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_07"));
}
TEST(RoboticsOracle, FuzzTopo08) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_08"));
}
TEST(RoboticsOracle, FuzzTopo09) {
    runOracleFuzzCase(robotics_oracle_loader::loadFuzzCase(kFixtureDir, "FuzzTopo_09"));
}
