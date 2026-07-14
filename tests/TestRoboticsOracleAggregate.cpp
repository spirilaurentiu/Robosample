// ============================================================================
//  TestRoboticsOracleAggregate.cpp -- the live-Simbody differential oracle,
//  aggregate/stress-case family (docs/specs/robotics-oracle-differential.md
//  §8.1), split out of TestRoboticsOracle.cpp (TEST-005). One `TEST` per
//  §8.1 stress case (depth chain, ill-conditioned Ball+Weld), run through
//  the shared `runOracleAggregateCase` staged comparison
//  (tests/support/RoboticsOracleRunners.hpp) -- AGGREGATE invariants +
//  end-products only (§6.1/§9, J1-B), never element-wise DI/PPlus/G, per
//  that header's own rationale.
// ============================================================================
#include <gtest/gtest.h>

#include "support/RoboticsOracleRunners.hpp"

using rtest::kFixtureDir;
using rtest::runOracleAggregateCase;

// ---- §8.1 stress cases (aggregate invariants only, §6.1/§9 J1-B) ----

// Depth-stress: a >=50-body linear Torsion chain -- targets stage-1
// kinematic product-chain drift (§8.1 mechanism 1).
TEST(RoboticsOracle, DepthChainStress) {
    runOracleAggregateCase(robotics_oracle_loader::loadAggregateCase(kFixtureDir, "DepthChain"));
}

// Conditioning-stress: a light Ball-jointed root carrying a heavy Weld
// (Rigid) child at a lever-arm offset, cond(D) ~ 1e6-1e8, min-eig(D) kept
// well above the 1e-12 lock (§8.1 mechanism 2).
TEST(RoboticsOracle, ConditioningStress) {
    runOracleAggregateCase(robotics_oracle_loader::loadAggregateCase(kFixtureDir, "ConditioningStress"));
}
