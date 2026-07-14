#pragma once
// ============================================================================
//  RobotFixture.hpp -- a gtest fixture base that owns the RobotModel +
//  RobotState + ConstraintSet + AnalyticForceBridge construction prologue one
//  canonical shape of test setup was hand-rolling per file (TESTS.md section 4):
//  build a model, allocate the state, seed a random-but-valid configuration,
//  realize it to Cartesian, then construct the analytic harmonic bridge whose
//  anchors are captured from that realized start configuration.
//
//  CONTRACT. `SetUp()` runs, in this exact order:
//    1. `buildModel(rng)` -- the derived fixture's builder (mirrors the lambda
//       builders in RobotBuilders.hpp: buildForest/buildChain/buildSingle/
//       buildBentTorsionChain all return a RobotModel).
//    2. `s.allocateFull(m)`.
//    3. `rtest::randomizeState(m, s, rng)` -- seeds q/u with the SAME RNG
//       stream `buildModel` was given, so a derived fixture that seeds its
//       fixture-scope Rng once and does nothing else draws a deterministic,
//       reproducible sample stream keyed only on that one seed.
//    4. `RobotEngine::realizePosition` + `RobotEngine::fillAtomPositionsFromBodies`
//       -- fills atomPosG()/X_GB() so the bridge constructor (which reads
//       atomPosG()) captures real anchors, not stale/zero positions.
//    5. `bridge.emplace(m, s, k)` -- the caller-chosen stiffness `k` (default 0,
//       i.e. KE-only; a derived fixture passes a nonzero k for PE-sensitive
//       suites), anchored at the just-realized start configuration.
//
//  USE ONLY where this is the exact prologue the converted test already ran
//  (TESTS.md section 5.3): same seed, same realize order, same k. A fixture that
//  reorders these steps or changes the seed/k draws a DIFFERENT sample stream
//  and is not a pure test motion.
// ============================================================================

#include <cstdint>
#include <gtest/gtest.h>
#include <optional>

#include "../AnalyticForceBridge.hpp"
#include "../RobotBuilders.hpp"
#include "../TestHelpers.hpp"
#include "Constraints.hpp"
#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"

namespace rtest {

// gtest fixture base for the build -> allocate -> randomize -> realize ->
// bridge prologue. Derived fixtures implement buildModel() and pass their
// fixed seed (and, for PE-sensitive suites, a nonzero stiffness) to the
// constructor so the drawn stream matches what the pre-fixture test body drew.
class RobotFixture : public ::testing::Test {
    public:
    explicit RobotFixture(std::uint64_t seed, robo::Real k = robo::Real(0))
        : rng(seed)
        , k_(k) {}

    protected:
    // Derived fixtures build their topology here. Called once, from SetUp(),
    // with the fixture's own Rng so topology randomization (e.g.
    // buildBentTorsionChain's per-body geometry) and state randomization
    // share one deterministic stream, mirroring the pre-fixture call sites.
    virtual RobotModel buildModel(Rng& rngIn) = 0;

    void SetUp() override {
        m = buildModel(rng);
        s.allocateFull(m);
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        bridge.emplace(m, s, k_);
    }

    RobotModel m;
    RobotState s;
    robo::ConstraintSet cs;
    std::optional<AnalyticForceBridge> bridge;
    Rng rng;

    private:
    robo::Real k_;
};

} // namespace rtest
