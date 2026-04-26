/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/**
 * IntegratorTests.cpp
 *
 * Pendulum-system integration tests for SimMath integrators.
 * Refactored to C++17 + Google Test.
 *
 * Each integrator variant is exercised via a value-parameterised fixture.
 * Shared mutable event-handler state is kept in a plain struct owned by the
 * fixture, eliminating all global/static data.
 */

#include <cmath>
#include <functional>
#include <gtest/gtest.h>
#include <memory>
#include <string_view>

#include "simmath/CPodesIntegrator.h"
#include "simmath/ExplicitEulerIntegrator.h"
#include "simmath/RungeKutta2Integrator.h"
#include "simmath/RungeKutta3Integrator.h"
#include "simmath/RungeKuttaFeldbergIntegrator.h"
#include "simmath/RungeKuttaMersonIntegrator.h"
#include "simmath/SemiExplicitEuler2Integrator.h"
#include "simmath/SemiExplicitEulerIntegrator.h"
#include "simmath/TimeStepper.h"
#include "simmath/VerletIntegrator.h"

#include "PendulumSystem.h"
#include "SimTKcommon.h"
#include "SimTKmath.h"

using namespace SimTK;

// -----------------------------------------------------------------------------
// Shared mutable state for all event handlers/reporters in one simulation run.
// Passed by pointer so handlers can mutate and query it without statics.
// -----------------------------------------------------------------------------
struct HandlerState {
    int periodicHandlerCount = 0;
    int zeroPositionCount = 0;
    Real zeroPositionLastTime = 0.0;
    bool hasAccelerated = false;
    int zeroVelocityCount = 0;
    Real zeroVelocityLastTime = 0.0;
    int periodicReporterCount = 0;
    bool onceOnlyOccurred = false;
    int discontinuousCount = 0;

    void reset() {
        *this = {};
    }
};

// -----------------------------------------------------------------------------
// Event handlers
// -----------------------------------------------------------------------------

class PeriodicHandler final : public PeriodicEventHandler {
    public:
    explicit PeriodicHandler(HandlerState& state)
        : PeriodicEventHandler(1.0)
        , state_(state) {
    }

    void handleEvent(State& s, Real /*accuracy*/, bool& /*terminate*/) const override {
        EXPECT_EQ(s.getTime(), getNextEventTime(s, true));
        ++state_.periodicHandlerCount;
    }

    private:
    HandlerState& state_;
};

class ZeroPositionHandler final : public TriggeredEventHandler {
    public:
    ZeroPositionHandler(PendulumSystem& pendulum, HandlerState& state)
        : TriggeredEventHandler(Stage::Velocity)
        , pendulum_(pendulum)
        , state_(state) {
    }

    Real getValue(const State& s) const override {
        return s.getQ(pendulum_.getGuts().getSubsysIndex())[0];
    }

    void handleEvent(State& s, Real /*accuracy*/, bool& /*terminate*/) const override {
        const Real x = s.getQ(pendulum_.getGuts().getSubsysIndex())[0];
        EXPECT_LT(std::abs(x), 0.01);
        EXPECT_GT(s.getTime(), state_.zeroPositionLastTime);

        ++state_.zeroPositionCount;
        state_.zeroPositionLastTime = s.getTime();

        // At t > 7, boost KE by ×1.5 (multiply speed by √1.5) once.
        if (s.getTime() > 7.0 && !state_.hasAccelerated) {
            state_.hasAccelerated = true;
            s.updU(pendulum_.getGuts().getSubsysIndex()) *= std::sqrt(1.5);
        }
    }

    private:
    PendulumSystem& pendulum_;
    HandlerState& state_;
};

class ZeroVelocityHandler final : public TriggeredEventHandler {
    public:
    ZeroVelocityHandler(PendulumSystem& pendulum, HandlerState& state)
        : TriggeredEventHandler(Stage::Velocity)
        , pendulum_(pendulum)
        , state_(state) {
        getTriggerInfo().setTriggerOnFallingSignTransition(false);
    }

    Real getValue(const State& s) const override {
        return s.getU(pendulum_.getGuts().getSubsysIndex())[0];
    }

    void handleEvent(State& s, Real /*accuracy*/, bool& /*terminate*/) const override {
        const Vector u = s.getU(pendulum_.getGuts().getSubsysIndex());
        EXPECT_LT(std::abs(u[0]), 0.01);
        EXPECT_GT(s.getTime(), state_.zeroVelocityLastTime);

        ++state_.zeroVelocityCount;
        state_.zeroVelocityLastTime = s.getTime();
    }

    private:
    PendulumSystem& pendulum_;
    HandlerState& state_;
};

// -----------------------------------------------------------------------------
// Event reporters
// -----------------------------------------------------------------------------

class PeriodicReporter final : public PeriodicEventReporter {
    public:
    PeriodicReporter(PendulumSystem& pendulum, HandlerState& state)
        : PeriodicEventReporter(1.0)
        , pendulum_(pendulum)
        , state_(state) {
    }

    void handleEvent(const State& s) const override {
        EXPECT_EQ(s.getTime(), getNextEventTime(s, true));
        ++state_.periodicReporterCount;

        // Verify energy conservation (within 5 %).
        const SubsystemIndex idx = pendulum_.getGuts().getSubsysIndex();
        const Vector q = s.getQ(idx);
        const Vector u = s.getU(idx);
        const Real mass = pendulum_.getMass(s);
        const Real gravity = pendulum_.getGravity(s);

        const Real energy = mass * (0.5 * (u[0] * u[0] + u[1] * u[1]) + gravity * (1.0 + q[1]));
        const Real expected = mass * gravity * (state_.hasAccelerated ? 1.5 : 1.0);

        EXPECT_NEAR(energy / expected, 1.0, 0.05);
    }

    private:
    PendulumSystem& pendulum_;
    HandlerState& state_;
};

class OnceOnlyReporter final : public ScheduledEventReporter {
    public:
    explicit OnceOnlyReporter(HandlerState& state)
        : state_(state) {
    }

    Real getNextEventTime(const State&, bool /*includeCurrentTime*/) const override {
        return 5.0;
    }

    void handleEvent(const State& /*s*/) const override {
        EXPECT_FALSE(state_.onceOnlyOccurred);
        state_.onceOnlyOccurred = true;
    }

    private:
    HandlerState& state_;
};

class DiscontinuousReporter final : public TriggeredEventReporter {
    public:
    explicit DiscontinuousReporter(HandlerState& state)
        : TriggeredEventReporter(Stage::Time)
        , state_(state) {
    }

    Real getValue(const State& s) const override {
        const Real step = std::fmod(std::floor(s.getTime()), 4.0);
        if (step == 0.0) {
            return 1.0;
        }
        if (step == 2.0) {
            return -1.0;
        }
        return 0.0;
    }

    void handleEvent(const State& s) const override {
        // Should only trigger when the value crosses zero (at t ≡ 1 mod 2).
        EXPECT_NEAR(std::fmod(s.getTime(), 2.0), 1.0, 0.01);
        ++state_.discontinuousCount;
    }

    private:
    HandlerState& state_;
};

// -----------------------------------------------------------------------------
// Integrator factory configuration
// -----------------------------------------------------------------------------
struct IntegratorConfig {
    std::string name;
    // Factory returns a freshly constructed integrator; called once per run.
    std::function<std::unique_ptr<Integrator>(PendulumSystem&)> factory;
    Real accuracy = 1e-4;
};

// -----------------------------------------------------------------------------
// Test fixture
// -----------------------------------------------------------------------------
class IntegratorTest : public ::testing::TestWithParam<IntegratorConfig> {
    protected:
    // Pointers to interval-configurable handlers (non-owning; system owns them).
    PeriodicHandler* periodicHandler_ = nullptr;
    PeriodicReporter* periodicReporter_ = nullptr;

    HandlerState state_;
    PendulumSystem sys_;

    void SetUp() override {
        // System takes raw-pointer ownership (SimTK convention).
        sys_.addEventHandler(new ZeroVelocityHandler(sys_, state_));
        sys_.addEventHandler(periodicHandler_ = new PeriodicHandler(state_));
        sys_.addEventHandler(new ZeroPositionHandler(sys_, state_));
        sys_.addEventReporter(periodicReporter_ = new PeriodicReporter(sys_, state_));
        sys_.addEventReporter(new OnceOnlyReporter(state_));
        sys_.addEventReporter(new DiscontinuousReporter(state_));
        sys_.realizeTopology();
    }

    // Run a full 0 -> 20.003 s simulation and verify all invariants.
    void runIntegrator(Integrator& integ, Real accuracy) {
        state_.reset();

        constexpr Real t0 = 0.0;
        constexpr Real tFinal = 20.003;
        const Real qi[] = {1.0, 0.0};
        const Real ui[] = {0.0, 0.0};

        sys_.setDefaultMass(10);
        sys_.setDefaultTimeAndState(t0, Vector(2, qi), Vector(2, ui));

        integ.setAccuracy(accuracy);
        integ.setConstraintTolerance(1e-4);
        integ.setFinalTime(tFinal);

        TimeStepper ts(sys_);
        ts.setIntegrator(integ);
        ts.initialize(sys_.getDefaultState());

        ASSERT_EQ(ts.getTime(), 0.0);

        // Fixed-size steps 1 -> 4 s.
        for (Real t = 1.0; t < 5.0; t += 1.0) {
            ts.stepTo(t);
            EXPECT_EQ(ts.getTime(), t);
        }
        EXPECT_FALSE(state_.onceOnlyOccurred);
        EXPECT_FALSE(state_.hasAccelerated);

        // Variable-size steps 5 -> 10 s.
        static Random::Uniform rng(0.0, 1.0);
        for (Real t = 5.0; t < 10.0; t += rng.getValue()) {
            EXPECT_EQ(state_.onceOnlyOccurred, ts.getTime() >= 5.0);
            ts.stepTo(t);
            EXPECT_EQ(ts.getTime(), t);
        }
        EXPECT_TRUE(state_.hasAccelerated);

        // One large step that should stop at tFinal.
        ts.stepTo(50.0);
        EXPECT_EQ(ts.getTime(), tFinal);
        EXPECT_EQ(integ.getTerminationReason(), Integrator::ReachedFinalTime);

        // Sanity-check event counts.
        EXPECT_GT(state_.zeroVelocityCount, 10);
        EXPECT_GT(state_.zeroPositionCount, 10);

        const auto periodicHandlerExpected =
            static_cast<int>(ts.getTime() / periodicHandler_->getEventInterval()) + 1;
        EXPECT_EQ(state_.periodicHandlerCount, periodicHandlerExpected);

        const auto periodicReporterExpected =
            static_cast<int>(ts.getTime() / periodicReporter_->getEventInterval()) + 1;
        EXPECT_EQ(state_.periodicReporterCount, periodicReporterExpected);

        EXPECT_EQ(state_.discontinuousCount, static_cast<int>(ts.getTime() / 2.0));
    }

    // Run both normal mode and single-step mode for one IntegratorConfig.
    void testConfig(const IntegratorConfig& cfg) {
        for (int i = 0; i < 4; ++i) {
            periodicHandler_->setEventInterval(i == 0 || i == 1 ? 0.01 : 2.0);
            periodicReporter_->setEventInterval(i == 0 || i == 2 ? 0.015 : 1.5);

            auto integ = cfg.factory(sys_);

            // Normal mode.
            runIntegrator(*integ, cfg.accuracy);

            // Single-step mode.
            integ->setReturnEveryInternalStep(true);
            runIntegrator(*integ, cfg.accuracy);
        }
    }
};

// -----------------------------------------------------------------------------
// Parameterised test body
// -----------------------------------------------------------------------------
TEST_P(IntegratorTest, RunsAndPassesAllChecks) {
    testConfig(GetParam());
}

// -----------------------------------------------------------------------------
// Integrator catalogue
// -----------------------------------------------------------------------------
static const IntegratorConfig kIntegratorConfigs[] = {
    {"VerletIntegrator",
     [](PendulumSystem& sys) {
         return std::make_unique<VerletIntegrator>(sys);
     }},
    {"RungeKutta2Integrator",
     [](PendulumSystem& sys) {
         return std::make_unique<RungeKutta2Integrator>(sys);
     }},
    {"RungeKutta3Integrator",
     [](PendulumSystem& sys) {
         return std::make_unique<RungeKutta3Integrator>(sys);
     }},
    {"RungeKuttaMersonIntegrator",
     [](PendulumSystem& sys) {
         return std::make_unique<RungeKuttaMersonIntegrator>(sys);
     }},
    {"RungeKuttaFeldbergIntegrator",
     [](PendulumSystem& sys) {
         return std::make_unique<RungeKuttaFeldbergIntegrator>(sys);
     }},
    {"SemiExplicitEuler2Integrator",
     [](PendulumSystem& sys) {
         return std::make_unique<SemiExplicitEuler2Integrator>(sys);
     },
     /*accuracy=*/1e-6},
    {"SemiExplicitEulerIntegrator",
     [](PendulumSystem& sys) {
         return std::make_unique<SemiExplicitEulerIntegrator>(sys, /*fixedStep=*/1e-4);
     }},
    {"ExplicitEulerIntegrator",
     [](PendulumSystem& sys) {
         return std::make_unique<ExplicitEulerIntegrator>(sys);
     },
     /*accuracy=*/1e-7},
};

INSTANTIATE_TEST_SUITE_P(AllIntegrators,
                         IntegratorTest,
                         ::testing::ValuesIn(kIntegratorConfigs),
                         [](const ::testing::TestParamInfo<IntegratorConfig>& info) {
                             return info.param.name;
                         });

// -----------------------------------------------------------------------------
// CPodesIntegrator — needs extra method-specific logic so kept as its own suite.
// -----------------------------------------------------------------------------
class CPodesIntegratorTest : public IntegratorTest {};

TEST_F(CPodesIntegratorTest, BDF_NormalAndSingleStep) {
    for (int i = 0; i < 4; ++i) {
        periodicHandler_->setEventInterval(i == 0 || i == 1 ? 0.01 : 2.0);
        periodicReporter_->setEventInterval(i == 0 || i == 2 ? 0.015 : 1.5);

        CPodesIntegrator integ(sys_, CPodes::BDF);
        runIntegrator(integ, 1e-4);

        integ.setReturnEveryInternalStep(true);
        runIntegrator(integ, 1e-4);
    }
}

TEST_F(CPodesIntegratorTest, Adams_NormalAndSingleStep) {
    for (int i = 0; i < 4; ++i) {
        periodicHandler_->setEventInterval(i == 0 || i == 1 ? 0.01 : 2.0);
        periodicReporter_->setEventInterval(i == 0 || i == 2 ? 0.015 : 1.5);

        CPodesIntegrator integ(sys_, CPodes::Adams);
        runIntegrator(integ, 1e-4);

        integ.setReturnEveryInternalStep(true);
        runIntegrator(integ, 1e-4);
    }
}

TEST_F(CPodesIntegratorTest, SetUseCPodesProjection_AfterInitThrows) {
    periodicHandler_->setEventInterval(0.01);
    periodicReporter_->setEventInterval(0.015);

    CPodesIntegrator integ(sys_, CPodes::Adams);
    runIntegrator(integ, 1e-4); // initialises the integrator internally

    // Calling setUseCPodesProjection() on an already-initialised integrator
    // must throw.
    EXPECT_THROW(integ.setUseCPodesProjection(), std::exception);
}

TEST_F(CPodesIntegratorTest, BDF_WithCPodesProjection) {
    for (int i = 0; i < 4; ++i) {
        periodicHandler_->setEventInterval(i == 0 || i == 1 ? 0.01 : 2.0);
        periodicReporter_->setEventInterval(i == 0 || i == 2 ? 0.015 : 1.5);

        CPodesIntegrator integ(sys_, CPodes::BDF);
        integ.setUseCPodesProjection();
        runIntegrator(integ, 1e-4);
    }
}