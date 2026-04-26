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
#include <gtest/gtest.h>

#include "PendulumSystem.h"
#include "SimTKmath.h"

using namespace SimTK;

class ZeroVelocityHandler : public TriggeredEventHandler {
    public:
    static int eventCount;
    static Real lastEventTime;
    ZeroVelocityHandler(PendulumSystem& pendulum)
        : TriggeredEventHandler(Stage::Velocity)
        , pendulum(pendulum) {
        getTriggerInfo().setTriggerOnFallingSignTransition(false);
    }
    [[nodiscard]] auto getValue(const State& state) const -> Real override {
        return state.getU(pendulum.getGuts().getSubsysIndex())[0];
    }
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const override {
        // This should be triggered when the pendulum reaches its farthest point in the
        // negative direction: q == -1, u == 0.
        const Real q = state.getQ(pendulum.getGuts().getSubsysIndex())[0];
        const Real u = state.getU(pendulum.getGuts().getSubsysIndex())[0];
        EXPECT_LT(std::abs(q + 1.0), 0.05);
        EXPECT_LT(std::abs(u), 0.01);
        EXPECT_GT(state.getTime(), lastEventTime);
        eventCount++;
        lastEventTime = state.getTime();
    }

    private:
    PendulumSystem& pendulum;
};

class PeriodicHandler : public ScheduledEventHandler {
    public:
    static int eventCount;
    static Real lastEventTime;
    [[nodiscard]] auto getNextEventTime(const State&, bool includeCurrentTime) const -> Real override {
        return lastEventTime + 1.5;
    }
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const override {
        // This should be triggered every 1.5 time units.
        EXPECT_EQ(state.getTime(), lastEventTime + 1.5);
        eventCount++;
        lastEventTime = state.getTime();
    }
};

class ZeroPositionReporter : public TriggeredEventReporter {
    public:
    static int eventCount;
    static Real lastEventTime;
    ZeroPositionReporter(PendulumSystem& pendulum)
        : TriggeredEventReporter(Stage::Velocity)
        , pendulum(pendulum) {
    }
    [[nodiscard]] auto getValue(const State& state) const -> Real override {
        return state.getQ(pendulum.getGuts().getSubsysIndex())[0];
    }
    void handleEvent(const State& state) const override {
        // This should be triggered when the pendulum crosses q == 0.
        const Real q = state.getQ(pendulum.getGuts().getSubsysIndex())[0];
        EXPECT_LT(std::abs(q), 0.01);
        EXPECT_GT(state.getTime(), lastEventTime);
        eventCount++;
        lastEventTime = state.getTime();
    }

    private:
    PendulumSystem& pendulum;
};

class PeriodicReporter : public ScheduledEventReporter {
    public:
    static int eventCount;
    static Real lastEventTime;
    [[nodiscard]] auto getNextEventTime(const State&, bool includeCurrentTime) const -> Real override {
        return lastEventTime * 2;
    }
    void handleEvent(const State& state) const override {
        // This should be triggered every 1.5 time units.
        EXPECT_EQ(state.getTime(), lastEventTime * 2);
        eventCount++;
        lastEventTime = state.getTime();
    }
};

int ZeroVelocityHandler::eventCount = 0;
Real ZeroVelocityHandler::lastEventTime = 0.0;
int PeriodicHandler::eventCount = 0;
Real PeriodicHandler::lastEventTime = 0.0;
int ZeroPositionReporter::eventCount = 0;
Real ZeroPositionReporter::lastEventTime = 0.0;
int PeriodicReporter::eventCount = 0;
Real PeriodicReporter::lastEventTime = 0.5;

TEST(TimeStepperTest, TimeStepper) {
    PendulumSystem sys;
    sys.addEventHandler(new ZeroVelocityHandler(sys));
    sys.addEventHandler(new PeriodicHandler());
    sys.addEventReporter(new ZeroPositionReporter(sys));
    sys.addEventReporter(new PeriodicReporter());
    sys.realizeTopology();

    RungeKuttaMersonIntegrator integ(sys);

    const Real t0 = 0;
    const Real qi[] = {1, 0}; // (x,y) = (1,0)
    const Real ui[] = {0, 0}; // v=0
    const Vector q0(2, qi);
    const Vector u0(2, ui);

    sys.setDefaultMass(10);
    sys.setDefaultTimeAndState(t0, q0, u0);

    integ.setAccuracy(1e-2);
    integ.setConstraintTolerance(1e-4);

    const Real tFinal = 20.003;
    const Real hReport = 1.;

    integ.setFinalTime(tFinal);

    TimeStepper ts(sys);
    ts.setIntegrator(integ);
    ts.initialize(sys.getDefaultState());
    EXPECT_EQ(ts.getTime(), 0.0);
    ts.stepTo(10.0);
    EXPECT_EQ(ts.getTime(), 10.0);
    ts.stepTo(50.0);
    EXPECT_EQ(ts.getTime(), tFinal);
    EXPECT_EQ(integ.getTerminationReason(), Integrator::ReachedFinalTime);
    EXPECT_GE(ZeroVelocityHandler::eventCount, 10);
    EXPECT_EQ(PeriodicHandler::eventCount, (int)(ts.getTime() / 1.5));
    EXPECT_GT(ZeroPositionReporter::eventCount, 10);
    EXPECT_EQ(PeriodicReporter::eventCount, (int)(std::log(ts.getTime()) / std::log(2.0)) + 1);
}
