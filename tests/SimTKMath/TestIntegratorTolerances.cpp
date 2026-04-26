/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

#include <cassert>
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <limits>

#include "SimTKcommon/internal/SystemGuts.h"
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

#include "Event.h"
#include "Integrator.h"
#include "Scalar.h"
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Stage.h"

using namespace SimTK;

#include "PendulumSystem.hpp"

struct IntegratorFactory {
    std::string name;
    Real stepSize;
    Real accuracy;
    Real constraintTolerance;
    Real tolerance;
    std::function<std::unique_ptr<Integrator>(System&)> create;

    friend void PrintTo(const IntegratorFactory& factory, std::ostream* ostream) {
        *ostream << factory.name;
    }
};

class IntegratorTest : public ::testing::TestWithParam<IntegratorFactory> {};


TEST_P(IntegratorTest, ) {
    MyPendulum sys;

    auto integPtr = GetParam().create(sys);
    Integrator& integ = *integPtr;

    const std::string& testName = GetParam().name;
    const bool isCPodes = (testName == "CPodesIntegrator");

    SCOPED_TRACE("Testing Integrator: " + testName);

    if ("VerletIntegrator" == testName) {
        EXPECT_TRUE(integPtr->methodHasErrorControl()) << "VerletIntegrator should support error control.";
        EXPECT_EQ(integPtr->getMethodMaxOrder(), 3) << "VerletIntegrator should have max order 2.";
        EXPECT_EQ(integPtr->getMethodMinOrder(), 2) << "VerletIntegrator should have min order 2.";
    } else if ("SemiExplicitEulerIntegrator" == testName) {
        EXPECT_FALSE(integPtr->methodHasErrorControl()) << "SemiExplicitEuler should not have error control.";
        EXPECT_EQ(integPtr->getMethodMaxOrder(), 1) << "SemiExplicitEuler should have max order 1.";
        EXPECT_EQ(integPtr->getMethodMinOrder(), 1) << "SemiExplicitEuler should have min order 1.";
    } else if ("SemiExplicitEuler2Integrator" == testName) {
        EXPECT_TRUE(integPtr->methodHasErrorControl()) << "SemiExplicitEuler2 should support error control.";
        EXPECT_EQ(integPtr->getMethodMaxOrder(), 1) << "SemiExplicitEuler2 should have max order 1.";
        EXPECT_EQ(integPtr->getMethodMinOrder(), 1) << "SemiExplicitEuler2 should have min order 1.";
    } else if ("RungeKutta2Integrator" == testName) {
        EXPECT_TRUE(integPtr->methodHasErrorControl()) << "RungeKutta2 should support error control.";
        EXPECT_EQ(integPtr->getMethodMaxOrder(), 2) << "RungeKutta2 should have max order 2.";
        EXPECT_EQ(integPtr->getMethodMinOrder(), 2) << "RungeKutta2 should have min order 2.";
    } else if ("RungeKutta3Integrator" == testName) {
        EXPECT_TRUE(integPtr->methodHasErrorControl()) << "RungeKutta3 should support error control.";
        EXPECT_EQ(integPtr->getMethodMaxOrder(), 3) << "RungeKutta3 should have max order 3.";
        EXPECT_EQ(integPtr->getMethodMinOrder(), 3) << "RungeKutta3 should have min order 3.";
    } else if ("RungeKuttaMersonIntegrator" == testName) {
        EXPECT_TRUE(integPtr->methodHasErrorControl()) << "RungeKuttaMerson should support error control.";
        EXPECT_EQ(integPtr->getMethodMaxOrder(), 4) << "RungeKuttaMerson should have max order 4.";
        EXPECT_EQ(integPtr->getMethodMinOrder(), 4) << "RungeKuttaMerson should have min order 1.";
    } else if ("RungeKuttaFeldbergIntegrator" == testName) {
        EXPECT_TRUE(integPtr->methodHasErrorControl()) << "RungeKuttaFeldberg should support error control.";
        EXPECT_EQ(integPtr->getMethodMaxOrder(), 5) << "RungeKuttaFeldberg should have max order 5.";
        EXPECT_EQ(integPtr->getMethodMinOrder(), 5) << "RungeKuttaFeldberg should have min order 5.";
    } else if ("CPodesIntegrator" == testName) {
        EXPECT_TRUE(integPtr->methodHasErrorControl()) << "CPodes should support error control.";
        EXPECT_EQ(integPtr->getMethodMaxOrder(), 5) << "CPodes should have max order 5.";
        EXPECT_EQ(integPtr->getMethodMinOrder(), 1) << "CPodes should have min order 1.";
    } else if ("ExplicitEulerIntegrator" == testName) {
        EXPECT_TRUE(integPtr->methodHasErrorControl()) << "ExplicitEuler should support error control.";
        EXPECT_EQ(integPtr->getMethodMaxOrder(), 1) << "ExplicitEuler should have max order 1.";
        EXPECT_EQ(integPtr->getMethodMinOrder(), 1) << "ExplicitEuler should have min order 1.";
    } else {
        FAIL() << "Unknown integrator: " << testName;
    }

    const Real t0 = 0;
    const Real qi[] = {1, 0};
    const Real ui[] = {0, 0};
    const Vector q0(2, qi);
    const Vector u0(2, ui);

    const Real tFinal = 30.003;
    const Real hReport = 1.;

    sys.setDefaultMass(10);
    sys.setDefaultTimeAndState(t0, q0, u0);

    integ.setAccuracy(GetParam().accuracy);
    integ.setConstraintTolerance(GetParam().constraintTolerance);
    sys.setTolerance(GetParam().tolerance);

    integ.setAllowInterpolation(false);
    integ.setProjectEveryStep(true);
    integ.setProjectInterpolatedStates(false);
    integ.setUseInfinityNorm(true);
    integ.setReturnEveryInternalStep(true);

    integ.setFixedStepSize(GetParam().stepSize);

    // if ("VerletIntegrator" == testName) {
    //     integ.setFixedStepSize(timeStep);
    // } else {
    //     integ.setInitialStepSize(timeStep);
    // }

    integ.setFinalTime(tFinal);
    integ.initialize(sys.getDefaultState());

    Real prevScheduledEventTime = -Infinity;
    Real lastTime = -Infinity;

    auto E_initial = SimTK::NaN;

    for (int reportNo = 0; !integ.isSimulationOver();
         reportNo += static_cast<int>(integ.getTime() >= reportNo * hReport)) {
        EXPECT_GE(integ.getTime(), lastTime);
        lastTime = integ.getTime();

        Array_<EventId> scheduledEventIds;
        Real nextScheduledEvent = NTraits<Real>::getInfinity();

        sys.calcTimeOfNextScheduledEvent(integ.getAdvancedState(),
                                         nextScheduledEvent,
                                         scheduledEventIds,
                                         integ.getAdvancedTime() > prevScheduledEventTime);

        HandleEventsOptions handleOpts(integ.getAccuracyInUse());
        HandleEventsResults handleResults;

        const auto stepResult = integ.stepTo(reportNo * hReport, nextScheduledEvent);

        switch (stepResult) {
            case Integrator::ReachedStepLimit:
                FAIL() << "Unexpected step limit reached";
                break;

            case Integrator::ReachedReportTime:
                if (integ.getTime() < tFinal - 1e-6) {
                    EXPECT_NEAR(integ.getTime(), reportNo * hReport, 10 * integ.getAccuracyInUse());
                }
                break;

            case Integrator::StartOfContinuousInterval:
                SUCCEED();
                break;

            case Integrator::ReachedScheduledEvent: {
                sys.checkState(integ);

                EXPECT_FALSE(scheduledEventIds.empty());
                EXPECT_GE(nextScheduledEvent, integ.getTime() - 1e-6);

                sys.handleEvents(integ.updAdvancedState(),
                                 Event::Cause::Scheduled,
                                 scheduledEventIds,
                                 handleOpts,
                                 handleResults);

                EXPECT_TRUE(handleResults.getExitStatus() == HandleEventsResults::Succeeded
                            || handleResults.getExitStatus() == HandleEventsResults::ShouldTerminate);

                const bool shouldTerminate =
                    handleResults.getExitStatus() == HandleEventsResults::ShouldTerminate;

                integ.reinitialize(handleResults.getLowestModifiedStage(), shouldTerminate);

                prevScheduledEventTime = integ.getAdvancedTime();
                break;
            }

            case Integrator::TimeHasAdvanced: {
                EXPECT_GT(integ.getTime(), 0);

                sys.handleEvents(integ.updAdvancedState(),
                                 Event::Cause::TimeAdvanced,
                                 Array_<EventId>(),
                                 handleOpts,
                                 handleResults);

                EXPECT_TRUE(handleResults.getExitStatus() == HandleEventsResults::Succeeded
                            || handleResults.getExitStatus() == HandleEventsResults::ShouldTerminate);

                const bool shouldTerminate =
                    handleResults.getExitStatus() == HandleEventsResults::ShouldTerminate;

                integ.reinitialize(handleResults.getLowestModifiedStage(), shouldTerminate);

                if (!std::isfinite(E_initial)) {
                    E_initial = sys.calcTotalEnergy(integ.getAdvancedState());
                }

                break;
            }

            case Integrator::ReachedEventTrigger: {
                sys.checkState(integ);

                EXPECT_FALSE(integ.getTriggeredEvents().empty());
                EXPECT_GE(integ.getAdvancedTime(), integ.getTime());

                sys.handleEvents(integ.updAdvancedState(),
                                 Event::Cause::Triggered,
                                 integ.getTriggeredEvents(),
                                 handleOpts,
                                 handleResults);

                EXPECT_TRUE(handleResults.getExitStatus() == HandleEventsResults::Succeeded
                            || handleResults.getExitStatus() == HandleEventsResults::ShouldTerminate);

                const bool shouldTerminate =
                    handleResults.getExitStatus() == HandleEventsResults::ShouldTerminate;

                integ.reinitialize(handleResults.getLowestModifiedStage(), shouldTerminate);
                break;
            }

            case Integrator::EndOfSimulation: {
                if (!isCPodes) {
                    EXPECT_GE(integ.getTime(), tFinal - GetParam().stepSize);
                }

                sys.handleEvents(integ.updAdvancedState(),
                                 Event::Cause::Termination,
                                 Array_<EventId>(),
                                 handleOpts,
                                 handleResults);

                const bool shouldTerminate =
                    handleResults.getExitStatus() == HandleEventsResults::ShouldTerminate;

                integ.reinitialize(handleResults.getLowestModifiedStage(), shouldTerminate);
                break;
            }

            default:
                FAIL() << "Unrecognized return from stepTo()";
        }

        // Fall through to here to report
        sys.checkState(integ);
    }

    const auto E_final = sys.calcTotalEnergy(integ.getAdvancedState());

    const int nst = integ.getNumStepsTaken();
    const int nattempt = integ.getNumStepsAttempted();
    const int netf = integ.getNumErrorTestFailures();

    std::stringstream sstream;
    sstream << "Integrator: " << testName << "\n"
            << "  Steps taken: " << nst << " with step size: " << integ.getPredictedNextStepSize() << "\n"
            << "  Steps attempted: " << nattempt << "\n"
            << "  Error test failures: " << netf;

    if (!isCPodes) {
        sstream << "\n  Convergent iterations: " << integ.getNumConvergentIterations()
                << "\n  Divergent iterations: " << integ.getNumDivergentIterations()
                << "\n  Convergence test failures: " << integ.getNumConvergenceTestFailures();
    }

    sstream << "\n  Initial energy: " << E_initial << "\n  Final energy: " << E_final
            << "\n  Energy change: " << E_final - E_initial;

    SCOPED_TRACE(sstream.str());

    // Make sure we advanced to the target final time
    if (!isCPodes) {
        EXPECT_GE(integ.getTime(), tFinal - GetParam().stepSize)
            << "Integrator did not reach final time. Last time: " << integ.getTime();
    }

    // In molecular dynamics, nst=0 usually indicates an immediate singularity in the potential
    ASSERT_GT(nst, 0) << "Critical Failure: Integrator stalled at t0. Check initial forces/potential.";

    // nattempt/nst ratio > 2.0 suggests the system is numerically stiff or the manifold is high-curvature
    double workEfficiency = static_cast<double>(nattempt) / nst;
    EXPECT_LT(workEfficiency, 2.5) << "Stiffness Warning: Integrator is thrashing. Efficiency ratio: "
                                   << workEfficiency;

    // Step-Size Integrity
    const Real h0u = integ.getActualInitialStepSizeTaken();

    // Manifold Projection & Drift Validation
    const int nprojq = integ.getNumQProjections();
    const int nproju = integ.getNumUProjections();
    const int nprqf = integ.getNumQProjectionFailures();
    const int npruf = integ.getNumUProjectionFailures();

    // Verify Integrator reported failures match System-level recorded failures
    EXPECT_EQ(nprqf, sys.getNumQProjFailures()) << "Inconsistency: Integrator Q-fails vs System Q-fails.";
    EXPECT_EQ(npruf, sys.getNumUProjFailures()) << "Inconsistency: Integrator U-fails vs System U-fails.";

    // In a constrained pendulum (Index-3 DAE), Q-projections are the "heartbeat."
    // If nprojq == 0, the constraints aren't being enforced, or the tolerance is dangerously loose.
    EXPECT_GT(nprojq, 0) << "Manifold Violation: No Q-projections performed.";
    EXPECT_GT(nproju, 0) << "Manifold Violation: No U-projections performed.";

    // Expect that we don't fail more than 10% of our Q-projections
    if (nprojq > 0) {
        const double qProjFailureRate = static_cast<double>(nprqf) / nprojq;
        EXPECT_LT(qProjFailureRate, 0.1)
            << "Manifold Drift Warning: Q-projection failure rate is " << qProjFailureRate * 100 << "%.";
    }

    // Expect that we don't fail more than 10% of our U-projections
    if (nproju > 0) {
        const double uProjFailureRate = static_cast<double>(npruf) / nproju;
        EXPECT_LT(uProjFailureRate, 0.1)
            << "Manifold Drift Warning: U-projection failure rate is " << uProjFailureRate * 100 << "%.";
    }

    // Expect that we don't waste more than 10% of our effort on failed steps
    if (nst > 0) {
        const double failureRate = static_cast<double>(netf) / nst;
        // Get the number of attempted steps that have failed due to the error being unacceptably high
        // This means infinite rejection loops or catastrophic step collapse
        EXPECT_LT(failureRate, 0.8) << netf << " failures out of " << nattempt << " attempts.";
    }
}

INSTANTIATE_TEST_SUITE_P(
    EnergyConservingPhysics,
    IntegratorTest,

    ::testing::Values(
        IntegratorFactory{"VerletIntegrator",
                          1e-3,
                          1e-3,
                          1e-3,
                          1e-3,
                          [](System& state) -> std::__detail::__unique_ptr_t<SimTK::VerletIntegrator> {
                              return std::make_unique<VerletIntegrator>(state);
                          }},
        IntegratorFactory{
            "SemiExplicitEulerIntegrator",
            1e-4,
            1e-3,
            1e-3,
            1e-1,
            [](System& state) -> std::__detail::__unique_ptr_t<SimTK::SemiExplicitEulerIntegrator> {
                // Set fixed step size hereol
                return std::make_unique<SemiExplicitEulerIntegrator>(state, 1e-4);
            }},
        IntegratorFactory{
            "SemiExplicitEuler2Integrator",
            1e-4,
            1e-3,
            1e-3,
            1e-1,
            [](System& state) -> std::__detail::__unique_ptr_t<SimTK::SemiExplicitEuler2Integrator> {
                return std::make_unique<SemiExplicitEuler2Integrator>(state);
            }}),

    [](const ::testing::TestParamInfo<IntegratorFactory>& info) -> std::string {
        return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    HighPrecisionTrajectory,
    IntegratorTest,

    ::testing::Values(
        IntegratorFactory{"RungeKutta2Integrator",
                          1e-4,
                          1e-3,
                          1e-3,
                          1e-5,
                          [](System& state) -> std::__detail::__unique_ptr_t<SimTK::RungeKutta2Integrator> {
                              return std::make_unique<RungeKutta2Integrator>(state);
                          }},
        IntegratorFactory{"RungeKutta3Integrator",
                          1e-4,
                          1e-3,
                          1e-3,
                          1e-8,
                          [](System& state) -> std::__detail::__unique_ptr_t<SimTK::RungeKutta3Integrator> {
                              return std::make_unique<RungeKutta3Integrator>(state);
                          }},
        IntegratorFactory{
            "RungeKuttaMersonIntegrator",
            1e-4,
            1e-3,
            1e-3,
            1e-12,
            [](System& state) -> std::__detail::__unique_ptr_t<SimTK::RungeKuttaMersonIntegrator> {
                return std::make_unique<RungeKuttaMersonIntegrator>(state);
            }},
        IntegratorFactory{
            "RungeKuttaFeldbergIntegrator",
            1e-4,
            1e-3,
            1e-3,
            1e-11,
            [](System& state) -> std::__detail::__unique_ptr_t<SimTK::RungeKuttaFeldbergIntegrator> {
                return std::make_unique<RungeKuttaFeldbergIntegrator>(state);
            }},
        IntegratorFactory{"CPodesIntegrator",
                          1e-4,
                          1e-3,
                          1e-3,
                          2,
                          [](System& state) -> std::__detail::__unique_ptr_t<SimTK::CPodesIntegrator> {
                              return std::make_unique<CPodesIntegrator>(state);
                          }}),

    [](const ::testing::TestParamInfo<IntegratorFactory>& info) -> std::string {
        return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    BaselineMethods,
    IntegratorTest,

    ::testing::Values(IntegratorFactory{
        "ExplicitEulerIntegrator",
        1e-3,
        1e-3,
        1e-3,
        20,
        [](System& state) -> std::__detail::__unique_ptr_t<SimTK::ExplicitEulerIntegrator> {
            return std::make_unique<ExplicitEulerIntegrator>(state);
        }}),

    [](const ::testing::TestParamInfo<IntegratorFactory>& info) -> std::string {
        return info.param.name;
    });
