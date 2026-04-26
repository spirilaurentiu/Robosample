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

#pragma once

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

// User-defined system to be integrated. This is a kind of SimTK::System.
class MyPendulum;
class MyPendulumGuts : public System::Guts {
    friend class MyPendulum;

    // TOPOLOGY STATE
    SubsystemIndex subsysIndex;

    // TOPOLOGY CACHE
    mutable DiscreteVariableIndex massIndex, lengthIndex, gravityIndex;
    mutable QIndex q0;
    mutable UIndex u0;
    mutable QErrIndex qerr0;
    mutable UErrIndex uerr0;
    mutable UDotErrIndex udoterr0;
    mutable EventTriggerByStageIndex trigger0;
    mutable CacheEntryIndex mgForceIndex; // a cache entry m*g calculated at Dynamics stage
    mutable EventId eventId0, eventId1, eventId2;

    mutable int qProj, qProjFail;
    mutable int uProj, uProjFail;

    mutable SimTK::Real E_start = SimTK::NaN;
    mutable SimTK::Real E_final = SimTK::NaN;
    mutable SimTK::Real E_min = std::numeric_limits<SimTK::Real>::max();
    mutable SimTK::Real E_max = std::numeric_limits<SimTK::Real>::lowest();
    mutable bool firstSample = true;

    mutable Real tolerance;

    public:
    // Index types set themselves invalid on construction.
    MyPendulumGuts() = default;

    inline auto getMyPendulum() const -> const MyPendulum&;

    /*virtual*/ auto cloneImpl() const -> MyPendulumGuts* override {
        return new MyPendulumGuts(*this);
    }

    /////////////////////////////////////////////////////////
    // Implementation of continuous DynamicSystem virtuals //
    /////////////////////////////////////////////////////////

    auto realizeTopologyImpl(State& state) const -> int override;
    auto realizeModelImpl(State& state) const -> int override;
    auto realizeInstanceImpl(const State& state) const -> int override;
    auto realizePositionImpl(const State& state) const -> int override;
    auto realizeVelocityImpl(const State& state) const -> int override;
    auto realizeDynamicsImpl(const State& state) const -> int override;
    auto realizeAccelerationImpl(const State& state) const -> int override;

    // qdot==u here so these are just copies
    void multiplyByNImpl(const State& state, const Vector& u, Vector& dq) const override {
        dq = u;
    }
    void multiplyByNTransposeImpl(const State& state, const Vector& fq, Vector& fu) const override {
        fu = fq;
    }
    void multiplyByNPInvImpl(const State& state, const Vector& dq, Vector& u) const override {
        u = dq;
    }
    void multiplyByNPInvTransposeImpl(const State& state, const Vector& fu, Vector& fq) const override {
        fq = fu;
    }

    // No prescribed motion.
    auto prescribeQImpl(State& state) const -> bool override {
        return false;
    }
    auto prescribeUImpl(State& state) const -> bool override {
        return false;
    }

    void projectQImpl(State& state,
                      Vector& qErrEst,
                      const ProjectOptions& opts,
                      ProjectResults& results) const override;
    void projectUImpl(State& state,
                      Vector& uErrEst,
                      const ProjectOptions& opts,
                      ProjectResults& results) const override;


    ////////////////////////////////////////////////
    // Implementation of discrete System virtuals //
    ////////////////////////////////////////////////

    auto calcEventTriggerInfoImpl(const State& state, Array_<EventTriggerInfo>& eti) const -> int override {
        eti.clear();
        eti.push_back(
            EventTriggerInfo(eventId0).setRequiredLocalizationTimeWindow(1).setTriggerOnRisingSignTransition(
                false));
        eti.push_back(EventTriggerInfo(eventId1));
        eti.push_back(EventTriggerInfo(eventId2).setTriggerOnFallingSignTransition(false));
        return 0;
    }

    auto calcTimeOfNextScheduledEventImpl(const State& state,
                                          Real& tNextEvent,
                                          Array_<EventId>& eventIds,
                                          bool includeCurrentTime) const -> int override {
        // Generate an event every 5.123 seconds.
        int nFives = (int)(state.getTime() / 5.123); // rounded down
        if (state.getTime() == 0) {
            nFives = 1; // don't start with the event
        }

        tNextEvent = nFives * Real(5.123);
        // Careful ...
        if (tNextEvent < state.getTime() || (tNextEvent == state.getTime() && !includeCurrentTime)) {
            tNextEvent += Real(5.123);
        }
        eventIds.push_back(eventId1); // event Id for scheduled pulse

        return 0;
    }

    // This should be called when the integrator returns indicating that
    // a discontinuity (event trigger) has been detected. The current
    // state is inconsistent in some way and we expect the event handlers
    // to correct that. Time will be the same before and after, but the
    // state may have changed discontinuously.
    void handleEventsImpl(State& state,
                          Event::Cause cause,
                          const Array_<EventId>& eventIds,
                          const HandleEventsOptions& options,
                          HandleEventsResults& results) const override {
        switch (cause) {
            case Event::Cause::Initialization: {
                E_start = calcTotalEnergy(state);
                break;
            }
            case Event::Cause::TimeAdvanced: {
                const auto curTotalEnergy = calcTotalEnergy(state);
                E_min = std::min(E_min, curTotalEnergy);
                E_max = std::max(E_max, curTotalEnergy);
                if (!std::isfinite(E_start)) {
                    E_start = curTotalEnergy;
                }
                break;
            }
            case Event::Cause::Triggered:
            case Event::Cause::Scheduled:
            case Event::Cause::Signaled:
            case Event::Cause::Termination:
            case Event::Cause::Invalid: {
                // Check energy conservation across the discontinuity
                E_final = calcTotalEnergy(state);
                const bool validEnergies = std::isfinite(E_start) && std::isfinite(E_final)
                                           && std::isfinite(E_min) && std::isfinite(E_max);

                if (validEnergies && !firstSample) {
                    const SimTK::Real E_ref = 0.5 * (E_max + E_min);
                    const SimTK::Real denom = std::max(std::abs(E_ref), 1e-6);

                    // Drift: endpoint deviation from the empirical center (measures bias within interval)
                    const SimTK::Real E_center = 0.5 * (E_max + E_min);
                    SimTK::Real drift =
                        std::max(std::abs(E_start - E_ref), std::abs(E_final - E_ref)) / denom;

                    // Amp: half-range of oscillation relative to center (measures energy fluctuation scale)
                    SimTK::Real amp = (E_max - E_min) / (2 * denom);

                    // Fallback for near-zero energy regimes where relative scaling is ill-conditioned
                    if (std::abs(E_ref) < 1e-6) {
                        amp = 0.5 * (E_max - E_min);
                        drift = std::max(std::abs(E_start - E_ref), std::abs(E_final - E_ref));
                    }

                    EXPECT_LT(drift, tolerance)
                        << "Energy drift too high at discontinuous event! E_start=" << E_start
                        << " E_final=" << E_final << ", E_min=" << E_min << " E_max=" << E_max;
                    EXPECT_LT(amp, tolerance)
                        << "Energy oscillation amplitude too high at discontinuous event! E_start=" << E_start
                        << " E_final=" << E_final << ", E_min=" << E_min << " E_max=" << E_max;

                    // std::cout << "drift=" << drift << " amp=" << amp << "\n";
                }

                // Triggered, Initialization etc - discontinuous events
                std::swap(state.updQ()[0], state.updQ()[1]); // invalidates Position stage
                state.updU() = 0;

                // Reset
                E_start = SimTK::NaN;
                E_min = std::numeric_limits<SimTK::Real>::max();
                E_max = std::numeric_limits<SimTK::Real>::lowest();
                firstSample = false;

                break;
            }
            default:
                FAIL() << "Unexpected event cause: " << static_cast<int>(cause);
        }

        results.setExitStatus(HandleEventsResults::Succeeded);
    }

    [[nodiscard]] auto getMass(const State& state) const -> Real {
        const AbstractValue& m = state.getDiscreteVariable(subsysIndex, massIndex);
        return Value<Real>::downcast(m).get();
    }
    [[nodiscard]] auto getDefaultMass() const -> Real {
        return getMass(getDefaultState());
    }

    [[nodiscard]] auto getLength(const State& state) const -> Real {
        const AbstractValue& d = state.getDiscreteVariable(subsysIndex, lengthIndex);
        return Value<Real>::downcast(d).get();
    }
    [[nodiscard]] auto getDefaultLength() const -> Real {
        return getLength(getDefaultState());
    }

    [[nodiscard]] auto getGravity(const State& state) const -> Real {
        const AbstractValue& g = state.getDiscreteVariable(subsysIndex, gravityIndex);
        return Value<Real>::downcast(g).get();
    }
    [[nodiscard]] auto getDefaultGravity() const -> Real {
        return getGravity(getDefaultState());
    }

    [[nodiscard]] auto calcKineticEnergy(const State& state) const -> Real {
        const Real m = getMass(state);
        const Vector& u = state.getU(SubsystemIndex(0));

        return 0.5 * m * (u[0] * u[0] + u[1] * u[1]);
    }

    [[nodiscard]] auto calcPotentialEnergy(const State& state) const -> Real {
        const Real m = getMass(state);
        const Real g = getGravity(state);
        const Vector& q = state.getQ(SubsystemIndex(0));

        const Real y = q[1];
        return m * g * y;
    }

    [[nodiscard]] auto calcTotalEnergy(const State& s) const -> Real {
        return calcKineticEnergy(s) + calcPotentialEnergy(s);
    }

    auto getNumQProj() const -> int {
        return qProj;
    }

    auto getNumQProjFailures() const -> int {
        return qProjFail;
    }

    auto getNumUProj() const -> int {
        return uProj;
    }

    auto getNumUProjFailures() const -> int {
        return uProjFail;
    }

    void setTolerance(Real tol) const {
        tolerance = tol;
    }
};

// This is the handle class for a MyPendulum System.
// It must not have any data members. Data, if needed, is
// in the corresponding "Guts" class.

class MyPendulum : public System {
    public:
    MyPendulum() {
        adoptSystemGuts(new MyPendulumGuts());
        DefaultSystemSubsystem defsub(*this);
        updGuts().subsysIndex = defsub.getMySubsystemIndex();

        setHasTimeAdvancedEvents(false);
        realizeTopology();
    }

    [[nodiscard]] auto getGuts() const -> const MyPendulumGuts& {
        return dynamic_cast<const MyPendulumGuts&>(getSystemGuts());
    }

    [[nodiscard]] auto updGuts() -> MyPendulumGuts& {
        return dynamic_cast<MyPendulumGuts&>(updSystemGuts());
    }

    // Instance variables are written to our defaultState.
    void setDefaultMass(Real mass) {
        const MyPendulumGuts& guts = getGuts();
        updDefaultState().updDiscreteVariable(guts.subsysIndex, guts.massIndex) = Value<Real>(mass);
    }

    void setDefaultLength(Real length) {
        const MyPendulumGuts& guts = getGuts();
        updDefaultState().updDiscreteVariable(guts.subsysIndex, guts.lengthIndex) = Value<Real>(length);
    }

    void setDefaultGravity(Real gravity) {
        const MyPendulumGuts& guts = getGuts();
        updDefaultState().updDiscreteVariable(guts.subsysIndex, guts.gravityIndex) = Value<Real>(gravity);
    }

    void setDefaultTimeAndState(Real t, const Vector& q, const Vector& u) {
        const MyPendulumGuts& guts = getGuts();
        updDefaultState().updU(guts.subsysIndex) = u;
        updDefaultState().updQ(guts.subsysIndex) = q;
        updDefaultState().updTime() = t;
    }

    [[nodiscard]] auto getMass(const State& s) const -> Real {
        return getGuts().getMass(s);
    }
    [[nodiscard]] auto getDefaultMass() const -> Real {
        return getGuts().getMass(getDefaultState());
    }

    [[nodiscard]] auto getLength(const State& s) const -> Real {
        return getGuts().getLength(s);
    }
    [[nodiscard]] auto getDefaultLength() const -> Real {
        return getGuts().getDefaultLength();
    }

    [[nodiscard]] auto getGravity(const State& s) const -> Real {
        return getGuts().getGravity(s);
    }
    [[nodiscard]] auto getDefaultGravity() const -> Real {
        return getGuts().getDefaultGravity();
    }

    [[nodiscard]] auto calcKineticEnergy(const State& s) const -> Real {
        return getGuts().calcKineticEnergy(s);
    }

    [[nodiscard]] auto calcPotentialEnergy(const State& s) const -> Real {
        return getGuts().calcPotentialEnergy(s);
    }

    [[nodiscard]] auto calcTotalEnergy(const State& s) const -> Real {
        return getGuts().calcTotalEnergy(s);
    }

    [[nodiscard]] auto getNumQProj() const -> int {
        return getGuts().getNumQProj();
    }

    [[nodiscard]] auto getNumQProjFailures() const -> int {
        return getGuts().getNumQProjFailures();
    }

    [[nodiscard]] auto getNumUProj() const -> int {
        return getGuts().getNumUProj();
    }

    [[nodiscard]] auto getNumUProjFailures() const -> int {
        return getGuts().getNumUProjFailures();
    }

    void checkState(const Integrator& integ) {
        const State& state = integ.getAdvancedState();
        realize(state);

        integ.getConstraintToleranceInUse();

        const bool isInterpolated = integ.isStateInterpolated();
        const Vector& q = state.getQ(SubsystemIndex(0));
        const Vector& u = state.getU(SubsystemIndex(0));
        const Real d = getLength(state);

        // Contextual information for failure logs
        SCOPED_TRACE(::testing::Message() << "Mode: " << (isInterpolated ? "INTERPOLATED" : "INTEGRATOR_STEP")
                                          << " Time: " << state.getTime());

        // Interpolation error in DAEs is typically O(h^k). We allow a slightly
        // wider margin (10x) for interpolated points because they haven't been
        // "projected" back to the constraint manifold yet.
        const auto baseTol = integ.getConstraintToleranceInUse();
        const Real constraintTol = isInterpolated ? baseTol * 10.0 : baseTol;

        // Position and velocity constraints
        const Real perr = q[0] * q[0] + q[1] * q[1] - d * d;
        const Real verr = q[0] * u[0] + q[1] * u[1];

        EXPECT_NEAR(perr, 0.0, constraintTol) << "Manifold Drift: Pendulum length error too high during "
                                              << (isInterpolated ? "interpolation." : "actual step.");

        EXPECT_NEAR(verr, 0.0, constraintTol);

        // If we are interpolated, the State's Error Estimate (YErr)
        // should technically be "unknown" or represent the interpolation error.
        // On an actual step, YErr must be strictly within the solver's accuracy.
        if (!isInterpolated) {
            for (int i = 0; i < state.getYErr().size(); ++i) {
                EXPECT_LE(std::abs(state.getYErr()[i]), integ.getAccuracyInUse() * 10)
                    << "Solver claimed success but local error estimate is too high at Y[" << i << "]";
            }
        }

        // If the state is NOT interpolated, the multipliers (lambda) should be
        // well-defined because the forces must balance to satisfy the constraints.
        // If interpolated, we check if the multipliers have "exploded," which
        // indicates a discontinuous interpolation of the force manifold.
        const auto& lambda = state.getMultipliers();
        for (int i = 0; i < lambda.size(); ++i) {
            if (!isInterpolated) {
                EXPECT_TRUE(std::isfinite(lambda[i])) << "Singularity detected at step!";
            } else {
                // We allow finite but perhaps noisy values during interpolation
                EXPECT_LT(std::abs(lambda[i]), 1e6) << "Interpolation hit a force singularity!";
            }
        }

        // For Actual Steps, UDotErr must be near zero (DAE index-1 satisfaction).
        // For Interpolated points, we accept that the second derivative is "dirty."
        if (!isInterpolated) {
            for (int i = 0; i < state.getUDotErr().size(); ++i) {
                EXPECT_NEAR(state.getUDotErr()[i], 0.0, constraintTol);
            }
        }
    }

    void setTolerance(Real tol) const {
        getGuts().setTolerance(tol);
    }

    void dump(const char* msg) const {
        const MyPendulumGuts& guts = getGuts();
        std::cout << std::string(msg) << ": MyPendulum default state dump:" << std::endl;
        std::cout << "  mass   =" << getDefaultMass() << std::endl;
        std::cout << "  length =" << getDefaultLength() << std::endl;
        std::cout << "  gravity=" << getDefaultGravity() << std::endl;
        std::cout << "  time=" << getDefaultState().getTime() << std::endl;
        std::cout << "  q=" << getDefaultState().getQ(guts.subsysIndex) << std::endl;
        std::cout << "  u=" << getDefaultState().getU(guts.subsysIndex) << std::endl;
    }
};

inline const MyPendulum& MyPendulumGuts::getMyPendulum() const {
    return static_cast<const MyPendulum&>(getSystem());
}

/*
 * This system is a 2d pendulum swinging in gravity. It is modeled as
 * a point mass free in the plane, plus a distance constraint to model
 * the rod.
 *
 *    y       | g               O
 *    ^       v                  \  d
 *    |                           \
 *    |                            * m
 *     ------> x
 *
 * Gravity acts in the y direction, the rod is length d, mass m, pivot
 * location is the ground origin (0,0).
 *
 * The DAE for a generic multibody system is:
 *       qdot = Qu
 *       M udot = f - ~A lambda
 *       A udot = b
 *       perr(t,q) = 0
 *       verr(t,q,u) = 0
 *
 * Let   r^2 = x^2  + y^2
 *       v^2 = x'^2 + y'^2
 * We will express the "rod length=d" constraint as
 *       (r^2 - d^2)/2 = 0    (perr)
 *           xx' + yy' = 0    (verr)
 *         xx'' + yy'' = -v^2 (aerr)
 *
 * So the matrix A = d perr/dq = [x y] and b = -v^2, and the
 * equations of motion are:
 *     [ m 0 x ] [ x'' ]   [  0  ]
 *     [ 0 m y ] [ y'' ] = [ -mg ]
 *     [ x y 0 ] [ L   ]   [-v^2 ]
 * where L (the Lagrange multiplier) is proportional to
 * the rod tension. You can solve this to get
 *     L   = (m*v^2 - mg*y)/(r^2)
 *     x'' = - x*L/m
 *     y'' = - y*L/m - g
 *
 */
int MyPendulumGuts::realizeTopologyImpl(State& s) const {
    // Instance variables mass, length, gravity
    massIndex = s.allocateDiscreteVariable(subsysIndex, Stage::Instance, new Value<Real>(1));
    lengthIndex = s.allocateDiscreteVariable(subsysIndex, Stage::Instance, new Value<Real>(1));
    gravityIndex =
        s.allocateDiscreteVariable(subsysIndex, Stage::Instance, new Value<Real>(13.7503716373294544));
    const Vector init(2, Real(0));
    q0 = s.allocateQ(subsysIndex, init);
    u0 = s.allocateU(subsysIndex, init);

    mgForceIndex = s.allocateCacheEntry(subsysIndex, Stage::Dynamics, new Value<Real>());
    System::Guts::realizeTopologyImpl(s);
    return 0;
}
int MyPendulumGuts::realizeModelImpl(State& s) const {
    System::Guts::realizeModelImpl(s);
    return 0;
}
int MyPendulumGuts::realizeInstanceImpl(const State& s) const {
    qerr0 = s.allocateQErr(subsysIndex, 1);
    uerr0 = s.allocateUErr(subsysIndex, 1);
    udoterr0 = s.allocateUDotErr(subsysIndex, 1); // and multiplier
    trigger0 = s.allocateEventTrigger(subsysIndex, Stage::Position, 3);
    eventId0 = getSystem().getDefaultSubsystem().createEventId(subsysIndex, s);
    eventId1 = getSystem().getDefaultSubsystem().createEventId(subsysIndex, s);
    eventId2 = getSystem().getDefaultSubsystem().createEventId(subsysIndex, s);
    System::Guts::realizeInstanceImpl(s);
    return 0;
}

auto MyPendulumGuts::realizePositionImpl(const State& state) const -> int {
    const Real d = getMyPendulum().getLength(state);
    const Vector& q = state.getQ(subsysIndex);

    // This is the perr() equation.
    state.updQErr(subsysIndex)[0] = (q[0] * q[0] + q[1] * q[1] - d * d) / 2;

    state.updEventTriggersByStage(subsysIndex, Stage::Position)[0] = (100 * q[0]) - q[1];

    // Make sure this boolean trigger *crosses* zero; it won't work right
    // if one end is actually zero. We'll use -.5 for false, .5 for true.
    state.updEventTriggersByStage(subsysIndex, Stage::Position)[1] =
        static_cast<double>(state.getTime() > /*1.49552*/ 1.49545 && state.getTime() < 12.28937) - 0.5;

    state.updEventTriggersByStage(subsysIndex, Stage::Position)[2] = state.getTime() - 1.495508;
    System::Guts::realizePositionImpl(state);
    return 0;
}

int MyPendulumGuts::realizeVelocityImpl(const State& s) const {
    const Vector& q = s.getQ(subsysIndex);
    const Vector& u = s.getU(subsysIndex);
    Vector& qdot = s.updQDot(subsysIndex);

    qdot[0] = u[0]; // qdot=u
    qdot[1] = u[1];

    // This is the verr() equation.
    s.updUErr(subsysIndex)[0] = q[0] * u[0] + q[1] * u[1];
    System::Guts::realizeVelocityImpl(s);
    return 0;
}

int MyPendulumGuts::realizeDynamicsImpl(const State& s) const {
    const Real m = getMyPendulum().getMass(s);
    const Real g = getMyPendulum().getGravity(s);

    Real& mg = Value<Real>::updDowncast(s.updCacheEntry(subsysIndex, mgForceIndex)).upd();
    // Calculate the force due to gravity.
    mg = m * g;
    System::Guts::realizeDynamicsImpl(s);
    return 0;
}

int MyPendulumGuts::realizeAccelerationImpl(const State& s) const {
    const Real m = getMyPendulum().getMass(s);
    const Real g = getMyPendulum().getGravity(s);
    // we're pretending we couldn't calculate this here!
    const Real mg = Value<Real>::updDowncast(s.updCacheEntry(subsysIndex, mgForceIndex)).get();

    const Vector& q = s.getQ(subsysIndex);
    const Vector& u = s.getU(subsysIndex);
    Vector& udot = s.updUDot(subsysIndex);
    Vector& qdotdot = s.updQDotDot(subsysIndex);

    const Real r2 = q[0] * q[0] + q[1] * q[1];
    const Real v2 = u[0] * u[0] + u[1] * u[1];
    const Real L = (m * v2 - mg * q[1]) / r2;
    udot[0] = -q[0] * L / m;
    udot[1] = -q[1] * L / m - g;
    qdotdot = udot; // N=identity for this problem
    s.updMultipliers(subsysIndex)[0] = L;
    s.updUDotErr(subsysIndex)[0] = q[0] * udot[0] + q[1] * udot[1] + v2;
    System::Guts::realizeAccelerationImpl(s);
    return 0;
}

/*
 * Here we want to remove any constraint errors from the current state,
 * and project out any component of the integrator's error estimate
 * perpendicular to the constraint manifold. We will do this sequentially
 * rather than handling position and velocity simultaneously.
 *
 * For this system we have P = d perr/dq = V = d verr/du = [x y].
 * Weighted, we have PW=tp*[x/wx y/wy] VW=tv*[x/wxd y/wyd].
 * With pinv(A)=~A*(A*~A)^-1, we have:
 *
 *    pinv(P)  = ~[            x             y] /  (    x ^2+     y ^2)
 *    pinv(PW) = ~(1/tp)*[(wx *wy ^2)*x (wx ^2*wy) *y] / ((wy *x)^2+(wx *y)^2)
 *    pinv(VW) = ~(1/tv)*[(wxd*wyd^2)*x (wxd^2*wyd)*y] / ((wyd*x)^2+(wxd*y)^2)
 *      (the latter assuming x,y already projected on position manifold)
 *
 * We want to solve
 *    |perr(q0 - dq)|_TRMS <= accuracy, such that dq=min_WLS(dq)
 *    PW(q0) dq = Tp * perr(q0); q = q0-dq
 * Then
 *    |verr(q,u0 - du)|_TRMS <= accuracy, du=min_WLS(du)
 *    VW(q) du = Tv * verr(q,u0); u = u0-du
 *
 *
 * To remove the corresponding error estimates:
 *    PW(q) qperp = PW(q) qerrest; qerrest -= qperp
 *    VW(q) uperp = VW(q) uerrest; uerrest -= uperp
 *
 *
 */
static Real wrms(const Vector& y, const Vector& w) {
    Real sumsq = 0;
    for (int i = 0; i < y.size(); ++i) {
        sumsq += square(y[i] * w[i]);
    }
    return std::sqrt(sumsq / y.size());
}

// qerrest is in/out
void MyPendulumGuts::projectQImpl(State& s,
                                  Vector& qerrest,
                                  const ProjectOptions& opts,
                                  ProjectResults& results) const

{
    const Real consAccuracy = opts.getRequiredAccuracy();
    const Real projLimit = opts.getProjectionLimit();
    const bool forceProj = opts.isOptionSet(ProjectOptions::ForceProjection);
    const Vector& uweights = s.getUWeights(subsysIndex);
    const Vector& ctols = s.getQErrWeights(subsysIndex);

    // Since qdot=u here we can use uweights directly as qweights.
    const Vec2& wq = Vec2::getAs(&uweights[0]);
    const Real& tp = ctols[0]; // inverse tolerances 1/ti

    const Vec2& q = Vec2::getAs(&s.getQ(subsysIndex)[0]); // set up aliases
    Real& ep = s.updQErr(subsysIndex)[0];                 // ep changes as we go

    results.setAnyChangeMade(false);

    // std::cout << "BEFORE wperr=" << tp*ep << std::endl;
    if (!forceProj && std::abs(tp * ep) <= consAccuracy) {
        results.setExitStatus(ProjectResults::Succeeded);
        return;
    }
    if (std::abs(tp * ep) > projLimit) {
        results.setProjectionLimitExceeded(true);
        results.setExitStatus(ProjectResults::FailedToConverge);
        ++qProjFail;
        return;
    }
    ++qProj;
    results.setAnyChangeMade(true);
    Real wqchg;
    do {
        // Position projection
        const Real r2 = ~q * q; // x^2+y^2
        const Real wqr2 = square(wq[1] * q[0]) + square(wq[0] * q[1]);
        const Row2 P(~q);
        const Row2 PW(tp * q[0] / wq[0], tp * q[1] / wq[1]);
        const Vec2 Pinv(q / r2);
        const Vec2 PWinv = Vec2(square(wq[1]) * wq[0] * q[0], square(wq[0]) * wq[1] * q[1]) / (tp * wqr2);
        const Vec2 dq = Pinv * (ep);        // std::cout << "dq=" << dq << std::endl;
        const Vec2 wdq = PWinv * (tp * ep); // std::cout << "wdq=" << wdq << std::endl;

        wqchg = std::sqrt(wdq.normSqr() / q.size()); // wrms norm

        s.updQ(subsysIndex)[0] -= wdq[0] / wq[0];
        s.updQ(subsysIndex)[1] -= wdq[1] / wq[1];
        realize(s, Stage::Position); // recalc QErr (ep)

        // std::cout << "AFTER q-=wdq/W wperr=" << tp*ep << " wqchg=" << wqchg << std::endl;
    } while (std::abs(tp * ep) > consAccuracy && wqchg >= 0.01 * consAccuracy);

    // std::cout << "...AFTER wperr=" << tp*ep << std::endl;

    // Now do error estimates.

    if (qerrest.size()) {
        Vec2& eq = Vec2::updAs(&qerrest[0]);

        // Recalc PW, PWInv:
        const Real wqr2 = square(wq[1] * q[0]) + square(wq[0] * q[1]);
        const Row2 PW = Row2(tp * q[0] / wq[0], tp * q[1] / wq[1]);
        const Vec2 PWinv = Vec2(wq[0] * square(wq[1]) * q[0], square(wq[0]) * wq[1] * q[1]) / (tp * wqr2);

        Vec2 qperp = PWinv * (PW * eq);

        // std::cout << "ERREST before=" << yerrest
        //      << " wrms=" << wrms(qerrest,qweights) << std::endl;
        // std::cout << "PW*eq=" << PW*eq << std::endl;
        eq -= qperp;

        // std::cout << "ERREST after=" << yerrest
        //      << " wrms=" << wrms(qerrest,qweights) << std::endl;
        // std::cout << "PW*eq=" << PW*eq << std::endl;
    }

    results.setExitStatus(ProjectResults::Succeeded);
}

void MyPendulumGuts::projectUImpl(State& s,
                                  Vector& uerrest,
                                  const ProjectOptions& opts,
                                  ProjectResults& results) const {
    const Real consAccuracy = opts.getRequiredAccuracy();
    const Real projLimit = opts.getProjectionLimit();
    const bool forceProj = opts.isOptionSet(ProjectOptions::ForceProjection);
    const Vector& uweights = s.getUWeights(subsysIndex);
    const Vector& ctols = s.getUErrWeights(subsysIndex);

    const Vec2& wu = Vec2::getAs(&uweights[0]);
    const Real& tv = ctols[0];

    const Vec2& q = Vec2::getAs(&s.getQ(subsysIndex)[0]); // set up aliases
    const Vec2& u = Vec2::getAs(&s.getU(subsysIndex)[0]);
    Real& ev = s.updUErr(subsysIndex)[0]; // ev changes as we go

    results.setAnyChangeMade(false);

    // std::cout << "BEFORE wverr=" << tv*ev << std::endl;
    if (!forceProj && std::abs(tv * ev) <= consAccuracy) {
        results.setExitStatus(ProjectResults::Succeeded);
        return;
    }
    if (std::abs(tv * ev) > projLimit) {
        results.setProjectionLimitExceeded(true);
        results.setExitStatus(ProjectResults::FailedToConverge);
        ++uProjFail;
        return;
    }

    ++uProj;
    results.setAnyChangeMade(true);

    // Do velocity projection at current values of q, which should have
    // been projected already.
    Real r2 = ~q * q; // x^2+y^2
    Real wur2 = square(wu[1] * q[0]) + square(wu[0] * q[1]);
    Row2 V(~q), VW(tv * q[0] / wu[0], tv * q[1] / wu[1]);
    Vec2 Vinv(q / r2);
    Vec2 VWinv = Vec2(square(wu[1]) * wu[0] * q[0], square(wu[0]) * wu[1] * q[1]) / (tv * wur2);
    realize(s, Stage::Velocity); // calculate UErr (ev)

    // std::cout << "BEFORE wverr=" << tv*ev << std::endl;
    Vec2 du = Vinv * (ev);        // std::cout << "du=" << du << std::endl;
    Vec2 wdu = VWinv * (tv * ev); // std::cout << "wdu=" << wdu << std::endl;

    s.updU(subsysIndex)[0] -= wdu[0] / wu[0];
    s.updU(subsysIndex)[1] -= wdu[1] / wu[1];

    realize(s, Stage::Velocity); // recalc UErr
    // std::cout << "AFTER u-=wdu wverr=" << tv*ev << std::endl;

    // std::cout << "...AFTER wverr=" << tv*ev << std::endl;

    // Now do error estimates.


    if (uerrest.size()) {
        Vec2& eu = Vec2::updAs(&uerrest[0]);
        Vec2 uperp = VWinv * (VW * eu);

        // std::cout << "ERREST before=" << uerrest
        //      << " wrms=" << wrms(uerrest,uweights) << std::endl;
        // std::cout << " VW*eu=" << VW*eu << std::endl;
        eu -= uperp;

        // std::cout << "ERREST after=" << yerrest
        //      << " wrms=" << wrms(uerrest,uweights) << std::endl;
        // std::cout << " VW*eu=" << VW*eu << std::endl;
    }

    results.setExitStatus(ProjectResults::Succeeded);
}