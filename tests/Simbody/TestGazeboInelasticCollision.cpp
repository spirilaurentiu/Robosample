/* -------------------------------------------------------------------------- *
 *                   Simbody(tm): Gazebo Inelastic Collision                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

/* This test is drawn from the Open Source Robotics Foundation Gazebo physics
regression test "CollisionTest". In the original test there is a cube and
a sphere, both of unit mass and same half-dimension 0.5. Initially they are 1
unit apart in x, with the block to the left of the sphere. The block is
accelerated to the right with a discrete applied force of 1000 applied for 1ms
for a total impulse of 1, which should produce a velocity of +1 that is steady
until the block hits the stationary sphere. The collision is supposed to be
fully inelastic (coefficient of restitution=0), so there should be no rebound
and the two objects should move off to the right together at half the impact
velocity.

The test runs with fixed 1ms steps and calculates the expected result by
integrating manually. The position and velocity changes during the first step
with the 1000 unit force active are
   x(.001) = 1/2 a t^2 = 1/2 1e3 1e-6 = .0005 length units
   v(.001) = a t       = 1e3 1e-3     = 1 velocity unit
After that first step, all subsequent steps prior to the collision have
   deltaX = v t       = 1 1e-3       = .001 length units
   deltaV = 0
We would like the collision to occur exactly at 1s, that is, between steps 1000
and 1001. The position x(1)=.0005 + 999*.001=.9995. So if we were going to apply
the entire impulse during the first step we would place the sphere initially so
that the block surface and sphere surface are separated by .9995. But read on.

NOTE: to avoid problems with first-order integrators (explanation below), we're
modifying the test here to apply 1/10th of the impulse over the first 10 steps
(.01 seconds) rather than all in the first step. Then
   x(.01) = 1/2 a t^2 = 1/2 1e2 1e-4 = .005 length units
   v(.01) = a t       = 1e2 1e-2     = 1 velocity unit
At 1s we'll have x(1)=.005 + 990*.001=.995. So we'll set the initial separation
to .995 to get a collision at 1s (between steps 1000 and 1001).

Why first order integrators cause trouble
-----------------------------------------
While any integrator can get the velocity correct in one step, no first order
integrator will be able to calculate x(.001) correctly in a single step.
  Explicit Euler:           v(.001)=v(0) + h*a(0)=1
                            x(.001)=x(0) + h*v(0)=0                 <--

  Semi-explicit Euler:      v(.001)=v(0) + h*a(0)=1
                            x(.001)=x(0) + h*v(.001)=.001.          <--

  Semi-explicit Euler 2:    v'(.0005)=v(0)+h/2 a(0)=.5
                            x'(.0005)=x(0)+h/2 v'(.0005)=.00025
                            v(.001)=v'(.0005)+h/2 a'(.0005)=1
                            x(.001)=x'(.0005)+h/2 v(.001)=.00075    <--

But any 2nd order integrator would give the correct result, for example
Explicit trapezoid rule:    v'(.001)=v(0)+h*a(0)=1
                            x(.001)=x(0)+h*(v(0)+v'(.001))/2 =.0005 <--

Explicit midpoint rule:     v'(.0005)=v(0)+h/2*a(0)=1/2
                            x(.001)=x(0)+h*v'(.0005)=.0005          <--

By spreading the impulse over several steps we substantially reduce the overall
error. For example, after ten steps we get
  ExplicitEuler:         x(.01)=.0045
  Semi-explicit Euler:   x(.01)=.0055
  Semi-explicit Euler 2: x(.01)=.00525
any of which we'll deem close enough to the right answer of .005.
*/

#include <gtest/gtest.h>

#include "Gazebo.hpp"

using namespace SimTK;

const Real Radius = 0.5;
const Real Mass = 1;

// Define an extremely stiff, lossy material.
const Real Stiffness = 1e8;
const Real Dissipation = 1000;

const Real MaxStepSize = Real(1 / 1000.); // 1 ms (1000 Hz)
const int DrawEveryN = 33;                // 33 ms frame update (30.3 Hz)
const Real SimTime = 2;

// Make this a whole number of viz frames
const int NSteps = DrawEveryN * (int(SimTime / MaxStepSize / DrawEveryN + 0.5));

const Real TotalImpulse = 1; // Applied only on the first step
const Real StepForce = TotalImpulse / MaxStepSize;

const Real IntegAccuracy = 1e-3;
const Real CheckAccuracy = 1e-3;

struct MyMultibodySystem {
    MyMultibodySystem()
        : m_matter(m_system)
        , m_tracker(m_system)
        , m_contact(m_system, m_tracker)
        , m_forces(m_system)
        , m_discrete(m_forces, m_matter) {
        ContactMaterial lossyMaterial(Stiffness,
                                      Dissipation,
                                      0,  // mu_static
                                      0,  // mu_dynamic
                                      0); // mu_viscous

        // no gravity
        Body::Rigid sphereBody(MassProperties(Mass, Vec3(0), UnitInertia::sphere(Radius)));
        sphereBody.addDecoration(Transform(), DecorativeSphere(Radius));
        sphereBody.addContactSurface(Transform(),
                                     ContactSurface(ContactGeometry::Sphere(Radius), lossyMaterial));

        // TODO: using inscribed sphere as contact shape within the block.
        Body::Rigid cubeBody(MassProperties(Mass, Vec3(0), UnitInertia::sphere(Radius)));
        cubeBody.addDecoration(Transform(), DecorativeBrick(Vec3(Radius)));
        cubeBody.addContactSurface(Transform(),
                                   ContactSurface(ContactGeometry::Sphere(Radius), // TODO!
                                                  lossyMaterial));

        MobilizedBody Ground = m_matter.Ground();

        m_cube = MobilizedBody::Slider(Ground, Transform(Vec3(0, 2, 0)), cubeBody, Transform(Vec3(0)));
        m_sphere =
            MobilizedBody::Slider(Ground, Transform(Vec3(2 - .005, 2, 0)), sphereBody, Transform(Vec3(0)));

        m_system.realizeTopology();
    }

    MultibodySystem m_system;
    SimbodyMatterSubsystem m_matter;
    ContactTrackerSubsystem m_tracker;
    CompliantContactSubsystem m_contact;
    GeneralForceSubsystem m_forces;
    Force::DiscreteForces m_discrete;
    MobilizedBody m_cube;
    MobilizedBody m_sphere;
};

void runOnce(const MyMultibodySystem& mbs, Integrator& integ, int nsteps) {
    // These variables are the manually calculated values for the cube's
    // x coordinate and x velocity in Ground. We'll calculate these assuming
    // a perfect inelastic collision occurring at t=1 and then compare with
    // the approximate solution produced by the integrator.
    Real x = 0;
    Real v = 0;

    for (int stepNum = 1; stepNum <= NSteps; ++stepNum) {
        // Get access to State being advanced by the integrator. Interpolation
        // must be off so that we're modifying the actual trajectory.
        State& state = integ.updAdvancedState();

        mbs.m_discrete.clearAllBodyForces(state);
        if (stepNum <= 10) {
            mbs.m_discrete.setOneBodyForce(state,
                                           mbs.m_cube,
                                           SpatialVec(Vec3(0), Vec3(StepForce / 10, 0, 0)));
        }

        // Advance time by MaxStepSize. Might take multiple internal steps to
        // get there, depending on difficulty and required accuracy.
        const Real tNext = stepNum * MaxStepSize;
        do {
            EXPECT_NO_THROW(integ.stepTo(tNext, tNext));
        } while (integ.getTime() < tNext);

        // From Gazebo test code in physics.cc:
        // integrate here to see when the collision should happen
        if (stepNum <= 10) {
            const Real acceleration = StepForce / 10 / Mass;

            // dx = 1/2 a t^2
            x += (v * MaxStepSize) + (acceleration * square(MaxStepSize) / 2);

            // dv = a t
            v += acceleration * MaxStepSize;
        } else {
            const Real impulse = StepForce * MaxStepSize;

            if (stepNum >= 1000) {
                // inelastic col. w/equal mass
                v = impulse / (2 * Mass);
            }
            x += v * MaxStepSize;
        }

        mbs.m_system.realize(state);

        // Allow some uncertainty between step 1000 and 1001.
        if (stepNum != 1000) {
            EXPECT_NEAR_CUSTOM_TOL_SIMTK(mbs.m_cube.getBodyOriginLocation(state)[0], x, CheckAccuracy);
            EXPECT_NEAR_CUSTOM_TOL_SIMTK(mbs.m_cube.getBodyOriginVelocity(state)[0], v, CheckAccuracy);
        }
    }
}

TEST(Simbody_GazeboInelasticCollision, RungeKuttaMersonIntegrator) {
    SCOPED_TRACE("RungeKuttaMersonIntegrator");

    MyMultibodySystem mbs;
    RungeKuttaMersonIntegrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ, NSteps);
}

TEST(Simbody_GazeboInelasticCollision, RungeKutta3Integrator) {
    SCOPED_TRACE("RungeKutta3Integrator");

    MyMultibodySystem mbs;
    RungeKutta3Integrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ, NSteps);
}

TEST(Simbody_GazeboInelasticCollision, RungeKutta2Integrator) {
    SCOPED_TRACE("RungeKutta2Integrator");

    MyMultibodySystem mbs;
    RungeKutta2Integrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ, NSteps);
}

TEST(Simbody_GazeboInelasticCollision, SemiExplicitEuler2Integrator) {
    SCOPED_TRACE("SemiExplicitEuler2Integrator");

    MyMultibodySystem mbs;
    SemiExplicitEuler2Integrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ, NSteps);
}

TEST(Simbody_GazeboInelasticCollision, ExplicitEulerIntegrator) {
    SCOPED_TRACE("ExplicitEulerIntegrator");

    MyMultibodySystem mbs;
    ExplicitEulerIntegrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ, NSteps);
}

// Fixed-step integrator that can't handle the compliant impact which demands a few smaller steps for
// stability
TEST(Simbody_GazeboInelasticCollision, SemiExplicitEulerIntegratorFixedStep) {
    SCOPED_TRACE("SemiExplicitEulerIntegratorFixedStep");

    MyMultibodySystem mbs;
    SemiExplicitEulerIntegrator integ(mbs.m_system, MaxStepSize);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ, 990); // stop before impact
}

// Fixed-step integrator that can't handle the compliant impact which demands a few smaller steps for
// stability
TEST(Simbody_GazeboInelasticCollision, ExplicitEulerIntegratorFixedStep) {
    SCOPED_TRACE("ExplicitEulerIntegratorFixedStep");

    MyMultibodySystem mbs;
    ExplicitEulerIntegrator integ(mbs.m_system, MaxStepSize);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ, 990); // stop before impact
}

// sherm 130617: this isn't accurate enough; I think it has a serious bug in its error estimator
TEST(Simbody_GazeboInelasticCollision, VerletIntegrator) {
    SCOPED_TRACE("VerletIntegrator");

    MyMultibodySystem mbs;
    VerletIntegrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ, NSteps);
}

// This has poor accuracy control also
TEST(Simbody_GazeboInelasticCollision, RungeKuttaFeldbergIntegrator) {
    SCOPED_TRACE("RungeKuttaFeldbergIntegrator");

    MyMultibodySystem mbs;
    RungeKuttaFeldbergIntegrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(IntegAccuracy);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ, NSteps);
}
