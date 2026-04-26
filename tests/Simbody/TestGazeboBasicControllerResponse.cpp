/* -------------------------------------------------------------------------- *
 *    Simbody(tm): Gazebo Basic Controller Response                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
 * Authors: Michael Sherman, John Hsu                                         *
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
regression test "PhysicsTest::TrikeyWheelResponse". From the Gazebo test:

    The trikey model has three wheels oriented in different directions.
    it caught a corner case in ODE where the inertia was being
    truncated unequally based on cartesian orientation of the link.
    Hence this test is added to ensure we get the same dynamics behavior
    regardless the orientation of the inertia matrices in world frame.

There is a controller on each wheel ("PID" but actually has only a P gain)
that is instructed to rotate each wheel to an arbitrary target position.

The test just requires that all three wheels behave identically; it doesn't
check whether they also behave correctly!
*/
#include <gtest/gtest.h>

#include "Gazebo.hpp"

// Control gain
const Real Kp1 = 10; // proportional to angle error

// Target angle
const Real Target1 = 1.3; // radians

const Real TriKeyRadius = 0.27;
const Real TriKeyHalfLength = 0.65 / 2;
const Real TriKeyMass = 32.7;
const Vec3 TriKeyCOM(-0.003048, 0.00254, 0.41415);
const Inertia TriKeyCentralInertia(1.747, 1.747, 1.192); // Central
const Inertia TriKeyInertia(TriKeyCentralInertia.shiftFromMassCenter(-TriKeyCOM, TriKeyMass));
const Vec3 TriKeyColor = Gray;

const Real WheelRadius = 0.101;
const Real WheelHalfLength = 0.025 / 2;
const Real WheelMass = 0.66725;
const Vec3 WheelCOM(0, 0, 0);
const Inertia WheelInertia(0.00160418,
                           0.00279814,
                           0.00160339, // xx,yy,zz
                           -3.598e-06,
                           1.6199e-05,
                           -4.1656e-05); // xy,xz,yz
const Vec3 WheelColors[] = {Red, Green, Blue};


const Real MaxStepSize = 1e-3; // 1 ms (1000 Hz)
const int DrawEveryN = 5;
const Real SimTime = 1.5;

// Make this a whole number of viz frames
const int NSteps = DrawEveryN * (int(SimTime / MaxStepSize / DrawEveryN + 0.5));

// Use this class to hold references into the Simbody system.
struct MyMultibodySystem {
    MyMultibodySystem(); // see below

    MultibodySystem m_system;
    SimbodyMatterSubsystem m_matter;
    GeneralForceSubsystem m_forces;
    Force::DiscreteForces m_discrete;
    MobilizedBody::Pin m_trikey_base;
    std::array<MobilizedBody::Pin, 3> m_wheel;
};

// Construct the multibody system. The dampers are built in here but the springs
// are applied during execution.
MyMultibodySystem::MyMultibodySystem()
    : m_matter(m_system)
    , m_forces(m_system)
    , m_discrete(m_forces, m_matter) {
    Force::Gravity(m_forces, m_matter, -ZAxis, 9.81);

    // Cylinder is along Z in Gazebo, Y in Simbody
    Rotation YtoZ(Pi / 2, XAxis);

    Body::Rigid triKeyInfo(MassProperties(TriKeyMass, Vec3(TriKeyCOM), TriKeyInertia));
    Body::Rigid wheelInfo(MassProperties(WheelMass, Vec3(WheelCOM), WheelInertia));

    MobilizedBody& Ground = m_matter.updGround(); // Nicer name for Ground.

    m_trikey_base = MobilizedBody::Pin(Ground, Vec3(0, 0, .5), triKeyInfo, Vec3(0));

    const Rotation ZtoMinusY(Pi / 2, XAxis);
    for (int i = 0; i < 3; ++i) {
        std::array<Real, 3> offs = {0., 1e-1, -1e-1};
        const Real angle = (Real(i) * 2 * Pi / 3) + offs[i]; // 0, 120, 240
        Rotation aboutZ(angle, ZAxis);
        const Transform X_IF(aboutZ * ZtoMinusY, aboutZ * Vec3(0, -.24, .1) + Vec3(offs[i], 0, 0));
        const Transform X_OM(YtoZ, Vec3(0));

        m_wheel[i] = MobilizedBody::Pin(m_trikey_base, X_IF, wheelInfo, X_OM);
    }

    m_system.realizeTopology();
}

auto calcTotalEnergy(const MyMultibodySystem& mbs, const State& state) -> Real {
    mbs.m_system.realize(state, Stage::Dynamics);
    // Calculate potential energy in controller.
    Real controllerPE = 0;

    for (int i = 0; i < 3; ++i) {
        const Real aerr = mbs.m_wheel[i].getAngle(state) - Target1;
        controllerPE += Kp1 * square(aerr) / 2; // 1/2 k x^2
    }
    return mbs.m_system.calcEnergy(state) + controllerPE;
}

// Run the system until it settles down, then check the answers.
void runOnce(const MyMultibodySystem& mbs, Integrator& integ) {
    for (int stepNum = 1; stepNum <= NSteps; ++stepNum) {
        // Get access to State being advanced by the integrator. Interpolation
        // must be off so that we're modifying the actual trajectory.
        State& state = integ.updAdvancedState();

        // Apply controller forces.
        for (int i = 0; i < 3; ++i) {
            const Real aerr = mbs.m_wheel[i].getAngle(state) - Target1;
            mbs.m_discrete.setOneMobilityForce(state, mbs.m_wheel[i], MobilizerUIndex(0), -Kp1 * aerr);
        }

        // Advance time by MaxStepSize. Might take multiple internal steps to
        // get there, depending on difficulty and required accuracy.
        const Real tNext = stepNum * MaxStepSize;
        do {
            EXPECT_NO_THROW(integ.stepTo(tNext, tNext));
        } while (integ.getTime() < tNext);
    }

    const State& state = integ.getAdvancedState();
    mbs.m_system.realize(state);

    // These should be very similar since they are all treated the same.
    // They might not be right, but they should match!
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(mbs.m_wheel[0].getAngle(state), mbs.m_wheel[1].getAngle(state), 1e-12);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(mbs.m_wheel[0].getAngle(state), mbs.m_wheel[2].getAngle(state), 1e-12);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(mbs.m_wheel[0].getRate(state), mbs.m_wheel[1].getRate(state), 1e-12);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(mbs.m_wheel[0].getRate(state), mbs.m_wheel[2].getRate(state), 1e-12);
}

TEST(Simbody_GazeboBasicControllerResponse, LowAccuracy) {
    MyMultibodySystem mbs;

    SemiExplicitEuler2Integrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(1e-3);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ);
}

TEST(Simbody_GazeboBasicControllerResponse, HighAccuracy) {
    MyMultibodySystem mbs;

    RungeKuttaMersonIntegrator rkm(mbs.m_system);
    rkm.setAllowInterpolation(false);
    rkm.setAccuracy(1e-6);
    rkm.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, rkm);
}
