/* -------------------------------------------------------------------------- *
 *                     Simbody(tm): Gazebo Reaction Force                     *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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
regression test "Joint_TEST::ForceTorque".

It is a small sphere pinned at its center to ground, rotating about X, then
an outboard brick pinned to the sphere and rotating about mutual Z. Oblique
gravity cause the brick to fall, hit the ground and come to rest. Then we check
against John Hsu's hand coded reaction force values.
*/

#include "Gazebo.hpp"

// Control gains
const Real Kp1 = 50000; // link1-2 joint stiffness
const Real Kp2 = 10000; // link2-3 joint stiffness
const Real Cd1 = 100;   // link1-2 joint damping
const Real Cd2 = 30;    // link2-3 joint damping

// Target angles
const Real Target1 = 0;
const Real Target2 = -Pi / 4;

const Real Mass1 = 100;
const Real Mass2 = 5;
const Real Mass3 = 1;

const Vec3 Cube(.5, .5, .5); // half-dimensions of cube
const Vec3 Box(.05, .1, .2); // half-dimensions of box
const Vec3 Centroid(.5, 0, .5);
const Vec3 COM1 = Centroid;
const Vec3 COM2 = Centroid;
const Vec3 COM3 = Centroid + Vec3(0, .5, 0);

// Simbody requires inertias to be expressed about body origin rather than COM.
const Inertia Inertia1 = Inertia(Vec3(1)).shiftFromMassCenter(-COM1, Mass1);
const Inertia Inertia2 = Inertia(Vec3(.05)).shiftFromMassCenter(-COM2, Mass2);
const Inertia Inertia3 = Inertia(Vec3(.001, .001, 0)).shiftFromMassCenter(-COM3, Mass3);

// Define a stiff, lossy material.
const Real Stiffness = 1e8;
const Real Dissipation = 10;
const Real Mu_s = 0.15;
const Real Mu_d = 0.1;
const Real Mu_v = 0;
const ContactMaterial lossyMaterial(Stiffness, Dissipation, Mu_s, Mu_d, Mu_v);

constexpr Real MaxStepSize = 1e-3; // 1 ms (1000 Hz)
constexpr int DrawEveryN = 33;     // 33 ms frame update (30.3 Hz)
constexpr Real SimTime = 5;

// Make this a whole number of viz frames
constexpr int NSteps = DrawEveryN * (int((SimTime / MaxStepSize / DrawEveryN) + 0.5));

// Use this class to hold references into the Simbody system.
struct MyMultibodySystem {
    MyMultibodySystem();

    MultibodySystem m_system;
    SimbodyMatterSubsystem m_matter;
    ContactTrackerSubsystem m_tracker;
    CompliantContactSubsystem m_contact;
    GeneralForceSubsystem m_forces;
    Force::Gravity m_gravity;
    MobilizedBody::Pin m_link1, m_link2;
};

MyMultibodySystem::MyMultibodySystem()
    : m_system()
    , m_matter(m_system)
    , m_tracker(m_system)
    , m_contact(m_system, m_tracker)
    , m_forces(m_system)
    , m_gravity(m_forces, m_matter, Vec3(-30, 10, -50)) {
    Body::Rigid link1Info(MassProperties(Mass1, COM1, Inertia1));
    Body::Rigid link2Info(MassProperties(Mass2, COM2, Inertia2));

    ContactGeometry::TriangleMesh cubeMesh(PolygonalMesh::createBrickMesh(Box, 3));
    link2Info.addContactSurface(COM2, ContactSurface(cubeMesh, lossyMaterial, 1));

    MobilizedBody& Ground = m_matter.updGround(); // Nicer name for Ground.

    // Add the Ground contact geometry. Contact half space has -XAxis normal
    // (right hand wall) so we have to rotate.
    const Rotation NegXToZ(Pi / 2, YAxis);
    Ground.updBody().addContactSurface(Transform(NegXToZ, Vec3(0)),
                                       ContactSurface(ContactGeometry::HalfSpace(), lossyMaterial));

    const Rotation ZtoX(Pi / 2, YAxis);
    m_link1 = MobilizedBody::Pin(Ground, Transform(ZtoX, COM1), link1Info, Transform(ZtoX, COM1));

    m_link2 = MobilizedBody::Pin(m_link1, Vec3(0, 0, 1.5), link2Info, Vec3(0));

    Force::MobilityLinearStop(m_forces, m_link1, MobilizerQIndex(0), 1000000., 1., -Pi / 2, Pi / 2);
    Force::MobilityLinearStop(m_forces, m_link2, MobilizerQIndex(0), 1000000., 1., -Infinity, Infinity);

    m_system.realizeTopology();
}

// Run the system until it settles down, then check the answers.
template <typename IntegratorType>
void runOnce(const MyMultibodySystem& mbs, IntegratorType& integ) {
    int stepNum = 0;
    while (stepNum++ < NSteps) {
        // Get access to State being advanced by the integrator. Interpolation must be off so that we're
        // modifying the actual trajectory.
        State& state = integ.updAdvancedState();

        // Advance time by MaxStepSize. Might take multiple internal steps to get there, depending on
        // difficulty and required accuracy.
        const Real tNext = stepNum * MaxStepSize;
        do {
            EXPECT_NO_THROW(integ.stepTo(tNext, tNext));
        } while (integ.getTime() < tNext);
    }

    const State& state = integ.getAdvancedState();
    mbs.m_system.realize(state);

    const ReactionPair reaction1 = getReactionPair(state, mbs.m_link1);
    const ReactionPair reaction2 = getReactionPair(state, mbs.m_link2);

    // Check the answers
    // Note (torque,force) ordering.
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(reaction1.reactionOnParentInParent,
                                 SpatialVec(Vec3(-750, 0, 450), Vec3(-600, 200, -1000)),
                                 0.5);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(reaction1.reactionOnChildInChild,
                                 SpatialVec(Vec3(750, 450, 0), Vec3(600, -1000, -200)),
                                 0.5);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(reaction2.reactionOnParentInParent,
                                 SpatialVec(Vec3(-250, -150, 0), Vec3(-300, 500, 100)),
                                 0.5);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(reaction2.reactionOnChildInChild,
                                 SpatialVec(Vec3(250, 150, 0), Vec3(300, -500, -100)),
                                 0.5);
}

TEST(Simbody_GazeboReactionForce, LowAccuracy) {
    MyMultibodySystem mbs;

    SemiExplicitEuler2Integrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(0.01);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ);
}

TEST(Simbody_GazeboReactionForce, HighAccuracy) {
    MyMultibodySystem mbs;

    RungeKuttaMersonIntegrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(1e-6);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ);
}
