/* -------------------------------------------------------------------------- *
 *      Simbody(tm): Gazebo Reaction Force With Applied Force  (Rigid)        *
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
regression test "Joint_TEST::GetForceTorqueWithAppliedForce". Here we are
using Simbody's unilateral rigid contact which should behave very similar to
the Gazebo reference.

See GazeboReactionForceWithAppliedForceCompliant.cpp for the same problem done
using compliant contact. (The rigid test case here was modified from the
compliant one so may have some dead code left over from that version.)

It is a stack of three cubes hinged together at their edges. The bottom block
is heavy and rests on the ground, the other two are light and have their
positions maintained by a pair of PD controllers like this:

              /
     link3  /   \         * = pin joint
           /     \
           \  1  /
             \  /                   z
        ------*  45 degrees         ^                     g = 0 0 -50
  link2 |     |                     |   y                       |
        |  6  |                     |  /                        |
   -----*------                     | /                         v
  |     | 0 degrees                 ---------> x
  | 100 |
  ------- link1
  contact
   GROUND

All the cube edges are of length 1. Masses are 100,6,1 as shown. The top block
has a COM that is offset into the +y direction by 0.5.

Expected reaction force results:
    pin1 on inboard:  tx=-25, ty= 175, fz=-300
    pin1 on outboard: tx= 25, ty=-175, fx=1, fz= 300,
    pin2 on inboard:  tx=-25, fz=-50
    pin2 in outboard: tx=25*pi/4, tz=-25*pi/4, fx=fz=50*pi/4

The reason for fx=1 in the second line is that with a gain of only 50000, the
first pin joint must be at an angle of 0.0035 radians to generate a torque of
175. That tips link2's frame in which (0,0,300)_G is reexpressed.
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
    MyMultibodySystem(); // see below

    MultibodySystem m_system;
    SimbodyMatterSubsystem m_matter;
    ContactTrackerSubsystem m_tracker;
    CompliantContactSubsystem m_contact;
    GeneralForceSubsystem m_forces;
    Force::DiscreteForces m_discrete;
    MobilizedBody m_link1;
    MobilizedBody::Pin m_link2, m_link3;
};

// Execute the simulation with a given integrator and accuracy, and verify that
// it produces the correct answers.
static void runOnce(const MyMultibodySystem& mbs, Integrator& integ, Real accuracy);

// Construct the multibody system. The dampers are built in here but the springs
// are applied during execution.
MyMultibodySystem::MyMultibodySystem()
    : m_matter(m_system)
    , m_tracker(m_system)
    , m_contact(m_system, m_tracker)
    , m_forces(m_system)
    , m_discrete(m_forces, m_matter) {
    Force::Gravity(m_forces, m_matter, -ZAxis, 50);

    Body::Rigid link1Info(MassProperties(Mass1, COM1, Inertia1));
    ContactGeometry::TriangleMesh cubeMesh(PolygonalMesh::createBrickMesh(Cube, 3));
    link1Info.addContactSurface(Centroid, ContactSurface(cubeMesh, lossyMaterial, 1));

    Body::Rigid link2Info(MassProperties(Mass2, COM2, Inertia2));
    Body::Rigid link3Info(MassProperties(Mass3, COM3, Inertia3));

    MobilizedBody& Ground = m_matter.updGround();

    // Add the Ground contact geometry. Contact half space has -XAxis normal (right hand wall) so we have to
    // rotate.
    const Rotation NegXToZ(Pi / 2, YAxis);
    Ground.updBody().addContactSurface(Transform(NegXToZ, Vec3(0)),
                                       ContactSurface(ContactGeometry::HalfSpace(), lossyMaterial));
    m_link1 = MobilizedBody::Free(Ground, Vec3(0), link1Info, Vec3(0));

    const double CoefRest = 0;
    for (int i = -1; i <= 1; i += 2) {
        for (int j = -1; j <= 1; j += 2) {
            for (int k = -1; k <= 1; k += 2) {
                const Vec3 point = Centroid + Vec3(i, j, k).elementwiseMultiply(Cube);
                auto* contact =
                    new PointPlaneContact(Ground, ZAxis, 0., m_link1, point, CoefRest, Mu_s, Mu_d, Mu_v);
                m_matter.adoptUnilateralContact(contact);
            }
        }
    }

    // Use this instead of the free joint to remove contact.
    // m_link1 = MobilizedBody::Weld(Ground,Vec3(0),
    //                              link1Info, Vec3(0));
    const Rotation ZtoY(-Pi / 2, XAxis);
    m_link2 = MobilizedBody::Pin(m_link1, Transform(ZtoY, 2 * Centroid), link2Info, Transform(ZtoY, Vec3(0)));
    m_link3 = MobilizedBody::Pin(m_link2, Transform(ZtoY, 2 * Centroid), link3Info, Transform(ZtoY, Vec3(0)));

    // It is more stable to build the springs into the mechanism like
    // this rather than apply them discretely.
    // Force::MobilityLinearSpring
    //   (m_forces, m_link2, MobilizerQIndex(0), Kp1, Target1);
    // Force::MobilityLinearSpring
    //   (m_forces, m_link3, MobilizerQIndex(0), Kp2, Target2);
    Force::MobilityLinearDamper(m_forces, m_link2, MobilizerUIndex(0), Cd1);
    Force::MobilityLinearDamper(m_forces, m_link3, MobilizerUIndex(0), Cd2);

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

        // Apply discrete spring forces.
        const Real a1err = mbs.m_link2.getAngle(state) - Target1;
        mbs.m_discrete.setOneMobilityForce(state, mbs.m_link2, MobilizerUIndex(0), -Kp1 * a1err);

        const Real a2err = mbs.m_link3.getAngle(state) - Target2;
        mbs.m_discrete.setOneMobilityForce(state, mbs.m_link3, MobilizerUIndex(0), -Kp2 * a2err);

        // Advance time by MaxStepSize. Might take multiple internal steps to get there, depending on
        // difficulty and required accuracy.
        const Real tNext = stepNum * MaxStepSize;
        do {
            EXPECT_NO_THROW(integ.stepTo(tNext));
        } while (integ.getTime() < tNext);
    }

    const State& state = integ.getAdvancedState();
    mbs.m_system.realize(state);

    const ReactionPair reaction2 = getReactionPair(state, mbs.m_link2);
    const ReactionPair reaction3 = getReactionPair(state, mbs.m_link3);

    // Check the answers
    // Note (torque,force) ordering.
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(reaction2.reactionOnParentInParent,
                                 SpatialVec(Vec3(-25, 175, 0), Vec3(0, 0, -300)),
                                 0.5);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(reaction2.reactionOnChildInChild,
                                 SpatialVec(Vec3(25, -175, 0), Vec3(-1, 0, 300)),
                                 0.5);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(reaction3.reactionOnParentInParent,
                                 SpatialVec(Vec3(-25, 0, 0), Vec3(0, 0, -50)),
                                 0.5);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(reaction3.reactionOnChildInChild,
                                 (Pi / 4) * SpatialVec(Vec3(25, 0, -25), Vec3(50, 0, 50)),
                                 0.5);
}

TEST(Simbody_GazeboReactionForces_Rigid, LowAccuracy) {
    MyMultibodySystem mbs;

    SemiExplicitEulerTimeStepper integ(mbs.m_system);
    integ.setPositionProjectionMethod(SemiExplicitEulerTimeStepper::Bilateral);
    integ.setAccuracy(1e-2);
    integ.setConstraintTolerance(.001);

    ImpulseSolver* solver = new PLUSImpulseSolver(integ.getDefaultFrictionTransitionVelocityInUse());
    integ.setImpulseSolver(solver);

    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ);
}

TEST(Simbody_GazeboReactionForces_Compliant, LowAccuracy) {
    MyMultibodySystem mbs;

    SemiExplicitEuler2Integrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(1e-2);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ);
}

TEST(Simbody_GazeboReactionForces_Compliant, HighAccuracy) {
    MyMultibodySystem mbs;

    RungeKuttaMersonIntegrator integ(mbs.m_system);
    integ.setAllowInterpolation(false);
    integ.setAccuracy(1e-6);
    integ.initialize(mbs.m_system.getDefaultState());

    runOnce(mbs, integ);
}
