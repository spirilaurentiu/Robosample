#include <gtest/gtest.h>

#include "Simbody.h"

using namespace SimTK;

#define EXPECT_SIMTK_SIZE(expected, actual, n)                            \
    EXPECT_TRUE(SimTK::Test::numericallyEqual((expected), (actual), (n))) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\nSize:      " << (n) << "\n"

#define EXPECT_NEAR_CUSTOM_TOL_SIMTK(expected, actual, tol)                              \
    EXPECT_TRUE(SimTK::Test::numericallyEqual((expected), (actual), /*scale*/ 1, (tol))) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\nTolerance: " << (tol) << "\n"

#define EXPECT_NEAR_DEFAULT_TOL_SIMTK(expected, actual) \
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((expected), (actual), 1e-6)

TEST(Simbody_MobilizedBody_KinematicsCalculationMethods,
     StationTransformsVelocitiesAndFrameConversionsConsistent) {
    // Create a system with two bodies.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Free b1(matter.Ground(), body);
    MobilizedBody::Free b2(matter.Ground(), body);

    // Set all the state variables to random values.
    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Gaussian random;

    for (int i = 0; i < state.getNY(); ++i) {
        state.updY()[i] = random.getValue();
    }

    system.realize(state, Stage::Acceleration);

    // Test the low level methods for transforming points and vectors.
    const Vec3 point(0.5, 1, -1.5);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.findStationLocationInGround(state, Vec3(0)),
                                  b1.getBodyOriginLocation(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(
        b1.findStationAtGroundPoint(state, b1.findStationLocationInGround(state, point)),
        point);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(
        b2.findStationAtGroundPoint(state, b1.findStationLocationInGround(state, point)),
        b1.findStationLocationInAnotherBody(state, point, b2));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(
        b2.findStationAtGroundPoint(state, b1.findStationLocationInGround(state, Vec3(0))).norm(),
        (b1.getBodyOriginLocation(state) - b2.getBodyOriginLocation(state)).norm());
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b2.findMassCenterLocationInGround(state),
                                  b2.findStationLocationInGround(state, b2.getBodyMassCenterStation(state)));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.expressVectorInGroundFrame(state, Vec3(0)), Vec3(0));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.expressVectorInGroundFrame(state, point),
                                  b1.getBodyRotation(state) * point);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(
        b1.expressGroundVectorInBodyFrame(state, b1.expressVectorInGroundFrame(state, point)),
        point);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(
        b2.expressGroundVectorInBodyFrame(state, b1.expressVectorInGroundFrame(state, point)),
        b1.expressVectorInAnotherBodyFrame(state, point, b2));

    // Test the routines for mapping locations, velocities, and accelerations.
    Vec3 r;
    Vec3 v;
    Vec3 a;
    b1.findStationLocationVelocityAndAccelerationInGround(state, point, r, v, a);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(v, b1.findStationVelocityInGround(state, point));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(a, b1.findStationAccelerationInGround(state, point));

    Vec3 r2;
    Vec3 v2;
    b1.findStationLocationAndVelocityInGround(state, point, r2, v2);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(r, r2);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(v, v2);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.findStationVelocityInGround(state, Vec3(0)),
                                  b1.getBodyOriginVelocity(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.findStationAccelerationInGround(state, Vec3(0)),
                                  b1.getBodyOriginAcceleration(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.findStationVelocityInGround(state, point),
                                  b1.findStationVelocityInAnotherBody(state, point, matter.Ground()));
}

TEST(Simbody_MobilizedBody_WeldMobilizerVsWeldConstraintEquivalence,
     WeldedSystemsRemainKinematicallyIdenticalDuringSimulation) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -1, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));

    // Create two pendulums, each with two welded bodies. One uses a Weld MobilizedBody,
    // and the other uses a Weld constraint.
    Transform inboard(Vec3(0.1, 0.5, -1));
    Transform outboard(Vec3(0.2, -0.2, 0));
    MobilizedBody::Ball p1(matter.updGround(), Vec3(0), body, Vec3(0, 1, 0));
    MobilizedBody::Ball p2(matter.updGround(), Vec3(0), body, Vec3(0, 1, 0));
    MobilizedBody::Weld c1(p1, inboard, body, outboard);
    MobilizedBody::Free c2(p2, inboard, body, outboard);
    Constraint::Weld constraint(p2, inboard, c2, outboard);

    // It is not a general test unless the Weld mobilizer has children!
    MobilizedBody::Pin wchild1(c1, inboard, body, outboard);
    MobilizedBody::Pin wchild2(c2, inboard, body, outboard);
    Force::MobilityLinearSpring(forces, wchild1, 0, 1000, 0);
    Force::MobilityLinearSpring(forces, wchild2, 0, 1000, 0);

    State state = system.realizeTopology();
    p1.setU(state, Vec3(1, 2, 3));
    p2.setU(state, Vec3(1, 2, 3));
    system.realize(state, Stage::Velocity);
    system.project(state, 1e-10);

    SimTK_TEST_EQ(c1.getBodyTransform(state), c2.getBodyTransform(state));
    SimTK_TEST_EQ(c1.getBodyVelocity(state), c2.getBodyVelocity(state));

    // Simulate it and see if both pendulums behave identically.
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5);
    system.realize(integ.getState(), Stage::Velocity);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(c1.getBodyTransform(integ.getState()),
                                 c2.getBodyTransform(integ.getState()),
                                 1e-10);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(c1.getBodyVelocity(integ.getState()),
                                 c2.getBodyVelocity(integ.getState()),
                                 1e-10);
}

/*
 * Create two pendulums with gimbal joints, one where the gimbals are faked with a series of pin joints that
 * should provide identical generalized coordinates and speeds (we're expecting the new gimbal definition in
 * which the generalized speeds are the Euler angle derivatives).
 */
TEST(Simbody_MobilizedBody_GimbalVsPinChainEquivalence,
     EulerAngleGimbalMatchesPinChainKinematicsAndDynamics) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.81, 0));
    Body::Rigid lumpy(MassProperties(3.1, Vec3(.1, .2, .3), UnitInertia(1.2, 1.1, 1.3, .01, .02, .03)));
    Body::Rigid massless(MassProperties(0, Vec3(0), UnitInertia(0)));

    Vec3 inboard(0.1, 0.5, -1);
    Vec3 outboard(0.2, -0.2, 0);
    MobilizedBody::Gimbal p1(matter.updGround(), inboard, lumpy, outboard);

    Rotation axisX(Pi / 2, YAxis);  // rotate z into x
    Rotation axisY(-Pi / 2, XAxis); // rotate z into y
    Rotation axisZ;                 // leave z where it is
    MobilizedBody::Pin dummy1(matter.updGround(),
                              Transform(axisX, inboard),
                              massless,
                              Transform(axisX, Vec3(0)));
    MobilizedBody::Pin dummy2(dummy1, Transform(axisY, Vec3(0)), massless, Transform(axisY, Vec3(0)));
    MobilizedBody::Pin p2(dummy2, Transform(axisZ, Vec3(0)), lumpy, Transform(axisZ, outboard));

    MobilizedBody::Weld c1(p1, inboard, lumpy, outboard);
    MobilizedBody::Weld c2(p2, inboard, lumpy, outboard);

    c1.addBodyDecoration(Transform(), DecorativeBrick(.1 * Vec3(1, 2, 3)).setColor(Cyan).setOpacity(0.2));
    c2.addBodyDecoration(Transform(),
                         DecorativeBrick(.1 * Vec3(1, 2, 3))
                             .setColor(Red)
                             .setRepresentation(DecorativeGeometry::DrawWireframe));

    State state = system.realizeTopology();
    p1.setQ(state, Vec3(.1, .2, .3));
    dummy1.setAngle(state, .1);
    dummy2.setAngle(state, .2);
    p2.setAngle(state, .3);

    p1.setU(state, Vec3(1, 2, 3));
    dummy1.setRate(state, 1);
    dummy2.setRate(state, 2);
    p2.setRate(state, 3);

    system.realize(state, Stage::Velocity);
    system.project(state, 1e-10);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(c1.getBodyTransform(state), c2.getBodyTransform(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(c1.getBodyVelocity(state), c2.getBodyVelocity(state));

    // Simulate it and see if both pendulums behave identically.
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper timeStepper(system, integ);

    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(5));

    system.realize(integ.getState(), Stage::Velocity);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(c1.getBodyTransform(integ.getState()),
                                 c2.getBodyTransform(integ.getState()),
                                 1e-10);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(c1.getBodyVelocity(integ.getState()),
                                 c2.getBodyVelocity(integ.getState()),
                                 1e-10);
}

/*
 * Create two pendulums with bushing joints, one where the bushings are faked with a Cartesian (3 sliders) and
 * a series of pin joints that should provide identical generalized coordinates and speeds, except that the
 * translational coordinates are the last three q's for the bushing but are interpreted as translating first,
 * then rotating.
 */
TEST(Simbody_MobilizedBody_BushingVsCompositePinTranslationRotationOrdering,
     EquivalentKinematicsWithReorderedGeneralizedCoordinates) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.81, 0));
    Body::Rigid lumpy(MassProperties(3.1, Vec3(.1, .2, .3), UnitInertia(1.2, 1.1, 1.3, .01, .02, .03)));
    Body::Rigid massless(MassProperties(0, Vec3(0), UnitInertia(0)));

    Vec3 inboard(0.1, 0.5, -1);
    Vec3 outboard(0.2, -0.2, 0);
    MobilizedBody::Bushing p1(matter.updGround(), inboard, lumpy, outboard);

    Rotation axisX(Pi / 2, YAxis);  // rotate z into x
    Rotation axisY(-Pi / 2, XAxis); // rotate z into y
    Rotation axisZ;                 // leave z where it is
    MobilizedBody::Translation dummy0(matter.updGround(), inboard, massless, Transform());
    MobilizedBody::Pin dummy1(dummy0, Transform(axisX, Vec3(0)), massless, Transform(axisX, Vec3(0)));
    MobilizedBody::Pin dummy2(dummy1, Transform(axisY, Vec3(0)), massless, Transform(axisY, Vec3(0)));
    MobilizedBody::Pin p2(dummy2, Transform(axisZ, Vec3(0)), lumpy, Transform(axisZ, outboard));

    MobilizedBody::Weld c1(p1, inboard, lumpy, outboard);
    MobilizedBody::Weld c2(p2, inboard, lumpy, outboard);

    c1.addBodyDecoration(Transform(), DecorativeBrick(.1 * Vec3(1, 2, 3)).setColor(Cyan).setOpacity(0.2));
    c2.addBodyDecoration(Transform(),
                         DecorativeBrick(.1 * Vec3(1, 2, 3))
                             .setColor(Red)
                             .setRepresentation(DecorativeGeometry::DrawWireframe));

    State state = system.realizeTopology();
    p1.setQ(state, Vec6(.1, .2, .3, 1, 2, 3));
    dummy0.setTranslation(state, Vec3(1, 2, 3));
    dummy1.setAngle(state, .1);
    dummy2.setAngle(state, .2);
    p2.setAngle(state, .3);

    p1.setU(state, Vec6(1, 2, 3, -.1, -.2, -.3));
    dummy0.setVelocity(state, Vec3(-.1, -.2, -.3));
    dummy1.setRate(state, 1);
    dummy2.setRate(state, 2);
    p2.setRate(state, 3);

    system.realize(state, Stage::Velocity);
    system.project(state, 1e-10);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(c1.getBodyTransform(state), c2.getBodyTransform(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(c1.getBodyVelocity(state), c2.getBodyVelocity(state));

    // Simulate it and see if both pendulums behave identically.
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper timeStepper(system, integ);

    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(5));

    system.realize(integ.getState(), Stage::Velocity);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(c1.getBodyTransform(integ.getState()),
                                 c2.getBodyTransform(integ.getState()),
                                 1e-10);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(c1.getBodyVelocity(integ.getState()),
                                 c2.getBodyVelocity(integ.getState()),
                                 1e-10);
}
