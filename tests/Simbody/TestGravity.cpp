#include <gtest/gtest.h>
#include <iostream>

#include "Simbody.h"
using std::cout;
using std::endl;

using namespace SimTK;


#define EXPECT_SIMTK_SIZE(expected, actual, n)                            \
    EXPECT_TRUE(SimTK::Test::numericallyEqual((expected), (actual), (n))) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\nSize:      " << (n) << "\n"

#define EXPECT_NEAR_CUSTOM_TOL_SIMTK(expected, actual, tol)                              \
    EXPECT_TRUE(SimTK::Test::numericallyEqual((expected), (actual), /*scale*/ 1, (tol))) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\nTolerance: " << (tol) << "\n"

#define EXPECT_NEAR_DEFAULT_TOL_SIMTK(expected, actual) \
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((expected), (actual), 1e-6)

#define EXPECT_NOT_NEAR_DEFAULT_TOL_SIMTK(expected, actual)                        \
    EXPECT_FALSE(SimTK::Test::numericallyEqual((expected), (actual), /*scale*/ 1)) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\n"


const Real Mass1 = 100, Mass2 = 5, Mass3 = 1;
const Vec3 Centroid(.5, 0, .5);
const Vec3 COM1 = Centroid, COM2 = Centroid, COM3 = Centroid + Vec3(0, .5, 0);
const UnitInertia Central(Vec3(.1), Vec3(.05));

// Simbody requires inertias to be expressed about body origin rather than COM.
const Inertia Inertia1 = Mass1 * Central.shiftFromCentroid(-COM1);
const Inertia Inertia2 = Mass2 * Central.shiftFromCentroid(-COM2);
const Inertia Inertia3 = Mass3 * Central.shiftFromCentroid(-COM3); // weird

Body::Rigid body1Info(MassProperties(Mass1, COM1, Inertia1));
Body::Rigid body2Info(MassProperties(Mass2, COM2, Inertia2));
Body::Rigid body3Info(MassProperties(Mass3, COM3, Inertia3));

// Make sure the constructors and default setters and getters work.
TEST(Simbody_Gravity_Construction, DefaultParametersAndGettersAreConsistent) {
    MultibodySystem mbs;
    SimbodyMatterSubsystem matter(mbs);
    GeneralForceSubsystem forces(mbs);

    Force::Gravity gravity1(forces, matter, -ZAxis, 50);
    EXPECT_EQ(gravity1.getDefaultMagnitude(), 50);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity1.getDefaultDownDirection(), UnitVec3(0, 0, -1));
    EXPECT_EQ(gravity1.getDefaultZeroHeight(), 0);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity1.getDefaultGravityVector(), Vec3(0, 0, -50));
    EXPECT_TRUE(gravity1.getDefaultBodyIsExcluded(MobodIndex(0)));
    EXPECT_FALSE(gravity1.getDefaultBodyIsExcluded(MobodIndex(1)));

    const Vec3 grav2(1, 2, 3);
    Force::Gravity gravity2(forces, matter, grav2);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity2.getDefaultMagnitude(), grav2.norm());
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity2.getDefaultDownDirection(), UnitVec3(grav2));
    EXPECT_EQ(gravity2.getDefaultZeroHeight(), 0);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity2.getDefaultGravityVector(), grav2);

    mbs.setUpDirection(XAxis);
    const Real mag = 16.75;
    Force::Gravity gravity3(forces, matter, mag);

    EXPECT_EQ(gravity3.getDefaultMagnitude(), mag);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity3.getDefaultDownDirection(), UnitVec3(-XAxis));
    EXPECT_EQ(gravity3.getDefaultZeroHeight(), 0);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity3.getDefaultGravityVector(), Vec3(-mag, 0, 0));

    // Using the vector constructor with a zero vector should pluck the
    // direction out of the System.
    Force::Gravity gravity4(forces, matter, Vec3(0));
    EXPECT_EQ(gravity4.getDefaultMagnitude(), 0);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity4.getDefaultDownDirection(), UnitVec3(-XAxis));
    EXPECT_EQ(gravity4.getDefaultZeroHeight(), 0);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity4.getDefaultGravityVector(), Vec3(0));

    // Make sure Ground can't be included.
    gravity1.setDefaultBodyIsExcluded(MobodIndex(0), false);
    EXPECT_TRUE(gravity1.getDefaultBodyIsExcluded(MobodIndex(0)));

    // Make sure we can exclude a body and that it doesn't leak over into
    // another body. Note that for defaults we don't yet know how may bodies
    // there might be so we can use arbitrary body numbers.
    EXPECT_FALSE(gravity1.getDefaultBodyIsExcluded(MobodIndex(13)));

    gravity1.setDefaultBodyIsExcluded(MobodIndex(13), true);
    EXPECT_TRUE(gravity1.getDefaultBodyIsExcluded(MobodIndex(13)));
    EXPECT_FALSE(gravity1.getDefaultBodyIsExcluded(MobodIndex(5)));

    gravity1.setDefaultGravityVector(Vec3(19, 20, 21));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity1.getDefaultGravityVector(), Vec3(19, 20, 21));

    gravity1.setDefaultMagnitude(3).setDefaultDownDirection(grav2);
    EXPECT_EQ(gravity1.getDefaultMagnitude(), 3);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity1.getDefaultDownDirection(), UnitVec3(grav2));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity1.getDefaultGravityVector(), Real(3) * UnitVec3(grav2));

    gravity1.setDefaultZeroHeight(5.25);
    EXPECT_EQ(gravity1.getDefaultZeroHeight(), 5.25);
}

// Make sure we can override default values in the State.
TEST(Simbody_Gravity_Parameters_StateOverridesAndDefaultPropagation,
     OverridesDefaultsAndPreservesIndependence) {
    MultibodySystem mbs;
    SimbodyMatterSubsystem matter(mbs);
    GeneralForceSubsystem forces(mbs);

    Force::Gravity gravity(forces, matter, -ZAxis, 49, -5);
    EXPECT_EQ(gravity.getDefaultZeroHeight(), -5);

    MobilizedBody::Weld mobod1(matter.Ground(), Vec3(0), body1Info, Vec3(0));
    MobilizedBody::Pin mobod2(mobod1, Vec3(0), body2Info, Vec3(0));
    MobilizedBody::Pin mobod3(mobod2, Vec3(0), body3Info, Vec3(0));
    MobilizedBody::Pin mobod4(mobod3, Vec3(0), body3Info, Vec3(0));
    MobilizedBody::Free mobod5(matter.Ground(), Vec3(0), body3Info, Vec3(0));

    State state = mbs.realizeTopology();

    // Make sure defaults made it to the state.
    EXPECT_TRUE(gravity.getBodyIsExcluded(state, MobodIndex(0)));
    for (MobodIndex i(1); i < matter.getNumBodies(); ++i) {
        EXPECT_FALSE(gravity.getBodyIsExcluded(state, i));
    }
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getGravityVector(state), Vec3(0, 0, -49));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getDownDirection(state), UnitVec3(0, 0, -1));
    EXPECT_EQ(gravity.getMagnitude(state), 49);
    EXPECT_EQ(gravity.getZeroHeight(state), -5);

    // Now make some changes in the state.

    // Shouldn't be able to include ground.
    gravity.setBodyIsExcluded(state, MobodIndex(0), false);
    EXPECT_TRUE(gravity.getBodyIsExcluded(state, MobodIndex(0)));

    gravity.setBodyIsExcluded(state, mobod3, true);
    EXPECT_TRUE(gravity.getBodyIsExcluded(state, mobod3));
    EXPECT_FALSE(gravity.getBodyIsExcluded(state, mobod2));
    EXPECT_FALSE(gravity.getBodyIsExcluded(state, mobod4));

    // That shouldn't have changed the default.
    EXPECT_FALSE(gravity.getDefaultBodyIsExcluded(mobod3));

    // When making changes in the state we are restricted to bodies that actually exist.
    EXPECT_THROW(gravity.setBodyIsExcluded(state, MobodIndex(13), true), std::exception);

    gravity.setGravityVector(state, Vec3(1, 2, 3));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getGravityVector(state), Vec3(1, 2, 3));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getDownDirection(state), UnitVec3(Vec3(1, 2, 3)));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getMagnitude(state), Vec3(1, 2, 3).norm());

    gravity.setDownDirection(state, UnitVec3(9, 10, 11));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getDownDirection(state), UnitVec3(9, 10, 11));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getMagnitude(state), Vec3(1, 2, 3).norm()); // no change
    EXPECT_THROW(gravity.setMagnitude(state, -5), std::exception);                    // must be >= 0

    gravity.setMagnitude(state, 5);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getDownDirection(state), UnitVec3(9, 10, 11)); // no change
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getMagnitude(state), 5);

    // Changing gravity vector to zero should leave direction unchanged.
    gravity.setGravityVector(state, Vec3(0));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getGravityVector(state), Vec3(0));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getDownDirection(state), UnitVec3(9, 10, 11));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getMagnitude(state), 0);

    gravity.setZeroHeight(state, 1.25);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getZeroHeight(state), 1.25);
}

TEST(Simbody_Gravity_ForcesDynamicsCachingAndStateDependence,
     ComputesCorrectForcesAndRespectsCachingAndExclusions) {
    MultibodySystem mbs;
    SimbodyMatterSubsystem matter(mbs);
    GeneralForceSubsystem forces(mbs);
    Force::DiscreteForces discrete(forces, matter);
    Force::Gravity gravity(forces, matter, -ZAxis, 50, 17);
    MobilizedBody::Weld mobod1(matter.Ground(), Vec3(0), body1Info, Vec3(0));

    const Rotation ZtoY(-Pi / 2, XAxis);
    MobilizedBody::Pin mobod2(mobod1, Transform(ZtoY, 2 * Centroid), body2Info, Transform(ZtoY, Vec3(0)));
    MobilizedBody::Pin mobod3(mobod2, Transform(ZtoY, 2 * Centroid), body3Info, Transform(ZtoY, Vec3(0)));

    State state = mbs.realizeTopology();
    state.updQ() = SimTK::Test::randVector(state.getNQ());
    state.updU() = SimTK::Test::randVector(state.getNU());

    mbs.realize(state, Stage::Dynamics);

    const Vector_<SpatialVec>& bodyForces = mbs.getRigidBodyForces(state, Stage::Dynamics);

    // Mobility forces should all be zero.
    const Vector& mobilityForces = mbs.getMobilityForces(state, Stage::Dynamics);
    for (int i = 0; i < state.getNU(); ++i) {
        EXPECT_EQ(mobilityForces[i], 0);
    }

    // Calculate body forces and torques and verify that they are correct.
    // (These are about the body origin, not the COM.)
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(bodyForces[0], SpatialVec(Vec3(0))); // Ground

    const Real g = gravity.getMagnitude(state);
    const UnitVec3 d = gravity.getDownDirection(state);
    Real h0 = gravity.getZeroHeight(state);

    Real pe = 0;
    for (MobodIndex i(1); i < matter.getNumBodies(); ++i) {
        const Mobod& mobod = matter.getMobilizedBody(i);
        const Real m = mobod.getBodyMass(state);
        const Vec3 p_BC = mobod.getBodyMassCenterStation(state);
        const Vec3 p_BC_G = mobod.expressVectorInGroundFrame(state, p_BC);
        const Vec3 F = m * g * d;
        const Vec3 M = p_BC_G % F;
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(bodyForces[i], SpatialVec(M, F));

        EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForce(state, i), SpatialVec(M, F));

        const Real h = -dot(mobod.findStationLocationInGround(state, p_BC), d) - h0;
        pe += m * g * h;
    }

    EXPECT_EQ(gravity.getPotentialEnergy(state), pe);
    EXPECT_EQ(mbs.calcPotentialEnergy(state), pe);

    // Change zero height and verify that the potential energy changes.
    h0 = 10.3;
    gravity.setZeroHeight(state, h0);
    mbs.realize(state, Stage::Dynamics);
    pe = 0;
    for (MobodIndex i(1); i < matter.getNumBodies(); ++i) {
        const Mobod& mobod = matter.getMobilizedBody(i);
        const Real m = mobod.getBodyMass(state);
        const Vec3 p_BC = mobod.getBodyMassCenterStation(state);
        const Real h = -dot(mobod.findStationLocationInGround(state, p_BC), d) - h0;
        pe += m * g * h;
    }
    EXPECT_EQ(gravity.getPotentialEnergy(state), pe);
    EXPECT_EQ(mbs.calcPotentialEnergy(state), pe);

    // Turn off a body and make sure it doesn't see gravity after that, and
    // that the other bodies are unchanged.
    gravity.setBodyIsExcluded(state, mobod2, true);
    mbs.realize(state, Stage::Dynamics);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(bodyForces[0], SpatialVec(Vec3(0))); // Ground
    pe = 0;
    for (MobodIndex i(1); i < matter.getNumBodies(); ++i) {
        const Mobod& mobod = matter.getMobilizedBody(i);
        const Real m = mobod.getBodyMass(state);
        const Vec3 p_BC = mobod.getBodyMassCenterStation(state);
        const Vec3 p_BC_G = mobod.expressVectorInGroundFrame(state, p_BC);
        const Vec3 F = m * g * d;
        const Vec3 M = p_BC_G % F;

        if (i != mobod2.getMobilizedBodyIndex()) {
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(bodyForces[i], SpatialVec(M, F));
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForce(state, i), SpatialVec(M, F));
            const Real h = -dot(mobod.findStationLocationInGround(state, p_BC), d) - h0;
            pe += m * g * h;
        } else {
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(bodyForces[i], SpatialVec(Vec3(0)));
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForce(state, i), SpatialVec(Vec3(0)));
        }
    }
    EXPECT_EQ(gravity.getPotentialEnergy(state), pe);
    EXPECT_EQ(mbs.calcPotentialEnergy(state), pe);

    // Test caching.
    const long long nevals1 = gravity.getNumEvaluations();
    EXPECT_TRUE(gravity.isForceCacheValid(state));

    // This should not require re-evaluation.
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForces(state)[0], SpatialVec(Vec3(0)));
    EXPECT_EQ(gravity.getNumEvaluations(), nevals1);

    // Force re-evaluation.
    gravity.invalidateForceCache(state);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForces(state)[0], SpatialVec(Vec3(0)));

    const long long nevals2 = gravity.getNumEvaluations();
    EXPECT_EQ(nevals2, nevals1 + 1);

    state.invalidateAllCacheAtOrAbove(Stage::Velocity); // shouldn't invalidate
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForces(state)[0], SpatialVec(Vec3(0)));
    EXPECT_EQ(gravity.getNumEvaluations(), nevals2);

    mbs.realize(state, Stage::Dynamics); // shouldn't reevaluate
    EXPECT_EQ(gravity.getNumEvaluations(), nevals2);

    // Bring mobod2 back in. This should only invalidate Dynamics stage, but
    // should nevertheless force recomputation of gravity.
    gravity.setBodyIsExcluded(state, mobod2, false);
    EXPECT_EQ(state.getSystemStage(), Stage(Stage::Dynamics - 1));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForces(state)[0], SpatialVec(Vec3(0)));

    const long long nevals3 = gravity.getNumEvaluations();
    EXPECT_EQ(nevals3, nevals2 + 1);

    // Make sure that setting gravity to zero works properly -- it is a
    // special case since the zeroes are precalculated. This should not
    // require a gravity evaluation.
    gravity.setMagnitude(state, 0);
    EXPECT_FALSE(gravity.isForceCacheValid(state));
    for (int i = 0; i < matter.getNumBodies(); ++i) {
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForces(state)[i], SpatialVec(Vec3(0)));
    }
    EXPECT_EQ(gravity.getNumEvaluations(), nevals3);

    // Setting to non-zero should invalidate, and then require just a single
    // evaluation to respond to multiple calls.
    gravity.setMagnitude(state, 9.8);
    EXPECT_FALSE(gravity.isForceCacheValid(state));
    for (int i = 1; i < matter.getNumBodies(); ++i) {
        EXPECT_NOT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForces(state)[i], SpatialVec(Vec3(0)));
    }
    const long long nevals4 = gravity.getNumEvaluations();
    EXPECT_EQ(nevals4, nevals3 + 1);

    // Turn off velocities for gravity compensation test to eliminate Coriolis
    // foces (shouldn't invalidate gravity forces).
    state.updU() = 0;
    mbs.realize(state, Stage::Acceleration);
    const Vector_<SpatialVec>& gfrc = gravity.getBodyForces(state);
    Vector f;
    matter.multiplyBySystemJacobianTranspose(state, gfrc, f);
    discrete.setAllMobilityForces(state, -f);
    mbs.realize(state, Stage::Acceleration);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(state.getUDot(), Vector(state.getNU(), Real(0)));
    EXPECT_EQ(gravity.getNumEvaluations(), nevals4); // all for free?

    // Sneaking a zero in by vector should behave just like setting the
    // magnitude to zero.
    gravity.setGravityVector(state, Vec3(0));
    EXPECT_FALSE(gravity.isForceCacheValid(state));
    for (int i = 0; i < matter.getNumBodies(); ++i) {
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(gravity.getBodyForces(state)[i], SpatialVec(Vec3(0)));
    }

    EXPECT_EQ(gravity.getNumEvaluations(), nevals4); // no eval needed
}
