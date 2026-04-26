#include <algorithm>
#include <gtest/gtest.h>

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

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


TEST(Simbody_HuntCrossleyForce_SphereHalfSpaceContact, NormalAndFrictionForcesMatchAnalyticalModel) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    GeneralForceSubsystem forces(system);
    const Vec3 gravity = Vec3(0, -9.8, 0);
    Force::UniformGravity(forces, matter, gravity, 0);
    const Real radius = 0.8;
    const Real k1 = 1.0;
    const Real k2 = 2.0;
    const Real stiffness1 = std::pow(k1, 2.0 / 3.0);
    const Real stiffness2 = std::pow(k2, 2.0 / 3.0);
    const Real dissipation1 = 0.5;
    const Real dissipation2 = 1.0;
    const Real us1 = 1.0;
    const Real us2 = 0.7;
    const Real ud1 = 0.5;
    const Real ud2 = 0.2;
    const Real uv1 = 0.1;
    const Real uv2 = 0.05;
    Random::Uniform random(0.0, 1.0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Translation sphere(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, sphere, ContactGeometry::Sphere(radius), Transform());
    contacts.addBody(setIndex,
                     matter.updGround(),
                     ContactGeometry::HalfSpace(),
                     Transform(Rotation(-0.5 * Pi, ZAxis), Vec3(0))); // y < 0
    HuntCrossleyForce hc(forces, contacts, setIndex);
    hc.setBodyParameters(ContactSurfaceIndex(0), k1, dissipation1, us1, ud1, uv1);
    hc.setBodyParameters(ContactSurfaceIndex(1), k2, dissipation2, us2, ud2, uv2);
    const Real vt = 0.001;
    hc.setTransitionVelocity(vt);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(vt, hc.getTransitionVelocity());
    State state = system.realizeTopology();

    // Position the sphere at a variety of positions and check the normal force.

    const Real stiffness = stiffness1 * stiffness2 / (stiffness1 + stiffness2);
    for (Real height = radius + 0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        system.realize(state, Stage::Dynamics);
        const Real depth = radius - height;
        Real f = 0;
        if (depth > 0) {
            f = (4.0 / 3.0) * stiffness * depth * std::sqrt(radius * stiffness * depth);
        }
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(
            system.getRigidBodyForces(state, Stage::Dynamics)[sphere.getMobilizedBodyIndex()][1],
            gravity + Vec3(0, f, 0));
    }

    // Now do it with a vertical velocity and see if the dissipation force is correct.
    const Real dissipation =
        (dissipation1 * stiffness2 + dissipation2 * stiffness1) / (stiffness1 + stiffness2);
    for (Real height = radius + 0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        const Real depth = radius - height;
        Real fh = 0;
        if (depth > 0) {
            fh = (4.0 / 3.0) * stiffness * depth * std::sqrt(radius * stiffness * depth);
        }
        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            sphere.setUToFitLinearVelocity(state, Vec3(0, -v, 0));
            system.realize(state, Stage::Dynamics);

            Real f = fh * (1.0 + 1.5 * dissipation * v);
            f = std::max<SimTK::Real>(f, 0);

            EXPECT_NEAR_DEFAULT_TOL_SIMTK(
                system.getRigidBodyForces(state, Stage::Dynamics)[sphere.getMobilizedBodyIndex()][1],
                gravity + Vec3(0, f, 0));
        }
    }

    // Do it with a horizontal velocity and see if the friction force is correct.
    const Real us = 2.0 * us1 * us2 / (us1 + us2);
    const Real ud = 2.0 * ud1 * ud2 / (ud1 + ud2);
    const Real uv = 2.0 * uv1 * uv2 / (uv1 + uv2);
    Vector_<SpatialVec> expectedForce(matter.getNumBodies());
    for (Real height = radius + 0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        const Real depth = radius - height;
        Real fh = 0;
        if (depth > 0) {
            fh = (4.0 / 3.0) * stiffness * depth * std::sqrt(radius * stiffness * depth);
        }
        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            sphere.setUToFitLinearVelocity(state, Vec3(v, 0, 0));
            system.realize(state, Stage::Dynamics);
            const Real vrel = std::abs(v / vt);
            Real ff = (v < 0 ? 1 : -1) * fh
                      * (std::min(vrel, 1.0) * (ud + 2 * (us - ud) / (1 + vrel * vrel)) + uv * std::fabs(v));
            const Vec3 totalForce = gravity + Vec3(ff, fh, 0);
            expectedForce = SpatialVec(Vec3(0), Vec3(0));
            Vec3 contactPointInSphere =
                sphere.findStationAtGroundPoint(state,
                                                Vec3(0, -stiffness1 * depth / (stiffness1 + stiffness2), 0));
            sphere.applyForceToBodyPoint(state, contactPointInSphere, totalForce, expectedForce);
            SpatialVec actualForce =
                system.getRigidBodyForces(state, Stage::Dynamics)[sphere.getMobilizedBodyIndex()];
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(actualForce[0], expectedForce[sphere.getMobilizedBodyIndex()][0]);
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(actualForce[1], expectedForce[sphere.getMobilizedBodyIndex()][1]);
        }
    }
}
