#include <algorithm>
#include <array>
#include <gtest/gtest.h>

#include "SimTKsimbody.h"

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


TEST(Simbody_ElasticFoundationForce, PyramidOnPlane) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    GeneralForceSubsystem forces(system);

    // Create a triangle mesh in the shape of a pyramid, with the square base having area 1 (split into two
    // triangles).
    std::vector<Vec3> vertices;
    vertices.emplace_back(0, 0, 0);
    vertices.emplace_back(1, 0, 0);
    vertices.emplace_back(1, 0, 1);
    vertices.emplace_back(0, 0, 1);
    vertices.emplace_back(0.5, 1, 0.5);

    const std::array<std::array<int, 3>, 6> faces = {
        {{0, 1, 2}, {0, 2, 3}, {1, 0, 4}, {2, 1, 4}, {3, 2, 4}, {0, 3, 4}}};
    std::vector<int> faceIndices;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) {
            faceIndices.push_back(faces[i][j]);
        }
    }

    // Create the mobilized bodies and configure the contact model.
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Translation mesh(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, mesh, ContactGeometry::TriangleMesh(vertices, faceIndices), Transform());
    contacts.addBody(setIndex,
                     matter.updGround(),
                     ContactGeometry::HalfSpace(),
                     Transform(Rotation(-0.5 * Pi, ZAxis), Vec3(0))); // y < 0

    const Real stiffness = 1e9;
    const Real dissipation = 0.01;
    const Real us = 0.1;
    const Real ud = 0.05;
    const Real uv = 0.01;
    const Real vt = 0.01;

    ElasticFoundationForce ef(forces, contacts, setIndex);
    ef.setBodyParameters(ContactSurfaceIndex(0), stiffness, dissipation, us, ud, uv);
    ef.setTransitionVelocity(vt);

    EXPECT_EQ(ef.getTransitionVelocity(), vt);
    State state = system.realizeTopology();

    // Position the pyramid at a variety of positions and check the normal force.
    for (Real depth = -0.1; depth < 0.1; depth += 0.01) {
        mesh.setQToFitTranslation(state, Vec3(0, -depth, 0));
        system.realize(state, Stage::Dynamics);

        Real f = 0;
        if (depth > 0) {
            f = stiffness * depth;
        }

        EXPECT_NEAR_DEFAULT_TOL_SIMTK(
            system.getRigidBodyForces(state, Stage::Dynamics)[mesh.getMobilizedBodyIndex()][1],
            Vec3(0, f, 0));
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(
            system.getRigidBodyForces(state, Stage::Dynamics)[matter.getGround().getMobilizedBodyIndex()][1],
            Vec3(0, -f, 0));
    }

    // Now do it with a vertical velocity and see if the dissipation force is correct.
    for (Real depth = -0.105; depth < 0.1; depth += 0.01) {
        mesh.setQToFitTranslation(state, Vec3(0, -depth, 0));
        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            mesh.setUToFitLinearVelocity(state, Vec3(0, -v, 0));
            system.realize(state, Stage::Dynamics);

            Real f = (depth > 0 ? stiffness * depth * (1 + dissipation * v) : 0);
            f = std::max<SimTK::Real>(f, 0);

            EXPECT_NEAR_DEFAULT_TOL_SIMTK(
                system.getRigidBodyForces(state, Stage::Dynamics)[mesh.getMobilizedBodyIndex()][1],
                Vec3(0, f, 0));
        }
    }

    // Do it with a horizontal velocity and see if the friction force is correct.
    Vector_<SpatialVec> expectedForce(matter.getNumBodies());
    for (Real depth = -0.105; depth < 0.1; depth += 0.01) {
        mesh.setQToFitTranslation(state, Vec3(0, -depth, 0));
        Real fh = 0;
        if (depth > 0) {
            fh = stiffness * depth;
        }

        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            mesh.setUToFitLinearVelocity(state, Vec3(v, 0, 0));
            system.realize(state, Stage::Dynamics);
            const Real vrel = std::abs(v / vt);
            Real ff = (v < 0 ? 1 : -1) * fh
                      * (std::min(vrel, 1.0) * (ud + 2 * (us - ud) / (1 + vrel * vrel)) + uv * std::fabs(v));

            const Vec3 totalForce = Vec3(ff, fh, 0);
            expectedForce = SpatialVec(Vec3(0), Vec3(0));

            const Vec3 contactPoint1 = mesh.findStationAtGroundPoint(state, Vec3(2.0 / 3.0, 0, 1.0 / 3.0));
            mesh.applyForceToBodyPoint(state, contactPoint1, 0.5 * totalForce, expectedForce);

            const Vec3 contactPoint2 = mesh.findStationAtGroundPoint(state, Vec3(1.0 / 3.0, 0, 2.0 / 3.0));
            mesh.applyForceToBodyPoint(state, contactPoint2, 0.5 * totalForce, expectedForce);

            const SpatialVec actualForce =
                system.getRigidBodyForces(state, Stage::Dynamics)[mesh.getMobilizedBodyIndex()];

            EXPECT_NEAR_DEFAULT_TOL_SIMTK(actualForce[0], expectedForce[mesh.getMobilizedBodyIndex()][0]);
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(actualForce[1], expectedForce[mesh.getMobilizedBodyIndex()][1]);
        }
    }
}

/**
 * @brief This test compares the numerical result of a sphere
 *        in contact with a plane, using the elastic foundation
 *        model.
 *        The analytical solution of this problem is given by
 *        the product of the stiffness with the volume of
 *        the sphere in the plane, i.e the volume of a spherical
 *        cap.
 *        The volume of a spherical cap is:
 *          Vcap = Pi*h*h/3.0*(3.0*r-h)
 *        where
 *          r is the radis of the sphere
 *          h the height of the cap. In our case, the penetration
 *          depth
 * @note If we want to go further, we can observe that doubling
 *       the penetration depth results in multiplying the normal
 *       effort by 4.
 *       This is different from Hertz theory, where doubling the
 *       penetration depth results in multiplying
 *       the normal effort by 2^(3/2)~2.68
 *
 */
TEST(Simbody_ElasticFoundationForce, SphereOnPlane) {
    // Material properties for sphere
    const Real stiffness = 1e9;
    const Real dissipation = 0.0;
    const Real us = 0.0;
    const Real ud = 0.0;
    const Real uv = 0.0;
    const Real vt = 0.0;

    // Sphere radius
    const Real radius = 1.0;

    // Define initial penetration
    const Real initialPenetration = 0.002;

    // Define the number of tests to perform
    const int maxLevel = 6;

    // Define some tolerances for each level in %
    static constexpr std::array<Real, 6> tolerances = {0.15, 0.07, 0.03, 0.02, 0.01, 0.02};

    for (int i = 0; i < maxLevel; ++i) {
        // For each level, penetration is double
        const Real penetration = initialPenetration * pow(2.0, (Real)i);

        // Creation of the classical problem
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        GeneralContactSubsystem contacts(system);
        GeneralForceSubsystem forces(system);
        const ContactSetIndex setIndex = contacts.createContactSet();

        // Creation a sphere with 6 levels of refinement
        const PolygonalMesh sphereMesh(PolygonalMesh::createSphereMesh(radius, 6));

        // Create the mobilized bodies and configure the contact model.
        const Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
        const MobilizedBody::Translation mesh(matter.updGround(), Transform(), body, Transform());
        contacts.addBody(setIndex, mesh, ContactGeometry::TriangleMesh(sphereMesh), Transform());
        contacts.addBody(setIndex,
                         matter.updGround(),
                         ContactGeometry::HalfSpace(),
                         Transform(Rotation(-0.5 * Pi, ZAxis),
                                   Vec3(0.0, penetration - radius, 0.0))); // y < penetration-radius

        ElasticFoundationForce ef(forces, contacts, setIndex);
        ef.setBodyParameters(ContactSurfaceIndex(0), stiffness, dissipation, us, ud, uv);
        ef.setTransitionVelocity(vt);

        const State state = system.realizeTopology();
        system.realize(state, Stage::Dynamics);

        const SpatialVec r = system.getRigidBodyForces(state, Stage::Dynamics)[mesh.getMobilizedBodyIndex()];
        const Real volumeSphericalCap = Pi * penetration * penetration / 3.0 * (3.0 * radius - penetration);
        const Real theoreticalResult = stiffness * volumeSphericalCap;
        const Real numericalResult = r[1][1];

        EXPECT_LT(std::abs(r[1][0]), 1e-6);
        EXPECT_LT(std::abs(r[1][2]), 1e-6);

        const Real relativeDifference = std::abs((numericalResult / theoreticalResult) - 1.0);
        EXPECT_LT(relativeDifference, tolerances[i]);
    }
}

TEST(Simbody_ElasticFoundationForce, SphereOnPlaneOldFormulation) {
    // Global stiffness of the contact: each material will have
    // twice this stiffness to obtain this global stiffness in the contact
    // 1/kG = 1/k1 + 1/k2
    const Real stiffness = 1e9;
    const Real dissipation = 0.0;
    const Real us = 0.0;
    const Real ud = 0.0;
    const Real uv = 0.0;
    const Real vt = 1.0e-2;

    // Sphere radius
    const Real radius = 1.0;

    // Define initial penetration
    const Real initialPenetration = 0.002;
    const int maxLevel = 6;

    // Define some tolerances for each level in %
    static constexpr std::array<Real, 6> tolerances = {0.15, 0.07, 0.03, 0.02, 0.01, 0.02};

    for (int i = 0; i < maxLevel; ++i) {
        // For each level, penetration is double
        const Real penetration = initialPenetration * pow(2.0, (Real)i);

        // Creation of the classical problem
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        ContactTrackerSubsystem tracker(system);
        CompliantContactSubsystem contactForces(system, tracker);
        contactForces.setTransitionVelocity(vt);
        matter.Ground().updBody().addContactSurface(
            Transform(Rotation(-0.5 * Pi, ZAxis),
                      Vec3(0.0, penetration - radius, 0.0)), // y < penetration-radius
            ContactSurface(ContactGeometry::HalfSpace(),
                           ContactMaterial(2.0 * stiffness, dissipation, us, ud, uv),
                           1.0));
        Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));

        body.addContactSurface(
            Transform(),
            ContactSurface(ContactGeometry::TriangleMesh(PolygonalMesh::createSphereMesh(radius, 6)),
                           ContactMaterial(2.0 * stiffness, dissipation, us, ud, uv),
                           1.0));

        const MobilizedBody::Translation mesh(matter.updGround(), Transform(), body, Transform());

        const State state = system.realizeTopology();
        system.realize(state, Stage::Dynamics);
        EXPECT_EQ(contactForces.getNumContactForces(state), 1);

        const ContactForce& force = contactForces.getContactForce(state, 0);
        const Vec3& frc = force.getForceOnSurface2()[1];
        EXPECT_LT(std::abs(frc[0]), 1e-6);
        EXPECT_LT(std::abs(frc[2]), 1e-6);

        const Real numericalResult = frc[1];
        const Real volumeSphericalCap = Pi * penetration * penetration / 3.0 * (3.0 * radius - penetration);
        const Real theoreticalResult = stiffness * volumeSphericalCap;
        const Real relativeDifference = std::abs((numericalResult / theoreticalResult) - 1.0);
        EXPECT_LT(relativeDifference, tolerances[i]);
    }
}
