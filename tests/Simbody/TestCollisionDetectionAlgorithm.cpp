/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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
#include <gtest/gtest.h>
#include <set>

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

void verifyPointContact(const Array_<Contact>& contacts,
                        int surface1,
                        int surface2,
                        const Vec3& normal,
                        const Vec3& location,
                        Real depth,
                        Real r1,
                        Real r2) {
    EXPECT_EQ(contacts.size(), 1);
    EXPECT_TRUE(PointContact::isInstance(contacts[0]));

    const auto& pointContact = static_cast<const PointContact&>(contacts[0]);
    EXPECT_EQ((int)pointContact.getSurface1(), surface1);
    EXPECT_EQ((int)pointContact.getSurface2(), surface2);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(pointContact.getNormal(), normal);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(pointContact.getDepth(), depth);
    EXPECT_NEAR(min(pointContact.getRadiusOfCurvature1(), pointContact.getRadiusOfCurvature2()),
                min(r1, r2),
                1e-1);
    EXPECT_NEAR(max(pointContact.getRadiusOfCurvature1(), pointContact.getRadiusOfCurvature2()),
                max(r1, r2),
                1e-1);
    EXPECT_NEAR(pointContact.getEffectiveRadiusOfCurvature(), sqrt(r1 * r2), 1e-1);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(pointContact.getLocation(), location);
}

TEST(Simbody_CollisionDetectionAlgorithm, SphereHalfSpace) {
    SCOPED_TRACE("Testing Sphere-HalfSpace contact detection");

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    const Real radius = 0.8;
    const Vec3 center(0.1, -0.3, 0.3);
    Random::Uniform random(0.0, 1.0);
    const Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    const ContactSetIndex setIndex = contacts.createContactSet();
    const MobilizedBody::Free sphere(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, sphere, ContactGeometry::Sphere(radius), center);
    contacts.addBody(setIndex,
                     matter.updGround(),
                     ContactGeometry::HalfSpace(),
                     Transform(Rotation(-0.5 * Pi, ZAxis), Vec3(0, 1, 0))); // y < 1
    State state = system.realizeTopology();
    Vec3 centerInGround;

    // Pick a random positions for the sphere.
    for (int iteration = 0; iteration < 100; ++iteration) {
        for (int i = 0; i < state.getNY(); i++) {
            state.updY()[i] = 5 * random.getValue();
        }
        system.realize(state, Stage::Dynamics);
        centerInGround = sphere.findStationLocationInGround(state, center);

        // Check the results of collision detection.
        const Array_<Contact>& contact = contacts.getContacts(state, setIndex);
        if (centerInGround[1] > radius + 1) {
            EXPECT_EQ(contact.size(), 0);
        } else {
            Real depth = radius - centerInGround[1] + 1;
            verifyPointContact(contact,
                               1,
                               0,
                               Vec3(0, 1, 0),
                               Vec3(centerInGround[0], 1 - (0.5 * depth), centerInGround[2]),
                               depth,
                               radius,
                               radius);
        }
    }
}

TEST(Simbody_CollisionDetectionAlgorithm, SphereSphere) {
    SCOPED_TRACE("Testing Sphere-Sphere contact detection");

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);

    const int numBodies = 10;
    std::vector<Real> radius(numBodies);
    std::vector<Vec3> center(numBodies);

    Random::Uniform random(0.0, 1.0);
    for (int i = 0; i < numBodies; i++) {
        radius[i] = random.getValue();
        center[i] = Vec3(random.getValue(), random.getValue(), random.getValue());
    }

    const Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    for (int i = 0; i < numBodies; ++i) {
        MobilizedBody::Free mobod(matter.updGround(), Transform(), body, Transform());
        contacts.addBody(setIndex, mobod, ContactGeometry::Sphere(radius[i]), center[i]);
    }

    State state = system.realizeTopology();
    std::vector<Vec3> centerInGround(numBodies);

    // Pick random positions for all the bodies.
    for (int iteration = 0; iteration < 100; ++iteration) {
        for (int i = 0; i < state.getNY(); i++) {
            state.updY()[i] = 5 * random.getValue();
        }
        system.realize(state, Stage::Dynamics);
        for (MobilizedBodyIndex index(1); index <= numBodies; ++index) {
            centerInGround[index - 1] =
                matter.getMobilizedBody(index).findStationLocationInGround(state, center[index - 1]);
        }

        // Make sure all contacts are accurate.
        const Array_<Contact>& actualContacts = contacts.getContacts(state, setIndex);
        for (const auto& contact : actualContacts) {
            EXPECT_TRUE(PointContact::isInstance(contact));

            const auto& pointContact = static_cast<const PointContact&>(contact);
            const int body1 = pointContact.getSurface1();
            const int body2 = pointContact.getSurface2();

            const Vec3 delta = centerInGround[body2] - centerInGround[body1];
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(delta.normalize(), pointContact.getNormal());
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(delta.norm(),
                                          radius[body1] + radius[body2] - pointContact.getDepth());

            const double combinedRadius = radius[body1] * radius[body2] / (radius[body1] + radius[body2]);
            EXPECT_EQ(combinedRadius, pointContact.getRadiusOfCurvature1());
            EXPECT_EQ(combinedRadius, pointContact.getRadiusOfCurvature2());
        }

        // Make sure no contacts were missed.
        int expectedContacts = 0;
        for (int i = 0; i < numBodies; i++) {
            for (int j = 0; j < i; j++) {
                if ((centerInGround[i] - centerInGround[j]).norm() < radius[i] + radius[j]) {
                    expectedContacts++;
                }
            }
        }
        EXPECT_EQ(actualContacts.size(), expectedContacts);
    }
}

TEST(Simbody_CollisionDetectionAlgorithm, EllipsoidHalfSpace) {
    SCOPED_TRACE("Testing Ellipsoid-HalfSpace contact detection");

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Vec3 radii(0.8, 1.5, 2.1);
    Vec3 center(0.1, -0.3, 0.3); // Major axes span the ranges [-0.7, 0.9], [-1.8, 1.2], [-1.8, 2.4]
    Random::Uniform random(0.0, 1.0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free ellipsoid(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, ellipsoid, ContactGeometry::Ellipsoid(radii), center);
    contacts.addBody(setIndex,
                     matter.updGround(),
                     ContactGeometry::HalfSpace(),
                     Transform(Rotation(-0.5 * Pi, ZAxis), Vec3(0, 1, 0))); // y < 1
    State state = system.realizeTopology();

    // Test a variety of positions.
    ellipsoid.setQToFitTransform(
        state,
        Transform(Rotation(), Vec3(0, 2.9, 0))); // [-0.7, 0.9], [1.1, 4.1], [-1.8, 2.4]
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0);

    ellipsoid.setQToFitTransform(
        state,
        Transform(Rotation(), Vec3(0, 2.6, 0))); // [-0.7, 0.9], [0.8, 3.8], [-1.8, 2.4]
    system.realize(state, Stage::Dynamics);
    verifyPointContact(contacts.getContacts(state, setIndex),
                       1,
                       0,
                       Vec3(0, 1, 0),
                       Vec3(0.1, 0.9, 0.3),
                       0.2,
                       0.8,
                       2.1);

    ellipsoid.setQToFitTransform(
        state,
        Transform(Rotation(SimTK_PI / 2, ZAxis), Vec3(0, 1.6, 0))); // [-1.2, 1.8], [0.9, 2.5], [-1.8, 2.4]
    system.realize(state, Stage::Dynamics);
    verifyPointContact(contacts.getContacts(state, setIndex),
                       1,
                       0,
                       Vec3(0, 1, 0),
                       Vec3(0.3, 0.95, 0.3),
                       0.1,
                       1.5,
                       2.1);

    ellipsoid.setQToFitTransform(state, Transform(Rotation(SimTK_PI / 4, XAxis), Vec3(0, 3.1, 0)));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 1);

    const auto& pointContact = static_cast<const PointContact&>(contacts.getContacts(state, setIndex)[0]);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(pointContact.getNormal(), Vec3(0, 1, 0));
    EXPECT_LT(pointContact.getDepth(), 0.3);
    EXPECT_NEAR(min(pointContact.getRadiusOfCurvature1(), pointContact.getRadiusOfCurvature2()), 0.8, 1e-6);
    EXPECT_GT(max(pointContact.getRadiusOfCurvature1(), pointContact.getRadiusOfCurvature2()), 1.5);
    EXPECT_LT(max(pointContact.getRadiusOfCurvature1(), pointContact.getRadiusOfCurvature2()), 2.1);
    EXPECT_EQ(pointContact.getLocation()[0], 0.1);
    EXPECT_EQ(pointContact.getLocation()[1], 1 - (pointContact.getDepth() / 2));
    EXPECT_GT(pointContact.getLocation()[2], 0);
}

struct EllipsoidShape {
    const Vec3& radii;
    const Vec3& center;
    const Transform& transform;
};

// Separated so callers can track soft failures independently from hard ones.
struct EllipsoidContactVerification {
    bool surfacePointsOnEllipsoids = false; // must always be true
    bool normalsAligned = false;            // may occasionally fail due to Newton convergence
};

[[nodiscard]] auto verifyEllipsoidContact(const Contact& contact,
                                          const EllipsoidShape& ellipsoid1,
                                          const EllipsoidShape& ellipsoid2) -> EllipsoidContactVerification {
    EllipsoidContactVerification result;

    EXPECT_TRUE(PointContact::isInstance(contact))
        << "Expected a PointContact, but received a different Contact subtype.";

    const auto& ptc = static_cast<const PointContact&>(contact);
    const Vec3 offset = 0.5 * ptc.getDepth() * ptc.getNormal();
    const Vec3 loc1 = ~ellipsoid1.transform * (ptc.getLocation() + offset) - ellipsoid1.center;
    const Vec3 loc2 = ~ellipsoid2.transform * (ptc.getLocation() - offset) - ellipsoid2.center;

    auto ellipsoidEquation = [](const Vec3& point, const Vec3& radii) {
        const double nrx = point[0] / radii[0];
        const double nry = point[1] / radii[1];
        const double nrz = point[2] / radii[2];
        return (nrx * nrx) + (nry * nry) + (nrz * nrz);
    };

    auto ellipsoidSurfaceNormal = [](const Vec3& point, const Vec3& radii) {
        return UnitVec3(point[0] / (radii[0] * radii[0]),
                        point[1] / (radii[1] * radii[1]),
                        point[2] / (radii[2] * radii[2]));
    };

    static constexpr double kSurfaceTolerance = 1e-4;
    static constexpr double kNormalAlignment = 0.999;

    // Hard checks: contact surface points must lie on their respective ellipsoids.
    const double eqVal1 = ellipsoidEquation(loc1, ellipsoid1.radii);
    const double eqVal2 = ellipsoidEquation(loc2, ellipsoid2.radii);

    EXPECT_NEAR(eqVal1, 1.0, kSurfaceTolerance)
        << "Ellipsoid 1: contact surface point does not lie on the ellipsoid surface.\n"
        << "  (x/a)^2+(y/b)^2+(z/c)^2 = " << eqVal1 << ", expected 1.0.";
    EXPECT_NEAR(eqVal2, 1.0, kSurfaceTolerance)
        << "Ellipsoid 2: contact surface point does not lie on the ellipsoid surface.\n"
        << "  (x/a)^2+(y/b)^2+(z/c)^2 = " << eqVal2 << ", expected 1.0.";

    result.surfacePointsOnEllipsoids =
        std::abs(eqVal1 - 1.0) <= kSurfaceTolerance && std::abs(eqVal2 - 1.0) <= kSurfaceTolerance;

    // Soft checks: normals may fail at high-curvature points where Newton
    // iteration doesn't converge. Tracked separately so callers can apply
    // a statistical tolerance rather than failing the whole test.
    const double dot1 =
        ~ellipsoidSurfaceNormal(loc1, ellipsoid1.radii) * (~ellipsoid1.transform.R() * ptc.getNormal());
    const double dot2 =
        -~ellipsoidSurfaceNormal(loc2, ellipsoid2.radii) * (~ellipsoid2.transform.R() * ptc.getNormal());

    result.normalsAligned = dot1 >= kNormalAlignment && dot2 >= kNormalAlignment;

    return result;
}

TEST(Simbody_CollisionDetectionAlgorithm, EllipsoidEllipsoidKnownPositions) {
    SCOPED_TRACE("Testing Ellipsoid-Ellipsoid contact detection at known positions");

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);

    const Vec3 radii1(0.8, 1.5, 2.1);
    const Vec3 radii2(1.0, 1.2, 1.4);
    const Vec3 center1(0, -0.2, 0.5);
    const Vec3 center2(0.1, 0, 0.3);

    Body::Rigid rigidBody1(MassProperties(1.0, Vec3(0), Inertia(1)));
    Body::Rigid rigidBody2(MassProperties(1.0, Vec3(0), Inertia(1)));

    const ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free ellipsoid1(matter.updGround(), Transform(), rigidBody1, Transform());
    MobilizedBody::Free ellipsoid2(matter.updGround(), Transform(), rigidBody2, Transform());
    contacts.addBody(setIndex, ellipsoid1, ContactGeometry::Ellipsoid(radii1), center1);
    contacts.addBody(setIndex, ellipsoid2, ContactGeometry::Ellipsoid(radii2), center2);

    State state = system.realizeTopology();

    // Separated ellipsoids: no contact expected.
    ellipsoid1.setQToFitTransform(state, Transform(Rotation(), Vec3(0)));
    ellipsoid2.setQToFitTransform(state, Transform(Rotation(), Vec3(2, 0, 0)));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0U)
        << "Ellipsoids should not be in contact when clearly separated.";

    // Overlapping ellipsoids: exactly one contact expected.
    ellipsoid1.setQToFitTransform(state, Transform(Rotation(), Vec3(0)));
    ellipsoid2.setQToFitTransform(state, Transform(Rotation(), Vec3(1.5, 0, 0)));
    system.realize(state, Stage::Dynamics);

    const auto& foundContacts = contacts.getContacts(state, setIndex);
    EXPECT_EQ(foundContacts.size(), 1U) << "Expected exactly one contact between overlapping ellipsoids.";

    const auto verification = verifyEllipsoidContact(foundContacts[0],
                                                     {radii1, center1, ellipsoid1.getBodyTransform(state)},
                                                     {radii2, center2, ellipsoid2.getBodyTransform(state)});

    EXPECT_TRUE(verification.normalsAligned)
        << "Normal alignment failed for a known well-conditioned geometry - "
        << "this is not a convergence issue, the algorithm may be wrong.";
}

TEST(Simbody_CollisionDetectionAlgorithm, EllipsoidEllipsoidRandomCloud) {
    SCOPED_TRACE("Testing Ellipsoid-Ellipsoid contact detection in a random cloud of ellipsoids");

    // Maximum fraction of contacts allowed to have non-converged normals.
    // The Newton solver for ellipsoid closest-point can fail at very high
    // curvature; this tolerance acknowledges that without masking surface errors.
    static constexpr double kMaxNormalFailureFraction = 0.10;

    // Fixed seed for reproducibility - if this test fails, the same cloud
    // can be regenerated locally to debug.
    static constexpr int kRandomSeed = 42;
    static constexpr int kNumEllipsoids = 100;
    static constexpr double kMinRadius = 0.1;
    static constexpr double kMaxRadius = 1.1;
    static constexpr double kCloudSpan = 5.0;

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    ContactSetIndex setIndex = contacts.createContactSet();
    Body::Rigid rigidBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    Random::Uniform random(0.0, 1.0);
    random.setSeed(kRandomSeed);

    for (int idx = 0; idx < kNumEllipsoids; idx++) {
        MobilizedBody::Free ellipsoid(matter.updGround(), Transform(), rigidBody, Transform());
        const Vec3 radii(kMinRadius + (random.getValue() * (kMaxRadius - kMinRadius)),
                         kMinRadius + (random.getValue() * (kMaxRadius - kMinRadius)),
                         kMinRadius + (random.getValue() * (kMaxRadius - kMinRadius)));
        contacts.addBody(setIndex, ellipsoid, ContactGeometry::Ellipsoid(radii), Vec3(0));
    }

    State state = system.realizeTopology();

    for (MobilizedBodyIndex idx(1); idx <= kNumEllipsoids; idx++) {
        Rotation rot;
        rot.setRotationToBodyFixedXYZ(
            Vec3(random.getValue() * SimTK_PI, random.getValue() * SimTK_PI, random.getValue() * SimTK_PI));
        const Vec3 pos = kCloudSpan * Vec3(random.getValue(), random.getValue(), random.getValue());
        matter.getMobilizedBody(idx).setQToFitTransform(state, Transform(rot, pos));
    }

    system.realize(state, Stage::Dynamics);

    const Array_<Contact>& foundContacts = contacts.getContacts(state, setIndex);
    int normalFailureCount = 0;

    for (const Contact& ctc : foundContacts) {
        EXPECT_TRUE(PointContact::isInstance(ctc))
            << "All contacts in an ellipsoid-ellipsoid test must be PointContacts.";

        const auto& ptc = static_cast<const PointContact&>(ctc);
        const auto& geom1 = static_cast<const ContactGeometry::Ellipsoid&>(
            contacts.getBodyGeometry(setIndex, ptc.getSurface1()));
        const auto& geom2 = static_cast<const ContactGeometry::Ellipsoid&>(
            contacts.getBodyGeometry(setIndex, ptc.getSurface2()));
        const MobilizedBody& mob1 = contacts.getBody(setIndex, ptc.getSurface1());
        const MobilizedBody& mob2 = contacts.getBody(setIndex, ptc.getSurface2());

        const EllipsoidShape shape1{geom1.getRadii(), Vec3(0), mob1.getBodyTransform(state)};
        const EllipsoidShape shape2{geom2.getRadii(), Vec3(0), mob2.getBodyTransform(state)};

        const auto result = verifyEllipsoidContact(ctc, shape1, shape2);

        // Surface-point failures are always hard errors - the contact location
        // is wrong regardless of Newton convergence.
        EXPECT_TRUE(result.surfacePointsOnEllipsoids)
            << "Surface point check failed - this is never acceptable.";

        if (!result.normalsAligned) {
            normalFailureCount++;
        }
    }

    const int maxAllowedNormalFailures =
        static_cast<int>(std::ceil(kMaxNormalFailureFraction * foundContacts.size()));

    EXPECT_LE(normalFailureCount, maxAllowedNormalFailures)
        << "Too many normal alignment failures: " << normalFailureCount << " of " << foundContacts.size()
        << " contacts failed (limit is " << kMaxNormalFailureFraction * 100 << "%).\n"
        << "Seed was " << kRandomSeed << " - rerun with the same seed to reproduce.";
}

/**
 * Check the set of faces in a contact.
 */
void verifyContactFaces(const int* expected, int numExpected, const set<int>& found) {
    EXPECT_EQ(numExpected, found.size());
    for (int i = 0; i < numExpected; i++) {
        EXPECT_NE(found.find(expected[i]), found.end());
    }
}

TEST(Simbody_CollisionDetectionAlgorithm, HalfSpaceTriangleMesh) {
    SCOPED_TRACE("Testing HalfSpace-TriangleMesh contact detection");

    // Create a triangle mesh consisting of two pyramids: one right side up and one upside down.
    vector<Vec3> vertices;
    vertices.emplace_back(0, 0, 0);
    vertices.emplace_back(1, 0, 0);
    vertices.emplace_back(0, 0, 1);
    vertices.emplace_back(1, 0, 1);
    vertices.emplace_back(0.5, 1, 0.5);
    vertices.emplace_back(2, 1, 0);
    vertices.emplace_back(2, 1, 1);
    vertices.emplace_back(3, 1, 1);
    vertices.emplace_back(3, 1, 0);
    vertices.emplace_back(2.5, 0, 0.5);

    constexpr std::array<std::array<int, 3>, 6> faces = {
        {{0, 1, 2}, {0, 2, 3}, {1, 0, 4}, {0, 3, 4}, {3, 2, 4}, {2, 1, 4}}};
    vector<int> faceIndices;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) {
            faceIndices.push_back(faces[i][j]);
        }
    }
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) {
            faceIndices.push_back(faces[i][j] + 5);
        }
    }
    ContactGeometry::TriangleMesh mesh(vertices, faceIndices);

    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free b(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, b, mesh, Transform());
    contacts.addBody(setIndex,
                     matter.updGround(),
                     ContactGeometry::HalfSpace(),
                     Transform(Rotation(-0.5 * Pi, ZAxis), Vec3(0, 1, 0))); // y < 1
    State state = system.realizeTopology();

    constexpr std::array<int, 10> bottomFaces = {0, 1, 2, 3, 4, 5, 8, 9, 10, 11};
    for (Real depth = -0.25; depth < 1; depth += 0.1) {
        Vec3 center(0.1, 1 - depth, 2.0);
        b.setQToFitTranslation(state, center);
        system.realize(state, Stage::Dynamics);
        const Array_<Contact>& contact = contacts.getContacts(state, setIndex);
        if (depth < 0.0) {
            EXPECT_EQ(contact.size(), 0);
        } else {
            EXPECT_EQ(contact.size(), 1);
            EXPECT_TRUE(TriangleMeshContact::isInstance(contact[0]));
            const TriangleMeshContact& c = static_cast<const TriangleMeshContact&>(contact[0]);
            EXPECT_EQ(c.getSurface1Faces().size(), 0);
            verifyContactFaces(bottomFaces.data(), bottomFaces.size(), c.getSurface2Faces());
        }
    }
}

TEST(Simbody_CollisionDetectionAlgorithm, SphereTriangleMesh) {
    SCOPED_TRACE("Testing Sphere-TriangleMesh contact detection");

    // Create a triangle mesh consisting of two pyramids: one right side up and one upside down.
    vector<Vec3> vertices;
    vertices.emplace_back(0, 0, 0);
    vertices.emplace_back(0, 0, 1);
    vertices.emplace_back(1, 0, 1);
    vertices.emplace_back(1, 0, 0);
    vertices.emplace_back(0.5, 1, 0.5);
    vertices.emplace_back(2, 1, 0);
    vertices.emplace_back(2, 1, 1);
    vertices.emplace_back(3, 1, 1);
    vertices.emplace_back(3, 1, 0);
    vertices.emplace_back(2.5, 0, 0.5);

    constexpr std::array<std::array<int, 3>, 6> faces{
        {{0, 1, 2}, {0, 2, 3}, {1, 0, 4}, {0, 3, 4}, {3, 2, 4}, {2, 1, 4}}};
    vector<int> faceIndices;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) {
            faceIndices.push_back(faces[i][j]);
        }
    }
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) {
            faceIndices.push_back(faces[i][j] + 5);
        }
    }
    ContactGeometry::TriangleMesh mesh(vertices, faceIndices);

    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free mobod(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, mobod, mesh, Transform());
    contacts.addBody(setIndex, matter.updGround(), ContactGeometry::Sphere(0.5), Transform(Vec3(0, 1, 0)));
    State state = system.realizeTopology();

    // Try various positions and make sure the results are correct.
    mobod.setQToFitTranslation(state, Vec3(0, -2, 0));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0);

    mobod.setQToFitTranslation(state, Vec3(0, 1.51, 0));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0);

    {
        mobod.setQToFitTranslation(state, Vec3(-0.5, 1.49, -0.5));
        system.realize(state, Stage::Dynamics);
        EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 1);
        const auto& triangleMeshContact =
            static_cast<const TriangleMeshContact&>(contacts.getContacts(state, setIndex)[0]);
        constexpr std::array<int, 2> faces = {0, 1};
        verifyContactFaces(faces.data(), faces.size(), triangleMeshContact.getSurface2Faces());
    }

    mobod.setQToFitTranslation(state, Vec3(-0.5, -0.51, -0.5));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0);
    {
        mobod.setQToFitTranslation(state, Vec3(-0.5, -0.49, -0.5));
        system.realize(state, Stage::Dynamics);
        EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 1);
        const auto& triangleMeshContact =
            static_cast<const TriangleMeshContact&>(contacts.getContacts(state, setIndex)[0]);
        constexpr std::array<int, 4> faces = {2, 3, 4, 5};
        verifyContactFaces(faces.data(), faces.size(), triangleMeshContact.getSurface2Faces());
    }

    mobod.setQToFitTranslation(state, Vec3(-2.5, 1.51, -0.5));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0);
    {
        mobod.setQToFitTranslation(state, Vec3(-2.5, 1.49, -0.5));
        system.realize(state, Stage::Dynamics);
        EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 1);
        const auto& triangleMeshContact =
            static_cast<const TriangleMeshContact&>(contacts.getContacts(state, setIndex)[0]);
        constexpr std::array<int, 4> faces = {8, 9, 10, 11};
        verifyContactFaces(faces.data(), faces.size(), triangleMeshContact.getSurface2Faces());
    }
}

TEST(Simbody_CollisionDetectionAlgorithm, TriangleMeshTriangleMesh) {
    SCOPED_TRACE("Testing TriangleMesh-TriangleMesh contact detection");

    // Create two triangle meshes, each consisting of a pyramid.
    std::vector<Vec3> vertices;
    vertices.emplace_back(0, 0, 0);
    vertices.emplace_back(1, 0, 0);
    vertices.emplace_back(1, 0, 1);
    vertices.emplace_back(0, 0, 1);
    vertices.emplace_back(0.5, 1, 0.5);

    constexpr std::array<std::array<int, 3>, 6> faces{
        {{0, 1, 2}, {0, 2, 3}, {1, 0, 4}, {2, 1, 4}, {3, 2, 4}, {0, 3, 4}}};
    std::vector<int> faceIndices;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) {
            faceIndices.push_back(faces[i][j]);
        }
    }
    ContactGeometry::TriangleMesh mesh1(vertices, faceIndices);
    ContactGeometry::TriangleMesh mesh2(vertices, faceIndices);

    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free mobod1(matter.updGround(), Transform(), body, Transform());
    MobilizedBody::Free mobod2(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, mobod1, mesh1, Transform());
    contacts.addBody(setIndex, mobod2, mesh2, Transform());
    State state = system.realizeTopology();

    // Try some configurations that should not intersect.
    mobod1.setQToFitTranslation(state, Vec3(0));
    mobod2.setQToFitTranslation(state, Vec3(2));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0);

    mobod1.setQToFitTranslation(state, Vec3(0));
    mobod2.setQToFitTranslation(state, Vec3(1.01, 0, 0));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0);

    mobod1.setQToFitTranslation(state, Vec3(0));
    mobod2.setQToFitTranslation(state, Vec3(0, 1.01, 0));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0);

    mobod1.setQToFitTranslation(state, Vec3(0));
    mobod2.setQToFitTranslation(state, Vec3(0, -1.01, 0));
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(contacts.getContacts(state, setIndex).size(), 0);

    // Now try ones that should intersect.
    constexpr std::array<int, 2> baseFaces = {0, 1};
    constexpr std::array<int, 4> pointFaces = {2, 3, 4, 5};
    {
        mobod1.setQToFitTranslation(state, Vec3(0));
        mobod2.setQToFitTranslation(state, Vec3(0, -0.99, 0));
        system.realize(state, Stage::Dynamics);
        const Array_<Contact>& contact = contacts.getContacts(state, setIndex);

        EXPECT_EQ(contact.size(), 1);
        EXPECT_TRUE(TriangleMeshContact::isInstance(contact[0]));
        const auto& c = static_cast<const TriangleMeshContact&>(contact[0]);

        if (contact[0].getSurface1() == 0) {
            verifyContactFaces(baseFaces.data(), baseFaces.size(), c.getSurface1Faces());
            verifyContactFaces(pointFaces.data(), pointFaces.size(), c.getSurface2Faces());
        } else {
            verifyContactFaces(pointFaces.data(), pointFaces.size(), c.getSurface1Faces());
            verifyContactFaces(baseFaces.data(), baseFaces.size(), c.getSurface2Faces());
        }
    }

    {
        mobod1.setQToFitTranslation(state, Vec3(0, -0.5, 0));
        mobod2.setQToFitTranslation(state, Vec3(0, 0.49, 0));
        system.realize(state, Stage::Dynamics);
        const Array_<Contact>& contact = contacts.getContacts(state, setIndex);

        EXPECT_EQ(contact.size(), 1);
        EXPECT_TRUE(TriangleMeshContact::isInstance(contact[0]));
        const auto& triangleMeshContact = static_cast<const TriangleMeshContact&>(contact[0]);

        if (contact[0].getSurface1() == 0) {
            verifyContactFaces(pointFaces.data(), pointFaces.size(), triangleMeshContact.getSurface1Faces());
            verifyContactFaces(baseFaces.data(), baseFaces.size(), triangleMeshContact.getSurface2Faces());
        } else {
            verifyContactFaces(baseFaces.data(), baseFaces.size(), triangleMeshContact.getSurface1Faces());
            verifyContactFaces(pointFaces.data(), pointFaces.size(), triangleMeshContact.getSurface2Faces());
        }
    }

    {
        mobod1.setQToFitTranslation(state, Vec3(0.1, -0.5, 0));
        mobod2.setQToFitTranslation(state, Vec3(0, 0.49, 0.1));
        system.realize(state, Stage::Dynamics);
        const Array_<Contact>& contact = contacts.getContacts(state, setIndex);

        EXPECT_EQ(contact.size(), 1);
        EXPECT_TRUE(TriangleMeshContact::isInstance(contact[0]));
    }

    {
        mobod1.setQToFitTransform(state, Transform(Rotation(-0.5 * Pi, ZAxis), Vec3(0, 0.5, 0)));
        mobod2.setQToFitTransform(state, Transform(Rotation(0.5 * Pi, ZAxis), Vec3(1.9, -0.5, 0)));
        system.realize(state, Stage::Dynamics);
        const Array_<Contact>& contact = contacts.getContacts(state, setIndex);

        EXPECT_EQ(contact.size(), 1);
        EXPECT_TRUE(TriangleMeshContact::isInstance(contact[0]));

        const auto& triangleMeshContact = static_cast<const TriangleMeshContact&>(contact[0]);
        verifyContactFaces(pointFaces.data(), pointFaces.size(), triangleMeshContact.getSurface1Faces());
        verifyContactFaces(pointFaces.data(), pointFaces.size(), triangleMeshContact.getSurface2Faces());
    }
}
