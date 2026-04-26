#include <gtest/gtest.h>

#include "SimTKsimbody.h"

using namespace SimTK;

// Define a tolerance for floating point comparisons if not already defined
const double TOL = 1e-10;

#define EXPECT_NEAR_CUSTOM_TOL_SIMTK(expected, actual, tol)                              \
    EXPECT_TRUE(SimTK::Test::numericallyEqual((expected), (actual), /*scale*/ 1, (tol))) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\nTolerance: " << (tol) << "\n"

#define EXPECT_NEAR_DEFAULT_TOL_SIMTK(expected, actual) \
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((expected), (actual), 1e-6)

TEST(Simbody_OrientedBoundingBox_ContainsPoint, AxisAlignedAtOrigin) {
    OrientedBoundingBox box(Vec3(0), Vec3(1, 2, 3));

    // Boundary and interior points
    EXPECT_TRUE(box.containsPoint(Vec3(0)));
    EXPECT_TRUE(box.containsPoint(Vec3(1, 2, 3)));
    EXPECT_TRUE(box.containsPoint(Vec3(0.5, 1.0, 1.5)));

    // Points outside each face (negative and positive directions)
    EXPECT_FALSE(box.containsPoint(Vec3(-0.5, 1.0, 1.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(0.5, -1.0, 1.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(0.5, 1.0, -1.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(1.01, 1.0, 1.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(0.5, 2.01, 1.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(0.5, 1.0, 3.01)));
}

TEST(Simbody_OrientedBoundingBox_ContainsPoint, AxisAlignedTranslated) {
    OrientedBoundingBox box(Vec3(3, 2, 1), Vec3(1, 2, 3));

    // Boundary and interior points
    EXPECT_TRUE(box.containsPoint(Vec3(3, 2, 1)));
    EXPECT_TRUE(box.containsPoint(Vec3(4, 4, 4)));
    EXPECT_TRUE(box.containsPoint(Vec3(3.5, 3.0, 2.5)));

    // Exterior points
    EXPECT_FALSE(box.containsPoint(Vec3(2.5, 3.0, 2.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(4.5, 1.0, 2.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(4.5, 3.0, -0.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(5.01, 3.0, 2.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(3.5, 4.01, 2.5)));
    EXPECT_FALSE(box.containsPoint(Vec3(3.5, 3.0, 4.01)));
}

TEST(Simbody_OrientedBoundingBox_ContainsPoint, RotatedOBB) {
    // Rotation of 45 degrees (Pi/4) around Z
    OrientedBoundingBox box(Rotation(0.25 * Pi, ZAxis), Vec3(1, 1, 1));

    EXPECT_TRUE(box.containsPoint(Vec3(0)));

    // Test point near the rotated corner/edge
    // Sqrt2 is typically defined in Simbody, or use std::sqrt(2.0)
    EXPECT_TRUE(box.containsPoint(Vec3(0, Sqrt2 - 1e-10, 0)));
    EXPECT_FALSE(box.containsPoint(Vec3(0, Sqrt2 + 1e-10, 0)));

    // Z-axis checks
    EXPECT_TRUE(box.containsPoint(Vec3(0, Sqrt2 - 1e-10, 1e-10)));
    EXPECT_FALSE(box.containsPoint(Vec3(0, Sqrt2 - 1e-10, -1e-10)));
}

/**
 * @brief Helper to verify that the generated corners match the expected corners,
 * regardless of the order in which they are returned.
 */
void verifyCorners(const Vec3 expected[8], const Vec3 found[8]) {
    for (int i = 0; i < 8; i++) {
        bool match = false;
        for (int j = 0; j < 8 && !match; j++) {
            if (std::abs(expected[i][0] - found[j][0]) < TOL && std::abs(expected[i][1] - found[j][1]) < TOL
                && std::abs(expected[i][2] - found[j][2]) < TOL) {
                match = true;
                break;
            }
        }
        EXPECT_TRUE(match) << "Could not find expected corner: " << expected[i];
    }
}

/**
 * @brief Tests the corner generation logic for OrientedBoundingBox (OBB).
 * This covers axis-aligned boxes at origin, translated boxes, and rotated boxes.
 */
TEST(Simbody_OrientedBoundingBox_GetCorners, CorrectlyCalculatesAllEightVertices) {
    Vec3 corners[8];

    // Case 1: Axis-aligned OBB at origin
    {
        SCOPED_TRACE("Axis-aligned Oriented Bounding Box at origin");
        OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)).getCorners(corners);
        Vec3 expected[] = {Vec3(0),
                           Vec3(1, 0, 0),
                           Vec3(0, 2, 0),
                           Vec3(1, 2, 0),
                           Vec3(0, 0, 3),
                           Vec3(1, 0, 3),
                           Vec3(0, 2, 3),
                           Vec3(1, 2, 3)};
        verifyCorners(expected, corners);
    }

    // Case 2: Axis-aligned OBB translated from origin
    {
        SCOPED_TRACE("Axis-aligned Oriented Bounding Box translated from origin");
        OrientedBoundingBox(Vec3(3, 2, 1), Vec3(1, 2, 3)).getCorners(corners);
        Vec3 expected[] = {Vec3(3, 2, 1),
                           Vec3(4, 2, 1),
                           Vec3(3, 4, 1),
                           Vec3(4, 4, 1),
                           Vec3(3, 2, 4),
                           Vec3(4, 2, 4),
                           Vec3(3, 4, 4),
                           Vec3(4, 4, 4)};
        verifyCorners(expected, corners);
    }

    // Case 3: OBB rotated 45 degrees around Z-axis
    {
        SCOPED_TRACE("Oriented Bounding Box rotated 45 degrees around Z-axis");
        OrientedBoundingBox(Rotation(0.25 * Pi, ZAxis), Vec3(1, 1, 1)).getCorners(corners);
        Real d = 0.5 * Sqrt2;
        Vec3 expected[] = {Vec3(0),
                           Vec3(-d, d, 0),
                           Vec3(d, d, 0),
                           Vec3(0, Sqrt2, 0),
                           Vec3(0, 0, 1),
                           Vec3(-d, d, 1),
                           Vec3(d, d, 1),
                           Vec3(0, Sqrt2, 1)};
        verifyCorners(expected, corners);
    }
}

/**
 * @brief Helper to verify if two boxes intersect.
 * Checks commutativity (A vs B and B vs A) as required by Simbody tests.
 */
void verifyBoxIntersection(bool shouldIntersect,
                           const OrientedBoundingBox& box1,
                           const OrientedBoundingBox& box2) {
    EXPECT_EQ(box1.intersectsBox(box2), shouldIntersect)
        << "Intersection check failed for box1.intersectsBox(box2)";
    EXPECT_EQ(box2.intersectsBox(box1), shouldIntersect)
        << "Intersection check failed for box2.intersectsBox(box1) (Commutativity)";
}

/**
 * @brief Tests intersection between boxes with identical (axis-aligned) orientations.
 */
TEST(Simbody_OrientedBoundingBox_IntersectsBox, IdenticalOrientation) {
    // Overlapping at origin
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(Vec3(0), Vec3(3, 2, 1)));

    // Just barely touching/overlapping
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(Vec3(1 - 1e-10, 2 - 1e-10, 3 - 1e-10), Vec3(3, 2, 1)));

    // Separated by epsilon in each axis (Positive direction)
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(Vec3(1 + 1e-10, 2 - 1e-10, 3 - 1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(Vec3(1 - 1e-10, 2 + 1e-10, 3 - 1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(Vec3(1 - 1e-10, 2 - 1e-10, 3 + 1e-10), Vec3(3, 2, 1)));

    // Just barely touching/overlapping (Negative direction)
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(Vec3(-3 + 1e-10, -2 + 1e-10, -1 + 1e-10), Vec3(3, 2, 1)));

    // Separated by epsilon in each axis (Negative direction)
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(Vec3(-3 - 1e-10, -2 + 1e-10, -1 + 1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(Vec3(-3 + 1e-10, -2 - 1e-10, -1 + 1e-10), Vec3(3, 2, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(Vec3(-3 + 1e-10, -2 + 1e-10, -1 - 1e-10), Vec3(3, 2, 1)));
}

/**
 * @brief Tests intersection between boxes rotated by 90 degrees relative to each other.
 */
TEST(Simbody_OrientedBoundingBox_IntersectsBox, Rotated90Degrees) {
    const Rotation R_Z90(0.5 * Pi, ZAxis);

    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(-0.5, 0, 0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(R_Z90, Vec3(3, 2, 1)));
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(-0.5, -2 + 1e-10, 0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(R_Z90, Vec3(3, 2, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(-0.5, -2 - 1e-10, 0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(R_Z90, Vec3(3, 2, 1)));
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(-0.5, 3 - 1e-10, 0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(R_Z90, Vec3(3, 2, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(-0.5, 3 + 1e-10, 0), Vec3(1, 2, 3)),
                          OrientedBoundingBox(R_Z90, Vec3(3, 2, 1)));
    verifyBoxIntersection(true,
                          OrientedBoundingBox(R_Z90, Vec3(1, 2, 3)),
                          OrientedBoundingBox(R_Z90, Vec3(3, 2, 1)));
}

/**
 * @brief Tests intersection between boxes rotated by 45 degrees.
 */
TEST(Simbody_OrientedBoundingBox_IntersectsBox, Rotated45Degrees) {
    const Real d = 0.5 * Sqrt2;
    const Rotation R_Z45(0.25 * Pi, ZAxis);
    const Rotation R_Xn45(-0.25 * Pi, XAxis);

    // Z-Axis rotation tests
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(-0.5, 0, 0), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Z45, Vec3(1, 1, 1)));
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(-0.5, -1 + 1e-10, 0), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Z45, Vec3(1, 1, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(-0.5, -1 - 1e-10, 0), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Z45, Vec3(1, 1, 1)));
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(-0.5, Sqrt2 - 1e-10, 0), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Z45, Vec3(1, 1, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(-0.5, Sqrt2 + 1e-10, 0), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Z45, Vec3(1, 1, 1)));

    // Translation tests against rotated box
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(d - 1e-10, 0, 0), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Z45, Vec3(1, 1, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(d + 1e-10, 0, 0), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Z45, Vec3(1, 1, 1)));
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(-1 - d + 1e-10, 0, 0), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Z45, Vec3(1, 1, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(-1 - d - 1e-10, 0, 0), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Z45, Vec3(1, 1, 1)));

    // X-Axis rotation tests
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(0, 0, d - 1e-10), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Xn45, Vec3(1, 1, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(0, 0, d + 1e-10), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Xn45, Vec3(1, 1, 1)));
    verifyBoxIntersection(true,
                          OrientedBoundingBox(Vec3(0, 0, -1 - d + 1e-10), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Xn45, Vec3(1, 1, 1)));
    verifyBoxIntersection(false,
                          OrientedBoundingBox(Vec3(0, 0, -1 - d - 1e-10), Vec3(1, 1, 1)),
                          OrientedBoundingBox(R_Xn45, Vec3(1, 1, 1)));
}

/**
 * @brief Regression test for a specific intersection bug involving
 * Euler angle transformations.
 */
TEST(Simbody_OrientedBoundingBox_IntersectsBox, KnownRegressionBug) {
    Rotation r;
    r.setRotationToBodyFixedXYZ(Vec3(Pi / 2, 0, -1.1));

    OrientedBoundingBox box1(Transform(r, Vec3(-0.95, 1.1, 2.5)), Vec3(2e-10, 1.118, 1));
    OrientedBoundingBox box2(Vec3(0, -50, -50), Vec3(100, 100, 100));

    verifyBoxIntersection(true, box1, box2);
}

void verifyRayIntersection(const OrientedBoundingBox& box,
                           const Vec3& origin,
                           const UnitVec3& direction,
                           bool shouldIntersect,
                           Real expectedDistance) {
    Real distance;
    EXPECT_TRUE(shouldIntersect == box.intersectsRay(origin, direction, distance));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(expectedDistance, distance);
}

void testIntersectsRay() {
    // Try rays starting inside the box.
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)),
                          Vec3(0.5, 0.5, 0.5),
                          UnitVec3(1, 0, 0),
                          true,
                          0);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)),
                          Vec3(0.5, 0.5, 0.5),
                          UnitVec3(0, 1, 0),
                          true,
                          0);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)),
                          Vec3(0.5, 0.5, 0.5),
                          UnitVec3(0, 0, 1),
                          true,
                          0);

    // Try rays that hit it on various sides.
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)),
                          Vec3(-1.5, 0.5, 0.5),
                          UnitVec3(1, 0, 0),
                          true,
                          1.5);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)),
                          Vec3(2.5, 0.5, 0.5),
                          UnitVec3(-1, 0, 0),
                          true,
                          1.5);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)),
                          Vec3(0.5, -0.5, 0.5),
                          UnitVec3(0, 1, 0),
                          true,
                          0.5);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)),
                          Vec3(0.5, 2.5, 0.5),
                          UnitVec3(0, -1, 0),
                          true,
                          0.5);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)),
                          Vec3(0.5, 0.5, -1.0),
                          UnitVec3(0, 0, 1),
                          true,
                          1.0);
    verifyRayIntersection(OrientedBoundingBox(Transform(), Vec3(1, 2, 3)),
                          Vec3(0.5, 0.5, 4.0),
                          UnitVec3(0, 0, -1),
                          true,
                          1.0);

    // Try a box at an angle.
    verifyRayIntersection(OrientedBoundingBox(Rotation(-0.25 * Pi, ZAxis), Vec3(2, 2, 2)),
                          Vec3(0, 0, 0.5),
                          UnitVec3(1, 0, 0),
                          true,
                          0);
    verifyRayIntersection(OrientedBoundingBox(Rotation(-0.25 * Pi, ZAxis), Vec3(2, 2, 2)),
                          Vec3(-1, 0, 0.5),
                          UnitVec3(1, 0, 0),
                          true,
                          1.0);
}

/**
 * @brief Stress test for OrientedBoundingBox creation from a point cloud.
 * * This test generates random volumes, fills them with random points,
 * and verifies that the resulting OBB:
 * 1. Strictly contains every point used to create it.
 * 2. Is reasonably efficient (the volume is within 50% of the generating volume).
 */
TEST(Simbody_OrientedBoundingBox_CreateFromPoints, RandomPointClouds) {
    // Simbody's random number generator
    Random::Uniform random(0, 1);

    for (int trial = 0; trial < 100; trial++) {
        // 1. Define a random volume (size, orientation, and position)
        Vec3 size(10 * random.getValue(), 10 * random.getValue(), 10 * random.getValue());

        Rotation rotation;
        rotation.setRotationToBodyFixedXYZ(Vec3(random.getValue(), random.getValue(), random.getValue()));

        Transform transform(rotation,
                            Vec3(10 * random.getValue(), 10 * random.getValue(), 10 * random.getValue()));

        // 2. Generate a random set of points within that volume
        int numPoints = static_cast<int>((50 * random.getValue()) + 1);
        Vector_<Vec3> points(numPoints);
        for (int i = 0; i < numPoints; i++) {
            Vec3 localPoint(size[0] * random.getValue(),
                            size[1] * random.getValue(),
                            size[2] * random.getValue());
            points[i] = transform * localPoint;
        }

        // 3. Construct the OBB from the point set
        OrientedBoundingBox box(points);

        // 4. Verification: All points must be inside the box
        for (int i = 0; i < numPoints; i++) {
            ASSERT_TRUE(box.containsPoint(points[i])) << "Trial " << trial << " failed: Point " << i << " ("
                                                      << points[i] << ") is outside the generated OBB.";
        }

        // 5. Verification: The box should be reasonably tight
        const Real expectedVolume = size[0] * size[1] * size[2];
        const Vec3 boxSize = box.getSize();
        const Real actualVolume = boxSize[0] * boxSize[1] * boxSize[2];

        // We use EXPECT_LT here so that even if one trial is slightly over
        // the volume limit, we continue to check other trials.
        EXPECT_LT(actualVolume, 1.5 * expectedVolume)
            << "Trial " << trial << " failed: OBB volume (" << actualVolume
            << ") exceeds 150% of generating volume (" << expectedVolume << ")";
    }
}

/**
 * @brief Tests findNearestPoint for points already inside the OBB.
 * The nearest point to an internal point should be the point itself.
 */
TEST(Simbody_OrientedBoundingBox_FindNearestPoint, PointsInsideBox) {
    Vec3 size(1, 1.5, 3);
    Transform trans(Rotation(0.3, XAxis), Vec3(1, 2, 0.5));
    OrientedBoundingBox box(trans, size);

    Random::Uniform random(0, 1);
    for (int i = 0; i < 100; i++) {
        // Generate a point in local space [0, size]
        Vec3 localP(random.getValue() * size[0], random.getValue() * size[1], random.getValue() * size[2]);

        // Transform to global space
        Vec3 globalP = trans * localP;
        Vec3 nearest = box.findNearestPoint(globalP);

        // Check equality within tolerance for all 3 coordinates
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(globalP[j], nearest[j], TOL);
        }
    }
}

/**
 * @brief Tests findNearestPoint for points outside the OBB.
 * Verifies that the nearest point is correctly clamped to the box boundaries.
 */
TEST(Simbody_OrientedBoundingBox_FindNearestPoint, PointsOutsideBox) {
    Vec3 size(1, 1.5, 3);
    Transform trans(Rotation(0.3, XAxis), Vec3(1, 2, 0.5));
    OrientedBoundingBox box(trans, size);

    // Helper lambda to check local coordinates
    auto checkNearestLocal = [&](const Vec3& globalTarget, const Vec3& localExpected) {
        Vec3 globalNearest = box.findNearestPoint(globalTarget);
        Vec3 localNearest = ~trans * globalNearest; // Move to box local frame
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(localNearest[i], localExpected[i], TOL);
        }
    };

    // Test 1: Point beyond the maximum corner
    checkNearestLocal(trans * (size + Vec3(1, 2, 3)), size);

    // Test 2: Point outside X and Y, but inside Z range
    checkNearestLocal(trans * Vec3(2, 3, 0.25), Vec3(1, 1.5, 0.25));

    // Test 3: Point beyond the minimum corner (origin)
    checkNearestLocal(trans * Vec3(-1, -1, -2), Vec3(0, 0, 0));

    // Test 4: Point outside X and Y min, but inside Z range
    checkNearestLocal(trans * Vec3(-1, -1, 0.5), Vec3(0, 0, 0.5));
}