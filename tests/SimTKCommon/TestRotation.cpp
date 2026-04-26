/**
 * TestRotation.cpp
 *
 * Google Test conversion of the original SimBody Rotation CTest suite.
 *
 * What the original authors tested
 * ---------------------------------
 * 1. Basic Rotation construction: default (identity), copy, assignment,
 *    and construction from a CoordinateAxis with a scalar angle.
 * 2. Single-axis rotations about X/Y/Z: the angle-about-axis constructor
 *    must match the equivalent unit-vector form, and the inverse
 *    (convertOneAxisRotationToOneAngle) must recover the original angle
 *    within machine precision when the angle is in [-π, π].
 * 3. Two-axis body/space rotation sequences: the compound constructor must
 *    match manual multiplication of individual rotations, and the inverse
 *    must recover both angles.
 * 4. Three-axis body/space rotation sequences: same as above for all 27
 *    axis triples, including previously-failing near-π edge cases and
 *    exhaustive near-singularity sweeps.
 * 5. Quaternion construction, normalization, and round-trip through
 *    Rotation: quaternion → Rotation → quaternion must round-trip cleanly,
 *    and both normalize() and normalizeThis() must produce unit norm.
 * 6. Construction of a nearby orthogonal rotation matrix from a raw Mat33.
 * 7. Special-cased body-fixed XYZ setter vs. the generic three-axis setter.
 * 8. Rotation from two given axes: the constructor must produce determinant
 *    +1 and both column axes must be oriented in the requested direction.
 * 9. Re-expression of a symmetric 3×3 matrix under rotation must match the
 *    full Mat33 triple-product, be invertible, and be unchanged by identity.
 * 10. Complex-valued SymMat Hermitian symmetry and multiplication correctness.
 */

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>

#include "SimTKcommon/Testing.h"

#include "SimTKcommon.h"
#include "Util.hpp"

using namespace SimTK;

// ============================================================
// Macro helpers
// ============================================================

#define EXPECT_SIMTK_EQ(actual, expected) \
    EXPECT_TRUE(AssertSimTKEqual(#actual, #expected, (actual), (expected)))

// ============================================================
// Internal helpers
// ============================================================

/// Returns true when the 1-angle inverse round-trip error is within
/// 10×SignificantReal, or when the original angle lies outside [-π, π]
/// (where uniqueness is not guaranteed).
static auto inverseAngle1IsValid(Real angle, Real theta) -> bool {
    const bool inRange = ((-SimTK_PI) <= angle && angle <= SimTK_PI);
    if (!inRange) {
        return true;
    }
    return std::fabs(angle - theta) < (10 * SignificantReal);
}

/// Returns true when both 2-angle inverse round-trips are within tolerance,
/// or when either input angle lies outside [-π, π].
static auto inverseAngle2IsValid(Real angle1, Real theta1, Real angle2, Real theta2) -> bool {
    const bool in1 = ((-SimTK_PI) <= angle1 && angle1 <= SimTK_PI);
    const bool in2 = ((-SimTK_PI) <= angle2 && angle2 <= SimTK_PI);
    if (!in1 || !in2) {
        return true;
    }
    const bool ok1 = std::fabs(angle1 - theta1) < (10 * SignificantReal);
    const bool ok2 = std::fabs(angle2 - theta2) < (10 * SignificantReal);
    return ok1 && ok2;
}

/// Validates the inverse for a two-axis (same first-and-last axis) sequence.
/// Near the singularity angle2 ≈ 0, only the sum angle1+angle3 is uniquely
/// recoverable.
static auto
inverseAngleTwoAxesIsValid(Real angle1, Real theta1, Real angle2, Real theta2, Real angle3, Real theta3)
    -> bool {
    const bool in1 = ((-SimTK_PI) <= angle1 && angle1 <= SimTK_PI);
    const bool in2 = (Real(0) <= angle2 && angle2 <= SimTK_PI);
    const bool in3 = ((-SimTK_PI) <= angle3 && angle3 <= SimTK_PI);
    if (!in1 || !in2 || !in3) {
        return true;
    }

    const bool ok1 = std::fabs(angle1 - theta1) < (10 * SignificantReal);
    const bool ok2 = std::fabs(angle2 - theta2) < (10 * SignificantReal);
    const bool ok3 = std::fabs(angle3 - theta3) < (10 * SignificantReal);
    if (ok1 && ok2 && ok3) {
        return true;
    }

    const Real singularity = Real(0);
    if (std::fabs(angle2 - singularity) <= SignificantReal) {
        const Real sumAngle = angle1 + angle3;
        const Real sumTheta = theta1 + theta3;
        const bool sumInRange = ((-SimTK_PI) <= sumAngle && sumAngle <= SimTK_PI);
        if (!sumInRange) {
            return true;
        }
        return (std::fabs(angle2 - theta2) < (10 * SignificantReal))
               && (std::fabs(sumAngle - sumTheta) < (10 * SignificantReal));
    }
    return false;
}

/// Validates the inverse for a three distinct-axis sequence.
/// Near the singularity angle2 ≈ π/2, only the sum angle1+angle3 is
/// uniquely recoverable.
static auto
inverseAngleThreeAxesIsValid(Real angle1, Real theta1, Real angle2, Real theta2, Real angle3, Real theta3)
    -> bool {
    const bool in1 = ((-SimTK_PI) <= angle1 && angle1 <= SimTK_PI);
    const bool in2 = ((-Real(0.5) * SimTK_PI) <= angle2 && angle2 <= (Real(0.5) * SimTK_PI));
    const bool in3 = ((-SimTK_PI) <= angle3 && angle3 <= SimTK_PI);
    if (!in1 || !in2 || !in3) {
        return true;
    }

    const bool ok1 = std::fabs(angle1 - theta1) < (10 * SignificantReal);
    const bool ok2 = std::fabs(angle2 - theta2) < (10 * SignificantReal);
    const bool ok3 = std::fabs(angle3 - theta3) < (10 * SignificantReal);
    if (ok1 && ok2 && ok3) {
        return true;
    }

    const Real singularity = Real(0.5) * SimTK_PI;
    if (std::fabs(angle2 - singularity) <= SignificantReal) {
        const Real sumAngle = angle1 + angle3;
        const Real sumTheta = theta1 + theta3;
        const bool sumInRange = ((-SimTK_PI) <= sumAngle && sumAngle <= SimTK_PI);
        if (!sumInRange) {
            return true;
        }
        return (std::fabs(angle2 - theta2) < (10 * SignificantReal))
               && (std::fabs(sumAngle - sumTheta) < (10 * SignificantReal));
    }
    return false;
}

// ============================================================
// Per-case rotation check helpers returning AssertionResult
// (used by both named and exhaustive tests)
// ============================================================

/// Verifies a single-axis rotation:
///  (a) Rotation(angle, axis) matches Rotation(angle, unitVector).
///  (b) Rebuilding from the recovered angle reproduces the matrix.
///  (c) The round-trip angle is within tolerance in the canonical range.
static auto checkOneAxisRotation(Real angle, const CoordinateAxis& axis) -> ::testing::AssertionResult {
    Rotation rotationSpecified;
    if (axis == XAxis) {
        rotationSpecified.setRotationFromAngleAboutX(angle);
    }
    if (axis == YAxis) {
        rotationSpecified.setRotationFromAngleAboutY(angle);
    }
    if (axis == ZAxis) {
        rotationSpecified.setRotationFromAngleAboutZ(angle);
    }

    const Real unitX = (axis == XAxis) ? Real(1) : Real(0);
    const Real unitY = (axis == YAxis) ? Real(1) : Real(0);
    const Real unitZ = (axis == ZAxis) ? Real(1) : Real(0);
    const UnitVec3 unitVector(unitX, unitY, unitZ);
    const Rotation testRotation(angle, unitVector);

    if (!rotationSpecified.areAllRotationElementsSameToMachinePrecision(testRotation)) {
        return ::testing::AssertionFailure() << "Rotation(angle, axis) != Rotation(angle, unitVector) "
                                             << "for angle=" << angle;
    }

    const Real theta = rotationSpecified.convertOneAxisRotationToOneAngle(axis);

    Rotation recovered;
    recovered.setRotationFromAngleAboutAxis(theta, axis);
    if (!rotationSpecified.areAllRotationElementsSameToMachinePrecision(recovered)) {
        return ::testing::AssertionFailure() << "Rotation rebuilt from recovered angle does not match "
                                             << "original for angle=" << angle;
    }

    if (!inverseAngle1IsValid(angle, theta)) {
        return ::testing::AssertionFailure()
               << "Round-trip angle error too large: angle=" << angle << " recovered=" << theta;
    }
    return ::testing::AssertionSuccess();
}

/// Verifies a two-axis rotation:
///  (a) The compound constructor matches manual multiplication.
///  (b) Rebuilding from recovered angles reproduces the matrix.
///  (c) Angle round-trips are within tolerance in the canonical range.
static auto checkTwoAxesRotation(BodyOrSpaceType bodyOrSpace,
                                 Real angle1,
                                 const CoordinateAxis& axis1,
                                 Real angle2,
                                 const CoordinateAxis& axis2) -> ::testing::AssertionResult {
    const Rotation rotationSpecified(bodyOrSpace, angle1, axis1, angle2, axis2);

    const Rotation AB(angle1, axis1);
    const Rotation BC(angle2, axis2);
    const Rotation testRotation = (bodyOrSpace == BodyRotationSequence) ? (AB * BC) : (BC * AB);

    if (!rotationSpecified.areAllRotationElementsSameToMachinePrecision(testRotation)) {
        return ::testing::AssertionFailure() << "Two-axis Rotation does not match manual composition";
    }

    const Vec2 angles = rotationSpecified.convertTwoAxesRotationToTwoAngles(bodyOrSpace, axis1, axis2);
    const Real theta1 = angles[0];
    const Real theta2 = angles[1];

    Rotation recovered;
    recovered.setRotationFromTwoAnglesTwoAxes(bodyOrSpace, theta1, axis1, theta2, axis2);
    if (!rotationSpecified.areAllRotationElementsSameToMachinePrecision(recovered)) {
        return ::testing::AssertionFailure() << "Two-axis rotation rebuilt from recovered angles does not "
                                             << "match original";
    }

    const bool angleOk = axis1.isSameAxis(axis2) ? inverseAngle1IsValid((angle1 + angle2), (theta1 + theta2))
                                                 : inverseAngle2IsValid(angle1, theta1, angle2, theta2);

    if (!angleOk) {
        return ::testing::AssertionFailure()
               << "Two-axis angle round-trip failed: "
               << "angle1=" << angle1 << " theta1=" << theta1 << " angle2=" << angle2 << " theta2=" << theta2;
    }
    return ::testing::AssertionSuccess();
}

/// Verifies a three-axis rotation:
///  (a) The compound constructor matches manual multiplication.
///  (b) Rebuilding from recovered angles reproduces the matrix.
///  (c) Angle round-trips are within tolerance, with special handling near
///      singularities and for degenerate axis combinations.
static auto checkThreeAxesRotation(BodyOrSpaceType bodyOrSpace,
                                   Real angle1,
                                   const CoordinateAxis& axis1,
                                   Real angle2,
                                   const CoordinateAxis& axis2,
                                   Real angle3,
                                   const CoordinateAxis& axis3) -> ::testing::AssertionResult {
    const Rotation rotationSpecified(bodyOrSpace, angle1, axis1, angle2, axis2, angle3, axis3);

    const Rotation AB(angle1, axis1);
    const Rotation BC(angle2, axis2);
    const Rotation CD(angle3, axis3);
    const Rotation testRotation = (bodyOrSpace == BodyRotationSequence) ? (AB * BC * CD) : (CD * BC * AB);

    if (!rotationSpecified.areAllRotationElementsSameToMachinePrecision(testRotation)) {
        return ::testing::AssertionFailure() << "Three-axis Rotation does not match manual composition";
    }

    const Vec3 angles =
        rotationSpecified.convertThreeAxesRotationToThreeAngles(bodyOrSpace, axis1, axis2, axis3);
    const Real theta1 = angles[0];
    const Real theta2 = angles[1];
    const Real theta3 = angles[2];

    Rotation recovered;
    recovered.setRotationFromThreeAnglesThreeAxes(bodyOrSpace, theta1, axis1, theta2, axis2, theta3, axis3);
    if (!rotationSpecified.areAllRotationElementsSameToMachinePrecision(recovered)) {
        return ::testing::AssertionFailure() << "Three-axis rotation rebuilt from recovered angles does not "
                                             << "match original";
    }

    bool angleOk = false;
    if (axis1.areAllSameAxes(axis2, axis3)) {
        angleOk = inverseAngle1IsValid((angle1 + angle2 + angle3), (theta1 + theta2 + theta3));
    } else if (axis1.isSameAxis(axis2)) {
        angleOk = inverseAngle2IsValid((angle1 + angle2), (theta1 + theta2), angle3, theta3);
    } else if (axis2.isSameAxis(axis3)) {
        angleOk = inverseAngle2IsValid(angle1, theta1, (angle2 + angle3), (theta2 + theta3));
    } else if (axis1.isSameAxis(axis3)) {
        angleOk = inverseAngleTwoAxesIsValid(angle1, theta1, angle2, theta2, angle3, theta3);
    } else {
        angleOk = inverseAngleThreeAxesIsValid(angle1, theta1, angle2, theta2, angle3, theta3);
    }

    if (!angleOk) {
        return ::testing::AssertionFailure() << "Three-axis angle round-trip failed";
    }
    return ::testing::AssertionSuccess();
}

// ============================================================
// Tests: Construction
// ============================================================

/// The default Rotation constructor must produce the 3×3 identity matrix.
TEST(SimTKCommon_Rotation_DefaultConstructor, ProducesIdentityMatrix) {
    const Rotation defaultRotation;
    Rotation identity;
    identity.setRotationToIdentityMatrix();

    EXPECT_TRUE(defaultRotation.areAllRotationElementsSameToMachinePrecision(identity));
}

/// Copying a Rotation must produce a matrix with identical elements.
TEST(SimTKCommon_Rotation_CopyConstructor, PreservesAllElements) {
    Rotation source;
    source.setRotationFromAngleAboutNonUnitVector(1.0, Vec3(0.2, 0.4, 0.6));
    const Rotation copy(source);

    EXPECT_TRUE(copy.areAllRotationElementsSameToMachinePrecision(source));
}

/// The assignment operator must produce a matrix with identical elements.
TEST(SimTKCommon_Rotation_AssignmentOperator, PreservesAllElements) {
    Rotation source;
    source.setRotationFromAngleAboutNonUnitVector(1.0, Vec3(0.2, 0.4, 0.6));
    const Rotation assigned = source;

    EXPECT_TRUE(assigned.areAllRotationElementsSameToMachinePrecision(source));
}

/// Rotation(angle, CoordinateAxis) must match Rotation(angle, unitVector)
/// for the same axis direction.
TEST(SimTKCommon_Rotation_CoordinateAxisConstructor, MatchesUnitVectorConstruction) {
    const Real angle = 0.1;
    const CoordinateAxis axis = XAxis;

    Rotation fromNonUnit;
    fromNonUnit.setRotationFromAngleAboutNonUnitVector(angle, Vec3(1.0, 0.0, 0.0));
    const Rotation fromCoordAxis(angle, axis);

    EXPECT_TRUE(fromCoordAxis.areAllRotationElementsSameToMachinePrecision(fromNonUnit));
}

/// The angle recovered by convertOneAxisRotationToOneAngle must match the
/// original construction angle within 10×SignificantReal.
TEST(SimTKCommon_Rotation_CoordinateAxisConstructor, AngleRoundTripIsWithinMachinePrecision) {
    const Real angle = 0.1;
    const CoordinateAxis axis = XAxis;
    const Rotation rotation(angle, axis);

    const Real recovered = rotation.convertOneAxisRotationToOneAngle(axis);

    EXPECT_LT(std::fabs(angle - recovered), (10 * SignificantReal));
}

// ============================================================
// Tests: Single-axis rotations
// ============================================================

/// Rotation about X, Y, Z with representative positive and negative angles
/// must both match the equivalent unit-vector construction and recover the
/// original angle.
TEST(SimTKCommon_Rotation_OneAxisRotation, MatchesAndRecovery_SelectedAngles) {
    EXPECT_TRUE(checkOneAxisRotation(0.2, XAxis));
    EXPECT_TRUE(checkOneAxisRotation(-0.2, XAxis));
    EXPECT_TRUE(checkOneAxisRotation(2.1, YAxis));
    EXPECT_TRUE(checkOneAxisRotation(-2.1, YAxis));
    EXPECT_TRUE(checkOneAxisRotation(3.1, ZAxis));
    EXPECT_TRUE(checkOneAxisRotation(-3.1, ZAxis));
}

/// Exhaustive sweep of all three axes over [-385°, +385°] in 0.5° steps.
TEST(SimTKCommon_Rotation_OneAxisRotation, ExhaustiveSweepAllAxesAndAngles) {
    const Real start = convertDegreesToRadians(-385.0);
    const Real stop = convertDegreesToRadians(385.0);
    const Real step = convertDegreesToRadians(0.5);

    for (int i = 0; i <= 2; ++i) {
        const CoordinateAxis axis = CoordinateAxis::getCoordinateAxis(i);
        for (Real angle = start; angle < stop; angle += step) {
            EXPECT_TRUE(checkOneAxisRotation(angle, axis)) << "axis=" << i << " angle=" << angle;
        }
    }
}

// ============================================================
// Tests: Two-axis rotations
// ============================================================

/// Selected two-axis body/space sequences covering all 9 axis-pair
/// combinations must match manual composition and recover both angles.
TEST(SimTKCommon_Rotation_TwoAxesRotation, BodyAndSpaceSelectedCases) {
    // XX
    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, 0.2, XAxis, 0.3, XAxis));
    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, 0.2, XAxis, 0.3, XAxis));
    // XY
    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, 1.2, XAxis, -1.3, YAxis));
    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, 1.2, XAxis, -1.3, YAxis));
    // XZ
    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, -3.1, XAxis, 1.2, ZAxis));
    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, -3.1, XAxis, 1.2, ZAxis));
    // YX
    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, 1.2, YAxis, 0.3, XAxis));
    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, 1.2, YAxis, 0.3, XAxis));
    // YY
    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, 2.2, YAxis, -1.3, YAxis));
    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, 2.2, YAxis, -1.3, YAxis));
    // YZ
    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, -3.1, YAxis, 1.2, ZAxis));
    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, -3.1, YAxis, 1.2, ZAxis));
    // ZX
    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, 1.2, ZAxis, 0.3, XAxis));
    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, 1.2, ZAxis, 0.3, XAxis));
    // ZY
    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, 2.2, ZAxis, -1.3, YAxis));
    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, 2.2, ZAxis, -1.3, YAxis));
    // ZZ
    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, -3.1, ZAxis, 1.2, ZAxis));
    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, -3.1, ZAxis, 1.2, ZAxis));
}

/// Exhaustive sweep of all 9 axis pairs over ±200° in 10° steps for both
/// body and space sequences.
TEST(SimTKCommon_Rotation_TwoAxesRotation, ExhaustiveSweepAllAxisPairsAndAngles) {
    const Real start = convertDegreesToRadians(-200.0);
    const Real stop = convertDegreesToRadians(200.0);
    const Real step = convertDegreesToRadians(10.0);

    for (int i = 0; i <= 2; ++i) {
        const CoordinateAxis axisi = CoordinateAxis::getCoordinateAxis(i);
        for (int j = 0; j <= 2; ++j) {
            const CoordinateAxis axisj = CoordinateAxis::getCoordinateAxis(j);
            for (Real ai = start; ai < stop; ai += step) {
                for (Real aj = start; aj < stop; aj += step) {
                    EXPECT_TRUE(checkTwoAxesRotation(BodyRotationSequence, ai, axisi, aj, axisj))
                        << "Body: i=" << i << " j=" << j << " ai=" << ai << " aj=" << aj;
                    EXPECT_TRUE(checkTwoAxesRotation(SpaceRotationSequence, ai, axisi, aj, axisj))
                        << "Space: i=" << i << " j=" << j << " ai=" << ai << " aj=" << aj;
                }
            }
        }
    }
}

// ============================================================
// Tests: Three-axis rotations
// ============================================================

/// Selected three-axis body/space sequences spanning all axis-triple
/// combinations starting with X, Y, or Z must match manual composition and
/// recover all three angles.
TEST(SimTKCommon_Rotation_ThreeAxesRotation, BodyAndSpaceSelectedCases) {
    // --- X first ---
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 0.2, XAxis, 0.3, XAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 0.2, XAxis, 0.3, XAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 1.2, XAxis, -1.3, XAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 1.2, XAxis, -1.3, XAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, -3.1, XAxis, 1.2, XAxis, 1.3, ZAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, -3.1, XAxis, 1.2, XAxis, 1.3, ZAxis));

    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 0.2, XAxis, 0.3, YAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 0.2, XAxis, 0.3, YAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 1.2, XAxis, -1.3, YAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 1.2, XAxis, -1.3, YAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, -3.1, XAxis, 1.2, YAxis, 1.3, ZAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, -3.1, XAxis, 1.2, YAxis, 1.3, ZAxis));

    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 0.2, XAxis, 0.3, ZAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 0.2, XAxis, 0.3, ZAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 1.2, XAxis, -1.3, ZAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 1.2, XAxis, -1.3, ZAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, -3.1, XAxis, 1.2, ZAxis, 1.3, ZAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, -3.1, XAxis, 1.2, ZAxis, 1.3, ZAxis));

    // --- Y first ---
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 0.2, YAxis, 0.3, XAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 0.2, YAxis, 0.3, XAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 1.2, YAxis, -1.3, XAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 1.2, YAxis, -1.3, XAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, -3.1, YAxis, 1.2, XAxis, 1.3, ZAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, -3.1, YAxis, 1.2, XAxis, 1.3, ZAxis));

    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 0.2, YAxis, 0.3, YAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 0.2, YAxis, 0.3, YAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 1.2, YAxis, -1.3, YAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 1.2, YAxis, -1.3, YAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, -3.1, YAxis, 1.2, YAxis, 1.3, ZAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, -3.1, YAxis, 1.2, YAxis, 1.3, ZAxis));

    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 0.2, YAxis, 0.3, ZAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 0.2, YAxis, 0.3, ZAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 1.2, YAxis, -1.3, ZAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 1.2, YAxis, -1.3, ZAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, -3.1, YAxis, 1.2, ZAxis, 1.3, ZAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, -3.1, YAxis, 1.2, ZAxis, 1.3, ZAxis));

    // --- Z first ---
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 0.2, ZAxis, 0.3, XAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 0.2, ZAxis, 0.3, XAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 1.2, ZAxis, -1.3, XAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 1.2, ZAxis, -1.3, XAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, -3.1, ZAxis, 1.2, XAxis, 1.3, ZAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, -3.1, ZAxis, 1.2, XAxis, 1.3, ZAxis));

    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 0.2, ZAxis, 0.3, YAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 0.2, ZAxis, 0.3, YAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 1.2, ZAxis, -1.3, YAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 1.2, ZAxis, -1.3, YAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, -3.1, ZAxis, 1.2, YAxis, 1.3, ZAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, -3.1, ZAxis, 1.2, YAxis, 1.3, ZAxis));

    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 0.2, ZAxis, 0.3, ZAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 0.2, ZAxis, 0.3, ZAxis, 0.4, XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, 1.2, ZAxis, -1.3, ZAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, 1.2, ZAxis, -1.3, ZAxis, -1.4, YAxis));
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence, -3.1, ZAxis, 1.2, ZAxis, 1.3, ZAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence, -3.1, ZAxis, 1.2, ZAxis, 1.3, ZAxis));
}

/// An XYX sequence with angles near ±π was a known failure before a bug fix.
/// Both body and space variants must now correctly match manual composition.
TEST(SimTKCommon_Rotation_ThreeAxesRotation, PreviouslyFailingXYXNearPi) {
    EXPECT_TRUE(checkThreeAxesRotation(BodyRotationSequence,
                                       -3.2288591161895095,
                                       XAxis,
                                       -SimTK_PI,
                                       YAxis,
                                       -SimTK_PI,
                                       XAxis));
    EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence,
                                       -3.2288591161895095,
                                       XAxis,
                                       -SimTK_PI,
                                       YAxis,
                                       -SimTK_PI,
                                       XAxis));
}

/// Exhaustive sweep of all 27 axis triples over ±200° in 40° steps for both
/// body and space sequences.
TEST(SimTKCommon_Rotation_ThreeAxesRotation, ExhaustiveSweepAllAxisTriplesAndAngles) {
    const Real start = convertDegreesToRadians(-200.0);
    const Real stop = convertDegreesToRadians(200.0);
    const Real step = convertDegreesToRadians(40.0);

    for (int i = 0; i <= 2; ++i) {
        const CoordinateAxis axisi = CoordinateAxis::getCoordinateAxis(i);
        for (int j = 0; j <= 2; ++j) {
            const CoordinateAxis axisj = CoordinateAxis::getCoordinateAxis(j);
            for (int k = 0; k <= 2; ++k) {
                const CoordinateAxis axisk = CoordinateAxis::getCoordinateAxis(k);
                for (Real ai = start; ai < stop; ai += step) {
                    for (Real aj = start; aj < stop; aj += step) {
                        for (Real ak = start; ak < stop; ak += step) {
                            EXPECT_TRUE(
                                checkThreeAxesRotation(BodyRotationSequence, ai, axisi, aj, axisj, ak, axisk))
                                << "Body: i=" << i << " j=" << j << " k=" << k << " ai=" << ai << " aj=" << aj
                                << " ak=" << ak;
                            EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence,
                                                               ai,
                                                               axisi,
                                                               aj,
                                                               axisj,
                                                               ak,
                                                               axisk))
                                << "Space: i=" << i << " j=" << j << " k=" << k << " ai=" << ai
                                << " aj=" << aj << " ak=" << ak;
                        }
                    }
                }
            }
        }
    }
}

/// For two-axis (same first-and-last axis) sequences, the middle angle sweeps
/// within ±1e-11° of the singularity at 0°.  All axis triples are covered
/// with outer angles varying over ±200° in 20° steps.
TEST(SimTKCommon_Rotation_ThreeAxesRotation, ExhaustiveSweepNearSingularityTwoAxesSequence) {
    const Real outerStart = convertDegreesToRadians(-200.0);
    const Real outerStop = convertDegreesToRadians(200.0);
    const Real outerStep = convertDegreesToRadians(20.0);

    const Real singStart = convertDegreesToRadians(0.0 - 1.0e-11);
    const Real singStop = convertDegreesToRadians(0.0 + 1.0e-11);
    const Real singStep = convertDegreesToRadians(1.0e-12);

    for (int i = 0; i <= 2; ++i) {
        const CoordinateAxis axisi = CoordinateAxis::getCoordinateAxis(i);
        for (int j = 0; j <= 2; ++j) {
            const CoordinateAxis axisj = CoordinateAxis::getCoordinateAxis(j);
            for (int k = 0; k <= 2; ++k) {
                const CoordinateAxis axisk = CoordinateAxis::getCoordinateAxis(k);
                for (Real ai = outerStart; ai < outerStop; ai += outerStep) {
                    for (Real aj = singStart; aj < singStop; aj += singStep) {
                        for (Real ak = outerStart; ak < outerStop; ak += outerStep) {
                            EXPECT_TRUE(
                                checkThreeAxesRotation(BodyRotationSequence, ai, axisi, aj, axisj, ak, axisk))
                                << "Body near sing0: i=" << i << " j=" << j << " k=" << k;
                            EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence,
                                                               ai,
                                                               axisi,
                                                               aj,
                                                               axisj,
                                                               ak,
                                                               axisk))
                                << "Space near sing0: i=" << i << " j=" << j << " k=" << k;
                        }
                    }
                }
            }
        }
    }
}

/// For three distinct-axis sequences, the middle angle sweeps within ±1e-11°
/// of the singularity at 90°.  All axis triples are covered with outer
/// angles varying over ±200° in 20° steps.
TEST(SimTKCommon_Rotation_ThreeAxesRotation, ExhaustiveSweepNearSingularityThreeAxesSequence) {
    const Real outerStart = convertDegreesToRadians(-200.0);
    const Real outerStop = convertDegreesToRadians(200.0);
    const Real outerStep = convertDegreesToRadians(20.0);

    const Real singStart = convertDegreesToRadians(90.0 - 1.0e-11);
    const Real singStop = convertDegreesToRadians(90.0 + 1.0e-11);
    const Real singStep = convertDegreesToRadians(1.0e-12);

    for (int i = 0; i <= 2; ++i) {
        const CoordinateAxis axisi = CoordinateAxis::getCoordinateAxis(i);
        for (int j = 0; j <= 2; ++j) {
            const CoordinateAxis axisj = CoordinateAxis::getCoordinateAxis(j);
            for (int k = 0; k <= 2; ++k) {
                const CoordinateAxis axisk = CoordinateAxis::getCoordinateAxis(k);
                for (Real ai = outerStart; ai < outerStop; ai += outerStep) {
                    for (Real aj = singStart; aj < singStop; aj += singStep) {
                        for (Real ak = outerStart; ak < outerStop; ak += outerStep) {
                            EXPECT_TRUE(
                                checkThreeAxesRotation(BodyRotationSequence, ai, axisi, aj, axisj, ak, axisk))
                                << "Body near sing90: i=" << i << " j=" << j << " k=" << k;
                            EXPECT_TRUE(checkThreeAxesRotation(SpaceRotationSequence,
                                                               ai,
                                                               axisi,
                                                               aj,
                                                               axisj,
                                                               ak,
                                                               axisk))
                                << "Space near sing90: i=" << i << " j=" << j << " k=" << k;
                        }
                    }
                }
            }
        }
    }
}

// ============================================================
// Tests: Quaternion
// ============================================================

/// Selected quaternion components must survive the round-trip
/// quaternion → Rotation → quaternion → Rotation with identical elements.
TEST(SimTKCommon_Rotation_Quaternion, RoundTripPreservesRotation) {
    auto checkQuaternion = [](Real e0, Real e1, Real e2, Real e3) -> ::testing::AssertionResult {
        const Quaternion q(e0, e1, e2, e3);
        const Rotation fromQ(q);
        const Quaternion qBack = fromQ.convertRotationToQuaternion();

        Rotation fromQBack;
        fromQBack.setRotationFromQuaternion(qBack);

        if (!fromQ.areAllRotationElementsSameToMachinePrecision(fromQBack)) {
            return ::testing::AssertionFailure() << "Quaternion round-trip failed for (" << e0 << ", " << e1
                                                 << ", " << e2 << ", " << e3 << ")";
        }
        return ::testing::AssertionSuccess();
    };

    EXPECT_TRUE(checkQuaternion(0.5, 0.1, 0.2, 0.3));
    EXPECT_TRUE(checkQuaternion(-0.5, 0.1, 0.2, -0.3));
}

/// Exhaustive grid of quaternion components in [-1, 1]^4 in 0.2 steps.
TEST(SimTKCommon_Rotation_Quaternion, ExhaustiveRoundTrip) {
    for (Real e0 = -1; e0 <= 1; e0 += 0.2) {
        for (Real e1 = -1; e1 <= 1; e1 += 0.2) {
            for (Real e2 = -1; e2 <= 1; e2 += 0.2) {
                for (Real e3 = -1; e3 <= 1; e3 += 0.2) {
                    const Quaternion q(e0, e1, e2, e3);
                    const Rotation fromQ(q);
                    const Quaternion qBack = fromQ.convertRotationToQuaternion();
                    Rotation fromQBack;
                    fromQBack.setRotationFromQuaternion(qBack);

                    EXPECT_TRUE(fromQ.areAllRotationElementsSameToMachinePrecision(fromQBack))
                        << "Failed for e0=" << e0 << " e1=" << e1 << " e2=" << e2 << " e3=" << e3;
                }
            }
        }
    }
}

/// An explicitly un-normalised Quaternion must not have unit norm.
/// Both normalize() and normalizeThis() must produce a unit quaternion.
TEST(SimTKCommon_Rotation_Quaternion, NormalizationProducesUnitNorm) {
    Quaternion unnorm(Vec4(1, 2, 3, 4), true); // skip normalisation intentionally

    EXPECT_FALSE(SimTK::Test::numericallyEqual(unnorm.norm(), Real(1), 1))
        << "Un-normalised quaternion should NOT have unit norm";

    const Quaternion fixedUp = unnorm.normalize();
    EXPECT_SIMTK_EQ(fixedUp.norm(), Real(1));

    unnorm.normalizeThis();
    EXPECT_SIMTK_EQ(unnorm.norm(), Real(1));
}

// ============================================================
// Tests: Nearby orthogonal rotation from Mat33
// ============================================================

/// Constructing a Rotation from the raw Mat33 of an existing Rotation must
/// reproduce the same (nearest orthogonal) rotation matrix.
TEST(SimTKCommon_Rotation_NearbyOrthogonal, ConstructionFromMat33MatchesOriginal) {
    Rotation source;
    source.setRotationFromAngleAboutNonUnitVector(1.0, Vec3(0.2, 0.4, 0.6));
    const Rotation nearby(source.asMat33());

    EXPECT_TRUE(nearby.areAllRotationElementsSameToMachinePrecision(source));
}

// ============================================================
// Tests: Rotation from two given axes
// ============================================================

/// Rotation(UnitVec3, CoordinateAxis, Vec3, CoordinateAxis) must produce a
/// proper rotation matrix (det = +1), and both column axes must point in the
/// positive direction of the given vectors.
TEST(SimTKCommon_Rotation_TwoGivenAxes, DeterminantIsOneAndAxesHaveCorrectOrientation) {
    const UnitVec3 vi(0.01, 0.02, 0.9);
    const Vec3 vj(-0.5, 0.5, 0.2);

    struct AxisPair {
        CoordinateAxis ai;
        CoordinateAxis aj;
    };

    const AxisPair pairs[] = {
        {XAxis, YAxis},
        {YAxis, XAxis},
        {YAxis, ZAxis},
        {ZAxis, YAxis},
        {ZAxis, XAxis},
        {XAxis, ZAxis},
    };

    for (const auto& [ai, aj] : pairs) {
        const Rotation R(UnitVec3(vi), ai, vj, aj);

        EXPECT_LE(std::fabs(det(R) - Real(1)), SignificantReal) << "det != 1 for ai=" << ai << " aj=" << aj;
        EXPECT_GT(dot(R(ai), vi), Real(0)) << "i-axis column not aligned for ai=" << ai << " aj=" << aj;
        EXPECT_GT(dot(R(aj), vj), Real(0)) << "j-axis column not aligned for ai=" << ai << " aj=" << aj;
    }
}

// ============================================================
// Tests: Body-fixed XYZ
// ============================================================

/// All three body-fixed XYZ setters must produce identical rotation matrices:
///  (a) generic three-axis setter with BodyRotationSequence,
///  (b) setRotationToBodyFixedXYZ(angles),
///  (c) setRotationToBodyFixedXYZ(cosines, sines).
TEST(SimTKCommon_Rotation_BodyFixedXYZ, AllSettersProduceSameResult) {
    const Real q0 = 0.123;
    const Real q1 = -0.234;
    const Real q2 = 0.787;

    Rotation R0;
    R0.setRotationFromThreeAnglesThreeAxes(BodyRotationSequence, q0, XAxis, q1, YAxis, q2, ZAxis);

    Rotation R1;
    R1.setRotationToBodyFixedXYZ(Vec3(q0, q1, q2));

    Rotation R2;
    R2.setRotationToBodyFixedXYZ(Vec3(std::cos(q0), std::cos(q1), std::cos(q2)),
                                 Vec3(std::sin(q0), std::sin(q1), std::sin(q2)));

    EXPECT_TRUE(R0.areAllRotationElementsSameToMachinePrecision(R1))
        << "Generic XYZ setter and setRotationToBodyFixedXYZ(angles) differ";
    EXPECT_TRUE(R0.areAllRotationElementsSameToMachinePrecision(R2))
        << "Generic XYZ setter and setRotationToBodyFixedXYZ(cos,sin) differ";
    EXPECT_TRUE(R1.areAllRotationElementsSameToMachinePrecision(R2))
        << "setRotationToBodyFixedXYZ(angles) and (cos,sin) differ";
}

// ============================================================
// Tests: Re-expression of a symmetric 3×3 matrix
// ============================================================

/// R_AB.reexpressSymMat33(S_BB) must equal R_AB * Mat33(S_BB) * ~R_AB
/// to within SignificantReal.
TEST(SimTKCommon_Rotation_ReexpressSymMat33, MatchesFullMatrixTripleProduct) {
    const Rotation R_AB(SimTK::Test::randRotation());
    const SymMat33 S_BB(SimTK::Test::randSymMat<3>());
    const Mat33 M_BB(S_BB);

    const SymMat33 S_AA = R_AB.reexpressSymMat33(S_BB);
    const Mat33 M_AA = (R_AB * M_BB) * ~R_AB;
    const Mat33 MS_AA(S_AA);

    EXPECT_LE((MS_AA - M_AA).norm(), SignificantReal);
}

/// Applying the InverseRotation to the re-expressed matrix must recover the
/// original symmetric matrix.
TEST(SimTKCommon_Rotation_ReexpressSymMat33, InverseRotationRecoverOriginalMatrix) {
    const Rotation R_AB(SimTK::Test::randRotation());
    const SymMat33 S_BB(SimTK::Test::randSymMat<3>());

    const SymMat33 S_AA = R_AB.reexpressSymMat33(S_BB);
    const SymMat33 isS_BB = (~R_AB).reexpressSymMat33(S_AA);

    EXPECT_LE((S_BB - isS_BB).norm(), SignificantReal);
}

/// Re-expression through the identity rotation must leave the matrix unchanged.
TEST(SimTKCommon_Rotation_ReexpressSymMat33, IdentityRotationLeavesMatrixUnchanged) {
    const SymMat33 S_BB(SimTK::Test::randSymMat<3>());
    const Rotation I;

    const SymMat33 result = I.reexpressSymMat33(S_BB);

    EXPECT_LE((result - S_BB).norm(), SignificantReal);
}

/// SymMat33 * SymMat33 must equal the product of the corresponding Mat33s.
TEST(SimTKCommon_Rotation_ReexpressSymMat33, SymMat33MultiplicationMatchesFullMat33) {
    const SymMat33 S1(SimTK::Test::randSymMat<3>());
    const SymMat33 S2(SimTK::Test::randSymMat<3>());
    const Mat33 M1(S1);
    const Mat33 M2(S2);

    const Mat33 fromSym(S1 * S2);
    const Mat33 fromFull(M1 * M2);

    EXPECT_LE((fromSym - fromFull).norm(), SignificantReal);
}

/// A random SymMat<3,Complex> must satisfy Hermitian symmetry:
/// elt(i,j) == conj(elt(j,i)) for all off-diagonal pairs.
TEST(SimTKCommon_Rotation_ReexpressSymMat33, ComplexSymMatHasHermitianSymmetry) {
    const SymMat<3, Complex> SC(SimTK::Test::randComplex(),
                                SimTK::Test::randComplex(),
                                SimTK::Test::randComplex(),
                                SimTK::Test::randComplex(),
                                SimTK::Test::randComplex(),
                                SimTK::Test::randComplex());

    EXPECT_TRUE(AssertSimTKEqual("SC.elt(1,0)", "conj(SC.elt(0,1))", SC.elt(1, 0), conj(SC.elt(0, 1))));
    EXPECT_TRUE(AssertSimTKEqual("SC.elt(2,0)", "conj(SC.elt(0,2))", SC.elt(2, 0), conj(SC.elt(0, 2))));
    EXPECT_TRUE(AssertSimTKEqual("SC.elt(1,2)", "conj(SC.elt(2,1))", SC.elt(1, 2), conj(SC.elt(2, 1))));
}

/// SymMat<3,Complex> * SymMat<3,Complex> must match the corresponding
/// Mat<3,3,Complex> product element-by-element.
TEST(SimTKCommon_Rotation_ReexpressSymMat33, ComplexSymMatMultiplicationMatchesFullComplexMat33) {
    const SymMat<3, Complex> SC1(SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex());
    const SymMat<3, Complex> SC2(SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex(),
                                 SimTK::Test::randComplex());

    const Mat<3, 3, Complex> MC1(SC1);
    const Mat<3, 3, Complex> MC2(SC2);
    const Mat<3, 3, Complex> fromSym(SC1 * SC2);
    const Mat<3, 3, Complex> fromFull(MC1 * MC2);

    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_TRUE(AssertSimTKEqual("fromSym(i,j)", "fromFull(i,j)", fromSym(i, j), fromFull(i, j)))
                << "Mismatch at (" << i << ", " << j << ")";
        }
    }
}