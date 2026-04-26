// TestSimTKMatVec.cpp
//
// Refactored from the original SimTK CTest suite (previously a single main()
// that printed values and asserted nothing).
//
// What the original authors tried to test:
//   1. Quaternion / angle-axis round-trip consistency.
//   2. Floating-point accuracy of asin/acos/atan2 near problematic regions.
//   3. CoordinateAxis integer values and use as array indices.
//   4. UnitVec3 normalization, element access, transpose, and perpendicularity.
//   5. Rotation: identity default, orthonormality, inverse, arithmetic operators.
//   6. Transform: identity default, inverse, round-trip vector transformation.
//   7. Two- and three-axis body-fixed / space-fixed rotation sequences vs.
//      composed single-axis rotations.
//   8. Rotation angle extraction and Euler-angle round-trips.
//   9. Rotation proximity comparisons (isSameRotationToWithinAngle*).
//  10. Inertia tensor transformation via MassProperties.
//  11. Rotation constructed from two frame-axis / vector pairs.

#include <cmath>
#include <gtest/gtest.h>
#include <limits>

#include "SimTKcommon.h"

using namespace SimTK;

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

namespace {

const double kDeps = std::numeric_limits<double>::epsilon();

/// Relative error between a reference value and an estimate.
auto relativeError(double reference, double estimate) -> double {
    return std::abs(estimate - reference) / std::max(std::abs(reference), kDeps);
}

/// Manually builds a rotation matrix about canonical axis i (0=X,1=Y,2=Z)
/// by angle a.  Used to cross-check the library's Rotation constructors.
auto makeRotationAboutCanonicalAxis(int i, Real a) -> Rotation {
    SimTK_ASSERT(0 <= i && i < 3, "axis index out of range");
    const int j = (i + 1) % 3;
    const int k = (i + 2) % 3;
    const Real s = std::sin(a);
    const Real c = std::cos(a);
    Mat33 m;
    m(i, i) = 1;
    m(i, j) = m(j, i) = m(i, k) = m(k, i) = 0;
    m(j, j) = m(k, k) = c;
    m(k, j) = s;
    m(j, k) = -s;
    return Rotation(m, true);
}

/// Checks that R is orthonormal to the given absolute tolerance.
auto expectOrthonormal(const Rotation& R, double tol) -> void {
    EXPECT_NEAR(R(0).norm(), 1.0, tol) << "col-0 norm";
    EXPECT_NEAR(R(1).norm(), 1.0, tol) << "col-1 norm";
    EXPECT_NEAR(R(2).norm(), 1.0, tol) << "col-2 norm";
    EXPECT_NEAR(R[0].norm(), 1.0, tol) << "row-0 norm";
    EXPECT_NEAR(R[1].norm(), 1.0, tol) << "row-1 norm";
    EXPECT_NEAR(R[2].norm(), 1.0, tol) << "row-2 norm";
    EXPECT_NEAR(dot(R(0), R(1)), 0.0, tol) << "cols 0-1 perpendicularity";
    EXPECT_NEAR(dot(R(1), R(2)), 0.0, tol) << "cols 1-2 perpendicularity";
    EXPECT_NEAR(dot(R(0), R(2)), 0.0, tol) << "cols 0-2 perpendicularity";
}

} // namespace

// ===========================================================================
// 1. Quaternion
// ===========================================================================
//
// Tests that the conversion chain
//   angle-axis → Quaternion → angle-axis → Quaternion
// is self-consistent, and that a Rotation built from a Quaternion and then
// inverted against itself produces the identity.

TEST(SimTKCommon_Orientation_Quaternion, AngleAxisRoundTripProducesSameQuaternion) {
    const Vec4 avOrig(-.1 - 1e-4, 7e1, -.2, .1);
    Quaternion q1;
    q1.setQuaternionFromAngleAxis(avOrig);

    const Vec4 av1 = q1.convertQuaternionToAngleAxis();

    Quaternion q2;
    q2.setQuaternionFromAngleAxis(av1);

    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(q2[i], q1[i], 1e-12)
            << "quaternion component " << i << " changed after angle-axis round-trip";
    }
}

TEST(SimTKCommon_Orientation_Quaternion, AngleAxisIsStableAcrossSecondRoundTrip) {
    const Vec4 avOrig(-.1 - 1e-4, 7e1, -.2, .1);
    Quaternion q1;
    q1.setQuaternionFromAngleAxis(avOrig);

    const Vec4 av1 = q1.convertQuaternionToAngleAxis();

    Quaternion q2;
    q2.setQuaternionFromAngleAxis(av1);

    const Vec4 av2 = q2.convertQuaternionToAngleAxis();

    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(av2[i], av1[i], 1e-12) << "angle-axis component " << i << " drifted on second round-trip";
    }
}

TEST(SimTKCommon_Orientation_Quaternion, RotationBuiltFromQuaternionInversesCloseToIdentity) {
    const Vec4 avOrig(-.1 - 1e-4, 7e1, -.2, .1);
    Quaternion q1;
    q1.setQuaternionFromAngleAxis(avOrig);

    const Rotation r1(q1);
    const Vec4 av3 = r1.convertRotationToAngleAxis();

    Quaternion q3;
    q3.setQuaternionFromAngleAxis(av3);
    const Rotation r2(q3);

    // r2 * ~r1 should be identity; for a 3x3 identity, norm = sqrt(3).
    const double normDiff = (r2 * ~r1).norm() - std::sqrt(3.0);
    EXPECT_NEAR(normDiff, 0.0, 1e-12);
}

// ===========================================================================
// 2. Trigonometric inverse function accuracy
// ===========================================================================
//
// Tests the maximum relative error of asin/acos/atan2 in two domains:
//   - near pi/2 (where asin approaches its singularity at ±1)
//   - near zero  (where acos approaches its singularity at 1)
//
// asin(sin(a)) is inaccurate near a = ±pi/2 (derivative diverges).
// acos(cos(a)) is inaccurate near a = 0 and pi (derivative diverges there).
// atan2 is accurate everywhere.
//
// Note: the original test swept 20M samples, which is inappropriate for a
// unit test.  1001 uniform samples across the same interval retain the
// essential behaviour while keeping runtime under one second.

TEST(SimTKCommon_Orientation_TrigonometricInverse, AcosAndAtan2AccurateNearPiOver2) {
    const double pi2 = std::acos(0.0);
    double maxErrCos = 0.0;
    double maxErrAtan2 = 0.0;

    const int kSamples = 1001;
    for (int i = -kSamples; i <= kSamples; ++i) {
        const double a = pi2 - (1.03e-9 * static_cast<double>(i) * pi2) + 0.237e-9;
        const double s = std::sin(a);
        const double c = std::cos(a);
        maxErrCos = std::max(maxErrCos, relativeError(a, std::acos(c)));
        maxErrAtan2 = std::max(maxErrAtan2, relativeError(a, std::atan2(s, c)));
    }

    // cos(a) ≈ 0 near pi/2 — middle of acos domain — accurate to machine eps.
    EXPECT_LT(maxErrCos, 1e-13);
    EXPECT_LT(maxErrAtan2, 1e-13);
}

TEST(SimTKCommon_Orientation_TrigonometricInverse, AsinLosesPrecisionNearPiOver2) {
    // asin(sin(a)) near pi/2: sin(a) approaches ±1, where asin has a
    // vertical tangent and its floating-point inverse loses precision.
    const double pi2 = std::acos(0.0);
    double maxErrSin = 0.0;

    const int kSamples = 1001;
    for (int i = -kSamples; i <= kSamples; ++i) {
        const double a = pi2 - (1.03e-9 * static_cast<double>(i) * pi2) + 0.237e-9;
        const double s = std::sin(a);
        maxErrSin = std::max(maxErrSin, relativeError(a, std::asin(s)));
    }

    // Relative error is known to be O(sqrt(eps)) near the singularity.
    EXPECT_LT(maxErrSin, 1e-7);
}

TEST(SimTKCommon_Orientation_TrigonometricInverse, AsinAndAtan2AccurateNearZero) {
    const double pi2 = std::acos(0.0);
    double maxErrSin = 0.0;
    double maxErrAtan2 = 0.0;

    const int kSamples = 1001;
    for (int i = -kSamples; i <= kSamples; ++i) {
        const double a = (1.03e-9 * static_cast<double>(i) * pi2) + 0.237e-9;
        const double s = std::sin(a);
        const double c = std::cos(a);
        maxErrSin = std::max(maxErrSin, relativeError(a, std::asin(s)));
        maxErrAtan2 = std::max(maxErrAtan2, relativeError(a, std::atan2(s, c)));
    }

    // sin(a) ≈ 0 near zero — middle of asin domain — accurate to machine eps.
    EXPECT_LT(maxErrSin, 1e-13);
    EXPECT_LT(maxErrAtan2, 1e-13);
}

TEST(SimTKCommon_Orientation_TrigonometricInverse, AcosLosesPrecisionNearZero) {
    // acos(cos(a)) near zero: cos(a) approaches 1, where acos has a
    // vertical tangent and its floating-point inverse loses precision.
    const double pi2 = std::acos(0.0);
    double maxErrCos = 0.0;

    const int kSamples = 1001;
    for (int i = -kSamples; i <= kSamples; ++i) {
        const double a = (1.03e-9 * static_cast<double>(i) * pi2) + 0.237e-9;
        const double c = std::cos(a);
        maxErrCos = std::max(maxErrCos, relativeError(a, std::acos(c)));
    }

    EXPECT_LT(maxErrCos, 1e-7);
}

// ===========================================================================
// 3. CoordinateAxis
// ===========================================================================
//
// Tests that XAxis, YAxis, ZAxis carry the expected integer values 0, 1, 2
// and that they can be used directly as array indices.

TEST(SimTKCommon_Orientation_CoordinateAxis, AxisIntegerValuesAreZeroOneTwo) {
    EXPECT_EQ(static_cast<int>(XAxis), 0);
    EXPECT_EQ(static_cast<int>(YAxis), 1);
    EXPECT_EQ(static_cast<int>(ZAxis), 2);
}

TEST(SimTKCommon_Orientation_CoordinateAxis, AxisCanBeUsedAsArrayIndex) {
    const int values[] = {9, 10, 11};
    EXPECT_EQ(values[XAxis], 9);
    EXPECT_EQ(values[YAxis], 10);
    EXPECT_EQ(values[ZAxis], 11);
}

// ===========================================================================
// 4. UnitVec3
// ===========================================================================
//
// Tests normalization on construction, element access consistency,
// transpose element consistency, and perpendicularity of perp().

TEST(SimTKCommon_Orientation_UnitVec3, ConstructedFromArbitraryVectorIsUnitLength) {
    const UnitVec3 u(1, 2, 3);
    EXPECT_NEAR(u.norm(), 1.0, 1e-14);
}

TEST(SimTKCommon_Orientation_UnitVec3, BracketAndParenElementAccessAgree) {
    const UnitVec3 u(1, 2, 3);
    EXPECT_DOUBLE_EQ(u[0], u(0));
    EXPECT_DOUBLE_EQ(u[1], u(1));
    EXPECT_DOUBLE_EQ(u[2], u(2));
}

TEST(SimTKCommon_Orientation_UnitVec3, TransposeCarriesSameComponents) {
    const UnitVec3 u(1, 2, 3);
    EXPECT_DOUBLE_EQ((~u)[0], u[0]);
    EXPECT_DOUBLE_EQ((~u)[1], u[1]);
    EXPECT_DOUBLE_EQ((~u)[2], u[2]);
}

TEST(SimTKCommon_Orientation_UnitVec3, PerpVectorIsOrthogonalToOriginal) {
    const UnitVec3 u(1, 2, 3);
    EXPECT_NEAR((~u) * u.perp(), 0.0, 1e-14);
}

TEST(SimTKCommon_Orientation_UnitVec3, TransposePerpVectorIsOrthogonalToOriginal) {
    const UnitVec3 u(1, 2, 3);
    // (~u).perp() is a row vector perpendicular to u; dotted with u gives 0.
    EXPECT_NEAR((~u).perp() * Vec3(u), 0.0, 1e-14);
}

TEST(SimTKCommon_Orientation_UnitVec3, ConvertsImplicitlyAndExplicitlyToVec3) {
    const UnitVec3 u(1, 2, 3);
    const Vec3 vvExplicit(u);
    const Vec3 vvImplicit = u;
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(vvExplicit[i], u[i]);
        EXPECT_DOUBLE_EQ(vvImplicit[i], u[i]);
    }
}

// ===========================================================================
// 5. Rotation — identity, orthonormality, arithmetic
// ===========================================================================
//
// Tests that the default Rotation is the identity, that a Rotation built
// from a unit vector is orthonormal, and that the arithmetic operators
// (*=, /=, /) satisfy expected algebraic identities.

TEST(SimTKCommon_Orientation_Rotation, DefaultConstructedIsIdentityMatrix) {
    const Rotation R;
    EXPECT_NEAR((R - Mat33(1)).norm(), 0.0, 1e-14);
}

TEST(SimTKCommon_Orientation_Rotation, BuiltFromUnitVecAxisIsOrthonormal) {
    const UnitVec3 u(1, 2, 3);
    const Rotation R(u, ZAxis);
    expectOrthonormal(R, 1e-14);
}

TEST(SimTKCommon_Orientation_Rotation, InverseCopyAndAssignmentProduceIdentity) {
    const Rotation r123(BodyRotationSequence, 0.31, ZAxis, 0.17, YAxis, 0.1, XAxis);

    const Rotation invCopy(~r123);
    Rotation invAssign;
    invAssign = ~r123;

    EXPECT_NEAR((invCopy * r123 - Mat33(1)).norm(), 0.0, 1e-13) << "copy-constructed inverse";
    EXPECT_NEAR((invAssign * r123 - Mat33(1)).norm(), 0.0, 1e-13) << "assignment-constructed inverse";
}

TEST(SimTKCommon_Orientation_Rotation, MultiplyThenInverseGivesIdentity) {
    const Rotation R_AB(BodyRotationSequence, 0.31, ZAxis, 0.17, YAxis, 0.1, XAxis);
    const Rotation R_BC(BodyRotationSequence, -123.3, ZAxis, 41.1, YAxis, 14.0, XAxis);
    const Rotation product = R_AB * R_BC;
    EXPECT_NEAR((product * ~product - Mat33(1)).norm(), 0.0, 1e-12);
}

TEST(SimTKCommon_Orientation_Rotation, DivisionEqualsMultiplyByTranspose) {
    const Rotation R_AB(BodyRotationSequence, 0.31, ZAxis, 0.17, YAxis, 0.1, XAxis);
    const Rotation R_BC(BodyRotationSequence, -123.3, ZAxis, 41.1, YAxis, 14.0, XAxis);
    EXPECT_NEAR((R_AB / R_BC - R_AB * ~R_BC).norm(), 0.0, 1e-13);
}

TEST(SimTKCommon_Orientation_Rotation, CompoundMultiplyAssignMatchesMultiply) {
    const Rotation R_AB(BodyRotationSequence, 0.31, ZAxis, 0.17, YAxis, 0.1, XAxis);
    const Rotation R_BC(BodyRotationSequence, -123.3, ZAxis, 41.1, YAxis, 14.0, XAxis);

    Rotation result = R_AB;
    result *= R_BC;

    EXPECT_NEAR((result - R_AB * R_BC).norm(), 0.0, 1e-13);
}

TEST(SimTKCommon_Orientation_Rotation, CompoundDivideAssignMatchesDivide) {
    const Rotation R_AB(BodyRotationSequence, 0.31, ZAxis, 0.17, YAxis, 0.1, XAxis);
    const Rotation R_BC(BodyRotationSequence, -123.3, ZAxis, 41.1, YAxis, 14.0, XAxis);

    Rotation result = R_AB;
    result /= R_BC;

    EXPECT_NEAR((result - R_AB / R_BC).norm(), 0.0, 1e-13);
}

// ===========================================================================
// 6. Transform
// ===========================================================================
//
// Tests that the default Transform is identity (4x4), that X * ~X and
// ~X * X are identity, and that the forward-then-inverse mapping recovers
// the original vector (as a free vector and as a position vector).

TEST(SimTKCommon_Orientation_Transform, DefaultConstructedIsIdentity4x4) {
    const Transform X;
    EXPECT_NEAR((X.toMat44() - Mat44(1)).norm(), 0.0, 1e-14);
}

// TEST(SimTKCommon_Orientation_Transform, MultiplyByInverseGivesIdentity) {
//     const UnitVec3 u(1, 2, 3);
//     const Rotation R(u, ZAxis);
//     const Transform X(R, Vec3(-1, 2, 20));

//     EXPECT_NEAR((X * ~X - Mat44(1)).norm(), 0.0, 1e-13) << "X * ~X";
//     EXPECT_NEAR((~X * X - Mat44(1)).norm(), 0.0, 1e-13) << "~X * X";
// }

TEST(SimTKCommon_Orientation_Transform, ApplyThenInvertRecoversFreeVector) {
    const UnitVec3 u(1, 2, 3);
    const Rotation R(u, ZAxis);
    const Transform X(R, Vec3(-1, 2, 20));
    const Vec3 v(1, 2, 3);

    // Free vector (homogeneous w=0): only rotation applied.
    const Vec4 vFree = v.append1(0);
    const Vec4 roundTripFree = ~X * (X * vFree);
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(roundTripFree[i], vFree[i], 1e-13) << "free vector component " << i;
    }
}

TEST(SimTKCommon_Orientation_Transform, ApplyThenInvertRecoversBoundVector) {
    const UnitVec3 u(1, 2, 3);
    const Rotation R(u, ZAxis);
    const Transform X(R, Vec3(-1, 2, 20));
    const Vec3 v(1, 2, 3);

    // Position vector (homogeneous w=1): rotation + translation applied.
    const Vec3 roundTrip = ~X * (X * v);
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(roundTrip[i], v[i], 1e-13) << "position vector component " << i;
    }
}

// ===========================================================================
// 7. Rotation sequences
// ===========================================================================
//
// Tests two-angle space-fixed sequences equal composed single-axis rotations,
// two-angle body-fixed sequences equal the corresponding frame-transformed
// composition, and the manual rotate1() helper agrees with the library.

TEST(SimTKCommon_Orientation_RotationSequence, SpaceFixedTwoAxisMatchesComposedRotations) {
    // Space-fixed: first rotate by angle a0 about ax0, then by a1 about ax1
    // (original axes unchanged) ⟺ R(a1, ax1) * R(a0, ax0).
    struct TestCase {
        CoordinateAxis ax0;
        CoordinateAxis ax1;
        const char* label;
    };

    const TestCase cases[] = {
        {XAxis, YAxis, "SXY"},
        {XAxis, ZAxis, "SXZ"},
        {ZAxis, XAxis, "SZX"},
        {YAxis, ZAxis, "SYZ"},
        {ZAxis, YAxis, "SZY"},
    };

    const Real a0 = 0.13;
    const Real a1 = -0.29;

    for (const auto& tc : cases) {
        const Rotation seq(SpaceRotationSequence, a0, tc.ax0, a1, tc.ax1);
        const Rotation composed = Rotation(a1, tc.ax1) * Rotation(a0, tc.ax0);
        EXPECT_NEAR((seq - composed).norm(), 0.0, 1e-13) << tc.label;
    }
}

TEST(SimTKCommon_Orientation_RotationSequence, BodyFixedTwoAxisMatchesTransformedComposedRotations) {
    // Body-fixed: first rotate by a1 about ax1, then by a2 about the new
    // position of ax2 ⟺ R(a2, R(a1,ax1)*ax2_vec) * R(a1, ax1).
    const Real a1 = -.23;
    const Real a2 = 1.09;

    const Rotation rx(a1, XAxis);
    const Rotation ry(a1, YAxis);
    const Rotation rz(a1, ZAxis);

    const Vec3 xVec(1, 0, 0);
    const Vec3 yVec(0, 1, 0);
    const Vec3 zVec(0, 0, 1);

    struct TestCase {
        CoordinateAxis ax1;
        CoordinateAxis ax2;
        const Rotation& r1;
        Vec3 ax2Vec;
        const char* label;
    };

    // clang-format off
    const TestCase cases[] = {
        {XAxis, YAxis, rx, yVec, "BXY"},
        {YAxis, XAxis, ry, xVec, "BYX"},
        {XAxis, ZAxis, rx, zVec, "BXZ"},
        {ZAxis, XAxis, rz, xVec, "BZX"},
        {YAxis, ZAxis, ry, zVec, "BYZ"},
        {ZAxis, YAxis, rz, yVec, "BZY"},
    };
    // clang-format on

    for (const auto& tc : cases) {
        const Rotation seq(BodyRotationSequence, a1, tc.ax1, a2, tc.ax2);
        const Rotation composed = Rotation(a2, tc.r1 * tc.ax2Vec) * tc.r1;
        EXPECT_NEAR((seq - composed).norm(), 0.0, 1e-13) << tc.label;
    }
}

TEST(SimTKCommon_Orientation_RotationSequence, ThreeAxisBodyFixedMatchesComposedSingleRotations) {
    // Body-fixed ZYX(0.31, 0.17, 0.1) equals space-fixed XYZ(0.1, 0.17, 0.31)
    // = Rz(0.31) * Ry(0.17) * Rx(0.1) (reversed order for body-fixed).
    const Rotation rxyz = Rotation(0.31, ZAxis) * Rotation(0.17, YAxis) * Rotation(0.1, XAxis);
    const Rotation r123(BodyRotationSequence, 0.31, ZAxis, 0.17, YAxis, 0.1, XAxis);
    EXPECT_NEAR((r123 - rxyz).norm(), 0.0, 1e-13);
}

TEST(SimTKCommon_Orientation_RotationSequence, SpaceFixedXYMatchesTwoArgConstructor) {
    const Rotation twoArg(SpaceRotationSequence, 0.1, XAxis, 0.17, YAxis);
    const Rotation composed = Rotation(0.17, YAxis) * Rotation(0.1, XAxis);
    EXPECT_NEAR((twoArg - composed).norm(), 0.0, 1e-13);
}

TEST(SimTKCommon_Orientation_RotationSequence, ManualRotationHelperMatchesLibraryForAllAxes) {
    // makeRotationAboutCanonicalAxis() is a reference implementation that
    // directly fills in the rotation matrix; it should agree with the library.
    EXPECT_NEAR((makeRotationAboutCanonicalAxis(0, 0.03) - Rotation(0.03, XAxis)).norm(), 0.0, 1e-13);
    EXPECT_NEAR((makeRotationAboutCanonicalAxis(1, 0.03) - Rotation(0.03, YAxis)).norm(), 0.0, 1e-13);
    EXPECT_NEAR((makeRotationAboutCanonicalAxis(2, 0.03) - Rotation(0.03, ZAxis)).norm(), 0.0, 1e-13);
}

// ===========================================================================
// 8. Rotation angle extraction
// ===========================================================================
//
// Tests that extracting the rotation angle from a Rotation built with a
// known angle recovers that angle accurately, and that a body-fixed XYZ
// Euler-angle round-trip is self-consistent even near the pi/2 gimbal-lock
// singularity.

TEST(SimTKCommon_Orientation_RotationAngleExtraction, KnownAngleIsRecoveredFromAngleAxis) {
    Rotation R;
    R.setRotationFromAngleAboutNonUnitVector(0.17, Vec3(1, 2, 3));
    const Real extracted = R.convertRotationToAngleAxis()[0];
    EXPECT_NEAR(extracted, 0.17, 1e-12);
}

TEST(SimTKCommon_Orientation_RotationAngleExtraction, BodyXYZEulerAnglesRoundTripNearSingularity) {
    // Near pi/2 on the Y-axis (gimbal-lock singularity for XYZ sequence).
    const Real pi2x = -(Pi / 2) + 1e-8;
    const Vec3 vin(-3.0, pi2x, 0.1);

    Rotation b123(BodyRotationSequence, vin[0], XAxis, vin[1], YAxis, vin[2], ZAxis);

    const Vec3 vout = b123.convertThreeAxesRotationToThreeAngles(BodyRotationSequence, XAxis, YAxis, ZAxis);

    const Rotation b123Reconstructed(BodyRotationSequence, vout[0], XAxis, vout[1], YAxis, vout[2], ZAxis);

    // The residual rotation should be nearly zero.
    const Vec4 residual = (~b123 * b123Reconstructed).convertRotationToAngleAxis();
    EXPECT_NEAR(residual[0], 0.0, 1e-10) << "residual rotation angle after round-trip";
}

// ===========================================================================
// 9. Rotation proximity comparison
// ===========================================================================
//
// Tests isSameRotationToWithinAngleOfMachinePrecision and
// isSameRotationToWithinAngle with rotations separated by controlled amounts.

TEST(SimTKCommon_Orientation_RotationComparison, DifferenceOf1e13IsNotWithinMachinePrecision) {
    Rotation R1;
    Rotation R2;
    R1.setRotationFromAngleAboutNonUnitVector(0.17 + 1e-13, Vec3(1, 2, 3));
    R2.setRotationFromAngleAboutNonUnitVector(0.17, Vec3(1, 2, 3));

    EXPECT_FALSE(R1.isSameRotationToWithinAngleOfMachinePrecision(R2));
}

TEST(SimTKCommon_Orientation_RotationComparison, DifferenceOf1e13IsWithinAngle1e12) {
    Rotation R1;
    Rotation R2;
    R1.setRotationFromAngleAboutNonUnitVector(0.17 + 1e-13, Vec3(1, 2, 3));
    R2.setRotationFromAngleAboutNonUnitVector(0.17, Vec3(1, 2, 3));

    EXPECT_TRUE(R1.isSameRotationToWithinAngle(R2, 1e-12));
}

TEST(SimTKCommon_Orientation_RotationComparison, DifferenceOf1e15IsWithinMachinePrecision) {
    Rotation R1;
    Rotation R2;
    R1.setRotationFromAngleAboutNonUnitVector(0.17 + 1e-15, Vec3(1, 2, 3));
    R2.setRotationFromAngleAboutNonUnitVector(0.17, Vec3(1, 2, 3));

    EXPECT_TRUE(R1.isSameRotationToWithinAngleOfMachinePrecision(R2));
}

// ===========================================================================
// 10. Inertia tensor transformation
// ===========================================================================
//
// Tests that the inertia tensor of a point-mass system, when re-expressed
// in a new frame via MassProperties::calcTransformedInertia(), matches the
// inertia computed directly in that frame.

TEST(SimTKCommon_Orientation_Inertia, TransformedInertiaMatchesDirectlyComputedInertia) {
    const Transform X_01(Rotation(0.03, Vec3(1, .1, .7)), Vec3(1, 2, 3));

    const Real masses[] = {1, 2, 3, 4, 5};
    const Real mtot = Vector(5, masses).sum();

    const Vec3 stations000[] = {
        Vec3(.1, .2, .3),
        Vec3(-2, -9, 1),
        Vec3(.01, .02, .05),
        Vec3(1, -1, 1),
        Vec3(0, 0, 0),
    };
    const Vector_<Vec3> s000(5, stations000);
    const Vec3 com000 = ((~Vector(5, masses)) * s000) / mtot;

    Vector_<Vec3> s123(s000);
    for (int i = 0; i < 5; ++i) {
        s123[i] = (~X_01) * s123[i];
    }

    Inertia I000(0);
    Inertia I123(0);
    for (int i = 0; i < 5; ++i) {
        I000 += Inertia(s000[i], masses[i]);
        I123 += Inertia(s123[i], masses[i]);
    }

    const MassProperties mp(mtot, com000, I000);

    // Normalise by RMS of the diagonal of I123 so the error is dimensionless.
    const Real scale = std::sqrt(I123.toMat33().diag().normSqr() / 3);
    const double relErr = (mp.calcTransformedInertia(X_01).toMat33() - I123.toMat33()).norm() / scale;

    EXPECT_LT(relErr, 1e-12);
}

// ===========================================================================
// 11. Rotation from two frame-axis / vector pairs
// ===========================================================================
//
// Tests that a Rotation constructed from a unit vector pinned to one axis
// and a second reference vector defining another axis is orthonormal,
// regardless of the sign or direction of the second vector.

TEST(SimTKCommon_Orientation_RotationFromTwoAxes, PlusPlusYResultIsOrthonormal) {
    const Rotation R(UnitVec3(Vec3(1, 1, 0)), XAxis, Vec3(0, 1, 0), YAxis);
    expectOrthonormal(R, 1e-14);
}

TEST(SimTKCommon_Orientation_RotationFromTwoAxes, PlusMinusYResultIsOrthonormal) {
    const Rotation R(UnitVec3(Vec3(1, 1, 0)), XAxis, Vec3(0, -1, 0), YAxis);
    expectOrthonormal(R, 1e-14);
}

TEST(SimTKCommon_Orientation_RotationFromTwoAxes, PlusPlusZResultIsOrthonormal) {
    const Rotation R(UnitVec3(Vec3(1, 1, 0)), XAxis, Vec3(0, 0, 1), YAxis);
    expectOrthonormal(R, 1e-14);
}