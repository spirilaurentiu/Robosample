// ============================================================================
//  TestOrientation.cpp -- port of Simbody SimTKcommon/tests/OrientationTest.cpp
//
//  The original is an ad-hoc, cout-based exploratory program (it prints rather
//  than asserts) built around SimTK-only API that the modern port does not
//  have: Quaternion::setQuaternionFromAngleAxis / convertQuaternionToAngleAxis,
//  Rotation::convertRotationToAngleAxis, the BodyRotationSequence /
//  SpaceRotationSequence Euler-angle constructors,
//  convertThreeAxesRotationToThreeAngles, setRotationFromAngleAboutNonUnitVector,
//  isSameRotationToWithinAngle, UnitVec3::perp(), Transform::asMat34/x()/y()/z(),
//  etc. None of that machinery exists in robot_math.hpp and none is needed by
//  the internal-coordinate dynamics.
//
//  WHAT IS WORTH PORTING is the one piece of real, checkable algebra buried in
//  the original: a space-fixed two-axis rotation sequence equals the reversed
//  product of the two single-axis rotations, e.g.
//        Rotation(SpaceSeq, a, X, b, Y) == Rotation(b,Y) * Rotation(a,X).
//  The original validates that against six hand-derived closed-form matrices
//  (aboutXThenOldY & friends). The modern port has no sequence constructor, but
//  it has Rotation(angle, axis) and Mat33 composition, so we assert the product
//  Rotation(b,axis2)*Rotation(a,axis1) against those same independent closed
//  forms. This is a third-source check of Rotation(angle,axis) AND Mat33::* that
//  the existing single-axis sweep (TestRotationConstruction) cannot give on its
//  own. The two-axis frame construction (chrisRot) and single-axis closed forms
//  from the original are already covered there, so they are not repeated.
// ============================================================================
#include <cmath>

#include "TestHelpers.hpp"
#include "robot_math.hpp"

using namespace rtest;
using robo::CoordinateAxis;
using robo::Mat33;
using robo::Rotation;
using robo::XAxis;
using robo::YAxis;
using robo::ZAxis;

namespace {
// Independent closed forms for "rotate about <first> by f, then about the OLD
// (space-fixed) <second> by s" == R(second, s) * R(first, f). Copied from the
// original OrientationTest helper matrices (aboutAThenOldB).
Mat33 aboutXThenOldY(Real x, Real y) {
    const Real s0 = std::sin(x), c0 = std::cos(x), s1 = std::sin(y), c1 = std::cos(y);
    return Mat33(c1, s0 * s1, c0 * s1, 0, c0, -s0, -s1, s0 * c1, c0 * c1);
}
Mat33 aboutZThenOldX(Real z, Real x) {
    const Real s0 = std::sin(z), c0 = std::cos(z), s1 = std::sin(x), c1 = std::cos(x);
    return Mat33(c0, -s0, 0, s0 * c1, c0 * c1, -s1, s0 * s1, c0 * s1, c1);
}
Mat33 aboutYThenOldZ(Real y, Real z) {
    const Real s0 = std::sin(y), c0 = std::cos(y), s1 = std::sin(z), c1 = std::cos(z);
    return Mat33(c0 * c1, -s1, s0 * c1, c0 * s1, c1, s0 * s1, -s0, 0, c0);
}
Mat33 aboutYThenOldX(Real y, Real x) {
    const Real s0 = std::sin(y), c0 = std::cos(y), s1 = std::sin(x), c1 = std::cos(x);
    return Mat33(c0, 0, s0, s0 * s1, c1, -c0 * s1, -s0 * c1, s1, c0 * c1);
}
Mat33 aboutXThenOldZ(Real x, Real z) {
    const Real s0 = std::sin(x), c0 = std::cos(x), s1 = std::sin(z), c1 = std::cos(z);
    return Mat33(c1, -c0 * s1, s0 * s1, s1, c0 * c1, -s0 * c1, 0, s0, c0);
}
Mat33 aboutZThenOldY(Real z, Real y) {
    const Real s0 = std::sin(z), c0 = std::cos(z), s1 = std::sin(y), c1 = std::cos(y);
    return Mat33(c0 * c1, -s0 * c1, s1, s0, c0, 0, -c0 * s1, s0 * s1, c1);
}

Rotation rot(Real angle, CoordinateAxis axis) {
    Rotation r;
    r.setRotationFromAngleAboutAxis(angle, axis);
    return r;
}
} // namespace

// --- space-fixed two-axis composition == reversed single-axis product -------
// For each ordered axis pair, R(second,s)*R(first,f) must equal the independent
// hand-derived closed form, over a small angle grid.
TEST(Orientation, SpaceFixedTwoAxisComposition) {
    const Real angles[] = {-1.27, -0.29, 0.0, 0.13, 0.9, 2.4};
    for (Real a : angles) {
        for (Real b : angles) {
            EXPECT_TRUE(NearMat33(rot(b, YAxis) * rot(a, XAxis), aboutXThenOldY(a, b), kAlg))
                << "X then old Y, a=" << a << " b=" << b;
            EXPECT_TRUE(NearMat33(rot(b, XAxis) * rot(a, ZAxis), aboutZThenOldX(a, b), kAlg))
                << "Z then old X, a=" << a << " b=" << b;
            EXPECT_TRUE(NearMat33(rot(b, ZAxis) * rot(a, YAxis), aboutYThenOldZ(a, b), kAlg))
                << "Y then old Z, a=" << a << " b=" << b;
            EXPECT_TRUE(NearMat33(rot(b, XAxis) * rot(a, YAxis), aboutYThenOldX(a, b), kAlg))
                << "Y then old X, a=" << a << " b=" << b;
            EXPECT_TRUE(NearMat33(rot(b, ZAxis) * rot(a, XAxis), aboutXThenOldZ(a, b), kAlg))
                << "X then old Z, a=" << a << " b=" << b;
            EXPECT_TRUE(NearMat33(rot(b, YAxis) * rot(a, ZAxis), aboutZThenOldY(a, b), kAlg))
                << "Z then old Y, a=" << a << " b=" << b;
        }
    }
}

// --- three-axis space-fixed sequence: R(c,Z) R(b,Y) R(a,X) is a proper rotation
// and equals the product of the three closed-form single-axis matrices. (The
// original builds r123 = rz*ry*rx and a BodyRotationSequence; we keep only the
// part expressible without a sequence ctor: associated product is orthonormal.)
TEST(Orientation, ThreeAxisSpaceFixedProductIsProperRotation) {
    const Real a = 0.1, b = 0.17, c = 0.31;
    const Rotation R = rot(c, ZAxis) * rot(b, YAxis) * rot(a, XAxis);
    EXPECT_TRUE(NearMat33(R * R.transpose(), Mat33::identity(), kAlg));
    // det(R) == +1
    const Real det = R(0, 0) * (R(1, 1) * R(2, 2) - R(1, 2) * R(2, 1))
                     - R(0, 1) * (R(1, 0) * R(2, 2) - R(1, 2) * R(2, 0))
                     + R(0, 2) * (R(1, 0) * R(2, 1) - R(1, 1) * R(2, 0));
    EXPECT_NEAR(det, 1.0, kAlg);
    // reverse-order body composition associativity: (Rz Ry) Rx == Rz (Ry Rx)
    EXPECT_TRUE(NearMat33((rot(c, ZAxis) * rot(b, YAxis)) * rot(a, XAxis),
                          rot(c, ZAxis) * (rot(b, YAxis) * rot(a, XAxis)),
                          kAlg));
}