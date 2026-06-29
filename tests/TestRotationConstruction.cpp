// ============================================================================
//  test_rotation_construction.cpp -- Part 1.4: Rotation from angles/axes
//  Exhaustive angle sweeps after Simbody RotationTest; two-axis Gram-Schmidt;
//  collinear fallback; orthonormality and det+1 everywhere.
// ============================================================================
#include <cmath>

#include "TestHelpers.hpp"
#include "robot_math.hpp"

using namespace rtest;
using robo::CoordinateAxis;
using robo::dot;
using robo::Mat33;
using robo::Rotation;
using robo::UnitVec3;
using robo::Vec3;
using robo::XAxis;
using robo::YAxis;
using robo::ZAxis;

namespace {
Real det3(const Mat33& m) {
    return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1))
           - m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0))
           + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
}
void expectProperRotation(const Mat33& R, Real tol) {
    EXPECT_TRUE(NearMat33(R * R.transpose(), Mat33::identity(), tol)) << "not orthonormal";
    EXPECT_NEAR(det3(R), 1.0, tol) << "det != +1 (improper)";
}
} // namespace

// --- single-axis exhaustive sweep, -385..+385 deg in 0.5 deg steps ----------
TEST(RotationConstruction, SingleAxisExhaustiveSweep) {
    const Real deg = M_PI / 180.0;
    for (int axis = 0; axis < 3; ++axis) {
        const CoordinateAxis ax = static_cast<CoordinateAxis>(axis);
        for (Real a = -385 * deg; a < 385 * deg; a += 0.5 * deg) {
            Rotation R;
            R.setRotationFromAngleAboutAxis(a, ax);
            expectProperRotation(R, kAlg);

            // analytic closed form
            const Real c = std::cos(a), s = std::sin(a);
            Mat33 ref = Mat33::identity();
            if (ax == XAxis) {
                ref = Mat33(1, 0, 0, 0, c, -s, 0, s, c);
            }
            if (ax == YAxis) {
                ref = Mat33(c, 0, s, 0, 1, 0, -s, 0, c);
            }
            if (ax == ZAxis) {
                ref = Mat33(c, -s, 0, s, c, 0, 0, 0, 1);
            }
            ASSERT_TRUE(NearMat33(R, ref, kAlg)) << "axis " << axis << " angle " << a;

            // R(a) R(-a) = I
            Rotation Rm;
            Rm.setRotationFromAngleAboutAxis(-a, ax);
            EXPECT_TRUE(NearMat33(R * Rm, Mat33::identity(), kAlg));
        }
    }
}

TEST(RotationConstruction, AboutZMatchesAboutAxisZ) {
    const Real deg = M_PI / 180.0;
    for (Real a = -200 * deg; a < 200 * deg; a += 7 * deg) {
        Rotation R1, R2;
        R1.setRotationFromAngleAboutZ(a);
        R2.setRotationFromAngleAboutAxis(a, ZAxis);
        EXPECT_TRUE(NearMat33(R1, R2, kTight));
    }
}

// --- one-axis construction places the axis and stays orthonormal ------------
TEST(RotationConstruction, OneAxisPlacesAxisOrthonormal) {
    Rng rng(1);
    for (int t = 0; t < 300; ++t) {
        const UnitVec3 dir(rng.unitVec3());
        for (int axis = 0; axis < 3; ++axis) {
            Rotation R;
            R.setRotationFromOneAxis(dir, static_cast<CoordinateAxis>(axis));
            expectProperRotation(R, kAlg);
            // column `axis` equals dir
            const Vec3 col(R(0, axis), R(1, axis), R(2, axis));
            EXPECT_TRUE(NearVec3(col, dir.asVec3(), kAlg));
        }
    }
    // axis-aligned and near-degenerate inputs (the anyUnitPerpendicular branches)
    for (const Vec3 d : {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(-1, 0, 0), Vec3(1e-9, 1, 0)}) {
        Rotation R;
        R.setRotationFromOneAxis(UnitVec3(d), XAxis);
        expectProperRotation(R, kAlg);
    }
}

// --- two-axis Gram-Schmidt construction, all 6 ordered axis pairs -----------
TEST(RotationConstruction, TwoAxisGramSchmidtAllPairs) {
    Rng rng(2);
    const int pairs[6][2] = {{0, 1}, {1, 2}, {2, 0}, {1, 0}, {2, 1}, {0, 2}};
    for (int t = 0; t < 200; ++t) {
        const UnitVec3 primary(rng.unitVec3());
        const Vec3 planeVec = rng.vec3();
        for (auto& pr : pairs) {
            if (pr[0] == pr[1]) {
                continue;
            }
            Rotation R;
            R.setRotationFromTwoAxes(primary,
                                     static_cast<CoordinateAxis>(pr[0]),
                                     planeVec,
                                     static_cast<CoordinateAxis>(pr[1]));
            expectProperRotation(R, kAlg);
            // primary axis is exactly `primary`
            const Vec3 colU(R(0, pr[0]), R(1, pr[0]), R(2, pr[0]));
            EXPECT_TRUE(NearVec3(colU, primary.asVec3(), kAlg));
            // plane axis is perpendicular to primary and in the (primary,planeVec) plane
            const Vec3 colV(R(0, pr[1]), R(1, pr[1]), R(2, pr[1]));
            EXPECT_NEAR(dot(colU, colV), 0.0, kAlg);
        }
    }
}

// --- collinear fallback: planeVec parallel to primary must not produce NaN --
TEST(RotationConstruction, TwoAxisCollinearFallback) {
    const UnitVec3 primary(Vec3(0, 0, 1));
    const Vec3 planeVec(0, 0, 5); // exactly parallel
    Rotation R;
    R.setRotationFromTwoAxes(primary, ZAxis, planeVec, XAxis);
    expectProperRotation(R, kAlg);
    for (int i = 0; i < 9; ++i) {
        EXPECT_TRUE(std::isfinite(R.elems[static_cast<std::size_t>(i)]));
    }
}

// --- directed two-axis case with hand-computed frame (after OrientationTest) -
// primary=(1,1,0)/sqrt2 on X, planeVec=(0,1,0) on Y. Gram-Schmidt gives
// colV=(-1,1,0)/sqrt2; cyclic (X,Y) => colW = colU x colV = (0,0,1).
TEST(RotationConstruction, TwoAxisDirectedFrame) {
    const Real r2 = std::sqrt(2.0);
    Rotation R;
    R.setRotationFromTwoAxes(UnitVec3(Vec3(1, 1, 0)), XAxis, Vec3(0, 1, 0), YAxis);
    EXPECT_TRUE(NearVec3(Vec3(R(0, 0), R(1, 0), R(2, 0)), Vec3(1 / r2, 1 / r2, 0), kAlg));
    EXPECT_TRUE(NearVec3(Vec3(R(0, 1), R(1, 1), R(2, 1)), Vec3(-1 / r2, 1 / r2, 0), kAlg));
    EXPECT_TRUE(NearVec3(Vec3(R(0, 2), R(1, 2), R(2, 2)), Vec3(0, 0, 1), kAlg));
    expectProperRotation(R, kAlg);
}

// --- transpose is the inverse rotation --------------------------------------
TEST(RotationConstruction, TransposeIsInverse) {
    Rng rng(3);
    for (int t = 0; t < 200; ++t) {
        const Rotation R = rng.rotation();
        EXPECT_TRUE(NearMat33(R * (~R), Mat33::identity(), kTight));
        EXPECT_TRUE(NearMat33((~R) * R, Mat33::identity(), kTight));
    }
}
