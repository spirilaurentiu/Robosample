// ============================================================================
//  test_geometry.cpp -- Part 1.9: calcDihedralAngle, pitch, safeLogSineSqr
// ============================================================================
#include <cmath>

#include "TestHelpers.hpp"
#include "engine_helpers.hpp"
#include "robot_math.hpp"

using namespace rtest;
using robo::calcDihedralAngle;
using robo::Quat;
using robo::Rotation;
using robo::Vec3;

namespace {
// Build 4 points realizing a given dihedral phi about a central bond on +x.
void dihedralPoints(Real phi, Vec3& a, Vec3& b, Vec3& c, Vec3& d) {
    b = Vec3(0, 0, 0);
    c = Vec3(1, 0, 0);                             // central bond along +x
    a = b + Vec3(0, 1, 0);                         // first arm in +y from b
    d = c + Vec3(0, std::cos(phi), std::sin(phi)); // arm rotated by phi about x
}
} // namespace

TEST(Geometry, DihedralKnownAngles) {
    // Construction: d is rotated by +phi about the b->c (+x) axis (right-hand).
    // calcDihedralAngle follows the standard (IUPAC) convention, which for this
    // parameterization returns -phi. We assert magnitude AND the code's sign
    // convention explicitly (so a future sign regression is caught).
    for (Real phi : {0.0, 0.3, 1.0, 2.0, -0.7, M_PI / 2}) {
        Vec3 a, b, c, d;
        dihedralPoints(phi, a, b, c, d);
        const Real got = calcDihedralAngle(a, b, c, d);
        EXPECT_NEAR(std::abs(got), std::abs(phi), kAlg) << "phi=" << phi;
        if (std::abs(phi) > 1e-6 && std::abs(std::abs(phi) - M_PI) > 1e-6) {
            EXPECT_NEAR(got, -phi, kAlg) << "code sign convention (IUPAC), phi=" << phi;
        }
    }
}

TEST(Geometry, DihedralRigidMotionInvariant) {
    Rng rng(1);
    Vec3 a, b, c, d;
    dihedralPoints(0.9, a, b, c, d);
    const Real ref = calcDihedralAngle(a, b, c, d);
    for (int t = 0; t < 50; ++t) {
        const Rotation R = rng.rotation();
        const Vec3 sh = rng.vec3();
        const Real got = calcDihedralAngle(R * a + sh, R * b + sh, R * c + sh, R * d + sh);
        EXPECT_NEAR(got, ref, kAlg);
    }
}

TEST(Geometry, DihedralChiralityFlipsSign) {
    Vec3 a, b, c, d;
    dihedralPoints(0.9, a, b, c, d);
    const Real ref = calcDihedralAngle(a, b, c, d);
    // reflect through the x-z plane (y -> -y): chirality flips
    auto refl = [](const Vec3& v) {
        return Vec3(v[0], -v[1], v[2]);
    };
    const Real got = calcDihedralAngle(refl(a), refl(b), refl(c), refl(d));
    EXPECT_NEAR(got, -ref, kAlg);
}

// --- pitch extraction: sinPitch = 2(qw*qy - qz*qx) --------------------------
TEST(Geometry, PitchFromRotation) {
    for (Real theta = -1.2; theta < 1.2; theta += 0.1) {
        // pure pitch = rotation about Y by theta
        Rotation R;
        R.setRotationFromAngleAboutAxis(theta, robo::YAxis);
        Real w, x, y, z;
        EngineHelpers::rotationToQuaternion(R, w, x, y, z);
        const Real sinPitch = 2 * (w * y - z * x);
        EXPECT_NEAR(std::asin(std::max(Real(-1), std::min(Real(1), sinPitch))), theta, kAlg)
            << "theta=" << theta;
    }
}

// --- safeLogSineSqr: floor at the pole, exact away from it -------------------
TEST(Geometry, SafeLogSineSqrFloor) {
    using EngineHelpers::safeLogSineSqr;
    // at pitch = pi/2, sin^2 = 1 -> log = 0
    EXPECT_NEAR(safeLogSineSqr(M_PI / 2), 0.0, kAlg);
    // away from pole: equals ln(sin^2)
    for (Real p : {0.3, 0.6, 1.0, 1.4}) {
        const Real s = std::sin(p);
        EXPECT_NEAR(safeLogSineSqr(p), std::log(s * s), kAlg);
    }
    // at the pole: floored, finite, equals ln(1e-12)
    EXPECT_NEAR(safeLogSineSqr(0.0), std::log(1e-12), kAlg);
    EXPECT_TRUE(std::isfinite(safeLogSineSqr(0.0)));
    EXPECT_TRUE(std::isfinite(safeLogSineSqr(M_PI)));
}
