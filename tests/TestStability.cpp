// ============================================================================
//  test_stability.cpp -- Part 1.10: orthonormality drift, exp-map norm,
//  scaling, NaN/Inf guards, determinism.
// ============================================================================
#include <cmath>

#include "EngineHelpers.hpp"
#include "TestHelpers.hpp"
#include "robot_math.hpp"

using namespace rtest;
using robo::Mat33;
using robo::Quat;
using robo::Rotation;
using robo::SymMat33;
using robo::Transform;
using robo::Vec3;

// --- composing many rotations does not lose orthonormality ------------------
TEST(Stability, OrthonormalityUnderManyProducts) {
    Rng rng(1);
    Rotation acc;
    for (int i = 0; i < 10000; ++i) {
        acc = Rotation(acc * rng.rotation());
    }
    EXPECT_TRUE(NearMat33(acc * acc.transpose(), Mat33::identity(), 1e-10))
        << "orthonormality drifted after 1e4 products";
}

// --- exp-map quaternion advance preserves unit norm and is reversible -------
TEST(Stability, ExpMapNormAndReversibility) {
    Rng rng(2);
    for (int t = 0; t < 200; ++t) {
        Quat q0 = rng.unitQuat();
        const Vec3 w = rng.vec3(-3, 3);
        const Real h = 0.01;
        Real q1[4];
        Real qq0[4] = {q0.elems[0], q0.elems[1], q0.elems[2], q0.elems[3]};
        EngineHelpers::advanceQuatExp(qq0, w, h, q1);
        const Real n1 = std::sqrt(q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2] + q1[3] * q1[3]);
        EXPECT_NEAR(n1, 1.0, 1e-12) << "exp-map must preserve |q|=1 by construction";

        // reversibility: advancing by -w returns to q0 (double-cover aware)
        Real q2[4];
        EngineHelpers::advanceQuatExp(q1, Vec3(0, 0, 0) - w, h, q2);
        EXPECT_LT(quatDist(Quat(q2[0], q2[1], q2[2], q2[3]), q0), 1e-12);
    }
}

// --- exp-map reduces to the linear map as h -> 0 ----------------------------
TEST(Stability, ExpMapReducesToLinearAtSmallH) {
    Rng rng(3);
    for (int t = 0; t < 100; ++t) {
        const Quat q0 = rng.unitQuat();
        const Vec3 w = rng.vec3(-2, 2);
        const Real h = 1e-7;
        Real qq0[4] = {q0.elems[0], q0.elems[1], q0.elems[2], q0.elems[3]};
        Real q1[4];
        EngineHelpers::advanceQuatExp(qq0, w, h, q1);
        // linear: q0 + h * angVelToQdot(q0, w)
        const Quat qd = Quat::angVelToQdot(q0, w);
        for (int r = 0; r < 4; ++r) {
            const Real lin = qq0[r] + h * qd.elems[static_cast<std::size_t>(r)];
            EXPECT_NEAR(q1[r], lin, 1e-12) << "component " << r;
        }
    }
}

// --- length scaling: translations scale, rotations unchanged ----------------
TEST(Stability, LengthScaling) {
    Rng rng(4);
    const Real lambda = 7.3;
    for (int t = 0; t < 100; ++t) {
        const Rotation R = rng.rotation();
        const Vec3 p = rng.vec3();
        const Transform X(R, p);
        const Transform Xs(R, p * lambda);
        // composing scaled transforms scales translation parts consistently
        const Vec3 pt = rng.vec3();
        const Vec3 unscaled = X * pt;
        const Vec3 scaledInput = Xs * (pt * lambda);
        EXPECT_TRUE(NearVec3(scaledInput, unscaled * lambda, kAlg));
    }
}

// --- NaN/Inf guards: normalize, exp-map small-w branch ----------------------
TEST(Stability, GuardsAgainstDegenerateInput) {
    // zero quaternion -> identity (guarded)
    Quat z(0, 0, 0, 0);
    z.normalize();
    EXPECT_TRUE(std::isfinite(z.norm()));
    EXPECT_NEAR(z.elems[0], 1.0, kTight);

    // exp-map with ~zero angular velocity uses the small-angle branch, stays unit
    Real q0[4] = {1, 0, 0, 0}, q1[4];
    EngineHelpers::advanceQuatExp(q0, Vec3(1e-15, 0, 0), 0.01, q1);
    const Real n = std::sqrt(q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2] + q1[3] * q1[3]);
    EXPECT_TRUE(std::isfinite(n));
}

// --- determinism: identical seed -> identical stream ------------------------
TEST(Stability, DeterministicGenerators) {
    Rng a(12345), b(12345);
    for (int t = 0; t < 1000; ++t) {
        const Vec3 va = a.vec3();
        const Vec3 vb = b.vec3();
        EXPECT_EQ(va[0], vb[0]);
        EXPECT_EQ(va[1], vb[1]);
        EXPECT_EQ(va[2], vb[2]);
    }
    // pure transforms are deterministic functions of their inputs
    Rng r(9);
    for (int t = 0; t < 100; ++t) {
        const Quat q = r.unitQuat();
        const Mat33 R1 = Rotation::fromQuaternion(q);
        const Mat33 R2 = Rotation::fromQuaternion(q);
        for (int i = 0; i < 9; ++i) {
            EXPECT_EQ(R1.elems[static_cast<std::size_t>(i)], R2.elems[static_cast<std::size_t>(i)]);
        }
    }
}
