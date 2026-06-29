// ============================================================================
//  test_transform.cpp -- Part 1.2: Transform composition and inversion
// ============================================================================
#include "TestHelpers.hpp"
#include "robot_math.hpp"

using namespace rtest;
using robo::Mat33;
using robo::Rotation;
using robo::Transform;
using robo::Vec3;

TEST(Transform, InverseIsIdentity) {
    Rng rng(1);
    for (int t = 0; t < 200; ++t) {
        const Transform X(rng.rotation(), rng.vec3());
        const Transform I1 = X.inverse() * X;
        const Transform I2 = X * X.inverse();
        EXPECT_TRUE(NearMat33(I1.R(), Mat33::identity(), kTight));
        EXPECT_TRUE(NearVec3(I1.p(), Vec3(0, 0, 0), kTight));
        EXPECT_TRUE(NearMat33(I2.R(), Mat33::identity(), kTight));
        EXPECT_TRUE(NearVec3(I2.p(), Vec3(0, 0, 0), kTight));
    }
}

// Directed inverse: inverse(R,p) MUST be (R^T, -R^T p). A wrong-but-cancelling
// inverse still gives ~X*X=I, so assert the blocks explicitly.
TEST(Transform, InverseBlocksExplicit) {
    Rng rng(2);
    for (int t = 0; t < 200; ++t) {
        const Rotation R = rng.rotation();
        const Vec3 p = rng.vec3();
        const Transform X(R, p);
        const Transform Xi = X.inverse();
        EXPECT_TRUE(NearMat33(Xi.R(), R.transpose(), kTight));
        EXPECT_TRUE(NearVec3(Xi.p(), Vec3(0, 0, 0) - (R.transpose() * p), kTight));
    }
}

TEST(Transform, CompositionAssociativeAndActsOnPoint) {
    Rng rng(3);
    for (int t = 0; t < 200; ++t) {
        const Transform X(rng.rotation(), rng.vec3());
        const Transform Y(rng.rotation(), rng.vec3());
        const Transform Z(rng.rotation(), rng.vec3());
        const Transform a = (X * Y) * Z;
        const Transform b = X * (Y * Z);
        EXPECT_TRUE(NearMat33(a.R(), b.R(), kAlg));
        EXPECT_TRUE(NearVec3(a.p(), b.p(), kAlg));

        // (X*Y) applied to a point == X applied to (Y applied to point)
        const Vec3 pt = rng.vec3();
        EXPECT_TRUE(NearVec3((X * Y) * pt, X * (Y * pt), kAlg));
        // action formula
        EXPECT_TRUE(NearVec3(X * pt, (X.R() * pt) + X.p(), kTight));
    }
}

TEST(Transform, OperatorTildeEqualsInverse) {
    Rng rng(4);
    for (int t = 0; t < 100; ++t) {
        const Transform X(rng.rotation(), rng.vec3());
        const Transform a = ~X;
        const Transform b = X.inverse();
        EXPECT_TRUE(NearMat33(a.R(), b.R(), kTight));
        EXPECT_TRUE(NearVec3(a.p(), b.p(), kTight));
    }
}
