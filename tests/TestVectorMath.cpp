// ============================================================================
//  test_vec_mat.cpp -- Part 1.1: Vec3 / Mat33 / SymMat33 / crossMat algebra
// ============================================================================
#include "TestHelpers.hpp"
#include "robot_math.hpp"

using namespace rtest;
using robo::crossMat;
using robo::dot;
using robo::Mat33;
using robo::SymMat33;
using robo::Vec3;

// --- cross product: right-hand axis identities + algebraic laws -------------
TEST(VecMat, CrossProductRightHandRule) {
    const Vec3 x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
    EXPECT_TRUE(NearVec3(x % y, z, kTight));
    EXPECT_TRUE(NearVec3(y % z, x, kTight));
    EXPECT_TRUE(NearVec3(z % x, y, kTight));
}

TEST(VecMat, CrossProductLaws) {
    Rng rng(1);
    for (int t = 0; t < 500; ++t) {
        const Vec3 a = rng.vec3(), b = rng.vec3();
        EXPECT_TRUE(NearVec3(a % b, Vec3(0, 0, 0) - (b % a), kTight)); // antisymmetry
        EXPECT_TRUE(NearVec3(a % a, Vec3(0, 0, 0), kTight));           // self -> 0
        EXPECT_NEAR(dot(a, a % b), 0.0, kTight);                       // perp to a
        EXPECT_NEAR(dot(b, a % b), 0.0, kTight);                       // perp to b
    }
}

// --- crossMat is the matrix form of the cross product (load-bearing) --------
TEST(VecMat, CrossMatEqualsCrossProduct) {
    Rng rng(2);
    for (int t = 0; t < 500; ++t) {
        const Vec3 v = rng.vec3(), x = rng.vec3();
        EXPECT_TRUE(NearVec3(crossMat(v) * x, v % x, kTight));
    }
}

TEST(VecMat, CrossMatIsSkewSymmetric) {
    Rng rng(3);
    for (int t = 0; t < 100; ++t) {
        const Vec3 v = rng.vec3();
        const Mat33 S = crossMat(v);
        EXPECT_TRUE(NearMat33(S, Mat33(0, 0, 0, 0, 0, 0, 0, 0, 0) - S.transpose(), kTight));
        EXPECT_TRUE(NearVec3(S * v, Vec3(0, 0, 0), kTight)); // S v = v x v = 0
    }
}

// --- Mat33 multiply / transpose -------------------------------------------
TEST(VecMat, MatMulAssociativeAndTranspose) {
    Rng rng(4);
    for (int t = 0; t < 100; ++t) {
        Mat33 A, B, C;
        for (int i = 0; i < 9; ++i) {
            A.elems[static_cast<std::size_t>(i)] = rng.uniform(-2, 2);
            B.elems[static_cast<std::size_t>(i)] = rng.uniform(-2, 2);
            C.elems[static_cast<std::size_t>(i)] = rng.uniform(-2, 2);
        }
        EXPECT_TRUE(NearMat33((A * B) * C, A * (B * C), kAlg));
        EXPECT_TRUE(NearMat33(A.transpose().transpose(), A, kTight));
        // (AB)^T = B^T A^T
        EXPECT_TRUE(NearMat33((A * B).transpose(), B.transpose() * A.transpose(), kAlg));
    }
}

TEST(VecMat, MatVecMatchesManual) {
    Rng rng(5);
    for (int t = 0; t < 100; ++t) {
        Mat33 A;
        for (int i = 0; i < 9; ++i) {
            A.elems[static_cast<std::size_t>(i)] = rng.uniform(-2, 2);
        }
        const Vec3 x = rng.vec3();
        const Vec3 y = A * x;
        for (int r = 0; r < 3; ++r) {
            const Real m = A(r, 0) * x[0] + A(r, 1) * x[1] + A(r, 2) * x[2];
            EXPECT_NEAR(y[r], m, kTight);
        }
    }
}

// --- SymMat33 storage <-> full round-trip and action -----------------------
TEST(VecMat, SymMatStorageRoundTrip) {
    Rng rng(6);
    for (int t = 0; t < 100; ++t) {
        const Real xx = rng.uniform(-2, 2), xy = rng.uniform(-2, 2), yy = rng.uniform(-2, 2);
        const Real xz = rng.uniform(-2, 2), yz = rng.uniform(-2, 2), zz = rng.uniform(-2, 2);
        const SymMat33 S(xx, xy, yy, xz, yz, zz);
        const Mat33 F = S.full();
        // full() must be symmetric and carry the right entries
        EXPECT_NEAR(F(0, 0), xx, kTight);
        EXPECT_NEAR(F(0, 1), xy, kTight);
        EXPECT_NEAR(F(1, 0), xy, kTight);
        EXPECT_NEAR(F(0, 2), xz, kTight);
        EXPECT_NEAR(F(2, 0), xz, kTight);
        EXPECT_NEAR(F(1, 2), yz, kTight);
        EXPECT_NEAR(F(2, 1), yz, kTight);
        EXPECT_NEAR(F(2, 2), zz, kTight);
        // round-trip fromSymmetricPart(full(S)) == S
        const SymMat33 back = SymMat33::fromSymmetricPart(F);
        for (int i = 0; i < 6; ++i) {
            EXPECT_NEAR(back.elems[static_cast<std::size_t>(i)],
                        S.elems[static_cast<std::size_t>(i)],
                        kTight);
        }
    }
}

TEST(VecMat, SymMatTimesVecEqualsFull) {
    Rng rng(7);
    for (int t = 0; t < 100; ++t) {
        const SymMat33 S(rng.uniform(-2, 2),
                         rng.uniform(-2, 2),
                         rng.uniform(-2, 2),
                         rng.uniform(-2, 2),
                         rng.uniform(-2, 2),
                         rng.uniform(-2, 2));
        const Vec3 x = rng.vec3();
        EXPECT_TRUE(NearVec3(S * x, S.full() * x, kTight));
    }
}

// --- vector arithmetic basics ----------------------------------------------
TEST(VecMat, VectorArithmetic) {
    Rng rng(8);
    for (int t = 0; t < 100; ++t) {
        const Vec3 a = rng.vec3(), b = rng.vec3();
        EXPECT_TRUE(NearVec3((a + b) - b, a, kTight));
        EXPECT_TRUE(NearVec3(a * 2.0, a + a, kTight));
        EXPECT_NEAR((a * 3.0).norm(), 3.0 * a.norm(), kAlg);
        EXPECT_NEAR(dot(a, b), a[0] * b[0] + a[1] * b[1] + a[2] * b[2], kTight);
    }
}
