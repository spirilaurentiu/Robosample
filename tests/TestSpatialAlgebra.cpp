// ============================================================================
//  test_spatial_algebra.cpp -- Part 1.6: SpatialVec, Phi/~Phi, inertias
//  Phi/~Phi assertions are the EXACT Simbody SpatialAlgebraTest::testPhiMatrix
//  goldens. ArticulatedInertia::shift == Phi P ~Phi.
// ============================================================================
#include "TestHelpers.hpp"
#include "robot_math.hpp"

using namespace rtest;
using robo::ArticulatedInertia;
using robo::dot;
using robo::Mat33;
using robo::PhiMatrix;
using robo::SpatialInertia;
using robo::SpatialVec;
using robo::SymMat33;
using robo::UnitInertia;
using robo::Vec3;

namespace {
SpatialVec randSV(Rng& rng) {
    return SpatialVec(rng.vec3(), rng.vec3());
}
} // namespace

// --- index ordering and the component pairing dot ---------------------------
TEST(SpatialAlgebra, IndexOrderingAndDot) {
    const Vec3 ang(1, 2, 3), lin(4, 5, 6);
    SpatialVec v(ang, lin);
    EXPECT_TRUE(NearVec3(v[0], ang, kTight)); // [0] = angular
    EXPECT_TRUE(NearVec3(v[1], lin, kTight)); // [1] = linear
    SpatialVec f(Vec3(7, 8, 9), Vec3(10, 11, 12));
    // dot = omega.tau + v.f
    EXPECT_NEAR(dot(v, f), dot(ang, f.angular) + dot(lin, f.linear), kTight);
}

// --- Phi golden (force/inertia shift): Phi(l)[a;b] = [a + l%b; b] ------------
TEST(SpatialAlgebra, PhiForceShiftGolden) {
    Rng rng(1);
    for (int t = 0; t < 300; ++t) {
        const Vec3 l = rng.vec3();
        const SpatialVec v = randSV(rng);
        const SpatialVec got = PhiMatrix(l) * v;
        EXPECT_TRUE(NearVec3(got[0], v[0] + (l % v[1]), kTight));
        EXPECT_TRUE(NearVec3(got[1], v[1], kTight));
    }
}

// --- ~Phi golden (velocity shift): ~Phi(l)[a;b] = [a; b + a%l] ---------------
TEST(SpatialAlgebra, PhiTransposeVelocityShiftGolden) {
    Rng rng(2);
    for (int t = 0; t < 300; ++t) {
        const Vec3 l = rng.vec3();
        const SpatialVec v = randSV(rng);
        const SpatialVec got = (~PhiMatrix(l)) * v;
        EXPECT_TRUE(NearVec3(got[0], v[0], kTight));
        EXPECT_TRUE(NearVec3(got[1], v[1] + (v[0] % l), kTight)); // == v[1] - l%v[0]
    }
}

// --- SpatialInertia * V agrees with ArticulatedInertia(SpatialInertia) * V ---
TEST(SpatialAlgebra, SpatialInertiaSingleSourceOfTruth) {
    Rng rng(3);
    for (int t = 0; t < 200; ++t) {
        const Real m = rng.uniform(0.5, 4);
        const Vec3 com = rng.vec3(-1, 1);
        // a valid (SPD) unit inertia: diagonal + small off-diagonal
        const UnitInertia G(rng.uniform(1, 3),
                            rng.uniform(0, 0.3),
                            rng.uniform(1, 3),
                            rng.uniform(0, 0.3),
                            rng.uniform(0, 0.3),
                            rng.uniform(1, 3));
        const SpatialInertia Mk(m, com, G);
        const SpatialVec v = randSV(rng);
        const SpatialVec a = Mk * v;
        const SpatialVec b = ArticulatedInertia(Mk) * v;
        EXPECT_TRUE(NearVec3(a[0], b[0], kTight));
        EXPECT_TRUE(NearVec3(a[1], b[1], kTight));
    }
}

// --- ArticulatedInertia block action matches the documented block form ------
TEST(SpatialAlgebra, ArticulatedInertiaBlockAction) {
    Rng rng(4);
    for (int t = 0; t < 200; ++t) {
        const SymMat33 J(rng.uniform(1, 3),
                         rng.uniform(0, 0.2),
                         rng.uniform(1, 3),
                         rng.uniform(0, 0.2),
                         rng.uniform(0, 0.2),
                         rng.uniform(1, 3));
        Mat33 F;
        for (int i = 0; i < 9; ++i) {
            F.elems[static_cast<std::size_t>(i)] = rng.uniform(-1, 1);
        }
        const SymMat33 M(rng.uniform(1, 3),
                         rng.uniform(0, 0.2),
                         rng.uniform(1, 3),
                         rng.uniform(0, 0.2),
                         rng.uniform(0, 0.2),
                         rng.uniform(1, 3));
        const ArticulatedInertia P(M, F, J); // (massBlock, momentBlock, inertiaBlock)
        const SpatialVec v = randSV(rng);
        const SpatialVec got = P * v;
        // out.angular = J*w + F*v ; out.linear = F^T*w + M*v
        EXPECT_TRUE(NearVec3(got[0], (J * v.angular) + (F * v.linear), kAlg));
        EXPECT_TRUE(NearVec3(got[1], F.transposeTimes(v.angular) + (M * v.linear), kAlg));
    }
}

// --- shift(offset) == Phi(offset) P ~Phi(offset)  (operator-applied) --------
TEST(SpatialAlgebra, ArticulatedInertiaShiftEqualsPhiTriple) {
    Rng rng(5);
    for (int t = 0; t < 200; ++t) {
        const Real m = rng.uniform(0.5, 4);
        const Vec3 com = rng.vec3(-0.5, 0.5);
        const UnitInertia G(rng.uniform(1, 3),
                            rng.uniform(0, 0.2),
                            rng.uniform(1, 3),
                            rng.uniform(0, 0.2),
                            rng.uniform(0, 0.2),
                            rng.uniform(1, 3));
        const ArticulatedInertia P(SpatialInertia(m, com, G));
        const Vec3 off = rng.vec3(-1, 1);
        const ArticulatedInertia Pshift = P.shift(off);

        for (int k = 0; k < 5; ++k) {
            const SpatialVec v = randSV(rng);
            const SpatialVec lhs = Pshift * v;
            const SpatialVec rhs = PhiMatrix(off) * (P * ((~PhiMatrix(off)) * v));
            EXPECT_TRUE(NearVec3(lhs[0], rhs[0], kAlg));
            EXPECT_TRUE(NearVec3(lhs[1], rhs[1], kAlg));
        }
    }
}

// --- additivity of articulated inertias -------------------------------------
TEST(SpatialAlgebra, ArticulatedInertiaAdditive) {
    Rng rng(6);
    const UnitInertia G(2, 0, 2, 0, 0, 2);
    const ArticulatedInertia A(SpatialInertia(1.0, Vec3(0.1, 0, 0), G));
    const ArticulatedInertia B(SpatialInertia(2.0, Vec3(0, 0.2, 0), G));
    ArticulatedInertia sum = A;
    sum += B;
    for (int t = 0; t < 50; ++t) {
        const SpatialVec v = randSV(rng);
        EXPECT_TRUE(NearVec3((sum * v)[0], (A * v)[0] + (B * v)[0], kAlg));
        EXPECT_TRUE(NearVec3((sum * v)[1], (A * v)[1] + (B * v)[1], kAlg));
    }
}
