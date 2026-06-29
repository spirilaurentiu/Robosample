// ============================================================================
//  TestInertia.cpp -- rigid-body mass properties: Inertia, MassProperties,
//  reexpress, mass scaling, and the cross-product / spatial-inertia block-form
//  goldens merged in from the former TestMassProperties.cpp (the two files shared
//  the same value types and fixtures; one file removes the duplication, keeping
//  the Simbody-golden checks as the MassProperties suite and the inertia algebra
//  as the InertiaTest suite).
// ============================================================================
#include "TestHelpers.hpp"
#include "robot_math.hpp"

using namespace rtest;
using robo::ArticulatedInertia;
using robo::crossMat;
using robo::Inertia;
using robo::MassProperties;
using robo::Mat33;
using robo::Rotation;
using robo::SpatialInertia;
using robo::SymMat33;
using robo::UnitInertia;
using robo::Vec3;

// --- point-mass inertia: I = m(|p|^2 I3 - p p^T) ----------------------------
TEST(InertiaTest, PointMassFormula) {
    Rng rng(1);
    for (int t = 0; t < 200; ++t) {
        const Vec3 p = rng.vec3(-2, 2);
        const Real m = rng.uniform(0.5, 5);
        const SymMat33 I = Inertia(p, m).asSymMat33();
        EXPECT_NEAR(I.elems[0], m * (p[1] * p[1] + p[2] * p[2]), kTight); // xx
        EXPECT_NEAR(I.elems[1], -m * p[0] * p[1], kTight);                // xy
        EXPECT_NEAR(I.elems[2], m * (p[0] * p[0] + p[2] * p[2]), kTight); // yy
        EXPECT_NEAR(I.elems[3], -m * p[0] * p[2], kTight);                // xz
        EXPECT_NEAR(I.elems[4], -m * p[1] * p[2], kTight);                // yz
        EXPECT_NEAR(I.elems[5], m * (p[0] * p[0] + p[1] * p[1]), kTight); // zz
    }
}

// --- additivity of point-mass inertias --------------------------------------
TEST(InertiaTest, Additive) {
    const Vec3 p1(1, 0, 0), p2(0, 1, 0);
    const Real m1 = 2, m2 = 3;
    Inertia sum = Inertia(p1, m1);
    sum += Inertia(p2, m2);
    const SymMat33 a = sum.asSymMat33();
    const SymMat33 b = (Inertia(p1, m1) + Inertia(p2, m2)).asSymMat33();
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(a.elems[static_cast<std::size_t>(i)], b.elems[static_cast<std::size_t>(i)], kTight);
    }
}

// --- reexpress is a similarity transform R^T I R ----------------------------
TEST(InertiaTest, ReexpressSimilarity) {
    Rng rng(2);
    for (int t = 0; t < 200; ++t) {
        const SymMat33 I(rng.uniform(1, 4),
                         rng.uniform(-0.3, 0.3),
                         rng.uniform(1, 4),
                         rng.uniform(-0.3, 0.3),
                         rng.uniform(-0.3, 0.3),
                         rng.uniform(1, 4));
        const Rotation R = rng.rotation();
        const SymMat33 Ir = I.reexpress(R);
        const Mat33 ref = R.transpose() * I.full() * R;
        EXPECT_TRUE(NearMat33(Ir.full(), ref, kAlg));
        // trace invariant
        EXPECT_NEAR(Ir.full()(0, 0) + Ir.full()(1, 1) + Ir.full()(2, 2),
                    I.full()(0, 0) + I.full()(1, 1) + I.full()(2, 2),
                    kAlg);
        // round-trip reexpress(R) then reexpress(R^T) == I
        const SymMat33 back = Ir.reexpress(Rotation(R.transpose()));
        for (int i = 0; i < 6; ++i) {
            EXPECT_NEAR(back.elems[static_cast<std::size_t>(i)], I.elems[static_cast<std::size_t>(i)], kAlg);
        }
    }
}

// --- MassProperties::reexpress: com' = R^T com ------------------------------
TEST(InertiaTest, MassPropertiesReexpress) {
    Rng rng(3);
    for (int t = 0; t < 100; ++t) {
        const Real m = rng.uniform(0.5, 4);
        const Vec3 com = rng.vec3(-1, 1);
        const UnitInertia G(2, 0, 2, 0, 0, 2);
        const MassProperties mp(m, com, G);
        const Rotation R = rng.rotation();
        const MassProperties mpr = mp.reexpress(R);
        EXPECT_TRUE(NearVec3(mpr.com, R.transpose() * com, kAlg));
        EXPECT_NEAR(mpr.mass, m, kTight); // mass invariant
    }
}

// --- mass scaling: ArticulatedInertia blocks scale linearly in mass ---------
TEST(InertiaTest, MassScalingLinear) {
    Rng rng(4);
    const Vec3 com(0.1, -0.2, 0.3);
    const UnitInertia G(2, 0.1, 3, 0.1, 0.1, 4);
    const Real m = 1.7, s = 3.5;
    const ArticulatedInertia P(SpatialInertia(m, com, G));
    const ArticulatedInertia Ps(SpatialInertia(s * m, com, G));
    for (int t = 0; t < 50; ++t) {
        const SpatialVec v(rng.vec3(), rng.vec3());
        const SpatialVec a = Ps * v;
        const SpatialVec b = (P * v) * s; // scaling mass scales the operator
        EXPECT_TRUE(NearVec3(a[0], b[0], kAlg));
        EXPECT_TRUE(NearVec3(a[1], b[1], kAlg));
    }
}

// --- zero-mass guard (virtual sites): unit inertia is zero, no divide-by-0 --
TEST(InertiaTest, ZeroMassGuard) {
    const MassProperties mp(0.0, Vec3(1, 2, 3), Inertia(5.0));
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(mp.unitInertia.elems[static_cast<std::size_t>(i)], 0.0, kTight);
    }
    EXPECT_NEAR(mp.mass, 0.0, kTight);
}

// ===========================================================================
//  suite MassProperties -- exact Simbody goldens for cross products, point
//  inertia, and the spatial / articulated inertia block form.
// ===========================================================================

// --- exact cross product & crossMat goldens (testCrossProduct, 3-D part) -----
// All operands are exactly representable in base 2, so these must hold to the
// last bit -- a transposition or sign error in % / crossMat is caught exactly.
TEST(MassProperties, CrossProductExactGolden) {
    const Vec3 w(1.25, 3, -2.5), v(-2.75, 2.125, 5);
    const Vec3 wxv = w % v;
    EXPECT_EQ(wxv[0], 20.3125);
    EXPECT_EQ(wxv[1], 0.625);
    EXPECT_EQ(wxv[2], 10.90625);
    const Vec3 vxw = v % w;
    EXPECT_EQ(vxw[0], -wxv[0]);
    EXPECT_EQ(vxw[1], -wxv[1]);
    EXPECT_EQ(vxw[2], -wxv[2]);
}

TEST(MassProperties, CrossMatExactGolden) {
    const Vec3 w(1.25, 3, -2.5);
    const Mat33 expected(0, 2.5, 3, -2.5, 0, -1.25, -3, 1.25, 0);
    EXPECT_TRUE(NearMat33(crossMat(w), expected, 0.0)); // exact
}

// crossMat(w) is skew: ~crossMat == -crossMat, so crossMat(-w) == ~crossMat(w).
TEST(MassProperties, CrossMatOfNegatedIsTranspose) {
    Rng rng(10);
    for (int t = 0; t < 200; ++t) {
        const Vec3 w = rng.vec3();
        EXPECT_TRUE(NearMat33(crossMat(Vec3(0, 0, 0) - w), crossMat(w).transpose(), kTight));
    }
}

// crossMat^2 * v == w%(w%v); crossMatSq == ~crossMat*crossMat == -crossMat^2.
TEST(MassProperties, CrossMatSquaredAction) {
    Rng rng(11);
    for (int t = 0; t < 300; ++t) {
        const Vec3 w = rng.vec3(), v = rng.vec3();
        const Mat33 wx = crossMat(w);
        EXPECT_TRUE(NearVec3(wx * (wx * v), w % (w % v), kTight));
        const Mat33 crossMatSq = wx.transpose() * wx;
        EXPECT_TRUE(NearVec3(crossMatSq * v, Vec3(0, 0, 0) - (w % (w % v)), kTight));
    }
}

// --- point-mass inertia == m * crossMatSq(p) (testInertia, shift-from-point) -
// Inertia(p,m) is the load-bearing path the engine uses to build body inertia
// from atom point masses: I_point = m (|p|^2 I3 - p p^T) = m * (~crossMat(p) crossMat(p)).
TEST(MassProperties, PointInertiaEqualsCrossMatSquared) {
    Rng rng(12);
    for (int t = 0; t < 300; ++t) {
        const Vec3 p = rng.vec3(-2, 2);
        const Real m = rng.uniform(0.3, 5);
        const Mat33 px = crossMat(p);
        const Mat33 crossMatSq = px.transpose() * px; // = |p|^2 I - p p^T
        const Mat33 got = Inertia(p, m).asSymMat33().full();
        EXPECT_TRUE(NearMat33(got, crossMatSq * m, kTight));
    }
}

// Point inertia is linear in the mass: m * I_unit == I(p,m).
TEST(MassProperties, PointInertiaLinearInMass) {
    Rng rng(13);
    for (int t = 0; t < 200; ++t) {
        const Vec3 p = rng.vec3(-2, 2);
        const Real m = rng.uniform(0.3, 5);
        const SymMat33 unit = Inertia(p, Real(1)).asSymMat33();
        const SymMat33 scaled = unit * m;
        const SymMat33 direct = Inertia(p, m).asSymMat33();
        for (int i = 0; i < 6; ++i) {
            EXPECT_NEAR(scaled.elems[static_cast<std::size_t>(i)],
                        direct.elems[static_cast<std::size_t>(i)],
                        kTight);
        }
    }
}

// --- SpatialInertia block form (testSpatialInertia) -------------------------
// The spatial inertia's single source of truth is ArticulatedInertia(Spatial-
// Inertia); the block goldens are asserted on its J/F/M blocks: (0,0)=m*G,
// (0,1)=m*crossMat(com), (1,1)=m*I.
TEST(MassProperties, SpatialInertiaBlockForm) {
    const Real mass = 1.125; // exact in base 2
    const Vec3 com(0.1, 0.2, 0.25);
    const UnitInertia gyration(1.8, 1.9, 2.1, 0.01, 0.03, 0.02);
    const ArticulatedInertia P(SpatialInertia(mass, com, gyration));

    EXPECT_TRUE(NearMat33(P.angAng.full(), gyration.full() * mass, kTight)); // (0,0) = m*G
    EXPECT_TRUE(NearMat33(P.angLin, crossMat(com) * mass, kTight));          // (0,1) = m*crossMat(com)
    EXPECT_TRUE(NearMat33(P.linLin.full(), Mat33::identity() * mass, kTight)); // (1,1) = m*I3
}

// --- ArticulatedInertia ctor block mapping (testArticulatedInertia) ---------
// Pins the (massBlock, momentBlock, inertiaBlock) -> (linLin=M, angLin=F,
// angAng=J) argument order directly.
TEST(MassProperties, ArticulatedInertiaBlockMapping) {
    Rng rng(14);
    const SymMat33 mass(rng.uniform(1, 3),
                        rng.uniform(-0.2, 0.2),
                        rng.uniform(1, 3),
                        rng.uniform(-0.2, 0.2),
                        rng.uniform(-0.2, 0.2),
                        rng.uniform(1, 3));
    Mat33 massMoment;
    for (int i = 0; i < 9; ++i) {
        massMoment.elems[static_cast<std::size_t>(i)] = rng.uniform(-1, 1);
    }
    const Vec3 q = rng.vec3();
    const SymMat33 inertia = SymMat33::fromSymmetricPart(crossMat(q).transpose() * crossMat(q));

    const ArticulatedInertia abi(mass, massMoment, inertia);
    EXPECT_TRUE(NearMat33(abi.angAng.full(), inertia.full(), kTight)); // (0,0) = inertia (J)
    EXPECT_TRUE(NearMat33(abi.angLin, massMoment, kTight));            // (0,1) = massMoment (F)
    EXPECT_TRUE(NearMat33(abi.linLin.full(), mass.full(), kTight));    // (1,1) = mass (M)
}
