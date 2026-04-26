#include <gtest/gtest.h>

#include "SimTKcommon/Testing.h"

#include "SimTKcommon.h"
#include "Util.hpp"

using namespace SimTK;

/**
 * DOCUMENTATION OF ORIGINAL TEST INTENT:
 * The authors sought to verify the mathematical integrity of the Mass Properties
 * subsystem. This includes:
 * 1. Algebra of 3D/2D cross products and their matrix representations.
 * 2. Inertia tensor properties: Trace preservation, scalar scaling, and rotation.
 * 3. The Parallel Axis Theorem (shifting) for both UnitInertia and Inertia.
 * 4. High-performance "Half-Cross" optimizations used in Articulated Body Inertia (ABI) updates.
 * 5. SpatialInertia and ArticulatedInertia consistency when shifted via Phi matrices vs. direct API.
 */

// --- CROSS PRODUCT TESTS ---

TEST(SimTKCommon_MassProperties_CrossProduct, Verify3DBehavior) {
    const Vec3 w(1.25, 3.0, -2.5), v(-2.75, 2.125, 5.0);
    const Vec<3, float> wf(1.25f, 3.0f, -2.5f), vf(-2.75f, 2.125f, 5.0f);

    const Vec3 wxv = (w % v);
    const Vec3 vxw = (v % w);
    const Vec<3, float> wxvf = (wf % vf);
    const Vec<3, float> vxwf = (vf % wf);

    EXPECT_TRUE(wxv == Vec3(20.3125, 0.625, 10.90625));
    EXPECT_TRUE(wxvf == (Vec<3, float>(20.3125f, 0.625f, 10.90625f)));
    EXPECT_TRUE(vxw == -wxv);
    EXPECT_TRUE(vxwf == -wxvf);

    // Layout agnosticism (Rows vs Vectors)
    EXPECT_TRUE(((~w) % v) == (~wxv));
    EXPECT_TRUE((w % (~v)) == (~wxv));
    EXPECT_TRUE(((~w) % (~v)) == (~wxv));

    const Mat33 expectedCrossMat(0, 2.5, 3, -2.5, 0, -1.25, -3, 1.25, 0);
    EXPECT_TRUE(crossMat(w) == expectedCrossMat);
    EXPECT_TRUE(crossMatSq(w) == ((~crossMat(w)) * crossMat(w)));
}

TEST(SimTKCommon_MassProperties_CrossProduct, VerifyMatrixInteractions) {
    const Vec3 v = SimTK::Test::randVec3();
    const Mat33 full33 = SimTK::Test::randMat33();
    const SymMat33 sym33 = SimTK::Test::randSymMat33();
    const Mat33 fsym33(sym33);

    EXPECT_TRUE(fsym33 == sym33);

    // v % M relations
    EXPECT_TRUE(AssertSimTKEqual("v % full33",
                                 "col-wise",
                                 (v % full33),
                                 Mat33((v % full33(0)), (v % full33(1)), (v % full33(2)))));
    EXPECT_TRUE(AssertSimTKEqual("v % full33", "crossMat*M", (v % full33), (crossMat(v) * full33)));
    EXPECT_TRUE(AssertSimTKEqual("v % sym33",
                                 "col-wise sym",
                                 (v % sym33),
                                 Mat33((v % fsym33(0)), (v % fsym33(1)), (v % fsym33(2)))));
    EXPECT_TRUE(AssertSimTKEqual("v % sym33", "crossMat*Msym", (v % sym33), (crossMat(v) * fsym33)));

    // M % v relations
    EXPECT_TRUE(AssertSimTKEqual("full33 % v", "M*crossMat", (full33 % v), (full33 * crossMat(v))));
    EXPECT_TRUE(AssertSimTKEqual("full33 % v", "-~(v % ~M)", (full33 % v), (-(~(v % (~full33))))));
    EXPECT_TRUE(AssertSimTKEqual("full33 % v", "~(-v % ~M)", (full33 % v), (~((-v) % (~full33)))));

    EXPECT_TRUE(AssertSimTKEqual("sym33 % v", "Msym*crossMat", (sym33 % v), (fsym33 * crossMat(v))));
    EXPECT_TRUE(AssertSimTKEqual("sym33 % v", "-~(v % sym33)", (sym33 % v), (-(~(v % sym33)))));
    EXPECT_TRUE(AssertSimTKEqual("sym33 % v", "~(-v % sym33)", (sym33 % v), (~((-v) % sym33))));

    EXPECT_TRUE(AssertSimTKEqual("det(sym33)", "det(fsym33)", det(sym33), det(fsym33)));
}

TEST(SimTKCommon_MassProperties_CrossProduct, Verify2DBehavior) {
    const Vec2 w2(1.25, 3.0), v2(-2.75, 2.125);
    const Vec<2, float> w2f(1.25f, 3.0f), v2f(-2.75f, 2.125f);

    const Real wxv2 = (w2 % v2);
    const float wxv2f = (w2f % v2f);

    EXPECT_EQ(wxv2, (Vec3(w2[0], w2[1], 0) % Vec3(v2[0], v2[1], 0))[2]);
    EXPECT_EQ(wxv2f, (Vec<3, float>(w2f[0], w2f[1], 0.0f) % Vec<3, float>(v2f[0], v2f[1], 0.0f))[2]);
    EXPECT_EQ((v2 % w2), -wxv2);

    EXPECT_TRUE(((~w2) % v2) == wxv2);
    EXPECT_TRUE((w2 % (~v2)) == wxv2);
    EXPECT_TRUE(((~w2) % (~v2)) == wxv2);

    EXPECT_TRUE(crossMat(w2) == Row2(-w2[1], w2[0]));
    EXPECT_EQ((crossMat(w2) * v2), (w2 % v2));
}

// --- INERTIA TESTS ---

TEST(SimTKCommon_MassProperties_Inertia, VerifyBasicProperties) {
    const Real mass = std::abs(SimTK::Test::randReal());
    const UnitInertia_<Real> G(Vec3(1, 2, 2.5), Vec3(0.1, 0.2, 0.3));
    const Real Gtrace = 5.5; // 1 + 2 + 2.5

    EXPECT_EQ(G.trace(), Gtrace);

    const Inertia_<Real> I = (mass * G);
    const SymMat33 sI = I.asSymMat33();
    const Mat33 mI = I.toMat33();

    EXPECT_TRUE(
        AssertSimTKEqual("sI", "mass * manual SymMat", sI, (mass * SymMat33(1, 0.1, 2, 0.2, 0.3, 2.5))));
    EXPECT_TRUE(mI.isExactlySymmetric());
    EXPECT_TRUE(mI == Mat33(sI));
    EXPECT_TRUE(sI == SymMat33(mI));
    EXPECT_TRUE(AssertSimTKEqual("I.trace()", "mass*Gtrace", I.trace(), (mass * Gtrace)));
}

TEST(SimTKCommon_MassProperties_Inertia, VerifyScalingAndMultiplication) {
    const Real mass = std::abs(SimTK::Test::randReal());
    const UnitInertia_<Real> G(Vec3(1, 2, 2.5), Vec3(0.1, 0.2, 0.3));
    const Inertia_<Real> I = (mass * G);
    const Real s = std::abs(SimTK::Test::randReal()) + 1.0;

    // Inertia * scalar
    EXPECT_TRUE(AssertSimTKEqual("(I*s).toMat33()", "I.toMat33()*s", (I * s).toMat33(), (I.toMat33() * s)));
    EXPECT_TRUE(AssertSimTKEqual("(s*I).toMat33()", "I.toMat33()*s", (s * I).toMat33(), (I.toMat33() * s)));
    EXPECT_TRUE(AssertSimTKEqual("(G*s).toMat33()", "G.toMat33()*s", (G * s).toMat33(), (G.toMat33() * s)));
    EXPECT_TRUE(AssertSimTKEqual("(s*G).toMat33()", "G.toMat33()*s", (s * G).toMat33(), (G.toMat33() * s)));

    EXPECT_TRUE(AssertSimTKEqual("(I*s).asSymMat33()",
                                 "I.asSymMat33()*s",
                                 (I * s).asSymMat33(),
                                 (I.asSymMat33() * s)));
    EXPECT_TRUE(AssertSimTKEqual("(s*I).asSymMat33()",
                                 "I.asSymMat33()*s",
                                 (s * I).asSymMat33(),
                                 (I.asSymMat33() * s)));

    // Inertia * vec
    const Vec3 w = SimTK::Test::randVec3();
    EXPECT_TRUE(AssertSimTKEqual("I*w", "I.toMat33()*w", (I * w), (I.toMat33() * w)));
    EXPECT_TRUE(AssertSimTKEqual("G*w", "G.toMat33()*w", (G * w), (G.toMat33() * w)));
}

TEST(SimTKCommon_MassProperties_Inertia, VerifyRotationAndShifting) {
    const Real mass = std::abs(SimTK::Test::randReal());
    const UnitInertia_<Real> G(Vec3(1, 2, 2.5), Vec3(0.1, 0.2, 0.3));
    const Inertia_<Real> I = (mass * G);
    const Rotation R = SimTK::Test::randRotation();

    const Inertia_<Real> mIR((~R * I.toMat33() * R));
    EXPECT_TRUE(AssertSimTKEqual("I.reexpress(R)", "mIR", I.reexpress(R).asSymMat33(), mIR.asSymMat33()));

    Inertia_<Real> J = I;
    J.reexpressInPlace(R);
    EXPECT_TRUE(AssertSimTKEqual("J.reexpressInPlace", "mIR", J.asSymMat33(), mIR.asSymMat33()));

    // Parallel Axis Shifting
    const Vec3 pLoc = SimTK::Test::randVec3();
    const SymMat33 psG = crossMatSq(pLoc);
    const UnitInertia_<Real> pG(psG);
    const Inertia_<Real> pI((mass * psG));

    EXPECT_TRUE(AssertSimTKEqual("G.shiftFromCentroid", "G+pG", G.shiftFromCentroid(pLoc), (G + pG)));
    EXPECT_TRUE(
        AssertSimTKEqual("I.shiftFromMassCenter", "I+pI", I.shiftFromMassCenter(pLoc, mass), (I + pI)));

    UnitInertia_<Real> Gshft(G);
    Gshft.shiftFromCentroidInPlace(pLoc);
    Inertia_<Real> Ishft(I);
    Ishft.shiftFromMassCenterInPlace(pLoc, mass);
    EXPECT_TRUE(AssertSimTKEqual("Gshft", "G+pG", Gshft, (G + pG)));
    EXPECT_TRUE(AssertSimTKEqual("Ishft", "I+pI", Ishft, (I + pI)));

    EXPECT_TRUE(AssertSimTKEqual("Gshft.shiftToCentroid", "G", Gshft.shiftToCentroid(pLoc), G));
    EXPECT_TRUE(AssertSimTKEqual("Ishft.shiftToMassCenter", "I", Ishft.shiftToMassCenter(pLoc, mass), I));
}

// --- OPTIMIZED KERNEL TESTS ---

TEST(SimTKCommon_MassProperties_HalfCross, VerifyKernels) {
    const Vec3 v = SimTK::Test::randVec3();
    const Mat33 F = SimTK::Test::randMat33();
    const Mat33 G = SimTK::Test::randMat33();
    const Mat33 vxF = (v % F);
    const Mat33 Gxv = (G % v);

    // Helpers to mimic the static inline templates in the original file
    auto halfCrossVF = [](const Vec3& v_in, const Mat33& F_in) -> SymMat33 {
        return SymMat33(v_in[1] * F_in(2, 0) - v_in[2] * F_in(1, 0),
                        v_in[2] * F_in(0, 0) - v_in[0] * F_in(2, 0),
                        v_in[2] * F_in(0, 1) - v_in[0] * F_in(2, 1),
                        v_in[0] * F_in(1, 0) - v_in[1] * F_in(0, 0),
                        v_in[0] * F_in(1, 1) - v_in[1] * F_in(0, 1),
                        v_in[0] * F_in(1, 2) - v_in[1] * F_in(0, 2));
    };

    auto halfCrossGV = [](const Mat33& G_in, const Vec3& v_in) -> SymMat33 {
        return SymMat33(v_in[2] * G_in(0, 1) - v_in[1] * G_in(0, 2),
                        v_in[2] * G_in(1, 1) - v_in[1] * G_in(1, 2),
                        v_in[0] * G_in(1, 2) - v_in[2] * G_in(1, 0),
                        v_in[2] * G_in(2, 1) - v_in[1] * G_in(2, 2),
                        v_in[0] * G_in(2, 2) - v_in[2] * G_in(2, 0),
                        v_in[1] * G_in(2, 0) - v_in[0] * G_in(2, 1));
    };

    auto halfCrossDiffVFG = [](const Vec3& v_in, const Mat33& F_in, const Mat33& G_in) -> SymMat33 {
        return SymMat33(v_in[1] * (F_in(2, 0) + G_in(0, 2)) - v_in[2] * (F_in(1, 0) + G_in(0, 1)),
                        v_in[2] * (F_in(0, 0) - G_in(1, 1)) - v_in[0] * F_in(2, 0) + v_in[1] * G_in(1, 2),
                        v_in[2] * (F_in(0, 1) + G_in(1, 0)) - v_in[0] * (F_in(2, 1) + G_in(1, 2)),
                        v_in[0] * F_in(1, 0) - v_in[2] * G_in(2, 1) - v_in[1] * (F_in(0, 0) - G_in(2, 2)),
                        v_in[0] * (F_in(1, 1) - G_in(2, 2)) - v_in[1] * F_in(0, 1) + v_in[2] * G_in(2, 0),
                        v_in[0] * (F_in(1, 2) + G_in(2, 1)) - v_in[1] * (F_in(0, 2) + G_in(2, 0)));
    };

    const SymMat33 hvxF(vxF(0, 0), vxF(1, 0), vxF(1, 1), vxF(2, 0), vxF(2, 1), vxF(2, 2));
    EXPECT_TRUE(AssertSimTKEqual("halfCross(v,F)", "hvxF", halfCrossVF(v, F), hvxF));

    const SymMat33 hGxv(Gxv(0, 0), Gxv(1, 0), Gxv(1, 1), Gxv(2, 0), Gxv(2, 1), Gxv(2, 2));
    EXPECT_TRUE(AssertSimTKEqual("halfCross(G,v)", "hGxv", halfCrossGV(G, v), hGxv));

    const SymMat33 hdiff = (hvxF - hGxv);
    EXPECT_TRUE(AssertSimTKEqual("halfCrossDiff", "hdiff", halfCrossDiffVFG(v, F, G), hdiff));
}

// --- SPATIAL / ARTICULATED TESTS ---

TEST(SimTKCommon_MassProperties_SpatialInertia, VerifyShifting) {
    const Real mass = 1.125;
    const Vec3 com(0.1, 0.2, 0.25);
    const UnitInertia gyration(1.8, 1.9, 2.1, 0.01, 0.03, 0.02);
    SpatialInertia si(mass, com, gyration);

    const SpatialMat msi = si.toSpatialMat();
    EXPECT_TRUE(msi(0, 0) == (mass * gyration.toMat33()));
    EXPECT_TRUE(msi(0, 1) == (mass * crossMat(com)));
    EXPECT_TRUE(msi(1, 0) == (mass * (~crossMat(com))));
    EXPECT_TRUE(msi(1, 1) == (mass * Mat33(1)));

    const Vec3 shiftVec(1, 2, 3);
    const PhiMatrix phi(-shiftVec);
    const SpatialMat msiShiftedByPhi = (phi * msi * (~phi));

    const SpatialInertia shiftSi = si.shift(shiftVec);
    EXPECT_TRUE(
        AssertSimTKEqual("shiftSi.toSpatialMat()", "Phi-shifted", shiftSi.toSpatialMat(), msiShiftedByPhi));
    EXPECT_TRUE(AssertSimTKEqual("Round trip", "Original", shiftSi.shift(-shiftVec).toSpatialMat(), msi));
}

TEST(SimTKCommon_MassProperties_ArticulatedInertia, VerifyMethodConsistency) {
    const ArticulatedInertia abi(SimTK::Test::randSymMat33(),
                                 SimTK::Test::randMat33(),
                                 crossMatSq(SimTK::Test::randVec3()));
    const SpatialMat mabi = abi.toSpatialMat();
    const Vec3 shiftVec(1, 2, 3);

    // Manual Spatial Shift
    const SpatialMat shiftMat(Mat33(1), crossMat(-shiftVec), Mat33(0), Mat33(1));
    const SpatialMat manual = (shiftMat * mabi * (~shiftMat));

    // Phi Shift
    const PhiMatrix phi(-shiftVec);
    const SpatialMat viaPhi = (phi * mabi * (~phi));
    EXPECT_TRUE(AssertSimTKEqual("Phi", "Manual", viaPhi, manual));

    // API Shift (Legacy: ABI shift(v) actually shifts by -v relative to SpatialInertia)
    const ArticulatedInertia shiftAbi = abi.shift(-shiftVec);
    EXPECT_TRUE(AssertSimTKEqual("ABI shift", "Manual", shiftAbi.toSpatialMat(), manual));
}

TEST(SimTKCommon_MassProperties_ArticulatedInertia, LargeScaleConsistency) {
    const int n_shifts = 1000; // Original main() test logic
    const ArticulatedInertia abi_in(SimTK::Test::randSymMat33(),
                                    SimTK::Test::randMat33(),
                                    crossMatSq(SimTK::Test::randVec3()));

    std::vector<Vec3> shifts(n_shifts);
    for (int i = 0; i < n_shifts; ++i) {
        shifts[i] = SimTK::Test::randVec3();
    }

    // 1. Manual result
    SpatialMat mabi_manual = abi_in.toSpatialMat();
    for (const auto& s : shifts) {
        SpatialMat S(Mat33(1));
        S(0, 1) = crossMat(s);
        mabi_manual = (S * mabi_manual * (~S));
    }

    // 2. Phi result
    SpatialMat mabi_phi = abi_in.toSpatialMat();
    for (const auto& s : shifts) {
        mabi_phi = (PhiMatrix(s) * mabi_phi * (~PhiMatrix(s)));
    }

    // 3. Fast result
    ArticulatedInertia abi_fast = abi_in;
    for (const auto& s : shifts) {
        abi_fast = abi_fast.shift(s);
    }

    const Real out1 = mabi_manual(1, 1)(2, 2);
    const Real out2 = mabi_phi(1, 1)(2, 2);
    const Real out3 = abi_fast.getMass()(2, 2);

    EXPECT_TRUE(AssertSimTKEqual("out1", "out2", out1, out2));
    EXPECT_TRUE(AssertSimTKEqual("out2", "out3", out2, out3));
}
