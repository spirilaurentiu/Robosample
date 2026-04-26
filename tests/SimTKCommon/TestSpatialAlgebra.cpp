#include <cmath>
#include <gtest/gtest.h>

#include "SimTKcommon/Testing.h"

#include "SimTKcommon.h"
#include "Util.hpp"

using namespace SimTK;

TEST(SimTKCommon_SpatialAlgebra_PhiMatrix, ForwardMultiplyMatchesAnalyticFormula) {
    const Vec3 p = SimTK::Test::randVec3();
    const PhiMatrix phi(p);
    const SpatialVec v1(SimTK::Test::randSpatialVec());

    EXPECT_TRUE(AssertSimTKEqual("phi * v1",
                                 "SpatialVec(v1[0] + p % v1[1], v1[1])",
                                 phi * v1,
                                 SpatialVec(v1[0] + (p % v1[1]), v1[1])));
}

TEST(SimTKCommon_SpatialAlgebra_PhiMatrix, TransposeMultiplyMatchesAnalyticFormula) {
    const Vec3 p = SimTK::Test::randVec3();
    const PhiMatrix phi(p);
    const SpatialVec v1(SimTK::Test::randSpatialVec());

    EXPECT_TRUE(AssertSimTKEqual("~phi * v1",
                                 "SpatialVec(v1[0], v1[1] - p % v1[0])",
                                 ~phi * v1,
                                 SpatialVec(v1[0], v1[1] - (p % v1[0]))));
}

TEST(SimTKCommon_SpatialAlgebra_PhiMatrix, ForwardVectorMultiplyMatchesDenseMatrixProduct) {
    const Vec3 p = SimTK::Test::randVec3();
    const PhiMatrix phi(p);
    const SpatialVec v1(SimTK::Test::randSpatialVec());

    EXPECT_TRUE(AssertSimTKEqual("phi * v1", "phi.toSpatialMat() * v1", phi * v1, phi.toSpatialMat() * v1));
}

TEST(SimTKCommon_SpatialAlgebra_PhiMatrix, ForwardMatrixMultiplyMatchesDenseMatrixProduct) {
    const Vec3 p = SimTK::Test::randVec3();
    const PhiMatrix phi(p);
    const SpatialMat m1(SimTK::Test::randSpatialMat());

    EXPECT_TRUE(AssertSimTKEqual("phi * m1", "phi.toSpatialMat() * m1", phi * m1, phi.toSpatialMat() * m1));
}

TEST(SimTKCommon_SpatialAlgebra_PhiMatrix, TransposeVectorMultiplyMatchesDenseMatrixProduct) {
    const Vec3 p = SimTK::Test::randVec3();
    const PhiMatrix phi(p);
    const SpatialVec v1(SimTK::Test::randSpatialVec());

    EXPECT_TRUE(
        AssertSimTKEqual("~phi * v1", "(~phi).toSpatialMat() * v1", ~phi * v1, (~phi).toSpatialMat() * v1));
}

TEST(SimTKCommon_SpatialAlgebra_PhiMatrix, MatrixTimesTransposeMatchesDenseMatrixProduct) {
    const Vec3 p = SimTK::Test::randVec3();
    const PhiMatrix phi(p);
    const SpatialMat m1(SimTK::Test::randSpatialMat());

    EXPECT_TRUE(
        AssertSimTKEqual("m1 * ~phi", "m1 * (~phi).toSpatialMat()", m1 * ~phi, m1 * (~phi).toSpatialMat()));
}

TEST(SimTKCommon_SpatialAlgebra_TypeTraits, BlockDimensionsAreTwo) {
    EXPECT_EQ(SpatialMat::NRows, 2);
    EXPECT_EQ(SpatialMat::NCols, 2);
}

TEST(SimTKCommon_SpatialAlgebra_TypeTraits, ArgDepthReflectsNestingHierarchy) {
    EXPECT_EQ(CNT<Real>::ArgDepth, 1);
    EXPECT_EQ(Mat33::ArgDepth, 2);
    EXPECT_EQ(SpatialVec::ArgDepth, 3);
    EXPECT_EQ(SpatialRow::ArgDepth, 3);
    EXPECT_EQ(SpatialMat::ArgDepth, 3);
    EXPECT_EQ(CNT<SpatialMat>::ArgDepth, SpatialMat::ArgDepth);
}

TEST(SimTKCommon_SpatialAlgebra_ScalarAssignment, VectorAssignedOneHasAllComponentsOne) {
    SpatialVec v;
    v = 1;

    EXPECT_TRUE(AssertSimTKEqual("v",
                                 "SpatialVec(Vec3(1,1,1), Vec3(1,1,1))",
                                 v,
                                 SpatialVec(Vec3(1, 1, 1), Vec3(1, 1, 1))));
}

// TEST(SimTKCommon_SpatialAlgebra_ScalarAssignment, RowAssignedOneHasAllComponentsOne) {
//     SpatialRow r;
//     r = 1;

//     EXPECT_TRUE(AssertSimTKEqual("r",
//                                  "SpatialRow(Vec3(1,1,1), Vec3(1,1,1))",
//                                  r,
//                                  SpatialRow(Vec3(1, 1, 1), Vec3(1, 1, 1))));
// }

TEST(SimTKCommon_SpatialAlgebra_ScalarAssignment, MatrixAssignedOneIsBlockIdentity) {
    SpatialMat m;
    m = 1;

    EXPECT_TRUE(AssertSimTKEqual("m", "SpatialMat(1)", m, SpatialMat(1)));
}

TEST(SimTKCommon_SpatialAlgebra_MemoryLayout, SpatialVecAndVec6ShareContiguousStorage) {
    SpatialVec sv(Vec3(1, 2, -3), Vec3(.1, -.2, .3));
    Vec6& sv6 = Vec6::updAs(&sv[0][0]);

    EXPECT_EQ(&sv[0][0], &sv6[0]);

    for (int i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ((&sv[0][0])[i], sv6[i])
            << "Scalar element [" << i << "] differs between SpatialVec and Vec6 alias.";
    }
}

TEST(SimTKCommon_SpatialAlgebra_MemoryLayout, NegatedAliasAndNegatedVec6ShareSameAddress) {
    SpatialVec sv(Vec3(1, 2, -3), Vec3(.1, -.2, .3));
    SpatialVec::TNeg& nsv = sv.updNegate();
    Vec6::TNeg& nsv6 = Vec6::TNeg::updAs(&nsv[0][0]);

    EXPECT_EQ(reinterpret_cast<const void*>(&nsv[0][0]), reinterpret_cast<const void*>(&nsv6[0]));
}

TEST(SimTKCommon_SpatialAlgebra_Normalization, NormalizedSpatialVecHasUnitNorm) {
    const SpatialVec sv(Vec3(1, 2, -3), Vec3(.1, -.2, .3));
    const SpatialVec normalized = sv.normalize();

    EXPECT_TRUE(AssertSimTKEqual("sv.normalize().norm()", "1.0", normalized.norm(), Real(1)));
}

TEST(SimTKCommon_SpatialAlgebra_Normalization, NormalizedNegatedAliasHasUnitNorm) {
    SpatialVec sv(Vec3(1, 2, -3), Vec3(.1, -.2, .3));
    SpatialVec::TNeg& nsv = sv.updNegate();

    EXPECT_TRUE(AssertSimTKEqual("nsv.normalize().norm()", "1.0", nsv.normalize().norm(), Real(1)));
}

TEST(SimTKCommon_SpatialAlgebra_Arithmetic, AddingVectorToItselfYieldsDoubleScalarMultiply) {
    SpatialVec v;
    v = 1;

    EXPECT_TRUE(AssertSimTKEqual("v + v", "2 * v", v + v, 2 * v));
}

TEST(SimTKCommon_SpatialAlgebra_Arithmetic, LeftAndRightScalarMultiplyAreEqual) {
    SpatialVec v;
    v = 1;

    EXPECT_TRUE(AssertSimTKEqual("v * 3.0", "3.0 * v", v * 3.0, 3.0 * v));
}

TEST(SimTKCommon_SpatialAlgebra_Arithmetic, DivisionByNegatorEqualsNegatedDivision) {
    SpatialVec v;
    v = 1;
    const negator<Real> neg3(3.0);

    EXPECT_TRUE(AssertSimTKEqual("v / negator<Real>(3.0)", "-(v / 3.0)", v / neg3, -(v / 3.0)));
}

TEST(SimTKCommon_SpatialAlgebra_ConformingMultiply, IdentityTimesVectorIsVector) {
    SpatialVec v;
    SpatialMat m;
    v = 1;
    m = 1;

    EXPECT_TRUE(AssertSimTKEqual("m * v", "v", m * v, v));
}

TEST(SimTKCommon_SpatialAlgebra_ConformingMultiply, RowTimesIdentityIsRow) {
    SpatialRow r;
    SpatialMat m;
    r = 1;
    m = 1;

    EXPECT_TRUE(AssertSimTKEqual("r * m", "r", r * m, r));
}

TEST(SimTKCommon_SpatialAlgebra_ConformingMultiply, IdentitySquaredIsIdentity) {
    SpatialMat m;
    m = 1;

    EXPECT_TRUE(AssertSimTKEqual("m * m", "m", m * m, m));
}

TEST(SimTKCommon_SpatialAlgebra_DotProduct, RowTimesColEqualsTransposedColTimesTransposedRow) {
    SpatialVec v;
    SpatialRow r;
    v = 1;
    r = 1;

    EXPECT_TRUE(AssertSimTKEqual("r * v", "~v * ~r", r * v, ~v * ~r));
}

TEST(SimTKCommon_SpatialAlgebra_DotProduct, AllOnesOperandsYieldSix) {
    SpatialVec v;
    SpatialRow r;
    v = 1;
    r = 1;
    EXPECT_TRUE(AssertSimTKEqual("r * v", "6.0", r * v, Real(6)));
}

TEST(SimTKCommon_SpatialAlgebra_OuterProduct, ColTimesRowEqualsTransposedColTimesTransposedRow) {
    SpatialVec v;
    SpatialRow r;
    v = 1;
    r = 1;

    EXPECT_TRUE(AssertSimTKEqual("v * r", "~r * ~v", v * r, ~r * ~v));
}

TEST(SimTKCommon_SpatialAlgebra_OuterProduct, AllOnesOperandsYieldAllOnesMatrix) {
    SpatialVec v;
    SpatialRow r;
    v = 1;
    r = 1;
    const Mat33 ones33(1, 1, 1, 1, 1, 1, 1, 1, 1);
    const SpatialMat expected(ones33, ones33, ones33, ones33);

    EXPECT_TRUE(AssertSimTKEqual("v * r", "all-ones SpatialMat", v * r, expected));
}

TEST(SimTKCommon_SpatialAlgebra_NonConformingMultiply, AllCombinationsCompileAndExecute) {
    SpatialMat m;
    SpatialRow r;
    SpatialVec v;
    m = 1;
    r = 1;
    v = 1;
    const Mat33 m33(1, 2, 3, 4, 5, 6, 7, 8, 9);
    const Mat22 m22(10, 20, 30, 40);

    [[maybe_unused]] const auto mm33 = m * m33;
    [[maybe_unused]] const auto tmm33 = ~m * m33;
    [[maybe_unused]] const auto m33m = m33 * m;
    [[maybe_unused]] const auto rm33 = r * m33;
    [[maybe_unused]] const auto tvm33 = ~v * m33;
    [[maybe_unused]] const auto mm22 = m * m22;

    SUCCEED();
}

TEST(SimTKCommon_SpatialAlgebra_ElementMultiply, SpatialRowTimesVec3YieldsRow2WithBlockDotProducts) {
    SpatialRow r;
    r = 1;

    const Row<2> result = r * Vec3(1, 2, 3);

    EXPECT_TRUE(AssertSimTKEqual("r * Vec3(1,2,3)", "Row<2>(6.0, 6.0)", result, Row<2>(6, 6)));
}

TEST(SimTKCommon_SpatialAlgebra_SymMat, AddingDiagonalMatrixShiftsOnlyDiagonalElements) {
    const SymMat<3> sy3 = SymMat<3>().setFromLower(Mat<3, 3>(2, 99, 99, 3, 4, 99, 5, 6, 7));
    const SymMat<3> sy3d(-5);
    const SymMat<3> result = sy3 + sy3d;

    EXPECT_TRUE(AssertSimTKEqual("result(0,0)", "-3.0", result(0, 0), Real(-3)));
    EXPECT_TRUE(AssertSimTKEqual("result(1,1)", "-1.0", result(1, 1), Real(-1)));
    EXPECT_TRUE(AssertSimTKEqual("result(2,2)", " 2.0", result(2, 2), Real(2)));

    EXPECT_TRUE(AssertSimTKEqual("result(1,0)", "3.0", result(1, 0), Real(3)));
    EXPECT_TRUE(AssertSimTKEqual("result(2,0)", "5.0", result(2, 0), Real(5)));
    EXPECT_TRUE(AssertSimTKEqual("result(2,1)", "6.0", result(2, 1), Real(6)));
}

TEST(SimTKCommon_SpatialAlgebra_MassProperties, ConstructionAndConversionProduceFiniteValues) {
    const MassProperties mprops(23, Vec3(1, 2, 3), UnitInertia::brick(.1, .2, .3));

    const Mat66 m66 = mprops.toMat66();
    [[maybe_unused]]
    const SpatialMat sm = mprops.toSpatialMat(); // verifies the call compiles

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            EXPECT_TRUE(std::isfinite(m66(i, j))) << "toMat66()(" << i << "," << j << ") is not finite.";
        }
    }
}