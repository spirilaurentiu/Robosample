// MatVecTest_gtest.cpp
//
// Google Test conversion of the original SimTK MatVecTest.cpp.
//
// WHAT THE AUTHORS ORIGINALLY TESTED
// ====================================
// The original test suite mixed five logical concerns:
//
//  1. testNegator        – Vec3/Row3 subtraction and negated-scalar scaling
//                          produce exact results (no floating-point round-off).
//  2. testElementwiseOps – elementwiseMultiply/Divide on Vec3, Row3, Mat22,
//                          and SymMat22 yield expected per-element results.
//  3. testSums           – rowSum/colSum/sum on real and complex Mat and SymMat;
//                          for Hermitian SymMat<Complex>, rowSum and colSum are
//                          conjugates rather than equal.
//  4. testMiscellaneous  – a large print-only function that demonstrated (but
//                          never asserted) at least ten distinct behaviours:
//                            a. imaginary-step differentiation recovers cos(x)
//                            b. dot()  is layout-agnostic (Vec vs Row)
//                            c. outer() is layout-agnostic
//                            d. cross() is layout-agnostic
//                            e. crossMat(vec) == crossMat(row) in 2-D and 3-D
//                            f. crossMat(v) * w == v % w
//                            g. negator<complex> aliased view sums to zero
//                            h. writing through a negator<conjugate> view
//                               updates the underlying storage correctly
//                            i. Vec::drop1 / append1 / insert1 reshape correctly
//                            j. abs() on scalar-element and nested-element Vec
//                            k. getSubVec / getSubRow / updSubMat address the
//                               correct elements
//                            l. diagonal of a rectangular matrix has length
//                               min(nrow, ncol)
//  5. testMatInverse     – a 20×20 random Matrix times its inverse is identity.
//
// CONVERSION NOTES
// =================
// * Every SimTK_TEST / SimTK_TEST_EQ assertion is preserved and converted to a
//   EXPECT_* counterpart; print-only statements have been given real assertions.
// * The naming convention is:
//     TEST(SimTKCommon_MatVec_<Suite>, <ExpectedBehaviour>)
// * Logically independent concerns are each their own TEST().
// * Modern C++17 features (structured bindings, inline variables, if-constexpr,
//   etc.) are used where they improve clarity.
// * The AssertSimTKEqual<> helper from Util.hpp is exposed via EXPECT_SIMTK_EQ.

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <type_traits>

#include "SimTKcommon/SmallMatrix.h"
#include "SimTKcommon/Testing.h"

#include "Util.hpp"

using namespace SimTK;

// =========================================================================
// Suite: Negator
// Authors verified that subtraction of a vector from itself is exactly zero
// and that negating before subtraction equals doubling the negative.
// =========================================================================

TEST(SimTKCommon_MatVec_Negator, Vec3SubtractionWithSelfIsZero) {
    const Vec3 v(1, 2, 3);
    EXPECT_EQ(v - v, Vec3(0));
}

TEST(SimTKCommon_MatVec_Negator, Vec3NegatedMinusSelfEqualsDoubleNegative) {
    const Vec3 v(1, 2, 3);
    EXPECT_EQ(-v - v, -2 * v);
}

TEST(SimTKCommon_MatVec_Negator, Row3SubtractionWithSelfIsZero) {
    const Row3 r(1, 2, 3);
    EXPECT_EQ(r - r, Row3(0));
}

TEST(SimTKCommon_MatVec_Negator, Row3NegatedMinusSelfEqualsDoubleNegative) {
    const Row3 r(1, 2, 3);
    EXPECT_EQ(-r - r, -2 * r);
}

// =========================================================================
// Suite: ElementwiseOps
// Authors verified that elementwiseMultiply and elementwiseDivide produce
// correct per-element results for Vec3, Row3, Mat22, and SymMat22.
// =========================================================================

TEST(SimTKCommon_MatVec_ElementwiseOps, Vec3MultiplyYieldsPerElementProduct) {
    const Vec3 v(1, 2, 3), w(7, 9, 2);
    EXPECT_TRUE(AssertSimTKEqual("v.elementwiseMultiply(w)",
                                 "Vec3(7, 18, 6)",
                                 v.elementwiseMultiply(w),
                                 Vec3(7, 18, 6)));
}

TEST(SimTKCommon_MatVec_ElementwiseOps, Vec3DivideYieldsPerElementQuotient) {
    const Vec3 v(1, 2, 3), w(7, 9, 2);
    EXPECT_TRUE(AssertSimTKEqual("v.elementwiseDivide(w)",
                                 "Vec3(Real(1) / 7, Real(2) / 9, Real(3) / 2)",
                                 v.elementwiseDivide(w),
                                 Vec3(Real(1) / 7, Real(2) / 9, Real(3) / 2)));
}

TEST(SimTKCommon_MatVec_ElementwiseOps, Row3MultiplyYieldsPerElementProduct) {
    const Row3 r(1, 2, 3), s(5, 4, 10);
    EXPECT_TRUE(AssertSimTKEqual("r.elementwiseMultiply(s)",
                                 "Row3(5, 8, 30)",
                                 r.elementwiseMultiply(s),
                                 Row3(5, 8, 30)));
}

TEST(SimTKCommon_MatVec_ElementwiseOps, Row3DivideYieldsPerElementQuotient) {
    const Row3 r(1, 2, 3), s(5, 4, 10);
    EXPECT_TRUE(AssertSimTKEqual("r.elementwiseDivide(s)",
                                 "Row3(Real(1) / 5, Real(2) / 4, Real(3) / 10)",
                                 r.elementwiseDivide(s),
                                 Row3(Real(1) / 5, Real(2) / 4, Real(3) / 10)));
}

TEST(SimTKCommon_MatVec_ElementwiseOps, Mat22MultiplyYieldsPerElementProduct) {
    const Mat22 m(1, 2, 3, 4), n(7, 9, 2, 3);
    EXPECT_TRUE(AssertSimTKEqual("m.elementwiseMultiply(n)",
                                 "Mat22(7, 18, 6, 12)",
                                 m.elementwiseMultiply(n),
                                 Mat22(7, 18, 6, 12)));
}

TEST(SimTKCommon_MatVec_ElementwiseOps, Mat22DivideYieldsPerElementQuotient) {
    const Mat22 m(1, 2, 3, 4), n(7, 9, 2, 3);
    EXPECT_TRUE(AssertSimTKEqual("m.elementwiseDivide(n)",
                                 "Mat22(Real(1) / 7, Real(2) / 9, Real(3) / 2, Real(4) / 3)",
                                 m.elementwiseDivide(n),
                                 Mat22(Real(1) / 7, Real(2) / 9, Real(3) / 2, Real(4) / 3)));
}

TEST(SimTKCommon_MatVec_ElementwiseOps, SymMat22MultiplyYieldsPerElementProduct) {
    const SymMat22 y(1, 2, 3), z(5, 4, 10);
    EXPECT_TRUE(AssertSimTKEqual("y.elementwiseMultiply(z)",
                                 "SymMat22(5, 8, 30)",
                                 y.elementwiseMultiply(z),
                                 SymMat22(5, 8, 30)));
}

TEST(SimTKCommon_MatVec_ElementwiseOps, SymMat22DivideYieldsPerElementQuotient) {
    const SymMat22 y(1, 2, 3), z(5, 4, 10);
    EXPECT_TRUE(AssertSimTKEqual("y.elementwiseDivide(z)",
                                 "SymMat22(Real(1) / 5, Real(2) / 4, Real(3) / 10)",
                                 y.elementwiseDivide(z),
                                 SymMat22(Real(1) / 5, Real(2) / 4, Real(3) / 10)));
}

// =========================================================================
// Suite: Sums
// Authors verified rowSum, colSum, and sum() for real and complex matrices.
// Key insight: for Hermitian SymMat<Complex>, rowSum ≠ colSum (they are
// conjugates), because the stored lower triangle is the conjugate of the
// implicit upper triangle.
// =========================================================================

TEST(SimTKCommon_MatVec_Sums, RealMat22ColSumSumsEachColumn) {
    const Mat22 m(1, 2, 3, 4);
    EXPECT_TRUE(AssertSimTKEqual("m.colSum()", "Row2(4, 6)", m.colSum(), Row2(4, 6)));
}

TEST(SimTKCommon_MatVec_Sums, RealMat22RowSumSumsEachRow) {
    const Mat22 m(1, 2, 3, 4);
    EXPECT_TRUE(AssertSimTKEqual("m.rowSum()", "Vec2(3, 7)", m.rowSum(), Vec2(3, 7)));
}

TEST(SimTKCommon_MatVec_Sums, RealMat22SumAliasesColSum) {
    const Mat22 m(1, 2, 3, 4);
    EXPECT_EQ(m.sum(), m.colSum());
}

TEST(SimTKCommon_MatVec_Sums, RealSymMat22RowAndColSumAreEqualForSymmetricReal) {
    // For real symmetric matrices the Hermitian transpose equals the ordinary
    // transpose, so row sums and column sums are identical.
    const SymMat22 y(3, 4, 5);
    EXPECT_TRUE(AssertSimTKEqual("y.colSum()", "Row2(7, 9)", y.colSum(), Row2(7, 9)));
    EXPECT_TRUE(AssertSimTKEqual("y.rowSum()", "Vec2(7, 9)", y.rowSum(), Vec2(7, 9)));
    EXPECT_EQ(y.sum(), y.colSum());
}

TEST(SimTKCommon_MatVec_Sums, FullMatExpansionOfSymMat22PreservesRowAndColSums) {
    const SymMat22 y(3, 4, 5);
    const Mat22 full(y);
    EXPECT_TRUE(AssertSimTKEqual("full.rowSum()", "Vec2(7, 9)", full.rowSum(), Vec2(7, 9)));
    EXPECT_TRUE(AssertSimTKEqual("full.colSum()", "Row2(7, 9)", full.colSum(), Row2(7, 9)));
}

TEST(SimTKCommon_MatVec_Sums, ComplexMat22ColSumSumsEachColumn) {
    using CRow2 = Row<2, Complex>;
    const Mat<2, 2, Complex> mc(1 + 2 * I, 3 + 4 * I, 5 + 6 * I, 7 + 8 * I);
    EXPECT_TRUE(AssertSimTKEqual("mc.colSum()",
                                 "Row<2, Complex>(6 + 8 * I, 10 + 12 * I)",
                                 mc.colSum(),
                                 CRow2(6 + 8 * I, 10 + 12 * I)));
}

TEST(SimTKCommon_MatVec_Sums, ComplexMat22RowSumSumsEachRow) {
    using CVec2 = Vec<2, Complex>;
    const Mat<2, 2, Complex> mc(1 + 2 * I, 3 + 4 * I, 5 + 6 * I, 7 + 8 * I);
    EXPECT_TRUE(AssertSimTKEqual("mc.rowSum()",
                                 "Vec<2, Complex>(4 + 6 * I, 12 + 14 * I)",
                                 mc.rowSum(),
                                 CVec2(4 + 6 * I, 12 + 14 * I)));
}

TEST(SimTKCommon_MatVec_Sums, ComplexMat22SumAliasesColSum) {
    const Mat<2, 2, Complex> mc(1 + 2 * I, 3 + 4 * I, 5 + 6 * I, 7 + 8 * I);
    EXPECT_EQ(mc.sum(), mc.colSum());
}

TEST(SimTKCommon_MatVec_Sums, HermitianSymMatComplexRowAndColSumAreConjugates) {
    // For a Hermitian complex SymMat, (colSum)[j] == conj((rowSum)[j]).
    using CRow2 = Row<2, Complex>;
    using CVec2 = Vec<2, Complex>;
    const SymMat<2, Complex> yc(1, 3 - 6 * I, 4);
    EXPECT_TRUE(AssertSimTKEqual("yc.colSum()",
                                 "Row<2, Complex>(4 - 6 * I, 7 + 6 * I)",
                                 yc.colSum(),
                                 CRow2(4 - 6 * I, 7 + 6 * I)));
    EXPECT_TRUE(AssertSimTKEqual("yc.rowSum()",
                                 "Vec<2, Complex>(4 + 6 * I, 7 - 6 * I)",
                                 yc.rowSum(),
                                 CVec2(4 + 6 * I, 7 - 6 * I)));
    EXPECT_EQ(yc.sum(), yc.colSum());
}

TEST(SimTKCommon_MatVec_Sums, FullMatExpansionOfComplexSymMatPreservesRowAndColSums) {
    const SymMat<2, Complex> yc(1, 3 - 6 * I, 4);
    const Mat<2, 2, Complex> full(yc);
    EXPECT_TRUE(AssertSimTKEqual("full.rowSum()",
                                 "Vec<2, Complex>(4 + 6 * I, 7 - 6 * I)",
                                 full.rowSum(),
                                 Vec<2, Complex>(4 + 6 * I, 7 - 6 * I)));
    EXPECT_TRUE(AssertSimTKEqual("full.colSum()",
                                 "Row<2, Complex>(4 - 6 * I, 7 + 6 * I)",
                                 full.colSum(),
                                 Row<2, Complex>(4 - 6 * I, 7 + 6 * I)));
}

// =========================================================================
// Suite: ComplexDotProduct
// Authors demonstrated (via prints only) that the global dot() function is
// layout-agnostic: dot(v,w) == dot(r,w) == dot(v,s) == dot(r,s) for any
// combination of Vec/Row arguments.  Unlike operator*, dot() does NOT
// conjugate either argument.
// =========================================================================

namespace {
// Shared complex test data for dot / outer / cross / crossMat suites.
inline const Complex kCplxData[] = {{1., 2.}, {3., 4.}, {5., 6.}, {7., 8.}, {9., 10.}, {10., 11.}};
} // anonymous namespace

TEST(SimTKCommon_MatVec_ComplexDotProduct, GlobalFunctionIsLayoutAgnostic) {
    const Vec<3, Complex> v(kCplxData);
    const Vec<3, Complex> w(&kCplxData[3]);
    const Row<3, Complex> r(kCplxData);
    const Row<3, Complex> s(&kCplxData[3]);

    const auto ref = dot(v, w);
    EXPECT_TRUE(AssertSimTKEqual("dot(r, w)", "Complex(38, 44)", dot(r, w), ref));
    EXPECT_TRUE(AssertSimTKEqual("dot(v, s)", "Complex(38, 44)", dot(v, s), ref));
    EXPECT_TRUE(AssertSimTKEqual("dot(r, s)", "Complex(38, 44)", dot(r, s), ref));
}

// =========================================================================
// Suite: ComplexOuterProduct
// Authors demonstrated that outer() is layout-agnostic.
// =========================================================================

TEST(SimTKCommon_MatVec_ComplexOuterProduct, GlobalFunctionIsLayoutAgnostic) {
    const Vec<3, Complex> v(kCplxData);
    const Vec<3, Complex> w(&kCplxData[3]);
    const Row<3, Complex> r(kCplxData);
    const Row<3, Complex> s(&kCplxData[3]);

    const auto ref = outer(v, w);
    EXPECT_TRUE(AssertSimTKEqual("outer(r, w)", "Mat<3, 3, Complex>(...)", outer(r, w), ref));
    EXPECT_TRUE(AssertSimTKEqual("outer(v, s)", "Mat<3, 3, Complex>(...)", outer(v, s), ref));
    EXPECT_TRUE(AssertSimTKEqual("outer(r, s)", "Mat<3, 3, Complex>(...)", outer(r, s), ref));
}

// =========================================================================
// Suite: ComplexCrossProduct
// Authors demonstrated that cross() is layout-agnostic.
// =========================================================================

TEST(SimTKCommon_MatVec_ComplexCrossProduct, GlobalFunctionIsLayoutAgnostic) {
    const Vec<3, Complex> v(kCplxData);
    const Vec<3, Complex> w(&kCplxData[3]);
    const Row<3, Complex> r(kCplxData);
    const Row<3, Complex> s(&kCplxData[3]);

    const auto cross_rw = cross(r, w);
    const auto cross_vs = cross(v, s);
    const auto cross_rs = cross(r, s);
    const auto cross_ref = cross(v, w);

    for (std::size_t i = 0; i < cross_ref.size(); ++i) {
        EXPECT_EQ(cross_rw[i], cross_ref[i])
            << "cross(r, w)[" << i << "] = " << cross_rw[i] << " but expected " << cross_ref[i];
        EXPECT_EQ(cross_vs[i], cross_ref[i])
            << "cross(v, s)[" << i << "] = " << cross_vs[i] << " but expected " << cross_ref[i];
        EXPECT_EQ(cross_rs[i], cross_ref[i])
            << "cross(r, s)[" << i << "] = " << cross_rs[i] << " but expected " << cross_ref[i];
    }
}

// =========================================================================
// Suite: CrossMat
// Authors demonstrated that:
//   a. crossMat(vec) produces the same skew-symmetric matrix as crossMat(row).
//   b. crossMat(v) * w == v % w (the cross-product operator).
// Both 2-D and 3-D cases were exercised.
// =========================================================================

TEST(SimTKCommon_MatVec_CrossMat, ThreeDMatrixSameFromVecAndRow) {
    const Vec<3, Complex> v(kCplxData);
    const Row<3, Complex> r(kCplxData);
    EXPECT_TRUE(AssertSimTKEqual("crossMat(v)",
                                 "Mat<3, 3, Complex>(...)",
                                 Mat<3, 3, Complex>(crossMat(v)),
                                 Mat<3, 3, Complex>(crossMat(r))));
}

TEST(SimTKCommon_MatVec_CrossMat, ThreeDMatrixTimesVecMatchesCrossOperator) {
    const Vec<3, Complex> v(kCplxData);
    const Vec<3, Complex> w(&kCplxData[3]);
    EXPECT_TRUE(AssertSimTKEqual("crossMat(v) * w", "Vec<3, Complex>(...)", crossMat(v) * w, v % w));
}

TEST(SimTKCommon_MatVec_CrossMat, TwoDMatrixSameFromVecAndRow) {
    const Vec<2, Complex> v2(kCplxData);
    const Row<2, Complex> r2(kCplxData);

    // crossMat of a 2-D vector returns a Row<2>; both must be equal.
    EXPECT_TRUE(AssertSimTKEqual("crossMat(v2)",
                                 "Row<2, Complex>(...)",
                                 Row<2, Complex>(crossMat(v2)),
                                 Row<2, Complex>(crossMat(r2))));
}

TEST(SimTKCommon_MatVec_CrossMat, TwoDMatrixTimesVecMatchesCrossOperator) {
    const Vec<2, Complex> v2(kCplxData);
    const Vec<2, Complex> w2(&kCplxData[2]);
    EXPECT_TRUE(AssertSimTKEqual("crossMat(v2) * w2", "Vec<2, Complex>(...)", crossMat(v2) * w2, v2 % w2));
}

// =========================================================================
// Suite: VecManipulation
// Authors exercised Vec::drop1, append1, and insert1 to verify that elements
// are removed or inserted at the correct index positions.
// =========================================================================

TEST(SimTKCommon_MatVec_VecManipulation, Drop1RemovesElementAtIndex0) {
    const Vec<3, float> v(39.f, 40.f, 50.f);
    EXPECT_EQ(v.drop1(0), (Vec<2, float>(40.f, 50.f)));
}

TEST(SimTKCommon_MatVec_VecManipulation, Drop1RemovesElementAtIndex1) {
    const Vec<3, float> v(39.f, 40.f, 50.f);
    EXPECT_EQ(v.drop1(1), (Vec<2, float>(39.f, 50.f)));
}

TEST(SimTKCommon_MatVec_VecManipulation, Drop1RemovesElementAtIndex2) {
    const Vec<3, float> v(39.f, 40.f, 50.f);
    EXPECT_EQ(v.drop1(2), (Vec<2, float>(39.f, 40.f)));
}

TEST(SimTKCommon_MatVec_VecManipulation, Append1AddsElementAtEnd) {
    const Vec<3, float> v(39.f, 40.f, 50.f);
    EXPECT_TRUE(AssertSimTKEqual("v.append1(3.3f)",
                                 "Vec<4, float>(39.f, 40.f, 50.f, 3.3f)",
                                 v.append1(3.3f),
                                 Vec<4, float>(39.f, 40.f, 50.f, 3.3f)));
}

TEST(SimTKCommon_MatVec_VecManipulation, Insert1AtIndex0PrependsSingleElement) {
    const Vec<3, float> v(39.f, 40.f, 50.f);
    EXPECT_TRUE(AssertSimTKEqual("v.insert1(0, 23.f)",
                                 "Vec<4, float>(23.f, 39.f, 40.f, 50.f)",
                                 v.insert1(0, 23.f),
                                 Vec<4, float>(23.f, 39.f, 40.f, 50.f)));
}

TEST(SimTKCommon_MatVec_VecManipulation, Insert1AtIndex1InsertsBetweenFirstAndSecond) {
    const Vec<3, float> v(39.f, 40.f, 50.f);
    EXPECT_TRUE(AssertSimTKEqual("v.insert1(1, 23.f)",
                                 "Vec<4, float>(39.f, 23.f, 40.f, 50.f)",
                                 v.insert1(1, 23.f),
                                 Vec<4, float>(39.f, 23.f, 40.f, 50.f)));
}

TEST(SimTKCommon_MatVec_VecManipulation, Insert1AtIndex2InsertsBetweenSecondAndThird) {
    const Vec<3, float> v(39.f, 40.f, 50.f);
    EXPECT_TRUE(AssertSimTKEqual("v.insert1(2, 23.f)",
                                 "Vec<4, float>(39.f, 40.f, 23.f, 50.f)",
                                 v.insert1(2, 23.f),
                                 Vec<4, float>(39.f, 40.f, 23.f, 50.f)));
}

TEST(SimTKCommon_MatVec_VecManipulation, Insert1AtIndex3AppendsElement) {
    const Vec<3, float> v(39.f, 40.f, 50.f);
    EXPECT_TRUE(AssertSimTKEqual("v.insert1(3, 23.f)",
                                 "Vec<4, float>(39.f, 40.f, 50.f, 23.f)",
                                 v.insert1(3, 23.f),
                                 Vec<4, float>(39.f, 40.f, 50.f, 23.f)));
}

TEST(SimTKCommon_MatVec_VecManipulation, Insert1ThenDrop1AtSameIndexIsIdentity) {
    // Inserting at position k then dropping at position k must round-trip.
    const Vec<3, float> v(39.f, 40.f, 50.f);
    EXPECT_TRUE(AssertSimTKEqual("v.insert1(2, 23.f).drop1(2)",
                                 "Vec<3, float>(39.f, 40.f, 50.f)",
                                 v.insert1(2, 23.f).drop1(2),
                                 Vec<3, float>(39.f, 40.f, 50.f)));
}

// =========================================================================
// Suite: AbsoluteValue
// Authors verified that abs() flips negative signs for both flat and nested
// Vec element types.
// =========================================================================

TEST(SimTKCommon_MatVec_AbsoluteValue, Vec3NegativeElementsBecomePositive) {
    const Vec3 v(1., -2., -3.);
    EXPECT_EQ(v.abs(), Vec3(1., 2., 3.));
}

TEST(SimTKCommon_MatVec_AbsoluteValue, NestedVec2OfVec3AllElementsBecomePositive) {
    const Vec3 pos(1., 2., 3.);
    const Vec<2, Vec3> nested(pos, -pos);
    const Vec<2, Vec3> expected(pos, pos);
    EXPECT_EQ(nested.abs(), expected);
}

// =========================================================================
// Suite: Normalize
// Authors printed normalize() results; the fundamental invariant is that the
// normalized vector has unit Euclidean norm.
// =========================================================================

TEST(SimTKCommon_MatVec_Normalize, NormalizedVec3HasUnitNorm) {
    const Vec3 v(1., 2., 3.);
    EXPECT_TRUE(AssertSimTKEqual("v.normalize().norm()", "Real(1)", v.normalize().norm(), Real(1)));
}

TEST(SimTKCommon_MatVec_Normalize, NormalizedComplexVec3HasUnitNorm) {
    const Vec<3, Complex> v(kCplxData);
    EXPECT_TRUE(AssertSimTKEqual("v.normalize().norm()", "Real(1)", v.normalize().norm(), Real(1)));
}

// =========================================================================
// Suite: NegatorViews
// Authors demonstrated that aliased negator<complex> and negator<conjugate>
// views share storage with the original Vec and that writes through those
// views update the underlying data consistently.
// =========================================================================

TEST(SimTKCommon_MatVec_NegatorViews, NegatorViewPlusOriginalIsExactlyZero) {
    // negCv2 is cv2 reinterpreted as a negated view; every element of
    // cv2 + negCv2 must therefore be exactly (0, 0).
    Vec<2, std::complex<float>, 1> cv2(std::complex<float>(1.f, 2.f), std::complex<float>(3.f, 4.f));
    auto& negCv2 = reinterpret_cast<Vec<2, negator<std::complex<float>>, 1>&>(cv2);

    using CVec2f = Vec<2, std::complex<float>>;
    const CVec2f result(cv2 + negCv2);
    EXPECT_TRUE(AssertSimTKEqual("cv2 + negCv2",
                                 "CVec2f(std::complex<float>(0.f), std::complex<float>(0.f))",
                                 result,
                                 CVec2f(std::complex<float>(0.f), std::complex<float>(0.f))));
}

TEST(SimTKCommon_MatVec_NegatorViews, NegatorViewElementEqualsNegativeOfOriginal) {
    Vec<2, std::complex<float>, 1> cv2(std::complex<float>(1.f, 2.f), std::complex<float>(3.f, 4.f));
    auto& negCv2 = reinterpret_cast<Vec<2, negator<std::complex<float>>, 1>&>(cv2);
    // negCv2[0] must appear as the negation of cv2[0].
    EXPECT_TRUE(AssertSimTKEqual("negCv2[0]",
                                 "-cv2[0]",
                                 std::complex<float>(negCv2[0]),
                                 -std::complex<float>(cv2[0])));
}

TEST(SimTKCommon_MatVec_NegatorViews, WritesThroughNegConjViewUpdateStoredValue) {
    // negator<conjugate<float>> stores s such that -(conj(s)) is the visible
    // value.  Assigning (8+9i) through the view should store -(conj(8+9i))
    // = -(8-9i) = (-8+9i) in the underlying cv2 element.
    Vec<2, std::complex<float>, 1> cv2(std::complex<float>(1.f, 2.f), std::complex<float>(3.f, 4.f));
    auto& negConjCv2 = reinterpret_cast<Vec<2, negator<conjugate<float>>, 1>&>(cv2);

    negConjCv2[0] = std::complex<float>(8.f, 9.f);

    EXPECT_TRUE(AssertSimTKEqual("cv2[0]",
                                 "std::complex<float>(-8.f, 9.f)",
                                 std::complex<float>(cv2[0]),
                                 std::complex<float>(-8.f, 9.f)));
}

// =========================================================================
// Suite: SubVecSubMat
// Authors verified that getSubVec, getSubRow, and updSubMat address the
// correct indices.
// =========================================================================

TEST(SimTKCommon_MatVec_SubVecSubMat, GetSubVec2From1ExtractsLastTwoElements) {
    // Vec3(1,-2,-3).getSubVec<2>(1) should yield the sub-range [-2,-3].
    const Vec3 v(1., -2., -3.);
    EXPECT_EQ((v.getSubVec<2>(1)), Vec2(-2., -3.));
}

TEST(SimTKCommon_MatVec_SubVecSubMat, UpdSubMat2x2WritesScalarToCorrectBlock) {
    // Build a 4×3 matrix with −1 on the diagonal, 0 elsewhere.
    Mat<4, 3> m;
    m = Real(1); // sets diagonal; off-diagonal stays 0
    Mat<4, 3> neg = -m;

    // Overwrite the 2×2 block at (row=2, col=1) with −27.
    neg.updSubMat<2, 2>(2, 1) = Real(-27);

    // All four cells in the block must be −27.
    EXPECT_EQ(neg(2, 1), Real(-27));
    EXPECT_EQ(neg(2, 2), Real(-27));
    EXPECT_EQ(neg(3, 1), Real(-27));
    EXPECT_EQ(neg(3, 2), Real(-27));

    // Cells outside the block must be untouched.
    EXPECT_EQ(neg(0, 0), Real(-1));
    EXPECT_EQ(neg(1, 1), Real(-1));
    EXPECT_EQ(neg(0, 1), Real(0));
    EXPECT_EQ(neg(1, 0), Real(0));
}

TEST(SimTKCommon_MatVec_SubVecSubMat, GetSubRowExtractsCorrectColumnsFromRow) {
    // Row 2 of neg4x3 is [0, 0, -1]; getSubRow<2>(1) picks columns 1 and 2.
    Mat<4, 3> m;
    m = Real(1);
    const Mat<4, 3> neg = -m;
    EXPECT_EQ((neg[2].getSubRow<2>(1)), (Row<2>(Real(0), Real(-1))));
}

// =========================================================================
// Suite: RectangularDiagonal
// Authors verified that the diagonal slice of a rectangular matrix has length
// equal to min(nrow, ncol), both for the matrix and its transpose.
// =========================================================================

TEST(SimTKCommon_MatVec_RectangularDiagonal, ThreeByTwoMatrixDiagonalLengthIsTwo) {
    // Mat<3,2> – min(3,2) = 2.
    using M32 = Mat<3, 2, Row3, 1, 2>;
    M32 H;
    EXPECT_EQ(H.diag().nrow(), 2);
}

TEST(SimTKCommon_MatVec_RectangularDiagonal, TransposeOfThreeByTwoAlsoDiagonalLengthTwo) {
    using M32 = Mat<3, 2, Row3, 1, 2>;
    typename M32::TransposeType Ht;
    EXPECT_EQ(Ht.diag().nrow(), 2);
}

// =========================================================================
// Suite: ComplexArithmetic
// Authors demonstrated (print-only) that imag(sin(x + i·h)) / h recovers
// cos(x) as h → 0, which validates correct complex-number arithmetic in the
// SimTK type system.  This is the classic "complex-step derivative" identity.
// =========================================================================

TEST(SimTKCommon_MatVec_ComplexArithmetic, ImaginaryStepDerivativeOfSinRecoversCos) {
    constexpr double x = 0.3;
    constexpr double h = 1e-20;
    const Complex approx_cos = std::sin(Complex(x, h)) / h;
    // imag(sin(x+ih))/h → cos(x) as h→0 (exact to machine precision for h=1e-20).
    EXPECT_NEAR(approx_cos.imag(), std::cos(x), 1e-10);
}

// =========================================================================
// Suite: MatrixInverse
// Authors verified that a 20×20 random real Matrix multiplied by its inverse
// produces the identity matrix within floating-point tolerance.
// =========================================================================

TEST(SimTKCommon_MatVec_MatrixInverse, TwentyByTwentyRandomMatrixTimesInverseIsIdentity) {
    constexpr int N = 20;
    const Matrix m = SimTK::Test::randMatrix(N, N);
    const Matrix mi = m.invert();

    Matrix id(N, N);
    id = Real(1); // sets diagonal to 1, off-diagonal stays 0

    const Matrix product = m * mi;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(product(i, j), expected, 1e-10) << "at (" << i << "," << j << ")";
        }
    }
}

// =========================================================================
// Entry point
// =========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}