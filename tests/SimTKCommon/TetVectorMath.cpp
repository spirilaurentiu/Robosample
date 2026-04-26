// TestVectorMath.cpp
//
// Google Test suite for element-wise math functions and reduction operations
// applied to SimTK vector/matrix container types.
//
// What the original authors tested
// ─────────────────────────────────
// 1. Unary element-wise math functions (abs, exp, log, sqrt, sin, cos, tan,
//    asin, acos, atan, sinh, cosh, tanh) applied to every SimTK container:
//      • Dynamic heap-allocated:  Vector, RowVector, Matrix
//      • Fixed stack-allocated:   Vec5, Row5, Mat<2,3>
//      • Symmetric matrix:        SymMat<2>
//    Each call must act independently on every stored element.
//
// 2. Reduction / aggregate functions (sum, min, max, mean, sort, median)
//    on the same container types, verifying:
//      • Scalar return for 1-D containers.
//      • Column-wise vector return for 2-D containers.
//      • sort(SymMat<N>) returns a non-symmetric Mat<N,N>.
//
// Input data (shared across all tests)
// ─────────────────────────────────────
//   Vectors / rows  : {-1, 2, -3, 4, -5}
//   Matrix (2×3)    : |  -1   2  -3 |
//                     |   4  -5   6 |
//   SymMat<2> with lower triangle {-1, 2, -3}:
//       (0,0) = -1,  (1,0) = (0,1) = 2,  (1,1) = -3
//
// NaN-producing functions on the above inputs (IEEE-754 real arithmetic)
// ────────────────────────────────────────────────────────────────────────
//   log  : log(x < 0) → NaN     [input elements: -1, -3, -5]
//   sqrt : sqrt(x < 0) → NaN    [input elements: -1, -3, -5]
//   asin : asin(|x| > 1) → NaN  [input elements:  2, -3,  4, -5]
//   acos : acos(|x| > 1) → NaN  [input elements:  2, -3,  4, -5]
//
// The NaN-aware helpers EqOrBothNaN / ExpectVecEq / ExpectMatEq /
// ExpectSymMatEq replicate the original testVector / testMatrix behaviour of
// treating (NaN, NaN) as a passing comparison.
//
// Note: the original source contains a copy-paste comment "Test the asin
// function" inside the acos section; this suite corrects it.

#include <cmath>
#include <gtest/gtest.h>
#include <string>

#include "SimTKcommon.h"
#include "Util.hpp"

using namespace SimTK;

// ─────────────────────────────────────────────────────────────────────────────
// Test-data factories and NaN-aware assertion helpers (internal linkage)
// ─────────────────────────────────────────────────────────────────────────────

namespace {

// ── Test-data factories ───────────────────────────────────────────────────────

[[nodiscard]] auto MakeVector() -> Vector {
    Vector v(5);
    v[0] = -1.0;
    v[1] = 2.0;
    v[2] = -3.0;
    v[3] = 4.0;
    v[4] = -5.0;
    return v;
}

[[nodiscard]] auto MakeMatrix() -> Matrix {
    Matrix m(2, 3);
    m[0][0] = -1.0;
    m[0][1] = 2.0;
    m[0][2] = -3.0;
    m[1][0] = 4.0;
    m[1][1] = -5.0;
    m[1][2] = 6.0;
    return m;
}

// SymMat<2> stores only the lower triangle as {diag0, offDiag10, diag1}.
// With data = {-1, 2, -3}:  (0,0)=-1,  (1,0)=(0,1)=2,  (1,1)=-3.
[[nodiscard]] auto MakeSymMat() -> SymMat<2> {
    const Real data[3] = {-1.0, 2.0, -3.0};
    return SymMat<2>(data);
}

// Fixed-size containers constructed once at namespace scope.
const Vec5 kVec(-1.0, 2.0, -3.0, 4.0, -5.0);
const Row5 kRow(-1.0, 2.0, -3.0, 4.0, -5.0);
const Mat<2, 3> kMat(-1.0, 2.0, -3.0, 4.0, -5.0, 6.0);

// ── NaN-aware assertion helpers ───────────────────────────────────────────────

// Returns success when both values are NaN, or when they are numerically equal
// within the default SimTK tolerance. Returns a descriptive failure otherwise.
[[nodiscard]] auto
EqOrBothNaN(Real actual, Real expected, const std::string& actual_expr, const std::string& expected_expr)
    -> ::testing::AssertionResult {
    if (SimTK::isNaN(expected)) {
        if (SimTK::isNaN(actual)) {
            return ::testing::AssertionSuccess();
        }
        return ::testing::AssertionFailure() << "Value of: " << actual_expr << "\n"
                                             << "  Actual: " << actual << "\n"
                                             << "Expected: NaN (from " << expected_expr << ")";
    }
    return AssertSimTKEqual(actual_expr, expected_expr, actual, expected);
}

// Element-wise comparison for any 1-D SimTK container supporting operator[].
// Uses ASSERT_EQ for the size check (to avoid out-of-bounds access on mismatch)
// and EXPECT_TRUE per element so all element failures are reported.
template <typename TActual, typename TExpected>
auto ExpectVecEq(const TActual& actual, const TExpected& expected, int n) -> void {
    ASSERT_EQ(static_cast<int>(actual.size()), n);
    for (int i = 0; i < n; ++i) {
        EXPECT_TRUE(EqOrBothNaN(actual[i],
                                expected[i],
                                ("actual[" + std::to_string(i) + "]"),
                                ("expected[" + std::to_string(i) + "]")))
            << "  at index i=" << i;
    }
}

// Element-wise comparison for any 2-D SimTK container supporting operator().
template <typename TActual, typename TExpected>
auto ExpectMatEq(const TActual& actual, const TExpected& expected, int rows, int cols) -> void {
    ASSERT_EQ(actual.nrow(), rows);
    ASSERT_EQ(actual.ncol(), cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            EXPECT_TRUE(EqOrBothNaN(actual(i, j),
                                    expected(i, j),
                                    ("actual(" + std::to_string(i) + "," + std::to_string(j) + ")"),
                                    ("expected(" + std::to_string(i) + "," + std::to_string(j) + ")")))
                << "  at (" << i << ", " << j << ")";
        }
    }
}

// Lower-triangle comparison for SymMat<N>; the upper triangle is implied by
// symmetry and need not be checked separately.
template <typename TActual, int N>
auto ExpectSymMatEq(const TActual& actual, const SymMat<N>& expected) -> void {
    ASSERT_EQ(actual.nrow(), N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            EXPECT_TRUE(EqOrBothNaN(actual(i, j),
                                    expected(i, j),
                                    ("actual(" + std::to_string(i) + "," + std::to_string(j) + ")"),
                                    ("expected(" + std::to_string(i) + "," + std::to_string(j) + ")")))
                << "  at (" << i << ", " << j << ")";
        }
    }
}

} // namespace

// ═════════════════════════════════════════════════════════════════════════════
// abs — |x| for every element; no NaN produced on any real input.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Abs, AppliesElementwiseToVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    // abs{-1, 2, -3, 4, -5} → {1, 2, 3, 4, 5}
    const Vec5 expected(1.0, 2.0, 3.0, 4.0, 5.0);

    ExpectVecEq(abs(vector), expected, 5);
    ExpectVecEq(abs(rowvector), expected, 5);
    ExpectVecEq(abs(kVec), expected, 5);
    ExpectVecEq(abs(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Abs, AppliesElementwiseToMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);

    ExpectMatEq(abs(matrix), expected, 2, 3);
    ExpectMatEq(abs(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Abs, AppliesElementwiseToSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    // Lower triangle {-1, 2, -3} → abs → {1, 2, 3}
    const Real data[3] = {1.0, 2.0, 3.0};
    const SymMat<2> expected(data);

    ExpectSymMatEq(abs(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// exp — defined for all real x; no NaN on these inputs.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Exp, AppliesElementwiseToVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::exp(-1.0), std::exp(2.0), std::exp(-3.0), std::exp(4.0), std::exp(-5.0));

    ExpectVecEq(exp(vector), expected, 5);
    ExpectVecEq(exp(rowvector), expected, 5);
    ExpectVecEq(exp(kVec), expected, 5);
    ExpectVecEq(exp(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Exp, AppliesElementwiseToMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::exp(-1.0),
                             std::exp(2.0),
                             std::exp(-3.0),
                             std::exp(4.0),
                             std::exp(-5.0),
                             std::exp(6.0));

    ExpectMatEq(exp(matrix), expected, 2, 3);
    ExpectMatEq(exp(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Exp, AppliesElementwiseToSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    // Lower triangle {-1, 2, -3} → {exp(-1), exp(2), exp(-3)}
    const Real data[3] = {std::exp(-1.0), std::exp(2.0), std::exp(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(exp(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// log — log(x < 0) yields NaN (IEEE-754).
// From {-1, 2, -3, 4, -5}: log(-1)=NaN, log(2)=valid, log(-3)=NaN,
//                           log(4)=valid, log(-5)=NaN.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Log, YieldsNaNForNegativeElementsInVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::log(-1.0), std::log(2.0), std::log(-3.0), std::log(4.0), std::log(-5.0));

    ExpectVecEq(log(vector), expected, 5);
    ExpectVecEq(log(rowvector), expected, 5);
    ExpectVecEq(log(kVec), expected, 5);
    ExpectVecEq(log(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Log, YieldsNaNForNegativeElementsInMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::log(-1.0),
                             std::log(2.0),
                             std::log(-3.0),
                             std::log(4.0),
                             std::log(-5.0),
                             std::log(6.0));

    ExpectMatEq(log(matrix), expected, 2, 3);
    ExpectMatEq(log(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Log, YieldsNaNForNegativeElementsInSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    // Lower triangle {-1, 2, -3}: log(-1)=NaN, log(2)=valid, log(-3)=NaN
    const Real data[3] = {std::log(-1.0), std::log(2.0), std::log(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(log(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// sqrt — sqrt(x < 0) yields NaN (IEEE-754).
// From {-1, 2, -3, 4, -5}: sqrt(-1)=NaN, sqrt(2)=valid, sqrt(-3)=NaN,
//                            sqrt(4)=2,    sqrt(-5)=NaN.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Sqrt, YieldsNaNForNegativeElementsInVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::sqrt(-1.0), std::sqrt(2.0), std::sqrt(-3.0), std::sqrt(4.0), std::sqrt(-5.0));

    ExpectVecEq(sqrt(vector), expected, 5);
    ExpectVecEq(sqrt(rowvector), expected, 5);
    ExpectVecEq(sqrt(kVec), expected, 5);
    ExpectVecEq(sqrt(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Sqrt, YieldsNaNForNegativeElementsInMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::sqrt(-1.0),
                             std::sqrt(2.0),
                             std::sqrt(-3.0),
                             std::sqrt(4.0),
                             std::sqrt(-5.0),
                             std::sqrt(6.0));

    ExpectMatEq(sqrt(matrix), expected, 2, 3);
    ExpectMatEq(sqrt(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Sqrt, YieldsNaNForNegativeElementsInSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    // Lower triangle {-1, 2, -3}: sqrt(-1)=NaN, sqrt(2)=valid, sqrt(-3)=NaN
    const Real data[3] = {std::sqrt(-1.0), std::sqrt(2.0), std::sqrt(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(sqrt(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// sin — defined for all real x; no NaN.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Sin, AppliesElementwiseToVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::sin(-1.0), std::sin(2.0), std::sin(-3.0), std::sin(4.0), std::sin(-5.0));

    ExpectVecEq(sin(vector), expected, 5);
    ExpectVecEq(sin(rowvector), expected, 5);
    ExpectVecEq(sin(kVec), expected, 5);
    ExpectVecEq(sin(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Sin, AppliesElementwiseToMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::sin(-1.0),
                             std::sin(2.0),
                             std::sin(-3.0),
                             std::sin(4.0),
                             std::sin(-5.0),
                             std::sin(6.0));

    ExpectMatEq(sin(matrix), expected, 2, 3);
    ExpectMatEq(sin(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Sin, AppliesElementwiseToSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Real data[3] = {std::sin(-1.0), std::sin(2.0), std::sin(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(sin(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// cos — defined for all real x; no NaN.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Cos, AppliesElementwiseToVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::cos(-1.0), std::cos(2.0), std::cos(-3.0), std::cos(4.0), std::cos(-5.0));

    ExpectVecEq(cos(vector), expected, 5);
    ExpectVecEq(cos(rowvector), expected, 5);
    ExpectVecEq(cos(kVec), expected, 5);
    ExpectVecEq(cos(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Cos, AppliesElementwiseToMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::cos(-1.0),
                             std::cos(2.0),
                             std::cos(-3.0),
                             std::cos(4.0),
                             std::cos(-5.0),
                             std::cos(6.0));

    ExpectMatEq(cos(matrix), expected, 2, 3);
    ExpectMatEq(cos(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Cos, AppliesElementwiseToSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Real data[3] = {std::cos(-1.0), std::cos(2.0), std::cos(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(cos(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// tan — singularities at ±π/2 + nπ; none of {-1,2,-3,4,-5} hit them.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Tan, AppliesElementwiseToVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::tan(-1.0), std::tan(2.0), std::tan(-3.0), std::tan(4.0), std::tan(-5.0));

    ExpectVecEq(tan(vector), expected, 5);
    ExpectVecEq(tan(rowvector), expected, 5);
    ExpectVecEq(tan(kVec), expected, 5);
    ExpectVecEq(tan(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Tan, AppliesElementwiseToMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::tan(-1.0),
                             std::tan(2.0),
                             std::tan(-3.0),
                             std::tan(4.0),
                             std::tan(-5.0),
                             std::tan(6.0));

    ExpectMatEq(tan(matrix), expected, 2, 3);
    ExpectMatEq(tan(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Tan, AppliesElementwiseToSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Real data[3] = {std::tan(-1.0), std::tan(2.0), std::tan(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(tan(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// asin — defined only on [-1, 1]; |x| > 1 → NaN.
// From {-1, 2, -3, 4, -5}: asin(-1) = -π/2 (valid); all others → NaN.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Asin, YieldsNaNForOutOfDomainInputsInVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::asin(-1.0), std::asin(2.0), std::asin(-3.0), std::asin(4.0), std::asin(-5.0));

    ExpectVecEq(asin(vector), expected, 5);
    ExpectVecEq(asin(rowvector), expected, 5);
    ExpectVecEq(asin(kVec), expected, 5);
    ExpectVecEq(asin(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Asin, YieldsNaNForOutOfDomainInputsInMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::asin(-1.0),
                             std::asin(2.0),
                             std::asin(-3.0),
                             std::asin(4.0),
                             std::asin(-5.0),
                             std::asin(6.0));

    ExpectMatEq(asin(matrix), expected, 2, 3);
    ExpectMatEq(asin(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Asin, YieldsNaNForOutOfDomainInputsInSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    // Lower triangle {-1, 2, -3}: asin(-1) = -π/2, asin(2) = NaN, asin(-3) = NaN
    const Real data[3] = {std::asin(-1.0), std::asin(2.0), std::asin(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(asin(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// acos — defined only on [-1, 1]; |x| > 1 → NaN.
// From {-1, 2, -3, 4, -5}: acos(-1) = π (valid); all others → NaN.
// Note: the original source labels this section "Test the asin function" —
//       a copy-paste bug preserved here as a comment, not as code.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Acos, YieldsNaNForOutOfDomainInputsInVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::acos(-1.0), std::acos(2.0), std::acos(-3.0), std::acos(4.0), std::acos(-5.0));

    ExpectVecEq(acos(vector), expected, 5);
    ExpectVecEq(acos(rowvector), expected, 5);
    ExpectVecEq(acos(kVec), expected, 5);
    ExpectVecEq(acos(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Acos, YieldsNaNForOutOfDomainInputsInMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::acos(-1.0),
                             std::acos(2.0),
                             std::acos(-3.0),
                             std::acos(4.0),
                             std::acos(-5.0),
                             std::acos(6.0));

    ExpectMatEq(acos(matrix), expected, 2, 3);
    ExpectMatEq(acos(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Acos, YieldsNaNForOutOfDomainInputsInSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    // Lower triangle {-1, 2, -3}: acos(-1) = π, acos(2) = NaN, acos(-3) = NaN
    const Real data[3] = {std::acos(-1.0), std::acos(2.0), std::acos(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(acos(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// atan — defined for all real x; no NaN.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Atan, AppliesElementwiseToVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::atan(-1.0), std::atan(2.0), std::atan(-3.0), std::atan(4.0), std::atan(-5.0));

    ExpectVecEq(atan(vector), expected, 5);
    ExpectVecEq(atan(rowvector), expected, 5);
    ExpectVecEq(atan(kVec), expected, 5);
    ExpectVecEq(atan(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Atan, AppliesElementwiseToMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::atan(-1.0),
                             std::atan(2.0),
                             std::atan(-3.0),
                             std::atan(4.0),
                             std::atan(-5.0),
                             std::atan(6.0));

    ExpectMatEq(atan(matrix), expected, 2, 3);
    ExpectMatEq(atan(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Atan, AppliesElementwiseToSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Real data[3] = {std::atan(-1.0), std::atan(2.0), std::atan(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(atan(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// sinh — defined for all real x; no NaN.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Sinh, AppliesElementwiseToVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::sinh(-1.0), std::sinh(2.0), std::sinh(-3.0), std::sinh(4.0), std::sinh(-5.0));

    ExpectVecEq(sinh(vector), expected, 5);
    ExpectVecEq(sinh(rowvector), expected, 5);
    ExpectVecEq(sinh(kVec), expected, 5);
    ExpectVecEq(sinh(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Sinh, AppliesElementwiseToMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::sinh(-1.0),
                             std::sinh(2.0),
                             std::sinh(-3.0),
                             std::sinh(4.0),
                             std::sinh(-5.0),
                             std::sinh(6.0));

    ExpectMatEq(sinh(matrix), expected, 2, 3);
    ExpectMatEq(sinh(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Sinh, AppliesElementwiseToSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Real data[3] = {std::sinh(-1.0), std::sinh(2.0), std::sinh(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(sinh(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// cosh — defined for all real x; no NaN.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Cosh, AppliesElementwiseToVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::cosh(-1.0), std::cosh(2.0), std::cosh(-3.0), std::cosh(4.0), std::cosh(-5.0));

    ExpectVecEq(cosh(vector), expected, 5);
    ExpectVecEq(cosh(rowvector), expected, 5);
    ExpectVecEq(cosh(kVec), expected, 5);
    ExpectVecEq(cosh(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Cosh, AppliesElementwiseToMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::cosh(-1.0),
                             std::cosh(2.0),
                             std::cosh(-3.0),
                             std::cosh(4.0),
                             std::cosh(-5.0),
                             std::cosh(6.0));

    ExpectMatEq(cosh(matrix), expected, 2, 3);
    ExpectMatEq(cosh(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Cosh, AppliesElementwiseToSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Real data[3] = {std::cosh(-1.0), std::cosh(2.0), std::cosh(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(cosh(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// tanh — defined for all real x; no NaN.
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Tanh, AppliesElementwiseToVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(std::tanh(-1.0), std::tanh(2.0), std::tanh(-3.0), std::tanh(4.0), std::tanh(-5.0));

    ExpectVecEq(tanh(vector), expected, 5);
    ExpectVecEq(tanh(rowvector), expected, 5);
    ExpectVecEq(tanh(kVec), expected, 5);
    ExpectVecEq(tanh(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Tanh, AppliesElementwiseToMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(std::tanh(-1.0),
                             std::tanh(2.0),
                             std::tanh(-3.0),
                             std::tanh(4.0),
                             std::tanh(-5.0),
                             std::tanh(6.0));

    ExpectMatEq(tanh(matrix), expected, 2, 3);
    ExpectMatEq(tanh(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Tanh, AppliesElementwiseToSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Real data[3] = {std::tanh(-1.0), std::tanh(2.0), std::tanh(-3.0)};
    const SymMat<2> expected(data);

    ExpectSymMatEq(tanh(symmat), expected);
}

// ═════════════════════════════════════════════════════════════════════════════
// sum
// For 1-D containers: scalar sum of all elements.
// For 2-D containers: column-wise sum → vector of length #cols.
//   col 0: -1+4=3,  col 1: 2+(-5)=-3,  col 2: -3+6=3
// For SymMat<2> columns (full matrix used):
//   col 0: -1+2=1,  col 1: 2+(-3)=-1
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Sum, ReturnsScalarSumForVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    // (-1) + 2 + (-3) + 4 + (-5) = -3
    EXPECT_TRUE(AssertSimTKEqual("sum(vector)", "-3.0", sum(vector), Real(-3.0)));
    EXPECT_TRUE(AssertSimTKEqual("sum(rowvector)", "-3.0", sum(rowvector), Real(-3.0)));
    EXPECT_TRUE(AssertSimTKEqual("sum(kVec)", "-3.0", sum(kVec), Real(-3.0)));
    EXPECT_TRUE(AssertSimTKEqual("sum(kRow)", "-3.0", sum(kRow), Real(-3.0)));
}

TEST(SimTKCommon_VectorMath_Sum, ReturnsColumnWiseSumForMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Vec3 expected(3.0, -3.0, 3.0);

    ExpectVecEq(sum(matrix), expected, 3);
    ExpectVecEq(sum(kMat), expected, 3);
}

TEST(SimTKCommon_VectorMath_Sum, ReturnsColumnWiseSumForSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Vec2 expected(1.0, -1.0);

    ExpectVecEq(sum(symmat), expected, 2);
}

// ═════════════════════════════════════════════════════════════════════════════
// min
// 1-D: minimum element; 2-D: column-wise minimum.
//   Scalar:  min{-1,2,-3,4,-5} = -5
//   Matrix cols: min(-1,4)=-1, min(2,-5)=-5, min(-3,6)=-3
//   SymMat cols: min(-1,2)=-1, min(2,-3)=-3
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Min, ReturnsScalarMinForVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    EXPECT_TRUE(AssertSimTKEqual("min(vector)", "-5.0", min(vector), Real(-5.0)));
    EXPECT_TRUE(AssertSimTKEqual("min(rowvector)", "-5.0", min(rowvector), Real(-5.0)));
    EXPECT_TRUE(AssertSimTKEqual("min(kVec)", "-5.0", min(kVec), Real(-5.0)));
    EXPECT_TRUE(AssertSimTKEqual("min(kRow)", "-5.0", min(kRow), Real(-5.0)));
}

TEST(SimTKCommon_VectorMath_Min, ReturnsColumnWiseMinForMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Vec3 expected(-1.0, -5.0, -3.0);

    ExpectVecEq(min(matrix), expected, 3);
    ExpectVecEq(min(kMat), expected, 3);
}

TEST(SimTKCommon_VectorMath_Min, ReturnsColumnWiseMinForSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Vec2 expected(-1.0, -3.0);

    ExpectVecEq(min(symmat), expected, 2);
}

// ═════════════════════════════════════════════════════════════════════════════
// max
// 1-D: maximum element; 2-D: column-wise maximum.
//   Scalar:  max{-1,2,-3,4,-5} = 4
//   Matrix cols: max(-1,4)=4, max(2,-5)=2, max(-3,6)=6
//   SymMat cols: max(-1,2)=2, max(2,-3)=2
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Max, ReturnsScalarMaxForVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    EXPECT_TRUE(AssertSimTKEqual("max(vector)", "4.0", max(vector), Real(4.0)));
    EXPECT_TRUE(AssertSimTKEqual("max(rowvector)", "4.0", max(rowvector), Real(4.0)));
    EXPECT_TRUE(AssertSimTKEqual("max(kVec)", "4.0", max(kVec), Real(4.0)));
    EXPECT_TRUE(AssertSimTKEqual("max(kRow)", "4.0", max(kRow), Real(4.0)));
}

TEST(SimTKCommon_VectorMath_Max, ReturnsColumnWiseMaxForMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Vec3 expected(4.0, 2.0, 6.0);

    ExpectVecEq(max(matrix), expected, 3);
    ExpectVecEq(max(kMat), expected, 3);
}

TEST(SimTKCommon_VectorMath_Max, ReturnsColumnWiseMaxForSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Vec2 expected(2.0, 2.0);

    ExpectVecEq(max(symmat), expected, 2);
}

// ═════════════════════════════════════════════════════════════════════════════
// mean
// 1-D: arithmetic mean of all elements; 2-D: column-wise mean.
//   Scalar:  (-1+2-3+4-5)/5 = -3/5 = -0.6
//   Matrix cols: (-1+4)/2=1.5, (2-5)/2=-1.5, (-3+6)/2=1.5
//   SymMat cols: (-1+2)/2=0.5, (2-3)/2=-0.5
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Mean, ReturnsScalarMeanForVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    EXPECT_TRUE(AssertSimTKEqual("mean(vector)", "-0.6", mean(vector), Real(-0.6)));
    EXPECT_TRUE(AssertSimTKEqual("mean(rowvector)", "-0.6", mean(rowvector), Real(-0.6)));
    EXPECT_TRUE(AssertSimTKEqual("mean(kVec)", "-0.6", mean(kVec), Real(-0.6)));
    EXPECT_TRUE(AssertSimTKEqual("mean(kRow)", "-0.6", mean(kRow), Real(-0.6)));
}

TEST(SimTKCommon_VectorMath_Mean, ReturnsColumnWiseMeanForMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Vec3 expected(1.5, -1.5, 1.5);

    ExpectVecEq(mean(matrix), expected, 3);
    ExpectVecEq(mean(kMat), expected, 3);
}

TEST(SimTKCommon_VectorMath_Mean, ReturnsColumnWiseMeanForSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Vec2 expected(0.5, -0.5);

    ExpectVecEq(mean(symmat), expected, 2);
}

// ═════════════════════════════════════════════════════════════════════════════
// sort
// 1-D: ascending sort of all elements.
// 2-D: column-wise ascending sort.
// sort(SymMat<2>) returns Mat<2,2> (symmetry is not preserved after sorting).
//
// {-1,2,-3,4,-5} sorted → {-5,-3,-1,2,4}
// Matrix col sorts: (-1,4)→(-1,4), (2,-5)→(-5,2), (-3,6)→(-3,6)
//   result rows: (-1,-5,-3) and (4,2,6)
// SymMat<2> full matrix: (-1,2 / 2,-3)
//   col sorts: (-1,2)→(-1,2), (2,-3)→(-3,2)
//   result rows: (-1,-3) and (2,2)
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Sort, SortsElementsAscendingForVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    const Vec5 expected(-5.0, -3.0, -1.0, 2.0, 4.0);

    ExpectVecEq(sort(vector), expected, 5);
    ExpectVecEq(sort(rowvector), expected, 5);
    ExpectVecEq(sort(kVec), expected, 5);
    ExpectVecEq(sort(kRow), expected, 5);
}

TEST(SimTKCommon_VectorMath_Sort, SortsColumnsAscendingForMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Mat<2, 3> expected(-1.0, -5.0, -3.0, 4.0, 2.0, 6.0);

    ExpectMatEq(sort(matrix), expected, 2, 3);
    ExpectMatEq(sort(kMat), expected, 2, 3);
}

TEST(SimTKCommon_VectorMath_Sort, SortsColumnsAscendingForSymMat) {
    // sort(SymMat<N>) returns Mat<N,N>, not SymMat<N>.
    const SymMat<2> symmat = MakeSymMat();
    const Mat<2, 2> expected(-1.0, -3.0, 2.0, 2.0);

    ExpectMatEq(sort(symmat), expected, 2, 2);
}

// ═════════════════════════════════════════════════════════════════════════════
// median
// 1-D odd length : middle element of sorted data.
// 1-D even length: mean of the two middle elements.
// 2-D            : column-wise median.
//
// Odd (n=5): sorted {-5,-3,-1,2,4} → median = -1
// Even (n=6): {6,1,5,2,4,3} sorted {1,2,3,4,5,6} → median = (3+4)/2 = 3.5
// Matrix col medians (n=2 per col): (-1+4)/2=1.5, (2-5)/2=-1.5, (-3+6)/2=1.5
// SymMat col medians: (-1+2)/2=0.5, (2-3)/2=-0.5
// ═════════════════════════════════════════════════════════════════════════════

TEST(SimTKCommon_VectorMath_Median, ReturnsMiddleElementForOddSizedVectorTypes) {
    const Vector vector = MakeVector();
    const RowVector rowvector = ~vector;
    EXPECT_TRUE(AssertSimTKEqual("median(vector)", "-1.0", median(vector), Real(-1.0)));
    EXPECT_TRUE(AssertSimTKEqual("median(rowvector)", "-1.0", median(rowvector), Real(-1.0)));
    EXPECT_TRUE(AssertSimTKEqual("median(kVec)", "-1.0", median(kVec), Real(-1.0)));
    EXPECT_TRUE(AssertSimTKEqual("median(kRow)", "-1.0", median(kRow), Real(-1.0)));
}

TEST(SimTKCommon_VectorMath_Median, ReturnsMeanOfMiddleTwoForEvenSizedVector) {
    // Sorted: {1,2,3,4,5,6} → middle pair (3,4) → mean = 3.5
    const Vec6 even_vec(6.0, 1.0, 5.0, 2.0, 4.0, 3.0);
    EXPECT_TRUE(AssertSimTKEqual("median(even_vec)", "3.5", median(even_vec), Real(3.5)));
}

TEST(SimTKCommon_VectorMath_Median, ReturnsColumnWiseMedianForMatrixTypes) {
    const Matrix matrix = MakeMatrix();
    const Vec3 expected(1.5, -1.5, 1.5);

    ExpectVecEq(median(matrix), expected, 3);
    ExpectVecEq(median(kMat), expected, 3);
}

TEST(SimTKCommon_VectorMath_Median, ReturnsColumnWiseMedianForSymMat) {
    const SymMat<2> symmat = MakeSymMat();
    const Vec2 expected(0.5, -0.5);

    ExpectVecEq(median(symmat), expected, 2);
}
