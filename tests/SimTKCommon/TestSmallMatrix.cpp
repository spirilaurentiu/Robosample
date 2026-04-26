/**
 * @file TestSmallMatrix.cpp
 *
 * Google Test port of the SimTK small-matrix CTest suite (TestSmallMatrix).
 *
 * A matrix operation of size N can be expected to achieve an accuracy of
 * about N*tol where tol is the expected accuracy of a scalar operation.
 *
 * Test suites
 * -----------
 * SimTKCommon_SmallMatrix_Inverse          mat * inv(mat) == I, sizes 1–10
 * SimTKCommon_SmallMatrix_DotProduct       Vec / Row / Spatial dot products
 * SimTKCommon_SmallMatrix_CrossProduct     cross-product algebra and matrix forms
 * SimTKCommon_SmallMatrix_SymMat           SymMat <-> Mat round-trip and multiply
 * SimTKCommon_SmallMatrix_NumericallyEqual numerical-equality helpers
 * SimTKCommon_SmallMatrix_UnitVec          UnitVec3 / CoordinateAxis algebra
 * SimTKCommon_SmallMatrix_AppendRowCol     structural matrix editing
 */

#include <complex>
#include <gtest/gtest.h>
#include <string>

#include "SimTKcommon/Testing.h"

#include "SimTKcommon.h"
#include "Util.hpp" // Provides AssertSimTKEqual

using namespace SimTK;

// ---------------------------------------------------------------------------
// Size-aware assertion helper (add to Util.hpp if not already present).
//
// Calls SimTK::SimTK::Test::numericallyEqual with an explicit scale factor `n`
// (matching SimTK_TEST_EQ_SIZE semantics) and produces a descriptive
// failure message.
// ---------------------------------------------------------------------------
template <typename T1, typename T2>
[[nodiscard]] auto AssertSimTKEqualSize(const std::string& actual_expr,
                                        const std::string& expected_expr,
                                        const T1& actual,
                                        const T2& expected,
                                        int n) -> ::testing::AssertionResult {
    if (SimTK::Test::numericallyEqual(actual, expected, n)) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << "Value of: " << actual_expr << "\n"
           << "  Actual: " << actual << "\n"
           << "Expected: " << expected_expr << "\n"
           << "  Which is: " << expected << "\n  (tolerance scale factor: " << n << ")";
}

// Convenience macros ---------------------------------------------------------

#define EXPECT_SIMTK_EQ(actual, expected) \
    EXPECT_TRUE(AssertSimTKEqual(#actual, #expected, (actual), (expected)))

#define EXPECT_SIMTK_EQ_SIZE(actual, expected, n) \
    EXPECT_TRUE(AssertSimTKEqualSize(#actual, #expected, (actual), (expected), (n)))

// ============================================================================
// SimTKCommon_SmallMatrix_Inverse
//
// For a random NxN matrix M, verifies that M * M^{-1} ≈ I within a tolerance
// scaled by N (matching the expected floating-point error growth for size-N
// operations).  Exercises three code paths: the plain matrix, its transpose
// (~mat), and its negation (-mat).
// ============================================================================

template <int N>
void ExpectInverseIsIdentity() {
    const Mat<N, N> identity(1);
    Mat<N, N> mat = SimTK::Test::randMat<N, N>();

    // Forward path: mat * inv(mat) ≈ I
    EXPECT_SIMTK_EQ_SIZE(mat * mat.invert(), identity, N);

    // Transposed path: ~mat * inv(~mat) ≈ I
    mat = ~mat;
    EXPECT_SIMTK_EQ_SIZE(mat * mat.invert(), identity, N);

    // Negated path: (-mat) * inv(-mat) ≈ I
    EXPECT_SIMTK_EQ_SIZE((-mat) * (-mat).invert(), identity, N);
}

TEST(SimTKCommon_SmallMatrix_Inverse, Size1RandomMatTimesInverseIsIdentity) {
    ExpectInverseIsIdentity<1>();
}
TEST(SimTKCommon_SmallMatrix_Inverse, Size2RandomMatTimesInverseIsIdentity) {
    ExpectInverseIsIdentity<2>();
}
TEST(SimTKCommon_SmallMatrix_Inverse, Size3RandomMatTimesInverseIsIdentity) {
    ExpectInverseIsIdentity<3>();
}
TEST(SimTKCommon_SmallMatrix_Inverse, Size5RandomMatTimesInverseIsIdentity) {
    ExpectInverseIsIdentity<5>();
}
TEST(SimTKCommon_SmallMatrix_Inverse, Size10RandomMatTimesInverseIsIdentity) {
    ExpectInverseIsIdentity<10>();
}

// ============================================================================
// SimTKCommon_SmallMatrix_DotProduct
//
// Verifies dot-product values for Vec3, Row3, and SpatialVec/SpatialRow.
// Checks that dot(v, −v) == −‖v‖², that mixed Vec/Row forms agree, and
// that the row*col multiplication operator gives the same result as the
// free-function dot().
// ============================================================================

TEST(SimTKCommon_SmallMatrix_DotProduct, VecOpposingSignsGiveNegativeResult) {
    const Vec3 v1(1.0, 2.0, 3.0);
    const Vec3 v2(-1.0, -2.0, -3.0);
    EXPECT_SIMTK_EQ(dot(v1, v2), -14.0);
}

TEST(SimTKCommon_SmallMatrix_DotProduct, RowOpposingSignsGiveNegativeResult) {
    const Row3 r1(0.1, 0.2, 0.3);
    const Row3 r2(-0.1, -0.2, -0.3);
    EXPECT_SIMTK_EQ(dot(r1, r2), -0.14);
}

TEST(SimTKCommon_SmallMatrix_DotProduct, MixedVecRowFormsGiveConsistentResult) {
    const Vec3 v1(1.0, 2.0, 3.0);
    const Vec3 v2(-1.0, -2.0, -3.0);
    const Row3 r1(0.1, 0.2, 0.3);
    const Row3 r2(-0.1, -0.2, -0.3);

    // dot(Vec, Row) and dot(Row, Vec) must agree
    EXPECT_SIMTK_EQ(dot(v1, r2), -1.4);
    EXPECT_SIMTK_EQ(dot(r1, v2), -1.4);

    // The Row * Vec multiplication operator must equal the free function
    EXPECT_SIMTK_EQ(r1 * v2, -1.4);
}

TEST(SimTKCommon_SmallMatrix_DotProduct, SpatialTypesProductIsCorrect) {
    const SpatialVec sv(Vec3(1.0, 2.0, 3.0), Vec3(4.0, 5.0, 6.0));
    const SpatialRow sr(Row3(1.0, 2.0, 3.0), Row3(4.0, 5.0, 6.0));
    EXPECT_SIMTK_EQ(sr * sv, 91.0);
}

// ============================================================================
// SimTKCommon_SmallMatrix_CrossProduct
//
// Verifies:
//   1. The component formula for w × v.
//   2. crossMat(w) * v == w % v; skew-symmetry: crossMat(-w) == ~crossMat(w).
//   3. Columnwise cross product:  v % M  (each column crossed with v).
//   4. Rowwise cross product:     M % w  (each row crossed with w).
//   5. Chain identity:  crossMat(v) * M * crossMat(v)  ==  v % M % v.
//
// Note: wp / vp (perpendicular vectors) and their cross-matrices were
// declared in the original suite but never asserted; they are omitted here.
// ============================================================================

TEST(SimTKCommon_SmallMatrix_CrossProduct, ComponentFormulaMatchesFreeFunction) {
    const Vec3 w = SimTK::Test::randVec3();
    const Vec3 v = SimTK::Test::randVec3();

    // Explicit component formula
    const Vec3 expected((w[1] * v[2]) - (w[2] * v[1]),
                        (w[2] * v[0]) - (w[0] * v[2]),
                        (w[0] * v[1]) - (w[1] * v[0]));
    EXPECT_SIMTK_EQ(w % v, expected);

    // free function cross() and operator% must agree for both Vec and Row
    EXPECT_SIMTK_EQ(w % v, cross(w, v));
    EXPECT_SIMTK_EQ(~w % ~v, -cross(~v, ~w));
}

TEST(SimTKCommon_SmallMatrix_CrossProduct, CrossMatEquivalentToOperatorAndSymmetry) {
    const Vec3 w = SimTK::Test::randVec3();
    const Vec3 v = SimTK::Test::randVec3();
    const Mat33 wx = crossMat(w);

    // crossMat(w) * v must equal w % v
    EXPECT_SIMTK_EQ(wx * v, w % v);

    // crossMat accepts a Row as well as a Vec
    EXPECT_SIMTK_EQ(crossMat(~w), wx);

    // Negating the vector transposes the matrix (skew-symmetry)
    EXPECT_SIMTK_EQ(crossMat(-w), ~wx);

    // crossMatSq(w) * v == −w × (w × v)
    EXPECT_SIMTK_EQ(crossMatSq(w) * v, -(w % (w % v)));
}

TEST(SimTKCommon_SmallMatrix_CrossProduct, ColumnwiseVectorMatrixCrossIsCorrect) {
    const Vec3 v = SimTK::Test::randVec3();
    const Mat34 m = SimTK::Test::randMat<3, 4>();
    const Mat<3, 1> m1 = SimTK::Test::randMat<3, 1>();
    const Mat33 vx = crossMat(v);

    // v % M must be equivalent to applying the cross product column-by-column
    const Mat34 c = v % m;
    EXPECT_SIMTK_EQ(c(0), v % m(0));
    EXPECT_SIMTK_EQ(c(1), v % m(1));
    EXPECT_SIMTK_EQ(c(2), v % m(2));
    EXPECT_SIMTK_EQ(c(3), v % m(3));

    // Equivalent to crossMat(v) * M
    EXPECT_SIMTK_EQ(c, vx * m);

    // A row vector is treated identically to a column vector here
    EXPECT_SIMTK_EQ(c, (~v) % m);

    // Single-column matrix variant
    const Mat<3, 1> c1 = v % m1;
    EXPECT_SIMTK_EQ(c1(0), v % m1(0));
}

TEST(SimTKCommon_SmallMatrix_CrossProduct, RowwiseMatrixVectorCrossIsCorrect) {
    const Vec3 w = SimTK::Test::randVec3();
    const Mat34 m = SimTK::Test::randMat<3, 4>();
    const Mat43 mt = ~m;
    const Mat33 wx = crossMat(w);

    // mt % w must be equivalent to applying the cross product row-by-row
    const Mat43 cr = mt % w;
    EXPECT_SIMTK_EQ(cr[0], mt[0] % w);
    EXPECT_SIMTK_EQ(cr[1], mt[1] % w);
    EXPECT_SIMTK_EQ(cr[2], mt[2] % w);
    EXPECT_SIMTK_EQ(cr[3], mt[3] % w);

    // Equivalent to M * crossMat(w)
    EXPECT_SIMTK_EQ(cr, mt * wx);

    // A row vector is treated identically to a column vector here
    EXPECT_SIMTK_EQ(cr, mt % (~w));
}

TEST(SimTKCommon_SmallMatrix_CrossProduct, ChainedCrossEqualsMatrixProduct) {
    const Vec3 v = SimTK::Test::randVec3();
    const Mat33 m33 = SimTK::Test::randMat33();
    const Mat33 vx = crossMat(v);

    // v × (M × v) via chained % must equal crossMat(v) * M * crossMat(v)
    EXPECT_SIMTK_EQ(vx * m33 * vx, v % m33 % v);
}

// ============================================================================
// SimTKCommon_SmallMatrix_SymMat
//
// Verifies that SymMat<N>:
//   • expands correctly to the equivalent full Mat<N,N> via explicit
//     conversion,
//   • round-trips through setFromSymmetric(),
//   • multiplies a Vec on the right (sm*v) and on the left (~v*sm)
//     identically to the equivalent full matrix.
//
// Sizes 2×2 and 3×3 may have specialised inline operators; 4×4 exercises
// the general code path.
//
// Complex matrices must respect Hermitian (conjugate) symmetry: diagonals
// must be real and corresponding off-diagonal entries are complex-conjugate
// pairs, NOT the same value, even though off-diagonal data is stored only
// once.
// ============================================================================

TEST(SimTKCommon_SmallMatrix_SymMat, TwoByTwoRealConversionRoundTripAndMultiply) {
    const Vec3 a = SimTK::Test::randVec3();
    const Vec<2> v = SimTK::Test::randVec<2>();

    const SymMat<2> sm(a[0], a[1], a[2]);
    const Mat<2, 2> m(a[0], a[1], a[1], a[2]);

    EXPECT_SIMTK_EQ((Mat<2, 2>(sm)), m);
    EXPECT_SIMTK_EQ(sm, SymMat<2>().setFromSymmetric(m));
    EXPECT_SIMTK_EQ(sm * v, m * v);
    EXPECT_SIMTK_EQ(~v * sm, ~v * m);
}

TEST(SimTKCommon_SmallMatrix_SymMat, ThreeByThreeRealConversionRoundTripAndMultiply) {
    const Vec<6> a3 = SimTK::Test::randVec<6>();
    const Vec<3> v3 = SimTK::Test::randVec<3>();

    const SymMat<3> sm3(a3[0], a3[1], a3[2], a3[3], a3[4], a3[5]);
    const Mat<3, 3> m3(a3[0], a3[1], a3[3], a3[1], a3[2], a3[4], a3[3], a3[4], a3[5]);

    EXPECT_SIMTK_EQ((Mat<3, 3>(sm3)), m3);
    EXPECT_SIMTK_EQ(sm3, SymMat<3>().setFromSymmetric(m3));
    EXPECT_SIMTK_EQ(sm3 * v3, m3 * v3);
    EXPECT_SIMTK_EQ(~v3 * sm3, ~v3 * m3);
}

TEST(SimTKCommon_SmallMatrix_SymMat, FourByFourRealConversionRoundTripAndMultiply) {
    // 4×4 exercises the general (non-specialised) code path.
    const Vec<10> a4 = SimTK::Test::randVec<10>();
    const Vec<4> v4 = SimTK::Test::randVec<4>();

    const SymMat<4> sm4(a4[0], a4[1], a4[2], a4[3], a4[4], a4[5], a4[6], a4[7], a4[8], a4[9]);
    const Mat<4, 4> m4(a4[0],
                       a4[1],
                       a4[3],
                       a4[6],
                       a4[1],
                       a4[2],
                       a4[4],
                       a4[7],
                       a4[3],
                       a4[4],
                       a4[5],
                       a4[8],
                       a4[6],
                       a4[7],
                       a4[8],
                       a4[9]);

    EXPECT_SIMTK_EQ((Mat<4, 4>(sm4)), m4);
    EXPECT_SIMTK_EQ(sm4, SymMat<4>().setFromSymmetric(m4));
    EXPECT_SIMTK_EQ(sm4 * v4, m4 * v4);
    EXPECT_SIMTK_EQ(~v4 * sm4, ~v4 * m4);
}

TEST(SimTKCommon_SmallMatrix_SymMat, TwoByTwoComplexHermitianConversionAndMultiply) {
    // Complex is tricky for symmetric (really Hermitian) matrices because
    // the diagonals must be real and the corresponding off-diagonals are
    // complex-conjugate pairs, NOT the same value, even though the
    // off-diagonal data is stored only once.
    const Vec<3, Complex> ac(SimTK::Test::randComplex(),
                             SimTK::Test::randComplex(),
                             SimTK::Test::randComplex());
    const Vec<2, Complex> vc(SimTK::Test::randComplex(), SimTK::Test::randComplex());

    const SymMat<2, Complex> smc(ac[0], ac[1], ac[2]);

    // The constructor must generate a conjugate element for the upper-right
    // entry in the expanded full Mat.
    const Mat<2, 2, Complex> mc(ac[0].real(), std::conj(ac[1]), ac[1], ac[2].real());

    const Mat<2, 2, Complex> sm2mc(smc);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            EXPECT_EQ(sm2mc(i, j), mc(i, j))
                << "sm2mc(" << i << "," << j << ") = " << sm2mc(i, j) << " but expected " << mc(i, j);
        }
    }

    EXPECT_SIMTK_EQ(smc, (SymMat<2, Complex>().setFromSymmetric(mc)));

    // Multiply: smc and ~smc (== smc for a Hermitian matrix) must agree
    // with the equivalent full matrix mc.
    for (std::size_t k = 0; k < 2; ++k) {
        EXPECT_EQ((smc * vc)[k], (mc * vc)[k]) << "(smc*vc)[" << k << "] mismatch";
        EXPECT_EQ((~vc * smc)[k], (~vc * mc)[k]) << "(~vc*smc)[" << k << "] mismatch";
        EXPECT_EQ((~smc * vc)[k], (~mc * vc)[k]) << "(~smc*vc)[" << k << "] mismatch";
        // For a Hermitian matrix ~smc * v must equal smc * v
        EXPECT_EQ((~smc * vc)[k], (smc * vc)[k]) << "(~smc*vc)[" << k << "] != (smc*vc)[" << k << "]";
    }
}

// ============================================================================
// SimTKCommon_SmallMatrix_NumericallyEqual
//
// Verifies isNumericallyEqual() for Mat, SymMat, Vec, and Row, covering:
//   • exact equality and self-equality
//   • perturbations within / beyond the default tolerance
//   • submatrix extraction that excludes a perturbed column
//   • tolerance override (tightening must make a nearly-equal pair fail)
//   • scalar overload
//   • mixed precision (double vs float): the looser (float) tolerance applies
//     unless both operands are double
//   • symmetry detection: isExactlySymmetric / isNumericallySymmetric
//   • complex Hermitian symmetry, including Debug-mode throw on bad input
// ============================================================================

// ── Mat<3,4,float> ──────────────────────────────────────────────────────────

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatMatEqualityWithinDefaultTolerance) {
    Mat<3, 4, float> fm1(1);
    Mat<3, 4, float> fm1e(1);
    // Half-tolerance perturbation: should still compare as equal
    fm1e(1, 1) += 0.5F * fm1.getDefaultTolerance();

    const Mat<3, 4, float> fmident34(1.0F, 0.0F, 0.0F, 0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 0.0F, 1.0F, 0.0F);

    EXPECT_EQ(fm1, fmident34);                 // exact equality via operator==
    EXPECT_TRUE(fm1.isNumericallyEqual(fm1));  // self-equality
    EXPECT_TRUE(fm1.isNumericallyEqual(fm1e)); // within half-tolerance
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatMatInequalityBeyondDefaultTolerance) {
    Mat<3, 4, float> fm1(1);
    Mat<3, 4, float> fm1n(1);
    Mat<3, 4, float> fm1nz(1);
    // Double-tolerance diagonal perturbation: should fail equality
    fm1n(2, 2) += 2.0F * fm1.getDefaultTolerance();
    // Nonzero off-diagonal where identity has zero: should fail equality
    fm1nz(1, 3) = 2.0F * fm1.getDefaultTolerance();

    EXPECT_FALSE(fm1.isNumericallyEqual(fm1n));
    EXPECT_FALSE(fm1.isNumericallyEqual(fm1nz));
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatMatSubmatrixExcludesPerturbedColumn) {
    Mat<3, 4, float> fm1(1);
    Mat<3, 4, float> fm1e(1);
    Mat<3, 4, float> fm1n(1);
    Mat<3, 4, float> fm1nz(1);
    fm1e(1, 1) += 0.5F * fm1.getDefaultTolerance();
    fm1n(2, 2) += 2.0F * fm1.getDefaultTolerance();
    // fm1nz perturbation lives in column 3, outside the 3×3 leading sub-block
    fm1nz(1, 3) = 2.0F * fm1.getDefaultTolerance();

    // The leading 3×3 sub-block of fm1nz is the identity → equal to fm1 and fm1e
    EXPECT_TRUE((fm1.getSubMat<3, 3>(0, 0).isNumericallyEqual(fm1nz.getSubMat<3, 3>(0, 0))));
    EXPECT_TRUE((fm1e.getSubMat<3, 3>(0, 0).isNumericallyEqual(fm1nz.getSubMat<3, 3>(0, 0))));
    // fm1n has its perturbation at (2,2), inside the sub-block → not equal
    EXPECT_FALSE((fm1n.getSubMat<3, 3>(0, 0).isNumericallyEqual(fm1nz.getSubMat<3, 3>(0, 0))));
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatMatTighterToleranceMakesPerturbedUnequal) {
    Mat<3, 4, float> fm1(1);
    Mat<3, 4, float> fm1e(1);
    fm1e(1, 1) += 0.5F * fm1.getDefaultTolerance();

    // fm1e passed at default tolerance; tightening to 0.3× must now reject it
    EXPECT_FALSE(fm1.isNumericallyEqual(fm1e, 0.3 * fm1.getDefaultTolerance()));
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatMatIsNumericallyEqualToScalarOne) {
    Mat<3, 4, float> fm1(1);
    Mat<3, 4, float> fm1e(1);
    Mat<3, 4, float> fm1n(1);
    fm1e(1, 1) += 0.5F * fm1.getDefaultTolerance();
    fm1n(2, 2) += 2.0F * fm1.getDefaultTolerance();

    EXPECT_TRUE(fm1.isNumericallyEqual(1.0F));
    EXPECT_TRUE(fm1e.isNumericallyEqual(1.0F));
    EXPECT_FALSE(fm1n.isNumericallyEqual(1.0F));
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, MixedPrecisionMatUsesLooserTolerance) {
    Mat<3, 4, float> fm1(1);
    Mat<3, 4, float> fm1e(1);
    Mat<3, 4, double> dm1(1);
    fm1e(1, 1) += 0.5F * fm1.getDefaultTolerance();
    const Mat<3, 4, double> dfm1e(fm1e);

    // double vs float: should use the looser (float) tolerance → equal
    EXPECT_TRUE(dm1.isNumericallyEqual(fm1));
    EXPECT_TRUE(dm1.isNumericallyEqual(fm1e));
    // double vs double (promoted from float): uses double tolerance → not equal
    EXPECT_FALSE(dm1.isNumericallyEqual(dfm1e));
    // Explicitly forcing float tolerance restores equality
    EXPECT_TRUE(dm1.isNumericallyEqual(dfm1e, fm1e.getDefaultTolerance()));
}

// ── SymMat<3,float> ─────────────────────────────────────────────────────────

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatSymMatEqualityWithinDefaultTolerance) {
    SymMat<3, float> fs1(1);
    SymMat<3, float> fs1e(1);
    fs1e(1, 1) += 0.5F * fs1.getDefaultTolerance();

    const SymMat<3, float> fsident3(1.0F, 0.0F, 1.0F, 0.0F, 0.0F, 1.0F);

    EXPECT_EQ(fs1, fsident3);
    EXPECT_TRUE(fs1.isNumericallyEqual(fs1));
    EXPECT_TRUE(fs1.isNumericallyEqual(fs1e));
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatSymMatInequalityBeyondDefaultTolerance) {
    SymMat<3, float> fs1(1);
    SymMat<3, float> fs1n(1);
    SymMat<3, float> fs1nz(1);
    fs1n(2, 2) += 2.0F * fs1.getDefaultTolerance();
    // (2,1) is a lower-triangle off-diagonal element
    fs1nz(2, 1) = 2.0F * fs1.getDefaultTolerance();

    EXPECT_FALSE(fs1.isNumericallyEqual(fs1n));
    EXPECT_FALSE(fs1.isNumericallyEqual(fs1nz));
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatSymMatTighterToleranceMakesPerturbedUnequal) {
    SymMat<3, float> fs1(1);
    SymMat<3, float> fs1e(1);
    fs1e(1, 1) += 0.5F * fs1.getDefaultTolerance();

    EXPECT_FALSE(fs1.isNumericallyEqual(fs1e, 0.3 * fs1.getDefaultTolerance()));
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatSymMatIsNumericallyEqualToScalarOne) {
    SymMat<3, float> fs1(1);
    SymMat<3, float> fs1e(1);
    SymMat<3, float> fs1n(1);
    fs1e(1, 1) += 0.5F * fs1.getDefaultTolerance();
    fs1n(2, 2) += 2.0F * fs1.getDefaultTolerance();

    EXPECT_TRUE(fs1.isNumericallyEqual(1.0F));
    EXPECT_TRUE(fs1e.isNumericallyEqual(1.0F));
    EXPECT_FALSE(fs1n.isNumericallyEqual(1.0F));
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, MixedPrecisionSymMatUsesLooserTolerance) {
    SymMat<3, float> fs1(1);
    SymMat<3, float> fs1e(1);
    SymMat<3, double> ds1(1);
    fs1e(1, 1) += 0.5F * fs1.getDefaultTolerance();
    const SymMat<3, double> dfs1e(fs1e);

    EXPECT_TRUE(ds1.isNumericallyEqual(fs1));
    EXPECT_TRUE(ds1.isNumericallyEqual(fs1e));
    EXPECT_FALSE(ds1.isNumericallyEqual(dfs1e));
    EXPECT_TRUE(ds1.isNumericallyEqual(dfs1e, fs1e.getDefaultTolerance()));
}

// ── Vec<3,float> and Row<3,float> ───────────────────────────────────────────

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatVecAndRowEqualityWithinTolerance) {
    Vec<3, float> fv1(1);
    Vec<3, float> fv1e(1);
    Vec<3, float> fv1n(1);
    Row<3, float> fr1(1);
    fv1e[1] += 0.5F * fv1.getDefaultTolerance();
    fv1n[0] += 2.0F * fv1.getDefaultTolerance();

    const Vec<3, float> fone(1.0F, 1.0F, 1.0F);

    EXPECT_EQ(fv1, fone);
    EXPECT_TRUE(fv1.isNumericallyEqual(fv1));
    EXPECT_TRUE(fv1.isNumericallyEqual(fv1e));
    EXPECT_FALSE(fv1.isNumericallyEqual(fv1n));

    EXPECT_EQ(fr1, ~fone);
    EXPECT_TRUE(fr1.isNumericallyEqual(~fv1));
    EXPECT_TRUE(fr1.isNumericallyEqual(~fv1e));
    EXPECT_FALSE(fr1.isNumericallyEqual(~fv1n));
}

// ── Symmetry detection ───────────────────────────────────────────────────────

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, NonSquareMatCannotBeSymmetric) {
    const Mat<2, 7, double> notSquare(0);
    EXPECT_FALSE(notSquare.isExactlySymmetric());
    EXPECT_FALSE(notSquare.isNumericallySymmetric());
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, FloatMatSymmetryDetection) {
    Mat<3, 3, float> f33(1);  // exactly symmetric (identity)
    Mat<3, 3, float> f33e(1); // numerically symmetric (tiny asymmetry)
    Mat<3, 3, float> f33n(1); // asymmetry too large to be considered symmetric
    f33e(1, 2) += 0.5F * f33.getDefaultTolerance();
    f33n(2, 0) += 2.0F * f33.getDefaultTolerance();

    EXPECT_TRUE(f33.isExactlySymmetric());
    EXPECT_FALSE(f33e.isExactlySymmetric());
    EXPECT_FALSE(f33n.isExactlySymmetric());

    EXPECT_TRUE(f33.isNumericallySymmetric());
    EXPECT_TRUE(f33e.isNumericallySymmetric());
    EXPECT_FALSE(f33n.isNumericallySymmetric());
}

TEST(SimTKCommon_SmallMatrix_NumericallyEqual, ComplexMatHermitianSymmetryDetection) {
    using Cplx = std::complex<double>;

    // Positional symmetry: M[i][j] == M[j][i], but NOT Hermitian because
    // the imaginary parts are not negated on the conjugate side.
    const Mat<2, 2, Cplx> mcp(Cplx(1.0, 2.0), Cplx(3.0, 4.0), Cplx(3.0, 4.0), Cplx(5.0, 6.0));

    // Hermitian: M[i][j] == conj(M[j][i]); diagonals are real.
    Mat<2, 2, Cplx> mch(Cplx(1.0, 0.0), Cplx(3.0, -4.0), Cplx(3.0, 4.0), Cplx(5.0, 0.0));

    EXPECT_FALSE(mcp.isExactlySymmetric());
    EXPECT_FALSE(mcp.isNumericallySymmetric());

    EXPECT_TRUE(mch.isExactlySymmetric());
    EXPECT_TRUE(mch.isNumericallySymmetric());

    // setFromSymmetric should succeed on an exactly Hermitian matrix
    SymMat<2, Cplx> symTest;
    EXPECT_NO_THROW(symTest.setFromSymmetric(mch));

    // Perturbation within tolerance: no longer exact, but still numerical
    mch(0, 1) += 0.5 * mch.getDefaultTolerance();
    EXPECT_FALSE(mch.isExactlySymmetric());
    EXPECT_TRUE(mch.isNumericallySymmetric());
    EXPECT_NO_THROW(symTest.setFromSymmetric(mch));

    // Perturbation beyond tolerance: no longer numerically symmetric
    mch(0, 1) += 5.0 * mch.getDefaultTolerance();
    EXPECT_FALSE(mch.isExactlySymmetric());
    EXPECT_FALSE(mch.isNumericallySymmetric());

    // In Debug mode setFromSymmetric must throw when mch is too far off
#ifndef NDEBUG
    EXPECT_ANY_THROW(symTest.setFromSymmetric(mch));
#endif
}

// ============================================================================
// SimTKCommon_SmallMatrix_UnitVec
//
// Verifies UnitVec3 / CoordinateAxis / CoordinateDirection:
//   • normalisation of a non-unit input vector
//   • construction from the six CoordinateAxis enum values
//   • implicit conversion to UnitVec3 inside a function call
//   • orthogonality dot products between canonical axes
//   • crossProductAxis() and crossProductSign() tables
//   • negation operator and direction algebra, including double-negation
// ============================================================================

namespace {
[[nodiscard]] auto IsXAxis(const UnitVec3& test) -> bool {
    return test == UnitVec3(1.0, 0.0, 0.0);
}
[[nodiscard]] auto IsNegZAxis(const UnitVec3& test) -> bool {
    return test == UnitVec3(0.0, 0.0, -1.0);
}
} // namespace

TEST(SimTKCommon_SmallMatrix_UnitVec, NormalisationOfDiagonalVector) {
    EXPECT_SIMTK_EQ(Vec3(UnitVec3(1.0, 1.0, 0.0)), Vec3((Sqrt2 / 2.0), (Sqrt2 / 2.0), 0.0));
}

TEST(SimTKCommon_SmallMatrix_UnitVec, ConstructionFromCoordinateAxisEnum) {
    EXPECT_TRUE(UnitVec3(XAxis) == UnitVec3(1.0, 0.0, 0.0));
    EXPECT_TRUE(UnitVec3(YAxis) == UnitVec3(0.0, 1.0, 0.0));
    EXPECT_TRUE(UnitVec3(ZAxis) == UnitVec3(0.0, 0.0, 1.0));
    EXPECT_TRUE(UnitVec3(NegXAxis) == UnitVec3(-1.0, 0.0, 0.0));
    EXPECT_TRUE(UnitVec3(NegYAxis) == UnitVec3(0.0, -1.0, 0.0));
    EXPECT_TRUE(UnitVec3(NegZAxis) == UnitVec3(0.0, 0.0, -1.0));
}

TEST(SimTKCommon_SmallMatrix_UnitVec, ImplicitConversionToUnitVec3InFunctionCall) {
    // CoordinateAxis and CoordinateDirection must implicitly convert to
    // UnitVec3 when passed to a function expecting UnitVec3.
    EXPECT_TRUE(IsXAxis(XAxis));
    EXPECT_FALSE(IsXAxis(YAxis));
    EXPECT_TRUE(IsNegZAxis(NegZAxis));
    EXPECT_FALSE(IsNegZAxis(NegYAxis));
    EXPECT_FALSE(IsNegZAxis(ZAxis));
}

TEST(SimTKCommon_SmallMatrix_UnitVec, AxisDotProducts) {
    EXPECT_EQ(XAxis.dotProduct(XAxis), 1);
    EXPECT_EQ(XAxis.dotProduct(YAxis), 0);
    EXPECT_EQ(XAxis.dotProduct(ZAxis), 0);
    EXPECT_EQ(ZAxis.dotProduct(YAxis), 0);
}

TEST(SimTKCommon_SmallMatrix_UnitVec, AxisCrossProductAxisAndSign) {
    // Parallel axes: sign must be zero
    EXPECT_EQ(XAxis.crossProductSign(XAxis), 0);

    // X × Y = +Z
    EXPECT_EQ(XAxis.crossProductAxis(YAxis), ZAxis);
    EXPECT_EQ(XAxis.crossProductSign(YAxis), 1);

    // X × Z = −Y
    EXPECT_EQ(XAxis.crossProductAxis(ZAxis), YAxis);
    EXPECT_EQ(XAxis.crossProductSign(ZAxis), -1);

    // Z × X = +Y
    EXPECT_EQ(ZAxis.crossProductAxis(XAxis), YAxis);
    EXPECT_EQ(ZAxis.crossProductSign(XAxis), 1);

    // Z × Y → NegX direction (out-parameter variant)
    int sign = 0;
    CoordinateAxis axis = ZAxis.crossProduct(YAxis, sign);
    EXPECT_TRUE(CoordinateDirection(axis, sign) == NegXAxis);

    // Explicit Negative() construction
    EXPECT_TRUE(CoordinateDirection(YAxis, CoordinateDirection::Negative()) == NegYAxis);

    // Negated-axis cross products
    EXPECT_EQ(NegXAxis.crossProductAxis(NegYAxis), ZAxis);
    EXPECT_EQ(NegYAxis.crossProductAxis(ZAxis), XAxis);
    EXPECT_EQ(NegYAxis.crossProductSign(ZAxis), -1);
    EXPECT_EQ(NegYAxis.crossProductSign(NegZAxis), 1);
}

TEST(SimTKCommon_SmallMatrix_UnitVec, NegationOperatorsAndDirectionAlgebra) {
    EXPECT_TRUE(-XAxis == NegXAxis);
    EXPECT_TRUE(-YAxis == NegYAxis);
    EXPECT_TRUE(-ZAxis == NegZAxis);
    EXPECT_TRUE(-ZAxis != NegYAxis);
    EXPECT_TRUE(ZAxis == -NegZAxis);

    // Double-negation laws
    EXPECT_TRUE(-(-NegXAxis) == -XAxis);
    EXPECT_TRUE(-(-NegYAxis) == -(-(-YAxis)));
}

// ============================================================================
// SimTKCommon_SmallMatrix_AppendRowCol
//
// Starting from a fixed 3×4 matrix, verifies every structural editing
// operation:
//   appendRow / insertRow at positions 0, 1, and last (== append)
//   appendCol / insertCol at positions 0, 2, and last (== append)
//   appendRowCol / insertRowCol at (0,0), (1,2), (3,4); corner element
//     comes from the new column vector, not from the new row vector
//   dropRow at all valid indices
//   dropCol at all valid indices
// ============================================================================

TEST(SimTKCommon_SmallMatrix_AppendRowCol, AppendAndInsertRowProducesCorrectMatrix) {
    const Mat34 m34(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

    const Row4 newRow(-1.0, -2.0, -3.0, -4.0);

    // Append row at the end
    const Mat44 appended = m34.appendRow(newRow);
    EXPECT_EQ(appended,
              Mat44(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, -1.0, -2.0, -3.0, -4.0));

    // Insert at position 0 (prepend)
    const Mat44 ins0 = m34.insertRow(0, newRow);
    EXPECT_EQ(ins0,
              Mat44(-1.0, -2.0, -3.0, -4.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0));

    // Insert at position 1 (middle)
    const Mat44 ins1 = m34.insertRow(1, newRow);
    EXPECT_EQ(ins1,
              Mat44(1.0, 2.0, 3.0, 4.0, -1.0, -2.0, -3.0, -4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0));

    // insertRow at last position must equal appendRow
    const Mat44 insLast = m34.insertRow(3, newRow);
    EXPECT_EQ(insLast, appended);
}

TEST(SimTKCommon_SmallMatrix_AppendRowCol, AppendAndInsertColProducesCorrectMatrix) {
    const Mat34 m34(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

    const Vec3 newCol(-6.0, -7.0, -8.0);

    // Append column at the end
    const Mat35 appended = m34.appendCol(newCol);
    EXPECT_EQ(appended,
              Mat35(1.0, 2.0, 3.0, 4.0, -6.0, 5.0, 6.0, 7.0, 8.0, -7.0, 9.0, 10.0, 11.0, 12.0, -8.0));

    // Insert at position 0 (prepend)
    const Mat35 ins0 = m34.insertCol(0, newCol);
    EXPECT_EQ(ins0, Mat35(-6.0, 1.0, 2.0, 3.0, 4.0, -7.0, 5.0, 6.0, 7.0, 8.0, -8.0, 9.0, 10.0, 11.0, 12.0));

    // Insert at position 2 (middle)
    const Mat35 ins2 = m34.insertCol(2, newCol);
    EXPECT_EQ(ins2, Mat35(1.0, 2.0, -6.0, 3.0, 4.0, 5.0, 6.0, -7.0, 7.0, 8.0, 9.0, 10.0, -8.0, 11.0, 12.0));

    // insertCol at last position must equal appendCol
    const Mat35 insLast = m34.insertCol(4, newCol);
    EXPECT_EQ(insLast, appended);
}

TEST(SimTKCommon_SmallMatrix_AppendRowCol, SimultaneousRowColInsertIsCorrect) {
    const Mat34 m34(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

    const Row5 newRow(-1.0, -2.0, -3.0, -4.0, -5.0);
    const Vec4 newCol(-6.0, -7.0, -8.0, -9.0);

    // appendRowCol: corner element comes from the column vector, not the row.
    const Mat45 appended = m34.appendRowCol(newRow, newCol);
    EXPECT_EQ(appended,
              Mat45(Row5(1.0, 2.0, 3.0, 4.0, -6.0),
                    Row5(5.0, 6.0, 7.0, 8.0, -7.0),
                    Row5(9.0, 10.0, 11.0, 12.0, -8.0),
                    Row5(-1.0, -2.0, -3.0, -4.0, -9.0)));

    // insertRowCol(0,0): both prepended; corner from newCol[0]
    const Mat45 ins00 = m34.insertRowCol(0, 0, newRow, newCol);
    EXPECT_EQ(ins00,
              Mat45(Row5(-6.0, -2.0, -3.0, -4.0, -5.0),
                    Row5(-7.0, 1.0, 2.0, 3.0, 4.0),
                    Row5(-8.0, 5.0, 6.0, 7.0, 8.0),
                    Row5(-9.0, 9.0, 10.0, 11.0, 12.0)));

    // insertRowCol(1,2): row at index 1, column at index 2
    const Mat45 ins12 = m34.insertRowCol(1, 2, newRow, newCol);
    EXPECT_EQ(ins12,
              Mat45(Row5(1.0, 2.0, -6.0, 3.0, 4.0),
                    Row5(-1.0, -2.0, -7.0, -4.0, -5.0),
                    Row5(5.0, 6.0, -8.0, 7.0, 8.0),
                    Row5(9.0, 10.0, -9.0, 11.0, 12.0)));

    // insertRowCol at (3,4) must equal appendRowCol
    const Mat45 insLast = m34.insertRowCol(3, 4, newRow, newCol);
    EXPECT_EQ(insLast, appended);
}

TEST(SimTKCommon_SmallMatrix_AppendRowCol, DropRowProducesCorrectMatrix) {
    const Mat34 m34(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

    EXPECT_EQ(m34.dropRow(0), Mat24(5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0));
    EXPECT_EQ(m34.dropRow(1), Mat24(1.0, 2.0, 3.0, 4.0, 9.0, 10.0, 11.0, 12.0));
    EXPECT_EQ(m34.dropRow(2), Mat24(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0));
}

TEST(SimTKCommon_SmallMatrix_AppendRowCol, DropColProducesCorrectMatrix) {
    const Mat34 m34(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

    EXPECT_EQ(m34.dropCol(0), Mat33(2.0, 3.0, 4.0, 6.0, 7.0, 8.0, 10.0, 11.0, 12.0));
    EXPECT_EQ(m34.dropCol(1), Mat33(1.0, 3.0, 4.0, 5.0, 7.0, 8.0, 9.0, 11.0, 12.0));
    EXPECT_EQ(m34.dropCol(2), Mat33(1.0, 2.0, 4.0, 5.0, 6.0, 8.0, 9.0, 10.0, 12.0));
    EXPECT_EQ(m34.dropCol(3), Mat33(1.0, 2.0, 3.0, 5.0, 6.0, 7.0, 9.0, 10.0, 11.0));
}
