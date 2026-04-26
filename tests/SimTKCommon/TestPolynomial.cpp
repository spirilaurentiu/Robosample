/**
 * @file TestPolynomial.cpp
 *
 * Google Test suite for SimTK::PolynomialRootFinder.
 *
 * What the original authors tested
 * ---------------------------------
 * 1. Quadratic polynomials with REAL coefficients
 *    - A handful of hand-chosen root pairs (0/0, 1/1, 0/5, 10/−5, trivial)
 *    - 1 000 randomly scaled polynomials with a repeated real root
 *    - 1 000 randomly scaled polynomials with two distinct real roots
 *    - 1 000 randomly scaled polynomials whose roots are a complex-conjugate
 *      pair (real coefficients ↔ conjugate roots)
 *    - That a zero leading coefficient throws ZeroLeadingCoefficient
 *
 * 2. Quadratic polynomials with COMPLEX coefficients
 *    - 1 000 randomly scaled polynomials with two arbitrary complex roots
 *    - That a zero leading coefficient throws ZeroLeadingCoefficient
 *
 * 3. Cubic polynomials with REAL coefficients – same families as quadratic:
 *    fixed cases, triple root, two-distinct, three-distinct, conjugate pair +
 *    real root, and the zero-leading-coefficient guard.
 *
 * 4. Cubic polynomials with COMPLEX coefficients
 *    - 1 000 randomly scaled polynomials with three arbitrary complex roots
 *    - Zero-leading-coefficient guard
 *
 * 5. Arbitrary-degree polynomials (degree 2–6)
 *    - 1 000 random real-coefficient polynomials verified by back-substitution
 *    - 1 000 random complex-coefficient polynomials verified by back-substitution
 *    - Zero-leading-coefficient guard for both overloads
 *
 * 6. Regression: a specific ill-conditioned 6th-degree polynomial from
 *    Ellipsoid::findNearestPoint() that silently returned *no* roots with the
 *    original Jenkins-Traub implementation (fixed in rpoly.cpp, 2013-04-10).
 */

#include <cmath>
#include <gtest/gtest.h>

#include "SimTKcommon.h"

using namespace SimTK;

namespace {

// ============================================================================
// Floating-point comparison helpers
// ============================================================================

/**
 * Returns true when both the real and imaginary relative errors between
 * @p expected and @p found are less than @p tol.  Absolute comparison is used
 * whenever the expected component is exactly zero.
 */
auto complexNearEqual(Complex expected, Complex found, Real tol) -> bool {
    Real diffr = std::abs(expected.real() - found.real());
    Real diffi = std::abs(expected.imag() - found.imag());
    Real scaler = std::max(std::abs(expected.real()), std::abs(found.real()));
    Real scalei = std::max(std::abs(expected.imag()), std::abs(found.imag()));
    if (expected.real() == 0.0) {
        scaler = 1.0;
    }
    if (expected.imag() == 0.0) {
        scalei = 1.0;
    }
    return (diffr < (tol * scaler)) && (diffi < (tol * scalei));
}

/**
 * High-accuracy comparison for quadratic roots (tolerance = sqrt(eps)).
 *
 * Additionally enforces that when the expected root is purely real the found
 * root must have an *exactly* zero imaginary part – the quadratic solver is
 * expected to be precise enough to never introduce spurious imaginary parts.
 */
auto quadraticRootsEqual(Complex expected, Complex found) -> bool {
    if (expected.imag() == 0.0 && found.imag() != 0.0) {
        return false;
    }
    return complexNearEqual(expected, found, std::sqrt(NTraits<Real>::getEps()));
}

// ============================================================================
// Polynomial-evaluation helpers (used for back-substitution)
// ============================================================================

auto evalPoly(const Vector_<Real>& coefficients, Complex value) -> Complex {
    Complex sum{0.0};
    for (int j = 0; j < coefficients.size(); ++j) {
        sum = (sum * value) + coefficients[j];
    }
    return sum;
}

auto evalPoly(const Vector_<Complex>& coefficients, Complex value) -> Complex {
    Complex sum{0.0};
    for (int j = 0; j < coefficients.size(); ++j) {
        sum = (sum * value) + coefficients[j];
    }
    return sum;
}

// ============================================================================
// Root-set comparison helpers that emit GTest assertions
// ============================================================================

/**
 * Asserts (order-independently) that the two-element root arrays @p expected
 * and @p found represent the same multiset, using quadratic-accuracy tolerance.
 */
auto expectQuadraticRootsMatch(Vec<2, Complex> expected, Vec<2, Complex> found) -> void {
    bool matched =
        (quadraticRootsEqual(expected[0], found[0]) && quadraticRootsEqual(expected[1], found[1]))
        || (quadraticRootsEqual(expected[0], found[1]) && quadraticRootsEqual(expected[1], found[0]));
    EXPECT_TRUE(matched) << "Expected roots: " << expected[0] << ",  " << expected[1] << "\n"
                         << "Found roots:    " << found[0] << ",  " << found[1];
}

/**
 * Asserts (order-independently) that the three-element root arrays @p expected
 * and @p found represent the same multiset within @p tol.
 */
auto expectCubicRootsMatch(Vec<3, Complex> expected, Vec<3, Complex> found, Real tol = 1e-4) -> void {
    auto eq = [tol](Complex a, Complex b) -> bool {
        return complexNearEqual(a, b, tol);
    };
    bool matched = false;
    if (eq(expected[0], found[0])) {
        matched = (eq(expected[1], found[1]) && eq(expected[2], found[2]))
                  || (eq(expected[1], found[2]) && eq(expected[2], found[1]));
    } else if (eq(expected[0], found[1])) {
        matched = (eq(expected[1], found[2]) && eq(expected[2], found[0]))
                  || (eq(expected[1], found[0]) && eq(expected[2], found[2]));
    } else if (eq(expected[0], found[2])) {
        matched = (eq(expected[1], found[0]) && eq(expected[2], found[1]))
                  || (eq(expected[1], found[1]) && eq(expected[2], found[0]));
    }
    EXPECT_TRUE(matched) << "Expected roots: " << expected[0] << ",  " << expected[1] << ",  " << expected[2]
                         << "\n"
                         << "Found roots:    " << found[0] << ",  " << found[1] << ",  " << found[2];
}

/**
 * Back-substitution check: asserts that every root in @p roots satisfies the
 * real-coefficient polynomial @p coefficients within @p tol.
 */
auto verifyRoots(const Vector_<Real>& coefficients, const Vector_<Complex>& roots, Real tol = 1e-2) -> void {
    for (int i = 0; i < roots.size(); ++i) {
        Complex residual = evalPoly(coefficients, roots[i]);
        EXPECT_TRUE(complexNearEqual(Complex{0.0}, residual, tol))
            << "Root[" << i << "] = " << roots[i] << " is not a root; residual = " << residual;
    }
}

/**
 * Overload for complex-coefficient polynomials.
 */
auto verifyRoots(const Vector_<Complex>& coefficients, const Vector_<Complex>& roots, Real tol = 1e-2)
    -> void {
    for (int i = 0; i < roots.size(); ++i) {
        Complex residual = evalPoly(coefficients, roots[i]);
        EXPECT_TRUE(complexNearEqual(Complex{0.0}, residual, tol))
            << "Root[" << i << "] = " << roots[i] << " is not a root; residual = " << residual;
    }
}

// ============================================================================
// Thin solver wrappers (call findRoots then delegate to comparison helpers)
// ============================================================================

auto solveAndCheckQuadratic(Vec3 coeff, Vec<2, Complex> expected) -> void {
    Vec<2, Complex> found;
    PolynomialRootFinder::findRoots(coeff, found);
    expectQuadraticRootsMatch(expected, found);
}

auto solveAndCheckQuadratic(Vec<3, Complex> coeff, Vec<2, Complex> expected) -> void {
    Vec<2, Complex> found;
    PolynomialRootFinder::findRoots(coeff, found);
    expectQuadraticRootsMatch(expected, found);
}

auto solveAndCheckCubic(Vec4 coeff, Vec<3, Complex> expected) -> void {
    Vec<3, Complex> found;
    PolynomialRootFinder::findRoots(coeff, found);
    expectCubicRootsMatch(expected, found);
}

auto solveAndCheckCubic(Vec<4, Complex> coeff, Vec<3, Complex> expected) -> void {
    Vec<3, Complex> found;
    PolynomialRootFinder::findRoots(coeff, found);
    expectCubicRootsMatch(expected, found);
}

// ============================================================================
// Polynomial-construction helpers
// (build the monic polynomial from known roots, apply a random scale, solve)
// ============================================================================

auto checkQuadraticFromRealRoots(Vec2 roots, Real scale) -> void {
    Vec3 coeff{1.0, -(roots[0] + roots[1]), roots[0] * roots[1]};
    coeff *= scale;
    solveAndCheckQuadratic(coeff, Vec<2, Complex>{roots[0], roots[1]});
}

auto checkQuadraticFromConjugateRoot(Complex root, Real scale) -> void {
    // For real-coefficient quadratic x^2 - 2·Re(r)·x + |r|^2
    Vec3 coeff{1.0, -(2.0 * root.real()), norm(root)};
    coeff *= scale;
    solveAndCheckQuadratic(coeff, Vec<2, Complex>{root, conj(root)});
}

auto checkQuadraticFromComplexRoots(Vec<2, Complex> roots, Complex scale) -> void {
    Vec<3, Complex> coeff{Complex{1.0}, -(roots[0] + roots[1]), roots[0] * roots[1]};
    coeff *= scale;
    solveAndCheckQuadratic(coeff, roots);
}

auto checkCubicFromRealRoots(Vec3 roots, Real scale) -> void {
    Vec4 coeff{1.0,
               -(roots[0] + roots[1] + roots[2]),
               (roots[0] * roots[1]) + (roots[1] * roots[2]) + (roots[2] * roots[0]),
               -(roots[0] * roots[1] * roots[2])};
    coeff *= scale;
    solveAndCheckCubic(coeff, Vec<3, Complex>{roots[0], roots[1], roots[2]});
}

auto checkCubicFromConjugateAndRealRoot(Complex root1, Real root2, Real scale) -> void {
    // Coefficients of (x - root1)(x - conj(root1))(x - root2) with real result
    Vec4 coeff{1.0,
               -(2.0 * root1.real()) - root2,
               norm(root1) + (2.0 * root1.real() * root2),
               -(norm(root1) * root2)};
    coeff *= scale;
    solveAndCheckCubic(coeff, Vec<3, Complex>{root1, conj(root1), root2});
}

auto checkCubicFromComplexRoots(Vec<3, Complex> roots, Complex scale) -> void {
    Vec<4, Complex> coeff{Complex{1.0},
                          -(roots[0] + roots[1] + roots[2]),
                          (roots[0] * roots[1]) + (roots[1] * roots[2]) + (roots[2] * roots[0]),
                          -(roots[0] * roots[1] * roots[2])};
    coeff *= scale;
    solveAndCheckCubic(coeff, roots);
}

} // namespace

// ============================================================================
// Quadratic – real coefficients
// ============================================================================

/**
 * Verifies PolynomialRootFinder for a small set of hand-chosen root pairs with
 * real coefficients, including the trivial all-zero case.
 */
TEST(SimTKCommon_Polynomial_QuadraticRealCoefficients, SolvesKnownFixedRoots) {
    // Trivial: x^2 = 0
    solveAndCheckQuadratic(Vec3{1.0, 0.0, 0.0}, Vec<2, Complex>{0.0, 0.0});
    checkQuadraticFromRealRoots(Vec2{0.0, 0.0}, 1.0);
    checkQuadraticFromRealRoots(Vec2{1.0, 1.0}, 1.0);
    checkQuadraticFromRealRoots(Vec2{0.0, 5.0}, 1.0);
    checkQuadraticFromRealRoots(Vec2{10.0, -5.0}, 1.0);
}

/**
 * Checks that the solver returns the correct double (repeated) real root for
 * 1 000 randomly scaled polynomials.
 *
 * The quadratic solver is expected to be accurate to sqrt(eps); the comparison
 * additionally requires the imaginary part of each found root to be exactly
 * zero when the expected root is real.
 */
TEST(SimTKCommon_Polynomial_QuadraticRealCoefficients, SolvesRepeatedRealRoot) {
    Random::Gaussian rng{0.0, 100.0};
    for (int i = 0; i < 1000; ++i) {
        Real root = rng.getValue();
        Real scale = rng.getValue();
        checkQuadraticFromRealRoots(Vec2{root, root}, scale);
    }
}

/**
 * Checks that the solver correctly recovers two distinct real roots for 1 000
 * randomly generated, randomly scaled quadratics.
 */
TEST(SimTKCommon_Polynomial_QuadraticRealCoefficients, SolvesTwoDistinctRealRoots) {
    Random::Gaussian rng{0.0, 100.0};
    for (int i = 0; i < 1000; ++i) {
        Real scale = rng.getValue();
        checkQuadraticFromRealRoots(Vec2{rng.getValue(), rng.getValue()}, scale);
    }
}

/**
 * Checks that when the two roots of a real-coefficient quadratic are a
 * complex-conjugate pair the solver returns exactly that pair (not merely two
 * complex numbers with small imaginary parts) for 1 000 random cases.
 */
TEST(SimTKCommon_Polynomial_QuadraticRealCoefficients, SolvesComplexConjugateRoots) {
    Random::Gaussian rng{0.0, 1000.0};
    for (int i = 0; i < 1000; ++i) {
        Real scale = rng.getValue();
        checkQuadraticFromConjugateRoot(Complex{rng.getValue(), rng.getValue()}, scale);
    }
}

/**
 * Verifies that passing a zero leading coefficient to the real-coefficient
 * quadratic overload throws ZeroLeadingCoefficient.
 */
TEST(SimTKCommon_Polynomial_QuadraticRealCoefficients, ThrowsOnZeroLeadingCoefficient) {
    Vec<2, Complex> found;
    EXPECT_THROW(PolynomialRootFinder::findRoots(Vec3{0.0, 1.0, 1.0}, found),
                 PolynomialRootFinder::ZeroLeadingCoefficient);
}

// ============================================================================
// Quadratic – complex coefficients
// ============================================================================

/**
 * Checks that the complex-coefficient quadratic solver correctly recovers two
 * arbitrary complex roots for 1 000 randomly scaled polynomials.
 */
TEST(SimTKCommon_Polynomial_QuadraticComplexCoefficients, SolvesTwoArbitraryComplexRoots) {
    Random::Gaussian rng{0.0, 100.0};
    for (int i = 0; i < 1000; ++i) {
        Complex scale{rng.getValue(), rng.getValue()};
        checkQuadraticFromComplexRoots(
            Vec<2, Complex>{Complex{rng.getValue(), rng.getValue()}, Complex{rng.getValue(), rng.getValue()}},
            scale);
    }
}

/**
 * Verifies that passing a zero leading coefficient to the complex-coefficient
 * quadratic overload throws ZeroLeadingCoefficient.
 */
TEST(SimTKCommon_Polynomial_QuadraticComplexCoefficients, ThrowsOnZeroLeadingCoefficient) {
    Vec<2, Complex> found;
    EXPECT_THROW(PolynomialRootFinder::findRoots(Vec<3, Complex>{0.0, 1.0, 1.0}, found),
                 PolynomialRootFinder::ZeroLeadingCoefficient);
}

// ============================================================================
// Cubic – real coefficients
// ============================================================================

/**
 * Verifies the cubic solver for a small set of hand-chosen root triples with
 * real coefficients, including the all-zero case.
 */
TEST(SimTKCommon_Polynomial_CubicRealCoefficients, SolvesKnownFixedRoots) {
    // Trivial: x^3 = 0
    solveAndCheckCubic(Vec4{1.0, 0.0, 0.0, 0.0}, Vec<3, Complex>{0.0, 0.0, 0.0});
    checkCubicFromRealRoots(Vec3{0.0, 0.0, 0.0}, 1.0);
    checkCubicFromRealRoots(Vec3{1.0, 1.0, 1.0}, 1.0);
    checkCubicFromRealRoots(Vec3{0.0, 5.0, 5.0}, 1.0);
    checkCubicFromRealRoots(Vec3{10.0, -5.0, 100.0}, 1.0);
}

/**
 * Checks that the solver handles a triple (repeated) real root for 1 000
 * randomly scaled cases.
 */
TEST(SimTKCommon_Polynomial_CubicRealCoefficients, SolvesTripleRepeatedRealRoot) {
    Random::Gaussian rng{0.0, 100.0};
    for (int i = 0; i < 1000; ++i) {
        Real root = rng.getValue();
        Real scale = rng.getValue();
        checkCubicFromRealRoots(Vec3{root, root, root}, scale);
    }
}

/**
 * Checks that the solver handles the case where two of three real roots are
 * equal (one repeated root + one distinct root) for 1 000 random cases.
 */
TEST(SimTKCommon_Polynomial_CubicRealCoefficients, SolvesOneRepeatedAndOneDistinctRealRoot) {
    Random::Gaussian rng{0.0, 100.0};
    for (int i = 0; i < 1000; ++i) {
        Real repeatedRoot = rng.getValue();
        Real scale = rng.getValue();
        checkCubicFromRealRoots(Vec3{repeatedRoot, repeatedRoot, rng.getValue()}, scale);
    }
}

/**
 * Checks that the solver correctly recovers three distinct real roots for
 * 1 000 randomly generated cubics.
 */
TEST(SimTKCommon_Polynomial_CubicRealCoefficients, SolvesThreeDistinctRealRoots) {
    Random::Gaussian rng{0.0, 100.0};
    for (int i = 0; i < 1000; ++i) {
        Real scale = rng.getValue();
        checkCubicFromRealRoots(Vec3{rng.getValue(), rng.getValue(), rng.getValue()}, scale);
    }
}

/**
 * Checks that the solver correctly recovers a complex-conjugate pair plus one
 * real root from a real-coefficient cubic for 1 000 random cases.
 */
TEST(SimTKCommon_Polynomial_CubicRealCoefficients, SolvesComplexConjugatePairPlusRealRoot) {
    Random::Gaussian rng{0.0, 1000.0};
    for (int i = 0; i < 1000; ++i) {
        Real scale = rng.getValue();
        checkCubicFromConjugateAndRealRoot(Complex{rng.getValue(), rng.getValue()}, rng.getValue(), scale);
    }
}

/**
 * Verifies that passing a zero leading coefficient to the real-coefficient
 * cubic overload throws ZeroLeadingCoefficient.
 */
TEST(SimTKCommon_Polynomial_CubicRealCoefficients, ThrowsOnZeroLeadingCoefficient) {
    Vec<3, Complex> found;
    EXPECT_THROW(PolynomialRootFinder::findRoots(Vec4{0.0, 0.0, 0.0, 0.0}, found),
                 PolynomialRootFinder::ZeroLeadingCoefficient);
}

// ============================================================================
// Cubic – complex coefficients
// ============================================================================

/**
 * Checks that the complex-coefficient cubic solver correctly recovers three
 * arbitrary complex roots for 1 000 randomly scaled polynomials.
 */
TEST(SimTKCommon_Polynomial_CubicComplexCoefficients, SolvesThreeArbitraryComplexRoots) {
    Random::Gaussian rng{0.0, 100.0};
    for (int i = 0; i < 1000; ++i) {
        Complex scale{rng.getValue(), rng.getValue()};
        checkCubicFromComplexRoots(Vec<3, Complex>{Complex{rng.getValue(), rng.getValue()},
                                                   Complex{rng.getValue(), rng.getValue()},
                                                   Complex{rng.getValue(), rng.getValue()}},
                                   scale);
    }
}

/**
 * Verifies that passing a zero leading coefficient to the complex-coefficient
 * cubic overload throws ZeroLeadingCoefficient.
 */
TEST(SimTKCommon_Polynomial_CubicComplexCoefficients, ThrowsOnZeroLeadingCoefficient) {
    Vec<3, Complex> found;
    EXPECT_THROW(PolynomialRootFinder::findRoots(Vec<4, Complex>{0.0, 0.0, 0.0, 0.0}, found),
                 PolynomialRootFinder::ZeroLeadingCoefficient);
}

// ============================================================================
// Arbitrary-degree polynomials
// ============================================================================

/**
 * Verifies the arbitrary-degree solver via back-substitution for 1 000 random
 * polynomials of degree 2–6 with real coefficients.
 *
 * Each root returned by findRoots() is substituted back into the polynomial;
 * the residual must be within 1 % of the polynomial's scale (tol = 1e-2).
 */
TEST(SimTKCommon_Polynomial_ArbitraryDegree, RealCoefficientRootsPassBackSubstitution) {
    Random::Uniform degreeRng{2.0, 7.0};
    Random::Gaussian valueRng{0.0, 100.0};
    for (int i = 0; i < 1000; ++i) {
        int length = degreeRng.getIntValue();
        Vector realCoeff(length);
        Vector_<Complex> roots(length - 1);
        for (int j = 0; j < length; ++j) {
            realCoeff[j] = valueRng.getValue();
        }
        PolynomialRootFinder::findRoots(realCoeff, roots);
        verifyRoots(realCoeff, roots);
    }
}

/**
 * Verifies the arbitrary-degree solver via back-substitution for 1 000 random
 * polynomials of degree 2–6 with complex coefficients.
 */
TEST(SimTKCommon_Polynomial_ArbitraryDegree, ComplexCoefficientRootsPassBackSubstitution) {
    Random::Uniform degreeRng{2.0, 7.0};
    Random::Gaussian valueRng{0.0, 100.0};
    for (int i = 0; i < 1000; ++i) {
        int length = degreeRng.getIntValue();
        Vector_<Complex> complexCoeff(length);
        Vector_<Complex> roots(length - 1);
        for (int j = 0; j < length; ++j) {
            complexCoeff[j] = Complex{valueRng.getValue(), valueRng.getValue()};
        }
        PolynomialRootFinder::findRoots(complexCoeff, roots);
        verifyRoots(complexCoeff, roots);
    }
}

/**
 * Verifies that a zero leading coefficient in a real-coefficient polynomial
 * throws ZeroLeadingCoefficient regardless of degree.
 */
TEST(SimTKCommon_Polynomial_ArbitraryDegree, RealCoefficientsThrowOnZeroLeadingCoefficient) {
    Random::Gaussian valueRng{0.0, 100.0};
    constexpr int kLength = 5;
    Vector realCoeff(kLength);
    for (int j = 0; j < kLength; ++j) {
        realCoeff[j] = valueRng.getValue();
    }
    realCoeff[0] = 0.0;
    Vector_<Complex> roots(kLength - 1);
    EXPECT_THROW(PolynomialRootFinder::findRoots(realCoeff, roots),
                 PolynomialRootFinder::ZeroLeadingCoefficient);
}

/**
 * Verifies that a zero leading coefficient in a complex-coefficient polynomial
 * throws ZeroLeadingCoefficient regardless of degree.
 */
TEST(SimTKCommon_Polynomial_ArbitraryDegree, ComplexCoefficientsThrowOnZeroLeadingCoefficient) {
    Random::Gaussian valueRng{0.0, 100.0};
    constexpr int kLength = 5;
    Vector_<Complex> complexCoeff(kLength);
    for (int j = 0; j < kLength; ++j) {
        complexCoeff[j] = Complex{valueRng.getValue(), valueRng.getValue()};
    }
    complexCoeff[0] = Complex{0.0};
    Vector_<Complex> roots(kLength - 1);
    EXPECT_THROW(PolynomialRootFinder::findRoots(complexCoeff, roots),
                 PolynomialRootFinder::ZeroLeadingCoefficient);
}

// ============================================================================
// Regression: ill-conditioned 6th-degree ellipsoid polynomial
// ============================================================================

/**
 * Regression test added by Sherm (2013-04-10).
 *
 * Three 6th-degree polynomials were generated by Ellipsoid::findNearestPoint().
 * "good1" and "good2" solved correctly with the original Jenkins-Traub code.
 * "bad" caused the solver to return *no* roots at all while appearing almost
 * identical to the other two polynomials.  A fix was applied to rpoly.cpp; this
 * test ensures the fix is never regressed.
 *
 * All roots are verified by back-substitution at tol = 1e-10 (much stricter
 * than the default 1e-2 used in the stochastic tests above).
 */
TEST(SimTKCommon_Polynomial_RegressionEllipsoid, JenkinsTraubSolverHandlesIllConditionedPolynomial) {
    Vector_<Complex> roots(6);

    // --- good1: solved correctly before the fix ---
    Vector good1(7);
    good1[0] = 1.0000000000000000;
    good1[1] = 0.093099999999999988;
    good1[2] = 0.00057211940879762883;
    good1[3] = 1.4343090324181468e-006;
    good1[4] = 1.6097307763625053e-009;
    good1[5] = 6.6189348690786845e-013;
    good1[6] = -3.0418139616145048e-018;
    PolynomialRootFinder::findRoots(good1, roots);
    verifyRoots(good1, roots, 1e-10);

    // --- good2: solved correctly before the fix ---
    Vector good2(7);
    good2[0] = 1.0000000000000000;
    good2[1] = 0.021700000000000004;
    good2[2] = 0.00013355532616269355;
    good2[3] = 1.8068473626980300e-007;
    good2[4] = 6.2414644021223684e-011;
    good2[5] = 6.4036661086410838e-015;
    good2[6] = -6.8627441483352105e-022;
    PolynomialRootFinder::findRoots(good2, roots);
    verifyRoots(good2, roots, 1e-10);

    // --- bad: broke the original Jenkins-Traub, returning no roots ---
    Vector bad(7);
    bad[0] = 1.0000000000000000;
    bad[1] = 0.021700000000000004;
    bad[2] = 2.9889970904696875e-005;
    bad[3] = 1.0901272298136685e-008;
    bad[4] = -4.4822782160985054e-012;
    bad[5] = -2.6193432740351220e-015;
    bad[6] = -3.0900602527225053e-019;
    PolynomialRootFinder::findRoots(bad, roots);
    verifyRoots(bad, roots, 1e-10);
}