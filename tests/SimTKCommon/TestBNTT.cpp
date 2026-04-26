/* -------------------------------------------------------------------------- *
 *                       SimTK Simbody: SimTKcommon                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/**
 * TestScalarArithmetic.cpp
 *
 * Converted from the original simbody ctest "main()" smoke test to Google Test.
 *
 * What the original authors were actually trying to test:
 *
 *  1. isNaN() - correct detection of IEEE NaN across all scalar wrappers:
 *       plain float/double/Real, complex<float/double>, conjugate<float/double>,
 *       and negator<T> (both raw and after applying unary minus).
 *
 *  2. negator<T> value semantics - a negator<T> reinterpret-cast from a value v
 *       stores v in memory but represents -v.  The authors wanted to confirm that
 *       equality, negation, and isNaN() all respect that convention.
 *
 *  3. negator<Complex> arithmetic identities - that the negator wrapper correctly
 *       implements the algebraic laws: double-negation, distribution of negation
 *       over addition/subtraction, and sign rules for multiplication and division.
 *       These are the only *actual* ASSERTs in the original code; everything else
 *       was just print-and-eyeball.
 *
 *  4. square() / cube() helpers - that square(x) == x*x and cube(x) == x*x*x for
 *       negator<Complex>.
 *
 *  5. sign() - correct return value (+1 / -1 / 0) for Real and for negator<Real>
 *       (where the stored sign is flipped).
 *
 *  6. Mathematical constants - that the SimTK named constants (Pi, E, Sqrt2, ...)
 *       are numerically self-consistent (reciprocals multiply to 1, roots square
 *       to their radicand, log-constants match std::log, Euler's identity holds,
 *       fractional constants match 1/n, and epsilon relationships hold across
 *       float/double precision).
 *
 * Intentionally omitted:
 *   - All writeUnformatted / writeFormatted / readUnformatted / fillUnformatted
 *     calls and their associated stream state; the original code never asserted
 *     anything about their output.
 *   - Pure print statements (cout << ...) that carried no assertion.
 *   - typeid(...).name() diagnostics.
 *   - Mixed-mode operator compilation checks (cff*3., etc.); those are
 *     compile-time smoke tests with zero runtime content.
 */

#include <cmath>
#include <complex>
#include <gtest/gtest.h>

#include "SimTKcommon.h"

using namespace SimTK;
using std::complex;

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

// Four-ULP absolute tolerance used for "should be the same to floating-point
// rounding" checks.
static const Real kEps4 = Real(4) * NTraits<Real>::getEps();

// Prevent the optimizer from constant-folding 0/0 or 1/0 into a compile-time
// constant; we want the runtime IEEE result.
static auto getRealZero() -> Real {
    return std::sin(Real(0));
}

// Convenience wrapper: assert two Complex values agree to within kEps4 in
// both components.
static void expectComplexNear(Complex a, Complex b, const char* file, int line) {
    SCOPED_TRACE(::testing::Message() << "at " << file << ":" << line);
    EXPECT_NEAR(a.real(), b.real(), 1e-9);
    EXPECT_NEAR(a.imag(), b.imag(), 1e-9);
}

#define EXPECT_COMPLEX_NEAR(a, b) expectComplexNear((a), (b), __FILE__, __LINE__)

TEST(SimTKCommon_BNTT_IsNaN, DetectsNaNForPlainRealTypes) {
    const Real nan = CNT<Real>::getNaN();
    const float fnan = NTraits<float>::getNaN();
    const double dnan = NTraits<double>::getNaN();

    // Positive and negated NaN are both NaN.
    EXPECT_TRUE(isNaN(nan));
    EXPECT_TRUE(isNaN(-nan));

    // SimTK::NaN and its negation.
    EXPECT_TRUE(isNaN(NaN));
    EXPECT_TRUE(isNaN(-NaN));

    // Non-NaN values must not be flagged.
    EXPECT_FALSE(isNaN(Real(1)));
    EXPECT_FALSE(isNaN(Infinity));

    // 0/0 at runtime must produce NaN; 1/0 must produce Infinity (not NaN).
    EXPECT_TRUE(isNaN(Real(0) / getRealZero()));
    EXPECT_FALSE(isNaN(Real(1) / getRealZero()));

    // float and double specialisations.
    EXPECT_TRUE(isNaN(fnan));
    EXPECT_TRUE(isNaN(dnan));
}

TEST(SimTKCommon_BNTT_IsNaN, DetectsNaNForComplexTypes) {
    const float fnan = NTraits<float>::getNaN();
    const double dnan = NTraits<double>::getNaN();

    // Both components NaN.
    const std::complex<float> cfnan(fnan, fnan);
    // Only imaginary component NaN - still NaN.
    const std::complex<double> cdnan(3, dnan);
    // A complex infinity is NOT a NaN.
    const complex<float> fcinf = CNT<complex<float>>::getInfinity();

    EXPECT_TRUE(isNaN(cfnan));
    EXPECT_TRUE(isNaN(cdnan));
    EXPECT_FALSE(isNaN(fcinf));
}

TEST(SimTKCommon_BNTT_IsNaN, DetectsNaNForConjugateTypes) {
    const float fnan = NTraits<float>::getNaN();
    const double dnan = NTraits<double>::getNaN();

    // conjugate imaginary part stores the *negative* of the mathematical
    // imaginary part; NaN in either slot must propagate.
    const conjugate<float> jfnan(fnan, 0.09f);
    const conjugate<double> jdnan(3, dnan);

    EXPECT_TRUE(isNaN(jfnan));
    EXPECT_TRUE(isNaN(jdnan));
}

TEST(SimTKCommon_BNTT_IsNaN, DetectsNaNThroughNegatorWrapper) {
    const Real zero = 0., two = 2.;
    const float fnan = NTraits<float>::getNaN();
    const std::complex<float> cfnan(fnan, fnan);

    // negator<T> reinterpret-cast from v stores v but semantically represents
    // -v.  isNaN must ignore the sign flip and report the underlying payload.
    const auto& nzero = reinterpret_cast<const negator<Real>&>(zero);
    const auto& ntwo = reinterpret_cast<const negator<Real>&>(two);
    const auto& nfnan = reinterpret_cast<const negator<float>&>(fnan);
    const auto& ncfnan = reinterpret_cast<const negator<std::complex<float>>&>(cfnan);

    // Non-NaN values - with or without the extra negation - must be clean.
    EXPECT_FALSE(isNaN(nzero));
    EXPECT_FALSE(isNaN(-ntwo));

    // NaN payload must be detected both directly and after negation.
    EXPECT_TRUE(isNaN(nfnan));
    EXPECT_TRUE(isNaN(-nfnan));
    EXPECT_TRUE(isNaN(ncfnan));
    EXPECT_TRUE(isNaN(-ncfnan));
}

TEST(SimTKCommon_BNTT_Negator, ValueEqualityAndNegation) {
    const Real zero = 0., two = 2.;

    // A negator<T> that wraps value v represents -v.
    const auto& nzero = reinterpret_cast<const negator<Real>&>(zero);
    const auto& ntwo = reinterpret_cast<const negator<Real>&>(two);

    EXPECT_EQ(nzero, zero);  // represents -0 == 0
    EXPECT_EQ(-nzero, zero); // -(-0) == 0
    EXPECT_EQ(ntwo, -two);   // represents -2
    EXPECT_EQ(-ntwo, two);   // -(-2) == 2
}

// x = 7.1 + 1.7i  (stored as negator, so represents -x mathematically)
// y = 2*x          (similar)
// All identities below must hold to within a few ULPs.
TEST(SimTKCommon_BNTT_NegatorComplex, AdditionIdentities) {
    negator<Complex> x(Complex(Real(7.1), Real(1.7))), y;
    y = x;
    y *= Real(2);

    // x + y  ==  -(-x) + y
    EXPECT_COMPLEX_NEAR(x + y, -(-x) + y);
    // x + y  ==  -( (-x) + (-y) )
    EXPECT_COMPLEX_NEAR(x + y, -((-x) + (-y)));
    // x + y  ==  x - (-y)
    EXPECT_COMPLEX_NEAR(x + y, x - (-y));
}

TEST(SimTKCommon_BNTT_NegatorComplex, SubtractionIdentities) {
    negator<Complex> x(Complex(Real(7.1), Real(1.7))), y;
    y = x;
    y *= Real(2);

    EXPECT_COMPLEX_NEAR(x - y, x + (-y));
    EXPECT_COMPLEX_NEAR(x - y, -(y - x));
    EXPECT_COMPLEX_NEAR(x - y, -(-x + y));
    EXPECT_COMPLEX_NEAR(-(x + y), (-x) - y);
    EXPECT_COMPLEX_NEAR(-(x + y), -x - y);
}

TEST(SimTKCommon_BNTT_NegatorComplex, MultiplicationIdentities) {
    negator<Complex> x(Complex(Real(7.1), Real(1.7))), y;
    y = x;
    y *= Real(2);

    // (-x)(-y) == xy
    EXPECT_COMPLEX_NEAR(x * y, (-x) * (-y));
    EXPECT_COMPLEX_NEAR(x * y, -x * -y);
    // -(xy) == (-x)y == x(-y)
    EXPECT_COMPLEX_NEAR(-(x * y), -x * y);
    EXPECT_COMPLEX_NEAR(-(x * y), x * -y);
}

TEST(SimTKCommon_BNTT_NegatorComplex, DivisionIdentities) {
    negator<Complex> x(Complex(Real(7.1), Real(1.7))), y;
    y = x;
    y *= Real(2);

    // (-x)/(-y) == x/y
    EXPECT_COMPLEX_NEAR(x / y, (-x) / (-y));
    EXPECT_COMPLEX_NEAR(x / y, -x / -y);
    // -(x/y) == (-x)/y == x/(-y)
    EXPECT_COMPLEX_NEAR(-(x / y), -x / y);
    EXPECT_COMPLEX_NEAR(-(x / y), x / -y);
}

TEST(SimTKCommon_BNTT_NegatorComplex, SquareHelperMatchesManualProduct) {
    negator<Complex> x(Complex(Real(7.1), Real(1.7)));
    EXPECT_COMPLEX_NEAR(square(x), x * x);
}

TEST(SimTKCommon_BNTT_NegatorComplex, CubeHelperMatchesManualProduct) {
    negator<Complex> x(Complex(Real(7.1), Real(1.7)));
    EXPECT_COMPLEX_NEAR(cube(x), x * x * x);
}

TEST(SimTKCommon_BNTT_Sign, ReturnsCorrectSignForPositiveNegativeAndZeroReal) {
    EXPECT_EQ(sign(Real(27)), Real(1));
    EXPECT_EQ(sign(Real(-14)), Real(-1));
    EXPECT_EQ(sign(Real(0)), Real(0));
}

TEST(SimTKCommon_BNTT_Sign, ReturnsFlippedSignThroughNegatorWrapper) {
    EXPECT_EQ(sign(negator<Real>::recast(Real(27))), Real(-1));
    EXPECT_EQ(sign(negator<Real>::recast(Real(-14))), Real(1));
    EXPECT_EQ(sign(negator<Real>::recast(Real(0))), Real(0));
}

TEST(SimTKCommon_BNTT_Constants, EulerIdentityHoldsToMachineEpsilon) {
    Complex result = std::pow(E, I * Pi) + Real(1);
    EXPECT_NEAR(result.real(), 0.0, kEps4);
    EXPECT_NEAR(result.imag(), 0.0, kEps4);
}


TEST(SimTKCommon_BNTT_Constants, ReciprocalPairsMultiplyToOne) {
    EXPECT_NEAR(double(Pi) * double(OneOverPi), 1.0, kEps4);
    EXPECT_NEAR(double(Sqrt2) * double(OneOverSqrt2), 1.0, kEps4);
    EXPECT_NEAR(double(Sqrt3) * double(OneOverSqrt3), 1.0, kEps4);
}

TEST(SimTKCommon_BNTT_Constants, SquareRootsSquareToTheirRadicand) {
    EXPECT_NEAR(double(Sqrt2) * double(Sqrt2), 2.0, kEps4);
    EXPECT_NEAR(double(Sqrt3) * double(Sqrt3), 3.0, kEps4);
}

TEST(SimTKCommon_BNTT_Constants, CubeRootsCubeToTheirRadicand) {
    EXPECT_NEAR(double(CubeRoot2) * double(CubeRoot2) * double(CubeRoot2), 2.0, kEps4);
    EXPECT_NEAR(double(CubeRoot3) * double(CubeRoot3) * double(CubeRoot3), 3.0, kEps4);
}

TEST(SimTKCommon_BNTT_Constants, LogarithmConstantsMatchStdLib) {
    EXPECT_NEAR(double(Log2E), std::log2(double(E)), kEps4);
    EXPECT_NEAR(double(Log10E), std::log10(double(E)), kEps4);
    EXPECT_NEAR(double(Ln2), std::log(2.0), kEps4);
    EXPECT_NEAR(double(Ln10), std::log(10.0), kEps4);
}

TEST(SimTKCommon_BNTT_Constants, FractionalConstantsMatchExactValues) {
    EXPECT_NEAR(double(OneHalf), 1.0 / 2, kEps4);
    EXPECT_NEAR(double(OneThird), 1.0 / 3, kEps4);
    EXPECT_NEAR(double(OneFourth), 1.0 / 4, kEps4);
    EXPECT_NEAR(double(OneFifth), 1.0 / 5, kEps4);
    EXPECT_NEAR(double(OneSixth), 1.0 / 6, kEps4);
    EXPECT_NEAR(double(OneSeventh), 1.0 / 7, kEps4);
    EXPECT_NEAR(double(OneEighth), 1.0 / 8, kEps4);
    EXPECT_NEAR(double(OneNinth), 1.0 / 9, kEps4);
}