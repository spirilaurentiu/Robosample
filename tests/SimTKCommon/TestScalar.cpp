// TestScalar_GTest.cpp
//
// Google Test port of the SimTK scalar utility function tests
// (originally TestScalar.cpp in the SimTK CTest suite).
//
// Covered functions:
//   isNaN           — detects NaN in scalar / complex / conjugate / negator types
//   isInf           — detects infinity (not NaN) in the same type family
//   isFinite        — all components must be finite
//   signBit         — raw sign-bit inspection (same bit as underlying for negator<T>)
//   sign            — -1 / 0 / +1 value (negated interpretation for negator<T>)
//   square / cube   — x² / x³ with negator and complex support
//   isNumericallyEqual — tolerance-based equality for float / double / negator
//   clamp / clampInPlace
//   stepUp / stepDown and their derivatives (dstepUp … d3stepDown)
//   stepAny / dstepAny / d2stepAny / d3stepAny

#include <complex>
#include <gtest/gtest.h>

#include "SimTKcommon.h"
#include "Util.hpp"

using namespace SimTK;

// Convenience wrapper: invokes AssertSimTKEqual (Util.hpp) and passes the
// result to EXPECT_TRUE so failures include variable names and values.
#define SIMTK_EXPECT_EQ(actual, expected) \
    EXPECT_TRUE(AssertSimTKEqual(#actual, #expected, (actual), (expected)))

// ============================================================
// isNaN
//
// isNaN(x) must return true when x contains any NaN component.
// For complex and conjugate types it is sufficient for either the
// real or the imaginary part to be NaN.
// ============================================================

// Scalar float and double NaN values are correctly detected.
TEST(SimTKCommon_Scalar_IsNaN, ScalarNaN_IsDetected) {
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();

    EXPECT_TRUE(isNaN(fltNaN));
    EXPECT_TRUE(isNaN(dblNaN));
}

// Negating a NaN value yields another NaN; isNaN must still return true.
TEST(SimTKCommon_Scalar_IsNaN, NegatedScalarNaN_IsStillNaN) {
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    const float nfltNaN = -fltNaN;
    const double ndblNaN = -dblNaN;

    EXPECT_TRUE(isNaN(nfltNaN));
    EXPECT_TRUE(isNaN(ndblNaN));
}

// Regular (finite, non-NaN) scalar values are not NaN.
TEST(SimTKCommon_Scalar_IsNaN, RegularScalar_IsNotNaN) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;

    EXPECT_FALSE(isNaN(fltRegular));
    EXPECT_FALSE(isNaN(dblRegular));
}

// Regular complex<T> and conjugate<T> values (finite real and imaginary
// parts) are not NaN.
TEST(SimTKCommon_Scalar_IsNaN, RegularComplexAndConjugate_AreNotNaN) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;

    const std::complex<float> cflt(fltRegular, (-2.0F * fltRegular));
    const std::complex<double> cdbl(dblRegular, (-2.0 * dblRegular));
    const conjugate<float> cjflt(fltRegular, (-2.0F * fltRegular));
    const conjugate<double> cjdbl(dblRegular, (-2.0 * dblRegular));

    EXPECT_FALSE(isNaN(cflt));
    EXPECT_FALSE(isNaN(cdbl));
    EXPECT_FALSE(isNaN(cjflt));
    EXPECT_FALSE(isNaN(cjdbl));
}

// negator<T> reinterprets the same memory with flipped sign semantics.
// Wrapping a regular scalar must not produce a NaN.
// This test also verifies that the negators themselves evaluate correctly.
TEST(SimTKCommon_Scalar_IsNaN, NegatorOfRegularScalar_IsNotNaN) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;

    const negator<float>& nflt = reinterpret_cast<const negator<float>&>(fltRegular);
    const negator<double>& ndbl = reinterpret_cast<const negator<double>&>(dblRegular);

    // Verify negator semantics before the isNaN check.
    SIMTK_EXPECT_EQ(nflt, -fltRegular);
    SIMTK_EXPECT_EQ(ndbl, -dblRegular);

    EXPECT_FALSE(isNaN(nflt));
    EXPECT_FALSE(isNaN(ndbl));
}

// negator<complex<T>> and negator<conjugate<T>> wrapping regular values
// are not NaN.  The negators are verified for correctness first.
TEST(SimTKCommon_Scalar_IsNaN, NegatorOfRegularComplexAndConjugate_AreNotNaN) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;

    std::complex<float> cflt(fltRegular, (-2.0F * fltRegular));
    std::complex<double> cdbl(dblRegular, (-2.0 * dblRegular));
    conjugate<float> cjflt(fltRegular, (-2.0F * fltRegular));
    conjugate<double> cjdbl(dblRegular, (-2.0 * dblRegular));

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    // Verify negator semantics.
    SIMTK_EXPECT_EQ(ncflt, -cflt);
    SIMTK_EXPECT_EQ(-ncflt, cflt);
    SIMTK_EXPECT_EQ(ncjflt, -cjflt);
    SIMTK_EXPECT_EQ(-ncjflt, cjflt);

    EXPECT_FALSE(isNaN(ncflt));
    EXPECT_FALSE(isNaN(ncdbl));
    EXPECT_FALSE(isNaN(ncjflt));
    EXPECT_FALSE(isNaN(ncjdbl));
}

// Should be NaN if either or both parts are NaN.
// Here: real is finite, imaginary is NaN.
TEST(SimTKCommon_Scalar_IsNaN, ComplexWithImaginaryPartNaN_IsNaN) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();

    std::complex<float> cflt(fltRegular, fltNaN);
    std::complex<double> cdbl(dblRegular, dblNaN);
    conjugate<float> cjflt(fltRegular, fltNaN);
    conjugate<double> cjdbl(dblRegular, dblNaN);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    EXPECT_TRUE(isNaN(cflt));
    EXPECT_TRUE(isNaN(cdbl));
    EXPECT_TRUE(isNaN(cjflt));
    EXPECT_TRUE(isNaN(cjdbl));
    EXPECT_TRUE(isNaN(ncflt));
    EXPECT_TRUE(isNaN(ncdbl));
    EXPECT_TRUE(isNaN(ncjflt));
    EXPECT_TRUE(isNaN(ncjdbl));
}

// Should be NaN if either or both parts are NaN.
// Here: both real and imaginary are NaN.
TEST(SimTKCommon_Scalar_IsNaN, ComplexWithBothPartsNaN_IsNaN) {
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();

    std::complex<float> cflt(fltNaN, fltNaN);
    std::complex<double> cdbl(dblNaN, dblNaN);
    conjugate<float> cjflt(fltNaN, fltNaN);
    conjugate<double> cjdbl(dblNaN, dblNaN);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    EXPECT_TRUE(isNaN(cflt));
    EXPECT_TRUE(isNaN(cdbl));
    EXPECT_TRUE(isNaN(cjflt));
    EXPECT_TRUE(isNaN(cjdbl));
    EXPECT_TRUE(isNaN(ncflt));
    EXPECT_TRUE(isNaN(ncdbl));
    EXPECT_TRUE(isNaN(ncjflt));
    EXPECT_TRUE(isNaN(ncjdbl));
}

// Real part only is NaN; imaginary part is finite.
TEST(SimTKCommon_Scalar_IsNaN, ComplexWithRealPartOnlyNaN_IsNaN) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();

    std::complex<float> cflt(fltNaN, fltRegular);
    std::complex<double> cdbl(dblNaN, dblRegular);
    conjugate<float> cjflt(fltNaN, fltRegular);
    conjugate<double> cjdbl(dblNaN, dblRegular);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    EXPECT_TRUE(isNaN(cflt));
    EXPECT_TRUE(isNaN(cdbl));
    EXPECT_TRUE(isNaN(cjflt));
    EXPECT_TRUE(isNaN(cjdbl));
    EXPECT_TRUE(isNaN(ncflt));
    EXPECT_TRUE(isNaN(ncdbl));
    EXPECT_TRUE(isNaN(ncjflt));
    EXPECT_TRUE(isNaN(ncjdbl));
}

// ============================================================
// isInf
//
// isInf(x) must return true when x contains at least one infinite
// component AND no NaN components.  For complex / conjugate types,
// either part being infinite is sufficient, but a NaN in any part
// disqualifies the whole.
// ============================================================

// Scalar +∞ and −∞ are both detected as infinite.
TEST(SimTKCommon_Scalar_IsInf, ScalarInfIncludingNegative_IsDetected) {
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    EXPECT_TRUE(isInf(fltInf));
    EXPECT_TRUE(isInf(dblInf));
    EXPECT_TRUE(isInf(mfltInf));
    EXPECT_TRUE(isInf(mdblInf));
}

// The negator of +∞ stores +∞ bits but evaluates as −∞; both are infinite.
TEST(SimTKCommon_Scalar_IsInf, NegatorOfInf_IsStillInf) {
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();

    const negator<float>& nfltInf = reinterpret_cast<const negator<float>&>(fltInf);
    const negator<double>& ndblInf = reinterpret_cast<const negator<double>&>(dblInf);

    EXPECT_EQ(nfltInf, -fltInf);
    EXPECT_EQ(ndblInf, -dblInf);

    EXPECT_TRUE(isInf(nfltInf));
    EXPECT_TRUE(isInf(ndblInf));
}

// Regular scalar, complex<T>, and conjugate<T> values are not infinite.
TEST(SimTKCommon_Scalar_IsInf, RegularValues_AreNotInf) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;

    const std::complex<float> cflt(fltRegular, (-2.0F * fltRegular));
    const std::complex<double> cdbl(dblRegular, (-2.0 * dblRegular));
    const conjugate<float> cjflt(fltRegular, (-2.0F * fltRegular));
    const conjugate<double> cjdbl(dblRegular, (-2.0 * dblRegular));

    EXPECT_FALSE(isInf(fltRegular));
    EXPECT_FALSE(isInf(dblRegular));
    EXPECT_FALSE(isInf(cflt));
    EXPECT_FALSE(isInf(cdbl));
    EXPECT_FALSE(isInf(cjflt));
    EXPECT_FALSE(isInf(cjdbl));
}

// negator<complex<T>> and negator<conjugate<T>> wrapping finite values
// are not infinite.  Negator semantics are verified first.
TEST(SimTKCommon_Scalar_IsInf, NegatorOfRegular_BehavesCorrectlyAndIsNotInf) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;

    std::complex<float> cflt(fltRegular, (-2.0F * fltRegular));
    std::complex<double> cdbl(dblRegular, (-2.0 * dblRegular));
    conjugate<float> cjflt(fltRegular, (-2.0F * fltRegular));
    conjugate<double> cjdbl(dblRegular, (-2.0 * dblRegular));

    const negator<float>& nflt = reinterpret_cast<const negator<float>&>(fltRegular);
    const negator<double>& ndbl = reinterpret_cast<const negator<double>&>(dblRegular);
    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    // Verify negator semantics.
    SIMTK_EXPECT_EQ(nflt, -fltRegular);
    SIMTK_EXPECT_EQ(ndbl, -dblRegular);
    SIMTK_EXPECT_EQ(ncflt, -cflt);
    SIMTK_EXPECT_EQ(-ncflt, cflt);
    SIMTK_EXPECT_EQ(ncjflt, -cjflt);
    SIMTK_EXPECT_EQ(-ncjflt, cjflt);

    EXPECT_FALSE(isInf(nflt));
    EXPECT_FALSE(isInf(ndbl));
    EXPECT_FALSE(isInf(ncflt));
    EXPECT_FALSE(isInf(ncdbl));
    EXPECT_FALSE(isInf(ncjflt));
    EXPECT_FALSE(isInf(ncjdbl));
}

// Should be Inf if either or both parts are Inf, as long as neither
// part is NaN.  Here: real is finite, imaginary is +∞.
TEST(SimTKCommon_Scalar_IsInf, ComplexWithImaginaryPartInf_IsInf) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();

    std::complex<float> cflt(fltRegular, fltInf);
    std::complex<double> cdbl(dblRegular, dblInf);
    conjugate<float> cjflt(fltRegular, fltInf);
    conjugate<double> cjdbl(dblRegular, dblInf);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    // Imaginary only is Inf.
    EXPECT_TRUE(isInf(cflt));
    EXPECT_TRUE(isInf(cdbl));
    EXPECT_TRUE(isInf(cjflt));
    EXPECT_TRUE(isInf(cjdbl));
    EXPECT_TRUE(isInf(ncflt));
    EXPECT_TRUE(isInf(ncdbl));
    EXPECT_TRUE(isInf(ncjflt));
    EXPECT_TRUE(isInf(ncjdbl));
}

// Both parts are Inf.
TEST(SimTKCommon_Scalar_IsInf, ComplexWithBothPartsInf_IsInf) {
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();

    std::complex<float> cflt(fltInf, fltInf);
    std::complex<double> cdbl(dblInf, dblInf);
    conjugate<float> cjflt(fltInf, fltInf);
    conjugate<double> cjdbl(dblInf, dblInf);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    EXPECT_TRUE(isInf(cflt));
    EXPECT_TRUE(isInf(cdbl));
    EXPECT_TRUE(isInf(cjflt));
    EXPECT_TRUE(isInf(cjdbl));
    EXPECT_TRUE(isInf(ncflt));
    EXPECT_TRUE(isInf(ncdbl));
    EXPECT_TRUE(isInf(ncjflt));
    EXPECT_TRUE(isInf(ncjdbl));
}

// Real part only is Inf (imaginary is finite), including negative infinity.
TEST(SimTKCommon_Scalar_IsInf, ComplexWithRealPartOnlyInf_IncludingNegative_IsInf) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    // Positive infinity real part.
    std::complex<float> cflt(fltInf, fltRegular);
    std::complex<double> cdbl(dblInf, dblRegular);
    conjugate<float> cjflt(fltInf, fltRegular);
    conjugate<double> cjdbl(dblInf, dblRegular);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    // Real part only is Inf.
    EXPECT_TRUE(isInf(cflt));
    EXPECT_TRUE(isInf(cdbl));
    EXPECT_TRUE(isInf(cjflt));
    EXPECT_TRUE(isInf(cjdbl));
    EXPECT_TRUE(isInf(ncflt));
    EXPECT_TRUE(isInf(ncdbl));
    EXPECT_TRUE(isInf(ncjflt));
    EXPECT_TRUE(isInf(ncjdbl));

    // Set real part to minus infinity — still infinite.
    cflt = std::complex<float>(mfltInf, cflt.imag());
    cdbl = std::complex<double>(mdblInf, cdbl.imag());
    cjflt = conjugate<float>(mfltInf, cjflt.imag());
    cjdbl = conjugate<double>(mdblInf, cjdbl.imag());

    EXPECT_TRUE(isInf(cflt));
    EXPECT_TRUE(isInf(cdbl));
    EXPECT_TRUE(isInf(cjflt));
    EXPECT_TRUE(isInf(cjdbl));
    EXPECT_TRUE(isInf(ncflt));
    EXPECT_TRUE(isInf(ncdbl));
    EXPECT_TRUE(isInf(ncjflt));
    EXPECT_TRUE(isInf(ncjdbl));
}

// A NaN in any part disqualifies isInf — even if the other part is finite.
TEST(SimTKCommon_Scalar_IsInf, ComplexWithNaNPart_IsNotInf) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();

    std::complex<float> cflt(fltNaN, fltRegular);
    std::complex<double> cdbl(dblNaN, dblRegular);
    conjugate<float> cjflt(fltNaN, fltRegular);
    conjugate<double> cjdbl(dblNaN, dblRegular);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    EXPECT_FALSE(isInf(cflt));
    EXPECT_FALSE(isInf(cdbl));
    EXPECT_FALSE(isInf(cjflt));
    EXPECT_FALSE(isInf(cjdbl));
    EXPECT_FALSE(isInf(ncflt));
    EXPECT_FALSE(isInf(ncdbl));
    EXPECT_FALSE(isInf(ncjflt));
    EXPECT_FALSE(isInf(ncjdbl));
}

// ============================================================
// isFinite
//
// isFinite(x) returns true only when all components of x are finite
// (neither NaN nor Inf).  For complex / conjugate types, both real
// and imaginary parts must be finite.
// ============================================================

// Regular (finite, non-NaN) scalar values are finite.
TEST(SimTKCommon_Scalar_IsFinite, RegularScalar_IsFinite) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;

    EXPECT_TRUE(isFinite(fltRegular));
    EXPECT_TRUE(isFinite(dblRegular));
}

// NaN and ±∞ are not finite.
TEST(SimTKCommon_Scalar_IsFinite, NaNAndInfinity_AreNotFinite) {
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    const float nfltNaN = -fltNaN;
    const double ndblNaN = -dblNaN;
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    EXPECT_FALSE(isFinite(fltNaN));
    EXPECT_FALSE(isFinite(dblNaN));
    EXPECT_FALSE(isFinite(nfltNaN));
    EXPECT_FALSE(isFinite(ndblNaN));
    EXPECT_FALSE(isFinite(fltInf));
    EXPECT_FALSE(isFinite(dblInf));
    EXPECT_FALSE(isFinite(mfltInf));
    EXPECT_FALSE(isFinite(mdblInf));
}

// complex<T> and conjugate<T> with both parts finite are finite.
TEST(SimTKCommon_Scalar_IsFinite, RegularComplexAndConjugate_AreFinite) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;

    const std::complex<float> cflt(fltRegular, (-2.0F * fltRegular));
    const std::complex<double> cdbl(dblRegular, (-2.0 * dblRegular));
    const conjugate<float> cjflt(fltRegular, (-2.0F * fltRegular));
    const conjugate<double> cjdbl(dblRegular, (-2.0 * dblRegular));

    EXPECT_TRUE(isFinite(cflt));
    EXPECT_TRUE(isFinite(cdbl));
    EXPECT_TRUE(isFinite(cjflt));
    EXPECT_TRUE(isFinite(cjdbl));
}

// negator<T> wrapping finite values is finite.  Negator semantics verified.
TEST(SimTKCommon_Scalar_IsFinite, NegatorOfRegular_BehavesCorrectlyAndIsFinite) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;

    std::complex<float> cflt(fltRegular, (-2.0F * fltRegular));
    std::complex<double> cdbl(dblRegular, (-2.0 * dblRegular));
    conjugate<float> cjflt(fltRegular, (-2.0F * fltRegular));
    conjugate<double> cjdbl(dblRegular, (-2.0 * dblRegular));

    const negator<float>& nflt = reinterpret_cast<const negator<float>&>(fltRegular);
    const negator<double>& ndbl = reinterpret_cast<const negator<double>&>(dblRegular);
    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    // Verify negator semantics.
    SIMTK_EXPECT_EQ(nflt, -fltRegular);
    SIMTK_EXPECT_EQ(ndbl, -dblRegular);
    SIMTK_EXPECT_EQ(ncflt, -cflt);
    SIMTK_EXPECT_EQ(-ncflt, cflt);
    SIMTK_EXPECT_EQ(ncjflt, -cjflt);
    SIMTK_EXPECT_EQ(-ncjflt, cjflt);

    EXPECT_TRUE(isFinite(nflt));
    EXPECT_TRUE(isFinite(ndbl));
    EXPECT_TRUE(isFinite(ncflt));
    EXPECT_TRUE(isFinite(ncdbl));
    EXPECT_TRUE(isFinite(ncjflt));
    EXPECT_TRUE(isFinite(ncjdbl));
}

// Should be finite only if both parts are finite.
// Here: imaginary part is non-finite (Inf or NaN), real part is finite.
TEST(SimTKCommon_Scalar_IsFinite, ComplexWithNonFiniteImaginaryPart_IsNotFinite) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    // Mix of non-finite imaginary values across the four types.
    std::complex<float> cflt(fltRegular, fltInf);
    std::complex<double> cdbl(dblRegular, mdblInf);
    conjugate<float> cjflt(fltRegular, fltNaN);
    conjugate<double> cjdbl(dblRegular, dblInf);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    // Imaginary part only is non-finite.
    EXPECT_FALSE(isFinite(cflt));
    EXPECT_FALSE(isFinite(cdbl));
    EXPECT_FALSE(isFinite(cjflt));
    EXPECT_FALSE(isFinite(cjdbl));
    EXPECT_FALSE(isFinite(ncflt));
    EXPECT_FALSE(isFinite(ncdbl));
    EXPECT_FALSE(isFinite(ncjflt));
    EXPECT_FALSE(isFinite(ncjdbl));
}

// Both parts are non-finite.
TEST(SimTKCommon_Scalar_IsFinite, ComplexWithBothPartsNonFinite_IsNotFinite) {
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    std::complex<float> cflt(fltInf, fltInf);
    std::complex<double> cdbl(mdblInf, mdblInf);
    conjugate<float> cjflt(fltNaN, fltNaN);
    conjugate<double> cjdbl(dblInf, dblInf);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    EXPECT_FALSE(isFinite(cflt));
    EXPECT_FALSE(isFinite(cdbl));
    EXPECT_FALSE(isFinite(cjflt));
    EXPECT_FALSE(isFinite(cjdbl));
    EXPECT_FALSE(isFinite(ncflt));
    EXPECT_FALSE(isFinite(ncdbl));
    EXPECT_FALSE(isFinite(ncjflt));
    EXPECT_FALSE(isFinite(ncjdbl));
}

// Real part only is non-finite; imaginary part is restored to regular.
TEST(SimTKCommon_Scalar_IsFinite, ComplexWithNonFiniteRealPart_IsNotFinite) {
    const float fltRegular = -12.34F;
    const double dblRegular = -12.34;
    const float fltNaN = NTraits<float>::getNaN();
    const double dblNaN = NTraits<double>::getNaN();
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    std::complex<float> cflt(fltInf, fltRegular);
    std::complex<double> cdbl(mdblInf, dblRegular);
    conjugate<float> cjflt(fltNaN, fltRegular);
    conjugate<double> cjdbl(dblInf, dblRegular);

    negator<std::complex<float>>& ncflt = reinterpret_cast<negator<std::complex<float>>&>(cflt);
    negator<std::complex<double>>& ncdbl = reinterpret_cast<negator<std::complex<double>>&>(cdbl);
    negator<conjugate<float>>& ncjflt = reinterpret_cast<negator<conjugate<float>>&>(cjflt);
    negator<conjugate<double>>& ncjdbl = reinterpret_cast<negator<conjugate<double>>&>(cjdbl);

    // Real part only is non-finite.
    EXPECT_FALSE(isFinite(cflt));
    EXPECT_FALSE(isFinite(cdbl));
    EXPECT_FALSE(isFinite(cjflt));
    EXPECT_FALSE(isFinite(cjdbl));
    EXPECT_FALSE(isFinite(ncflt));
    EXPECT_FALSE(isFinite(ncdbl));
    EXPECT_FALSE(isFinite(ncjflt));
    EXPECT_FALSE(isFinite(ncjdbl));
}

// ============================================================
// signBit
//
// signBit(x) inspects the raw IEEE sign bit.
// Unsigned types have no sign bit → always false.
// For negator<T>, signBit returns the SAME result as for the
// underlying T (the bit is not flipped — only its interpretation
// changes).
// ============================================================

// Unsigned integer types have no sign bit, so signBit is always false
// regardless of the stored value.
TEST(SimTKCommon_Scalar_SignBit, UnsignedTypes_NeverHaveSignBit) {
    const unsigned char ucm = 0xFFU;
    const unsigned char ucz = 0U;
    const unsigned char ucp = 27U;
    const unsigned short usm = 0xFFFFU;
    const unsigned short usz = 0U;
    const unsigned short usp = 2342U;
    const unsigned int uim = 0xFFFFFFFFU;
    const unsigned int uiz = 0U;
    const unsigned int uip = 2342344U;
    const auto ulm = static_cast<unsigned long>(-23423L);
    const unsigned long ulz = 0UL;
    const unsigned long ulp = 234234UL;
    const auto ullm = static_cast<unsigned long long>(-234234234LL);
    const unsigned long long ullz = 0ULL;
    const unsigned long long ullp = 234234234ULL;

    EXPECT_FALSE(signBit(ucm));
    EXPECT_FALSE(signBit(ucz));
    EXPECT_FALSE(signBit(ucp));
    EXPECT_FALSE(signBit(usm));
    EXPECT_FALSE(signBit(usz));
    EXPECT_FALSE(signBit(usp));
    EXPECT_FALSE(signBit(uim));
    EXPECT_FALSE(signBit(uiz));
    EXPECT_FALSE(signBit(uip));
    EXPECT_FALSE(signBit(ulm));
    EXPECT_FALSE(signBit(ulz));
    EXPECT_FALSE(signBit(ulp));
    EXPECT_FALSE(signBit(ullm));
    EXPECT_FALSE(signBit(ullz));
    EXPECT_FALSE(signBit(ullp));
}

// Signed integer types: signBit is true only for strictly negative values.
// Note: signBit(char) is not provided by SimTK; use signed char.
TEST(SimTKCommon_Scalar_SignBit, SignedIntegerTypes_SignBitReflectsSign) {
    const signed char cm = -23;
    const signed char cz = 0;
    const signed char cp = 99;
    const short sm = -1234;
    const short sz = 0;
    const short sp = 23423;
    const int im = -2342343;
    const int iz = 0;
    const int ip = 29472383;
    const long lm = -43488L;
    const long lz = 0L;
    const long lp = 3454545L;
    const long long llm = -2342342343433LL;
    const long long llz = 0LL;
    const long long llp = 874578478478574LL;

    EXPECT_TRUE(signBit(cm));
    EXPECT_FALSE(signBit(cz));
    EXPECT_FALSE(signBit(cp));
    EXPECT_TRUE(signBit(sm));
    EXPECT_FALSE(signBit(sz));
    EXPECT_FALSE(signBit(sp));
    EXPECT_TRUE(signBit(im));
    EXPECT_FALSE(signBit(iz));
    EXPECT_FALSE(signBit(ip));
    EXPECT_TRUE(signBit(lm));
    EXPECT_FALSE(signBit(lz));
    EXPECT_FALSE(signBit(lp));
    EXPECT_TRUE(signBit(llm));
    EXPECT_FALSE(signBit(llz));
    EXPECT_FALSE(signBit(llp));
}

// Float and double: signBit correctly reflects the IEEE sign bit.
TEST(SimTKCommon_Scalar_SignBit, FloatAndDouble_SignBitReflectsSign) {
    const float fm = -12398.34F;
    const float fz = 0.0F;
    const float fp = 4354.331F;
    const double dm = -234234.454;
    const double dz = 0.0;
    const double dp = 345345.2342;

    EXPECT_TRUE(signBit(fm));
    EXPECT_FALSE(signBit(fz));
    EXPECT_FALSE(signBit(fp));
    EXPECT_TRUE(signBit(dm));
    EXPECT_FALSE(signBit(dz));
    EXPECT_FALSE(signBit(dp));
}

// Note: signBit of a negator<float/double> must return the *same* result
// as for the underlying float/double — the negator changes the sign
// interpretation, not the stored bit that signBit inspects.
TEST(SimTKCommon_Scalar_SignBit, NegatorFloatAndDouble_SignBitSameAsUnderlying) {
    const float fm = -12398.34F;
    const float fz = 0.0F;
    const float fp = 4354.331F;
    const double dm = -234234.454;
    const double dz = 0.0;
    const double dp = 345345.2342;
    // -0 (may or may not be produced depending on compiler optimisation).
    const float mfz = -fz;
    const double mdz = -dz;

    const negator<float>& nfm = reinterpret_cast<const negator<float>&>(fm);
    const negator<float>& nfz = reinterpret_cast<const negator<float>&>(fz);
    const negator<float>& nfp = reinterpret_cast<const negator<float>&>(fp);
    const negator<float>& nmfz = reinterpret_cast<const negator<float>&>(mfz);
    const negator<double>& ndm = reinterpret_cast<const negator<double>&>(dm);
    const negator<double>& ndz = reinterpret_cast<const negator<double>&>(dz);
    const negator<double>& ndp = reinterpret_cast<const negator<double>&>(dp);
    const negator<double>& nmdz = reinterpret_cast<const negator<double>&>(mdz);

    EXPECT_TRUE(signBit(nfm));
    EXPECT_FALSE(signBit(nfz));
    EXPECT_FALSE(signBit(nfp));
    EXPECT_TRUE(signBit(ndm));
    EXPECT_FALSE(signBit(ndz));
    EXPECT_FALSE(signBit(ndp));
    // The sign bit of the negator of -0 matches that of -0 itself.
    EXPECT_EQ(signBit(nmfz), signBit(mfz));
    EXPECT_EQ(signBit(nmdz), signBit(mdz));
}

// +∞ has signBit = false; -∞ has signBit = true.
TEST(SimTKCommon_Scalar_SignBit, Infinity_SignBitReflectsSign) {
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    EXPECT_FALSE(signBit(fltInf));
    EXPECT_FALSE(signBit(dblInf));
    EXPECT_TRUE(signBit(mfltInf));
    EXPECT_TRUE(signBit(mdblInf));
}

// ============================================================
// sign
//
// sign(x) returns -1, 0, or +1 as the appropriate scalar type.
// For unsigned types, -1 is never returned.
// For negator<T>, sign returns the *opposite* of the underlying T
// because the stored value is interpreted as negated.
// ============================================================

// Unsigned integer types: sign is always 0 (zero) or 1 (positive).
TEST(SimTKCommon_Scalar_Sign, UnsignedTypes_SignNeverNegative) {
    const unsigned char ucm = 0xFFU;
    const unsigned char ucz = 0U;
    const unsigned char ucp = 27U;
    const unsigned short usm = 0xFFFFU;
    const unsigned short usz = 0U;
    const unsigned short usp = 2342U;
    const unsigned int uim = 0xFFFFFFFFU;
    const unsigned int uiz = 0U;
    const unsigned int uip = 2342344U;
    const auto ulm = static_cast<unsigned long>(-23423L);
    const unsigned long ulz = 0UL;
    const unsigned long ulp = 234234UL;
    const auto ullm = static_cast<unsigned long long>(-234234234LL);
    const unsigned long long ullz = 0ULL;
    const unsigned long long ullp = 234234234ULL;

    EXPECT_EQ(sign(ucm), 1);
    EXPECT_EQ(sign(ucz), 0);
    EXPECT_EQ(sign(ucp), 1);
    EXPECT_EQ(sign(usm), 1);
    EXPECT_EQ(sign(usz), 0);
    EXPECT_EQ(sign(usp), 1);
    EXPECT_EQ(sign(uim), 1);
    EXPECT_EQ(sign(uiz), 0);
    EXPECT_EQ(sign(uip), 1);
    EXPECT_EQ(sign(ulm), 1);
    EXPECT_EQ(sign(ulz), 0);
    EXPECT_EQ(sign(ulp), 1);
    EXPECT_EQ(sign(ullm), 1);
    EXPECT_EQ(sign(ullz), 0);
    EXPECT_EQ(sign(ullp), 1);
}

// Signed integer types: sign returns -1, 0, or +1 correctly.
// Note: sign(char) is not provided; use signed char.
TEST(SimTKCommon_Scalar_Sign, SignedIntegerTypes_SignReflectsSign) {
    const signed char cm = -23;
    const signed char cz = 0;
    const signed char cp = 99;
    const short sm = -1234;
    const short sz = 0;
    const short sp = 23423;
    const int im = -2342343;
    const int iz = 0;
    const int ip = 29472383;
    const long lm = -43488L;
    const long lz = 0L;
    const long lp = 3454545L;
    const long long llm = -2342342343433LL;
    const long long llz = 0LL;
    const long long llp = 874578478478574LL;

    EXPECT_EQ(sign(cm), -1);
    EXPECT_EQ(sign(cz), 0);
    EXPECT_EQ(sign(cp), 1);
    EXPECT_EQ(sign(sm), -1);
    EXPECT_EQ(sign(sz), 0);
    EXPECT_EQ(sign(sp), 1);
    EXPECT_EQ(sign(im), -1);
    EXPECT_EQ(sign(iz), 0);
    EXPECT_EQ(sign(ip), 1);
    EXPECT_EQ(sign(lm), -1);
    EXPECT_EQ(sign(lz), 0);
    EXPECT_EQ(sign(lp), 1);
    EXPECT_EQ(sign(llm), -1);
    EXPECT_EQ(sign(llz), 0);
    EXPECT_EQ(sign(llp), 1);
}

// Float and double: sign returns -1, 0, or +1 based on value.
TEST(SimTKCommon_Scalar_Sign, FloatAndDouble_SignReflectsSign) {
    const float fm = -12398.34F;
    const float fz = 0.0F;
    const float fp = 4354.331F;
    const double dm = -234234.454;
    const double dz = 0.0;
    const double dp = 345345.2342;

    EXPECT_EQ(sign(fm), -1);
    EXPECT_EQ(sign(fz), 0);
    EXPECT_EQ(sign(fp), 1);
    EXPECT_EQ(sign(dm), -1);
    EXPECT_EQ(sign(dz), 0);
    EXPECT_EQ(sign(dp), 1);
}

// Negative zero (-0.0) has sign = 0 regardless of the sign bit.
TEST(SimTKCommon_Scalar_Sign, NegativeZero_SignIsZero) {
    const float fz = 0.0F;
    const double dz = 0.0;
    // -0 (may or may not be produced depending on compiler).
    const float mfz = -fz;
    const double mdz = -dz;

    EXPECT_EQ(sign(mfz), 0);
    EXPECT_EQ(sign(mdz), 0);
}

// Note: sign of negator<float/double> must be the *opposite* of the
// sign of the underlying float/double, because the negator's value is
// the negation of the stored bits.
TEST(SimTKCommon_Scalar_Sign, NegatorFloatAndDouble_SignIsOpposite) {
    const float fm = -12398.34F;
    const float fz = 0.0F;
    const float fp = 4354.331F;
    const double dm = -234234.454;
    const double dz = 0.0;
    const double dp = 345345.2342;
    const float mfz = -fz;
    const double mdz = -dz;

    const auto& nfm = reinterpret_cast<const negator<float>&>(fm);
    const auto& nfz = reinterpret_cast<const negator<float>&>(fz);
    const auto& nfp = reinterpret_cast<const negator<float>&>(fp);
    const auto& nmfz = reinterpret_cast<const negator<float>&>(mfz);
    const auto& ndm = reinterpret_cast<const negator<double>&>(dm);
    const auto& ndz = reinterpret_cast<const negator<double>&>(dz);
    const auto& ndp = reinterpret_cast<const negator<double>&>(dp);
    const auto& nmdz = reinterpret_cast<const negator<double>&>(mdz);

    // The sign of a negator is opposite to the sign of the stored value.
    EXPECT_EQ(sign(nfm), 1);
    EXPECT_EQ(sign(nfz), 0);
    EXPECT_EQ(sign(nfp), -1);
    EXPECT_EQ(sign(ndm), 1);
    EXPECT_EQ(sign(ndz), 0);
    EXPECT_EQ(sign(ndp), -1);
    // -0 sign is 0 regardless of negator.
    EXPECT_EQ(sign(nmfz), 0);
    EXPECT_EQ(sign(nmdz), 0);
}

// +∞ → sign = +1; -∞ → sign = -1.  negator(+∞) evaluates as -∞ → sign = -1.
TEST(SimTKCommon_Scalar_Sign, Infinity_SignReflectsSign) {
    const float fltInf = NTraits<float>::getInfinity();
    const double dblInf = NTraits<double>::getInfinity();
    const float mfltInf = -fltInf;
    const double mdblInf = -dblInf;

    const auto& nfltInf = reinterpret_cast<const negator<float>&>(fltInf);
    const auto& ndblInf = reinterpret_cast<const negator<double>&>(dblInf);

    EXPECT_EQ(sign(fltInf), 1);
    EXPECT_EQ(sign(dblInf), 1);
    EXPECT_EQ(sign(mfltInf), -1);
    EXPECT_EQ(sign(mdblInf), -1);
    // The negator of +inf stores +inf bits but evaluates as -inf.
    EXPECT_EQ(sign(nfltInf), -1);
    EXPECT_EQ(sign(ndblInf), -1);
}

// ============================================================
// square / cube
//
// square(x) = x*x, cube(x) = x*x*x.
// For negator<T>: square(neg) = square(val) since (-x)²  = x²;
//                 cube(neg)   = -cube(val)  since (-x)^3  = -x³.
// ============================================================

// Basic float square and cube.
TEST(SimTKCommon_Scalar_SquareAndCube, ScalarFloat_SquareAndCube) {
    const float fval = -23.33F;

    SIMTK_EXPECT_EQ(square(fval), (fval * fval));
    SIMTK_EXPECT_EQ(cube(fval), (fval * fval * fval));
}

// Basic double square and cube.
TEST(SimTKCommon_Scalar_SquareAndCube, ScalarDouble_SquareAndCube) {
    const double dval = -234443.441;

    SIMTK_EXPECT_EQ(square(dval), (dval * dval));
    SIMTK_EXPECT_EQ(cube(dval), (dval * dval * dval));
}

// For scalar negators: square cancels the sign; cube preserves the negation.
TEST(SimTKCommon_Scalar_SquareAndCube, NegatorScalar_SquareAndCube) {
    const float fval = -23.33F;
    const double dval = -234443.441;

    const auto& nfval = reinterpret_cast<const negator<float>&>(fval);
    const auto& ndval = reinterpret_cast<const negator<double>&>(dval);

    SIMTK_EXPECT_EQ(square(nfval), (nfval * nfval));
    SIMTK_EXPECT_EQ(square(nfval), (fval * fval));
    SIMTK_EXPECT_EQ(square(ndval), (ndval * ndval));
    SIMTK_EXPECT_EQ(square(ndval), (dval * dval));

    SIMTK_EXPECT_EQ(cube(nfval), (nfval * nfval * nfval));
    SIMTK_EXPECT_EQ(cube(nfval), (-fval * fval * fval));
    SIMTK_EXPECT_EQ(cube(ndval), (ndval * ndval * ndval));
    SIMTK_EXPECT_EQ(cube(ndval), (-dval * dval * dval));
}

// Verify that conjugate<T> obeys the same square/cube rules as
// complex<T>.  A manual complex is constructed with the same bits
// to cross-check (sign change only → should be exact).
TEST(SimTKCommon_Scalar_SquareAndCube, ComplexAndConjugate_SquareAndCube) {
    const std::complex<float> fc(-234.343F, 45345e7F);
    const std::complex<double> dc(-234.343, 45345e7);
    const conjugate<float> fcj(-19.1e3F, -454.234F);
    const conjugate<double> dcj(-19.1e3, -454.234);

    // Manual conjugates constructed to match fcj/dcj bit-for-bit.
    const std::complex<float> fcmj(fcj.real(), fcj.imag());
    const std::complex<double> dcmj(dcj.real(), dcj.imag());

    // Sign-change only; equality should be exact.
    EXPECT_TRUE(fcj == fcmj);
    EXPECT_TRUE(dcj == dcmj);

    SIMTK_EXPECT_EQ((fcj * fcj), (fcmj * fcmj));
    SIMTK_EXPECT_EQ((dcj * dcj), (dcmj * dcmj));
    SIMTK_EXPECT_EQ((fcj * fcj * fcj), (fcmj * fcmj * fcmj));
    SIMTK_EXPECT_EQ((dcj * dcj * dcj), (dcmj * dcmj * dcmj));

    SIMTK_EXPECT_EQ(square(fc), (fc * fc));
    SIMTK_EXPECT_EQ(cube(fc), (fc * fc * fc));
    SIMTK_EXPECT_EQ(square(dc), (dc * dc));
    SIMTK_EXPECT_EQ(cube(dc), (dc * dc * dc));
    SIMTK_EXPECT_EQ(square(fcj), (fcj * fcj));
    SIMTK_EXPECT_EQ(cube(fcj), (fcj * fcj * fcj));
    SIMTK_EXPECT_EQ(square(dcj), (dcj * dcj));
    SIMTK_EXPECT_EQ(cube(dcj), (dcj * dcj * dcj));
}

// Tests involving negators of complex<T> and conjugate<T>:
// square cancels the sign; cube negates the result.
TEST(SimTKCommon_Scalar_SquareAndCube, NegatorOfComplexAndConjugate_SquareAndCube) {
    std::complex<float> fc(-234.343F, 45345e7F);
    std::complex<double> dc(-234.343, 45345e7);
    conjugate<float> fcj(-19.1e3F, -454.234F);
    conjugate<double> dcj(-19.1e3, -454.234);

    auto& nfc = reinterpret_cast<negator<std::complex<float>>&>(fc);
    auto& ndc = reinterpret_cast<negator<std::complex<double>>&>(dc);
    auto& nfcj = reinterpret_cast<negator<conjugate<float>>&>(fcj);
    auto& ndcj = reinterpret_cast<negator<conjugate<double>>&>(dcj);

    // Change of sign should be exact.
    EXPECT_TRUE(nfc == -fc);
    EXPECT_TRUE(ndc == -dc);
    EXPECT_TRUE(nfcj == -fcj);
    EXPECT_TRUE(ndcj == -dcj);

    // square: (-x)^2 = x^2, so the sign cancels.
    SIMTK_EXPECT_EQ(square(nfc), (nfc * nfc));
    SIMTK_EXPECT_EQ(square(nfc), (fc * fc));
    SIMTK_EXPECT_EQ(square(ndc), (ndc * ndc));
    SIMTK_EXPECT_EQ(square(ndc), (dc * dc));
    SIMTK_EXPECT_EQ(square(nfcj), (nfcj * nfcj));
    SIMTK_EXPECT_EQ(square(nfcj), (fcj * fcj));
    SIMTK_EXPECT_EQ(square(ndcj), (ndcj * ndcj));
    SIMTK_EXPECT_EQ(square(ndcj), (dcj * dcj));

    // cube: (-x)^3 = -x^3, so the result is negated.
    SIMTK_EXPECT_EQ(cube(nfc), (nfc * nfc * nfc));
    SIMTK_EXPECT_EQ(cube(nfc), (-fc * fc * fc));
    SIMTK_EXPECT_EQ(cube(ndc), (ndc * ndc * ndc));
    SIMTK_EXPECT_EQ(cube(ndc), (-dc * dc * dc));
    SIMTK_EXPECT_EQ(cube(nfcj), (nfcj * nfcj * nfcj));
    SIMTK_EXPECT_EQ(cube(nfcj), (-fcj * fcj * fcj));
    SIMTK_EXPECT_EQ(cube(ndcj), (ndcj * ndcj * ndcj));
    SIMTK_EXPECT_EQ(cube(ndcj), (-dcj * dcj * dcj));
}

// ============================================================
// isNumericallyEqual
// ============================================================

// float: value equal to itself, and to a value within default tolerance.
// A value outside default tolerance must be not-equal.
TEST(SimTKCommon_Scalar_IsNumericallyEqual, Float_WithinDefaultTolerance) {
    // f  — base value
    // fn — differs by 1e-5 (outside float default tolerance ~1e-6)
    // fe — differs by 1e-9 (inside float default tolerance)
    const float f = 1.234F;
    const float fn = 1.234F + 1e-5F;
    const float fe = 1.234F + 1e-9F;

    EXPECT_TRUE(isNumericallyEqual(f, f));
    EXPECT_TRUE(isNumericallyEqual(f, fe));  // within default tol
    EXPECT_FALSE(isNumericallyEqual(f, fn)); // outside default tol
}

// Providing an explicit tolerance overrides the default.
TEST(SimTKCommon_Scalar_IsNumericallyEqual, Float_WithCustomTolerance) {
    const float f = 1.234F;
    const float fn = 1.234F + 1e-5F;

    EXPECT_TRUE(isNumericallyEqual(f, fn, 1e-4F));  // loose tol → equal
    EXPECT_FALSE(isNumericallyEqual(f, fn, 1e-6F)); // tight tol → not equal
}

// CNT<float>::isNumericallyEqual must behave identically to the free function.
TEST(SimTKCommon_Scalar_IsNumericallyEqual, CNTFloat_WithinTolerance) {
    const float f = 1.234F;
    const float fn = 1.234F + 1e-5F;
    const float fe = 1.234F + 1e-9F;

    EXPECT_TRUE(CNT<float>::isNumericallyEqual(f, f));
    EXPECT_TRUE(CNT<float>::isNumericallyEqual(f, fe));
    EXPECT_FALSE(CNT<float>::isNumericallyEqual(f, fn));
    EXPECT_TRUE(CNT<float>::isNumericallyEqual(f, fn, 1e-4F));
    EXPECT_FALSE(CNT<float>::isNumericallyEqual(f, fn, 1e-6F));
}

// negator<float> stores f but evaluates as -f.
// isNumericallyEqual(nf, -f) must be true; isNumericallyEqual(nf, f) false.
TEST(SimTKCommon_Scalar_IsNumericallyEqual, NegatorFloat_EqualToNegatedValue) {
    const float f = 1.234F;

    const negator<float>& nf = negator<float>::recast(f);

    EXPECT_TRUE(nf.isNumericallyEqual(nf)); // self-equal
    EXPECT_TRUE(nf.isNumericallyEqual(-f)); // equal to negation of f
    EXPECT_FALSE(nf.isNumericallyEqual(f)); // not equal to f itself
}

// float compared against an integer: 1000*f ≈ 1234 within float tolerance.
TEST(SimTKCommon_Scalar_IsNumericallyEqual, FloatWithIntegerComparison) {
    const float f = 1.234F;
    const float fn = 1.234F + 1e-5F;
    const float fe = 1.234F + 1e-9F;

    EXPECT_TRUE(isNumericallyEqual(1000.0F * f, 1234));
    EXPECT_TRUE(isNumericallyEqual(1234, 1000.0F * f));
    EXPECT_TRUE(isNumericallyEqual(1000.0F * fe, 1234));
    EXPECT_TRUE(isNumericallyEqual(1234, 1000.0F * fe));
    EXPECT_FALSE(isNumericallyEqual(1000.0F * fn, 1234));
    EXPECT_FALSE(isNumericallyEqual(1234, 1000.0F * fn));
}

// double: value equal to itself and to values within double default tolerance.
TEST(SimTKCommon_Scalar_IsNumericallyEqual, Double_WithinDefaultTolerance) {
    // d  — base value
    // dn — differs by 1e-12 (outside double default tolerance ~1e-14)
    // de — differs by 1e-15 (inside double default tolerance)
    const double d = 1.234;
    const double dn = 1.234 + 1e-12;
    const double de = 1.234 + 1e-15;

    EXPECT_TRUE(isNumericallyEqual(d, d));
    EXPECT_TRUE(isNumericallyEqual(d, de));
    EXPECT_FALSE(isNumericallyEqual(d, dn));
}

// double compared against an integer: 1000*d ≈ 1234 within double tolerance.
TEST(SimTKCommon_Scalar_IsNumericallyEqual, DoubleWithIntegerComparison) {
    const double d = 1.234;
    const double dn = 1.234 + 1e-12;
    const double de = 1.234 + 1e-15;

    EXPECT_TRUE(isNumericallyEqual(1000.0 * d, 1234));
    EXPECT_TRUE(isNumericallyEqual(1234, 1000.0 * d));
    EXPECT_TRUE(isNumericallyEqual(1000.0 * de, 1234));
    EXPECT_TRUE(isNumericallyEqual(1234, 1000.0 * de));
    EXPECT_FALSE(isNumericallyEqual(1000.0 * dn, 1234));
    EXPECT_FALSE(isNumericallyEqual(1234, 1000.0 * dn));
}

// Mixed float/double: when one operand is float the comparison uses the
// (wider) float tolerance, not the tighter double tolerance.
// isNumericallyEqual(fe, de) must be true (float tol is sufficient), but
// isNumericallyEqual((double)fe, de) must be false (double tol is too tight).
TEST(SimTKCommon_Scalar_IsNumericallyEqual, MixedFloatDouble_UsesFloatTolerance) {
    const float fe = 1.234F + 1e-9F;
    const double de = 1.234 + 1e-15;

    // Mixed — uses float (wider) tolerance.
    EXPECT_TRUE(isNumericallyEqual(fe, de));
    // Both double — uses double (tighter) tolerance; fe cast to double
    // diverges more than double tol from de.
    EXPECT_FALSE(isNumericallyEqual(static_cast<double>(fe), de));
}

// ============================================================
// clamp / clampInPlace
//
// clamp(low, value, high) returns:
//   low   if value < low
//   high  if value > high
//   value otherwise
// clampInPlace modifies the variable in place and also returns the
// new value.
// ============================================================

// Integer clamping: exact boundary values and interior values.
TEST(SimTKCommon_Scalar_Clamp, Int_Clamp) {
    const int i4 = 4;

    EXPECT_EQ(clamp(4, i4, 4), 4);    // exactly on both bounds
    EXPECT_EQ(clamp(0, i4, 9), 4);    // below low → low
    EXPECT_EQ(clamp(5, i4, 9), 5);    // interior
    EXPECT_EQ(clamp(-7, i4, -5), -5); // below low → high? value=4 > high=-5 → high
}

// Double clamping.
TEST(SimTKCommon_Scalar_Clamp, Double_Clamp) {
    const double d325 = 3.25;

    EXPECT_EQ(clamp(3.25, d325, 3.25), 3.25);
    EXPECT_EQ(clamp(0., d325, 9.), 3.25);
    EXPECT_EQ(clamp(5., d325, 9.), 5.0);
    EXPECT_EQ(clamp(-7., d325, -5.), -5.0);
}

// Float clamping (note: fn325 is negative).
TEST(SimTKCommon_Scalar_Clamp, Float_Clamp) {
    const float fn325 = -3.25F;

    EXPECT_EQ(clamp(-3.25F, fn325, -3.25F), -3.25F);
    EXPECT_EQ(clamp(-9.F, fn325, 0.F), -3.25F);
    EXPECT_EQ(clamp(-9.F, fn325, -5.F), -5.0F);
    EXPECT_EQ(clamp(5.F, fn325, 7.F), 5.0F);
}

// Integer bounds supplied together with a double value.
TEST(SimTKCommon_Scalar_Clamp, IntegerBoundsOnDouble) {
    const double d325 = 3.25;

    EXPECT_EQ(clamp(0, d325, 9), 3.25);
    EXPECT_EQ(clamp(5, d325, 9), 5.0);
    EXPECT_EQ(clamp(-7, d325, -5), -5.0);
}

// Integer bounds supplied together with a float value.
TEST(SimTKCommon_Scalar_Clamp, IntegerBoundsOnFloat) {
    const float fn325 = -3.25F;

    EXPECT_EQ(clamp(-9, fn325, 0), -3.25F);
    EXPECT_EQ(clamp(-9, fn325, -5), -5.0F);
    EXPECT_EQ(clamp(5, fn325, 7), 5.0F);
}

// Mixed-type bounds: one bound is double / float, the other is int.
TEST(SimTKCommon_Scalar_Clamp, MixedTypeBounds) {
    const double d325 = 3.25;
    const float fn325 = -3.25F;

    // double value, mixed bounds
    EXPECT_EQ(clamp(0., d325, 9), 3.25);
    EXPECT_EQ(clamp(5., d325, 9), 5.0);
    EXPECT_EQ(clamp(-7., d325, -5), -5.0);

    EXPECT_EQ(clamp(0, d325, 9.), 3.25);
    EXPECT_EQ(clamp(5, d325, 9.), 5.0);
    EXPECT_EQ(clamp(-7, d325, -5.), -5.0);

    // float value, mixed bounds
    EXPECT_EQ(clamp(-9.F, fn325, 0), -3.25F);
    EXPECT_EQ(clamp(-9.F, fn325, -5), -5.0F);
    EXPECT_EQ(clamp(5.F, fn325, 7), 5.0F);

    EXPECT_EQ(clamp(-9, fn325, 0.F), -3.25F);
    EXPECT_EQ(clamp(-9, fn325, -5.F), -5.0F);
    EXPECT_EQ(clamp(5, fn325, 7.F), 5.0F);
}

// clampInPlace mutates its argument, returns the clamped value, and the
// variable holds the new value afterwards.
TEST(SimTKCommon_Scalar_Clamp, ClampInPlace) {
    int i = 4;
    double d = 3.25;
    float f = -3.25F;

    // int: 4 is above high=3 → clamped to 3
    EXPECT_EQ(clampInPlace(-2, i, 3), 3);
    EXPECT_EQ(i, 3);

    // double: 3.25 is above high=3.0 → clamped to 3
    EXPECT_EQ(clampInPlace(-2., d, 3.), 3.0);
    EXPECT_EQ(d, 3.0);

    // float: -3.25 is below low=-2 → clamped to -2
    EXPECT_EQ(clampInPlace(-2, f, 3), -2.0F);
    EXPECT_EQ(f, -2.0F);
}

// Char and related small-integer types.
TEST(SimTKCommon_Scalar_Clamp, CharTypes) {
    const char c = 'j';
    const unsigned char uc = 3U;
    const signed char sc = -2;

    EXPECT_EQ(clamp('a', c, 'e'), 'e');
    EXPECT_EQ(clamp('a', c, 'z'), 'j');
    EXPECT_EQ(clamp(static_cast<unsigned char>(4), uc, static_cast<unsigned char>(5)),
              static_cast<unsigned char>(4));
    EXPECT_EQ(clamp(static_cast<signed char>(-7), sc, static_cast<signed char>(-1)),
              static_cast<signed char>(-2));
}

// short and large unsigned int.
TEST(SimTKCommon_Scalar_Clamp, ShortAndLargeUnsignedInt) {
    const short s = -32000;
    const unsigned short us = 17U;
    const unsigned int ui = 4023456789U;

    EXPECT_EQ(clamp(static_cast<short>(-29000), s, static_cast<short>(400)), -29000);
    EXPECT_EQ(clamp(static_cast<unsigned short>(4), us, static_cast<unsigned short>(15)), 15U);
    EXPECT_EQ(clamp(100000000U, ui, 4010000000U), 4010000000U);
}

// long, unsigned long, long long, unsigned long long.
TEST(SimTKCommon_Scalar_Clamp, LongAndLongLongTypes) {
    const long l = -234234L;
    const unsigned long ul = 293493849UL;
    const long long ll = -123456789123LL;
    const unsigned long long ull = 123456789123ULL;

    EXPECT_EQ(clamp(-1000000L, l, -200000L), -234234L);
    EXPECT_EQ(clamp(1000000UL, ul, 4000000000UL), 293493849UL);
    // value=-123456789123, low=-100000000000 > value → low wins
    EXPECT_EQ(clamp(-100000000000LL, ll, 27LL), -100000000000LL);
    // value=-123456789123, low=-1000000000000 < value → value in range
    EXPECT_EQ(clamp(-1000000000000LL, ll, 27LL), -123456789123LL);
    // ull is unused as a clamp target, use it for a basic in-range test
    (void)ull;
}

// ============================================================
// stepUp / stepDown and their derivatives
//
// stepUp(x):   smooth S-curve from 0 (at x=0) to 1 (at x=1)
// stepDown(x): smooth S-curve from 1 (at x=0) to 0 (at x=1)
// Both are antisymmetric around x=0.5; the midpoint is exactly 0.5.
//
// dstepUp  / dstepDown  — first derivative  (zero at boundaries)
// d2stepUp / d2stepDown — second derivative (zero at boundaries)
// d3stepUp / d3stepDown — third derivative  (non-zero in interior)
//
// stepAny(y0, yr, x0, ooxr, x): generalised step function
//   y goes from y0 to y0+yr as x goes from x0 to x0+xr,
//   where ooxr = 1/xr.
// ============================================================

// Boundary and midpoint values for double stepUp / stepDown.
TEST(SimTKCommon_Scalar_ScalarStep, StepUpAndDown_Double_BoundaryValues) {
    EXPECT_EQ(stepUp(0.), 0.0);
    EXPECT_EQ(stepUp(0.5), 0.5);
    EXPECT_EQ(stepUp(1.), 1.0);
    EXPECT_EQ(stepDown(0.), 1.0);
    EXPECT_EQ(stepDown(0.5), 0.5);
    EXPECT_EQ(stepDown(1.), 0.0);
}

// stepUp is monotonically increasing; stepDown decreasing.
TEST(SimTKCommon_Scalar_ScalarStep, StepUpAndDown_Double_MonotonicInInterior) {
    // stepUp(.3) should be strictly between 0 and 0.5.
    EXPECT_LT(0.0, stepUp(0.3));
    EXPECT_LT(stepUp(0.3), 0.5);
    // stepUp(.7) should be strictly between 0.5 and 1.
    EXPECT_LT(0.5, stepUp(0.7));
    EXPECT_LT(stepUp(0.7), 1.0);

    // stepDown(.3) should be strictly between 0.5 and 1.
    EXPECT_LT(0.5, stepDown(0.3));
    EXPECT_LT(stepDown(0.3), 1.0);
    // stepDown(.7) should be strictly between 0 and 0.5.
    EXPECT_LT(0.0, stepDown(0.7));
    EXPECT_LT(stepDown(0.7), 0.5);
}

// Boundary and midpoint values for float stepUp / stepDown.
TEST(SimTKCommon_Scalar_ScalarStep, StepUpAndDown_Float_BoundaryValues) {
    EXPECT_EQ(stepUp(0.F), 0.0F);
    EXPECT_EQ(stepUp(0.5F), 0.5F);
    EXPECT_EQ(stepUp(1.F), 1.0F);
    EXPECT_EQ(stepDown(0.F), 1.0F);
    EXPECT_EQ(stepDown(0.5F), 0.5F);
    EXPECT_EQ(stepDown(1.F), 0.0F);
}

// Monotonicity for float.
TEST(SimTKCommon_Scalar_ScalarStep, StepUpAndDown_Float_MonotonicInInterior) {
    EXPECT_LT(0.0F, stepUp(0.3F));
    EXPECT_LT(stepUp(0.3F), 0.5F);
    EXPECT_LT(0.5F, stepUp(0.7F));
    EXPECT_LT(stepUp(0.7F), 1.0F);

    EXPECT_LT(0.5F, stepDown(0.3F));
    EXPECT_LT(stepDown(0.3F), 1.0F);
    EXPECT_LT(0.0F, stepDown(0.7F));
    EXPECT_LT(stepDown(0.7F), 0.5F);
}

// int arguments are treated as double (only endpoints are meaningful).
TEST(SimTKCommon_Scalar_ScalarStep, StepUpAndDown_Int_BoundaryValues) {
    EXPECT_EQ(stepUp(0), 0.0);
    EXPECT_EQ(stepUp(1), 1.0);
    EXPECT_EQ(stepDown(0), 1.0);
    EXPECT_EQ(stepDown(1), 0.0);
}

// First and second derivatives are zero at the boundary (double).
TEST(SimTKCommon_Scalar_ScalarStep, StepDerivatives_Double_BoundaryValues) {
    EXPECT_EQ(dstepUp(0.), 0.0);
    EXPECT_EQ(dstepUp(1.), 0.0);
    EXPECT_GT(dstepUp(0.5), 0.0); // strictly positive in interior
    EXPECT_EQ(dstepDown(0.), 0.0);
    EXPECT_EQ(dstepDown(1.), 0.0);
    EXPECT_LT(dstepDown(0.5), 0.0); // strictly negative in interior

    EXPECT_EQ(d2stepUp(0.), 0.0);
    EXPECT_EQ(d2stepUp(1.), 0.0);
    EXPECT_EQ(d2stepDown(0.), 0.0);
    EXPECT_EQ(d2stepDown(1.), 0.0);
}

// First and second derivatives are zero at the boundary (float).
TEST(SimTKCommon_Scalar_ScalarStep, StepDerivatives_Float_BoundaryValues) {
    EXPECT_EQ(dstepUp(0.F), 0.0F);
    EXPECT_EQ(dstepUp(1.F), 0.0F);
    EXPECT_GT(dstepUp(0.5F), 0.0F);
    EXPECT_EQ(dstepDown(0.F), 0.0F);
    EXPECT_EQ(dstepDown(1.F), 0.0F);
    EXPECT_LT(dstepDown(0.5F), 0.0F);

    EXPECT_EQ(d2stepUp(0.F), 0.0F);
    EXPECT_EQ(d2stepUp(1.F), 0.0F);
    EXPECT_EQ(d2stepDown(0.F), 0.0F);
    EXPECT_EQ(d2stepDown(1.F), 0.0F);
}

// Central finite-difference validation of the analytic derivatives (double).
// Central difference estimates should give around 10 decimal places in double.
TEST(SimTKCommon_Scalar_ScalarStep, StepDerivatives_Double_FiniteDifference) {
    constexpr double h = 1e-6;
    constexpr double tol = 1e-8;

    const double dupEst = (stepUp(0.799 + h) - stepUp(0.799 - h)) / (2.0 * h);
    const double ddnEst = (stepDown(0.799 + h) - stepDown(0.799 - h)) / (2.0 * h);
    const double d2upEst = (dstepUp(0.723 + h) - dstepUp(0.723 - h)) / (2.0 * h);
    const double d2dnEst = (dstepDown(0.723 + h) - dstepDown(0.723 - h)) / (2.0 * h);
    const double d3upEst = (d2stepUp(0.123 + h) - d2stepUp(0.123 - h)) / (2.0 * h);
    const double d3dnEst = (d2stepDown(0.123 + h) - d2stepDown(0.123 - h)) / (2.0 * h);

    EXPECT_NEAR(dstepUp(0.799), dupEst, tol);
    EXPECT_NEAR(dstepDown(0.799), ddnEst, tol);
    EXPECT_NEAR(d2stepUp(0.723), d2upEst, tol);
    EXPECT_NEAR(d2stepDown(0.723), d2dnEst, tol);
    EXPECT_NEAR(d3stepUp(0.123), d3upEst, tol);
    EXPECT_NEAR(d3stepDown(0.123), d3dnEst, tol);
}

// Central finite-difference validation of analytic derivatives (float).
// Float only gives ~4 decimal places; step h is coarser, tol is wider.
TEST(SimTKCommon_Scalar_ScalarStep, StepDerivatives_Float_FiniteDifference) {
    constexpr float h = 1e-3F;
    constexpr float tol = 1e-3F;

    const float fdupEst = (stepUp(0.699F + h) - stepUp(0.699F - h)) / (2.0F * h);
    const float fddnEst = (stepDown(0.699F + h) - stepDown(0.699F - h)) / (2.0F * h);
    const float fd2upEst = (dstepUp(0.623F + h) - dstepUp(0.623F - h)) / (2.0F * h);
    const float fd2dnEst = (dstepDown(0.623F + h) - dstepDown(0.623F - h)) / (2.0F * h);
    const float fd3upEst = (d2stepUp(0.211F + h) - d2stepUp(0.211F - h)) / (2.0F * h);
    const float fd3dnEst = (d2stepDown(0.211F + h) - d2stepDown(0.211F - h)) / (2.0F * h);

    EXPECT_NEAR(dstepUp(0.699F), fdupEst, tol);
    EXPECT_NEAR(dstepDown(0.699F), fddnEst, tol);
    EXPECT_NEAR(d2stepUp(0.623F), fd2upEst, tol);
    EXPECT_NEAR(d2stepDown(0.623F), fd2dnEst, tol);
    EXPECT_NEAR(d3stepUp(0.211F), fd3upEst, tol);
    EXPECT_NEAR(d3stepDown(0.211F), fd3dnEst, tol);
}

// stepAny: generalised step with arbitrary output range and x-domain (double).
// y = stepAny(y0, yr, x0, 1/xr, x)
// y goes from y0 to y0+yr as x goes from x0 to x0+xr.
TEST(SimTKCommon_Scalar_ScalarStep, StepAny_Double_BoundaryValues) {
    // Simple case: y goes from -1 to 1 as x goes from 0 to 1.
    EXPECT_EQ(stepAny(-1, 2, 0, 1, 0.), -1.0);
    EXPECT_EQ(stepAny(-1, 2, 0, 1, 0.5), 0.0);
    EXPECT_EQ(stepAny(-1, 2, 0, 1, 1.), 1.0);

    // Arbitrary range: y goes from -7 down to -14 as x goes from -3.1 to 429.3.
    const double x0 = -3.1;
    const double x1 = 429.3;
    const double y0 = -7.0;
    const double y1 = -14.0;
    const double xr = (x1 - x0);
    const double ooxr = 1.0 / xr;
    const double yr = (y1 - y0);

    SIMTK_EXPECT_EQ(stepAny(y0, yr, x0, ooxr, x0), y0);
    SIMTK_EXPECT_EQ(stepAny(y0, yr, x0, ooxr, x1), y1);
    SIMTK_EXPECT_EQ(stepAny(y0, yr, x0, ooxr, (x0 + (xr / 2.0))), (y0 + (yr / 2.0)));
}

// stepAny: boundary values for float.
TEST(SimTKCommon_Scalar_ScalarStep, StepAny_Float_BoundaryValues) {
    // Simple case (float).
    EXPECT_EQ(stepAny(-1, 2, 0, 1, 0.F), -1.0F);
    EXPECT_EQ(stepAny(-1, 2, 0, 1, 0.5F), 0.0F);
    EXPECT_EQ(stepAny(-1, 2, 0, 1, 1.F), 1.0F);

    const float fx0 = -3.1F;
    const float fx1 = 429.3F;
    const float fy0 = -7.0F;
    const float fy1 = -14.0F;
    const float fxr = (fx1 - fx0);
    const float fooxr = 1.0F / fxr;
    const float fyr = (fy1 - fy0);

    SIMTK_EXPECT_EQ(stepAny(fy0, fyr, fx0, fooxr, fx0), fy0);
    SIMTK_EXPECT_EQ(stepAny(fy0, fyr, fx0, fooxr, fx1), fy1);
    SIMTK_EXPECT_EQ(stepAny(fy0, fyr, fx0, fooxr, (fx0 + (fxr / 2.0F))), (fy0 + (fyr / 2.0F)));
}

// Central finite-difference validation of dstepAny / d2stepAny / d3stepAny
// (double).
TEST(SimTKCommon_Scalar_ScalarStep, StepAny_Double_Derivatives_FiniteDifference) {
    constexpr double h = 1e-6;
    constexpr double tol = 1e-8;

    const double x0 = -3.1;
    const double x1 = 429.3;
    const double y0 = -7.0;
    const double y1 = -14.0;
    const double xr = (x1 - x0);
    const double ooxr = 1.0 / xr;
    const double yr = (y1 - y0);

    const double danyEst =
        (stepAny(y0, yr, x0, ooxr, 0.799 + h) - stepAny(y0, yr, x0, ooxr, 0.799 - h)) / (2.0 * h);
    const double d2anyEst =
        (dstepAny(yr, x0, ooxr, 0.723 + h) - dstepAny(yr, x0, ooxr, 0.723 - h)) / (2.0 * h);
    const double d3anyEst =
        (d2stepAny(yr, x0, ooxr, 0.123 + h) - d2stepAny(yr, x0, ooxr, 0.123 - h)) / (2.0 * h);

    EXPECT_NEAR(dstepAny(yr, x0, ooxr, 0.799), danyEst, tol);
    EXPECT_NEAR(d2stepAny(yr, x0, ooxr, 0.723), d2anyEst, tol);
    EXPECT_NEAR(d3stepAny(yr, x0, ooxr, 0.123), d3anyEst, tol);
}

// Central finite-difference validation of dstepAny / d2stepAny / d3stepAny
// (float).
TEST(SimTKCommon_Scalar_ScalarStep, StepAny_Float_Derivatives_FiniteDifference) {
    constexpr float h = 1e-3F;
    constexpr float tol = 1e-3F;

    const float fx0 = -3.1F;
    const float fx1 = 429.3F;
    const float fy0 = -7.0F;
    const float fy1 = -14.0F;
    const float fxr = (fx1 - fx0);
    const float fooxr = 1.0F / fxr;
    const float fyr = (fy1 - fy0);

    const float fdanyEst =
        (stepAny(fy0, fyr, fx0, fooxr, 0.799F + h) - stepAny(fy0, fyr, fx0, fooxr, 0.799F - h)) / (2.0F * h);
    const float fd2anyEst =
        (dstepAny(fyr, fx0, fooxr, 0.723F + h) - dstepAny(fyr, fx0, fooxr, 0.723F - h)) / (2.0F * h);
    const float fd3anyEst =
        (d2stepAny(fyr, fx0, fooxr, 0.123F + h) - d2stepAny(fyr, fx0, fooxr, 0.123F - h)) / (2.0F * h);

    EXPECT_NEAR(dstepAny(fyr, fx0, fooxr, 0.799F), fdanyEst, tol);
    EXPECT_NEAR(d2stepAny(fyr, fx0, fooxr, 0.723F), fd2anyEst, tol);
    EXPECT_NEAR(d3stepAny(fyr, fx0, fooxr, 0.123F), fd3anyEst, tol);
}
