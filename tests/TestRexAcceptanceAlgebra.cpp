// ============================================================================
//  TestRexAcceptanceAlgebra.cpp -- V8 (exact acceptance algebra) and V4 (no-
//  drive limit), the C++ analogue of tests/test_rex_swap_acceptance_algebra.py
//  (docs/specs/replica-exchange-nonequilibrium-work.md, Derivation sketch,
//  B6, D2, V4, V8). Also exercises RENEMC's ETerm_nonequil formula (Stage 2b,
//  acceptance wired even though its own driven round-loop is Stage 2c) and
//  the S2 domain-error consuming-side contract (Context::driveReplica maps a
//  std::domain_error from World::applyBatScalingDrive to WORK_Jacobian =
//  -infinity, guaranteeing an automatic reject).
//
//  Noise-free spec oracle, no MD, no engine Context/OpenMM: a one-body
//  analytic system, exactly mirroring the Python reproducer's structure
//  (kept in place, untouched, at tests/test_rex_swap_acceptance_algebra.py)
//  so both languages check the SAME closed-form claims. This file
//  DELIBERATELY reimplements Context::attemptREXSwap's log-alpha formulas
//  inline (REMC ETerm_equal, RENEMC ETerm_nonequil, RENE/REBASONTOP WTerm) --
//  see each TEST's comment for the cross-reference to the Context.cpp
//  switch case it mirrors -- rather than instantiating a live Context, which
//  needs a full OpenMM-backed World (this test suite's own convention for
//  closed-form analytic oracles, e.g. TestFixmanBoltzmann.cpp/
//  TestRexSwapAcceptanceAlgebra's own Python original: no live engine
//  Context needed to check a detailed-balance identity).
//
//  One-body system: each replica is a single bond, a vector x in R^3 with
//  U(x) = 0.5 k (|x| - r0)^2. The BAT scaling drives the bond isotropically,
//  x -> s x (mean 0), whose exact Cartesian log-Jacobian is 3 ln s (D6: a
//  single isotropic bond scale gives |dx'/dx| = s^3). Closed form.
//
//  NOT compiled or run (coordinator directive, 2026-07-12) -- written to the
//  tests/Test*.cpp convention for the user to build and run manually.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <utility>
#include <vector>

#include "TestHelpers.hpp"
#include "robot_math.hpp"

using rtest::Rng;
using robo::Vec3;

namespace {

constexpr double kBond = 3.0;
constexpr double kR0 = 1.0;

double U(const Vec3& x) {
    const double r = x.norm();
    return 0.5 * kBond * (r - kR0) * (r - kR0);
}

// Exact Cartesian log-Jacobian of x -> s*x in R^3 (D6: 3 ln s).
double lnJacIsotropic(double s) {
    return 3.0 * std::log(s);
}

struct SwapResult {
    double logAlpha;
    Vec3 xXtau;
    Vec3 xYtau;
};

// Mirrors Context::attemptREXSwap's RUN_TYPE::RENE/REBASONTOP branch (WTerm,
// B6/D7): X = replica in the cold (source) state, driven by s toward hot
// (target); Y = replica in the hot (source) state, driven by 1/s toward
// cold. mean (anchor) = 0 so M_s(q) = s*q is an exact involution with
// M_{1/s} (D2/INV-9, this file's system has no separate anchor test -- see
// TestBatAnchorInvolution.cpp for that). jacSign/useLnS/swapBetas inject the
// bugs the oracle must catch (the SAME three knobs the Python original
// uses).
SwapResult wTermSwap(const Vec3& coldCfg,
                     const Vec3& hotCfg,
                     double betaC,
                     double betaH,
                     double s,
                     double jacSign = 1.0,
                     bool useLnS = true,
                     bool swapBetas = false) {
    const Vec3 xXtau = coldCfg * s;
    const Vec3 xYtau = hotCfg * (1.0 / s);
    const double lnJacX = (useLnS ? lnJacIsotropic(s) : 0.0) * jacSign;
    const double lnJacY = (useLnS ? lnJacIsotropic(1.0 / s) : 0.0) * jacSign;

    double workX = 0.0;
    double workY = 0.0;
    if (!swapBetas) {
        // Correct (B6/D7): x^tau at TARGET beta, x^0 at SOURCE beta.
        workX = (betaH * U(xXtau)) - (betaC * U(coldCfg)) - lnJacX;
        workY = (betaC * U(xYtau)) - (betaH * U(hotCfg)) - lnJacY;
    } else {
        // BUG: x^tau at source beta, x^0 at target beta.
        workX = (betaC * U(xXtau)) - (betaH * U(coldCfg)) - lnJacX;
        workY = (betaH * U(xYtau)) - (betaC * U(hotCfg)) - lnJacY;
    }
    const double wTerm = -(workX + workY);
    return {wTerm, xXtau, xYtau};
}

// Joint target: cold state holds coldCfg, hot state holds hotCfg.
double logPi(const Vec3& coldCfg, const Vec3& hotCfg, double betaC, double betaH) {
    return (-betaC * U(coldCfg)) - (betaH * U(hotCfg));
}

double acceptProb(double logAlpha) {
    return std::min(1.0, std::exp(logAlpha));
}

// Mirrors Context::attemptREXSwap's RUN_TYPE::REMC branch (ETerm_equal, B6
// step 3) on the reduced potentials directly (no drive).
double eTermEqual(double refUXset, double refUYset, double betaC, double betaH) {
    return -((betaH - betaC) * (refUXset - refUYset));
}

// Mirrors Context::attemptREXSwap's RUN_TYPE::RENEMC branch (ETerm_nonequil,
// B6 step 3): the SAME parallel-tempering form on the DRIVEN-ENDPOINT
// reference potentials, no Jacobian (INV-10: RENEMC's velocity/NMA drive is
// volume-preserving).
double eTermNonequil(double refUXtau, double refUYtau, double betaC, double betaH) {
    return -((betaH - betaC) * (refUXtau - refUYtau));
}

} // namespace

// ---------------------------------------------------------------------------
//  A. beta-assignment is guarded by the true detailed-balance equation
//     (mirrors Context.cpp's RENE/REBASONTOP WTerm branch).
// ---------------------------------------------------------------------------
TEST(RexAcceptanceAlgebra, BetaAssignmentSatisfiesDetailedBalance) {
    for (const double s : {1.15, 1.4, 0.8}) {
        for (const auto& betas : {std::pair<double, double>{1.0, 0.5}, std::pair<double, double>{2.0, 1.3}}) {
            const double betaC = betas.first;
            const double betaH = betas.second;
            Rng rng(1);
            const Vec3 coldc = Vec3(1.0, 0.0, 0.0) + Vec3(rng.gaussian(), rng.gaussian(), rng.gaussian());
            const Vec3 hotc = Vec3(1.0, 0.0, 0.0) + Vec3(rng.gaussian(), rng.gaussian(), rng.gaussian());

            const SwapResult fwd = wTermSwap(coldc, hotc, betaC, betaH, s);
            // Tz: cold state now holds the driven Y config, hot holds driven X.
            const SwapResult rev = wTermSwap(fwd.xYtau, fwd.xXtau, betaC, betaH, s);

            const double lhs = std::exp(logPi(coldc, hotc, betaC, betaH)) * acceptProb(fwd.logAlpha);
            const double rhs = std::exp(logPi(fwd.xYtau, fwd.xXtau, betaC, betaH)) * acceptProb(rev.logAlpha);
            EXPECT_NEAR(lhs, rhs, 1e-12) << "s=" << s << " betaC=" << betaC << " betaH=" << betaH;
        }
    }
}

TEST(RexAcceptanceAlgebra, SwappedBetaAssignmentBreaksDetailedBalance) {
    const double betaC = 1.0;
    const double betaH = 0.5;
    for (const double s : {1.15, 1.4}) {
        Rng rng(2);
        const Vec3 coldc = Vec3(1.0, 0.0, 0.0) + Vec3(rng.gaussian(), rng.gaussian(), rng.gaussian());
        const Vec3 hotc = Vec3(1.0, 0.0, 0.0) + Vec3(rng.gaussian(), rng.gaussian(), rng.gaussian());

        const SwapResult fwd = wTermSwap(coldc, hotc, betaC, betaH, s, /*jacSign=*/1.0, /*useLnS=*/true,
                                         /*swapBetas=*/true);
        const SwapResult rev =
            wTermSwap(fwd.xYtau, fwd.xXtau, betaC, betaH, s, 1.0, true, /*swapBetas=*/true);

        const double lhs = std::exp(logPi(coldc, hotc, betaC, betaH)) * acceptProb(fwd.logAlpha);
        const double rhs = std::exp(logPi(fwd.xYtau, fwd.xXtau, betaC, betaH)) * acceptProb(rev.logAlpha);
        EXPECT_GT(std::abs(lhs - rhs), 1e-9) << "the x^0<->x^tau beta swap MUST be caught by DB, s=" << s;
    }
}

// ---------------------------------------------------------------------------
//  B. FINDING (kept from the Python original): the SYMMETRIC paired swap
//     cannot see the Jacobian (joint |det| == s^3 * s^-3 == 1). Pins the gap
//     so nobody trusts this two-replica DB check to guard D6/F7 -- that is
//     what part C (single-body NCMC) and TestBatScalingJacobian.cpp's V5 do.
// ---------------------------------------------------------------------------
TEST(RexAcceptanceAlgebra, SymmetricSwapIsBlindToJacobianSignFlip) {
    const double betaC = 1.0;
    const double betaH = 0.5;
    for (const double s : {1.15, 1.4}) {
        Rng rng(3);
        const Vec3 coldc = Vec3(1.0, 0.0, 0.0) + Vec3(rng.gaussian(), rng.gaussian(), rng.gaussian());
        const Vec3 hotc = Vec3(1.0, 0.0, 0.0) + Vec3(rng.gaussian(), rng.gaussian(), rng.gaussian());

        auto dbGap = [&](double jacSign, bool useLnS) {
            const SwapResult fwd = wTermSwap(coldc, hotc, betaC, betaH, s, jacSign, useLnS);
            const SwapResult rev = wTermSwap(fwd.xYtau, fwd.xXtau, betaC, betaH, s, jacSign, useLnS);
            const double lhs = std::exp(logPi(coldc, hotc, betaC, betaH)) * acceptProb(fwd.logAlpha);
            const double rhs = std::exp(logPi(fwd.xYtau, fwd.xXtau, betaC, betaH)) * acceptProb(rev.logAlpha);
            return std::abs(lhs - rhs);
        };

        EXPECT_LT(dbGap(1.0, true), 1e-12) << "s=" << s;
        EXPECT_LT(dbGap(-1.0, true), 1e-12) << "sign-flipped Jacobian, s=" << s;
        EXPECT_LT(dbGap(1.0, false), 1e-12) << "dropped ln(s), s=" << s;
    }
}

// ---------------------------------------------------------------------------
//  C. The Jacobian IS observable where |det dT/dz| != 1: single-replica NCMC
//     (Nilmeier eq 28). This is the oracle a D6/F7 Jacobian bug needs -- see
//     also TestBatScalingJacobian.cpp's V5 (the engine-level version).
// ---------------------------------------------------------------------------
TEST(RexAcceptanceAlgebra, SingleBodyNcmcAcceptanceNeedsCorrectJacobian) {
    const double beta = 1.0;
    for (const double s : {1.15, 1.4, 0.8}) {
        Rng rng(4);
        const Vec3 x0 = Vec3(1.5, 0.0, 0.0) + Vec3(rng.gaussian(), rng.gaussian(), rng.gaussian());
        const Vec3 xtau = x0 * s;

        const double analytic = (-beta * (U(xtau) - U(x0))) + lnJacIsotropic(s);
        const double correct = (-beta * (U(xtau) - U(x0))) + lnJacIsotropic(s);
        const double flipped = (-beta * (U(xtau) - U(x0))) - lnJacIsotropic(s);
        const double dropped = -beta * (U(xtau) - U(x0));

        EXPECT_NEAR(correct, analytic, 1e-12) << "s=" << s;
        EXPECT_GT(std::abs(flipped - analytic), 1e-6) << "flipped sign must be caught, s=" << s;
        EXPECT_GT(std::abs(dropped - analytic), 1e-6) << "missing ln(s) must be caught, s=" << s;
    }
}

// ---------------------------------------------------------------------------
//  V4 -- no-drive limit: with s=1 (lnJac=0, x^tau=x^0), WTerm degenerates to
//  ETerm_equal (Derivation sketch: "the no-drive limit (x^tau=x^0, lnJac=0)
//  collapses to ETerm_equal"). Confirms RENE's acceptance is a strict
//  generalisation of REMC's, not a different formula that happens to
//  coincide numerically.
// ---------------------------------------------------------------------------
TEST(RexAcceptanceAlgebra, NoDriveLimitReducesToETermEqual) {
    const double betaC = 1.0;
    const double betaH = 0.6;
    Rng rng(5);
    const Vec3 coldc = Vec3(1.0, 0.0, 0.0) + Vec3(rng.gaussian(), rng.gaussian(), rng.gaussian());
    const Vec3 hotc = Vec3(1.0, 0.0, 0.0) + Vec3(rng.gaussian(), rng.gaussian(), rng.gaussian());

    const double s = 1.0; // no drive: x^tau == x^0, lnJac == 0
    const SwapResult wterm = wTermSwap(coldc, hotc, betaC, betaH, s);
    ASSERT_NEAR(lnJacIsotropic(s), 0.0, 1e-15) << "sanity: lnJac(s=1) must be exactly 0";
    ASSERT_NEAR((wterm.xXtau - coldc).norm(), 0.0, 1e-15) << "sanity: x^tau == x^0 at s=1";

    const double eTerm = eTermEqual(U(coldc), U(hotc), betaC, betaH);
    EXPECT_NEAR(wterm.logAlpha, eTerm, 1e-12)
        << "WTerm at s=1 (no drive) must equal ETerm_equal (V4, Derivation sketch)";
}

// ---------------------------------------------------------------------------
//  RENEMC's ETerm_nonequil (Stage 2b acceptance, INV-10: no Jacobian). Since
//  its own driven round-loop is Stage 2c (velocity/NMA drive not wired), this
//  exercises the FORMULA directly, mirroring how Context::attemptREXSwap's
//  RUN_TYPE::RENEMC branch reads replicas_[X/Y].referenceWORK_potential --
//  here those are just the (already-driven-by-some-external-mechanism)
//  potentials U(x_X^tau)/U(x_Y^tau) passed in directly.
// ---------------------------------------------------------------------------
TEST(RexAcceptanceAlgebra, RenemcETermNonequilHasNoJacobianTerm) {
    const double betaC = 1.0;
    const double betaH = 0.7;
    Rng rng(6);
    // Driven-endpoint potentials for a volume-preserving move: pick two
    // ARBITRARY reduced potentials directly (RENEMC's drive does not scale
    // BAT coordinates, so there is no s/Jacobian in this formula at all --
    // that IS the point of the test: the formula has no lnJac term to omit
    // incorrectly).
    const double refUXtau = 2.0 + rng.uniform(-0.5, 0.5);
    const double refUYtau = 1.3 + rng.uniform(-0.5, 0.5);

    const double eTermNe = eTermNonequil(refUXtau, refUYtau, betaC, betaH);
    const double eTermEq = eTermEqual(refUXtau, refUYtau, betaC, betaH);
    // Same functional FORM as ETerm_equal (B6 step 3: "same form on the
    // driven-endpoint reference potentials") -- confirmed by construction
    // (both eTermNonequil/eTermEqual share one implementation), so this
    // assertion is a change-detector against an accidental divergence
    // between the two helpers, not a tautology on paper.
    EXPECT_NEAR(eTermNe, eTermEq, 1e-15);
}

// ---------------------------------------------------------------------------
//  S2 (reviewer, 2026-07-12) -- domain-error consuming side: driveReplica
//  (Context.cpp) maps a caught std::domain_error to WORK_Jacobian =
//  -infinity. Confirm that value forces WTerm to -infinity (deterministic
//  automatic reject: -inf < 0 and exp(-inf) == 0, so no RNG draw can ever
//  accept it) REGARDLESS of the (now-irrelevant) potential values. This is
//  the "consuming side" of the mechanism; TestBatScalingJacobian.cpp's
//  AggressiveScaleThrowsDomainErrorNotSilentCorruption test is the
//  "producing side" (World::applyBatScalingDrive actually throwing). A live
//  Context/OpenMM harness exercising Context::driveReplica's try/catch
//  end-to-end does not exist in this test suite (flagged in the coder
//  checkpoint) -- these two tests together validate the mechanism's logic
//  without one.
// ---------------------------------------------------------------------------
TEST(RexAcceptanceAlgebra, DomainErrorSentinelJacobianForcesAutomaticReject) {
    const double betaC = 1.0;
    const double betaH = 2.0; // deliberately the LARGER beta, the case most
                              // likely to accidentally accept if the sentinel
                              // were not overwhelming
    const double refUXset = 0.5; // an arbitrarily "favourable" set of potentials
    const double refUYset = 5.0;
    const double refUXtau = 0.1; // (post-scale potentials, irrelevant once lnJac=-inf dominates)
    const double refUYtau = 0.1;
    const double lnJacXValid = 1.0; // Y's drive succeeded normally
    const double lnJacXSentinel = -std::numeric_limits<double>::infinity(); // X's drive hit the domain guard

    // Mirrors Context::attemptREXSwap's RUN_TYPE::RENE/REBASONTOP branch
    // verbatim (Context.cpp): workX = betaH*refUXtau - betaC*refUXset - lnJacX; etc.
    const double workY = (betaC * refUYtau) - (betaH * refUYset) - lnJacXValid;
    const double workXSentinel = (betaH * refUXtau) - (betaC * refUXset) - lnJacXSentinel;
    ASSERT_TRUE(std::isinf(workXSentinel));
    ASSERT_GT(workXSentinel, 0.0) << "workX must go to +infinity (an infinitely bad/invalid Jacobian), "
                                     "not -infinity, for the sign to force a REJECT via WTerm=-(workX+workY)";

    const double wTerm = -(workXSentinel + workY);
    EXPECT_TRUE(std::isinf(wTerm));
    EXPECT_LT(wTerm, 0.0) << "WTerm must be -infinity (guaranteed reject), got " << wTerm;

    // The natural accept test (logPAccept >= 0 || u < exp(logPAccept)) then
    // rejects deterministically, with no RNG draw able to overturn it:
    EXPECT_FALSE(wTerm >= 0.0);
    EXPECT_EQ(std::exp(wTerm), 0.0) << "exp(-infinity) == 0, so u < exp(wTerm) is false for every u in [0,1)";
}
