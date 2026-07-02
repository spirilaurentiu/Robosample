// ============================================================================
//  TestEnsembleValidation.cpp -- distribution-level canonical-ensemble checks
//  (Shirts 2013, "Simple quantitative tests to validate sampling from
//  thermodynamic ensembles", eqs. 19-20): KE and PE INDEPENDENTLY obey
//        P(E|β₂)/P(E|β₁) = const · exp(−(β₂−β₁)E).
//  These are stronger than the first-moment equipartition check (⟨2KE⟩/n=kT):
//  they validate the whole energy distribution.
//
//   * KEisMaxwellBoltzmann -- the momentum draw gives KE ~ Gamma(n_dof/2, kT):
//     density ∝ E^{n/2−1} e^{−E/kT}. A wrong DOF count or a broken sqrt(M⁻¹) draw
//     shows up as a distribution-shape mismatch, not just a mean shift.
//   * PEobeysEnsembleSlope -- sampling a bound system at two temperatures, the
//     log-ratio of PE histograms is linear in E with slope −(β₂−β₁).
// ============================================================================
#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <vector>

#include "AnalyticForceBridge.hpp"
#include "Constraints.hpp"
#include "HmcDriver.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "StatTest.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::attachAtoms;
using rtest::BodySpec;
using rtest::buildBentTorsionChain;
using rtest::buildForest;
using rtest::HmcDriver;
using rtest::Rng;
using rtest::stat::chiSquareCritical;
using rtest::stat::chiSquareStatistic;
using rtest::stat::expectedFromWeights;
using rtest::stat::Histogram;

namespace {

constexpr double kB = 0.0083144626;
constexpr double T0 = 300.0;
constexpr double kT300 = kB * T0;

bool slowEnabled() {
    return std::getenv("ROBOSAMPLE_SLOW_TESTS") != nullptr;
}

// Sub-sampled per-bin Gamma(shape, scaleKT) weight over a [0, hiBound] histogram
// (lo == 0 assumed, matching every Histogram this file builds for a KE domain).
// Bin-centre evaluation is biased for this curved density and inflates chi^2 at
// large N (KEisMaxwellBoltzmann's comment); sub-sampling each bin fixes that.
// Generalized over `shape` (vs. KEisMaxwellBoltzmann's fixed nDof/2) so T0.2 can
// evaluate the SAME histogram against the correct shape and the wrong-DOF
// controls without re-histogramming.
std::vector<double> gammaSubsampledWeights(const Histogram& h, double hiBound, double shape, double scaleKT) {
    const double binW = hiBound / h.nbins();
    std::vector<double> w(static_cast<std::size_t>(h.nbins()), 0.0);
    for (int b = 0; b < h.nbins(); ++b) {
        const double lo = h.center(b) - 0.5 * binW;
        const int sub = 16;
        double acc = 0.0;
        for (int k = 0; k < sub; ++k) {
            const double E = lo + (k + 0.5) * binW / sub;
            if (E > 0) {
                acc += std::pow(E, shape - 1.0) * std::exp(-E / scaleKT);
            }
        }
        w[static_cast<std::size_t>(b)] = acc / sub;
    }
    return w;
}

RobotModel freeBody() {
    BodySpec b;
    b.parent = 0;
    b.joint = JointType::Free;
    b.mass = Real(2.0);
    b.inertia_B = UnitInertia(Real(0.4), Real(0.55), Real(0.7));
    RobotModel m = buildForest({b});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.1, 0, 0), Vec3(0, 0.1, 0)}, {12.0, 1.0, 1.0});
    return m;
}

} // namespace

// ---------------------------------------------------------------------------
//  KE from the momentum draw is Gamma(n_dof/2, kT): density ∝ E^{n/2−1} e^{−E/kT}.
//  Free body: n_dof = 6, so density ∝ E² e^{−E/kT}.
// ---------------------------------------------------------------------------
TEST(EnsembleValidation, KEisMaxwellBoltzmann) {
    RobotModel m = freeBody();
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1;
    ConstraintSet cs;
    AnalyticForceBridge bridge(m, s, Real(0));
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, Real(0.01), 1, 0xE5);

    const int nDof = m.nu; // 6
    const double meanKE = 0.5 * nDof * kT300;
    Histogram h(0.0, 8.0 * meanKE, 30);
    const long N = 400000;
    for (long i = 0; i < N; ++i) {
        drv.seedMomenta();
        RobotEngine::realizeVelocity(m, s);
        h.add(static_cast<double>(RobotEngine::calcKineticEnergy(m, s)));
    }
    // Gamma(n/2, kT) per-bin weight = ∫_bin E^{n/2−1} e^{−E/kT} dE, computed by
    // sub-sampling the bin (center evaluation is biased for this curved density and
    // inflates chi^2 at large N).
    const double binW = (8.0 * meanKE) / h.nbins();
    std::vector<double> w(static_cast<std::size_t>(h.nbins()), 0.0);
    for (int b = 0; b < h.nbins(); ++b) {
        const double lo = h.center(b) - 0.5 * binW;
        const int sub = 16;
        double acc = 0.0;
        for (int k = 0; k < sub; ++k) {
            const double E = lo + (k + 0.5) * binW / sub;
            if (E > 0) {
                acc += std::pow(E, nDof / 2.0 - 1.0) * std::exp(-E / kT300);
            }
        }
        w[static_cast<std::size_t>(b)] = acc / sub;
    }
    const double chi = chiSquareStatistic(h.counts(), expectedFromWeights(w, h.total()));
    EXPECT_LT(chi, chiSquareCritical(h.nbins() - 1, 1e-4))
        << "KE distribution is not Gamma(n/2,kT): chi2=" << chi << " (n_dof=" << nDof << ")";
}

// ---------------------------------------------------------------------------
//  T0.2 -- KE ~ Gamma(n_dof/2, kT) on the SAME dense-M multi-body fixture T0.1
//  uses (RobotBuilders.hpp::buildBentTorsionChain, seed 0x70 -- literally the
//  same construction as TestEquipartition.cpp::MultiBodyBentChain), n_dof = 14.
//
//  The discriminating part is the CONTROL: the identical histogram must REJECT
//  Gamma((n_dof-1)/2,kT) and Gamma((n_dof+1)/2,kT) -- a plausible-but-wrong DOF
//  count -- at the same alpha. Without this, a test that only asserts the
//  correct-shape PASS cannot distinguish "the M^{-1/2} draw and n_dof count are
//  right" from "the shape test has no power here" (foundations Sec.5.2, tier0
//  Verification plan). This is the real gate; the PASS half alone is not.
//
//  Gated behind ROBOSAMPLE_SLOW_TESTS: alpha=1e-4 goodness-of-fit at N large
//  enough to also have power against a +/-1 DOF shift needs more draws than the
//  always-on equipartition mean checks, and each draw here costs an O(n_dof)
//  articulated-body sweep on a 14-DOF dense-M fixture (heavier than freeBody's
//  6-DOF always-on KEisMaxwellBoltzmann above).
// ---------------------------------------------------------------------------
TEST(EnsembleValidation, KEGammaMultiBodyBentChainWithWrongDofControl) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    const int nTorsions = 8;
    Rng rng(0x70); // SAME seed as TestEquipartition.cpp::MultiBodyBentChain
    RobotModel m = buildBentTorsionChain(nTorsions, rng);
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1; // unit quaternion for the Free root
    ConstraintSet cs;
    AnalyticForceBridge bridge(m, s, Real(0));
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, Real(0.01), 1, 0xE580);

    const int nDof = m.nu; // 6 + 8 = 14
    ASSERT_EQ(nDof, 14) << "fixture drifted -- update the n_dof-based comments/expectations above";
    const double meanKE = 0.5 * nDof * kT300;
    const double hiBound = 8.0 * meanKE;
    Histogram h(0.0, hiBound, 30);
    const long N = 1'200'000;
    for (long i = 0; i < N; ++i) {
        drv.seedMomenta();
        RobotEngine::realizeVelocity(m, s);
        h.add(static_cast<double>(RobotEngine::calcKineticEnergy(m, s)));
    }

    const double crit = chiSquareCritical(h.nbins() - 1, 1e-4);

    const double chiCorrect =
        chiSquareStatistic(h.counts(), expectedFromWeights(gammaSubsampledWeights(h, hiBound, nDof / 2.0, kT300), h.total()));
    EXPECT_LT(chiCorrect, crit) << "KE distribution is not Gamma(n_dof/2,kT): chi2=" << chiCorrect << " crit=" << crit
                                << " (n_dof=" << nDof << ")";

    const double chiMinusOne = chiSquareStatistic(
        h.counts(), expectedFromWeights(gammaSubsampledWeights(h, hiBound, (nDof - 1) / 2.0, kT300), h.total()));
    EXPECT_GT(chiMinusOne, crit) << "CONTROL FAILED: histogram should REJECT Gamma((n_dof-1)/2,kT)=Gamma("
                                 << (nDof - 1) / 2.0 << ",kT) but chi2=" << chiMinusOne << " <= crit=" << crit
                                 << " -- the test cannot distinguish a correct DOF count from an off-by-one";

    const double chiPlusOne = chiSquareStatistic(
        h.counts(), expectedFromWeights(gammaSubsampledWeights(h, hiBound, (nDof + 1) / 2.0, kT300), h.total()));
    EXPECT_GT(chiPlusOne, crit) << "CONTROL FAILED: histogram should REJECT Gamma((n_dof+1)/2,kT)=Gamma("
                                << (nDof + 1) / 2.0 << ",kT) but chi2=" << chiPlusOne << " <= crit=" << crit
                                << " -- the test cannot distinguish a correct DOF count from an off-by-one";
}

// ---------------------------------------------------------------------------
//  PE log-ratio between two temperatures is linear with slope −(β₂−β₁)
//  (Shirts eq. 20). A bound Free body (harmonic AnalyticForceBridge) sampled at
//  T1 and T2; the slope of ln[P(E|β₂)/P(E|β₁)] over the overlap region matches.
// ---------------------------------------------------------------------------
TEST(EnsembleValidation, PEobeysEnsembleSlope) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    const double T1 = 280.0, T2 = 320.0;
    const double b1 = 1.0 / (kB * T1), b2 = 1.0 / (kB * T2);

    auto samplePE = [&](double T, std::uint64_t seed) {
        RobotModel m = freeBody();
        RobotState s;
        s.allocateFull(m);
        s.q()[0] = 1;
        RobotEngine::realizePosition(m, s);
        ConstraintSet cs;
        AnalyticForceBridge bridge(m, s, /*k*/ Real(150)); // harmonic well -> bound PE
        HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kB * T, Real(2e-3), 10, seed);
        Histogram h(0.0, 60.0, 40);
        const long N = 1'500'000;
        for (long i = 0; i < N; ++i) {
            drv.move();
            RobotEngine::realizePosition(m, s);
            RobotEngine::fillAtomPositionsFromBodies(m, s);
            bridge.evaluate(s);
            h.add(static_cast<double>(bridge.calcPotentialEnergy(s)));
        }
        return h;
    };

    const Histogram h1 = samplePE(T1, 0xA1);
    const Histogram h2 = samplePE(T2, 0xA2);

    // Weighted least-squares slope of ln(n2/n1) vs E over bins both well populated.
    double Sw = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    for (int b = 0; b < h1.nbins(); ++b) {
        const long c1 = h1.counts()[static_cast<std::size_t>(b)];
        const long c2 = h2.counts()[static_cast<std::size_t>(b)];
        if (c1 < 500 || c2 < 500) {
            continue; // restrict to well-populated overlap bins (tail bias is large)
        }
        const double E = h1.center(b);
        const double y = std::log(static_cast<double>(c2) / static_cast<double>(c1));
        const double wgt = 1.0 / (1.0 / c1 + 1.0 / c2); // inverse variance of y
        Sw += wgt;
        Sx += wgt * E;
        Sy += wgt * y;
        Sxx += wgt * E * E;
        Sxy += wgt * E * y;
    }
    ASSERT_GT(Sw, 0.0);
    const double slope = (Sw * Sxy - Sx * Sy) / (Sw * Sxx - Sx * Sx);
    const double expected = -(b2 - b1);
    // The histogram-ratio slope is a known-noisy estimator (Shirts §2.4 notes large
    // statistical error); assert the correct sign and magnitude within 25%.
    EXPECT_NEAR(slope, expected, 0.25 * std::abs(expected))
        << "PE ensemble-validation slope " << slope << " != -(β₂−β₁)=" << expected;
}
