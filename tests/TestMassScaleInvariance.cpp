// ============================================================================
//  TestMassScaleInvariance.cpp -- fictitious mass scaling enters ONLY the proposal
//  (momentum draw, KE, Fixman ln det M) and must cancel in dH, so it raises the
//  stable timestep ~sqrt(scale) with ZERO configurational bias.
//
//  TWO CLAIMS (World::setMassScaleByJoint / RobotModel::bodyMassScale):
//   1. STABLE-dt SCALING. Uniformly scaling the spatial inertia by s slows every
//      mode frequency by sqrt(s), so the largest reversible Verlet step grows by
//      ~sqrt(s). Checked with RobotEngine::checkReversibility (fast).
//   2. CONFIGURATIONAL INVARIANCE. Scaling M by s multiplies det M by s^{nu}, a
//      q-INDEPENDENT constant, so the metric weight sqrt(det M) keeps the same
//      SHAPE and the configurational marginal is unchanged. Verified by sampling
//      the phi2 marginal of the two-torsion chain at s=1 and s=4 and showing both
//      match the same (unscaled) sqrt(det M) reference (slow tier).
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
#include "RobotIntegrator.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "StatTest.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::attachAtoms;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::HmcDriver;
using rtest::randomizeState;
using rtest::Rng;
using rtest::stat::chiSquareCritical;
using rtest::stat::chiSquareStatistic;
using rtest::stat::expectedFromWeights;
using rtest::stat::Histogram;

namespace {

constexpr double kT300 = 0.0083144626 * 300.0;

bool slowEnabled() {
    return std::getenv("ROBOSAMPLE_SLOW_TESTS") != nullptr;
}

RobotModel twoTorsionChain() {
    const Rotation bend(Real(M_PI_2), XAxis);
    BodySpec b1;
    b1.parent = 0;
    b1.joint = JointType::Torsion;
    b1.mass = Real(1.5);
    b1.com_B = Vec3(Real(0.12), 0, 0);
    b1.inertia_B = UnitInertia(Real(0.30), Real(0.45), Real(0.55));
    BodySpec b2;
    b2.parent = 1;
    b2.joint = JointType::Torsion;
    b2.X_PF = Transform(bend, Vec3(Real(0.15), 0, 0));
    b2.mass = Real(1.0);
    b2.com_B = Vec3(Real(0.14), 0, Real(0.05));
    b2.inertia_B = UnitInertia(Real(0.25), Real(0.40), Real(0.50));
    return buildForest({b1, b2});
}

void setUniformMassScale(RobotModel& m, double scale) {
    m.bodyMassScale.assign(static_cast<std::size_t>(m.numBodies), static_cast<Real>(scale));
}

// A shell of Free-rooted "solvent" bodies on Ground -- the explicit-solvent case
// whose high-frequency librational modes set the stable timestep and which
// mass-scaling is meant to tame (World.hpp: setMassScaleByJoint(Free, 16)).
RobotModel freeShell(int n) {
    std::vector<BodySpec> specs;
    for (int i = 0; i < n; ++i) {
        BodySpec b;
        b.parent = 0;
        b.joint = JointType::Free;
        b.mass = Real(1.0);
        b.inertia_B = UnitInertia(Real(0.4), Real(0.5), Real(0.6));
        specs.push_back(b);
    }
    RobotModel m = buildForest(specs);
    for (int b = 1; b < m.numBodies; ++b) {
        attachAtoms(m, b, {Vec3(0, 0, 0), Vec3(0.10, 0, 0), Vec3(0, 0.10, 0)}, {16.0, 1.0, 1.0});
    }
    return m;
}

double wrapPi(double a) {
    return std::atan2(std::sin(a), std::cos(a));
}

double detMAt(RobotModel& m, RobotState& s, double phi1, double phi2) {
    s.q()[0] = static_cast<Real>(phi1);
    s.q()[1] = static_cast<Real>(phi2);
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    return std::exp(static_cast<double>(RobotEngine::calcLogDetM(m, s)));
}

// Largest reversible step found by geometric scan: the biggest h whose round-trip
// residual stays under tol. A harmonic AnalyticForceBridge (k>0) sets the mode
// frequencies the step must resolve.
Real findStableH(RobotModel& m, Real k, Real tol) {
    RobotState s;
    s.allocateFull(m);
    Rng rng(0x7E57);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    AnalyticForceBridge bridge(m, s, k);
    ConstraintSet cs;
    // seed a thermal velocity so the round trip is non-trivial.
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    std::vector<Real> g(static_cast<std::size_t>(m.nu)), u(static_cast<std::size_t>(m.nu));
    for (int i = 0; i < m.nu; ++i) {
        g[i] = rng.gaussian();
    }
    RobotEngine::multiplyBySqrtMInv(m, s, g.data(), u.data());
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = std::sqrt(static_cast<Real>(kT300)) * u[i];
    }
    // seed the derivative chain (qdot0/udot0/qdotdot0) the first verletStep reads.
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    bridge.evaluate(s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    RobotEngine::calcUDot(m, s);
    RobotEngine::calcQDot(m, s, s.qdot());
    RobotEngine::calcQDotDot(m, s);

    // The symmetric step is reversible up to a small floor (set by the implicit-
    // trapezoid corrector tol), then the residual JUMPS to O(1) at the instability
    // onset h_c ~ 1/omega_max ~ sqrt(scale). `tol` separates the floor from blow-up.
    // The round-trip residual grows monotonically as res(h, s) = F(h / sqrt(s))
    // (the only dynamical scale is omega ~ 1/sqrt(s)). The largest h whose residual
    // stays under a fixed level c is therefore h_c(s) = sqrt(s) * F^{-1}(c), so the
    // ratio h_c(4)/h_c(1) = 2 for ANY level c the curve crosses in range -- the
    // sqrt(scale) law. `tol` is chosen so both s=1 and s=4 cross within the scan.
    Real best = 0;
    for (Real h = Real(5e-5); h < Real(0.5); h *= Real(1.12)) {
        const Real res = RobotEngine::checkReversibility(m, s, bridge, cs, /*nSteps=*/100, h);
        if (res < tol) {
            best = h;
        }
    }
    return best;
}

Histogram samplePhi2(double scale, bool useFixman, std::uint64_t seed, long nMoves, int stride) {
    RobotModel m = twoTorsionChain();
    setUniformMassScale(m, scale);
    RobotState s;
    s.allocateFull(m);
    for (int i = 0; i < m.nq; ++i) {
        s.q()[i] = 0;
    }
    ConstraintSet cs;
    AnalyticForceBridge bridge(m, s, Real(0));
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, Real(0.02), 12, seed);
    drv.useFixman = useFixman;
    Histogram h(-M_PI, M_PI, 24);
    for (long i = 0; i < nMoves; ++i) {
        drv.move();
        if (i % stride == 0) {
            h.add(wrapPi(static_cast<double>(s.q()[1])));
        }
    }
    return h;
}

} // namespace

// ---------------------------------------------------------------------------
//  STABLE-dt SCALING: the largest reversible step grows ~sqrt(scale).
// ---------------------------------------------------------------------------
TEST(MassScaleInvariance, StableStepGrowsAsSqrtScale) {
    RobotModel m1 = twoTorsionChain();
    RobotModel m4 = twoTorsionChain();
    setUniformMassScale(m4, 4.0);

    // A residual level both s=1 and s=4 cross within the scan (their crossings sit
    // near h~0.1 and h~0.2, both < 0.5), so the sqrt(scale)=2 ratio is observable.
    const Real tol = Real(1e-5);
    const Real h1 = findStableH(m1, /*k*/ Real(60), tol);
    const Real h4 = findStableH(m4, Real(60), tol);
    ASSERT_GT(h1, Real(0)) << "no stable step found at scale 1";
    ASSERT_GT(h4, Real(0)) << "no stable step found at scale 4";

    const double ratio = static_cast<double>(h4 / h1);
    // sqrt(4) = 2; the geometric scan (1.15x ladder) resolves the ratio coarsely.
    EXPECT_GT(ratio, 1.5) << "stable step did not grow with mass scale (ratio " << ratio << ")";
    EXPECT_LT(ratio, 2.7) << "stable step grew too much (ratio " << ratio << ")";
}

// ---------------------------------------------------------------------------
//  The same sqrt(scale) law on a Free-rooted solvent shell -- the explicit-solvent
//  case the feature targets (stiff librational modes set dt; mass-scaling tames them).
// ---------------------------------------------------------------------------
TEST(MassScaleInvariance, StableStepGrowsAsSqrtScale_FreeShell) {
    RobotModel m1 = freeShell(4);
    RobotModel m4 = freeShell(4);
    setUniformMassScale(m4, 4.0);

    const Real tol = Real(1e-5);
    const Real h1 = findStableH(m1, /*k*/ Real(60), tol);
    const Real h4 = findStableH(m4, Real(60), tol);
    ASSERT_GT(h1, Real(0)) << "no stable step found for the Free shell at scale 1";
    ASSERT_GT(h4, Real(0)) << "no stable step found for the Free shell at scale 4";

    const double ratio = static_cast<double>(h4 / h1);
    EXPECT_GT(ratio, 1.5) << "Free-shell stable step did not grow ~sqrt(scale) (ratio " << ratio << ")";
    EXPECT_LT(ratio, 2.7) << "Free-shell stable step grew too much (ratio " << ratio << ")";
}

// ---------------------------------------------------------------------------
//  CONFIGURATIONAL INVARIANCE: the phi2 marginal at scale 1 and scale 4 both match
//  the same unscaled sqrt(det M) reference -- mass scaling leaves the configuration
//  distribution unchanged (Fixman OFF here, so the marginal is sqrt(det M)).
// ---------------------------------------------------------------------------
TEST(MassScaleInvariance, ConfigurationalMarginalUnchanged) {
    if (!slowEnabled()) {
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    const Histogram h1 = samplePhi2(/*scale*/ 1.0, /*useFixman*/ false, 0xA1, 1'000'000, 4);
    const Histogram h4 = samplePhi2(/*scale*/ 4.0, /*useFixman*/ false, 0xA4, 1'000'000, 4);

    // Unscaled sqrt(det M) reference (scale only rescales det M by a constant).
    RobotModel m = twoTorsionChain();
    RobotState s;
    s.allocateFull(m);
    std::vector<double> w(static_cast<std::size_t>(h1.nbins()), 0.0);
    for (int b = 0; b < h1.nbins(); ++b) {
        w[static_cast<std::size_t>(b)] = std::sqrt(detMAt(m, s, 0.0, h1.center(b)));
    }

    const double crit = chiSquareCritical(h1.nbins() - 1, 1e-4);
    const double chi1 = chiSquareStatistic(h1.counts(), expectedFromWeights(w, h1.total()));
    const double chi4 = chiSquareStatistic(h4.counts(), expectedFromWeights(w, h4.total()));
    EXPECT_LT(chi1, crit) << "scale=1 marginal off sqrt(det M): chi2=" << chi1;
    EXPECT_LT(chi4, crit) << "scale=4 marginal off sqrt(det M): chi2=" << chi4 << " (mass scale biased config!)";
}
