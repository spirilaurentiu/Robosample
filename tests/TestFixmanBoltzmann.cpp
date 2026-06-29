// ============================================================================
//  TestFixmanBoltzmann.cpp -- the Fixman keystone (Pear-Weiner benchmark).
//
//  THEORY. Generalized-coordinate HMC draws the momentum p ~ N(0, M(q)), so
//  integrating p out of the canonical measure exp(-beta H) leaves a configurational
//  marginal weighted by the metric volume:
//        rho(phi) ~ sqrt(det M(phi)) * exp(-beta U(phi)).
//  With U = 0 (free internal rotation) the bare torsional sampler is therefore
//  NOT flat -- it is rho_OFF(phi) ~ sqrt(det M(phi)) = sqrt(g(phi)). The Fixman
//  compensating potential
//        U_F(phi) = (1/2) RT ln det M(phi)            (World::calcFixman)
//  contributes exp(-beta U_F) = det M^{-1/2}, which cancels the metric factor and
//  flattens the marginal: rho_ON(phi) ~ const. This is exactly the Pear-Weiner /
//  Fixman statement (Fixman 1974; Patriciu, Chirikjian, Pappu 2004) and the claim
//  the calcFixman comment relies on.
//
//  SELF-CONTAINED REFERENCE. Rather than depend on a literature constant (the
//  often-quoted 1.438119 mean), this test computes its OWN g(phi) = det M(phi) on a
//  grid using the same engine operator the sampler uses (calcLogDetM after
//  realizeArticulatedBodyInertias). The reference curve and the sampled histogram
//  thus come from one source of truth; a divergence is a real Fixman bug, not a
//  mismatched external number.
//
//  FIXTURE. A two-torsion chain with a 90-degree bend between the two torsion axes
//  (the canonical fixed-bond-angle chain). The mass-metric coupling M_12 between
//  the non-parallel axes makes det M depend on the relative dihedral phi2, so
//  g(phi2) genuinely varies (the test asserts the variation is non-trivial, else it
//  would be vacuous). phi1 only rigidly reorients the whole chain in space and
//  leaves the internal metric unchanged, so phi1 is uniform and phi2 carries the
//  Fixman signal -- the role of "the torsion" in the Pear-Weiner setup.
//
//  Full statistics run under ROBOSAMPLE_SLOW_TESTS; a small smoke is always on.
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
using rtest::BodySpec;
using rtest::buildForest;
using rtest::HmcDriver;
using rtest::stat::chiSquareCritical;
using rtest::stat::chiSquareStatistic;
using rtest::stat::expectedFromWeights;
using rtest::stat::Histogram;
using rtest::stat::uniformExpected;

namespace {

constexpr double kT300 = 0.0083144626 * 300.0;

bool slowEnabled() {
    return std::getenv("ROBOSAMPLE_SLOW_TESTS") != nullptr;
}

// Two-torsion chain, 90-degree bend between the torsion axes. Off-axis body mass
// (com offset) makes the inter-axis metric coupling -- and hence det M(phi2) --
// non-trivial.
RobotModel twoTorsionChain() {
    const Rotation bend(Real(M_PI_2), XAxis); // tilt body2's torsion axis 90 deg
    BodySpec b1;
    b1.parent = 0;
    b1.joint = JointType::Torsion;
    b1.X_PF = Transform();
    b1.X_BM = Transform();
    b1.mass = Real(1.5);
    b1.com_B = Vec3(Real(0.12), Real(0.0), Real(0.0));
    b1.inertia_B = UnitInertia(Real(0.30), Real(0.45), Real(0.55));

    BodySpec b2;
    b2.parent = 1;
    b2.joint = JointType::Torsion;
    b2.X_PF = Transform(bend, Vec3(Real(0.15), Real(0.0), Real(0.0)));
    b2.X_BM = Transform();
    b2.mass = Real(1.0);
    b2.com_B = Vec3(Real(0.14), Real(0.0), Real(0.05));
    b2.inertia_B = UnitInertia(Real(0.25), Real(0.40), Real(0.50));

    return buildForest({b1, b2});
}

void initState(const RobotModel& m, RobotState& s) {
    s.allocateFull(m);
    for (int i = 0; i < m.nq; ++i) {
        s.q()[i] = 0;
    }
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = 0;
    }
}

double wrapPi(double a) {
    return std::atan2(std::sin(a), std::cos(a));
}

// det M as a function of the two torsions, via the engine's O(n) log-det operator.
double detMAt(RobotModel& m, RobotState& s, double phi1, double phi2) {
    s.q()[0] = static_cast<Real>(phi1);
    s.q()[1] = static_cast<Real>(phi2);
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    return std::exp(static_cast<double>(RobotEngine::calcLogDetM(m, s)));
}

// Reference per-bin weight for the phi2 marginal: sqrt(det M(phi2)), marginalized
// over a small phi1 grid (so the reference is correct even if det M depended on
// phi1, which it does not for this fixture).
std::vector<double> sqrtDetMWeights(RobotModel& m, RobotState& s, const Histogram& h) {
    const int nPhi1 = 8;
    std::vector<double> w(static_cast<std::size_t>(h.nbins()), 0.0);
    for (int b = 0; b < h.nbins(); ++b) {
        const double phi2 = h.center(b);
        double acc = 0.0;
        for (int k = 0; k < nPhi1; ++k) {
            const double phi1 = -M_PI + (k + 0.5) * (2.0 * M_PI / nPhi1);
            acc += std::sqrt(detMAt(m, s, phi1, phi2));
        }
        w[static_cast<std::size_t>(b)] = acc / nPhi1;
    }
    return w;
}

struct ChainRun {
    Histogram phi2{-M_PI, M_PI, 24};
    long accepted = 0;
    long moves = 0;
};

ChainRun sampleChain(bool useFixman, std::uint64_t seed, long nMoves, int stride) {
    RobotModel m = twoTorsionChain();
    RobotState s;
    initState(m, s);
    ConstraintSet cs;
    AnalyticForceBridge bridge(m, s, /*k*/ Real(0)); // U = 0
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, /*h*/ Real(0.02), /*mdSteps*/ 12, seed);
    drv.useFixman = useFixman;

    ChainRun r;
    for (long i = 0; i < nMoves; ++i) {
        r.accepted += drv.move() ? 1 : 0;
        ++r.moves;
        if (i % stride == 0) {
            r.phi2.add(wrapPi(static_cast<double>(s.q()[1])));
        }
    }
    return r;
}

} // namespace

// ---------------------------------------------------------------------------
//  SMOKE: the chain runs, the metric is genuinely phi2-dependent (else the whole
//  test is vacuous), and phi2 explores its range.
// ---------------------------------------------------------------------------
TEST(FixmanBoltzmann, Smoke) {
    RobotModel m = twoTorsionChain();
    RobotState s;
    initState(m, s);

    // det M(phi2) must vary appreciably for the Fixman test to mean anything.
    double lo = 1e300, hi = 0.0;
    for (int k = 0; k < 32; ++k) {
        const double phi2 = -M_PI + (k + 0.5) * (2.0 * M_PI / 32);
        const double g = detMAt(m, s, 0.0, phi2);
        lo = std::min(lo, g);
        hi = std::max(hi, g);
    }
    EXPECT_GT(hi / lo, 1.05) << "det M(phi2) is nearly flat (ratio " << hi / lo
                             << "); fixture would make the Fixman test vacuous";

    ChainRun r = sampleChain(/*useFixman*/ false, 0xF0, /*moves*/ 4000, /*stride*/ 2);
    EXPECT_GT(r.accepted, 0);
    int populated = 0;
    for (long c : r.phi2.counts()) {
        populated += (c > 0) ? 1 : 0;
    }
    EXPECT_GE(populated, 8) << "phi2 did not explore";
}

// ---------------------------------------------------------------------------
//  Fixman OFF: the phi2 marginal follows sqrt(det M), and REJECTS flat.
// ---------------------------------------------------------------------------
TEST(FixmanBoltzmann, OffMatchesSqrtDetM) {
    if (!slowEnabled()) {
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    ChainRun r = sampleChain(/*useFixman*/ false, 0xB0117, /*moves*/ 1'000'000, /*stride*/ 4);
    const long N = r.phi2.total();

    RobotModel m = twoTorsionChain();
    RobotState s;
    initState(m, s);
    const std::vector<double> wRef = sqrtDetMWeights(m, s, r.phi2);

    const double crit = chiSquareCritical(r.phi2.nbins() - 1, 1e-4);
    const double chiRef = chiSquareStatistic(r.phi2.counts(), expectedFromWeights(wRef, N));
    const double chiFlat = chiSquareStatistic(r.phi2.counts(), uniformExpected(r.phi2.nbins(), N));

    EXPECT_LT(chiRef, crit) << "OFF marginal should match sqrt(det M): chi2=" << chiRef << " crit=" << crit;
    EXPECT_GT(chiFlat, crit) << "OFF marginal should NOT be flat: chi2=" << chiFlat;
}

// ---------------------------------------------------------------------------
//  Fixman ON: the phi2 marginal is flat, and REJECTS the sqrt(det M) shape.
// ---------------------------------------------------------------------------
TEST(FixmanBoltzmann, OnIsFlat) {
    if (!slowEnabled()) {
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    ChainRun r = sampleChain(/*useFixman*/ true, 0xF1A7, /*moves*/ 1'000'000, /*stride*/ 4);
    const long N = r.phi2.total();

    RobotModel m = twoTorsionChain();
    RobotState s;
    initState(m, s);
    const std::vector<double> wRef = sqrtDetMWeights(m, s, r.phi2);

    const double crit = chiSquareCritical(r.phi2.nbins() - 1, 1e-4);
    const double chiFlat = chiSquareStatistic(r.phi2.counts(), uniformExpected(r.phi2.nbins(), N));
    const double chiRef = chiSquareStatistic(r.phi2.counts(), expectedFromWeights(wRef, N));

    EXPECT_LT(chiFlat, crit) << "ON marginal should be flat: chi2=" << chiFlat << " crit=" << crit;
    EXPECT_GT(chiRef, crit) << "ON marginal should NOT follow sqrt(det M): chi2=" << chiRef;
}
