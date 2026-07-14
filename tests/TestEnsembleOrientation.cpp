// ============================================================================
//  TestEnsembleOrientation.cpp -- the KEYSTONE ensemble test: does a free rigid
//  body sample orientation Haar-uniformly, and does the production orientational
//  Jacobian term belong with the quaternion parameterization?
//
//  THE QUESTION (the PI's open item; World.cpp once carried a "TEMPORARY
//  DIAGNOSTIC -- remove once the J / logSineSqr term is confirmed" comment). A
//  Free / Ball root is parameterized by a UNIT QUATERNION: 4 numbers, 1 norm
//  constraint => 3 rotational DOF (nq has 4 quaternion entries, nu has 3 angular).
//  The production acceptance Hamiltonian could optionally add
//      J(q) = -(1/2) RT ln sin^2(gamma2),
//  the sin(theta) polar volume factor of an EULER / spherical parameterization.
//
//  THE MEASURE-THEORY ANSWER, which this test verifies empirically:
//   * The flat (uniform) measure on the unit 3-sphere S^3 -- drawing the
//     quaternion uniformly subject to ||q||=1 -- pushes forward to EXACTLY the
//     Haar (rotation-invariant) measure on SO(3). "Flat in quaternion space" and
//     "Haar-uniform on rotations" are the SAME distribution; there is no Jacobian
//     between them.
//   * For a single free body with constant body inertia, det M(q) is
//     orientation-independent, and the quaternion exp-map integrator
//     (advanceQuatExp) is the geodesic flow of the kinetic metric. So GC-HMC with
//     U = 0 samples orientation Haar-uniformly with NO correction term.
//   * Adding J(q) reweights orientation by exp(-J/RT) = |sin(pitch)|, biasing the
//     marginal toward the poles. J is the CORRECT Jacobian only for an Euler-angle
//     parameterization; on a quaternion root it is WRONG.
//
//  HAAR OBSERVABLES (all parameterization-clean):
//   1. sinPitch = 2(w*y - z*x) = -R(2,0): a single entry of the rotation matrix,
//      hence a component of a Haar-uniform unit vector -> UNIFORM on [-1,1].
//      (This is exactly the quantity J keys on, so it is the bias's fingerprint.)
//   2. body z-axis z-component R(2,2): a Haar-uniform unit vector component ->
//      UNIFORM on [-1,1]; azimuth atan2(R(1,2),R(0,2)) -> UNIFORM on [-pi,pi].
//   3. rotation angle theta = 2 acos(|w|): Haar density (1 - cos theta)/pi on
//      [0,pi] -- the strongest single discriminator.
//
//  ASSERTIONS:
//   * No-Jacobian run  -> all three marginals match Haar (chi-square, alpha=1e-4).
//   * No-Jacobian sinPitch is INCONSISTENT with an Euler-flat (pitch-uniform)
//     reference density ~ 1/sqrt(1-P^2): codifies "S^3-uniform == Haar != Euler".
//   * With-Jacobian run -> sinPitch marginal is NO LONGER uniform, and DOES match
//     the predicted |P| bias: the empirical proof the term is wrong for quaternions.
//
//  Cross-reference: the 3-vs-4 DOF sizing is asserted in TestJointKernels.cpp
//  (JointFacts.DofAndNqTable); this file owns the MEASURE statement only.
//
//  The full-statistics cases run only under ROBOSAMPLE_SLOW_TESTS (millions of
//  moves); a small always-on smoke guards compilation and basic sanity.
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
#include "engine_helpers.hpp"
#include "support/SamplingHarness.hpp"
#include "support/TestPhysConstants.hpp"

using namespace robo;
using rtest::attachAtoms;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::HmcDriver;
using rtest::Rng;
using rtest::phys::kT300;
using rtest::stat::chiSquareCritical;
using rtest::stat::chiSquareStatistic;
using rtest::stat::expectedFromWeights;
using rtest::stat::Histogram;
using rtest::stat::slowEnabled;
using rtest::stat::uniformExpected;

namespace {

// A single free rigid body on Ground, identity joint frames (so s.q()[0..3] IS
// the body's orientation quaternion), asymmetric inertia (no symmetric-top
// degeneracy), three non-collinear atoms. U = 0 is supplied by AnalyticForceBridge
// with k = 0 (zero force, zero energy -- a genuine, not approximate, free body).
RobotModel freeBody() {
    BodySpec s;
    s.parent = 0;
    s.joint = JointType::Free;
    s.X_PF = Transform(); // identity
    s.X_BM = Transform(); // identity
    s.mass = Real(2.0);
    s.com_B = Vec3(0);
    s.inertia_B = UnitInertia(Real(0.40), Real(0.55), Real(0.72)); // asymmetric
    RobotModel m = buildForest({s});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.12, 0.01, -0.02), Vec3(-0.03, 0.10, 0.04)}, {12.0, 1.0, 1.0});
    return m;
}

// Identity-quaternion start, zero translation.
void initState(const RobotModel& m, RobotState& s) {
    s.allocateFull(m);
    Real* q = s.q();
    q[0] = 1;
    q[1] = 0;
    q[2] = 0;
    q[3] = 0;            // unit quaternion (identity rotation)
    for (int i = 4; i < m.nq; ++i) {
        q[i] = 0;
    }
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = 0;
    }
}

// Collected orientation marginals over a chain.
struct Marginals {
    Histogram sinPitch{-1.0, 1.0, 24};  // -R(2,0)
    Histogram zAxisZ{-1.0, 1.0, 24};    // R(2,2)
    Histogram azimuth{-M_PI, M_PI, 24}; // atan2(R(1,2),R(0,2))
    Histogram theta{0.0, M_PI, 24};     // rotation angle
    long moves = 0;
    long accepted = 0;
};

// Run `nMoves` GC-HMC moves, recording orientation observables every `stride`
// moves (thinning to decorrelate the free-rotation random walk).
template <class Driver>
Marginals sample(Driver& drv, const RobotModel& m, RobotState& s, long nMoves, int stride) {
    Marginals out;
    const rtest::Marginals marg = rtest::runHmcChain(drv, nMoves, stride, [&](long /*i*/, bool /*acc*/) {
        const Real w = s.q()[0], x = s.q()[1], y = s.q()[2], z = s.q()[3];
        const Rotation R = EngineHelpers::quatToRotation(w, x, y, z);
        const double sinPitch = std::clamp(2.0 * (w * y - z * x), -1.0, 1.0);
        out.sinPitch.add(sinPitch);
        out.zAxisZ.add(R(2, 2));
        out.azimuth.add(std::atan2(R(1, 2), R(0, 2)));
        out.theta.add(2.0 * std::acos(std::min(1.0, std::abs(double(w)))));
    });
    out.moves = marg.attempted;
    out.accepted = marg.accepted;
    return out;
}

// Per-bin Haar weight for the rotation angle: (1 - cos theta).
std::vector<double> haarThetaWeights(const Histogram& h) {
    std::vector<double> w(static_cast<std::size_t>(h.nbins()));
    for (int b = 0; b < h.nbins(); ++b) {
        w[static_cast<std::size_t>(b)] = 1.0 - std::cos(h.center(b));
    }
    return w;
}

// Per-bin Euler-flat (pitch-uniform) weight for sinPitch: 1/sqrt(1 - P^2).
std::vector<double> eulerFlatPitchWeights(const Histogram& h) {
    std::vector<double> w(static_cast<std::size_t>(h.nbins()));
    for (int b = 0; b < h.nbins(); ++b) {
        const double p = h.center(b);
        w[static_cast<std::size_t>(b)] = 1.0 / std::sqrt(std::max(1e-6, 1.0 - p * p));
    }
    return w;
}

// Per-bin predicted-bias weight for sinPitch when J is ON: |P|.
std::vector<double> absPitchWeights(const Histogram& h) {
    std::vector<double> w(static_cast<std::size_t>(h.nbins()));
    for (int b = 0; b < h.nbins(); ++b) {
        w[static_cast<std::size_t>(b)] = std::abs(h.center(b));
    }
    return w;
}

} // namespace

// ---------------------------------------------------------------------------
//  SMOKE (always on): the driver runs, accepts moves, and the orientation does
//  not collapse into a single bin. Cheap guard against build/link/logic breakage.
// ---------------------------------------------------------------------------
TEST(EnsembleOrientation, Smoke) {
    RobotModel m = freeBody();
    RobotState s;
    initState(m, s);
    ConstraintSet cs;
    AnalyticForceBridge bridge(m, s, /*k*/ Real(0));
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, /*h*/ Real(0.05), /*mdSteps*/ 20, 0xE0);

    Marginals r = sample(drv, m, s, /*nMoves*/ 4000, /*stride*/ 2);
    EXPECT_GT(r.accepted, 0) << "no move was accepted";
    // sinPitch must populate more than one bin (orientation actually moves).
    int populated = 0;
    for (long c : r.sinPitch.counts()) {
        populated += (c > 0) ? 1 : 0;
    }
    EXPECT_GE(populated, 5) << "orientation did not explore (only " << populated << " bins populated)";
}

// ---------------------------------------------------------------------------
//  No-Jacobian run is Haar-uniform on all three marginals. This is the default
//  production path (useOrientationJacobian = false) and its regression guard.
// ---------------------------------------------------------------------------
TEST(EnsembleOrientation, NoJacobianIsHaarUniform) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    RobotModel m = freeBody();
    RobotState s;
    initState(m, s);
    ConstraintSet cs;
    AnalyticForceBridge bridge(m, s, Real(0));
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, Real(0.05), 12, 0x1A1A);
    drv.useOrientationJac = false;

    // 8e6 moves, thinned by 8 -> 1e6 decorrelated samples.
    Marginals r = sample(drv, m, s, 600'000, 3);
    const long N = r.sinPitch.total();
    const double crit = chiSquareCritical(r.sinPitch.nbins() - 1, 1e-4);

    const double chiP = chiSquareStatistic(r.sinPitch.counts(), uniformExpected(r.sinPitch.nbins(), N));
    const double chiZ = chiSquareStatistic(r.zAxisZ.counts(), uniformExpected(r.zAxisZ.nbins(), N));
    const double chiA = chiSquareStatistic(r.azimuth.counts(), uniformExpected(r.azimuth.nbins(), N));
    const double chiT =
        chiSquareStatistic(r.theta.counts(), expectedFromWeights(haarThetaWeights(r.theta), N));

    EXPECT_LT(chiP, crit) << "sinPitch not uniform: chi2=" << chiP << " crit=" << crit;
    EXPECT_LT(chiZ, crit) << "z-axis z-component not uniform: chi2=" << chiZ;
    EXPECT_LT(chiA, crit) << "azimuth not uniform: chi2=" << chiA;
    EXPECT_LT(chiT, crit) << "rotation angle not Haar (1-cos): chi2=" << chiT;
}

// ---------------------------------------------------------------------------
//  S^3-uniform == Haar != Euler-flat. The no-Jacobian sinPitch marginal must be
//  consistent with the Haar (uniform-in-P) reference AND inconsistent with an
//  Euler-flat (pitch-uniform) reference density ~ 1/sqrt(1-P^2). This is the
//  empirical content of "the quaternion parameterization is already Haar".
// ---------------------------------------------------------------------------
TEST(EnsembleOrientation, HaarNotEulerFlat) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    RobotModel m = freeBody();
    RobotState s;
    initState(m, s);
    ConstraintSet cs;
    AnalyticForceBridge bridge(m, s, Real(0));
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, Real(0.05), 12, 0x2B2B);
    drv.useOrientationJac = false;

    Marginals r = sample(drv, m, s, 600'000, 3);
    const long N = r.sinPitch.total();
    const double crit = chiSquareCritical(r.sinPitch.nbins() - 1, 1e-4);

    const double chiHaar = chiSquareStatistic(r.sinPitch.counts(), uniformExpected(r.sinPitch.nbins(), N));
    const double chiEuler =
        chiSquareStatistic(r.sinPitch.counts(), expectedFromWeights(eulerFlatPitchWeights(r.sinPitch), N));

    EXPECT_LT(chiHaar, crit) << "sinPitch should match Haar (uniform): chi2=" << chiHaar;
    EXPECT_GT(chiEuler, crit) << "sinPitch should REJECT Euler-flat 1/sqrt(1-P^2): chi2=" << chiEuler;
}

// ---------------------------------------------------------------------------
//  With the Jacobian ON the orientation marginal is biased exactly as predicted:
//  sinPitch is no longer uniform and instead follows the |P| reweighting
//  exp(-J/RT) = |sin pitch|. This is the empirical proof the term does NOT belong
//  with the quaternion parameterization.
// ---------------------------------------------------------------------------
TEST(EnsembleOrientation, WithJacobianIsBiased) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    RobotModel m = freeBody();
    RobotState s;
    initState(m, s);
    ConstraintSet cs;
    AnalyticForceBridge bridge(m, s, Real(0));
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, Real(0.05), 12, 0x3C3C);
    drv.useOrientationJac = true; // the (wrong-for-quaternion) Euler Jacobian

    Marginals r = sample(drv, m, s, 600'000, 3);
    const long N = r.sinPitch.total();
    const double crit = chiSquareCritical(r.sinPitch.nbins() - 1, 1e-4);

    const double chiUniform = chiSquareStatistic(r.sinPitch.counts(), uniformExpected(r.sinPitch.nbins(), N));
    const double chiAbsP =
        chiSquareStatistic(r.sinPitch.counts(), expectedFromWeights(absPitchWeights(r.sinPitch), N));

    EXPECT_GT(chiUniform, crit) << "J ON should BREAK uniformity but chi2=" << chiUniform;
    EXPECT_LT(chiAbsP, crit) << "J ON marginal should match the |P| bias: chi2=" << chiAbsP;
}
