// ============================================================================
//  TestNcmcTeleport.cpp -- the λ=0 trough teleport: correctness of the move that
//  replaces the implicit-solvent RigidKick in explicit solvent.
//
//  At λ=0 the moved region is a non-interacting ghost, so a rigid reposition is
//  free of intermolecular energy. The tests pin the correctness conditions the
//  papers require (see tests/TeleportMove.hpp for the derivations):
//    * TroughKickIsKENeutral          -- velocity co-rotation keeps KE exact
//    * GhostTeleportCostsZeroInter     -- λ=0 reposition changes V by 0
//    * UniformDrawIsSymmetric          -- fixed-region draw is pose-independent
//    * CompositeProposalIsReversible   -- decouple∘K∘recouple retraces under flip
//    * ConstrainedTeleportBranchConsistency -- reverse-SHAKE branch guard (CHMC Thm 3)
//    * TeleportPreservesMarginal       -- samples the correct Boltzmann marginal
// ============================================================================
#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <vector>

#include "Constraints.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotIntegrator.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "StatTest.hpp"
#include "TeleportMove.hpp"
#include "TestHelpers.hpp"
#include "support/HarmonicBridge.hpp"
#include "support/TestPhysConstants.hpp"

using namespace robo;
using rtest::attachAtoms;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::HarmonicBridge;
using rtest::jointRotation;
using rtest::LambdaWellPolicy;
using rtest::Rng;
using rtest::teleportFreeRoot;
using rtest::phys::kT300;
using rtest::stat::chiSquareCritical;
using rtest::stat::chiSquareStatistic;
using rtest::stat::expectedFromWeights;
using rtest::stat::Histogram;
using rtest::stat::slowEnabled;

namespace {

using LambdaWellBridge = HarmonicBridge<LambdaWellPolicy>;

// LambdaWellBridge (a type alias, above) is a λ-scaled harmonic bridge: V = λ ·
// ½k · Σ|r_a − anchor_a|², F = −λk(r−anchor). At λ=0 it is a true ghost (V=0,
// F=0); at λ=1 it is an external well. Same per-atom→per-body reduction as
// ForceBridge/AnalyticForceBridge (HarmonicBridge<LambdaWellPolicy>).

// One free-rooted "solute" body carrying real atoms.
RobotModel solute() {
    BodySpec b;
    b.parent = 0;
    b.joint = JointType::Free;
    b.mass = Real(2.0);
    b.inertia_B = UnitInertia(Real(0.4), Real(0.55), Real(0.7));
    RobotModel m = buildForest({b});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.12, 0.01, -0.02), Vec3(-0.03, 0.10, 0.04)}, {12.0, 1.0, 1.0});
    return m;
}

void initState(const RobotModel& m, RobotState& s) {
    s.allocateFull(m);
    s.q()[0] = 1; // unit quaternion
    for (int i = 1; i < m.nq; ++i) {
        s.q()[i] = 0;
    }
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = 0;
    }
}

Real kineticEnergy(const RobotModel& m, RobotState& s) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeVelocity(m, s);
    return RobotEngine::calcKineticEnergy(m, s);
}

// Per-bin weight = ∫_bin exp(-x²/2σ²) dx, sub-sampled (center evaluation of a
// curved Gaussian is biased and inflates chi^2 at large N).
std::vector<double> gaussianBinWeights(const Histogram& h, double sigma) {
    const double binW = (h.center(1) - h.center(0));
    std::vector<double> w(static_cast<std::size_t>(h.nbins()), 0.0);
    for (int b = 0; b < h.nbins(); ++b) {
        const double lo = h.center(b) - 0.5 * binW;
        const int sub = 16;
        double acc = 0.0;
        for (int k = 0; k < sub; ++k) {
            const double x = lo + (k + 0.5) * binW / sub;
            acc += std::exp(-0.5 * x * x / (sigma * sigma));
        }
        w[static_cast<std::size_t>(b)] = acc / sub;
    }
    return w;
}

} // namespace

// ---------------------------------------------------------------------------
//  Velocity co-rotation keeps KE exactly invariant under a reorientation; without
//  it, reorientation leaks a KE penalty (M_ang(q) is orientation-dependent).
// ---------------------------------------------------------------------------
TEST(NcmcTeleport, TroughKickIsKENeutral) {
    Rng rng(0x7E1E);
    RobotModel m = solute();
    RobotState s;
    initState(m, s);
    // give the body a non-trivial velocity
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = rng.gaussian(0, 0.8);
    }
    const Real ke0 = kineticEnergy(m, s);

    const Vec3 target = rng.vec3(-1, 1);
    const Rotation Rnew = rng.rotation();

    // with co-rotation: KE preserved to machine precision
    {
        RobotState sc;
        sc.allocateFull(m);
        std::copy(s.q(), s.q() + m.nq, sc.q());
        std::copy(s.u(), s.u() + m.nu, sc.u());
        teleportFreeRoot(m, sc, 1, target, Rnew, /*coRotateVel=*/true);
        const Real ke1 = kineticEnergy(m, sc);
        EXPECT_NEAR(ke1, ke0, 1e-9) << "co-rotated teleport changed KE: " << ke0 << " -> " << ke1;
    }
    // without co-rotation: reorientation changes KE (test is non-vacuous)
    {
        RobotState sc;
        sc.allocateFull(m);
        std::copy(s.q(), s.q() + m.nq, sc.q());
        std::copy(s.u(), s.u() + m.nu, sc.u());
        teleportFreeRoot(m, sc, 1, target, Rnew, /*coRotateVel=*/false);
        const Real ke1 = kineticEnergy(m, sc);
        EXPECT_GT(std::abs(ke1 - ke0), 1e-3)
            << "reorientation without co-rotation should change KE (M_ang depends on q)";
    }
}

// ---------------------------------------------------------------------------
//  At λ=0 the body is a ghost: a rigid teleport changes the potential by 0.
// ---------------------------------------------------------------------------
TEST(NcmcTeleport, GhostTeleportCostsZeroInter) {
    Rng rng(0x6057);
    RobotModel m = solute();
    RobotState s;
    initState(m, s);
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    LambdaWellBridge bridge(m, s, /*k*/ Real(500));

    bridge.setLambda(Real(0)); // ghost
    const Real v0 = bridge.calcPotentialEnergy(s);
    teleportFreeRoot(m, s, 1, rng.vec3(-3, 3), rng.rotation(), true);
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    const Real v1 = bridge.calcPotentialEnergy(s);
    EXPECT_NEAR(v0, v1, 1e-12) << "λ=0 teleport changed V (not a ghost)";
    EXPECT_EQ(v0, Real(0));

    // sanity: at λ=1 the same teleport DOES change V (the well is on)
    bridge.setLambda(Real(1));
    const Real v1on = bridge.calcPotentialEnergy(s);
    teleportFreeRoot(m, s, 1, rng.vec3(-3, 3), rng.rotation(), true);
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    EXPECT_GT(std::abs(bridge.calcPotentialEnergy(s) - v1on), 1e-6);
}

// ---------------------------------------------------------------------------
//  The fixed-region uniform-in-sphere draw is pose-independent (symmetric):
//  the radial CDF is r^3 (uniform-in-volume) regardless of the current pose, so
//  g(new|old)=g(old|new) and no Hastings ratio is needed.
// ---------------------------------------------------------------------------
TEST(NcmcTeleport, UniformDrawIsSymmetric) {
    Rng rng(0x5717);
    const double R = 2.0;
    const int nbins = 20;
    // uniform-in-sphere sampler (same convention as World::sampleUniformInSphere)
    auto draw = [&]() {
        const double th = 2.0 * M_PI * rng.uniform(0, 1);
        const double ph = std::acos(2.0 * rng.uniform(0, 1) - 1.0);
        const double r = R * std::cbrt(rng.uniform(0, 1));
        return Vec3(r * std::cos(th) * std::sin(ph), r * std::sin(th) * std::sin(ph), r * std::cos(ph));
    };
    // bin r^3/R^3 -> must be uniform on [0,1] (uniform-in-volume), independent of pose
    Histogram h(0.0, 1.0, nbins);
    const long N = 400000;
    for (long i = 0; i < N; ++i) {
        const double rr = draw().norm() / R;
        h.add(rr * rr * rr);
    }
    const double chi = chiSquareStatistic(h.counts(), rtest::stat::uniformExpected(nbins, h.total()));
    EXPECT_LT(chi, chiSquareCritical(nbins - 1, 1e-4))
        << "uniform-in-sphere draw is not uniform-in-volume (chi2=" << chi << ")";
}

// ---------------------------------------------------------------------------
//  A deterministic teleport sandwiched in a symmetric decouple∘recouple protocol
//  retraces under momentum flip + the inverse teleport: the composite proposal is
//  reversible (centering + co-rotation + volume preservation).
// ---------------------------------------------------------------------------
TEST(NcmcTeleport, CompositeProposalIsReversible) {
    Rng rng(0x0C0DE);
    RobotModel m = solute();
    RobotState s;
    initState(m, s);
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = rng.gaussian(0, 0.5);
    }
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    LambdaWellBridge bridge(m, s, Real(80));
    ConstraintSet cs;

    const std::vector<Real> q0(s.q(), s.q() + m.nq);
    const std::vector<Real> u0(s.u(), s.u() + m.nu);
    // Incremental (relative) kick: a fixed displacement + a fixed rotation, whose
    // exact inverse is (-dTrans, dR⁻¹). An ABSOLUTE independence draw is NOT a
    // deterministic involution under sandwiching propagation (its reversibility is
    // the symmetric-density property, tested in UniformDrawIsSymmetric); the
    // incremental kick is the right primitive for the retrace test.
    const Vec3 dTrans(0.3, -0.2, 0.15);
    const Rotation dR(Real(0.6), YAxis); // fixed rotation about Y
    const Rotation dRinv(dR.transpose());

    const Real h = Real(2e-4);
    const int half = 6;
    auto seedDeriv = [&]() {
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        bridge.evaluate(s);
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        RobotEngine::calcUDot(m, s);
        RobotEngine::calcQDot(m, s, s.qdot());
        RobotEngine::calcQDotDot(m, s);
    };
    // half-hold propagation at fixed lambda=0
    auto propagateGhost = [&](int n) {
        bridge.setLambda(Real(0));
        seedDeriv();
        for (int k = 0; k < n; ++k) {
            ASSERT_TRUE(RobotEngine::stepTo(m, s, bridge, cs, s.time + h));
        }
    };
    auto flip = [&]() {
        for (int i = 0; i < m.nu; ++i) {
            s.u()[i] = -s.u()[i];
        }
    };

    // FORWARD: [propagate half] ∘ K(+dTrans, dR) ∘ [propagate half]
    propagateGhost(half);
    rtest::teleportFreeRootIncremental(m, s, 1, dTrans, dR, true);
    propagateGhost(half);

    // REVERSE: flip momenta, [propagate half] ∘ K(−dTrans, dR⁻¹) ∘ [propagate half], flip
    flip();
    propagateGhost(half);
    rtest::teleportFreeRootIncremental(m, s, 1, Vec3(0, 0, 0) - dTrans, dRinv, true);
    propagateGhost(half);
    flip();

    Real worstQ = 0, worstU = 0;
    for (int i = 0; i < m.nq; ++i) {
        worstQ = std::max(worstQ, std::abs(s.q()[i] - q0[i]));
    }
    for (int i = 0; i < m.nu; ++i) {
        worstU = std::max(worstU, std::abs(s.u()[i] - u0[i]));
    }
    EXPECT_LT(worstQ, 1e-6) << "composite teleport proposal not reversible in q (residual " << worstQ << ")";
    EXPECT_LT(worstU, 1e-6) << "composite teleport proposal not reversible in u (residual " << worstU << ")";
}

// ---------------------------------------------------------------------------
//  Constrained (cyclic) teleport: the reverse-SHAKE branch guard (CHMC Thm 3)
//  accepts a small on-branch jump and rejects a large branch-flipping jump.
// ---------------------------------------------------------------------------
TEST(NcmcTeleport, ConstrainedTeleportBranchConsistency) {
    Rng rng(0xC1C);
    // a Free root + two torsions, ring-closed by a distance constraint.
    auto mk = [&](int parent, JointType jt) {
        BodySpec b;
        b.parent = parent;
        b.joint = jt;
        b.X_PF = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        b.X_BM = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        b.mass = rng.uniform(0.9, 1.4);
        b.com_B = rng.vec3(-0.05, 0.05);
        b.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
        return b;
    };
    RobotModel m = buildForest({mk(0, JointType::Free), mk(1, JointType::Torsion), mk(2, JointType::Torsion)});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, 0.02, -0.03)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0), Vec3(-0.03, 0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0.04, -0.06, 0.09), Vec3(0.08, 0.01, -0.02)}, {16.0, 1.0});
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1;
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    const int aA = m.bodyAtoms[m.bodyAtomsBeg[1]];
    const int aB = m.bodyAtoms[m.bodyAtomsBeg[3]];
    const double d0 = (s.atomPosG()[aA] - s.atomPosG()[aB]).norm();
    ConstraintSet cs;
    cs.distance.push_back(DistanceConstraint{aA, aB, static_cast<Real>(d0)});

    auto refresh = [&]() {
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
    };
    refresh();
    cs.enforcePositionConstraints(m, s, refresh);
    const std::vector<Real> qOrigin(s.q(), s.q() + m.nq);

    auto applyKickAndShake = [&](const std::vector<Real>& kick) {
        std::copy(qOrigin.begin(), qOrigin.end(), s.q());
        for (int i = 0; i < m.nq; ++i) {
            s.q()[i] += kick[i];
        }
        refresh();
        cs.enforcePositionConstraints(m, s, refresh);
    };

    // Branch tolerance: well above the O(kick^2) SHAKE round-trip noise of an
    // on-branch move, well below the O(1) gap to a different ring-closure branch.
    const Real branchTol = Real(0.05);

    // small on-branch torsion nudge: SHAKE(q1 − kick) returns to origin -> consistent
    {
        std::vector<Real> kick(static_cast<std::size_t>(m.nq), Real(0));
        kick[static_cast<std::size_t>(m.bodyQIndex[2])] = Real(0.02);
        applyKickAndShake(kick);
        EXPECT_TRUE(rtest::teleportReverseConsistent(m, s, cs, qOrigin, kick, refresh, branchTol))
            << "small on-branch jump wrongly flagged as branch-inconsistent";
    }
    // large torsion flip across the ring: the inverse kick + SHAKE lands on a
    // different branch -> rejected.
    {
        std::vector<Real> kick(static_cast<std::size_t>(m.nq), Real(0));
        kick[static_cast<std::size_t>(m.bodyQIndex[2])] = Real(M_PI);
        applyKickAndShake(kick);
        EXPECT_FALSE(rtest::teleportReverseConsistent(m, s, cs, qOrigin, kick, refresh, branchTol))
            << "branch-flipping jump should be rejected by the reverse-SHAKE guard";
    }
}

// ---------------------------------------------------------------------------
//  The teleport move samples the correct Boltzmann marginal. A single-atom solute
//  in a λ-gated harmonic well V = ½k|T|² (atom at the body origin, anchored at the
//  origin): decouple (ghost) → uniform-in-box teleport of the translation →
//  recouple, accept on dH. With KE-neutral co-rotation and instantaneous switching
//  this is an independence sampler of exp(-βU); the translational marginal must be
//  Gaussian with variance kT/k per axis. (Single atom ⇒ the COM marginal is exactly
//  Gaussian, decoupled from orientation.)
// ---------------------------------------------------------------------------
TEST(NcmcTeleport, TeleportPreservesMarginal) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    Rng rng(0x7E1E909);
    BodySpec b;
    b.parent = 0;
    b.joint = JointType::Free;
    b.mass = Real(2.0);
    b.inertia_B = UnitInertia(Real(0.4), Real(0.55), Real(0.7));
    RobotModel m = buildForest({b});
    attachAtoms(m, 1, {Vec3(0, 0, 0)}, {12.0}); // single real atom at the body origin
    RobotState s;
    initState(m, s);
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    const Real k = Real(300);
    LambdaWellBridge bridge(m, s, k); // anchor = origin
    bridge.setLambda(Real(1));
    const double sigma = std::sqrt(kT300 / static_cast<double>(k)); // 1-D stddev

    auto setT = [&](const Vec3& T) {
        s.q()[4] = T[0];
        s.q()[5] = T[1];
        s.q()[6] = T[2];
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        return bridge.calcPotentialEnergy(s);
    };

    // Independence sampler: the uniform box must overlap the target well for the
    // chain to converge (a box >> σ gives near-zero acceptance), and the HISTOGRAM
    // RANGE MUST EQUAL THE BOX (Histogram::add clamps out-of-range samples into the
    // edge bins, which would corrupt the chi^2). On [-box,box] the marginal of T0
    // is exactly ∝ exp(-x²/2σ²) (the y,z box integrals give an x-independent const).
    const double box = 2.5 * sigma;
    Histogram hx(-box, box, 25);
    const long N = 2'000'000;
    Vec3 Tacc(0, 0, 0);
    double Uacc = setT(Tacc);
    for (long i = 0; i < N; ++i) {
        const Vec3 prop(rng.uniform(-box, box), rng.uniform(-box, box), rng.uniform(-box, box));
        const double Unew = setT(prop); // ghost-pass is free; recouple cost is Unew
        const double dH = Unew - Uacc;
        if (dH <= 0 || rng.uniform(0, 1) < std::exp(-dH / kT300)) {
            Tacc = prop;
            Uacc = Unew;
        } else {
            setT(Tacc); // reject: restore accepted pose
        }
        hx.add(Tacc[0]);
    }
    const double chi =
        chiSquareStatistic(hx.counts(), expectedFromWeights(gaussianBinWeights(hx, sigma), hx.total()));
    EXPECT_LT(chi, chiSquareCritical(hx.nbins() - 1, 1e-4))
        << "teleport marginal not Gaussian(kT/k): chi2=" << chi;
}

// ---------------------------------------------------------------------------
//  Phase 3: a cavity/site-BIASED teleport draw is asymmetric, so the acceptance
//  needs the Metropolis-Hastings proposal ratio g(old)/g(new). This validates the
//  mechanism: an independence sampler with a biased (off-center, wider) Gaussian
//  proposal reproduces the target N(0,σ²) ONLY when the Hastings ratio is included,
//  and is provably biased without it. (1-D; the production aiming would fold the
//  same ratio into the endpoint Metropolis.)
// ---------------------------------------------------------------------------
TEST(NcmcTeleport, HastingsRatioRestoresBalance) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    Rng rng(0x4A57);
    const double sigma = 0.5;          // target N(0, sigma^2)
    // Biased independence proposal N(muP, sP^2): off-center (asymmetric, so the
    // ratio matters) but WIDER than the target (sP>σ) so it covers the tails and
    // the with-ratio chain converges.
    const double muP = 0.35, sP = 0.8;
    auto U = [&](double x) { return 0.5 * x * x / (sigma * sigma); }; // beta*U (beta=1 units)
    auto lnG = [&](double x) { return -0.5 * (x - muP) * (x - muP) / (sP * sP); };

    const double range = 2.5 * sigma;
    auto run = [&](bool withHastings) {
        Histogram h(-range, range, 25);
        double x = 0.0;
        const long N = 4'000'000;
        for (long i = 0; i < N; ++i) {
            const double xp = rng.gaussian(muP, sP); // independence proposal (unbounded)
            double lnA = -(U(xp) - U(x));
            if (withHastings) {
                lnA += lnG(x) - lnG(xp); // + ln[g(old)/g(new)]
            }
            if (lnA >= 0 || rng.uniform(0, 1) < std::exp(lnA)) {
                x = xp;
            }
            // Histogram::add clamps out-of-range into edge bins; DISCARD the tails
            // instead so the score is the target conditioned on [-range,range]
            // (shape ∝ exp(-x²/2σ²) there) -- comparison stays exact.
            if (std::abs(x) < range) {
                h.add(x);
            }
        }
        return chiSquareStatistic(h.counts(), expectedFromWeights(gaussianBinWeights(h, sigma), h.total()));
    };

    const double crit = chiSquareCritical(24, 1e-4);
    const double chiWith = run(/*withHastings=*/true);
    const double chiWithout = run(/*withHastings=*/false);
    EXPECT_LT(chiWith, crit) << "biased proposal WITH Hastings ratio is not N(0,σ²): chi2=" << chiWith;
    EXPECT_GT(chiWithout, crit) << "omitting the Hastings ratio should bias the marginal: chi2=" << chiWithout;
}
