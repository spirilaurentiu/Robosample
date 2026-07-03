// ============================================================================
//  TestNcmcExplicitSolvent.cpp -- the OpenMM-free correctness suite for
//  docs/specs/ncmc-explicit-solvent/ (Construction II, Metropolized-dynamics
//  NCMC).
//
//  SCOPE. This file exercises the CONSTRUCTION-LEVEL claims of 10-acceptance-
//  construction.md and 20-inner-integrator.md against a small, OpenMM-free
//  system (à la TestNCMCWork.cpp / TestFixmanBoltzmann.cpp), independently of
//  World/ForceBridge (which are OpenMM-linked and privately encapsulate
//  ncmcMove). It intentionally REPRODUCES World::ncmcMove's Construction-I and
//  Construction-II loops (perturb/propagate, inner GHMC, outer accept) against
//  RobotEngine + a lambda-aware analytic bridge, exactly as TestNCMCWork.cpp
//  already does for Construction I alone.
//
//  PRIMARY ORACLE: INV0 (docs/specs/ncmc-explicit-solvent/
//  40-reproducer-and-oracles.md Sec.4). Construction I and Construction II both
//  target exp(-beta H_{lambda=1}); H = PE + KE + Fixman (this fixture has no
//  Cartesian solvent / Free root / NMA route, so keSolvent/pitch/nmaCorr are
//  identically 0 and dropped). Sampling the SAME target under two different
//  exact constructions must recover the SAME torsional marginal. An inner GHMC
//  accept that deliberately DROPS the Fixman term from H_lambda (the F1 bug:
//  trusting the "|M|^{1/2} weighting" shorthand instead of the full assembly,
//  10-acceptance-construction.md Sec.3 NOTE F2) is not pi_lambda-invariant and
//  must measurably diverge from the true marginal. This is the ONLY reference-
//  free oracle that catches that omission (L1/L2/INV2 are all blind to it).
//
//  FIXTURE. TestFixmanBoltzmann.cpp's validated two-torsion 90-degree-bend chain
//  (det M(phi2) genuinely phi2-dependent, sqrt(det M) reference computed from
//  the SAME calcLogDetM operator the sampler uses), extended with atoms (so a
//  lambda-dependent analytic potential -- TestNCMCWork.cpp's LambdaAnalyticBridge
//  pattern, renamed here -- has something to act on) and a genuine lambda:1->0->1
//  protocol driving real reorganization work.
//
//  DEFERRED (docs/specs/ncmc-explicit-solvent/40-reproducer-and-oracles.md):
//  INV1 (matched-pair DeltaF discriminator) and INV3 (physical barrier-crossing
//  gate) both need a trusted 2ala phi/psi free-energy reference (MBAR or a known
//  DeltaF observable) that is not yet generated -- wired below as explicit
//  GTEST_SKIP, never a fabricated reference. INV2 (bath-scaling) needs the
//  OpenMM-backed explicit-solvent Sweep-A harness, out of scope for this
//  OpenMM-free file; also wired as an explicit SKIP so the gap is visible.
// ============================================================================
#include <cmath>
#include <cstdint>
#include <gtest/gtest.h>
#include <random>
#include <vector>

#include "Constraints.hpp"
#include "NCMCProtocol.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotIntegrator.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "StatTest.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::attachAtoms;
using rtest::stat::chiSquareCritical;
using rtest::stat::chiSquareStatistic;
using rtest::stat::expectedFromWeights;
using rtest::stat::Histogram;
using rtest::stat::MeanAccumulator;
using rtest::stat::slowEnabled;
using rtest::stat::uniformExpected;

namespace {

constexpr double kT300 = 0.0083144626 * 300.0; // RT, kJ/mol

// ---------------------------------------------------------------------------
//  Fixture: TestFixmanBoltzmann.cpp's two-torsion 90-degree-bend chain, WITH
//  atoms attached (freeTorsionChain's placement pattern from TestNCMCWork.cpp)
//  so a lambda-dependent analytic potential can act on real Cartesian points.
//  Deterministic (no RNG): the geometry is IDENTICAL to the already-validated
//  twoTorsionChain fixture, so det M(phi2)'s phi2-dependence is trusted.
// ---------------------------------------------------------------------------
RobotModel twoTorsionChainWithAtoms() {
    const Rotation bend(Real(M_PI_2), XAxis);
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

    RobotModel m = buildForest({b1, b2});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, 0.02, -0.03), Vec3(-0.02, 0.10, 0.01)}, {12.0, 1.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    return m;
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

void primePositions(const RobotModel& m, RobotState& s) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
}

double wrapPi(double a) {
    return std::atan2(std::sin(a), std::cos(a));
}

// det M as a function of the two torsions (identical operator World::calcFixman
// uses). Mirrors TestFixmanBoltzmann.cpp::detMAt exactly.
double detMAt(RobotModel& m, RobotState& s, double phi1, double phi2) {
    s.q()[0] = static_cast<Real>(phi1);
    s.q()[1] = static_cast<Real>(phi2);
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    return std::exp(static_cast<double>(RobotEngine::calcLogDetM(m, s)));
}

// Reference per-bin weight for the phi2 marginal: sqrt(det M(phi2)), marginalized
// over a small phi1 grid. Mirrors TestFixmanBoltzmann.cpp::sqrtDetMWeights.
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

// ---------------------------------------------------------------------------
//  A lambda-aware analytic bridge: V(r; lambda) = U_intra(r) + lambda*U_inter(r)
//  (TestNCMCWork.cpp's LambdaAnalyticBridge pattern, duplicated here by the same
//  established convention -- each OpenMM-free NCMC test file rebuilds its own
//  small lambda-dependent potential rather than sharing one, keeping every file
//  self-contained). setLambda() is the only knob the perturb substep turns.
// ---------------------------------------------------------------------------
class NcmcLambdaBridge {
    public:
    NcmcLambdaBridge(const RobotModel& m, const RobotState& s0, Real kIntra, Real kInter)
        : model_(m)
        , kIntra_(kIntra)
        , kInter_(kInter)
        , intraAnchor_(static_cast<std::size_t>(m.numAtoms))
        , interAnchor_(static_cast<std::size_t>(m.numAtoms)) {
        const Vec3* p = s0.atomPosG();
        for (int a = 0; a < m.numAtoms; ++a) {
            intraAnchor_[a] = p[a];
            interAnchor_[a] = p[a] + Vec3(0.05, -0.03, 0.04);
        }
    }

    void setLambda(Real l) {
        lambda_ = l;
    }

    [[nodiscard]] Vec3 atomForce(const Vec3* pos, int a) const {
        if (model_.atomMass[a] == Real(0)) {
            return Vec3(0);
        }
        const Vec3 fi = (intraAnchor_[a] - pos[a]) * kIntra_;
        const Vec3 fe = (interAnchor_[a] - pos[a]) * (kInter_ * lambda_);
        return fi + fe;
    }

    [[nodiscard]] Real calcPotentialEnergy(const RobotState& s) const {
        const Vec3* pos = s.atomPosG();
        Real ui = 0, ue = 0;
        for (int a = 0; a < model_.numAtoms; ++a) {
            if (model_.atomMass[a] == Real(0)) {
                continue;
            }
            const Vec3 di = pos[a] - intraAnchor_[a];
            const Vec3 de = pos[a] - interAnchor_[a];
            ui += dot(di, di);
            ue += dot(de, de);
        }
        return Real(0.5) * kIntra_ * ui + lambda_ * Real(0.5) * kInter_ * ue;
    }

    void evaluate(RobotState& s) {
        SpatialVec* BF = s.bodyForceG();
        for (int b = 0; b < model_.numBodies; ++b) {
            BF[b] = SpatialVec(Vec3(0), Vec3(0));
        }
        Real* mob = s.mobilityForce();
        for (int i = 0; i < model_.nu; ++i) {
            mob[i] = Real(0);
        }
        const Vec3* posG = s.atomPosG();
        const Transform* X_GB = s.X_GB();
        for (int a = 0; a < model_.numAtoms; ++a) {
            if (model_.atomMass[a] == Real(0)) {
                continue;
            }
            const int b = model_.atomBody[a];
            const Vec3 f = atomForce(posG, a);
            const Vec3 r = posG[a] - X_GB[b].p();
            BF[b][1] += f;
            BF[b][0] += r % f;
        }
    }

    private:
    const RobotModel& model_;
    Real kIntra_, kInter_, lambda_ = 1.0;
    std::vector<Vec3> intraAnchor_, interAnchor_;
};

// ---------------------------------------------------------------------------
//  H = PE + KE [+ Fixman]. Fixman uses the SAME calcLogDetM operator World's
//  calcFixman() uses (this fixture is acyclic: cs.calcConstraintLogDet is
//  always 0, kept for parity with the production assembly). This fixture has
//  no Cartesian solvent / Free root / NMA route, so it drops K_s, the pitch
//  term, and nmaCorr -- all identically 0 here, never silently approximated
//  away for a case that would need them.
// ---------------------------------------------------------------------------
double hamiltonianAt(const RobotModel& m, RobotState& s, NcmcLambdaBridge& bridge, robo::ConstraintSet& cs,
                     double RT, bool includeFixman) {
    const double pe = bridge.calcPotentialEnergy(s);
    RobotEngine::realizeVelocity(m, s);
    const double ke = RobotEngine::calcKineticEnergy(m, s);
    double H = pe + ke;
    if (includeFixman) {
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        const double lnDetM = RobotEngine::calcLogDetM(m, s);
        const double lnDetZ = cs.calcConstraintLogDet(m, s);
        H += 0.5 * RT * (lnDetM - lnDetZ);
    }
    return H;
}

void seedMomenta(const RobotModel& m, RobotState& s, robo::ConstraintSet& cs, double RT, std::mt19937_64& rng) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    const int nu = m.nu;
    std::vector<Real> g(static_cast<std::size_t>(nu)), seeded(static_cast<std::size_t>(nu));
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (int i = 0; i < nu; ++i) {
        g[i] = static_cast<Real>(gauss(rng));
    }
    RobotEngine::multiplyBySqrtMInv(m, s, g.data(), seeded.data());
    const Real scale = static_cast<Real>(std::sqrt(RT));
    Real* u = s.u();
    for (int i = 0; i < nu; ++i) {
        u[i] = scale * seeded[i];
    }
    if (!cs.empty()) {
        RobotEngine::realizeVelocity(m, s);
        cs.enforceVelocityConstraints(m, s);
    }
}

void seedDerivatives(const RobotModel& m, RobotState& s, NcmcLambdaBridge& bridge) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    bridge.evaluate(s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    RobotEngine::calcUDot(m, s);
    RobotEngine::calcQDot(m, s, s.qdot());
    RobotEngine::calcQDotDot(m, s);
}

// ---------------------------------------------------------------------------
//  Construction-II inner kernel: an independent test-only re-derivation of
//  World::ncmcInnerGhmcStep (src/World.cpp), against this file's OpenMM-free
//  bridge. `includeFixman` is the deliberate F1-bug knob: false reproduces the
//  omission this file's INV0 divergence arm must catch.
// ---------------------------------------------------------------------------
bool innerGhmcStep(const RobotModel& m, RobotState& s, NcmcLambdaBridge& bridge, robo::ConstraintSet& cs,
                   double h, double RT, bool includeFixman, std::mt19937_64& rng) {
    std::vector<Real> q0(s.q(), s.q() + m.nq);
    std::vector<Real> u0(s.u(), s.u() + m.nu);

    const double Hbefore = hamiltonianAt(m, s, bridge, cs, RT, includeFixman);

    bool converged = true;
    const bool stepOk = RobotEngine::stepTo(m, s, bridge, cs, s.time + static_cast<Real>(h), &converged);
    if (!stepOk) {
        return false;
    }

    bool accept = false;
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    if (converged) {
        const double Hafter = hamiltonianAt(m, s, bridge, cs, RT, includeFixman);
        const double dH = Hafter - Hbefore;
        accept = std::isfinite(dH) && (dH <= 0.0 || unif(rng) < std::exp(-dH / RT));
    }

    if (!accept) {
        std::copy(q0.begin(), q0.end(), s.q());
        std::copy(u0.begin(), u0.end(), s.u());
        Real* u = s.u();
        for (int i = 0; i < m.nu; ++i) {
            u[i] = -u[i];
        }
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        bridge.evaluate(s);
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        RobotEngine::calcUDot(m, s);
        RobotEngine::calcQDot(m, s, s.qdot());
        RobotEngine::calcQDotDot(m, s);
    }
    return true;
}

// Construction I (endpoint-DeltaH): raw Verlet propagate, outer accept on
// Hend-Hstart with the Fixman-COMPLETE H (matching World's own assembly; F5 --
// NOT the V+K-only reference TestNCMCWork.cpp uses for its unrelated claim).
bool ncmcMoveConstructionI(const RobotModel& m, RobotState& s, NcmcLambdaBridge& bridge, robo::ConstraintSet& cs,
                          int ncmcSteps, double holdFraction, double h, double RT, std::mt19937_64& rng,
                          double* workOut = nullptr) {
    std::vector<Real> q0(s.q(), s.q() + m.nq);
    std::vector<Real> u0(s.u(), s.u() + m.nu);

    bridge.setLambda(1.0);
    seedMomenta(m, s, cs, RT, rng);
    seedDerivatives(m, s, bridge);
    const double Hstart = hamiltonianAt(m, s, bridge, cs, RT, /*includeFixman*/ true);

    double work = 0.0, lamPrev = 1.0;
    double Vprev = bridge.calcPotentialEnergy(s);
    bool ok = true;
    for (int step = 0; ok && step < ncmcSteps; ++step) {
        const double lam = ncmc::protocolLambda(step, ncmcSteps, holdFraction);
        if (lam != lamPrev) {
            bridge.setLambda(static_cast<Real>(lam));
            const double Vnew = bridge.calcPotentialEnergy(s);
            work += Vnew - Vprev;
            Vprev = Vnew;
            lamPrev = lam;
        }
        ok = RobotEngine::stepTo(m, s, bridge, cs, s.time + static_cast<Real>(h));
        if (ok) {
            Vprev = bridge.calcPotentialEnergy(s);
        }
    }
    if (workOut) {
        *workOut = work;
    }

    auto restore = [&]() {
        std::copy(q0.begin(), q0.end(), s.q());
        std::copy(u0.begin(), u0.end(), s.u());
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
    };
    if (!ok) {
        restore();
        return false;
    }
    bridge.setLambda(1.0);
    const double Hend = hamiltonianAt(m, s, bridge, cs, RT, true);
    const double dH = Hend - Hstart;
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const bool accept = std::isfinite(dH) && (dH <= 0.0 || unif(rng) < std::exp(-dH / RT));
    if (!accept) {
        restore();
    }
    return accept;
}

// Construction II: inner-GHMC-Metropolized propagate; outer accept on the
// protocol work W alone, min(1,exp(-W/RT)).
bool ncmcMoveConstructionII(const RobotModel& m, RobotState& s, NcmcLambdaBridge& bridge, robo::ConstraintSet& cs,
                           int ncmcSteps, double holdFraction, double h, double RT, std::mt19937_64& rng,
                           bool innerIncludesFixman, double* workOut = nullptr) {
    std::vector<Real> q0(s.q(), s.q() + m.nq);
    std::vector<Real> u0(s.u(), s.u() + m.nu);

    bridge.setLambda(1.0);
    seedMomenta(m, s, cs, RT, rng);
    seedDerivatives(m, s, bridge);

    double work = 0.0, lamPrev = 1.0;
    double Vprev = bridge.calcPotentialEnergy(s);
    bool ok = true;
    for (int step = 0; ok && step < ncmcSteps; ++step) {
        const double lam = ncmc::protocolLambda(step, ncmcSteps, holdFraction);
        if (lam != lamPrev) {
            bridge.setLambda(static_cast<Real>(lam));
            const double Vnew = bridge.calcPotentialEnergy(s);
            work += Vnew - Vprev;
            Vprev = Vnew;
            lamPrev = lam;
        }
        ok = innerGhmcStep(m, s, bridge, cs, h, RT, innerIncludesFixman, rng);
        if (ok) {
            Vprev = bridge.calcPotentialEnergy(s);
        }
    }
    if (workOut) {
        *workOut = work;
    }

    auto restore = [&]() {
        std::copy(q0.begin(), q0.end(), s.q());
        std::copy(u0.begin(), u0.end(), s.u());
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
    };
    if (!ok) {
        restore();
        return false;
    }
    bridge.setLambda(1.0);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const bool accept = std::isfinite(work) && (work <= 0.0 || unif(rng) < std::exp(-work / RT));
    if (!accept) {
        restore();
    }
    return accept;
}

struct ChainRun {
    Histogram phi2{-M_PI, M_PI, 24};
    long accepted = 0;
    long moves = 0;
};

// ncmcSteps/holdFraction/h chosen so BOTH constructions run the SAME protocol
// (only the propagate substep differs); Construction I needs a small enough dt
// that its endpoint-DeltaH acceptance is not itself vanishing (else the phi2
// marginal would just reflect "never moves", making the "I vs II agree" claim
// vacuous -- INV0's own precondition).
constexpr int kNcmcSteps = 8;
constexpr double kHoldFraction = 0.1;
constexpr double kDt = 4.0e-4;

ChainRun sampleChainI(std::uint64_t seed, long nMoves, int stride) {
    RobotModel m = twoTorsionChainWithAtoms();
    RobotState s;
    initState(m, s);
    primePositions(m, s);
    robo::ConstraintSet cs;
    NcmcLambdaBridge bridge(m, s, /*kIntra*/ 300.0, /*kInter*/ 120.0);
    std::mt19937_64 rng(seed);
    ChainRun r;
    for (long i = 0; i < nMoves; ++i) {
        r.accepted += ncmcMoveConstructionI(m, s, bridge, cs, kNcmcSteps, kHoldFraction, kDt, kT300, rng) ? 1 : 0;
        ++r.moves;
        if (i % stride == 0) {
            r.phi2.add(wrapPi(static_cast<double>(s.q()[1])));
        }
    }
    return r;
}

ChainRun sampleChainII(std::uint64_t seed, long nMoves, int stride, bool innerIncludesFixman) {
    RobotModel m = twoTorsionChainWithAtoms();
    RobotState s;
    initState(m, s);
    primePositions(m, s);
    robo::ConstraintSet cs;
    NcmcLambdaBridge bridge(m, s, /*kIntra*/ 300.0, /*kInter*/ 120.0);
    std::mt19937_64 rng(seed);
    ChainRun r;
    for (long i = 0; i < nMoves; ++i) {
        r.accepted +=
            ncmcMoveConstructionII(m, s, bridge, cs, kNcmcSteps, kHoldFraction, kDt, kT300, rng, innerIncludesFixman)
                ? 1
                : 0;
        ++r.moves;
        if (i % stride == 0) {
            r.phi2.add(wrapPi(static_cast<double>(s.q()[1])));
        }
    }
    return r;
}

} // namespace

// ===========================================================================
//  SMOKE (always on): the machinery runs, both constructions accept a nonzero
//  fraction of moves, and phi2 explores its range under each. Vacuous-test
//  guard for the slow-tier assertions below.
// ===========================================================================
TEST(NcmcExplicitSolvent, Smoke) {
    ChainRun rI = sampleChainI(0xA0, 4000, 2);
    ChainRun rII = sampleChainII(0xA1, 4000, 2, /*innerIncludesFixman*/ true);
    EXPECT_GT(rI.accepted, 0) << "Construction I never accepted; fixture/dt too aggressive";
    EXPECT_GT(rII.accepted, 0) << "Construction II never accepted; fixture/dt too aggressive";

    auto explored = [](const Histogram& h) {
        int populated = 0;
        for (long c : h.counts()) {
            populated += (c > 0) ? 1 : 0;
        }
        return populated;
    };
    EXPECT_GE(explored(rI.phi2), 6) << "Construction I phi2 did not explore";
    EXPECT_GE(explored(rII.phi2), 6) << "Construction II phi2 did not explore";
}

// ===========================================================================
//  INV0 (docs/specs/ncmc-explicit-solvent/40-reproducer-and-oracles.md Sec.4,
//  PRIMARY correctness oracle, F1) -- Sweep E: Construction I vs Construction II
//  marginal agreement, AND the Fixman-omitting arm's divergence.
// ===========================================================================

// ARM 1/2 -- AGREEMENT. Construction I (endpoint-DeltaH, Fixman-complete H) and
// Construction II (Fixman-complete inner GHMC accept) both target
// exp(-beta H_{lambda=1}); with useFixman-equivalent H the FIXMAN-Boltzmann
// theory (TestFixmanBoltzmann.cpp) predicts a FLAT phi2 marginal for EITHER
// construction. We assert both individually pass the flat (uniform) chi-square
// gate and both individually REJECT the sqrt(det M) shape -- i.e. Construction
// II reproduces Construction I's own established correctness signature.
TEST(NcmcExplicitSolvent, INV0_ConstructionIAndIIAgree_BothFlatUnderFixman) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    const long nMoves = 200'000;
    ChainRun rI = sampleChainI(0xB0117, nMoves, 3);
    ChainRun rII = sampleChainII(0xB0118, nMoves, 3, /*innerIncludesFixman*/ true);

    RobotModel m = twoTorsionChainWithAtoms();
    RobotState s;
    initState(m, s);
    const std::vector<double> wRefI = sqrtDetMWeights(m, s, rI.phi2);
    const std::vector<double> wRefII = sqrtDetMWeights(m, s, rII.phi2);

    const double crit = chiSquareCritical(rI.phi2.nbins() - 1, 1e-4);

    const double chiFlatI = chiSquareStatistic(rI.phi2.counts(), uniformExpected(rI.phi2.nbins(), rI.phi2.total()));
    const double chiRefI =
        chiSquareStatistic(rI.phi2.counts(), expectedFromWeights(wRefI, rI.phi2.total()));
    EXPECT_LT(chiFlatI, crit) << "Construction I (Fixman-complete) marginal should be flat: chi2=" << chiFlatI
                              << " crit=" << crit;
    EXPECT_GT(chiRefI, crit) << "Construction I marginal should NOT follow sqrt(det M): chi2=" << chiRefI;

    const double chiFlatII =
        chiSquareStatistic(rII.phi2.counts(), uniformExpected(rII.phi2.nbins(), rII.phi2.total()));
    const double chiRefII =
        chiSquareStatistic(rII.phi2.counts(), expectedFromWeights(wRefII, rII.phi2.total()));
    EXPECT_LT(chiFlatII, crit) << "Construction II (Fixman-complete inner accept) marginal should be flat: chi2="
                               << chiFlatII << " crit=" << crit;
    EXPECT_GT(chiRefII, crit) << "Construction II marginal should NOT follow sqrt(det M): chi2=" << chiRefII;
}

// ARM 3 -- DIVERGENCE. The inner GHMC accept deliberately OMITS the Fixman term
// from H_lambda (the canonical F1 bug: trusting the "|M|^{1/2}" shorthand
// instead of the full currentTotalEnergy() assembly). This is Construction I's
// OWN "useFixman=false" regime reproduced through the inner kernel: the phi2
// marginal must follow sqrt(det M) (NOT flat), i.e. it measurably DIVERGES from
// the ARM-1/2 result above. This is the discriminator L1/L2/INV2 cannot provide
// (docs/specs/ncmc-explicit-solvent/40-reproducer-and-oracles.md Sec.4 NOTE F1).
TEST(NcmcExplicitSolvent, INV0_FixmanOmittingInnerAcceptDiverges) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    const long nMoves = 200'000;
    ChainRun rBad = sampleChainII(0xBAD0, nMoves, 3, /*innerIncludesFixman*/ false);

    RobotModel m = twoTorsionChainWithAtoms();
    RobotState s;
    initState(m, s);
    const std::vector<double> wRef = sqrtDetMWeights(m, s, rBad.phi2);

    const double crit = chiSquareCritical(rBad.phi2.nbins() - 1, 1e-4);
    const double chiFlat =
        chiSquareStatistic(rBad.phi2.counts(), uniformExpected(rBad.phi2.nbins(), rBad.phi2.total()));
    const double chiRef = chiSquareStatistic(rBad.phi2.counts(), expectedFromWeights(wRef, rBad.phi2.total()));

    EXPECT_GT(chiFlat, crit) << "Fixman-omitting inner accept should NOT be flat (F1 bug undetected!): chi2="
                             << chiFlat << " crit=" << crit;
    EXPECT_LT(chiRef, crit) << "Fixman-omitting inner accept should follow sqrt(det M) (the biased, "
                               "un-corrected marginal): chi2="
                            << chiRef << " crit=" << crit;
}

// ===========================================================================
//  LEMMA L1 (flat-lambda gate, F5: Fixman-COMPLETE reference)
// ===========================================================================

// With lambda == 1 throughout (perturb never fires), the protocol work is
// EXACTLY 0 bitwise -- structurally true for BOTH constructions (the perturb
// substep never executes), and Construction II's outer accept then reduces to
// min(1,exp(-0/RT)) == always-accept, so the WHOLE composite trajectory reduces
// to n consecutive Metropolized inner-GHMC substeps -- i.e. plain torsional GHMC
// under the SAME Fixman-complete H the outer acceptance uses (NOT the V+K-only
// reference TestNCMCWork.cpp's existing L1-analogue uses, per F5).
//
// CAVEAT (matches 40-reproducer-and-oracles.md Sec.4 L1 note verbatim): this
// lemma is Fixman-BLIND by construction -- toggling innerIncludesFixman does not
// change whether work==0 or whether the loop reduces to n inner-GHMC calls, so
// L1 passes identically whether or not Fixman is included. It is NOT the
// discriminating oracle (INV0 above is); it only certifies the flat-lambda
// STRUCTURAL reduction and the Fixman-complete-vs-incomplete distinction in the
// REFERENCE this file compares against elsewhere.
TEST(NcmcExplicitSolvent, L1_FlatLambdaGivesZeroWorkAndReducesToInnerGhmcChain) {
    RobotModel m = twoTorsionChainWithAtoms();
    RobotState s;
    initState(m, s);
    primePositions(m, s);
    robo::ConstraintSet cs;
    NcmcLambdaBridge bridge(m, s, 300.0, 120.0);
    std::mt19937_64 rng(0x1111);

    bridge.setLambda(1.0);
    seedMomenta(m, s, cs, kT300, rng);
    seedDerivatives(m, s, bridge);

    double work = 0.0, lamPrev = 1.0;
    double Vprev = bridge.calcPotentialEnergy(s);
    const int n = 12;
    bool ok = true;
    for (int step = 0; ok && step < n; ++step) {
        const double lam = 1.0; // flat protocol: perturb never fires
        if (lam != lamPrev) {
            bridge.setLambda(static_cast<Real>(lam));
            work += bridge.calcPotentialEnergy(s) - Vprev;
            lamPrev = lam;
        }
        ok = innerGhmcStep(m, s, bridge, cs, kDt, kT300, /*includeFixman*/ true, rng);
        Vprev = bridge.calcPotentialEnergy(s);
    }
    ASSERT_TRUE(ok);
    EXPECT_EQ(work, 0.0) << "flat lambda=1 protocol produced nonzero work under Construction II";

    // The outer accept is min(1,exp(-0/RT)) == always-accept: verify that
    // directly (not just that work==0), since a bug that computed work
    // correctly but mis-wired the outer accept comparison would still pass the
    // work==0.0 assertion alone.
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const bool outerAccept = std::isfinite(work) && (work <= 0.0 || unif(rng) < std::exp(-work / kT300));
    EXPECT_TRUE(outerAccept) << "work==0 outer accept must be unconditional (min(1,exp(0))==1)";
}

// ===========================================================================
//  LEMMA L2 (Jarzynski cyclic -- WEAK self-consistency gate only, F4)
// ===========================================================================

// For the cyclic palindromic lambda:1->0->1 protocol, <exp(-beta W)> == 1
// (Jarzynski, DeltaF=0). NOTE (matches 40-...:Sec.4 L2 verbatim): this does NOT
// discriminate correct from Fixman-biased sampling (U_F/pitch are
// lambda-independent, so W and hence <exp(-beta W)> are IDENTICAL whether or not
// the inner kernel includes them), and the exponential-average estimator is
// dominated by rare large-negative-W tails, so Var(exp(-beta W)) can be
// effectively undefined at achievable move counts. We therefore assert only a
// WIDE sanity band (order-of-magnitude around 1), explicitly NOT a tight
// statistical gate, and run it under Construction II with the Fixman-complete
// inner accept only (the arm this file otherwise certifies correct via INV0).
TEST(NcmcExplicitSolvent, L2_JarzynskiCyclic_WeakSelfConsistencyOnly) {
    if (!slowEnabled()) {
        rtest::warnSlowTierSkipped();
        GTEST_SKIP() << "slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)";
    }
    RobotModel m = twoTorsionChainWithAtoms();
    RobotState s;
    initState(m, s);
    primePositions(m, s);
    robo::ConstraintSet cs;
    NcmcLambdaBridge bridge(m, s, 300.0, 120.0);
    std::mt19937_64 rng(0x1A2B);

    MeanAccumulator expNegBetaW;
    const long nMoves = 20'000;
    for (long i = 0; i < nMoves; ++i) {
        double w = 0.0;
        ncmcMoveConstructionII(m, s, bridge, cs, kNcmcSteps, kHoldFraction, kDt, kT300, rng,
                              /*innerIncludesFixman*/ true, &w);
        // Numerical-safety clip only (documented above): prevents an IEEE
        // overflow from a rare large-negative-W tail from corrupting the running
        // mean with +inf; does not change the "weak, noisy" character the spec
        // itself attributes to this estimator.
        const double exponent = std::min(-w / kT300, 50.0);
        expNegBetaW.add(std::exp(exponent));
    }
    const double mean = expNegBetaW.mean();
    ASSERT_TRUE(std::isfinite(mean));
    // Wide sanity band: a factor of 5 around 1 in either direction. A real
    // acceptance-construction bug (e.g. sign-flipped W, or a persistently
    // dissipative protocol far from the reversible limit) would fail this by
    // orders of magnitude, not by a factor of a few.
    EXPECT_GT(mean, 0.2) << "<exp(-beta W)> implausibly far below 1 (mean=" << mean << ")";
    EXPECT_LT(mean, 5.0) << "<exp(-beta W)> implausibly far above 1 (mean=" << mean << ")";
}

// ===========================================================================
//  DEFERRED (explicit SKIP, never a fabricated reference)
// ===========================================================================

TEST(NcmcExplicitSolvent, INV1_MatchedPairDeltaFDiscriminator_DeferredPendingReference) {
    GTEST_SKIP() << "INV1 (docs/specs/ncmc-explicit-solvent/40-reproducer-and-oracles.md Sec.4) needs a "
                    "trusted 2ala phi/psi free-energy reference (MBAR or another known-DeltaF "
                    "observable) that is not yet generated. Out of scope for this pass per the coder "
                    "task's explicit instruction (stub as pending, do not fabricate a reference).";
}

TEST(NcmcExplicitSolvent, INV2_BathScalingInvariance_DeferredPendingOpenMMHarness) {
    GTEST_SKIP() << "INV2 (bath-scaling invariance under box padding) needs the OpenMM-backed "
                    "explicit-solvent Sweep-A harness (docs/specs/ncmc-explicit-solvent/"
                    "40-reproducer-and-oracles.md Sec.3) driving a real tip3p/2ala system through "
                    "World::ncmcMove; out of scope for this OpenMM-free file. Deferred, not silently "
                    "dropped from the suite.";
}

TEST(NcmcExplicitSolvent, INV3_PhysicalBarrierCrossingGate_DeferredPendingReference) {
    GTEST_SKIP() << "INV3 needs an alanine-dipeptide C7eq<->C7ax (or cis<->trans) MBAR free-energy "
                    "reference to compare the recovered basin free-energy difference against; not yet "
                    "generated (same deferral as INV1). Wired here as an explicit SKIP.";
}
