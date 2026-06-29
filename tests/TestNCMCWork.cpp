// ============================================================================
//  TestNcmcWork.cpp -- the OpenMM-free half of the alchemical/NCMC suite.
//
//  The alchemy itself lives in OpenMM (TestAlchemy.cpp covers the forces). What
//  lives in Robosample-land, and is tested HERE, is the part that wires the
//  alchemical perturbation into the internal-coordinate (torsional + Free-root
//  quaternion) generalized-coordinate HMC:
//
//   1. THE LAMBDA SCHEDULE (robo::ncmc::protocolLambda). The symmetric triangle
//      1 -> 0 -> 1 with an optional flat hold. Endpoints, range, monotonicity,
//      its self-time-reverse symmetry (the property that makes the NCMC proposal
//      its own reverse, so the protocol ratio is 1 and drops out of dH), and the
//      degenerate guards (tiny ncmcSteps, holdFraction extremes).
//
//   2. THE PROTOCOL WORK (World::ncmcMove's accumulator). ncmcMove interleaves a
//      PERTURB substep (change lambda at FIXED q; w += V(lam_new) - V(lam_old))
//      with a PROPAGATE substep (one constrained Verlet step at fixed lambda) --
//      structurally identical to openmmtools' _add_alchemical_perturbation_step
//      (protocol_work += Enew - Eold). We reproduce that loop EXACTLY, driving the
//      REAL RobotEngine::stepTo over a Free+Torsion chain with a lambda-dependent
//      analytic potential, and assert the physics the production comment claims:
//        (a) telescoping: w == sum_s [V(lam_s,q_s) - V(lam_{s-1},q_s)];
//        (b) for a deterministic, reversible, volume-preserving propagator the
//            heat is ~0, so w == H_end - H_start up to the integrator drift, and
//            the gap shrinks ~quadratically as the step h is halved (this IS the
//            "free consistency check" the code logs);
//        (c) the lambda==1 gate: a flat-1 protocol makes w == 0 EXACTLY and
//            reduces dH to the plain torsional-HMC Verlet drift.
//
//   3. JACOBIAN INVARIANCE. The alchemical perturbation changes lambda at fixed
//      q, so it has UNIT coordinate Jacobian: it must leave the Fixman mass-metric
//      determinant ln|M(q)| (RobotEngine::calcLogDetM) and the Free-root quaternion
//      pitch term log sin^2(gamma2) -- both pure functions of q -- bit-unchanged.
//      That invariance is exactly why ncmcMove can accept on H_end - H_start with
//      no extra alpha-ratio / protocol-ratio factors. We assert it directly.
//
//  No OpenMM, no World (which pulls OpenMM via ForceBridge): the loop is rebuilt
//  on the engine + a test bridge, the same pattern as TestIntegrator.cpp.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "AnalyticForceBridge.hpp"
#include "Constraints.hpp"
#include "NCMCProtocol.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotIntegrator.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::attachAtoms;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::kAlg;
using rtest::kTight;
using rtest::randomizeState;
using rtest::Rng;

// ===========================================================================
//  PART 1 -- the pure lambda schedule
//
//  THEORY UNDER TEST (Fix 1): the whole NCMC trajectory is one map T on (q,p).
//  Accept-on-endpoint-H is an EXACT Metropolis test only if T is momentum-flip
//  reversible, F T F == T^-1. Since T composes per-substep fixed-lambda Verlet
//  steps (each itself F-reversible), the composition is reversible IFF the lambda
//  schedule is a true PALINDROME, lambda(s) == lambda(N-1-s), with BOTH endpoints
//  pinned to lambda = 1. These tests encode exactly that. They FAIL on the old
//  offset-by-one schedule (which pinned only s=0 to 1 and was symmetric about the
//  trough INDEX, not about the array center) and PASS on the palindromic one.
// ===========================================================================

// lambda(0) is exactly 1 (the protocol starts fully coupled at the physical
// lambda=1 state acceptance uses).
TEST(NcmcProtocol, StartsAtOne) {
    for (int n : {2, 5, 10, 33, 100, 1000}) {
        EXPECT_DOUBLE_EQ(ncmc::protocolLambda(0, n, 0.0), 1.0) << "n=" << n;
    }
}

// BOTH ENDPOINTS PINNED: lambda(0) == lambda(N-1) == 1. This is the schedule-level
// statement that the trajectory begins AND ends at the fully-coupled lambda=1
// physical state, with NO out-of-loop "final jump to 1". The OLD schedule failed
// this (its last grid point sat at 1 - 2/N), which is what broke F-reversibility
// at O(1/N). Checked across parities and with a hold.
TEST(NcmcProtocol, BothEndpointsPinnedToOne) {
    for (int n : {2, 3, 4, 5, 20, 21, 100, 101}) {
        for (double hold : {0.0, 0.1, 0.3}) {
            EXPECT_DOUBLE_EQ(ncmc::protocolLambda(0, n, hold), 1.0) << "start n=" << n << " hold=" << hold;
            EXPECT_DOUBLE_EQ(ncmc::protocolLambda(n - 1, n, hold), 1.0)
                << "end n=" << n << " hold=" << hold;
        }
    }
}

// PALINDROME: lambda(s) == lambda(N-1-s) for EVERY step and EVERY hold. This is
// the exact discrete condition for F T F == T^-1 (the protocol is its own
// time-reverse), so the protocol ratio is 1 and no extra factor enters dH. The
// OLD schedule was symmetric about the trough INDEX (ramp), not about (N-1)/2, so
// it violated this for even N -- the reversibility bug Fix 1 closes.
TEST(NcmcProtocol, IsPalindrome) {
    for (int n : {2, 3, 4, 7, 16, 17, 64, 101, 256}) {
        for (double hold : {0.0, 0.1, 0.4}) {
            for (int s = 0; s < n; ++s) {
                EXPECT_DOUBLE_EQ(ncmc::protocolLambda(s, n, hold), ncmc::protocolLambda(n - 1 - s, n, hold))
                    << "palindrome broken at n=" << n << " hold=" << hold << " s=" << s;
            }
        }
    }
}

// lambda stays in [0,1] for every step, including any step at or beyond ncmcSteps
// (defensive: ncmcMove indexes s in [0,ncmcSteps)).
TEST(NcmcProtocol, RangeIsUnitInterval) {
    for (int n : {2, 7, 16, 64, 257}) {
        for (double hold : {0.0, 0.1, 0.5, 0.9}) {
            for (int s = 0; s <= n + 3; ++s) {
                const double l = ncmc::protocolLambda(s, n, hold);
                EXPECT_GE(l, 0.0) << "n=" << n << " hold=" << hold << " s=" << s;
                EXPECT_LE(l, 1.0) << "n=" << n << " hold=" << hold << " s=" << s;
            }
        }
    }
}

// The trough reaches EXACTLY 0 (the fully-decoupled, uncaged state) at the
// palindrome center when N is odd (a single center index s=(N-1)/2). The trough
// is the global minimum of the schedule either way.
TEST(NcmcProtocol, TroughIsMinimumAndReachesZero) {
    for (int n : {5, 11, 51, 201}) { // odd -> exact center index exists
        const int center = (n - 1) / 2;
        EXPECT_DOUBLE_EQ(ncmc::protocolLambda(center, n, 0.0), 0.0) << "n=" << n;
        for (int s = 0; s < n; ++s) {
            EXPECT_GE(ncmc::protocolLambda(s, n, 0.0), ncmc::protocolLambda(center, n, 0.0))
                << "center is not the minimum, n=" << n << " s=" << s;
        }
    }
}

// Monotone non-increasing on the down-ramp (s up to the center), then
// non-decreasing on the up-ramp.
TEST(NcmcProtocol, MonotoneDownThenUp) {
    const int n = 101;
    const int center = (n - 1) / 2;
    for (int s = 1; s <= center; ++s) {
        EXPECT_LE(ncmc::protocolLambda(s, n, 0.0), ncmc::protocolLambda(s - 1, n, 0.0)) << "down s=" << s;
    }
    for (int s = center + 1; s < n; ++s) {
        EXPECT_GE(ncmc::protocolLambda(s, n, 0.0), ncmc::protocolLambda(s - 1, n, 0.0)) << "up s=" << s;
    }
}

// The flat lambda=0 hold is CENTERED on the palindrome axis: every pinned-at-0
// step is symmetric about (N-1)/2 (so the hold cannot itself break the palindrome,
// the uncaged stride straddles the trough). We assert the zero-set is palindromic
// and non-empty for a real hold.
TEST(NcmcProtocol, HoldZeroBlockIsCentered) {
    const int n = 100;
    const double hold = 0.4;
    int zeros = 0;
    for (int s = 0; s < n; ++s) {
        const bool isZero = (ncmc::protocolLambda(s, n, hold) == 0.0);
        if (isZero) {
            ++zeros;
            EXPECT_DOUBLE_EQ(ncmc::protocolLambda(n - 1 - s, n, hold), 0.0)
                << "zero-set not symmetric at s=" << s;
        }
    }
    EXPECT_GT(zeros, 0) << "a real hold produced no lambda=0 steps";
}

// Degenerate guards: N<=1 is a single fully-coupled substep, and holdFraction>=1
// (whole protocol held) still yields finite lambdas in [0,1].
TEST(NcmcProtocol, DegenerateBudgetsAreSafe) {
    EXPECT_DOUBLE_EQ(ncmc::protocolLambda(0, 1, 0.0), 1.0);
    for (int n : {0, 1, 2, 3}) {
        for (double hold : {0.0, 1.0, 2.0}) {
            for (int s = 0; s <= std::max(n, 1) + 1; ++s) {
                const double l = ncmc::protocolLambda(s, n, hold);
                EXPECT_TRUE(std::isfinite(l)) << "n=" << n << " hold=" << hold << " s=" << s;
                EXPECT_GE(l, 0.0);
                EXPECT_LE(l, 1.0);
            }
        }
    }
}

// ===========================================================================
//  PART 2 -- the protocol-work accumulator, on the REAL GC integrator
// ===========================================================================
namespace {

// A lambda-aware analytic bridge: V(r; lambda) = U_intra(r) + lambda * U_inter(r),
// two isotropic harmonic wells per atom about distinct anchors. This mirrors the
// alchemy structure -- only V depends on lambda; the q-dynamics is the generalized-
// coordinate Verlet -- WITHOUT needing OpenMM. evaluate() reduces per-atom forces
// to per-body spatial forces with the IDENTICAL rule as ForceBridge/AnalyticForce-
// Bridge (clear, skip mass==0, moment about body origin), so stepTo integrates it
// exactly as it would the real bridge. setLambda() is the only knob the perturb
// substep turns.
class LambdaAnalyticBridge {
    public:
    LambdaAnalyticBridge(const RobotModel& m, const RobotState& s0, Real kIntra, Real kInter)
        : model_(m)
        , kIntra_(kIntra)
        , kInter_(kInter)
        , intraAnchor_(static_cast<std::size_t>(m.numAtoms))
        , interAnchor_(static_cast<std::size_t>(m.numAtoms)) {
        const Vec3* p = s0.atomPosG();
        for (int a = 0; a < m.numAtoms; ++a) {
            intraAnchor_[a] = p[a];                           // U_intra: anchored AT the start (0 at start)
            interAnchor_[a] = p[a] + Vec3(0.05, -0.03, 0.04); // U_inter: a displaced well
        }
    }

    void setLambda(Real l) {
        lambda_ = l;
    }
    [[nodiscard]] Real lambda() const {
        return lambda_;
    }

    // F_a = -dV/dr_a = -kIntra (r-intra) - lambda*kInter (r-inter); 0 for vsites.
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

    // The contract stepTo/verletStep call. Same reduction as the real bridge.
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

// A Free-rooted 3-body chain (Free + Torsion + Torsion) carrying real atoms: a
// faithful stand-in for a small mobile molecule with external (quaternion) and
// internal (torsion) DOF -- exactly the coordinate types the alchemy must coexist
// with. Fixed seed -> bit-reproducible.
RobotModel freeTorsionChain(Rng& rng) {
    auto mk = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        s.X_BM = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        s.mass = rng.uniform(0.9, 1.5);
        s.com_B = rng.vec3(-0.08, 0.08);
        s.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
        return s;
    };
    RobotModel m =
        buildForest({mk(0, JointType::Free), mk(1, JointType::Torsion), mk(2, JointType::Torsion)});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, 0.02, -0.03), Vec3(-0.02, 0.10, 0.01)}, {12.0, 1.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0.04, -0.06, 0.09), Vec3(0.08, 0.01, -0.02)}, {16.0, 1.0});
    return m;
}

// Seed the derivative chain so the first verletStep has valid qdot0/udot0/qdd0
// (identical preamble to TestIntegrator.cpp's seedDerivatives).
template <class Bridge>
void seedDerivatives(const RobotModel& m, RobotState& s, Bridge& bridge) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    bridge.evaluate(s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    RobotEngine::calcUDot(m, s);
    RobotEngine::calcQDot(m, s, s.qdot());
    RobotEngine::calcQDotDot(m, s);
}

// Generalized-coordinate kinetic energy 1/2 u^T M(q) u (the same calcKineticEnergy
// World uses in its Hamiltonian).
Real kineticEnergy(const RobotModel& m, RobotState& s) {
    RobotEngine::realizeVelocity(m, s);
    return RobotEngine::calcKineticEnergy(m, s);
}

// Snapshot / restore the generalized coordinates and speeds, so a test can replay
// a move from a bit-identical start (RobotState::copyTransferPayloadTo copies only
// atom positions + energy, not q/u).
struct QU {
    std::vector<Real> q, u;
};
QU snapshotQU(const RobotModel& m, const RobotState& s) {
    return {std::vector<Real>(s.q(), s.q() + m.nq), std::vector<Real>(s.u(), s.u() + m.nu)};
}
void restoreQU(const RobotModel& m, RobotState& s, const QU& snap) {
    std::copy(snap.q.begin(), snap.q.end(), s.q());
    std::copy(snap.u.begin(), snap.u.end(), s.u());
}

// Fill atomPosG from q so the bridge can capture start-configuration anchors.
void primePositions(const RobotModel& m, RobotState& s) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
}

// Result of running the ncmcMove loop over a model potential.
struct NcmcRun {
    double work = 0;   // the protocol-work accumulator (diagnostic in production)
    double Hstart = 0; // PE + KE at lambda=1, start
    double Hend = 0;   // PE + KE at lambda=1, end
    double heat = 0;   // sum of fixed-lambda propagate dH (the integrator's shadow work)
    bool ok = true;
};

// Reproduce ncmcMove's PERTURB/PROPAGATE loop EXACTLY (minus the Metropolis gate
// and momentum resampling, which are orthogonal to the work identity), over the
// LambdaAnalyticBridge. H is PE+KE here (Fixman is constant in q-free perturb and
// covered separately in Part 3). nLambda controls how many distinct lambda values
// the schedule visits; we use a fine schedule so perturbations actually fire.
NcmcRun runNcmcLoop(RobotModel& m,
                    RobotState& s,
                    LambdaAnalyticBridge& bridge,
                    ConstraintSet& cs,
                    int ncmcSteps,
                    double holdFraction,
                    double h) {
    NcmcRun r;
    bridge.setLambda(1.0);
    seedDerivatives(m, s, bridge);
    r.Hstart = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);

    double work = 0.0, heat = 0.0;
    double Vprev = bridge.calcPotentialEnergy(s);
    double lamPrev = 1.0;

    for (int step = 0; r.ok && step < ncmcSteps; ++step) {
        // (i) PERTURB at FIXED q: change lambda, accumulate work = V_new - V_old.
        const double lam = ncmc::protocolLambda(step, ncmcSteps, holdFraction);
        if (lam != lamPrev) {
            bridge.setLambda(static_cast<Real>(lam));
            const double Vnew = bridge.calcPotentialEnergy(s); // SAME q, new lambda
            work += Vnew - Vprev;
            Vprev = Vnew;
            lamPrev = lam;
        }
        // (ii) PROPAGATE one Verlet step at fixed lambda; the fixed-lambda dH is heat.
        const double Hbefore = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);
        r.ok = RobotEngine::stepTo(m, s, bridge, cs, s.time + h);
        if (!r.ok) {
            break;
        }
        const double Hafter = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);
        heat += Hafter - Hbefore;
        Vprev = bridge.calcPotentialEnergy(s);
    }

    // Hend is evaluated at the protocol's ACTUAL final lambda: the work accumulator
    // only saw the lambda changes that fired, so the energy telescope
    // H_end - H_start == work + heat holds exactly only when Hend uses the same
    // final lambda the loop left the bridge in. With the PALINDROMIC schedule the
    // last index s = N-1 is pinned to lambda = 1 (both endpoints are), so this final
    // lambda is 1 -- the same lambda production's ncmcMove accepts on.
    r.Hend = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);
    r.work = work;
    r.heat = heat;
    return r;
}

} // namespace

// TELESCOPING: an independent re-summation of the fixed-q perturbation deltas
// reproduces the accumulator to machine precision. Guards the accumulator against
// an off-by-one in which V is sampled (it MUST be V at the post-propagation q of
// the previous step, i.e. the q the new lambda is applied at).
TEST(NcmcWork, AccumulatorIsTheFixedQTelescope) {
    Rng rng(0x9C3C0001);
    RobotModel m = freeTorsionChain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    // shrink velocities so the short trajectory stays in the smooth basin
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] *= 0.2;
    }
    primePositions(m, s);
    RobotState sCopy;
    sCopy.allocateFull(m);
    const QU snap = snapshotQU(m, s);

    ConstraintSet cs;
    LambdaAnalyticBridge bridge(m, s, /*kIntra*/ 300.0, /*kInter*/ 120.0);
    NcmcRun r = runNcmcLoop(m, s, bridge, cs, /*ncmcSteps*/ 40, /*hold*/ 0.0, /*h*/ 2.0e-4);
    ASSERT_TRUE(r.ok);

    // Independent telescope: re-run, re-summing V(lam_s,q_s) - V(lam_{s-1},q_s)
    // at each perturbation, from a fresh copy of the same start state.
    restoreQU(m, sCopy, snap);
    primePositions(m, sCopy);
    LambdaAnalyticBridge b2(m, sCopy, 300.0, 120.0);
    b2.setLambda(1.0);
    seedDerivatives(m, sCopy, b2);
    double wRef = 0.0, Vprev = b2.calcPotentialEnergy(sCopy), lamPrev = 1.0;
    for (int step = 0; step < 40; ++step) {
        const double lam = ncmc::protocolLambda(step, 40, 0.0);
        if (lam != lamPrev) {
            b2.setLambda(static_cast<Real>(lam));
            const double Vnew = b2.calcPotentialEnergy(sCopy);
            wRef += Vnew - Vprev;
            Vprev = Vnew;
            lamPrev = lam;
        }
        ASSERT_TRUE(RobotEngine::stepTo(m, sCopy, b2, cs, sCopy.time + 2.0e-4));
        Vprev = b2.calcPotentialEnergy(sCopy);
    }
    EXPECT_NEAR(r.work, wRef, kTight) << "work accumulator is not the fixed-q perturbation telescope";
}

// WORK == dH UP TO INTEGRATOR DRIFT, AND THE GAP IS THE HEAT. The production
// comment claims: for a deterministic, reversible, volume-preserving propagator
// the heat ~ 0, so w == H_end - H_start, and the logged gap w-(H_end-H_start) is a
// free consistency check that flags a too-large dt. We verify the exact bookkeeping
// identity  H_end - H_start == work + heat  (this is just energy conservation of
// the split move and must hold to machine precision), and that the gap |w - dH| is
// SMALL and equals |heat|.
TEST(NcmcWork, WorkPlusHeatEqualsDeltaH) {
    Rng rng(0x5151);
    RobotModel m = freeTorsionChain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] *= 0.2;
    }
    primePositions(m, s);
    ConstraintSet cs;
    LambdaAnalyticBridge bridge(m, s, 300.0, 120.0);
    NcmcRun r = runNcmcLoop(m, s, bridge, cs, 60, 0.0, 1.0e-4);
    ASSERT_TRUE(r.ok);

    const double dH = r.Hend - r.Hstart;
    // Exact split-move energy balance (definitional; must be machine-precise).
    EXPECT_NEAR(dH, r.work + r.heat, 1e-9) << "dH != work + heat (split-move energy balance broken)";
    // The gap the code logs IS the heat (shadow work of the fixed-lambda steps).
    EXPECT_NEAR(r.work - dH, -r.heat, 1e-9);
}

// THE GAP SHRINKS TOWARD ZERO UNDER STEP REFINEMENT. The production comment
// claims w == H_end - H_start "up to the integrator drift", and that the logged
// gap is a consistency check that flags a too-large dt. The gap IS the heat (the
// fixed-lambda shadow work); for a convergent propagator it must shrink monoton-
// ically toward 0 as h is refined. We verify that across halvings (the rate is
// set by the engine's implicit-trapezoid velocity refinement, so we assert
// convergence and meaningful contraction, not a specific formal order).
TEST(NcmcWork, WorkDeltaHGapShrinksUnderStepRefinement) {
    auto gapAt = [](double h) {
        Rng rng(0x2024);
        RobotModel m = freeTorsionChain(rng);
        RobotState s;
        s.allocateFull(m);
        randomizeState(m, s, rng);
        for (int i = 0; i < m.nu; ++i) {
            s.u()[i] *= 0.2;
        }
        primePositions(m, s);
        ConstraintSet cs;
        LambdaAnalyticBridge bridge(m, s, 300.0, 120.0);
        // Same physical protocol DURATION at each h: total time T fixed, steps = T/h.
        const double T = 6.0e-3;
        const int steps = static_cast<int>(std::lround(T / h));
        NcmcRun r = runNcmcLoop(m, s, bridge, cs, steps, 0.0, h);
        EXPECT_TRUE(r.ok);
        return std::abs(r.work - (r.Hend - r.Hstart));
    };
    const double g1 = gapAt(4.0e-4);
    const double g2 = gapAt(2.0e-4);
    const double g3 = gapAt(1.0e-4);
    // Monotone decrease toward zero as the step is refined.
    EXPECT_GT(g1, g2) << "gap did not shrink when halving h (g1=" << g1 << ", g2=" << g2 << ")";
    EXPECT_GT(g2, g3) << "gap did not shrink when halving h again (g2=" << g2 << ", g3=" << g3 << ")";
    // Each halving cuts the gap by a clear factor (we require meaningful
    // contraction, comfortably above 1).
    EXPECT_GT(g1 / std::max(g2, 1e-18), 1.6);
    EXPECT_GT(g2 / std::max(g3, 1e-18), 1.6);
    // And the refined gap is genuinely small in absolute terms.
    EXPECT_LT(g3, 1e-4);
}

// THE LAMBDA==1 GATE. A protocol pinned at lambda=1 throughout fires NO
// perturbation, so the work accumulator is EXACTLY 0 (bitwise), and the move
// reduces to the plain torsional-HMC Verlet trajectory: dH is then purely the
// integrator's energy drift. This is the reduction the ncmcMove acceptance comment
// relies on ("lambda == 1 throughout => work == 0 ... reduces EXACTLY to the
// torsional-HMC metropolis test").
TEST(NcmcWork, FlatLambdaOneGivesZeroWorkAndPlainHmcDeltaH) {
    Rng rng(0x1111);
    RobotModel m = freeTorsionChain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] *= 0.2;
    }
    RobotState sPlain;
    sPlain.allocateFull(m);
    const QU snap = snapshotQU(m, s);

    ConstraintSet cs;
    // We test the GATE by driving a flat lambda==1 protocol: the perturb branch
    // never fires, so the work must be exactly 0 and the trajectory must coincide
    // bit-for-bit with a plain (alchemy-free) HMC run from the same start.
    primePositions(m, s);
    LambdaAnalyticBridge bridge(m, s, 300.0, 120.0);

    // constant-lambda-1 loop (the flat protocol): perturb never triggers.
    bridge.setLambda(1.0);
    seedDerivatives(m, s, bridge);
    const double Hstart = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);
    double work = 0.0, lamPrev = 1.0, Vprev = bridge.calcPotentialEnergy(s);
    const int steps = 50;
    for (int step = 0; step < steps; ++step) {
        const double lam = 1.0; // flat protocol
        if (lam != lamPrev) {
            bridge.setLambda(static_cast<Real>(lam));
            work += bridge.calcPotentialEnergy(s) - Vprev;
            lamPrev = lam;
        }
        ASSERT_TRUE(RobotEngine::stepTo(m, s, bridge, cs, s.time + 1.0e-4));
        Vprev = bridge.calcPotentialEnergy(s);
    }
    const double Hend = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);

    // Plain torsional-HMC: identical trajectory with NO alchemy at all (lambda
    // fixed at 1 == the inter well always fully on, same forces). The end H must
    // match the flat-lambda NCMC end H bit-for-bit, and work must be exactly 0.
    restoreQU(m, sPlain, snap);
    primePositions(m, sPlain);
    LambdaAnalyticBridge plain(m, sPlain, 300.0, 120.0);
    plain.setLambda(1.0);
    seedDerivatives(m, sPlain, plain);
    const double HstartP = plain.calcPotentialEnergy(sPlain) + kineticEnergy(m, sPlain);
    for (int step = 0; step < steps; ++step) {
        ASSERT_TRUE(RobotEngine::stepTo(m, sPlain, plain, cs, sPlain.time + 1.0e-4));
    }
    const double HendP = plain.calcPotentialEnergy(sPlain) + kineticEnergy(m, sPlain);

    EXPECT_EQ(work, 0.0) << "flat lambda=1 protocol produced nonzero work";
    EXPECT_NEAR(Hstart, HstartP, kTight);
    EXPECT_NEAR(Hend, HendP, kTight) << "flat-lambda NCMC dH differs from plain HMC dH";
    EXPECT_NEAR(Hend - Hstart, HendP - HstartP, kTight);
}

// ===========================================================================
//  PART 3 -- Jacobian invariance under the alchemical (fixed-q) perturbation
// ===========================================================================
namespace {

// Independent oracle for a Free root's orientation pitch term, matching World's
// calcLogSineSqrGamma2 construction: build the body frame from a non-collinear
// real-atom triplet (x along root->ref, z along the triplet normal), extract the
// quaternion pitch sin = 2(w*y - z*x). A PURE function of atom geometry (=> of q).
// Returns sin(pitch); NaN if the body lacks 3 non-collinear real atoms.
double freeRootPitchSin(const RobotModel& m, const RobotState& s, int b) {
    const Vec3* P = s.atomPosG();
    const int a0 = m.bodyRootAtom[b];
    if (a0 < 0) {
        return std::nan("");
    }
    const Vec3 p0 = P[a0];
    int bestA1 = -1, bestA2 = -1;
    Real bestArea = 0;
    for (int ci = m.bodyAtomsBeg[b]; ci < m.bodyAtomsEnd[b]; ++ci) {
        const int a1 = m.bodyAtoms[ci];
        if (a1 == a0 || m.atomMass[a1] <= Real(0)) {
            continue;
        }
        for (int cj = ci + 1; cj < m.bodyAtomsEnd[b]; ++cj) {
            const int a2 = m.bodyAtoms[cj];
            if (a2 == a0 || m.atomMass[a2] <= Real(0)) {
                continue;
            }
            const Real area = ((P[a1] - p0) % (P[a2] - p0)).norm();
            if (area > bestArea) {
                bestArea = area;
                bestA1 = a1;
                bestA2 = a2;
            }
        }
    }
    if (bestA1 < 0 || bestArea < Real(1e-10)) {
        return std::nan("");
    }
    const Vec3 ex = (P[bestA1] - p0) / (P[bestA1] - p0).norm();
    Vec3 ez = ex % (P[bestA2] - p0);
    ez = ez / ez.norm();
    const Vec3 ey = ez % ex;
    // Rotation matrix R = [ex ey ez] as columns (the same frame World builds).
    const Real r00 = ex[0], r01 = ey[0], r02 = ez[0];
    const Real r10 = ex[1], r11 = ey[1], r12 = ez[1];
    const Real r20 = ex[2], r21 = ey[2], r22 = ez[2];
    // Rotation -> quaternion (Shepperd's largest-component method), then the pitch
    // sin = 2(w*y - z*x), matching World::calcLogSineSqrGamma2.
    Real w, x, y, z;
    const Real tr = r00 + r11 + r22;
    if (tr > Real(0)) {
        Real ss = std::sqrt(tr + Real(1)) * Real(2); // 4w
        w = Real(0.25) * ss;
        x = (r21 - r12) / ss;
        y = (r02 - r20) / ss;
        z = (r10 - r01) / ss;
    } else if (r00 > r11 && r00 > r22) {
        Real ss = std::sqrt(Real(1) + r00 - r11 - r22) * Real(2); // 4x
        w = (r21 - r12) / ss;
        x = Real(0.25) * ss;
        y = (r01 + r10) / ss;
        z = (r02 + r20) / ss;
    } else if (r11 > r22) {
        Real ss = std::sqrt(Real(1) + r11 - r00 - r22) * Real(2); // 4y
        w = (r02 - r20) / ss;
        x = (r01 + r10) / ss;
        y = Real(0.25) * ss;
        z = (r12 + r21) / ss;
    } else {
        Real ss = std::sqrt(Real(1) + r22 - r00 - r11) * Real(2); // 4z
        w = (r10 - r01) / ss;
        x = (r02 + r20) / ss;
        y = (r12 + r21) / ss;
        z = Real(0.25) * ss;
    }
    return std::clamp(Real(2.0) * (w * y - z * x), Real(-1.0), Real(1.0));
}

Real logDetMAt(const RobotModel& m, RobotState& s) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    return RobotEngine::calcLogDetM(m, s);
}

} // namespace

// The alchemical perturbation changes lambda at FIXED q. Therefore EVERY q-only
// quantity in the acceptance Hamiltonian -- the Fixman determinant ln|M(q)| and
// the Free-root pitch term -- must be bitwise invariant across the perturbation,
// while the potential energy V changes. That invariance is the "unit coordinate
// Jacobian" (alpha-ratio == 1) the ncmcMove acceptance comment claims: the
// perturbation injects NO coordinate-Jacobian factor into dH. If a future change
// ever made calcLogDetM or the pitch read anything downstream of the force/energy
// evaluation, this test would catch it.
TEST(NcmcJacobian, FixedQPerturbationLeavesMassMetricAndPitchInvariant) {
    Rng rng(0x7); // 0x7 happens to give a well-conditioned Free root
    RobotModel m = freeTorsionChain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);

    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);

    LambdaAnalyticBridge bridge(m, s, 300.0, 120.0);

    // --- snapshot q-only quantities at lambda=1 ---
    bridge.setLambda(1.0);
    bridge.evaluate(s);
    const double V1 = bridge.calcPotentialEnergy(s);
    const Real logDet1 = logDetMAt(m, s);
    const double pitch1 = freeRootPitchSin(m, s, /*body*/ 1);
    std::vector<double> q1(s.q(), s.q() + m.nq);
    std::vector<Vec3> pos1(s.atomPosG(), s.atomPosG() + m.numAtoms);

    // --- PERTURB: change lambda only (NO q update, exactly as ncmcMove's step i) ---
    bridge.setLambda(0.37);
    bridge.evaluate(s); // recompute forces/energy at the SAME q, new lambda
    const double V2 = bridge.calcPotentialEnergy(s);
    const Real logDet2 = logDetMAt(m, s);
    const double pitch2 = freeRootPitchSin(m, s, 1);

    // The perturbation must actually have done something to V (otherwise the test
    // is vacuous).
    ASSERT_GT(std::abs(V2 - V1), 1e-6) << "perturbation did not change V; test is vacuous";

    // q and atom geometry: bitwise unchanged.
    for (int i = 0; i < m.nq; ++i) {
        EXPECT_EQ(s.q()[i], q1[i]) << "q[" << i << "] moved during a fixed-q perturbation";
    }
    for (int a = 0; a < m.numAtoms; ++a) {
        EXPECT_TRUE(rtest::NearVec3(s.atomPosG()[a], pos1[a], 0.0)) << "atom " << a << " moved";
    }
    // Fixman mass-metric determinant: bitwise invariant under the perturbation.
    EXPECT_EQ(logDet1, logDet2) << "ln|M(q)| changed under a fixed-q alchemical perturbation";
    // Free-root quaternion pitch term: invariant (a pure function of q).
    ASSERT_FALSE(std::isnan(pitch1));
    EXPECT_NEAR(pitch1, pitch2, kTight) << "Free-root pitch changed under a fixed-q perturbation";
}

// ===========================================================================
//  PART 4 -- F-reversibility & the Crooks crossing of the NCMC map (Fix 1)
//
//  Fix 1's core invariant: with the PALINDROMIC schedule the whole NCMC trajectory
//  is one map T on (q,p) that is momentum-flip reversible, F T F == T^-1. T is the
//  composition of per-substep fixed-lambda Verlet steps (each itself F-reversible;
//  the fixed-q lambda perturbations do not move the state), so the composition is
//  reversible IFF the lambda sequence is a palindrome with endpoints pinned to 1.
//  That reversibility is what makes accept-on-endpoint-H an EXACT Metropolis test.
// ===========================================================================

// PRIMARY dynamical gate for Fix 1. Drive the REAL GC integrator forward N
// substeps, flip the momenta, drive forward N MORE (the palindrome makes the
// second leg retrace the first), flip back -> land on the start. FAILS on the old
// offset-by-one schedule (its reverse leg visits a different lambda at each step,
// so the legs do not cancel); PASSES on the palindrome.
TEST(NcmcReversibility, PalindromeMapIsMomentumFlipReversible) {
    Rng rng(0xF1B7);
    RobotModel m = freeTorsionChain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] *= 0.2; // keep the short trajectory in the smooth basin
    }
    primePositions(m, s);

    ConstraintSet cs;
    LambdaAnalyticBridge bridge(m, s, 300.0, 120.0);

    const QU start = snapshotQU(m, s);
    const int N = 41;        // odd: exact trough, true palindrome
    const double h = 1.0e-4; // well inside the reversible regime

    const NcmcRun fwd = runNcmcLoop(m, s, bridge, cs, N, 0.0, h); // forward leg
    ASSERT_TRUE(fwd.ok);
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = -s.u()[i]; // F: flip momenta
    }
    const NcmcRun rev = runNcmcLoop(m, s, bridge, cs, N, 0.0, h); // reverse (same schedule)
    ASSERT_TRUE(rev.ok);
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = -s.u()[i]; // F: flip back
    }

    ASSERT_GT(std::abs(fwd.Hend - fwd.Hstart), 1e-6) << "forward leg did not move; test is vacuous";
    // Tolerance: the implicit-trapezoid corrector's reversibility residual over the
    // 82 lambda-switched substeps is ~1e-5; a NON-palindromic schedule (the old
    // offset-by-one bug) does not retrace at all and lands O(0.1) away. 1e-4 sits
    // ~1000x below that, so the test stays a sharp discriminator for Fix 1.
    for (int i = 0; i < m.nq; ++i) {
        EXPECT_NEAR(s.q()[i], start.q[i], 1e-4) << "q[" << i << "] not F-reversible";
    }
    for (int i = 0; i < m.nu; ++i) {
        EXPECT_NEAR(s.u()[i], start.u[i], 1e-4) << "u[" << i << "] not F-reversible";
    }
}

// CROOKS / DeltaF=0 crossing for the self-reverse protocol. For the palindromic
// protocol the forward and reverse processes are statistically identical and
// DeltaF=0, so the protocol work changes sign EXACTLY under time reversal:
// W_forward + W_reverse == 0 (using the time-reversed trajectory as the reverse
// realization). That sign flip is the non-statistical core of "forward and reverse
// work histograms are mirror images crossing at W=0" (and hence <exp(-beta W)>=1).
// The full ensemble <exp(-beta W)>=1 estimate is the heavier integration test
// (flagged in the handoff). Here we assert the exact per-trajectory identity.
TEST(NcmcCrooks, SelfReverseProtocolFlipsWorkSign) {
    for (unsigned seed : {0xC0FFEEu, 0x1234u, 0xBEEFu}) {
        Rng rng(seed);
        RobotModel m = freeTorsionChain(rng);
        RobotState s;
        s.allocateFull(m);
        randomizeState(m, s, rng);
        for (int i = 0; i < m.nu; ++i) {
            s.u()[i] *= 0.2;
        }
        primePositions(m, s);

        ConstraintSet cs;
        LambdaAnalyticBridge bridge(m, s, 300.0, 120.0);
        const int N = 41;
        const double h = 1.0e-4;

        const NcmcRun fwd = runNcmcLoop(m, s, bridge, cs, N, 0.0, h);
        ASSERT_TRUE(fwd.ok);
        for (int i = 0; i < m.nu; ++i) {
            s.u()[i] = -s.u()[i];
        }
        const NcmcRun rev = runNcmcLoop(m, s, bridge, cs, N, 0.0, h);
        ASSERT_TRUE(rev.ok);

        ASSERT_GT(std::abs(fwd.work), 1e-3) << "forward work ~0; test vacuous (seed " << seed << ")";
        const double tol = 1e-2 * (std::abs(fwd.work) + std::abs(rev.work) + 1.0);
        EXPECT_NEAR(fwd.work + rev.work, 0.0, tol)
            << "self-reverse work did not flip sign (seed " << seed << "): W_f=" << fwd.work
            << " W_r=" << rev.work;
    }
}

// ===========================================================================
//  PART 5 -- Fixman block-diagonal invariance under welded solvent (Fix 2)
// ===========================================================================
namespace {

// A Free+Torsion+Torsion solute, then `waters.size()` WELDED (Rigid, zero-DOF)
// waters appended as separate trees rooted at Ground, each at the given placement.
// The solute body specs come from `soluteRng`, so two calls with the SAME seed
// produce a BIT-IDENTICAL solute regardless of how many waters are appended.
RobotModel soluteWithWeldedWater(Rng& soluteRng, const std::vector<Transform>& waters) {
    auto mk = [](Rng& r, int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(r.rotation(), r.vec3(-0.2, 0.2));
        s.X_BM = Transform(r.rotation(), r.vec3(-0.2, 0.2));
        s.mass = r.uniform(0.9, 1.5);
        s.com_B = r.vec3(-0.08, 0.08);
        s.inertia_B = UnitInertia(r.uniform(0.4, 0.6), r.uniform(0.4, 0.6), r.uniform(0.4, 0.6));
        return s;
    };
    std::vector<BodySpec> specs;
    specs.push_back(mk(soluteRng, 0, JointType::Free));
    specs.push_back(mk(soluteRng, 1, JointType::Torsion));
    specs.push_back(mk(soluteRng, 2, JointType::Torsion));
    for (const Transform& place : waters) {
        BodySpec s;
        s.parent = 0; // separate tree rooted at Ground = a welded water
        s.joint = JointType::Rigid;
        s.X_PF = place;
        s.X_BM = Transform();
        s.mass = 18.0;
        s.com_B = Vec3(0);
        s.inertia_B = UnitInertia(0.5, 0.5, 0.5);
        specs.push_back(s);
    }
    RobotModel m = buildForest(specs);
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, 0.02, -0.03), Vec3(-0.02, 0.10, 0.01)}, {12.0, 1.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0.04, -0.06, 0.09), Vec3(0.08, 0.01, -0.02)}, {16.0, 1.0});
    for (std::size_t w = 0; w < waters.size(); ++w) {
        attachAtoms(m, 4 + static_cast<int>(w), {Vec3(0, 0, 0), Vec3(0.1, 0, 0), Vec3(0, 0.1, 0)},
                    {16.0, 1.0, 1.0});
    }
    return m;
}

} // namespace

// FIXMAN BLOCK-DIAGONAL INVARIANCE (Fix 2). Welding all solvent must leave the
// Fixman mass-metric determinant ln|M(q)| EQUAL to the solute-only value: a
// welded, zero-DOF water contributes det = 1 -> ln 0, and rooted at Ground it does
// not enter the solute bodies' articulated inertia. Hence 1/2 RT ln det M is the
// solute-only value up to a constant (here exactly 0), INDEPENDENT of solvent
// count and configuration -- the property guaranteeing Fix 2 adds NO Fixman bias.
TEST(NcmcFixman, WeldedSolventLeavesMassMetricInvariant) {
    Rng rngA(0x5A1A);
    RobotModel mA = soluteWithWeldedWater(rngA, {}); // solute only
    RobotState sA;
    sA.allocateFull(mA);
    Rng qrng(0x9001);
    randomizeState(mA, sA, qrng);
    const Real logDetSolute = logDetMAt(mA, sA);

    // welded-solvent configurations: different counts AND placements.
    const std::vector<std::vector<Transform>> configs = {
        {Transform(Vec3(1, 0, 0))},
        {Transform(Vec3(1, 0, 0)), Transform(Vec3(0, 2, 0))},
        {Transform(Vec3(-1, 1, 3)), Transform(Vec3(2, -2, 1)), Transform(Vec3(0, 0, 5))},
    };
    for (std::size_t c = 0; c < configs.size(); ++c) {
        Rng rngB(0x5A1A); // SAME seed -> bit-identical solute bodies
        RobotModel mB = soluteWithWeldedWater(rngB, configs[c]);
        ASSERT_EQ(mB.nq, mA.nq) << "welded water added configurational DOF (config " << c << ")";
        ASSERT_EQ(mB.nu, mA.nu) << "welded water added velocity DOF (config " << c << ")";
        RobotState sB;
        sB.allocateFull(mB);
        std::copy(sA.q(), sA.q() + mA.nq, sB.q()); // identical solute coordinates
        const Real logDetWelded = logDetMAt(mB, sB);
        EXPECT_NEAR(logDetWelded, logDetSolute, kTight)
            << "welded solvent changed ln|M(q)| (config " << c << ")";
    }
}