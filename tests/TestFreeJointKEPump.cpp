// ============================================================================
//  TestFreeJointKEPump.cpp -- DIAGNOSIS reproducer for the documented Free/Ball
//  root-joint kinetic-energy pump (RobotIntegrator.hpp:216-227,
//  docs/specs/two-robot-contact/30-friction-diagnosis.md,
//  docs/specs/two-robot-contact/50-validation.md LEMMA DR1/DR2).
//
//  QUESTION. Is the "friction/dampening" seen by the two-robot contact-world
//  campaign the documented numerical KE pump (gamma-pump), or benign O(dt^2)
//  physical discretization (gamma-phys)?
//
//  This file does NOT touch RobotIntegrator.hpp or the exp-map propagator. It
//  runs the REAL, SHIPPED RobotEngine::verletStep (which already uses the
//  exp-map quaternion advance) through an NVE protocol (no momentum refresh,
//  pure Hamiltonian flow) with an OpenMM-free AnalyticForceBridge, mirroring
//  the fixtures/harness already established in tests/TestIntegrator.cpp
//  (singleFree/chain3/freeBall, buildSeeded, runTrajectory, lsSlope) --
//  extended here with (a) an explicit dt-scaling sweep and (b) per-body KE
//  localization, per DR1/DR2.
//
//  A SECOND, DELIBERATE ablation is included: `ablatedVerletStepLinearTaylorQuat`
//  is a test-only, byte-for-byte-except-one-block COPY of the production
//  verletStep (RobotIntegrator.hpp:102-433) with ONLY the exp-map quaternion
//  overwrite loop (RobotIntegrator.hpp:232-261) REMOVED. Since the raw
//  linear-Taylor q update (q1 = q0 + h*qdot0 + h^2/2*qddot0, computed
//  unconditionally just above that loop) is left standing, and
//  normalizeQuaternions() still runs unconditionally afterward, this
//  reproduces exactly the historical "linear-Taylor-plus-renormalize" quaternion
//  update the RobotIntegrator.hpp comment blames for the "~+1300 kJ/mol/traj"
//  pump -- giving a clean, in-repo POSITIVE CONTROL for the same NVE protocol,
//  run through the IDENTICAL fixtures/seeds. This ablation is EXPLORATORY /
//  documentation-only (it exercises a path nothing in production takes) and is
//  never used to gate the build; only the SHIPPED-propagator tests do that.
//  This mirrors the existing test-side-reproduction convention in
//  tests/TestNcmcExplicitSolvent.cpp ("intentionally REPRODUCES World::
//  ncmcMove's Construction-I and Construction-II").
//
//  UNITS CAVEAT. AnalyticForceBridge's harmonic constant k and the hand-built
//  masses/inertias are NOT calibrated against a real force field (no OpenMM
//  dependency, by design -- "OpenMM-free if possible so it runs anywhere").
//  `h` is in the engine's internal time unit (the corrector's own diagnostic
//  message labels it "ps"); the dt sweep below spans a factor-of-8 range
//  (kBaseH/2 .. 4*kBaseH) matching the REQUESTED 0.25/0.5/1/2 fs ratios, not
//  literal wall-clock femtoseconds. The discriminating structure (O(dt^2)
//  scaling, no secular term, no per-body localization) is unit-independent.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

#include "AnalyticForceBridge.hpp"
#include "Constraints.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotIntegrator.hpp" // the templated, SHIPPED verletStep/stepTo (exp-map)
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::attachAtoms;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::randomizeState;
using rtest::Rng;

namespace {

// ---------------------------------------------------------------------------
//  Fixtures (mirror tests/TestIntegrator.cpp singleFree/chain3/freeBall).
// ---------------------------------------------------------------------------

// A single Free root -- the cleanest isolation named in the spec.
RobotModel singleFree(Rng& rng) {
    BodySpec s;
    s.parent = 0;
    s.joint = JointType::Free;
    s.X_PF = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
    s.X_BM = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
    s.mass = rng.uniform(0.8, 1.6);
    s.com_B = rng.vec3(-0.1, 0.1);
    s.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
    RobotModel m = buildForest({s});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.10, 0.02, -0.03), Vec3(-0.04, 0.07, 0.0)}, {12.0, 1.0, 1.0});
    return m;
}

// Free root -> Ball child: BOTH bodies carry a quaternion DOF, so a root-joint
// quaternion-kinematics defect that "propagates through the articulated
// recursion" (RobotIntegrator.hpp:223) has a second body to show up in.
RobotModel freeBallChain(Rng& rng) {
    auto mk = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        s.X_BM = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        s.mass = rng.uniform(0.8, 1.6);
        s.com_B = rng.vec3(-0.1, 0.1);
        s.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
        return s;
    };
    RobotModel m = buildForest({mk(0, JointType::Free), mk(1, JointType::Ball)});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.10, 0.02, -0.03)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0), Vec3(-0.03, 0, 0.06)}, {14.0, 1.0});
    return m;
}

void seedDerivatives(const RobotModel& m, RobotState& s, AnalyticForceBridge& bridge) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    bridge.evaluate(s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    RobotEngine::calcUDot(m, s);
    RobotEngine::calcQDot(m, s, s.qdot());
    RobotEngine::calcQDotDot(m, s);
}

constexpr Real kBoostRT = Real(2.5);

// Production momentum draw (u = sqrt(RT) * sqrt(M^-1) * gaussian), displaced
// off the anchor pose (so U != 0). Identical convention to TestIntegrator.cpp.
void buildSeeded(const RobotModel& m,
                 RobotState& s,
                 Rng& rng,
                 Real k,
                 std::vector<AnalyticForceBridge>& bridgeStore) {
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    bridgeStore.emplace_back(m, s, k); // anchors = this pose

    randomizeState(m, s, rng); // displace off the anchors
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);

    std::vector<Real> w(static_cast<std::size_t>(m.nu)), u(static_cast<std::size_t>(m.nu));
    for (int i = 0; i < m.nu; ++i) {
        w[static_cast<std::size_t>(i)] = rng.gaussian();
    }
    RobotEngine::multiplyBySqrtMInv(m, s, w.data(), u.data());
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = std::sqrt(kBoostRT) * u[static_cast<std::size_t>(i)];
    }

    seedDerivatives(m, s, bridgeStore.back());
}

Real lsSlope(const std::vector<Real>& y) {
    const Real n = static_cast<Real>(y.size());
    Real sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < y.size(); ++i) {
        const Real x = static_cast<Real>(i);
        sx += x;
        sy += y[i];
        sxx += x * x;
        sxy += x * y[i];
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

Real bodyKE(const RobotModel& m, const RobotState& s, int b) {
    (void)m;
    const SpatialInertia* Mk = s.Mk_G();
    const SpatialVec* V = s.V_GB();
    return Real(0.5) * dot(V[b], Mk[b] * V[b]);
}

// ---------------------------------------------------------------------------
//  NVE trajectory driver over the SHIPPED verletStep. Records total energy
//  and per-body KE at every step. `ok` is false iff a step was rejected
//  (non-finite force / corrector runaway) -- DR1/DR2 are not meaningful past
//  a rejected step, so callers ASSERT_TRUE(d.ok).
// ---------------------------------------------------------------------------
struct NveTrace {
    bool ok = true;
    std::vector<Real> totalE;         // KE+U per step
    std::vector<std::vector<Real>> keByBody; // [body][step], body in [1, numBodies)
};

NveTrace runNve(const RobotModel& m,
                RobotState& s,
                AnalyticForceBridge& bridge,
                const ConstraintSet& cset,
                Real h,
                int nSteps) {
    NveTrace tr;
    tr.keByBody.resize(static_cast<std::size_t>(m.numBodies));
    for (int n = 0; n < nSteps; ++n) {
        if (!RobotEngine::verletStep(m, s, bridge, cset, h)) {
            tr.ok = false;
            return tr;
        }
        tr.totalE.push_back(RobotEngine::calcKineticEnergy(m, s) + bridge.calcPotentialEnergy(s));
        for (int b = 1; b < m.numBodies; ++b) {
            tr.keByBody[static_cast<std::size_t>(b)].push_back(bodyKE(m, s, b));
        }
    }
    return tr;
}

// DR1 summary for one trace: bounded max drift, and a SCALE-FREE secularity
// ratio (|slope|*N / maxAbs) -- this isolates "is there a trend" from "how
// big is the O(dt^2) oscillation", so it is meaningful across the whole dt
// sweep without a dt-dependent magic number.
struct Dr1Summary {
    Real maxAbs = 0;
    Real secularRatio = 0; // |slope|*N / maxAbs; ~0 for pure oscillation, ~1 for a pure ramp
};
Dr1Summary dr1(const std::vector<Real>& e, Real e0) {
    Dr1Summary r;
    for (Real v : e) {
        r.maxAbs = std::max(r.maxAbs, std::abs(v - e0));
    }
    const Real slope = std::abs(lsSlope(e)) * static_cast<Real>(e.size());
    r.secularRatio = slope / (r.maxAbs + Real(1e-300));
    return r;
}

// ===========================================================================
//  Ablated positive control: byte-for-byte copy of RobotIntegrator.hpp's
//  verletStep (lines 102-433 at the time of writing) with ONLY the exp-map
//  quaternion-overwrite loop removed, so the position update for Free/Ball/
//  FreeLine quaternion blocks falls back to the raw linear-Taylor result
//  (q1 = q0 + h*qdot0 + h^2/2*qddot0) computed by the SAME loop that already
//  handles every scalar DOF, renormalized afterward by the SAME
//  normalizeQuaternions() call the production step already makes. This is
//  exactly the "linear-Taylor-plus-renormalize" scheme the production comment
//  (RobotIntegrator.hpp:214-227) names as the historical KE-pump source.
//  NEVER call this outside this diagnostic file; it exists solely to confirm
//  or refute that the documented defect is real when the exp-map workaround
//  is absent, on the exact same fixtures/seeds as the shipped-propagator
//  tests above. If RobotIntegrator.hpp's verletStep changes, re-sync this by
//  diffing against the current source.
// ===========================================================================
template <class Bridge>
bool ablatedVerletStepLinearTaylorQuat(const RobotModel& m,
                                       RobotState& s,
                                       Bridge& bridge,
                                       const robo::ConstraintSet& cset,
                                       Real h,
                                       bool* correctorConverged = nullptr) {
    if (correctorConverged) {
        *correctorConverged = false;
    }
    const int nq = m.nq, nu = m.nu;
    Real* q = s.q();
    Real* u = s.u();
    Real* qdot = s.qdot();
    Real* udot = s.udot();
    Real* qdd = s.qdotdot();

    std::vector<Real> q0(q, q + nq);
    std::vector<Real> u0(u, u + nu);
    std::vector<Real> qdot0(qdot, qdot + nq);
    std::vector<Real> udot0(udot, udot + nu);
    std::vector<Real> qdd0(qdd, qdd + nq);

    // no Cartesian solvent in this diagnostic's fixtures -- solvAtoms empty.
    const std::vector<int>& solvAtoms = s.cartSolventAtoms();
    const std::vector<Real>& solvInvM = s.cartSolventInvMass();
    const int nSolv = static_cast<int>(solvAtoms.size());
    Vec3* posG = s.atomPosG();
    Vec3* velG = s.atomVelG();
    Vec3* frcG = s.atomForceG();
    std::vector<Vec3> xs0(nSolv), vs0(nSolv), fs0(nSolv);
    for (int j = 0; j < nSolv; ++j) {
        const int a = solvAtoms[j];
        xs0[j] = posG[a];
        vs0[j] = velG[a];
        fs0[j] = frcG[a];
    }

    auto forcesFinite = [&]() -> bool {
        const SpatialVec* bf = s.bodyForceG();
        for (int b = 1; b < m.numBodies; ++b) {
            if (!std::isfinite(bf[b][0][0]) || !std::isfinite(bf[b][0][1]) || !std::isfinite(bf[b][0][2])
                || !std::isfinite(bf[b][1][0]) || !std::isfinite(bf[b][1][1])
                || !std::isfinite(bf[b][1][2])) {
                return false;
            }
        }
        return true;
    };

    Real uSeedNorm2 = 0;
    for (int i = 0; i < nu; ++i) {
        uSeedNorm2 += u0[i] * u0[i];
    }
    const Real uCap2 = (uSeedNorm2 + Real(1e-30)) * Real(1e6);
    auto velocitiesSane = [&]() -> bool {
        Real n2 = 0;
        for (int i = 0; i < nu; ++i) {
            if (!std::isfinite(u[i])) {
                return false;
            }
            n2 += u[i] * u[i];
        }
        return n2 <= uCap2;
    };

    auto restorePreStep = [&]() {
        std::copy(q0.begin(), q0.end(), q);
        std::copy(u0.begin(), u0.end(), u);
        std::copy(qdot0.begin(), qdot0.end(), qdot);
        std::copy(udot0.begin(), udot0.end(), udot);
        std::copy(qdd0.begin(), qdd0.end(), qdd);
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        for (int j = 0; j < nSolv; ++j) {
            const int a = solvAtoms[j];
            posG[a] = xs0[j];
            velG[a] = vs0[j];
            frcG[a] = fs0[j];
        }
    };

    // ---- position drift: linear-Taylor for EVERY q slot, quaternion blocks
    //      INCLUDED (this is the ablation -- production overwrites the
    //      quaternion blocks with the exp-map here; that overwrite loop is
    //      deliberately NOT reproduced below). ----
    for (int i = 0; i < nq; ++i) {
        q[i] = q0[i] + (h * qdot0[i]) + ((h * h / 2) * qdd0[i]);
    }
    RobotEngine::normalizeQuaternions(m, s); // renormalize the raw Taylor quaternion (historical scheme)

    {
        const Real h2half = Real(0.5) * h * h;
        for (int j = 0; j < nSolv; ++j) {
            const int a = solvAtoms[j];
            const Real im = solvInvM[j];
            posG[a] = Vec3(xs0[j][0] + h * vs0[j][0] + h2half * im * fs0[j][0],
                           xs0[j][1] + h * vs0[j][1] + h2half * im * fs0[j][1],
                           xs0[j][2] + h * vs0[j][2] + h2half * im * fs0[j][2]);
        }
    }

    auto refreshPos = [&]() {
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
    };
    refreshPos();
    cset.enforcePositionConstraints(m, s, refreshPos);

    for (int i = 0; i < nu; ++i) {
        u[i] = u0[i] + (h * udot0[i]);
    }

    auto evalDerivs = [&]() -> bool {
        RobotEngine::realizePosition(m, s);
        bridge.evaluate(s);
        if (!forcesFinite()) {
            return false;
        }
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        RobotEngine::calcUDot(m, s);
        for (int i = 0; i < nu; ++i) {
            if (!std::isfinite(udot[i])) {
                return false;
            }
        }
        RobotEngine::calcQDot(m, s, qdot);
        RobotEngine::calcQDotDot(m, s);
        return true;
    };
    if (!evalDerivs()) {
        restorePreStep();
        return false;
    }

    const Real tol = Real(1e-4);
    Real prevChange = std::numeric_limits<Real>::infinity();
    Real lastChange = std::numeric_limits<Real>::infinity();
    int usedIters = 0;
    bool converged = false;

    for (int iter = 0; iter < 10; ++iter) {
        ++usedIters;
        Real num = 0;
        Real den = 0;

        for (int i = 0; i < nu; ++i) {
            const Real un = u0[i] + ((h / 2) * (udot0[i] + udot[i]));
            const Real d = un - u[i];
            num += d * d;
            den += u[i] * u[i];
            u[i] = un;
        }

        if (!velocitiesSane()) {
            restorePreStep();
            return false;
        }
        if (!evalDerivs()) {
            restorePreStep();
            return false;
        }

        const Real change = std::sqrt(num) / (std::sqrt(den) + Real(1e-30));
        lastChange = change;
        if (!std::isfinite(change) || change > Real(1e6)) {
            restorePreStep();
            return false;
        }

        if (change <= tol) {
            converged = true;
            break;
        }

        if (iter > 1 && change > prevChange) {
            break;
        }

        prevChange = change;
    }

    if (correctorConverged) {
        *correctorConverged = converged;
    }
    (void)usedIters;
    (void)lastChange;

    cset.enforceVelocityConstraints(m, s);
    RobotEngine::realizeVelocity(m, s);

    for (int j = 0; j < nSolv; ++j) {
        const int a = solvAtoms[j];
        const Real hHalfInvM = Real(0.5) * h * solvInvM[j];
        velG[a] = Vec3(vs0[j][0] + hHalfInvM * (fs0[j][0] + frcG[a][0]),
                       vs0[j][1] + hHalfInvM * (fs0[j][1] + frcG[a][1]),
                       vs0[j][2] + hHalfInvM * (fs0[j][2] + frcG[a][2]));
    }

    s.time += h;
    return true;
}

NveTrace runNveAblated(const RobotModel& m,
                       RobotState& s,
                       AnalyticForceBridge& bridge,
                       const ConstraintSet& cset,
                       Real h,
                       int nSteps) {
    NveTrace tr;
    tr.keByBody.resize(static_cast<std::size_t>(m.numBodies));
    for (int n = 0; n < nSteps; ++n) {
        if (!ablatedVerletStepLinearTaylorQuat(m, s, bridge, cset, h)) {
            tr.ok = false;
            return tr;
        }
        tr.totalE.push_back(RobotEngine::calcKineticEnergy(m, s) + bridge.calcPotentialEnergy(s));
        for (int b = 1; b < m.numBodies; ++b) {
            tr.keByBody[static_cast<std::size_t>(b)].push_back(bodyKE(m, s, b));
        }
    }
    return tr;
}

// dt sweep base unit: matches tests/TestIntegrator.cpp's kSafeH (a step well
// inside the reversible regime for this harness).
constexpr Real kBaseH = Real(5e-4);
// matched total physical time across the sweep. Long enough to cover many
// oscillation periods at every dt in the sweep (matches the window used by
// the FreeBallChain DR2 test below, 4000 steps at kBaseH) -- a short window
// (e.g. tests/TestIntegrator.cpp's SecondOrderEnergyError 600-step window,
// fine for an O(dt^2)-amplitude-only check) makes a least-squares secular-
// trend fit pick up partial-period curvature as a false "trend": that is a
// window artifact, not a signature of a real pump (a genuine dt-independent
// pump would make the trend/oscillation RATIO grow as dt shrinks, not stay
// flat/shrink the way a windowing artifact does -- see the dt-monotonicity
// assertion below).
constexpr Real kMatchedTime = Real(4000) * kBaseH;

} // namespace

// ---------------------------------------------------------------------------
//  DR1 (LEMMA, docs/specs/two-robot-contact/50-validation.md 5): total-energy
//  drift over the SHIPPED (exp-map) propagator is bounded, non-secular, and
//  O(dt^2) -- across a dt sweep spanning the requested 0.25x/0.5x/1x/2x/4x
//  ratios, at MATCHED physical time, on a single Free body (the cleanest
//  isolation named in the diagnosis task).
// ---------------------------------------------------------------------------
TEST(FreeJointKEPump, DR1_ShippedPropagator_SingleFree_NoSecularTerm_DtSweep) {
    const ConstraintSet cset;
    const Real dts[] = {kBaseH / 2, kBaseH, 2 * kBaseH, 4 * kBaseH};

    Real prevMaxAbs = -1;
    std::vector<Real> secularRatios;
    for (Real h : dts) {
        Rng rng(0x6001);
        RobotModel m = singleFree(rng);
        RobotState s;
        s.allocateFull(m);
        std::vector<AnalyticForceBridge> store;
        store.reserve(1);
        Rng r(0x6001 ^ 0xABCDu); // identical start at every dt
        buildSeeded(m, s, r, /*k=*/60.0, store);

        const int nSteps = static_cast<int>(std::lround(static_cast<double>(kMatchedTime / h)));
        const Real e0 = RobotEngine::calcKineticEnergy(m, s) + store.back().calcPotentialEnergy(s);
        const NveTrace tr = runNve(m, s, store.back(), cset, h, nSteps);
        ASSERT_TRUE(tr.ok) << "h=" << h << ": a step was rejected -- DR1 not measurable at this dt "
                              "(that is itself a gamma-corr signal, not gamma-pump)";

        const Dr1Summary d = dr1(tr.totalE, e0);
        // no secular (monotone) term: the trend must be a small fraction of the
        // trajectory's own oscillation amplitude, at EVERY dt in the sweep.
        EXPECT_LT(d.secularRatio, Real(0.15))
            << "h=" << h << ": secular-vs-oscillation ratio " << d.secularRatio
            << " -- a secular term surviving across the dt sweep is gamma-pump, not gamma-phys";
        std::fprintf(stderr, "[DR1 singleFree] h=%.6g nSteps=%d maxAbsDrift=%.6g secularRatio=%.4g\n",
                     (double)h, nSteps, (double)d.maxAbs, (double)d.secularRatio);
        secularRatios.push_back(d.secularRatio);

        if (prevMaxAbs >= 0) {
            const Real ratio = d.maxAbs / (prevMaxAbs + Real(1e-300));
            // O(dt^2): each doubling of h should ~4x the drift amplitude.
            EXPECT_GE(ratio, Real(2.5))
                << "h=" << h << ": drift did not grow like O(dt^2) (ratio=" << ratio << ")";
            EXPECT_LE(ratio, Real(6.0))
                << "h=" << h << ": drift grew FASTER than O(dt^2) (ratio=" << ratio
                << ") -- consistent with a dt-independent pump term, not pure discretization";
        }
        prevMaxAbs = d.maxAbs;
    }

    // The crux of LEMMA DR1 ("drift RATE -> 0 as dt -> 0"): a genuine
    // dt-independent secular pump has a roughly CONSTANT total injection at
    // matched physical time (accumulated over more, smaller steps), while the
    // O(dt^2) truncation-error oscillation amplitude SHRINKS -- so their ratio
    // would grow toward 1 as dt->0. Pure O(dt^2) discretization (gamma-phys)
    // instead keeps the trend/oscillation ratio flat or shrinking as dt->0,
    // since both terms come from the same truncation error. Assert the
    // smallest-dt ratio has NOT grown relative to the largest-dt ratio.
    ASSERT_EQ(secularRatios.size(), std::size(dts));
    EXPECT_LE(secularRatios.front(), secularRatios.back() * Real(2.0))
        << "secular-ratio GREW as dt shrank (h/2: " << secularRatios.front() << ", 4h: "
        << secularRatios.back() << ") -- the signature of a dt-independent pump surviving the dt->0 limit";
}

// ---------------------------------------------------------------------------
//  DR2 (LEMMA): per-body KE localization on a Free->Ball chain (both bodies
//  carry a quaternion DOF -- the maximal stress case for "propagated through
//  the articulated recursion", RobotIntegrator.hpp:223). No individual body's
//  KE trace may show a secular (monotone) trend under the shipped propagator;
//  an elastic/thermalizing exchange integrates to ~0 (LEMMA DR2).
// ---------------------------------------------------------------------------
TEST(FreeJointKEPump, DR2_ShippedPropagator_PerBodyKENotSecular_FreeBallChain) {
    const ConstraintSet cset;
    Rng rng(0x6101);
    RobotModel m = freeBallChain(rng);
    RobotState s;
    s.allocateFull(m);
    std::vector<AnalyticForceBridge> store;
    store.reserve(1);
    Rng r(0x6111);
    buildSeeded(m, s, r, /*k=*/60.0, store);

    const Real e0 = RobotEngine::calcKineticEnergy(m, s) + store.back().calcPotentialEnergy(s);
    const int nSteps = 4000;
    const NveTrace tr = runNve(m, s, store.back(), cset, kBaseH, nSteps);
    ASSERT_TRUE(tr.ok) << "a safe-step trajectory was rejected";

    const Dr1Summary total = dr1(tr.totalE, e0);
    EXPECT_LT(total.secularRatio, Real(0.2)) << "total-energy secular ratio " << total.secularRatio;

    Real maxBodySecular = 0;
    int maxBody = -1;
    for (int b = 1; b < m.numBodies; ++b) {
        const std::vector<Real>& ke = tr.keByBody[static_cast<std::size_t>(b)];
        // per-body KE is not centered at 0 (it's a strictly-positive quantity),
        // so measure the SECULAR component against the trajectory's own KE
        // fluctuation amplitude (max - min), not against a KE0 offset.
        Real keMin = ke[0], keMax = ke[0];
        for (Real v : ke) {
            keMin = std::min(keMin, v);
            keMax = std::max(keMax, v);
        }
        const Real amp = (keMax - keMin) + Real(1e-300);
        const Real slope = std::abs(lsSlope(ke)) * static_cast<Real>(ke.size());
        const Real secularRatio = slope / amp;
        std::fprintf(stderr, "[DR2 freeBallChain] body=%d keMin=%.6g keMax=%.6g secularRatio=%.4g\n", b,
                     (double)keMin, (double)keMax, (double)secularRatio);
        if (secularRatio > maxBodySecular) {
            maxBodySecular = secularRatio;
            maxBody = b;
        }
        EXPECT_LT(secularRatio, Real(0.35))
            << "body " << b
            << ": per-body KE shows a secular (monotone) trend -- a defect (DR2), not elastic exchange";
    }
    std::fprintf(stderr, "[DR2 freeBallChain] most-secular body=%d ratio=%.4g (root=body 1)\n", maxBody,
                 (double)maxBodySecular);
}

// ---------------------------------------------------------------------------
//  Positive control (exploratory, documentation-only -- NOT a build gate on
//  the shipped propagator). Runs the SAME single-Free-body fixture/seed/dt
//  sweep through the ablated linear-Taylor-plus-renormalize quaternion update
//  to determine whether the documented historical defect is empirically real
//  when the exp-map workaround is absent. This directly answers the "is the
//  ~1300 kJ/mol/traj pump real" question the diagnosis was asked to convict or
//  clear, WITHOUT touching production code.
// ---------------------------------------------------------------------------
TEST(FreeJointKEPump, Ablated_LinearTaylorQuat_PositiveControl_SingleFree) {
    const ConstraintSet cset;
    const Real dts[] = {kBaseH / 2, kBaseH, 2 * kBaseH, 4 * kBaseH};

    for (Real h : dts) {
        Rng rng(0x6001);
        RobotModel m = singleFree(rng);
        RobotState s;
        s.allocateFull(m);
        std::vector<AnalyticForceBridge> store;
        store.reserve(1);
        Rng r(0x6001 ^ 0xABCDu); // IDENTICAL start to the shipped-propagator DR1 test
        buildSeeded(m, s, r, /*k=*/60.0, store);

        const int nSteps = static_cast<int>(std::lround(static_cast<double>(kMatchedTime / h)));
        const Real e0 = RobotEngine::calcKineticEnergy(m, s) + store.back().calcPotentialEnergy(s);
        const NveTrace tr = runNveAblated(m, s, store.back(), cset, h, nSteps);
        if (!tr.ok) {
            std::fprintf(stderr, "[ablated singleFree] h=%.6g: step REJECTED (runaway)\n", (double)h);
            continue;
        }
        const Dr1Summary d = dr1(tr.totalE, e0);
        std::fprintf(stderr,
                     "[ablated singleFree] h=%.6g nSteps=%d maxAbsDrift=%.6g secularRatio=%.4g "
                     "(shipped-propagator bound: secularRatio<0.2)\n",
                     (double)h, nSteps, (double)d.maxAbs, (double)d.secularRatio);
    }
    SUCCEED() << "exploratory positive control -- see stderr for the measured ablated-path drift; "
                 "this test documents the measurement, it does not gate the shipped propagator "
                 "(see FreeJointKEPump.DR1_ShippedPropagator_* for the gate)";
}

// ---------------------------------------------------------------------------
//  Second positive control (exploratory, documentation-only): the single-body,
//  small-rotation-per-step probe above shows the shipped and ablated paths are
//  numerically indistinguishable (both O(dt^2), same magnitude) -- expected,
//  since exp-map and linear-Taylor agree to O(h^2) and only differ starting at
//  O(h^3), AND a lone Free body has no child to propagate an error into ("the
//  articulated recursion", RobotIntegrator.hpp:223). This probe widens BOTH
//  knobs the comment's claim depends on: (a) a Free->Ball CHAIN (a body for
//  the error to propagate into), and (b) a much LARGER per-step rotation
//  (bigger u, i.e. bigger theta = 0.5|w|h -- where exp-map's exact curvature
//  vs. linear-Taylor's flat extrapolation should start to separate), swept out
//  to the dt at which the shipped propagator's OWN corrector starts rejecting
//  steps (the practical ceiling), to see whether the ablated path degrades
//  FASTER (a real, if latent, defect) or in step with the shipped path
//  (indistinguishable in this harness).
// ---------------------------------------------------------------------------
TEST(FreeJointKEPump, Ablated_LinearTaylorQuat_PositiveControl_FreeBallChain_LargeRotationSweep) {
    const ConstraintSet cset;
    // large-rotation boost: pushes theta = 0.5|w|h well past the small-angle
    // regime at the coarser end of the sweep (unlike kBoostRT=2.5 used
    // elsewhere in this file, which keeps every dt in the safe/small-theta
    // regime by construction).
    constexpr Real kBigBoostRT = Real(80);
    const Real dts[] = {kBaseH, 2 * kBaseH, 4 * kBaseH, 8 * kBaseH, 16 * kBaseH, 32 * kBaseH};

    auto buildBigSeeded = [](const RobotModel& m, RobotState& s, Rng& rng, Real k,
                             std::vector<AnalyticForceBridge>& store) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        store.emplace_back(m, s, k);
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        std::vector<Real> w(static_cast<std::size_t>(m.nu)), u(static_cast<std::size_t>(m.nu));
        for (int i = 0; i < m.nu; ++i) {
            w[static_cast<std::size_t>(i)] = rng.gaussian();
        }
        RobotEngine::multiplyBySqrtMInv(m, s, w.data(), u.data());
        for (int i = 0; i < m.nu; ++i) {
            s.u()[i] = std::sqrt(kBigBoostRT) * u[static_cast<std::size_t>(i)];
        }
        seedDerivatives(m, s, store.back());
    };

    for (Real h : dts) {
        // shipped (exp-map)
        {
            Rng rng(0x6201);
            RobotModel m = freeBallChain(rng);
            RobotState s;
            s.allocateFull(m);
            std::vector<AnalyticForceBridge> store;
            store.reserve(1);
            Rng r(0x6211);
            buildBigSeeded(m, s, r, /*k=*/60.0, store);
            const Real e0 = RobotEngine::calcKineticEnergy(m, s) + store.back().calcPotentialEnergy(s);
            const int nSteps = 400;
            const NveTrace tr = runNve(m, s, store.back(), cset, h, nSteps);
            if (!tr.ok) {
                std::fprintf(stderr, "[shipped freeBallChain BIG] h=%.6g: REJECTED\n", (double)h);
            } else {
                const Dr1Summary d = dr1(tr.totalE, e0);
                std::fprintf(stderr, "[shipped  freeBallChain BIG] h=%.6g maxAbsDrift=%.6g secularRatio=%.4g\n",
                             (double)h, (double)d.maxAbs, (double)d.secularRatio);
            }
        }
        // ablated (linear-Taylor + renormalize)
        {
            Rng rng(0x6201);
            RobotModel m = freeBallChain(rng);
            RobotState s;
            s.allocateFull(m);
            std::vector<AnalyticForceBridge> store;
            store.reserve(1);
            Rng r(0x6211); // IDENTICAL start to the shipped arm above
            buildBigSeeded(m, s, r, /*k=*/60.0, store);
            const Real e0 = RobotEngine::calcKineticEnergy(m, s) + store.back().calcPotentialEnergy(s);
            const int nSteps = 400;
            const NveTrace tr = runNveAblated(m, s, store.back(), cset, h, nSteps);
            if (!tr.ok) {
                std::fprintf(stderr, "[ablated freeBallChain BIG] h=%.6g: REJECTED\n", (double)h);
            } else {
                const Dr1Summary d = dr1(tr.totalE, e0);
                std::fprintf(stderr, "[ablated freeBallChain BIG] h=%.6g maxAbsDrift=%.6g secularRatio=%.4g\n",
                             (double)h, (double)d.maxAbs, (double)d.secularRatio);
            }
        }
    }
    SUCCEED() << "exploratory large-rotation positive control -- see stderr for shipped-vs-ablated "
                 "drift at each dt; documentation only, not a build gate";
}
