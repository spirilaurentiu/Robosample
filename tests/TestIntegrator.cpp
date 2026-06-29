// ============================================================================
//  TestIntegratorSmoke.cpp -- Phase 1 keystone proof: drive the REAL, now-
//  templated RobotEngine::verletStep with the OpenMM-free AnalyticForceBridge.
//
//  This is the payoff of templating the integrator: the production velocity-
//  Verlet step runs end-to-end against an analytic potential, no OpenMM, no
//  faking. It is a smoke test (full energy-drift / reversibility goldens are
//  Phase 5); here we assert the step instantiates, runs, returns success, keeps
//  the state finite, and -- since AnalyticForceBridge is a genuine conservative
//  harmonic system integrated by a symplectic step -- that total energy KE+U
//  does not drift more than a few percent over a short trajectory at small h.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

#include "AnalyticForceBridge.hpp"
#include "Constraints.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotIntegrator.hpp" // the templated verletStep/stepTo/checkReversibility
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

RobotModel chain(Rng& rng) {
    auto mk = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        s.mass = rng.uniform(0.8, 1.6);
        s.com_B = rng.vec3(-0.1, 0.1);
        s.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
        return s;
    };
    RobotModel m = buildForest({mk(0, JointType::Free), mk(1, JointType::Torsion)});
    attachAtoms(m, 1, {Vec3(0.0, 0.0, 0.0), Vec3(0.10, 0.02, -0.03)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    return m;
}

// Seed the derivative chain so the first verletStep has valid qdot0/udot0/qdd0.
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

// ---------------------------------------------------------------------------
//  Phase 5 fixtures: a single Free body, and a Free-rooted 3-body chain. Both
//  carry real atoms so the harmonic AnalyticForceBridge has something to pull on.
// ---------------------------------------------------------------------------
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

RobotModel chain3(Rng& rng) {
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
    RobotModel m =
        buildForest({mk(0, JointType::Free), mk(1, JointType::Torsion), mk(2, JointType::Torsion)});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.10, 0.02, -0.03)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0), Vec3(-0.03, 0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0, 0.08, -0.04), Vec3(0.04, 0, 0.05)}, {16.0, 1.0});
    return m;
}

// A Free root with a Ball child: two quaternion bodies, for the |q|=1 test.
RobotModel freeBall(Rng& rng) {
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

// A forest exercising the joints whose headers warn they are unvalidated beyond
// the kernel finite-difference checks (JointKernels.hpp): BendStretch and
// SphericalCoords (internal 1-3/1-4 bond-angle + length DOF) as children of a
// Free root, plus a FreeLine root (5-DOF: quaternion minus roll + 2 translation).
// Energy conservation here is the dynamics-level proof those joints' H_FM/HDot
// assemble into a symplectic step -- the "energy conservation before trusting
// these joints" the kernel header explicitly demands.
RobotModel underValidatedJoints(Rng& rng) {
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
    // tree 1: Free -> BendStretch -> SphericalCoords ; tree 2: FreeLine root.
    RobotModel m = buildForest({mk(0, JointType::Free),
                                mk(1, JointType::BendStretch),
                                mk(2, JointType::SphericalCoords),
                                mk(0, JointType::FreeLine)});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.10, 0.02, -0.03)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0), Vec3(-0.03, 0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0, 0.08, -0.04), Vec3(0.04, 0, 0.05)}, {16.0, 1.0});
    attachAtoms(m, 4, {Vec3(0.06, -0.05, 0.02), Vec3(-0.02, 0.07, 0.0)}, {15.0, 1.0});
    return m;
}

// thermal scale for the momentum draw (arbitrary positive constant; only sets
// the energy scale, not correctness).
constexpr Real kBoostRT = Real(2.5);

// Build a ready-to-integrate state: anchors captured at a first random pose, then
// the body is displaced to a SECOND random pose (so U != 0), and u is seeded the
// production way -- u = sqrt(RT) * sqrt(M^-1) * gaussian (multiplyBySqrtMInv).
// Returns the bridge by value-construction into `bridgeStore`. After this the
// derivative chain is seeded and the first verletStep is valid.
void buildSeeded(const RobotModel& m,
                 RobotState& s,
                 Rng& rng,
                 Real k,
                 std::vector<AnalyticForceBridge>& bridgeStore) {
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    bridgeStore.emplace_back(m, s, k); // anchors = this pose
    AnalyticForceBridge& bridge = bridgeStore.back();

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

    seedDerivatives(m, s, bridge);
}

// least-squares slope of y vs index 0..n-1 (for the monotone-drift check).
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

// max |E_t - E_0| / |E_0| and the linear-fit total drift slope*N/|E_0| over a
// trajectory; returns false if any step is rejected.
struct DriftResult {
    bool ok = true;
    Real maxRel = 0;
    Real slopeTotal = 0;
    Real maxAbs = 0;
};
DriftResult runTrajectory(const RobotModel& m,
                          RobotState& s,
                          AnalyticForceBridge& bridge,
                          const ConstraintSet& cset,
                          Real h,
                          int nSteps) {
    DriftResult r;
    const Real E0 = RobotEngine::calcKineticEnergy(m, s) + bridge.calcPotentialEnergy(s);
    std::vector<Real> energies;
    energies.reserve(static_cast<std::size_t>(nSteps));
    for (int n = 0; n < nSteps; ++n) {
        if (!RobotEngine::verletStep(m, s, bridge, cset, h)) {
            r.ok = false;
            return r;
        }
        const Real E = RobotEngine::calcKineticEnergy(m, s) + bridge.calcPotentialEnergy(s);
        energies.push_back(E);
        r.maxRel = std::max(r.maxRel, std::abs(E - E0) / (std::abs(E0) + Real(1e-12)));
        r.maxAbs = std::max(r.maxAbs, std::abs(E - E0));
    }
    r.slopeTotal = std::abs(lsSlope(energies)) * static_cast<Real>(nSteps) / (std::abs(E0) + Real(1e-12));
    return r;
}

// ---- Phase 5 tolerances ----
constexpr Real kEdrift = Real(1e-4);   // relative energy drift at the safe step
constexpr Real kRevTight = Real(1e-7); // round-trip reversibility residual at safe h
constexpr Real kSafeH = Real(5e-4);    // a step well inside the reversible regime

} // namespace

// ---------------------------------------------------------------------------
//  The templated verletStep instantiates with a non-ForceBridge type, runs, and
//  returns success on a benign step.
// ---------------------------------------------------------------------------
TEST(IntegratorSmoke, VerletStepRunsWithAnalyticBridge) {
    Rng rng(0x5709);
    RobotModel m = chain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    AnalyticForceBridge bridge(m, s, /*k=*/120.0);
    const ConstraintSet cset; // empty -> unconstrained

    // small displacement + gentle velocities so the step is benign
    randomizeState(m, s, rng);
    Real* u = s.u();
    for (int i = 0; i < m.nu; ++i) {
        u[i] *= 0.05;
    }
    seedDerivatives(m, s, bridge);

    const bool ok = RobotEngine::verletStep(m, s, bridge, cset, /*h=*/1e-3);
    EXPECT_TRUE(ok) << "templated verletStep rejected a benign step";
    for (int i = 0; i < m.nq; ++i) {
        EXPECT_TRUE(std::isfinite(s.q()[i]));
    }
    for (int i = 0; i < m.nu; ++i) {
        EXPECT_TRUE(std::isfinite(s.u()[i]));
    }
}

// ---------------------------------------------------------------------------
//  stepTo (also templated) drives a short trajectory; total energy KE+U stays
//  bounded (symplectic Verlet on a smooth harmonic potential).
// ---------------------------------------------------------------------------
TEST(IntegratorSmoke, EnergyBoundedOverShortTrajectory) {
    Rng rng(0x570A);
    RobotModel m = chain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    AnalyticForceBridge bridge(m, s, 120.0);
    const ConstraintSet cset;

    // displace off the anchors and give a modest velocity
    randomizeState(m, s, rng);
    Real* u = s.u();
    for (int i = 0; i < m.nu; ++i) {
        u[i] *= 0.1;
    }
    seedDerivatives(m, s, bridge);

    const Real h = 5e-4;
    const Real E0 = RobotEngine::calcKineticEnergy(m, s) + bridge.calcPotentialEnergy(s);

    Real maxRelDrift = 0;
    int steps = 0;
    for (int n = 0; n < 200; ++n) {
        const bool ok = RobotEngine::verletStep(m, s, bridge, cset, h);
        ASSERT_TRUE(ok) << "step " << n << " was rejected";
        ++steps;
        const Real E = RobotEngine::calcKineticEnergy(m, s) + bridge.calcPotentialEnergy(s);
        maxRelDrift = std::max(maxRelDrift, std::abs(E - E0) / (std::abs(E0) + Real(1e-12)));
    }
    EXPECT_EQ(steps, 200);
    // symplectic Verlet on a smooth potential: bounded, small drift (not secular
    // growth). Loose bound -- this is a smoke test, not the Phase 5 energy golden.
    EXPECT_LT(maxRelDrift, 0.05) << "energy drifted " << maxRelDrift << " over 200 steps";
}

// ===========================================================================
//  Phase 5 -- verlet / reversibility / stepTo goldens. Production momentum draw
//  (multiplyBySqrtMInv of a Gaussian), no constraints (empty ConstraintSet).
// ===========================================================================

// ---------------------------------------------------------------------------
//  P5.1: total energy KE+U stays bounded with NO secular (monotone) drift over a
//        long trajectory at a safe step, for both the single Free body and the
//        3-chain. FAIL sensitivity: at 10x the step, the drift blows past kEdrift.
// ---------------------------------------------------------------------------
TEST(Integrator, EnergyConservationBounded) {
    const ConstraintSet cset; // unconstrained

    struct Case {
        const char* name;
        RobotModel (*build)(Rng&);
        unsigned seed;
    };
    const Case cases[] = {
        {"Free", &singleFree, 0x5101},
        {"Chain3", &chain3, 0x5102},
    };

    for (const Case& c : cases) {
        Rng rng(c.seed);
        RobotModel m = c.build(rng);

        // safe step: bounded, non-secular drift.
        {
            RobotState s;
            s.allocateFull(m);
            std::vector<AnalyticForceBridge> store;
            store.reserve(1);
            Rng r(c.seed ^ 0xABCDu);
            buildSeeded(m, s, r, /*k=*/60.0, store);
            const DriftResult d = runTrajectory(m, s, store.back(), cset, kSafeH, 2000);
            ASSERT_TRUE(d.ok) << c.name << ": a safe-step trajectory was rejected";
            EXPECT_LT(d.maxRel, kEdrift) << c.name << ": energy drift " << d.maxRel;
            EXPECT_LT(d.slopeTotal, kEdrift)
                << c.name << ": secular (monotone) drift slope*N/|E0| = " << d.slopeTotal;
        }

        // FAIL sensitivity: 10x the step must drive the drift past kEdrift, proving
        // the test can actually SEE a bad step (it is not vacuously loose).
        {
            RobotState s;
            s.allocateFull(m);
            std::vector<AnalyticForceBridge> store;
            store.reserve(1);
            Rng r(c.seed ^ 0xABCDu); // SAME start as the safe run
            buildSeeded(m, s, r, 60.0, store);
            const DriftResult d = runTrajectory(m, s, store.back(), cset, 10 * kSafeH, 2000);
            // either it was rejected (also a clear "bad step" signal) or the drift
            // exceeded the budget.
            const bool sawBadStep = (!d.ok) || (d.maxRel > kEdrift);
            EXPECT_TRUE(sawBadStep) << c.name << ": 10x step stayed under kEdrift (maxRel=" << d.maxRel
                                    << ", ok=" << d.ok << ") -- test is insensitive";
        }
    }
}

// ---------------------------------------------------------------------------
//  P5.2: velocity-Verlet is globally second order -- the energy drift at h is
//        ~4x the drift at h/2 over the SAME physical time from the SAME start.
// ---------------------------------------------------------------------------
TEST(Integrator, SecondOrderEnergyError) {
    const ConstraintSet cset;

    struct Case {
        const char* name;
        RobotModel (*build)(Rng&);
        unsigned seed;
    };
    const Case cases[] = {
        {"Free", &singleFree, 0x5201},
        {"Chain3", &chain3, 0x5202},
    };

    for (const Case& c : cases) {
        Rng rng(c.seed);
        RobotModel m = c.build(rng);

        auto driftAt = [&](Real h, int nSteps) -> Real {
            RobotState s;
            s.allocateFull(m);
            std::vector<AnalyticForceBridge> store;
            store.reserve(1);
            Rng r(c.seed ^ 0x55AAu); // fixed seed: identical start for h and h/2
            buildSeeded(m, s, r, 60.0, store);
            const DriftResult d = runTrajectory(m, s, store.back(), cset, h, nSteps);
            return d.ok ? d.maxAbs : std::numeric_limits<Real>::infinity();
        };

        // same physical time: 600 steps at h, 1200 at h/2.
        const Real driftH = driftAt(kSafeH, 600);
        const Real driftHalf = driftAt(kSafeH / 2, 1200);
        ASSERT_GT(driftHalf, 0.0) << c.name;
        const Real ratio = driftH / driftHalf;
        EXPECT_GE(ratio, 3.0) << c.name << ": O(h^2) ratio too small (" << ratio << ")";
        EXPECT_LE(ratio, 5.0) << c.name << ": O(h^2) ratio too large (" << ratio << ")";
    }
}

// ---------------------------------------------------------------------------
//  P5.3: the engine's own checkReversibility (which no other test invokes) gives
//        a tiny round-trip residual at a safe h. FAIL sensitivity: at a large h
//        the residual is O(1) -- documenting the header's configuration-dependent
//        safe-h claim. (Uses the 3-chain: a single Free body's exp-map advance is
//        reversible even at large h, so the chain is what exposes the breakdown.)
// ---------------------------------------------------------------------------
TEST(Integrator, TimeReversibilityResidualSmall) {
    const ConstraintSet cset;

    Rng rng(0x5301);
    RobotModel m = chain3(rng);
    RobotState s;
    s.allocateFull(m);
    std::vector<AnalyticForceBridge> store;
    store.reserve(1);
    Rng r(0x5311);
    buildSeeded(m, s, r, 60.0, store);

    // checkReversibility is non-destructive, so both calls run from the same state.
    const Real safe = RobotEngine::checkReversibility(m, s, store.back(), cset, /*nSteps=*/200, kSafeH);
    EXPECT_LT(safe, kRevTight) << "reversible step gave residual " << safe;

    const Real big = RobotEngine::checkReversibility(m, s, store.back(), cset, 200, /*h=*/5e-2);
    EXPECT_GT(big, 0.1) << "large-h residual " << big << " is not O(1) -- safe-h claim not exercised";
}

// ---------------------------------------------------------------------------
//  P5.4: every quaternion q-block of every Free/Ball body stays unit-norm over a
//        long trajectory. The exp-map advance preserves |q|=1 by construction;
//        this guards against a regression to linear-Taylor-without-renormalize.
// ---------------------------------------------------------------------------
TEST(Integrator, QuaternionStaysUnitOverTrajectory) {
    const ConstraintSet cset;

    Rng rng(0x5401);
    RobotModel m = freeBall(rng); // a Free body and a Ball body
    RobotState s;
    s.allocateFull(m);
    std::vector<AnalyticForceBridge> store;
    store.reserve(1);
    Rng r(0x5411);
    buildSeeded(m, s, r, 60.0, store);
    AnalyticForceBridge& bridge = store.back();

    bool anyQuat = false;
    for (int b = 1; b < m.numBodies; ++b) {
        anyQuat = anyQuat || m.isQuaternionBody(b);
    }
    ASSERT_TRUE(anyQuat) << "fixture has no quaternion body -- test is vacuous";

    for (int n = 0; n < 1500; ++n) {
        ASSERT_TRUE(RobotEngine::verletStep(m, s, bridge, cset, kSafeH)) << "step " << n << " rejected";
        for (int b = 1; b < m.numBodies; ++b) {
            if (!m.isQuaternionBody(b)) {
                continue;
            }
            const int q = m.bodyQIndex[b];
            const Real norm = std::sqrt(s.q()[q] * s.q()[q] + s.q()[q + 1] * s.q()[q + 1]
                                        + s.q()[q + 2] * s.q()[q + 2] + s.q()[q + 3] * s.q()[q + 3]);
            ASSERT_NEAR(norm, 1.0, rtest::kTight)
                << "body " << b << " quaternion lost unit norm at step " << n;
        }
    }
}

// ---------------------------------------------------------------------------
//  P5.6: energy conservation + reversibility for the joints whose headers warn
//        they are unvalidated beyond the kernel FD checks: BendStretch,
//        SphericalCoords (1-3/1-4 bond-angle + length DOF) and FreeLine. This is
//        the dynamics-level proof the kernel header demands -- it ensures these
//        DOF types are integrated without secular energy pumping, the practical
//        failure the NCMC gap= diagnostic watches for.
// ---------------------------------------------------------------------------
TEST(Integrator, UnderValidatedJointsConserveEnergy) {
    const ConstraintSet cset;
    Rng rng(0x5601);
    RobotModel m = underValidatedJoints(rng);

    // safe step: bounded, non-secular drift.
    {
        RobotState s;
        s.allocateFull(m);
        std::vector<AnalyticForceBridge> store;
        store.reserve(1);
        Rng r(0x5611);
        buildSeeded(m, s, r, 60.0, store);
        const DriftResult d = runTrajectory(m, s, store.back(), cset, kSafeH, 2000);
        ASSERT_TRUE(d.ok) << "a safe-step trajectory over BendStretch/SphericalCoords/FreeLine was rejected";
        EXPECT_LT(d.maxRel, kEdrift) << "energy drift " << d.maxRel;
        EXPECT_LT(d.slopeTotal, kEdrift) << "secular (monotone) energy pumping slope*N/|E0| = " << d.slopeTotal;
    }

    // O(h^2): halving the step quarters the drift over the same physical time.
    {
        auto driftAt = [&](Real h, int nSteps) -> Real {
            RobotState s;
            s.allocateFull(m);
            std::vector<AnalyticForceBridge> store;
            store.reserve(1);
            Rng r(0x5622);
            buildSeeded(m, s, r, 60.0, store);
            const DriftResult d = runTrajectory(m, s, store.back(), cset, h, nSteps);
            return d.ok ? d.maxAbs : std::numeric_limits<Real>::infinity();
        };
        const Real ratio = driftAt(kSafeH, 600) / driftAt(kSafeH / 2, 1200);
        EXPECT_GE(ratio, 3.0) << "O(h^2) ratio too small (" << ratio << ")";
        EXPECT_LE(ratio, 5.0) << "O(h^2) ratio too large (" << ratio << ")";
    }

    // round-trip reversibility residual is tiny at a safe step.
    {
        RobotState s;
        s.allocateFull(m);
        std::vector<AnalyticForceBridge> store;
        store.reserve(1);
        Rng r(0x5633);
        buildSeeded(m, s, r, 60.0, store);
        const Real res = RobotEngine::checkReversibility(m, s, store.back(), cset, /*nSteps=*/200, kSafeH);
        EXPECT_LT(res, kRevTight) << "reversibility residual " << res << " on under-validated joints";
    }
}

// ---------------------------------------------------------------------------
//  P5.5: stepTo is deterministic -- two identical states advanced over the same
//        schedule produce BITWISE-identical q,u (mirrors Stability.Deterministic-
//        Generators). Steps are genuinely taken (asserted), not rejected no-ops.
// ---------------------------------------------------------------------------
TEST(Integrator, StepToReachesEndTimeDeterministically) {
    const ConstraintSet cset;

    Rng rng(0x5501);
    RobotModel m = chain3(rng);
    RobotState s1;
    s1.allocateFull(m);
    std::vector<AnalyticForceBridge> store;
    store.reserve(1);
    Rng r(0x5511);
    buildSeeded(m, s1, r, 60.0, store);
    AnalyticForceBridge& bridge = store.back();

    // exact clone of the seeded state
    RobotState s2;
    s2.allocateFull(m);
    std::copy(s1.q(), s1.q() + m.nq, s2.q());
    std::copy(s1.u(), s1.u() + m.nu, s2.u());
    std::copy(s1.qdot(), s1.qdot() + m.nq, s2.qdot());
    std::copy(s1.udot(), s1.udot() + m.nu, s2.udot());
    std::copy(s1.qdotdot(), s1.qdotdot() + m.nq, s2.qdotdot());
    s2.time = s1.time;

    int taken1 = 0, taken2 = 0;
    for (int n = 0; n < 200; ++n) {
        if (RobotEngine::stepTo(m, s1, bridge, cset, s1.time + kSafeH)) {
            ++taken1;
        }
    }
    for (int n = 0; n < 200; ++n) {
        if (RobotEngine::stepTo(m, s2, bridge, cset, s2.time + kSafeH)) {
            ++taken2;
        }
    }
    ASSERT_EQ(taken1, 200) << "stepTo rejected steps -- determinism test would be vacuous";
    ASSERT_EQ(taken2, 200);

    for (int i = 0; i < m.nq; ++i) {
        EXPECT_EQ(s1.q()[i], s2.q()[i]) << "q[" << i << "] diverged -- stepTo not deterministic";
    }
    for (int i = 0; i < m.nu; ++i) {
        EXPECT_EQ(s1.u()[i], s2.u()[i]) << "u[" << i << "] diverged -- stepTo not deterministic";
    }
}

// A test-only bridge that injects a non-finite body force, to drive the
// reject-and-restore path of verletStep without OpenMM. Same evaluate() contract.
struct InfForceBridge {
    const RobotModel& model;
    explicit InfForceBridge(const RobotModel& m)
        : model(m) {
    }
    void evaluate(RobotState& s) {
        SpatialVec* BF = s.bodyForceG();
        for (int b = 0; b < model.numBodies; ++b) {
            BF[b] = SpatialVec(Vec3(0), Vec3(0));
        }
        Real* mob = s.mobilityForce();
        for (int i = 0; i < model.nu; ++i) {
            mob[i] = 0;
        }
        BF[1].linear[0] = std::numeric_limits<Real>::infinity(); // poison one body force
    }
};

// ---------------------------------------------------------------------------
//  P5.6: a bad step (non-finite force, OR an absurdly large h) makes verletStep
//        return false AND leave q,u EXACTLY equal to the pre-step snapshot. The
//        restorePreStep contract: the caller must never be handed a NaN to reject.
// ---------------------------------------------------------------------------
TEST(Integrator, BadStepReturnsFalseAndRestores) {
    const ConstraintSet cset;

    // ---- (a) non-finite force ----
    {
        Rng rng(0x5601);
        RobotModel m = chain3(rng);
        RobotState s;
        s.allocateFull(m);
        std::vector<AnalyticForceBridge> store;
        store.reserve(1);
        Rng r(0x5611);
        buildSeeded(m, s, r, 60.0, store); // valid derivative seed via the finite bridge

        std::vector<Real> q0(s.q(), s.q() + m.nq), u0(s.u(), s.u() + m.nu);
        InfForceBridge bad(m);
        const bool ok = RobotEngine::verletStep(m, s, bad, cset, kSafeH);
        EXPECT_FALSE(ok) << "verletStep accepted a step with a non-finite force";
        for (int i = 0; i < m.nq; ++i) {
            EXPECT_EQ(s.q()[i], q0[static_cast<std::size_t>(i)])
                << "q[" << i << "] not restored after bad step";
        }
        for (int i = 0; i < m.nu; ++i) {
            EXPECT_EQ(s.u()[i], u0[static_cast<std::size_t>(i)])
                << "u[" << i << "] not restored after bad step";
        }
    }

    // ---- (b) absurd step size (CFL runaway -> velocitiesSane fails) ----
    {
        Rng rng(0x5602);
        RobotModel m = chain3(rng);
        RobotState s;
        s.allocateFull(m);
        std::vector<AnalyticForceBridge> store;
        store.reserve(1);
        Rng r(0x5612);
        buildSeeded(m, s, r, 60.0, store);

        std::vector<Real> q0(s.q(), s.q() + m.nq), u0(s.u(), s.u() + m.nu);
        const bool ok = RobotEngine::verletStep(m, s, store.back(), cset, /*h=*/5.0);
        EXPECT_FALSE(ok) << "verletStep accepted an absurd-h runaway";
        for (int i = 0; i < m.nq; ++i) {
            EXPECT_EQ(s.q()[i], q0[static_cast<std::size_t>(i)])
                << "q[" << i << "] not restored after runaway";
        }
        for (int i = 0; i < m.nu; ++i) {
            EXPECT_EQ(s.u()[i], u0[static_cast<std::size_t>(i)])
                << "u[" << i << "] not restored after runaway";
        }
    }
}