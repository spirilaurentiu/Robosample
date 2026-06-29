// ============================================================================
//  TestAnalyticForce.cpp -- Phase 1: self-tests for the AnalyticForceBridge
//  harness (the keystone that makes integrator/energy/ensemble tests OpenMM-free).
//
//  The harness is only trustworthy as an oracle if (a) its force really is the
//  negative gradient of its stated potential, and (b) its per-atom -> per-body
//  reduction is byte-identical to the production ForceBridge reduction. These
//  tests pin both, each with an explicit must-FAIL control.
// ============================================================================
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "AnalyticForceBridge.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
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

RobotModel twoBodyWithAtoms(Rng& rng) {
    auto mk = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        s.mass = rng.uniform(0.5, 2.0);
        s.com_B = rng.vec3(-0.2, 0.2);
        s.inertia_B = UnitInertia(rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7));
        return s;
    };
    RobotModel m = buildForest({mk(0, JointType::Free), mk(1, JointType::Torsion)});
    attachAtoms(m, 1, {Vec3(0.0, 0.0, 0.0), Vec3(0.10, -0.03, 0.05)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.0, 0.12, 0.0), Vec3(-0.04, 0.0, 0.07)}, {16.0, 1.0});
    return m;
}

} // namespace

// ---------------------------------------------------------------------------
//  ForceMatchesNegGradU: per-atom force == -dU/dr_a by central FD on the
//  Cartesian atom position. U is a pure function of the flat atom-position
//  array (calcPotentialEnergyAt), so each component is perturbed independently.
//  This certifies the harness IS conservative -- the property energy tests need.
// ---------------------------------------------------------------------------
TEST(AnalyticForce, ForceMatchesNegGradU) {
    Rng rng(0xA1F0);
    RobotModel m = twoBodyWithAtoms(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    AnalyticForceBridge bridge(m, s, /*k=*/250.0); // anchors = these start positions

    // move OFF the anchors so forces are nonzero
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    std::vector<Vec3> base(s.atomPosG(), s.atomPosG() + m.numAtoms);
    const Real h = 1e-6;
    for (int a = 0; a < m.numAtoms; ++a) {
        if (m.atomMass[a] == Real(0)) {
            continue;
        }
        const Vec3 F = bridge.atomForceAt(base.data(), a);
        for (int i = 0; i < 3; ++i) {
            std::vector<Vec3> p = base;
            p[static_cast<std::size_t>(a)][i] += h;
            const Real up = bridge.calcPotentialEnergyAt(p.data(), m.numAtoms);
            p[static_cast<std::size_t>(a)][i] -= 2 * h;
            const Real um = bridge.calcPotentialEnergyAt(p.data(), m.numAtoms);
            const Real dUdr = (up - um) / (2 * h);
            EXPECT_NEAR(-dUdr, F[i], rtest::kFD) << "atom " << a << " comp " << i;
        }
    }
}

// ---------------------------------------------------------------------------
//  ReductionMatchesRealBridge: bodyForceG from evaluate() equals the hand
//  reduction (sum f ; sum (r_a - X_GB[b].p()) % f) to kTight.
//  MUST-FAIL guard: the moment about X_BM.p() (a constant body-frame point)
//  instead of the body origin must NOT reproduce the bridge's angular row.
// ---------------------------------------------------------------------------
TEST(AnalyticForce, ReductionMatchesRealBridge) {
    Rng rng(0xA1F1);
    RobotModel m = twoBodyWithAtoms(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    AnalyticForceBridge bridge(m, s, 180.0);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    bridge.evaluate(s);

    const SpatialVec* BF = s.bodyForceG();
    const Vec3* posG = s.atomPosG();
    const Transform* X_GB = s.X_GB();

    std::vector<Vec3> lin(static_cast<std::size_t>(m.numBodies), Vec3(0));
    std::vector<Vec3> ang(static_cast<std::size_t>(m.numBodies), Vec3(0));
    std::vector<Vec3> angWrong(static_cast<std::size_t>(m.numBodies), Vec3(0));
    for (int a = 0; a < m.numAtoms; ++a) {
        if (m.atomMass[a] == Real(0)) {
            continue;
        }
        const int b = m.atomBody[a];
        const Vec3 f = bridge.atomForce(s, a);
        lin[static_cast<std::size_t>(b)] += f;
        ang[static_cast<std::size_t>(b)] += (posG[a] - X_GB[b].p()) % f;        // correct: body origin
        angWrong[static_cast<std::size_t>(b)] += (posG[a] - m.X_BM[b].p()) % f; // wrong lever
    }

    bool anyAngDiffersFromWrong = false;
    for (int b = 1; b < m.numBodies; ++b) {
        EXPECT_TRUE(rtest::NearVec3(BF[b][1], lin[static_cast<std::size_t>(b)], rtest::kTight))
            << "linear body " << b;
        EXPECT_TRUE(rtest::NearVec3(BF[b][0], ang[static_cast<std::size_t>(b)], rtest::kTight))
            << "angular body " << b;
        if (!rtest::NearVec3(BF[b][0], angWrong[static_cast<std::size_t>(b)], 1e-9)) {
            anyAngDiffersFromWrong = true;
        }
    }
    // must-fail control: the wrong-origin moment must differ on at least one body
    // (a robot whose every atom sits at X_BM would be degenerate; ours does not).
    EXPECT_TRUE(anyAngDiffersFromWrong)
        << "angular row matched the WRONG lever origin (X_BM.p) -- convention not pinned";
}

// ---------------------------------------------------------------------------
//  VirtualSiteSkipped: an atomMass==0 site contributes ZERO to bodyForceG.
//  MUST-FAIL guard: the force the site WOULD carry if (wrongly) counted is
//  nonzero, so including it would change the body's linear row.
// ---------------------------------------------------------------------------
TEST(AnalyticForce, VirtualSiteSkipped) {
    Rng rng(0xA1F2);
    BodySpec root;
    root.parent = 0;
    root.joint = JointType::Free;
    root.X_PF = Transform(rng.rotation(), rng.vec3());
    root.X_BM = Transform(rng.rotation(), rng.vec3());
    root.mass = 1.0;
    RobotModel m = buildForest({root});
    // one real atom (mass 12) and one VIRTUAL site (mass 0)
    attachAtoms(m, 1, {Vec3(0.05, 0.0, 0.0), Vec3(0.0, 0.07, 0.0)}, {12.0, 0.0});

    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    AnalyticForceBridge bridge(m, s, 200.0);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    bridge.evaluate(s);

    const SpatialVec* BF = s.bodyForceG();
    const Vec3* posG = s.atomPosG();
    const Transform* X_GB = s.X_GB();

    // expected: only the REAL atom (index 0) contributes
    const Vec3 fReal = bridge.atomForce(s, 0);
    EXPECT_TRUE(rtest::NearVec3(BF[1][1], fReal, rtest::kTight));
    EXPECT_TRUE(rtest::NearVec3(BF[1][0], (posG[0] - X_GB[1].p()) % fReal, rtest::kTight));

    // the site genuinely carries no force
    EXPECT_TRUE(rtest::NearVec3(bridge.atomForce(s, 1), Vec3(0), rtest::kTight));

    // must-fail control: the site IS displaced from its anchor, so the harmonic
    // force it WOULD carry is nonzero -- including it would change BF[1].linear.
    const Vec3 fSiteIfCounted = (bridge.anchor(1) - posG[1]) * bridge.stiffness();
    EXPECT_GT(fSiteIfCounted.norm(), 1e-6) << "site not displaced -> test is vacuous; reseed";
    EXPECT_FALSE(rtest::NearVec3(BF[1][1], fReal + fSiteIfCounted, rtest::kTight))
        << "body force includes the virtual-site contribution -- skip was lost";
}