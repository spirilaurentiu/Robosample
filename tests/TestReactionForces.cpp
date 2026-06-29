// ============================================================================
//  TestReactionForces.cpp -- the newly-ported mobilizer reaction-force operator
//  RobotEngine::calcMobilizerReactionForces / findMobilizerReactionOnBodyAtM-
//  InGround (the SimTK calcMobilizerReactionForces functionality that was
//  previously documented "NOT PORTED -- no operator in Robosample").
//
//  The reaction on body b is the spatial force its inboard mobilizer transmits
//  to it, in Ground, from a rigid Newton-Euler inward sweep on the TRUE body
//  accelerations:
//      reac_b@Bo = Mk_b A_GB_b + gyro_b - F_ext_b + sum_c Phi[c] reac_c@Bo
//  with F_ext_b = bodyForceG[b] + sum_j H_j mobilityForce[uOff+j]. The public
//  output is reported at the outboard M frame (Simbody's convention).
//
//  Oracle strategy: no OpenMM. Forces come from the analytic harmonic
//  AnalyticForceBridge, the full dynamics chain is realized (position, velocity,
//  ABI, calcUDot), and the operator is checked against (a) the defining rigid
//  recursion recomputed independently, (b) the physical "free mobilizer
//  transmits no constraint force" law, and (c) a whole-tree Newton-Euler balance.
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

// A forest rooted by `root`, with two torsion children carrying atoms, so the
// bridge has something to pull on and the joints transmit nontrivial force.
RobotModel chain(Rng& rng, JointType root) {
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
    RobotModel m = buildForest({mk(0, root), mk(1, JointType::Torsion), mk(2, JointType::Torsion)});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.10, 0.02, -0.03)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0), Vec3(-0.03, 0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0, 0.08, -0.04), Vec3(0.04, 0, 0.05)}, {16.0, 1.0});
    return m;
}

// Realize the full dynamics chain so A_GB / gyro / Mk_G / bodyForceG / mobilityForce
// are all valid (the operator's precondition).
void realizeDynamics(const RobotModel& m, RobotState& s, AnalyticForceBridge& bridge) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    bridge.evaluate(s); // fills bodyForceG (mobilityForce zeroed inside)
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    RobotEngine::calcUDot(m, s); // A_GB, udot
}

// shift a spatial force [t;f] applied at point P (Ground) to the global origin:
// torque_about_O = t + P x f ; force unchanged.
SpatialVec shiftToOrigin(const SpatialVec& F, const Vec3& P) {
    return SpatialVec(F.angular + (P % F.linear), F.linear);
}

} // namespace

// ---------------------------------------------------------------------------
//  1. ReactionMatchesRigidRecursion: the operator reproduces the defining rigid
//     Newton-Euler inward sweep recomputed here independently, body for body, to
//     kTight. Pins the algebra of the port.
// ---------------------------------------------------------------------------
TEST(ReactionForces, ReactionMatchesRigidRecursion) {
    Rng rng(0x7A01);
    RobotModel m = chain(rng, JointType::Free);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 10; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        AnalyticForceBridge bridge(m, s, 150.0);
        randomizeState(m, s, rng);
        realizeDynamics(m, s, bridge);

        std::vector<SpatialVec> reacBo(static_cast<std::size_t>(m.numBodies));
        std::vector<SpatialVec> reacM(static_cast<std::size_t>(m.numBodies));
        RobotEngine::calcMobilizerReactionForces(m, s, reacBo.data(), reacM.data());

        // independent recomputation of the same recursion
        const SpatialInertia* Mk = s.Mk_G();
        const SpatialVec* A = s.A_GB();
        const SpatialVec* gyro = s.gyro();
        const SpatialVec* bF = s.bodyForceG();
        const Real* mF = s.mobilityForce();
        const SpatialVec* H = s.H();
        const PhiMatrix* Phi = s.Phi();
        const Transform* X_GB = s.X_GB();

        std::vector<SpatialVec> ref(static_cast<std::size_t>(m.numBodies), SpatialVec(Vec3(0), Vec3(0)));
        for (int b = m.numBodies - 1; b >= 1; --b) {
            const int uOff = m.bodyUIndex[b], dof = m.bodyNU[b];
            SpatialVec fExt = bF[b];
            for (int j = 0; j < dof; ++j) {
                fExt += H[uOff + j] * mF[uOff + j];
            }
            SpatialVec r = (Mk[b] * A[b]) + gyro[b] - fExt;
            for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
                const int c = m.bodyChildren[ci];
                r += Phi[c] * ref[static_cast<std::size_t>(c)];
            }
            ref[static_cast<std::size_t>(b)] = r;
            EXPECT_TRUE(rtest::NearVec3(reacBo[b].angular, r.angular, rtest::kTight))
                << "ang body " << b << " rep " << rep;
            EXPECT_TRUE(rtest::NearVec3(reacBo[b].linear, r.linear, rtest::kTight))
                << "lin body " << b << " rep " << rep;

            // M-frame output must be the Bo reaction shifted by the rigid Bo->Mo
            // offset: [t;f]@Bo -> [t - (R_GB X_BM.p) x f ; f]@Mo. Pins the frame
            // convention independently of the recursion algebra.
            const Vec3 p_BoMo_G = X_GB[b].R() * m.X_BM[b].p();
            const SpatialVec rM(r.angular - (p_BoMo_G % r.linear), r.linear);
            EXPECT_TRUE(rtest::NearVec3(reacM[b].angular, rM.angular, rtest::kTight))
                << "M-shift ang body " << b << " rep " << rep;
            EXPECT_TRUE(rtest::NearVec3(reacM[b].linear, rM.linear, rtest::kTight))
                << "M-shift lin body " << b << " rep " << rep;
        }
    }
}

// ---------------------------------------------------------------------------
//  2. FreeMobilizerTransmitsNoConstraintForce: a Free (6-dof) root applies no
//     constraint force across its joint, so its reaction reported at M is zero
//     (only the mobility force, which is zero here, could appear). The classic
//     Simbody reaction sanity check. A Torsion root, by contrast, transmits a
//     clearly nonzero reaction (it constrains 5 of 6 spatial dof).
// ---------------------------------------------------------------------------
TEST(ReactionForces, FreeMobilizerTransmitsNoConstraintForce) {
    Rng rng(0x7A02);
    for (int rep = 0; rep < 10; ++rep) {
        RobotModel mFree = chain(rng, JointType::Free);
        RobotState s;
        s.allocateFull(mFree);
        randomizeState(mFree, s, rng);
        RobotEngine::realizePosition(mFree, s);
        AnalyticForceBridge bridge(mFree, s, 150.0);
        randomizeState(mFree, s, rng);
        realizeDynamics(mFree, s, bridge);

        const SpatialVec reacM = RobotEngine::findMobilizerReactionOnBodyAtMInGround(mFree, s, 1);
        EXPECT_LT(reacM.angular.norm(), rtest::kLoose) << "free-root reaction torque nonzero, rep " << rep;
        EXPECT_LT(reacM.linear.norm(), rtest::kLoose) << "free-root reaction force nonzero, rep " << rep;
    }

    // contrast: a Torsion root genuinely transmits a reaction.
    {
        RobotModel mTor = chain(rng, JointType::Torsion);
        RobotState s;
        s.allocateFull(mTor);
        randomizeState(mTor, s, rng);
        RobotEngine::realizePosition(mTor, s);
        AnalyticForceBridge bridge(mTor, s, 150.0);
        randomizeState(mTor, s, rng);
        realizeDynamics(mTor, s, bridge);

        const SpatialVec reacM = RobotEngine::findMobilizerReactionOnBodyAtMInGround(mTor, s, 1);
        const Real mag = reacM.angular.norm() + reacM.linear.norm();
        EXPECT_GT(mag, 1.0) << "a 1-dof Torsion root should transmit a real reaction";
    }
}

// ---------------------------------------------------------------------------
//  3. WholeTreeNewtonEulerBalance: the only inboard connection of the robot to
//     Ground is the root mobilizer, so the root reaction (at Bo, shifted to the
//     global origin) must equal the sum over all bodies of (Mk A + gyro - F_ext)
//     shifted to the same origin. A global statics/dynamics consistency law that
//     is INDEPENDENT of the per-body recursion's internal Phi transmission.
// ---------------------------------------------------------------------------
TEST(ReactionForces, WholeTreeNewtonEulerBalance) {
    Rng rng(0x7A03);
    RobotModel m = chain(rng, JointType::Free);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 10; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        AnalyticForceBridge bridge(m, s, 150.0);
        randomizeState(m, s, rng);
        realizeDynamics(m, s, bridge);

        std::vector<SpatialVec> reacBo(static_cast<std::size_t>(m.numBodies));
        RobotEngine::calcMobilizerReactionForces(m, s, reacBo.data(), nullptr);

        const SpatialInertia* Mk = s.Mk_G();
        const SpatialVec* A = s.A_GB();
        const SpatialVec* gyro = s.gyro();
        const SpatialVec* bF = s.bodyForceG();
        const Real* mF = s.mobilityForce();
        const SpatialVec* H = s.H();
        const Transform* X_GB = s.X_GB();

        SpatialVec sumNE(Vec3(0), Vec3(0));
        for (int b = 1; b < m.numBodies; ++b) {
            const int uOff = m.bodyUIndex[b], dof = m.bodyNU[b];
            SpatialVec fExt = bF[b];
            for (int j = 0; j < dof; ++j) {
                fExt += H[uOff + j] * mF[uOff + j];
            }
            const SpatialVec ne = (Mk[b] * A[b]) + gyro[b] - fExt; // at Bo
            sumNE += shiftToOrigin(ne, X_GB[b].p());
        }
        const SpatialVec rootReacO = shiftToOrigin(reacBo[1], X_GB[1].p());

        EXPECT_TRUE(rtest::NearVec3(sumNE.angular, rootReacO.angular, rtest::kLoose))
            << "torque balance, rep " << rep;
        EXPECT_TRUE(rtest::NearVec3(sumNE.linear, rootReacO.linear, rtest::kLoose))
            << "force balance, rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  4. WeldChildTransmitsFullLoad: a Rigid (0-dof) body welded to a MOVABLE parent
//     has no freedom of its own, so it must transmit the entire inertial-minus-
//     applied load of its subtree across the weld. (A Rigid body welded straight
//     to Ground has a fixed X_GB and cannot be displaced from its anchors, so it
//     would carry no load -- the weld must sit below a movable joint to be
//     exercised.) Its reaction at Bo equals Mk A + gyro - bodyForce (zero dof ->
//     no mobility force, no children), and is nonzero once the parent moves it.
// ---------------------------------------------------------------------------
TEST(ReactionForces, WeldChildTransmitsFullLoad) {
    Rng rng(0x7A04);
    auto mk = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        s.X_BM = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        s.mass = rng.uniform(0.8, 1.6);
        s.com_B = rng.vec3(-0.1, 0.1);
        s.inertia_B = UnitInertia(0.5, 0.5, 0.5);
        return s;
    };
    // Free root (body 1) -> Rigid welded child (body 2).
    RobotModel m = buildForest({mk(0, JointType::Free), mk(1, JointType::Rigid)});
    attachAtoms(m, 1, {Vec3(0.0, 0.0, 0.0), Vec3(0.09, -0.02, 0.05)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.0, 0.0), Vec3(0.0, 0.07, -0.02)}, {13.0, 14.0});
    const int weld = 2;

    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    AnalyticForceBridge bridge(m, s, 150.0);
    randomizeState(m, s, rng); // move the Free root -> displaces the welded child's atoms
    realizeDynamics(m, s, bridge);

    std::vector<SpatialVec> reacBo(static_cast<std::size_t>(m.numBodies));
    RobotEngine::calcMobilizerReactionForces(m, s, reacBo.data(), nullptr);

    // welded body: 0 dof (no mobility force), no children -> reaction at Bo is
    // exactly Mk A + gyro - bodyForce.
    const SpatialVec expected = (s.Mk_G()[weld] * s.A_GB()[weld]) + s.gyro()[weld] - s.bodyForceG()[weld];
    EXPECT_TRUE(rtest::NearVec3(reacBo[weld].angular, expected.angular, rtest::kTight));
    EXPECT_TRUE(rtest::NearVec3(reacBo[weld].linear, expected.linear, rtest::kTight));

    // and the weld genuinely carries a load (its parent displaced its atoms).
    EXPECT_GT(reacBo[weld].linear.norm() + reacBo[weld].angular.norm(), 1e-6)
        << "welded child under load should transmit a nonzero reaction";
}

// ---------------------------------------------------------------------------
//  6. MobilityForceEntersReaction: with a NONZERO applied mobility (generalized
//     joint) force, the reaction must include its spatial image sum_j H_j tau_j
//     in F_ext. We set body forces and mobility forces directly (no bridge),
//     run calcUDot, and check the operator against the full recursion -- which
//     now genuinely depends on the mobility-force term.
//     FAIL guard: dropping the H*mobilityForce mapping makes this go red.
// ---------------------------------------------------------------------------
TEST(ReactionForces, MobilityForceEntersReaction) {
    Rng rng(0x7A06);
    RobotModel m = chain(rng, JointType::Free);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 10; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);

        // applied forces set DIRECTLY: random per-body spatial forces + random
        // per-dof mobility forces (the bridge would zero the latter).
        SpatialVec* bF = s.bodyForceG();
        for (int b = 0; b < m.numBodies; ++b) {
            bF[b] = SpatialVec(rng.vec3(-2, 2), rng.vec3(-2, 2));
        }
        Real* mF = s.mobilityForce();
        for (int i = 0; i < m.nu; ++i) {
            mF[i] = rng.gaussian();
        }

        RobotEngine::calcUDot(m, s); // A_GB consistent with these forces

        std::vector<SpatialVec> reacBo(static_cast<std::size_t>(m.numBodies));
        RobotEngine::calcMobilizerReactionForces(m, s, reacBo.data(), nullptr);

        const SpatialInertia* Mk = s.Mk_G();
        const SpatialVec* A = s.A_GB();
        const SpatialVec* gyro = s.gyro();
        const SpatialVec* H = s.H();
        const PhiMatrix* Phi = s.Phi();

        // recursion WITH the mobility-force term; and a CONTROL that omits it.
        std::vector<SpatialVec> ref(static_cast<std::size_t>(m.numBodies), SpatialVec(Vec3(0), Vec3(0)));
        std::vector<SpatialVec> refNoMob(static_cast<std::size_t>(m.numBodies), SpatialVec(Vec3(0), Vec3(0)));
        bool mobMatters = false;
        for (int b = m.numBodies - 1; b >= 1; --b) {
            const int uOff = m.bodyUIndex[b], dof = m.bodyNU[b];
            SpatialVec fExt = bF[b], fExtNoMob = bF[b];
            for (int j = 0; j < dof; ++j) {
                fExt += H[uOff + j] * mF[uOff + j];
            }
            SpatialVec r = (Mk[b] * A[b]) + gyro[b] - fExt;
            SpatialVec rNo = (Mk[b] * A[b]) + gyro[b] - fExtNoMob;
            for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
                const int c = m.bodyChildren[ci];
                r += Phi[c] * ref[static_cast<std::size_t>(c)];
                rNo += Phi[c] * refNoMob[static_cast<std::size_t>(c)];
            }
            ref[static_cast<std::size_t>(b)] = r;
            refNoMob[static_cast<std::size_t>(b)] = rNo;
            if ((r.angular - rNo.angular).norm() + (r.linear - rNo.linear).norm() > 1e-6) {
                mobMatters = true;
            }

            EXPECT_TRUE(rtest::NearVec3(reacBo[b].angular, r.angular, rtest::kTight))
                << "ang b" << b << " rep " << rep;
            EXPECT_TRUE(rtest::NearVec3(reacBo[b].linear, r.linear, rtest::kTight))
                << "lin b" << b << " rep " << rep;
        }
        // the mobility term must actually MATTER here, else the test is vacuous.
        EXPECT_TRUE(mobMatters) << "mobility force did not affect the reaction -- term untested, rep " << rep;
    }
}
// ---------------------------------------------------------------------------
//  7. GroundHasNoReaction: body 0 (Ground) has no inboard mobilizer; its output
//     slot is left zero in both the Bo and M frames.
// ---------------------------------------------------------------------------
TEST(ReactionForces, GroundHasNoReaction) {
    Rng rng(0x7A05);
    RobotModel m = chain(rng, JointType::Free);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    AnalyticForceBridge bridge(m, s, 150.0);
    randomizeState(m, s, rng);
    realizeDynamics(m, s, bridge);

    std::vector<SpatialVec> reacBo(static_cast<std::size_t>(m.numBodies));
    std::vector<SpatialVec> reacM(static_cast<std::size_t>(m.numBodies));
    RobotEngine::calcMobilizerReactionForces(m, s, reacBo.data(), reacM.data());

    EXPECT_EQ(reacBo[0].angular.norm(), 0.0);
    EXPECT_EQ(reacBo[0].linear.norm(), 0.0);
    EXPECT_EQ(reacM[0].angular.norm(), 0.0);
    EXPECT_EQ(reacM[0].linear.norm(), 0.0);
}