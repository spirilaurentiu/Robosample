// ============================================================================
//  TestConstraints.cpp -- Phase 3: the PUBLIC loop-closure paths of
//  robo::ConstraintSet on a real hand-built robot forest -- SHAKE
//  (enforcePositionConstraints), RATTLE (enforceVelocityConstraints), the
//  Jacobian-transpose force map (mapAtomForcesToGeneralizedForces), the
//  loop-closure Fixman term (calcConstraintLogDet), and manifold preservation
//  under the real templated verletStep + AnalyticForceBridge.
//
//  The private dense solvers (solveSmallSpd / solveCoupling) already have a
//  direct test in TestConstraintSolver.cpp via ConstraintTestAccess; here we
//  exercise the production-facing entry points end to end.
//
//  FIXTURE. A single ring is the macrocycle the header is written for:
//      Ground -> b1 (Free) -> b2 (Torsion) -> b3 (Torsion)
//  with atoms on each body and ONE DistanceConstraint joining an atom on b1 to
//  an atom on b3 (different bodies, so the bond's path runs through the shared
//  torsion u -- the loop the spanning tree carries one too many DOF for). d0 is
//  taken from a chosen reference configuration so C(q_ref) == 0 exactly, then
//  the state is opened so |C| >> tol and the projection has real work to do.
//
//  ORACLE. The dense coupling matrix A = G M^-1 G^T is assembled with the SAME
//  primitives the engine uses -- mapAtomForcesToGeneralizedForces on the
//  (+vecAB at A, -vecAB at B) unit pattern gives a row of G^T, multiplyByMInv
//  gives the matching M^-1 G^T column, and A[r][c] = (G^T_r) . (M^-1 G^T_c).
//  That is byte-for-byte how solveCoupling builds A, so ln det A from
//  robo_linalg::logDetSymPD is an independent check of calcConstraintLogDet.
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
#include "RobotIntegrator.hpp"
#include "RobotLinearAlgebra.hpp"
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

// ---- named tolerances local to the constraint paths ------------------------
constexpr Real kShakeTol = Real(1e-10); // SHAKE convergence target (== solver default)
constexpr Real kConstr = Real(1e-9);    // post-RATTLE bond-rate residual
constexpr Real kManifold = Real(1e-8);  // trajectory drift of |C| / |G u|

// The ring atoms inside the loop fixture. b1 has atoms {0,1}; b2 {2,3}; b3 {4,5}.
// The ring-closing bond joins atom A on b1 to atom B on b3.
constexpr int kAtomA = 0; // on body 1 (Free)
constexpr int kAtomB = 4; // on body 3 (Torsion)

// ---------------------------------------------------------------------------
//  Build Ground -> Free -> Torsion -> Torsion with two atoms per moving body.
//  No constraint is added here; the caller attaches the DistanceConstraint with
//  a d0 it reads off a reference configuration (so C(q_ref) == 0 exactly).
// ---------------------------------------------------------------------------
RobotModel buildRingChain(Rng& rng) {
    auto mk = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3(-0.3, 0.3));
        s.X_BM = Transform(rng.rotation(), rng.vec3(-0.3, 0.3));
        s.mass = rng.uniform(0.8, 1.6);
        s.com_B = rng.vec3(-0.1, 0.1);
        s.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
        return s;
    };
    RobotModel m =
        buildForest({mk(0, JointType::Free), mk(1, JointType::Torsion), mk(2, JointType::Torsion)});
    attachAtoms(m, 1, {Vec3(0.0, 0.0, 0.0), Vec3(0.11, -0.02, 0.04)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.10, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0.0, 0.07, -0.05), Vec3(0.04, 0.0, 0.05)}, {16.0, 1.0});
    return m;
}

// refresh closure the SHAKE solver needs: realizePosition refills X_GB/H/atomPosG
// (C is nonlinear in q), realizeArticulatedBodyInertias refills P/G/DI for the
// multiplyByMInv inside the projection.
void refreshFor(const RobotModel& m, RobotState& s) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
}

// C(q) = |r_AB|^2 - d0^2 for one distance bond, read off the current atom posns.
Real evalC(const RobotModel& m, const RobotState& s, const DistanceConstraint& bond) {
    (void)m;
    const Vec3* p = s.atomPosG();
    const Vec3 vecAB = p[bond.atomA] - p[bond.atomB];
    return dot(vecAB, vecAB) - bond.restLength * bond.restLength;
}

// G u = relative velocity along the bond = dot(vecAB, vA - vB) -- the exact rhs
// enforceVelocityConstraints projects to zero (== half d/dt C along u).
Real evalGu(const RobotModel& m, const RobotState& s, const DistanceConstraint& bond) {
    const Transform* X_GB = s.X_GB();
    const Vec3* posG = s.atomPosG();
    const SpatialVec* V_GB = s.V_GB();
    const Vec3 vecAB = posG[bond.atomA] - posG[bond.atomB];
    const int bA = m.atomBody[bond.atomA];
    const int bB = m.atomBody[bond.atomB];
    const Vec3 stA = posG[bond.atomA] - X_GB[bA].p();
    const Vec3 stB = posG[bond.atomB] - X_GB[bB].p();
    const Vec3 vA = V_GB[bA].linear + (V_GB[bA].angular % stA);
    const Vec3 vB = V_GB[bB].linear + (V_GB[bB].angular % stB);
    return dot(vecAB, vA - vB);
}

// Assemble one row of G^T (mapAtomForcesToGeneralizedForces on the +/-vecAB
// pattern) and its M^-1 G^T partner -- the engine's exact assembly. Fills both
// jacT and minvJacT (each length nu). Requires position + ABI realized.
void assembleConstraintRow(const RobotModel& m,
                           RobotState& s,
                           const DistanceConstraint& bond,
                           std::vector<Real>& jacT,
                           std::vector<Real>& minvJacT) {
    const Vec3* posG = s.atomPosG();
    const Vec3 vecAB = posG[bond.atomA] - posG[bond.atomB];
    std::vector<Vec3> atomForce(static_cast<std::size_t>(m.numAtoms), Vec3(0));
    atomForce[static_cast<std::size_t>(bond.atomA)] = vecAB;
    atomForce[static_cast<std::size_t>(bond.atomB)] = Vec3(0) - vecAB;
    jacT.assign(static_cast<std::size_t>(m.nu), 0);
    minvJacT.assign(static_cast<std::size_t>(m.nu), 0);
    ConstraintSet::mapAtomForcesToGeneralizedForces(m, s, atomForce.data(), jacT.data());
    RobotEngine::multiplyByMInv(m, s, jacT.data(), minvJacT.data());
}

// Dense A = G M^-1 G^T (numC x numC, row-major) assembled exactly as solveCoupling
// does. Requires position + ABI realized.
std::vector<Real> denseCoupling(const RobotModel& m, RobotState& s, const ConstraintSet& cset) {
    const int numC = cset.numConstraints();
    std::vector<std::vector<Real>> jacT(static_cast<std::size_t>(numC));
    std::vector<std::vector<Real>> minvJacT(static_cast<std::size_t>(numC));
    for (int c = 0; c < numC; ++c) {
        assembleConstraintRow(m,
                              s,
                              cset.distance[static_cast<std::size_t>(c)],
                              jacT[static_cast<std::size_t>(c)],
                              minvJacT[static_cast<std::size_t>(c)]);
    }
    std::vector<Real> A(static_cast<std::size_t>(numC * numC), 0);
    for (int r = 0; r < numC; ++r) {
        for (int c = 0; c < numC; ++c) {
            Real acc = 0;
            for (int i = 0; i < m.nu; ++i) {
                acc += jacT[static_cast<std::size_t>(r)][static_cast<std::size_t>(i)]
                       * minvJacT[static_cast<std::size_t>(c)][static_cast<std::size_t>(i)];
            }
            A[static_cast<std::size_t>(r * numC + c)] = acc;
        }
    }
    return A;
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

} // namespace

// ===========================================================================
//  1. MapAtomForcesEqualsJacobianTranspose
//     mapAtomForcesToGeneralizedForces IS J^T. With atomForce = (+vecAB at A,
//     -vecAB at B), tau = (J_A - J_B)^T vecAB = (1/2) dC/dq. We verify it as a
//     central finite difference of C taken in u-space (perturb each u direction
//     by applyVelSpaceIncrementToQ, refresh atoms, recompute C). MUST-FAIL: the
//     flipped A/B sign breaks it -- pinning the "verify once against a finite-
//     difference of C" comment in the code.
// ===========================================================================
TEST(Constraints, MapAtomForcesEqualsJacobianTranspose) {
    Rng rng(0xC03A);
    RobotModel m = buildRingChain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    refreshFor(m, s);

    DistanceConstraint bond;
    bond.atomA = kAtomA;
    bond.atomB = kAtomB;
    bond.restLength = 0; // C = |rAB|^2 here; d0 irrelevant to dC/dq.

    // tau = J^T (vecAB at A, -vecAB at B): the engine's assembly.
    std::vector<Real> jacT, minvJacT;
    assembleConstraintRow(m, s, bond, jacT, minvJacT);

    // Central FD of C along each u-direction. applyVelSpaceIncrementToQ maps a
    // u-space increment onto q (quaternion DOF go through N(q)), so this is the
    // directional derivative dC/du_i, which must equal 2*tau[i].
    const Real h = Real(1e-6);
    for (int i = 0; i < m.nu; ++i) {
        std::vector<Real> qSave(s.q(), s.q() + m.nq);

        std::vector<Real> dvel(static_cast<std::size_t>(m.nu), 0);
        dvel[static_cast<std::size_t>(i)] = h;
        ConstraintSet::applyVelSpaceIncrementToQ(m, s, dvel.data());
        RobotEngine::realizePosition(m, s);
        const Real cPlus = evalC(m, s, bond);

        std::copy(qSave.begin(), qSave.end(), s.q());
        dvel[static_cast<std::size_t>(i)] = -h;
        ConstraintSet::applyVelSpaceIncrementToQ(m, s, dvel.data());
        RobotEngine::realizePosition(m, s);
        const Real cMinus = evalC(m, s, bond);

        std::copy(qSave.begin(), qSave.end(), s.q());
        RobotEngine::realizePosition(m, s);

        const Real dCdu = (cPlus - cMinus) / (2 * h);
        EXPECT_NEAR(jacT[static_cast<std::size_t>(i)], Real(0.5) * dCdu, rtest::kFD) << "u-direction " << i;
    }

    // MUST-FAIL guard: flip the A/B sign of the atom force. tau negates, so it
    // can no longer equal +1/2 dC/du -- at least one component must disagree.
    std::vector<Vec3> flipped(static_cast<std::size_t>(m.numAtoms), Vec3(0));
    const Vec3* p = s.atomPosG();
    const Vec3 vecAB = p[bond.atomA] - p[bond.atomB];
    flipped[static_cast<std::size_t>(bond.atomA)] = Vec3(0) - vecAB; // wrong sign
    flipped[static_cast<std::size_t>(bond.atomB)] = vecAB;           // wrong sign
    std::vector<Real> tauFlipped(static_cast<std::size_t>(m.nu), 0);
    ConstraintSet::mapAtomForcesToGeneralizedForces(m, s, flipped.data(), tauFlipped.data());

    bool anyDisagree = false;
    const Real h2 = Real(1e-6);
    for (int i = 0; i < m.nu; ++i) {
        std::vector<Real> qSave(s.q(), s.q() + m.nq);
        std::vector<Real> dvel(static_cast<std::size_t>(m.nu), 0);
        dvel[static_cast<std::size_t>(i)] = h2;
        ConstraintSet::applyVelSpaceIncrementToQ(m, s, dvel.data());
        RobotEngine::realizePosition(m, s);
        const Real cPlus = evalC(m, s, bond);
        std::copy(qSave.begin(), qSave.end(), s.q());
        dvel[static_cast<std::size_t>(i)] = -h2;
        ConstraintSet::applyVelSpaceIncrementToQ(m, s, dvel.data());
        RobotEngine::realizePosition(m, s);
        const Real cMinus = evalC(m, s, bond);
        std::copy(qSave.begin(), qSave.end(), s.q());
        RobotEngine::realizePosition(m, s);
        const Real halfDC = Real(0.5) * (cPlus - cMinus) / (2 * h2);
        if (std::abs(tauFlipped[static_cast<std::size_t>(i)] - halfDC) > 1e-5) {
            anyDisagree = true;
        }
    }
    EXPECT_TRUE(anyDisagree)
        << "flipped-sign atom force still matched 1/2 dC/du -- the A/B convention is not pinned";
}

// ===========================================================================
//  2. ShakeDrivesPositionConstraintToZero
//     Open the ring (|C| >> tol), SHAKE it shut. PASS: returns iter < maxIter and
//     |C| <= tol. PASS: a body NOT in the loop is unmoved beyond the shared-u
//     coupling (forest locality). MUST-FAIL: maxIter == 0 leaves |C| unchanged.
// ===========================================================================
TEST(Constraints, ShakeDrivesPositionConstraintToZero) {
    Rng rng(0xC03B);
    // The ring chain plus an independent second robot (a lone Free body with
    // atoms) that shares NO u with the loop -- forest locality must leave it put.
    auto mk = [&](int parent, JointType jt) {
        BodySpec sp;
        sp.parent = parent;
        sp.joint = jt;
        sp.X_PF = Transform(rng.rotation(), rng.vec3(-0.3, 0.3));
        sp.X_BM = Transform(rng.rotation(), rng.vec3(-0.3, 0.3));
        sp.mass = rng.uniform(0.8, 1.6);
        sp.com_B = rng.vec3(-0.1, 0.1);
        sp.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
        return sp;
    };
    RobotModel m = buildForest({
        mk(0, JointType::Free),    // 1  loop body
        mk(1, JointType::Torsion), // 2  loop body
        mk(2, JointType::Torsion), // 3  loop body
        mk(0, JointType::Free),    // 4  SECOND robot (not in loop)
    });
    attachAtoms(m, 1, {Vec3(0.0, 0.0, 0.0), Vec3(0.11, -0.02, 0.04)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.10, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0.0, 0.07, -0.05), Vec3(0.04, 0.0, 0.05)}, {16.0, 1.0});
    const int farAtomBeg = m.numAtoms;
    attachAtoms(m, 4, {Vec3(0.02, 0.0, 0.0), Vec3(0.0, 0.08, 0.0)}, {15.0, 1.0});
    const int farAtomEnd = m.numAtoms;

    RobotState s;
    s.allocateFull(m);

    // Reference config: set d0 to the actual current A-B distance so C == 0,
    // then OPEN the ring by re-randomizing.
    randomizeState(m, s, rng);
    refreshFor(m, s);
    const Vec3* p0 = s.atomPosG();
    const Real d0 = (p0[kAtomA] - p0[kAtomB]).norm();

    ConstraintSet cset;
    cset.distance.push_back(DistanceConstraint{kAtomA, kAtomB, d0});

    // Open the loop: a fresh random state breaks |rAB| away from d0.
    randomizeState(m, s, rng);
    refreshFor(m, s);
    const Real cOpen = evalC(m, s, cset.distance[0]);
    ASSERT_GT(std::abs(cOpen), 1e-3) << "setup failed to open the ring; reseed";

    // snapshot the far robot's atoms BEFORE projection
    std::vector<Vec3> farBefore;
    for (int a = farAtomBeg; a < farAtomEnd; ++a) {
        farBefore.push_back(s.atomPosG()[a]);
    }

    // MUST-FAIL guard: maxIter == 0 does no projection -> violation unchanged.
    {
        RobotState sGuard;
        sGuard.allocateFull(m);
        std::copy(s.q(), s.q() + m.nq, sGuard.q());
        std::copy(s.u(), s.u() + m.nu, sGuard.u());
        refreshFor(m, sGuard);
        const Real cBefore = evalC(m, sGuard, cset.distance[0]);
        const int it0 = cset.enforcePositionConstraints(
            m,
            sGuard,
            [&]() {
                refreshFor(m, sGuard);
            },
            kShakeTol,
            /*maxIter=*/0);
        const Real cAfter = evalC(m, sGuard, cset.distance[0]);
        EXPECT_EQ(it0, 0);
        EXPECT_NEAR(cBefore, cAfter, 1e-14)
            << "maxIter=0 changed C -- setup, not projection, is satisfying it";
    }

    // Real SHAKE: converges in < maxIter, drives |C| to <= tol.
    const int iter = cset.enforcePositionConstraints(
        m,
        s,
        [&]() {
            refreshFor(m, s);
        },
        kShakeTol,
        50);
    refreshFor(m, s);
    EXPECT_LT(iter, 50) << "SHAKE did not converge within the iteration cap";
    EXPECT_LE(std::abs(evalC(m, s, cset.distance[0])), kShakeTol) << "ring still open after SHAKE";

    // Forest locality: the second robot shares no u with the loop, so M^-1 G^T
    // is identically zero on its block and its atoms cannot move.
    for (int k = 0; k < farAtomEnd - farAtomBeg; ++k) {
        EXPECT_TRUE(rtest::NearVec3(s.atomPosG()[farAtomBeg + k],
                                    farBefore[static_cast<std::size_t>(k)],
                                    rtest::kTight))
            << "off-loop atom " << (farAtomBeg + k) << " moved -- forest locality broken";
    }
}

// ===========================================================================
//  3. RattleProjectsVelocityOntoConstraintSurface
//     After realizeVelocity, G u (the bond-length rate) is nonzero; RATTLE kills
//     it. PASS: |G u| <= tol. PASS: idempotent (second call is a no-op on u).
//     PASS: only the along-bond (minvJacobianT) component is removed; the
//     orthogonal complement of u is preserved.
// ===========================================================================
TEST(Constraints, RattleProjectsVelocityOntoConstraintSurface) {
    Rng rng(0xC03C);
    RobotModel m = buildRingChain(rng);
    RobotState s;
    s.allocateFull(m);

    randomizeState(m, s, rng);
    refreshFor(m, s);
    const Vec3* p0 = s.atomPosG();
    const Real d0 = (p0[kAtomA] - p0[kAtomB]).norm();
    ConstraintSet cset;
    cset.distance.push_back(DistanceConstraint{kAtomA, kAtomB, d0});

    // generic velocity -> nonzero bond-length rate
    randomizeState(m, s, rng);
    refreshFor(m, s);
    RobotEngine::realizeVelocity(m, s);
    const Real guBefore = evalGu(m, s, cset.distance[0]);
    ASSERT_GT(std::abs(guBefore), 1e-3) << "velocity already on the surface; reseed";

    // the along-bond direction in u-space: M^-1 G^T (the only direction RATTLE
    // is allowed to touch). Captured BEFORE the projection.
    std::vector<Real> jacT, minvJacT;
    assembleConstraintRow(m, s, cset.distance[0], jacT, minvJacT);
    std::vector<Real> uBefore(s.u(), s.u() + m.nu);

    // RATTLE
    cset.enforceVelocityConstraints(m, s);
    RobotEngine::realizeVelocity(m, s);
    EXPECT_LE(std::abs(evalGu(m, s, cset.distance[0])), kConstr) << "bond-length rate not killed";

    // Idempotent: a second projection leaves u unchanged.
    std::vector<Real> uAfter1(s.u(), s.u() + m.nu);
    cset.enforceVelocityConstraints(m, s);
    RobotEngine::realizeVelocity(m, s);
    for (int i = 0; i < m.nu; ++i) {
        EXPECT_NEAR(s.u()[i], uAfter1[static_cast<std::size_t>(i)], kConstr)
            << "RATTLE not idempotent at u[" << i << "]";
    }

    // Energy-only: the update du = u_after - u_before is a pure multiple of
    // minvJacobianT (a single ring -> a single direction). Any component of u
    // orthogonal to minvJacobianT is therefore preserved exactly.
    Real mm = 0;
    for (int i = 0; i < m.nu; ++i) {
        mm += minvJacT[static_cast<std::size_t>(i)] * minvJacT[static_cast<std::size_t>(i)];
    }
    ASSERT_GT(mm, 0.0);
    // du must be parallel to minvJacT: du - (du.mhat) mhat == 0.
    Real duDotM = 0;
    for (int i = 0; i < m.nu; ++i) {
        duDotM += (uAfter1[static_cast<std::size_t>(i)] - uBefore[static_cast<std::size_t>(i)])
                  * minvJacT[static_cast<std::size_t>(i)];
    }
    const Real coeff = duDotM / mm;
    for (int i = 0; i < m.nu; ++i) {
        const Real du = uAfter1[static_cast<std::size_t>(i)] - uBefore[static_cast<std::size_t>(i)];
        EXPECT_NEAR(du, coeff * minvJacT[static_cast<std::size_t>(i)], rtest::kLoose)
            << "RATTLE moved u off the minvJacobianT direction at u[" << i << "]";
    }

    // Concretely: pick a test direction orthogonal to minvJacT and confirm its
    // projection of u is unchanged to kLoose.
    std::vector<Real> e(static_cast<std::size_t>(m.nu), 0);
    e[0] = 1.0; // arbitrary probe
    Real eDotM = 0;
    for (int i = 0; i < m.nu; ++i) {
        eDotM += e[static_cast<std::size_t>(i)] * minvJacT[static_cast<std::size_t>(i)];
    }
    for (int i = 0; i < m.nu; ++i) {
        e[static_cast<std::size_t>(i)] -=
            (eDotM / mm) * minvJacT[static_cast<std::size_t>(i)]; // make e _|_ minvJacT
    }
    Real projBefore = 0, projAfter = 0;
    for (int i = 0; i < m.nu; ++i) {
        projBefore += e[static_cast<std::size_t>(i)] * uBefore[static_cast<std::size_t>(i)];
        projAfter += e[static_cast<std::size_t>(i)] * uAfter1[static_cast<std::size_t>(i)];
    }
    EXPECT_NEAR(projBefore, projAfter, rtest::kLoose)
        << "u component orthogonal to the bond direction was not preserved";
}

// ===========================================================================
//  4. ConstraintManifoldPreservedUnderIntegration
//     Integrate the ring with the real templated verletStep + AnalyticForceBridge
//     + the ConstraintSet. PASS: max_t |C| and max_t |G u| stay <= 1e-8 across
//     2000 steps (the secular drift the header warns about never opens the ring).
//     PASS (O(h^2)): halving h shrinks the energy drift by ~4x.
// ===========================================================================
TEST(Constraints, ConstraintManifoldPreservedUnderIntegration) {
    Rng rng(0xC03D);
    RobotModel m = buildRingChain(rng);
    RobotState s;
    s.allocateFull(m);

    // Build a closed, gently-moving start state and freeze the ring at d0.
    randomizeState(m, s, rng);
    Real* u0 = s.u();
    for (int i = 0; i < m.nu; ++i) {
        u0[i] *= 0.05; // gentle velocities so the step is benign
    }
    refreshFor(m, s);
    const Vec3* p0 = s.atomPosG();
    const Real d0 = (p0[kAtomA] - p0[kAtomB]).norm();

    ConstraintSet cset;
    cset.distance.push_back(DistanceConstraint{kAtomA, kAtomB, d0});

    AnalyticForceBridge bridge(m, s, /*k=*/80.0);

    // Project the start state ONTO the manifold first (SHAKE then RATTLE), so the
    // trajectory begins exactly on C = 0, G u = 0.
    cset.enforcePositionConstraints(
        m,
        s,
        [&]() {
            refreshFor(m, s);
        },
        kShakeTol,
        50);
    refreshFor(m, s);
    RobotEngine::realizeVelocity(m, s);
    cset.enforceVelocityConstraints(m, s);
    seedDerivatives(m, s, bridge);

    const Real h = Real(5e-4);
    Real maxC = 0, maxGu = 0;
    for (int n = 0; n < 2000; ++n) {
        const bool ok = RobotEngine::verletStep(m, s, bridge, cset, h);
        ASSERT_TRUE(ok) << "step " << n << " rejected";
        RobotEngine::realizeVelocity(m, s);
        maxC = std::max(maxC, std::abs(evalC(m, s, cset.distance[0])));
        maxGu = std::max(maxGu, std::abs(evalGu(m, s, cset.distance[0])));
    }
    EXPECT_LE(maxC, kManifold) << "position constraint drifted: max|C| = " << maxC;
    EXPECT_LE(maxGu, kManifold) << "velocity constraint drifted: max|G u| = " << maxGu;

    // O(h^2): integrate the SAME closed start at h and h/2 over the SAME physical
    // time; the RATTLE energy drift is quadratic (Andersen Appendix B), so the
    // ratio drift(h)/drift(h/2) lands near 4. The two runs MUST share an identical
    // start configuration -- otherwise the ratio compares two unrelated drifts and
    // is only statistically ~4. A fixed local seed gives both runs the same q0,u0.
    auto driftAt = [&](Real step, int nSteps) -> Real {
        Rng seedRng(0x0F1E2D3C); // identical start state for every call
        RobotState st;
        st.allocateFull(m);
        randomizeState(m, st, seedRng);
        Real* uu = st.u();
        for (int i = 0; i < m.nu; ++i) {
            uu[i] *= 0.05;
        }
        refreshFor(m, st);
        const Vec3* pp = st.atomPosG();
        const Real dd = (pp[kAtomA] - pp[kAtomB]).norm();
        ConstraintSet cs;
        cs.distance.push_back(DistanceConstraint{kAtomA, kAtomB, dd});
        AnalyticForceBridge br(m, st, 80.0);
        cs.enforcePositionConstraints(
            m,
            st,
            [&]() {
                refreshFor(m, st);
            },
            kShakeTol,
            50);
        refreshFor(m, st);
        RobotEngine::realizeVelocity(m, st);
        cs.enforceVelocityConstraints(m, st);
        seedDerivatives(m, st, br);
        const Real E0 = RobotEngine::calcKineticEnergy(m, st) + br.calcPotentialEnergy(st);
        Real worst = 0;
        for (int n = 0; n < nSteps; ++n) {
            const bool ok = RobotEngine::verletStep(m, st, br, cs, step);
            if (!ok) {
                return std::numeric_limits<Real>::infinity();
            }
            const Real E = RobotEngine::calcKineticEnergy(m, st) + br.calcPotentialEnergy(st);
            worst = std::max(worst, std::abs(E - E0));
        }
        return worst;
    };
    // same total time: 600 steps at h, 1200 steps at h/2.
    const Real driftH = driftAt(h, 600);
    const Real driftHalf = driftAt(h / 2, 1200);
    ASSERT_GT(driftHalf, 0.0);
    const Real ratio = driftH / driftHalf;
    EXPECT_GE(ratio, 3.0) << "energy drift not shrinking ~4x under h/2 (ratio " << ratio << ")";
    EXPECT_LE(ratio, 5.0) << "energy drift shrinking too fast for O(h^2) (ratio " << ratio << ")";
}

// ===========================================================================
//  5. ConstraintLogDetMatchesDenseLnDet
//     calcConstraintLogDet returns ln det(G M^-1 G^T) (up to the documented
//     constant that cancels in any dH). Compare against the dense A built from
//     the same assembly via robo_linalg::logDetSymPD. PASS: single ring (1x1)
//     and a two-ring (2x2) case agree to 1e-7. PASS: acyclic -> exactly 0.0.
//     MUST-FAIL: the cyclic return must be != 0 (catches a regression that
//     early-returns 0 for cyclic systems -- the bug the header was written for).
// ===========================================================================
TEST(Constraints, ConstraintLogDetMatchesDenseLnDet) {
    Rng rng(0xC03E);
    RobotModel m = buildRingChain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);

    RobotEngine::realizePosition(m, s);
    const Vec3* p0 = s.atomPosG();
    const Real d0 = (p0[kAtomA] - p0[kAtomB]).norm();

    // ---- single ring: A is 1x1, ln det A == ln(A00) ----
    {
        ConstraintSet cset;
        cset.distance.push_back(DistanceConstraint{kAtomA, kAtomB, d0});

        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        const std::vector<Real> A = denseCoupling(m, s, cset); // 1x1
        const Real lnDetDense = robo_linalg::logDetSymPD(A.data(), 1);

        const Real lnDet = cset.calcConstraintLogDet(m, s);
        EXPECT_NEAR(lnDet, lnDetDense, 1e-7) << "single-ring ln det mismatch";

        // MUST-FAIL guard: a cyclic system must NOT return 0 (the early-return bug).
        EXPECT_GT(std::abs(lnDet), 1e-12)
            << "cyclic calcConstraintLogDet returned ~0 -- loop-closure Fixman term lost";
    }

    // ---- two rings: A is 2x2, exercises the dense logdet path ----
    {
        // A second independent ring-closing bond on the SAME chain: join atom 1
        // (on b1) to atom 5 (on b3). Distinct from the first bond, so A is a
        // genuine 2x2 (the two constraints share u, off-diagonal nonzero).
        RobotEngine::realizePosition(m, s);
        const Vec3* p = s.atomPosG();
        const Real d1 = (p[1] - p[5]).norm();

        ConstraintSet cset2;
        cset2.distance.push_back(DistanceConstraint{kAtomA, kAtomB, d0});
        cset2.distance.push_back(DistanceConstraint{1, 5, d1});

        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        const std::vector<Real> A = denseCoupling(m, s, cset2); // 2x2
        const Real lnDetDense = robo_linalg::logDetSymPD(A.data(), 2);

        const Real lnDet = cset2.calcConstraintLogDet(m, s);
        EXPECT_NEAR(lnDet, lnDetDense, 1e-7) << "two-ring ln det mismatch";
    }

    // ---- acyclic: empty distance set returns EXACTLY 0.0 ----
    {
        const ConstraintSet empty;
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        const Real lnDet = empty.calcConstraintLogDet(m, s);
        EXPECT_EQ(lnDet, 0.0) << "acyclic calcConstraintLogDet must be an exact no-op";
    }
}