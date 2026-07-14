// ============================================================================
//  TestCyclicBoltzmann.cpp -- the loop-closure Fixman term that World::calcFixman
//  admits "was never tested for ring Boltzmann correctness."
//
//  For a cyclic molecule the spanning tree carries one too many flexible torsions;
//  the ring is closed by a RATTLE distance constraint that removes one DOF, and the
//  flexible mass-metric determinant gains a factor:
//        |M_{N_f}| = |M_tree| / det(G M^-1 G^T),
//  so the Fixman potential gains -(1/2) RT ln det(G M^-1 G^T)
//  (ConstraintSet::calcConstraintLogDet). This file validates that previously
//  untested term two ways:
//
//   1. ORACLE (deterministic, rigorous). calcConstraintLogDet is cross-checked
//      against an INDEPENDENT assembly of ln det(G M^-1 G^T): build each constraint
//      row G_i with mapAtomForcesToGeneralizedForces (atomForce = vecAB at A,
//      -vecAB at B, the same convention the production assembly uses), form
//      A_ij = G_i . (M^-1 G_j) via multiplyByMInv, and take ln|det A| directly.
//      Done for a 1-loop and a 2-loop ring, so both the scalar and the coupled
//      multi-loop determinant paths are pinned.
//
//   2. MANIFOLD SMOKE. A constrained ring sampled with HMC (SHAKE position +
//      RATTLE velocity projection) stays on the loop-closure manifold C(q)=0 and
//      its dihedral explores -- the dynamics the Fixman term corrects actually runs.
//
//  NOTE / REPORTED LIMITATION. A full ring-Boltzmann brute-force reference (the
//  unconstrained-Cartesian marginal the Fixman term should reproduce) needs an
//  independent Cartesian sampler on the constraint manifold and is left for future
//  work; the oracle above pins the term that was actually missing, deterministically.
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
#include "TestHelpers.hpp"
#include "support/SamplingHarness.hpp"
#include "support/TestPhysConstants.hpp"

using namespace robo;
using rtest::attachAtoms;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::HmcDriver;
using rtest::Rng;
using rtest::phys::kT300;

namespace {

// A Free-rooted torsion chain whose first and last bodies carry atoms that we will
// tie together with loop-closure distance constraints to form a ring.
RobotModel torsionChain(int nTorsions, Rng& rng) {
    std::vector<BodySpec> specs;
    {
        BodySpec root;
        root.parent = 0;
        root.joint = JointType::Free;
        root.mass = rng.uniform(1.0, 1.5);
        root.com_B = rng.vec3(-0.05, 0.05);
        root.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
        specs.push_back(root);
    }
    for (int i = 0; i < nTorsions; ++i) {
        BodySpec b;
        b.parent = static_cast<int>(specs.size());
        b.joint = JointType::Torsion;
        b.X_PF = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        b.X_BM = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        b.mass = rng.uniform(0.9, 1.4);
        b.com_B = rng.vec3(-0.08, 0.08);
        b.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
        specs.push_back(b);
    }
    RobotModel m = buildForest(specs);
    for (int b = 1; b < m.numBodies; ++b) {
        attachAtoms(m, b, {Vec3(0, 0, 0), Vec3(0.10, 0.02, -0.03)}, {12.0, 1.0});
    }
    return m;
}

// Independent ln|det(G M^-1 G^T)| assembled from scratch, matching the production
// G convention (atomForce = vecAB at A, -vecAB at B). Requires realizePosition +
// realizeArticulatedBodyInertias current.
double oracleLogDetGMInvGt(const RobotModel& m, RobotState& s, const ConstraintSet& cs) {
    const int nc = cs.numConstraints();
    const int nu = m.nu;
    const Vec3* pos = s.atomPosG();

    std::vector<std::vector<Real>> G(static_cast<std::size_t>(nc), std::vector<Real>(nu, 0));
    std::vector<std::vector<Real>> MinvG(static_cast<std::size_t>(nc), std::vector<Real>(nu, 0));
    std::vector<Vec3> atomForce(static_cast<std::size_t>(m.numAtoms), Vec3(0));

    for (int c = 0; c < nc; ++c) {
        const DistanceConstraint& d = cs.distance[static_cast<std::size_t>(c)];
        const Vec3 vecAB = pos[d.atomA] - pos[d.atomB];
        for (Vec3& f : atomForce) {
            f = Vec3(0);
        }
        atomForce[d.atomA] = vecAB;
        atomForce[d.atomB] = Vec3(0) - vecAB;
        ConstraintSet::mapAtomForcesToGeneralizedForces(m, s, atomForce.data(), G[c].data());
        RobotEngine::multiplyByMInv(m, s, G[c].data(), MinvG[c].data());
    }

    // A_ij = G_i . (M^-1 G_j)
    std::vector<std::vector<double>> A(static_cast<std::size_t>(nc), std::vector<double>(nc, 0.0));
    for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nc; ++j) {
            double acc = 0.0;
            for (int k = 0; k < nu; ++k) {
                acc += static_cast<double>(G[i][k]) * static_cast<double>(MinvG[j][k]);
            }
            A[i][j] = acc;
        }
    }

    // ln|det A| via LU with partial pivoting (nc is tiny).
    double logdet = 0.0;
    for (int col = 0; col < nc; ++col) {
        int piv = col;
        for (int r = col + 1; r < nc; ++r) {
            if (std::abs(A[r][col]) > std::abs(A[piv][col])) {
                piv = r;
            }
        }
        std::swap(A[piv], A[col]);
        const double d = A[col][col];
        logdet += std::log(std::abs(d));
        for (int r = col + 1; r < nc; ++r) {
            const double f = A[r][col] / d;
            for (int cc = col; cc < nc; ++cc) {
                A[r][cc] -= f * A[col][cc];
            }
        }
    }
    return logdet;
}

// Close `pairs` loops: tie atomA(bodyA) to atomA(bodyB), restLength = current dist.
ConstraintSet closeLoops(const RobotModel& m, RobotState& s, const std::vector<std::pair<int, int>>& pairs) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    ConstraintSet cs;
    for (auto [ba, bb] : pairs) {
        const int a = m.bodyAtoms[m.bodyAtomsBeg[ba]];
        const int b = m.bodyAtoms[m.bodyAtomsBeg[bb]];
        const double d0 = (s.atomPosG()[a] - s.atomPosG()[b]).norm();
        cs.distance.push_back(DistanceConstraint{a, b, static_cast<Real>(d0)});
    }
    return cs;
}

} // namespace

// ---------------------------------------------------------------------------
//  ORACLE, single loop: calcConstraintLogDet == ln det(G M^-1 G^T) (a scalar here).
// ---------------------------------------------------------------------------
TEST(CyclicBoltzmann, LoopLogDetMatchesOracleOneLoop) {
    Rng rng(0xC1);
    RobotModel m = torsionChain(/*nTorsions*/ 4, rng);
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1; // unit quaternion on the Free root
    rtest::randomizeState(m, s, rng);
    s.q()[0] = std::abs(s.q()[0]) + Real(0.1); // keep a valid quaternion sign
    // tie body 1 to the last body to close one ring.
    ConstraintSet cs = closeLoops(m, s, {{1, m.numBodies - 1}});

    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    const double term = static_cast<double>(cs.calcConstraintLogDet(m, s));
    const double oracle = oracleLogDetGMInvGt(m, s, cs);
    EXPECT_NEAR(term, oracle, 1e-9) << "loop-closure ln det(G M^-1 G^T) disagrees with the oracle";
}

// ---------------------------------------------------------------------------
//  ORACLE, two loops: pins the coupled 2x2 determinant path.
// ---------------------------------------------------------------------------
TEST(CyclicBoltzmann, LoopLogDetMatchesOracleTwoLoops) {
    Rng rng(0xC2);
    RobotModel m = torsionChain(6, rng);
    RobotState s;
    s.allocateFull(m);
    rtest::randomizeState(m, s, rng);
    s.q()[0] = std::abs(s.q()[0]) + Real(0.1);
    ConstraintSet cs = closeLoops(m, s, {{1, 4}, {3, m.numBodies - 1}});

    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    const double term = static_cast<double>(cs.calcConstraintLogDet(m, s));
    const double oracle = oracleLogDetGMInvGt(m, s, cs);
    EXPECT_NEAR(term, oracle, 1e-9) << "two-loop ln det(G M^-1 G^T) disagrees with the oracle";
}

// ---------------------------------------------------------------------------
//  MANIFOLD SMOKE: a constrained ring sampled with the loop Fixman term ON stays
//  on the loop-closure manifold and explores; the loop term is non-zero (so the
//  test is not vacuous).
// ---------------------------------------------------------------------------
TEST(CyclicBoltzmann, ConstrainedRingStaysOnManifold) {
    Rng rng(0xC3);
    RobotModel m = torsionChain(4, rng);
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1;
    ConstraintSet cs = closeLoops(m, s, {{1, m.numBodies - 1}});

    // SHAKE the start onto C(q)=0.
    auto refresh = [&]() {
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
    };
    refresh();
    cs.enforcePositionConstraints(m, s, refresh);

    RobotEngine::realizeArticulatedBodyInertias(m, s);
    const double loopTerm = static_cast<double>(cs.calcConstraintLogDet(m, s));
    EXPECT_NE(loopTerm, 0.0) << "loop-closure term is zero -- ring test is vacuous";

    AnalyticForceBridge bridge(m, s, /*k*/ Real(0));
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, Real(0.005), 8, 0xC30D);
    drv.useFixman = true; // includes the loop-closure term

    const DistanceConstraint& d = cs.distance[0];
    double worst = 0.0;
    const rtest::Marginals marg = rtest::runHmcChain(drv, 3000, /*stride*/ 1, [&](long /*i*/, bool /*acc*/) {
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        const double dist = (s.atomPosG()[d.atomA] - s.atomPosG()[d.atomB]).norm();
        worst = std::max(worst, std::abs(dist - static_cast<double>(d.restLength)));
    });
    const long accepted = marg.accepted;
    EXPECT_GT(accepted, 0) << "no move accepted on the constrained ring";
    EXPECT_LT(worst, 1e-4) << "ring drifted off the loop-closure manifold (worst |dist-d0|=" << worst << ")";
}
