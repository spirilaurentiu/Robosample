// ============================================================================
//  TestTwoRobotContact.cpp -- the correctness gate for docs/specs/
//  two-robot-contact/, PRIMARILY 10-mixed-integrator-correctness.md and
//  50-validation.md Sec.2 (INV-WELD, INV-FIX, INV-DRAW, INV-REV, INV-KE).
//
//  SCOPE. Two disjoint construction styles, matched to what each invariant
//  actually exercises:
//
//   * INV-WELD needs the REAL World::setCartesianSolvent guard (the one new
//     code path the spec adds), so it drives a genuine World through
//     World::buildModel(SystemTopology, Selection, rootMobilities) -- which,
//     like World::reinitialize, needs no OpenMM force evaluation, only the
//     topology/mobility inputs (mirrors TestAlchemy.cpp's hand-built
//     SystemTopology, no prmtop file). setCartesianSolvent/calcFixman/
//     drawSolventVelocities/calcSolventKE are private on World, so INV-WELD
//     is the only test here that goes through World at all.
//
//   * INV-FIX/INV-DRAW/INV-REV/INV-KE are mechanics claims (M1-M4 of
//     10-mixed-integrator-correctness.md) about RobotEngine + RobotState, so
//     -- exactly TestNcmcExplicitSolvent.cpp's established convention --
//     they drive RobotModel/RobotState/RobotEngine directly against a small
//     hand-built two-robot fixture (rtest::buildForest/attachAtoms) and an
//     analytic force bridge, REPRODUCING World::calcFixman/drawSolvent-
//     Velocities/calcSolventKE's documented formulas (World.hpp:408-446,
//     World.cpp:851-880) rather than reaching into World's private surface.
//
//  FIXTURE. R = TestFixmanBoltzmann.cpp's validated two-torsion 90-degree-
//  bend chain (duplicated here per the established per-file convention --
//  see TestNcmcExplicitSolvent.cpp's own header comment), atoms 0-4. E = one
//  ADDITIONAL 0-DOF Rigid body welded straight to Ground (a separate branch
//  of the SAME forest, atom-disjoint from R), atoms 5-6 -- the two-robot
//  contact world of docs/specs/two-robot-contact/00-diagnosis-and-
//  reformulation.md Sec.2.
// ============================================================================
#include <cmath>
#include <cstdint>
#include <gtest/gtest.h>
#include <random>
#include <stdexcept>
#include <vector>

#include "Constraints.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotIntegrator.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "StatTest.hpp"
#include "TestHelpers.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"
#include "support/HarmonicBridge.hpp"
#include "support/TestPhysConstants.hpp"

using namespace robo;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::attachAtoms;
using rtest::HarmonicBridge;
using rtest::SingleWellPolicy;
using rtest::phys::kT300;
using rtest::stat::MeanAccumulator;

using TwoRobotBridge = HarmonicBridge<SingleWellPolicy>;

namespace {

// ---------------------------------------------------------------------------
//  R (twoTorsionChain, TestFixmanBoltzmann.cpp/TestNcmcExplicitSolvent.cpp's
//  validated fixture) + E (one 0-DOF Rigid body welded to Ground, its own
//  branch of the same forest). Atoms 0,1,2 -> R body 1; atoms 3,4 -> R body
//  2; atoms 5,6 -> E's welded body. nu == 2 (R's two torsions only -- E
//  contributes 0, per P1).
// ---------------------------------------------------------------------------
RobotModel twoRobotFixture() {
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

    // E: a SEPARATE branch off Ground (parent = 0), Rigid == Weld, 0 DOF --
    // the P1 construction rule (docs/specs/two-robot-contact/
    // 10-mixed-integrator-correctness.md Sec.2). Placed away from R so the
    // two robots do not spatially overlap in this synthetic fixture.
    BodySpec e;
    e.parent = 0;
    e.joint = JointType::Rigid;
    e.X_PF = Transform(Vec3(Real(1.0), Real(0.0), Real(0.0)));
    e.X_BM = Transform();
    e.mass = Real(17.0);
    e.com_B = Vec3(Real(0.02), Real(0.0), Real(0.0));
    e.inertia_B = UnitInertia(Real(0.3), Real(0.3), Real(0.3));

    RobotModel m = buildForest({b1, b2, e});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, 0.02, -0.03), Vec3(-0.02, 0.10, 0.01)}, {12.0, 1.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0.0, 0.0, 0.0), Vec3(0.06, -0.02, 0.01)}, {16.0, 1.0});
    return m;
}

// E's atom indices, read back from the fixture's body-atom CSR (not
// hardcoded) so a fixture edit cannot silently desync the tests below.
std::vector<int> eAtoms(const RobotModel& m) {
    const int eBody = 3;
    std::vector<int> atoms;
    for (int ci = m.bodyAtomsBeg[eBody]; ci < m.bodyAtomsEnd[eBody]; ++ci) {
        atoms.push_back(m.bodyAtoms[ci]);
    }
    return atoms;
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

// ---------------------------------------------------------------------------
//  Small dense Gauss-Jordan inverse (partial pivoting), n small (here n ==
//  model.nu == 2). Used ONLY to reconstruct M_phi = (L L^T)^-1 from the SAME
//  multiplyBySqrtMInv operator the momentum draw uses, independently of
//  calcKineticEnergy's own articulated-recursion formula (INV-KE).
// ---------------------------------------------------------------------------
void invertSmall(const std::vector<double>& Ain, int n, std::vector<double>& Aout) {
    std::vector<double> A = Ain;
    Aout.assign(static_cast<std::size_t>(n * n), 0.0);
    for (int i = 0; i < n; ++i) {
        Aout[static_cast<std::size_t>(i * n + i)] = 1.0;
    }
    for (int col = 0; col < n; ++col) {
        int piv = col;
        double best = std::abs(A[static_cast<std::size_t>(col * n + col)]);
        for (int r = col + 1; r < n; ++r) {
            const double v = std::abs(A[static_cast<std::size_t>(r * n + col)]);
            if (v > best) {
                best = v;
                piv = r;
            }
        }
        if (piv != col) {
            for (int k = 0; k < n; ++k) {
                std::swap(A[static_cast<std::size_t>(col * n + k)], A[static_cast<std::size_t>(piv * n + k)]);
                std::swap(Aout[static_cast<std::size_t>(col * n + k)], Aout[static_cast<std::size_t>(piv * n + k)]);
            }
        }
        const double d = A[static_cast<std::size_t>(col * n + col)];
        for (int k = 0; k < n; ++k) {
            A[static_cast<std::size_t>(col * n + k)] /= d;
            Aout[static_cast<std::size_t>(col * n + k)] /= d;
        }
        for (int r = 0; r < n; ++r) {
            if (r == col) {
                continue;
            }
            const double f = A[static_cast<std::size_t>(r * n + col)];
            if (f == 0.0) {
                continue;
            }
            for (int k = 0; k < n; ++k) {
                A[static_cast<std::size_t>(r * n + k)] -= f * A[static_cast<std::size_t>(col * n + k)];
                Aout[static_cast<std::size_t>(r * n + k)] -= f * Aout[static_cast<std::size_t>(col * n + k)];
            }
        }
    }
}

// L = sqrt(M_phi^-1) as an explicit nu x nu row-major matrix: column j is
// multiplyBySqrtMInv(e_j) -- the SAME operator World::reinitialize's momentum
// draw calls (World.cpp:1686). PRECONDITION: realizeArticulatedBodyInertias
// already current for the state's q.
std::vector<double> sqrtMInvMatrix(const RobotModel& m, RobotState& s) {
    const int nu = m.nu;
    std::vector<double> L(static_cast<std::size_t>(nu * nu), 0.0);
    for (int j = 0; j < nu; ++j) {
        std::vector<Real> e(static_cast<std::size_t>(nu), Real(0)), col(static_cast<std::size_t>(nu));
        e[static_cast<std::size_t>(j)] = Real(1);
        RobotEngine::multiplyBySqrtMInv(m, s, e.data(), col.data());
        for (int i = 0; i < nu; ++i) {
            L[static_cast<std::size_t>(i * nu + j)] = static_cast<double>(col[static_cast<std::size_t>(i)]);
        }
    }
    return L;
}

// TwoRobotBridge (a type alias, above) is a lambda-free analytic bridge:
// independent per-atom harmonic anchors (AnalyticForceBridge.hpp's pattern),
// EXTENDED (HarmonicBridge<SingleWellPolicy> constructed with cacheAtomForces
// = true, below) to also cache the raw per-atom Cartesian force into
// atomForceG() when state.wantsAtomForces() -- the Cartesian-solvent Verlet
// block (RobotIntegrator.hpp:264-276, 419-429) reads that cache for E's
// atoms, exactly matching ForceBridge::getForcesFromOpenMM's real contract
// (ForceBridge.hpp:77-108).

void seedDerivatives(const RobotModel& m, RobotState& s, TwoRobotBridge& bridge) {
    RobotEngine::realizePosition(m, s);
    bridge.evaluate(s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    RobotEngine::calcUDot(m, s);
    RobotEngine::calcQDot(m, s, s.qdot());
    RobotEngine::calcQDotDot(m, s);
}

// Mirrors World::calcFixman()'s formula EXACTLY (World.cpp:851-880):
//   U_F = 1/2 RT (lnDetM - lnDetZ - lnDetMCartesian_), lnDetZ == 0 (acyclic).
// `cartesianRefAtomCount` is the knob INV-FIX turns: how many LEADING atoms
// enter the Cartesian reference sum. World's actual construction (M2a) sums
// ALL atoms (== m.numAtoms here); the pure mixed-manifold U_F^mixed (M2) sums
// only R's atoms. calcLogDetM itself is UNCHANGED by this knob (E's 0-DOF
// body is skipped by the dof==0 guard regardless), so this isolates exactly
// the Cartesian-reference term M2a's derivation is about.
double fixmanLike(const RobotModel& m, RobotState& s, int cartesianRefAtomCount) {
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    const double lnDetM = static_cast<double>(RobotEngine::calcLogDetM(m, s));
    double lnDetMCartesian = 0.0;
    for (int a = 0; a < cartesianRefAtomCount; ++a) {
        if (m.atomMass[static_cast<std::size_t>(a)] > Real(0)) {
            lnDetMCartesian += 3.0 * std::log(static_cast<double>(m.atomMass[static_cast<std::size_t>(a)]));
        }
    }
    return 0.5 * kT300 * (lnDetM - lnDetMCartesian);
}

} // namespace

// ===========================================================================
//  INV-WELD (50-validation.md Sec.2, PRIMARY) -- the ONE new code path
//  (10-mixed-integrator-correctness.md Sec.6 touch list): World::
//  setCartesianSolvent SHALL throw on a Free-rooted or partially-masked E,
//  and pass a correctly welded E. Drives the REAL World through buildModel
//  on a hand-built two-molecule SystemTopology (TestAlchemy.cpp's pattern --
//  no prmtop file, no OpenMM force evaluation needed for buildModel/
//  setCartesianSolvent).
// ===========================================================================
namespace {

// Two molecules, two atoms each, one intramolecular bond each (default Weld
// -- Selection::bondMobility left empty), no ring closures. Molecule 0 == R
// (root mobility supplied per test), molecule 1 == E.
SystemTopology twoMoleculeTopology() {
    SystemTopology sys;
    sys.numAtoms = 4;
    sys.numMolecules = 2;
    sys.atomsBegin = {0, 2};
    sys.atomsEnd = {2, 4};
    sys.atomsMass = {12.0, 1.0, 16.0, 1.0};
    sys.atomsX = {0.0, 0.10, 1.00, 1.10};
    sys.atomsY = {0.0, 0.0, 0.0, 0.0};
    sys.atomsZ = {0.0, 0.0, 0.0, 0.0};
    sys.numBonds = 2;
    sys.bondsI = {0, 2};
    sys.bondsJ = {1, 3};
    sys.bondsRingClosing = {false, false};
    return sys;
}

} // namespace

TEST(TwoRobotContact, INV_WELD_WeldedEPasses) {
    SystemTopology sys = twoMoleculeTopology();
    Selection sel; // empty -> both bonds default Rigid/Weld
    World world(0, /*cartesian=*/false, /*seed=*/1u);
    world.buildModel(sys, sel, {JointType::Free, JointType::Rigid});

    EXPECT_NO_THROW(world.setCartesianSolvent({2, 3}));

    const RobotModel& m = world.model();
    for (int a : {2, 3}) {
        const int b = m.atomBody[a];
        EXPECT_EQ(m.bodyNU[b], 0) << "E's atoms must sit in a 0-DOF body (P1)";
    }
}

TEST(TwoRobotContact, INV_WELD_FreeRootedEThrows) {
    SystemTopology sys = twoMoleculeTopology();
    Selection sel;
    World world(0, /*cartesian=*/false, /*seed=*/2u);
    // Wrong-arm: E's root is Free (6 DOF) -- nothing in the engine forbids
    // constructing this, so the guard is the only thing standing between it
    // and a silent double-counted-KE / double-drawn-momentum corruption.
    world.buildModel(sys, sel, {JointType::Free, JointType::Free});

    EXPECT_THROW(world.setCartesianSolvent({2, 3}), std::runtime_error);
}

TEST(TwoRobotContact, INV_WELD_PartiallyMaskedBodyThrows) {
    SystemTopology sys = twoMoleculeTopology();
    Selection sel;
    World world(0, /*cartesian=*/false, /*seed=*/3u);
    world.buildModel(sys, sel, {JointType::Free, JointType::Rigid});

    // Atom 3 (E's second atom, same 0-DOF body as atom 2) is NOT flagged: the
    // mask covers E's body only partially.
    EXPECT_THROW(world.setCartesianSolvent({2}), std::runtime_error);
}

// ===========================================================================
//  INV-FIX (10-mixed-integrator-correctness.md Sec.3, M2/M2a) -- Delta U_F
//  over a move of BOTH phi and x_E is bitwise identical whether or not E's
//  atoms sit in the Cartesian reference lnDetMCartesian_ (a canceling
//  constant); moving x_E alone leaves U_F exactly unchanged (E contributes NO
//  configuration-dependent Fixman term).
// ===========================================================================
TEST(TwoRobotContact, INV_FIX_DeltaUFIndependentOfEInCartesianReference) {
    RobotModel m = twoRobotFixture();
    RobotState s;
    initState(m, s);
    primePositions(m, s);
    const std::vector<int> eAt = eAtoms(m);
    s.setCartSolvent(eAt, m.atomMass.data());

    ASSERT_EQ(m.numAtoms, 7);
    const int numRAtoms = eAt.front(); // R's atoms are the leading contiguous block
    ASSERT_EQ(numRAtoms, 5);

    // ---- move x_E ONLY (phi fixed): U_F must be EXACTLY unchanged ------------
    RobotEngine::realizePosition(m, s);
    const double ufBefore = fixmanLike(m, s, m.numAtoms);
    Vec3* posG = s.atomPosG();
    for (int a : eAt) {
        posG[a] = posG[a] + Vec3(0.2, -0.15, 0.05);
    }
    const double ufAfterMoveXEOnly = fixmanLike(m, s, m.numAtoms);
    EXPECT_DOUBLE_EQ(ufBefore, ufAfterMoveXEOnly)
        << "moving x_E alone must leave U_F unchanged (M2: E contributes NO "
           "configuration-dependent Fixman term)";

    // ---- move BOTH phi and x_E: Delta U_F identical with/without E in the
    //      Cartesian reference (M2a's canceling -1/2 RT ln|M_E| constant) -----
    const double ufWithE_q0 = fixmanLike(m, s, m.numAtoms);
    const double ufRonly_q0 = fixmanLike(m, s, numRAtoms);

    s.q()[0] = Real(0.9);
    s.q()[1] = Real(-1.1);
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s); // skips E (cartSolvent-flagged)
    for (int a : eAt) {
        posG[a] = posG[a] + Vec3(-0.3, 0.4, -0.1); // move x_E again, alongside phi
    }
    const double ufWithE_q1 = fixmanLike(m, s, m.numAtoms);
    const double ufRonly_q1 = fixmanLike(m, s, numRAtoms);

    const double dUF_withE = ufWithE_q1 - ufWithE_q0;
    const double dUF_Ronly = ufRonly_q1 - ufRonly_q0;
    EXPECT_NEAR(dUF_withE, dUF_Ronly, 1e-10 * std::max(1.0, std::abs(dUF_Ronly)))
        << "Delta U_F over a joint (phi,x_E) move must be identical whether or not E's "
           "atoms sit in the Cartesian reference: dUF_withE=" << dUF_withE
        << " dUF_Ronly=" << dUF_Ronly;
}

// ===========================================================================
//  INV-DRAW (10-...md Sec.4, M3) -- the sampled covariance of (u, v_E)
//  matches RT*diag(M_phi^-1, M_E^-1), cross-block ~0.
// ===========================================================================
TEST(TwoRobotContact, INV_DRAW_BlockDiagonalCovarianceMatchesRTMinv) {
    RobotModel m = twoRobotFixture();
    RobotState s;
    initState(m, s);
    s.q()[0] = Real(0.3);
    s.q()[1] = Real(-0.4); // non-trivial configuration: M_phi(q) genuinely non-identity
    primePositions(m, s);
    const std::vector<int> eAt = eAtoms(m);
    s.setCartSolvent(eAt, m.atomMass.data());

    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);

    const int nu = m.nu;
    const int nE = static_cast<int>(eAt.size());
    const int dim = nu + 3 * nE;

    const std::vector<double> L = sqrtMInvMatrix(m, s); // same operator the draw uses (M3)

    std::vector<double> predicted(static_cast<std::size_t>(dim * dim), 0.0);
    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j < nu; ++j) {
            double acc = 0.0;
            for (int k = 0; k < nu; ++k) {
                acc += L[static_cast<std::size_t>(i * nu + k)] * L[static_cast<std::size_t>(j * nu + k)];
            }
            predicted[static_cast<std::size_t>(i * dim + j)] = kT300 * acc;
        }
    }
    const std::vector<Real>& invM = s.cartSolventInvMass();
    for (int j = 0; j < nE; ++j) {
        for (int c = 0; c < 3; ++c) {
            const int idx = nu + 3 * j + c;
            predicted[static_cast<std::size_t>(idx * dim + idx)] = kT300 * static_cast<double>(invM[static_cast<std::size_t>(j)]);
        }
    }
    // Cross-block (u vs v_E) predicted entries stay 0.0 (default) -- M1's zero
    // cross-Jacobian claim.

    std::mt19937_64 rng(0xD8A7);
    std::normal_distribution<double> gauss(0.0, 1.0);
    std::vector<MeanAccumulator> cov(static_cast<std::size_t>(dim * dim));
    const long nDraws = 300000;
    std::vector<Real> g(static_cast<std::size_t>(nu)), seeded(static_cast<std::size_t>(nu));
    std::vector<double> x(static_cast<std::size_t>(dim));
    const double scaleU = std::sqrt(kT300);
    for (long n = 0; n < nDraws; ++n) {
        for (int i = 0; i < nu; ++i) {
            g[static_cast<std::size_t>(i)] = static_cast<Real>(gauss(rng));
        }
        RobotEngine::multiplyBySqrtMInv(m, s, g.data(), seeded.data());
        for (int i = 0; i < nu; ++i) {
            x[static_cast<std::size_t>(i)] = scaleU * static_cast<double>(seeded[static_cast<std::size_t>(i)]);
        }
        for (int j = 0; j < nE; ++j) {
            const double sigma = std::sqrt(kT300 * static_cast<double>(invM[static_cast<std::size_t>(j)]));
            for (int c = 0; c < 3; ++c) {
                x[static_cast<std::size_t>(nu + 3 * j + c)] = sigma * gauss(rng);
            }
        }
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                cov[static_cast<std::size_t>(i * dim + j)].add(x[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(j)]);
            }
        }
    }

    // TOLERANCE PHILOSOPHY (StatTest.hpp): |observed-expected| <= 4*stderr.
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            const MeanAccumulator& acc = cov[static_cast<std::size_t>(i * dim + j)];
            const double expected = predicted[static_cast<std::size_t>(i * dim + j)];
            const double tol = std::max(4.0 * acc.stderrMean(), 1e-8);
            EXPECT_NEAR(acc.mean(), expected, tol)
                << "Cov(x_" << i << ",x_" << j << ") sample=" << acc.mean() << " expected=" << expected
                << " stderr=" << acc.stderrMean();
        }
    }
}

// ===========================================================================
//  INV-KE (50-validation.md Sec.2) -- ke (calcKineticEnergy, used in H)
//  matches 1/2 u^T M_phi u for M_phi INDEPENDENTLY reconstructed from
//  multiplyBySqrtMInv (used in the draw); E's mass metric is self-consistent
//  between cartSolventInvMass_ (the draw/Verlet path) and model.atomMass (the
//  reference).
// ===========================================================================
TEST(TwoRobotContact, INV_KE_MatchesIndependentlyReconstructedMetric) {
    RobotModel m = twoRobotFixture();
    RobotState s;
    initState(m, s);
    s.q()[0] = Real(0.5);
    s.q()[1] = Real(0.2);
    primePositions(m, s);
    const std::vector<int> eAt = eAtoms(m);
    s.setCartSolvent(eAt, m.atomMass.data());

    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);

    const int nu = m.nu;
    const std::vector<double> L = sqrtMInvMatrix(m, s);
    std::vector<double> LLt(static_cast<std::size_t>(nu * nu), 0.0);
    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j < nu; ++j) {
            double acc = 0.0;
            for (int k = 0; k < nu; ++k) {
                acc += L[static_cast<std::size_t>(i * nu + k)] * L[static_cast<std::size_t>(j * nu + k)];
            }
            LLt[static_cast<std::size_t>(i * nu + j)] = acc;
        }
    }
    std::vector<double> Mphi;
    invertSmall(LLt, nu, Mphi); // M_phi = (L L^T)^-1 -- the metric the draw implies

    std::mt19937_64 rng(0x77A1);
    std::normal_distribution<double> gauss(0.0, 1.0);
    std::vector<double> uD(static_cast<std::size_t>(nu));
    Real* u = s.u();
    for (int i = 0; i < nu; ++i) {
        uD[static_cast<std::size_t>(i)] = gauss(rng);
        u[i] = static_cast<Real>(uD[static_cast<std::size_t>(i)]);
    }
    RobotEngine::realizeVelocity(m, s);
    const double keActual = static_cast<double>(RobotEngine::calcKineticEnergy(m, s));

    double keRef = 0.0;
    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j < nu; ++j) {
            keRef += uD[static_cast<std::size_t>(i)] * Mphi[static_cast<std::size_t>(i * nu + j)] * uD[static_cast<std::size_t>(j)];
        }
    }
    keRef *= 0.5;

    EXPECT_NEAR(keActual, keRef, 1e-9 * std::max(1.0, std::abs(keRef)))
        << "calcKineticEnergy (used in H) disagrees with the metric multiplyBySqrtMInv "
           "implies (used in the momentum draw) -- draw/acceptance metric mismatch. "
           "keActual=" << keActual << " keRef=" << keRef;

    // E block: 1/2 sum m_a |v_a|^2 must be self-consistent between
    // cartSolventInvMass_ (drives calcSolventKE/the Cartesian Verlet) and
    // model.atomMass (the reference table setCartSolvent derives it from).
    Vec3* velG = s.atomVelG();
    const std::vector<Real>& invM = s.cartSolventInvMass();
    double keSolventViaInvMass = 0.0, keSolventViaAtomMass = 0.0;
    for (std::size_t j = 0; j < eAt.size(); ++j) {
        const int a = eAt[j];
        const Vec3 v(gauss(rng), gauss(rng), gauss(rng));
        velG[a] = v;
        const double mFromInv = (invM[j] > Real(0)) ? (1.0 / static_cast<double>(invM[j])) : 0.0;
        const double mDirect = static_cast<double>(m.atomMass[static_cast<std::size_t>(a)]);
        keSolventViaInvMass += 0.5 * mFromInv * static_cast<double>(dot(v, v));
        keSolventViaAtomMass += 0.5 * mDirect * static_cast<double>(dot(v, v));
    }
    EXPECT_NEAR(keSolventViaInvMass, keSolventViaAtomMass, 1e-12 * std::max(1.0, keSolventViaAtomMass))
        << "cartSolventInvMass_ disagrees with model.atomMass -- E's mass metric is "
           "inconsistent between the draw/Verlet path and the reference.";
}

// ===========================================================================
//  INV-REV (10-...md Sec.5, M4-risk) -- checkReversibility on the JOINT (R+E)
//  map at the two-robot working dt: round-trip residual <= 1e-6. Must be run
//  WITH E present (not R alone), per the spec.
// ===========================================================================
TEST(TwoRobotContact, INV_REV_JointRoundTripResidualWithinTolerance) {
    RobotModel m = twoRobotFixture();
    RobotState s;
    initState(m, s);
    s.q()[0] = Real(0.25);
    s.q()[1] = Real(-0.35);
    primePositions(m, s);
    const std::vector<int> eAt = eAtoms(m);
    s.setCartSolvent(eAt, m.atomMass.data());

    robo::ConstraintSet cs; // acyclic fixture: empty
    TwoRobotBridge bridge(m, s, /*k*/ Real(300.0), /*cacheAtomForces*/ true);
    seedDerivatives(m, s, bridge);

    // Seed nonzero momenta on BOTH blocks (a zero-momentum state would make
    // checkReversibility vacuously trivial).
    std::mt19937_64 rng(0x9C31);
    std::normal_distribution<double> gauss(0.0, 1.0);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    const int nu = m.nu;
    std::vector<Real> gvec(static_cast<std::size_t>(nu)), seeded(static_cast<std::size_t>(nu));
    for (int i = 0; i < nu; ++i) {
        gvec[static_cast<std::size_t>(i)] = static_cast<Real>(gauss(rng));
    }
    RobotEngine::multiplyBySqrtMInv(m, s, gvec.data(), seeded.data());
    const Real scaleU = static_cast<Real>(std::sqrt(kT300));
    Real* u = s.u();
    for (int i = 0; i < nu; ++i) {
        u[i] = scaleU * seeded[static_cast<std::size_t>(i)];
    }
    Vec3* velG = s.atomVelG();
    const std::vector<Real>& invM = s.cartSolventInvMass();
    for (std::size_t j = 0; j < eAt.size(); ++j) {
        const Real sigma = static_cast<Real>(std::sqrt(kT300 * static_cast<double>(invM[j])));
        velG[eAt[j]] = Vec3(sigma * static_cast<Real>(gauss(rng)),
                            sigma * static_cast<Real>(gauss(rng)),
                            sigma * static_cast<Real>(gauss(rng)));
    }
    seedDerivatives(m, s, bridge);

    const Real h = Real(4.0e-4);
    const Real residual = RobotEngine::checkReversibility(m, s, bridge, cs, /*nSteps*/ 6, h);

    ASSERT_TRUE(std::isfinite(residual)) << "reversibility probe diverged (dt too large for this fixture)";
    EXPECT_LE(residual, Real(1e-6)) << "joint (R+E) round-trip residual=" << residual << " exceeds 1e-6";
}
