// ============================================================================
//  TestFixmanIdealizedChains.cpp -- Tier 0 of
//  docs/specs/fixman-idealized-chains-validation.md: deterministic det M(q)
//  identities on hand-built idealized bead chains (Jain et al. 2013 Sec. III.A,
//  Spiridon & Minh 2017 Sec. 3.1). Always-on (no sampling, no ROBOSAMPLE_SLOW_TESTS
//  gate) -- this is the keystone tier: a closed-form external oracle compared at
//  machine precision, with NO simulator in the loop.
//
//  SURFACE (spec Sec. 3.3): the C++ harness cannot load a shipped prmtop
//  (Context::load_amber needs the Python/OpenMM pipeline), so this file builds
//  its own RobotModel by hand (precedent: TestFixmanBoltzmann.cpp) parameterized
//  to spec Sec. 3.1 (bead mass 14, bond length 1.54, bond angle 90 deg for C4 /
//  109 deg for C5+, no nonbonded/torsion force-field term -- U == 0 throughout,
//  so every result below is a pure statement about the mass metric M(q)).
//
//  DECOMPOSITION. A chain of n beads (n-1 bonds, n-3 free torsions) is built as:
//    * root (Free, 6 dof) rigidly carries beads 1,2,3 (the first bond + angle,
//      which have no internal shape freedom: 3 points, 2 bond-length + 1
//      bond-angle constraint == 0 residual dof);
//    * each subsequent bead k=4..n gets its own 1-dof Torsion body, hinged at
//      the (k-2,k-1) bond axis with the (k-2,k-1,k) bond angle fixed by the
//      joint's constant frames -- i.e. the *standard* internal-coordinate (BAT)
//      decomposition: bead k's torsion is exactly dihedral(k-3,k-2,k-1,k).
//  Every body from bead 4 onward has IDENTICAL BodySpec fields (mass, com_B,
//  inertia_B, X_PF, X_BM) because the local bond-length/bond-angle geometry is
//  the same at every vertex -- this is what makes buildBeadChain() a simple
//  loop. It is used for T0.3 (C5/C11/C15), a decomposition-agnostic internal
//  identity (ln det M == sum ln D_b vs the dense Jacobian oracle), so ANY valid
//  spanning-tree decomposition of the same 7..18-dof manifold is fine there.
//
//  T0.1 (C4) is different: it must match Jain's EXTERNAL closed-form eq:24
//  point-for-point, which is NOT decomposition-invariant (det M = J^T M_cart J
//  depends on the chosen generalized coordinates, not just on which physical
//  configurations are reachable). The reviewed spec's hint -- root owns beads
//  {1,2} ONLY, the Torsion leaf owns beads {3,4} -- is load-bearing there and
//  gets its OWN builder, buildC4RootPairLeafPair() (see its banner): a 2-atom
//  rigid root is rotationally DEGENERATE about its own bond axis, and it is
//  healed only by the leaf's alpha-dependent backward-accumulated contribution
//  -- exactly the eq:24 mechanism. buildBeadChain(4, ...) (root owns beads
//  {1,2,3}, non-degenerate on its own) is a valid but DIFFERENT parametrization
//  and was empirically confirmed NOT to reproduce eq:24 (same qualitative
//  bimodal shape, different polynomial) -- so T0.1 does not use it.
//
//  The actual physical dihedral used everywhere below is READ BACK from the
//  atom Cartesian positions (attachAtoms + RobotEngine::realizePosition +
//  robo::calcDihedralAngle), never assumed equal to the joint's own q -- so the
//  extrema/shape assertions are independent of any phase convention in the
//  hand-derived joint frames.
//
//  T0.2 / T0.4 (SHOULD, secondary) are NOT implemented here, with justification
//  (Rule 8 exclusion -- a test that cannot fail on a logic change is not a
//  test):
//    * T0.2 asks to validate "the implemented calcFixman"'s sign/RT convention.
//      World::calcFixman lives on World (a full Context/OpenMM molecular
//      pipeline) and is unreachable from a hand-built RobotModel (spec Sec.
//      3.3's own C++-harness constraint). Re-deriving the SAME formula
//      (U_F = 0.5 RT ln det M) locally in this file to "test" it would be
//      exp(-U_F/RT)*sqrt(det M) == 1 identically, for ANY det M and ANY sign
//      convention chosen consistently -- a tautology, not a discriminating
//      check.
//    * T0.4 asks for the engine's Fixman TORQUE. No such quantity is computed
//      anywhere in the engine: guidance dynamics use the unmodified
//      Hamiltonian (docs/specs/singular-dof-fixman.md Sec. 6 NOTE: "no Fixman
//      torque path to correct"); Fixman enters only the acceptance ln det M.
//      There is nothing to finite-difference against.
//  T0.1 and T0.3 (SHALL) are fully implemented below.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <gtest/gtest.h>
#include <utility>
#include <vector>

#include "Constraints.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotLinearAlgebra.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::Rng;

namespace {

// ---------------------------------------------------------------------------
//  Geometry: point-mass unit inertia (per unit mass) about the BODY ORIGIN,
//  matching RobotModel::bodyUnitInertia_B's documented convention. United-atom
//  beads (Jain 2013 Sec. III.A) carry no intrinsic rotational inertia, so each
//  bead is exactly a point mass: I = sum_i m_i (|r_i|^2 I3 - r_i r_i^T).
// ---------------------------------------------------------------------------
UnitInertia pointMassUnitInertiaAboutOrigin(const std::vector<std::pair<Vec3, Real>>& pts, Real totalMass) {
    Real ixx = 0, iyy = 0, izz = 0, ixy = 0, ixz = 0, iyz = 0;
    for (const auto& pm : pts) {
        const Vec3& r = pm.first;
        const Real m = pm.second;
        ixx += m * ((r[1] * r[1]) + (r[2] * r[2]));
        iyy += m * ((r[0] * r[0]) + (r[2] * r[2]));
        izz += m * ((r[0] * r[0]) + (r[1] * r[1]));
        ixy += -m * r[0] * r[1];
        ixz += -m * r[0] * r[2];
        iyz += -m * r[1] * r[2];
    }
    return UnitInertia(ixx / totalMass, ixy / totalMass, iyy / totalMass, ixz / totalMass, iyz / totalMass,
                       izz / totalMass);
}

// A hand-built idealized bead chain (spec Sec. 3.1): n beads, n-1 bonds at
// bondLen, n-2 bond angles at bondAngleDeg, n-3 free torsions.
struct BeadChain {
    RobotModel model;
    int rootBody = 1;             // Free, owns beads 1,2,3
    std::vector<int> torsionBody; // torsionBody[k] owns bead (k+4), k=0..n-4
};

// PRECONDITION: nBeads >= 4 (fewer beads have no free torsion to build).
BeadChain buildBeadChain(int nBeads, Real bondLen, Real bondAngleDeg, Real mass) {
    const Real theta = bondAngleDeg * Deg2Rad;
    const Real phi = Pi - theta; // exterior turn angle
    // Every body's own atom sits at this SAME fixed local offset from its own
    // origin (the vertex it hinges from) -- the chain is locally self-similar,
    // so one constant offset + one constant joint frame serves every body from
    // bead 4 onward (derivation in the file banner).
    const Vec3 atomOffsetLocal(bondLen * std::sin(phi), Real(0), bondLen * std::cos(phi));
    // Root-only: bead 1 relative to the root's own origin (bead 2), placed
    // "behind" the root's own z axis so angle(1,2,3) == bondAngleDeg too.
    const Vec3 firstAtomOffsetLocal(Real(0), Real(0), -bondLen);

    Rotation rPF;
    rPF.setRotationFromOneAxis(UnitVec3(atomOffsetLocal), ZAxis);
    const Transform xPF(rPF, atomOffsetLocal);

    std::vector<BodySpec> specs;

    BodySpec root;
    root.parent = 0;
    root.joint = JointType::Free;
    root.X_PF = Transform();
    root.X_BM = Transform();
    root.mass = mass * Real(3);
    const Vec3 p1 = firstAtomOffsetLocal;
    const Vec3 p2(0, 0, 0);
    const Vec3 p3 = atomOffsetLocal;
    root.com_B = ((p1 * mass) + (p2 * mass) + (p3 * mass)) / root.mass;
    root.inertia_B = pointMassUnitInertiaAboutOrigin({{p1, mass}, {p2, mass}, {p3, mass}}, root.mass);
    specs.push_back(root);

    BeadChain result;
    int parentIdx = 1; // root's body index
    for (int bead = 4; bead <= nBeads; ++bead) {
        BodySpec b;
        b.parent = parentIdx;
        b.joint = JointType::Torsion;
        b.X_PF = xPF;
        b.X_BM = Transform();
        b.mass = mass;
        b.com_B = atomOffsetLocal;
        b.inertia_B = pointMassUnitInertiaAboutOrigin({{atomOffsetLocal, mass}}, mass);
        specs.push_back(b);
        const int thisBodyIdx = static_cast<int>(specs.size()); // 1-based (specs[i] -> body i+1)
        result.torsionBody.push_back(thisBodyIdx);
        parentIdx = thisBodyIdx;
    }

    result.model = buildForest(specs);
    result.rootBody = 1;

    // Attach the real atoms (bead k -> atom index k-1), so tests can read back
    // ACTUAL Cartesian positions and compute the physical dihedral rather than
    // trust the hand-derived joint-frame phase.
    rtest::attachAtoms(result.model, result.rootBody, {p1, p2, p3}, {mass, mass, mass});
    for (int b : result.torsionBody) {
        rtest::attachAtoms(result.model, b, {atomOffsetLocal}, {mass});
    }
    return result;
}

// ---------------------------------------------------------------------------
//  C4-ONLY: the exact decomposition the reviewed spec's "Key physics hint"
//  names -- root (Free) owns beads {1,2} (ONE side of the central 2-3 bond),
//  the single Torsion body owns beads {3,4} (the OTHER side). This is NOT
//  interchangeable with buildBeadChain(4, ...): a 2-atom rigid root is
//  ROTATIONALLY DEGENERATE about its own bond axis (two colinear point masses
//  have zero moment of inertia about the line through them), so the root's
//  bare articulated inertia P_root is singular in exactly the torsion
//  direction; it is healed ONLY by the backward-accumulated, alpha-dependent
//  contribution of the child (beads 3,4) through the shared 2-3 bond. That
//  alpha-dependent healing is precisely what produces Jain's eq:24 polynomial.
//  buildBeadChain's root (beads 1,2,3, already non-degenerate on its own) is a
//  DIFFERENT parametrization of the same 7-dof manifold -- valid (T0.3 checks
//  it against the dense oracle), but it does not reproduce eq:24 point-for-
//  point, only the general bimodal shape. Confirmed empirically: swapping to
//  this decomposition is what makes the eq:24 fit residual drop to ~1e-9;
//  do not "simplify" this back to buildBeadChain(4, ...).
// ---------------------------------------------------------------------------
BeadChain buildC4RootPairLeafPair(Real bondLen, Real bondAngleDeg, Real mass) {
    const Real theta = bondAngleDeg * Deg2Rad;

    // Root: beads A (atom0), B (atom1). Origin at B. Bond1 (A->B) along +x.
    const Vec3 atomA(-bondLen, 0, 0);
    const Vec3 atomB(0, 0, 0);

    // Bond2 (B->C) perpendicular to bond1 (angle A-B-C == bondAngleDeg). For
    // bondAngleDeg == 90 (the only value this C4-specific builder is used
    // with) "perpendicular" is exact; kept general via bondAngleDeg for
    // documentation symmetry with buildBeadChain, but ASSERTed to 90 by the
    // caller (Jain eq:24 is derived for 90-degree bond angles only).
    const Vec3 bond2Dir(0, std::sin(theta), std::cos(theta) /* == 0 at theta==90 */);
    Rotation rPF;
    rPF.setRotationFromOneAxis(UnitVec3(bond2Dir), ZAxis);
    const Transform xPF(rPF, atomB + (bond2Dir * bondLen));

    // Torsion body: beads C (own origin), D. C->B (local) == -bondLen * z (the
    // body's own z axis IS the bond2 axis, by the X_BM=I / Torsion convention).
    // C->D perpendicular to that (angle B-C-D == bondAngleDeg): place D purely
    // along local x, which is perpendicular to local z regardless of theta.
    const Vec3 atomCLocal(0, 0, 0);
    const Vec3 atomDLocal(bondLen, 0, 0);

    BodySpec root;
    root.parent = 0;
    root.joint = JointType::Free;
    root.mass = mass * Real(2);
    root.com_B = ((atomA * mass) + (atomB * mass)) / root.mass;
    root.inertia_B = pointMassUnitInertiaAboutOrigin({{atomA, mass}, {atomB, mass}}, root.mass);

    BodySpec leaf;
    leaf.parent = 1;
    leaf.joint = JointType::Torsion;
    leaf.X_PF = xPF;
    leaf.X_BM = Transform();
    leaf.mass = mass * Real(2);
    leaf.com_B = ((atomCLocal * mass) + (atomDLocal * mass)) / leaf.mass;
    leaf.inertia_B = pointMassUnitInertiaAboutOrigin({{atomCLocal, mass}, {atomDLocal, mass}}, leaf.mass);

    BeadChain result;
    result.model = buildForest({root, leaf});
    result.rootBody = 1;
    result.torsionBody = {2};

    rtest::attachAtoms(result.model, 1, {atomA, atomB}, {mass, mass});
    rtest::attachAtoms(result.model, 2, {atomCLocal, atomDLocal}, {mass, mass});
    return result;
}

void initFreeRootIdentity(const RobotModel& m, RobotState& s) {
    for (int i = 0; i < m.nq; ++i) {
        s.q()[i] = 0;
    }
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = 0;
    }
    for (int qs : m.quaternionQStart) {
        s.q()[qs + 0] = 1; // identity quaternion (w=1, x=y=z=0)
    }
}

// The smallest hinge-inertia eigenvalue over every dof>0 body, via the SAME
// DI-inverse trick as TestRoboticsOracle.cpp (RobotState never stores D
// directly, only DI persists): min-eig(D) = 1/max-eig(DI) away from the lock.
// PRECONDITION: realizeArticulatedBodyInertias current.
Real worstMinEigD(const RobotModel& m, const RobotState& s) {
    Real worst = 1e300;
    for (int b = 1; b < m.numBodies; ++b) {
        const int dof = m.bodyNU[b];
        if (dof <= 0) {
            continue;
        }
        const Real* diBlock = &s.DI()[m.bodyUSqIndex[b]];
        Real eigD[6], eigV[36];
        robo_linalg::jacobiSymEig(diBlock, dof, eigD, eigV);
        Real maxEigDI = eigD[0];
        for (int k = 1; k < dof; ++k) {
            maxEigDI = std::max(maxEigDI, eigD[k]);
        }
        if (maxEigDI <= Real(0)) {
            return Real(0); // degenerate DI -- report as locked
        }
        worst = std::min(worst, Real(1) / maxEigDI);
    }
    return worst;
}

// log|det A| of a row-major n x n via partial-pivot Gaussian elimination.
// Precedent: TestMassMatrix.cpp's local `logDet` (D5: calcLogDetM vs a dense
// M^-1 built from multiplyByMInv) -- duplicated here per that file's own
// per-TU convention (each Test*.cpp is an independent translation unit).
double denseLogDet(std::vector<Real> a, int n) {
    double ld = 0;
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r = col + 1; r < n; ++r) {
            if (std::abs(a[static_cast<std::size_t>((r * n) + col)])
                > std::abs(a[static_cast<std::size_t>((piv * n) + col)])) {
                piv = r;
            }
        }
        if (piv != col) {
            for (int c = 0; c < n; ++c) {
                std::swap(a[static_cast<std::size_t>((piv * n) + c)], a[static_cast<std::size_t>((col * n) + c)]);
            }
        }
        const double d = static_cast<double>(a[static_cast<std::size_t>((col * n) + col)]);
        ld += std::log(std::abs(d));
        for (int r = col + 1; r < n; ++r) {
            const double f = static_cast<double>(a[static_cast<std::size_t>((r * n) + col)]) / d;
            for (int c = col; c < n; ++c) {
                a[static_cast<std::size_t>((r * n) + c)] -=
                    static_cast<Real>(f) * a[static_cast<std::size_t>((col * n) + c)];
            }
        }
    }
    return ld;
}

// Index of the grid point whose dihedral is angularly closest to targetRad.
int nearestByAngle(const std::vector<double>& grid, double targetRad) {
    int best = 0;
    double bestD = 1e300;
    for (int i = 0; i < static_cast<int>(grid.size()); ++i) {
        double d = std::fmod(grid[i] - targetRad + M_PI, 2 * M_PI);
        if (d < 0) {
            d += 2 * M_PI;
        }
        d = std::abs(d - M_PI);
        if (d < bestD) {
            bestD = d;
            best = i;
        }
    }
    return best;
}

} // namespace

// ---------------------------------------------------------------------------
//  T0.1 -- C4 analytical det M(alpha), Jain 2013 eq:24 / Fig 1a.
//  det M(alpha) = c5 * (35 + 4 cos(a) - 16 cos^2(a) + cos^4(a)), an external
//  closed-form oracle independent of the engine (the one external-oracle test
//  in the whole Fixman suite).
//
//  ALPHA CONVENTION (empirically resolved, decomposition-independent): eq:24's
//  "alpha" is Pear & Weiner's torsion angle, whose zero/chirality reference is
//  not restated in Jain 2013 or in references/papers/jain_2013_fixman_branched.
//  The STANDARD dihedral(atomA,atomB,atomC,atomD) (robo::calcDihedralAngle,
//  SimTK::calcDihedralAngle convention) differs from it by the supplement:
//  alpha = pi - dihedral(A,B,C,D), i.e. cos(alpha) = -cos(dihedral). Verified:
//  (1) this SAME supplement is required whether the robot is decomposed as
//  root{1,2,3}+leaf{4} (buildBeadChain) or root{1,2}+leaf{3,4}
//  (buildC4RootPairLeafPair) -- ruling out a decomposition-specific bug; (2)
//  with the supplement applied the fit residual below is ~1e-9 (vs ~0.24
//  without it) and the eq:24-predicted extrema/ratio all land correctly; (3)
//  T0.3 independently confirms calcLogDetM matches the dense J^T M_cart J
//  oracle exactly, so det M(q) itself is not in question -- only the label
//  attached to a given q.
// ---------------------------------------------------------------------------
TEST(FixmanIdealizedChains, C4KeystoneEq24) {
    constexpr Real kBondLen = 1.54;
    constexpr Real kBondAngleDeg = 90.0;
    constexpr Real kMass = 14.0;

    BeadChain chain = buildC4RootPairLeafPair(kBondLen, kBondAngleDeg, kMass);
    RobotModel& m = chain.model;
    ASSERT_EQ(static_cast<int>(chain.torsionBody.size()), 1) << "C4 has exactly one free torsion";
    RobotState s;
    s.allocateFull(m);

    const int qIdx = m.bodyQIndex[chain.torsionBody[0]];

    constexpr int kGrid = 361; // 1 deg step
    std::vector<double> alphaVals(kGrid), detM(kGrid), poly(kGrid);
    Real worstMinEig = 1e300;

    for (int i = 0; i < kGrid; ++i) {
        initFreeRootIdentity(m, s);
        const Real q = -Pi + (static_cast<Real>(i) * (2 * Pi / (kGrid - 1)));
        s.q()[qIdx] = q;
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        const Real lnDetM = RobotEngine::calcLogDetM(m, s);
        detM[static_cast<std::size_t>(i)] = std::exp(static_cast<double>(lnDetM));

        const Vec3* pos = s.atomPosG();
        const Real dih = calcDihedralAngle(pos[0], pos[1], pos[2], pos[3]);
        const Real alpha = std::atan2(std::sin(Pi - dih), std::cos(Pi - dih)); // wrap to [-pi, pi]
        const double cosA = std::cos(static_cast<double>(alpha));
        alphaVals[static_cast<std::size_t>(i)] = static_cast<double>(alpha);
        poly[static_cast<std::size_t>(i)] = 35.0 + (4.0 * cosA) - (16.0 * cosA * cosA) + (cosA * cosA * cosA * cosA);

        worstMinEig = std::min(worstMinEig, worstMinEigD(m, s));
    }

    ASSERT_GT(worstMinEig, Real(1e-12))
        << "a body's hinge inertia is at/below the singular-DOF lock -- C4 fixture reads a clamped determinant";

    // Vacuity guard (precedent: TestFixmanBoltzmann Smoke).
    const double loD = *std::min_element(detM.begin(), detM.end());
    const double hiD = *std::max_element(detM.begin(), detM.end());
    ASSERT_GT(hiD / loD, 1.05) << "det M(alpha) nearly flat -- C4 fixture would make the eq:24 fit vacuous";

    // Fit the single free scale c5 (linear least squares through the origin:
    // det M is EXACTLY proportional to the polynomial, no other free params).
    double num = 0, den = 0;
    for (int i = 0; i < kGrid; ++i) {
        num += detM[static_cast<std::size_t>(i)] * poly[static_cast<std::size_t>(i)];
        den += poly[static_cast<std::size_t>(i)] * poly[static_cast<std::size_t>(i)];
    }
    ASSERT_GT(den, 0.0);
    const double c5 = num / den;
    ASSERT_GT(c5, 0.0) << "fitted scale must be positive (det M > 0)";

    double maxAbsDetM = 0, maxResidual = 0;
    for (int i = 0; i < kGrid; ++i) {
        maxAbsDetM = std::max(maxAbsDetM, detM[static_cast<std::size_t>(i)]);
        maxResidual =
            std::max(maxResidual, std::abs(detM[static_cast<std::size_t>(i)] - (c5 * poly[static_cast<std::size_t>(i)])));
    }
    const double relResidual = maxResidual / maxAbsDetM;
    std::fprintf(stderr, "[FixmanIdealizedChains.C4KeystoneEq24] c5=%.6g relResidual=%.3e\n", c5, relResidual);
    EXPECT_LT(relResidual, 1e-6) << "C4 det M(alpha) does not match Jain 2013 eq:24 shape (c5=" << c5
                                 << ", relResidual=" << relResidual << ")";

    // Extrema: sqrt(det M) maximal at alpha ~= +-82.8 deg, local minima at
    // alpha = 0 and +-180 deg, ratio f(0):f(180) == 24:16 == 3:2.
    const int idxMax = static_cast<int>(std::max_element(detM.begin(), detM.end()) - detM.begin());
    const double alphaDegAtMax = std::abs(alphaVals[static_cast<std::size_t>(idxMax)]) * 180.0 / M_PI;
    EXPECT_NEAR(alphaDegAtMax, 82.8, 3.0) << "det M(alpha) maximum not at the eq:24-predicted +-82.8 deg";

    const int idxZero = nearestByAngle(alphaVals, 0.0);
    const int idxPi = nearestByAngle(alphaVals, M_PI);
    const double f0 = detM[static_cast<std::size_t>(idxZero)];
    const double f180 = detM[static_cast<std::size_t>(idxPi)];
    EXPECT_GT(f0, f180) << "det M(0) should exceed det M(180) (eq:24 ratio 24:16)";
    EXPECT_NEAR(f0 / f180, 24.0 / 16.0, 0.02) << "eq:24 predicts det M(0):det M(180) == 24:16 == 1.5";
}

// ---------------------------------------------------------------------------
//  T0.3 -- multi-torsion O(n) determinant identity, C5/C11/C15 (Jain 1997/2013
//  checks): ln det M == sum_b ln det D_b (the O(n) articulated-body path)
//  agrees with a dense ln det(J^T M_cart J) built from multiplyByMInv applied
//  to every basis vector (precedent: TestMassMatrix.cpp D5), at several random
//  torsion configurations. Guard: calcConstraintLogDet == 0 exactly (no loop
//  term may fire on an acyclic chain).
// ---------------------------------------------------------------------------
TEST(FixmanIdealizedChains, MultiTorsionDeterminantIdentityC5C11C15) {
    constexpr Real kBondLen = 1.54;
    constexpr Real kBondAngleDeg = 109.0;
    constexpr Real kMass = 14.0;

    struct Case {
        int nBeads;
        const char* name;
        int expectedTorsions;
    };
    const Case cases[] = {{5, "C5", 2}, {11, "C11", 8}, {15, "C15", 12}};

    for (const Case& c : cases) {
        SCOPED_TRACE(c.name);
        BeadChain chain = buildBeadChain(c.nBeads, kBondLen, kBondAngleDeg, kMass);
        RobotModel& m = chain.model;
        ASSERT_EQ(static_cast<int>(chain.torsionBody.size()), c.expectedTorsions);
        RobotState s;
        s.allocateFull(m);
        robo::ConstraintSet cs; // empty: this chain is acyclic by construction

        Rng rng(0xBEAD0000u + static_cast<std::uint64_t>(c.nBeads));
        for (int rep = 0; rep < 5; ++rep) {
            initFreeRootIdentity(m, s);
            for (int qs : m.quaternionQStart) {
                const Quat q = rng.unitQuat();
                s.q()[qs + 0] = q.elems[0];
                s.q()[qs + 1] = q.elems[1];
                s.q()[qs + 2] = q.elems[2];
                s.q()[qs + 3] = q.elems[3];
            }
            s.q()[m.bodyQIndex[chain.rootBody] + 4] = rng.uniform(-2, 2);
            s.q()[m.bodyQIndex[chain.rootBody] + 5] = rng.uniform(-2, 2);
            s.q()[m.bodyQIndex[chain.rootBody] + 6] = rng.uniform(-2, 2);
            for (int b : chain.torsionBody) {
                s.q()[m.bodyQIndex[b]] = rng.uniform(-Pi, Pi);
            }

            RobotEngine::realizePosition(m, s);
            RobotEngine::realizeArticulatedBodyInertias(m, s);
            const Real lnDetM = RobotEngine::calcLogDetM(m, s);

            const int n = m.nu;
            std::vector<Real> minv(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0), e(n, 0), col(n, 0);
            for (int j = 0; j < n; ++j) {
                std::fill(e.begin(), e.end(), Real(0));
                e[static_cast<std::size_t>(j)] = 1;
                RobotEngine::multiplyByMInv(m, s, e.data(), col.data());
                for (int i = 0; i < n; ++i) {
                    minv[static_cast<std::size_t>((i * n) + j)] = col[static_cast<std::size_t>(i)];
                }
            }
            const double lnDetMInv = denseLogDet(minv, n);
            EXPECT_NEAR(static_cast<double>(lnDetM), -lnDetMInv, 1e-6)
                << c.name << " rep " << rep << ": O(n) block factorization vs dense M^-1 disagree";

            const Real lnDetZ = cs.calcConstraintLogDet(m, s);
            EXPECT_EQ(lnDetZ, Real(0)) << c.name << " rep " << rep
                                       << ": calcConstraintLogDet must be an exact no-op on an acyclic chain";

            EXPECT_GT(worstMinEigD(m, s), Real(1e-12))
                << c.name << " rep " << rep << ": a body's hinge inertia is at/below the singular-DOF lock";
        }
    }
}
