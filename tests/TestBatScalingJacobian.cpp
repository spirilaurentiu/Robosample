// ============================================================================
//  TestBatScalingJacobian.cpp -- V5 (engine-level Jacobian gate), the highest-
//  risk check in docs/specs/replica-exchange-nonequilibrium-work.md (D6/F7):
//  a wrong Jacobian sign silently biases the entire nonequilibrium sampler.
//
//  Calls the ENGINE directly: World::buildModel (real z-matrix population,
//  src/World.cpp) + World::previewBatScaling (robo::applyBatScaling /
//  robo::calcBatVolumeLogJac, src/BatScaling.cpp) -- the SAME code path
//  World::applyBatScalingDrive uses in production. NOT the pybind layer, NOT
//  numpy (tests/bat_jacobian_scaling_check.py remains the separate,
//  standalone analytic half; untouched by this file).
//
//  Three cases, all on hand-built two/four-atom SystemTopology fixtures
//  (TestTwoRobotContact.cpp/TestAlchemy.cpp's established "no prmtop file"
//  convention -- World::buildModel needs only the topology/mobility inputs):
//
//   * SingleSliderBondMatches3LnS (V5a): N_scaled=1 (r only), the D6-cited
//     numeric example lnJac = 3 ln(s) exactly.
//   * ThreeBodyChainMatchesClosedForm (V5c): a >=3-BODY chain (three
//     consecutive BendStretch bonds) scaling an UPSTREAM bond AND an
//     upstream angle together -- the discriminating case (spec V5 note): a
//     geometric-BAT/mobilizer-coordinate mismatch in the z-matrix ancestor
//     lookup (zJ/zK reading ALREADY-CASCADED positions) would show up ONLY
//     here, not in a single-DOF case. N_scaled=5: the first body is
//     root-adjacent (its grandparent is Ground, so its angle DOF has no zK
//     and is skipped, BatScaling.hpp's documented limitation -- contributes
//     r-only); the other two each contribute r+theta.
//   * Both cases include an F7 REGRESSION GUARD: the original's inverted-sign,
//     double-counted composition (J_ini + J_scale - J_fin, D6) is computed
//     INLINE and asserted to differ from the engine's corrected lnJac by far
//     more than the finite-difference tolerance, for every s != 1 -- so a
//     future regression to the buggy formula fails loudly.
//   * AggressiveScaleThrowsDomainError (S2, reviewer 2026-07-12): an s large
//     enough to push a scaled angle outside (0, pi) MUST throw
//     std::domain_error (BatScaling.cpp's fail-loud domain guard), not
//     silently return a corrupted lnJac. Context::driveReplica (Context.cpp,
//     Stage 2b) catches exactly this exception type and converts it to an
//     automatic swap reject -- see tests/TestRexAcceptanceAlgebra.cpp's
//     WORK_Jacobian=-infinity consuming-side test for the other half of that
//     mechanism (this file only exercises the THROW; a live Context/OpenMM
//     harness to exercise Context::driveReplica's try/catch end-to-end does
//     not exist in this test suite -- flagged in the coder checkpoint).
//
//  Fixture geometry is chosen so every measured angle is EXACTLY 90 degrees
//  (orthogonal zigzag: bond0 along +X, bond1 along +Y, bond2 along +Z), which
//  is verifiable by hand (dot product of consecutive bond directions is
//  exactly 0) and stays safely inside (0, pi) for every s used except the
//  deliberately-aggressive S2 case.
//
//  NOT compiled or run (coordinator directive, 2026-07-12: "drop compiling
//  and running entirely") -- written to the tests/Test*.cpp convention
//  (auto-discovered by CMakeLists.txt's robosample_add_test glob) for the
//  user to build and run manually.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BatScaling.hpp"
#include "RobotModel.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"
#include "robot_math.hpp"

using robo::Vec3;

namespace {

// ---------------------------------------------------------------------------
//  Fixtures (hand-built SystemTopology, no prmtop file -- World::buildModel
//  only reads numAtoms/numMolecules/atomsBegin/atomsMass/numBonds/bondsI/
//  bondsJ/bondsRingClosing; atomsX/Y/Z are set for completeness but are not
//  read by buildModel itself -- the actual geometry fed to previewBatScaling
//  is the explicit x0 vector each test builds).
// ---------------------------------------------------------------------------
SystemTopology twoAtomSliderTopology() {
    SystemTopology sys;
    sys.numAtoms = 2;
    sys.numMolecules = 1;
    sys.atomsBegin = {0};
    sys.atomsEnd = {2};
    sys.atomsMass = {12.0, 12.0};
    sys.atomsX = {0.0, 0.15};
    sys.atomsY = {0.0, 0.0};
    sys.atomsZ = {0.0, 0.0};
    sys.numBonds = 1;
    sys.bondsI = {0};
    sys.bondsJ = {1};
    sys.bondsRingClosing = {false};
    return sys;
}

// Orthogonal zigzag: atom0->atom1 along +X, atom1->atom2 along +Y,
// atom2->atom3 along +Z, all bonds 0.15 nm. angle(atom2,atom1,atom0) and
// angle(atom3,atom2,atom1) are EXACTLY 90 degrees (perpendicular consecutive
// bond directions, dot product 0 by construction).
SystemTopology fourAtomBendStretchChainTopology() {
    SystemTopology sys;
    sys.numAtoms = 4;
    sys.numMolecules = 1;
    sys.atomsBegin = {0};
    sys.atomsEnd = {4};
    sys.atomsMass = {12.0, 12.0, 12.0, 12.0};
    sys.atomsX = {0.0, 0.15, 0.15, 0.15};
    sys.atomsY = {0.0, 0.0, 0.15, 0.15};
    sys.atomsZ = {0.0, 0.0, 0.0, 0.15};
    sys.numBonds = 3;
    sys.bondsI = {0, 1, 2};
    sys.bondsJ = {1, 2, 3};
    sys.bondsRingClosing = {false, false, false};
    return sys;
}

std::vector<Vec3> fourAtomChainPositions() {
    return {Vec3(0.0, 0.0, 0.0), Vec3(0.15, 0.0, 0.0), Vec3(0.15, 0.15, 0.0), Vec3(0.15, 0.15, 0.15)};
}

// ---------------------------------------------------------------------------
//  Linear algebra: log|det(A)| of an n x n matrix via LU decomposition with
//  partial pivoting (standard algorithm; A is copied, not mutated in place
//  from the caller's view). Returns {NaN, 0} if a pivot is (numerically)
//  singular.
// ---------------------------------------------------------------------------
struct LuLogDet {
    double logAbsDet;
    int sign; // +1 / -1, 0 if singular
};

LuLogDet luLogDet(std::vector<double> a, int n) {
    double logAbsDet = 0.0;
    int sign = 1;
    for (int col = 0; col < n; ++col) {
        int piv = col;
        double best = std::abs(a[static_cast<std::size_t>(col) * n + col]);
        for (int r = col + 1; r < n; ++r) {
            const double v = std::abs(a[static_cast<std::size_t>(r) * n + col]);
            if (v > best) {
                best = v;
                piv = r;
            }
        }
        if (best < 1e-14) {
            return {std::nan(""), 0};
        }
        if (piv != col) {
            for (int c = 0; c < n; ++c) {
                std::swap(a[static_cast<std::size_t>(col) * n + c], a[static_cast<std::size_t>(piv) * n + c]);
            }
            sign = -sign;
        }
        const double pivotVal = a[static_cast<std::size_t>(col) * n + col];
        logAbsDet += std::log(std::abs(pivotVal));
        sign *= (pivotVal < 0.0) ? -1 : 1;
        for (int r = col + 1; r < n; ++r) {
            const double factor = a[static_cast<std::size_t>(r) * n + col] / pivotVal;
            for (int c = col; c < n; ++c) {
                a[static_cast<std::size_t>(r) * n + c] -= factor * a[static_cast<std::size_t>(col) * n + c];
            }
        }
    }
    return {logAbsDet, sign};
}

// Central-difference log|det dx'/dx| of World::previewBatScaling (FIXED s,
// anchor), over the FULL 3*numAtoms Cartesian input -- the whole map, not
// just the atoms a single body's own scaling touches, so cross-terms from an
// upstream body's rigid cascade into a downstream body are exercised too.
double fdLogAbsDetJac(const World& w,
                      const std::vector<Vec3>& x0,
                      double s,
                      const std::unordered_map<int, double>& anchorR = {},
                      const std::unordered_map<int, double>& anchorTheta = {},
                      double h = 1e-6) {
    const int n = static_cast<int>(x0.size()) * 3;
    std::vector<double> jac(static_cast<std::size_t>(n) * n, 0.0);
    for (int k = 0; k < n; ++k) {
        std::vector<Vec3> xp = x0;
        std::vector<Vec3> xm = x0;
        xp[static_cast<std::size_t>(k / 3)][k % 3] += h;
        xm[static_cast<std::size_t>(k / 3)][k % 3] -= h;
        const auto rp = w.previewBatScaling(xp, s, anchorR, anchorTheta);
        const auto rm = w.previewBatScaling(xm, s, anchorR, anchorTheta);
        for (int i = 0; i < n; ++i) {
            const double vp = rp.atomPos[static_cast<std::size_t>(i / 3)][i % 3];
            const double vm = rm.atomPos[static_cast<std::size_t>(i / 3)][i % 3];
            jac[static_cast<std::size_t>(i) * n + k] = (vp - vm) / (2.0 * h);
        }
    }
    const LuLogDet d = luLogDet(jac, n);
    EXPECT_GT(d.sign, 0) << "map orientation reversed -- something is badly wrong";
    return d.logAbsDet;
}

double angleAt(const Vec3& a, const Vec3& vertex, const Vec3& c) {
    const Vec3 v1 = a - vertex;
    const Vec3 v2 = c - vertex;
    double cosTheta = robo::dot(v1, v2) / (v1.norm() * v2.norm());
    cosTheta = std::min(1.0, std::max(-1.0, cosTheta));
    return std::acos(cosTheta);
}

} // namespace

// ---------------------------------------------------------------------------
//  V5a -- single Slider bond (N_scaled=1, r only): lnJac == 3 ln(s) exactly
//  (D6's own numeric counter-example to the F7 bug).
// ---------------------------------------------------------------------------
TEST(BatScalingJacobian, SingleSliderBondMatches3LnS) {
    SystemTopology sys = twoAtomSliderTopology();
    Selection sel;
    sel.bondMobility = {JointType::Slider};
    World world(0, /*cartesian=*/false, /*seed=*/11u);
    world.buildModel(sys, sel, {JointType::Free});

    const std::vector<Vec3> x0 = {Vec3(0.0, 0.0, 0.0), Vec3(0.15, 0.0, 0.0)};

    for (const double s : {1.3, 0.8}) {
        const World::BatScalingResult result = world.previewBatScaling(x0, s);
        EXPECT_EQ(result.nScaled, 1) << "s=" << s;

        const double want = 3.0 * std::log(s);
        EXPECT_NEAR(result.lnJac, want, 1e-9) << "analytic formula, s=" << s;

        const double fd = fdLogAbsDetJac(world, x0, s);
        EXPECT_NEAR(fd, want, 5e-4) << "ENGINE finite-difference lnJac != 3 ln(s) -- V5 GATE, s=" << s;
        EXPECT_NEAR(fd, result.lnJac, 5e-4) << "engine FD vs analytic formula, s=" << s;

        // F7 regression guard: J_ini + J_scale - J_fin (original composition,
        // D6), bond ratio double-counted, MUST NOT match the corrected lnJac.
        const double r0 = (x0[1] - x0[0]).norm();
        const double r1 = (result.atomPos[1] - result.atomPos[0]).norm();
        const double jIni = 2.0 * std::log(r0);
        const double jFin = 2.0 * std::log(r1);
        const double jScale = 2.0 * std::log(r1 / r0); // double-counted (D6/F7)
        const double buggy = (jIni - jFin) + jScale;
        EXPECT_GT(std::abs(buggy - result.lnJac), 0.05)
            << "F7 regression guard: buggy composition (" << buggy << ") must NOT match the "
            << "corrected lnJac (" << result.lnJac << "), s=" << s;
    }
}

// ---------------------------------------------------------------------------
//  V5c -- >=3-body chain (three consecutive BendStretch bonds), scaling an
//  upstream bond AND an upstream angle together: the discriminating case.
// ---------------------------------------------------------------------------
TEST(BatScalingJacobian, ThreeBodyChainMatchesClosedForm) {
    SystemTopology sys = fourAtomBendStretchChainTopology();
    Selection sel;
    sel.bondMobility = {JointType::BendStretch, JointType::BendStretch, JointType::BendStretch};
    World world(0, /*cartesian=*/false, /*seed=*/12u);
    world.buildModel(sys, sel, {JointType::Free});

    const std::vector<Vec3> x0 = fourAtomChainPositions();

    for (const double s : {1.1, 0.9}) {
        const World::BatScalingResult result = world.previewBatScaling(x0, s);
        // body(atom1): root-adjacent (its parent's parent is Ground) -> no
        // zK -> r-only, contributes 1 (BatScaling.hpp's documented
        // limitation). body(atom2)/body(atom3): valid zK -> r+theta, 2 each.
        EXPECT_EQ(result.nScaled, 5) << "s=" << s;

        const double fd = fdLogAbsDetJac(world, x0, s, {}, {}, 1e-6);
        EXPECT_NEAR(fd, result.lnJac, 1.5e-3)
            << "engine FD != analytic on a >=3-body chain (upstream bond+angle) -- V5 GATE, s=" << s;

        // F7 regression guard on the FULL multi-body composition: J_ini +
        // J_scale - J_fin, bond ratios double-counted, angle ratios single
        // (D6: "value-ratios log(r'/r)+log(theta'/theta) over every Q, bond
        // double-counted"). J(x0)/J(x') via the SAME engine function
        // previewBatScaling itself calls internally (robo::
        // calcBatVolumeLogJac) -- reusing it here checks the COMPOSITION,
        // not the per-body element formula (already confirmed by the FD
        // check above).
        const double j0 = robo::calcBatVolumeLogJac(world.model(), x0);
        const double j1 = robo::calcBatVolumeLogJac(world.model(), result.atomPos);

        double jScale = 0.0;
        {
            // body(atom1): r-only, zJ=atom0.
            const double r0 = (x0[1] - x0[0]).norm();
            const double r1 = (result.atomPos[1] - result.atomPos[0]).norm();
            jScale += 2.0 * std::log(r1 / r0);
        }
        {
            // body(atom2): r (zJ=atom1) + theta (zJ=atom1, zK=atom0).
            const double r0 = (x0[2] - x0[1]).norm();
            const double r1 = (result.atomPos[2] - result.atomPos[1]).norm();
            jScale += 2.0 * std::log(r1 / r0);
            const double th0 = angleAt(x0[2], x0[1], x0[0]);
            const double th1 = angleAt(result.atomPos[2], result.atomPos[1], result.atomPos[0]);
            jScale += std::log(th1 / th0);
        }
        {
            // body(atom3): r (zJ=atom2) + theta (zJ=atom2, zK=atom1).
            const double r0 = (x0[3] - x0[2]).norm();
            const double r1 = (result.atomPos[3] - result.atomPos[2]).norm();
            jScale += 2.0 * std::log(r1 / r0);
            const double th0 = angleAt(x0[3], x0[2], x0[1]);
            const double th1 = angleAt(result.atomPos[3], result.atomPos[2], result.atomPos[1]);
            jScale += std::log(th1 / th0);
        }
        const double buggy = (j0 - j1) + jScale;
        EXPECT_GT(std::abs(buggy - result.lnJac), 0.05)
            << "F7 regression guard: buggy composition (" << buggy << ") must NOT match the "
            << "corrected lnJac (" << result.lnJac << "), s=" << s;
    }
}

// ---------------------------------------------------------------------------
//  S2 (reviewer, 2026-07-12) -- an aggressive scale that pushes a scaled
//  angle outside (0, pi) MUST throw std::domain_error, not silently corrupt
//  lnJac. Context::driveReplica (Stage 2b) catches exactly this and converts
//  it to an automatic swap reject (see TestRexAcceptanceAlgebra.cpp's
//  WORK_Jacobian=-infinity consuming-side test).
// ---------------------------------------------------------------------------
TEST(BatScalingJacobian, AggressiveScaleThrowsDomainErrorNotSilentCorruption) {
    SystemTopology sys = fourAtomBendStretchChainTopology();
    Selection sel;
    sel.bondMobility = {JointType::BendStretch, JointType::BendStretch, JointType::BendStretch};
    World world(0, /*cartesian=*/false, /*seed=*/13u);
    world.buildModel(sys, sel, {JointType::Free});

    const std::vector<Vec3> x0 = fourAtomChainPositions();

    // theta0 for body(atom2)/body(atom3) is EXACTLY pi/2 by construction
    // (orthogonal zigzag). With the default (zero) anchor, s=3.0 gives
    // theta1 = 3*(pi/2) = 4.712 rad ~= 270 deg, well outside (0, pi).
    EXPECT_THROW(world.previewBatScaling(x0, 3.0), std::domain_error);

    // Sanity: a MILD scale (s=1.1, already exercised above) does NOT throw
    // -- the guard is a genuine domain check, not an over-eager one.
    EXPECT_NO_THROW(world.previewBatScaling(x0, 1.1));
}
