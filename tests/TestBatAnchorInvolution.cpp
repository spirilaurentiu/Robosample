// ============================================================================
//  TestBatAnchorInvolution.cpp -- V10 (anchor involution), guards INV-9/D2
//  (docs/specs/replica-exchange-nonequilibrium-work.md).
//
//  CLAIM (D2, algebraic): the affine BAT-scaling map q' = s*q - (s-1)*mu is
//  invertible with M_{1/s,mu} = (M_{s,mu})^-1 EXACTLY when both legs share
//  the SAME mu:
//      M_{1/s,mu}(M_{s,mu}(q)) = (1/s)*(s*q-(s-1)*mu) - ((1/s)-1)*mu
//                              = q - ((s-1)/s)*mu + ((s-1)/s)*mu = q.
//  This cancellation is INDEPENDENT of s and holds per-DOF; by induction over
//  the kinematic tree (an ancestor body's OWN involution restores its
//  descendants' z-matrix reference atoms exactly, so a descendant's
//  RELATIVE (r, theta) round-trips too), it holds for the FULL multi-body
//  map robo::applyBatScaling (World::previewBatScaling) as well -- not just
//  a single isolated body. Using a DIFFERENT mu on the reverse leg (the
//  per-thermodynamic-state running mean the original used, D2's REFUTED
//  design) breaks the cancellation and the round-trip fails to close.
//
//  Calls the ENGINE directly (World::buildModel + World::previewBatScaling),
//  the SAME >=3-body BendStretch chain fixture as TestBatScalingJacobian.cpp
//  (kept independent/duplicated per this test suite's established
//  per-file-fixture convention, e.g. TestTwoRobotContact.cpp's own header
//  comment).
//
//  NOT compiled or run (coordinator directive, 2026-07-12) -- written to the
//  tests/Test*.cpp convention for the user to build and run manually.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <unordered_map>
#include <vector>

#include "RobotModel.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"
#include "robot_math.hpp"

using robo::Vec3;

namespace {

// Orthogonal zigzag fixture (see TestBatScalingJacobian.cpp for the by-hand
// geometric verification: angle(atom2,atom1,atom0) == angle(atom3,atom2,
// atom1) == exactly pi/2).
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

double maxAbsDiff(const std::vector<Vec3>& a, const std::vector<Vec3>& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (int d = 0; d < 3; ++d) {
            m = std::max(m, std::abs(a[i][d] - b[i][d]));
        }
    }
    return m;
}

} // namespace

TEST(BatAnchorInvolution, SharedAnchorRoundTrips) {
    SystemTopology sys = fourAtomBendStretchChainTopology();
    Selection sel;
    sel.bondMobility = {JointType::BendStretch, JointType::BendStretch, JointType::BendStretch};
    World world(0, /*cartesian=*/false, /*seed=*/21u);
    world.buildModel(sys, sel, {JointType::Free});

    const std::vector<Vec3> x0 = fourAtomChainPositions();

    const double s = 1.15;
    // Nontrivial, SHARED anchor (INV-9): away from the current (r0, theta0)
    // values (0.15 nm, pi/2 rad) so the round-trip is a real algebraic
    // cancellation, not a degenerate mu==q no-op. Kept modest so theta1
    // stays safely inside (0, pi) for s=1.15 (theta0=pi/2=1.5708,
    // theta1 = 1.15*1.5708 - 0.15*0.10 = 1.7914 rad ~= 102.7 deg, and the
    // reverse leg's intermediate values are algebraically bounded the same
    // way -- see the module docstring's derivation).
    const std::unordered_map<int, double> anchorR = {{1, 0.02}, {2, 0.03}, {3, 0.04}};
    const std::unordered_map<int, double> anchorTheta = {{2, 0.10}, {3, -0.05}};

    const World::BatScalingResult forward = world.previewBatScaling(x0, s, anchorR, anchorTheta);
    const World::BatScalingResult back = world.previewBatScaling(forward.atomPos, 1.0 / s, anchorR, anchorTheta);

    const double maxErr = maxAbsDiff(back.atomPos, x0);
    EXPECT_LT(maxErr, 1e-9) << "M_{1/s,mu}(M_{s,mu}(q)) != q for the SHARED anchor (max err " << maxErr << ")";
}

TEST(BatAnchorInvolution, DifferentPerPartnerAnchorBreaksRoundTrip) {
    SystemTopology sys = fourAtomBendStretchChainTopology();
    Selection sel;
    sel.bondMobility = {JointType::BendStretch, JointType::BendStretch, JointType::BendStretch};
    World world(0, /*cartesian=*/false, /*seed=*/22u);
    world.buildModel(sys, sel, {JointType::Free});

    const std::vector<Vec3> x0 = fourAtomChainPositions();

    const double s = 1.15;
    const std::unordered_map<int, double> anchorR = {{1, 0.02}, {2, 0.03}, {3, 0.04}};
    const std::unordered_map<int, double> anchorTheta = {{2, 0.10}, {3, -0.05}};

    const World::BatScalingResult forward = world.previewBatScaling(x0, s, anchorR, anchorTheta);

    // The D2-REFUTED design: the reverse leg uses a DIFFERENT (e.g.
    // per-thermodynamic-state) anchor instead of the shared, frozen one
    // (INV-9). The round-trip SHALL break -- this is the discriminating
    // failure V10 exists to catch (a future regression to per-state means
    // would otherwise pass silently). Per D2's own algebra,
    // M_{1/s,mu'}(M_{s,mu}(q)) = q + ((s-1)/s)*(mu'-mu); with s=1.15,
    // (s-1)/s ~= 0.13, so a delta=0.05 perturbation gives a residual of
    // order 0.13*0.05 ~= 6.5e-3 per affected body -- comfortably above the
    // 1e-3 threshold below even allowing for sign cancellation across the
    // 3-body chain.
    std::unordered_map<int, double> anchorROther = anchorR;
    std::unordered_map<int, double> anchorThetaOther = anchorTheta;
    for (auto& kv : anchorROther) {
        kv.second += 0.05;
    }
    for (auto& kv : anchorThetaOther) {
        kv.second += 0.05;
    }

    const World::BatScalingResult backWrong =
        world.previewBatScaling(forward.atomPos, 1.0 / s, anchorROther, anchorThetaOther);

    const double maxErrWrong = maxAbsDiff(backWrong.atomPos, x0);
    EXPECT_GT(maxErrWrong, 1e-3) << "using a DIFFERENT anchor on the reverse leg SHOULD break the "
                                    "round-trip (INV-9) -- got a suspiciously small error, the "
                                    "anchor may not be wired into the map correctly";
}
