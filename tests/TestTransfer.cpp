// ============================================================================
//  TestTransfer.cpp -- Phase 6: the internal-q -> Cartesian transfer, the
//  currency every downstream consumer reads.
//
//  realizePosition fills, per atom,
//      atomStationG[a] = R_GB * station_B[a]
//      atomPosG[a]     = X_GB[b].p() + atomStationG[a]
//  and fillAtomPositionsFromBodies recomputes atomPosG[a] = X_GB[b].p() +
//  X_GB[b].R() * station_B[a] from already-realized body transforms. These arrays
//  are what the forcefield, the DCD writer, the statistics layer, and SHAKE all
//  consume, yet nothing pins them directly. This file does.
//
//  Fixtures:
//   * freeTree  -- a SINGLE Free-rooted tree (Free -> Torsion -> Ball). Every
//                  body descends from the Free root, so a rigid motion applied at
//                  the root moves the WHOLE robot -- the equivariance fixture.
//   * weldForest -- a Free-rooted chain PLUS an independent Rigid (welded-to-
//                  Ground) robot, for the rigid-within-body and frozen-weld checks.
// ============================================================================
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"
#include "engine_helpers.hpp" // rotationToQuaternion (build a root q' for the rigid map)

using namespace robo;
using rtest::attachAtoms;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::randomizeState;
using rtest::Rng;

namespace {

// A single Free-rooted tree: every moving body hangs off body 1, so a rigid map
// applied at the root transports the entire robot.
RobotModel freeTree(Rng& rng) {
    auto F = [&](int parent, JointType jt) {
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
    RobotModel m = buildForest({F(0, JointType::Free), F(1, JointType::Torsion), F(2, JointType::Ball)});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, -0.02, 0.04), Vec3(-0.05, 0.07, 0.0)}, {12.0, 1.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.10, 0), Vec3(-0.03, 0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0, 0.07, -0.05), Vec3(0.04, 0, 0.05)}, {16.0, 1.0});
    return m;
}

// A Free-rooted chain plus an independent Rigid (welded-to-Ground) robot.
// Bodies: 1=Free, 2=Torsion (descendants); 3=Rigid welded straight to Ground.
RobotModel weldForest(Rng& rng) {
    auto F = [&](int parent, JointType jt) {
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
    RobotModel m = buildForest({F(0, JointType::Free), F(1, JointType::Torsion), F(0, JointType::Rigid)});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, -0.02, 0.04)}, {12.0, 1.0});     // atoms 0,1
    attachAtoms(m, 2, {Vec3(0.05, 0.10, 0), Vec3(-0.03, 0, 0.06)}, {14.0, 1.0});  // atoms 2,3
    attachAtoms(m, 3, {Vec3(0.02, 0.0, 0.0), Vec3(0.0, 0.06, 0.0)}, {15.0, 1.0}); // atoms 4,5 (welded)
    return m;
}

} // namespace

// ---------------------------------------------------------------------------
//  1. StationFormulaExact: after realizePosition, both filled arrays match the
//     defining formula exactly (kTight).
// ---------------------------------------------------------------------------
TEST(Transfer, StationFormulaExact) {
    Rng rng(0x6101);
    RobotModel m = freeTree(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 10; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);

        const Transform* X_GB = s.X_GB();
        const Vec3* posG = s.atomPosG();
        const Vec3* stG = s.atomStationG();
        for (int a = 0; a < m.numAtoms; ++a) {
            const int b = m.atomBody[a];
            const Vec3 station = X_GB[b].R() * m.atomStation_B[a];
            EXPECT_TRUE(rtest::NearVec3(stG[a], station, rtest::kTight))
                << "atomStationG wrong, atom " << a << " rep " << rep;
            EXPECT_TRUE(rtest::NearVec3(posG[a], X_GB[b].p() + station, rtest::kTight))
                << "atomPosG wrong, atom " << a << " rep " << rep;
        }
    }
}

// ---------------------------------------------------------------------------
//  2. FillAtomPositionsRoundTrips: with body transforms unchanged,
//     fillAtomPositionsFromBodies reproduces the realizePosition atom fill
//     BIT-FOR-BIT (the two share the X_GB.p + X_GB.R*station expression).
// ---------------------------------------------------------------------------
TEST(Transfer, FillAtomPositionsRoundTrips) {
    Rng rng(0x6102);
    RobotModel m = freeTree(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 10; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);

        // snapshot the realizePosition fill
        std::vector<Vec3> ref(s.atomPosG(), s.atomPosG() + m.numAtoms);

        // clobber, then refill from the (unchanged) body transforms
        for (int a = 0; a < m.numAtoms; ++a) {
            s.atomPosG()[a] = Vec3(-1e9, -1e9, -1e9);
        }
        RobotEngine::fillAtomPositionsFromBodies(m, s);

        const Vec3* posG = s.atomPosG();
        for (int a = 0; a < m.numAtoms; ++a) {
            // bit-for-bit: the two code paths must produce identical doubles
            EXPECT_EQ(posG[a][0], ref[static_cast<std::size_t>(a)][0]) << "atom " << a << " x, rep " << rep;
            EXPECT_EQ(posG[a][1], ref[static_cast<std::size_t>(a)][1]) << "atom " << a << " y, rep " << rep;
            EXPECT_EQ(posG[a][2], ref[static_cast<std::size_t>(a)][2]) << "atom " << a << " z, rep " << rep;
        }
    }
}

// ---------------------------------------------------------------------------
//  3. AtomsRigidWithinBody: two atoms on the same body keep a constant distance
//     under arbitrary state (kTight). Couples with RigidRootIsFrozen: the two
//     atoms on the welded Rigid root never move at all.
// ---------------------------------------------------------------------------
TEST(Transfer, AtomsRigidWithinBody) {
    Rng rng(0x6103);
    RobotModel m = weldForest(rng);
    RobotState s;
    s.allocateFull(m);

    // atom layout: body1 -> {0,1}, body2 -> {2,3}, body3(welded) -> {4,5}
    const int weldBeg = 4, weldEnd = 6;

    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    // reference intra-body distances and the frozen weld positions
    const Real d01 = (s.atomPosG()[0] - s.atomPosG()[1]).norm(); // same body (1)
    const Real d23 = (s.atomPosG()[2] - s.atomPosG()[3]).norm(); // same body (2)
    const Real d45 = (s.atomPosG()[4] - s.atomPosG()[5]).norm(); // same body (3, welded)
    std::vector<Vec3> weldRef(s.atomPosG() + weldBeg, s.atomPosG() + weldEnd);

    for (int rep = 0; rep < 25; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);

        EXPECT_NEAR((s.atomPosG()[0] - s.atomPosG()[1]).norm(), d01, rtest::kTight)
            << "body1 not rigid, rep " << rep;
        EXPECT_NEAR((s.atomPosG()[2] - s.atomPosG()[3]).norm(), d23, rtest::kTight)
            << "body2 not rigid, rep " << rep;
        EXPECT_NEAR((s.atomPosG()[4] - s.atomPosG()[5]).norm(), d45, rtest::kTight)
            << "weld not rigid, rep " << rep;

        // the welded-to-Ground body has zero dof: its atoms are frozen for ALL state.
        for (int a = weldBeg; a < weldEnd; ++a) {
            EXPECT_TRUE(rtest::NearVec3(s.atomPosG()[a],
                                        weldRef[static_cast<std::size_t>(a - weldBeg)],
                                        rtest::kTight))
                << "welded-root atom " << a << " moved, rep " << rep;
        }
    }
}

// ---------------------------------------------------------------------------
//  4. TransferIsRigidMotionEquivariant: apply a random rigid map g to the whole
//     Free-rooted robot (by composing g onto the root's across-joint transform),
//     and every atom position transforms by exactly g (kAlg).
//
//     Construction: X_GB[1] = X_PF[1] * X_FM(q_root) * ~X_BM[1]. To get
//     X_GB'[1] = g * X_GB[1] (which propagates g to every descendant since
//     X_GB[child] = X_GB[1] * (joint-only chain)), solve
//         X_FM(q') = ~X_PF[1] * g * X_PF[1] * X_FM(q_root)
//     and write q' = (quat(X_FM'.R), X_FM'.p) back into the Free q-block.
// ---------------------------------------------------------------------------
TEST(Transfer, TransferIsRigidMotionEquivariant) {
    Rng rng(0x6104);
    RobotModel m = freeTree(rng);
    RobotState s;
    s.allocateFull(m);

    for (int rep = 0; rep < 20; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        std::vector<Vec3> before(s.atomPosG(), s.atomPosG() + m.numAtoms);

        // random rigid map g
        const Transform g(rng.rotation(), rng.vec3(-1.0, 1.0));

        // X_FM'(root) = ~X_PF[1] * g * X_PF[1] * X_FM(root)
        const Transform XFMnew = (~m.X_PF[1]) * g * m.X_PF[1] * s.X_FM()[1];
        Real qw, qx, qy, qz;
        EngineHelpers::rotationToQuaternion(XFMnew.R(), qw, qx, qy, qz);
        const int qOff = m.bodyQIndex[1];
        s.q()[qOff + 0] = qw;
        s.q()[qOff + 1] = qx;
        s.q()[qOff + 2] = qy;
        s.q()[qOff + 3] = qz;
        s.q()[qOff + 4] = XFMnew.p()[0];
        s.q()[qOff + 5] = XFMnew.p()[1];
        s.q()[qOff + 6] = XFMnew.p()[2];

        RobotEngine::realizePosition(m, s);

        const Vec3* posG = s.atomPosG();
        for (int a = 0; a < m.numAtoms; ++a) {
            EXPECT_TRUE(rtest::NearVec3(posG[a], g * before[static_cast<std::size_t>(a)], rtest::kAlg))
                << "atom " << a << " not equivariant under the root rigid map, rep " << rep;
        }
    }
}