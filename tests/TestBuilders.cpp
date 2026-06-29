// ============================================================================
//  TestBuilders.cpp -- Phase 0.2: the atom-attachment builder extension.
//
//  buildForest hands back atom-free robots (numAtoms == 0). Constraints
//  (SHAKE/RATTLE), the force bridge, and the internal->Cartesian transfer all
//  read model.atomBody / atomStation_B / atomPosG, so attachAtoms (RobotBuilders
//  .hpp) is the precondition every one of those downstream test files relies on.
//
//  This file pins attachAtoms two ways:
//   * PURE (no engine): the arrays and the body->atoms CSR it builds are
//     internally consistent (every atom appears exactly once, under its own
//     body; the root atom is set; virtual sites are stored verbatim).
//   * ENGINE: after RobotState::allocateFull + realizePosition, the per-atom
//     Ground positions obey the station formula
//         atomPosG[a] == X_GB[b].p() + X_GB[b].R() * atomStation_B[a]
//     and atomStationG[a] == R_GB * atomStation_B[a], to machine precision.
//     This is the kTight golden the Phase 3/6 constraint+transfer tests assume.
//
//  NOT PORTED / DEFERRED: the z-matrix-driven atom placement (model.zI/zJ/zK/zL)
//  is the production path; these hand-built robots set stations directly, which
//  is what the multibody/constraint tests need. The z-matrix transfer is a
//  separate (OpenMM-pipeline) concern and is not exercised here.
// ============================================================================
#include <array>
#include <gtest/gtest.h>
#include <vector>

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

// Ground -> Free(body 1) -> Torsion(body 2). Random static frames.
RobotModel twoBody(Rng& rng) {
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
    return buildForest({mk(0, JointType::Free), mk(1, JointType::Torsion)});
}

} // namespace

// ---------------------------------------------------------------------------
//  P1 (pure): attachAtoms populates the per-atom arrays and sets the root atom.
// ---------------------------------------------------------------------------
TEST(Builders, AttachAtomsPopulatesArrays) {
    Rng rng(0xA70);
    RobotModel m = twoBody(rng);
    EXPECT_EQ(m.numAtoms, 0); // buildForest is atom-free

    const std::vector<Vec3> st1 = {Vec3(0, 0, 0), Vec3(0.1, 0, 0)};
    const std::vector<Real> ms1 = {12.0, 1.0};
    const std::vector<Vec3> st2 = {Vec3(0, 0.15, 0)};
    const std::vector<Real> ms2 = {16.0};

    attachAtoms(m, 1, st1, ms1);
    attachAtoms(m, 2, st2, ms2);

    ASSERT_EQ(m.numAtoms, 3);
    // atom 0,1 -> body 1 ; atom 2 -> body 2 (append order)
    EXPECT_EQ(m.atomBody[0], 1);
    EXPECT_EQ(m.atomBody[1], 1);
    EXPECT_EQ(m.atomBody[2], 2);
    EXPECT_TRUE(rtest::NearVec3(m.atomStation_B[0], st1[0], rtest::kTight));
    EXPECT_TRUE(rtest::NearVec3(m.atomStation_B[1], st1[1], rtest::kTight));
    EXPECT_TRUE(rtest::NearVec3(m.atomStation_B[2], st2[0], rtest::kTight));
    EXPECT_EQ(m.atomMass[0], 12.0);
    EXPECT_EQ(m.atomMass[2], 16.0);
    // first atom on each body becomes that body's root atom
    EXPECT_EQ(m.bodyRootAtom[1], 0);
    EXPECT_EQ(m.bodyRootAtom[2], 2);
    // compound index == array position (the atom index identity)
    EXPECT_EQ(m.atomCompoundIndex[2], 2);
}

// ---------------------------------------------------------------------------
//  P2 (pure): the body->atoms CSR partitions every atom exactly once, each
//  under its own body. (Constraint/transfer tests walk a body's atoms via CSR.)
// ---------------------------------------------------------------------------
TEST(Builders, AtomCsrPartitionsAllAtoms) {
    Rng rng(0xA71);
    RobotModel m = twoBody(rng);
    attachAtoms(m, 2, {Vec3(0, 0, 0), Vec3(0.1, 0.1, 0)}, {14.0, 1.0});
    attachAtoms(m, 1, {Vec3(0, 0, 0)}, {12.0}); // attach to a DIFFERENT body second

    // CSR ranges are non-overlapping, cover [0,numAtoms), and every listed atom
    // sits under the body whose range contains it.
    std::vector<int> seen(static_cast<std::size_t>(m.numAtoms), 0);
    int total = 0;
    for (int b = 0; b < m.numBodies; ++b) {
        ASSERT_LE(m.bodyAtomsBeg[b], m.bodyAtomsEnd[b]);
        for (int i = m.bodyAtomsBeg[b]; i < m.bodyAtomsEnd[b]; ++i) {
            const int a = m.bodyAtoms[static_cast<std::size_t>(i)];
            EXPECT_EQ(m.atomBody[a], b) << "atom " << a << " listed under wrong body";
            seen[static_cast<std::size_t>(a)]++;
            total++;
        }
    }
    EXPECT_EQ(total, m.numAtoms);
    for (int a = 0; a < m.numAtoms; ++a) {
        EXPECT_EQ(seen[static_cast<std::size_t>(a)], 1) << "atom " << a << " not covered exactly once";
    }
}

// ---------------------------------------------------------------------------
//  P3 (pure): a virtual site (mass == 0) is stored verbatim, not dropped. The
//  force reduction skips it, but the model must still carry it (it has a station
//  and a body), so a regression that filters mass==0 at attach time is caught.
// ---------------------------------------------------------------------------
TEST(Builders, VirtualSiteMassZeroPreserved) {
    Rng rng(0xA72);
    RobotModel m = twoBody(rng);
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.05, 0, 0)}, {16.0, 0.0}); // 2nd is an M-site
    ASSERT_EQ(m.numAtoms, 2);
    EXPECT_EQ(m.atomMass[1], 0.0);
    EXPECT_TRUE(rtest::NearVec3(m.atomStation_B[1], Vec3(0.05, 0, 0), rtest::kTight));
    EXPECT_EQ(m.atomBody[1], 1);
}

// ---------------------------------------------------------------------------
//  P4 (engine): the station formula golden -- the kTight invariant every
//  downstream constraint/transfer test depends on.
//      atomPosG[a]     == X_GB[b].p() + X_GB[b].R() * atomStation_B[a]
//      atomStationG[a] == X_GB[b].R() * atomStation_B[a]
//  attachAtoms MUST precede allocateFull (the state sizes atom arrays from
//  model.numAtoms); this test pins that ordering contract by construction.
// ---------------------------------------------------------------------------
TEST(Builders, StationFormulaExactAfterRealize) {
    Rng rng(0xA73);
    RobotModel m = twoBody(rng);
    attachAtoms(m, 1, {Vec3(0.0, 0.0, 0.0), Vec3(0.12, -0.03, 0.05)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.0, 0.10, 0.0), Vec3(-0.04, 0.0, 0.07)}, {16.0, 1.0});

    RobotState s;
    s.allocateFull(m); // AFTER attachAtoms -> atom arrays sized for numAtoms
    for (int rep = 0; rep < 20; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        const Transform* X_GB = s.X_GB();
        const Vec3* posG = s.atomPosG();
        const Vec3* stG = s.atomStationG();
        for (int a = 0; a < m.numAtoms; ++a) {
            const int b = m.atomBody[a];
            const Vec3 st = X_GB[b].R() * m.atomStation_B[a];
            EXPECT_TRUE(rtest::NearVec3(stG[a], st, rtest::kTight)) << "stationG atom " << a;
            EXPECT_TRUE(rtest::NearVec3(posG[a], X_GB[b].p() + st, rtest::kTight))
                << "posG atom " << a << " rep " << rep;
        }
    }
}