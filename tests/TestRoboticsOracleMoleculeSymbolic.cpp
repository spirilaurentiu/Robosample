// ============================================================================
//  TestRoboticsOracleMoleculeSymbolic.cpp -- Scope B (molecule) of the
//  robotics oracle (docs/specs/robotics-oracle-differential.md §4.2/§4.4/
//  §5.1), SYMBOLIC-differential family, split out of
//  TestRoboticsOracleMolecule.cpp (TEST-005).
//
//  SCOPE. Proves and locks in the PORT-side half of Scope B's Stage 0/1 --
//  SystemTopology (dumped by the port's OWN context.load_amber, the same path
//  loader_differential uses) -> Context::buildFlexibilities -> World::
//  buildModel -> RobotModel, plus the §5.1 T1/T2 topological preflight -- for
//  the rigid/regular (examples/10ala.*) and cyclic (examples/1APQ.*) topology
//  classes. No live-Simbody numeric comparison here (no `World`-vs-clone state
//  diff) -- that is the NUMERIC family, TestRoboticsOracleMoleculeNumeric.cpp,
//  which also carries the full debugged history (RESOLVED/FINDING/ADJUDICATED
//  notes) behind the tolerances and guards the numeric differential needs.
//
//  §5.1 T1/T2 preflight note: T1/T2 are computed HERE, independently, from
//  the public SystemTopology + RobotModel::atomBody (World::model(),
//  public) rather than from the port's private ConstraintSet
//  (include/World.hpp has no public accessor for it, and the HARD
//  CONSTRAINT forbids adding one to disasm include/). This still fully
//  re-validates the atom->body (DSU) reduction T1/T2 target (§5.1's named
//  "un-guarded gap"), using the IDENTICAL inter-body ring-closing filter
//  formula src/World.cpp:670-680 uses (same inputs: bondsRingClosing +
//  atomBody), but does not additionally catch a bug isolated to that one
//  filter step in isolation from atomBody -- documented, not silently
//  narrowed (Rule 11).
// ============================================================================
#include <algorithm>
#include <vector>
#include <gtest/gtest.h>

#include "Context.hpp"
#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "RoboticsOracleMoleculeLoader.hpp"
#include "TestHelpers.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"

using namespace robo;

namespace {

#ifndef ROBOTICS_ORACLE_MOLECULES_FIXTURE_DIR
#error "ROBOTICS_ORACLE_MOLECULES_FIXTURE_DIR must be defined by the build"
#endif
const std::string kFixtureDir = ROBOTICS_ORACLE_MOLECULES_FIXTURE_DIR;

// Trivial union-find, local to this preflight (mirrors World.cpp's own
// unexported DSU -- not shared, Rule 3: this file's need is tiny and
// World.cpp's DSU is anonymous-namespace/private to that TU).
struct DSU {
    std::vector<int> parent;
    explicit DSU(int n) : parent(static_cast<std::size_t>(n)) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int x) {
        while (parent[static_cast<std::size_t>(x)] != x) {
            parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    }
    void join(int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) {
            parent[static_cast<std::size_t>(a)] = b;
        }
    }
};

// §5.1 T1/T2 topological preflight, run on the BODY-adjacency graph (nodes =
// bodies 1..numBodies-1, Ground excluded -- it is never a bond endpoint) --
// independently rebuilt from `sys.bondsI/J/bondsRingClosing` + `model.
// atomBody` (both public), NOT from `model.bodyParent` (§5.1: bodyParent is
// a forest by construction from the BFS, so checking it can never fail --
// a Rule-8 vacuous test).
struct Preflight {
    int numBodyNodes = 0;   // bodies 1..numBodies-1
    int treeEdges = 0;      // inter-body, non-ring bonds
    int ringInterBodyEdges = 0; // inter-body, ring-closing bonds (T2's L candidates)
    int numComponents = 0;
    bool t1AcyclicOnTreeEdges = false; // treeEdges == numBodyNodes - numComponents
};

Preflight runPreflight(const SystemTopology& sys, const RobotModel& model) {
    Preflight pf;
    pf.numBodyNodes = model.numBodies - 1; // exclude Ground
    DSU dsu(model.numBodies); // indexed 0..numBodies-1; Ground unused

    for (int k = 0; k < sys.numBonds; ++k) {
        const int bi = model.atomBody[sys.bondsI[k]];
        const int bj = model.atomBody[sys.bondsJ[k]];
        if (bi == bj) {
            continue; // intra-body (rigid-merged or ring collapsed within one body)
        }
        if (sys.bondsRingClosing[k]) {
            ++pf.ringInterBodyEdges;
        } else {
            ++pf.treeEdges;
            dsu.join(bi, bj);
        }
    }

    std::vector<bool> seenRoot(static_cast<std::size_t>(model.numBodies), false);
    for (int b = 1; b < model.numBodies; ++b) {
        seenRoot[static_cast<std::size_t>(dsu.find(b))] = true;
    }
    pf.numComponents = static_cast<int>(std::count(seenRoot.begin(), seenRoot.end(), true));
    pf.t1AcyclicOnTreeEdges = (pf.treeEdges == pf.numBodyNodes - pf.numComponents);
    return pf;
}

} // namespace

// ---------------------------------------------------------------------------
//  Rigid class (§4.4 item 1): examples/10ala.*, all bonds Rigid (no explicit
//  flexibilities => Context::addRoboticWorld(Selection{}) defaults every
//  bond to JointType::Rigid, src/World.cpp:380/413) -- the whole 112-atom
//  peptide collapses to ONE rigid body. §8 fail-loud assert: numBodies==2
//  (Ground + the one rigid body).
// ---------------------------------------------------------------------------
TEST(RoboticsOracleMolecule, RigidWeldRoot) {
    const SystemTopology sys = robotics_oracle_loader::loadPortSystemTopology(kFixtureDir, "10ala");
    Context ctx("RoboticsOracleMolecule_RigidWeldRoot", 0);
    ctx.systemTopology = sys;
    // Port default root mobility (python/robosample/context.py) is Rigid
    // (Weld) for a no-box system -- verify the fixture actually carries that,
    // not just assume it (Rule 11).
    ASSERT_EQ(sys.numMolecules, 1);
    ASSERT_EQ(static_cast<int>(sys.rootMobilities[0]), static_cast<int>(JointType::Rigid))
        << "fixture's default root mobility changed -- Weld-root assumption stale";

    World& world = ctx.addRoboticWorld(Selection{});
    const RobotModel& m = world.model();

    // §8 mandated fail-loud assert: rigid 10ala -> numBodies==2.
    EXPECT_EQ(m.numBodies, 2) << "rigid 10ala must collapse to Ground + one rigid body";
    EXPECT_EQ(m.nu, 0) << "Weld root: 0 total DOF";
    EXPECT_EQ(m.nq, 0);
}

// Free-root variant of the rigid class (§4.4 item 1): same rigid-body
// decomposition, but the root attaches to Ground with 6 DOF instead of 0.
TEST(RoboticsOracleMolecule, RigidFreeRoot) {
    const SystemTopology sys = robotics_oracle_loader::loadPortSystemTopology(kFixtureDir, "10ala");
    Context ctx("RoboticsOracleMolecule_RigidFreeRoot", 0);
    ctx.systemTopology = sys;
    ctx.systemTopology.rootMobilities[0] = JointType::Free;

    World& world = ctx.addRoboticWorld(Selection{});
    const RobotModel& m = world.model();

    EXPECT_EQ(m.numBodies, 2);
    EXPECT_EQ(m.nu, 6) << "Free root: 6 total DOF (classic free-rigid-body case)";
    EXPECT_EQ(m.nq, 7);
}

// ---------------------------------------------------------------------------
//  Regular class (§4.4 item 2): examples/10ala.*, default-flexible mobility
//  (Context::buildFlexibilities(nullopt, Torsion, false) -- every eligible
//  non-terminal, non-ring bond becomes Torsion; ring-closing bonds always
//  stay Rigid and are never tree edges). §8 fail-loud assert: 0 active
//  ring-closure constraints (10ala/alanine has no ring side chain).
// ---------------------------------------------------------------------------
TEST(RoboticsOracleMolecule, RegularFlexible) {
    const SystemTopology sys = robotics_oracle_loader::loadPortSystemTopology(kFixtureDir, "10ala");
    ASSERT_EQ(std::count(sys.bondsRingClosing.begin(), sys.bondsRingClosing.end(), true), 0)
        << "10ala has no ring side chain -- atom-graph ring count must be 0 (deliberate control, §8)";

    Context ctx("RoboticsOracleMolecule_RegularFlexible", 0);
    ctx.systemTopology = sys;
    const Selection sel = ctx.buildFlexibilities(std::nullopt, JointType::Torsion, false);

    World& world = ctx.addRoboticWorld(sel);
    const RobotModel& m = world.model();

    EXPECT_GT(m.numBodies, 2) << "flexible selection must produce more than one rigid body";
    EXPECT_GT(m.nu, 0) << "flexible selection must produce nonzero internal DOF";

    const Preflight pf = runPreflight(sys, m);
    EXPECT_TRUE(pf.t1AcyclicOnTreeEdges)
        << "T1: body-tree (from bond list + atomBody, NOT bodyParent) must be a spanning forest -- "
           "treeEdges=" << pf.treeEdges << " nodes=" << pf.numBodyNodes << " components=" << pf.numComponents;
    // §8 mandated fail-loud assert: regular 10ala -> 0 active ring-closure
    // constraints (T2: L == 0, no ring side chain to close).
    EXPECT_EQ(pf.ringInterBodyEdges, 0) << "T2: regular 10ala must have 0 active (inter-body) ring closures";
}

// §5.1 T1 discriminating regression: inject a real ring bond with
// bondsRingClosing FALSELY set to false (simulating a mis-reduced loop-
// closure that leaked into the tree-edge filter) and assert T1 fails. This
// is what makes T1 non-vacuous (§5.1: "a loop-closing bond mis-reduced into
// a tree edge... must fail loud").
TEST(RoboticsOracleMolecule, T1FailsOnMisflaggedRingBond) {
    SystemTopology sys = robotics_oracle_loader::loadPortSystemTopology(kFixtureDir, "1APQ");
    ASSERT_GT(std::count(sys.bondsRingClosing.begin(), sys.bondsRingClosing.end(), true), 0)
        << "1APQ must carry real ring-closing bonds for this regression to be meaningful";

    Context ctx("RoboticsOracleMolecule_T1Regression", 0);
    ctx.systemTopology = sys;
    const Selection sel = ctx.buildFlexibilities(std::nullopt, JointType::Torsion, false);
    World& world = ctx.addRoboticWorld(sel);
    const RobotModel& m = world.model();

    // Sanity: the UNMODIFIED topology's preflight is clean (T1 holds).
    ASSERT_TRUE(runPreflight(sys, m).t1AcyclicOnTreeEdges);

    // Flip the FIRST inter-body ring-closing bond's flag to false: the
    // preflight now sees it as a tree edge alongside the rest of the real
    // spanning tree it closes a loop on top of -- a cycle.
    bool flipped = false;
    for (int k = 0; k < sys.numBonds && !flipped; ++k) {
        if (sys.bondsRingClosing[k] && m.atomBody[sys.bondsI[k]] != m.atomBody[sys.bondsJ[k]]) {
            sys.bondsRingClosing[k] = false;
            flipped = true;
        }
    }
    ASSERT_TRUE(flipped) << "1APQ must have >=1 inter-body ring-closing bond to flip";

    const Preflight pfBroken = runPreflight(sys, m);
    EXPECT_FALSE(pfBroken.t1AcyclicOnTreeEdges)
        << "T1 must fail loud when a real ring-closing bond is mis-flagged as a tree edge";
}

// ---------------------------------------------------------------------------
//  Cyclic class (§4.4 item 3): examples/1APQ.* (prolines + disulfide CYX
//  cross-links), default-flexible mobility (§4.2: "default for cyclic").
//  §8 mandated report: BOTH the atom-graph #bondsRingClosing AND the
//  body-graph active ring-closure count (T2's L) -- they may differ (§5.1
//  T2: a ring collapsed entirely inside one rigid body contributes to the
//  atom-graph cyclomatic number but no body-graph edge), so this test
//  asserts ONLY #bondsRingClosing >= L (never equality) plus the §8 minimum
//  L >= 1 (1APQ's disulfides are inter-body regardless of mobility).
// ---------------------------------------------------------------------------
TEST(RoboticsOracleMolecule, Cyclic1APQ) {
    const SystemTopology sys = robotics_oracle_loader::loadPortSystemTopology(kFixtureDir, "1APQ");
    const int bondsRingClosingCount =
        static_cast<int>(std::count(sys.bondsRingClosing.begin(), sys.bondsRingClosing.end(), true));
    // Sanity check on the fixture itself: for a single connected molecule,
    // atom-graph cyclomatic number == numBonds - numAtoms + numMolecules.
    ASSERT_EQ(bondsRingClosingCount, sys.numBonds - sys.numAtoms + sys.numMolecules)
        << "1APQ atom-graph ring count vs Euler cyclomatic-number cross-check (engine-independent)";

    Context ctx("RoboticsOracleMolecule_Cyclic1APQ", 0);
    ctx.systemTopology = sys;
    const Selection sel = ctx.buildFlexibilities(std::nullopt, JointType::Torsion, false);
    World& world = ctx.addRoboticWorld(sel);
    const RobotModel& m = world.model();

    const Preflight pf = runPreflight(sys, m);
    EXPECT_TRUE(pf.t1AcyclicOnTreeEdges) << "T1: 1APQ's spanning tree (bonds minus ring-closing) must be acyclic";

    // §8 mandated fail-loud assert: 1APQ -> >=1 active constraint.
    EXPECT_GE(pf.ringInterBodyEdges, 1) << "T2: 1APQ must have >=1 active (inter-body) ring closure";
    // §5.1 T2: never assert equality to the raw atom-graph count -- only
    // #bondsRingClosing >= L (a ring collapsed intra-body contributes to the
    // atom-graph count but not the body graph).
    EXPECT_GE(bondsRingClosingCount, pf.ringInterBodyEdges)
        << "T2: #bondsRingClosing must be >= the body-graph active ring-closure count";

    // Report both counts (§8: "report both #bondsRingClosing and
    // numConstraints() for 1APQ").
    std::cout << "[ RoboticsOracleMolecule.Cyclic1APQ ] #bondsRingClosing=" << bondsRingClosingCount
              << " bodyGraphActiveRingClosures(L)=" << pf.ringInterBodyEdges << " numBodies=" << m.numBodies
              << " nu=" << m.nu << std::endl;
}
