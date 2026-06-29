// ============================================================================
//  TestGibbsWorlds.cpp -- does the internal<->Cartesian transfer (the inter-world
//  hand-off) move the measure?
//
//  In a Gibbs-across-worlds scheme a configuration is passed between a torsional
//  world (internal q) and a Cartesian world (atom positions). The hand-off must be
//  a PURE function of the configuration: realizing geometry from q must carry no
//  hidden state and must be idempotent, otherwise repeated hand-offs would drift
//  the distribution (the recomputeGeometry concern called out in World.cpp).
//
//  This file pins the engine-level transfer:
//   1. IDEMPOTENT / STATELESS: realizePosition + fillAtomPositionsFromBodies, run
//      twice from the same q, produces BITWISE-identical atom positions and body
//      transforms X_GB. A second realization adds nothing -- no accumulation.
//   2. DETERMINISTIC ACROSS INSTANCES: two independently allocated states with the
//      same q realize to bitwise-identical geometry (no per-instance state leaks
//      into the map).
//
//  REPORTED LIMITATION. The full cross-world Gibbs invariance -- a known 1-D
//  marginal preserved under [torsional o Cartesian] block alternation -- lives in
//  World::recomputeGeometry / setAtomsLocationsInGround (the per-block root-frame
//  reset to identity with q=0). That path is World-only (it pulls OpenMM via
//  ForceBridge) and cannot be exercised at engine level; it needs a World-level
//  harness. The idempotence proven here is the necessary engine-level precondition
//  for that invariance.
// ============================================================================
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

RobotModel freeTorsionChain(Rng& rng) {
    auto mk = [&](int parent, JointType jt) {
        BodySpec b;
        b.parent = parent;
        b.joint = jt;
        b.X_PF = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        b.X_BM = Transform(rng.rotation(), rng.vec3(-0.2, 0.2));
        b.mass = rng.uniform(0.9, 1.5);
        b.com_B = rng.vec3(-0.08, 0.08);
        b.inertia_B = UnitInertia(rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6), rng.uniform(0.4, 0.6));
        return b;
    };
    RobotModel m = buildForest({mk(0, JointType::Free), mk(1, JointType::Torsion), mk(2, JointType::Torsion)});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, 0.02, -0.03), Vec3(-0.02, 0.10, 0.01)}, {12.0, 1.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0.04, -0.06, 0.09), Vec3(0.08, 0.01, -0.02)}, {16.0, 1.0});
    return m;
}

void realizeGeometry(const RobotModel& m, RobotState& s) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
}

} // namespace

// ---------------------------------------------------------------------------
//  The geometry map is idempotent: realizing twice from the same q leaves atom
//  positions and body transforms bitwise unchanged (no accumulated state).
// ---------------------------------------------------------------------------
TEST(GibbsWorlds, TransferIsIdempotent) {
    Rng rng(0x6101);
    RobotModel m = freeTorsionChain(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);

    realizeGeometry(m, s);
    std::vector<Vec3> pos1(s.atomPosG(), s.atomPosG() + m.numAtoms);
    std::vector<Transform> xgb1(s.X_GB(), s.X_GB() + m.numBodies);

    realizeGeometry(m, s); // second pass from the SAME q
    for (int a = 0; a < m.numAtoms; ++a) {
        EXPECT_TRUE(rtest::NearVec3(s.atomPosG()[a], pos1[static_cast<std::size_t>(a)], 0.0))
            << "atom " << a << " moved on a second realization (transfer not idempotent)";
    }
    for (int b = 0; b < m.numBodies; ++b) {
        EXPECT_TRUE(rtest::NearVec3(s.X_GB()[b].p(), xgb1[static_cast<std::size_t>(b)].p(), 0.0))
            << "body " << b << " origin moved on a second realization";
        EXPECT_TRUE(rtest::NearMat33(s.X_GB()[b].R(), xgb1[static_cast<std::size_t>(b)].R(), 0.0))
            << "body " << b << " rotation moved on a second realization";
    }
}

// ---------------------------------------------------------------------------
//  The map is a pure function of q: a second, independently allocated state with
//  the same q realizes to bitwise-identical geometry (no per-instance state leak).
// ---------------------------------------------------------------------------
TEST(GibbsWorlds, TransferIsDeterministicAcrossInstances) {
    Rng rng(0x6202);
    RobotModel m = freeTorsionChain(rng);
    RobotState s1;
    s1.allocateFull(m);
    randomizeState(m, s1, rng);

    RobotState s2;
    s2.allocateFull(m);
    std::copy(s1.q(), s1.q() + m.nq, s2.q());
    std::copy(s1.u(), s1.u() + m.nu, s2.u());

    realizeGeometry(m, s1);
    realizeGeometry(m, s2);
    for (int a = 0; a < m.numAtoms; ++a) {
        EXPECT_TRUE(rtest::NearVec3(s1.atomPosG()[a], s2.atomPosG()[a], 0.0))
            << "atom " << a << " differs between two states with identical q";
    }
}
