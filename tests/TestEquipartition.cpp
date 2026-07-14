// ============================================================================
//  TestEquipartition.cpp -- the cheapest ensemble sanity check: does the momentum
//  draw satisfy equipartition, <2 KE> / n_dof = kT, and is n_dof counted correctly
//  when loop-closure constraints remove degrees of freedom?
//
//  THEORY. GC-HMC seeds u = sqrt(RT) * M^{-1/2} g, g ~ N(0,I), so the kinetic
//  energy KE = 1/2 u^T M u has
//        2 KE = u^T M u = RT * g^T M^{-1/2} M M^{-1/2} g = RT * g^T g,
//  hence <2 KE> = RT * n_dof exactly. When the system has n_C loop closures, the
//  RATTLE projection (enforceVelocityConstraints) removes n_C velocity DOF in the
//  M-metric, so the EFFECTIVE count drops to n_dof = nu - n_C and equipartition
//  must read <2 KE>/(nu - n_C) = kT. That second case is the real point: it
//  verifies the constraint DOF bookkeeping the Fixman loop term also depends on.
//
//  This is pure momentum-space statistics (no trajectory), so it is fast and lives
//  in the always-on gate. Tolerance: 4 * stderr of the mean.
// ============================================================================
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "AnalyticForceBridge.hpp"
#include "Constraints.hpp"
#include "HmcDriver.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "StatTest.hpp"
#include "TestHelpers.hpp"
#include "support/TestPhysConstants.hpp"

using namespace robo;
using rtest::attachAtoms;
using rtest::BodySpec;
using rtest::buildBentTorsionChain;
using rtest::buildForest;
using rtest::HmcDriver;
using rtest::Rng;
using rtest::phys::kB;
using rtest::phys::kT300;
using rtest::stat::MeanAccumulator;

namespace {

// Measure <2 KE> over many momentum draws and assert it equals kT * nDofExpected
// within 4 standard errors.
void checkEquipartition(RobotModel& m, RobotState& s, ConstraintSet& cs, int nDofExpected, const char* tag) {
    AnalyticForceBridge bridge(m, s, /*k*/ Real(0)); // KE-only; U irrelevant here
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, kT300, /*h*/ Real(0.01), /*mdSteps*/ 1, 0x9001);

    MeanAccumulator acc;
    const long N = 200000;
    for (long i = 0; i < N; ++i) {
        drv.seedMomenta();
        RobotEngine::realizeVelocity(m, s);
        const double twoKE = 2.0 * static_cast<double>(RobotEngine::calcKineticEnergy(m, s));
        acc.add(twoKE / nDofExpected);
    }
    const double obs = acc.mean();
    const double tol = 4.0 * acc.stderrMean();
    EXPECT_NEAR(obs, kT300, tol) << tag << ": <2KE>/nDof=" << obs << " kT=" << kT300 << " nDof=" << nDofExpected
                                 << " (4*stderr=" << tol << ")";
}

// Draw N momenta at temperature T and accumulate the raw (un-normalized) KE
// samples -- the T-ratio invariant (foundations Sec.5.3, DOF-free) needs the
// mean KE at two different temperatures on the SAME fixture, not the
// equipartition-normalized quantity checkEquipartition asserts.
MeanAccumulator measureMeanKE(RobotModel& m, RobotState& s, ConstraintSet& cs, double T, long N, std::uint64_t seed) {
    const Real RT = static_cast<Real>(kB * T);
    AnalyticForceBridge bridge(m, s, /*k*/ Real(0)); // KE-only
    HmcDriver<AnalyticForceBridge> drv(m, s, bridge, cs, RT, /*h*/ Real(0.01), /*mdSteps*/ 1, seed);

    MeanAccumulator acc;
    for (long i = 0; i < N; ++i) {
        drv.seedMomenta();
        RobotEngine::realizeVelocity(m, s);
        acc.add(static_cast<double>(RobotEngine::calcKineticEnergy(m, s)));
    }
    return acc;
}

// n_dof = 6 (Free root) + nTorsions. Fixed at module scope so the fixture is
// identical between the equipartition test and the T-ratio test.
constexpr int kNTorsions = 8;
constexpr int kMultiBodyNDof = 6 + kNTorsions;

} // namespace

// A free rigid body: 6 DOF (3 rotation + 3 translation).
TEST(Equipartition, FreeBody) {
    BodySpec b;
    b.parent = 0;
    b.joint = JointType::Free;
    b.mass = Real(2.0);
    b.inertia_B = UnitInertia(Real(0.4), Real(0.55), Real(0.7));
    RobotModel m = buildForest({b});
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.1, 0, 0), Vec3(0, 0.1, 0)}, {12.0, 1.0, 1.0});
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1; // unit quaternion
    ConstraintSet cs;
    checkEquipartition(m, s, cs, /*nDof*/ 6, "Free body");
}

// A Free + Torsion + Torsion chain: 6 + 1 + 1 = 8 DOF.
TEST(Equipartition, FreeTorsionChain) {
    rtest::Rng rng(0x55);
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
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1;
    ConstraintSet cs;
    checkEquipartition(m, s, cs, /*nDof*/ 8, "Free+Torsion+Torsion");
}

// A constrained ring: a Free+Torsion+Torsion chain (8 DOF) with ONE loop-closure
// distance constraint. The RATTLE projection removes one velocity DOF, so the
// effective count is 7 -- and equipartition must read kT against 7, not 8. This is
// the constraint DOF-counting check.
TEST(Equipartition, ConstrainedRing) {
    rtest::Rng rng(0x56);
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
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, 0.02, -0.03)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0.04, -0.06, 0.09), Vec3(0.08, 0.01, -0.02)}, {16.0, 1.0});
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1;

    // Place atoms, then close a loop with a Rod constraint set to the CURRENT
    // distance (so the constraint Jacobian G is well defined at this geometry).
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    const int atomA = m.bodyAtoms[m.bodyAtomsBeg[1]]; // first atom of body 1
    const int atomB = m.bodyAtoms[m.bodyAtomsBeg[3]]; // first atom of body 3
    const double d0 = (s.atomPosG()[atomA] - s.atomPosG()[atomB]).norm();
    ConstraintSet cs;
    cs.distance.push_back(DistanceConstraint{atomA, atomB, static_cast<Real>(d0)});

    checkEquipartition(m, s, cs, /*nDof*/ 8 - 1, "Constrained ring (nu - 1)");
}

// ---------------------------------------------------------------------------
//  T0.3 -- the n_C RATTLE removes (equipartition denominator nu - n_C, asserted
//  by ConstrainedRing above) MUST be the SAME accounting the Fixman loop-closure
//  term depends on (foundations Sec.2; World.cpp:860 uses ConstraintSet::
//  calcConstraintLogDet, exactly what HmcDriver::fixmanTerm's lnDetZ is here).
//  Rebuilds the identical ring fixture (same seed/topology as ConstrainedRing) and
//  checks: (a) n_C == cs.distance.size() == 1, the exact quantity subtracted from
//  nu above; (b) that same constraint set's calcConstraintLogDet is genuinely
//  nonzero for this closed ring. If (b) silently returned 0 (the early-return-
//  cyclic bug TestConstraints.cpp::ConstraintLogDetMatchesDenseLnDet guards against),
//  the ConstrainedRing equipartition test would still numerically pass on nu-1 with
//  NO matching Fixman correction -- a constraint-bookkeeping desync this test closes.
// ---------------------------------------------------------------------------
TEST(Equipartition, ConstrainedRingFixmanCoupling) {
    rtest::Rng rng(0x56);
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
    attachAtoms(m, 1, {Vec3(0, 0, 0), Vec3(0.11, 0.02, -0.03)}, {12.0, 1.0});
    attachAtoms(m, 2, {Vec3(0.05, 0.11, 0.0), Vec3(-0.03, 0.0, 0.06)}, {14.0, 1.0});
    attachAtoms(m, 3, {Vec3(0.04, -0.06, 0.09), Vec3(0.08, 0.01, -0.02)}, {16.0, 1.0});
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1;

    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    const int atomA = m.bodyAtoms[m.bodyAtomsBeg[1]]; // first atom of body 1
    const int atomB = m.bodyAtoms[m.bodyAtomsBeg[3]]; // first atom of body 3
    const double d0 = (s.atomPosG()[atomA] - s.atomPosG()[atomB]).norm();
    ConstraintSet cs;
    cs.distance.push_back(DistanceConstraint{atomA, atomB, static_cast<Real>(d0)});

    const int nC = static_cast<int>(cs.distance.size());
    ASSERT_EQ(nC, 1) << "ring fixture must have exactly one loop-closure constraint";
    const int nDofEquipartition = m.nu - nC; // exactly what ConstrainedRing above asserts (8-1=7)
    EXPECT_EQ(nDofEquipartition, 7);

    // The SAME lnDetZ the Fixman loop term (World.cpp:860, HmcDriver::fixmanTerm)
    // uses to correct for these n_C removed DOF.
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    const Real lnDetZ = cs.calcConstraintLogDet(m, s);
    EXPECT_GT(std::abs(lnDetZ), 1e-8)
        << "loop-closure Fixman term (lnDetZ) is ~0 for a genuinely closed ring -- the n_C=" << nC
        << " DOF the equipartition test removes would then have no matching Fixman correction "
           "(constraint-bookkeeping desync)";
}

// ---------------------------------------------------------------------------
//  T0.1 -- a realistic multi-body internal chain where M(phi) is dense and
//  strongly configuration-dependent (unlike the 6-8 DOF toys above): a Free
//  6-DOF root plus kNTorsions Torsion joints in a bent (non-coaxial) chain, so
//  n_dof = 6 + kNTorsions (RobotModel::nu). Same claim as checkEquipartition
//  above, on a fixture that actually exercises the dense articulated-body path.
// ---------------------------------------------------------------------------
TEST(Equipartition, MultiBodyBentChain) {
    Rng rng(0x70);
    RobotModel m = buildBentTorsionChain(kNTorsions, rng);
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1; // unit quaternion for the Free root
    ConstraintSet cs;
    checkEquipartition(m, s, cs, kMultiBodyNDof, "Free+8Torsion bent chain");
}

// ---------------------------------------------------------------------------
//  T0.1 -- the T-ratio invariant (foundations Sec.5.3, DOF-free): <KE(T1)>/
//  <KE(T2)> = T1/T2 holds even though M(phi) is configuration-dependent -- the
//  n_dof count cancels out of the ratio. This is a WEAKER but INDEPENDENT check
//  from the absolute equipartition assertion above: it would still pass with a
//  globally-wrong n_dof (both means scale by the same wrong constant and cancel),
//  but it catches a wrong kT/temperature scaling or a T-dependent bug in the
//  sqrt(RT)*M^{-1/2} draw that a single-temperature check cannot see.
// ---------------------------------------------------------------------------
TEST(Equipartition, MultiBodyBentChain_TRatio) {
    Rng rngM(0x70); // SAME seed as MultiBodyBentChain -> identical fixture geometry
    RobotModel m = buildBentTorsionChain(kNTorsions, rngM);
    RobotState s;
    s.allocateFull(m);
    s.q()[0] = 1;
    ConstraintSet cs;

    const double T1 = 250.0, T2 = 400.0;
    const long N = 200000;
    const MeanAccumulator acc1 = measureMeanKE(m, s, cs, T1, N, 0x7001);
    const MeanAccumulator acc2 = measureMeanKE(m, s, cs, T2, N, 0x7002);

    const double ratio = acc1.mean() / acc2.mean();
    const double expected = T1 / T2;
    // Propagate the ratio's stderr (independent samples, first-order delta method):
    //   Var(A/B) ~= (A/B)^2 * ( (sA/A)^2 + (sB/B)^2 ).
    const double relErr1 = acc1.stderrMean() / acc1.mean();
    const double relErr2 = acc2.stderrMean() / acc2.mean();
    const double stderrRatio = ratio * std::sqrt(relErr1 * relErr1 + relErr2 * relErr2);
    const double tol = 4.0 * stderrRatio;
    EXPECT_NEAR(ratio, expected, tol) << "<KE(T1)>/<KE(T2)>=" << ratio << " != T1/T2=" << expected
                                      << " (4*stderr=" << tol << ")";
}
