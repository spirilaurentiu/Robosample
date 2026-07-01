// ============================================================================
//  TestKineticEnergy.cpp -- PHYSICAL, closed-form validation of
//  RobotEngine::calcKineticEnergy against textbook rigid-body mechanics.
//
//  WHY this exists (it is not redundant with TestMassMatrix):
//    The existing KE tests (TestMassMatrix D4 KineticEnergyMatchesSqrtMInv,
//    P4.3 KineticEnergyMatchesDenseForm) are CONSISTENCY checks: their oracles
//    (buildDenseM, multiplyBySqrtMInv) are themselves assembled from the SAME
//    Mk_G, SpatialInertia::operator* and spatialDot primitives that
//    calcKineticEnergy uses. A sign/reference-point/factor error shared by the
//    primitive and the oracle would pass both. None of them tie KE to a value
//    computed WITHOUT the spatial machinery.
//
//  This file's oracle is independent: it computes the per-body kinetic energy
//  from elementary Vec3/Mat33 algebra only -- the standard split of rigid-body
//  KE about the body origin Bo,
//
//      KE_b = 1/2 m |v_Bo|^2  +  m w . (c x v_Bo)  +  1/2 m w . (R I_B R^T) w ,
//
//  where V_GB[b] = (w ; v_Bo) is the body spatial velocity in Ground (angular
//  first), c = R * com_B is the Bo->COM offset in Ground, and I_B is the UNIT
//  inertia about Bo in the body frame (so mass*I_B is the physical inertia).
//  This is the exact expansion of Simbody's authoritative per-body form
//      1/2 V_GB . (Mk_G V_GB)   (RigidBodyNode::calcKineticEnergy),
//  which Robosample reproduces verbatim in RobotEngine::calcKineticEnergy, and
//  the SUM over non-ground bodies mirrors SimbodyMatterSubsystemRep::
//  calcKineticEnergy. If the two agree at random configurations AND velocities,
//  the whole chain (Mk_G construction incl. reexpress handedness, V_GB, the
//  angular/linear ordering, the m[c]x coupling block, the 1/2, ground exclusion)
//  is correct against physics, not just against itself.
// ============================================================================
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::randomizeState;
using rtest::Rng;

namespace {

// Independent per-body spatial kinetic energy from textbook mechanics. Uses ONLY
// Vec3/Mat33/SymMat33 elementary ops and the engine's position/velocity
// kinematics (X_GB, V_GB) -- never SpatialInertia::operator* or spatialDot, the
// two primitives under test. `keCoupling` is returned separately so a test can
// assert the m[c]x block is actually exercised.
Real bodyKEOracle(const RobotModel& m, const RobotState& s, int b, Real* keCoupling = nullptr) {
    const SpatialVec V = s.V_GB()[b];
    const Vec3& w = V.angular;      // body angular velocity in Ground
    const Vec3& v = V.linear;       // velocity of body origin Bo in Ground
    const Rotation& R = s.X_GB()[b].R();
    const Real mass = m.bodyMass[b];

    // c = Bo -> COM in Ground (matches p_BBc_G = R_GB * com_B in realizePosition).
    const Vec3 c = R * m.bodyCom_B[b];

    // I_B (unit inertia about Bo, body frame) re-expressed to Ground as R I_B R^T,
    // built with plain Mat33 multiplies -- deliberately NOT UnitInertia::reexpress,
    // so a handedness bug in reexpress cannot hide by appearing on both sides.
    const Mat33 Ib = m.bodyUnitInertia_B[b].full();
    const Mat33 Ig = (R * Ib) * R.transpose();
    const Vec3 Igw = Ig * w;

    const Real keLin = Real(0.5) * mass * dot(v, v);
    const Real keCpl = mass * dot(w, c % v);        // w . (c x v)
    const Real keRot = Real(0.5) * mass * dot(w, Igw);
    if (keCoupling != nullptr) {
        *keCoupling = keCpl;
    }
    return keLin + keCpl + keRot;
}

// Single free body with a chosen mass / COM / inertia and identity joint frames.
RobotModel makeFreeBody(Real mass, const Vec3& comB, const UnitInertia& inertiaB) {
    BodySpec s;
    s.parent = 0;
    s.joint = JointType::Free;
    s.X_PF = Transform();
    s.X_BM = Transform();
    s.mass = mass;
    s.com_B = comB;
    s.inertia_B = inertiaB;
    return buildForest({s});
}

// A mixed forest: several robots, offset COMs, anisotropic inertias -- so the
// per-body sum is non-trivial and every joint's H participates in V_GB.
RobotModel makeMixedForest(Rng& rng) {
    auto F = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        s.mass = rng.uniform(0.5, 3.0);
        s.com_B = rng.vec3(-0.3, 0.3);   // offset COM -> exercises the coupling block
        s.inertia_B = UnitInertia(rng.uniform(0.3, 0.8), rng.uniform(0.3, 0.8), rng.uniform(0.3, 0.8));
        return s;
    };
    std::vector<BodySpec> specs;
    specs.push_back(F(0, JointType::Free));      // 1
    specs.push_back(F(1, JointType::Ball));       // 2
    specs.push_back(F(2, JointType::Torsion));    // 3
    specs.push_back(F(0, JointType::Cartesian));  // 4 (second robot)
    specs.push_back(F(4, JointType::Cylinder));   // 5
    return buildForest(specs);
}

void zeroConfig(const RobotModel& m, RobotState& s) {
    for (int i = 0; i < m.nq; ++i) {
        s.q()[i] = Real(0);
    }
    // Unit quaternions where present (all-zero quat is not a valid rotation).
    for (int qs : m.quaternionQStart) {
        s.q()[qs + 0] = Real(1);
        s.q()[qs + 1] = Real(0);
        s.q()[qs + 2] = Real(0);
        s.q()[qs + 3] = Real(0);
    }
    for (int i = 0; i < m.nu; ++i) {
        s.u()[i] = Real(0);
    }
}

} // namespace

// ---------------------------------------------------------------------------
//  Golden: pure translation. A 3-dof Cartesian body at identity has V_GB =
//  (0 ; u), so KE must be exactly 1/2 m |v|^2 -- a hand-checkable number that
//  pins the mass, the linear inertia block, and the absence of any spurious
//  angular contribution. m=2, v=(3,4,0) -> KE = 1/2 * 2 * 25 = 25.
// ---------------------------------------------------------------------------
TEST(KineticEnergy, TranslationOnlyIsHalfMVSquared) {
    RobotModel m = [] {
        BodySpec s;
        s.parent = 0;
        s.joint = JointType::Cartesian;
        s.X_PF = Transform();
        s.X_BM = Transform();
        s.mass = Real(2);
        s.com_B = Vec3(0);
        s.inertia_B = UnitInertia(Real(0.4), Real(0.5), Real(0.6));
        return buildForest({s});
    }();
    RobotState s;
    s.allocateFull(m);
    zeroConfig(m, s);
    RobotEngine::realizePosition(m, s);

    const Real vel[3] = {Real(3), Real(4), Real(0)};
    for (int i = 0; i < 3; ++i) {
        s.u()[i] = vel[i];
    }
    RobotEngine::realizeVelocity(m, s);

    // Document/verify the kinematic mapping this golden relies on.
    const SpatialVec V = s.V_GB()[1];
    EXPECT_NEAR(std::sqrt(dot(V.angular, V.angular)), 0.0, rtest::kLoose) << "Cartesian must not rotate";
    EXPECT_NEAR(V.linear[0], vel[0], rtest::kLoose);
    EXPECT_NEAR(V.linear[1], vel[1], rtest::kLoose);
    EXPECT_NEAR(V.linear[2], vel[2], rtest::kLoose);

    const Real ke = RobotEngine::calcKineticEnergy(m, s);
    EXPECT_NEAR(ke, 25.0, rtest::kLoose);                    // 1/2 * 2 * (9+16)
    EXPECT_NEAR(ke, bodyKEOracle(m, s, 1), rtest::kLoose);
}

// ---------------------------------------------------------------------------
//  Free body, COM at the origin (no coupling block): calcKineticEnergy equals
//  the textbook 1/2 m|v|^2 + 1/2 m w.(R I R^T)w at random orientation and
//  random spatial velocity. Isolates the rotational-inertia reexpression and
//  the angular/linear split.
// ---------------------------------------------------------------------------
TEST(KineticEnergy, FreeBodyNoOffsetMatchesClosedForm) {
    Rng rng(0x4B02);
    RobotModel m = makeFreeBody(Real(1.7), Vec3(0), UnitInertia(Real(0.35), Real(0.55), Real(0.8)));
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 25; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);

        const Real ke = RobotEngine::calcKineticEnergy(m, s);
        EXPECT_NEAR(ke, bodyKEOracle(m, s, 1), rtest::kLoose) << "rep " << rep;
        EXPECT_GT(ke, 0.0) << "rep " << rep;   // KE is a positive-definite form
    }
}

// ---------------------------------------------------------------------------
//  Free body with an OFFSET COM: now the m[c]x coupling block is live. This is
//  the term that distinguishes an inertia about the body origin (Simbody's and
//  Robosample's convention) from one about the COM; getting the reference point
//  wrong changes KE here but not in the no-offset test. We also assert the
//  coupling term is materially non-zero in at least one rep, so the test really
//  exercises it rather than passing vacuously.
// ---------------------------------------------------------------------------
TEST(KineticEnergy, FreeBodyOffsetCOMMatchesClosedForm) {
    Rng rng(0x4B03);
    RobotModel m = makeFreeBody(Real(2.4), Vec3(Real(0.2), Real(-0.15), Real(0.3)),
                                UnitInertia(Real(0.4), Real(0.6), Real(0.5)));
    RobotState s;
    s.allocateFull(m);
    Real maxCoupling = 0;
    for (int rep = 0; rep < 25; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);

        Real keCpl = 0;
        const Real oracle = bodyKEOracle(m, s, 1, &keCpl);
        maxCoupling = std::max(maxCoupling, std::abs(keCpl));
        EXPECT_NEAR(RobotEngine::calcKineticEnergy(m, s), oracle, rtest::kLoose) << "rep " << rep;
    }
    EXPECT_GT(maxCoupling, 1e-3) << "coupling block never exercised -- test is vacuous";
}

// ---------------------------------------------------------------------------
//  Multi-body forest: calcKineticEnergy equals the SUM over non-ground bodies
//  of the independent per-body oracle. This mirrors, line for line, Simbody's
//  SimbodyMatterSubsystemRep::calcKineticEnergy (skip level 0 = Ground, sum
//  each node's 1/2 V.Mk V), validating the assembly and the ground exclusion.
// ---------------------------------------------------------------------------
TEST(KineticEnergy, MixedForestMatchesPerBodySum) {
    Rng rng(0x4B04);
    RobotModel m = makeMixedForest(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 20; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);

        Real oracle = 0;
        for (int b = 1; b < m.numBodies; ++b) { // b = 0 is Ground: excluded, exactly as Simbody
            oracle += bodyKEOracle(m, s, b);
        }
        const Real ke = RobotEngine::calcKineticEnergy(m, s);
        EXPECT_NEAR(ke, oracle, rtest::kLoose) << "rep " << rep;
        EXPECT_GT(ke, 0.0) << "rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  Quadratic-form sanity, independent of the oracle: KE is a quadratic form in
//  the spatial velocity, so KE(-u) == KE(u) and KE(a*u) == a^2 KE(u). These
//  hold for the physical form regardless of configuration and catch any stray
//  linear-in-velocity term.
// ---------------------------------------------------------------------------
TEST(KineticEnergy, IsQuadraticInVelocity) {
    Rng rng(0x4B05);
    RobotModel m = makeMixedForest(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        std::vector<Real> u(s.u(), s.u() + m.nu);

        RobotEngine::realizeVelocity(m, s);
        const Real ke = RobotEngine::calcKineticEnergy(m, s);

        for (int i = 0; i < m.nu; ++i) {
            s.u()[i] = -u[i];
        }
        RobotEngine::realizeVelocity(m, s);
        EXPECT_NEAR(RobotEngine::calcKineticEnergy(m, s), ke, rtest::kLoose) << "rep " << rep;

        const Real a = Real(1.5);
        for (int i = 0; i < m.nu; ++i) {
            s.u()[i] = a * u[i];
        }
        RobotEngine::realizeVelocity(m, s);
        EXPECT_NEAR(RobotEngine::calcKineticEnergy(m, s), a * a * ke, rtest::kLoose) << "rep " << rep;
    }
}
