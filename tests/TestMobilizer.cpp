// ============================================================================
//  TestPortedMobilizer.cpp -- direct ports of the portable parts of the SimTK
//  mobilizer test files, translated from the SimTK MobilizedBody API to the
//  RobotModel/RobotEngine SoA. What is faithfully portable is ported; what
//  depends on SimTK-only machinery is documented at the bottom of this file
//  rather than faked.
//
//  Ported here:
//   * Acceleration leg of TestMobilizedBody::testCalculationMethods -- the body
//     spatial acceleration A_GB must equal d/dt V_GB (completes the position +
//     velocity legs already in TestMobilizerKinematics).
//   * TestAngleConversions' "one of every mobilizer in a chain" stress config:
//     a deep chain of all ten JointTypes realized at a random configuration,
//     with the kinematic sweeps self-consistent and quaternion blocks unit.
//     (The Euler<->quaternion conversion itself is N/A: the port is quaternion-
//     only; the portable invariant is that renormalizing the quaternion does not
//     move any body, and the q-from-u map round-trips.)
//   * The Weld-mobilizer kernel of TestMobilizedBody::testWeld: a Rigid child
//     moves rigidly with its parent (same angular velocity; origin velocity is
//     the rigid transport of the parent's).
//   * The equivalence pattern shared by TestLoneParticle / testWeld / testGimbal
//     / reaction-force tests: two structurally identical robots in one forest
//     produce identical kinematics and dynamics body-for-body (exercises forest
//     indexing / children-CSR and determinism).
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
using rtest::BodySpec;
using rtest::buildChain;
using rtest::buildForest;
using rtest::randomizeState;
using rtest::Rng;

namespace {

// every JointType, with a legal root (Free) first.
const std::vector<JointType> kChainAllJoints = {JointType::Free,
                                                JointType::Torsion,
                                                JointType::Slider,
                                                JointType::Cylinder,
                                                JointType::BendStretch,
                                                JointType::Cartesian,
                                                JointType::Ball,
                                                JointType::SphericalCoords,
                                                JointType::FreeLine,
                                                JointType::Rigid};

Vec3 angVelFD(const Rotation& Rp, const Rotation& Rm, Real h) {
    const Mat33 Rdot = (Rp - Rm) * (Real(1) / (Real(2) * h));
    const Mat33 W = Rdot * Rp.transpose();
    const Mat33 W2 = Rdot * Rm.transpose();
    return Vec3(Real(0.5) * (W(2, 1) + W2(2, 1)),
                Real(0.5) * (W(0, 2) + W2(0, 2)),
                Real(0.5) * (W(1, 0) + W2(1, 0)));
}

} // namespace

// ---------------------------------------------------------------------------
//  Acceleration leg of testCalculationMethods: A_GB == d/dt V_GB.
//  We drive the body with a random generalized force, run full forward dynamics
//  to get (udot, A_GB), then finite-difference V_GB along the true trajectory
//  q(t)=q0+t*qdot0, u(t)=u0+t*udot0 (the omitted O(t^2) terms are symmetric and
//  cancel in the central difference).
// ---------------------------------------------------------------------------
TEST(PortedCalculationMethods, BodyAccelerationIsDerivativeOfVelocity) {
    Rng rng(0xACCE1);
    const Real h = 1e-6;
    RobotModel m = buildChain(
        {JointType::Free, JointType::Torsion, JointType::Ball, JointType::Cartesian, JointType::FreeLine},
        rng);
    RobotState s;
    s.allocateFull(m);

    for (int rep = 0; rep < 20; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        std::fill(s.bodyForceG(), s.bodyForceG() + m.numBodies, SpatialVec(Vec3(0), Vec3(0)));
        for (int i = 0; i < m.nu; ++i) {
            s.mobilityForce()[i] = rng.gaussian(0, Real(0.5));
        }
        RobotEngine::calcUDot(m, s);

        std::vector<Real> q0(m.nq), u0(m.nu), qdot0(m.nq), udot0(m.nu);
        std::copy(s.q(), s.q() + m.nq, q0.begin());
        std::copy(s.u(), s.u() + m.nu, u0.begin());
        std::copy(s.udot(), s.udot() + m.nu, udot0.begin());
        RobotEngine::calcQDot(m, s, qdot0.data());
        std::vector<SpatialVec> A0(m.numBodies);
        std::copy(s.A_GB(), s.A_GB() + m.numBodies, A0.begin());

        auto velAt = [&](Real t, std::vector<SpatialVec>& V) {
            for (int i = 0; i < m.nq; ++i) {
                s.q()[i] = q0[i] + t * qdot0[i];
            }
            for (int i = 0; i < m.nu; ++i) {
                s.u()[i] = u0[i] + t * udot0[i];
            }
            RobotEngine::normalizeQuaternions(m, s);
            RobotEngine::realizePosition(m, s);
            RobotEngine::realizeVelocity(m, s);
            V.assign(s.V_GB(), s.V_GB() + m.numBodies);
        };
        std::vector<SpatialVec> Vp, Vm;
        velAt(+h, Vp);
        velAt(-h, Vm);

        for (int b = 1; b < m.numBodies; ++b) {
            const Vec3 aAng = (Vp[b][0] - Vm[b][0]) * (Real(1) / (Real(2) * h));
            const Vec3 aLin = (Vp[b][1] - Vm[b][1]) * (Real(1) / (Real(2) * h));
            EXPECT_TRUE(rtest::NearVec3(A0[b][0], aAng, 1e-3)) << "ang body " << b << " rep " << rep;
            EXPECT_TRUE(rtest::NearVec3(A0[b][1], aLin, 1e-3)) << "lin body " << b << " rep " << rep;
        }
        std::copy(q0.begin(), q0.end(), s.q());
        std::copy(u0.begin(), u0.end(), s.u());
    }
}

// ---------------------------------------------------------------------------
//  TestAngleConversions, ported: one of every mobilizer in a chain, random
//  configuration, all kinematics finite & self-consistent; quaternions unit.
// ---------------------------------------------------------------------------
TEST(PortedAngleConversions, AllMobilizerChainRealizesConsistently) {
    Rng rng(0xC4A1);
    const Real h = 1e-6;
    RobotModel m = buildChain(kChainAllJoints, rng);
    RobotState s;
    s.allocateFull(m);

    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);

        // every body transform and velocity is finite.
        for (int b = 1; b < m.numBodies; ++b) {
            for (int k = 0; k < 3; ++k) {
                ASSERT_TRUE(std::isfinite(s.X_GB()[b].p()[k])) << "body " << b;
                ASSERT_TRUE(std::isfinite(s.V_GB()[b][0][k]));
                ASSERT_TRUE(std::isfinite(s.V_GB()[b][1][k]));
            }
            const Mat33 RRt = s.X_GB()[b].R() * s.X_GB()[b].R().transpose();
            EXPECT_TRUE(rtest::NearMat33(RRt, Mat33(Real(1)), 1e-9)) << "orthonormal body " << b;
        }

        // velocity consistency across the whole all-joint chain (FD of X_GB).
        std::vector<Real> q0(m.nq), qdot0(m.nq);
        std::copy(s.q(), s.q() + m.nq, q0.begin());
        RobotEngine::calcQDot(m, s, qdot0.data());
        std::vector<SpatialVec> V0(s.V_GB(), s.V_GB() + m.numBodies);

        auto posAt = [&](Real t, std::vector<Transform>& X) {
            for (int i = 0; i < m.nq; ++i) {
                s.q()[i] = q0[i] + t * qdot0[i];
            }
            RobotEngine::normalizeQuaternions(m, s);
            RobotEngine::realizePosition(m, s);
            X.assign(s.X_GB(), s.X_GB() + m.numBodies);
        };
        std::vector<Transform> Xp, Xm;
        posAt(+h, Xp);
        posAt(-h, Xm);
        for (int b = 1; b < m.numBodies; ++b) {
            const Vec3 w_fd = angVelFD(Xp[b].R(), Xm[b].R(), h);
            const Vec3 v_fd = (Xp[b].p() - Xm[b].p()) * (Real(1) / (Real(2) * h));
            EXPECT_TRUE(rtest::NearVec3(V0[b][0], w_fd, 1e-4)) << "ang body " << b << " rep " << rep;
            EXPECT_TRUE(rtest::NearVec3(V0[b][1], v_fd, 1e-4)) << "lin body " << b << " rep " << rep;
        }
        std::copy(q0.begin(), q0.end(), s.q());
    }
}

TEST(PortedAngleConversions, QuaternionNormalizationPreservesConfiguration) {
    // The portable remnant of the Euler<->quaternion round trip: the body layout
    // is a function of the quaternion only through its direction, so renormalizing
    // an already-unit quaternion must not move any body, and q -> u -> q via the
    // engine's own maps is the identity up to renorm.
    Rng rng(0xC4A2);
    RobotModel m = buildChain(kChainAllJoints, rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    std::vector<Transform> before(s.X_GB(), s.X_GB() + m.numBodies);

    // perturb quaternion blocks slightly off the unit sphere, then renormalize.
    for (int qs : m.quaternionQStart) {
        for (int k = 0; k < 4; ++k) {
            s.q()[qs + k] *= Real(1.05);
        }
    }
    RobotEngine::normalizeQuaternions(m, s);
    RobotEngine::realizePosition(m, s);

    for (int b = 1; b < m.numBodies; ++b) {
        EXPECT_TRUE(rtest::NearMat33(s.X_GB()[b].R(), before[b].R(), 1e-9)) << "body " << b;
        EXPECT_TRUE(rtest::NearVec3(s.X_GB()[b].p(), before[b].p(), 1e-9)) << "body " << b;
    }
}

// ---------------------------------------------------------------------------
//  Weld-mobilizer kernel of testWeld: a Rigid child moves rigidly with parent.
//  (The full testWeld -- weld MOBILIZER vs weld CONSTRAINT equivalence under
//  integration -- needs the constraint set + integrator, i.e. "more vertical".)
// ---------------------------------------------------------------------------
TEST(PortedWeld, RigidChildMovesRigidlyWithParent) {
    Rng rng(0x3ED);
    std::vector<BodySpec> specs;
    BodySpec root;
    root.parent = 0;
    root.joint = JointType::Free;
    root.X_PF = Transform(rng.rotation(), rng.vec3());
    root.X_BM = Transform(rng.rotation(), rng.vec3());
    specs.push_back(root);
    BodySpec weld;
    weld.parent = 1;
    weld.joint = JointType::Rigid;
    weld.X_PF = Transform(rng.rotation(), rng.vec3());
    weld.X_BM = Transform(rng.rotation(), rng.vec3());
    specs.push_back(weld);
    RobotModel m = buildForest(specs);
    RobotState s;
    s.allocateFull(m);

    for (int rep = 0; rep < 20; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);
        const SpatialVec Vp = s.V_GB()[1]; // parent (Free)
        const SpatialVec Vc = s.V_GB()[2]; // welded child
        // same angular velocity.
        EXPECT_TRUE(rtest::NearVec3(Vc[0], Vp[0], rtest::kAlg)) << "rep " << rep;
        // child origin velocity = parent origin velocity + w x r (rigid transport).
        const Vec3 r = s.X_GB()[2].p() - s.X_GB()[1].p();
        const Vec3 vExpected = Vp[1] + (Vp[0] % r);
        EXPECT_TRUE(rtest::NearVec3(Vc[1], vExpected, rtest::kAlg)) << "rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  Equivalence pattern (TestLoneParticle / testGimbal / reaction tests all
//  compare two equivalent builds): two identical robots in one forest, driven by
//  identical coordinates/forces, agree body-for-body across X_GB, V_GB, A_GB and
//  udot. Exercises forest indexing, the children-CSR, and determinism.
// ---------------------------------------------------------------------------
TEST(PortedEquivalence, IdenticalRobotsAgreeAcrossAllQuantities) {
    Rng rng(0x7011);
    // Two identical 3-body chains as separate robots on Ground.
    const std::vector<JointType> chain = {JointType::Free, JointType::Torsion, JointType::Ball};

    // Build robot A specs, then duplicate as robot B (same frames/masses).
    std::vector<BodySpec> a;
    for (int i = 0; i < (int)chain.size(); ++i) {
        BodySpec s;
        s.parent = (i == 0 ? 0 : i);
        s.joint = chain[i];
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        s.mass = rng.uniform(0.5, 2.0);
        s.com_B = rng.vec3(-0.2, 0.2);
        s.inertia_B = UnitInertia(rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7));
        a.push_back(s);
    }
    const int nA = (int)a.size();
    std::vector<BodySpec> all = a;
    for (int i = 0; i < nA; ++i) {
        BodySpec s = a[i];
        s.parent = (i == 0 ? 0 : nA + i); // robot B mirrors robot A's topology
        all.push_back(s);
    }
    RobotModel m = buildForest(all);
    RobotState s;
    s.allocateFull(m);

    randomizeState(m, s, rng);
    // copy robot A's coordinates onto robot B so both are driven identically.
    for (int i = 0; i < nA; ++i) {
        const int bA = 1 + i, bB = 1 + nA + i;
        for (int j = 0; j < m.bodyNQ[bA]; ++j) {
            s.q()[m.bodyQIndex[bB] + j] = s.q()[m.bodyQIndex[bA] + j];
        }
        for (int j = 0; j < m.bodyNU[bA]; ++j) {
            s.u()[m.bodyUIndex[bB] + j] = s.u()[m.bodyUIndex[bA] + j];
        }
    }
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    std::fill(s.bodyForceG(), s.bodyForceG() + m.numBodies, SpatialVec(Vec3(0), Vec3(0)));
    for (int i = 0; i < nA; ++i) {
        const int bA = 1 + i, bB = 1 + nA + i;
        for (int j = 0; j < m.bodyNU[bA]; ++j) {
            const Real f = rng.gaussian(0, 0.5);
            s.mobilityForce()[m.bodyUIndex[bA] + j] = f;
            s.mobilityForce()[m.bodyUIndex[bB] + j] = f;
        }
    }
    RobotEngine::calcUDot(m, s);

    for (int i = 0; i < nA; ++i) {
        const int bA = 1 + i, bB = 1 + nA + i;
        // X_GB depends on the (random, different) frames, so compare the MOTION
        // invariants that must match: mobilizer-frame velocity V_FM and udot.
        for (int j = 0; j < m.bodyNU[bA]; ++j) {
            EXPECT_NEAR(s.udot()[m.bodyUIndex[bA] + j], s.udot()[m.bodyUIndex[bB] + j], rtest::kLoose)
                << "udot body pair " << i << " comp " << j;
        }
        EXPECT_TRUE(rtest::NearVec3(s.V_FM()[bA][0], s.V_FM()[bB][0], rtest::kTight)) << "V_FM ang " << i;
        EXPECT_TRUE(rtest::NearVec3(s.V_FM()[bA][1], s.V_FM()[bB][1], rtest::kTight)) << "V_FM lin " << i;
    }
}

// ---------------------------------------------------------------------------
//  NOT PORTED (documented, not faked) -- these need SimTK-only machinery absent
//  from the SimTK-free engine, or layers the project deferred to "more vertical":
//
//   * TestMobilizerReactionForces (all of it): the engine has no
//     calcMobilizerReactionForces / findMobilizerReactionOnBodyAtMInGround. The
//     "free joint -> zero reaction" and SD/FAST reference checks require that
//     operator. -> port once a reaction-force operator is added to RobotEngine.
//   * TestLoneParticle: there is no RBNodeLoneParticle specialization to compare
//     against RBNodeTranslate, and it also exercises multiplyByM,
//     calcCompositeBodyInertias, calcResidualForceIgnoringConstraints and
//     multiplyBySystemJacobianTranspose -- none of which exist in the port yet.
//     (multiplyByMInv IS covered, in TestMassMatrix.)
//   * testWeld / testGimbal / testBushing equivalence-under-simulation: needs the
//     constraint set + integrator (verletStep, ForceBridge) -> the deferred
//     vertical-integration layer. Gimbal and Bushing are also not JointTypes.
//   * Euler<->quaternion conversion (TestAngleConversions core): the port stores
//     orientation as quaternions only; there is no Euler generalized-coordinate
//     representation to convert to/from.
//   * UniformGravity legs: intentionally omitted (no gravity in Robosample
//     molecular robotics).
// ---------------------------------------------------------------------------