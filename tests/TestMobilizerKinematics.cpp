// ============================================================================
//  TestMobilizerKinematics.cpp -- COMPOSITION up to the "robot". We hand-build
//  forests of robots (trees rooted on the shared Ground) and exercise the engine
//  kinematics sweeps (realizePosition / realizeVelocity), no World/forcefield.
//
//  This is the port of Simbody's TestMobilizedBody::testCalculationMethods (two
//  Free bodies on Ground, random state, station/vector/velocity self-consistency)
//  -- translated from the SimTK MobilizedBody API to the RobotModel/RobotEngine
//  SoA: a "station location in Ground" is X_GB[b]*p, a "vector in Ground" is
//  R_GB*v, and so on. We add the structural facts that matter for a FOREST: each
//  robot is kinematically independent through the shared Ground, and a Rigid root
//  is frozen.
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
using rtest::buildForest;
using rtest::randomizeState;
using rtest::Rng;

namespace {

// A forest of several robots on the shared Ground, covering a spread of joints:
//   robot 1: Free root  -> Torsion child   (a 2-body articulated robot)
//   robot 2: Ball root
//   robot 3: Cartesian root -> Slider child
//   robot 4: Rigid (welded) root            (frozen)
RobotModel makeForest(Rng& rng) {
    auto F = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        s.mass = rng.uniform(0.5, 3.0);
        s.com_B = rng.vec3(-0.3, 0.3);
        s.inertia_B = UnitInertia(rng.uniform(0.3, 0.8), rng.uniform(0.3, 0.8), rng.uniform(0.3, 0.8));
        return s;
    };
    std::vector<BodySpec> specs;
    specs.push_back(F(0, JointType::Free));         // body 1 (robot 1 root)
    specs.push_back(F(1, JointType::Torsion));      // body 2 (robot 1 child)
    specs.push_back(F(0, JointType::Ball));         // body 3 (robot 2)
    specs.push_back(F(0, JointType::Cartesian));    // body 4 (robot 3 root)
    specs.push_back(F(4, JointType::Slider));       // body 5 (robot 3 child)
    specs.push_back(F(0, JointType::Rigid));        // body 6 (robot 4, frozen)
    return buildForest(specs);
}

void realizeFull(const RobotModel& m, RobotState& s) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeVelocity(m, s);
}

} // namespace

// ---------------------------------------------------------------------------
//  R1: station / vector transforms (port of testCalculationMethods).
// ---------------------------------------------------------------------------
TEST(RobotKinematics, StationAndVectorTransforms) {
    Rng rng(0xC0FFEE);
    for (int rep = 0; rep < 20; ++rep) {
        RobotModel m = makeForest(rng);
        RobotState s; s.allocateFull(m);
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        const Transform* X_GB = s.X_GB();

        const Vec3 point(0.5, 1.0, -1.5);
        for (int b = 1; b < m.numBodies; ++b) {
            // station(0) in ground == body origin.
            EXPECT_TRUE(rtest::NearVec3(X_GB[b] * Vec3(0), X_GB[b].p(), rtest::kAlg)) << "body " << b;
            // ground -> body -> ground round trip.
            const Vec3 g = X_GB[b] * point;
            const Vec3 backToBody = (~X_GB[b]) * g;
            EXPECT_TRUE(rtest::NearVec3(backToBody, point, rtest::kAlg)) << "round trip body " << b;
            // express a vector in ground = R_GB * v; round trip.
            const Vec3 v(-0.2, 0.7, 0.4);
            const Vec3 vG = X_GB[b].R() * v;
            const Vec3 vBack = X_GB[b].R().transpose() * vG;
            EXPECT_TRUE(rtest::NearVec3(vBack, v, rtest::kAlg)) << "vec round trip body " << b;
        }
        // cross-body: a ground point expressed in body b2 equals (b1 station -> b2).
        const int b1 = 1, b2 = 3;
        const Vec3 g = X_GB[b1] * point;
        const Vec3 inB2_direct = (~X_GB[b2]) * g;
        // distance between two body origins is frame-independent.
        const Real d = (X_GB[b1].p() - X_GB[b2].p()).norm();
        const Vec3 b2OriginInB1 = (~X_GB[b1]) * X_GB[b2].p();
        EXPECT_NEAR(b2OriginInB1.norm(), d, rtest::kAlg);
        (void)inB2_direct;
    }
}

// ---------------------------------------------------------------------------
//  R2: velocity consistency. The spatial velocity V_GB from realizeVelocity must
//  equal d/dt X_GB along the trajectory q(t) = q0 + t*qdot (qdot from calcQDot).
//  Couples realizePosition + calcQDot + realizeVelocity across the whole forest.
// ---------------------------------------------------------------------------
TEST(RobotKinematics, VGBisDerivativeOfXGB) {
    Rng rng(0xD00D);
    const Real h = 1e-6;
    RobotModel m = makeForest(rng);
    RobotState s; s.allocateFull(m);

    for (int rep = 0; rep < 25; ++rep) {
        randomizeState(m, s, rng);
        realizeFull(m, s);

        // snapshot q0, u, and the analytic V_GB.
        std::vector<Real> q0(m.nq), u(m.nu), qdot(m.nq);
        std::copy(s.q(), s.q() + m.nq, q0.begin());
        std::copy(s.u(), s.u() + m.nu, u.begin());
        RobotEngine::calcQDot(m, s, qdot.data());
        std::vector<SpatialVec> Vsnap(m.numBodies);
        std::copy(s.V_GB(), s.V_GB() + m.numBodies, Vsnap.begin());

        auto setAndRealize = [&](Real t) {
            for (int i = 0; i < m.nq; ++i) s.q()[i] = q0[i] + t * qdot[i];
            RobotEngine::normalizeQuaternions(m, s);
            RobotEngine::realizePosition(m, s);
        };
        setAndRealize(+h);
        std::vector<Transform> Xp(m.numBodies);
        std::copy(s.X_GB(), s.X_GB() + m.numBodies, Xp.begin());
        setAndRealize(-h);
        std::vector<Transform> Xm(m.numBodies);
        std::copy(s.X_GB(), s.X_GB() + m.numBodies, Xm.begin());

        for (int b = 1; b < m.numBodies; ++b) {
            // angular vel in Ground from FD of R_GB.
            const Mat33 Rdot = (Xp[b].R() - Xm[b].R()) * (Real(1) / (Real(2) * h));
            const Mat33 W = Rdot * Xp[b].R().transpose();
            const Mat33 W2 = Rdot * Xm[b].R().transpose();
            const Vec3 w_fd(Real(0.5) * (W(2, 1) + W2(2, 1)),
                            Real(0.5) * (W(0, 2) + W2(0, 2)),
                            Real(0.5) * (W(1, 0) + W2(1, 0)));
            const Vec3 v_fd = (Xp[b].p() - Xm[b].p()) * (Real(1) / (Real(2) * h));
            EXPECT_TRUE(rtest::NearVec3(Vsnap[b][0], w_fd, 1e-4)) << "ang body " << b << " rep " << rep;
            EXPECT_TRUE(rtest::NearVec3(Vsnap[b][1], v_fd, 1e-4)) << "lin body " << b << " rep " << rep;
        }
        // restore for the next rep
        std::copy(q0.begin(), q0.end(), s.q());
        std::copy(u.begin(), u.end(), s.u());
    }
}

// ---------------------------------------------------------------------------
//  R3: forest independence. Perturbing one robot's coordinates must not move any
//  body of another robot -- the trees couple ONLY through the shared Ground.
// ---------------------------------------------------------------------------
TEST(RobotForest, RobotsAreKinematicallyIndependent) {
    Rng rng(0x1234);
    RobotModel m = makeForest(rng);
    RobotState s; s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    // robot 3 occupies bodies {4,5}; record their transforms.
    std::vector<Transform> before = {s.X_GB()[4], s.X_GB()[5]};

    // perturb robot 1's coordinates (bodies {1,2} -> q indices of those bodies).
    for (int b : {1, 2}) {
        for (int j = 0; j < m.bodyNQ[b]; ++j) s.q()[m.bodyQIndex[b] + j] += 0.3;
    }
    RobotEngine::normalizeQuaternions(m, s);
    RobotEngine::realizePosition(m, s);

    EXPECT_TRUE(rtest::NearMat33(s.X_GB()[4].R(), before[0].R(), rtest::kTight));
    EXPECT_TRUE(rtest::NearVec3(s.X_GB()[4].p(), before[0].p(), rtest::kTight));
    EXPECT_TRUE(rtest::NearMat33(s.X_GB()[5].R(), before[1].R(), rtest::kTight));
    EXPECT_TRUE(rtest::NearVec3(s.X_GB()[5].p(), before[1].p(), rtest::kTight));
}

// ---------------------------------------------------------------------------
//  R4: a Rigid (welded) root is frozen -- zero dof, and its ground transform is
//  fixed at X_PF * ~X_BM regardless of the rest of the state.
// ---------------------------------------------------------------------------
TEST(RobotForest, RigidRootIsFrozen) {
    Rng rng(0x9001);
    RobotModel m = makeForest(rng);
    RobotState s; s.allocateFull(m);

    const int rigidBody = 6;
    EXPECT_EQ(m.bodyNQ[rigidBody], 0);
    EXPECT_EQ(m.bodyNU[rigidBody], 0);

    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    const Transform fixed = s.X_GB()[rigidBody];
    // X_GB of a Ground-rooted weld = X_PF * X_FM(identity) * ~X_BM.
    const Transform expected = m.X_PF[rigidBody] * (~m.X_BM[rigidBody]);
    EXPECT_TRUE(rtest::NearMat33(fixed.R(), expected.R(), rtest::kAlg));
    EXPECT_TRUE(rtest::NearVec3(fixed.p(), expected.p(), rtest::kAlg));

    // re-randomize everything else; the weld stays put.
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);
    EXPECT_TRUE(rtest::NearMat33(s.X_GB()[rigidBody].R(), fixed.R(), rtest::kTight));
    EXPECT_TRUE(rtest::NearVec3(s.X_GB()[rigidBody].p(), fixed.p(), rtest::kTight));
}
