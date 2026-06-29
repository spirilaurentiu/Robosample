// ============================================================================
//  TestReverseMobilizer.cpp -- the portable core of Simbody's
//  TestReverseMobilizers, translated to the SimTK-free engine.
//
//  WHAT REVERSE MOBILIZERS ARE, AND WHY THE LITERAL TEST IS N/A HERE:
//  Simbody's MobilizedBody::Reverse builds a joint whose F (inboard) and M
//  (outboard) frames are swapped, producing a reversed H/HDot. molmodel never
//  uses this when mapping a molecule to a multibody tree (CompoundSystem.cpp uses
//  "Reverse" exactly zero times -- the tree is always built forward, outward from
//  the root), and the Robosample port has no reverse flag in RobotModel. So the
//  forward-system-vs-reverse-system comparison the original test performs cannot
//  be reproduced: there is no reverse mobilizer to compare against.
//
//  WHAT IS PORTABLE (and ported here): the physical invariants the original test
//  rests on, which the forward engine MUST satisfy --
//   1. the spatial-velocity reversal map is an involution (reversing twice is
//      identity) and equals the rigid spatial-transform inverse;
//   2. the engine's cross-mobilizer velocity V_FM is consistent with the body
//      velocities V_GB it was assembled from (independently re-derived from the
//      parent/child spatial velocities and the static joint frames) -- across all
//      mobilizer types;
//   3. reversing the engine's forward V_FM yields the correct relative velocity
//      seen from the other end (M->F), i.e. forward and reverse descriptions of
//      the SAME relative motion agree -- which is exactly the equality the
//      original fwd/rev test checks, minus the reverse-mobilizer machinery.
//
//  NOT PORTED (documented at bottom): reaction-force reversal (no reaction-force
//  operator in the port) and acceleration/HDot reversal (the HDot-reversal terms
//  are the reverse-mobilizer-specific machinery that does not exist here; forward
//  A_GB correctness is already covered by TestPortedMobilizer). Gravity omitted.
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

// Reverse a spatial velocity across a rigid transform: given X_AB (frame B in
// frame A) and V_AB (spatial velocity of B measured in A, expressed in A),
// return V_BA (spatial velocity of A measured in B, expressed in B). This is the
// exact map used in Simbody's TestReverseMobilizers to relate fwd/rev speeds.
SpatialVec reverseSpatialVel(const Transform& X_AB, const SpatialVec& V_AB) {
    const Transform X_BA = ~X_AB;
    return SpatialVec((X_BA.R() * V_AB[0]) * Real(-1),
                      (X_BA.R() * (V_AB[1] + (X_AB.p() % V_AB[0]))) * Real(-1));
}

const std::array<JointType, 8> kChildJoints = {JointType::Torsion,
                                               JointType::Slider,
                                               JointType::Cylinder,
                                               JointType::Cartesian,
                                               JointType::BendStretch,
                                               JointType::Ball,
                                               JointType::FreeLine,
                                               JointType::Free};

const char* jn(JointType jt) {
    switch (jt) {
        case JointType::Torsion:
            return "Torsion";
        case JointType::Slider:
            return "Slider";
        case JointType::Cylinder:
            return "Cylinder";
        case JointType::Cartesian:
            return "Cartesian";
        case JointType::BendStretch:
            return "BendStretch";
        case JointType::Ball:
            return "Ball";
        case JointType::FreeLine:
            return "FreeLine";
        case JointType::Free:
            return "Free";
        default:
            return "?";
    }
}

// Build Ground -> Free(A) -> jt(B), random frames/masses.
RobotModel twoBody(JointType jt, Rng& rng) {
    auto mk = [&](int parent, JointType j) {
        BodySpec s;
        s.parent = parent;
        s.joint = j;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        s.mass = rng.uniform(0.5, 2.0);
        s.com_B = rng.vec3(-0.2, 0.2);
        s.inertia_B = UnitInertia(rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7));
        return s;
    };
    return buildForest({mk(0, JointType::Free), mk(1, jt)});
}

} // namespace

// ---------------------------------------------------------------------------
//  1. The spatial-velocity reversal map is an involution and matches the rigid
//     spatial-transform inverse. Pure algebra (no engine).
// ---------------------------------------------------------------------------
TEST(ReverseMobilizer, SpatialVelocityReversalIsInvolution) {
    Rng rng(0x5EED1234);
    for (int t = 0; t < 500; ++t) {
        const Transform X_AB(rng.rotation(), rng.vec3());
        const SpatialVec V_AB(rng.vec3(), rng.vec3());
        const SpatialVec V_BA = reverseSpatialVel(X_AB, V_AB);
        // reversing again (now across X_BA) must return the original.
        const SpatialVec V_AB2 = reverseSpatialVel(~X_AB, V_BA);
        EXPECT_TRUE(rtest::NearVec3(V_AB2[0], V_AB[0], rtest::kAlg)) << "t=" << t;
        EXPECT_TRUE(rtest::NearVec3(V_AB2[1], V_AB[1], rtest::kAlg)) << "t=" << t;
    }
}

// ---------------------------------------------------------------------------
//  2. Cross-mobilizer velocity consistency: the engine's V_FM equals the
//     relative spatial velocity of M w.r.t. F, independently re-derived from the
//     parent/child body velocities and the static joint frames. Across every
//     mobilizer type. (V_FM was FD-validated against d/dt X_FM elsewhere; here we
//     pin it to the assembled body velocities -- the relation the fwd/rev test
//     exploits.)
// ---------------------------------------------------------------------------
TEST(ReverseMobilizer, EngineVFMMatchesRelativeVelocityFromBodies) {
    Rng rng(0x2222);
    for (JointType jt : kChildJoints) {
        RobotModel m = twoBody(jt, rng);
        RobotState s;
        s.allocateFull(m);
        for (int rep = 0; rep < 30; ++rep) {
            randomizeState(m, s, rng);
            RobotEngine::realizePosition(m, s);
            RobotEngine::realizeVelocity(m, s);

            const int b = 2, p = 1;
            const Transform X_GF = s.X_GB()[p] * m.X_PF[b];
            const Transform X_GM = s.X_GB()[b] * m.X_BM[b];
            const Rotation R_GF = X_GF.R();
            const Vec3 p_F = X_GF.p(), p_M = X_GM.p();

            const Vec3 w_p = s.V_GB()[p][0], v_p = s.V_GB()[p][1], oP = s.X_GB()[p].p();
            const Vec3 w_c = s.V_GB()[b][0], v_c = s.V_GB()[b][1], oC = s.X_GB()[b].p();

            // relative spatial velocity of M measured in F, expressed in F.
            const Vec3 relW_G = w_c - w_p;
            const Vec3 vM_G = v_c + (w_c % (p_M - oC));   // M-origin velocity (on child)
            const Vec3 vFpt_G = v_p + (w_p % (p_M - oP)); // F-point coincident w/ M-origin
            const Vec3 vfm_ang = R_GF.transpose() * relW_G;
            const Vec3 vfm_lin = R_GF.transpose() * (vM_G - vFpt_G);

            EXPECT_TRUE(rtest::NearVec3(s.V_FM()[b][0], vfm_ang, 1e-10)) << jn(jt) << " ang rep " << rep;
            EXPECT_TRUE(rtest::NearVec3(s.V_FM()[b][1], vfm_lin, 1e-10)) << jn(jt) << " lin rep " << rep;
        }
    }
}

// ---------------------------------------------------------------------------
//  3. Forward/reverse agreement: reversing the engine's forward V_FM yields the
//     correct relative velocity seen from M (i.e. F-in-M), independently derived
//     from body velocities. This is the faithful analogue of "forward mobilizer
//     == reverse mobilizer": both describe the SAME relative motion.
// ---------------------------------------------------------------------------
TEST(ReverseMobilizer, ReversingForwardVFMGivesReverseRelativeVelocity) {
    Rng rng(0x3333);
    for (JointType jt : kChildJoints) {
        RobotModel m = twoBody(jt, rng);
        RobotState s;
        s.allocateFull(m);
        for (int rep = 0; rep < 30; ++rep) {
            randomizeState(m, s, rng);
            RobotEngine::realizePosition(m, s);
            RobotEngine::realizeVelocity(m, s);

            const int b = 2, p = 1;
            const Transform X_GF = s.X_GB()[p] * m.X_PF[b];
            const Transform X_GM = s.X_GB()[b] * m.X_BM[b];
            const Rotation R_GM = X_GM.R();
            const Vec3 p_F = X_GF.p(), p_M = X_GM.p();

            const Vec3 w_p = s.V_GB()[p][0], v_p = s.V_GB()[p][1], oP = s.X_GB()[p].p();
            const Vec3 w_c = s.V_GB()[b][0], v_c = s.V_GB()[b][1], oC = s.X_GB()[b].p();

            // reverse via the reversal map applied to the engine's forward V_FM,
            // using X_FM = ~X_GF * X_GM (= the engine's stored X_FM).
            const Transform X_FM = (~X_GF) * X_GM;
            const SpatialVec V_MF_reversed = reverseSpatialVel(X_FM, s.V_FM()[b]);

            // independent reference: relative spatial velocity of F measured in M.
            const Vec3 relW_G = w_p - w_c;
            const Vec3 vF_G = v_p + (w_p % (p_F - oP));   // F-origin velocity (on parent)
            const Vec3 vMpt_G = v_c + (w_c % (p_F - oC)); // M-point coincident w/ F-origin
            const Vec3 vmf_ang = R_GM.transpose() * relW_G;
            const Vec3 vmf_lin = R_GM.transpose() * (vF_G - vMpt_G);

            EXPECT_TRUE(rtest::NearVec3(V_MF_reversed[0], vmf_ang, 1e-9)) << jn(jt) << " ang rep " << rep;
            EXPECT_TRUE(rtest::NearVec3(V_MF_reversed[1], vmf_lin, 1e-9)) << jn(jt) << " lin rep " << rep;
        }
    }
}

// ---------------------------------------------------------------------------
//  NOT PORTED (documented, not faked):
//   * Reaction-force reversal (reacBA == rev reaction): no calcMobilizerReaction-
//     Forces operator in the port -> port once that operator exists.
//   * Acceleration / HDot reversal: the reverse-mobilizer HDot terms are the
//     feature-specific machinery absent here; forward A_GB correctness is covered
//     in TestPortedMobilizer (BodyAccelerationIsDerivativeOfVelocity).
//   * The full forward-system-vs-reverse-system construction: requires a
//     MobilizedBody::Reverse equivalent in RobotModel, which molmodel never uses
//     (CompoundSystem builds forward only) and the port does not implement.
// ---------------------------------------------------------------------------