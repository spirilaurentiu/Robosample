// ============================================================================
//  RobotEngine_reaction.cpp - mobilizer reaction forces (port of
//  calcMobilizerReactionForces).
//
//  One of four cohesive splits of the former RobotEngine.cpp (SPLIT-R3, pure
//  code motion); see RobotEngine_kinematics.cpp's header comment for the
//  Simbody provenance / validation notes shared by all four
//  RobotEngine_*.cpp translation units.
// ============================================================================

#include "RobotEngine.hpp"

#include <vector>

#include "RobotEngine_internal.hpp"
#include "robot_math.hpp"

using robo::ArticulatedInertia;
using robo::Mat33;
using robo::PhiMatrix;
using robo::Quaternion;
using robo::Real;
using robo::Rotation;
using robo::SpatialInertia;
using robo::SpatialVec;
using robo::SymMat33;
using robo::Transform;
using robo::Vec3;
using robo::Vec4;

// ============================================================================
//  MOBILIZER REACTION FORCES   (port of calcMobilizerReactionForces)
//
//  The reaction transmitted to body b across its inboard mobilizer, in Ground,
//  is rigid Newton-Euler with the TRUE accelerations A_GB plus what the children
//  pass back inward:
//      reac_b@Bo = Mk_b A_GB_b + gyro_b - F_ext_b + sum_c Phi[c] reac_c@Bo
//  F_ext_b is the applied spatial BODY force on b only (bodyForceG from the
//  bridge). Matching Simbody (SimbodyMatterSubsystemRep::calcMobilizerReactionForces
//  / calcMobilizerReactionForcesUsingFreebodyMethod, SimbodyMatterSubsystemRep.cpp
//  :6061-6062,6067-6168): "any generalized forces applied at the mobilities end
//  up included in the reaction forces" -- i.e. applied mobility (generalized
//  joint) forces are NOT subtracted out here; the reported reaction is the one
//  actually transmitted across the joint given whatever generalized force was
//  applied. NOTE: this is behavior-neutral today because the force bridge
//  zeroes mobilityForce every step (include/ForceBridge.hpp:73-75), so mobF is
//  identically 0 and dropping the H*mobF term changes nothing numerically. It
//  becomes live once a Fixman/biasing generalized torque is introduced --  at
//  that point this convention (reaction includes actuation) is the intended one.
//  Phi[c] (offset = parent->child origin in Ground) shifts a child's force from
//  the child origin to b's origin -- identical to calcUDot's pass-1 transmission.
// ============================================================================
void RobotEngine::calcMobilizerReactionForces(const RobotModel& m,
                                              const RobotState& s,
                                              SpatialVec* reactionAtBoInG,
                                              SpatialVec* reactionAtMInG) {
    const SpatialInertia* Mk = s.Mk_G();
    const SpatialVec* A_GB = s.A_GB();
    const SpatialVec* gyro = s.gyro();
    const SpatialVec* bodyF = s.bodyForceG();
    const PhiMatrix* Phi = s.Phi();
    const Transform* X_GB = s.X_GB();

    // reaction at body origin, accumulated inward. Local scratch so the operator
    // is side-effect free on the cache (callers may not want a dedicated slot).
    std::vector<SpatialVec> reacBo(static_cast<std::size_t>(m.numBodies), SpatialVec(Vec3(0), Vec3(0)));

    for (int b = m.numBodies - 1; b >= 1; --b) {
        // applied spatial force on b: bridge body force only (see Simbody
        // convention note above -- generalized/mobility forces are NOT
        // subtracted here).
        const SpatialVec& fExt = bodyF[b];

        // rigid Newton-Euler residual at Bo, in Ground.
        SpatialVec reac = (Mk[b] * A_GB[b]) + gyro[b] - fExt;

        // add what the outboard children transmit inward (shift child Bo -> b Bo).
        for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
            const int c = m.bodyChildren[ci];
            reac += Phi[c] * reacBo[static_cast<std::size_t>(c)];
        }

        reacBo[static_cast<std::size_t>(b)] = reac;
        if (reactionAtBoInG != nullptr) {
            reactionAtBoInG[b] = reac;
        }
        if (reactionAtMInG != nullptr) {
            // shift the spatial force from Bo to the outboard frame origin Mo:
            //   p_BoMo_G = R_GB * X_BM.p ;  [t;f]@Bo -> [t - p x f ; f]@Mo.
            const Vec3 p_BoMo_G = X_GB[b].R() * m.X_BM[b].p();
            reactionAtMInG[b] = SpatialVec(reac.angular - (p_BoMo_G % reac.linear), reac.linear);
        }
    }

    if (reactionAtBoInG != nullptr) {
        reactionAtBoInG[0] = SpatialVec(Vec3(0), Vec3(0)); // Ground: no inboard joint
    }
    if (reactionAtMInG != nullptr) {
        reactionAtMInG[0] = SpatialVec(Vec3(0), Vec3(0));
    }
}

auto RobotEngine::findMobilizerReactionOnBodyAtMInGround(const RobotModel& m, const RobotState& s, int body)
    -> SpatialVec {
    std::vector<SpatialVec> atM(static_cast<std::size_t>(m.numBodies), SpatialVec(Vec3(0), Vec3(0)));
    calcMobilizerReactionForces(m, s, nullptr, atM.data());
    if (body < 0 || body >= m.numBodies) {
        return SpatialVec(Vec3(0), Vec3(0));
    }
    return atM[static_cast<std::size_t>(body)];
}
