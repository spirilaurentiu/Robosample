#pragma once

// ============================================================================
//  ReactionReporter -- per-frame, per-body reaction-force snapshot row
//  (docs/specs/reaction-force-monitoring.md Sec.1/3/4). Relocated verbatim
//  from World.hpp (SPLIT-W3, pure code motion): `force`/`torque`, about the
//  body origin Bo, in Ground, is the SUM of whichever term(s)
//  World::enableReactionReporter enabled -- the OpenMM net applied force
//  (`bodyForceG`, read directly off `ForceBridge::evaluate`) and/or the
//  static (u=0) mobilizer reaction (`RobotEngine::calcMobilizerReactionForces`).
//  Both terms are spatial forces about the SAME point (Bo, in Ground), so
//  summing them is a valid spatial-force addition, not an apples-to-oranges
//  combination; which term(s) are included is a per-reporter-world CHOICE,
//  not encoded in the row itself (docs/specs/reaction-force-monitoring.md
//  Sec.1.2). RAW (native units, no normalization -- normalization, if
//  wanted, is a render-only transform applied downstream of this CSV, in
//  vmd/arrows.tcl, never here). `bodyIdx` is the RobotModel body index;
//  `atomIdx` is the VMD/DCD atom index (prmtop order,
//  SystemTopology::atomsPrmtopIndex) of a REPRESENTATIVE atom of that body,
//  for placement in vmd/arrows.tcl.
// ============================================================================

#include "robot_math.hpp"

/**
 * @brief One per-frame, per-body reaction-force snapshot row
 *        (docs/specs/reaction-force-monitoring.md Sec.1/3/4).
 *
 * @c force and @c torque are the spatial force about the body origin Bo, in
 * Ground, and are the SUM of whichever term(s) World::enableReactionReporter
 * enabled: the OpenMM net applied force (@c bodyForceG) and/or the static (u=0)
 * mobilizer reaction (RobotEngine::calcMobilizerReactionForces). Both terms are
 * spatial forces about the SAME point (Bo, in Ground), so summing them is a
 * valid spatial-force addition, not an apples-to-oranges combination; which
 * term(s) are included is a per-reporter-world choice, not encoded in the row
 * itself (Sec.1.2).
 *
 * @note Values are RAW (native units, no normalization). Normalization, if
 *       wanted, is a render-only transform applied downstream of this CSV, in
 *       vmd/arrows.tcl, never here.
 */
struct ReactionSample {
    /** @brief RobotModel body index. */
    int bodyIdx;
    /**
     * @brief VMD/DCD atom index (prmtop order, SystemTopology::atomsPrmtopIndex)
     *        of a representative atom of the body, for placement in
     *        vmd/arrows.tcl.
     */
    int atomIdx;
    /** @brief Net force term(s) about Bo, in Ground (native units). */
    robo::Vec3 force;
    /** @brief Net torque term(s) about Bo, in Ground (native units). */
    robo::Vec3 torque;
};
