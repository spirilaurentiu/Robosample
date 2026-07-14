#pragma once

// ============================================================================
//  ForceReducer -- THE single per-body force->wrench reduction (INV-1/INV-2,
//  see docs/specs/refactor/SPLIT-DEDUP-FORCEREDUCER.md). Hoisted out of
//  ForceBridge::getForcesFromOpenMM (the host path) so the CUDA reduceForces
//  kernel (src/OpenMMContext.cpp) has exactly one host-side definition to be
//  pinned against (tests/TestForceReducer.cpp). OpenMM-free: no OpenMM types
//  appear in this header or its translation unit.
//
//  Convention (restated in full):
//    bodyForceG[b].linear  = sum_a f_a
//    bodyForceG[b].angular = sum_a (r_a - origin_b) x f_a
//  summed over body b's REAL atoms (mass != 0); virtual-site slots are
//  skipped -- OpenMM has already projected their force onto real parent
//  atoms, so re-reducing the leftover slot would double-count it.
// ============================================================================

#include "robot_math.hpp"

/**
 * @brief Reduces per-atom Ground-frame forces into per-body spatial wrenches
 *        about each body's Ground origin, the sole host implementation of the
 *        force->wrench convention (INV-1).
 *
 * The reduction borrows the caller's per-atom arrays read-only and accumulates
 * (`+=`) into @p bodyForceG. It owns nothing and allocates nothing. Real atoms
 * (nonzero mass) belonging to body `b` contribute; massless slots are skipped
 * (INV-2, see @note).
 *
 * @param[in]     atomForceG  Per-atom Cartesian force in Ground, `[numAtoms]`,
 *                            kJ/mol/nm. Borrowed, read-only.
 * @param[in]     atomPosG    Per-atom Cartesian position in Ground, `[numAtoms]`,
 *                            nm. Borrowed, read-only.
 * @param[in]     atomMass    Per-atom mass in Daltons, `[numAtoms]`. A zero entry
 *                            marks a virtual site and is skipped (INV-2).
 * @param[in]     atomBody    Body index of each atom, `[numAtoms]`, each in
 *                            `[0, numBodies)`. Borrowed, read-only.
 * @param[in]     X_GB        Body-to-Ground transforms, `[numBodies]`. Only
 *                            `X_GB[b].p()` (the body origin in Ground) is read.
 * @param[in]     numAtoms    Length of the per-atom arrays.
 * @param[in]     numBodies   Length of @p X_GB and @p bodyForceG.
 * @param[in,out] bodyForceG  Per-body wrench `(angular, linear)`, `[numBodies]`,
 *                            owned and sized by the caller. Accumulated into.
 *
 * @pre  @p bodyForceG is caller-sized to `numBodies` and zero-initialized before
 *       the call; this routine accumulates and never zeros. A non-zeroed buffer
 *       yields wrenches summed on top of stale data.
 * @pre  Every `atomBody[a]` lies in `[0, numBodies)`; out-of-range indices are
 *       not checked and read/write out of bounds.
 *
 * @post For each body `b`, over its real (mass != 0) atoms `a`, with body origin
 *       `origin_b = X_GB[b].p()`:
 *       `bodyForceG[b].linear  += sum_a f_a` and
 *       `bodyForceG[b].angular += sum_a (r_a - origin_b) x f_a`,
 *       i.e. `(torque about the body origin, net force)` in the Ground frame
 *       (INV-1). `SpatialVec` index `[0]` is angular, `[1]` is linear.
 *
 * @note INV-2: a massless slot (`atomMass[a] == 0`, e.g. an OPC/TIP4P M-site) is
 *       skipped by contract, not as an incidental loop guard. OpenMM's
 *       `getState(Forces)` has already projected the site's force onto its real
 *       parent atoms; reducing the leftover force still sitting in the site's
 *       own slot would double-count it and inject a non-conservative kick.
 *
 * @see tests/TestForceReducer.cpp - pins this reducer and the on-device
 *      `reduceForces` kernel against each other on identical input (INV-1/INV-2).
 */
void reduceAtomForcesToBodies(const robo::Vec3* atomForceG,
                              const robo::Vec3* atomPosG,
                              const robo::Real* atomMass,
                              const int* atomBody,
                              const robo::Transform* X_GB,
                              int numAtoms,
                              int numBodies,
                              robo::SpatialVec* bodyForceG);
