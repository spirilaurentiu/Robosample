#pragma once

/**
 * @file JointKernels.hpp
 * @brief The single per-joint-type kernel set: cross-mobilizer transform,
 *        velocity Jacobian and its rate, the q/u <-> qdot/qddot coupling, and
 *        the exact quaternion advance, dispatched by @c JointType.
 *
 * These are stateless free functions; each reads its inputs and writes into a
 * caller-provided output buffer, holding no state. Dispatch on @ref JointType is
 * centralized here (adding a joint touches this file), so the kinematic
 * recursion, the velocity Jacobian, and the position drift stay mutually
 * consistent. The supported joint set is exactly the @ref JointType enum:
 * Rigid, Torsion, Slider, Cylinder, BendStretch, Cartesian, Ball,
 * SphericalCoords, FreeLine, Free.
 */

#include "RobotModel.hpp"
#include "robot_math.hpp"

namespace robo {

/**
 * @brief Cross-mobilizer transform @c X_FM(q) of one body: child frame M in its
 *        parent joint frame F.
 * @param[in] jt   Joint type selecting the kinematic map.
 * @param[in] q    Full generalized-coordinate array.
 * @param[in] qOff Offset of this body's coordinates within @p q.
 * @return @c X_FM for the given @p jt and coordinates. Rigid returns the
 *         identity transform (0 dof).
 */
auto jointX_FM(JointType jt, const Real* q, int qOff) -> Transform;

/**
 * @brief Columns of the cross-mobilizer velocity Jacobian @c H_FM (in F).
 * @param[in]  jt    Joint type.
 * @param[in]  X_FM  This body's cross-mobilizer transform; read only for the
 *                   three q-dependent joints (BendStretch, SphericalCoords,
 *                   FreeLine), ignored by the others (their @c H_FM is constant
 *                   in F).
 * @param[out] Hcol  Jacobian columns written to @c Hcol[0..dof-1].
 * @post @c V_FM == sum_i Hcol[i] * u[i]: @c H_FM maps generalized speeds to the
 *       cross-mobilizer spatial velocity in F.
 */
void jointH_FM(JointType jt, const Transform& X_FM, SpatialVec* Hcol);

/**
 * @brief Time derivative columns @c dH_FM/dt (in F) of the velocity Jacobian.
 * @param[in]  jt   Joint type.
 * @param[in]  X_FM This body's cross-mobilizer transform.
 * @param[in]  V_FM Cross-mobilizer spatial velocity in F (@c V_FM[0] = w_FM,
 *                  @c V_FM[1] = v_FM).
 * @param[out] Hdot Rate columns written to @c Hdot[0..dof-1].
 * @post Identically zero for the seven joints whose @c H_FM is constant in F;
 *       the engine calls this only for BendStretch, SphericalCoords, and
 *       FreeLine, whose M axes rotate. These feed the mobilizer-bias (Coriolis)
 *       acceleration.
 * @note The @c HDot / Coriolis path has historically hidden silent energy-pumping
 *       bugs, so the BendStretch, SphericalCoords, and FreeLine cases are pinned by
 *       three independent oracles (OQ-2 closed): the external Simbody differential
 *       @c TestRoboticsOracle.{BendStretch,SphericalCoords,FreeLine} compares the
 *       HDot-dependent @c udot / @c A_GB at 1e-8; @c JointJacobianDot.HDotMatchesDerivativeOfHFM
 *       checks @c HDot_FM == d/dt @c jointH_FM by finite difference; and
 *       @c Integrator.UnderValidatedJointsConserveEnergy /
 *       @c BendStretchIsolatedConservesEnergy are the energy-conservation oracle.
 */
void jointHDot_FM(JointType jt, const Transform& X_FM, const SpatialVec& V_FM, SpatialVec* Hdot);

/**
 * @brief Per-body @c qdot from generalized speeds @c u (quaternion coupling for
 *        Ball / Free / FreeLine).
 * @param[in]  jt      Joint type.
 * @param[in]  q       Full coordinate array.
 * @param[in]  qOff    This body's coordinate offset.
 * @param[in]  u       Full generalized-speed array.
 * @param[in]  uOff    This body's speed offset.
 * @param[in]  dof     This body's speed count.
 * @param[in]  X_FM    This body's cross-mobilizer transform.
 * @param[out] qdotOut Coordinate rates for this body.
 * @post For quaternion-bearing joints @c qdot is the parent-frame quaternion map
 *       of the angular speeds (@ref quaternionDotFromAngVel); for every other
 *       joint @c qdot == u.
 */
void jointQDot(JointType jt, const Real* q, int qOff, const Real* u, int uOff, int dof,
               const Transform& X_FM, Real* qdotOut);

/**
 * @brief Per-body @c qddot from accelerations @c udot (quaternion second
 *        derivative for the quaternion-bearing joints).
 * @param[in]  jt     Joint type.
 * @param[in]  q      Full coordinate array.
 * @param[in]  qOff   This body's coordinate offset.
 * @param[in]  u      Full generalized-speed array.
 * @param[in]  udot   Full generalized-acceleration array.
 * @param[in]  uOff   This body's speed/acceleration offset.
 * @param[in]  dof    This body's speed count.
 * @param[in]  X_FM   This body's cross-mobilizer transform.
 * @param[out] qddOut Coordinate second derivatives for this body.
 * @post For non-quaternion joints @c qddot == udot.
 */
void jointQDotDot(JointType jt,
                  const Real* q,
                  int qOff,
                  const Real* u,
                  const Real* udot,
                  int uOff,
                  int dof,
                  const Transform& X_FM,
                  Real* qddOut);

/**
 * @brief Advance a unit quaternion by an exponential-map step under a constant
 *        angular velocity.
 * @param[in]  q0 Initial unit quaternion (4 components).
 * @param[in]  wF Angular velocity expressed in the PARENT frame F, consistent
 *                with the engine's @c qdot = 1/2 (0,wF) (x) q left-multiply map.
 * @param[in]  h  Time step.
 * @param[out] q1 Advanced quaternion (4 components).
 * @post @p q1 is unit by construction (the exponential map composes two unit
 *       quaternions), so the overflow-to-zero failure of a linear
 *       @c q += h*qdot drift cannot occur (INV-9).
 * @post Reversible: negating @p wF maps @p q1 back to @p q0, which HMC
 *       reversibility relies on. As @c h -> 0 it reduces to the linear drift
 *       @c q0 + h*qdot0.
 */
void advanceQuatExp(const Real* q0, const Vec3& wF, Real h, Real* q1);

/**
 * @brief Position drift of one quaternion body's orientation over a step, using
 *        its midpoint angular velocity.
 * @param[in]  jt    Joint type (quaternion-bearing: Ball / Free / FreeLine).
 * @param[in]  q0    Initial coordinate array.
 * @param[in]  qOff  This body's coordinate offset.
 * @param[in]  u0    Initial generalized speeds.
 * @param[in]  udot0 Initial accelerations (for the midpoint speed).
 * @param[in]  uOff  This body's speed offset.
 * @param[in]  X_FM  This body's cross-mobilizer transform.
 * @param[in]  h     Time step.
 * @param[out] qOut  The four advanced quaternion components; the caller writes
 *                   them into @c q[qOff..qOff+3].
 * @post @p qOut is the exact @ref advanceQuatExp advance under the midpoint
 *       angular velocity, hence unit and reversible (INV-9).
 */
void jointDriftQuat(JointType jt,
                    const Real* q0,
                    int qOff,
                    const Real* u0,
                    const Real* udot0,
                    int uOff,
                    const Transform& X_FM,
                    Real h,
                    Real* qOut);

} // namespace robo