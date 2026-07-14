#pragma once

/**
 * @file RobotEngine.hpp
 * @brief The SimTK-free articulated-body dynamics: kinematics, articulated-body
 *        inertia recursion, forward dynamics, the matrix-free mass operators, the
 *        mass-matrix log-determinant, and reaction forces.
 *
 * A namespace of static methods over an immutable @c RobotModel and a mutable
 * @c RobotState (a structure-of-arrays driver, not an object). The engine owns no
 * state: it reads the model and writes cache fields on the state. Mobilizer
 * identity is a @c switch on @c model.bodyJoint[b], replacing virtual dispatch.
 *
 * Realization stages run in a fixed order (INV-4): position -> velocity ->
 * articulated-body inertias -> forward dynamics (@c udot). Each method below
 * documents which cache stages it requires valid on entry and which it makes
 * valid on exit; these entry/exit contracts are the write-side counterpart of
 * the @c RobotState accessor contracts and SHALL agree with them.
 * @note Tree order: @c model.bodyParent[b] < b, so an outward sweep is forward
 *       iteration and an inward sweep is reverse iteration; Ground is body 0.
 *       Every method is a behavior-preserving port of a named Simbody routine.
 */

#include "ForceBridge.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "robot_math.hpp"

// Forward declaration breaks the Constraints <-> RobotEngine include cycle.
// verletStep takes ConstraintSet only by reference, so an incomplete type is
// enough here; the full definition is pulled in by RobotEngine.cpp.
namespace robo {
class ConstraintSet;
}

class RobotEngine {
    public:
    // ---- KINEMATICS --------------------------------------------------------
    /**
     * @brief Realize the position stage: body transforms, the velocity Jacobian,
     *        and per-atom Ground positions from the generalized coordinates.
     * @param[in]     m Immutable model.
     * @param[in,out] s State; reads @c q, writes the position caches
     *                  (@c X_FM, @c X_PB, @c X_GB, @c Phi, @c H_FM, @c H_PB_G) and
     *                  the per-atom Ground positions/stations.
     * @pre @c q is set for @p s.
     * @post The position stage is valid (INV-4); this is the first stage every
     *       downstream recursion depends on.
     */
    static void realizePosition(const RobotModel& m, RobotState& s);

    /**
     * @brief Realize the velocity stage: cross-mobilizer and body spatial
     *        velocities, Coriolis and gyroscopic terms.
     * @param[in]     m Immutable model.
     * @param[in,out] s State; reads @c u and the position caches, writes
     *                  @c V_FM, @c V_GB, the Coriolis (@c AC_GB precursor), and
     *                  gyroscopic terms.
     * @pre @ref realizePosition is current and @c u is set for @p s.
     * @post The velocity stage is valid (INV-4).
     */
    static void realizeVelocity(const RobotModel& m, RobotState& s);

    /**
     * @brief Coordinate rates @c qdot from generalized speeds @c u (with the
     *        quaternion coupling for Ball/Free/FreeLine).
     * @param[in]  m       Immutable model.
     * @param[in]  s       State; reads @c q, @c u, and the position caches.
     * @param[out] qdotOut Coordinate rates, length @c m.nq.
     * @pre @ref realizePosition is current for @p s.
     */
    static void calcQDot(const RobotModel& m, const RobotState& s, robo::Real* qdotOut);
    /**
     * @brief Coordinate second derivatives @c qdotdot from accelerations @c udot.
     * @param[in]     m Immutable model.
     * @param[in,out] s State; reads @c q, @c u, @c udot and writes @c qdotdot.
     * @pre @ref realizePosition is current and @c udot is set (see @ref calcUDot).
     */
    static void calcQDotDot(const RobotModel& m, RobotState& s);

    // ---- ARTICULATED-BODY DYNAMICS ----------------------------------------
    /**
     * @brief Factorize the articulated-body inertias (position-only part).
     * @param[in]     m Immutable model.
     * @param[in,out] s State; writes the hinge factorization caches @c P,
     *                  @c PPlus, @c D, @c DI, @c G.
     * @pre @ref realizePosition is current.
     * @post The @c q-only inertia factorization is valid; this part alone is
     *       enough for the matrix-free mass operators and @ref calcLogDetM.
     * @note @c DI is the @ref robo::detail::invertDense pseudo-inverse of each
     *       hinge block @c D; a null-locked direction is frozen (INV-5).
     */
    static void factorizeArticulatedInertias(const RobotModel& m, RobotState& s); // P,PPlus,D,DI,G (q-only)
    /**
     * @brief Seed the velocity-dependent centrifugal (Coriolis) articulated force.
     * @param[in]     m Immutable model.
     * @param[in,out] s State; writes @c abcf = P*a_mob + gyro.
     * @pre @ref realizeVelocity and @ref factorizeArticulatedInertias are current.
     */
    static void seedArticulatedCentrifugal(const RobotModel& m, RobotState& s);   // abcf = P*a_mob+gyro (u)
    /**
     * @brief Realize the articulated-body inertia stage
     *        (@ref factorizeArticulatedInertias then @ref seedArticulatedCentrifugal).
     * @param[in]     m Immutable model.
     * @param[in,out] s State; writes @c P, @c PPlus, @c D, @c DI, @c G, @c abcf.
     * @pre @ref realizePosition and @ref realizeVelocity are current.
     * @post The articulated-body inertia stage is valid (INV-4); required before
     *       @ref calcUDot.
     */
    static void realizeArticulatedBodyInertias(const RobotModel& m, RobotState& s); // factorize + seed

    /**
     * @brief Forward dynamics: solve accelerations @c udot from the applied
     *        forces.
     * @param[in]     m Immutable model.
     * @param[in,out] s State; reads the inertia factorization, @c bodyForceG and
     *                  @c mobilityForce, and writes @c udot and the body spatial
     *                  accelerations @c A_GB.
     * @pre @ref realizeArticulatedBodyInertias is current and @c bodyForceG /
     *      @c mobilityForce are set (from the force bridge and any Fixman term).
     * @post @c udot and @c A_GB are valid, completing the realization sequence
     *       (INV-4).
     */
    static void calcUDot(const RobotModel& m, RobotState& s);

    // ---- MASS-MATRIX OPERATORS (no M ever formed) -------------------------
    /**
     * @brief Apply the inverse mass matrix, @c out = M^-1 * in, matrix-free.
     * @param[in]     m   Immutable model.
     * @param[in,out] s   State; uses inertia caches (and velocity scratch).
     * @param[in]     in  Input generalized vector, length @c m.nu.
     * @param[out]    out Result @c M^-1 * in, length @c m.nu.
     * @pre @ref factorizeArticulatedInertias is current (reads @c P, @c DI, @c G).
     * @note The global @c nu x nu mass matrix is never formed (INV-5).
     */
    static void multiplyByMInv(const RobotModel& m, RobotState& s, const robo::Real* in, robo::Real* out);

    /**
     * @brief Apply the inverse mass-matrix square root,
     *        @c out = sqrt(M)^-1 * in, matrix-free.
     * @param[in]     m   Immutable model.
     * @param[in,out] s   State; uses inertia caches (and velocity scratch).
     * @param[in]     in  Input generalized vector, length @c m.nu.
     * @param[out]    out Result @c sqrt(M)^-1 * in, length @c m.nu.
     * @pre @ref factorizeArticulatedInertias is current.
     * @post @c out satisfies the mass metric so that HMC velocity seeding
     *       @c u = sqrt(boostRT) * multiplyBySqrtMInv(gaussian) has covariance
     *       @c boostRT * M^-1; this is the same metric @ref calcKineticEnergy and
     *       @ref calcLogDetM use (INV-5).
     */
    static void multiplyBySqrtMInv(const RobotModel& m, RobotState& s, const robo::Real* in, robo::Real* out);

    /**
     * @brief Apply the mass-matrix square root, @c out = sqrt(M) * in,
     *        matrix-free.
     * @param[in]     m   Immutable model.
     * @param[in,out] s   State; uses local scratch and does NOT disturb
     *                    @c s.V_GB() (safe to call alongside
     *                    @ref calcKineticEnergy).
     * @param[in]     in  Input generalized vector, length @c m.nu.
     * @param[out]    out Result @c sqrt(M) * in, length @c m.nu.
     * @pre @ref factorizeArticulatedInertias is current.
     * @post Exact algebraic inverse of @ref multiplyBySqrtMInv on the unlocked
     *       subspace; the NMA Route-B acceptance uses it to recover the
     *       white-noise coordinates @c w = sqrt(M) u / sqrt(RT) (INV-5).
     */
    static void multiplyBySqrtM(const RobotModel& m, RobotState& s, const robo::Real* in, robo::Real* out);

    // ---- MASS-MATRIX LOG-DETERMINANT (Fixman) -----------------------------
    /**
     * @brief The log-determinant of the mass matrix,
     *        @c ln|M| = sum_b ln det(D_b), the O(n) articulated-body determinant.
     * @param[in] m Immutable model.
     * @param[in] s Realized state.
     * @return @c ln|M|, the Fixman kinetic term.
     * @pre @ref realizeArticulatedBodyInertias is current (reads @c P and @c H).
     * @note Matrix-free: only each per-body hinge inertia @c D_b = H^T P H
     *       (@c dof <= 6) is formed; its pseudo-log-determinant
     *       (@ref robo::detail::pseudoLogDet) omits null-locked directions, so
     *       @c ln|M| stays consistent with the forward-dynamics pseudo-inverse
     *       (INV-5).
     */
    static robo::Real calcLogDetM(const RobotModel& m, const RobotState& s);


    /**
     * @brief The kinetic energy @c 1/2 u^T M u.
     * @param[in] m Immutable model.
     * @param[in] s Realized state (reads @c u and the inertia caches).
     * @return The kinetic energy under the same mass metric as the momentum
     *         seeding operators (INV-5).
     * @pre @ref realizeArticulatedBodyInertias is current.
     */
    static robo::Real calcKineticEnergy(const RobotModel& m, const RobotState& s);

    // ---- MOBILIZER REACTION FORCES ----------------------------------------
    /**
     * @brief Reaction forces the inboard mobilizer transmits to each body,
     *        expressed in Ground.
     * @param[in]  m               Immutable model.
     * @param[in]  s               Fully realized state.
     * @param[out] reactionAtBoInG Reaction at each body origin @c Bo in Ground,
     *                             length @c numBodies (slot 0 for Ground is left
     *                             zero); pass @c nullptr to skip.
     * @param[out] reactionAtMInG  The same reaction shifted to each outboard
     *                             mobilizer frame origin @c Mo in Ground
     *                             (@c [t;f] at Bo -> @c [t - p_BoMo_G x f; f] at
     *                             Mo); pass @c nullptr to skip.
     * @pre @ref realizePosition, @ref realizeVelocity,
     *      @ref realizeArticulatedBodyInertias, and @ref calcUDot are all current.
     *      At least one output pointer is non-null.
     * @note Exact rigid Newton-Euler inward transmission (no articulated
     *       inertia). Following Simbody, applied mobility (generalized joint)
     *       forces are NOT subtracted, so they are included in the reported
     *       reaction; only @c bodyForceG counts as the external force.
     */
    static void calcMobilizerReactionForces(const RobotModel& m,
                                            const RobotState& s,
                                            robo::SpatialVec* reactionAtBoInG,
                                            robo::SpatialVec* reactionAtMInG);

    /**
     * @brief The reaction on one body at its M-frame origin, in Ground.
     * @param[in] m    Immutable model.
     * @param[in] s    Fully realized state (same preconditions as
     *                 @ref calcMobilizerReactionForces).
     * @param[in] body Body index.
     * @return The reaction spatial force at @c Mo in Ground.
     */
    static robo::SpatialVec
    findMobilizerReactionOnBodyAtMInGround(const RobotModel& m, const RobotState& s, int body);

    // ---- INTEGRATOR -------------------------------------------------------
    /**
     * @brief Advance the state by exactly one fixed-step leapfrog HMC step in
     *        generalized coordinates.
     * @tparam Bridge Force provider; the integrator touches it only through
     *         @c bridge.evaluate(s), which fills @c bodyForceG / @c mobilityForce
     *         and reports whether the forces are finite. Any type meeting that
     *         interface may instantiate this (the production @c ForceBridge and
     *         the test @c AnalyticForceBridge both do); OpenMM stays out of the
     *         header.
     * @param[in]     m      Immutable model.
     * @param[in,out] s      State advanced by one step @p h.
     * @param[in,out] bridge Force provider (see @p Bridge).
     * @param[in]     cset   Loop-closure constraints; the SHAKE position
     *                       projection and the trapezoid velocity correction use
     *                       the same @c G as the Fixman log-det (INV-6).
     * @param[in]     h      Step size.
     * @param[out]    correctorConverged Optional; if non-null, receives whether
     *                       the implicit-trapezoid velocity corrector reached its
     *                       fixed point. The step is taken unconditionally either
     *                       way; a caller that needs F-reversibility (e.g. inner
     *                       NCMC) SHALL treat a @c false readback as a reject.
     * @return @c true on success; @c false if a non-finite force or coordinate
     *         appeared during the step, in which case the pre-step @c q / @c u are
     *         restored so the caller sees a finite, unmodified state to reject
     *         from.
     * @post One reversible leapfrog step: position drift (quaternion DOF via the
     *       exact exp-map, INV-9) with quaternion renormalization, SHAKE position
     *       projection (INV-6), then the implicit-trapezoid velocity correction.
     *       The realization stages run in the canonical order (INV-4). A
     *       behavior-preserving reorder of the internal helper calls does not
     *       change this contract.
     */
    template <class Bridge>
    static bool verletStep(const RobotModel& m,
                           RobotState& s,
                           Bridge& bridge,
                           const robo::ConstraintSet& cset,
                           robo::Real h,
                           bool* correctorConverged = nullptr);

    /**
     * @brief Drive @ref verletStep from the current time to @p tEnd at fixed step.
     * @tparam Bridge Force provider (see @ref verletStep).
     * @param[in]     m      Immutable model.
     * @param[in,out] s      State advanced to @p tEnd.
     * @param[in,out] bridge Force provider.
     * @param[in]     cset   Loop-closure constraints.
     * @param[in]     tEnd   Target time.
     * @param[out]    correctorConverged Optional; see @ref verletStep.
     * @return @c true when every step succeeded; @c false on the first failed
     *         step (with that step's state restored, per @ref verletStep).
     */
    template <class Bridge>
    static bool stepTo(const RobotModel& m,
                       RobotState& s,
                       Bridge& bridge,
                       const robo::ConstraintSet& cset,
                       robo::Real tEnd,
                       bool* correctorConverged = nullptr);

    /**
     * @brief Reversibility diagnostic: forward-integrate, flip momenta,
     *        back-integrate, and report the relative round-trip residual.
     * @tparam Bridge Force provider (see @ref verletStep).
     * @param[in]     m      Immutable model.
     * @param[in,out] s      State; restored to its input value before return
     *                       (non-destructive).
     * @param[in,out] bridge Force provider.
     * @param[in]     cset   Loop-closure constraints.
     * @param[in]     nSteps Steps in each leg.
     * @param[in]     h      Step size to certify.
     * @return @c ||(q,u)_returned - (q,u)_start|| / ||(q,u)_start||: near machine
     *         epsilon for a reversible step, @c O(1) when @p h is too large.
     *         Quaternion DOF use a double-cover-aware distance (INV-9).
     * @note Certifies @p h only for the configuration it is called from (safe
     *       @p h depends on @c M(q) and force stiffness); a startup or periodic
     *       smoke test, not a whole-run guarantee.
     */
    template <class Bridge>
    static robo::Real checkReversibility(const RobotModel& m,
                                         RobotState& s,
                                         Bridge& bridge,
                                         const robo::ConstraintSet& cset,
                                         int nSteps,
                                         robo::Real h);

    // ---- TRANSFER: internal q/u -> Cartesian (accept path) ----------------
    /**
     * @brief Recompute per-atom Ground positions from the realized body frames,
     *        @c atomPosG[a] = X_GB[body].p + X_GB[body].R * station.
     * @param[in]     m Immutable model.
     * @param[in,out] s State; reads @c X_GB, writes @c atomPosG.
     * @pre @ref realizePosition is current.
     */
    static void fillAtomPositionsFromBodies(const RobotModel& m, RobotState& s);

    // ---- TRANSFER: Cartesian targets -> internal state --------------------
    /**
     * @brief Renormalize every quaternion coordinate block to unit length.
     * @param[in]     m Immutable model (supplies @c quaternionQStart).
     * @param[in,out] s State; renormalizes the listed quaternion @c q-blocks in
     *                  place.
     * @post Each quaternion coordinate block is unit; this is the
     *       no-constraint Verlet projection step (INV-9).
     */
    static void normalizeQuaternions(const RobotModel& m, RobotState& s);

    // NOTE: the per-joint kernels (X_FM, H_FM, qdot, qddot) are transcribed
    // inline inside realizePosition / realizeVelocity / calcQDot[Dot] rather
    // than as separate static helpers, so there are no extra symbols to define.
};