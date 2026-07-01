#pragma once

// ============================================================================
//  RobotEngine — the SimTK-free dynamics, living in Robosample. Operates on a
//  (const RobotModel&, RobotState&) pair. NO SimTK subsystem/handle types, NO
//  vtables: mobilizer identity is a `switch (model.bodyJoint[b])`, replacing
//  Simbody's RigidBodyNode virtual dispatch.
//
//  EVERY method here is a *faithful port* of a named Simbody routine. The point
//  is behavior preservation. NOTE: the migration is complete and Simbody has
//  been removed from the build, so the old "build beside Simbody and diff each
//  operator" gate no longer exists (see the VALIDATION note in RobotEngine.cpp
//  for the analytic / finite-difference / statistical oracle-of-record that
//  replaced it). These routines are transcriptions, not inventions, and are
//  validated against physics and against each other, not against a live Simbody.
//
//  Tree order invariant: model.bodyParent[b] < b. Outward sweep = forward
//  iteration over bodies; inward sweep = reverse iteration. Ground == body 0.
//
//  OpenMP-future: all sweeps are index loops over contiguous SoA; level-parallel
//  variants can later use model.bodyLevel without changing the math. No OpenMP
//  now.
// ============================================================================

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
    // Port of RigidBodyNodeSpec::realizePosition (RigidBodyNodeSpec.h) +
    // per-joint calcAcrossJointVelocityJacobian (H_FM) + calcBodyTransforms.
    // Outward sweep. Fills X_FM, X_PB, X_GB, Phi, H_FM, H_PB_G, and then the
    // per-atom Ground positions/stations (port of
    // DuMMForceFieldSubsystemRep::realizeSubsystemPositionImpl).
    static void realizePosition(const RobotModel& m, RobotState& s);

    // Port of RigidBodyNodeSpec::realizeVelocity. Fills V_FM, V_GB, Coriolis
    // (AC_GB precursor), gyroscopic terms. Outward sweep.
    static void realizeVelocity(const RobotModel& m, RobotState& s);

    // q-dot from u: per-joint calcQDot / multiplyByN (quaternion coupling for
    // Ball/Free lives here). Port of RigidBodyNodeSpec::calcQDot.
    static void calcQDot(const RobotModel& m, const RobotState& s, robo::Real* qdotOut);
    static void calcQDotDot(const RobotModel& m, RobotState& s);

    // ---- ARTICULATED-BODY DYNAMICS ----------------------------------------
    // Port of realizeArticulatedBodyInertiasInward (RigidBodyNodeSpec.cpp).
    // Inward sweep. Fills P, PPlus, D, DI, G.
    static void realizeArticulatedBodyInertias(const RobotModel& m, RobotState& s);

    // Forward dynamics: forces -> udot. Port of calcUDotPass1Inward (Z, eps,
    // zPlus) + calcUDotPass2Outward (udot, A_GB). Requires bodyForceG +
    // mobilityForce already set (ForceBridge + optional Fixman).
    static void calcUDot(const RobotModel& m, RobotState& s);

    // ---- MASS-MATRIX OPERATORS (no M ever formed) -------------------------
    // Port of multiplyByMInvPass1Inward + multiplyByMInvPass2Outward.
    static void multiplyByMInv(const RobotModel& m, RobotState& s, const robo::Real* in, robo::Real* out);

    // Port of multiplyBySqrtMInvPassOutward (+ inward companion). Used for HMC
    // velocity seeding: u = sqrt(boostRT) * sqrtMInv * gaussian.
    static void multiplyBySqrtMInv(const RobotModel& m, RobotState& s, const robo::Real* in, robo::Real* out);

    // sqrt(M): exact algebraic inverse of multiplyBySqrtMInv (see RobotEngine.cpp).
    // Used by the NMA Route B acceptance to recover the white-noise coordinates
    // w = M^(1/2) u / sqrt(RT) of a trajectory-end velocity. Uses local scratch,
    // so it does NOT disturb s.V_GB() (safe to call alongside calcKineticEnergy).
    static void multiplyBySqrtM(const RobotModel& m, RobotState& s, const robo::Real* in, robo::Real* out);

    // ---- MASS-MATRIX LOG-DETERMINANT (Fixman) -----------------------------
    // ln|M_phi| = sum_b ln det(D_b), the O(n) articulated-body determinant
    // (port of RigidBodyNodeSpec::calcDetMPass2Outward, which accumulates
    //  detM += log(det(D))). Matrix-free: the global nu x nu mass matrix is
    //  NEVER formed; only the per-body dof x dof hinge inertia D_b = ~H P H
    //  (dof <= 6) is built and its determinant taken. REQUIRES
    //  realizeArticulatedBodyInertias to have run first (reads P and H).
    static robo::Real calcLogDetM(const RobotModel& m, const RobotState& s);


    // Port of SimbodyMatterSubsystemRep::calcKineticEnergy ( 1/2 u^T M u ).
    static robo::Real calcKineticEnergy(const RobotModel& m, const RobotState& s);

    // ---- MOBILIZER REACTION FORCES ----------------------------------------
    // Port of SimbodyMatterSubsystemRep::calcMobilizerReactionForces. The
    // reaction on body b is the spatial force its INBOARD mobilizer transmits to
    // it, expressed in Ground. It is obtained by a rigid Newton-Euler inward
    // sweep using the TRUE body accelerations A_GB (so calcUDot must be current):
    //   reac_b@Bo = Mk_b A_GB_b + gyro_b - F_ext_b + sum_children Phi[c] reac_c@Bo
    // where F_ext_b = bodyForceG[b] ONLY. Applied mobility (generalized joint)
    // forces are NOT subtracted -- matching Simbody, they end up included in the
    // reported reaction (SimbodyMatterSubsystemRep.cpp:6061-6062); see the
    // definition's doc comment in RobotEngine.cpp for the full note.
    // No articulated inertia is used; this is exact rigid force transmission.
    //
    // reactionAtBoInG[b] (length numBodies, [0]=Ground unused) receives the
    // reaction reported AT THE BODY ORIGIN Bo in Ground. Pass nullptr to skip.
    // reactionAtMInG[b] receives the SAME reaction shifted to the outboard
    // mobilizer frame origin Mo in Ground (Simbody's findMobilizerReactionOn-
    // BodyAtMInGround convention): [t;f] at Bo -> [t - p_BoMo_G x f; f] at Mo.
    // Pass nullptr to skip. At least one output must be non-null.
    //
    // PRECONDITION: realizePosition, realizeVelocity, realizeArticulatedBody-
    // Inertias, and calcUDot all current (reads Mk_G, A_GB, gyro, Phi, H,
    // bodyForceG, mobilityForce). Ground (body 0) has no inboard joint; its slot
    // is left zero.
    static void calcMobilizerReactionForces(const RobotModel& m,
                                            const RobotState& s,
                                            robo::SpatialVec* reactionAtBoInG,
                                            robo::SpatialVec* reactionAtMInG);

    // Convenience: the reaction on body b at its M frame origin, in Ground.
    static robo::SpatialVec
    findMobilizerReactionOnBodyAtMInGround(const RobotModel& m, const RobotState& s, int body);

    // ---- INTEGRATOR -------------------------------------------------------
    // Port of VerletIntegrator.cpp::attemptDAEStep, FIXED step. Sequence:
    //   q1 = q0 + h*qdot0 + (h^2/2)*qdotdot0
    //   u1_est = u0 + h*udot0
    //   realizePosition; projectQ (== quaternion renormalization ONLY, no
    //     constraints — see MIGRATION_DESIGN.md §1/§4.3)
    //   realizeVelocity; forces; calcUDot
    //   implicit-trapezoid refine u (<=10 iters): u1 = u0 + (h/2)(udot0+udot1)
    // Advances state by exactly one step h. Returns false if a non-finite force
    // or coordinate appeared during the step (e.g. a hard steric overlap drove an
    // LJ force to Inf); in that case the pre-step q/u are restored so the caller
    // sees a finite, unmodified state to reject from -- never a NaN it must chase.
    template <class Bridge>
    static bool verletStep(const RobotModel& m,
                           RobotState& s,
                           Bridge& bridge,
                           const robo::ConstraintSet& cset,
                           robo::Real h);

    // Drives verletStep until t_end. Port of TimeStepper::stepTo for the fixed-
    // step velocity-Verlet path Robosample uses. Returns success.
    template <class Bridge>
    static bool stepTo(const RobotModel& m,
                       RobotState& s,
                       Bridge& bridge,
                       const robo::ConstraintSet& cset,
                       robo::Real tEnd);

    // Reversibility diagnostic: integrate nSteps forward at step h, flip the
    // momenta, integrate nSteps back, flip again, and return the RELATIVE
    // round-trip residual ||(q,u)_returned - (q,u)_start|| / ||(q,u)_start||.
    // ~machine-eps for a reversible step; O(1) when h is too large. Quaternion
    // DOF use a double-cover-aware distance. NON-DESTRUCTIVE (restores s).
    // NOTE: certifies h only for the configuration it is called from -- the safe
    // h is configuration dependent (M(q), force stiffness), so this is a startup/
    // periodic smoke test, not a whole-run guarantee. The per-step corrector
    // throw in verletStep is the ongoing guard. See THEORY 5.5.
    template <class Bridge>
    static robo::Real checkReversibility(const RobotModel& m,
                                         RobotState& s,
                                         Bridge& bridge,
                                         const robo::ConstraintSet& cset,
                                         int nSteps,
                                         robo::Real h);

    // ---- TRANSFER: internal q/u -> Cartesian (accept path) ----------------
    // Port of the accept-branch loop in HMCSampler::sampleIteration:
    //   atomPosG[a] = X_GB[body].p + X_GB[body].R * station_B[a]
    // Assumes realizePosition has run.
    static void fillAtomPositionsFromBodies(const RobotModel& m, RobotState& s);

    // ---- TRANSFER: Cartesian targets -> internal state --------------------
    // ---- quaternion projection (no-constraint Verlet projection) ----------
    // Renormalizes each quaternion q-block listed in model.quaternionQStart.
    static void normalizeQuaternions(const RobotModel& m, RobotState& s);

    // NOTE: the per-joint kernels (X_FM, H_FM, qdot, qddot) are transcribed
    // inline inside realizePosition / realizeVelocity / calcQDot[Dot] rather
    // than as separate static helpers, so there are no extra symbols to define.
};