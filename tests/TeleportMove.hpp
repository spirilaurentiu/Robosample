#pragma once
// ============================================================================
//  TeleportMove.hpp -- the λ=0 trough teleport, at the engine level (the math the
//  production World::ncmcMove restructure mirrors).
//
//  At the λ=0 trough of an NCMC protocol the moved region is a non-interacting
//  ghost, so a rigid reposition costs zero intermolecular energy and can move the
//  solute arbitrarily far. The subtlety (vs the docking RigidKick, which kicks
//  BEFORE momenta are drawn) is that the kick here is mid-trajectory: momenta are
//  live, so the kick must preserve the canonical measure.
//
//  CORRECTNESS (papers cross-check):
//   * Velocity co-rotation. Translation is KE-neutral (M's linear block is m·I,
//     rotation/position-independent). Reorientation is NOT: M_ang(q)=R I_b Rᵀ. A
//     rigid rotation of the WHOLE phase point by ΔR = R_new R_oldᵀ -- rotate the
//     orientation AND the angular+linear generalized speeds -- is a canonical
//     transformation: ½(ΔR u)ᵀ(ΔR M ΔRᵀ)(ΔR u) = ½ uᵀ M u, KE exactly preserved.
//   * Symmetric independence draw (uniform-in-fixed-region translation + Haar
//     rotation) ⇒ g(new|old)=g(old|new) ⇒ folds into the bare endpoint Metropolis
//     (Duane 1987, "correct distribution for any reversible proposal").
//   * For CONSTRAINED (cyclic) regions, a long jump + SHAKE can land on a different
//     Lagrange-multiplier branch than the reverse projection (Brubaker, Salzmann,
//     Urtasun 2012, CHMC Thm 3). reverseShakeReturnsToOrigin() is the branch-
//     consistency guard the move must apply before accepting.
//
//  Operates on the joint quaternion block (Free/Ball root): the orientation it
//  manipulates is R_FM(q); for a root with identity static frames that is the body
//  orientation in Ground, and u's angular/linear parts are in the parent (Ground)
//  frame -- consistent with rotating both by ΔR.
// ============================================================================

#include <algorithm>
#include <cmath>
#include <vector>

#include "Constraints.hpp"
#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "engine_helpers.hpp"
#include "robot_math.hpp"

namespace rtest {

// Orientation rotation R_FM(q) read from body b's quaternion block.
inline robo::Rotation jointRotation(const RobotModel& m, const RobotState& s, int b) {
    const int qOff = m.bodyQIndex[b];
    return EngineHelpers::quatToRotation(s.q()[qOff + 0], s.q()[qOff + 1], s.q()[qOff + 2], s.q()[qOff + 3]);
}

// Rotate a Ground-frame 3-vector by R (left multiply).
inline robo::Vec3 applyR(const robo::Rotation& R, const robo::Vec3& v) {
    return robo::Vec3(R(0, 0) * v[0] + R(0, 1) * v[1] + R(0, 2) * v[2],
                      R(1, 0) * v[0] + R(1, 1) * v[1] + R(1, 2) * v[2],
                      R(2, 0) * v[0] + R(2, 1) * v[1] + R(2, 2) * v[2]);
}

// Rigid teleport of a Free root body b: set translation q to newTransG, set
// orientation to Rnew. With coRotateVel, rotate the body's angular AND linear
// generalized speeds by ΔR = Rnew R_oldᵀ (the canonical, KE-preserving choice).
// Returns ΔR. nu layout for a Free root: [w(3); v(3)]; nq: [quat(4); trans(3)].
inline robo::Rotation teleportFreeRoot(const RobotModel& m,
                                       RobotState& s,
                                       int b,
                                       const robo::Vec3& newTransG,
                                       const robo::Rotation& Rnew,
                                       bool coRotateVel) {
    using robo::Rotation;
    const int qOff = m.bodyQIndex[b];
    const int uOff = m.bodyUIndex[b];

    const Rotation Rold = jointRotation(m, s, b);
    const Rotation dR(Rnew * Rold.transpose()); // ΔR = Rnew R_oldᵀ

    // set orientation quaternion
    robo::Real w, x, y, z;
    EngineHelpers::rotationToQuaternion(Rnew, w, x, y, z);
    s.q()[qOff + 0] = w;
    s.q()[qOff + 1] = x;
    s.q()[qOff + 2] = y;
    s.q()[qOff + 3] = z;
    // set translation
    s.q()[qOff + 4] = newTransG[0];
    s.q()[qOff + 5] = newTransG[1];
    s.q()[qOff + 6] = newTransG[2];

    if (coRotateVel) {
        robo::Real* u = s.u();
        const robo::Vec3 wAng(u[uOff + 0], u[uOff + 1], u[uOff + 2]);
        const robo::Vec3 vLin(u[uOff + 3], u[uOff + 4], u[uOff + 5]);
        const robo::Vec3 wAngR = applyR(dR, wAng);
        const robo::Vec3 vLinR = applyR(dR, vLin);
        u[uOff + 0] = wAngR[0];
        u[uOff + 1] = wAngR[1];
        u[uOff + 2] = wAngR[2];
        u[uOff + 3] = vLinR[0];
        u[uOff + 4] = vLinR[1];
        u[uOff + 5] = vLinR[2];
    }
    return dR;
}

// Incremental (relative) rigid kick: compose dR onto the current orientation and
// add dTransG to the current translation, co-rotating u by dR. Unlike the absolute
// teleportFreeRoot (an independence draw), the incremental kick is a deterministic
// INVOLUTION-compatible map: K(-dT, dR⁻¹) ∘ K(dT, dR) = identity even when
// sandwiched between propagation legs, so it is the right primitive for the
// composite-reversibility check (and for a short-range relative teleport).
inline void teleportFreeRootIncremental(const RobotModel& m,
                                        RobotState& s,
                                        int b,
                                        const robo::Vec3& dTransG,
                                        const robo::Rotation& dR,
                                        bool coRotateVel) {
    using robo::Rotation;
    const int qOff = m.bodyQIndex[b];
    const int uOff = m.bodyUIndex[b];
    const Rotation Rnew(dR * jointRotation(m, s, b));
    robo::Real w, x, y, z;
    EngineHelpers::rotationToQuaternion(Rnew, w, x, y, z);
    s.q()[qOff + 0] = w;
    s.q()[qOff + 1] = x;
    s.q()[qOff + 2] = y;
    s.q()[qOff + 3] = z;
    s.q()[qOff + 4] += dTransG[0];
    s.q()[qOff + 5] += dTransG[1];
    s.q()[qOff + 6] += dTransG[2];
    if (coRotateVel) {
        robo::Real* u = s.u();
        const robo::Vec3 wR = applyR(dR, robo::Vec3(u[uOff + 0], u[uOff + 1], u[uOff + 2]));
        const robo::Vec3 vR = applyR(dR, robo::Vec3(u[uOff + 3], u[uOff + 4], u[uOff + 5]));
        u[uOff + 0] = wR[0];
        u[uOff + 1] = wR[1];
        u[uOff + 2] = wR[2];
        u[uOff + 3] = vR[0];
        u[uOff + 4] = vR[1];
        u[uOff + 5] = vR[2];
    }
}

// Branch-consistency guard for a CONSTRAINED region (CHMC Thm 3). Detailed balance
// of a teleport-then-SHAKE move requires the REVERSE move (apply the inverse kick,
// then SHAKE) to return to the originating manifold point; otherwise the jump
// crossed to a different Lagrange-multiplier branch and must be rejected.
//
// Precondition: s holds q1 = SHAKE(qOrigin + kick) (the forward-proposed point);
// `kick` is the RAW q-space displacement applied before the forward SHAKE.
// Returns true iff SHAKE(q1 − kick) == qOrigin within tol. Restores q1 on exit.
template <typename RefreshFn>
inline bool teleportReverseConsistent(const RobotModel& m,
                                      RobotState& s,
                                      const robo::ConstraintSet& cs,
                                      const std::vector<robo::Real>& qOrigin,
                                      const std::vector<robo::Real>& kick,
                                      RefreshFn&& refresh,
                                      robo::Real tol = robo::Real(1e-5)) {
    if (cs.empty()) {
        return true; // acyclic: no manifold, guaranteed no-op
    }
    std::vector<robo::Real> q1(s.q(), s.q() + m.nq);
    for (int i = 0; i < m.nq; ++i) {
        s.q()[i] = q1[i] - kick[i]; // inverse kick from the proposed point
    }
    refresh();
    cs.enforcePositionConstraints(m, s, refresh);
    robo::Real worst = 0;
    for (int i = 0; i < m.nq; ++i) {
        worst = std::max(worst, std::abs(s.q()[i] - qOrigin[i]));
    }
    std::copy(q1.begin(), q1.end(), s.q()); // restore the forward-proposed q
    refresh();
    return worst <= tol;
}

} // namespace rtest
