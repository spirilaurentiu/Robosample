#include "JointKernels.hpp"

#include <cmath>

namespace robo {

// -------- per-joint cross-mobilizer transform X_FM(q) ----------------------
// Faithful to RigidBodyNodeSpec_{Torsion,Slider,Cylinder,Translation,Ball,
// BendStretch,SphericalCoords,FreeLine,Free}.h / RigidBodyNode_Weld.
auto jointX_FM(JointType jt, const Real* q, int qOff) -> Transform {
    switch (jt) {
        case JointType::Rigid:
            return {}; // identity
        case JointType::Torsion: {
            Rotation R;
            R.setRotationFromAngleAboutZ(q[qOff]);
            return Transform(R, Vec3(0));
        }
        case JointType::Slider:
            // 1 dof: translation along F's x axis (the bond axis stays on x).
            return Transform(Rotation(), Vec3(q[qOff], 0, 0));
        case JointType::Cylinder: {
            // 2 dof: rotation about z (q0) + translation along z (q1).
            Rotation R;
            R.setRotationFromAngleAboutZ(q[qOff]);
            return Transform(R, Vec3(0, 0, q[qOff + 1]));
        }
        case JointType::BendStretch: {
            // 2 dof: rotation about z (q0) + translation along M's x (q1),
            // re-expressed in F. RigidBodyNodeSpec_Derived.cpp::RBNodeBendStretch.
            Rotation R;
            R.setRotationFromAngleAboutZ(q[qOff]);
            return Transform(R, R * Vec3(q[qOff + 1], 0, 0));
        }
        case JointType::Cartesian:
            return Transform(Rotation(), Vec3(q[qOff], q[qOff + 1], q[qOff + 2]));
        case JointType::Ball: {
            // 3 dof rotation; q = [q0 q1 q2 q3] unit quaternion, no translation.
            const Vec4 qv(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            Quaternion quat(qv); // normalizes
            Rotation R;
            R.setRotationFromQuaternion(quat);
            return Transform(R, Vec3(0));
        }
        case JointType::SphericalCoords: {
            // 3 dof BAT: body-fixed Z-Y sequence (azimuth q0 about Fz, zenith q1
            // about the new My) + radius q2 along the RADIAL axis. We implement
            // the molmodel default (BondMobility::Spherical): radial axis = Mz,
            // zero offsets, no negation. RigidBodyNodeSpec_SphericalCoords.h.
            // (OrthoSpherical / radial-X would set axisT = Mx; not represented by
            // the single SphericalCoords value -- see RobotModel note.)
            Rotation Rz;
            Rz.setRotationFromAngleAboutAxis(q[qOff], robo::ZAxis);
            Rotation Ry;
            Ry.setRotationFromAngleAboutAxis(q[qOff + 1], robo::YAxis);
            const Rotation R = Rz * Ry; // body-fixed 3-2
            const Vec3 Mz_F = R * Vec3(0, 0, 1);
            return Transform(R, Mz_F * q[qOff + 2]);
        }
        case JointType::FreeLine: {
            // 5 dof: orientation as unit quaternion q[0..3] (R_FM), translation
            // q[4..6]. The missing 6th dofs (Torsion about the body's own line, M's
            // z) lives only in the velocity Jacobian, not here.
            const Vec4 qv(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            Quaternion quat(qv);
            Rotation R;
            R.setRotationFromQuaternion(quat);
            return Transform(R, Vec3(q[qOff + 4], q[qOff + 5], q[qOff + 6]));
        }
        case JointType::Free: {
            // q = [q0 q1 q2 q3 x y z], quaternion (normalized) + translation.
            const Vec4 qv(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            Quaternion quat(qv); // normalizes
            Rotation R;
            R.setRotationFromQuaternion(quat);
            return Transform(R, Vec3(q[qOff + 4], q[qOff + 5], q[qOff + 6]));
        }
        default:
            return Transform();
    }
}

// -------- per-joint H_FM columns (in F), written into Hcol[0..dof-1] --------
// X_FM is needed only by the three q-dependent joints (BendStretch,
// SphericalCoords, FreeLine); the other seven ignore it (their H_FM is constant
// in F). Mirrors each node's calcAcrossJointVelocityJacobian.
void jointH_FM(JointType jt, const Transform& X_FM, SpatialVec* Hcol) {
    switch (jt) {
        case JointType::Rigid:
            break;
        case JointType::Torsion:
            Hcol[0] = SpatialVec(Vec3(0, 0, 1), Vec3(0));
            break;
        case JointType::Slider:
            Hcol[0] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            break;
        case JointType::Cylinder:
            Hcol[0] = SpatialVec(Vec3(0, 0, 1), Vec3(0)); // rotation about z
            Hcol[1] = SpatialVec(Vec3(0), Vec3(0, 0, 1)); // translation along z
            break;
        case JointType::Cartesian:
            Hcol[0] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            Hcol[1] = SpatialVec(Vec3(0), Vec3(0, 1, 0));
            Hcol[2] = SpatialVec(Vec3(0), Vec3(0, 0, 1));
            break;
        case JointType::Ball:
            Hcol[0] = SpatialVec(Vec3(1, 0, 0), Vec3(0));
            Hcol[1] = SpatialVec(Vec3(0, 1, 0), Vec3(0));
            Hcol[2] = SpatialVec(Vec3(0, 0, 1), Vec3(0));
            break;
        case JointType::BendStretch: {
            // u0: rotation about Fz; u1: translation along M's x (in F).
            const Vec3 p_FM = X_FM.p();
            const Vec3 Mx_F = X_FM.R() * Vec3(1, 0, 0);
            Hcol[0] = SpatialVec(Vec3(0, 0, 1), Vec3(0, 0, 1) % p_FM);
            Hcol[1] = SpatialVec(Vec3(0), Mx_F);
            break;
        }
        case JointType::SphericalCoords: {
            // Z-radial default: u0 about Fz, u1 about My, u2 along Mz.
            const Rotation& R = X_FM.R();
            const Vec3 p_FM = X_FM.p();
            const Vec3 sFz(0, 0, 1);
            const Vec3 sMy = R * Vec3(0, 1, 0);
            const Vec3 sMt = R * Vec3(0, 0, 1); // radial axis = Mz
            Hcol[0] = SpatialVec(sFz, sFz % p_FM);
            Hcol[1] = SpatialVec(sMy, sMy % p_FM);
            Hcol[2] = SpatialVec(Vec3(0), sMt);
            break;
        }
        case JointType::FreeLine: {
            // u0,u1: (x,y) angular velocity of M in F, expressed in M (so the
            // columns are M's x/y axes in F); u2..4: translation in F.
            const Rotation& R = X_FM.R();
            const Vec3 Mx_F = R * Vec3(1, 0, 0);
            const Vec3 My_F = R * Vec3(0, 1, 0);
            Hcol[0] = SpatialVec(Mx_F, Vec3(0));
            Hcol[1] = SpatialVec(My_F, Vec3(0));
            Hcol[2] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            Hcol[3] = SpatialVec(Vec3(0), Vec3(0, 1, 0));
            Hcol[4] = SpatialVec(Vec3(0), Vec3(0, 0, 1));
            break;
        }
        case JointType::Free:
            Hcol[0] = SpatialVec(Vec3(1, 0, 0), Vec3(0));
            Hcol[1] = SpatialVec(Vec3(0, 1, 0), Vec3(0));
            Hcol[2] = SpatialVec(Vec3(0, 0, 1), Vec3(0));
            Hcol[3] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            Hcol[4] = SpatialVec(Vec3(0), Vec3(0, 1, 0));
            Hcol[5] = SpatialVec(Vec3(0), Vec3(0, 0, 1));
            break;
        default:
            break;
    }
}

// -------- per-joint dH_FM/dt columns (in F), written into Hdot[0..dof-1] -----
// Zero for the seven joints with H_FM constant in F (jointHasConstantHFM); the
// engine never calls this for them. For BendStretch / SphericalCoords / FreeLine
// these are q,u-dependent (the M axes rotate), ported verbatim from each node's
// calcAcrossJointVelocityJacobianDot. V_FM is the cross-mobilizer spatial
// velocity in F: V_FM[0] = w_FM, V_FM[1] = v_FM.
//
// VALIDATION: these three feed the mobilizer-bias (Coriolis) acceleration, the
// same hot path whose HDot/Coriolis terms have historically hidden silent
// energy-pumping bugs (see realizeVelocity). They are golden-tested by three
// independent oracles (OQ-2 closed, 2026-07-13):
//   1. Simbody differential -- TestRoboticsOracle {BendStretch, SphericalCoords,
//      FreeLine} (+ FuzzStates) compare udot and A_GB (both HDot-dependent)
//      against Simbody fixtures at 1e-8.
//   2. FD kinematic consistency -- JointJacobianDot.HDotMatchesDerivativeOfHFM
//      checks HDot_FM == d/dt jointH_FM along the q-flow, all three types.
//   3. Energy conservation -- Integrator.UnderValidatedJointsConserveEnergy and
//      Integrator.BendStretchIsolatedConservesEnergy integrate a free body/chain.
void jointHDot_FM(JointType jt, const Transform& X_FM, const SpatialVec& V_FM, SpatialVec* Hdot) {
    const Rotation& R = X_FM.R();
    const Vec3 p_FM = X_FM.p();
    const Vec3& w_FM = V_FM[0];
    const Vec3& v_FM = V_FM[1];
    switch (jt) {
        case JointType::BendStretch:
            Hdot[0] = SpatialVec(Vec3(0), Vec3(0, 0, 1) % v_FM);
            Hdot[1] = SpatialVec(Vec3(0), w_FM % (R * Vec3(1, 0, 0)));
            break;
        case JointType::SphericalCoords: {
            const Vec3 sFz(0, 0, 1);
            const Vec3 sMy = R * Vec3(0, 1, 0);
            const Vec3 sMt = R * Vec3(0, 0, 1);
            const Vec3 dsMy = w_FM % sMy;
            const Vec3 dsMt = w_FM % sMt;
            Hdot[0] = SpatialVec(Vec3(0), sFz % v_FM);
            Hdot[1] = SpatialVec(dsMy, (dsMy % p_FM) + (sMy % v_FM));
            Hdot[2] = SpatialVec(Vec3(0), dsMt);
            break;
        }
        case JointType::FreeLine:
            Hdot[0] = SpatialVec(w_FM % (R * Vec3(1, 0, 0)), Vec3(0));
            Hdot[1] = SpatialVec(w_FM % (R * Vec3(0, 1, 0)), Vec3(0));
            Hdot[2] = SpatialVec(Vec3(0), Vec3(0));
            Hdot[3] = SpatialVec(Vec3(0), Vec3(0));
            Hdot[4] = SpatialVec(Vec3(0), Vec3(0));
            break;
        default:
            break;
    }
}

// -------- per-body qdot from u (quaternion coupling for Ball/Free/FreeLine) --
// Centralizes the per-joint-type switch previously duplicated in
// RobotEngine::calcQDot (SPLIT-R4). X_FM is the body's own (already-indexed)
// cross-mobilizer transform; dof/qOff/uOff are the body's own offsets. The
// Weld/Torsion/Slider/Cylinder/Translation/BendStretch/SphericalCoords default
// case is qdot == u.
void jointQDot(JointType jt, const Real* q, int qOff, const Real* u, int uOff, int dof,
               const Transform& X_FM, Real* qdotOut) {
    switch (jt) {
        case JointType::Ball: {
            // 3 dof rotation: quaternion qdot = N(q) * w_FM. No translation.
            const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            const Vec3 w_FM(u[uOff], u[uOff + 1], u[uOff + 2]);
            const Vec4 qd = Rotation::convertAngVelToQuaternionDot(quat, w_FM);
            qdotOut[qOff] = qd[0];
            qdotOut[qOff + 1] = qd[1];
            qdotOut[qOff + 2] = qd[2];
            qdotOut[qOff + 3] = qd[3];
            break;
        }
        case JointType::FreeLine: {
            // 2 rotational speeds are (x,y) of w_FM expressed in M, so
            // w_FM (in F) = R_FM * (u0, u1, 0); qdot_quat = N(q) * w_FM.
            // The 3 translational speeds (u2..4) are qdot of x,y,z directly.
            const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            const Rotation& R_FM = X_FM.R();
            const Vec3 w_FM = R_FM * Vec3(u[uOff], u[uOff + 1], 0);
            const Vec4 qd = Rotation::convertAngVelToQuaternionDot(quat, w_FM);
            qdotOut[qOff] = qd[0];
            qdotOut[qOff + 1] = qd[1];
            qdotOut[qOff + 2] = qd[2];
            qdotOut[qOff + 3] = qd[3];
            qdotOut[qOff + 4] = u[uOff + 2];
            qdotOut[qOff + 5] = u[uOff + 3];
            qdotOut[qOff + 6] = u[uOff + 4];
            break;
        }
        case JointType::Free: {
            // quaternion qdot = N(q)*w ; translation qdot = v.
            const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            const Vec3 w_FM(u[uOff], u[uOff + 1], u[uOff + 2]);
            const Vec4 qd = Rotation::convertAngVelToQuaternionDot(quat, w_FM);
            qdotOut[qOff] = qd[0];
            qdotOut[qOff + 1] = qd[1];
            qdotOut[qOff + 2] = qd[2];
            qdotOut[qOff + 3] = qd[3];
            qdotOut[qOff + 4] = u[uOff + 3];
            qdotOut[qOff + 5] = u[uOff + 4];
            qdotOut[qOff + 6] = u[uOff + 5];
            break;
        }
        default: // Weld/Torsion/Slider/Cylinder/Translation/BendStretch/SphericalCoords: qdot == u
            for (int j = 0; j < dof; ++j) {
                qdotOut[qOff + j] = u[uOff + j];
            }
            break;
    }
}

// -------- per-body qddot from udot (quaternion second derivative) -----------
// Centralizes the per-joint-type branch previously duplicated in
// RobotEngine::calcQDotDot (SPLIT-R4).
void jointQDotDot(JointType jt,
                  const Real* q,
                  int qOff,
                  const Real* u,
                  const Real* udot,
                  int uOff,
                  int dof,
                  const Transform& X_FM,
                  Real* qddOut) {
    if (jt == JointType::Ball) {
        const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
        const Vec3 w(u[uOff], u[uOff + 1], u[uOff + 2]);
        const Vec3 wd(udot[uOff], udot[uOff + 1], udot[uOff + 2]);
        const Vec4 qdd4 = Rotation::convertAngVelDotToQuaternionDotDot(quat, w, wd);
        for (int k = 0; k < 4; ++k) {
            qddOut[qOff + k] = qdd4[k];
        }
    } else if (jt == JointType::FreeLine) {
        // w_FM = R_FM*(u0,u1,0); wdot_FM = R_FM*(udot0,udot1,0) (the
        // R_FM_dot*(u0,u1,0) term is w_FM x w_FM = 0, see derivation).
        const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
        const Rotation& R_FM = X_FM.R();
        const Vec3 w = R_FM * Vec3(u[uOff], u[uOff + 1], 0);
        const Vec3 wd = R_FM * Vec3(udot[uOff], udot[uOff + 1], 0);
        const Vec4 qdd4 = Rotation::convertAngVelDotToQuaternionDotDot(quat, w, wd);
        for (int k = 0; k < 4; ++k) {
            qddOut[qOff + k] = qdd4[k];
        }
        qddOut[qOff + 4] = udot[uOff + 2];
        qddOut[qOff + 5] = udot[uOff + 3];
        qddOut[qOff + 6] = udot[uOff + 4];
    } else if (jt == JointType::Free) {
        const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
        const Vec3 w(u[uOff], u[uOff + 1], u[uOff + 2]);
        const Vec3 wd(udot[uOff], udot[uOff + 1], udot[uOff + 2]);
        const Vec4 qdd4 = Rotation::convertAngVelDotToQuaternionDotDot(quat, w, wd);
        for (int k = 0; k < 4; ++k) {
            qddOut[qOff + k] = qdd4[k];
        }
        qddOut[qOff + 4] = udot[uOff + 3];
        qddOut[qOff + 5] = udot[uOff + 4];
        qddOut[qOff + 6] = udot[uOff + 5];
    } else {
        for (int j = 0; j < dof; ++j) {
            qddOut[qOff + j] = udot[uOff + j];
        }
    }
}

// -------- exact unit-quaternion advance (exponential map) -------------------
// Advance a unit quaternion under a constant angular velocity w_F expressed in
// the PARENT (F) frame, consistent with the engine's qdot = 1/2 (0,w_F) (x) q
// (left multiply; matches convertAngVelToQuaternionDot / N(q)). With theta =
// 1/2 |w| h: Dq = (cos theta, sin theta * w_hat) and q1 = Dq (x) q0. Both
// operands are unit, so q1 is unit BY CONSTRUCTION -- the |q|^2 -> Inf ->
// q/sqrt(Inf) = 0 overflow path of the linear "q += h*qdot" drift cannot occur.
// Reversible: negating w gives Dq^-1, so q0 = Dq^-1 (x) q1 (HMC needs this).
// As h -> 0 it reduces to q0 + h * (1/2 (0,w) (x) q0) = q0 + h * qdot0.
void advanceQuatExp(const Real* q0, const Vec3& wF, Real h, Real* q1) {
    const Real wn = std::sqrt(wF[0] * wF[0] + wF[1] * wF[1] + wF[2] * wF[2]);
    Real a, b, c, d; // Dq = (a, b, c, d)
    if (wn > Real(1e-12)) {
        const Real theta = Real(0.5) * wn * h;
        const Real ssc = std::sin(theta) / wn; // sin(theta) / |w|
        a = std::cos(theta);
        b = ssc * wF[0];
        c = ssc * wF[1];
        d = ssc * wF[2];
    } else {
        a = Real(1); // small-angle limit: Dq ~ (1, 1/2 h w)
        b = Real(0.5) * h * wF[0];
        c = Real(0.5) * h * wF[1];
        d = Real(0.5) * h * wF[2];
    }
    const Real w = q0[0], x = q0[1], y = q0[2], z = q0[3];
    // Hamilton product q1 = Dq (x) q0  (LEFT multiply), consistent with the
    // parent-frame map qdot = 1/2 (0,w) (x) q. (C-3 fix: the body previously
    // computed the right product q0 (x) Dq, which pairs with the body-frame map.)
    q1[0] = a * w - b * x - c * y - d * z;
    q1[1] = a * x + b * w + c * z - d * y;
    q1[2] = a * y + c * w + d * x - b * z;
    q1[3] = a * z + d * w + b * y - c * x;
}

// -------- per-body quaternion position drift (verletStep's driftPositions) --
// wHalf construction + advanceQuatExp call for one quaternion body's midpoint
// angular velocity (SPLIT-R4). qOut receives the 4 new quaternion components;
// the caller (driftPositions) writes them into q[qOff..qOff+3] and keeps the
// surrounding scalar Taylor drift / normalizeQuaternions call unchanged.
void jointDriftQuat(JointType jt,
                    const Real* q0,
                    int qOff,
                    const Real* u0,
                    const Real* udot0,
                    int uOff,
                    const Transform& X_FM,
                    Real h,
                    Real* qOut) {
    // Half-step angular velocity in F. Ball/Free: the first 3 u ARE w_FM.
    // FreeLine: the 2 rotational speeds are (x,y) of w_FM expressed in M, so
    // w_FM = R_FM*(u0,u1,0) (R_FM from the step-start X_FM, still valid here).
    Vec3 wHalf;
    if (jt == JointType::FreeLine) {
        const Rotation& R_FM = X_FM.R();
        wHalf = R_FM
                * Vec3(u0[uOff + 0] + (0.5 * h * udot0[uOff + 0]),
                       u0[uOff + 1] + (0.5 * h * udot0[uOff + 1]),
                       0);
    } else {
        wHalf = Vec3(u0[uOff + 0] + (0.5 * h * udot0[uOff + 0]),
                     u0[uOff + 1] + (0.5 * h * udot0[uOff + 1]),
                     u0[uOff + 2] + (0.5 * h * udot0[uOff + 2])); // w_FM in F, step start
    }
    advanceQuatExp(&q0[qOff], wHalf, h, qOut); // overrides the 4 quaternion slots
}

} // namespace robo