#pragma once
// ============================================================================
//  engine_helpers.hpp -- small orientation / exp-map helpers used by the engine.
//
//  SINGLE SOURCE OF TRUTH. advanceQuatExp was a file-local static in
//  src/RobotEngine.cpp; rotationToQuaternion / safeLogSineSqr / quatToRotation
//  were file-local statics in src/World.cpp. They are hoisted here verbatim so
//  the engine TUs and the unit-test TU bind the *same* symbol, replacing the
//  hand-synced copy tests/EngineHelpers.hpp (now deleted, both sides #include
//  this).
//
//   - advanceQuatExp        : exact exp-map unit-quaternion advance (parent frame)
//   - rotationToQuaternion  : Shepperd R -> quaternion (w,x,y,z)
//   - quatToRotation        : unit quaternion -> R (row-major)
//   - safeLogSineSqr        : floored ln sin^2(pitch) for the torsional J(q)
// ============================================================================

#include <cmath>

#include "robot_math.hpp"

namespace EngineHelpers {

using robo::Mat33;
using robo::Real;
using robo::Rotation;
using robo::Vec3;

// -------- exact unit-quaternion advance (exponential map) ------------------
// Advance a unit quaternion under a constant angular velocity w_F expressed in
// the PARENT (F) frame, consistent with the engine's qdot = 1/2 (0,w_F) (x) q
// (left multiply; matches convertAngVelToQuaternionDot / N(q)). With theta =
// 1/2 |w| h: Dq = (cos theta, sin theta * w_hat) and q1 = Dq (x) q0. Both
// operands are unit, so q1 is unit BY CONSTRUCTION -- the |q|^2 -> Inf ->
// q/sqrt(Inf) = 0 overflow path of the linear "q += h*qdot" drift cannot occur.
// Reversible: negating w gives Dq^-1, so q0 = Dq^-1 (x) q1 (HMC needs this).
// As h -> 0 it reduces to q0 + h * (1/2 (0,w) (x) q0) = q0 + h * qdot0.
inline void advanceQuatExp(const Real* q0, const Vec3& wF, Real h, Real* q1) {
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

// Rotation (row-major Mat33) -> unit quaternion (w,x,y,z). Shepperd's method.
inline void rotationToQuaternion(const Rotation& R, Real& w, Real& x, Real& y, Real& z) {
    const Real m00 = R(0, 0), m11 = R(1, 1), m22 = R(2, 2);
    const Real tr = m00 + m11 + m22;
    if (tr > 0) {
        Real s = std::sqrt(tr + 1.0) * 2.0;
        w = 0.25 * s;
        x = (R(2, 1) - R(1, 2)) / s;
        y = (R(0, 2) - R(2, 0)) / s;
        z = (R(1, 0) - R(0, 1)) / s;
    } else if (m00 > m11 && m00 > m22) {
        Real s = std::sqrt(1.0 + m00 - m11 - m22) * 2.0;
        w = (R(2, 1) - R(1, 2)) / s;
        x = 0.25 * s;
        y = (R(0, 1) + R(1, 0)) / s;
        z = (R(0, 2) + R(2, 0)) / s;
    } else if (m11 > m22) {
        Real s = std::sqrt(1.0 + m11 - m00 - m22) * 2.0;
        w = (R(0, 2) - R(2, 0)) / s;
        x = (R(0, 1) + R(1, 0)) / s;
        y = 0.25 * s;
        z = (R(1, 2) + R(2, 1)) / s;
    } else {
        Real s = std::sqrt(1.0 + m22 - m00 - m11) * 2.0;
        w = (R(1, 0) - R(0, 1)) / s;
        x = (R(0, 2) + R(2, 0)) / s;
        y = (R(1, 2) + R(2, 1)) / s;
        z = 0.25 * s;
    }
}

// log(sin^2(pitch)) with a floor that keeps it finite near the singularity.
inline Real safeLogSineSqr(Real pitch) {
    const Real s = std::sin(pitch);
    Real s2 = s * s;
    constexpr Real kFloor = 1e-12;
    if (s2 < kFloor) {
        s2 = kFloor;
    }
    return std::log(s2);
}

// Unit quaternion (w,x,y,z) -> rotation matrix (row-major).
inline Rotation quatToRotation(Real qw, Real qx, Real qy, Real qz) {
    const Real xx = qx * qx, yy = qy * qy, zz = qz * qz;
    const Real xy = qx * qy, xz = qx * qz, yz = qy * qz;
    const Real wx = qw * qx, wy = qw * qy, wz = qw * qz;
    return Rotation(Mat33(1 - 2 * (yy + zz),
                          2 * (xy - wz),
                          2 * (xz + wy),
                          2 * (xy + wz),
                          1 - 2 * (xx + zz),
                          2 * (yz - wx),
                          2 * (xz - wy),
                          2 * (yz + wx),
                          1 - 2 * (xx + yy)));
}

} // namespace EngineHelpers