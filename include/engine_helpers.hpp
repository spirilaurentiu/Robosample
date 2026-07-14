#pragma once

// ============================================================================
//  engine_helpers.hpp -- small orientation / exp-map helpers used by the
//  engine.
//
//  SINGLE SOURCE OF TRUTH. rotationToQuaternion / quatToRotation /
//  safeLogSineSqr were file-local statics in the anonymous namespace of
//  src/World.cpp; they are hoisted here verbatim (SPLIT-W5) so every engine
//  TU that needs them (FixmanCorrection.cpp, DockingMove.cpp, NcmcMove.cpp,
//  the still-present World.cpp) and the unit-test TU (tests/EngineHelpers.hpp)
//  bind the SAME symbol.
//
//   - rotationToQuaternion : Shepperd R -> quaternion (w,x,y,z)
//   - quatToRotation       : unit quaternion -> R (row-major)
//   - safeLogSineSqr       : floored ln sin^2(pitch) for the torsional J(q)
// ============================================================================

#include <cmath>

#include "robot_math.hpp"

namespace EngineHelpers {

using robo::Real;
using robo::Rotation;

/**
 * @brief Converts a rotation (row-major Mat33) to a unit quaternion (w,x,y,z)
 *        by Shepperd's method.
 *
 * @param[in]  R  Rotation matrix, row-major.
 * @param[out] w  Scalar quaternion component.
 * @param[out] x  Quaternion i component.
 * @param[out] y  Quaternion j component.
 * @param[out] z  Quaternion k component.
 */
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

/**
 * @brief Returns log(sin^2(@p pitch)) with a floor that keeps it finite near
 *        the singularity.
 *
 * @param[in] pitch  Pitch angle in radians.
 * @return log(sin^2(pitch)), with sin^2 floored at 1e-12 so the result stays
 *         finite as @p pitch approaches a multiple of pi.
 */
inline Real safeLogSineSqr(Real pitch) {
    const Real s = std::sin(pitch);
    Real s2 = s * s;
    constexpr Real kFloor = 1e-12;
    if (s2 < kFloor) {
        s2 = kFloor;
    }
    return std::log(s2);
}

/**
 * @brief Converts a unit quaternion (w,x,y,z) to a row-major rotation matrix.
 *
 * @param[in] qw  Scalar quaternion component.
 * @param[in] qx  Quaternion i component.
 * @param[in] qy  Quaternion j component.
 * @param[in] qz  Quaternion k component.
 * @return The corresponding rotation matrix (row-major).
 */
inline auto quatToRotation(Real qw, Real qx, Real qy, Real qz) -> Rotation {
    const Real xx = qx * qx;
    const Real yy = qy * qy;
    const Real zz = qz * qz;
    const Real xy = qx * qy;
    const Real xz = qx * qz;
    const Real yz = qy * qz;
    const Real wx = qw * qx;
    const Real wy = qw * qy;
    const Real wz = qw * qz;
    return {robo::Mat33(1 - (2 * (yy + zz)),
                        2 * (xy - wz),
                        2 * (xz + wy),
                        2 * (xy + wz),
                        1 - (2 * (xx + zz)),
                        2 * (yz - wx),
                        2 * (xz - wy),
                        2 * (yz + wx),
                        1 - (2 * (xx + yy)))};
}

} // namespace EngineHelpers
