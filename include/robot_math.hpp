#pragma once
// ============================================================================
//  robot_math_v2.hpp -- in-house POD value types, SimTK-free.
//
//  TARGET: x86-64-v3 (AVX2 + FMA + BMI). Build with
//      -O3 -march=x86-64-v3 -ffp-contract=fast -funroll-loops
//
//  Scalar PODs for the per-body articulated recursion (sequential parent<->child,
//  not vectorisable across a branch; the win is flat topological-order layout +
//  FMA contraction, not intrinsics). The per-ATOM bulk SIMD work (position
//  broadcast, force gather) lives in the engine's hot loops over RobotState
//  arrays, NOT here -- value types stay branch-free PODs. (The earlier AtomSoA /
//  transformAtomsInPlace helpers were never referenced by the engine, so they
//  are removed; when the per-atom kernels are migrated they belong in
//  ForceBridge / RobotState next to the CSR body->atom map, gathering rather
//  than scattering.)
//
//  PRECISION: double throughout the dynamics. Only the OpenMM transfer narrows.
//  Highest numeric risk: the quaternion kinematics -- golden-test them.
//
//  clang-tidy: std::array (no C arrays), trailing return types, parameter names
//  >= 3 chars, explicit parentheses around '*' vs '+/-', no redundant `inline`
//  (constexpr / in-class members are implicitly inline).
// ============================================================================

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace robo {

using Real = double;

inline constexpr Real Pi = Real(3.14159265358979323846);
inline constexpr Real Deg2Rad = Pi / Real(180);

// coordinate axes (unscoped so `XAxis` etc. read like the SimTK names)
enum CoordinateAxis : int {
    XAxis = 0,
    YAxis = 1,
    ZAxis = 2
};

// ---------------------------------------------------------------------------
//  Vec3
// ---------------------------------------------------------------------------
struct alignas(8) Vec3 {
    std::array<Real, 3> elems;

    Vec3() = default;
    constexpr explicit Vec3(Real fill)
        : elems{fill, fill, fill} {
    }
    constexpr Vec3(Real valX, Real valY, Real valZ)
        : elems{valX, valY, valZ} {
    }

    constexpr auto operator[](int idx) -> Real& {
        return elems[static_cast<std::size_t>(idx)];
    }
    constexpr auto operator[](int idx) const -> Real {
        return elems[static_cast<std::size_t>(idx)];
    }

    constexpr auto operator+(const Vec3& rhs) const -> Vec3 {
        return {elems[0] + rhs.elems[0], elems[1] + rhs.elems[1], elems[2] + rhs.elems[2]};
    }
    constexpr auto operator-(const Vec3& rhs) const -> Vec3 {
        return {elems[0] - rhs.elems[0], elems[1] - rhs.elems[1], elems[2] - rhs.elems[2]};
    }
    constexpr auto operator*(Real scalar) const -> Vec3 {
        return {elems[0] * scalar, elems[1] * scalar, elems[2] * scalar};
    }
    constexpr auto operator/(Real scalar) const -> Vec3 {
        return {elems[0] / scalar, elems[1] / scalar, elems[2] / scalar};
    }
    constexpr auto operator+=(const Vec3& rhs) -> Vec3& {
        elems[0] += rhs.elems[0];
        elems[1] += rhs.elems[1];
        elems[2] += rhs.elems[2];
        return *this;
    }
    constexpr auto operator-=(const Vec3& rhs) -> Vec3& {
        elems[0] -= rhs.elems[0];
        elems[1] -= rhs.elems[1];
        elems[2] -= rhs.elems[2];
        return *this;
    }
    constexpr auto operator%(const Vec3& rhs) const -> Vec3 { // cross product
        return {(elems[1] * rhs.elems[2]) - (elems[2] * rhs.elems[1]),
                (elems[2] * rhs.elems[0]) - (elems[0] * rhs.elems[2]),
                (elems[0] * rhs.elems[1]) - (elems[1] * rhs.elems[0])};
    }
    [[nodiscard]] auto norm() const -> Real {
        return std::sqrt(normSqr());
    }
    [[nodiscard]] constexpr auto normSqr() const -> Real {
        return (elems[0] * elems[0]) + (elems[1] * elems[1]) + (elems[2] * elems[2]);
    }
};
constexpr auto dot(const Vec3& lhs, const Vec3& rhs) -> Real {
    return (lhs.elems[0] * rhs.elems[0]) + (lhs.elems[1] * rhs.elems[1]) + (lhs.elems[2] * rhs.elems[2]);
}
constexpr auto operator*(Real scalar, const Vec3& vec) -> Vec3 {
    return vec * scalar;
}

// ---------------------------------------------------------------------------
//  Mat33 (row-major). Non-aggregate: diagonal ctor + 9-element ctor.
// ---------------------------------------------------------------------------
struct Mat33 {
    std::array<Real, 9> elems; // row-major: elems[(row*3)+col]

    Mat33() = default;
    explicit constexpr Mat33(Real diag)
        : elems{diag, 0, 0, 0, diag, 0, 0, 0, diag} {
    }
    constexpr Mat33(Real m00, Real m01, Real m02, Real m10, Real m11, Real m12, Real m20, Real m21, Real m22)
        : elems{m00, m01, m02, m10, m11, m12, m20, m21, m22} {
    }

    static constexpr auto identity() -> Mat33 {
        return Mat33(Real(1));
    }
    static constexpr auto zero() -> Mat33 {
        return Mat33(Real(0));
    }

    constexpr auto operator()(int row, int col) -> Real& {
        return elems[static_cast<std::size_t>((row * 3) + col)];
    }
    constexpr auto operator()(int row, int col) const -> Real {
        return elems[static_cast<std::size_t>((row * 3) + col)];
    }

    constexpr auto operator*(const Vec3& vec) const -> Vec3 {
        return {(elems[0] * vec[0]) + (elems[1] * vec[1]) + (elems[2] * vec[2]),
                (elems[3] * vec[0]) + (elems[4] * vec[1]) + (elems[5] * vec[2]),
                (elems[6] * vec[0]) + (elems[7] * vec[1]) + (elems[8] * vec[2])};
    }
    [[nodiscard]] constexpr auto transposeTimes(const Vec3& vec) const -> Vec3 {
        return {(elems[0] * vec[0]) + (elems[3] * vec[1]) + (elems[6] * vec[2]),
                (elems[1] * vec[0]) + (elems[4] * vec[1]) + (elems[7] * vec[2]),
                (elems[2] * vec[0]) + (elems[5] * vec[1]) + (elems[8] * vec[2])};
    }
    constexpr auto operator+(const Mat33& rhs) const -> Mat33 {
        return Mat33(elems[0] + rhs.elems[0],
                     elems[1] + rhs.elems[1],
                     elems[2] + rhs.elems[2],
                     elems[3] + rhs.elems[3],
                     elems[4] + rhs.elems[4],
                     elems[5] + rhs.elems[5],
                     elems[6] + rhs.elems[6],
                     elems[7] + rhs.elems[7],
                     elems[8] + rhs.elems[8]);
    }
    constexpr auto operator-(const Mat33& rhs) const -> Mat33 {
        return Mat33(elems[0] - rhs.elems[0],
                     elems[1] - rhs.elems[1],
                     elems[2] - rhs.elems[2],
                     elems[3] - rhs.elems[3],
                     elems[4] - rhs.elems[4],
                     elems[5] - rhs.elems[5],
                     elems[6] - rhs.elems[6],
                     elems[7] - rhs.elems[7],
                     elems[8] - rhs.elems[8]);
    }
    constexpr auto operator*(const Mat33& rhs) const -> Mat33 {
        Mat33 out(Real(0));
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                out(row, col) = ((*this)(row, 0) * rhs(0, col)) + ((*this)(row, 1) * rhs(1, col))
                                + ((*this)(row, 2) * rhs(2, col));
            }
        }
        return out;
    }
    constexpr auto operator*(Real scalar) const -> Mat33 {
        return Mat33(elems[0] * scalar,
                     elems[1] * scalar,
                     elems[2] * scalar,
                     elems[3] * scalar,
                     elems[4] * scalar,
                     elems[5] * scalar,
                     elems[6] * scalar,
                     elems[7] * scalar,
                     elems[8] * scalar);
    }
    [[nodiscard]] constexpr auto transpose() const -> Mat33 {
        return Mat33(elems[0],
                     elems[3],
                     elems[6],
                     elems[1],
                     elems[4],
                     elems[7],
                     elems[2],
                     elems[5],
                     elems[8]);
    }
};

constexpr auto crossMat(const Vec3& vec) -> Mat33 { // crossMat(v)*x == v % x
    return Mat33(0, -vec[2], vec[1], vec[2], 0, -vec[0], -vec[1], vec[0], 0);
}

// ---------------------------------------------------------------------------
//  UnitVec3 (normalized on construction)
// ---------------------------------------------------------------------------
struct UnitVec3 {
    Vec3 dir{0, 0, 1};

    UnitVec3() = default;
    explicit UnitVec3(const Vec3& vec) {
        const Real len = vec.norm();
        dir = (len > Real(1e-300)) ? (vec / len) : Vec3(0, 0, 1);
    }
    [[nodiscard]] auto asVec3() const -> const Vec3& {
        return dir;
    }
    operator const Vec3&() const {
        return dir;
    } // NOLINT(google-explicit-constructor): a UnitVec3 IS-A Vec3
    constexpr auto operator[](int idx) const -> Real {
        return dir[idx];
    }
};

// ---------------------------------------------------------------------------
//  Vec4 (raw quaternion storage / qdot result)
// ---------------------------------------------------------------------------
struct alignas(8) Vec4 {
    std::array<Real, 4> elems;

    Vec4() = default;
    constexpr Vec4(Real valW, Real valX, Real valY, Real valZ)
        : elems{valW, valX, valY, valZ} {
    }
    constexpr auto operator[](int idx) -> Real& {
        return elems[static_cast<std::size_t>(idx)];
    }
    constexpr auto operator[](int idx) const -> Real {
        return elems[static_cast<std::size_t>(idx)];
    }
};

// ---------------------------------------------------------------------------
//  Quaternion (w, x, y, z). q represents R_FM; angVel is w_FM in F (== u[0..2]).
//  qdot = N(q) * w. GOLDEN-TEST against the SimTK reference.
// ---------------------------------------------------------------------------
struct Quat {
    std::array<Real, 4> elems; // w, x, y, z

    Quat() = default;
    constexpr Quat(Real valW, Real valX, Real valY, Real valZ)
        : elems{valW, valX, valY, valZ} {
    }
    explicit Quat(const Vec4& raw)
        : elems{raw[0], raw[1], raw[2], raw[3]} {
        normalize();
    }

    [[nodiscard]] auto norm() const -> Real {
        return std::sqrt((elems[0] * elems[0]) + (elems[1] * elems[1]) + (elems[2] * elems[2])
                         + (elems[3] * elems[3]));
    }
    auto normalize() -> void {
        const Real len = norm();
        if (len > Real(1e-12)) {
            const Real inv = Real(1) / len;
            for (Real& comp : elems) {
                comp *= inv;
            }
        } else {
            elems = {1, 0, 0, 0};
        }
    }
    static auto angVelToQdot(const Quat& quat, const Vec3& angVel) -> Quat {
        const Real qw = quat.elems[0];
        const Real qx = quat.elems[1];
        const Real qy = quat.elems[2];
        const Real qz = quat.elems[3];
        return Quat(Real(0.5) * ((-qx * angVel[0]) - (qy * angVel[1]) - (qz * angVel[2])),
                    Real(0.5) * ((qw * angVel[0]) - (qz * angVel[1]) + (qy * angVel[2])),
                    Real(0.5) * ((qz * angVel[0]) + (qw * angVel[1]) - (qx * angVel[2])),
                    Real(0.5) * ((-qy * angVel[0]) + (qx * angVel[1]) + (qw * angVel[2])));
    }
};
using Quaternion = Quat; // SimTK-compatible name

// ---------------------------------------------------------------------------
//  SymMat33. Storage order matches the engine's 6-arg ctor: (xx, xy, yy, xz, yz, zz).
//      full = [[xx xy xz],[xy yy yz],[xz yz zz]]
// ---------------------------------------------------------------------------
struct Rotation; // fwd: needed by SymMat33::reexpress

struct SymMat33 {
    std::array<Real, 6> elems; // xx, xy, yy, xz, yz, zz

    SymMat33() = default;
    // SimTK-compatible 6-arg order (xx, xy, yy, xz, yz, zz)
    constexpr SymMat33(Real exx, Real exy, Real eyy, Real exz, Real eyz, Real ezz)
        : elems{exx, exy, eyy, exz, eyz, ezz} {
    }
    // diagonal moments (UnitInertia(xx, yy, zz))
    constexpr SymMat33(Real exx, Real eyy, Real ezz)
        : elems{exx, 0, eyy, 0, 0, ezz} {
    }

    static constexpr auto zero() -> SymMat33 {
        return SymMat33(Real(0), Real(0), Real(0), Real(0), Real(0), Real(0));
    }
    static constexpr auto diagonal(Real value) -> SymMat33 {
        return SymMat33(value, value, value);
    }

    constexpr auto operator*(const Vec3& vec) const -> Vec3 {
        return {(elems[0] * vec[0]) + (elems[1] * vec[1]) + (elems[3] * vec[2]),
                (elems[1] * vec[0]) + (elems[2] * vec[1]) + (elems[4] * vec[2]),
                (elems[3] * vec[0]) + (elems[4] * vec[1]) + (elems[5] * vec[2])};
    }
    constexpr auto operator*(Real scalar) const -> SymMat33 {
        return SymMat33(elems[0] * scalar,
                        elems[1] * scalar,
                        elems[2] * scalar,
                        elems[3] * scalar,
                        elems[4] * scalar,
                        elems[5] * scalar);
    }
    constexpr auto operator+(const SymMat33& rhs) const -> SymMat33 {
        return SymMat33(elems[0] + rhs.elems[0],
                        elems[1] + rhs.elems[1],
                        elems[2] + rhs.elems[2],
                        elems[3] + rhs.elems[3],
                        elems[4] + rhs.elems[4],
                        elems[5] + rhs.elems[5]);
    }
    constexpr auto operator-(const SymMat33& rhs) const -> SymMat33 {
        return SymMat33(elems[0] - rhs.elems[0],
                        elems[1] - rhs.elems[1],
                        elems[2] - rhs.elems[2],
                        elems[3] - rhs.elems[3],
                        elems[4] - rhs.elems[4],
                        elems[5] - rhs.elems[5]);
    }
    constexpr auto operator+=(const SymMat33& rhs) -> SymMat33& {
        for (std::size_t idx = 0; idx < 6; ++idx) {
            elems[idx] += rhs.elems[idx];
        }
        return *this;
    }
    [[nodiscard]] constexpr auto full() const -> Mat33 {
        return Mat33(elems[0],
                     elems[1],
                     elems[3],
                     elems[1],
                     elems[2],
                     elems[4],
                     elems[3],
                     elems[4],
                     elems[5]);
    }
    static constexpr auto fromSymmetricPart(const Mat33& mat) -> SymMat33 {
        return SymMat33(mat(0, 0),
                        (mat(1, 0) + mat(0, 1)) / 2,
                        mat(1, 1),
                        (mat(2, 0) + mat(0, 2)) / 2,
                        (mat(2, 1) + mat(1, 2)) / 2,
                        mat(2, 2));
    }
    // re-express in a rotated frame: ~rot * (*this) * rot (SimTK convention). Defined after Rotation.
    [[nodiscard]] auto reexpress(const Rotation& rot) const -> SymMat33;
};
using UnitInertia = SymMat33; // inertia per unit mass about a point

// ---------------------------------------------------------------------------
//  Rotation: orthonormal Mat33 with assorted constructors / setters.
// ---------------------------------------------------------------------------
struct Rotation : Mat33 {
    Rotation()
        : Mat33(Real(1)) {
    }
    Rotation(const Mat33& mat)
        : Mat33(mat) {
    } // NOLINT(google-explicit-constructor): Rotation IS-A Mat33
    Rotation(Real angle, CoordinateAxis axis) {
        setRotationFromAngleAboutAxis(angle, axis);
    }
    Rotation(const UnitVec3& primaryDir,
             CoordinateAxis primaryAxis,
             const Vec3& planeVec,
             CoordinateAxis planeAxis) {
        setRotationFromTwoAxes(primaryDir, primaryAxis, planeVec, planeAxis);
    }

    auto setColumn(int col, const Vec3& vec) -> void {
        elems[static_cast<std::size_t>(col)] = vec[0];
        elems[static_cast<std::size_t>(3 + col)] = vec[1];
        elems[static_cast<std::size_t>(6 + col)] = vec[2];
    }

    auto setRotationFromAngleAboutAxis(Real angle, CoordinateAxis axis) -> void {
        const Real cosA = std::cos(angle);
        const Real sinA = std::sin(angle);
        switch (axis) {
            case XAxis:
                elems = {1, 0, 0, 0, cosA, -sinA, 0, sinA, cosA};
                break;
            case YAxis:
                elems = {cosA, 0, sinA, 0, 1, 0, -sinA, 0, cosA};
                break;
            case ZAxis:
            default:
                elems = {cosA, -sinA, 0, sinA, cosA, 0, 0, 0, 1};
                break;
        }
    }
    auto setRotationFromAngleAboutZ(Real angle) -> void {
        setRotationFromAngleAboutAxis(angle, ZAxis);
    }
    auto setRotationFromQuaternion(const Quat& quat) -> void {
        *this = fromQuaternion(quat);
    }

    // place `dir` on coordinate axis `axis`; complete with an arbitrary orthonormal pair.
    auto setRotationFromOneAxis(const UnitVec3& dir, CoordinateAxis axis) -> void {
        const int idx1 = static_cast<int>(axis);
        const int idx2 = (idx1 + 1) % 3;
        const int idx3 = (idx1 + 2) % 3;
        const Vec3 col1 = dir.asVec3();
        const Vec3 col2 = anyUnitPerpendicular(col1);
        const Vec3 col3 = col1 % col2; // (idx1,idx2,idx3) cyclic -> col1 x col2 = col3
        setColumn(idx1, col1);
        setColumn(idx2, col2);
        setColumn(idx3, col3);
    }

    // primaryAxis aligned to primaryDir; planeAxis as close to planeVec as possible.
    auto setRotationFromTwoAxes(const UnitVec3& primaryDir,
                                CoordinateAxis primaryAxis,
                                const Vec3& planeVec,
                                CoordinateAxis planeAxis) -> void {
        const int idxU = static_cast<int>(primaryAxis);
        const int idxV = static_cast<int>(planeAxis);
        const int idxW = 3 - idxU - idxV;
        const Vec3 colU = primaryDir.asVec3();
        const Vec3 perp = planeVec - (colU * dot(planeVec, colU)); // Gram-Schmidt toward planeVec
        const Real perpNorm = perp.norm();
        const Vec3 colV = (perpNorm > Real(1e-12)) ? (perp / perpNorm) : anyUnitPerpendicular(colU);
        const bool cyclic = (idxV == ((idxU + 1) % 3));
        const Vec3 colW = cyclic ? (colU % colV) : (colV % colU);
        setColumn(idxU, colU);
        setColumn(idxV, colV);
        setColumn(idxW, colW);
    }

    [[nodiscard]] auto operator~() const -> Rotation {
        return Rotation(transpose());
    }

    static auto fromQuaternion(const Quat& source) -> Rotation {
        Quat normalized = source;
        normalized.normalize();
        const Real qw = normalized.elems[0];
        const Real qx = normalized.elems[1];
        const Real qy = normalized.elems[2];
        const Real qz = normalized.elems[3];
        Rotation rot;
        rot.elems = {1 - (2 * ((qy * qy) + (qz * qz))),
                     2 * ((qx * qy) - (qw * qz)),
                     2 * ((qx * qz) + (qw * qy)),
                     2 * ((qx * qy) + (qw * qz)),
                     1 - (2 * ((qx * qx) + (qz * qz))),
                     2 * ((qy * qz) - (qw * qx)),
                     2 * ((qx * qz) - (qw * qy)),
                     2 * ((qy * qz) + (qw * qx)),
                     1 - (2 * ((qx * qx) + (qy * qy)))};
        return rot;
    }

    // qdot = N(quat) * angVel (SimTK Rotation::convertAngVelToQuaternionDot). quat treated as Vec4.
    static auto convertAngVelToQuaternionDot(const Vec4& quat, const Vec3& angVel) -> Vec4 {
        const Real ew = quat[0];
        const Real ex = quat[1];
        const Real ey = quat[2];
        const Real ez = quat[3];
        return Vec4(Real(0.5) * ((-ex * angVel[0]) - (ey * angVel[1]) - (ez * angVel[2])),
                    Real(0.5) * ((ew * angVel[0]) - (ez * angVel[1]) + (ey * angVel[2])),
                    Real(0.5) * ((ez * angVel[0]) + (ew * angVel[1]) - (ex * angVel[2])),
                    Real(0.5) * ((-ey * angVel[0]) + (ex * angVel[1]) + (ew * angVel[2])));
    }
    static auto convertAngVelToQuaternionDot(const Quat& quat, const Vec3& angVel) -> Vec4 {
        return convertAngVelToQuaternionDot(Vec4(quat.elems[0], quat.elems[1], quat.elems[2], quat.elems[3]),
                                            angVel);
    }
    // qddot = N*angVelDot + Ndot*angVel ; Ndot built from qdot = N*angVel (N linear in q).
    static auto
    convertAngVelDotToQuaternionDotDot(const Vec4& quat, const Vec3& angVel, const Vec3& angVelDot) -> Vec4 {
        const Vec4 termA = convertAngVelToQuaternionDot(quat, angVelDot);
        const Vec4 qdot = convertAngVelToQuaternionDot(quat, angVel);
        const Vec4 termB = convertAngVelToQuaternionDot(qdot, angVel);
        return Vec4(termA[0] + termB[0], termA[1] + termB[1], termA[2] + termB[2], termA[3] + termB[3]);
    }

    private:
    static auto anyUnitPerpendicular(const Vec3& vec) -> Vec3 {
        const Real absX = std::abs(vec[0]);
        const Real absY = std::abs(vec[1]);
        const Real absZ = std::abs(vec[2]);
        const Vec3 reference = (absX <= absY && absX <= absZ) ? Vec3(1, 0, 0)
                               : (absY <= absZ)               ? Vec3(0, 1, 0)
                                                              : Vec3(0, 0, 1);
        const Vec3 perp = vec % reference;
        const Real len = perp.norm();
        return (len > Real(1e-12)) ? (perp / len) : Vec3(1, 0, 0);
    }
};

inline auto SymMat33::reexpress(const Rotation& rot) const -> SymMat33 {
    const Mat33 rotated = (rot.transpose() * full()) * rot; // ~R * I * R
    return SymMat33::fromSymmetricPart(rotated);
}

// free dihedral angle (radians), SimTK::calcDihedralAngle(a, b, c, d)
inline auto calcDihedralAngle(const Vec3& atomA, const Vec3& atomB, const Vec3& atomC, const Vec3& atomD)
    -> Real {
    const Vec3 edge1 = atomB - atomA;
    const Vec3 edge2 = atomC - atomB;
    const Vec3 edge3 = atomD - atomC;
    const Vec3 normal1 = edge1 % edge2;
    const Vec3 normal2 = edge2 % edge3;
    const Vec3 frame = normal1 % (edge2 / (edge2.norm() + Real(1e-300)));
    return std::atan2(dot(frame, normal2), dot(normal1, normal2));
}

// ---------------------------------------------------------------------------
//  Transform: rotation + translation, X_AB (frame B expressed in A).
// ---------------------------------------------------------------------------
struct Transform {
    Rotation rot;
    Vec3 trans{0, 0, 0};

    Transform()
        : rot(Mat33::identity()) {
    }
    Transform(const Rotation& rotation, const Vec3& translation)
        : rot(rotation)
        , trans(translation) {
    }
    explicit Transform(const Rotation& rotation)
        : rot(rotation) {
    }
    explicit Transform(const Vec3& translation)
        : rot(Mat33::identity())
        , trans(translation) {
    }

    [[nodiscard]] auto p() const -> const Vec3& {
        return trans;
    }
    [[nodiscard]] auto p() -> Vec3& {
        return trans;
    }
    [[nodiscard]] auto R() const -> const Rotation& {
        return rot;
    }
    [[nodiscard]] auto R() -> Rotation& {
        return rot;
    }

    auto operator*(const Transform& rhs) const -> Transform {
        return Transform{Rotation(rot * rhs.rot), (rot * rhs.trans) + trans};
    }
    auto operator*(const Vec3& point) const -> Vec3 {
        return (rot * point) + trans;
    }
    [[nodiscard]] auto inverse() const -> Transform {
        const Mat33 rotT = rot.transpose();
        return Transform{Rotation(rotT), Vec3(0) - (rotT * trans)};
    }
    [[nodiscard]] auto operator~() const -> Transform {
        return inverse();
    }
};

// ---------------------------------------------------------------------------
//  SpatialVec : 6-vector [angular; linear]. Index [0]=angular, [1]=linear.
// ---------------------------------------------------------------------------
struct SpatialVec {
    Vec3 angular;
    Vec3 linear;

    SpatialVec() = default;
    SpatialVec(const Vec3& angularPart, const Vec3& linearPart)
        : angular(angularPart)
        , linear(linearPart) {
    }

    auto operator[](int idx) -> Vec3& {
        return idx == 0 ? angular : linear;
    }
    auto operator[](int idx) const -> const Vec3& {
        return idx == 0 ? angular : linear;
    }

    auto operator+(const SpatialVec& rhs) const -> SpatialVec {
        return {angular + rhs.angular, linear + rhs.linear};
    }
    auto operator-(const SpatialVec& rhs) const -> SpatialVec {
        return {angular - rhs.angular, linear - rhs.linear};
    }
    auto operator*(Real scalar) const -> SpatialVec {
        return {angular * scalar, linear * scalar};
    }
    auto operator+=(const SpatialVec& rhs) -> SpatialVec& {
        angular += rhs.angular;
        linear += rhs.linear;
        return *this;
    }
};
inline auto dot(const SpatialVec& lhs, const SpatialVec& rhs) -> Real {
    return dot(lhs.angular, rhs.angular) + dot(lhs.linear, rhs.linear);
}
inline auto operator*(Real scalar, const SpatialVec& vec) -> SpatialVec {
    return vec * scalar;
}

// ---------------------------------------------------------------------------
//  Inertia: full 2nd mass moment about a point. SimTK Inertia(p, m) / Inertia(0).
// ---------------------------------------------------------------------------
struct Inertia {
    SymMat33 moments{SymMat33::zero()};

    Inertia() = default;
    explicit Inertia(Real uniformMoment)
        : moments(SymMat33::diagonal(uniformMoment)) {
    }
    Inertia(const Vec3& point, Real pointMass)
        : moments(pointMassInertia(point, pointMass)) {
    }

    auto operator+=(const Inertia& rhs) -> Inertia& {
        moments += rhs.moments;
        return *this;
    }
    auto operator+(const Inertia& rhs) const -> Inertia {
        Inertia out;
        out.moments = moments + rhs.moments;
        return out;
    }
    [[nodiscard]] auto asSymMat33() const -> const SymMat33& {
        return moments;
    }

    private:
    static auto pointMassInertia(const Vec3& point, Real pointMass) -> SymMat33 {
        const Real posX = point[0];
        const Real posY = point[1];
        const Real posZ = point[2];
        // I = m * (|p|^2 I3 - p p^T); arg order (xx, xy, yy, xz, yz, zz)
        return SymMat33(pointMass * ((posY * posY) + (posZ * posZ)),
                        pointMass * (-(posX * posY)),
                        pointMass * ((posX * posX) + (posZ * posZ)),
                        pointMass * (-(posX * posZ)),
                        pointMass * (-(posY * posZ)),
                        pointMass * ((posX * posX) + (posY * posY)));
    }
};

// ---------------------------------------------------------------------------
//  SpatialInertia: rigid-body spatial inertia about the body origin.
//  Mk_G = SpatialInertia(mass, p_BBc_G, G_Bo_G).  Block form (about origin):
//      [[mass*G, mass*crossMat(com)], [~., mass*I]].
// ---------------------------------------------------------------------------
struct SpatialInertia {
    Real mass{0};
    Vec3 com{0, 0, 0};                         // p_BBc (origin -> COM)
    UnitInertia unitInertia{SymMat33::zero()}; // G about origin

    SpatialInertia() = default;
    SpatialInertia(Real bodyMass, const Vec3& massCenter, const UnitInertia& gyration)
        : mass(bodyMass)
        , com(massCenter)
        , unitInertia(gyration) {
    }

    [[nodiscard]] auto getMass() const -> Real {
        return mass;
    }
    [[nodiscard]] auto getMassCenter() const -> const Vec3& {
        return com;
    }
    [[nodiscard]] auto getUnitInertia() const -> const UnitInertia& {
        return unitInertia;
    }

    // spatial momentum = M_spatial * spatial velocity. Defined after ArticulatedInertia.
    [[nodiscard]] auto operator*(const SpatialVec& vel) const -> SpatialVec;
};

// ---------------------------------------------------------------------------
//  ArticulatedInertia: symmetric 6x6 in 3 blocks
//      P = [[angAng, angLin],[~angLin, linLin]]
//      out.angular = angAng*w + angLin*v ;  out.linear = ~angLin*w + linLin*v
//  Engine ctor order: (mass=linLin, massMoment=angLin, inertia=angAng).
// ---------------------------------------------------------------------------
struct ArticulatedInertia {
    SymMat33 angAng{SymMat33::zero()}; // J
    Mat33 angLin{Mat33::zero()};       // F
    SymMat33 linLin{SymMat33::zero()}; // M

    ArticulatedInertia() = default;
    ArticulatedInertia(const SymMat33& massBlock, const Mat33& momentBlock, const SymMat33& inertiaBlock)
        : angAng(inertiaBlock)
        , angLin(momentBlock)
        , linLin(massBlock) {
    }
    explicit ArticulatedInertia(const SpatialInertia& spatial)
        : angAng(spatial.unitInertia * spatial.mass)
        , angLin(crossMat(spatial.com) * spatial.mass)
        , linLin(SymMat33::diagonal(spatial.mass)) {
    }

    auto operator+=(const ArticulatedInertia& rhs) -> ArticulatedInertia& {
        angAng += rhs.angAng;
        angLin = angLin + rhs.angLin;
        linLin += rhs.linLin;
        return *this;
    }
    auto operator-(const ArticulatedInertia& rhs) const -> ArticulatedInertia {
        ArticulatedInertia out;
        out.angAng = angAng - rhs.angAng;
        out.angLin = angLin - rhs.angLin;
        out.linLin = linLin - rhs.linLin;
        return out;
    }
    auto operator*(const SpatialVec& vec) const -> SpatialVec {
        return {(angAng * vec.angular) + (angLin * vec.linear),
                angLin.transposeTimes(vec.angular) + (linLin * vec.linear)};
    }

    // shift to a point `offset` away: Phi(offset) * P * ~Phi(offset).
    [[nodiscard]] auto shift(const Vec3& offset) const -> ArticulatedInertia {
        const Mat33 skew = crossMat(offset);
        const Mat33 massFull = linLin.full();
        const Mat33 momentFull = angLin;
        const Mat33 newMoment = momentFull + (skew * massFull);
        const Mat33 term1 = skew * momentFull.transpose(); // sx * ~F
        const Mat33 term2 = momentFull * skew;             // F * sx
        const Mat33 term3 = (skew * massFull) * skew;      // sx * M * sx
        const Mat33 newInertiaFull = (angAng.full() + term1) - (term2 + term3);
        ArticulatedInertia out;
        out.linLin = linLin;
        out.angLin = newMoment;
        out.angAng = SymMat33::fromSymmetricPart(newInertiaFull);
        return out;
    }
};

inline auto SpatialInertia::operator*(const SpatialVec& vel) const -> SpatialVec {
    return ArticulatedInertia(*this) * vel; // identical block math, single source of truth
}

// ---------------------------------------------------------------------------
//  PhiMatrix: rigid shift with offset l() (parent origin -> child).
//      Phi(l)  * [a; b] = [a + l x b; b]
//      ~Phi(l) * [a; b] = [a; b + a x l]
// ---------------------------------------------------------------------------
struct PhiMatrixTranspose;

struct PhiMatrix {
    Vec3 offset{0, 0, 0}; // l()

    PhiMatrix() = default;
    explicit PhiMatrix(const Vec3& location)
        : offset(location) {
    }

    [[nodiscard]] auto l() const -> const Vec3& {
        return offset;
    }
    auto operator*(const SpatialVec& vec) const -> SpatialVec {
        return {vec.angular + (offset % vec.linear), vec.linear};
    }
    [[nodiscard]] auto operator~() const -> PhiMatrixTranspose;
};

struct PhiMatrixTranspose {
    Vec3 offset{0, 0, 0};
    auto operator*(const SpatialVec& vec) const -> SpatialVec {
        return {vec.angular, vec.linear + (vec.angular % offset)};
    }
};

inline auto PhiMatrix::operator~() const -> PhiMatrixTranspose {
    return PhiMatrixTranspose{offset};
}

// ---------------------------------------------------------------------------
//  MassProperties: mass + COM + unit inertia.
// ---------------------------------------------------------------------------
struct MassProperties {
    Real mass{0};
    Vec3 com{0, 0, 0};
    UnitInertia unitInertia{SymMat33::zero()};

    MassProperties() = default;
    MassProperties(Real bodyMass, const Vec3& massCenter, const UnitInertia& gyration)
        : mass(bodyMass)
        , com(massCenter)
        , unitInertia(gyration) {
    }
    // from a FULL inertia: store unit inertia = inertia / mass.
    MassProperties(Real bodyMass, const Vec3& massCenter, const Inertia& inertia)
        : mass(bodyMass)
        , com(massCenter)
        , unitInertia(bodyMass > Real(0) ? inertia.asSymMat33() * (Real(1) / bodyMass) : SymMat33::zero()) {
    }

    [[nodiscard]] auto getMass() const -> Real {
        return mass;
    }
    [[nodiscard]] auto getMassCenter() const -> const Vec3& {
        return com;
    }
    [[nodiscard]] auto getUnitInertia() const -> const UnitInertia& {
        return unitInertia;
    }

    // SimTK MassProperties::reexpress(R_BC): com' = ~R_BC*com, inertia' = unitInertia.reexpress(R_BC).
    [[nodiscard]] auto reexpress(const Rotation& rotation) const -> MassProperties {
        return MassProperties{mass, rotation.transpose() * com, unitInertia.reexpress(rotation)};
    }
    [[nodiscard]] auto toSpatialInertia() const -> SpatialInertia {
        return SpatialInertia{mass, com, unitInertia};
    }
};

// slab-safety (the MemoryArena frees without running destructors)
static_assert(sizeof(Vec3) == 3 * sizeof(Real), "Vec3 layout");
static_assert(sizeof(SpatialVec) == 6 * sizeof(Real), "SpatialVec layout");
static_assert(std::is_trivially_destructible<Vec3>::value, "Vec3 trivially destructible");
static_assert(std::is_trivially_destructible<SpatialVec>::value, "SpatialVec trivially destructible");
static_assert(std::is_trivially_destructible<Transform>::value, "Transform trivially destructible");
static_assert(std::is_trivially_destructible<ArticulatedInertia>::value,
              "ArticulatedInertia trivially destructible");
static_assert(std::is_trivially_destructible<MassProperties>::value, "MassProperties trivially destructible");

} // namespace robo
