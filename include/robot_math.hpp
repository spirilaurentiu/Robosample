#pragma once
/**
 * @file robot_math.hpp
 * @brief SimTK-free POD value types of spatial algebra: vectors, rotations,
 *        transforms, quaternions, and the spatial/articulated inertia operators
 *        the articulated-body recursion runs on.
 *
 * Every type here is a trivially-destructible value: copied by value, holds no
 * heap, and introduces no aliasing or ownership concern. The elementary
 * arithmetic operators (@c +, @c -, scalar @c *, @c [], @c +=) implement the
 * obvious componentwise algebra of the type that declares them and inherit that
 * type's contract; this file documents each type's meaning, storage convention,
 * and the semantically-loaded operations (quaternion kinematic map, frame
 * re-expression, spatial-inertia products) where a caller relies on more than
 * componentwise arithmetic.
 *
 * @note Precision: @ref robo::Real is @c double throughout the dynamics; only
 *       the transfer to OpenMM narrows. The quaternion kinematics carry the
 *       highest numeric risk and are pinned by golden tests.
 * @note Layout is fixed by @c static_assert (see end of file): the arena frees
 *       slabs without running destructors, so these types SHALL stay trivially
 *       destructible with the asserted sizes.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace robo {

/// The scalar type of the whole dynamics; @c double.
using Real = double;

inline constexpr Real Pi = Real(3.14159265358979323846);
inline constexpr Real Deg2Rad = Pi / Real(180);

/// Coordinate axis selector; unscoped so @c XAxis etc. read like the SimTK names.
enum CoordinateAxis : int {
    XAxis = 0,
    YAxis = 1,
    ZAxis = 2
};

/**
 * @brief A 3-vector over @ref Real.
 * @note @c operator% is the cross product (not modulo): @c a % b. The free
 *       function @ref dot(const Vec3&, const Vec3&) is the inner product;
 *       @ref normSqr avoids the square root when only a comparison is needed.
 */
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

/**
 * @brief A dense 3x3 matrix, row-major (@c elems[row*3+col]).
 * @note @c operator* multiplies @c M*v; @ref transposeTimes computes @c M^T*v
 *       without forming the transpose. @c operator()(row,col) is the element
 *       accessor.
 */
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

/**
 * @brief Skew-symmetric cross-product matrix of @p vec.
 * @param[in] vec Source vector.
 * @return The matrix @c S with @c S*x == vec % x for every @c x.
 */
constexpr auto crossMat(const Vec3& vec) -> Mat33 {
    return Mat33(0, -vec[2], vec[1], vec[2], 0, -vec[0], -vec[1], vec[0], 0);
}

/**
 * @brief A unit-length 3-vector, normalized at construction.
 * @note Constructing from a vector shorter than 1e-300 yields the fallback
 *       direction @c (0,0,1) rather than a NaN. Converts implicitly to
 *       @ref Vec3 (a @c UnitVec3 is-a direction), so it may be passed anywhere a
 *       @c Vec3 is expected.
 */
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

/**
 * @brief A 4-vector, used as raw quaternion storage @c (w,x,y,z) and as a
 *        @c qdot / @c qddot result before it is wrapped in a normalized @ref Quat.
 * @note Unlike @ref Quat, @c Vec4 is not normalized and carries no orientation
 *       contract; it is plain storage.
 */
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

/**
 * @brief Parent-frame quaternion kinematic map @f$ \dot q = N(q)\,\omega_{FM} @f$.
 * @param[in] qw Scalar component of @c q (the orientation @c R_FM of the body
 *               frame M in its parent F, as built by @ref Rotation::fromQuaternion).
 * @param[in] qx Vector component of @c q, x.
 * @param[in] qy Vector component of @c q, y.
 * @param[in] qz Vector component of @c q, z.
 * @param[in] w  Generalized angular velocity @c w_FM, expressed in the parent
 *               frame F.
 * @return @c qdot as a @ref Vec4, the left Hamilton product
 *         @f$ \tfrac12 (0,\omega)\otimes q @f$.
 * @warning The @p w argument SHALL be expressed in the parent frame F to match
 *          the standard @c R_FM convention of @ref Rotation::fromQuaternion.
 *          Feeding a body-frame angular velocity (the right-product map) advances
 *          the orientation with a wrong-handed velocity and pumps kinetic energy
 *          (INV-9).
 * @note Single source of truth: @ref Quat::angVelToQdot and
 *       @ref Rotation::convertAngVelToQuaternionDot both delegate here, so the
 *       two quaternion-derivative paths cannot diverge.
 */
inline auto quaternionDotFromAngVel(Real qw, Real qx, Real qy, Real qz, const Vec3& w) -> Vec4 {
    return Vec4(Real(0.5) * ((-qx * w[0]) - (qy * w[1]) - (qz * w[2])),
                Real(0.5) * ((qw * w[0]) + (qz * w[1]) - (qy * w[2])),
                Real(0.5) * ((-qz * w[0]) + (qw * w[1]) + (qx * w[2])),
                Real(0.5) * ((qy * w[0]) - (qx * w[1]) + (qw * w[2])));
}

/**
 * @brief A unit quaternion @c (w,x,y,z) representing the body orientation
 *        @c R_FM (frame M in its parent F).
 * @note Constructing from a @ref Vec4 normalizes; @ref normalize renormalizes in
 *       place and falls back to the identity @c (1,0,0,0) when the norm is below
 *       1e-12. @ref angVelToQdot takes a parent-frame angular velocity @c w_FM
 *       (the leading generalized speeds of a free/quaternion joint) and returns
 *       @c qdot via @ref quaternionDotFromAngVel. Represents an orientation up to
 *       the double cover @c q ~ -q (INV-9).
 */
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
        const Vec4 d =
            quaternionDotFromAngVel(quat.elems[0], quat.elems[1], quat.elems[2], quat.elems[3], angVel);
        return Quat(d[0], d[1], d[2], d[3]);
    }
};
using Quaternion = Quat; // SimTK-compatible name

/**
 * @brief A symmetric 3x3 matrix stored as 6 elements in order
 *        @c (xx, xy, yy, xz, yz, zz).
 *
 * The full matrix is @c [[xx xy xz],[xy yy yz],[xz yz zz]]. Aliased as
 * @ref UnitInertia (inertia per unit mass about a point).
 * @note @ref fromSymmetricPart symmetrizes an arbitrary @ref Mat33;
 *       @ref reexpress rotates the tensor into a new frame as @c ~R * (*this) * R
 *       (SimTK convention).
 */
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

/**
 * @brief An orthonormal @ref Mat33 representing a rigid-body rotation.
 *
 * Columns are the rotated axes; @c operator~ returns the transpose, which for an
 * orthonormal matrix is the inverse rotation. Constructors build a rotation from
 * an angle about a coordinate axis, from a quaternion (@ref fromQuaternion,
 * which normalizes its argument), or from one/two direction axes.
 * @note @ref fromQuaternion realizes the standard @c R_FM map paired with
 *       @ref quaternionDotFromAngVel (INV-9). @ref setRotationFromTwoAxes aligns
 *       @c primaryAxis exactly and places @c planeAxis as close to @c planeVec as
 *       orthonormality allows (Gram-Schmidt).
 */
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

    // qdot = N(quat) * angVel, parent-frame map (SimTK calcUnnormalizedNForQuaternion).
    // Delegates to the shared quaternionDotFromAngVel (single source of truth).
    static auto convertAngVelToQuaternionDot(const Vec4& quat, const Vec3& angVel) -> Vec4 {
        return quaternionDotFromAngVel(quat[0], quat[1], quat[2], quat[3], angVel);
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

/**
 * @brief Dihedral (torsion) angle of the four points @p atomA -> @p atomD.
 * @param[in] atomA First point.
 * @param[in] atomB Second point (first bond axis end).
 * @param[in] atomC Third point (second bond axis end).
 * @param[in] atomD Fourth point.
 * @return The signed dihedral angle in radians, in @c (-Pi, Pi], measured about
 *         the @c B->C axis with the SimTK sign convention.
 */
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

/**
 * @brief A rigid-body transform @c X_AB: frame B expressed in frame A
 *        (rotation @ref R plus translation @ref p).
 *
 * @c operator*(Transform) composes frames (@c X_AB * X_BC == X_AC);
 * @c operator*(Vec3) maps a point from B into A (@c R*point + p);
 * @ref inverse (also @c operator~) returns @c X_BA.
 */
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

/**
 * @brief A spatial 6-vector split into an angular and a linear @ref Vec3.
 * @note @c operator[] indexes @c [0] == angular, @c [1] == linear (SimTK
 *       convention). Used for spatial velocities, accelerations, and wrenches
 *       (angular = moment/torque, linear = force). The free @ref dot sums both
 *       halves.
 */
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

/**
 * @brief The full second mass moment (mass-weighted inertia tensor) about a point.
 *
 * Built either as an isotropic diagonal (@c Inertia(moment)) or as the inertia of
 * a point mass at an offset (@c Inertia(point, mass) == @c m(|p|^2 I - p p^T));
 * additive via @c operator+. Distinct from @ref MassProperties (which stores the
 * per-unit-mass @ref UnitInertia) and from @ref SpatialInertia (the 6x6 operator);
 * see the OQ-4 note in the module findings.
 */
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

/**
 * @brief Rigid-body spatial inertia about the body origin, as mass, mass center
 *        offset @c com (origin->COM), and @ref UnitInertia about the origin.
 *
 * Acts as the 6x6 spatial mass operator: @c operator*(SpatialVec) maps a spatial
 * velocity to spatial momentum (block form @c [[m*G, m*crossMat(com)],[~., m*I]],
 * delegated to @ref ArticulatedInertia for a single source of truth). Distinct
 * from @ref MassProperties and @ref Inertia; see the OQ-4 note in the findings.
 */
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

/**
 * @brief A symmetric spatial 6x6 inertia in three 3x3 blocks
 *        @c P = [[angAng, angLin],[~angLin, linLin]].
 *
 * Acting on a spatial vector: @c out.angular = angAng*w + angLin*v,
 * @c out.linear = ~angLin*w + linLin*v. The three-block constructor takes
 * @c (massBlock=linLin, momentBlock=angLin, inertiaBlock=angAng), matching the
 * engine's call order. @ref shift rigidly translates the inertia to a point
 * @c offset away (@c Phi(offset) P ~Phi(offset)); the articulated-body recursion
 * uses @c operator+= / @c operator- to assemble and remove child contributions.
 */
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

/**
 * @brief The rigid-body shift operator @c Phi(l) with offset @c l (parent origin
 *        -> child origin).
 *
 * @c Phi(l) * [a; b] == [a + l x b; b] shifts a spatial force/velocity across a
 * rigid offset; its transpose @ref PhiMatrixTranspose gives
 * @c ~Phi(l) * [a; b] == [a; b + a x l]. @c operator~ returns the transpose.
 */
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

/**
 * @brief Body mass properties: mass, mass center @c com, and @ref UnitInertia
 *        (inertia per unit mass) about the mass center.
 *
 * Constructing from a full @ref Inertia divides by mass to store the unit
 * inertia (and yields a zero tensor for a massless body). @ref reexpress rotates
 * the properties into a new frame (@c com' = ~R*com, unit inertia re-expressed);
 * @ref toSpatialInertia converts to the 6x6 @ref SpatialInertia operator. This is
 * the storage/description type; @ref SpatialInertia and @ref Inertia are the
 * operator and full-tensor forms (OQ-4, see findings).
 */
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