#pragma once

/**
 * @file BatScaling.hpp
 * @brief Deterministic BAT (bond/angle) position-scaling drive and its Cartesian
 *        log-Jacobian, used by the nonequilibrium replica-exchange move (RENE).
 *
 * A scale factor @c s displaces the scaled bond/angle coordinates of every
 * eligible body about a shared anchor @c mu, moving each body's subtree as one
 * rigid assembly so that its own @c (r,theta) relative to its parent become
 * @f$ (s r - (s-1)\mu_r,\; s\theta - (s-1)\mu_\theta) @f$. Only Slider/Cylinder
 * (bond length @c r) and BendStretch/SphericalCoords (@c r and bond angle
 * @c theta) expose a scaled DOF; all other joints are left unmoved. The map and
 * its log-Jacobian read the identical @c (r,theta) geometry from the z-matrix, so
 * they stay mutually consistent (INV-7).
 *
 * @note This module is the drive and Jacobian only; it performs no acceptance,
 *       no work accumulation, and no run-type dispatch.
 */

#include <unordered_map>
#include <vector>

#include "RobotModel.hpp"
#include "robot_math.hpp"

namespace robo {

/// @brief Which BAT coordinates a joint exposes to scaling: bond length @c r
///        and/or bond angle @c theta. Torsion/azimuth and root DOFs are never
///        scaled.
struct ScaledDofFlags {
    bool r = false;
    bool theta = false;
};

/**
 * @brief The scaled BAT coordinates of a joint type.
 * @param[in] jt Joint type.
 * @return @c {r=true} for Slider and Cylinder (bond length only; the rotation is
 *         a torsion and excluded); @c {r=true,theta=true} for BendStretch and
 *         SphericalCoords; @c {false,false} for all other joints.
 */
[[nodiscard]] constexpr ScaledDofFlags batScaledDofs(JointType jt) {
    switch (jt) {
        case JointType::Slider:
            return {true, false}; // r only
        case JointType::Cylinder:
            return {true, false}; // r only (rotation DOF is TORSION, excluded)
        case JointType::BendStretch:
            return {true, true}; // theta, r
        case JointType::SphericalCoords:
            return {true, true}; // zenith theta, radius r (azimuth excluded)
        default:
            return {false, false};
    }
}

/**
 * @brief Topology-only count of scaled bond+angle DOFs over the whole model.
 * @param[in] model Immutable model.
 * @return The sum over bodies of scaled-DOF contributions (Slider/Cylinder = 1,
 *         BendStretch/SphericalCoords = 2, others 0).
 * @warning This is the STATIC count. The count of DOFs actually displaced by
 *          @ref applyBatScaling may be lower in a degenerate geometry; the
 *          log-Jacobian SHALL use @ref applyBatScaling's @c outNScaled (what the
 *          map did), not this static count.
 */
[[nodiscard]] int countScaledDofsStatic(const RobotModel& model);

/**
 * @brief Cartesian volume log-term @c J(x) = sum over scaled bodies of
 *        @c [2 ln r + ln sin(theta)], evaluated on a given position snapshot.
 * @param[in] model   Immutable model (supplies the z-matrix).
 * @param[in] atomPos Cartesian positions to evaluate on; typically a pre- or
 *                    post-scale copy, not required to be the live state.
 * @return The sum of @c [2 ln r + ln sin(theta)] over scaled bodies, using only
 *         the coordinate(s) each body's joint exposes. A body whose z-row is
 *         missing or whose @c zJ/zK ancestor is absent (root-adjacent)
 *         contributes 0.
 */
[[nodiscard]] double calcBatVolumeLogJac(const RobotModel& model, const std::vector<Vec3>& atomPos);

/**
 * @brief Apply the deterministic BAT-scaling map by factor @p s and return both
 *        the scaled positions and the log-Jacobian the map actually realized.
 * @param[in]  model       Immutable model (topology and z-matrix).
 * @param[in]  atomPosIn   Input Cartesian positions; not mutated.
 * @param[in]  s           Scale factor applied to the scaled BAT coordinates.
 * @param[in]  anchorR     Per-body bond-length anchor @c mu_r, keyed by the
 *                         scaled body's placed (@c zI) atom index; a missing key
 *                         defaults to 0.
 * @param[in]  anchorTheta Per-body bond-angle anchor @c mu_theta, same keying.
 * @param[out] outNScaled  Count of DOFs the map actually displaced (may be below
 *                         @ref countScaledDofsStatic in degenerate geometry).
 * @param[out] outLnJac    The Cartesian log-Jacobian @c (J(x')-J(x0)) +
 *                         @c outNScaled*ln(s), consistent with @p outNScaled.
 * @return A scaled copy of the positions; @p atomPosIn is unchanged.
 * @post @p outLnJac is the log-Jacobian of the geometry actually produced, so it
 *       matches the map (INV-7); @p outNScaled is the count used to form it.
 * @warning The paired map is an exact involution
 *          (@c M_{1/s,mu} o M_{s,mu} == identity) ONLY when both legs share the
 *          SAME anchor @c mu. A caller SHALL pass the SAME frozen
 *          @p anchorR / @p anchorTheta to both partners of a swap pair in one
 *          round (INV-7); see @ref BatAnchorStats.
 * @note @p outLnJac depends on @c mu because the produced geometry
 *       @c x' = f(x0, s, mu) does; only the map's local determinant @c s^N is
 *       anchor-independent.
 * @note Limitation: a scaled body whose parent is itself Ground-adjacent has no
 *       @c zK reference atom, so its angle DOF is silently skipped (omitted from
 *       both @p outNScaled and @p outLnJac, kept self-consistent). Production
 *       configurations SHOULD avoid stacking a scaled joint directly on a
 *       root-adjacent body. See the module findings.
 */
std::vector<Vec3> applyBatScaling(const RobotModel& model,
                                  const std::vector<Vec3>& atomPosIn,
                                  double s,
                                  const std::unordered_map<int, double>& anchorR,
                                  const std::unordered_map<int, double>& anchorTheta,
                                  int& outNScaled,
                                  double& outLnJac);

/**
 * @brief Accumulator for the shared, state-independent BAT scaling anchor @c mu
 *        that @ref applyBatScaling requires (INV-7).
 *
 * Holds one running mean per scaled BAT coordinate (keyed by the scaled body's
 * @c zI atom index, the key @ref applyBatScaling reads), updated from
 * equilibrium samples across all thermodynamic states. The anchor SHALL be a
 * single value shared by both partners of a swap pair and frozen across their
 * forward and reverse drives; otherwise the paired map is not an exact
 * involution. Owned by @c Context, not per thermodynamic state.
 */
class BatAnchorStats {
    public:
    /**
     * @brief Fold one configuration's scaled @c (r, theta) into the running
     *        means (Welford online update).
     * @param[in] model   Immutable model.
     * @param[in] atomPos Configuration to sample; one @c r (and @c theta if that
     *                    body scales an angle) per scaled body.
     */
    void accumulate(const RobotModel& model, const std::vector<Vec3>& atomPos);

    /// @brief Immutable per-coordinate anchor means keyed by scaled-body @c zI.
    struct Snapshot {
        std::unordered_map<int, double> meanR;
        std::unordered_map<int, double> meanTheta;
    };
    /**
     * @brief Freeze the current running means into an immutable snapshot.
     * @return The current per-coordinate means.
     * @warning Take exactly ONE snapshot per round and pass that same object to
     *          both partners' drives that round, so the paired map is an exact
     *          involution (INV-7).
     */
    [[nodiscard]] Snapshot snapshot() const;

    /// @brief Discard all accumulated statistics.
    void reset();

    private:
    struct RunningMean {
        long long n = 0;
        double mean = 0.0;
    };
    std::unordered_map<int, RunningMean> meanR_;
    std::unordered_map<int, RunningMean> meanTheta_;
};

} // namespace robo
