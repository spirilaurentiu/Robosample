#pragma once
/**
 * @file Constraints.hpp
 * @brief Holonomic loop-closure distance constraints for internal-coordinate
 *        dynamics: SHAKE position projection, RATTLE velocity projection, and
 *        the loop-closure Fixman log-determinant.
 *
 * A ring-closing bond is modelled as a fixed-distance constraint. The three
 * consumers - SHAKE, RATTLE, and the Fixman log-det - assemble the constraint
 * Jacobian @c G through one shared helper (@c assembleConstraintRow), so the
 * correction matches the projection exactly (INV-6). Projection lives in the
 * nu-dimensional generalized space under the mass metric
 * @f$ \delta q = -M^{-1} G^\top (G M^{-1} G^\top)^{-1} C @f$; the coupling matrix
 * @f$ G M^{-1} G^\top @f$ has size (number of constraints)^2 and is solved
 * directly.
 * @note @c ConstraintSet is owned by @c World by value and operates on a borrowed
 *       @c RobotState. Both SHAKE and RATTLE are required: either projection
 *       alone leaves a secular drift that pumps energy and opens the ring.
 */

#include <array>
#include <cmath>
#include <vector>

#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "robot_math.hpp"

namespace robo {

/// @brief One ring-closing bond as a fixed-distance (rod) constraint between two
///        atoms at rest length @c restLength (nm).
struct DistanceConstraint {
    int atomA{-1};
    int atomB{-1};
    Real restLength{0}; // d0 [nm]
};

/**
 * @brief The set of loop-closure distance constraints of one molecular system,
 *        with the SHAKE/RATTLE projections and the Fixman log-det over them.
 */
class ConstraintSet {
    public:
    std::vector<DistanceConstraint> distance;

    /// @return The number of loop-closure distance constraints.
    [[nodiscard]] auto numConstraints() const -> int {
        return static_cast<int>(distance.size());
    }
    /// @return @c true when there are no constraints (an acyclic system).
    [[nodiscard]] auto empty() const -> bool {
        return distance.empty();
    }

    /**
     * @brief Project per-atom Cartesian forces to the generalized force
     *        @c tau = J^T f by a rigid inward transmission sweep (kinematic, no
     *        inertia).
     * @param[in]  model     Immutable model.
     * @param[in]  state     Realized state (positions/kinematics current).
     * @param[in]  atomForce Per-atom Cartesian forces, length @c model.numAtoms.
     * @param[out] tauOut    Generalized force, length @c model.nu.
     * @note The sign follows the ForceBridge force-on-body convention.
     */
    static auto mapAtomForcesToGeneralizedForces(const RobotModel& model,
                                                 const RobotState& state,
                                                 const Vec3* atomForce,
                                                 Real* tauOut) -> void;

    /**
     * @brief RATTLE: project the generalized velocities so that @c G u = 0.
     * @param[in]     model Immutable model.
     * @param[in,out] state Velocities @c u are projected in place.
     * @pre @c realizeVelocity() is current for @p state.
     * @post @c G u == 0: no relative velocity along any constrained bond.
     */
    auto enforceVelocityConstraints(const RobotModel& model, RobotState& state) const -> void;

    /**
     * @brief Loop-closure contribution to the Fixman potential,
     *        @c ln det(G M^-1 G^T).
     * @param[in]     model Immutable model.
     * @param[in,out] state Realized state; used as scratch by @c multiplyByMInv,
     *                      not left with observable coordinate changes.
     * @return @c ln det(G M^-1 G^T) up to an additive constant that depends only
     *         on the (fixed) constraint count and cancels in any Metropolis
     *         @c dH. Returns exactly 0 when there are no loop-closure constraints
     *         (an acyclic molecule, whose flexible-coordinate determinant is the
     *         articulated-body tree determinant from @c RobotEngine::calcLogDetM
     *         alone).
     * @pre @c realizeArticulatedBodyInertias() is current for @p state (needed by
     *      @c multiplyByMInv).
     * @note @c G is assembled by the same @c assembleConstraintRow helper as
     *       SHAKE and RATTLE (INV-6), so this correction matches the projected
     *       subspace exactly. It supplies the ring term
     *       @c |M_flex| = |M_tree| / |G M^-1 G^T| that the tree determinant omits
     *       for cyclic systems (Spiridon & Minh 2017, JCTC 13:4649, Eqs. 2-3).
     */
    auto calcConstraintLogDet(const RobotModel& model, RobotState& state) const -> Real;

    /**
     * @brief SHAKE: Newton-iterate the generalized coordinates until every
     *        constraint @c C(q) = |rAB|^2 - d0^2 is satisfied to @p tolerance.
     * @tparam RefreshFn Callable @c void() that rebuilds atom positions and the
     *         kinematic Jacobian from the current @c q; invoked once per
     *         iteration because @c C is nonlinear.
     * @param[in]     model             Immutable model.
     * @param[in,out] state             Coordinates @c q are projected in place.
     * @param[in]     refreshKinematics Re-realizes positions/@c H from @c q.
     * @param[in]     tolerance         Convergence bound on the worst @c |C|
     *                                  (default 1e-10).
     * @param[in]     maxIter           Iteration cap (default 50).
     * @return The number of Newton iterations taken (0 when there are no
     *         constraints; @p maxIter if it did not converge).
     * @post On convergence every constrained bond length equals its rest length
     *       to @p tolerance.
     */
    template <typename RefreshFn>
    auto enforcePositionConstraints(const RobotModel& model,
                                    RobotState& state,
                                    RefreshFn&& refreshKinematics,
                                    Real tolerance = Real(1e-10),
                                    int maxIter = 50) const -> int {
        const int numC = numConstraints();
        if (numC == 0) {
            return 0;
        }
        const int numU = model.nu;

        std::vector<std::vector<Real>> jacobianT(static_cast<std::size_t>(numC), std::vector<Real>(numU, 0));
        std::vector<std::vector<Real>> minvJacobianT(static_cast<std::size_t>(numC),
                                                     std::vector<Real>(numU, 0));
        std::vector<Real> violation(static_cast<std::size_t>(numC), 0);
        std::vector<Vec3> atomForce(static_cast<std::size_t>(model.numAtoms), Vec3(0));

        int iter = 0;
        for (; iter < maxIter; ++iter) {
            refreshKinematics();

            Real worst = 0;
            for (int con = 0; con < numC; ++con) {
                const DistanceConstraint& bond = distance[con];
                const Vec3 vecAB = assembleConstraintRow(model, state, bond, atomForce,
                                                          jacobianT[con].data(), minvJacobianT[con].data());
                violation[con] = dot(vecAB, vecAB) - (bond.restLength * bond.restLength);
                worst = std::max(worst, std::abs(violation[con]));
            }
            if (worst <= tolerance) {
                break;
            }

            std::vector<Real> halfViolation(static_cast<std::size_t>(numC), 0);
            for (int con = 0; con < numC; ++con) {
                halfViolation[con] = Real(0.5) * violation[con]; // dC/dq = 2 G
            }
            std::vector<Real> lambda = solveCoupling(jacobianT, minvJacobianT, halfViolation, numC, numU);

            std::vector<Real> deltaVel(static_cast<std::size_t>(numU), 0);
            for (int con = 0; con < numC; ++con) {
                for (int idx = 0; idx < numU; ++idx) {
                    deltaVel[idx] -= lambda[con] * minvJacobianT[con][idx];
                }
            }
            applyVelSpaceIncrementToQ(model, state, deltaVel.data());
        }
        return iter;
    }

    /**
     * @brief Advance the generalized coordinates by a velocity-space increment,
     *        @c q <- q (+) deltaVel.
     * @param[in]     model    Immutable model.
     * @param[in,out] state    Coordinates @c q advanced in place.
     * @param[in]     deltaVel Increment in generalized-speed space, length
     *                         @c model.nu.
     * @post Scalar DOFs add directly; quaternion DOFs advance by the
     *       parent-frame map @c Delta(quat) = N(quat) * Delta(omega) and are
     *       renormalized, so quaternion coordinates stay unit (INV-9).
     */
    static auto applyVelSpaceIncrementToQ(const RobotModel& model, RobotState& state, const Real* deltaVel)
        -> void;

    private:
    // Test-only access to the private static solvers below (ConstraintTestAccess);
    // a friend declaration changes neither layout nor codegen.
    friend struct ConstraintTestAccess;

    /**
     * @brief Assemble one distance constraint's Jacobian row and its
     *        mass-weighted image (internal; the shared @c G assembly, INV-6).
     * @param[in]     model            Immutable model.
     * @param[in,out] state            Realized state (read; @c multiplyByMInv uses
     *                                 velocity scratch).
     * @param[in]     bond             The distance constraint.
     * @param[in,out] atomForceScratch Per-atom force scratch buffer.
     * @param[out]    jacobianTRow     @c tau_c = J^T vecAB, length @c model.nu.
     * @param[out]    minvJacobianTRow @c M^-1 tau_c, length @c model.nu.
     * @return @c vecAB (the bond vector A->B) so callers reuse it for their
     *         violation/rhs term without recomputing.
     * @note Called by RATTLE, SHAKE, and @c calcConstraintLogDet so all three use
     *       a byte-identical @c G (INV-6).
     */
    static auto assembleConstraintRow(const RobotModel& model,
                                      RobotState& state,
                                      const DistanceConstraint& bond,
                                      std::vector<Vec3>& atomForceScratch,
                                      Real* jacobianTRow,
                                      Real* minvJacobianTRow) -> Vec3;

    /**
     * @brief Solve the coupling system @c (G M^-1 G^T) x = rhs (internal).
     * @param[in]  jacobianT     Constraint Jacobian rows.
     * @param[in]  minvJacobianT Mass-weighted Jacobian rows.
     * @param[in]  rhs           Right-hand side, length @p numC.
     * @param[in]  numC          Constraint count.
     * @param[in]  numU          Generalized-speed count.
     * @param[out] logAbsDetOut  If non-null, receives @c ln|det(G M^-1 G^T)| as a
     *                           by-product of the factorization.
     * @return The solution vector @c x, length @p numC.
     */
    static auto solveCoupling(const std::vector<std::vector<Real>>& jacobianT,
                              const std::vector<std::vector<Real>>& minvJacobianT,
                              const std::vector<Real>& rhs,
                              int numC,
                              int numU,
                              Real* logAbsDetOut = nullptr) -> std::vector<Real>;

    /**
     * @brief Dense SPD solve with partial pivoting for a handful of constraints
     *        (internal).
     * @param[in]  matA         The coupling matrix @c G M^-1 G^T (SPD).
     * @param[in]  vecB         Right-hand side.
     * @param[out] logAbsDetOut If non-null, set to @c sum_k ln|pivot_k| ==
     *                          @c ln|det matA|.
     * @return The solution vector.
     */
    static auto solveSmallSpd(std::vector<std::vector<Real>> matA,
                              std::vector<Real> vecB,
                              Real* logAbsDetOut = nullptr) -> std::vector<Real>;
};

} // namespace robo