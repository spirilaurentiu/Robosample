#pragma once
// ============================================================================
//  Constraints.hpp -- holonomic loop-closure constraints for internal-coordinate
//  dynamics, the SimTK-free port of RBDistanceConstraint + LengthConstraints
//  (Simbody's Newton-Raphson loop solver).
//
//  Both position (SHAKE) and velocity (RATTLE) projection are required: SHAKE
//  pulls q back onto C(q)=0 after the position update; RATTLE removes the
//  relative-velocity-along-the-bond component after the velocity update. Either
//  alone leaves a secular drift that pumps energy and opens the ring -- this is
//  exactly Simbody's projectQ / projectU pair.
//
//  In internal coordinates atoms move only through q, so the projection lives in
//  the nu-dimensional generalized space using the mass metric:
//      dq = - M^-1 G^T (G M^-1 G^T)^-1 C.
//  G_c^T (an nu-vector) is the generalized force from Cartesian forces +rAB at A
//  and -rAB at B -- one inward rigid force-transmission sweep (forward
//  kinematics transposed). M^-1 G_c^T is one call to multiplyByMInv. For a few
//  rings the coupling matrix G M^-1 G^T is tiny and solved directly.
// ============================================================================

#include <array>
#include <cmath>
#include <vector>

#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "robot_math.hpp"

namespace robo {

// One ring-closing bond modelled as a constant-distance (Rod) constraint.
struct DistanceConstraint {
    int atomA{-1};
    int atomB{-1};
    Real restLength{0}; // d0 [nm]
};

class ConstraintSet {
    public:
    std::vector<DistanceConstraint> distance;

    [[nodiscard]] auto numConstraints() const -> int {
        return static_cast<int>(distance.size());
    }
    [[nodiscard]] auto empty() const -> bool {
        return distance.empty();
    }

    // tau = J^T * (spatial forces from per-atom Cartesian forces): rigid inward
    // transmission, purely kinematic (no inertia). Sign matches ForceBridge's
    // force-on-body convention (verify once against a finite-difference of C).
    static auto mapAtomForcesToGeneralizedForces(const RobotModel& model,
                                                 const RobotState& state,
                                                 const Vec3* atomForce,
                                                 Real* tauOut) -> void;

    // RATTLE: project u so that G u = 0. Requires realizeVelocity() current.
    auto enforceVelocityConstraints(const RobotModel& model, RobotState& state) const -> void;

    // Loop-closure contribution to the Fixman potential: ln det( G M^-1 G^T ).
    //
    // Why this exists (Spiridon & Minh 2017, JCTC 13:4649, Eqs. 2-3):
    //   The marginal density a constrained move samples is rho(phi_f) ~
    //   |M_{N_f}|^{1/2} exp(-beta U)  (their Eq. 2), and the Fixman potential
    //   U_F = (1/2) RT ln( |M_{N_f}| / |M_3N| )  (their Eq. 3) is what flattens it
    //   back to the unconstrained Cartesian Boltzmann marginal. |M_{N_f}| is the
    //   mass-metric determinant of the *flexible* coordinates of the constrained
    //   system.
    //
    //   For an ACYCLIC molecule the flexible coordinates are exactly the tree's
    //   torsions, and |M_{N_f}| is the articulated-body tree determinant computed
    //   by RobotEngine::calcLogDetM (Jain et al. O(n) algorithm, refs 9 & 23 in
    //   the paper -- derived for *branched* molecules, i.e. trees). No loop
    //   closures, so this function returns 0 and the tree determinant is the whole
    //   story. This is the regime the paper validated (C4 chain, butane, alanine
    //   dipeptide -- all acyclic).
    //
    //   For a CYCLIC molecule (a macrocycle: a ring closed by a RATTLE distance
    //   constraint) the spanning tree carries one *too many* flexible torsions --
    //   the loop-closure constraint removes one DOF per ring. Integrating the
    //   RATTLE-projected momenta (constrained to G M^-1 p = 0) over that reduced
    //   subspace contributes a factor det(G M^-1 G^T)^{-1/2} to the marginal, so
    //   the correct flexible-coordinate determinant is
    //         |M_{N_f}| = |M_tree| / |G M^-1 G^T|.
    //   Hence the full Fixman gains a -(1/2) RT ln det(G M^-1 G^T) term. The Jain
    //   tree algorithm does not supply it (it knows nothing about loop closures),
    //   and the paper never tested ring Boltzmann correctness (its macrocycle
    //   results, Sec. 3.5, measure only sampling *efficiency*), so this term has
    //   simply been absent for cyclic systems.
    //
    // Implementation note: G here is assembled exactly as in
    //   enforceVelocityConstraints (atomForce = vecAB = (1/2) dC/dr for the
    //   squared-distance constraint C = |rAB|^2 - d0^2). That (1/2) and the
    //   squared-vs-linear constraint convention put a *constant* multiplicative
    //   factor in det(.) -- a constant per fixed constraint count -- which cancels
    //   identically in any Metropolis dH. Only the configuration-dependent part of
    //   ln det(G M^-1 G^T) survives, and that is what this returns (up to that
    //   constant). Requires realizeArticulatedBodyInertias() current (needed by
    //   multiplyByMInv). Returns 0 when there are no loop-closure constraints.
    auto calcConstraintLogDet(const RobotModel& model, RobotState& state) const -> Real;

    // SHAKE: Newton-iterate q so that C(q) = |rAB|^2 - d0^2 = 0. refreshKinematics
    // rebuilds atom positions + H from q each iteration (C is nonlinear).
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
            const Vec3* posGround = state.atomPosG();

            Real worst = 0;
            for (int con = 0; con < numC; ++con) {
                const DistanceConstraint& bond = distance[con];
                const Vec3 vecAB = posGround[bond.atomA] - posGround[bond.atomB];
                violation[con] = dot(vecAB, vecAB) - (bond.restLength * bond.restLength);
                worst = std::max(worst, std::abs(violation[con]));

                for (Vec3& force : atomForce) {
                    force = Vec3(0);
                }
                atomForce[bond.atomA] = vecAB;
                atomForce[bond.atomB] = Vec3(0) - vecAB;
                mapAtomForcesToGeneralizedForces(model, state, atomForce.data(), jacobianT[con].data());
                RobotEngine::multiplyByMInv(model, state, jacobianT[con].data(), minvJacobianT[con].data());
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

    // q <- q (+) dvel : plain add for scalar DOFs; quaternion DOFs use
    // Delta(quat) = N(quat) * Delta(omega), then renormalize.
    static auto applyVelSpaceIncrementToQ(const RobotModel& model, RobotState& state, const Real* deltaVel)
        -> void;

    private:
    // Test-only access to the private static solvers below. A friend declaration
    // changes neither layout nor codegen; it lets the unit tests pin solveSmallSpd
    // / solveCoupling directly in addition to the public SHAKE/RATTLE/Fixman paths.
    friend struct ConstraintTestAccess;

    // Build A = G M^-1 G^T (numC x numC, symmetric) and solve A x = rhs.
    // If logAbsDetOut != nullptr it receives ln|det A| as a free by-product of the
    // factorization (used by calcConstraintLogDet for the loop-closure Fixman term).
    static auto solveCoupling(const std::vector<std::vector<Real>>& jacobianT,
                              const std::vector<std::vector<Real>>& minvJacobianT,
                              const std::vector<Real>& rhs,
                              int numC,
                              int numU,
                              Real* logAbsDetOut = nullptr) -> std::vector<Real>;

    // Dense solve with partial pivoting (numC is a handful). If logAbsDetOut !=
    // nullptr it is set to sum_k ln|pivot_k| = ln|det matA| (row swaps only flip
    // the sign of det; matA = G M^-1 G^T is SPD so |det| = product of |pivots|).
    static auto solveSmallSpd(std::vector<std::vector<Real>> matA,
                              std::vector<Real> vecB,
                              Real* logAbsDetOut = nullptr) -> std::vector<Real>;
};

} // namespace robo