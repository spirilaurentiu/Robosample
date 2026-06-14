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
    // Build A = G M^-1 G^T (numC x numC, symmetric) and solve A x = rhs.
    static auto solveCoupling(const std::vector<std::vector<Real>>& jacobianT,
                              const std::vector<std::vector<Real>>& minvJacobianT,
                              const std::vector<Real>& rhs,
                              int numC,
                              int numU) -> std::vector<Real>;

    // Dense solve with partial pivoting (numC is a handful).
    static auto solveSmallSpd(std::vector<std::vector<Real>> matA, std::vector<Real> vecB)
        -> std::vector<Real>;
};

} // namespace robo