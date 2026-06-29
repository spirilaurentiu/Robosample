#pragma once
// ============================================================================
//  ConstraintTestAccess.hpp -- Phase 0.3: test-only window onto the two private
//  static solvers of robo::ConstraintSet (solveSmallSpd / solveCoupling).
//
//  ConstraintSet grants `friend struct ConstraintTestAccess;` (include/
//  Constraints.hpp); this struct forwards to the private statics so the unit
//  tests can pin the dense SPD solve and its ln|det| by-product DIRECTLY, in
//  addition to exercising them through the public SHAKE/RATTLE/Fixman paths.
//  No production behaviour changes -- a friend declaration alters neither layout
//  nor codegen.
// ============================================================================

#include <vector>

#include "Constraints.hpp"
#include "robot_math.hpp"

namespace robo {

struct ConstraintTestAccess {
    using Real = robo::Real;

    // Solve A x = b for a small dense SPD A (Gauss-Jordan, partial pivot); if
    // logAbsDetOut != nullptr it receives sum_k ln|pivot_k| = ln|det A|.
    static auto solveSmallSpd(std::vector<std::vector<Real>> matA,
                              std::vector<Real> vecB,
                              Real* logAbsDetOut = nullptr) -> std::vector<Real> {
        return ConstraintSet::solveSmallSpd(std::move(matA), std::move(vecB), logAbsDetOut);
    }

    // Assemble A = G M^-1 G^T from (G^T rows, M^-1 G^T rows) and solve A x = rhs;
    // logAbsDetOut receives ln|det A| from the same factorization.
    static auto solveCoupling(const std::vector<std::vector<Real>>& jacobianT,
                              const std::vector<std::vector<Real>>& minvJacobianT,
                              const std::vector<Real>& rhs,
                              int numC,
                              int numU,
                              Real* logAbsDetOut = nullptr) -> std::vector<Real> {
        return ConstraintSet::solveCoupling(jacobianT, minvJacobianT, rhs, numC, numU, logAbsDetOut);
    }
};

} // namespace robo