// ============================================================================
//  TestConstraintSolver.cpp -- Phase 0.3: the private dense coupling solvers of
//  ConstraintSet, reached via the friend ConstraintTestAccess.
//
//  solveSmallSpd(A, b, &lnDet) solves A x = b for a small SPD A (Gauss-Jordan,
//  partial pivot) and returns ln|det A| as a by-product. solveCoupling assembles
//  A = G M^-1 G^T from (G^T rows, M^-1 G^T rows) and forwards to solveSmallSpd.
//  Both are load-bearing for SHAKE/RATTLE (the Lagrange-multiplier solve) and
//  for calcConstraintLogDet (the loop-closure Fixman term), so they get a direct
//  test here in addition to the public-path coverage planned for TestConstraints.
//
//  Cross-link: the ln|det A| solveSmallSpd reports must match the independent
//  Cholesky log-determinant robo_linalg::logDetSymPD (the freshly hoisted Phase
//  0.1 symbol) for SPD inputs.
// ============================================================================
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "ConstraintTestAccess.hpp"
#include "RobotLinearAlgebra.hpp"
#include "TestHelpers.hpp"

using robo::ConstraintTestAccess;
using robo::Real;
using rtest::Rng;

namespace {

// dense A (vector-of-rows) times x.
std::vector<Real> matVecRows(const std::vector<std::vector<Real>>& A, const std::vector<Real>& x) {
    const int n = static_cast<int>(x.size());
    std::vector<Real> y(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        Real s = 0;
        for (int j = 0; j < n; ++j) {
            s += A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] * x[static_cast<std::size_t>(j)];
        }
        y[static_cast<std::size_t>(i)] = s;
    }
    return y;
}

Real maxResidual(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        d = std::max(d, std::abs(a[i] - b[i]));
    }
    return d;
}

} // namespace

// ---------------------------------------------------------------------------
//  1. 2x2 SPD golden: known solution and known ln det.
// ---------------------------------------------------------------------------
TEST(ConstraintSolver, SolveSmallSpd2x2Golden) {
    // A = [[4,1],[1,3]], det = 11.  b = [1,2].
    std::vector<std::vector<Real>> A = {{4, 1}, {1, 3}};
    std::vector<Real> b = {1, 2};
    Real lnDet = 0;
    const std::vector<Real> x = ConstraintTestAccess::solveSmallSpd(A, b, &lnDet);

    // A x == b
    EXPECT_LT(maxResidual(matVecRows(A, x), b), rtest::kLoose);
    // closed-form: x = A^-1 b = (1/11)[[3,-1],[-1,4]] [1;2] = [1/11; 7/11]
    EXPECT_NEAR(x[0], 1.0 / 11.0, rtest::kLoose);
    EXPECT_NEAR(x[1], 7.0 / 11.0, rtest::kLoose);
    EXPECT_NEAR(lnDet, std::log(11.0), rtest::kLoose);
}

// ---------------------------------------------------------------------------
//  2. 1x1 trivial: x = b/A, ln det = ln A.
// ---------------------------------------------------------------------------
TEST(ConstraintSolver, SolveSmallSpd1x1) {
    std::vector<std::vector<Real>> A = {{7.0}};
    std::vector<Real> b = {3.5};
    Real lnDet = 0;
    const std::vector<Real> x = ConstraintTestAccess::solveSmallSpd(A, b, &lnDet);
    EXPECT_NEAR(x[0], 0.5, rtest::kTight);
    EXPECT_NEAR(lnDet, std::log(7.0), rtest::kTight);
}

// ---------------------------------------------------------------------------
//  3. Degenerate pivot stays FINITE (must not NaN): a zero diagonal is floored,
//     the log stays finite, and the corresponding solution component locks to 0.
//     (Mirrors the loop-closure case where a constraint row is rank-deficient.)
// ---------------------------------------------------------------------------
TEST(ConstraintSolver, SolveSmallSpdDegeneratePivotFinite) {
    std::vector<std::vector<Real>> A = {{0.0, 0.0}, {0.0, 2.0}};
    std::vector<Real> b = {1.0, 4.0};
    Real lnDet = 0;
    const std::vector<Real> x = ConstraintTestAccess::solveSmallSpd(A, b, &lnDet);
    EXPECT_TRUE(std::isfinite(lnDet)) << "degenerate pivot must keep ln|det| finite";
    for (Real xi : x) {
        EXPECT_TRUE(std::isfinite(xi));
    }
    EXPECT_EQ(x[0], 0.0) << "null direction must lock to 0, not blow up";
    EXPECT_NEAR(x[1], 2.0, rtest::kLoose); // the well-conditioned row still solves
}

// ---------------------------------------------------------------------------
//  4. solveSmallSpd ln|det| == independent Cholesky log-det (robo_linalg), for
//     random SPD matrices of size 1..6. Ties the constraint solver's determinant
//     by-product to the Phase 0.1 hoisted kernel.
// ---------------------------------------------------------------------------
TEST(ConstraintSolver, LnDetMatchesCholeskyLogDet) {
    Rng rng(0xC0501);
    for (int n : {1, 2, 3, 6}) {
        for (int t = 0; t < 50; ++t) {
            Real flat[36];
            rng.spd(n, flat, 0.4, 4.0); // row-major SPD
            std::vector<std::vector<Real>> A(static_cast<std::size_t>(n),
                                             std::vector<Real>(static_cast<std::size_t>(n)));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = flat[i * n + j];
                }
            }
            std::vector<Real> b(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) {
                b[static_cast<std::size_t>(i)] = rng.gaussian();
            }
            Real lnDet = 0;
            const std::vector<Real> x = ConstraintTestAccess::solveSmallSpd(A, b, &lnDet);
            EXPECT_LT(maxResidual(matVecRows(A, x), b), 1e-7) << "n=" << n;
            EXPECT_NEAR(lnDet, robo_linalg::logDetSymPD(flat, n), rtest::kLoose) << "n=" << n;
        }
    }
}

// ---------------------------------------------------------------------------
//  5. solveCoupling assembles A = G M^-1 G^T correctly: with G^T = I (jacobianT
//     identity rows) and M^-1 G^T = S (minvJacobianT = S rows), the coupling
//     matrix is exactly S, so the solve must match solveSmallSpd(S, rhs). Pins
//     the assembly AND proves solveCoupling is reachable.
// ---------------------------------------------------------------------------
TEST(ConstraintSolver, SolveCouplingAssemblesGMInvGt) {
    Rng rng(0xC0502);
    for (int n : {1, 2, 3}) {
        Real flat[9];
        rng.spd(n, flat, 0.5, 3.0);
        // jacobianT = identity (n rows of length numU=n); minvJacobianT = S rows.
        std::vector<std::vector<Real>> jacT(static_cast<std::size_t>(n),
                                            std::vector<Real>(static_cast<std::size_t>(n), 0));
        std::vector<std::vector<Real>> minvJacT(static_cast<std::size_t>(n),
                                                std::vector<Real>(static_cast<std::size_t>(n), 0));
        std::vector<std::vector<Real>> S(static_cast<std::size_t>(n),
                                         std::vector<Real>(static_cast<std::size_t>(n)));
        for (int i = 0; i < n; ++i) {
            jacT[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1.0;
            for (int j = 0; j < n; ++j) {
                const Real v = flat[i * n + j];
                minvJacT[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = v;
                S[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = v;
            }
        }
        std::vector<Real> rhs(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            rhs[static_cast<std::size_t>(i)] = rng.gaussian();
        }

        Real lnDetCoupling = 0, lnDetDirect = 0;
        const std::vector<Real> viaCoupling =
            ConstraintTestAccess::solveCoupling(jacT, minvJacT, rhs, n, n, &lnDetCoupling);
        const std::vector<Real> viaDirect = ConstraintTestAccess::solveSmallSpd(S, rhs, &lnDetDirect);

        EXPECT_LT(maxResidual(viaCoupling, viaDirect), rtest::kLoose) << "n=" << n;
        EXPECT_NEAR(lnDetCoupling, lnDetDirect, rtest::kLoose) << "n=" << n;
    }
}