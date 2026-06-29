// ============================================================================
//  test_linalg.cpp -- Part 1.8: jacobiSymEig, invertDense, symSqrt(Inv),
//  logDetSymPD. Exercises the REAL helpers (extracted into RobotLinearAlgebra.hpp).
//  Self-consistency + crafted singular/conditioned matrices.
// ============================================================================
#include <cmath>

#include "RobotLinearAlgebra.hpp"
#include "TestHelpers.hpp"

using namespace rtest;
namespace L = robo_linalg;

// --- jacobiSymEig: V diag(d) V^T = A and V^T V = I --------------------------
TEST(LinAlg, JacobiEigReconstructs) {
    Rng rng(1);
    for (int n : {1, 3, 6}) {
        for (int t = 0; t < 100; ++t) {
            Real A[36];
            rng.spd(n, A, 0.3, 5.0);
            Real d[6], V[36];
            L::jacobiSymEig(A, n, d, V);
            // reconstruct
            Real recon[36];
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    Real acc = 0;
                    for (int k = 0; k < n; ++k) {
                        acc += V[i * n + k] * d[k] * V[j * n + k];
                    }
                    recon[i * n + j] = acc;
                }
            }
            EXPECT_LT(matDiff(recon, A, n), kLoose) << "n=" << n;
            // V^T V = I
            Real VtV[36], I[36];
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    Real acc = 0;
                    for (int k = 0; k < n; ++k) {
                        acc += V[k * n + i] * V[k * n + j];
                    }
                    VtV[i * n + j] = acc;
                }
            }
            identity(n, I);
            EXPECT_LT(matDiff(VtV, I, n), kLoose) << "eigenvectors not orthonormal, n=" << n;
        }
    }
}

// --- invertDense: A A^-1 = I for non-degenerate -----------------------------
TEST(LinAlg, InvertDenseInverts) {
    Rng rng(2);
    for (int n : {1, 3, 6}) {
        for (int t = 0; t < 100; ++t) {
            Real A[36];
            rng.spd(n, A, 0.5, 5.0);
            Real Ai[36], prod[36], I[36];
            L::invertDense(A, n, Ai);
            matMat(A, Ai, n, prod);
            identity(n, I);
            EXPECT_LT(matDiff(prod, I, n), kLoose) << "n=" << n;
        }
    }
}

// --- invertDense LOCKS a null direction (returns 0, not 1e14) ---------------
TEST(LinAlg, InvertDenseLocksNullDirection) {
    // 1-DOF: a tiny "inertia" on its own axis must lock to 0.
    Real d = 1e-33, di;
    L::invertDense(&d, 1, &di);
    EXPECT_EQ(di, 0.0) << "near-null 1-DOF must lock to 0, not 1/eps";

    // 3x3 with one ~null eigenvalue: that direction's inverse contributes 0.
    Rng rng(3);
    Real Q[9];
    rng.spd(3, Q, 1.0, 1.0); // orthonormal-ish basis with eqvals 1 (well-cond)
    // build A = sum lambda_k v_k v_k^T with lambda = {2, 3, 1e-20}
    Real d3[3], V[9];
    L::jacobiSymEig(Q, 3, d3, V);
    const Real lam[3] = {2.0, 3.0, 1e-20};
    Real A[9];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Real acc = 0;
            for (int k = 0; k < 3; ++k) {
                acc += V[i * 3 + k] * lam[k] * V[j * 3 + k];
            }
            A[i * 3 + j] = acc;
        }
    }
    Real Ai[9];
    L::invertDense(A, 3, Ai);
    // Ai must be finite and bounded (the null direction was dropped, not 1e20)
    for (int i = 0; i < 9; ++i) {
        EXPECT_TRUE(std::isfinite(Ai[i]));
        EXPECT_LT(std::abs(Ai[i]), 1e6) << "null direction was regularized to a huge inverse";
    }
}

// --- tolerance straddle: stiff-but-nonsingular must NOT lock (see 12.6) ------
TEST(LinAlg, InvertDenseToleranceStraddle) {
    // eigenvalue 1e-6 is small but >> 1e-12 lock tol: must invert to ~1e6.
    Real d = 1e-6, di;
    L::invertDense(&d, 1, &di);
    EXPECT_NEAR(di, 1e6, 1.0);
    // eigenvalue 1e-13 < 1e-12: must lock to 0.
    Real d2 = 1e-13, di2;
    L::invertDense(&d2, 1, &di2);
    EXPECT_EQ(di2, 0.0);
}

// --- symSqrt: S*S = A ; symSqrtInv*symSqrt = I ------------------------------
TEST(LinAlg, SymSqrtSquaresToA) {
    Rng rng(4);
    for (int n : {1, 3, 6}) {
        for (int t = 0; t < 100; ++t) {
            Real A[36];
            rng.spd(n, A, 0.3, 5.0);
            Real S[36], SS[36];
            L::symSqrt(A, n, S);
            matMat(S, S, n, SS);
            EXPECT_LT(matDiff(SS, A, n), 1e-11) << "n=" << n;

            Real Sinv[36], prod[36], I[36];
            L::symSqrtInv(A, n, Sinv);
            matMat(Sinv, S, n, prod);
            identity(n, I);
            EXPECT_LT(matDiff(prod, I, n), 1e-11) << "symSqrtInv*symSqrt != I, n=" << n;
        }
    }
}

// --- logDetSymPD: = log det (via eigenvalue product) ------------------------
TEST(LinAlg, LogDetMatchesEigenProduct) {
    Rng rng(5);
    for (int n : {1, 3, 6}) {
        for (int t = 0; t < 100; ++t) {
            Real A[36];
            rng.spd(n, A, 0.4, 4.0);
            const Real ld = L::logDetSymPD(A, n);
            Real d[6], V[36];
            L::jacobiSymEig(A, n, d, V);
            Real ref = 0;
            for (int k = 0; k < n; ++k) {
                ref += std::log(d[k]);
            }
            EXPECT_NEAR(ld, ref, kLoose) << "n=" << n;
        }
    }
}

// --- logDetSymPD pivot clamp: near-singular -> finite, not NaN --------------
TEST(LinAlg, LogDetPivotClampFinite) {
    // a 2x2 that is (numerically) singular
    Real A[4] = {1e-200, 0, 0, 1e-200};
    const Real ld = L::logDetSymPD(A, 2);
    EXPECT_TRUE(std::isfinite(ld)) << "pivot clamp must keep log-det finite";
}
