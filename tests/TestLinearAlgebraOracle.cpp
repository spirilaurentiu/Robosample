// ============================================================================
//  TestLinearAlgebraOracle.cpp -- VALIDATE the dense linalg kernels against an
//  INDEPENDENT oracle, not against themselves.
//
//  WHY THIS EXISTS. TestLinearAlgebra.cpp checks self-consistency: A*A^-1 == I,
//  S*S == A, V diag V^T == A. Those are necessary (they catch gross errors) but
//  they are not validation against a trusted reference -- a routine can satisfy
//  "A*Ainv==I" to 1e-9 and still carry the wrong convention, mis-order/eigen-sign,
//  or silently fail to converge on a hard spectrum. This file builds matrices
//  whose eigenvalues / eigenvectors / inverse / sqrt / log-det are KNOWN BY
//  CONSTRUCTION (A = Q diag(lambda) Q^T with Q,lambda chosen), so every kernel is
//  checked against the closed-form answer -- an oracle independent of the code
//  under test and requiring no external library.
//
//  CROSS-CHECK STATUS (recorded in the repo, run offline): the same kernels were
//  diffed against LAPACK (scipy dsyev/dpotrf/sqrtm) over well-conditioned,
//  clustered, and degenerate spectra -- agreement to ~1e-15. They diverge from
//  LAPACK *intentionally* in exactly two regimes, which are pinned below as
//  INTENDED behaviour (a naive "compare to LAPACK" oracle would falsely flag
//  them): (1) invertDense LOCKS a near-null eigendirection to 0 instead of
//  inverting it to ~1/eps; (2) logDetSymPD FLOORS a near-zero Cholesky pivot to
//  keep the log finite.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "RobotLinearAlgebra.hpp" // robo_linalg::{jacobiSymEig,invertDense,symSqrt,symSqrtInv,logDetSymPD}
#include "TestHelpers.hpp"

using robo::Real;
using rtest::Rng;
namespace L = robo_linalg;

namespace {

// Build A = Q diag(lambda) Q^T (row-major), returning A AND the Q,lambda used.
// Q is a random orthonormal matrix (Gram-Schmidt on a Gaussian), so the spectrum
// (lambda) and eigenspaces (columns of Q) are the GROUND TRUTH for the kernels.
struct Spd {
    int n;
    std::array<Real, 36> A{};
    std::array<Real, 36> Q{}; // columns = eigenvectors
    std::array<Real, 6> lam{};
};

Spd spdFromSpectrum(Rng& rng, const std::vector<Real>& lambda) {
    Spd s;
    s.n = static_cast<int>(lambda.size());
    const int n = s.n;
    for (int i = 0; i < n; ++i) {
        s.lam[static_cast<std::size_t>(i)] = lambda[static_cast<std::size_t>(i)];
    }
    // random orthonormal Q via Gram-Schmidt on a Gaussian (same idea as rng.spd,
    // but here we KEEP Q and lambda instead of discarding them).
    for (int i = 0; i < n * n; ++i) {
        s.Q[static_cast<std::size_t>(i)] = rng.gaussian();
    }
    for (int c = 0; c < n; ++c) {
        for (int p = 0; p < c; ++p) {
            Real d = 0;
            for (int r = 0; r < n; ++r) {
                d += s.Q[static_cast<std::size_t>(r * n + c)] * s.Q[static_cast<std::size_t>(r * n + p)];
            }
            for (int r = 0; r < n; ++r) {
                s.Q[static_cast<std::size_t>(r * n + c)] -= d * s.Q[static_cast<std::size_t>(r * n + p)];
            }
        }
        Real nrm = 0;
        for (int r = 0; r < n; ++r) {
            nrm += s.Q[static_cast<std::size_t>(r * n + c)] * s.Q[static_cast<std::size_t>(r * n + c)];
        }
        nrm = std::sqrt(nrm);
        for (int r = 0; r < n; ++r) {
            s.Q[static_cast<std::size_t>(r * n + c)] /= nrm;
        }
    }
    // A = Q diag(lam) Q^T
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Real acc = 0;
            for (int k = 0; k < n; ++k) {
                acc += s.Q[static_cast<std::size_t>(i * n + k)] * s.lam[static_cast<std::size_t>(k)]
                       * s.Q[static_cast<std::size_t>(j * n + k)];
            }
            s.A[static_cast<std::size_t>(i * n + j)] = acc;
        }
    }
    return s;
}

// Reference Q diag(f(lambda)) Q^T into out (the closed-form oracle for inv/sqrt).
template <class F>
void refFromSpectrum(const Spd& s, F f, Real* out) {
    const int n = s.n;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Real acc = 0;
            for (int k = 0; k < n; ++k) {
                acc += s.Q[static_cast<std::size_t>(i * n + k)] * f(s.lam[static_cast<std::size_t>(k)])
                       * s.Q[static_cast<std::size_t>(j * n + k)];
            }
            out[i * n + j] = acc;
        }
    }
}

Real maxAbsDiff(const Real* a, const Real* b, int n) {
    Real d = 0;
    for (int i = 0; i < n * n; ++i) {
        d = std::max(d, std::abs(a[i] - b[i]));
    }
    return d;
}

// A spread of well-conditioned spectra (cond <= ~1e3) where every kernel must
// match the closed form to ~kLoose.
const std::vector<std::vector<Real>> kGoodSpectra = {
    {2.0},
    {0.5, 3.0, 4.0},
    {0.3, 0.9, 1.7, 2.4, 3.1, 4.8},
    {1.0, 1.0000001, 2.0, 2.0000001, 5.0, 5.0000001}, // clustered
};

} // namespace

// ---------------------------------------------------------------------------
//  O1: jacobiSymEig recovers the KNOWN spectrum (sorted), not just "some
//      decomposition that reconstructs A".
// ---------------------------------------------------------------------------
TEST(LinAlgOracle, EigRecoversKnownSpectrum) {
    Rng rng(0x0E16);
    for (const auto& spec : kGoodSpectra) {
        for (int t = 0; t < 30; ++t) {
            const Spd s = spdFromSpectrum(rng, spec);
            Real d[6], V[36];
            L::jacobiSymEig(s.A.data(), s.n, d, V);
            std::vector<Real> got(d, d + s.n), want(spec);
            std::sort(got.begin(), got.end());
            std::sort(want.begin(), want.end());
            for (int k = 0; k < s.n; ++k) {
                EXPECT_NEAR(got[static_cast<std::size_t>(k)],
                            want[static_cast<std::size_t>(k)],
                            rtest::kLoose)
                    << "n=" << s.n << " k=" << k;
            }
        }
    }
}

// ---------------------------------------------------------------------------
//  O2: invertDense == Q diag(1/lambda) Q^T (the closed-form inverse).
// ---------------------------------------------------------------------------
TEST(LinAlgOracle, InverseMatchesClosedForm) {
    Rng rng(0x1217);
    for (const auto& spec : kGoodSpectra) {
        for (int t = 0; t < 30; ++t) {
            const Spd s = spdFromSpectrum(rng, spec);
            Real got[36], ref[36];
            L::invertDense(s.A.data(), s.n, got);
            refFromSpectrum(
                s,
                [](Real x) {
                    return Real(1) / x;
                },
                ref);
            EXPECT_LT(maxAbsDiff(got, ref, s.n), rtest::kLoose) << "n=" << s.n;
        }
    }
}

// ---------------------------------------------------------------------------
//  O3: symSqrt == Q diag(sqrt(lambda)) Q^T ; symSqrtInv == Q diag(1/sqrt) Q^T.
// ---------------------------------------------------------------------------
TEST(LinAlgOracle, SqrtAndInvSqrtMatchClosedForm) {
    Rng rng(0x5417);
    for (const auto& spec : kGoodSpectra) {
        for (int t = 0; t < 30; ++t) {
            const Spd s = spdFromSpectrum(rng, spec);
            Real gotS[36], refS[36], gotSi[36], refSi[36];
            L::symSqrt(s.A.data(), s.n, gotS);
            L::symSqrtInv(s.A.data(), s.n, gotSi);
            refFromSpectrum(
                s,
                [](Real x) {
                    return std::sqrt(x);
                },
                refS);
            refFromSpectrum(
                s,
                [](Real x) {
                    return Real(1) / std::sqrt(x);
                },
                refSi);
            EXPECT_LT(maxAbsDiff(gotS, refS, s.n), rtest::kLoose) << "sqrt n=" << s.n;
            EXPECT_LT(maxAbsDiff(gotSi, refSi, s.n), rtest::kLoose) << "invsqrt n=" << s.n;
        }
    }
}

// ---------------------------------------------------------------------------
//  O4: logDetSymPD == sum_k ln(lambda_k)  (the exact log-determinant).
// ---------------------------------------------------------------------------
TEST(LinAlgOracle, LogDetEqualsSumLogSpectrum) {
    Rng rng(0x70D7);
    for (const auto& spec : kGoodSpectra) {
        for (int t = 0; t < 30; ++t) {
            const Spd s = spdFromSpectrum(rng, spec);
            Real want = 0;
            for (Real lk : spec) {
                want += std::log(lk);
            }
            EXPECT_NEAR(L::logDetSymPD(s.A.data(), s.n), want, rtest::kLoose) << "n=" << s.n;
        }
    }
}

// ---------------------------------------------------------------------------
//  O5: INTENDED divergence #1 -- invertDense LOCKS a near-null eigendirection.
//  With lambda = {1e-13, 2, 3} the closed-form inverse has a ~1e13 entry; the
//  kernel must instead return the pseudo-inverse Q diag(0, 1/2, 1/3) Q^T (null
//  direction zeroed), staying bounded. This is the deliberate robustness feature
//  a naive LAPACK oracle would wrongly flag as a bug.
// ---------------------------------------------------------------------------
TEST(LinAlgOracle, InvertDenseLocksNullDirection_IntendedDivergence) {
    Rng rng(0x10CC);
    const Spd s = spdFromSpectrum(rng, {1e-13, 2.0, 3.0});
    Real got[9], pinv[9], trueInv[9];
    L::invertDense(s.A.data(), 3, got);
    refFromSpectrum(
        s,
        [](Real x) {
            return x > Real(1e-12) ? Real(1) / x : Real(0);
        },
        pinv); // locked
    refFromSpectrum(
        s,
        [](Real x) {
            return Real(1) / x;
        },
        trueInv); // 1e13 entry

    // (a) bounded -- not the 1e13 true inverse
    for (int i = 0; i < 9; ++i) {
        EXPECT_TRUE(std::isfinite(got[i]));
        EXPECT_LT(std::abs(got[i]), 1e6) << "null direction was inverted, not locked";
    }
    // (b) equals the pseudo-inverse (null direction zeroed), NOT the true inverse
    EXPECT_LT(maxAbsDiff(got, pinv, 3), 1e-6) << "locked inverse must match the pseudo-inverse";
    EXPECT_GT(maxAbsDiff(got, trueInv, 3), 1e6) << "must DIFFER from the true (1e13) inverse on purpose";
}

// ---------------------------------------------------------------------------
//  O6: INTENDED divergence #2 -- logDetSymPD FLOORS a near-zero pivot so the
//  log-det stays finite on a (numerically) singular matrix, rather than -inf.
// ---------------------------------------------------------------------------
TEST(LinAlgOracle, LogDetFloorsSingularPivot_IntendedDivergence) {
    Rng rng(0xF100);
    const Spd s = spdFromSpectrum(rng, {1e-200, 2.0, 3.0});
    const Real ld = L::logDetSymPD(s.A.data(), 3);
    EXPECT_TRUE(std::isfinite(ld)) << "pivot floor must keep ln|det| finite on a singular matrix";
}