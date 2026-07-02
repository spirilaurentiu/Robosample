#pragma once
// ============================================================================
//  robo_linalg.hpp -- small dense linear-algebra kernels for the articulated-
//  body dynamics (the dof x dof hinge block D = ~H P H, dof <= 6).
//
//  TEST-LOCAL REFERENCE IMPLEMENTATION -- NOT the single source of truth and
//  NOT bound to the same symbols as the engine. `src/RobotEngine.cpp` does
//  NOT #include this header; it keeps its own file-local copies of these
//  kernels in its anonymous namespace (jacobiSymEig, invertDense, symSqrt,
//  symSqrtInv, and -- since docs/specs/singular-dof-fixman.md -- pseudoLogDet
//  in place of logDetSymPD). This header is a separate, hand-maintained copy
//  used ONLY by the test TUs as an independent oracle; despite an earlier
//  banner here claiming otherwise, it was never deduplicated with the engine
//  and no such dedup copy was ever deleted.
//
//  CONVENTION DRIFT (deliberate, tracked): `invertDense`/`symSqrt`/
//  `symSqrtInv` below already match the engine's current CC1/CC4 null-space
//  lock (a near-null eigendirection of D is treated as exactly null, locked
//  to 0, not regularized). `logDetSymPD` below does NOT -- it still floors a
//  near-zero Cholesky pivot at 1e-300 (the PRE-FIX convention `calcLogDetM`
//  used before it switched to `pseudoLogDet`'s null-space lock). This is kept
//  ON PURPOSE as the "old convention" reference for tests that need to
//  reconstruct what the pre-fix engine would have returned (e.g.
//  TestMassMatrix.cpp's LogDetMExcludesStructuralPhantomNullDirection) and as
//  the oracle for `Constraints::solveSmallSpd`/`calcConstraintLogDet`
//  (Constraints.cpp), which the spec's S2 explicitly did NOT fold into the
//  shared hinge-inertia lock (`G M^-1 G^T` is a dimensionally distinct
//  quantity). Do not "fix" `logDetSymPD` here to match `pseudoLogDet` without
//  re-deriving which callers depend on which convention (see TestLinearAlgebraOracle.cpp
//  O6 and TestConstraints.cpp/TestConstraintSolver.cpp).
//
//  Consumers (test TUs only): invertDense -> D^-1 (=> M^-1 seed, udot);
//  symSqrt / symSqrtInv -> sqrt(D^-1)/sqrt(D) (the momentum draw); logDetSymPD
//  -> ln det(D) (the OLD-convention Fixman tree term / constraint log-det).
//  jacobiSymEig underlies all of them.
// ============================================================================

#include <algorithm>
#include <cmath>

#include "robot_math.hpp"

namespace robo_linalg {

using robo::Real;

inline robo::Real logDetSymPD(const robo::Real* A, int n) {
    if (n <= 0) {
        return robo::Real(0);
    }
    robo::Real L[36]; // n <= 6
    for (int i = 0; i < n * n; ++i) {
        L[i] = robo::Real(0);
    }
    robo::Real logDet = robo::Real(0);
    for (int j = 0; j < n; ++j) {
        robo::Real sum = A[j * n + j];
        for (int k = 0; k < j; ++k) {
            sum -= L[j * n + k] * L[j * n + k];
        }
        // Guard against a non-finite / non-PD block (near-singular hinge): clamp
        // the pivot so the caller gets a large-but-finite Fixman rather than NaN,
        // which the Metropolis step can then reject cleanly.
        if (!(sum > robo::Real(1e-300))) {
            sum = robo::Real(1e-300);
        }
        const robo::Real Ljj = std::sqrt(sum);
        L[j * n + j] = Ljj;
        logDet += robo::Real(2) * std::log(Ljj);
        for (int i = j + 1; i < n; ++i) {
            robo::Real s = A[i * n + j];
            for (int k = 0; k < j; ++k) {
                s -= L[i * n + k] * L[j * n + k];
            }
            L[i * n + j] = s / Ljj;
        }
    }
    return logDet;
}

// Symmetric eigensolver (cyclic Jacobi) for a small dense symmetric matrix
// (n <= 6 here). Row-major input A[n*n]; outputs eigenvalues d[n] and
// eigenvectors as COLUMNS of V (V[i*n+k] = component i of eigenvector k).
// Self-contained: no SimTK::Eigen / LAPACK dependency (that symbol is not in
// the link line and produced an undefined-symbol ImportError).
inline void jacobiSymEig(const Real* Ain, int n, Real* d, Real* V) {
    Real a[36];
    for (int i = 0; i < n * n; ++i) {
        a[i] = Ain[i];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            V[i * n + j] = (i == j) ? Real(1) : Real(0);
        }
    }
    for (int sweep = 0; sweep < 100; ++sweep) {
        Real off = 0;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                off += a[p * n + q] * a[p * n + q];
            }
        }
        if (off <= Real(1e-30)) {
            break;
        }
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const Real apq = a[p * n + q];
                if (apq == Real(0)) {
                    continue;
                }
                const Real app = a[p * n + p];
                const Real aqq = a[q * n + q];
                const Real theta = (aqq - app) / (2 * apq);
                const Real t =
                    (theta >= 0 ? Real(1) : Real(-1)) / (std::abs(theta) + std::sqrt(theta * theta + 1));
                const Real cs = Real(1) / std::sqrt(t * t + 1);
                const Real sn = t * cs;
                // A <- J^T A J  (rotate columns then rows of the (p,q) plane).
                for (int i = 0; i < n; ++i) {
                    const Real aip = a[i * n + p];
                    const Real aiq = a[i * n + q];
                    a[i * n + p] = cs * aip - sn * aiq;
                    a[i * n + q] = sn * aip + cs * aiq;
                }
                for (int i = 0; i < n; ++i) {
                    const Real api = a[p * n + i];
                    const Real aqi = a[q * n + i];
                    a[p * n + i] = cs * api - sn * aqi;
                    a[q * n + i] = sn * api + cs * aqi;
                }
                // V <- V J  (accumulate eigenvectors as columns).
                for (int i = 0; i < n; ++i) {
                    const Real vip = V[i * n + p];
                    const Real viq = V[i * n + q];
                    V[i * n + p] = cs * vip - sn * viq;
                    V[i * n + q] = sn * vip + cs * viq;
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        d[i] = a[i * n + i];
    }
}

// Symmetric (Moore-Penrose) pseudo-inverse of the dof x dof articulated hinge
// matrix D = ~H P H, row-major in/out. D is symmetric and, for a well-posed
// hinge, SPD; this routine returns its inverse.
//
// CRITICAL ROBUSTNESS PROPERTY: a *degenerate* hinge -- a DOF with ~zero
// articulated inertia -- must be LOCKED, not regularized into a huge finite
// inverse. The canonical case is a single atom sitting exactly on its own
// torsion axis (D ~ 1e-33). It arises naturally once a ring-closing bond is
// removed and a ring atom (e.g. a proline/ring nitrogen) is left as a lone
// body whose root atom lies on the bond it rotates about; rotating a point
// mass on its own axis is unobservable, so the DOF has no dynamics.
//
// The previous code clamped such a D to +-1e-14 and returned 1/eps ~ 1e14.
// That inverse is enormous but FINITE, so it slips past the non-finite guard,
// multiplies the (tiny, noisy) force residual into a ~1e14 acceleration, and
// blows up on the next step -- exactly the observed failure. The correct
// action is a pseudo-inverse: invert the well-conditioned directions and set
// the inverse of any ~null direction to 0. DI=0 => G=0 => PPlus=P (full
// inertia passes through, i.e. the joint behaves as a weld) => udot=0 and the
// HMC seed velocity is 0 too. The DOF is frozen at its initial value, which is
// the physically exact treatment of a null coordinate.
inline void invertDense(const Real* A, int n, Real* Ainv) {
    const Real lockTol = Real(1e-12); // |eigenvalue| <= this is treated as null

    // Dominant case: 1-DOF Torsion. Closed form, with a lock instead of a clamp.
    if (n == 1) {
        const Real d = A[0];
        Ainv[0] = (std::abs(d) > lockTol) ? (Real(1) / d) : Real(0);
        return;
    }

    // n in {3,6}: symmetric pseudo-inverse via the in-house Jacobi eigensolver
    // (self-contained; no LAPACK/Eigen link dependency). For a non-degenerate
    // hinge every eigenvalue is well above tol, so this equals the true inverse;
    // any null eigendirection contributes 0 (that direction is locked).
    Real d[6];
    Real V[36];
    jacobiSymEig(A, n, d, V);
    Real scale = 0;
    for (int k = 0; k < n; ++k) {
        scale = std::max(scale, std::abs(d[k]));
    }
    const Real tol = std::max(lockTol, scale * Real(1e-12));
    for (int i = 0; i < n * n; ++i) {
        Ainv[i] = Real(0);
    }
    for (int k = 0; k < n; ++k) {
        if (std::abs(d[k]) <= tol) {
            continue; // null direction -> pseudo-inverse contributes 0 (locked)
        }
        const Real inv = Real(1) / d[k];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Ainv[i * n + j] += inv * V[i * n + k] * V[j * n + k];
            }
        }
    }
}

// dof==1 is a scalar sqrt; otherwise S = V diag(sqrt(lambda)) ~V, exactly as
// Simbody's multiplyBySqrtMInvPassOutward does (sqrtDI = U diag(sqrt) ~U).
inline void symSqrt(const Real* A, int n, Real* S) {
    if (n == 1) {
        S[0] = std::sqrt(std::max(Real(0), A[0]));
        return;
    }
    Real d[6];
    Real V[36];
    jacobiSymEig(A, n, d, V);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Real acc = 0;
            for (int k = 0; k < n; ++k) {
                acc += V[i * n + k] * std::sqrt(std::max(Real(0), d[k])) * V[j * n + k];
            }
            S[i * n + j] = acc;
        }
    }
}

// Inverse symmetric square root of a dof x dof SPD block: S = V diag(1/sqrt(lambda)) ~V.
// This is the EXACT inverse of symSqrt (same eigenvectors, reciprocal sqrt eigenvalues),
// so symSqrtInv(A) * symSqrt(A) == I. Used by multiplyBySqrtM (the forward sqrt sweep).
// Verified numerically to ~2e-13 over random SPD blocks (dof 1..6).
inline void symSqrtInv(const Real* A, int n, Real* S) {
    if (n == 1) {
        S[0] = Real(1) / std::sqrt(std::max(Real(1e-300), A[0]));
        return;
    }
    Real d[6];
    Real V[36];
    jacobiSymEig(A, n, d, V);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Real acc = 0;
            for (int k = 0; k < n; ++k) {
                acc += V[i * n + k] * (Real(1) / std::sqrt(std::max(Real(1e-300), d[k]))) * V[j * n + k];
            }
            S[i * n + j] = acc;
        }
    }
}

} // namespace robo_linalg