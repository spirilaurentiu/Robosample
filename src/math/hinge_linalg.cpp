#include "math/hinge_linalg.hpp"

namespace robo::detail {

robo::Real nullLockTol(robo::Real scale) {
    return std::max(kNullLockAbs, scale * kNullLockAbs);
}

void eigDecompAndTol(const robo::Real* A, int n, robo::Real* d, robo::Real* V, robo::Real& tol) {
    if (n == 1) {
        d[0] = A[0];
        V[0] = robo::Real(1);
        tol = kNullLockAbs;
        return;
    }
    jacobiSymEig(A, n, d, V);
    robo::Real scale = 0;
    for (int k = 0; k < n; ++k) {
        scale = std::max(scale, std::abs(d[k]));
    }
    tol = nullLockTol(scale);
}

robo::Real pseudoLogDet(const robo::Real* A, int n) {
    if (n <= 0) {
        return robo::Real(0);
    }
    robo::Real d[6];
    robo::Real V[36];
    robo::Real tol;
    eigDecompAndTol(A, n, d, V, tol);
    robo::Real logDet = robo::Real(0);
    for (int k = 0; k < n; ++k) {
        if (std::abs(d[k]) <= tol) {
            continue; // null direction: contributes 0 (pseudo-determinant), not ln(residue)
        }
        logDet += std::log(std::abs(d[k]));
    }
    return logDet;
}

int invertDense(const Real* A, int n, Real* Ainv) {
    // Dominant case: 1-DOF Torsion. Closed form, with a lock instead of a clamp.
    if (n == 1) {
        const Real d = A[0];
        const bool locked = std::abs(d) <= kNullLockAbs;
        Ainv[0] = locked ? Real(0) : (Real(1) / d);
        return locked ? 1 : 0;
    }

    // n in {3,6}: symmetric pseudo-inverse via the in-house Jacobi eigensolver
    // (self-contained; no LAPACK/Eigen link dependency). For a non-degenerate
    // hinge every eigenvalue is well above tol, so this equals the true inverse;
    // any null eigendirection contributes 0 (that direction is locked).
    Real d[6];
    Real V[36];
    Real tol;
    eigDecompAndTol(A, n, d, V, tol);
    for (int i = 0; i < n * n; ++i) {
        Ainv[i] = Real(0);
    }
    int numLocked = 0;
    for (int k = 0; k < n; ++k) {
        if (std::abs(d[k]) <= tol) {
            ++numLocked; // null direction -> pseudo-inverse contributes 0 (locked)
            continue;
        }
        const Real inv = Real(1) / d[k];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Ainv[i * n + j] += inv * V[i * n + k] * V[j * n + k];
            }
        }
    }
    return numLocked;
}

// Symmetric eigensolver (cyclic Jacobi) for a small dense symmetric matrix
// (n <= 6 here). Row-major input A[n*n]; outputs eigenvalues d[n] and
// eigenvectors as COLUMNS of V (V[i*n+k] = component i of eigenvector k).
// Self-contained: no SimTK::Eigen / LAPACK dependency (that symbol is not in
// the link line and produced an undefined-symbol ImportError).
void jacobiSymEig(const Real* Ain, int n, Real* d, Real* V) {
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

// Symmetric matrix square root of a dof x dof SPD block (for multiplyBySqrtMInv).
// dof==1 is a scalar sqrt; otherwise S = V diag(sqrt(lambda)) ~V, exactly as
// Simbody's multiplyBySqrtMInvPassOutward does (sqrtDI = U diag(sqrt) ~U).
void symSqrt(const Real* A, int n, Real* S) {
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

// Inverse symmetric square root of a dof x dof block A == DI (invertDense's
// pseudo-inverse of D): S = V diag(1/sqrt(lambda)) ~V. On every UNLOCKED
// direction this is the exact inverse of symSqrt(DI) (same eigenvectors,
// reciprocal sqrt eigenvalues), so symSqrtInv(DI)*symSqrt(DI) == I there --
// verified numerically to ~2e-13 over random SPD blocks (dof 1..6). Used by
// multiplyBySqrtM (the forward sqrt sweep, NMA route B).
//
// NULL-DIRECTION CONVENTION (S3, docs/specs/singular-dof-fixman.md): DI's
// eigenvalue is EXACTLY 0 on any direction invertDense locked (a structural
// phantom). The OLD code floored that 0 to 1e-300 and returned
// 1/sqrt(1e-300) ~ 3e149 -- an enormous but technically-finite regularizer
// that is NOT the pseudo-inverse convention every other null-space lock in
// this file uses, and is only harmless today because the caller's `in` on
// that direction is always bit-exact 0 (u_phantom is frozen by the SAME
// lock), so huge*0==0. That correctness depends on an invariant this
// function cannot see or enforce. Locking the direction to 0 HERE instead
// (same shared tol as invertDense/pseudoLogDet) makes symSqrtInv robust on
// its own terms -- multiplyBySqrtM's derivation only needs sqrt(D)==0 on a
// direction that mobility WILL NEVER be seeded on, and 0 is the correct
// "no information on this gauge direction" answer, not the largest float
// that still avoids inf.
void symSqrtInv(const Real* A, int n, Real* S) {
    if (n == 1) {
        const bool locked =
            A[0] <= kNullLockAbs; // DI >= 0 always; <=0 means locked (or non-PSD, caught downstream)
        S[0] = locked ? Real(0) : (Real(1) / std::sqrt(A[0]));
        return;
    }
    Real d[6];
    Real V[36];
    Real tol;
    eigDecompAndTol(A, n, d, V, tol);
    for (int i = 0; i < n * n; ++i) {
        S[i] = Real(0);
    }
    for (int k = 0; k < n; ++k) {
        if (std::abs(d[k]) <= tol) {
            continue; // null direction -> contributes 0 (locked; matches invertDense/pseudoLogDet)
        }
        const Real inv = Real(1) / std::sqrt(d[k]);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                S[i * n + j] += inv * V[i * n + k] * V[j * n + k];
            }
        }
    }
}

} // namespace robo::detail
