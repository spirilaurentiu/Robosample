#pragma once

/**
 * @file hinge_linalg.hpp
 * @brief Dense symmetric linear-algebra kernels for a single joint's @c n <= 6
 *        hinge-inertia block: pseudo-inverse, pseudo-log-determinant, and
 *        symmetric (inverse) square root, all sharing one null-space lock.
 *
 * These are the pure numeric primitives the articulated-body solver applies to
 * each body's dof x dof hinge inertia @f$ D_b = H_b^\top P_b H_b @f$. They own
 * no state, allocate nothing beyond fixed-size locals (capped at @c n == 6), and
 * carry no dependency on engine types or on SimTK/LAPACK. All matrices are
 * row-major and dense.
 *
 * The single correctness link across the four consumers is the shared
 * null-space lock (see @c kNullLockAbs): @c invertDense (forward dynamics),
 * @c pseudoLogDet (Fixman kinetic log-det, INV-5), and @c symSqrtInv (the
 * mass-metric square root operators, INV-5) all classify a direction as null on
 * the identical test, so a gauge direction removed by the dynamics contributes
 * factor 1 (0 to the log-det) and 0 to the metric square root -- never a
 * floating-point residue treated as physically present.
 */

#include <algorithm>
#include <cmath>

#include "robot_math.hpp"

namespace robo::detail {

/**
 * @brief Absolute floor of the shared null-space lock tolerance.
 *
 * A direction of a hinge-inertia block is classified as null (a structural or
 * gauge phantom, e.g. a lone single-atom Torsion body whose atom lies on its own
 * rotation axis) when its eigenvalue magnitude is at or below the lock returned
 * by @ref nullLockTol, whose value never drops below this constant. Flooring the
 * relative tolerance here keeps an all-tiny block locked rather than rescaled to
 * nothing. This constant is the sole source of truth for null-direction
 * classification across @c invertDense, @c pseudoLogDet, and @c symSqrtInv;
 * they SHALL agree on which directions are null so their results stay mutually
 * consistent (INV-5/INV-6, docs/specs/singular-dof-fixman.md CC1/CC4/N1).
 */
constexpr robo::Real kNullLockAbs = robo::Real(1e-12);

/**
 * @brief Null-space lock tolerance for a block whose largest eigenvalue
 *        magnitude is @p scale.
 * @param[in] scale Largest absolute eigenvalue of the block (its magnitude).
 * @return @c max(kNullLockAbs, scale * kNullLockAbs): a relative lock floored at
 *         @c kNullLockAbs.
 */
robo::Real nullLockTol(robo::Real scale);

/**
 * @brief Eigen-decompose a symmetric block and return the null-space lock
 *        tolerance for it, so every caller locks on identical @c (d, V, tol).
 * @param[in]  A   Row-major symmetric @p n x @p n block.
 * @param[in]  n   Block dimension (@c 1 <= n <= 6).
 * @param[out] d   Eigenvalues, @p n entries; unordered (see @ref jacobiSymEig).
 * @param[out] V   Eigenvectors as columns, @p n x @p n row-major
 *                 (@c V[i*n+k] is component @c i of eigenvector @c k).
 * @param[out] tol Null-space lock tolerance for this block.
 * @note @c n == 1 is closed form: the eigenvalue is the scalar itself, @c V[0]
 *       is 1, and @p tol is the bare @c kNullLockAbs.
 */
void eigDecompAndTol(const robo::Real* A, int n, robo::Real* d, robo::Real* V, robo::Real& tol);

/**
 * @brief Pseudo-log-determinant @f$ \sum \ln|\lambda_k| @f$ of a symmetric block,
 *        with null directions omitted.
 * @param[in] A Row-major symmetric @p n x @p n block.
 * @param[in] n Block dimension; @c 0 <= n <= 6.
 * @return The sum of @c ln of the absolute non-null eigenvalues. A null direction
 *         (per the shared lock) contributes 0, i.e. factor 1 to the determinant,
 *         never @c ln of a residue. Returns 0 for @c n == 0 (a Weld body:
 *         @c det == 1).
 * @pre @p A is symmetric and, for a well-posed hinge, positive definite. The
 *      magnitude is taken defensively against roundoff sign noise at the lock
 *      boundary; a genuinely negative eigenvalue is a non-PD input (a caller bug)
 *      surfaced downstream as NaN/Inf, not corrected here.
 * @note This is the null-lock counterpart of @c invertDense: a direction the
 *       dynamics pseudo-inverse removes contributes exactly 0 here, so
 *       @c calcLogDetM stays consistent with the forward dynamics (INV-5, CC1).
 * @see invertDense, symSqrtInv
 */
robo::Real pseudoLogDet(const robo::Real* A, int n);

/**
 * @brief Symmetric (Moore-Penrose) pseudo-inverse of a dof x dof hinge block
 *        @f$ D = H^\top P H @f$.
 * @param[in]  A    Row-major symmetric @p n x @p n block; SPD for a well-posed hinge.
 * @param[in]  n    Block dimension (@c 1 <= n <= 6).
 * @param[out] Ainv Row-major @p n x @p n pseudo-inverse.
 * @return Number of locked (null) directions: 0 for a well-conditioned @p A.
 *         The forward-dynamics caller uses a nonzero count to fire the CC1/CC4
 *         fail-loud gate on any body whose lock was not anticipated.
 * @post Every well-conditioned direction is inverted exactly (@c Ainv equals the
 *       true inverse when no direction locks); every null direction inverts to 0,
 *       so @c Ainv is the pseudo-inverse, not a large finite regularization.
 * @warning A degenerate DOF (near-zero articulated inertia, e.g. a point mass on
 *          its own torsion axis) SHALL be locked to 0, not clamped to a huge
 *          finite inverse. Zeroing that direction freezes the DOF at its initial
 *          value (a weld), the physically exact treatment of a null coordinate; a
 *          finite clamp would instead produce a ~1e14 acceleration and diverge.
 * @see pseudoLogDet, symSqrtInv
 */
int invertDense(const Real* A, int n, Real* Ainv);

/**
 * @brief Symmetric eigensolver (cyclic Jacobi) for a small dense symmetric matrix.
 * @param[in]  Ain Row-major symmetric @p n x @p n matrix (@c n <= 6).
 * @param[in]  n   Matrix dimension (@c n <= 6).
 * @param[out] d   Eigenvalues, @p n entries.
 * @param[out] V   Eigenvectors as columns, @p n x @p n row-major
 *                 (@c V[i*n+k] is component @c i of eigenvector @c k).
 * @post @c V is orthonormal and @f$ V \, \mathrm{diag}(d) \, V^\top = Ain @f$.
 * @warning The eigenvalues are returned in no particular order; a caller that
 *          needs them sorted (e.g. NMA mode ordering) SHALL sort them itself.
 */
void jacobiSymEig(const Real* Ain, int n, Real* d, Real* V);

/**
 * @brief Symmetric matrix square root @f$ S = V \, \mathrm{diag}(\sqrt{\lambda}) \, V^\top @f$
 *        of a dof x dof SPD block.
 * @param[in]  A Row-major symmetric positive-semidefinite @p n x @p n block.
 * @param[in]  n Block dimension (@c 1 <= n <= 6).
 * @param[out] S Row-major @p n x @p n symmetric square root, @c S*S == A.
 * @post Negative eigenvalues (roundoff noise) are floored to 0 before the square
 *       root, so @p S is real for any input; a materially negative eigenvalue is a
 *       non-PSD caller error, not corrected here.
 */
void symSqrt(const Real* A, int n, Real* S);

/**
 * @brief Inverse symmetric square root
 *        @f$ S = V \, \mathrm{diag}(1/\sqrt{\lambda}) \, V^\top @f$ of a dof x dof
 *        block, with null directions locked to 0.
 * @param[in]  A Row-major symmetric block, conventionally @c DI (the
 *               @c invertDense pseudo-inverse of a hinge block @c D).
 * @param[in]  n Block dimension (@c 1 <= n <= 6).
 * @param[out] S Row-major @p n x @p n inverse square root.
 * @post On every unlocked direction @p S is the exact inverse of @ref symSqrt of
 *       the same block (@c symSqrtInv(DI) * symSqrt(DI) == I there). On any
 *       direction locked by the shared null-space test @p S is 0.
 * @note The zero on a locked direction is the intended "no information on this
 *       gauge direction" value; it uses the same lock as @c invertDense and
 *       @c pseudoLogDet, so the mass-metric operators stay consistent with the
 *       forward dynamics (INV-5, docs/specs/singular-dof-fixman.md S3).
 * @see symSqrt, invertDense
 */
void symSqrtInv(const Real* A, int n, Real* S);

} // namespace robo::detail
