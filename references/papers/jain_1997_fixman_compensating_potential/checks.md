# Checks / fixtures — Jain 1997, Compensating Mass Matrix Potential

This is a derivation/algorithm paper; it states no benchmark energy tables. The
testable content is a set of exact structural identities that any port of the
compensating-potential algorithm must satisfy, plus consistency checks against a
brute-force reference implementation.

## Exact identities (regression assertions)

- **Determinant reduction.** For the innovations factorization,
  $\det\{I+H\phi K\}=1$ (block lower triangular with identity diagonal blocks),
  therefore
  `det(M) == det(D) == prod_i det(D(i))`  (within float tolerance).
  Test: form the dense mass matrix $\mathcal{M}=H\phi M\phi^*H^*$ and the
  block-diagonal $D$ from Eq. (3.6); assert `abs(logdet(M) - sum_i logdet(D(i))) < 1e-8 * |logdet(M)|`.

- **Compensating potential consistency.**
  `Vc == 0.5*ln det(M) == 0.5*sum_i ln det(D(i))`  (Eq. 2.9 vs Eq. 4.20).
  Given the same configuration $\boldsymbol{\theta}$, the direct log-det of the
  dense mass matrix and the sum-of-block-log-dets must agree.

- **Compensating torque consistency.** For every hinge $i$:
  `Tc(i)  (Eq. 4.18, -h(i)^T F[Q11+Q22])`
  `     == Trace{P(i) Y(i) H(i)}  (Eq. 4.17)`
  `     == 0.5 * Trace{ M^{-1} * dM/dtheta_i }  (Eq. 2.14, brute force)`.
  Reference `dM/dtheta_i` by finite-differencing the dense mass matrix; assert
  all three agree to FD tolerance.

- **Gradient-of-potential check.** $T_c = \nabla_\theta \mathcal{V}_c$: central
  finite difference of $\mathcal{V}_c(\theta)=\tfrac12\sum_i\ln\det D(i)$ w.r.t.
  each $\theta(i)$ must equal the analytic $T_c(i)$ from Eq. (4.18).

- **Base-cluster torque is zero.** The compensating torque on the 6 rigid-body
  degrees of freedom of the base cluster is identically zero:
  `Tc[base 6 dof] == 0`.

- **Axial map self-consistency.** For any 3-vector $v$: `F[skew(v)] == v`, and for
  any $A in R^{3x3}$: `Trace{A * skew(v)} == -v^T F[A]` (Eq. 4.18b).

## Structural / complexity facts

- Serial chain of $n$ single-dof rotational hinges plus a 6-dof base cluster has
  `N = n + 5` total degrees of freedom.
- Articulated body inertia blocks: `P(k) in R^{6x6}`, `M(k) in R^{6x6}`.
- For a single-dof rotational hinge, `D(k)` is a scalar (`1x1`);
  `ln det D(k) = ln D(k)`.
- The compensating-torque computation (Step 3: Eq. 4.16 + Eq. 4.18) adds an
  $O(\mathcal{N})$ pass; total algorithm cost stays $O(\mathcal{N})$
  (vs $O(\mathcal{N}^3)$ for the naive dense $\tfrac12\operatorname{Trace}\{M^{-1}M_{\theta_i}\}$).

## Qualitative validation targets (from cited prior work)

- Compensating potential is negligible for rigid **bond-stretch** constraints but
  significant for rigid **bond-angle** constraints — a port can sanity-check that
  $T_c$ is near zero when only bond lengths are frozen and non-negligible when
  angles are frozen.
- For **n-butane**, including the compensating potential brings the number of
  dihedral transitions in constrained MD into agreement with unconstrained
  Cartesian MD (qualitative reproduction target, not a numeric fixture).
