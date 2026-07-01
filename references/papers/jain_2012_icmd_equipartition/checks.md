# Checks / Fixtures — Jain 2012, Equipartition Principle for ICMD

This is a theory paper with no benchmark result tables. The testable content is
a set of exact statistical/algebraic identities an implementation must satisfy.
Treat these as regression fixtures (all `*` = transpose; kT = thermal energy).

## Statistical (ensemble) fixtures

- **Cartesian equipartition (eq:7):** given a canonical ensemble at temperature
  T, expect `⟨p_i² / 2m_i⟩ = kT/2` for every Cartesian momentum coordinate i.
- **Modal equipartition (eq:25):** given modal coordinates ν sampled per the
  procedure, expect `⟨ν_i ν_j⟩ = kT δ_{ij}` — i.e. off-diagonal covariances
  zero, and `⟨ν_i²⟩ = kT`, so mean modal kinetic energy `⟨ν_i²/2⟩ = kT/2`.
- **Matrix form (eq:26):** expect the modal covariance `⟨ν ν*⟩ = kT · I_N`.
- **NO physical-coordinate equipartition:** the physical ICMD analog is
  `⟨θ̇ p*⟩ = kT · I_N` (eq:14), but `⟨θ̇_i θ̇_j⟩` and `⟨p_i p_j⟩` are NOT
  diagonal in general (mass matrix is dense) — a correct implementation must
  reproduce the dense cross-covariance, not a diagonal one.

## Algebraic-identity fixtures (unit tests on the operators)

- **Factorization (eq:15/16):** for the computed factor m(θ), expect
  `m(θ) m*(θ) = M(θ)` and, with l = m^{-1}, `l*(θ) l(θ) = M^{-1}(θ)`.
- **Round-trip transform (eq:27/33):** for any θ̇, expect
  `l*(θ) · (m*(θ) θ̇) = θ̇` (physical → modal → physical is identity);
  equivalently Table-1-left followed by Table-1-right returns the input to
  machine precision.
- **Kinetic-energy invariance (eq:8/18):** expect
  `½ ν*ν = ½ θ̇* M(θ) θ̇ = ½ p* M^{-1}(θ) p` for the same state
  (ν = m*θ̇, p = M θ̇). The modal norm equals the physical kinetic energy.
- **Mass matrix from Jacobian (eq:29):** expect `M(θ) = J*(θ) M J(θ)` where M is
  the diagonal Cartesian atom-mass matrix.
- **BAT factor (eq:30):** for a full BAT model (N=3n, J square invertible),
  expect the operator factor m(θ) to equal `J*(θ) M^{1/2}` (up to an orthogonal
  factor; both satisfy m m* = M).
- **Constrained mapping validity (eq:38/39):** with `P = M^{-1} H φ M^{1/2}`,
  expect `m* P P* m = I_N`, hence `⟨ν ν*⟩ = kT I_N` after θ̇ = P ẋ.
- **BAT initialization correctness (eq:36):** initializing Cartesian velocities
  ẋ with `⟨M^{1/2} ẋ ẋ* M^{1/2}⟩ = kT I_{3n}` then setting `θ̇ = J^{-1} ẋ`
  must yield `⟨ν ν*⟩ = kT I_{3n}`.

## Complexity fixtures

- Naive evaluation of the eq:33 transforms via explicit m/l matrices: **O(N²)**.
- Table 1 recursive base-to-tip algorithm (per transform): **O(N)** — linear in
  the number of DOF, no explicit formation of M, m, or l.

## Recommended velocity-initialization procedure (constrained & unconstrained)

1. Draw each modal velocity ν_i i.i.d. Gaussian with variance kT (mean 0), so
   `⟨ν_i²⟩ = kT` (eq:26 / equipartition).
2. Convert ν → θ̇ via the O(N) recursion (Table 1, right column / eq:33 right).
   Expect the resulting θ̇ to reproduce, in ensemble, `⟨ν ν*⟩ = kT I_N` and the
   correct dense `⟨θ̇ θ̇*⟩ = kT · M^{-1}(θ)` covariance.
This avoids the mass-matrix inverse required by the eq:39 Cartesian route.
