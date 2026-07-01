# Checks / fixtures - Forrest & Suter 1994

## Force-field constants (exact, testable directly)

Ryckaert-Bellemans torsional potential `V_φ(φ) = C Σ_{n=0}^5 a_n cos^n φ`:
- given `C = 9.0 kJ/mol`, `a = [1, 1.31, -1.414, -0.3297, 2.828, -3.3943]`
- expect `V_φ(0) = C·(a0+a1+a2+a3+a4+a5) = 9.0 × (1 + 1.31 - 1.414 - 0.3297 + 2.828 - 3.3943) = 9.0 × 0.0000 = 0.0 kJ/mol`
  (coefficients sum to ≈0, so trans minimum `φ=0` has `V_φ≈0`). CHECK arithmetic: sum = 0.0000.
- given `φ = π` (cis): `cos π = -1`, `V_φ(π) = C·(a0 - a1 + a2 - a3 + a4 - a5) = 9.0 × (1 - 1.31 - 1.414 + 0.3297 + 2.828 + 3.3943) = 9.0 × 4.828 = 43.45 kJ/mol`.

Lennard-Jones `V_LJ(r) = 4ε[(σ/r)^12 - (σ/r)^6]`:
- given `ε = 410 J/mol`, `σ = 3.94 Å`
- expect minimum at `r_min = 2^{1/6} σ = 4.423 Å` with depth `V_LJ(r_min) = -ε = -410 J/mol`.
- expect `V_LJ(σ) = 0` at `r = 3.94 Å`.

## System / state point

- 20 polybead chains, each 24 CH2 units (C24), cubic box length 25 Å, PBC.
- Temperature `T = 480 K`.
- Density ≈ 0.68 g/cm3 (≈ polyethylene at 480 K, 1 atm).
- Fixed geometry: bond length `lb = 1.53 Å`, bond angle `θ = 112°`.
- Coordinates per chain: `Nb + 3 = 27` (3 Cartesian origin + 3 Euler + 21 torsion)
  vs `3 Nb = 72` for full Cartesian. DOF reduction factor ≈ 2.67.
- CH2 group mass = 14 g/mol (or `mb = 1` in reduced units).

## Moment of inertia

- Outermost torsional angle `φ24` rotates only the last bead, giving a constant
  instantaneous moment of inertia `I_24 = mb lb^2 sin^2 θ`.
  Given `mb=1`, `lb=1.53 Å`, `θ=112°`: `I_24 = 1 × 1.53^2 × sin^2(112°) = 2.3409 × 0.8597 = 2.012` (reduced units).
- Ordering: instantaneous moments of inertia descend from Euler angles (largest,
  rotate whole molecule) through innermost torsions to `φ24` (smallest, constant).

## Time-scale mapping (with Ĩ_k = ⟨I_k⟩)

- Each MC step ≈ 150 fs of MD trajectory; each MD step ≈ 15 fs.
- Mean effective time-step over 24 angular DOF: 14.9 fs.
- Mean effective time-step from 3 chain-origin (translational) DOF: 14.7 fs.
  (Agreement of 14.9 vs 14.7 fs validates the unique-time-scale claim.)
- With the alternative `Ĩ_k = I_24` for all k, angular DOF move on different time
  scales and translations match only `φ24` (slower approach to equilibrium).

## Sampling-efficiency optima

- Fixed `δt̃_MD = 1.5×10^-3`, vary `N_MD`: optimal trajectory length `N_MD ≈ 75`
  (flat minimum) for both `τ_s` and `τ_0.75`.
- `N_MD = 10` at its optimal step-size gives mean acceptance ≈ 65%.
- Optimal performance for any N_MD has `0.3 < ⟨Pacc⟩ < 0.7`.
- Larger `N_MD ≥ 100` gives more efficient sampling at lower acceptance.
- Contrast: Lennard-Jones fluid optimal `N_MD ≈ 10` (much shorter than polymer).
- Radius of gyration used as diffusion yardstick: `⟨s^2⟩ ≈ 40 Å^2`.

## Correctness identities (regression checks for an implementation)

- Standard HMC fluctuation identity: `⟨exp(-ΔH̃/kB T)⟩ = 1` (plain HMC).
- GC-HMC identity (eq 18): `⟨ Πc · exp(-ΔH̃/kB T) ⟩ = 1`, where
  `Πc = Π_{c=1}^{Nc} sin(φ2^(c)(t̃+Δt̃))/sin(φ2^(c)(t̃))`.
- Equivalent standard form with `ΔH̃* = ΔH̃ - kB T ln Πc`: `⟨exp(-ΔH̃*/kB T)⟩ = 1`.
- Mean acceptance (Gupta et al.): `⟨Pacc⟩ ≈ erfc( ½ ⟨ΔH̃*/kB T⟩^{1/2} )`.
- Leap-frog phase-space volume conservation: `[dq0][dπ̃0] = [dqn][dπ̃n]`.
