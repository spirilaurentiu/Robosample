# Checks — Wu & Brooks 2011, SGLD

## Sanity / limit fixtures (analytic)

- **LD limit:** given `λ = 0`, expect `λ_lf = λ_hf = 1`, `T_lf = T_hf = T`, `T̃ = T̃_0`,
  `χ_lf = 1`, `w_SGLD = 1` (constant), `T_sg = T`, and `Θ_SGLD = Θ_LD` (eqs 13, 24, 27).
- **Self-guiding temperature sign:** given `λ > 0` expect `T̃ > T̃_0` and `T_sg > T`; given `λ < 0`
  expect `T̃ < T̃_0` and `T_sg < T` (eq:27).
- **Energy-conservation factor:** given the definition of ξ (eq:12), expect total guiding-force
  power `Σ_i g_i · ṙ_i = 0`.
- **Frequency response of evolving average (eq:3):** for `q(t)=sin(2π ϖ t)`, given `2π ϖ t_L << 1`
  expect `q̃(t) ≈ q(t)` (amplitude -> 1); given `2π ϖ t_L >> 1` expect amplitude ∝ `1/ϖ` -> 0. In
  Fig 1: at `ϖ t_L = 0.1` average amplitude ≈ that of q; at `ϖ t_L = 10` amplitude ≈ 0.
- **Reweighted average (eq:26):** if `w_SGLD ≡ 1` (λ=0), expect `⟨A⟩ = ⟨A⟩_SGLD`.

## Skewed double well (single argon atom)

- Potential (kcal/mol): `ε_p = 500(x² + z²) + y²(y-2)² + 0.25 y`. Two minima along y near y=0 Å and
  y=2 Å.
- Params: T = 80 K, `t_L = 0.2 ps`, `δt = 1 fs`, run length 100 ns, collision frequency γ = 10/ps.
- Fixture: given `λ = 1`, expect `T_sg = 100.7 K` and ~10x more well-to-well transitions than
  `λ = 0` (`T_sg = 80 K`).
- Reweighting: given SGLD energy / y-coordinate distributions at various λ, after applying `w_SGLD`
  (eq:25/A9) expect convergence to the `λ = 0` canonical distribution (Figs 4b, 5b).

## Argon fluid

- 500 argon atoms, Lennard-Jones `ε = 119.8 K`, `σ = 3.405 Å`, cubic periodic box
  28.53 × 28.53 × 28.53 Å³.
- Params: `δt = 1 fs`, run length 10 ns, T = 100 K, collision frequency γ = 1/ps.
- IPS dispersion potential (r ≤ R):
  `ε_disp = -C_ij/r⁶ - (C_ij/R⁶)(1341/3064 + (77/141)(r/R)² + (61/141)(r/R)⁴ + (56/141)(r/R)⁸)`;
  0 for r > R.
- IPS repulsion potential (r ≤ R):
  `ε_rep = A_ij/r¹² + (A_ij/R¹²)(23/3620 + (8/151)(r/R)² + (66/151)(r/R)⁶ + (100/151)(r/R)¹⁰)`;
  0 for r > R.
- Fixture: reweighted energy distributions converge for `λ ≤ 1`; fail to converge for `λ > 1`
  (poor sampling / exponential weight blow-up).
- Qualitative: comparing LD at 100 K vs 140 K shows little distribution overlap, whereas SGLD stays
  near the 100 K distribution while increasing the diffusion constant.

## Alanine dipeptide

- CHARMM all-atom force field, distance-dependent dielectric `4r`, nonbonded cutoff 100 Å,
  `δt = 2 fs`, SHAKE on bond lengths, run length 200 ns, frames saved every 2 ps, `t_L = 0.2 ps`,
  T = 300 K, collision frequency γ = 10/ps.
- Dihedrals: φ = CT–N–Cα–C, ψ = N–Cα–C–NT. Transition counted when (φ,ψ) moves from within 40° of
  `(-90°, -70°)` to within 40° of `(-90°, 170°)`.
- Fixture: given SGLD guiding factors λ = 0.2, 0.5, 1.0 (T = 300 K), expect self-guiding
  temperatures `T_sg = 346 K, 458 K, 1067 K` respectively.
- Reweighting: φ–ψ distributions at λ = 0.7 and λ = 1 after applying `w_SGLD` recover the LD (λ=0)
  distribution (peak heights / baseline), noisier at larger λ.
