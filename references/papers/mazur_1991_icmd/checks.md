# Checks / fixtures — Mazur 1991, explicit internal-coordinate MD

Test system: oligopeptide (Ala)_n, α-helical start, single molecule in vacuum
(solvent neglected). Potentials: mixed ECEPP + CHARMM. Integrator: Beeman's method
with 4th-order velocity prediction. Linear system eq:18 solved each step by Cholesky
("Kholezki") decomposition. Heated to ~300 K, equilibrated 4 ps before production.
Trajectories: 110 steps, first 10 ignored, statistics over the remaining 100 steps.

## Model definitions

| model | description | atoms | internal DOF (N) | step time (s, EC1051) |
|---|---|---|---|---|
| 1 | all H explicit; all bond lengths fixed; phase angles fixed | 93 | 133 | 6.0 |
| 2 | model 1 + all valence angles fixed | 93 | 42 | 1.1 |
| 3 | model 2 + united-atom Ala side chains (no methyl torsions) + C-terminus hydroxyl torsion fixed | (united) | 32 | 3.5 |

Note: step time for model 1 is dominated by solving the 133-equation linear system;
for models 2 and 3 the energy-gradient (force) calculation dominates.

## Energy-conservation fixtures (time step h = 0.5 fs, over 100 steps)

| model | ⟨E⟩ (kcal/mol) | δ_E (relative RMS energy error) |
|---|---|---|
| 1 | 23.5  | ~0.4 × 10⁻⁶ |
| 2 | −22.5 | ~0.4 × 10⁻⁵ |
| 3 | −26.8 | ~0.8 × 10⁻⁶ |

<!-- CHECK: OCR of the δ_E exponents is garbled ("0.4 x lo-', 0.4 x 10 P5, 0.8 x 10mm6").
Interpreted as 0.4e-6, 0.4e-5, 0.8e-6 respectively; magnitudes uncertain, treat as
order-of-magnitude fixtures only. -->

## Time-step scaling fixtures

- As h increases, δ_E increases, with the ratios between the three models staying
  approximately constant.
- Given the common protein-MD step h = 1 fs, model 1 reaches δ_E ≈ 10⁻⁴.
  <!-- CHECK: OCR "lo--'" read as 1e-4; exact exponent uncertain. -->
- The SAME accuracy level (that model 1 reaches at 1 fs) is reached by:
  - model 2 at h ≈ 9 fs,
  - model 3 at h ≈ 13 fs.
  => freezing fast DOF allows ~9-13× larger time step for equal energy conservation.

## Qualitative correctness checks (asserted, not numeric)

- Forward/backward time-reversal integration consistent for each model.
- Conservation of total energy, linear momentum, and angular momentum verified.
- For the fully unfixed molecule, trajectories from eqs. 16 and 18 match a
  traditional Cartesian-coordinate MD simulation.

## Constants / relations for a port

- Temperature: T = 2⟨K⟩/(N k_B) with N = number of internal DOF (see table), not 3×atoms.
- Only varying parts of U are computed; intra-rigid-body interactions are excluded,
  which is why ⟨E⟩ differs so much between models.
