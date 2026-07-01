# Checks / fixtures - Spiridon & Minh 2017

## Idealized 4-bead serial chain C4 (Sec 3.1)
- Model: 3 bonds fixed at 1.54 Å; angles fixed at 90°; harmonic restraints of 83.66 kcal/mol on bonds and 43.46 kcal/mol on angles.
- No nonbonded or torsional force-field terms.
- Expected result: **uniform distribution in torsion space** (analytic).
- Fixture: given CDHMC **with** Fixman potential (or MIXED), expect a uniform torsion-angle histogram. Given CDHMC **without** Fixman potential, expect a periodic distortion (WRONG target) - confirms Fixman potential is required. Fixman torque is NOT required for correctness.
- Sim setup: 10 independent runs, 10^6 samples each; sample after each accept/reject step.

## Butane - Shirts' test (Sec 3.2, Table 1)
- Force field: GAFF (AmberTools 12), AM1-BCC charges. Temperatures 300 K and 450 K.
- 10 independent runs, 2×10^5 accept/reject steps each per temperature.
- Expected slope of log-histogram-ratio: `-(β_2 - β_1) = 0.13363`.
- Fitted slopes m (weights = inverse empirical variance):

| Type | DOF | m | σ | |m - 0.13363| |
|---|---|---|---|---|
| CDHMC | Torsional | 0.13693 | 3.86e-3 | 3.30e-3 |
| CDHMC | Rigid Body | 0.13524 | 5.03e-3 | 1.61e-3 |
| MIXED | Torsional | 0.13357 | 1.61e-3 | 0.06e-3 |
| MIXED | Rigid Body | 0.12259 | 4.03e-2 | 1.10e-2 |

- Fixture: all fitted slopes within ~1σ of 0.13363 -> passes Shirts' test.

## Alanine dipeptide free energy (Sec 3.3)
- Vacuum, AMBER ff12SB, 300 K, 7 independent runs, 10^7 steps each.
- UDHMC: 100 MD steps/trajectory, 1.5 fs. MIXED: 10 CDHMC : 1 UDHMC moves.
- Free energies relative to C7eq (kJ/mol), UDHMC vs MIXED (mean ± std over 7 runs):
  - C5: 0.18 ± 0.13 vs 0.27 ± 0.27
  - PPII: 2.15 ± 0.15 vs 2.29 ± 0.30
  - α_L: 8.39 ± 0.26 vs 8.29 ± 0.57
- Fixture: UDHMC and MIXED estimates agree within 1 std -> same distribution.

## Conformational region boundaries (Table 2, backbone dihedrals, degrees)

| region | φmin | φmax | ψmin | ψmax |
|---|---|---|---|---|
| C5 | -180 | -95 | 105 | 180 |
| PPII | -96 | -45 | 105 | 180 |
| C7eq | -96 | -45 | -25 | 104 |
| αL | 35 | 85 | -180 | 25 |

<!-- CHECK: Table 2 header in raw reads "φmin φmax ψmax ψmax" (duplicated ψmax); third column is ψmin. -->

## Macrocycle acceptance rates (Sec 3.5, Fig 5)
- Systems (PDB CCD): 1R6, AA0, AB0, ACZ, ADN. GAFF, AM1-BCC. 10 MD steps/move.
- UDHMC: ~35% acceptance at 2 fs; **exactly 0%** at 3 fs.
- CDHMC (torsional): monotonically decreasing with step size, still >0.1 at 10 fs; all systems >0.7 acceptance at 4 fs.
- 3 independent runs, 2×10^5 trials each.

## MIXED sampling efficiency (macrocycles, Sec 3.5)
- Runs at 625 K, 13×10^6 MD steps. UDHMC: snapshot every 20 trials (65,000 snapshots). MIXED: snapshot every 15 CDHMC + 5 UDHMC trials = every 650 MD steps (200,000 snapshots).
- Hierarchical clustering: weighted linkage, distance cutoff 0.25 Å (SciPy 0.19.0).
- Fixture: MIXED clusters ≈ UNION, UDHMC clusters ≈ INTERSECTION -> MIXED conformations nearly a superset of UDHMC.

## Misc parameters / constants
- Constraint tolerance η = 10^-4.
- Fixman torque omission gives ~25% computational speedup (ref 9).
- Empirical transition matrices (Table 3) and MFPTs (Table 4): MIXED has larger off-diagonal / shorter MFPT than UDHMC (per equal integrator-step count) -> MIXED transitions faster. Transitions counted every 2 steps (UDHMC) or 11 steps (MIXED) for equal MD-step count.
