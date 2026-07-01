# Checks / fixtures - Larsen 2014 GNEIMO CASP refinement

This is an application/benchmark paper. The numbers below are protocol parameters
and reported refinement results. They are regression targets for a GNEIMO-REXMD
refinement port, not exact reproducible unit-test outputs (stochastic MD).

## Protocol parameters (deterministic config fixtures)

- given force field, expect **AMBER99SB**.
- given implicit solvent, expect **GB/SA OBC** with solute dielectric **1.5**,
  solvent dielectric **78.3**, nonpolar probe radius **1.4 Angstrom**.
- given nonbonded cutoff, expect switch-off at **20 Angstrom**.
- given thermostat, expect **Nose-Hoover** (constant temperature).
- given integrator, expect **Lobatto**, time step **5 fs** (stable up to 10 fs).
- given REXMD setup, expect **32 replicas**, temperature range **310-415 K**,
  Metropolis temperature swap attempted every **5 ps**.
- given per-target sim time, expect **15-100 ns per replica** (32 replicas).

## Metric definitions (fixtures for scoring code)

- given GDT_TS, expect distance-to-native cutoff set **{8, 4, 2, 1} Angstrom**,
  averaging aligned C-alpha counts over the four cutoffs.
- given native-contact definition, expect contact = C-alpha(i)-C-alpha(j) distance
  **< 8 Angstrom** AND residues **> 4 apart** in sequence; a snapshot pair counts as
  the same contact if within **0.5 Angstrom** of the native distance.

## Aggregate refinement results

- given 23 TR refinement targets, expect avg improvement **+4.9 GDT_TS**,
  **+0.04 TM-score**, **0.52 Angstrom RMSD**; refinement in **19/23** (GDT/TM),
  **21/23** (RMSD-based, per Conclusions). Max single-target GDT gain **14.0**,
  max TM gain **0.13**.
- given 7 T0 prediction targets, expect avg improvement **+4.5 GDT**, **+0.04 TM**,
  **0.7 Angstrom RMSD**; all within 4 Angstrom RMSD of crystal except **T0488**.
- given the whole 30-target set, expect refinement of up to **1.3 Angstrom** RMSD.

## Table I - refinement (TR) targets, C-alpha, {start / best GNEIMO / best CASP}

Format per target: GDT_TS | TM-Score | RMSD (Angstrom).

- TR429: 31.5 / 45.7 / 39.8 ; 0.46 / 0.59 / 0.53 ; 6.82 / 5.76 / 6.62
- TR435: 80.2 / 87.9 / 83.4 ; 0.86 / 0.91 / 0.89 ; 2.14 / 1.65 / 1.88
- TR453: 86.6 / 91.5 / 86.6 ; 0.87 / 0.92 / 0.88 ; 1.51 / 1.10 / 1.48
- TR454: 58.5 / 71.0 / 60.2 ; 0.79 / 0.87 / 0.81 ; 3.26 / 2.47 / 3.09
- TR461: 89.4 / 91.2 / 90.4 ; 0.93 / 0.94 / 0.94 ; 1.63 / 1.55 / 1.60
- TR462: 63.8 / 67.1 / 69.1 ; 0.80 / 0.81 / 0.83 ; 2.55 / 2.55 / 2.28
- TR464: 75.4 / 83.3 / 81.2 ; 0.76 / 0.81 / 0.82 ; 2.77 / 2.45 / 2.28
- TR469: 76.6 / 80.3 / 89.3 ; 0.74 / 0.79 / 0.85 ; 2.13 / 1.89 / 1.68
- TR476: 36.5 / 45.8 / 42.5 ; 0.42 / 0.50 / 0.47 ; 6.92 / 6.31 / 5.42
- TR488: 85.3 / 86.8 / 90.5 ; 0.88 / 0.89 / 0.92 ; 2.13 / 1.91 / 1.57
- TR517: 68.5 / 72.8 / 69.4 ; 0.77 / 0.80 / 0.78 ; 4.60 / 3.59 / 3.95
- TR530: 82.4 / 90.7 / 88.5 ; 0.84 / 0.90 / 0.88 ; 2.00 / 1.33 / 1.63
- TR557: 63.4 / 68.0 / 66.6 ; 0.73 / 0.76 / 0.78 ; 4.10 / 3.37 / 3.30
- TR568: 50.8 / 53.9 / 56.2 ; 0.55 / 0.57 / 0.60 ; 6.26 / 5.60 / 4.26
- TR569: 68.4 / 72.2 / 77.8 ; 0.71 / 0.73 / 0.81 ; 3.05 / 2.94 / 1.98
- TR574: 57.3 / 66.4 / 58.6 ; 0.64 / 0.72 / 0.65 ; 3.52 / 2.90 / 3.37
- TR576: 61.3 / 61.3 / 66.4 ; 0.72 / 0.72 / 0.76 ; 6.67 / 6.67 / 3.86
- TR592: 89.8 / 93.5 / 93.4 ; 0.92 / 0.94 / 0.95 ; 1.22 / 1.09 / 0.95
- TR594: 85.3 / 85.5 / 85.8 ; 0.90 / 0.91 / 0.91 ; 1.83 / 1.62 / 1.64
- TR606: 67.1 / 67.1 / 75.9 ; 0.73 / 0.73 / 0.81 ; 4.87 / 3.95 / 2.91
- TR614: 71.9 / 71.9 / 80.2 ; 0.76 / 0.76 / 0.84 ; 5.36 / 4.41 / 2.78
- TR622: 66.7 / 66.7 / 73.5 ; 0.74 / 0.74 / 0.78 ; 6.54 / 6.17 / 3.25
- TR624: 50.0 / 59.3 / 63.4 ; 0.49 / 0.58 / 0.63 ; 5.21 / 3.95 / 3.86
- avg:   68.1 / 73.0 / 73.4 ; 0.74 / 0.78 / 0.79 ; 3.79 / 3.27 / 2.85
- improvement (best GNEIMO / best CASP): GDT +4.9 / +5.3 ; TM +0.04 / +0.05 ;
  RMSD 0.52 / 0.93.

## Table II - structure-prediction (T0) targets, C-alpha, {start / best GNEIMO / best CASP}

Format per target: GDT_TS | TM-Score | RMSD (Angstrom).

- T0387: 85.7 / 89.3 / 95.5 ; 0.88 / 0.90 / 0.94 ; 1.95 / 1.43 / 1.01
- T0453: 80.3 / 83.0 / 87.1 ; 0.82 / 0.84 / 0.89 ; 3.75 / 1.85 / 1.47
- T0469: 82.0 / 86.5 / 73.4 ; 0.79 / 0.83 / 0.74 ; 1.93 / 1.88 / 2.47
- T0472: 89.5 / 89.7 / 61.8 ; 0.93 / 0.93 / 0.76 ; 1.21 / 1.06 / 2.68
- T0488: 71.8 / 75.3 / 86.0 ; 0.75 / 0.79 / 0.87 ; 4.60 / 3.54 / 1.97
- T0492: 82.3 / 92.7 / 85.8 ; 0.82 / 0.91 / 0.87 ; 1.67 / 1.16 / 1.70
- T0554: 67.1 / 73.7 / 32.3 ; 0.80 / 0.84 / 0.44 ; 3.40 / 2.65 / 8.31
