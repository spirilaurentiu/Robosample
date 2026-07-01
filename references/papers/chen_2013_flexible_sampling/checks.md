# Checks / fixtures

Benchmark values from a study of commercial software (MOE, MacroModel). These are NOT reproducible unit-test fixtures for a from-scratch implementation; they are reference protocol settings and expected performance ranges. Use them as sanity targets for a low-mode sampler and as documentation of best-practice hyperparameters.

## Compound sets

| set | # compounds | rotatable bonds (opr_nrot) |
|---|---|---|
| Drug-like | 253 | 1-13 |
| Flexible (non-macrocyclic) | 50 | 12-17 (filter 12-20) |
| Macrocycle (ring >= 9 atoms) | 30 | 9-30 (10 compounds have >= 20) |

## Table 1 - default protocol settings

| software | method | force field | solvation | Duplicate RMS (A) | DE (kcal/mol) | Max-Iterations | RotSteps |
|---|---|---|---|---|---|---|---|
| MOE | LowModeMD | MMFF94x | Diel | 0.25 | 7 | 10,000 | NA |
| MOE | Stochastic Search | MMFF94x | Diel | 0.25 | 7 | 10,000 | NA |
| MacroModel | MT/LMOD | OPLS2005 | GB | 0.50 | 5 | 1000 | 100 |
| MacroModel | LMOD | OPLS2005 | GB | 0.50 | 5 | 1000 | 100 |
| MacroModel | LLMOD | OPLS2005 | GB | 0.50 | 5 | 1000 | 100 |
| MacroModel | MT/LLMOD | OPLS2005 | GB | 0.50 | 5 | 1000 | 100 |
| MacroModel | MD/LLMOD | OPLS2005 | GB | 0.75 | 10 | 5000 (LLMOD steps) + 5000 MD SA cycles | NA |

## Table 6 - best-performing "enhanced" protocols

| software | method | solvation | Duplicate RMS (A) | DE (kcal/mol) | Max-Iterations | RotSteps |
|---|---|---|---|---|---|---|
| MOE | LowModeMD | GB | 0.25 | 15 | 10,000 | NA |
| MacroModel | MT/LMOD | GB | 0.25 | 15 | 10,000 | 400 |
| MacroModel | MT/LLMOD | GB | 0.25 | 15 | 10,000 | 400 |

- MacroModel energy minimization: Polak-Ribiere conjugate gradient, convergence gradient <= 0.05 kJ/mol/A, 3000 minimization iterations.
- Global-minimum "found" criterion: conformer within 0.5 kcal/mol AND 0.5 A of the aggregated reference minimum.
- 3D pharmacophore distance bins: 2 A wide, [0,2) [2,4) ... [20, inf).

## Reproduction of bioactive structures (%BioConf_Rep, mean over 3 runs)

Given: DE, solvation, force field, RotSteps -> expect %BioConf_Rep at RMSD thresholds [0.5 / 1.0 / 1.5 / 2.0 A] and NbConfs.

### Drug-like set (Table 2)
- MOE LowModeMD, DE=7, Diel, MMFF94x: 38 / 77 / 87 / 94; NbConfs 156.
- MOE LowModeMD, DE=7, GB: 46 / 91 / 98 / 100; NbConfs 304.
- MOE LowModeMD, DE=15, GB (enhanced): 46 / 94 / 100 / 100; NbConfs 755.
- MOE Stochastic, DE=15, GB: 45 / 94 / 98 / 100; NbConfs 705.
- MacroModel MT/LMOD, DE=5, GB, OPLS2005 (default): 49 / 89 / 97 / 99; NbConfs 101.
- MacroModel MT/LMOD, DE=5, Diel: 41 / 73 / 81 / 87; NbConfs 55.
- MacroModel MT/LMOD, DE=15, GB, RotSteps=400 (enhanced): 65 / 97 / 100 / 100; NbConfs 1032.
- MacroModel LMOD, DE=5, GB (default): 43 / 78 / 87 / 95; NbConfs 75.
- MacroModel LLMOD, DE=5, GB (default): 39 / 67 / 79 / 90; NbConfs 41.
- MacroModel MT/LLMOD, DE=15, GB, RotSteps=400 (enhanced): 65 / 97 / 100 / 100; NbConfs 1112.

### Flexible set (Table 3)
- MOE LowModeMD, DE=7, Diel (default): 0 / 9 / 35 / 63; NbConfs 408.
- MOE LowModeMD, DE=7, GB: 5 / 47 / 86 / 95; NbConfs 1853.
- MOE LowModeMD, DE=15, GB (enhanced): 5 / 68 / 97 / 98; NbConfs 5448.
- MOE LowModeMD, DE=20, GB: 5 / 70 / 96 / 98; NbConfs 5712.
- MOE Stochastic, DE=7, Diel (default): 2 / 7 / 29 / 59; NbConfs 236.
- MOE Stochastic, DE=15, GB: 8 / 59 / 93 / 95; NbConfs 3986.
- MacroModel MT/LMOD, DE=5, GB (default): 1 / 20 / 59 / 89; NbConfs 205.
- MacroModel MT/LMOD, DE=15, GB, RotSteps=100: 5 / 45 / 87 / 98; NbConfs 1159.
- MacroModel MT/LMOD, DE=15, GB, RotSteps=200: 9 / 55 / 91 / 100; NbConfs 2281.
- MacroModel MT/LMOD, DE=15, GB, RotSteps=400 (enhanced): 7 / 65 / 95 / 100; NbConfs 4452.
- MacroModel LMOD, DE=5, GB (default): 1 / 7 / 28 / 63; NbConfs 143.
- MacroModel MT/LLMOD, DE=15, GB, RotSteps=400: 11 / 65 / 97 / 100; NbConfs 4950.
- MD/LLMOD, DE=10, GB (default): 0 / 50 / 78 / 94; NbConfs 1093.

### Macrocycle set (Table 4)
- MOE LowModeMD, DE=7, Diel (default): 29 / 49 / 68 / 80; NbConfs 61.
- MOE LowModeMD, DE=15, GB (enhanced): 34 / 72 / 84 / 89; NbConfs 1675.
- MOE Stochastic, DE=15, GB: 30 / 53 / 60 / 72; NbConfs 479.
- MacroModel MT/LMOD, DE=5, GB (default): 23 / 49 / 66 / 79; NbConfs 48.
- MacroModel MT/LMOD, DE=15, GB, RotSteps=100: 32 / 64 / 87 / 93; NbConfs 855.
- MacroModel MT/LMOD, DE=15, GB, RotSteps=400 (enhanced): 39 / 79 / 96 / 97; NbConfs 2415.
- MacroModel LMOD, DE=5, GB (default): 19 / 40 / 58 / 77; NbConfs 40.
- MacroModel MT/LLMOD, DE=15, GB, RotSteps=400: 40 / 77 / 92 / 94; NbConfs 2123.

### Macrocycle set, MD/LLMOD (Table 5)
- Default (DE=10, GB, DupRMS=0.75, OPLS2005, 5000 MD cycles, eigenvectors on Global_min, 5000 LLMOD steps): 37 / 77 / 93 / 97; NbConfs 297.
- Eigenvectors Initial_only: 27 / 67 / 90 / 97; NbConfs 269.
- Eigenvectors Every_step: 43 / 73 / 93 / 97; NbConfs 362.
- MD cycles = 0 (MD stage off): 37 / 60 / 83 / 90; NbConfs 176  (demonstrates MD stage is essential: 1.0 A drops 77 -> 60).
- Diel instead of GB: 30 / 63 / 87 / 93; NbConfs 364.
- DE=30, GB, DupRMS=0.25, 10,000 LLMOD steps: 43 / 77 / 100 / 100; NbConfs 2469.

## Global-minimum location (%GlobMin_found, three runs, enhanced protocols, GB) - Table 7

- Drug-like, LowModeMD (MMFF94x, DE=15, RotSteps NA): 98 / 98 / 97.
- Drug-like, MT/LMOD (OPLS2005, DE=15, RotSteps=400): 98 / 97 / 97.
- Drug-like, MT/LLMOD (RotSteps=400): 97 / 97 / 97.
- Drug-like, MT/LMOD default (DE=5, RotSteps=100): 90 / 91 / 92.
- Flexible, LowModeMD (DE=15): 70 / 66 / 66.
- Flexible, Stochastic (DE=15): 52 / 42 / 56.
- Flexible, MT/LMOD (DE=15, RotSteps=400): 64 / 64 / 74.
- Flexible, MT/LLMOD (RotSteps=400): 76 / 68 / 78.
- Macrocycle, LowModeMD (DE=15): 57 / 53 / 57.
- Macrocycle, MT/LMOD (DE=15, RotSteps=400): 70 / 67 / 60.
- Macrocycle, MT/LLMOD (RotSteps=400): 60 / 67 / 70.

## 3D pharmacophore coverage (aggregated over set)

- Flexible: enhanced LowModeMD 28,185; default LowModeMD 19,021; enhanced Stochastic 26,754; enhanced MT/LMOD ~28,000.
- Macrocycle: enhanced MT/LMOD 34,283; enhanced Stochastic 33,403; enhanced LowModeMD 33,098; default MD/LLMOD 29,524.

## Cost
- Enhanced MT/LMOD ~6x, enhanced LowModeMD ~7x slower than default counterparts.
- Doubling RotSteps ~doubles NbConfs and compute time.
