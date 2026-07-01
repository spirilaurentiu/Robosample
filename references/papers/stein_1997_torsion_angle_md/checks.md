# Checks / Fixtures - Stein, Rice & Brünger 1997 (Torsion-Angle MD)

These are benchmark numbers usable as regression fixtures. Robosample's relevant
overlap is the torsion-angle (internal-coordinate) MD kinematics (Eqs. 11-13) and
the restraint/energy functional forms (Eqs. 4-8); the success-rate/timing tables
are algorithm-level references, not unit tests.

## Protocol parameters (Table 2, protein; nucleic acid in parentheses)

| Stage | 1 (high-T TAMD) | 2 (slow-cool TAMD) | 3 (Cartesian MD) | 4 (minimization) |
|---|---|---|---|---|
| Temperature | 50,000 K (20,000 K) | 50,000 -> 1000 K (20,000->) | 1000 -> 300 K | - |
| Time step | 0.015 ps | 0.015 ps | 0.003 ps | - |
| Duration Delta t | 15 ps (60 ps) | 15 ps (60 ps) | 6 ps | - |
| w_NOE | 150 | 150 | 150 | 50 |
| w_dihedral | 100 (5) | 100 (5) | 100 | 300 |
| w_vdw | 0.1 | 0.1 -> 1.0 | 1.0 | 1.0 |

- Extended-strand init: atoms placed along x at 0.1 A intervals; y,z random in [0,1].
- Acceptance criterion: no NOE violation > 0.5 A; no dihedral violation > 5 deg;
  RMS bond deviation <= 0.02 A; RMS angle deviation <= 2.0 deg.
- NOE soft-asymptote switch point: R = d_upper + 0.5 A (E_NOE differentiable there).

## Success rate and computational efficiency (Table 6, HP-735)

| Method | metric | BPTI | Protein G | IL-8 | Villin 14T |
|---|---|---|---|---|---|
| TAMD | success rate | 98.0% | 87.7% | 89.3% | 84.7% |
| TAMD | struc. calc time | 921 s | 959 s | 2992 s | 1942 s |
| TAMD | comp. efficiency | 940 s | 1093 s | 3352 s | 2892 s |
| SA | success rate | 78.1% | 78.1% | 69.4% | 36.8% |
| SA | struc. calc time | 523 s | 563 s | 5597 s | 4612 s |
| SA | comp. efficiency | 670 s | 720 s | 8060 s | 12546 s |
| DGSA | success rate | 64.0% | 100.0% | 24.3% | 32.5% |
| DGSA | struc. calc time | 557 s | 558 s | 1631 s | 1351 s |
| DGSA | comp. efficiency | 870 s | 558 s | 6690 s | 4164 s |

## DNA dodecamer (CGCGAATTCGCG) results (Table 7)

| quantity | TAMD | Cartesian (SA) | DGSA | Original | Re-refined original |
|---|---|---|---|---|---|
| Success rate | 52.0% | 0.0% | 0.0% | - | - |
| <N_NOE> (>0.5 A outside bounds) | 0.0 | 0.9 | 19.5 | 0.0 | 0.0 |
| <Delta_NOE> | 0.050 A | 0.089 A | 0.26 A | 0.054 A | 0.051 A |
| <Delta_dihedral> | 0.45 deg | 11.89 deg | 14.52 deg | 1.68 deg | 0.54 deg |
| <Delta_bonds> | 0.012 A | 0.019 A | 0.021 A | 0.0086 A | 0.012 A |
| <Delta_angles> | 0.99 deg | 1.80 deg | 1.82 deg | 5.29 deg | 1.00 deg |
| <E_vdw> (LJ, Eq.4) | -269.2 kcal/mol | 7820 kcal/mol | 458.5 kcal/mol | -352.3 kcal/mol | -253.4 kcal/mol |
| RMS dev from original | 2.66 +/- 0.27 A | 5.23 +/- 0.61 A | 6.23 +/- 0.56 A | 0 | 2.26 +/- 0.04 A |

- TAMD ensemble atomic RMS to original 2.67 A; re-refined original moves to 1.25 A
  from the TAMD structure.

## Geometry and energy statistics (Table 5, selected; TAMD column)

| system | <Delta_NOE> | <Delta_bonds> | <Delta_angles> | <E_vdw> (LJ) |
|---|---|---|---|---|
| BPTI | 0.047 +/- 0.0025 A | 0.0034 +/- 0.00015 A | 0.67 +/- 0.013 deg | -50.3 kcal/mol |
| Protein G | 0.017 +/- 0.00075 A | 0.0018 +/- 0.00003 A | 0.48 +/- 0.036 deg | -103.5 kcal/mol |
| Interleukin-8 | 0.025 +/- 0.0026 A | 0.0025 +/- 0.00036 A | 0.54 +/- 0.012 deg | -112.9 kcal/mol |
| Villin 14T | 0.025 +/- 0.0047 A | 0.0027 +/- 0.00044 A | 0.48 +/- 0.033 deg | -75.19 kcal/mol |

## RMS from average structure (Table 3, backbone O,C,Ca,N)

| method | BPTI | Protein G | IL-8 (three ensembles) | Villin 14T |
|---|---|---|---|---|
| TAMD | 0.34 +/- 0.07 A | 0.26 +/- 0.05 A | 1.68/1.81/1.72 A | 1.38 +/- 0.27 A |
| SA | 0.43 +/- 0.15 A | 0.35 +/- 0.05 A | - | 1.42 +/- 0.36 A |
| DGSA | 0.43 +/- 0.14 A | 0.30 +/- 0.04 A | - | 1.54 +/- 0.17 A |

## System sizes (Table 1)

| | BPTI | Protein G | IL-8 | Villin 14T | DNA |
|---|---|---|---|---|---|
| Residues | 58 | 56 | 144 (dimer) | 126 | 24 |
| Total NOE distance restraints | 712 | 789 | 1720 | 1320 | 228 |
| Total H-bond restraints | 0 | 68 | 124 | 86 | 30 |
| Total dihedral restraints | 0 | 105 | 362 | 120 | 136 |

## Kinematics sanity checks (Eqs. 11-13)

- Rigid bond: |h_ij| constant => body j has exactly one rotational DOF relative to
  body i (about h_hat_ij). A unit test can verify omega_j - omega_i is parallel to
  h_hat_ij for any q̇_ij (Eq. 11).
- Velocity recurrence (Eq. 13): setting q̇_ij = 0 and omega_i = 0 must give
  r_dot_j = r_dot_i (pure translation of a rigid pair).
