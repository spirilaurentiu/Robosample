# Checks / Fixtures - ICFF (Katritch, Totrov, Abagyan 2003)

Regression fixtures pulled from the paper. Source force field: MMFF94s. Test set: Halgren's diverse organic-molecule conformational set (refs 33, 66), excluding ring-conformation pairs.

## Global headline accuracy (Abstract / Conclusions)
- given ICFF projection of MMFF94s, expect equilibrium conformational energy-difference RMSD ~ **0.64 kcal**.
- given ICFF projection of MMFF94s, expect detailed torsion energy-profile RMSD ~ **0.37 kcal** (at 5-kcal MMFF94s cutoff).
- given rigid MMFF94s (no implicit flexibility), expect equilibrium energy RMSD > **1.2 kcal** and torsion-profile RMSD ~ **3.4-3.5 kcal**. (ICFF reduces 3.5 -> 0.37 kcal.)
- given full source MMFF94s vs QM, expect conformational energy RMSD ~ **0.3-0.5 kcal**.

## Fitting / algorithm parameters (fixtures for the parameterizer)
- Torsion grid scan: torsion angle from **0° to 360°**, increment **12°**, giving **30 conformations** per torsion fragment.
- Torsion fragment = the two bond-flanking atoms + all their immediate neighbors = **4 to 8 atoms** total.
- Fourier series order: **sixfold** (k = 0..6) reproduces profiles to **~0.01 kcal**. Fivefold -> ~0.08 kcal; fourfold -> ~0.2 kcal.
- Weighted least-squares fit (weight $[E-\min E+1]^{-1}$) improves fit RMSD to **~0.001 kcal** below the 5-kcal cutoff.
- Fourier coefficients with magnitude < **0.01** can be nullified without noticeable accuracy loss.
- Torsion-constraint restraint: $C_r = 10000$ kcal, keeps angle deviation from target **< 0.1°**.
- Soft-repulsion blend factor: optimal RMSD for **0.50 < C < 0.60** across four objective functions; consensus **C = 0.55**.
- ICFF is ~**50% faster** than ECEPP torsion force field (fewer pairwise interactions).
- Parameterization speed: ~**0.1 s** per average drug-like compound.

## Table 1 - Equilibrium energy comparisons vs MMFF94s Cartesian (and vs QM MP4SDQ/TZP)
Columns: ICFF | MMFF94rigid | ICFF(vs QM) | MMFF94s(vs QM) | MM3(vs QM).

| quantity | ICFF | MMFF94rigid | ICFF vs QM | MMFF94s vs QM | MM3 vs QM |
|---|---|---|---|---|---|
| number of comparisons | 120 | 118 | 120 | 121 | 108 |
| conformers with no local minimum | 1 | 3 | 1 | 0 | 13 |
| mean geometry RMSD (Å) | 0.04 | 0.06 | 0.04 | 0.04 | - |
| max geometry RMSD (Å) | 0.27 | 0.74 | 0.48 | 0.31 | - |
| number favoring wrong conformer | 7 | 8 | 6 | 6 | 18 |
| max energy deviation (kcal) | 4.01 | 8.12 | 4.2 | 1.56 | 3.47 |
| count in 1-2 kcal | 12 | 14 | 14 | 1 | 22 |
| count in 2-3 kcal | 0 | 11 | 2 | 0 | 5 |
| count in 3-5 kcal | 1 | 1 | 1 | 0 | 0 |
| energy RMSD (kcal) | **0.64** | **1.21** | **0.74** | **0.31** | **1.03** |
| energy RMSD when E_MMFF94s < 4 kcal | 0.46 | 0.92 | 0.57 | 0.30 | 0.91 |

- Consistency check: RMSD(ICFF vs QM)=0.74 ≈ quadrature of RMSD(ICFF vs MMFF94s)=0.64 and RMSD(MMFF94s vs QM)=0.31, i.e. $\sqrt{0.64^2+0.31^2}\approx0.71$ (independent error contributions).
- Added-error estimate: given Cartesian FF RMSD 1.03 and ICFF contribution 0.64, ICFF projection adds ~18%: $\sqrt{1.03^2+0.64^2}/1.03 \approx 1.18$.

## Table 2 - Torsion term accuracy vs MMFF94s Cartesian local energy
Local energy = torsion term only (ICFF) vs torsion+bend+stretch+OOP+1-4vdW (MMFF94s rigid). Columns: E_cut (kcal) | N comparisons | ICFF RMSD (kcal) | MMFF94s-rigid-local RMSD (kcal).

| E_cut | N | ICFF RMSD | MMFF94s rigid RMSD |
|---|---|---|---|
| 1 | 1010 | 0.12 | 0.23 |
| 2 | 2125 | 0.20 | 0.42 |
| 5 | 2651 | 0.26 | 0.58 |
| 7 | 2841 | 0.30 | 0.66 |
| 10 | 2956 | 0.33 | 0.80 |
| 20 | 3091 | 0.42 | 1.23 |
| All | 3190 | 0.53 | 1.53 |

(110 dihedral angles profiled. Note: leftmost E_cut column values 1,2,5,7,10,20,All reconstructed from monotonic N; original OCR fused columns.)
<!-- CHECK: Table 2 E_cut column was garbled in OCR (rows read "1010, 2125, 2651..."); the first number is N (comparisons), not E_cut. E_cut values inferred as 1,2,5,7,10,20,All from Tables 3/4 pattern. Verify against original PDF. -->

## Table 3 - Full-energy torsion profiles: "1-5,1-6" repulsion function accuracy vs MMFF94s Cartesian
Columns: E_cut (kcal) | N | ICFF (soft, C=0.55) RMSD | ICFF no-1-5,6 RMSD | ICFF with MMFF94 hard vdW RMSD.

| E_cut | N | ICFF soft | ICFF no1-56 | ICFF hard vdW |
|---|---|---|---|---|
| 1 | 796 | 0.21 | 0.33 | 0.87 |
| 3 | 1869 | 0.31 | 0.54 | 1.67 |
| 5 | 2414 | 0.37 | 0.62 | 3.49 |
| 7 | 2652 | 0.45 | 0.73 | 12.37 |
| 10 | 2851 | 0.56 | 0.87 | 17.73 |
| 20 | 3079 | 1.02 | 1.22 | 61.13 |
| All | 3190 | 1.40 | 1.31 | 60.59 |

- Hard MMFF94 vdW in rigid geometry blows up (up to 61 kcal RMSD); soft term (eq. 2) keeps it ~0.37 kcal at 5-kcal cutoff.
- Example clash (vinyl formate, C=C-O-C, 1-6 O-H contact): hard MMFF94s vdW reaches **24 kcal**; relaxed/flexible geometry never exceeds **1 kcal**.

## Table 4 - ICFF accuracy vs covalent geometry generation method
Covalent geometry from: local MMFF94s min | full MMFF94s min | MP2 QM opt. Columns: E_cut | N | ICFF(local) | ICFF//MMFF94s_full | ICFF//MP2.

| E_cut | N | ICFF (local) | ICFF//MMFF94s_full | ICFF//MP2 |
|---|---|---|---|---|
| 1 | 796 | 0.21 | 0.23 | 0.23 |
| 3 | 1869 | 0.31 | 0.34 | 0.36 |
| 5 | 2414 | 0.37 | 0.43 | 0.49 |
| 7 | 2652 | 0.45 | 0.51 | 0.56 |
| 10 | 2851 | 0.56 | 0.66 | 0.63 |
| 20 | 3079 | 1.02 | 1.10 | 0.99 |
| All | 3190 | 1.40 | 1.55 | 1.23 |

- "Local energy" geometry gives best low-energy accuracy; difference between methods < 25% across cutoffs (low sensitivity to bonded-geometry generation).

## Eq. 2 boundary-condition unit tests (for the soft repulsion implementation)
- given $R_{IJ}=R_{IJ}^*$: expect $E_{\text{vw}(1\text{-}5,1\text{-}6)} = \varepsilon_{IJ}$ and $dE/dR = 0$.
- given $R_{IJ}=R_{IJ}^0$: expect $E_{\text{vw}(1\text{-}5,1\text{-}6)} = 0$ (matches source vdW zero crossing).
- given $R_{IJ} \ge R_{IJ}^*$: use original attractive vdW branch (eq. 2 not applied).
- given $R_{IJ}=0$: expect finite (non-infinite) energy.
