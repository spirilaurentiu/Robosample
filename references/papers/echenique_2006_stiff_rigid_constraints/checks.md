# Checks / fixtures - Echenique 2006 (stiff vs rigid constraints)

System: model dipeptide HCO-L-Ala-NH2. n=16 atoms, N=48 DOF, M=2 soft internals (φ,ψ),
L=40 hard internals. Levels of theory: MP2/6-31++G(d,p) and HF/6-31++G(d,p).
All energy functions referenced to zero in the grid. Units: kcal/mol unless RT stated. RT≈0.6 kcal/mol at 300 K.

## System configuration
- Ramachandran grid: 12×12, φ,ψ ∈ [-165°, 165°] step 30°, giving 144 conformations.
- Displaced grids for finite-difference derivatives: +2° in φ and +2° in ψ.
- Convergence criteria (GAMESS): OPTTOL = 1e-5 (gradient), CONV = 1e-6 (SCF).
- Absolute PES minimum in the grid (MP2/6-31++G(d,p)): **-416.0733418995 hartree**.

## Table 4 - Max variation / Average / Std over the 12×12 grid (kcal/mol)
Format: quantity -> (Max, Ave, Std) at [MP2] / [HF].
- V_Σ:      (21.64, 6.76, 3.88) / (23.62, 6.92, 4.35)
- F_s:      (21.43, 6.47, 3.93) / (23.78, 7.17, 4.38)
- F_r:      (21.09, 6.46, 3.82) / (23.09, 6.76, 4.31)
- -T S_s^k: (0.24, 0.09, 0.05) / (0.23, 0.09, 0.04)   [mass-metric G term, smallest]
- -T S_s^c: (1.67, 0.98, 0.32) / (1.34, 0.63, 0.30)   [Hessian term, largest]
- -T S_r^k: (0.81, 0.37, 0.12) / (0.75, 0.38, 0.12)   [reduced g term]
- V_F:      (1.68, 0.89, 0.30) / (1.35, 0.55, 0.27)   [Fixman potential]

## Table 5 - Pearson correlation r_12 with V_Σ (grid)
MP2/6-31++G(d,p):
- V_Σ vs -T S_s^c : 0.1572
- V_Σ vs -T S_s^k : -0.0008
- V_Σ vs -T S_r^k : -0.3831
- V_Σ vs V_F      : 0.3334
HF/6-31++G(d,p):
- V_Σ vs -T S_s^c : 0.0682
- V_Σ vs -T S_s^k : 0.0897
- V_Σ vs -T S_r^k : -0.3544
- V_Σ vs V_F      : 0.2404
MP2 vs HF cross-correlation (same term, two theory levels):
- -T S_s^c : 0.9136
- -T S_s^k : 0.9808
- -T S_r^k : 0.9316
- V_F      : 0.9217

## Table 6 - Statistical distance d_12 (in RT), N_res, slope b_12, r_12 - FULL grid (144 confs)
Given (correction term dropped, V_1 reference, V_2 approx): expect (d_12, N_res, b_12, r_12).
MP2/6-31++G(d,p):
- (-TS_s^k -TS_s^c; F_s; V_Σ):              (0.74 RT, 1.82, 0.98, 0.9967)
- (-TS_s^c;        F_s; V_Σ-TS_s^k):        (0.74 RT, 1.83, 0.98, 0.9967)
- (-TS_s^k;        F_s; V_Σ-TS_s^c):        (0.11 RT, 80.45, 1.00, 0.9999)
- (-TS_r^k;        F_r; V_Σ):               (0.29 RT, 11.62, 1.01, 0.9995)
- (V_F;            F_s; F_r):               (0.67 RT, 2.24, 0.97, 0.9972)
HF/6-31++G(d,p):
- (-TS_s^k -TS_s^c; F_s; V_Σ):              (0.73 RT, 1.90, 0.99, 0.9975)
- (-TS_s^c;        F_s; V_Σ-TS_s^k):        (0.71 RT, 2.00, 0.99, 0.9976)
- (-TS_s^k;        F_s; V_Σ-TS_s^c):        (0.10 RT, 90.99, 1.00, 0.9999)
- (-TS_r^k;        F_r; V_Σ):               (0.26 RT, 14.83, 1.01, 0.9997)
- (V_F;            F_s; F_r):               (0.61 RT, 2.69, 0.98, 0.9982)
MP2 vs HF (same potential, two theory levels):
- V_Σ vs V_Σ : (1.25 RT, 0.64, 1.12, 0.9925)
- F_s  vs F_s : (1.18 RT, 0.72, 1.11, 0.9934)
- F_r  vs F_r : (1.18 RT, 0.72, 1.12, 0.9932)

## Table 7 - Secondary-structure Ramachandran angles (degrees) [ref 88]
- α-helix:            φ=-57,  ψ=-47  (raw ψ blank; standard value -47)  <!-- CHECK: ψ missing in OCR table for α-helix; conventional value ≈ -47 -->
- 3_10-helix:         φ=-49,  ψ=-26
- π-helix:            φ=-57,  ψ=-70
- polyproline II:     φ=-79,  ψ=149
- parallel β-sheet:   φ=-119, ψ=113
- antiparallel β-sheet: φ=-139, ψ=135  (minimum-energy reference in Fig 4)

## Table 8 - Statistical distance d_12 (RT), N_res, b_12, r_12 - 6 secondary-structure confs
MP2/6-31++G(d,p):
- (-TS_s^k -TS_s^c; F_s; V_Σ):       (0.22 RT, 19.72, 0.99, 0.9990)
- (-TS_s^c;        F_s; V_Σ-TS_s^k): (0.26 RT, 14.07, 0.98, 0.9985)
- (-TS_s^k;        F_s; V_Σ-TS_s^c): (0.06 RT, 298.13, 1.01, 0.9999)
- (-TS_r^k;        F_r; V_Σ):        (0.20 RT, 25.64, 0.99, 0.9992)
- (V_F;            F_s; F_r):        (0.34 RT, 8.73, 0.99, 0.9977)
HF/6-31++G(d,p):
- (-TS_s^k -TS_s^c; F_s; V_Σ):       (0.14 RT, 47.94, 1.00, 0.9997)
- (-TS_s^c;        F_s; V_Σ-TS_s^k): (0.15 RT, 46.12, 1.00, 0.9997)
- (-TS_s^k;        F_s; V_Σ-TS_s^c): (0.05 RT, 380.30, 1.00, 0.9999)
- (-TS_r^k;        F_r; V_Σ):        (0.15 RT, 41.85, 0.99, 0.9997)
- (V_F;            F_s; F_r):        (0.18 RT, 30.12, 1.01, 0.9996)
MP2 vs HF:
- V_Σ: (0.77 RT, 1.68, 1.28, 0.9929)
- F_s: (0.77 RT, 1.69, 1.26, 0.9928)
- F_r: (0.71 RT, 1.96, 1.28, 0.9939)

## Derived cross-checks (verify implementation)
- N_res = (RT/d_12)^2  (eq 33). E.g. d_12=0.11 RT -> N_res=(1/0.11)^2=82.6 (table lists 80.45; close, small rounding). d_12=0.29 RT -> 11.89 (table 11.62). d_12=0.74 RT -> 1.826 (table 1.82). PASS as consistency check within rounding.
- d_12 = √2 · σ_2 · (1 - r_12^2)^{1/2}  (eq 34).
- σ_2 (std of V_2 over working set): ≈ 4 kcal/mol for the full grid, ≈ 2 kcal/mol for the 6 secondary-structure elements.

## Qualitative conclusions (as regression expectations)
- Ordering of correction-term magnitude: |−TS_s^c| (Hessian) > |−TS_r^k| (reduced g) > |−TS_s^k| (mass-metric G).
- Drop-tolerance residue limits (whole Ramachandran space, MP2): −TS_s^k up to ~80 residues; −TS_r^k up to ~12; −TS_s^c and V_F only up to ~2.
- Restricting to the low-energy secondary-structure region multiplies these residue limits by ≈ 4.
- HF gives results very similar to MP2 (~1/10 CPU cost); differences are mostly a linear scaling (b_12 > 1).
