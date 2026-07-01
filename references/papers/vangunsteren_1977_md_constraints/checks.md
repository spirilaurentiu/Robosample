# Checks / fixtures - van Gunsteren & Berendsen 1977

Concrete numbers an implementation of the Gear/Verlet predictor-corrector and the SHAKE-constraint scheme can be regression-tested against.

## 1. Gear corrector vector `a` (N-representation, 2nd-order ODE)

From Gear (1971), p. 154. Index $i = 0 \dots k-1$. Note $a_2 = 1$ for all $k$.

| k | a0 | a1 | a2 | a3 | a4 | a5 | a6 | a7 |
|---|----|----|----|----|----|----|----|----|
| 4 | 1/6 | 5/6 | 1 | 1/3 | | | | |
| 5 | 19/120 | 3/4 | 1 | 1/2 | 1/12 | | | |
| 6 | 3/20 | 251/360 | 1 | 11/18 | 1/6 | 1/60 | | |
| 7 | 863/6048 | 665/1008 | 1 | 25/36 | 35/144 | 1/24 | 1/360 | |
| 8 | 275/2016 | 19087/30240 | 1 | 137/180 | 5/16 | 17/240 | 1/120 | 1/2520 |

- **given** k and the table row, **expect** these exact rational corrector coefficients in eq 2.9.

## 2. Predictor matrix A (N-representation) = Pascal triangle, k<=8

Element (row $r$, col $c$), 0-based: $A_{rc} = \binom{c}{r}$ (0 for $c<r$). For k=8:

```
1 1 1 1  1  1  1  1
0 1 2 3  4  5  6  7
0 0 1 3  6 10 15 21
0 0 0 1  4 10 20 35
0 0 0 0  1  5 15 35
0 0 0 0  0  1  6 21
0 0 0 0  0  0  1  7
0 0 0 0  0  0  0  1
```
- **given** k, **expect** A = upper-triangular Pascal (binomial) matrix.

## 3. F-representation corrector vector b (eq A.7)

- **given** the Gear `a` (table above), **expect** $b_0=a_0$, $b_1=a_1$, $b_2=a_2=1$, and $b_i=0$ for $i>2$.

## 4. F-representation predictor matrices B = T A T^{-1}

Exact values printed in the paper (rows are matrix rows; top-left is B[0][0]).

### k=4
```
1  1   4/3  -1/3
0  1   3    -1
0  0   2    -1
0  0   1     0
```

### k=5
```
1  1  19/12  -5/6   1/4
0  1  23/6   -8/3   5/6
0  0   3     -3     1
0  0   1      0     0
0  0   0      1     0
```

### k=6
```
1  1  323/180  -22/15   53/60   -19/90
0  1   55/12   -59/12   37/12   -3/4
0  0    4       -6       4       -1
0  0    1        0       0        0
0  0    0        1       0        0
0  0    0        0       1        0
```

### k=7
```
1  1  1427/720   -133/60    241/120   -173/180   3/16
0  1  1901/360   -1387/180  109/15    -637/180   251/360
0  0    5         -10        10        -5         1
0  0    1          0          0         0         0
0  0    0          1          0         0         0
0  0    0          0          1         0         0
0  0    0          0          0         1         0
```
<!-- CHECK: last two printed rows of the k=7 B in the source had an OCR-doubled trailing column; reconstructed as an identity shift so the last k-3 rows are "zero except a single 1" per eq A.7 discussion. -->

### k=8
```
1  1  2713/1260  -15487/5040  1172/315   -6737/2520  263/252   -863/5040
0  1  4277/720   -2641/240    4991/360   -3649/360   959/240   -95/144
0  0    6         -15          20         -15         6         -1
0  0    1          0            0          0          0          0
0  0    0          1            0          0          0          0
0  0    0          0            1          0          0          0
0  0    0          0            0          1          0          0
0  0    0          0            0          0          1          0
```

- **structural invariant:** the last $(k-3)$ rows of B are all zeros except a single 1 (identity-shift), and the last $(k-3)$ entries of b are 0. Row 2 (0-based) of B is the finite-difference / binomial stencil with alternating signs: e.g. k=8 row 2 = [6, -15, 20, -15, 6, -1] (binomials of $(1-x)^? $ pattern) preceded by 0,0.

## 5. T^{-1} (N<-F) fixtures

The paper prints $\mathbf{T}^{-1}$ for k=4..8; representative:

### k=4
```
1  0   0     0
0  1   0     0
0  0   1     0
0  0  1/3  -1/3
```
- **given** k=4, **expect** T^{-1} bottom row [0, 0, 1/3, -1/3].

## 6. BPTI system parameters (regression sanity)

- Protein: BPTI, 58 residues, single molecule in vacuum.
- Extended (united) atoms: 458 (= 454 in BPTI + 4 H-bonded waters).
- Cartesian degrees of freedom: 1374 (= 3 x 458).
- Possible bond-length constraints: 468. Possible angle constraints: 626.
- Interaction function: Gelin-Karplus (bond stretch, angle bend, dihedral, H-bond, VdW, Coulomb).

## 7. Integration / accuracy benchmarks (CDC Cyber 74-16, 100-step runs)

Time-step unit = $4.889\times10^{-14}$ s.

| quantity | value |
|---|---|
| Equilibration time step | $9.778\times10^{-16}$ s (= h=0.02 in figure units) |
| Velocity scaling (150 K -> 300 K) | multiply all velocities by 1.74 |
| Verlet CPU/step | 1.2 s; storage = 3 vectors |
| Gear CPU/step (k=4) | 1.3 s; storage = k+1 vectors |
| Gear CPU/step (k=8) | 1.5 s |
| Verlet + SHAKE CPU/step | ~1.6 s; +1 stored vector |
| Stability limit (no constraints) | algorithms become unstable beyond h ~ 0.02 |
| Crossover Gear vs Verlet accuracy | h ~ 0.03 (below: Gear better; above: Verlet better) |
| Optimal k (small step, no constraints) | k = 7 |
| E_tot fluctuation vs E_kin fluctuation (h<0.02) | E_tot fluctuation > 100x smaller than E_kin fluctuation |
| Drift/step at h=0.01, Gear | ~$10^{-6}$ kcal/mole |
| Drift/step at h=0.01, Verlet | ~$10^{-3}$ kcal/mole |
| Adequate SHAKE tol (for h>0.01) | $10^{-6}$ (energy fluctuation changes <1% on further decrease) |
| Constraint dynamics: allowed time step | ~4x larger than non-constraint for equal accuracy |
| Constraint dynamics per-step cost | ~1.3x slower than non-constraint |
| Net speedup from bond-length constraints | ~factor 3 |
| Best k for constraint dynamics | k=4 (larger k gives no improvement) |

- **qualitative expectation:** bond-angle constraints give no speedup (SHAKE cost grows, no time-step increase allowed) — do not implement them for macromolecules.
