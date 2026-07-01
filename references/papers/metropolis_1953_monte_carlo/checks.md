# Checks / Fixtures - Metropolis et al. 1953

## Simulation setup constants

- Given the actual calculation: expect N = 224 particles (indices 0..223),
  arranged as a trigonal lattice, 14 per row x 16 per column, in a unit square
  of unit area.
- Given the lattice spacing: expect d = 1/14 ≈ 0.0714286 (six nearest neighbors
  at ~equal distance d).
- Given the step size choice in the runs: expect maximum displacement
  $\alpha = d - d_0$; this made about half the moves forbidden.
- Given a cycle: it consists of moving every particle once. Averaging/extrapolation
  done every 16 cycles. Production runs: ~48 to 64 cycles per state point.
- Given the MANIAC: ~3 minutes per cycle; ~4-5 hours per pressure-curve point.
- A consistency check reported: runs with 56 particles and with 224 particles
  agreed within statistical error. Average error ~3 percent.

## Eq. (11a/11b) parameterization (verify numerically)

- Given $d = 1/14$, expect $d_0 = d(1 - 2^{\nu-8})$.
- Given the constant in Eq. (11b): expect
  $(A/A_0) = 1/[0.98974329 (1 - 2^{\nu-8})^2]$.
- Spot check $\nu = 5$: $2^{\nu-8} = 2^{-3} = 0.125$, so $(1-0.125)^2 = 0.765625$,
  and $A/A_0 = 1/(0.98974329 \times 0.765625) = 1.31966$. Matches Table I row
  $\nu=5$. <!-- CHECK: verified by hand -->
- Spot check $\nu = 2$: $2^{-6} = 0.015625$, $(0.984375)^2 = 0.969194$,
  $A/A_0 = 1/(0.98974329 \times 0.969194) = 1.04269$. Matches Table I.

## Table I - Results (rigid spheres, 2D)

Columns: $\nu$, $A/A_0$, $X_1 = (PA/NkT)-1$ from this paper, $X_2 = $ free volume
theory, $X_3 = $ four-term virial expansion, and $PA_0/NkT$ from this paper.

| $\nu$ | $A/A_0$ | $X_1$ (MC) | $X_2$ (free vol) | $X_3$ (virial) | $PA_0/NkT$ |
|------|---------|--------|-------|--------|--------------|
| 2    | 1.04269 | 49.17  | 47.35 | 9.77   | 48.11        |
| 4    | 1.14957 | 13.95  | 13.85 | 7.55   | 13.01        |
| 5    | 1.31966 | 6.43   | 6.72  | 5.35   | 5.63         |
| 5.5  | 1.4909  | 4.41   | 4.53  | 4.02   | 3.63         |
| 6    | 1.7962  | 2.929  | 2.939 | 2.680  | 2.187        |
| 6.25 | 2.04616 | 2.186  | 2.323 | 2.065  | 1.557        |
| 6.5  | 2.41751 | 1.486  | 1.802 | 1.514  | 1.028        |
| 7    | 4.04145 | 0.6766 | 0.990 | 0.667  | 0.4149       |

<!-- CHECK: raw table printed X3 header as "X_2" (typo) and value "2,929" for
     nu=6 X_1 (comma is a decimal point -> 2.929); corrected here from context. -->

Notes for use as fixtures:
- $X_1$ is the Monte Carlo result; a correct MC EOS implementation should
  reproduce these within ~3 percent statistical error.
- $X_2$ (free volume, Wood 1952) and $X_1$ agree well at small area
  ($A/A_0 < 1.8$); deviation begins with a sudden break near $\nu=6$
  ($A/A_0 \sim 1.8$).
- $X_3$ (four-term virial) and $X_1$ agree well for $A/A_0 > 2.5$.

## Radial distribution / extrapolation example ($\nu = 5$)

- Given $\nu=5$, $A/A_0 = 1.31966$, $K = 1.5$: three straight-line least-square
  extrapolations of the first 16 $N_m$ gave $N_{1/2}^{(1)} = 6367$,
  $N_{1/2}^{(2)} = 6160$, $N_{1/2}^{(3)} = 6377$.
- Their average: $\bar{N}_{1/2} = 6301$.
- Resultant $(PA/NkT)-1 = 64\,\bar{N}_{1/2}/[N^2(K^2-1)]$. With $N=224$, $K=1.5$:
  $64 \times 6301 / (224^2 \times (2.25-1)) = 403264 / (50176 \times 1.25)
  = 403264/62720 = 6.43$. Matches reported 6.43. <!-- CHECK: verified by hand -->

## Virial coefficients (Section V)

- Given $C_1 = \pi/3^{1/2}$, expect $C_1 = 1.813799$.
- Given cluster ratio: expect $A_{4,6}/A_{4,4} = 0.752 \ (\pm 0.002)$.
- Given $C_4 = 8\pi^3(0.585)/135$: expect $C_4 \approx 3.38$ ($\pm \sim 5$ percent).
  Numeric: $8 \times 31.00628 \times 0.585/135 = 145.109/135 = 1.075$?
  <!-- CHECK: direct evaluation gives ~1.075, but Eq.(14) reports C4=3.38;
       the 0.585 is the bracketed cluster-integral combination, not the full C4,
       and additional numerical factors/normalizations of the A_{i,k} enter.
       Trust the tabulated final coefficients in Eq. (14) as the authoritative
       fixture, not this back-of-envelope. -->
- Final four-term virial EOS (Eq. 14), authoritative coefficients:
  $C_1 = 1.813799$, $C_2 = 2.57269$, $C_3 = 3.179$, $C_4 = 3.38$:
  $(PA/NkT)-1 = 1.813799(A_0/A) + 2.57269(A_0/A)^2 + 3.179(A_0/A)^3 + 3.38(A_0/A)^4$.

## Algorithm regression fixture (Metropolis acceptance)

- Given a symmetric proposal and $\Delta E \le 0$: expect acceptance probability
  = 1 (always accept).
- Given $\Delta E > 0$: expect acceptance probability $= \exp(-\Delta E/kT)$;
  accept iff drawn $\xi_3 \sim U(0,1)$ satisfies $\xi_3 < \exp(-\Delta E/kT)$.
- Given a hard-sphere move creating an overlap: expect rejection (return particle
  to old position) since $\Delta E = +\infty$.
- Stationary distribution the chain must reproduce: $\nu_r \propto \exp(-E_r/kT)$
  (Boltzmann/canonical).
