# Checks / fixtures — Vaidehi, Jain, Goddard 1996 (NEIMO constant-temperature)

All runs: $T_B = 300$ K, 400 ps, torsional dof only (bonds and angles frozen). The reported temperature-fluctuation statistic is $((\mathcal{N}/2)\langle\delta T_{\rm calc}^2\rangle)^{1/2}$, which should equal 300 K at equilibrium (eq:44). Averages taken over the last 100 ps.

## Degrees-of-freedom bookkeeping (structural fixtures)

- polyethylene $C_pH_{2p+2}$: given atoms → $3N = 9p+6$, torsional dof $\mathcal{N} = p-1$.
- pe50 = $CH_3(CH_2)_{48}CH_3$: 99 hinges → total dof $\mathcal{N} = 99+6 = 105$; 99 independent dof (linear+angular momentum conserved).
- pe20 / pe30 / pe40 / pe50 total dof (incl. 6 base): **45 / 65 / 85 / 105**.
- PVDF50 = $CF_3(CH_2CF_2)_{49}CH_3$: 99 torsional dof.
- PVDF66 = $CF_3(CH_2CF_2)_{65}CH_3$: 132 carbons, 131 torsional dof, total $\mathcal{N} = 137$; cubic cell $a=b=c=18.0$ Å.
- (PVC20)₄ = $CH_3(CCl\,CH_2)_{19}CH_2Cl$ ×4 chains: 488 atoms, 180 torsional dof (incl. 4 base bodies); cell $a=22.8$, $b=23.1$, $c=12.1$ Å, $\alpha=94.82°$, $\beta=89.17°$, $\gamma=84.92°$.

## Analytic-formula fixtures

- $\tau_s$ lower bound (eq:38): given $\delta$, expect $\tau_s \ge 1.6\delta$. Examples the paper checks:
  - $\delta = 0.01$ ps → $\tau_s > 0.016$ ps; run at $\tau_s = 0.007$ ps **blew up**.
  - $\delta = 0.005$ ps → $\tau_s > 0.008$ ps; run at $\tau_s = 0.007$ ps **blew up**.
  - $\delta = 30$ fs → $\tau_s \ge 0.048$ ps.
  - $\delta = 44$ fs → $\tau_s = 0.07$ ps still consistent.
- Total-time bound (eq:41): given $\mathcal{N}, \tau_s$, expect $t_{\rm total} = 40\pi\sqrt{\mathcal{N}}\,\tau_s$. Example: $\tau_s = 10$ ps with $\mathcal{N}=99$ → $t_{\rm total} > 12600$ ps required (400 ps run at $\tau_s=10$ ps gave poor/blown-up results).
- eq:37→38 constant: $10/(2\pi) = 1.59 \approx 1.6$.

## Table 1 — pe50, $\delta = 5$ fs, $T_B = 300$ K

| $\tau_s$ (ps) | $\langle T_{\rm calc}\rangle$ (K) | $((\mathcal{N}/2)\langle\delta T_{\rm calc}^2\rangle)^{1/2}$ (K) |
|---|---|---|
| 0.007 | blew up | — |
| 0.01 | 303.16 | 240.15 |
| 0.03 | 302.13 | 156.43 |
| 0.05 | 301.68 | 238.96 |
| 0.07 | 302.74 | 295.31 |
| 0.1 | 302.32 | 285.88 |
| 0.3 | 302.95 | 304.35 |
| 0.5 | 302.13 | 367.47 |
| 0.8 | 302.46 | 354.77 |
| 1.0 | 304.96 | 923.79 |
| 10.0 | blew up | — |

Recommended for pe50 @ $\delta=5$ fs: $0.05 \le \tau_s \le 0.10$ ps (eq:45). $\langle T_{\rm calc}\rangle$ closest to 300 at $\tau_s=0.05$; fluctuation statistic closest to 300 at $\tau_s=0.07$ and 0.3.

## Table 2 — pe50, $\delta = 10$ fs, $T_B = 300$ K

| $\tau_s$ (ps) | $\langle T_{\rm calc}\rangle$ (K) | $((\mathcal{N}/2)\langle\delta T_{\rm calc}^2\rangle)^{1/2}$ (K) |
|---|---|---|
| 0.01 | blew up | — |
| 0.03 | 303.93 | 142.39 |
| 0.05 | 305.71 | 223.08 |
| 0.07 | 305.29 | 276.19 |
| 0.1 | 304.99 | 335.41 |
| 0.3 | 303.55 | 705.89 |
| 0.5 | 308.10 | 1591.40 |
| 0.8 | 307.58 | 1868.88 |
| 1.0 | 306.93 | 2829.92 |
| 10.0 | blew up | — |

$\langle T_{\rm calc}\rangle$ too high by 3–6 K (large-$\delta$ integration error). Recommended $0.05 \le \tau_s \le 0.07$ ps (eq:46).

## Table 3 — PVDF50, $\delta = 10$ fs, $T_B = 300$ K

| $\tau_s$ (ps) | $\langle T_{\rm calc}\rangle$ (K) | $((\mathcal{N}/2)\langle\delta T_{\rm calc}^2\rangle)^{1/2}$ (K) |
|---|---|---|
| 0.01 | blew up | — |
| 0.03 | 301.67 | 120.91 |
| 0.05 | 301.54 | 311.23 |
| 0.07 | 302.10 | 181.99 |
| 0.1 | 301.66 | 160.00 |
| 0.3 | 302.45 | 389.79 |
| 0.5 | 301.02 | 548.20 |
| 0.8 | 302.84 | 662.94 |
| 1.0 | 308.60 | 1134.51 |
| 10.0 | 379.21 | 28501.37 |

Both stats closest to 300 K at $\tau_s = 0.1$ ps. Recommended $0.05 < \tau_s < 0.10$ (eq:47).

## Overall recommendation and time-step claims

- Recommended thermostat time constant (combining eq:45–47): **$0.05 < \tau_s < 0.07$ ps** (eq:48). Applications use $\tau_s = 0.05$ ps ($0.07$ ps for $\delta = 25$–30 fs).
- pe50 torsion frequencies (velocity autocorrelation): 200–600 cm⁻¹ → periods 166–55 fs → suggested $\delta \approx 6$–16 fs.
- NEIMO–Hoover stable time steps: **20–30 fs** (isolated chains 25–30 fs), vs Cartesian–Hoover limited to **$\le 2$–3 fs** → ~**10× larger** time step for equal total-energy conservation.
- Amorphous (PVC20)₄ under PBC: $\delta = 20$ fs gives stable dynamics.
- Verlet predictor/corrector convergence criterion: velocity change < **0.001 MD units** (time unit = 0.0488 ps); converges in 1–2 iterations.
- Energy minimization pre-run: rms force < 0.1 (kcal/mol)/Å (PVDF66, PVC); pe50 minimized to 0.09 (kcal/mol)/Å.
- Scaling: NEIMO cost $O(\mathcal{N})$ vs $O(\mathcal{N}^3)$ for explicit $\mathcal{M}^{-1}$; for 1001 atoms (p=333) explicit method inverts a 332×332 matrix each step.
