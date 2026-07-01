# Checks / fixtures - Shirts 2012, Ensemble validation

Regression fixtures for a port of the ensemble-consistency test
(`checkensemble`, hosted at simtk.org/home/checkensemble).

## Analytic identities (exact, no simulation)

- D-dimensional harmonic oscillator, equal spring K, x_0=0:
  - Partition function `Q(beta) = (2*pi/(beta*K))^(D/2)`.
  - Free energy `A(beta) = -(D/(2*beta)) * ln(2*pi/(beta*K))`.
  - Average energy `<E> = D/(2*beta)`.
  - Given D=20, K=1, beta=0.6: expect `<E> = 20/(2*0.6) = 16.667`.
  - Given D=20, K=1, beta=1.3: expect `<E> = 20/2.6 = 7.692` (paper "7.7 k_B T").
  - Given D=20, K=1, beta=0.7: expect `<E> = 20/1.4 = 14.286` (paper "14.3 k_B T").
- Canonical log-ratio slope: for any pair, expect fitted slope of ln[P(E|b2)/P(E|b1)] vs E to equal `-(b2 - b1)`.
- Kinetic energy: mean `<E_kin> = 3N/(2*beta) = (k_B T/2)*DOF`; variance `sigma^2 = 3N/(2*beta^2)`.
- NPT toy: `Delta(P,beta) = (beta*P)^-2 * sqrt(2*pi/(a^2*beta))`; Gaussian width of x at fixed V is `sigma = V/(a*sqrt(beta))`; average rejection-sampling efficiency `exp(-beta/2)`.

## Toy harmonic oscillator, canonical (Table 1)

Setup: K=1, D=20, beta=1.3 and 0.7 (so true `beta_2 - beta_1 = 0.6`);
500,000 samples per temperature, 200 repetitions. Add noise dE = nu*|N(0,1)|.
Fitted slope +/- error (sigma-from-true in parentheses):

| nu | linear (analytic) | nonlinear (analytic) | max-likelihood (analytic) |
|---|---|---|---|
| 0.0    | 0.6006 +/- 0.0012 (0.5) | 0.6028 +/- 0.0013 (2.3) | 0.6016 +/- 0.0012 (1.3) |
| 0.005  | 0.5955 +/- 0.0012 (3.7) | 0.6019 +/- 0.0012 (1.6) | 0.5973 +/- 0.0012 (2.3) |
| 0.0075 | 0.5927 +/- 0.0012 (6.0) | 0.5969 +/- 0.0012 (2.5) | 0.5939 +/- 0.0012 (5.2) |
| 0.01   | 0.5924 +/- 0.0012 (6.3) | 0.5916 +/- 0.0012 (6.9) | 0.5936 +/- 0.0012 (5.5) |
| 0.02   | 0.5841 +/- 0.0012 (13.3)| 0.5899 +/- 0.0012 (8.3) | 0.5850 +/- 0.0012 (13.0)|

- Expect: at nu=0 all methods recover ~0.600. Deviations become >3 sigma consistently for nu >= 0.0075 (< 1% of k_B T).
- Bootstrap (200 samples) and 200-replicate sample std match the analytic errors for linear and ML; nonlinear analytic error underestimates (its sample/bootstrap error is ~0.003 vs analytic ~0.0012).

## Optimal temperature-gap sweep (Table 2)

Setup: HO, nu=0.01, (beta_1+beta_2)/2 = 1 fixed, ML analytic error, 500k samples each.
Expect estimated slope ~ true (b2-b1), with sigma-deviation peaking at intermediate gap:

| b2 | b1 | true b2-b1 | estimated | sigma dev |
|---|---|---|---|---|
| 1.05 | 0.95 | 0.1 | 0.0993 +/- 0.0006 | 1.1 |
| 1.20 | 0.80 | 0.4 | 0.3960 +/- 0.0009 | 4.7 |
| 1.30 | 0.70 | 0.6 | 0.5936 +/- 0.0012 | 5.4 (peak) |
| 1.50 | 0.50 | 1.0 | 0.9907 +/- 0.0027 | 3.5 |
| 1.70 | 0.30 | 1.4 | 1.3916 +/- 0.0100 | 0.8 |

- Rule of thumb validation: at peak (b1=0.7,b2=1.3) distribution-center gap 6.6 k_B T ~ sum of std devs 6.9 k_B T.

## NPT harmonic-oscillator toy (Section 2.7.1)

250,000 samples per paired distribution, maximum-likelihood.
- Enthalpy: beta_1=2/3, beta_2=2, P_1=P_2=1. Expect `beta_2-beta_1 = 1.3341 +/- 0.0040` (true 4/3=1.3333; 0.2 quantiles off).
- Volume: beta_1=beta_2=1.0, P_1=1.3, P_2=0.7. Expect `beta(P_2-P_1) = -0.6013 +/- 0.0025` (true -0.6; 0.53 quantiles).
- Joint (E,V): beta_1=0.6, beta_2=0.8, P_1=0.8, P_2=1.2. Expect slope `(beta_2-beta_1) = 0.20035 +/- 0.00318` (true 0.2, 0.1 quantile) and `(beta_2 P_2 - beta_1 P_1) = -0.48129 +/- 0.00185` (true -0.48, 0.7 quantile).

## Lennard-Jones argon MD (system params, Section 3.2)

- Rowley-Nicholson-Parsonage argon: sigma = 0.3405 nm, epsilon = 119.8 K, k_B = 0.996072 kJ/mol.
- rho = 0.85 rho_c -> box length 3.5328256 nm; T = 0.85 T_c = 135.0226 K; 300 LJ particles.
- Reduced time unit sigma*(M/epsilon)^(1/2) = 0.1245 ps; 48 fs step (~0.386 reduced) is stability limit.
- Default temperatures T_low=132.915 K, T_high=137.138 K (~0.7x ideal gap; C_V ~= 8.5 kJ/(mol K)).

## Thermostat validation, argon (Table 3)

True slope 0.027865 (k_B T)^-1, equivalent Delta_T = 4.223 K. Estimated Delta_T (sigma dev):

| thermostat | total | potential | kinetic |
|---|---|---|---|
| None (NVE) | N/A (E const) | 4.388 +/- 0.115 (1.4) | 3.048 +/- 0.112 (10.5) |
| Berendsen  | 9.369 +/- 0.122 (42.2) | 4.606 +/- 0.086 (4.5) | 29.034 +/- 0.364 (68.3) |
| Stochastic | 4.172 +/- 0.066 (0.8) | 4.098 +/- 0.081 (1.6) | 4.251 +/- 0.091 (0.3) |
| Nose-Hoover| 4.197 +/- 0.067 (0.4) | 4.220 +/- 0.082 (0.03)| 4.186 +/- 0.090 (0.4) |
| Andersen   | 4.212 +/- 0.066 (0.2) | 4.226 +/- 0.081 (0.03)| 4.226 +/- 0.090 (0.03)|
| Bussi-Parrinello | 4.167 +/- 0.066 (0.8) | 4.272 +/- 0.082 (0.6) | 4.155 +/- 0.089 (0.8) |

- Expect: all correct thermostats within ~1 sigma; Berendsen wildly off (kinetic 68 sigma). NVE kinetic deviates (10.5 sigma) but NVE potential does not.

## Step-size effect, argon Bussi-Parrinello (Table 4)

True T_low=132.915 K, true Delta_T=4.223 K. Estimated Delta_T (sigma dev):

| dt (fs) | total | potential | kinetic |
|---|---|---|---|
| 8  | 4.230 +/- 0.047 (0.2) | 4.237 +/- 0.058 (0.2) | 4.186 +/- 0.063 (0.6) |
| 16 | 4.183 +/- 0.032 (1.2) | 4.253 +/- 0.040 (0.8) | 4.106 +/- 0.043 (2.7) |
| 24 | 4.058 +/- 0.026 (6.4) | 4.140 +/- 0.032 (2.6) | 4.023 +/- 0.035 (5.8) |
| 32 | 3.967 +/- 0.030 (8.6) | 4.199 +/- 0.028 (0.9) | 4.054 +/- 0.022 (7.6) |
| 40 | 3.988 +/- 0.020 (11.6)| 4.178 +/- 0.026 (1.7) | 3.877 +/- 0.027 (12.9)|
| 40 (E_kin ave, leapfrog) | 4.266 +/- 0.021 (2.6) | 4.275 +/- 0.026 (2.0) | 4.296 +/- 0.029 (2.6) |

- Expect: kinetic/total deviate as dt grows; potential stays near-correct. Half-step (leapfrog) KE estimator markedly better than full-step (velocity-Verlet) KE at 40 fs.

## Abrupt cutoff effect, argon (Table 5)

True Delta_T=4.223. Estimated T_low (K) sigma dev, and Delta_T for total/pot/kin:

| r_c (LJ sigma) | E drift (kJ/mol/ns) | T_low (sigma) | total | potential | kinetic |
|---|---|---|---|---|---|
| 2.0 | 9400 | 133.952 +/- 0.045 (23.0) | 4.102 (2.1) | 4.018 (2.4) | 4.122 (1.4) |
| 2.5 | 1140 | 133.043 +/- 0.045 (2.9)  | 4.206 (0.3) | 4.177 (0.8) | 4.176 (0.7) |
| 3.0 | 239  | 132.941 +/- 0.045 (0.6)  | 4.232 (0.2) | 4.291 (1.1) | 4.213 (0.1) |
| 3.5 | 104  | 132.930 +/- 0.045 (0.3)  | 4.226 (0.1) | 4.302 (1.4) | 4.135 (1.2) |
| 4.0 | 78   | 132.929 +/- 0.045 (0.3)  | 4.302 (1.6) | 4.307 (1.6) | 4.192 (0.4) |

- Expect: only cutoffs < 3 sigma show clear violations; average-KE temperature deviation (23 sigma at r_c=2) more sensitive here than the ensemble test.

## Pressure control, argon (Table 7)

P_ave=90 bar, T_ave=125 K; volume test Delta_P=120 bar (Berendsen 4 bar); enthalpy Delta_T=7.138 K.
Estimated slopes (sigma dev):

| barostat | enthalpy | volume | joint E-V (E slope) | joint (V slope) |
|---|---|---|---|---|
| Berendsen | 4.176 +/- 0.121 (19.8) | 79.5 +/- 4.4 (17.1) | 0.69 +/- 0.14 (7.6) | -318.661 +/- 7.322 (43.9) |
| Parrinello-Rahman | 7.022 +/- 0.033 (3.5) | 114.58 +/- 0.57 (9.5) | 7.168 +/- 0.036 (0.8) | 110.971 +/- 0.529 (7.4) |
| MTTK | 7.105 +/- 0.029 (1.2) | 115.51 +/- 0.50 (9.0) | 7.152 +/- 0.031 (0.5) | 111.312 +/- 0.457 (7.8) |

- Expect: PR and MTTK give correct enthalpy but volume Delta_P off by ~5 bar (~5%, 7-9 sigma). Berendsen fails all three.

## Water TIP3P (Sections 3.5, Tables 8-9)

- System: 900 TIP3P, velocity Verlet, PME (cutoff 1.0 nm, order 6, tol 1e-6), LJ switch 0.8-0.9 nm, SETTLE, dt=2 fs, 20 ns.
- NVT: T = 298 and 301 K, Delta_T=3 K -> true slope 0.004023 (k_B T)^-1.
- NVT slopes (sigma dev): Berendsen 51.6 +/- 1.1 (44.2, FAIL); Stochastic 2.998 +/- 0.059 (0.04); Nose-Hoover 2.921 (1.4); Andersen 3.028 (0.4); Bussi-Parrinello 2.955 (0.8).
- NPT (PR & MTTK): true Delta_T=3, true Delta_P=350 bar. MTTK volume 335.7 +/- 3.9 (3.7 sigma); Parrinello-Rahman volume 309.3 +/- 3.7 (11.1 sigma). Berendsen (Delta_T=1, Delta_P=30) fails.

## Statistical primitives to unit-test
- Bin probability variance: `var(p_k) = p_k*(1-p_k)/N`.
- Log-ratio variance in bin k: `1/n_{k,1} - 1/N_1 + 1/n_{k,2} - 1/N_2`.
- Fermi function: `f(x) = 1/(1+exp(-x))` (eq 8 arg convention).
- Gap rules: `DeltaT/T = sqrt(2*k_B/C_V)`; `|DeltaP| = 2*k_B*T/sigma_V = sqrt(2*k_B*T/(V*kappa_T))`.
- Flag ensemble inconsistency when fitted slope is consistently > 2-3 sigma from true across repeats.
