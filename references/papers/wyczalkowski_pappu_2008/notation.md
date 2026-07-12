# Notation - Wyczalkowski & Pappu (2008)

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $\Gamma$ | full system configuration (all coordinates) | phase-space point | integrated over in $\int d\Gamma$ |
| $\lambda$ | Kirkwood coupling parameter | dimensionless | $0 \le \lambda \le 1$; $\lambda=0$ pure solvent, $\lambda=1$ full solute |
| $\lambda_C, \lambda_{LJ}$ | Coulomb / Lennard-Jones coupling | dimensionless | scaled together, $\lambda_C = \lambda_{LJ} = \lambda$ in this work |
| $\beta$ | inverse temperature | $1/\text{energy}$ | $\beta = (k_B T)^{-1}$ |
| $T$ | temperature | K | 298 K in the numerical study |
| $U_i(\Gamma), U(\Gamma,\lambda_i)$ | potential energy of replica $i$ | energy (kcal/mol) | $U_0, U_1$ are the two swap endpoints |
| $F_i$ | Helmholtz free energy of replica $i$ | energy | eq:2 |
| $\rho_i(\Gamma)$ | equilibrium (Boltzmann) density of ensemble $i$ | probability density | normalized, eq:3 |
| $\delta F$ | free-energy change $\lambda_0 \to \lambda_1$ | energy | true (unknown) value |
| $\Delta F$ | total free energy of hydration | kcal/mol | sum of neighboring $\delta F$, eq:dF_sum |
| $\Delta G$ | Gibbs free energy of hydration | kcal/mol | literature/experiment; $\approx \Delta F$ here |
| $W^F, W^R$ | forward / reverse switching work | energy | $W^R = -W^F$; eq:4a,4b |
| $W_D^F, W_D^R$ | forward / reverse dissipated work | energy | $W_D^F = -W_D^R$ for a given $\Gamma$; eq:7a,7b |
| $W_D$ | dissipated-work value (argument of distributions) | energy | dummy variable |
| $P_F, P_R$ | ideal forward/reverse dissipated-work distributions | probability density | Crooks pair, eq:1 |
| $P_F^*, P_R^*$ | observed (finite-sample) work distributions | normalized histogram | eq:8 |
| $\delta F_{FEP}^F, \delta F_{FEP}^R$ | forward / reverse FEP estimators | energy | eq:5a (sign $-\beta^{-1}$), eq:5b (sign $+\beta^{-1}$) |
| $\epsilon_H$ | hysteresis error | energy (kcal/mol) | $\delta F_{FEP}^F - \delta F_{FEP}^R$, eq:6 |
| $(\epsilon_H)_i$ | hysteresis error between replicas $i, i+1$ | kcal/mol | |
| $\epsilon_{rms}$ | root-mean-square hysteresis error | kcal/mol | eq:eps_rms |
| $\epsilon_{FT}^*$ | fluctuation-theorem error term | energy | function of $W_D$; zero when Crooks holds, eq:8 |
| $P_{swap}$ | Metropolis swap acceptance probability | probability | used for actual RE moves, eq:10 |
| $p_{swap}$ | Fermi swap probability | probability | theoretical analysis, $f(\beta \Delta U_{swap})$ |
| $\langle p_{swap} \rangle$ | average Fermi swap probability | probability | overlap integral, eq:17b |
| $\Delta U_{swap}$ | swap energy change | energy | $= W^F + W^R = W_D^F + W_D^R$, eq:11 |
| $\gamma, \gamma'$ | native / swapped configuration pair | $(\Gamma_0, \Gamma_1)$ / $(\Gamma_1, \Gamma_0)$ | |
| $\rho_N, \rho_N'$ | native / swapped joint density | probability density | eq:12a,12b; $\rho_N'(\gamma) = \rho_N(\gamma')$ |
| $\rho_\epsilon$ | sampling error in $\rho_N'$ | probability density | small, eq:14 |
| $f(x)$ | Fermi (logistic) function | dimensionless | $1/(1+e^x)$, eq:15 |
| $C_\lambda$ | variance of $\partial U/\partial \lambda$ | energy$^2$ | $\langle (\partial U/\partial\lambda)^2\rangle_0 - \langle \partial U/\partial\lambda\rangle_0^2$, eq:18 |
| $\delta_\lambda$ (also $\delta$) | $\lambda$ spacing between neighbor replicas | dimensionless | small-parameter of the linearization |
| $\partial U/\partial \lambda$ | conjugate force to $\lambda$ | energy | $= \partial U/\partial\lambda_{LJ} + \partial U/\partial\lambda_C$ |
| $V_0, W_0$ | $\partial U/\partial\lambda$, $\partial^2 U/\partial\lambda^2$ at $\lambda_0$ | energy | Taylor coefficients, appendix A4 |
| $M$ | total number of replicas | integer | 21 in the study |
| $Q_\delta$ | partition function at $\lambda_0 + \delta$ | | expansion in appendix A4 |
| $\langle\cdot\rangle_i$ | equilibrium average over ensemble $i$ | | ideal |
| $\langle\cdot\rangle_i^*$ | finite-simulation (estimated) average | | starred |
| $r$ | interatomic distance | Angstrom | |
| $q_i, q_j$ | atomic charges | e | |
| $\alpha_C$ | Coulomb soft-core parameter | Angstrom | 1.5 Angstrom |
| $\alpha_{LJ}$ | LJ soft-core parameter | Angstrom | 0.5 Angstrom |
| $b$ (a) | LJ soft-core exponent | dimensionless | $= 4$ |
| $k$ | LJ well-depth scaling constant | dimensionless | $= 1$ |
| $\sigma, \epsilon$ | LJ diameter / well depth | Angstrom / energy | |

## Sign / convention notes

- $\beta = 1/(k_B T)$; all "work" and "free energy" quantities are energies.
- The forward FEP estimator (eq:5a) has a leading $-\beta^{-1}$; the reverse
  (eq:5b) has a leading $+\beta^{-1}$. Both estimate $\delta F$ for the same
  process $\lambda_0 \to \lambda_1$.
- Dissipated work is antisymmetric per configuration: $W_D^F(\Gamma) =
  -W_D^R(\Gamma)$. This drives the FT and the cancellation in eq:11c.
- Second law: dissipated work is on average nonnegative,
  $\langle W_D \rangle \ge 0$; the identity eq:19 requires sampling the rare
  $W_D^F < 0$ tail.
- Two swap probabilities coexist: $P_{swap}$ (Metropolis, eq:10) drives actual
  RE moves; $p_{swap}$ (Fermi, eq:15) is used for analysis. Both give the same
  Boltzmann distribution of swapped/unswapped states (eq:13).
- Helmholtz vs Gibbs: the study is N-V-T (Helmholtz $\Delta F$); literature is
  N-P-T (Gibbs $\Delta G$); the two are treated as equal here because the box
  at $\lambda=0$ corresponds to 1 atm and P-V work is reversible.
