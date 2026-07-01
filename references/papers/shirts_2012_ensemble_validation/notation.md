# Notation - Shirts 2012, Ensemble validation

Reduced/unspecified units throughout the toy models (k_B = 1 unless a physical
system is stated). For molecular systems energies are kJ/mol, temperatures K,
pressures bar, volumes nm^3.

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| beta | inverse temperature | 1/(k_B T); reduced (k_B=1) in toys, else (kJ/mol)^-1 | beta = 1/(k_B T) |
| T | temperature | K (or reduced) | |
| k_B | Boltzmann constant | e.g. 0.996072 kJ/mol per argon reduced unit; 1 in toys | |
| E | total energy (extensive) | kJ/mol or reduced | the fit variable in canonical test |
| E_kin, E_pot | kinetic / potential energy | same as E | separable in NVT |
| H | Hamiltonian / instantaneous enthalpy E+PV | energy units | in NPT tests H is the fit variable |
| P | pressure | bar (or reduced) | |
| V | volume | nm^3 (or reduced) | fit variable in volume test |
| N, N_i | particle number(s) | integer(s) | fit variable in grand canonical |
| mu, mu_i | chemical potential(s) | energy | -mu plays role of P in GC |
| Omega(E) | density of states | = exp(S/k_B); function of E only, NOT beta | cancels in ratio |
| Omega'(H,P) | DOS at fixed H=E+PV | function of P, not beta | |
| Q(beta) | canonical partition function | integral Omega(E)exp(-beta E)dE; function of beta only | |
| Q_kin, Q_pot | kinetic / potential partition functions | Q_kin = prod_i (m_i/(pi beta))^{3/2} | mass cancels in ratio |
| Delta(beta,P) | isothermal-isobaric partition function | | Delta = Q_kin * Delta_pot |
| Xi(beta,mu) | grand partition function | | |
| A | Helmholtz free energy | = -beta^-1 ln Q | absolute value arbitrary; only DeltaA meaningful |
| G | Gibbs free energy (free enthalpy) | replaces A in NPT | set G_1+G_2=0 |
| S | entropy | Omega=exp(S/k_B) | |
| alpha_0 | fit intercept | = beta_2 A_2 - beta_1 A_1 = ln(Q_1/Q_2) | contains free energies |
| alpha_1 | fit slope | = -(beta_2 - beta_1) | KNOWN from user temperatures |
| alpha (vector) | fit params | length M+1; alpha.X := alpha_0 + sum_j alpha_j X_j | NOT plain dot product |
| Delta_beta | beta_2 - beta_1 | | slope of canonical/enthalpy test |
| Delta_P | P_1 - P_2 (eq 15) | bar | sign convention as printed |
| Delta_T, Delta_G | T and G differences | | |
| beta_ave | (beta_1 + beta_2)/2 | fixed at user average | energy-scale choice |
| P_ave | (P_1 + P_2)/2 | fixed at user average | |
| f(x) | Fermi / logistic function | eq:8 form [1+exp(-x)]^-1; appendix form [1+exp(x)]^-1 | args flipped so consistent |
| D | dimensionality of HO | integer | use D=20, not D=1 (atypical DOS) |
| K | harmonic spring constant | ; K=(a/V)^2 in NPT toy | |
| x_{i,0} | HO equilibrium position | set to 0 in tests | |
| a | constant in NPT-toy spring K=(a/V)^2 | | |
| 3N | kinetic degrees of freedom | integer | replace with actual DOF if constrained / COM removed |
| C_V, C_P | heat capacities | kJ/(mol K) | kinetic contribution 3Nk_B/2 (ideal gas) |
| sigma_E, sigma_V | std of energy / volume distribution | = T sqrt(C_V k_B), etc. | used for gap selection |
| kappa_T | isothermal compressibility | = -(1/V)(dV/dP)_T | |
| p_{k,i} | empirical bin probability | = n_{k,i}/N_i | var(p_k)=p_k(1-p_k)/N |
| n_{k,i}, N_i | bin counts / total samples | integers | sim i in {1,2} |
| r_k | histogram ratio in bin k | = p_{k,2}/p_{k,1} | |
| X | design matrix | (M+1) x N; col1 ones | linear fit |
| W | weight matrix | diagonal, W_ii = 1/variance_i | inverse-variance weights |
| J | Jacobian of model wrt alpha | at minimum | nonlinear fit |
| M (appendix) | ln(N_1/N_2) offset | scalar | logistic intercept offset |
| nu | noise amplitude added to toy energies | dE = nu*|N(0,1)| | error-injection test |
| sigma deviation | # of std devs of fitted slope from true | dimensionless | flag if consistently >2-3 sigma |

## Conventions / gotchas
- **Slope sign:** canonical/enthalpy log-ratio P2/P1 has slope -(beta_2-beta_1)=alpha_1 on E (or H). Volume test slope is -beta(P_2-P_1) on V.
- **Two free parameters only:** although A_1,A_2,beta_1,beta_2 appear, only DeltaA and Delta_beta are free (fix beta_ave and G/A sum).
- **DOF:** always subtract removed COM DOF and constraint DOF from 3N before using kinetic-energy formulas.
- **Kinetic-energy estimator matters:** leapfrog half-step-averaged KE is more accurate than velocity-Verlet full-step KE; the latter shows larger deviations at long time steps.
- **Choose T_1 != T_2:** if beta_1=beta_2 no information; if too far apart, tail small-sample noise dominates. Optimal gap ~ sum of the two distribution std devs (~2 sigma), 1.5-2x for molecular systems.
