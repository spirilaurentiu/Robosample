# Notation - Kofke 2002

Reduced Lennard-Jones units throughout the simulations: LJ diameter $\sigma = 1$
and well depth $\epsilon = 1$. Energies, temperatures and densities are in LJ
units. $\beta = 1/kT$.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $\beta_i$ | reciprocal temperature of replica $i$, $1/kT_i$ | scalar, 1/energy | subscript 0 = cold, 1 = hot; $\beta_0 > \beta_1$ |
| $T_i$ | temperature of replica $i$ | scalar, LJ units | $T_0 < T_1$ |
| $U_i$ | potential energy of system $i$ before the swap | scalar, energy | evaluated pre-exchange |
| $k$ | Boltzmann constant | energy/temperature | often set to 1 in reduced units |
| $U_m$ | lowest possible energy | scalar, energy | $U_m = U_r - C/\beta_r$ in the constant-$C$ model |
| $\bar{p}_{\rm acc}$ | average swap-acceptance probability | dimensionless, [0,1] | fraction of all swap trials accepted |
| $p_i(U)$ | energy distribution at temperature $\beta_i$ | probability density, 1/energy | normalized: $\int p_i\,dU = 1$ |
| $\Omega(U)$ | density of states | 1/energy | independent of $\beta$ |
| $Q(\beta_i)$ | canonical partition function | dimensionless | $Q = \int \Omega(U) e^{-\beta_i U}dU$ |
| $S$ | entropy | energy/temperature | $S = k\ln\Omega$ (bridge equation) |
| $\Delta S$ | hot-minus-cold entropy difference | energy/temperature | $\Delta S > 0$; extensive in $N$ |
| $C$ | extensive constant-volume heat capacity in units of $k$ | dimensionless | $C \equiv C_V/k$ |
| $C_V$ | extensive constant-volume heat capacity | energy/temperature | $C_V = N c_V$ |
| $c_V$ | molar constant-volume heat capacity | per-molecule | |
| $N$ | number of molecules (system size) | integer | 32, 64, 108 in the study |
| $\rho$ | number density | LJ units | fixed at 0.80 |
| $B$ | temperature ratio $\beta_1/\beta_0 = T_0/T_1$ | dimensionless, $<1$ | asymptotic-formula variable |
| $\kappa$ | ratio integration variable (eq:sub, eq:10) | dimensionless | range $[0, \beta_1/\beta_0]$ |
| $\gamma$ | sum variable (eq:sub) | dimensionless | integrated out |
| $\Gamma$ | gamma function | - | |
| $U_r, T_r, \beta_r$ | arbitrary reference state | - | reference for the constant-$C$ density of states |
