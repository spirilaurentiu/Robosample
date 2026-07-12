# Notation - Nilmeier 2011 NCMC

Reduced (dimensionless) units throughout: energies are expressed as reduced potentials
`u = beta * (energy terms)`, so `beta` is absorbed where noted. The Hamiltonian is presumed
invariant under momentum inversion (`pi(x,lambda) = pi(tilde x, lambda)`).

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| `x` | microstate | point in phase space Gamma; may include coords, momenta, box dims | full state, not just coords |
| `tilde x` | momentum-inverted microstate | same as x | velocities negated; coords unchanged |
| `Gamma` | phase space | domain of integration | - |
| `lambda` | thermodynamic parameters | e.g. {beta, H, p} | defines a thermodynamic state |
| `pi_lambda(x)` | Boltzmann density at state lambda | probability density | Eq. 1 |
| `pi(x,lambda)` | expanded-ensemble joint density | probability density | Eq. 3; reduces to pi_lambda if single state |
| `u_lambda(x)` | reduced potential | dimensionless (already x beta) | ensemble-dependent, Eq. 2 |
| `Z_lambda` | partition function | dimensionless | normalizer of pi_lambda |
| `beta` | inverse temperature | 1/energy = 1/(k_B T) | reduced units |
| `H(x)` | Hamiltonian | energy | may include external biasing potential |
| `V(x)` | system volume | volume | NpT ensemble |
| `p` | external pressure | pressure | NpT ensemble |
| `omega_lambda` | expanded-ensemble weight of state lambda | > 0, dimensionless | externally imposed |
| `G` | set of allowed lambda values | discrete or continuous | single value => single state |
| `Lambda` | protocol | ordered list {alpha_1,K_1,...,alpha_T,K_T} | forward protocol |
| `tilde Lambda` | reverse protocol | {K_T,alpha_T,...,K_0,alpha_0} | order reversed |
| `P(Lambda|x_0,lambda_0)` | protocol selection probability | probability | depends only on initial state |
| `alpha_t(x,y)` | perturbation kernel at step t | conditional prob. density | drives selected DOF |
| `K_t(x,y)` | propagation kernel at step t | conditional prob. density | MCMC or MD at fixed lambda_t |
| `T` | number of protocol / switching steps | integer | switching length (free parameter) |
| `x_t^*` | perturbed config (after alpha_t) | microstate | pre-propagation |
| `x_t` | config after propagation (after K_t) | microstate | end of step t |
| `X` | forward trajectory (x_0,...,x_T) | sequence | - |
| `tilde X` | time-reversed trajectory | sequence | momenta inverted |
| `A(X|Lambda)` | acceptance probability of trajectory X | in [0,1] | Eq. 12 |
| `Pi(X|x_0,Lambda)` | pathwise generation probability | probability | product of kernels, Eqs. 7-8 |
| `rho(X,Lambda|x_0,lambda_0)` | joint traj+protocol probability | probability | = Pi * P |
| `Delta u(X|Lambda)` | reduced energy difference | dimensionless | u_T(x_T) - u_0(x_0) |
| `Delta S(X|Lambda)` | conditional path-action difference | dimensionless | Eq. 10; = 0 for symplectic; = -q for reversible MCMC |
| `w(X|Lambda)` | nonequilibrium work | dimensionless (reduced) | sum over perturbation steps, Eq. 16 |
| `q(X|Lambda)` | heat | dimensionless (reduced) | sum over propagation steps, Eq. 17; w+q=Delta u |
| `phi, phi'` | torsion angle (current, proposed) | radians | von Mises proposal Eq. 13 |
| `kappa` | von Mises concentration / force constant | dimensionless, > 0 | Eq. 13 |
| `I_0(kappa)` | modified Bessel function order 0 | - | normalizer of von Mises |
| `J(phi)`, `J_r`, `J(x)` | Jacobian of coordinate change | dimensionless | = 1 for torsion rotation; = (r_new/r_old)^2 for radial |
| `s` (box) | box scaling factor | dimensionless | [(V+dV)/V]^(1/3); ratio s^(3N) |
| `N` (box) | number of molecular centers | integer | box-scaling Jacobian exponent |
| `r` | dimer interparticle distance / extension | length (Angstrom) | order parameter |
| `r_0` | compact-minimum distance | length | = r_WCA |
| `r_WCA` | WCA cutoff (LJ minimum) | length | = 2^(1/6) * sigma |
| `h` | double-well barrier height | energy | = 5 k_B T |
| `s` (bond) | double-well width parameter | length | = r_WCA / 2 |
| `sigma` | LJ / WCA diameter | length | = 3.4 Angstrom |
| `epsilon` | LJ / WCA well depth | energy | = 120 k_B T |
| `m` | particle mass | mass | = 39.9 amu (argon-like) |
| `rho sigma^3` | reduced number density | dimensionless | = 0.96 (solvated) |
| `k_B T / epsilon` | reduced temperature | dimensionless | = 0.824 |
| `tau` (time unit) | reduced time unit | time | = sqrt(sigma^2 m / epsilon) |
| `Delta t` | integration timestep | time | = 0.002 tau |
| `gamma` (friction) | collision frequency / friction | 1/time | = tau^-1 (GHMC); Langevin friction |
| `F_t(x)` | systematic force | force | = -dH_t/dx |
| `xi_t, xi_t'` | Gaussian noise variates | mean 0, variance beta^-1 | one (Ermak-Yeh) or two (BBK) per DOF per step |
| `tilde xi_t` | reverse noise variate | same | generates reverse transition |
| `v_t, v_t^*, v_t'` | velocities (post, pre, auxiliary) | velocity | BBK integrator |
| `tau` (correlation) | integrated autocorrelation time | iterations | Eq. 21 |
| `g` | statistical inefficiency | iterations | = 1 + 2 tau |
| `gamma` (acceptance) | average NCMC acceptance probability | in [0,1] | symmetric between states |
| `tau_MD, tau_NCMC, tau_eff` | correlation times | iterations | Eqs. 22, 54 |
| `E` | efficiency gain vs MD alone | dimensionless | Eq. 23 |
| `mu` | nonunit eigenvalue of 2x2 transition matrix | in (-1,1) | tau = -1/ln(mu) |
| `K` (umbrella) | umbrella force constant | energy/length^2 | = k_B T / eta^2 |
| `eta` | umbrella width | length | = 0.3 Angstrom |
| `theta` | Heaviside step function | - | 1 for arg >= 0, else 0 |
