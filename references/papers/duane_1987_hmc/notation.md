# Notation - Duane et al. 1987, Hybrid Monte Carlo

Conventions: natural/reduced units throughout. The Boltzmann factor is written $\exp[-S(\phi)]$ with $S$ dimensionless (temperature absorbed into $S$; equivalently $\beta = 1$). Momenta have **unit mass** (kinetic term $\tfrac12\pi^2$). Field $\phi$ plays the role of "position", momentum $\pi$ is refreshed each iteration from a unit Gaussian. For a molecular port: $S(\phi) \to \beta U(x)$, $\phi \to$ coordinates, $\pi \to$ momenta with mass matrix (here identity).

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $\phi$ | field configuration ("position") | real array over lattice sites | one component per d.o.f. |
| $\pi$ | conjugate momentum | real array, same shape as $\phi$ | mass = 1; $\tfrac12\pi^2$ kinetic |
| $S(\phi)$ | action (dimensionless energy) | scalar | Boltzmann weight $\exp[-S]$, so $\beta$ absorbed |
| $S'(\phi)$ | guidance action | scalar | may differ from $S$ (guidance $\neq$ acceptance) |
| $H(\phi,\pi)$ | acceptance Hamiltonian | scalar | $=\tfrac12\pi^2 + S(\phi)$; defines equilibrium dist |
| $H'(\phi,\pi)$ | guidance Hamiltonian | scalar | $=\tfrac12\pi^2 + S'(\phi)$; drives EOM |
| $\delta H$ | energy change (see eq:12) | scalar | $\delta H = H_{\text{old}} - H_{\text{new}}$; accept $\min(1,e^{\delta H})$ |
| $\tau$ | fictitious "computer" / MD time | scalar | trajectory parameter, not physical time |
| $\tau_0$ | total trajectory length | scalar | $= n\,\delta\tau$ |
| $\delta\tau$ | leapfrog step size | scalar | area-preserving for any $\delta\tau$; "safe" hybrid $\lesssim0.01$-$0.02$ |
| $n$ | number of MD (leapfrog) steps | integer | $n = \tau_0/\delta\tau$ |
| $Z$ | partition function | scalar | normalization of $P_S$ |
| $P_S(\phi)$ | target distribution | prob. density | $=(1/Z)\exp[-S(\phi)]$ |
| $P_G(\pi)$ | momentum heatbath dist | prob. density | $\propto\exp(-\pi^2/2)$, mean 0, unit variance |
| $P_C$ | candidate/proposal probability | prob. | generic in eq:6 |
| $P_H$ | phase-space proposal kernel | prob. | delta fn (deterministic), must be reversible |
| $P_A$ | acceptance probability | prob. in $[0,1]$ | Metropolis, eqs 6, 12 |
| $P_M$ | marginal field transition | prob. | eq:13 |
| $\Omega$ | observable / operator | scalar | $\langle\Omega\rangle$ expectation |
| $T$ | number of MC samples | integer | estimator error $\mathcal{O}(1/\sqrt T)$ |
| $\mathcal{M}$ | fermion matrix | complex matrix | pseudofermion kernel; $(\mathcal{M}^\dagger\mathcal{M})^{-1}$ via CG |
| $\chi,\chi^*$ | pseudofermion (bosonic) fields | complex array | refreshed $\chi=\mathcal{M}^\dagger\eta$, $\eta$ Gaussian; fixed during MD |
| $\bar\psi,\psi$ | Grassmann fermion fields | Grassmann | integrated out to $\chi$ |
| $\beta$ | lattice coupling constant | scalar | test QED; distinct from inverse temperature |
| $U_{\mu\nu}(x)$ | plaquette (link product) | U(1)/SU(N) element | Wilson gauge action |
| $m$ | fermion (electron) mass | lattice units | in $(-D^2+m^2)^{-1}$ |
| $L$ | linear lattice size | integer | volume $= L^4$ |
| $a$ | acceptance rate | fraction in $[0,1]$ | effective step size $= a\,\delta\tau$ |
