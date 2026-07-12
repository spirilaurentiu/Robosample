# Equations - Wyczalkowski & Pappu (2008)

Convention: $\beta = (k_B T)^{-1}$. $\Gamma$ is a full system configuration.
$\lambda$ is the Kirkwood coupling parameter, $0 \le \lambda \le 1$. Subscripts
$0, 1$ denote the two Hamiltonians $U_0 = U(\cdot, \lambda_0)$ and
$U_1 = U(\cdot, \lambda_1)$. $\langle \cdot \rangle_i$ is an equilibrium average
over ensemble $i$; a superscript $*$ marks a finite-sample estimate.

<!-- eq:1 -->
$$ \exp(\beta W_D) = \frac{P_F(\beta W_D)}{P_R(-\beta W_D)}. $$
- **what:** Crooks fluctuation theorem: ratio of forward to reverse
  dissipated-work distributions equals $\exp(\beta W_D)$.
- **symbols:** $W_D$ - dissipated work value (energy); $P_F$ - forward
  dissipated-work distribution; $P_R$ - reverse dissipated-work distribution;
  $\beta = 1/(k_B T)$.

<!-- eq:2 -->
$$ F_i = -\beta^{-1} \ln \left\{ \int d\Gamma \exp[-\beta U_i(\Gamma)] \right\}. $$
- **what:** Helmholtz free energy of replica $i$ from its configurational
  partition function.
- **symbols:** $F_i$ - free energy of replica $i$; $U_i(\Gamma) = U(\Gamma,
  \lambda_i)$ - potential of replica $i$; $\Gamma$ - configuration.

<!-- eq:3 -->
$$ \rho_i(\Gamma) = \exp\{\beta[F_i - U_i(\Gamma)]\}. $$
- **what:** normalized equilibrium (Boltzmann) probability density of
  configuration $\Gamma$ in ensemble $i$.
- **symbols:** $\rho_i$ - equilibrium density of ensemble $i$; $F_i$ - free
  energy (eq:2); $U_i$ - potential.

<!-- eq:4a -->
$$ W^{F}(\Gamma) = U_{1}(\Gamma) - U_{0}(\Gamma). $$
- **what:** forward switching work for a configuration $\Gamma$ (instantaneous
  $\lambda_0 \to \lambda_1$).
- **symbols:** $W^F$ - forward work; $U_0, U_1$ - potentials at the two
  $\lambda$ endpoints.

<!-- eq:4b -->
$$ W^{R}(\Gamma) = U_{0}(\Gamma) - U_{1}(\Gamma). $$
- **what:** reverse switching work ($\lambda_1 \to \lambda_0$);
  $W^R = -W^F$.
- **symbols:** $W^R$ - reverse work.

<!-- eq:5a -->
$$ \delta F_{FEP}^F = -\beta^{-1} \ln \langle \exp(-\beta W^F) \rangle_0. $$
- **what:** forward free-energy perturbation (FEP) estimator of $\delta F$,
  averaging forward work over the $U_0$ ensemble.
- **symbols:** $\delta F_{FEP}^F$ - forward FEP estimate; $\langle \cdot
  \rangle_0$ - average over the equilibrium ensemble of $U_0$.

<!-- eq:5b -->
$$ \delta F_{FEP}^{R} = +\beta^{-1} \ln \langle \exp(-\beta W^{R}) \rangle_{1}. $$
- **what:** reverse FEP estimator of $\delta F$ (same process
  $\lambda_0 \to \lambda_1$), averaging reverse work over the $U_1$ ensemble.
- **symbols:** $\delta F_{FEP}^R$ - reverse FEP estimate; $\langle \cdot
  \rangle_1$ - average over the equilibrium ensemble of $U_1$. Note the leading
  sign is $+\beta^{-1}$.

<!-- eq:6 -->
$$ \epsilon_H \equiv \delta F_{FEP}^F - \delta F_{FEP}^R. $$
- **what:** hysteresis error: difference between forward and reverse FEP
  estimates; sampling-quality measure between two replicas.
- **symbols:** $\epsilon_H$ - hysteresis error (energy).

<!-- eq:7a -->
$$ W_D^F(\Gamma) = W^F(\Gamma) - \delta F. $$
- **what:** forward dissipated work (work minus true free-energy change).
- **symbols:** $W_D^F$ - forward dissipated work; $\delta F$ - true free-energy
  change $\lambda_0 \to \lambda_1$.

<!-- eq:7b -->
$$ W_D^R(\Gamma) = W^R(\Gamma) + \delta F. $$
- **what:** reverse dissipated work; note $W_D^F = -W_D^R$ for a given $\Gamma$.
- **symbols:** $W_D^R$ - reverse dissipated work.

<!-- eq:8 -->
$$ \exp[\beta W_D + \beta \epsilon_{FT}^*(W_D)] = \frac{P_F^*(\beta W_D)}{P_R^*(-\beta W_D)}. $$
- **what:** Crooks theorem with a finite-sampling error term; recovers eq:1 when
  $\epsilon_{FT}^* = 0$.
- **symbols:** $\epsilon_{FT}^*(W_D)$ - fluctuation error (function of $W_D$);
  $P_F^*, P_R^*$ - observed (finite-sample) dissipated-work distributions.

<!-- eq:9 -->
$$ \epsilon_H = -\beta^{-1} \ln \langle \exp(-\beta \epsilon_{FT}^*) \rangle_0^*. $$
- **what:** relates hysteresis error to the fluctuation error; smaller deviation
  from Crooks means smaller $\epsilon_H$.
- **symbols:** $\langle \cdot \rangle_0^*$ - finite-simulation average over the
  $U_0$ ensemble.

<!-- eq:10 -->
$$ P_{swap} = \min[1, \exp(-\beta \Delta U_{swap})]. $$
- **what:** Metropolis acceptance probability for a Hamiltonian replica exchange
  swap.
- **symbols:** $P_{swap}$ - Metropolis swap acceptance probability;
  $\Delta U_{swap}$ - swap energy change (eq:11).

<!-- eq:11a -->
$$ \Delta U_{swap} = U_0(\Gamma_1) + U_1(\Gamma_0) - U_0(\Gamma_0) - U_1(\Gamma_1). $$
- **what:** energy change on swapping configurations $\Gamma_0 \leftrightarrow
  \Gamma_1$ between replicas 0 and 1.
- **symbols:** $\Gamma_0$ - config drawn from $U_0$ ensemble; $\Gamma_1$ -
  config drawn from $U_1$ ensemble.

<!-- eq:11b -->
$$ \Delta U_{swap} = W^F + W^R. $$
- **what:** swap energy equals sum of forward and reverse work
  ($W^F$ evaluated at $\Gamma_0$, $W^R$ at $\Gamma_1$).
- **symbols:** $W^F = W^F(\Gamma_0)$; $W^R = W^R(\Gamma_1)$.

<!-- eq:11c -->
$$ \Delta U_{swap} = W_D^F + W_D^R. $$
- **what:** swap energy equals sum of forward and reverse dissipated work (the
  $\delta F$ terms cancel).
- **symbols:** $W_D^F = W_D^F(\Gamma_0)$; $W_D^R = W_D^R(\Gamma_1)$.

<!-- eq:12a -->
$$ \rho_N(\gamma) = \rho_0(\Gamma_0)\rho_1(\Gamma_1). $$
- **what:** native joint probability of a configuration pair
  $\gamma = (\Gamma_0, \Gamma_1)$.
- **symbols:** $\rho_N$ - native joint density; $\gamma = (\Gamma_0, \Gamma_1)$.

<!-- eq:12b -->
$$ \rho_N'(\gamma) = \rho_0(\Gamma_1)\rho_1(\Gamma_0) = \rho_N(\gamma'). $$
- **what:** swapped joint probability; equals the native probability of the
  swapped pair $\gamma' = (\Gamma_1, \Gamma_0)$.
- **symbols:** $\rho_N'$ - swapped joint density; $\gamma' = (\Gamma_1,
  \Gamma_0)$.

<!-- eq:13 -->
$$ \frac{\rho_N'}{\rho_N} = \exp(-\beta \Delta U_{swap}). $$
- **what:** interreplica equilibrium relationship: relative probability of
  swapped vs native configurations.
- **symbols:** derived from eq:12, eq:3, eq:11a.

<!-- eq:14 -->
$$ \epsilon_H \simeq -\beta^{-1} \int d\Gamma_0 d\Gamma_1\, \rho_\epsilon. $$
- **what:** hysteresis error equals ($-\beta^{-1}$ times) the integrated error
  in the estimated swapped distribution; $\epsilon_H \to 0$ as $\rho_\epsilon
  \to 0$.
- **symbols:** $\rho_\epsilon(\Gamma_0, \Gamma_1)$ - sampling error in the
  swapped distribution $\rho_N'$.

<!-- eq:15 -->
$$ f(x) = \frac{1}{1 + \exp(x)}. $$
- **what:** Fermi (logistic) function; Fermi swap probability is
  $p_{swap} = f(\beta \Delta U_{swap})$.
- **symbols:** $f$ - Fermi function; $p_{swap}$ - Fermi swap probability
  (theory), distinct from the Metropolis $P_{swap}$ (eq:10).

<!-- eq:16a -->
$$ \langle p_{swap} \rangle \equiv \langle \langle f(\beta \Delta U_{swap}) \rangle_0 \rangle_1. $$
- **what:** average Fermi swap probability as a double equilibrium average over
  both independent ensembles.
- **symbols:** $\langle p_{swap} \rangle$ - average swap probability.

<!-- eq:16b -->
$$ \langle p_{swap} \rangle = \int d\Gamma_0 d\Gamma_1\, \rho_N\, f(\beta \Delta U_{swap}). $$
- **what:** same average written as an integral over the native joint density.
- **symbols:** $\rho_N$ - native joint density (eq:12a).

<!-- eq:17a -->
$$ \langle p_{swap} \rangle = \left\langle \left\langle \frac{\rho_N'}{\rho_N + \rho_N'} \right\rangle_0 \right\rangle_1. $$
- **what:** average swap probability in terms of native/swapped densities.
- **symbols:** $\rho_N, \rho_N'$ - native and swapped joint densities.

<!-- eq:17b -->
$$ \langle p_{swap} \rangle = \int d\Gamma_0 d\Gamma_1\, \frac{\rho_N \rho_N'}{\rho_N + \rho_N'}. $$
- **what:** average swap probability equals the overlap integral of $\rho_N$ and
  $\rho_N'$; large value means large ensemble overlap.
- **symbols:** integrand is a normalized probability of a configuration pair.

<!-- eq:18 -->
$$ \langle p_{swap} \rangle \simeq \frac{1}{2} - \frac{\beta^2 \delta_{\lambda}^2}{4} C_{\lambda}, \qquad C_{\lambda} \equiv \operatorname{var}\left(\frac{\partial U}{\partial \lambda}\right) = \langle (\partial U/\partial \lambda)^2 \rangle_0 - \langle \partial U/\partial \lambda \rangle_0^2. $$
- **what:** linearized average swap probability for small $\lambda$ separation;
  $C_\lambda$ sets the decline rate. At $\delta_\lambda=0$, $\langle p_{swap}
  \rangle = 1/2$.
- **symbols:** $\delta_\lambda$ - difference in $\lambda$ between neighboring
  replicas; $C_\lambda$ - variance of $\partial U/\partial \lambda$ in the
  $U_0$ ensemble.

<!-- eq:19 -->
$$ \langle \exp(-\beta W_D^F) \rangle_0 = 1. $$
- **what:** exact identity (from eq:5a); requires sampling rare $W_D^F < 0$
  configurations for convergence.
- **symbols:** $\langle \cdot \rangle_0$ - average over $U_0$ equilibrium
  ensemble.

<!-- eq:20a -->
$$ \frac{\rho_0(\Gamma_0)}{\rho_1(\Gamma_0)} = \exp[\beta W_D^F(\Gamma_0)]. $$
- **what:** density ratio at a configuration equals the exponential of forward
  dissipated work.
- **symbols:** $\Gamma_0$ - configuration from the $U_0$ ensemble.

<!-- eq:20b -->
$$ \frac{\rho_1(\Gamma_1)}{\rho_0(\Gamma_1)} = \exp[\beta W_D^R(\Gamma_1)]. $$
- **what:** reverse density ratio equals the exponential of reverse dissipated
  work.
- **symbols:** $\Gamma_1$ - configuration from the $U_1$ ensemble.

<!-- eq:dF_sum -->
$$ \Delta F \equiv \sum_{i}^{M-1} (\delta F)_{i}. $$
- **what:** total hydration free energy as the sum of neighboring pairwise
  free-energy changes (each from BAR).
- **symbols:** $\Delta F$ - total free energy of hydration; $M$ - total number
  of replicas; $(\delta F)_i$ - free-energy change between replicas $i$ and
  $i+1$.

<!-- eq:eps_rms -->
$$ \epsilon_{\rm rms} \equiv \sqrt{\sum_{i}^{M-1} (\epsilon_H)_i^2 / M}. $$
- **what:** root-mean-square hysteresis error over the $\lambda$ schedule.
- **symbols:** $(\epsilon_H)_i$ - hysteresis error between replicas $i$ and
  $i+1$; $M$ - number of replicas.

<!-- eq:B1 -->
$$ U_C(r, \lambda_C) = \lambda_C \frac{q_i q_j}{\alpha_C(1 - \lambda_C) + r}. $$
- **what:** soft-core scaled Coulomb potential between two atoms; the
  $\alpha_C(1-\lambda_C)$ term imposes a minimum effective separation at small
  $\lambda_C$.
- **symbols:** $r$ - interatomic distance; $q_i, q_j$ - atomic charges;
  $\lambda_C$ - Coulomb coupling; $\alpha_C = 1.5$ Angstrom - soft-core
  parameter.

<!-- eq:B2 -->
$$ U_{LJ}(r, \lambda_{LJ}) = B\,A\,(A - 1). $$
- **what:** general Lennard-Jones form; unscaled recovers the 12-6 potential
  with $A = (\sigma/r)^6$, $B = 4\epsilon$.
- **symbols:** $A$ - dimensionless radial factor (eq:B3a); $B$ - well-depth
  prefactor (eq:B3b); $\sigma, \epsilon$ - LJ parameters.

<!-- eq:B3a -->
$$ A(r, \lambda_{LJ}) = 1 \left/ \left[ \alpha_{LJ} (1 - \lambda_{LJ})^b + \left( \frac{r}{\sigma} \right)^6 \right] \right. . $$
- **what:** exponential soft-core radial factor; softens the $r^{-6}$ singularity
  at small $\lambda_{LJ}$.
- **symbols:** $\alpha_{LJ} = 0.5$ Angstrom - soft-core parameter; $b = 4$ -
  soft-core exponent (called $a$ in the prose, $a=4$); $\sigma$ - LJ diameter.
  <!-- CHECK: source writes the exponent as "b" in eq B3a but states "a=4" in prose; treat b=a=4. -->

<!-- eq:B3b -->
$$ B(\lambda_{LJ}) = 4\epsilon\, \frac{1 - e^{-k\lambda_{LJ}}}{1 - e^{-k}}. $$
- **what:** exponential scaling of the LJ well depth; $B \to 0$ at
  $\lambda_{LJ}=0$ and $B \to 4\epsilon$ at $\lambda_{LJ}=1$.
- **symbols:** $\epsilon$ - LJ well depth; $k = 1$ - scaling constant.

## Derivations note

Appendix A (A1-A4, and the intermediate FEP/partition-function expansion steps)
contains proof algebra only; those steps are kept in `paper.md` under the
`Derivation (not implemented)` heading and are not re-emitted here. Only the
resulting closed forms (eq:1, eq:9, eq:14, eq:18) are implementable.
