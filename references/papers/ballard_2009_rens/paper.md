# Replica Exchange with Nonequilibrium Switches (RENS)

Andrew J. Ballard, Christopher Jarzynski. PNAS 2009, 106(30). DOI: 10.1073/pnas.0900406106.

## Abstract

A replica exchange (parallel tempering) method in which attempted configuration
swaps are generated using nonequilibrium work simulations. By effectively
increasing phase space overlap, this approach mitigates the need for many
replicas. It achieves the computational efficiency of ordinary replica exchange
using fewer replicas.

Keywords: parallel tempering, nonequilibrium work simulations, molecular
dynamics, Monte Carlo, rough energy landscapes.

## Introduction

Replica exchange (parallel tempering) simulates M independent copies of the
system, typically ordered by increasing temperature, and performs "swaps" in
which an attempted exchange of configurations between adjacent replicas is
accepted or rejected according to a Metropolis-like scheme. Low-temperature
replicas gain access to the broader phase space explored by high-temperature
replicas.

The weakness of replica exchange is its scaling with system size N: the number
of replicas needed typically grows as $N^{1/2}$. This scaling is rooted in a
phase space overlap requirement. To achieve a reasonable frequency of accepted
swaps, neighboring replicas should overlap in phase space; with increasing
system size more replicas are needed. The difficulty is acute for large
molecules in explicit solvent, where poor overlap is due mostly to the large
number of solvent molecules.

RENS uses nonequilibrium simulations to increase overlap between replicas. When
attempting a swap between replicas A and B, first perform a finite-time "work
simulation" in replica A, dragging the system toward the region of phase space
characteristic of equilibrium ensemble B, and vice versa in replica B. Then
attempt to swap the generated configurations using a work-based acceptance
criterion (Eq. 2). The acceptance/rejection step ensures detailed balance is
preserved and the replicas' equilibrium states are undisturbed. The CPU time
devoted to the work simulations is an added cost, but the return is an increased
acceptance probability.

The paper uses REM for the usual replica exchange method and RENS for the method
based on nonequilibrium work simulations.

## Description of Method

Let $\mathcal{R}_1, \mathcal{R}_2, \ldots \mathcal{R}_M$ denote the collection of
replicas, and let $H_i(x)$ and $T_i$ denote the Hamiltonian (energy function)
and temperature of the canonical ensemble simulated in $\mathcal{R}_i$. Here $x$
denotes a point in configuration space or phase space. Often the $H_i$'s are
identical and only the temperatures differ, or vice versa, but this need not be
assumed. Define a reduced Hamiltonian $h_i(x) = H_i(x)/k_B T_i$, so that the
equilibrium distribution in $\mathcal{R}_i$ takes the form
$p_i^{\text{eq}} \propto \exp(-h_i)$.

In ordinary REM, if 2 replicas $\mathcal{R}_A$ and $\mathcal{R}_B$ are found in
configurations $x$ and $y$ at the time of an attempted swap, then the swap
$x \leftrightarrow y$ is accepted with probability
$P_{\text{acc}} = \min\{1, e^{-\Delta h}\}$, where $\Delta h$ is defined in
Eq. 1. If there is little overlap between $p_A^{\text{eq}}$ and
$p_B^{\text{eq}}$, then typically $P_{\text{acc}} \ll 1$ and the swap is most
likely rejected.

RENS proceeds instead as follows. The interval from $t_0$ to $t_1$ corresponds
to independent equilibrium sampling in each of the M replicas, using the reduced
Hamiltonians $h_1, \ldots h_M$. At time $t_1$, a swap between replicas A and B is
attempted. In lieu of an instantaneous exchange, first perform a pair of work
simulations. In $\mathcal{R}_A$ the system evolves from $t_1$ to $t_2$ as the
reduced Hamiltonian is parametrically switched from $h_A$ to $h_B$ (Eq. 3.1).
Let $x'$ denote the final configuration and $w_A$ the reduced work (Eq. 4.1). In
$\mathcal{R}_B$ simulate the reverse process, switching from $h_B$ to $h_A$, and
define $y'$ and $w_B$ analogously. Then attempt a swap with acceptance
probability given by Eq. 2, where $w = w_A + w_B$. If accepted, $y'$ is copied
into replica A and $x'$ into replica B; if rejected, the configurations revert
to $x$ and $y$. Subsequently equilibrium sampling continues.

Each replica alternates between sampling intervals at fixed $h_i$ and work
intervals. RENS satisfies detailed balance in the following sense: in each
replica $\mathcal{R}_i$, if the data generated during work intervals is discarded
and the remaining sampling intervals are stitched together, the result is a long
trajectory that samples the distribution $p_i^{\text{eq}}$. The acceptance
criterion compensates for the fact that the system is driven out of equilibrium
during the work simulations.

The method traces its validity to Crooks' extension of detailed balance to
nonequilibrium trajectories. The definition of reduced work depends on the
dynamics chosen. For discrete-time Monte Carlo dynamics, RENS is equivalent to
the annealed swapping method of Opps and Schofield and closely related to the
cool-walking algorithm of Brown and Head-Gordon. The derivation below is for
deterministic, reversible molecular dynamics; the SI Appendix extends it to
stochastic evolution.

## Derivation (not implemented)

For a pair of replicas $\mathcal{R}_A$ and $\mathcal{R}_B$, introduce a
parametrized Hamiltonian $h(x; \lambda)$ interpolating from $h(x;0) = h_A(x)$ to
$h(x;1) = h_B(x)$. Specify a switching protocol $\lambda_t^A$, with
$\lambda_0^A = 0$ and $\lambda_\tau^A = 1$. In $\mathcal{R}_A$, starting from
$x_0 = x$, generate a trajectory $\gamma_A$ during which the system evolves as
$\lambda$ is varied from 0 to 1:

$$ \gamma_A: x = x_0 \xrightarrow{\lambda \to 1} x_{\tau} = x'. $$

Simultaneously, in $\mathcal{R}_B$ generate a trajectory $\gamma_B$, varying
$\lambda$ from 1 to 0 using the time-reversed protocol
$\lambda_t^B = \lambda_{\tau-t}^A$:

$$ \gamma_B : y' = y_\tau \xleftarrow{0 \leftarrow \lambda} y_0 = y. $$

These trajectories are generated by deterministic equations of motion symmetric
under time-reversal: for any trajectory $\gamma_A = (x_0 \to x_\tau)$ evolving
under $\lambda_t^A$, the time-reversed trajectory
$\bar{\gamma}_B = (\bar{x}_0 \leftarrow \bar{x}_\tau)$ evolves under
$\lambda_t^B$, where $\bar{x}$ denotes inversion of momenta,
$\mathbf{p} \to -\mathbf{p}$. This assumption holds for Hamiltonian dynamics,
Nose-Hoover dynamics, and others, provided the Hamiltonian is time-reversal
symmetric, $h(x; \lambda) = h(\bar{x}; \lambda)$.

The reduced work is defined in Eqs. 4.1 and 4.2, with Jacobians
$J_A = |\partial x_\tau / \partial x_0|$ and
$J_B = |\partial y_\tau / \partial y_0|$ associated with propagating the system.
Eq. 4 is analogous to the first law of thermodynamics, with $\ln J$ representing
a heat term associated with increase of system entropy.

To analyze the method, think of an expanded phase space containing 2 copies of
the system. We wish to sample $p_{AB}^{\text{eq}}(x,y) \propto e^{-h_A(x)-h_B(y)}$.
The work simulations plus attempted swap represent a trial Monte Carlo move
$(x,y) \to (y',x')$. Let $P(y',x'|x,y)$ be the corresponding transition
probability. To establish detailed balance, show that the net probability of the
transition $(x,y) \to (y',x')$ equals that of the reverse transition (Eq. 5).

Treat the final microstate as a function of the initial microstate,
$x_\tau = M_A(x_0)$, $y_\tau = M_B(y_0)$ (Eq. 6), obtained by integrating the
deterministic equations of motion. The work-simulation arrival probabilities are
Dirac deltas $\pi_A(x'|x) = \delta(x' - M_A(x))$ (Eq. 7), with joint probability
$\pi(x',y'|x,y) = \pi_A(x'|x)\pi_B(y'|y)$ (Eq. 8) and swap acceptance
$\alpha(x',y'|x,y) = \min\{1, e^{-w_A - w_B}\}$ (Eq. 9). The transition
probability is $P = \pi\alpha$.

Time-reversal symmetry relates $M_A$ and $M_B$: if $x' = M_A(x)$ then
$\bar{x} = M_B(\bar{x}')$, implying $\pi_A(x'|x) = \pi_B(\bar{x}|\bar{x}')/J_A(x)$
(Eq. 10). Identifying $q = \ln J$ as reduced heat, this is equivalent to
equation 9 of Crooks (ref 18). The reduced work is odd under time-reversal,
$w_A(x \to x') = -w_B(\bar{x}' \to \bar{x})$ (Eq. 11). Combining Eqs. 4 and
8-11 gives Eq. 12, which is equivalent to Eq. 5, so the scheme preserves
equilibrium in each replica.

## Illustrative Dynamics

A simple dynamical scheme illustrates the method. Suppose $\mathcal{R}_A$ and
$\mathcal{R}_B$ are described by different temperatures $T_A < T_B$ but the same
$H$, and construct Hamilton's equations augmented by a term proportional to
$\dot\lambda$ (Eq. 13), with $s_\lambda = (1/2T_\lambda)(dT_\lambda/d\lambda)$.
Here $T_\lambda$ interpolates from $T_0 = T_A$ to $T_1 = T_B$. The extra term
provides rudimentary temperature control during the work simulations. As
$\lambda$ is varied from 0 to 1 in replica A, the momenta are scaled up,
effectively heating the system; in replica B, the system is cooled. An
equivalent rescaling of momenta is standard practice in REM, where it is done
instantaneously rather than over a trajectory (Sugita and Okamoto, equation 12
of ref 24).

These dynamics do not preserve phase space volume:
$\nabla \cdot \dot{x} = N \dot\lambda s_\lambda \neq 0$, where $N$ is the number
of degrees of freedom. The Jacobian for a work simulation in $\mathcal{R}_A$ is
Eq. 14, and in $\mathcal{R}_B$ we have $J_B = J_A^{-1}$. (The fact that $J_A$ and
$J_B$ do not depend on initial conditions is specific to these dynamics.)

By heating the system in $\mathcal{R}_A$ and cooling it in $\mathcal{R}_B$, the
scaling term $\dot\lambda s_\lambda \mathbf{p}_i$ increases the probability of
accepting the configuration swap. For a system of ideal gas particles, evolution
under Eq. 13 exactly transforms a Maxwell-Boltzmann distribution from $T_A$ to
$T_B$, or vice versa. In this special case, $w_A + w_B = 0$ and $P_{\text{acc}} = 1$.

Even with stochastic equations of motion (Langevin dynamics or the Andersen
thermostat, see SI Appendix), it is useful to include the scaling term in
Eq. 13, as it dynamically adjusts the momenta in response to the changing
temperature $T_\lambda$.

## Efficiency Considerations and Numerical Results

The method is valid for arbitrary switching time $\tau$. Two limiting cases are
instructive. In the sudden limit $\tau = 0$, we have $w = \Delta h$ (Eqs. 1, 3,
4), as noted by Wyczalkowski and Pappu. In this case RENS reduces to REM, and if
there is little overlap between $p_A^{\text{eq}}$ and $p_B^{\text{eq}}$, then
$P_{\text{acc}} \ll 1$. In the opposite quasi-static limit $\tau \to \infty$, a
properly thermalized system evolves reversibly, and the reduced work equals the
reduced free energy difference (Eq. 15), where $f_i = -\ln \int dx\, e^{-h_i}$;
hence $w = 0$ and $P_{\text{acc}} = 1$. The acceptance probability can be
manipulated by adjusting the switching time $\tau$. Generically, the more slowly
the work simulation is performed, the greater the probability of accepting the
swap.

Note: if the scaling term of Eq. 13 is included, then only the potential energy
contributes to $w$, exactly as with REM when momenta are rescaled (equations 12
and 15 of Sugita and Okamoto, ref 24).

### Sample cost analysis

Consider M replicas, with primary interest in one of them (the primary replica);
the remaining replicas serve only to enhance sampling in the primary. The output
trajectory is obtained by concatenating the sampling intervals generated in the
primary replica after discarding the work intervals. Let $\bar\tau_{eq}$ denote
the average duration of a sampling interval. Let $t_c$ denote a characteristic
correlation time of the output trajectory, and define $X$ (Eq. 16).

The sample cost $t^*$ (Eq. 17) measures the total computational cost, summed over
all M replicas, of producing a single statistically independent sample in the
primary replica. The factor $(1+X)$ accounts for the overhead of the work
intervals. The smaller $t^*$, the more efficient the use of resources. The
correlation time $t_c$ generally decreases with increasing $M$ or $\tau$, but in
Eq. 17 this competes with the overhead factors $M$ and $1+X$.

### Model system

A model system of $n_p = 10$ particles moving independently in a rough potential
$U(x)$ (adapted from Frenkel and Smit, Chapter 14). $M = 2$ replicas, at
$T_A = 0.30$ and $T_B = 2.0$ (arbitrary units). In the primary replica at
$T_A = 0.30$, sampling is hindered by barriers separating local minima; at
$T_B = 2.0$ the particles jump from well to well. MD simulations were performed
using Eq. 13 combined with an Andersen thermostat.

While REM performs well for a single particle, at $n_p = 10$ it encounters
difficulties due to poor phase space overlap (Kofke's argument). In
$\mathcal{R}_A$ each particle is near a local minimum; in $\mathcal{R}_B$ they
are distributed more uniformly. A swap is accepted only if all particles in
$\mathcal{R}_B$ are near minima, unlikely for $n_p \gg 1$. With RENS, the work
simulation in $\mathcal{R}_B$ increases swap acceptability by shepherding
particles closer to the minima of $U(x)$.

Implementation: during a sampling interval, a work simulation was initiated at
random with attempt rate $r = 0.166$. Once initiated, the work interval lasted
for the prescribed switching time $\tau$, after which the replicas reverted to
sampling. The average sampling interval was $\bar\tau_{eq} = 1/r \approx 6.0$,
roughly 3 times the relaxation rate within one of the local wells.

25 test runs were performed with $\tau$ ranging from 0 to 100. Empirical
occupation probabilities for the 4 wells (following a tagged particle) agreed
with the equilibrium distribution within statistical error. The swap acceptance
frequency and average reduced work are plotted against the fraction of
simulation time in work intervals $f_{sw}$ (Eq. 18). With increasing $f_{sw}$
(or $\tau$) the system approaches the reversible limit $w = 0$,
$P_{\text{acc}} = 1$.

The correlation time was evaluated using block-averaging:
$t_c = (1/\sigma^2)\int_{-\infty}^{+\infty} dt\, c(t)$, where $\sigma^2$ and
$c(t)$ are the variance and autocorrelation of $n_4(t)$ (number of particles in
the fourth well). At $f_{sw} = 0$ (REM), the sample cost is high, $t^* > 4000$.
As $f_{sw}$ increases, $t^*$ drops significantly, reaching a broad minimum
$t^* \sim 450-500$ for $f_{sw} \sim 0.2-0.6$. For $f_{sw} > 0.6$, diminishing
returns set in.

The computational cost of the acceptance/rejection step and the possible
exchange of configurations were neglected, and identical costs per unit
simulation time were assumed for work and sampling intervals. To drop the latter
assumption, replace $X$ by $\alpha X$ in Eqs. 17 and 18, where $\alpha$ is the
observed CPU cost of a work simulation relative to a sampling trajectory of equal
duration. In the test runs, $\alpha = 2.9$. (With particle-particle
interactions, $\alpha$ would be closer to unity.)

REM test runs at $\tau = 0$ were performed for $M = 2, 3, \ldots 11$, with
$T_1 = 0.30$ and $T_M = 2.0$, intermediate replicas spaced evenly in $T^{-1}$.
The smallest sample cost, $t^* = 706$, was achieved with $M = 4$ replicas. This
is comparable to the optimal sample cost achieved with RENS using $M = 2$. Thus
RENS matches the efficiency of REM with fewer replicas.

## Discussion

When applying REM, the phase space overlap requirement dictates a minimum number
of replicas $M^*$. With RENS, work simulations increase phase space overlap,
allowing fewer replicas $M < M^*$. Reasons to exploit this flexibility:

- (i) With a cluster of P processors, RENS allows one replica per processor even
  if $P < M^*$.
- (ii) Replica exchange can be pictured as a diffusion process; $\sim M^2$
  successful swaps are needed for a trajectory to transit $\mathcal{R}_1$ to
  $\mathcal{R}_M$. Fewer replicas reduce interprocessor communication cost.
- (iii) REM is often implemented synchronously (swaps only after every replica
  completes a fixed duration), limited by the slowest processor. RENS lends
  itself to asynchronous implementation: a master process initiates work
  simulations in a randomly chosen replica pair while the rest continue sampling.
- (iv) Efficiency can be improved adaptively during a production run by adjusting
  the durations of the work simulations. If a low $P_{\text{acc}}$ between
  $\mathcal{R}_n$ and $\mathcal{R}_{n+1}$ is a bottleneck, increase the switching
  time $\tau_{n,n+1}$ for that pair.
- (v) Data generated during work simulations can be scavenged for equilibrium
  information by statistical reweighting (Hummer and Szabo, ref 31), analogous to
  Frenkel's waste-recycling of rejected Monte Carlo moves.

RENS can be enhanced by combination with solute tempering, generalized effective
potentials, large time steps, or artificial flow fields during the work
simulations.
