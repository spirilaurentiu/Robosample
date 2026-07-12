# Satisfying the fluctuation theorem in free-energy calculations with Hamiltonian replica exchange

Matthew A. Wyczalkowski and Rohit V. Pappu (2008). Phys. Rev. E 77, 026104.
DOI: 10.1103/PhysRevE.77.026104

## Abstract

An error measure, the hysteresis error, is developed from the Crooks fluctuation
theorem to evaluate sampling quality in free-energy calculations. Theory and
numerical free energy of hydration calculations show that Hamiltonian replica
exchange provides a direct route for minimizing the hysteresis error. Replica
exchange swap probabilities yield the rate at which the hysteresis error falls
with simulation length; this can be used to decrease bias and statistical errors
in free-energy calculations based on multicanonical simulations.

## I. Introduction and overview

The free energy of solvation in the canonical ensemble is the free-energy change
$\Delta F$ for transferring a solute from the gas phase to a fixed position in
the solvent. Following Kirkwood, a single coupling parameter $\lambda$ with
$0 \le \lambda \le 1$ modulates solute-solvent interactions in the potential.
The limits $\lambda=0$ and $\lambda=1$ correspond to pure solvent and solvent
plus fully grown solute. One runs a series of independent canonical simulations,
each at a distinct $\lambda$, generating work values used to estimate the
free-energy change across the $\lambda$ schedule.

Standard multicanonical free-energy calculations suffer slow convergence and
inaccurate $\Delta F$. Errors split into statistical and bias (finite-sampling)
errors. Statistical error decreases as the inverse square root of simulation
length. Bias error is an error of the mean itself and changes with simulation
length. Following Zuckerman and Woolf, bias errors arise because free-energy
estimates are nonlinear averages and work distributions have long, rarely
sampled tails that dominate the average. Rare events dominate free-energy
estimates, so the average drifts with simulation length. The magnitude of the
bias error is hard to quantify directly since it requires the true free-energy
difference.

Crooks derived a fluctuation theorem valid for stochastic, microscopically
reversible dynamics, relating forward and reverse dissipated-work distributions
(eq:1). If each canonical simulation samples its equilibrium ensemble adequately,
the observed dissipated-work distributions satisfy eq:1. This work develops the
hysteresis error $\epsilon_H$, which quantifies how far observed work
distributions deviate from the Crooks theorem. Hamiltonian replica exchange, a
multicanonical equilibration technique, reduces $\epsilon_H$. The average swap
probability relates to the overlap between equilibrium ensembles and to the rate
at which $\epsilon_H$ falls, giving a route to an optimized $\lambda$ schedule.

## II. Theory

### A. Background

The free energy of replica $i$ at temperature $T$, with potential
$U_i(\Gamma) = U(\Gamma, \lambda_i)$ over configuration $\Gamma$, is given by
eq:2. At equilibrium the probability of observing $\Gamma$ is eq:3.

To calculate the free-energy change $\delta F$ for switching the Hamiltonian
from $U_0$ to $U_1$, run simulations at $\lambda_0$ and $\lambda_1$ and compute
forward and reverse work (eq:4a, eq:4b). Configurations are drawn from the
equilibrium ensemble of $U_0$ (forward) and $U_1$ (reverse). The free-energy
perturbation (FEP) method gives two independent estimators for $\delta F$
(eq:5a, eq:5b). Both are associated with switching $\lambda_0 \to \lambda_1$.
They have different convergence rates, so with finite sampling they generally
differ.

The Bennett acceptance ratio (BAR) uses both $W^F$ and $W^R$ distributions,
is generally more accurate, and is used for the numerical free-energy estimates
here. BAR is not used for the theoretical development.

### B. Hysteresis error

The hysteresis error $\epsilon_H$ is the difference between forward and reverse
$\delta F_{FEP}$ estimates (eq:6). It has contributions from both statistical
and bias error; the bias errors of the two estimators are typically opposite in
direction. In averages over multiple short simulations, the dominant
contribution to the average hysteresis error is the sum of the forward and
reverse FEP bias. $\epsilon_H$ is taken as a measure of sampling quality; the
goal is to minimize its magnitude between all pairs of neighboring replicas.

Switching $\lambda_0 \to \lambda_1$ performs nonequilibrium work. The dissipated
work is the difference between the work performed and the free-energy change
(eq:7a, eq:7b). Crooks equates $W_D^F$ and $W_D^R$ to the entropy production
caused by changing $\lambda_0 \to \lambda_1$ and $\lambda_1 \to \lambda_0$,
respectively, for a given configuration.

The distributions $P_F(W_D)$ and $P_R(W_D)$ give the probability of realizing a
specific dissipated-work value forward and reverse. They are related by the
fluctuation theorem (eq:1), rederived in the Derivation below for instantaneous
switching. With finite sampling, eq:1 is not satisfied exactly; introducing an
arbitrary error term $\epsilon_{FT}^*$ and observed distributions gives eq:8.
The Crooks theorem is recovered ($\epsilon_{FT}^*=0$) when observed
distributions match ideal ones. The hysteresis error and fluctuation error are
related by eq:9. The more closely a simulation obeys eq:1, the smaller
$\epsilon_H$.

Figure 1 note: replica exchange is visualized by plotting independent
configurations $\Gamma_0$ and $\Gamma_1$ on orthogonal axes, with the joint
equilibrium ensemble $\rho_N$ as an isocontour. A swap is a reflection of the
configuration pair $\gamma$ about the $\Gamma_0 = \Gamma_1$ diagonal axis. Swaps
are feasible only for configuration pairs belonging to both $\rho_N$ and its
swapped image $\rho_N'$; that overlap region ($\rho_{swap}$) is where the
integrand of eq:17b is large and its integral is the average swap probability.

### C. Replica exchange

In Hamiltonian replica exchange, Monte Carlo moves exchange configurations
$\Gamma$ (equivalently parameters $\lambda$) between two replicas with
probability eq:10, where $\Delta U_{swap}$ is given by eq:11a-11c. $\Gamma_0$
and $\Gamma_1$ are configurations drawn at random from the equilibrium ensembles
of $U_0$ and $U_1$. Write $\gamma = (\Gamma_0, \Gamma_1)$ for a pair and
$\gamma' = (\Gamma_1, \Gamma_0)$ for the swapped pair.

Since $\Gamma_0$ and $\Gamma_1$ are independent, the native joint probability is
$\rho_N(\gamma)$ (eq:12a) and the swapped joint probability is $\rho_N'(\gamma)$
(eq:12b). At equilibrium, the relative probability of swapped versus native
configurations is eq:13 (the interreplica equilibrium relationship), derived
from eq:12, eq:3, and eq:11a.

In an infinitely long simulation eq:13 holds exactly. With finite simulations,
inadequate sampling gives inaccurate probability estimates. Replica exchange
directly populates swapped configurations, improving the statistics of
$\rho_N'$, so interreplica equilibrium is reached more quickly. The degree to
which eq:13 is satisfied determines the magnitude of the hysteresis error.
Introducing a small sampling error $\rho_\epsilon$ in $\rho_N'$ and integrating
over all configuration pairs gives eq:14, so $\epsilon_H$ is minimized when the
estimated swapped-configuration probabilities are consistent with equilibrium.

### D. Swap probability

The Metropolis function (eq:10) is not analytical, so for interpretation the
Fermi swap probability $p_{swap} = f(\beta \Delta U_{swap})$ is used, with
$f(x)$ given by eq:15. $p_{swap}$ denotes the Fermi swap probability and
$P_{swap}$ the Metropolis swap probability; theory uses $p_{swap}$, while
actual replica exchange moves are accepted or rejected with $P_{swap}$. Both
yield the Boltzmann distribution of swapped and unswapped configurations
(eq:13).

The average Fermi swap probability for two independently evolving systems is
eq:16a-16b, which can be rewritten as an overlap integral eq:17a-17b. A large
average swap probability means large overlap between the two replicas'
equilibrium distributions; a low value means the replicas adopt distinct
configurations. Expanding eq:16a in a Taylor series about
$\lambda = \lambda_0 + \delta_\lambda$ gives, to leading order, the linearized
swap probability eq:18, where $C_\lambda = \mathrm{var}(\partial U/\partial
\lambda)$. $C_\lambda$ determines the rate at which the average swap probability
declines as $\delta_\lambda$ increases; the linear analysis is accurate only for
small $\delta_\lambda$.

### E. Swap probability and the hysteresis error convergence rate

The average swap probability measures how quickly the hysteresis error
decreases. Since forward and reverse FEP estimators do not converge at equal
rates, the slower governs convergence of $\epsilon_H$. Rewriting eq:5a gives
eq:19. For this to hold, one must sample configurations where $W_D^F < 0$;
since dissipated work is on average positive (second law), such configurations
are rare. The convergence rate of $\delta F_{FEP}^F$ is governed by the
probability of observing negative dissipated forward work; likewise
$\delta F_{FEP}^R$ by observations of $W_D^R < 0$. Graphically (eq:20a, eq:20b),
$W_D^F < 0$ corresponds to sampling from $\rho_0$ where $\rho_1 > \rho_0$, and
$W_D^R < 0$ requires $\rho_0 > \rho_1$ when sampled from $\rho_1$.

$\Delta U_{swap} = W_D^F + W_D^R$ is negative whenever $\rho_N' > \rho_N$ (by
eq:13). The larger the swap domain (given by the average swap probability), the
more frequently negative $W_D^F$ and $W_D^R$ are observed and the more quickly
the hysteresis error converges.

## III. Methods

The system consists of 21 replicas, each with a different $\lambda$, simulated
independently to obtain equilibrium statistics. $\lambda$ controls the nonbonded
interactions between an acetamide (ACE) solute and water molecules. The
Lennard-Jones and Coulomb interactions between water and ACE are scaled by
$\lambda_{LJ}$ and $\lambda_C$; both are scaled simultaneously,
$\lambda_{LJ} = \lambda_C = \lambda$. $\lambda$ varies across the 21 replicas
from 0 to 1 in increments of 0.05. The functional forms of the scaled Coulomb
(eq:B1) and Lennard-Jones (eq:B2, eq:B3a, eq:B3b) potentials are given in the
appendix.

Each replica has 343 water molecules and one ACE molecule, rigid and fixed in
the central box. All simulations were at constant temperature (298 K) and volume
(21.8 Angstrom cubic box) using Metropolis Monte Carlo. OPLS-AA force field and
TIP4P four-site water model were used. Minimum image boundary conditions and
spherical cutoffs were used: 10.5 Angstrom for electrostatics (group based) and
10 Angstrom for van der Waals (atom based). No long-range corrections. Simulated
with MCCCS Towhee.

Initial configurations for all replicas were identical (endpoint of a
preequilibration run). Each replica ran $2\times10^6$ cycles (a cycle = 343 MC
moves; each move combines rotations and translations of a randomly chosen water).
The initial $10^5$ cycles were discarded for equilibration. Average acceptance
rate was 31%.

The replica exchange simulation consists of simulation rounds (each replica
evolves independently) separated by swap rounds. The simulation round length was
drawn from a normal distribution, mean 500 and standard deviation 50 cycles.
500 cycles is the approximate energy autocorrelation "time." The swap round is
$21^2$ swap attempts between randomly selected replica pairs. Allowing swaps
beyond neighboring replicas increases efficiency.

The native $[U_i(\Gamma_i)]$ and foreign $[U_{j\neq i}(\Gamma_i)]$ potential
energies, and $dU/d\lambda_C$ and $dU/d\lambda_{LJ}$ (with $dU/d\lambda =
dU/d\lambda_{LJ} + dU/d\lambda_C$), were saved every ten cycles and
post-processed to obtain free energies, hysteresis error, swap probabilities,
and $C_\lambda$, regardless of whether actual swaps took place. The total free
energy of hydration $\Delta F$ is the sum of neighboring free-energy changes
$(\delta F)_i$ from BAR (eq:dF_sum). $\epsilon_{rms}$ is the root mean square of
the neighboring hysteresis errors (eq:eps_rms).

Statistical errors for $\Delta F$ were estimated by bootstrap: from $N$
observations, draw $n^*$ observations at random with replacement to form one
$\Delta F^*$; repeat 10000 times; the standard deviation of $\Delta F^*$ is the
error. $n^*$ is the expected number of independent observations, here
$n^*=1900$, assuming one independent observation per two internal-energy
autocorrelation "times."

## IV. Results

The hydration free energies for acetamide are consistent with other researchers.
Calculations were in the canonical ensemble, giving Helmholtz $\Delta F$, while
experimental and other computational values are Gibbs $\Delta G$. Since the box
volume at $\lambda=0$ corresponds to 1 atm, the distinction is negligible;
N-P-T test calculations confirm this.

The methods to calculate $\Delta G$ in N-P-T do not differ from $\Delta F$ in
N-V-T; the swap probability (eq:10) does not change since pressure-volume work
is reversible and does not contribute to dissipated work.

The root-mean-square hysteresis error is lowered by an order of magnitude when
replica exchange is coupled to multicanonical sampling. The bootstrap
statistical error remains unaffected: fluctuations in $\delta F$ originate in
fluctuations of the underlying work distribution (eq:1), and as long as both
simulations sample the work distribution adequately they have similar
statistical error. Low statistical errors can also be caused by inadequate
sampling of the work distributions.

A simulation with replica exchange achieves the same magnitude of
$\epsilon_{rms}$ about five times more quickly (4-8 times shorter) than one
without.

The swap probability shows downward spikes for $\lambda$ where the hysteresis
error is large, with a positive spike in $C_\lambda$ in the same region
(consistent with eq:18). $C_\lambda$ can be obtained from a single simulation,
whereas swap probability requires two. Evaluation of $C_\lambda(\lambda)$ using
a coarse $\lambda$ schedule can identify low-swap-probability regions and
construct optimal schedules.

## V. Discussion

### A. Physical interpretation of the $C_\lambda$ profile

Water occupancy around the growing solute decreases rapidly near
$\lambda \sim 0.15$. Expulsion and rearrangement of water during cavitation
shifts the equilibrium ensemble, giving a pronounced spike in $C_\lambda$.
Smaller shifts near $\lambda=1$ reflect electrostatic effects. $C_\lambda$
profiles serve as probes for detecting large shifts in equilibrium ensembles;
regions where ensembles change most rapidly contribute the largest free-energy
errors.

### B. Optimal $\lambda$ schedule for free-energy calculations

With the number of replicas and simulation length fixed, the root-mean-square
hysteresis error is decreased by optimizing the $\lambda$ schedule. The swap
probability gives the rate at which the average hysteresis error falls between
two replicas; in an optimized simulation it is uniform across all replica pairs.
Reasonable schedules can be built using the linearized swap probability (eq:18):
perform preliminary simulations to obtain $C_\lambda$ along a coarse schedule
(these need not be as long as production runs since $C_\lambda$ converges more
quickly than $\delta F$), then adjust the schedule so the linear swap
probability is uniform, or shift replicas from where $C_\lambda$ is small to
where it is large. Both are approximate and break down when the linear-response
assumption fails; they may be applied iteratively. The aim is to place replicas
close together where the $C_\lambda$ profile shows spikes.

### C. Replica exchange

Replica exchange lets a replica access a distant part of its equilibrium
ensemble in one step but is no substitute for conformational exploration within
a replica. The hysteresis error does not report on intrareplica sampling
quality. As an extreme case, frozen replicas undergoing only swap moves attain
the distribution of eq:13 with a modest number of swaps and zero hysteresis
error, yet the intrareplica probability distribution is not obtained. In
practice most Monte Carlo moves must be within a replica. The optimal frequency
of swap moves is an open question; preliminary simulations suggest more frequent
swaps reduce the hysteresis error more quickly.

Compared with $\lambda$-dynamical methods, where $\lambda$ evolves dynamically
according to the conjugate force $\partial U/\partial \lambda$, replica exchange
uses a Markov chain to distribute configurations across $\lambda$ values
according to a Boltzmann distribution. Sampling difficulties in both cases are
associated with phase changes and large variance of $\partial U/\partial
\lambda$ ($C_\lambda$), giving low swap probabilities or hard-to-traverse
regions. Replica exchange is readily parallelizable.

## VI. Summary and conclusion

Swapping configurations between replicas is a nonequilibrium work process
described by the Crooks fluctuation theorem. The hysteresis error $\epsilon_H$
measures how closely a simulation reproduces the work distributions between a
pair of replicas. It reports on the combined bias of forward and reverse FEP
and measures how completely individual replicas sample their equilibrium
ensemble. The root-mean-square hysteresis error may be decreased by longer
simulations, replica exchange, an improved $\lambda$ schedule, or all of these.
The average swap probability, and the related $C_\lambda$, determines the rate
at which the hysteresis error decreases; a uniform average swap probability
makes the hysteresis error fall evenly between all replica pairs.

## Derivation (not implemented): Appendix A

### A.1 Fluctuation theorem derivation

For instantaneous switching $\lambda_0 \to \lambda_1$ (forward) and
$\lambda_1 \to \lambda_0$ (reverse), expanding the ratio $\rho_0/\rho_1$ with
eq:3 gives eq:20a (= eq:A1a):

$$ \frac{\rho_0(\Gamma)}{\rho_1(\Gamma)} = \exp[\beta(F_0 - F_1) - \beta(U_0 - U_1)] = \exp(-\beta \delta F + \beta W^F) = \exp(\beta W_D^F), $$

and similarly $\rho_1(\Gamma)/\rho_0(\Gamma) = \exp[\beta W_D^R(\Gamma)]$.
Integrating $\rho_0$ over configurations with a fixed dissipated-work value
$W_D$ and using $W_D^F(\Gamma) = -W_D^R(\Gamma)$, together with the definitions
of the forward and reverse dissipated-work distributions

$$ P^F(W_D) = \int d\Gamma\, \rho_0(\Gamma)\, \delta[\beta W_D - \beta W_D^F(\Gamma)], $$
$$ P^R(W_D) = \int d\Gamma\, \rho_1(\Gamma)\, \delta[\beta W_D - \beta W_D^R(\Gamma)], $$

yields $\exp(-\beta W_D) P^F(\beta W_D) = P^R(-\beta W_D)$, equivalent to eq:1.

### A.2 Fluctuation theorem and hysteresis error

Rewriting eq:8, inserting the $\delta F_{FEP}^R$ definition (eq:5b) into eq:6,
expanding reverse work with eq:7b, and using the $\delta F_{FEP}^F$ estimate for
$\delta F$ gives $\epsilon_H = -\beta^{-1} \ln [\langle \exp(-\beta W_D^R)
\rangle_1^*]$. Expanding the estimated average as an integral over $\beta W_D^R$
with $P_R^*$ the normalized histogram, changing the dummy variable to
$-\beta W_D$, and using eq:8 reduces to eq:9.

### A.3 Interreplica equilibrium and hysteresis error

Rewriting eq:13 with a small error $\rho_\epsilon$ in $\rho_N'$,
$\rho_N' + \rho_\epsilon = \rho_N \exp(-\beta \Delta U_{swap})$, integrating over
all configuration pairs and expanding $\Delta U_{swap}$ with eq:11b, the
$\rho_N'$ term integrates to 1. Taking the logarithm, dividing by $\beta$, and
using $\ln(1+x) \simeq x$ for small $x$ with the definition of $\epsilon_H$
(eq:6) gives eq:14.

### A.4 Linearized average swap probability

Define $\mu \equiv \beta \Delta U_{swap}$ for replicas whose $\lambda$ differ by
$\delta$. Expanding $U_\delta$ as a Taylor series about $\lambda_0$ with
$V_0 \equiv \partial U/\partial \lambda |_{\lambda_0}$ and
$W_0 \equiv \partial^2 U/\partial \lambda^2 |_{\lambda_0}$, and using the
expansions of $\exp(x)$ and $1/(1+x)$, the Fermi swap probability is
$p_{swap} = \tfrac12 - \tfrac14 \mu + O(\mu^3)$. Averaging over configuration
pairs and evaluating $\langle V_0 \rangle_\delta$ and $\langle W_0
\rangle_\delta$ via the partition-function expansion $Q_\delta = Q_0[1 - \beta
\delta \langle V_0 \rangle_0 + O(\delta^2)]$ yields eq:18 (= eq:A10):

$$ \langle\langle p_{swap}\rangle_0\rangle_{\delta} = \frac{1}{2} - \frac{\beta^2 \delta^2}{4} (\langle V_0^2\rangle_0 - \langle V_0\rangle_0^2) + O(\delta^3), $$

valid for small $\delta$.

## Appendix B: $U_{LJ}$ and $U_C$ functional forms

Functional forms of the Coulomb and Lennard-Jones potentials were developed for
this work under three criteria: (i) overlapping solute-solvent configurations
observable for $\lambda=0$, with swaps permitted at reasonable frequency for
small $\lambda$ and falling off quickly thereafter; (ii) $\partial U/\partial
\lambda$ not always zero at $\lambda=0$ (to avoid TI complications); (iii) since
$\lambda_{LJ} = \lambda_C$, Lennard-Jones repulsion must dominate Coulombic
attraction at very small atomic separations.

Coulomb scaling: a modified linear soft-core scaling (eq:B1), with $\alpha_C$
controlling the soft-core term ($\alpha_C = 1.5$ Angstrom). Lennard-Jones
scaling: the general form (eq:B2) with the exponential soft core (eq:B3a,
eq:B3b), with $a=4$, $k=1$, $\alpha_{LJ}=0.5$ Angstrom.
