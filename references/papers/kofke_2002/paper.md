# On the acceptance probability of replica-exchange Monte Carlo trials

David A. Kofke, Department of Chemical Engineering, University at Buffalo, SUNY.
J. Chem. Phys. 117, 6911-6914 (2002). DOI: 10.1063/1.1507776.
Erratum: J. Chem. Phys. 120, 10852 (2004). DOI: 10.1063/1.1738103.

NOTE: This body incorporates the published erratum. Eq. (12) and Eq. (13) below
are the corrected forms; see the flags in `equations.md`.

## Abstract

An analysis is presented of the average probability of accepting an exchange
trial in the parallel-tempering Monte Carlo molecular simulation method.
Arguments are given that this quantity should be related to the entropy
difference between the phases, and results from simulations of a simple
Lennard-Jones system support this argument qualitatively. Another analysis
based on the energy distributions of a replica pair yields an exact expression
for the trial-move acceptance probability in terms of the overlap of these
distributions. A more detailed expression uses an approximation of constant
heat capacity, and an asymptotic form good for large system sizes is reported.
The detailed analyses are in quantitative agreement with the simulation data.
Treating the energy distributions as Gaussians is an inappropriate way to
analyze the acceptance probability.

## I. Introduction

Replica-exchange Monte Carlo (also known as parallel tempering) enhances
sampling in systems with many disconnected low-energy regions of configuration
space. Several independent realizations of the system are sampled
simultaneously, each differing in the temperature that governs sampling. As the
simulation proceeds, configurations from a pair of systems (subscript 0 and 1)
are exchanged occasionally, and acceptance of the new state uses the
probability

<!-- eq:1 -->
$$ \min[1, \exp(-(\beta_0 - \beta_1)(U_1 - U_0))] $$

where the argument of the exponential is the difference in the sum of $\beta U$
for the two systems, before and after the trial. Here $\beta_i$ is the
reciprocal temperature $1/kT$ of system $i$, with $k$ the Boltzmann constant,
and $U_i$ is the potential energy of system $i$ before the exchange.

The higher-temperature systems surmount barriers between low-energy states
better, and supply the low-temperature systems with trial configurations
covering a broader range of configuration space. Unlike simulated annealing,
the high-temperature system need not quench and climb out of the well; with
parallel tempering, quenching (acceptance of the exchange trial) occurs only
when the high-temperature system already resides in a low-energy region. All
replicas remain equilibrated at all times.

If replica temperatures differ by too much, the high-temperature system samples
configurations largely unacceptable to the low-temperature one, so the
acceptance rate becomes small. The remedy is to introduce intermediate-
temperature stages. A related problem arises as system size $N$ increases: the
average energy difference grows proportionally with $N$, while the width of the
energy distribution grows only as $\sqrt{N}$, so the distributions become
relatively narrower and farther apart and may lose their overlap. Again the
remedy is intermediate stages.

## II. Entropy model

The energy-distribution picture implies an inappropriate symmetry between the
systems. A consideration in configuration space is more telling. The
configurations important to the low-temperature system form a subset of those
important to the high-temperature system (Fig. 1(b) shows them as disconnected
subregions). The high-temperature trajectory moves about its region and
occasionally happens into one of the low-temperature subregions. An exchange
trial performed at that time is accepted; otherwise it is rejected.

Thus the acceptance rate for the trials equals the likelihood that the system
sampling the high temperature resides in one of the low-temperature important
regions. This likelihood is approximated as the fraction those regions occupy
in the high-temperature important region. The size of an important region
correlates with the entropy and is proportional to $\exp(S/k)$. Thus the
average acceptance probability $\bar{p}_{\rm acc}$ (the fraction of all trials
accepted) goes as

<!-- eq:2 -->
$$ \bar{p}_{\rm acc} \sim \exp(-\Delta S/k) $$

where $\Delta S > 0$ is the entropy of the high-temperature system minus that
of the low-temperature one.

Two conclusions follow. First, because $\Delta S$ is extensive, the acceptance
rate decreases exponentially with system size $N$. Second, the heat capacity
modulates the decrease as the temperature interval is widened. Assuming a
constant heat capacity over the interval, and using $T\Delta S = N c_V \Delta T$,

<!-- eq:3 -->
$$ \bar{p}_{\rm acc} \sim \left(\frac{T_0}{T_1}\right)^{N c_V/k} $$

where $T_0 < T_1$ and $c_V$ is the molar heat capacity at constant volume. The
connection to heat capacity is consistent with the distribution-overlap view,
since the distance between distributions, their widths, and their overlap are
all related to the heat capacity.

As a test, replica-exchange Monte Carlo simulations of a Lennard-Jones (LJ)
liquid were performed for three system sizes (32, 64, and 108 particles) and
several temperature differences, recording the acceptance rate. All simulations
used number density $\rho = 0.80$, temperatures from $T = 1.0$ to $2.0$ in
increments of $0.2$ (LJ units, $\sigma = \epsilon = 1$). Only one pair of
systems was studied per simulation. Each simulation ran between $5\times10^5$
and $2\times10^7$ elementary Monte Carlo trials. In each trial, either an atom
displacement or a replica exchange was attempted, chosen at random with
displacements attempted 100 times more often than exchanges. Entropy
differences were computed by thermodynamic integration of the internal energy.
The molar heat capacity at $T = 1.5$ was $C_V/Nk = 0.78$.

For a fixed system size, the acceptance rate correlates strongly with the
entropy. There is an additional system-size dependence not captured by the
extensive entropy: the decay of the acceptance probability with increasing
system size is attenuated, so a size increase does not lead to the full
exponential decay implied by Eq. (2). A more detailed model is required.

## III. Energy-distribution model

Define energy distributions $p_i(U) = p(U;\beta_i)$ such that $p_i(U)dU$ is the
likelihood to observe a system of temperature $\beta_i$ with energy in the
range $U$ to $U+dU$. The acceptance probability sums over all system energies
$U_0, U_1$ the probability that the two systems have those energies times the
probability that a trial exchanging them is accepted (Eq. 1):

<!-- eq:4 -->
$$ \bar{p}_{acc} = \int_{U_m}^{\infty} dU_1\, p_1(U_1) \int_{U_m}^{\infty} dU_0\, p_0(U_0) \times \min[1, \exp(-(\beta_0 - \beta_1)(U_1 - U_0))] $$

where $U_m$ is the lowest possible energy. The energy distributions decompose as

<!-- eq:5 -->
$$ p_i(U) = \frac{1}{Q(\beta_i)} \Omega(U) \exp(-\beta_i U) $$

where $\Omega(U)$ is the density of states, independent of $\beta$, and
$Q(\beta_i)$ is the canonical-ensemble partition function. Requiring
$\beta_0 > \beta_1$ ($T_0 < T_1$) lets the min function be evaluated by whether
$U_1$ is greater or less than $U_0$. Substituting Eq. (5) into Eq. (4) gives

<!-- eq:6 -->
$$ \bar{p}_{\text{acc}} = 2 \int_{U_m}^{\infty} dU_1 \int_{U_m}^{U_1} dU_0 \frac{\Omega(U_1)\Omega(U_0)}{Q(\beta_1)Q(\beta_0)} e^{-\beta_0 U_1} e^{-\beta_1 U_0} $$

which is expressed back in terms of the distribution functions as

<!-- eq:7 -->
$$ \bar{p}_{\text{acc}} = 2 \int_{U_m}^{\infty} dU_0\, p_0(U_0) \int_{U_m}^{U_0} dU_1\, p_1(U_1) $$

This integral quantifies the overlap of the distribution functions: it ranges
from zero to unity for cases from no overlap to complete overlap
($p_0(U) = p_1(U)$). The inner integral is for the high-temperature
distribution $p_1$, which peaks at larger $U$ than $p_0$. If $p_0$ peaks and
returns to zero before $p_1$ begins to rise, there is no contribution.

Assume the heat capacity is constant across the temperature range between
$\beta_0$ and $\beta_1$. Then the density of states across this range is

<!-- eq:8 -->
$$ \Omega(U) = \left(1 + \frac{1}{C}\beta_r(U - U_r)\right)^C \Omega(U_r) $$

where the $r$ subscript indicates an arbitrary reference state with energy
$U_r$ and temperature $T_r$, and $C \equiv C_V/k$ is the extensive,
constant-volume heat capacity in units of the Boltzmann constant.

### Derivation (not implemented)

Eq. (8) is derived by eliminating $\beta$ between the constant-$C$ expressions
for the entropy $S(\beta) = S_r + kC\ln(\beta_r/\beta)$ and energy
$U(\beta) = U_r + C(1/\beta - 1/\beta_r)$, along with the bridge equation
$S = k\ln\Omega$.

For this density of states the normalized energy distribution is

<!-- eq:9 -->
$$ p_i(U) = \frac{\beta_i}{C\,\Gamma(C)} [\beta_i(U - U_m)]^C \exp[-\beta_i(U - U_m)] $$

with $U_m = U_r - C/\beta_r$ identified as the minimum possible energy for this
model of the density of states.

Combining Eqs. (7) and (9) and using the variable substitutions

<!-- eq:sub -->
$$ \kappa \equiv \frac{\beta_1(U_1 - U_m)}{\beta_0(U_0 - U_m)}, \qquad \gamma = \beta_1(U_1 - U_m) + \beta_0(U_0 - U_m) $$

gives

<!-- eq:10 -->
$$ \bar{p}_{acc} = 4 \frac{\Gamma(2C)}{[\Gamma(C)]^2} \frac{2C+1}{C} \int_0^{\beta_1/\beta_0} d\kappa \frac{\kappa^C}{(1+\kappa)^{2(C+1)}} $$

With $\beta_1 < \beta_0$, the range of the integral is zero to unity at most. As
the temperature difference grows, the upper limit decreases and the acceptance
probability diminishes. The acceptance probability depends on the temperatures
only through their ratio, connecting to the entropy difference, because for
constant heat capacity

<!-- eq:11 -->
$$ \Delta S/k = -C \ln(\beta_1/\beta_0) $$

Figure 2 shows that Eq. (10), evaluated numerically with entropy from Eq. (11),
describes the full behavior of the average acceptance probability with
quantitative accuracy.

An asymptotic analysis reveals the additional size dependence beyond the
extensive entropy. Expressed partially in terms of the entropy difference via
Eq. (11), the result (for $\beta_1/\beta_0$ not too close to unity) is

<!-- eq:12 -->
$$ \bar{p}_{\rm acc} \sim \frac{\exp(-\Delta S/k)}{(\pi C)^{1/2}} \left[ \frac{4}{(1+B)^2} \right]^{C+1} \frac{1+B}{1-B} (1+O(C^{-1/2})), \quad C \to \infty $$

where $B = \beta_1/\beta_0 < 1$. NOTE: this is the erratum-corrected form; the
original 2002 text had $(1+B)/(1-2B+B^2) = (1+B)/(1-B)^2$ in place of
$(1+B)/(1-B)$. The bracketed term with the $C+1$ exponent exceeds unity and is
largely responsible for attenuating the decay of the acceptance probability
with increasing system size.

Treating the distributions as Gaussians is inappropriate. The concern is the
overlap of the distributions and thus their tail behavior. The Gaussian form is
inconsistent with the exponential decay on the high-energy side and the
algebraic decay on the other side, so it fails to capture the dependence of the
acceptance probability on the tail overlap. A Gaussian analysis leads via
Eq. (7) to

<!-- eq:13 -->
$$ \bar{p}_{\text{acc}} = \text{erfc}\left[ \left(\tfrac{1}{2}C\right)^{1/2} \frac{1 - \beta_1/\beta_0}{\left(1 + (\beta_1/\beta_0)^2\right)^{1/2}} \right] $$

NOTE: this is the erratum-corrected form; the original 2002 text omitted the
$(1+(\beta_1/\beta_0)^2)^{1/2}$ denominator, giving
$\text{erfc}[(\tfrac{1}{2}C)^{1/2}(1-\beta_1/\beta_0)]$. The corrected formula
gives more credence to a Gaussian approximation of the energy distribution.

## IV. Conclusion

The primary results are: identifying the acceptance probability with the
entropy (Eq. 2), quantifying its dependence on the overlap (Eq. 7), and
developing, within a constant-heat-capacity approximation, exact (Eq. 10) and
asymptotic (Eq. 12) formulas. This provides an analytical justification for the
empirical observation that the acceptance probability can be made uniform
across a multireplica partition by selecting temperature intervals such that
all adjacent temperatures are in a fixed ratio.

The analysis considers parallel tempering in its original form, where replicas
differ in temperature only. Generalizations such as hyperparallel tempering
(replicas differing in chemical potential or Hamiltonian) can also be analyzed,
but the connection to entropy is less certain because the subset relation of
Fig. 1(b) is not assured. The distribution-function analysis extends
straightforwardly, and is expected to show that the acceptance probability
depends on the ratio of the relevant field variable and the appropriate
susceptibility.
