# Hybrid Monte Carlo

Simon Duane, A.D. Kennedy, Brian J. Pendleton, Duncan Roweth. *Physics Letters B* 195(2), 1987, 216-222.

## Abstract

We present a new method for the numerical simulation of lattice field theory. A hybrid (molecular dynamics / Langevin) algorithm is used to guide a Monte Carlo simulation. There are no discretization errors even for large step sizes. The method is especially efficient for systems such as quantum chromodynamics which contain fermionic degrees of freedom. Detailed results are presented for four-dimensional compact quantum electrodynamics including the dynamical effects of electrons.

## Motivation

Computer simulation of lattice field theory with fermionic degrees of freedom is notoriously difficult. Because of the Grassmann nature of the fermions one cannot use standard methods to perform Monte Carlo calculations. Instead one integrates out the fermion fields leading to an effective action (the logarithm of the fermion determinant) which is highly nonlocal in the remaining bosonic degrees of freedom. Local updates therefore require calculations that depend on the state of the whole system.

A popular remedy is to replace the Monte Carlo algorithm by an equation of motion describing the evolution of the system in a new fictitious time variable $\tau$. The solution for asymptotically large times gives the desired probability distribution. Numerical integration schemes avoid the non-locality problem by evolving the whole system in parallel through many small steps, but replace one approximation (incomplete equilibration) with another (truncation errors from numerical integration).

The purpose of this paper is to describe a new algorithm which combines the ease of calculation of an equation-of-motion method and the absence of truncation error of an exact Monte Carlo. The algorithm involves **parallel updates of fields at all lattice sites followed by an accept/reject decision for the whole configuration**, yet has *no* truncation errors. The molecular-dynamics step size in the standard hybrid stochastic method is chosen as large as possible while keeping the Monte Carlo acceptance rate satisfactorily high.

## Expectation values and Markov chains

The fundamental objective of a quantum field theory is to calculate the expectation value of some operator $\langle\Omega\rangle$, where the field $\phi$ has dynamics governed by the action $S(\phi)$:

<!-- eq:1 -->
$$\langle \Omega \rangle = \frac{1}{Z} \int [d\phi]\, \exp[-S(\phi)]\, \Omega(\phi),$$

with partition function

<!-- eq:2 -->
$$Z = \int [d\phi]\, \exp[-S(\phi)].$$

The Monte Carlo method computes $\langle\Omega\rangle$ by generating field configurations at random with probability $P_S(\phi) = (1/Z)\exp[-S(\phi)]$, and then measuring

<!-- eq:3 -->
$$\bar{\Omega} = \frac{1}{T} \sum_{t=1}^{T} \Omega(\phi_t)$$

on a sequence $\{\phi_t\}$ of such configurations. As $T \to \infty$,

<!-- eq:4 -->
$$\bar{\Omega} = \langle \Omega \rangle + \mathcal{O}(1/\sqrt{T}).$$

The most useful technique for generating a sequence of configurations with the desired distribution is to construct a Markov process, which generates a new configuration $\phi'$ from its predecessor $\phi$ with probability $P_M(\phi \mapsto \phi')$. Any Markov process converges to a unique fixed-point distribution $P_S$ provided it is ergodic and satisfies detailed balance:

<!-- eq:5 -->
$$P_S(\phi)\, P_M(\phi \mapsto \phi') = P_S(\phi')\, P_M(\phi' \mapsto \phi).$$

It is convenient to construct the Markov process in two parts. First choose a candidate configuration $\phi'$ with probability $P_C(\phi \mapsto \phi')$ by some as-yet-unspecified procedure, then accept $\phi'$ with probability $P_A(\phi \mapsto \phi')$ or reject it and keep $\phi$. One choice of $P_A$ that enables detailed balance to be satisfied for any $P_C$ is a generalization of the Metropolis algorithm:

<!-- eq:6 -->
$$P_{A}(\phi \mapsto \phi') = \min\!\left(1,\ \frac{P_{S}(\phi')\, P_{C}(\phi' \mapsto \phi)}{P_{S}(\phi)\, P_{C}(\phi \mapsto \phi')}\right).$$

We require a method for choosing candidate configurations that can be computed efficiently and whose reverse probability $P_C(\phi' \mapsto \phi)$ is easy to obtain. Since we update all fields simultaneously, we insist that the acceptance rate $P_A$ is large and does not depend too strongly on system size, and we want to minimize correlation between successive configurations.

## Hamiltonian guidance dynamics

An elegant method follows the hybrid molecular dynamics / Langevin algorithm. We introduce a new "computer" time parameter $\tau$ and a Hamiltonian dynamics for the field $\phi(\tau)$. Introducing conjugate momenta $\pi(\tau)$, we impose the Hamiltonian by fiat:

<!-- eq:7 -->
$$H'(\phi,\pi) \equiv \tfrac{1}{2}\pi^2 + S'(\phi),$$

where $S'$ is some arbitrary (guidance) action. This gives the equations of motion:

<!-- eq:8 -->
$$\dot{\phi} = \delta H'/\delta \pi = \pi, \qquad \dot{\pi} = -\delta H'/\delta \phi = -\delta S'/\delta \phi.$$

The Hamiltonian itself is a constant of the motion. A special case is

<!-- eq:9 -->
$$H(\phi, \pi) \equiv \tfrac{1}{2}\pi^2 + S(\phi),$$

where $S$ is the action of eq. (1).

Our procedure for generating a new configuration $\phi'$ is:

1. Select initial momenta $\pi$ at random from a Gaussian distribution of mean zero and unit variance,

<!-- eq:10 -->
$$P_{G}(\pi) \propto \exp(-\pi^2/2),$$

2. Let the system evolve deterministically through $(\phi, \pi)$-phase space for a fixed time $\tau_0$ according to Hamilton's equations (8). If the trajectory $(\phi(\tau), \pi(\tau))$ solves Hamilton's equations, this evolution defines a mapping on phase space $(\phi, \pi) = (\phi(0), \pi(0)) \mapsto (\phi(\tau_0), \pi(\tau_0))$.

The probability $P_H$ of choosing the candidate phase-space configuration $(\phi', \pi')$ is thus a delta function (because the dynamics is completely deterministic; this is not necessary, as long as $P_H$ is reversible the algorithm remains valid):

<!-- eq:11 -->
$$P_{H}\big((\phi, \pi) \mapsto (\phi', \pi')\big) = \delta\big[(\phi', \pi') - (\phi(\tau_{0}), \pi(\tau_{0}))\big].$$

Finally we accept this candidate with probability (a slight generalization of eq. (6), since $P_A$ now depends on the conjugate momenta $\pi$ as well as $\phi$):

<!-- eq:12 -->
$$P_{A}\big((\phi,\pi)\mapsto(\phi',\pi')\big) = \min\!\big(1,\ \exp(\delta H)\big),$$

where $\delta H = -\big[H(\phi', \pi') - H(\phi, \pi)\big] = H(\phi,\pi) - H(\phi',\pi')$, $H$ being the special Hamiltonian of eq. (9). <!-- CHECK: OCR garbled "5H=-H(O',re')-H(O,re)"; sign fixed so that eq (15) identity holds (min(1,exp(-ΔH)) with ΔH=H_new-H_old). -->

The transition probability restricted to the $\phi$-field alone is

<!-- eq:13 -->
$$P_{M}(\phi \mapsto \phi') = \int [d\pi]\, [d\pi']\, P_{G}(\pi)\, P_{H}\big((\phi, \pi) \mapsto (\phi', \pi')\big)\, P_{A}\big((\phi, \pi) \mapsto (\phi', \pi')\big).$$

### Derivation (not implemented): detailed balance proof

We wish to show that (13) satisfies detailed balance (5) exactly. For this we require the dynamics to be reversible:

<!-- eq:14 -->
$$P_{H}\big((\phi,\pi) \mapsto (\phi',\pi')\big) = P_{H}\big((\phi',-\pi') \mapsto (\phi,-\pi)\big),$$

which holds for the Hamiltonian dynamics introduced above. From the identity

<!-- eq:15 -->
$$\exp[-H(\phi, \pi)] \min(1, \exp(-\delta H)) = \min\big(\exp[-H(\phi, \pi)],\ \exp[-H(\phi', \pi')]\big) = \exp[-H(\phi', \pi')] \min(\exp(\delta H), 1),$$

and observing that $P_S(\phi)P_G(\pi) \propto \exp[-H(\phi,\pi)]$ and $H$ is invariant under $\pi \mapsto -\pi$, we obtain

<!-- eq:16 -->
$$P_{S}(\phi)P_{G}(\pi)P_{A}\big((\phi,\pi)\mapsto(\phi',\pi')\big) = P_{S}(\phi')P_{G}(\pi')P_{A}\big((\phi',\pi')\mapsto(\phi,\pi)\big) = P_{S}(\phi')P_{G}(-\pi')P_{A}\big((\phi',-\pi')\mapsto(\phi,-\pi)\big).$$

Multiplying by $P_H$, integrating over $\pi$ and $\pi'$ and using the reversibility condition (14):

<!-- eq:17 -->
$$\int [d\pi][d\pi']\, P_{S}(\phi)P_{G}(\pi)P_{H}\big((\phi,\pi)\mapsto(\phi',\pi')\big)P_{A}\big((\phi,\pi)\mapsto(\phi',\pi')\big) = \int [d\pi][d\pi']\, P_{S}(\phi')P_{G}(-\pi')P_{H}\big((\phi',-\pi')\mapsto(\phi,-\pi)\big)P_{A}\big((\phi',-\pi')\mapsto(\phi,-\pi)\big),$$

which yields the detailed balance equation (5) from (13) and the invariance of the measure $[d\pi][d\pi'] = [d(-\pi)][d(-\pi')]$.

In the case $H' = H$ the dynamics conserves energy, $\delta H = 0$, so $P_A = 1$: this limit is the usual hybrid algorithm. But we have proved a more general result: the algorithm generates $\phi$-field configurations with the correct distribution for **any** $\delta H$ (i.e. for any guidance Hamiltonian $H' \neq H$).

## Leapfrog integrator

In practice we can only approximately integrate the equations of motion in discrete steps of duration $\delta\tau$. A simple scheme that ensures exact reversibility (14) and generates an area-preserving map on phase space for any value of $\delta\tau$ is the **leapfrog** algorithm. An initial half-step

<!-- eq:18 -->
$$\pi(\delta \tau/2) = \pi(0) - [\delta S(0)/\delta \phi]\, \delta \tau/2$$

is followed by $n = \tau_0/\delta\tau$ steps in $\phi$ and $n-1$ steps in $\pi$ of the form

<!-- eq:19 -->
$$\phi(\tau + \delta \tau) = \phi(\tau) + \pi(\tau + \delta \tau/2)\cdot \delta \tau,$$

<!-- eq:20 -->
$$\pi(\tau + \delta \tau/2) = \pi(\tau - \delta \tau/2) - [\delta S(\tau)/\delta \phi]\cdot \delta \tau,$$

and a final half-step

<!-- eq:21 -->
$$\pi(\tau_0) = \pi(\tau_0 - \delta \tau/2) - [\delta S(\tau_0)/\delta \phi]\cdot \delta \tau/2.$$

The half-steps differ from exact integration of Hamilton's equations by errors of order $\delta\tau^2$, whereas the intermediate steps have errors of order $\delta\tau^3$.

## Dynamical fermions

Dynamical fermions may be included as in the hybrid or molecular dynamics approaches. The Grassmann fields $\bar\psi$ and $\psi$ are replaced by bosonic fields $\chi^*$ and $\chi$ with non-local interactions:

<!-- eq:22 -->
$$P_S(\phi) = \frac{1}{Z} \int [d\bar{\psi}][d\psi]\, \exp[-S(\phi) - \bar{\psi}\mathcal{M}\psi] = \frac{1}{Z'} \det(\mathcal{M})\, \exp[-S(\phi)] = \frac{1}{Z''} \int [d\chi^*][d\chi]\, \exp[-S(\phi) - \chi^{*}(\mathcal{M}^{\dagger}\mathcal{M})^{-1}\chi].$$

The quadratic form $\chi^*(\mathcal{M}^\dagger\mathcal{M})^{-1}\chi$ is used to ensure convergence of the bosonic Gaussian integrals. The Hamiltonian and equations of motion become

<!-- eq:23 -->
$$H(\phi, \pi) = \tfrac{1}{2}\pi^2 + S(\phi) + \chi^{*}(\mathcal{M}^{\dagger}\mathcal{M})^{-1}\chi,$$

<!-- eq:24 -->
$$\dot{\phi} = \pi, \qquad \dot{\pi} = -\delta S/\delta \phi + \chi^{*}(\mathcal{M}^{\dagger}\mathcal{M})^{-1}\big[\mathcal{M}^{\dagger}\, \delta\mathcal{M}/\delta \phi + (\delta \mathcal{M}^{\dagger}/\delta \phi)\, \mathcal{M}\big](\mathcal{M}^{\dagger}\mathcal{M})^{-1}\chi.$$

The $\chi$ field is held fixed during the molecular dynamics steps and is updated by an exact heatbath, $\chi = \mathcal{M}^\dagger \eta$ for Gaussian noise $\eta$, in between. A single conjugate-gradient inversion is required for each molecular dynamics step.

## Numerical test: compact QED

We tested the algorithm on compact quantum electrodynamics, a lattice gauge theory model of interacting photons and electrons, using the standard Wilson action for the gauge field and the staggered fermion formulation for the electrons. The Hamiltonian is

<!-- eq:25 -->
$$H = \frac{1}{2} \sum_{x\mu} \pi_{\mu}^{2}(x) + S(U, \bar{\psi}, \psi),$$

where the action $S$ is

<!-- eq:26 -->
$$S(U, \bar{\psi}, \psi) = \beta \sum_{\substack{x\nu\mu \\ \nu>\mu}} \left[ 1 - \operatorname{Re} U_{\mu\nu}(x) \right] + \sum_{x\nu} \bar{\psi}(x)\, (-D^2 + m^2)^{-1}\, \psi(y).$$

The phases $U_{\mu\nu}(x)$ are the usual products of link variables around the elementary plaquettes of the lattice.

### Results (see checks.md)

Initial tests on the pure gauge theory (no fermions) measured the average plaquette on an $8^4$ lattice as a function of step size $\delta\tau$ for both the standard hybrid algorithm and the new hybrid Monte Carlo algorithm. Results agree with previous high-statistics data. The optimal step size for the new algorithm on this lattice with coupling $\beta = 0.97$ is roughly 0.1. The effective step size ($a\,\delta\tau$ where $a$ is the acceptance rate) is then 0.06, several times larger than the maximum "safe" hybrid step size of 0.01-0.02. The new algorithm does not require extrapolation to $\delta\tau = 0$.

Simulations on lattices of size $4^4$, $8^4$, $12^4$ show that the step size must be scaled as $\delta\tau \propto 1/L$ to keep a constant acceptance rate, where $L$ is the linear lattice size. This is consistent with the acceptance rate being proportional to $\exp(-L^2\delta\tau^2)$.

## Two Hamiltonians: acceptance vs guidance

In the hybrid Monte Carlo method the Hamiltonian serves two distinct roles:

1. The **acceptance** Hamiltonian $H$ defines the equilibrium distribution, entering the acceptance probability (12).
2. The **guidance** Hamiltonian $H'$ appears in the equations of motion (via $S'$, eq. 7-8).

Nothing requires these two Hamiltonians to be equal; the generalization $H \neq H'$ introduces scope for optimization of the acceptance rate. Discretization errors in the integration are in some cases equivalent to performing an exact computation for a theory with an unknown action $S''$ differing from $S$ by renormalizations. By adjusting parameters in the guidance Hamiltonian there should be an optimal $H'$ for which the equivalent action $S''$ for the discretized dynamics approaches the desired action $S$ in the acceptance Hamiltonian, maximizing the acceptance rate. Numerically, the acceptance rate as a function of $\beta - \beta'$ (couplings in $H$ and $H'$) shows a peak for $\beta \neq \beta'$, improving over $H = H'$.

The method is easily generalized to quantum chromodynamics with dynamical quark fields.
