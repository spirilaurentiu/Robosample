# Nonequilibrium candidate Monte Carlo is an efficient tool for equilibrium simulation

Jerome P. Nilmeier, Gavin E. Crooks, David D. L. Minh, John D. Chodera.
PNAS 108(45):E1009-E1018, 2011. DOI: 10.1073/pnas.1106094108.

> NOTE (published Corrections, PNAS 2012, doi:10.1073/pnas.1207617109): Equations 39,
> 41, and 42 appeared incorrectly in the original article. This file uses the corrected
> forms. In the corrected forms the systematic-force term carries a `+` sign (not `-`),
> and the reverse-noise Eq. 42 carries a leading `-1/sqrt(2)` factor. See `equations.md`.

## Abstract

Metropolis Monte Carlo simulation is a powerful tool for studying the equilibrium
properties of matter. In complex condensed-phase systems it is difficult to design Monte
Carlo moves with high acceptance probabilities that also rapidly sample uncorrelated
configurations. This paper introduces a class of moves based on nonequilibrium dynamics:
candidate configurations are generated through a finite-time process in which a system is
actively driven out of equilibrium, then accepted with criteria that preserve the
equilibrium distribution. The acceptance rule is similar to the Metropolis acceptance
probability but related to the nonequilibrium work rather than the instantaneous energy
difference. The method applies to sampling from a single thermodynamic state or a mixture
of thermodynamic states, and allows both coordinates and thermodynamic parameters to be
driven in nonequilibrium proposals. Driving some degrees of freedom while allowing others
to evolve naturally can greatly enhance acceptance probabilities and reduce structural
correlation times.

## Introduction and motivation

The goal is to construct efficient Markov chain Monte Carlo (MCMC) moves that both have
high acceptance rates and allow rapid transit through configuration space. The Metropolis
Monte Carlo sampling procedure is generalized by using nonequilibrium processes to
generate candidates. Moves that are efficient for an isolated part of a system but lead to
near-universal rejection in standard MC of dense mixtures can be converted to
nonequilibrium processes that generate candidates with higher acceptance. The acceptance
criterion is related to the nonequilibrium work rather than the potential energy
difference used in traditional MC moves.

Designing efficient moves requires balancing rapid traversal of phase space against
reasonable acceptance probabilities. An experimenter may possess physical insight about
one component (e.g. a biomolecule) permitting moves that would be efficient in the absence
of other components (e.g. solvent), but that encounter unfavorable interactions in their
presence. As an illustrative example, consider a bistable dimer: a pair of particles with
a potential having minima in compact or extended configurations separated by a high
barrier. In vacuum, an effective standard MC move instantaneously changes the
interparticle distance from compact to extended (or the reverse). When the dimer is
immersed in a dense solvent, this move meets near-universal rejection because solvent
molecules overlap with proposed configurations.

One remedy is to use a nonequilibrium process to generate candidates, allowing
unperturbed degrees of freedom to relax and maintaining a reasonable acceptance rate.
Applying the appropriate acceptance criterion for the final configuration preserves the
equilibrium distribution. For the bistable dimer in dense solvent, the extension (or
contraction) is carried out over a finite number of increments interleaved with standard
Metropolis Monte Carlo or molecular dynamics steps that let the solvent reorganize.

The idea of using nonequilibrium driven processes as MC moves has precedents in the
statistical and chemical literature (work-bias Monte Carlo; constant-pH schemes; driving
a subset of degrees of freedom with approximate acceptance criteria). This paper unifies
these ideas into a framework, nonequilibrium candidate Monte Carlo (NCMC), applicable to
single thermodynamic states (NVT, NpT, muVT) and to mixtures of thermodynamic states
(expanded ensembles). Nonequilibrium proposals may drive a subset of degrees of freedom,
the thermodynamic parameters, or both.

## Equilibrium and Expanded Thermodynamic Ensembles

For physical systems in equilibrium, the probability of observing a microstate is given by
the Boltzmann distribution (Eq. 1). Here $x \in \Gamma$ denotes a microstate (which may
include coordinates, momenta, and other dynamical variables such as simulation box
dimensions), $\lambda$ denotes a set of thermodynamic parameters defining a thermodynamic
state, and $Z_\lambda$ is the partition function.

The reduced potential $u_\lambda(x)$ depends on the thermodynamic ensemble. In an
isothermal-isobaric (NpT) ensemble it takes the form of Eq. 2, depending on the
Hamiltonian $H(x)$ (which may include an external biasing potential and is presumed
invariant under momentum inversion) and the system volume $V(x)$. The controllable
thermodynamic parameters $\lambda \equiv \{\beta, H, p\}$ include inverse temperature
$\beta$, the Hamiltonian, and external pressure $p$. Other conjugate variables can be
included or excluded to generate alternative ensembles.

To sample multiple thermodynamic states within one simulation, an expanded ensemble
(Eq. 3) specifies a joint distribution for $(x,\lambda)$ in a weighted mixture, where
$\omega_\lambda > 0$ is an externally imposed weight for state $\lambda$. The set of
allowed $\lambda \in \mathcal{G}$ may be discrete or continuous. If $\mathcal{G}$ has a
single value, one thermodynamic state is sampled and $\pi(x,\lambda) = \pi_\lambda(x)$.
States may correspond to temperatures (simulated tempering), alchemical states (simulated
scaling), or protonation states (constant-pH).

## Nonequilibrium Candidate Monte Carlo

At the start of an iteration, the current sample $(x^{(n)}, \lambda^{(n)})$, assumed drawn
from $\pi(x,\lambda)$, initializes a trajectory $(x_0,\lambda_0) = (x^{(n)},
\lambda^{(n)})$. A candidate $(x_T, \lambda_T)$ is proposed through a nonequilibrium
process in which a set of degrees of freedom and/or thermodynamic parameters may be driven
according to some protocol, selected with a probability dependent only on $(x_0,
\lambda_0)$. Even to sample from a single thermodynamic state $\lambda$, one may use a
protocol that transiently drives the thermodynamic parameters away from $\lambda$ and back
again. Finally an acceptance probability is computed and used to decide whether the next
sample is the candidate $(x_T, \lambda_T)$ or the momentum reversal of the initial sample
$(\tilde{x}^{(n)}, \lambda^{(n)})$.

An NCMC move begins by selecting a protocol $\Lambda$ with probability $P(\Lambda | x_0,
\lambda_0)$, such that there exists a reverse protocol $\tilde{\Lambda}$ with
$P(\Lambda | \tilde{x}_T, \lambda_T) > 0$. A protocol $\Lambda$ specifies $T$ perturbation
kernels $\alpha_t(x,y)$ and propagation kernels $K_t(x,y)$ arranged in an alternating
pattern $\Lambda \equiv \{\alpha_1, K_1, \alpha_2, K_2, \ldots, \alpha_T, K_T\}$. Both
$\alpha_t$ and $K_t$ are conditional probabilities of $y \in \Gamma$ given $x \in \Gamma$,
and must satisfy: if $p(x,y) > 0$ then $p(y,x) > 0$, for $p$ being either $\alpha_t$ or
$K_t$.

Each perturbation kernel $\alpha_t$ drives some or all degrees of freedom $x$
stochastically or deterministically (e.g. driving a torsion angle, an interatomic
distance, or the simulation cell volume). Each propagation kernel $K_t$ propagates some or
all coordinates at fixed $\lambda_t$ according to some form of MCMC or MD (Metropolis
Monte Carlo, velocity Verlet deterministic dynamics, or overdamped Langevin stochastic
dynamics), possibly depending on the time index $t$. Interleaving perturbation and
propagation lets energetically unfavorable interactions introduced by perturbation relax
during propagation, increasing acceptance relative to instantaneous perturbations.

A trajectory $X \equiv (x_0, x_1, \ldots, x_T)$ is generated from $x_0$ under protocol
$\Lambda$ via the scheme of Eq. 4: application of $\alpha_t$ to $x_{t-1}$ generates a
perturbed configuration $x_t^*$, then propagated by $K_t$ to $x_t$.

The reverse protocol $\tilde{\Lambda} \equiv \{K_T, \alpha_T, \ldots, K_0, \alpha_0\}$
reverses the order of perturbation and propagation, generating the time-reversed
trajectory $\tilde{X} \equiv \{\tilde{x}_T, \ldots, \tilde{x}_0\}$ where $\tilde{x}$
denotes $x$ with inverted momenta (Eq. 5).

To preserve the stationary distribution $\pi(x,\lambda)$, a strict pathwise form of
detailed balance is enforced (Eq. 6). The pathwise probabilities $\Pi$ are products of the
kernels (Eqs. 7, 8). Summing Eq. 6 over all trajectories starting with $x_0$ and ending
with $x_T$ recovers standard detailed balance (see Appendix).

Defining the ratio of proposal kernels (Eq. 9) and the ratio of propagation kernels as the
exponentiated forward/backward conditional path action difference (Eq. 10), and using
momentum invariance $\pi(x,\lambda) = \pi(\tilde{x},\lambda)$, gives the ratio of
acceptance probabilities (Eq. 11), the main result. Here $\Delta u(X|\Lambda) \equiv
u_T(x_T) - u_0(x_0)$ is the energy difference. Eq. 11 is general with respect to the choice
of perturbation and propagation protocol.

Many acceptance probabilities satisfy Eq. 11, including the Metropolis-Hastings criterion
(Eq. 12). After generating $(x_T, \lambda_T)$ and evaluating $A(X|\Lambda)$, draw a uniform
variate $U$. If $A(X|\Lambda) > U$ the candidate becomes the next value, $(x^{(n+1)},
\lambda^{(n+1)}) = (x_T, \lambda_T)$. Otherwise it is rejected, a momentum flip is
performed, and the next value becomes $(\tilde{x}_0, \lambda_0)$. Alternately, flip the
momentum upon acceptance ($(\tilde{x}_T, \lambda_T)$) and preserve it upon rejection
($(x_0, \lambda_0)$). The momentum flip cannot be ignored; it is necessary to preserve the
equilibrium distribution (see Appendix).

NCMC need not be used exclusively to sample $\pi(x,\lambda)$; it can be mixed with other
MCMC moves or with MD. For example, one may reinitialize velocities from the
Maxwell-Boltzmann distribution after each NCMC step; this is a Gibbs sampling MCMC move
using the marginal distribution for velocities.

## Perturbation Kernels

Many choices are available for $\alpha_t(x,y)$. Judicious selection allows nonequilibrium
proposals that carry a component of the system from one high-probability region to another
with high acceptance.

**Stochastically Driven Degrees of Freedom.** To drive a torsion angle $\phi$ (subtended
by four bonded atoms) stochastically to a new angle $\phi'$ (holding angles and bonds
fixed), draw $\phi'$ from the von Mises circular distribution centered on $\phi$ (Eq. 13),
with $I_0(\kappa)$ the modified Bessel function of order zero and $\kappa > 0$ a
dimensionless force constant. Because the perturbation is made in a non-Cartesian
coordinate, a Jacobian $J(\phi)$ must be included to compute $\alpha(x,y)$ in Cartesian
coordinates; the resulting kernel ratio (Eq. 14) equals 1 because the rotation about a
bond vector preserves Cartesian phase-space volume, so $J(\phi') = J(\phi) = 1$.

**Deterministically Driven Degrees of Freedom.** Instead of stochastic perturbation, drive
the torsion in small fixed increments $\Delta\phi$. Define an invertible map $\mathcal{M}$
taking $x \to y$ such that $y = \mathcal{M}x$ differs from $x$ only in rotating $\phi$ by
$\Delta\phi$. Choose $\Delta\phi$ from a distribution where $\pm\Delta\phi$ are equally
probable, and drive $\phi(x)$ from $\phi_0$ to $\phi_T = \phi_0 + \Delta\phi$ over $T$
steps in equal increments, so $\phi(x_t)$ is constrained to $\phi_t \equiv (1-t/T)\phi_0 +
(t/T)\phi_T$. Then $\alpha_t(x,y) = \delta(y - \mathcal{M}x) J(x)$, where $J(x)$ is the
Cartesian phase-space compression factor, again unity for rotation about a torsion, and by
invertibility the ratio $\alpha_t(\tilde{y},\tilde{x})/\alpha_t(x,y) = 1$.

**Simulation Box Scaling.** A barostat scales molecular centers and box geometry by a
factor $s = [(V(x) + \Delta V)/V(x)]^{1/3}$ with $\Delta V$ chosen uniformly from
$[V - \Delta V_0, V + \Delta V_0]$, applied as a factor of $s^{1/T}$ over $T$ steps. The
perturbation kernel is a delta function; the ratio of perturbation kernels is
$\alpha(\tilde{X}|\tilde{\Lambda})/\alpha(X|\Lambda) = s^{3N}$ where $N$ is the number of
molecular centers.

**Thermodynamic Perturbation.** In many driven processes there is no direct perturbation
to coordinates, so $\alpha_t(x,y) = \delta(x-y)$ and
$\alpha(\tilde{X}|\tilde{\Lambda})/\alpha(X|\Lambda) = 1$. Only the thermodynamic
parameters $\lambda$ vary in time, carrying the system out of equilibrium through $K_t$. If
$u_t$ is a linear interpolation $u_t(x) = (1-t/T)u_0(x) + (t/T)u_T(x)$, the protocol
probability is symmetric, and MC is used for $K_t$, this recovers Neal's method.

## Propagation Kernels

The choice of propagation kernels is broad. With strong driving in $\alpha$, one may use a
time-independent kernel $K_t(x,y) \equiv K(x,y)$ sampling from a stationary distribution
$\pi(x)$. Alternatively a strongly time-dependent $K_t$ can transiently drive the system
out of equilibrium.

**Reversible Markov Chain Monte Carlo.** Propagate some or all degrees of freedom (e.g.
those not affected by $\alpha_t$) by a method satisfying detailed balance in $\pi_t$
(Eq. 15), where $\pi_t(x) \equiv Z_t^{-1} e^{-u_t(x)}$. Metropolis and hybrid Monte Carlo
(HMC) algorithms obey detailed balance.

By analogy with Crooks, define a work $w$ and heat $q$ for the nonequilibrium process
(Eqs. 16, 17) such that $w(X|\Lambda) + q(X|\Lambda) = \Delta u(X|\Lambda)$, a restatement
of the first law. The conditional path action difference is then written in terms of the
heat (Eq. 18), leading to an acceptance ratio like standard MC but with the work
$w(X|\Lambda)$ replacing the instantaneous potential energy difference (Eq. 19).

**Deterministic Dynamics.** When an isolated system is propagated by a symplectic
integrator (reversible, deterministic, phase-space-volume preserving), the propagation
kernels satisfy $K_t(x,y) = K_t(\tilde{y}, \tilde{x})$. Hence $\Delta\mathcal{S}(X|\Lambda)
= 0$ and the acceptance ratio reduces to Eq. 20. The equivalence of work and energy
difference for volume-preserving integrators is known from fluctuation-theorem
calculations. Symplectic integrators include velocity Verlet, and remain symplectic under
constraints (e.g. RATTLE) provided constraints are iterated to convergence each timestep.

**Stochastic Dynamics.** Stochastic integrators such as velocity Verlet discretizations of
Langevin dynamics sample a modified distribution differing from $\pi_t$ in a
timestep-dependent manner. Computing the relative action $\Delta\mathcal{S}(X|\Lambda)$ is
straightforward, and the NCMC acceptance criterion ensures NCMC-sampled configurations are
distributed according to the desired equilibrium ensemble. The Appendix computes
$\Delta\mathcal{S}(X|\Lambda)$ for the overdamped Langevin (Ermak-Yeh) integrator and the
Brünger-Brooks-Karplus (BBK) Langevin integrator.

## Illustrative Application: Bistable Dimer in a WCA Fluid

Simulations of a bistable dimer were run in vacuum and in a dense fluid. The dimer is a
pair of bonded particles interacting via a double-well potential with minima at $r = r_0$
(compact) and $r = 2r_0$ (extended), separated by a $5 k_B T$ barrier. In solvated
simulations the dimer is immersed in a dense bath (reduced density $\rho\sigma^3 = 0.96$)
of particles interacting via the Weeks-Chandler-Andersen (WCA) soft repulsive potential.
Each simulation iteration consisted of velocity reassignment from Maxwell-Boltzmann, 500
steps of generalized hybrid Monte Carlo (GHMC) dynamics (a Metropolis-corrected form of
Langevin dynamics, referred to as MD), optionally followed by an instantaneous MC move or
an NCMC move.

The rate of generating uncorrelated samples is quantified by the correlation time $\tau$
for the dimer extension $r(t)$. This is the asymptotic decay time of the correlation
function $C(t) = \langle r(0) r(t)\rangle$ (Eq. 21) with $C_0 = \langle r^2\rangle$ and
$C_\infty = \langle r\rangle^2$. The correlation time relates to the statistical
inefficiency $g = 1 + 2\tau$, the number of iterations needed to generate an effectively
uncorrelated sample.

For MD in vacuum, slow hopping gives $\tau = 59.2$ iterations. Introducing instantaneous
MC extension/contraction moves ($\Delta r = \pm r_0$) reduces $\tau \approx 0.0$. When the
dimer is immersed in dense WCA fluid, 500 MD steps alone give extremely slow barrier
crossings, $g \approx 600$ iterations per uncorrelated sample. Instantaneous MC moves do
not significantly reduce $\tau$ in solvent. Performing the same expansion/contraction over
2,048-step NCMC moves reduces the correlation time to $\tau = 4.0$ iterations. Each
iteration requires a fivefold increase in effort (500 MD + 2,048 NCMC switching = 2,548
force evaluations versus 500 for MD alone), but a 67-fold reduction in correlation time
yields an order-of-magnitude gain in overall efficiency.

The NCMC switching length is a free parameter. Instantaneous MC proposals of $\pm r_0$ are
accepted with probability $\approx 10^{-27}$. Dividing the move into 1-8 steps gives
little increase; 16-1,024 steps give a superlinear boost; useful levels are reached around
2,000 steps (12% acceptance at 2,048 steps, 38% at 8,192 steps).

Under assumptions relevant to the bistable dimer, the NCMC acceptance probability links to
$\tau_{\rm eff}$ (Eq. 22), where the NCMC correlation time is estimated from the average
acceptance probability $\gamma$ via $\tau_{\rm NCMC} \approx -1/\ln(1 - 2\gamma)$ (see
Appendix). The effective correlation time $\tau_{\rm eff}$ is only diminished when the NCMC
acceptance is large enough that $\tau_{\rm NCMC} \approx \tau_{\rm MD}$ (about 256
switching steps or more).

Efficiency relative to MD alone is the efficiency gain $E$ (Eq. 23). There is a slight
efficiency loss at short switching times (minimum 86.9% of MD-alone efficiency at 128
steps), then a rapid gain plateauing at approximately 13x for 2,048-4,096-step NCMC
proposals. Longer switching times (8,192) do not achieve as high a gain because the
correlation-time reduction does not offset the additional cost.

## Epilogue

NCMC uses nonequilibrium proposals within MCMC to enhance acceptance rates and improve
statistical efficiency. A straightforward approach is to borrow Metropolis Monte Carlo
proposals reasonable for one component of a system in isolation and convert them to
nonequilibrium proposals. Selecting efficient nonequilibrium proposals resembles choosing
good reaction coordinates: drive the system along slow collective coordinates where
orthogonal degrees of freedom relax quickly. Switching trajectories contain potentially
useful information; it is straightforward to incorporate information from rejected NCMC
proposals in the estimation of equilibrium averages.

## Materials and Methods

The dimer consists of two particles interacting via a double-well bonded potential in the
interatomic distance $r$ (Eq. 24), with $h = 5 k_B T$, $r_0 = r_{\rm WCA}$, and
$s = r_{\rm WCA}/2$, where $r_{\rm WCA} \equiv 2^{1/6}\sigma$. "Vacuum" simulations contain
only these two particles; "solvated" simulations add a dense bath interacting via the WCA
nonbonded potential (Eq. 25), with mass $m = 39.9$ amu, $\sigma = 3.4$ Angstrom, and
$\varepsilon = 120\, k_B T$. The nonbonded WCA interaction is excluded between the two
bonded particles. The solvated system contains 216 WCA particles at reduced density
$\rho\sigma^3 = 0.96$. For all simulations the reduced temperature is $k_B T/\varepsilon =
0.824$. Simulations used a custom Python code with the GPU-accelerated OpenMM package and
the PyOpenMM wrapper.

To ensure differences were not due to changes in the integrator's stationary distribution,
GHMC was used for all simulations. GHMC is based on a velocity Verlet discretization of
Langevin dynamics (equivalent in the small-timestep limit) but includes an
acceptance/rejection step correcting finite-timestep errors so the stationary distribution
is exact. The timestep was $0.002\tau$ where $\tau = \sqrt{\sigma^2 m/\varepsilon}$, and
the collision rate was $\tau^{-1}$. The GHMC acceptance probability is $99.929 \pm
0.001\%$.

For instantaneous Monte Carlo moves, a perturbation $\Delta r$ to the interatomic distance
was chosen by Eq. 26. The dimer was contracted or expanded about the bond midpoint to
generate $x_{\rm new}$ with extension $r_{\rm new}$ from $x_{\rm old}$ with extension
$r_{\rm old}$, and accepted or rejected with the Metropolis-Hastings criterion (Eq. 27),
where the Jacobian ratio $J_r(x_{\rm old}, x_{\rm new}) = (r_{\rm new}/r_{\rm old})^2$
accounts for expansion/contraction of phase space.

For $T$-step NCMC moves, proposals were made by selecting a new velocity vector from
Maxwell-Boltzmann, integrating $T$ steps of velocity Verlet dynamics for all bath atoms as
the dimer extension was driven from $r_{\rm old}$ to $r_{\rm new}$ in equal steps of size
$\Delta r/T$, and accepting or rejecting based on the modified Metropolis criterion for
symplectic integrators (Eq. 28). The Jacobian ratio is again $(r_{\rm new}/r_{\rm
old})^2$. The corresponding perturbation kernel is Eq. 29. The propagation kernel
$K_t(x,y)$ is velocity Verlet dynamics with the dimer atoms held fixed in space.

The mean acceptance probability for each switching time is estimated by the sample mean
(Eq. 30). For numerical stability, logarithms $a_n \equiv \ln A(X_n)$ were stored, and
$\ln\langle A\rangle_\tau$ was estimated by Eq. 31 with $b \equiv \max_n a_n$. Integrated
autocorrelation times were estimated using the rapid scheme of ref. 42.

The acceptance probabilities plotted in Fig. 4 were estimated from 10,000 iterations of
2,048-step NCMC, with 500 GHMC steps between each NCMC trial. Statistical error was
estimated by 1,000 bootstrap trials.

The reference distribution $\mathcal{P}(r)$ was computed analytically for the vacuum
system (Eq. 32). For the solvated system it was estimated from an umbrella sampling
simulation with a modified bonded potential removing the barrier (Eq. 33), where
$r_{\min} = r_0$, $r_{\max} = 2.05 r_0$, $K = k_B T/\eta^2$ with $\eta = 0.3$ Angstrom, and
$\theta$ is the Heaviside function. The true solvated distribution $\mathcal{P}(r)$ was
estimated by reweighting (Eq. 34).

## Appendix (derivations, not implemented)

### Derivation (not implemented): Proof that NCMC preserves the equilibrium distribution

Following the GHMC proof, define the expected acceptance rate for NCMC moves from
$(x,\lambda)$ (Eq. 35). Given a variate $(x^{(n)},\lambda^{(n)})$ drawn from
$\pi(x,\lambda)$, the next value's density has a rejection contribution (Eq. 36, momentum
flip $\to (\tilde{x}^{(n)}, \lambda^{(n)})$) and an acceptance contribution (Eq. 37,
$\to (x_T, \lambda_T)$), where $\rho(X,\Lambda|x_0,\lambda_0) \equiv
\Pi(X|x_0,\Lambda)P(\Lambda|x_0,\lambda_0)$ and the pathwise detailed balance condition
(Eq. 6) is used. Summing Eqs. 36 and 37 shows the equilibrium distribution is preserved
(Eq. 38). Maintaining momentum on rejection and flipping it on acceptance preserves the
distribution equally.

### Derivation (not implemented): Ermak-Yeh overdamped Langevin acceptance criterion

The Ermak-Yeh overdamped (Brownian) integrator propagates only coordinates $x$ (Eq. 39,
corrected). $F_t(x) \equiv -(\partial/\partial x) H_t(x)$ is the systematic force, $\gamma$
an effective collision frequency (inverse time), $m$ the particle mass. The noise
$\xi_t$ is a normal variate with zero mean and variance $\beta^{-1}$ (Eq. 40). For each
step $x_t^* \to x_t$ there is a reverse noise $\tilde{\xi}_t$ generating $x_t \to x_t^*$
(Eqs. 41, 42, corrected). The conditional path action difference (Eq. 43) follows because
the Jacobians $|\partial x_t/\partial\xi_t|$ cancel and tildes drop (no momenta).

### Derivation (not implemented): BBK Langevin acceptance criterion

The velocity Verlet discretization of the Brünger-Brooks-Karplus (BBK) integrator (Eq. 44)
uses an auxiliary velocity $v_t'$ and requires two random variates $\xi_t, \xi_t'$ per
degree of freedom per timestep, with joint distribution Eq. 45. The reverse noise
variables generating $(r_t, -v_t) \to (r_t^*, -v_t^*)$ are Eq. 46. The Jacobian $J(\xi_t,
\xi_t')$ (Eq. 47) is independent of the noise variates, so the conditional path action
difference (Eq. 48) is the difference of squared noise magnitudes; the Jacobian ratio
cancels.

### Derivation (not implemented): Effective correlation time for mixed MD/NCMC sampling

Assume two long-lived conformational states of equal population with extensions $r_c$ and
$r_e$. The MD iteration transition is the column-stochastic matrix $T_{\rm MD}$ (Eq. 49).
The autocorrelation function for a $2\times 2$ column-stochastic $T$ (Eq. 50) gives
$C(n\Delta t) = (C_0 - C_\infty)\mu^n + C_\infty$ with $C_0 = (1/2)(r_c^2 + r_e^2)$,
$C_\infty = (1/4)(r_c + r_e)^2$, and $\mu$ the nonunit eigenvalue. Fitting Eq. 51 gives
$\tau = -1/\ln\mu$, hence $\tau_{\rm MD} = -1/\ln\mu_{\rm MD}$ with $\mu_{\rm MD} = 1 -
2\alpha$. The NCMC transition matrix is Eq. 52 with $\tau_{\rm NCMC} = -1/\ln\mu_{\rm
NCMC}$, $\mu_{\rm NCMC} = 1 - 2\gamma$. The effective transition matrix $T_{\rm eff} =
T_{\rm MD} T_{\rm NCMC}$ (Eq. 53) gives $\mu_{\rm eff} = 1 - 2[(\alpha+\gamma) -
2\alpha\gamma]$. Substituting $\alpha = (1 - e^{-1/\tau_{\rm MD}})/2$ and $\gamma =
(1 - e^{-1/\tau_{\rm NCMC}})/2$ yields $\tau_{\rm eff} = \tau_{\rm MD}\tau_{\rm NCMC}/
(\tau_{\rm MD} + \tau_{\rm NCMC})$ (Eq. 54). Check: for 2,048-step NCMC, $\tau_{\rm eff}
\approx 4.0$ from $\tau_{\rm MD} = 299.8$ and $\gamma = 12.1\%$; the measured value is
$\tau_{\rm eff} = 4.0$.
