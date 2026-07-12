# Nonequilibrium Measurements of Free Energy Differences for Microscopically Reversible Markovian Systems

Gavin E. Crooks (1998), J. Stat. Phys. 90(5/6):1481-1487. DOI 10.1023/A:1023208217925.

## Abstract

An equality (the Jarzynski equality) relates the free energy difference between
two equilibrium ensembles of a system to an ensemble average of the work
required to switch between these two configurations. This paper shows that the
result follows from two assumptions only: the system's dynamics is Markovian and
microscopically reversible (obeys detailed balance).

## 1. Introduction

Consider a classical system in contact with a constant-temperature heat bath
where some controllable degree of freedom is manipulated. Let $\lambda$ be a
parameter specifying the current value of that degree of freedom. $\lambda$ is
switched between an initial and final value over a finite length of time. The
free energy difference $\Delta F$ between the equilibrium ensembles corresponding
to the initial and final values of $\lambda$ relates to an average of the work
$W$ expended during switching:

<!-- eq:1 -->
$$ \overline{e^{-\beta W}} = e^{-\beta \Delta F} $$

where $\beta = 1/k_B T$, $k_B$ is Boltzmann's constant, and $T$ is the heat-bath
temperature. The overbar indicates an ensemble average over all possible paths
through phase space, given an equilibrium initial state and a fixed protocol of
the control parameter. It can be approximated by an average over many independent
measurements of the switching work. This is a nonequilibrium measurement: the
system need not be in equilibrium during switching.

Equation (1) generalizes several standard free-energy relations.

In the limit of an infinitely long (reversible) switching process, Eq. (1)
becomes thermodynamic integration: the system is always in equilibrium, no energy
is dissipated, and $\Delta F = W$. Finite switching time causes systematic
overestimation of $\Delta F$.

In the limit of infinitely fast switching, Eq. (1) reduces to the
Zwanzig thermodynamic-perturbation relation:

<!-- eq:2 -->
$$ \langle e^{-\beta W} \rangle_0 = e^{-\beta \Delta F} $$

The angled brackets indicate an equilibrium average with $\lambda$ fixed at its
initial value. Here $W$ is the work to instantaneously change the control
parameter while holding the system configuration fixed. To reduce statistical
error the control parameter is changed in small steps, equilibrating at each
step; finite relaxation time again produces systematic error.

Both thermodynamic integration and perturbation suffer systematic error because
the simulated systems are never truly in equilibrium. Equation (1) requires only
that the initial ensemble is in equilibrium, reducing systematic error at the
cost of increased statistical error.

This relation was previously derived for a Hamiltonian system weakly coupled to a
heat bath, and via a master-equation approach. This paper shows the equality
follows directly if the dynamics are Markovian (memoryless) and microscopically
reversible (detailed balance). Detailed balance ensures time reversibility and
that equilibrium distributions are canonical. Both conditions hold for many
computer simulations.

## 2. Notation and Assumptions

The system is held at constant temperature $T$. Its state is specified by two
parameters. The internal state at time $t$ is labeled $i_t$. The energy of each
state is affected by an externally controlled parameter $\lambda_t$. The energy
at time $t$ is $E(i_t, \lambda_t)$.

Example 1: a classical gas in a cylinder. $i_t$ is the vector of instantaneous
momenta and positions of all particles; $\lambda$ controls the volume via a
piston. Time evolution could be Langevin dynamics (deterministic equations of
motion plus frictional and stochastic forces).

Example 2: an Ising ferromagnet. $i_t$ indexes the up/down positions of all
spins; $\lambda_t$ is the external magnetic field. Time evolution is single-spin-
flip Metropolis Monte Carlo.

For convenience assume discrete time and discrete phase space (true for any
digital simulation). The results generalize to continuous time and phase space.

In a canonical ensemble the equilibrium probability of state $A$ given a fixed
control parameter $\lambda$ is

<!-- eq:3 -->
$$ P(A \mid \lambda) = \frac{e^{-\beta E(A, \lambda)}}{\sum_{i} e^{-\beta E(i, \lambda)}} = \exp[\,\beta F(\beta, \lambda) - \beta E(A, \lambda)\,] $$

The sum is over all states. $F(\beta, \lambda) = -\beta^{-1} \ln \sum_i e^{-\beta E(i, \lambda)}$ is the Helmholtz free energy.

We consider the evolution of the system as the control parameter moves through a
fixed sequence $\{\lambda_0, \lambda_1, \dots, \lambda_\tau\}$. A path through
phase space is written

<!-- eq:4 -->
$$ i_0 \xrightarrow{\lambda_1} i_1 \xrightarrow{\lambda_2} i_2 \xrightarrow{\lambda_3} \cdots \xrightarrow{\lambda_\tau} i_\tau $$

At $t = 0$ the system is in state $i_0$ and the control parameter is $\lambda_0$.
Each time step occurs in two substeps. First the control parameter moves to a new
value $\lambda_{t+1}$; this takes work $E(i_t, \lambda_{t+1}) - E(i_t, \lambda_t)$.
Then the state evolves, at constant $\lambda_{t+1}$, to state $i_{t+1}$; during
this evolution the system exchanges heat $E(i_{t+1}, \lambda_{t+1}) - E(i_t, \lambda_{t+1})$
with the reservoir. This repeats for $\tau$ time steps.

The total work $W$, total heat $Q$, and total energy change $\Delta E$ are

<!-- eq:5 -->
$$ W = \sum_{t=0}^{\tau-1} [\,E(i_t, \lambda_{t+1}) - E(i_t, \lambda_t)\,] $$

<!-- eq:6 -->
$$ Q = \sum_{t=1}^{\tau} [\,E(i_t, \lambda_t) - E(i_{t-1}, \lambda_t)\,] $$

<!-- eq:7 -->
$$ \Delta E = Q + W = E(i_\tau, \lambda_\tau) - E(i_0, \lambda_0) $$

The reversible work $W_r = \Delta F = F(\beta, \lambda_\tau) - F(\beta, \lambda_0)$
is the free energy difference between the two equilibrium ensembles. The
dissipative work $W_d = W - W_r$ is the difference between actual and reversible
work. If work is expended in changing the free energy, the change in entropy of
the universe is $\beta W_d$ (in units of $k_B$). Work and dissipative work depend
on the phase-space path; the reversible work depends only on the initial and
final ensembles.

Reversing the direction of time, the reverse-time path corresponding to Eq. (4)
is $i_0 \xleftarrow{\lambda_1} i_1 \xleftarrow{\lambda_2} i_2 \cdots \xleftarrow{\lambda_\tau} i_\tau$.
The order in which states are visited is reversed, as is the order in which
$\lambda$ changes. The forward path begins with a change in $\lambda$; the reverse
path begins with a change in the internal state. Work, heat, energy change, and
free energy change are defined in the forward direction; in the reverse direction
they are the negatives of the forward values.

If the evolution is Markovian, the transition probability $P(i_t \xrightarrow{\lambda} i_{t+1})$
depends only on the state at time $t$, not on prior history. The probability of a
path (given the initial state $i_0$ and the control parameter at all times) factors
into single-step parts:

<!-- eq:path-factor -->
$$ P(i_0 \xrightarrow{\lambda_1} i_1 \xrightarrow{\lambda_2} i_2 \cdots \xrightarrow{\lambda_\tau} i_\tau) = P(i_0 \xrightarrow{\lambda_1} i_1)\, P(i_1 \xrightarrow{\lambda_2} i_2) \cdots P(i_{\tau-1} \xrightarrow{\lambda_\tau} i_\tau) $$

The single time steps are microscopically reversible and obey detailed balance
for every fixed value of the control parameter $\lambda$:

<!-- eq:8 -->
$$ \frac{P(A \xrightarrow{\lambda} B)}{P(A \xleftarrow{\lambda} B)} = \frac{P(B \mid \lambda)}{P(A \mid \lambda)} = \frac{e^{-\beta E(B, \lambda)}}{e^{-\beta E(A, \lambda)}} $$

An analogous detailed-balance condition holds for a multi-step process in which
an arbitrary amount of work is performed. Given Markovian dynamics and single-step
detailed balance, the ratio of the probability of a forward path to the
corresponding time-reversed path is

<!-- eq:9 -->
$$ \frac{P(i_0 \xrightarrow{\lambda_1} i_1 \cdots \xrightarrow{\lambda_\tau} i_\tau)}{P(i_0 \xleftarrow{\lambda_1} i_1 \cdots \xleftarrow{\lambda_\tau} i_\tau)} = \frac{e^{-\beta E(i_1, \lambda_1)} e^{-\beta E(i_2, \lambda_2)} \cdots e^{-\beta E(i_\tau, \lambda_\tau)}}{e^{-\beta E(i_0, \lambda_1)} e^{-\beta E(i_1, \lambda_2)} \cdots e^{-\beta E(i_{\tau-1}, \lambda_\tau)}} = e^{-\beta Q} $$

Here $Q$ is the energy exchanged with the heat bath along the forward path
(Eq. 6), and $-\beta Q$ is the corresponding entropy change of the bath (in units
of $k_B$). This states that detailed balance continues to hold for Markovian,
microscopically reversible systems regardless of how much work is performed.

Specifying that both forward and reverse paths start from equilibrium
distributions gives

<!-- eq:10 -->
$$ \frac{P(i_0 \mid \lambda_0)\, P(i_0 \xrightarrow{\lambda_1} i_1 \cdots \xrightarrow{\lambda_\tau} i_\tau)}{P(i_\tau \mid \lambda_\tau)\, P(i_0 \xleftarrow{\lambda_1} i_1 \cdots \xleftarrow{\lambda_\tau} i_\tau)} = e^{\beta \Delta E - \beta \Delta F}\, e^{-\beta Q} = e^{\beta W_d} $$

## 3. Derivation

### Derivation (not implemented)

These consequences of detailed balance and the Markov condition give a simple
proof that $\overline{e^{-\beta W}} = e^{-\beta \Delta F}$ (Eq. 1). The overbar
averages over all phase-space paths, given a fixed control-parameter protocol and
a canonical equilibrium initial distribution:

<!-- eq:11 -->
$$ \overline{e^{-\beta W}} = \sum_{i_0, i_1, \dots, i_\tau} P(i_0 \mid \lambda_0)\, P(i_0 \xrightarrow{\lambda_1} i_1 \cdots \xrightarrow{\lambda_\tau} i_\tau)\, e^{-\beta W} $$

The forward-path average is converted to a reverse-path average using Eq. (10):

<!-- eq:12 -->
$$ \overline{e^{-\beta W}} = \sum_{i_0, i_1, \dots, i_\tau} P(i_\tau \mid \lambda_\tau)\, P(i_0 \xleftarrow{\lambda_1} i_1 \cdots \xleftarrow{\lambda_\tau} i_\tau)\, e^{\beta W_d - \beta W} = e^{-\beta \Delta F} $$

The last step holds because the reversible work $\Delta F = W_r = W - W_d$ is path
independent and because the reverse-path probabilities are normalized.

## Generalization to the isothermal-isobaric ensemble

The result generalizes to other isothermal ensembles. For a classical gas at
constant temperature and pressure $p$, work performed by slow-growth particle
insertion gives $\Delta G = G(\beta, p, \lambda_\tau) - G(\beta, p, \lambda_0)$
equal to the excess chemical potential of the inserted particle. The detailed
balance condition becomes

<!-- eq:13 -->
$$ \frac{P(A \xrightarrow{\lambda} B)}{P(A \xleftarrow{\lambda} B)} = \frac{e^{-\beta E(B, \lambda) - \beta p V(B, \lambda)}}{e^{-\beta E(A, \lambda) - \beta p V(A, \lambda)}} = e^{-\beta Q - \beta p \Delta V} $$

where $V$ is the volume and $\Delta V$ is a baric equivalent of the thermal heat.
The ratio of equilibrium probabilities is
$P(i_\tau \mid \lambda_\tau)/P(i_0 \mid \lambda_0) = \exp(\beta \Delta E + \beta p \Delta V - \beta \Delta G)$.
Equation (10) is unmodified and holds in both the canonical and isothermal-
isobaric ensembles, so $\overline{e^{-\beta W}} = e^{-\beta \Delta G}$.

<!-- CHECK: raw OCR eq (13) had e^{-beta p V(B,lambda)} in the denominator; corrected to V(A,lambda) so the ratio yields e^{-beta p Delta V} with Delta V = V(B)-V(A). -->

## 4. Summary

Free energy differences of isothermal systems relate directly to a nonequilibrium
exponential average of the work required to switch between ensembles (Eq. 1). The
proof (Eqs. 11-12) rests on the assumptions that the system is Markovian and
microscopically reversible (Eqs. 8-10). Thermodynamic integration and perturbation
are limiting cases. Explicit use of Eq. (1) reduces systematic error but may be
limited by increased statistical error.
