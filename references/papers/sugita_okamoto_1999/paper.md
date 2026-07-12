# Replica-exchange molecular dynamics method for protein folding

Yuji Sugita, Yuko Okamoto. Chemical Physics Letters 314 (1999) 141-151.

## Abstract

We have developed a formulation for a molecular dynamics algorithm for the
replica-exchange method. The effectiveness of the method for the protein-folding
problem is tested with the penta-peptide Met-enkephalin. The method can overcome
the multiple-minima problem by exchanging non-interacting replicas of the system
at several temperatures. From only one simulation run, one can obtain probability
distributions in the canonical ensemble for a wide temperature range using
multiple-histogram reweighting techniques, which allows the calculation of any
thermodynamic quantity as a function of temperature in that range.

## 1. Introduction

In protein-folding simulations, it is usually difficult to obtain accurate
canonical distributions at low temperatures by conventional simulation methods
because simulations at low temperatures tend to get trapped in one of a huge
number of local minimum-energy states. One way to overcome this multiple-minima
problem is to perform a simulation based on non-Boltzmann probability weight
factors so that a random walk in energy space may be realized. Random walks allow
the simulation to pass any energy barrier and to sample a much wider phase space
than by conventional methods.

Well-known generalized-ensemble methods include the *multicanonical algorithm*
(a free 1D random walk in energy space) and *simulated tempering* (a free random
walk in temperature space, which induces a random walk in energy space). These
methods perform random walks in energy space due to non-Boltzmann weight factors
and are given the generic name *generalized-ensemble algorithm*.

The generalized-ensemble method is powerful but, in the above two methods, the
probability weight factors are not a priori known and have to be determined by
iterations of short trial simulations. This process can be non-trivial and tedious.
In the present work we develop a molecular dynamics (MD) algorithm based on a new
generalized-ensemble algorithm, the *replica-exchange method* (REM). (The method
is also referred to as *replica Monte Carlo*, *multiple Markov chain method*, and
*parallel tempering*.) In this method the weight factor is essentially known and
there is no complication in its determination. The Monte Carlo (MC) algorithm (and
MD algorithm in dihedral space) in this generalized ensemble has been applied to
oligopeptide systems. Details for the MD algorithm in Cartesian coordinates had
yet to be worked out, and this is the purpose of the present Letter. The
performance of the new algorithm is tested with the penta-peptide Met-enkephalin
in gas phase.

## 2. Methods

Consider a system of $N$ atoms of mass $m_k$ ($k = 1, \dots, N$) with coordinate
vectors $q = \{q_1, \dots, q_N\}$ and momentum vectors $p = \{p_1, \dots, p_N\}$.
The Hamiltonian $H(q,p)$ is the sum of the kinetic energy $K(p)$ and the potential
energy $E(q)$:

$$ H(q,p) = K(p) + E(q), $$

where

$$ K(p) = \sum_{k=1}^{N} \frac{p_k^2}{2 m_k}. $$

In the canonical ensemble at temperature $T$, each state $x \equiv (q,p)$ with
Hamiltonian $H(q,p)$ is weighted by the Boltzmann factor
$W_B(x;T) = e^{-\beta H(q,p)}$, where $\beta = 1/k_B T$. The average kinetic energy
at temperature $T$ is $\langle K(p)\rangle_T = \frac{3}{2} N k_B T$.

In the original version of the replica-exchange method (REM), the MC algorithm was
used. Here we describe the method in the context of the MD algorithm.

The generalized ensemble for REM consists of $M$ non-interacting copies (replicas)
of the original system in the canonical ensemble at $M$ different temperatures
$T_m$ ($m = 1, \dots, M$). We arrange the replicas so that there is always exactly
one replica at each temperature. There is then a one-to-one correspondence between
replicas and temperatures; the label $i$ for replicas is a permutation of the
label $m$ for temperatures, and vice versa, via the permutation functions
$i = i(m) \equiv f(m)$ and $m = m(i) \equiv f^{-1}(i)$.

Let $X = (x_1^{[i(1)]}, \dots, x_M^{[i(M)]}) = (x_{m(1)}^{[1]}, \dots,
x_{m(M)}^{[M]})$ stand for a 'state' in this generalized ensemble. The superscript
labels the replica and the subscript labels the temperature in $x_m^{[i]}$. The
state $X$ is specified by the $M$ sets of coordinates $q^{[i]}$ and momenta
$p^{[i]}$ of $N$ atoms in replica $i$ at temperature $T_m$: $x_m^{[i]} =
(q^{[i]}, p^{[i]})_m$.

Because the replicas are non-interacting, the weight factor for the state $X$ is
the product of Boltzmann factors for each replica (see eq:7).

We now consider exchanging a pair of replicas $i$ and $j$ at temperatures $T_m$
and $T_n$. In the original MC implementation of REM, only the coordinates $q$
(and the potential energy $E(q)$) had to be taken into account. Here, in the MD
implementation, we also have to deal with the momenta $p$. We propose a momentum
assignment (see eq:12) in which we rescale uniformly the velocities of all atoms
in a replica by the square root of the ratio of the two temperatures, so that the
temperature condition on the average kinetic energy is satisfied. We believe this
is the simplest and most natural choice.

In order for the exchange process to converge towards an equilibrium distribution,
it is sufficient to impose the detailed balance condition on the transition
probability $w(X \to X')$ (see eq:13). Working through the algebra gives the
acceptance ratio in terms of a single quantity $\Delta$ (see eq:14, eq:15), which
depends only on the potential energies and inverse temperatures. This can be
satisfied by the usual Metropolis criterion (see eq:17). Note that this is exactly
the same criterion originally derived for the MC algorithm.

### Derivation (not implemented)

From the Boltzmann weight (eq:7), the detailed balance condition (eq:13), and the
momentum rescaling (eq:12), the ratio of transition probabilities is:

$$ \frac{w(X \to X')}{w(X' \to X)} = \exp\Big\{-\beta_m [K(p^{[j]'}) + E(q^{[j]})]
- \beta_n [K(p^{[i]'}) + E(q^{[i]})] + \beta_m [K(p^{[i]}) + E(q^{[i]})]
+ \beta_n [K(p^{[j]}) + E(q^{[j]})] \Big\}. $$

Substituting the momentum rescaling $K(p^{[i]'}) = (T_n/T_m) K(p^{[i]})$ and
$K(p^{[j]'}) = (T_m/T_n) K(p^{[j]})$, and using $\beta_m T_m = \beta_n T_n
= 1/k_B$, the kinetic-energy contributions cancel exactly:
$-\beta_m (T_m/T_n) K(p^{[j]}) + \beta_n K(p^{[j]}) = 0$ and
$-\beta_n (T_n/T_m) K(p^{[i]}) + \beta_m K(p^{[i]}) = 0$. What remains is
$\exp(-\Delta)$ with $\Delta = (\beta_n - \beta_m)(E(q^{[i]}) - E(q^{[j]}))$
(eq:15). The momentum rescaling is thus exactly what makes the kinetic terms
drop out, leaving the same potential-energy criterion as the MC version.

### Simulation protocol

Without loss of generality assume $\beta_1 < \beta_2 < \dots < \beta_M$. A REM
simulation alternates two steps:

1. Each replica in the canonical ensemble at its fixed temperature is simulated
   *simultaneously* and *independently* for a certain number of MC or MD steps.
2. A pair of replicas at neighboring temperatures, say $x_m^{[i]}$ and
   $x_{m+1}^{[j]}$, are exchanged with the probability $w(x_m^{[i]}|x_{m+1}^{[j]})$
   (eq:17).

Here Step (1) uses the MD algorithm. In Step (2) only pairs of replicas at
neighboring temperatures are exchanged, because the acceptance ratio decreases
exponentially with the difference of the two $\beta$s. Whenever a replica exchange
is accepted, the permutation functions are updated.

The major advantage of REM over multicanonical algorithm and simulated tempering
is that the weight factor is a priori known (eq:7), while in the latter algorithms
the determination of the weight factors can be tedious and time-consuming. For
optimal performance one still has to choose an appropriate temperature distribution.

### Expectation values and reweighting

The canonical expectation value of a physical quantity $A$ at temperature $T_m$ is
computed by the usual arithmetic mean (eq:18, eq:19). For expectation values at
any intermediate temperature we use the multiple-histogram reweighting techniques
(an extension of which is also referred to as WHAM), eq:20-eq:22. Suppose we have
made $R$ independent simulation runs at $R$ different temperatures. Let $N_m(E)$
and $n_m$ be the energy histogram and the total number of samples obtained in the
$m$th run. In REM, $n_m = N_{\text{sim}}$. Here $g_m = 1 + 2\tau_m$ and $\tau_m$
is the integrated autocorrelation time at temperature $T_m$. $P(E;\beta)$ and
$f_m$ are solved self-consistently by iteration.

## 3. Results and discussion

The algorithm was tested for the penta-peptide Met-enkephalin (sequence
Tyr-Gly-Gly-Phe-Met) in gas phase. The N and C termini were blocked with acetyl
and N-methyl groups. Force-field parameters were taken from the all-atom version
of AMBER; the dielectric constant was set to 1. The temperature during the MD
simulations was controlled by the constraint method. The unit timestep was 0.5 fs.
An MD simulation of $2 \times 10^6$ timesteps (1.0 ns) was made for each replica,
starting from an extended conformation. Before taking data, regular canonical MD
simulations were run for 100 ps at each temperature, followed by a
replica-exchange simulation of 100 ps for thermalization.

Eight temperatures were used ($M = 8$): 700, 585, 489, 409, 342, 286, 239, and
200 K, distributed exponentially following the annealing schedule of simulated
annealing simulations. This choice already gave an optimal temperature
distribution. Replica exchange was tried every 10 fs, and data were stored just
before the replica exchange for later analyses. Thus $N_{\text{sim}} = 10^5$ for
each replica.

A replica-exchange simulation is particularly suitable for parallel computers.
Because one can minimize the amount of information exchanged among nodes, it is
best to assign each replica to a node (exchanging pairs of temperature values
among nodes is much faster than exchanging coordinates and momenta). This means we
keep track of the permutation function $m(i;t) = f^{-1}(i;t)$ during the
simulation. After every 10 fs of parallel MD, four pairs of replicas corresponding
to neighboring temperatures were exchanged, and the pairing was alternated between
the two possible choices.

For expectation values at various temperatures, Eqs. 20-22 were used with
$R = M = 8$, taking into account all replica runs. For biomolecular systems the
integrated autocorrelation times are approximately equal, so $g = \text{const}$
was set in eq:21 for simplicity. Also $n_m = N_{\text{sim}} = 10^5$.

Three points were checked to verify proper performance: (a) whether temperatures
were optimally distributed, (b) whether the number of replicas was sufficient, and
(c) whether the highest temperature was sufficiently high. Points (a) and (b) are
checked by examining the acceptance ratios of replica exchange for adjacent
temperature pairs. Optimal distribution implies all acceptance ratios are equal
(a free random walk in temperature space); a sufficient number of replicas
requires acceptance ratios not too small (greater than about 0.1). The measured
ratios were uniform (about 15%) and larger than 10%, so criteria (a) and (b) were
met.

The distributions of the pair of dihedral angles $(\phi,\psi)$ of Gly-2 at
$T = 200$ K show that the regular canonical simulation is localized with one
dominant peak, whereas the replica-exchange simulation has several peaks. The
replica-exchange run samples a much broader configuration space at low
temperatures. The average potential energy at 200 K of the conformation
corresponding to the highest canonical peak is about 2 kcal/mol higher than for
the replica-exchange simulation (-141 versus ca. -143 kcal/mol). At $T = 700$ K
the results are similar, implying regular canonical simulation gives accurate
thermodynamic quantities at high temperatures.

A single replica-exchange simulation can give any thermodynamic quantity as a
function of temperature via multiple-histogram reweighting (Eqs. 20-22). About
10-100 iterations were necessary for convergence. The average total potential
energy as a function of temperature shows that canonical simulations at low
temperatures got trapped in energy local minima, starting near 300 K and below,
which is an experimentally relevant temperature.

## 4. Conclusions

We have presented a formulation of an MD algorithm for the replica-exchange method.
In this method the weight factor is essentially known and there is no complication
in its determination, while in other generalized-ensemble algorithms the
determination of the non-Boltzmann weight factor can be tedious and time-consuming.
The effectiveness of the method was tested with Met-enkephalin. From a single
simulation run one can obtain various thermodynamic quantities as a function of
temperature for a wide range. The method is particularly useful for studying the
protein-folding problem where information about a wide conformational space (from a
random-coil state to the native folded state) is required.
