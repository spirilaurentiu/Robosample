# Replica exchange and expanded ensemble simulations as Gibbs sampling: Simple improvements for enhanced mixing

John D. Chodera, Michael R. Shirts. J. Chem. Phys. 135, 194110 (2011). doi:10.1063/1.3660669

## Abstract

Replica exchange and expanded ensemble algorithms are shown to be special cases of Gibbs sampling within a Markov chain Monte Carlo framework. Gibbs sampling alternately updates different random variables from their conditional distributions. While updating conformational degrees of freedom by Metropolis MC or MD generates correlated samples, judicious updating of the thermodynamic state indices (temperature, alchemical coupling, etc.) can substantially increase mixing while still sampling the correct distributions. State-update methods in common use can mix suboptimally; the paper presents simple, inexpensive alternatives that increase mixing of the overall Markov chain, reducing simulation time to a target precision. Demonstrated on an alchemical expanded-ensemble simulation, parallel tempering, and 2D replica-exchange umbrella sampling.

## Key idea (routing summary)

Both expanded-ensemble (single walker over states) and replica-exchange (K coupled walkers, a permutation of states) simulations sample a joint distribution over configuration and thermodynamic-state variables. They are Gibbs samplers that alternate: (1) update configuration(s) given fixed state(s) via MD/MC, and (2) update the discrete state index / permutation given fixed configuration(s). The novel contribution is better step (2): replacing the traditional single neighbor-exchange move with **independence sampling** (draw the new state directly from the exact conditional), **Metropolized independence sampling**, or **restricted-range sampling**; and for replica exchange, running many random pair-swaps per iteration to approach independence sampling of the permutation. These state updates are essentially free (no extra energy evaluations for temperature/pressure/pH exchange, and reuse MBAR energies otherwise) yet accelerate both state-space and configuration-space mixing.

## II. Theory

### A. Thermodynamic states and ensembles

A thermodynamic state is parameterized by a vector $\lambda\equiv\{\beta,H,p,\mu,\dots\}$. The reduced potential is defined in Eq. 1. Configuration $x\in\Omega$; unnormalized density $q(x)>0$; partition function Eq. 2; normalized density Eq. 3. Boltzmann statistics give $q(x)=e^{-u(x)}$ (Eq. 4). A set of $K$ states $\lambda_k$, $k=1,\dots,K$, each has reduced potential $u_k(x)$, densities $q_k,\pi_k$ (also written $\pi(k,x)$), and partition function $Z_k$. Non-Boltzmann statistics (multicanonical, Tsallis) are allowed via alternative $q_k$. All states must specify the same set of thermodynamic parameters (values may differ) so every configuration has finite nonzero density in all $K$ states.

### B. Gibbs sampling

To sample a joint $\pi(x,y)$ that is hard to draw from directly, alternately sample from the conditionals $\pi(x|y)$ and $\pi(y|x)$ (Eq. gibbs-update). Samples may be correlated; the scheme still targets the correct joint. Choice of which variable to update can be stochastic (obeys detailed balance) or deterministic (obeys the weaker balance condition); both preserve $\pi(x,y)$ but have different correlation structure. Expanded ensemble and replica exchange are Gibbs sampling on $\pi(x,k)$. The state variable $k$ is discrete here but could be continuous.

### C. Expanded ensembles

Single walker samples pairs $(x,k)$ from $\pi(x,k)\propto\exp[-u_k(x)+g_k]$ (Eq. 5), $g_k$ a state log-weight (often $g_k=-\ln Z_k$ for equal per-state probability, determined iteratively e.g. Wang-Landau). Gibbs conditionals: $\pi(x|k)$ (Eq. 6) sampled by MD/MC; $\pi(k|x)$ (Eq. 7) has a simple closed form enabling many state-update choices.

### D. Replica exchange ensembles

$K$ simulations, one per state; sampler state $(X,S)$ with $X=\{x_1,\dots,x_K\}$ and $S$ a permutation of state indices. Joint $\pi(X,S)$ (Eq. 8), conditionals $\pi(X|S)$ (Eq. 9, per-replica MD/MC) and $\pi(S|X)$ (Eq. 10). The permutation-conditional normalizer sums over all $K!$ permutations (the matrix permanent, #P-complete), so exact independent sampling of $S$ is impractical; effective decorrelation via a short MCMC swap loop is used instead.

## III. Algorithms

### A. Expanded ensemble state updates

Target conditional $\pi(k|x)$ (Eq. 7). Any proposal/acceptance sampling it in the long run is valid. One may choose stochastically to update $k$ or $x$ (detailed balance) or alternate $N_k$ and $N_x$ steps (weaker balance, still correct). History-dependent proposals (adaptive weights) break equilibrium unless handled carefully.

1. **Neighbor exchange** (Marinari-Parisi): propose $i\pm1$ each with prob 1/2 (Eq. 11), accept by Eq. 12. A torus variant identifies $i+nK\equiv i$. A boundary-corrected proposal (Eq. 13) with acceptance Eq. 14 avoids wasted out-of-range proposals.

2. **Independence sampling:** propose $j\sim\pi(k|x)$ (Eq. 15) and always accept. Implement by uniform $r\in[0,1)$ and CDF search. Generates uncorrelated state indices.

3. **Metropolized independence sampling:** propose from Eq. 16 (never propose current state), accept by Eq. 17. Provably faster mixing in $\pi(k|x)$ than plain independence sampling (Peskun-type argument) because it always tries to move away from the current state.

4. **Restricted range sampling:** for each state $i$ define a symmetric proposal set $\mathcal{S}_i$ (with $i\in\mathcal{S}_j\iff j\in\mathcal{S}_i$); propose within $\mathcal{S}_i$ by Eq. 18, accept by Eq. 19. Useful when evaluating $q_k(x)$ for all $k$ is costly. Reduces to independence sampling when $\mathcal{S}_i=\{1,\dots,K\}$; equals Metropolized independence when $\mathcal{S}_i$ excludes $i$. Care: some set schemes preserve detailed balance but are non-ergodic (e.g. odd/even partitions for $K=6$).

For $\mathcal{S}_i=\{i-n,\dots,i+n\}$, $n\ll K$, one only computes reduced potentials for states in $\{\min(1,i-2n),\dots,\max(K,i+2n)\}$ to evaluate both sums in the acceptance ratio.

### Derivation (not implemented): detailed balance of restricted range sampling

With $Z(\mathcal{S}_i)=\sum_{k\in\mathcal{S}_i}e^{g_k-u_k(x)}$ and $\mathcal{S}_{\text{all}}=\{1,\dots,K\}$, the probability of starting in $i\in\mathcal{S}_j$ and transitioning to $j\in\mathcal{S}_i$ ($j\ne i$) is

$$\pi(i|x)\,\alpha(j|x,i)\,P_{\text{accept}}(j|x,i) = \left[\frac{e^{g_i-u_i(x)}}{Z(\mathcal{S}_{\text{all}})}\right]\left[\frac{e^{g_j-u_j(x)}}{Z(\mathcal{S}_i)}\right]\left[\min\left(1,\frac{Z(\mathcal{S}_i)}{Z(\mathcal{S}_j)}\right)\right]$$

$$= \left[\frac{e^{g_j-u_j(x)}e^{g_i-u_i(x)}}{Z(\mathcal{S}_{\text{all}})}\right]\left[\min\left(Z^{-1}(\mathcal{S}_i),Z^{-1}(\mathcal{S}_j)\right)\right]$$

$$= \left[\frac{e^{g_j-u_j(x)}}{Z(\mathcal{S}_{\text{all}})}\right]\left[\frac{e^{g_i-u_i(x)}}{Z(\mathcal{S}_j)}\right]\left[\min\left(1,\frac{Z(\mathcal{S}_j)}{Z(\mathcal{S}_i)}\right)\right] = \pi(j|x)\,\alpha(i|x,j)\,P_{\text{accept}}(i|x,j),$$

which is exactly detailed balance, so the scheme samples $\pi(k|x)$.

5. **Other schemes:** any move sampling $\pi(k|x)$ is valid; compositions allowed (e.g. repeating neighbor exchange several times per iteration).

### B. Replica exchange state updates

1. **Neighbor exchange:** attempt exchanging either the pairs $\{(1,2),(3,4),\dots\}$ or $\{(2,3),(4,5),\dots\}$ with equal probability; each pair swap accepted by Eq. 24.

2. **Independence sampling of the permutation:** exact independent draws from $\pi(S|X)$ require the $K!$ permanent (intractable). Instead run a short MCMC swap loop: repeatedly pick a random distinct pair $(i,j)$, apply the Eq. 24 swap criterion, updating labels after each accepted swap. With precomputed $\mathbf{U}=(u_{ij})$, $u_{ij}=u_i(x_j)$, these updates are cheap; $K^3$ to $K^5$ swaps per iteration effectively decorrelate the permutation. If all $u_{ij}$ equal, $\sim K\ln K$ swaps suffice (Aldous-Diaconis). For parallel tempering $\mathbf{U}$ is trivial from the $K$ potential energies; for alchemical/Hamiltonian exchange it needs $K^2$ energy evaluations, but these are exactly the energies MBAR needs, so no extra work. Multiple passes of neighbor exchange also approach independence sampling.

3. **Other schemes:** any sampler preserving $\pi(S|X)$ is valid; e.g. local exchanges within a compute node.

### C. Metrics of efficiency

Three surrogate "mixing times" from the projected state-index trajectory $\mathbf{s}=\{s_0,s_1,\dots\}$:

1. **$\tau_2$** from the second eigenvalue of the empirical (symmetrized, row-stochastic) transition matrix $\mathbf{T}$ (Eq. 25); relaxation time Eq. 26. $\mu_2=1$ signals a decomposable (non-mixing) chain. $\tau_2$ is a lower bound on the observed state and configuration correlation times.
2. **$\tau_{ac}$** the integrated autocorrelation time of the state index (statistical inefficiency $=2\tau_{ac}+1$).
3. **$\tau_{end}$** average end-to-end (round-trip) transit time of the state index between $k=1$ and $k=K$.

## IV. Model illustration

1D double-well $U(x)=10(x-1)^2(x+1)^2$ (Eq. 27), simulated tempering with $K$ geometric temperatures $\beta_k=10^{-(k-1)/(K-1)}$ over $k_BT\in[1,10]$ (Eq. 29), exact weights $g_k$ by Eq. 28. Each iteration = 1 state update + 100 Metropolis MC steps (Gaussian, std 0.1). Independence sampling greatly reduces state-index correlation and, because $k$ and $x$ are coupled, also substantially reduces the configuration correlation $\tau_x$ (see checks.md). Increasing $K$ raises correlation times for neighbor exchange but keeps them small for independence sampling.

## V. Applications

Three systems (numbers in checks.md):

- **A. Expanded-ensemble alchemical decoupling** of UA methane ($K=6$) and a large LJ sphere ($K=18$) in TIP3P water (GROMACS, softcore LJ Eq. 30). Independence and Metropolized independence sampling give statistically significant speedups (up to ~3.5x) over neighbor exchange, most pronounced for infrequent state updates. Repeating any correct move many times converges to independence sampling.
- **B. Parallel tempering** of blocked alanine dipeptide, OBC GBSA implicit solvent, OpenMM, $K^3$ random pair swaps. State-space mixing accelerated ~1-2 orders of magnitude; torsional correlation ~2-10x.
- **C. 2D replica-exchange umbrella sampling** of alanine dipeptide, $K=101$ on a $10\times10$ toroidal grid, periodic von Mises bias (Eq. 31, $\kappa=(2\pi/30)^{-2}\beta^{-1}$, neighbors $3\sigma$ apart). State relaxation reduced 2-6x, structural correlation 4-5x.

## VI. Discussion / recommendations

- For temperature/pressure/pH exchange (potential independent of state), independence sampling is effectively free; adopt it.
- If state updates are costly (simulated scaling, Hamiltonian exchange), either reuse MBAR energies, or use restricted-range updates for improved mixing at minimal extra evaluations.
- More frequent state updates are always better than less frequent; frequent neighbor exchange can beat infrequent independence sampling.
- If parallel state updates are expensive/complex, simply repeating the existing swap scheme several times per iteration approaches independence sampling with little code change.
- The empirical transition matrix and its dominant eigenvalues are useful diagnostics for equilibration, convergence, and poor state overlap.
- These methods are not a cure-all: first-order (phase-transition-like) systems can still mix exponentially slowly; optimal state selection remains an open problem.
- **Continuous tempering limit:** as $K\to\infty$ between fixed limits, the discrete index becomes a continuous parameter $\lambda$ sampling $\pi(x,\lambda)\propto\exp[-\lambda h(x)+g(\lambda)]$ (Eq. continuous-tempering), with continuous log-weight $g(\lambda)$.

## Note on MD sampling of $\pi(x|k)$

Caution: MD used to sample $\pi(x|k)\propto e^{-u_k(x)}$ deviates from the target in a timestep-dependent way unless wrapped in a Monte Carlo acceptance step (hybrid MC, metropolized Langevin, generalized hybrid MC). Uncorrected deviation can bias the joint $(x,k)$ chain.
