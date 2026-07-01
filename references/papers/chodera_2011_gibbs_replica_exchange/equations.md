# Equations - Chodera & Shirts 2011, Gibbs sampling for replica exchange / expanded ensemble

Reduced units throughout: energies are in units of $k_B T$ (via the reduced potential $u_k$).
State index $k \in \{1,\dots,K\}$; configuration $x$.

<!-- eq:1 -->
$$ u(x) = \beta \left[ H(x) + p V(x) + \sum_{i} \mu_i n_i(x) + \cdots \right] $$
- **what:** Reduced (dimensionless) potential of a physical system for a given thermodynamic state.
- **symbols:** $u(x)$ - reduced potential (dimensionless); $\beta = 1/(k_B T)$ - inverse temperature; $H(x)$ - Hamiltonian/energy (energy units); $p$ - pressure; $V(x)$ - volume; $\mu_i$ - chemical potential of component $i$; $n_i(x)$ - number of molecules of component $i$; $x$ - configuration.

<!-- eq:2 -->
$$ Z \equiv \int_{\Omega} dx \, q(x) $$
- **what:** Partition function (normalization constant) of an unnormalized density over configuration space $\Omega$.
- **symbols:** $Z$ - partition function; $\Omega$ - accessible configuration space; $q(x)$ - unnormalized probability density.

<!-- eq:3 -->
$$ \pi(x) = Z^{-1} q(x) $$
- **what:** Normalized probability density.
- **symbols:** $\pi(x)$ - normalized density; $Z$ - partition function; $q(x)$ - unnormalized density.

<!-- eq:4 -->
$$ q(x) \equiv e^{-u(x)} $$
- **what:** Boltzmann unnormalized density in reduced-potential form.
- **symbols:** $q(x)$ - unnormalized density; $u(x)$ - reduced potential (Eq. 1).

<!-- eq:gibbs-update -->
$$ x^{(n+1)}\,|\,y^{(n)} \sim \pi(x\,|\,y^{(n)}), \qquad y^{(n+1)}\,|\,x^{(n+1)} \sim \pi(y\,|\,x^{(n+1)}) $$
- **what:** Gibbs sampler alternating conditional updates of two blocks of variables; generates samples of the joint $\pi(x,y)$.
- **symbols:** $x,y$ - two random-variable blocks; $\pi(\cdot|\cdot)$ - conditional distribution; superscript $(n)$ - iteration index; $\sim$ - "sampled from".

<!-- eq:5 -->
$$ \pi(x,k) \propto \exp[-u_k(x) + g_k] $$
- **what:** Expanded-ensemble joint distribution of configuration $x$ and state index $k$ (single walker).
- **symbols:** $\pi(x,k)$ - joint density; $u_k(x)$ - reduced potential of state $k$; $g_k$ - state-dependent log-weight (bias) factor; choosing $g_k=-\ln Z_k$ gives equal probability per state.

<!-- eq:6 -->
$$ \pi(x\,|\,k) = \frac{q_k(x)}{\int_{\Omega} dx \, q_k(x)} = \frac{e^{-u_k(x)}}{\int_{\Omega} dx \, e^{-u_k(x)}} $$
- **what:** Conditional distribution of configuration given state index (sampled by MD/MC).
- **symbols:** $q_k(x)=e^{-u_k(x)}$ - unnormalized density of state $k$; $u_k$ - reduced potential of state $k$.

<!-- eq:7 -->
$$ \pi(k\,|\,x) = \frac{e^{g_k} q_k(x)}{\sum_{k'=1}^{K} e^{g_{k'}} q_{k'}(x)} = \frac{e^{g_k - u_k(x)}}{\sum_{k'=1}^{K} e^{g_{k'} - u_{k'}(x)}} $$
- **what:** Conditional distribution of state index given configuration (the discrete state-update target for expanded ensemble). Central to all expanded-ensemble state moves.
- **symbols:** $K$ - number of thermodynamic states; $g_k$ - log-weight; $u_k(x)$ - reduced potential; $k'$ - summation index over states.

<!-- eq:8 -->
$$ \pi(X,S) \propto \prod_{i=1}^{K} q_{s_i}(x_i) \propto \exp\left[-\sum_{i=1}^{K} u_{s_i}(x_i)\right] $$
- **what:** Replica-exchange joint distribution over configuration vector $X$ and state-index permutation $S$.
- **symbols:** $X=\{x_1,\dots,x_K\}$ - vector of replica configurations; $S=\{s_1,\dots,s_K\}$ - permutation of state indices $\{1,\dots,K\}$; $s_i$ - state index assigned to replica $i$; $u_{s_i}(x_i)$ - reduced potential of state $s_i$ evaluated at config $x_i$.

<!-- eq:9 -->
$$ \pi(X\,|\,S) = \prod_{i=1}^{K} \left[ \frac{e^{-u_{s_i}(x_i)}}{\int_{\Omega} dx \, e^{-u_{s_i}(x_i)}} \right] $$
- **what:** Conditional over configurations given fixed permutation (independent per-replica MD/MC updates).
- **symbols:** as Eq. 8.

<!-- eq:10 -->
$$ \pi(S\,|\,X) = \frac{\exp\left[-\sum_{i=1}^{K} u_{s_i}(x_i)\right]}{\sum_{S' \in \mathcal{S}_K} \exp\left[-\sum_{i=1}^{K} u_{s_i'}(x_i)\right]} $$
- **what:** Conditional over permutations given fixed configurations; denominator sums over all $K!$ permutations (the matrix permanent - #P-complete).
- **symbols:** $\mathcal{S}_K$ - symmetric group of all permutations of $K$ indices; $S'=\{s_1',\dots,s_K'\}$ - a permutation.

<!-- eq:11 -->
$$ \alpha(j\,|\,x,i) = \begin{cases} \tfrac{1}{2} & \text{if } j = i - 1\\ \tfrac{1}{2} & \text{if } j = i + 1\\ 0 & \text{else} \end{cases} $$
- **what:** Neighbor-exchange proposal probability for the expanded-ensemble state index (Marinari-Parisi).
- **symbols:** $\alpha(j|x,i)$ - probability of proposing new state $j$ from current state $i$ (config $x$ fixed).

<!-- eq:12 -->
$$ P_{\text{accept}}(j\,|\,x,i) = \begin{cases} 0 & \text{if } j \notin \{1,\dots,K\} \\ \min\left\{1,\; \dfrac{e^{g_j - u_j(x)}}{e^{g_i - u_i(x)}}\right\} & \text{else} \end{cases} $$
- **what:** Metropolis acceptance for neighbor-exchange state move (Eq. 11 proposal).
- **symbols:** $g_j,g_i$ - log-weights; $u_j(x),u_i(x)$ - reduced potentials at states $j,i$; out-of-range proposals rejected.

<!-- eq:13 -->
$$ \alpha(j\,|\,x,i) = \begin{cases} \tfrac{1}{2} & \text{if } i \in \{2,\dots,K-1\},\; |j-i| = 1\\ 1 & \text{if } i = 1,\; j = i+1 \le K\\ 1 & \text{if } i = K,\; j = i-1 \ge 1\\ 0 & \text{else} \end{cases} $$
- **what:** Boundary-corrected neighbor proposal that avoids proposing out-of-range states.
- **symbols:** as Eq. 11. <!-- CHECK: source printed "k" in first case; interpreted as current index i per context -->

<!-- eq:14 -->
$$ P_{\text{accept}}(j\,|\,x,i) = \min\left\{1,\; \frac{1}{2}\,\frac{e^{g_j - u_j(x)}}{e^{g_i - u_i(x)}}\right\} $$
- **what:** Acceptance for the two boundary moves ($i=1$ or $i=K$) of Eq. 13, including the proposal-ratio factor $\tfrac12$.
- **symbols:** as Eq. 12.

<!-- eq:15 -->
$$ \alpha(j\,|\,x,i) = \pi(j\,|\,x) $$
- **what:** Independence sampling: propose new state directly from the exact conditional $\pi(k|x)$ (Eq. 7) and always accept. Implement by drawing uniform $r\in[0,1)$ and picking smallest $k$ with $r < \sum_{i=1}^{k}\pi(i|x)$.
- **symbols:** $\pi(j|x)$ - conditional state distribution (Eq. 7). <!-- CHECK: source wrote alpha=pi(i|x); should be pi(j|x), the proposed state -->

<!-- eq:16 -->
$$ \alpha(j\,|\,x,i) = \begin{cases} \dfrac{\pi(j\,|\,x)}{1-\pi(i\,|\,x)} & j \neq i\\[2mm] 0 & j = i \end{cases} $$
- **what:** Metropolized independence sampling proposal: propose only states other than current, weighted by conditional excluding current state.
- **symbols:** $\pi(\cdot|x)$ - conditional state distribution (Eq. 7); denominator $1-\pi(i|x)$ renormalizes over $j\ne i$. <!-- CHECK: source printed denominator 1-pi(j|x,i); corrected to 1-pi(i|x) so alpha sums to 1 over j != i -->

<!-- eq:17 -->
$$ P_{\text{accept}}(j\,|\,x,i) = \min\left\{1,\; \frac{1 - \pi(i\,|\,x)}{1 - \pi(j\,|\,x)}\right\} $$
- **what:** Acceptance for Metropolized independence sampling; provably faster mixing in $\pi(k|x)$ than plain independence sampling.
- **symbols:** $\pi(i|x),\pi(j|x)$ - conditional probabilities of current and proposed states.

<!-- eq:18 -->
$$ \alpha(j\,|\,x,i) = \begin{cases} \dfrac{e^{g_j - u_j(x)}}{\sum\limits_{k \in \mathcal{S}_i} e^{g_k - u_k(x)}} & j \in \mathcal{S}_i \\ 0 & j \notin \mathcal{S}_i \end{cases} $$
- **what:** Restricted-range sampling proposal: sample proposed state from the conditional restricted to a neighbor set $\mathcal{S}_i$.
- **symbols:** $\mathcal{S}_i$ - proposal set for state $i$ (require $i\in\mathcal{S}_j \iff j\in\mathcal{S}_i$); e.g. $\mathcal{S}_i=\{i-n,\dots,i+n\}$ with $n\ll K$.

<!-- eq:19 -->
$$ P_{\text{accept}}(j\,|\,x,i) = \min\left(1,\; \frac{\sum_{k \in \mathcal{S}_i} e^{g_k - u_k(x)}}{\sum_{k' \in \mathcal{S}_j} e^{g_{k'} - u_{k'}(x)}}\right) $$
- **what:** Acceptance for restricted-range sampling; ratio of restricted partition functions $Z(\mathcal{S}_i)/Z(\mathcal{S}_j)$.
- **symbols:** $Z(\mathcal{S}_i) = \sum_{k\in\mathcal{S}_i} e^{g_k - u_k(x)}$ - restricted partition function. <!-- CHECK: source denominator exponent printed g_{k'}-g_{k'}(x); corrected to g_{k'}-u_{k'}(x) per Z(S_j) definition -->

<!-- eq:23 -->
$$ Z(\mathcal{S}_i) = \sum_{k \in \mathcal{S}_i} e^{g_k - u_k(x)}, \qquad \mathcal{S}_{\text{all}} = \{1, \dots, K\} $$
- **what:** Definition of the restricted partition function used in Eqs. 18-19 (detailed-balance proof, Eqs. 20-23).
- **symbols:** $\mathcal{S}_i$ - restricted proposal set; $\mathcal{S}_{\text{all}}$ - full state set.

<!-- eq:24 -->
$$ P_{\text{accept}}(x_i, i, x_j, j) = \min\left\{1,\; \frac{e^{-[u_i(x_j) + u_j(x_i)]}}{e^{-[u_i(x_i) + u_j(x_j)]}}\right\} $$
- **what:** Replica-exchange pair-swap Metropolis acceptance for swapping states $i,j$ between configs $x_i,x_j$. Used both in neighbor exchange and in the MCMC-permutation ("independence") swap loop.
- **symbols:** $u_i(x_j)$ - reduced potential of state $i$ evaluated at config $x_j$; the four terms are the cross vs. current energy assignments.

<!-- eq:25 -->
$$ T_{ij} \approx \frac{N_{ij} + N_{ji}}{\sum_{k=1}^{K} [N_{ik} + N_{ki}]} $$
- **what:** Empirical (symmetrized, row-stochastic) state-transition matrix estimator used to diagnose mixing.
- **symbols:** $T_{ij}$ - estimated probability of being in state $j$ one update after state $i$; $N_{ij}$ - observed count of transitions $i\to j$ across the pooled replica trajectories; symmetrization $N_{ij}+N_{ji}$ enforces real eigenvalues.

<!-- eq:26 -->
$$ \tau_2 = \frac{\tau}{1 - \mu_2} $$
- **what:** Relaxation time from second-largest eigenvalue of $\mathbf{T}$; lower bound on state-index correlation time.
- **symbols:** $\tau_2$ - relaxation time; $\tau$ - effective time between exchange attempts; $\mu_2$ - second-largest eigenvalue of $\mathbf{T}$ (with $1=\mu_1\ge\mu_2\ge\cdots\ge\mu_K$).

<!-- eq:27 -->
$$ U(x) = 10\,(x-1)^{2}(x+1)^{2} $$
- **what:** 1D double-well model potential used in the illustration.
- **symbols:** $U(x)$ - potential energy; $x$ - scalar coordinate; minima at $x=\pm1$, barrier at $x=0$.

<!-- eq:28 -->
$$ g_k = -\ln \int_{-\infty}^{+\infty} dx \, e^{-\beta_k U(x)} $$
- **what:** Exact log-weights giving each of the $K$ states equal probability (numerically computed for the 1D model).
- **symbols:** $g_k$ - log-weight of state $k$; $\beta_k$ - inverse temperature of state $k$; $U(x)$ - Eq. 27 potential.

<!-- eq:29 -->
$$ \beta_k = 10^{-(k-1)/(K-1)} \quad \text{for } k = 1, \dots, K $$
- **what:** Geometric spacing of the $K$ inverse temperatures over $k_B T \in [1,10]$ in the 1D model.
- **symbols:** $\beta_k$ - inverse temperature of state $k$; $K$ - number of temperature states ($k_BT$ ranges 1 at $k=1$ to 10 at $k=K$).

<!-- eq:30 -->
$$ U_{ij}(r;\lambda) = 4\epsilon_{ij}\,\lambda\, f(r;\lambda)\,[1 - f(r;\lambda)], \qquad f(r;\lambda) \equiv [\alpha(1-\lambda) + (r/\sigma_{ij})^6]^{-1} $$
- **what:** Softcore Lennard-Jones potential used for alchemical decoupling in the expanded-ensemble applications.
- **symbols:** $U_{ij}$ - pair energy; $r$ - pair distance; $\lambda$ - alchemical coupling ($\lambda=1$ full interaction, $\lambda=0$ decoupled); $\epsilon_{ij}$ - LJ well depth; $\sigma_{ij}$ - LJ diameter; $\alpha$ - softcore parameter. <!-- CHECK: softcore alpha value not given numerically in text -->

<!-- eq:31 -->
$$ U_k'(x) \equiv -\kappa \left[ \cos\left(\phi - \phi_k^0\right) + \cos\left(\psi - \psi_k^0\right) \right] $$
- **what:** Periodic (von Mises) umbrella bias potential restraining torsions near reference values in 2D REX umbrella sampling.
- **symbols:** $U_k'$ - bias potential of state $k$ (energy units); $\kappa$ - force constant (energy); $\phi,\psi$ - torsion angles; $(\phi_k^0,\psi_k^0)$ - reference torsions of state $k$. Approx. Gaussian width $\sigma \equiv (\beta\kappa)^{-1/2}$. <!-- CHECK: text wrote sigma=(beta kappa)^{1/2}; corrected to (beta kappa)^{-1/2} so larger kappa gives narrower well -->

<!-- eq:continuous-tempering -->
$$ \pi(x,\lambda) \propto \exp\left[-\lambda h(x) + g(\lambda)\right] $$
- **what:** Continuous-tempering limit of independence sampling: state index $k$ becomes continuous parameter $\lambda$ multiplying a conjugate variable $h(x)$.
- **symbols:** $\lambda$ - continuous thermodynamic parameter; $h(x)$ - conjugate configuration-dependent variable; $g(\lambda)$ - continuous log-weight function.
