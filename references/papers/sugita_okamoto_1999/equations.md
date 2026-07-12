# Equations - Sugita & Okamoto 1999 (Replica-Exchange MD)

<!-- eq:1 -->
$$ H(q,p) = K(p) + E(q) $$
- **what:** total Hamiltonian is kinetic plus potential energy.
- **symbols:** $q = \{q_1,\dots,q_N\}$ - atomic coordinates; $p = \{p_1,\dots,p_N\}$ - atomic momenta; $K$ - kinetic energy; $E$ - potential energy.

<!-- eq:2 -->
$$ K(p) = \sum_{k=1}^{N} \frac{p_k^2}{2 m_k} $$
- **what:** kinetic energy of $N$ atoms.
- **symbols:** $p_k$ - momentum of atom $k$; $m_k$ - mass of atom $k$; $N$ - number of atoms.

<!-- eq:3 -->
$$ W_B(x;T) = e^{-\beta H(q,p)} $$
- **what:** canonical Boltzmann weight of a phase-space state $x=(q,p)$ at temperature $T$.
- **symbols:** $x \equiv (q,p)$ - phase-space state; $\beta = 1/k_B T$ - inverse temperature; $k_B$ - Boltzmann constant.

<!-- eq:4 -->
$$ \langle K(p) \rangle_T = \left\langle \sum_{k=1}^N \frac{p_k^2}{2 m_k} \right\rangle_T = \frac{3}{2} N k_B T $$
- **what:** equipartition; average kinetic energy at temperature $T$. This is the condition the momentum rescaling (eq:12) must preserve after a swap.
- **symbols:** $\langle\cdot\rangle_T$ - canonical average at temperature $T$; $N$ - number of atoms.

<!-- eq:5 -->
$$ i = i(m) \equiv f(m), \qquad m = m(i) \equiv f^{-1}(i) $$
- **what:** one-to-one map between replica label $i$ and temperature label $m$. $f$ is a permutation; each temperature holds exactly one replica.
- **symbols:** $i$ - replica label ($1,\dots,M$); $m$ - temperature label ($1,\dots,M$); $f$ - permutation function; $f^{-1}$ - its inverse; $M$ - number of replicas/temperatures.

<!-- eq:6 -->
$$ x_m^{[i]} = (q^{[i]}, p^{[i]})_m $$
- **what:** microstate of replica $i$ currently at temperature $T_m$ (its coordinates and momenta).
- **symbols:** superscript $[i]$ - replica index; subscript $m$ - temperature index; $q^{[i]}, p^{[i]}$ - coordinates and momenta of replica $i$.

<!-- eq:7 -->
$$ W_{\text{REM}}(X) = \exp\left\{-\sum_{i=1}^{M} \beta_{m(i)} H(q^{[i]}, p^{[i]})\right\} = \exp\left\{-\sum_{m=1}^{M} \beta_m H(q^{[i(m)]}, p^{[i(m)]})\right\} $$
- **what:** weight of a full generalized-ensemble state $X$; product of Boltzmann factors over non-interacting replicas. A priori known (main advantage of REM).
- **symbols:** $X$ - full state of all $M$ replicas; $\beta_m = 1/k_B T_m$ - inverse temperature of temperature slot $m$; $i(m), m(i)$ - permutation functions (eq:5).

<!-- eq:8 -->
$$ X = (\dots, x_m^{[i]}, \dots, x_n^{[j]}, \dots) \to X' = (\dots, x_m^{[j]'}, \dots, x_n^{[i]'}, \dots) $$
- **what:** a replica-exchange move: replicas $i$ and $j$ (at temperatures $T_m$ and $T_n$) swap temperature slots.
- **symbols:** primed states carry rescaled momenta (eq:12); $m,n$ - the two temperature slots being exchanged.

<!-- eq:12 -->
$$ p^{[i]'} \equiv \sqrt{\frac{T_n}{T_m}}\, p^{[i]}, \qquad p^{[j]'} \equiv \sqrt{\frac{T_m}{T_n}}\, p^{[j]} $$
- **what:** momentum (velocity) rescaling on a swap. Replica $i$ moves from $T_m$ to $T_n$, so its momenta are scaled by $\sqrt{T_n/T_m}$; replica $j$ moves from $T_n$ to $T_m$. Preserves the equipartition condition (eq:4). This is the key MD-specific ingredient.
- **symbols:** $T_m, T_n$ - the two temperatures being exchanged; unprimed = before swap, primed = after swap.

<!-- eq:13 -->
$$ W_{\text{REM}}(X)\, w(X \to X') = W_{\text{REM}}(X')\, w(X' \to X) $$
- **what:** detailed balance condition imposed on the exchange transition probability $w$.
- **symbols:** $w(X\to X')$ - probability of accepting the swap $X\to X'$.

<!-- eq:14 -->
$$ \frac{w(X \to X')}{w(X' \to X)} = \exp(-\Delta) $$
- **what:** after substituting eq:7 and the momentum rescaling eq:12, the kinetic-energy terms cancel and the transition-probability ratio reduces to $\exp(-\Delta)$. (Full expansion in paper.md Derivation.)
- **symbols:** $\Delta$ - reduced swap cost (eq:15).

<!-- eq:15 -->
$$ \Delta \equiv (\beta_n - \beta_m)\,(E(q^{[i]}) - E(q^{[j]})) $$
- **what:** the swap acceptance argument. Depends only on potential energies of the two replicas and the two inverse temperatures; momenta drop out because of eq:12. Feed this into eq:17.
- **symbols:** $\beta_m = 1/k_B T_m$, $\beta_n = 1/k_B T_n$ - inverse temperatures of slots $m,n$; $E(q^{[i]}), E(q^{[j]})$ - potential energies of replicas $i,j$; $i=f(m)$, $j=f(n)$ before the exchange.

<!-- eq:17 -->
$$ w(X \to X') \equiv w(x_m^{[i]} \mid x_n^{[j]}) = \begin{cases} 1, & \Delta \le 0, \\ \exp(-\Delta), & \Delta > 0, \end{cases} $$
- **what:** Metropolis acceptance criterion for a replica exchange; identical to the MC-REM criterion. In practice only neighboring temperature pairs ($n=m+1$) are swapped.
- **symbols:** $w(x_m^{[i]}|x_n^{[j]})$ - probability of exchanging the pair at temperatures $T_m,T_n$; $\Delta$ from eq:15.

<!-- eq:18 -->
$$ \langle A \rangle_{T_m} = \frac{1}{N_{\text{sim}}} \sum_{t=1}^{N_{\text{sim}}} \sum_{i=1}^{M} A\!\left[ x_{f^{-1}(i;t)}^{[i]}(t) \right] \delta_{f^{-1}(i;t),\,m} $$
- **what:** canonical expectation of observable $A$ at temperature $T_m$ (replica view): average over measurements, selecting the replica currently at slot $m$.
- **symbols:** $N_{\text{sim}}$ - total measurements per replica; $f^{-1}(i;t)$ - temperature slot of replica $i$ at measurement $t$; $\delta_{k,l}$ - Kronecker delta; $A[\cdot]$ - observable evaluated on a microstate.

<!-- eq:19 -->
$$ \langle A \rangle_{T_m} = \frac{1}{N_{\text{sim}}} \sum_{t=1}^{N_{\text{sim}}} A\!\left( x_m^{[f(m;t)]}(t) \right) $$
- **what:** same expectation, temperature-exchange view: at each measurement, evaluate $A$ on whichever replica occupies temperature slot $m$.
- **symbols:** $f(m;t)$ - replica occupying temperature slot $m$ at measurement $t$.

<!-- eq:20 -->
$$ \langle A \rangle_T = \frac{\sum_E A(E)\, P(E;\beta)}{\sum_E P(E;\beta)} $$
- **what:** multiple-histogram (WHAM) reweighting: expectation of $A$ at any intermediate temperature $T = 1/k_B\beta$ from energy histograms.
- **symbols:** $A(E)$ - observable binned by energy; $P(E;\beta)$ - reweighted density of states (eq:21); sums run over energy bins $E$.

<!-- eq:21 -->
$$ P(E;\beta) = \frac{\sum_{m=1}^{R} g_m^{-1}\, N_m(E)\, e^{-\beta E}}{\sum_{m=1}^{R} n_m\, g_m^{-1}\, e^{f_m - \beta_m E}} $$
- **what:** WHAM estimate of the (reweighted) density of states from $R$ runs' histograms. Solved self-consistently with eq:22.
- **symbols:** $R$ - number of independent runs (temperatures); $N_m(E)$ - energy histogram of run $m$; $n_m$ - total samples in run $m$ (in REM, $n_m = N_{\text{sim}}$); $g_m = 1 + 2\tau_m$ - statistical inefficiency; $\tau_m$ - integrated autocorrelation time; $f_m$ - dimensionless free energy of run $m$ (eq:22); $\beta_m$ - inverse temperature of run $m$.

<!-- eq:22 -->
$$ e^{-f_m} = \sum_{E} P(E; \beta_m) $$
- **what:** self-consistency (normalization) condition fixing the free energies $f_m$. Iterate eq:21 and eq:22 to convergence (10-100 iterations for Met-enkephalin).
- **symbols:** $f_m$ - dimensionless free energy of run $m$; $P(E;\beta_m)$ - density of states at $\beta_m$ from eq:21.

<!-- eq:g_m -->
$$ g_m = 1 + 2\tau_m $$
- **what:** statistical inefficiency of run $m$ in terms of the integrated autocorrelation time. For biomolecular systems the $\tau_m$ are approximately equal, so $g_m = \text{const}$ may be used.
- **symbols:** $\tau_m$ - integrated autocorrelation time at temperature $T_m$.
