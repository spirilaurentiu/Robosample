# Checks / fixtures - Chodera & Shirts 2011

Concrete numbers an implementation can regression-test against. Mixing times in ps; "relative speedup" is relative to neighbor exchange (=1.0). Uncertainties are one standard error of the mean.

## Algebraic / structural invariants (implementable directly)

- **Independence-sampling CDF selection:** draw uniform $r\in[0,1)$; select smallest $k$ with $r<\sum_{i=1}^{k}\pi(i|x)$. Always accepted.
- **Restricted range reduces to independence sampling** when $\mathcal{S}_i=\{1,\dots,K\}$ for all $i$ (all proposals accepted).
- **Metropolized independence sampling == restricted range** with $\mathcal{S}_i=\{1,\dots,K\}\setminus\{i\}$ (current state excluded, $\alpha(i|x,i)=0$).
- **Transition matrix eigenvalue ordering:** $1=\mu_1\ge\mu_2\ge\cdots\ge\mu_K$. $\mu_2=1\Rightarrow$ decomposable chain (disjoint state sets, no observed transitions between them).
- **Statistical inefficiency** of a Markovian series $=2\tau_{ac}+1$ (samples per uncorrelated sample).
- **Aldous-Diaconis:** if all $u_{ij}$ equal, number of random pair swaps to mix a $K$-permutation is $\sim K\ln K$.
- **Swap-loop empirical guidance:** $K^3$ to $K^5$ random pair swaps per state-update iteration were sufficient to effectively decorrelate the permutation for the molecular systems studied.
- **Permutation denominator** $\pi(S|X)$ normalizer is the matrix permanent of $\mathbf{A}=(a_{ij})$, $a_{ij}=e^{-u_i(x_j)}$ (#P-complete; intractable for large $K$).
- **$\tau_2$ scales linearly with update interval:** infrequent (5 ps) mixing times were exactly $50\times$ = (5 ps / 0.1 ps) the frequent-update $\tau_2$.

## 1D double-well simulated tempering (Sec IV)

Fixture: potential $U(x)=10(x-1)^2(x+1)^2$, $K$ temperatures geometric $\beta_k=10^{-(k-1)/(K-1)}$ over $k_BT\in[1,10]$. Start $(x_0,k_0)=(-1,1)$. Each iteration: 1 state update + 100 Metropolis MC steps (Gaussian proposal, mean 0, std 0.1). $10^6$ iterations.

- $K=16$, **neighbor exchange:** $\tau_x = 24.1 \pm 0.9$ iterations.
- $K=16$, **independence sampling:** $\tau_k = 0.243 \pm 0.004$, $\tau_x = 9.6 \pm 0.2$ iterations.
- Increasing $K$ increases both $\tau_k$ and $\tau_x$ for neighbor exchange; independence sampling keeps both small as $K$ grows.

## UA methane, expanded ensemble alchemical (Table I)

$K=6$ softcore-LJ states, $\lambda_k=\{0.0,0.3,0.6,0.7,0.8,1.0\}$. Perfect log-weights $g_k=\{0.0,0.32,-0.46,-1.67,-2.83,-3.66\}$ (units $k_BT$). $\sigma$(methane)=0.373 nm, $\epsilon=1.230096$ kJ/mol, methane-water $\sigma_{ij}=0.3428$ nm. 893 TIP3P waters, 298 K, 2 fs step. Coordination number $N$ = O atoms within 0.3 nm.

| scheme | interval | $\tau_2$ | $\tau_{ac}$ | $\tau_{end}$ | $\tau_N$ | speedup $\tau_2$ |
|---|---|---|---|---|---|---|
| Neighbor | 0.1 ps | 1.693±0.008 | 6.7±0.4 | 11.9±0.2 | 5.9±0.4 | 1.0 |
| Independence | 0.1 ps | 0.771±0.004 | 6.2±0.2 | 7.2±0.1 | 5.4±0.2 | 2.20±0.02 |
| Metropolized indep | 0.1 ps | 0.645±0.003 | 4.6±0.2 | 6.6±0.1 | 4.4±0.2 | 2.62±0.02 |
| Neighbor | 5 ps | 85.8±2.3 | 177.7±17.6 | 330.0±16.1 | 105.3±12.1 | 1.0 |
| Independence | 5 ps | 39.0±0.9 | 69.2±6.1 | 141.1±4.7 | 49.1±3.8 | 2.20±0.08 |
| Metropolized indep | 5 ps | 31.8±0.4 | 51.4±1.9 | 115.7±3.4 | 37.4±1.4 | 2.70±0.08 |

- 1000 state moves per 0.1 ps: all three schemes converge to independence sampling (speedups ~1.0, $\tau_2\approx0.77$).

## Large LJ sphere, expanded ensemble (Table II)

$K=18$, $\lambda=[0,0.15,0.3,0.45,0.55,0.6,0.64,0.66,0.68,0.70,0.72,0.75,0.78,0.81,0.84,0.87,0.90,1.0]$. $\sigma_{ii}=1.09$ nm, $\epsilon_{ii}=1.230096$ kJ/mol, sphere-water $\sigma_{ij}=0.561$ nm (~5x methane volume). $g_k=\{0.0,1.74,2.96,3.39,2.84,2.01,0.73,-0.34,-1.75,-3.35,-4.96,-7.19,-9.11,-10.70,-11.98,-12.98,-13.72,-14.65\}$. $N$ = O atoms within 0.5 nm.

| scheme | interval | $\tau_2$ | $\tau_{ac}$ | $\tau_{end}$ | $\tau_N$ | speedup $\tau_2$ |
|---|---|---|---|---|---|---|
| Neighbor | 0.1 ps | 9.51±0.01 | 65.8±4.2 | 126.3±4.2 | 58.1±4.3 | 1.0 |
| Independence | 0.1 ps | 2.586±0.009 | 42.9±2.4 | 88.4±2.7 | 41.5±2.0 | 3.68±0.01 |
| Metropolized indep | 0.1 ps | 2.181±0.006 | 48.6±4.0 | 88.3±3.0 | 46.7±3.4 | 4.36±0.01 |
| Neighbor | 1 ps | 95.0±0.2 | 211.1±58.9 | 507.6±19.3 | 167.6±16.0 | 1.0 |
| Independence | 1 ps | 25.8±0.1 | 67.3±3.6 | 196.0±5.8 | 63.1±3.3 | 3.69±0.02 |
| Metropolized indep | 1 ps | 21.6±0.1 | 66.8±2.4 | 169.2±4.7 | 62.1±2.5 | 4.40±0.02 |

## Parallel tempering, alanine dipeptide implicit solvent (Table III)

Ace-Ala-Nme, AMBER parm96, OBC GBSA (igb=2), OpenMM. 2000 iterations, 500 MD steps/iter, 2 fs step. Replica-mix: $K^3$ random pair swaps. Velocities reassigned Maxwell-Boltzmann each iteration.

| scheme | $\tau_2$ | $\tau_{ac}$ | $\tau_{end}$ | $\tau_{\cos\phi}$ | $\tau_{\sin\phi}$ | $\tau_{\cos\psi}$ | $\tau_{\sin\psi}$ |
|---|---|---|---|---|---|---|---|
| Neighbor | 91.8±0.6 | 80±2 | 360±30 | 25±2 | 110±9 | 25±2 | 66±6 |
| Independence | 2.62±0.01 | 1.60±0.06 | 28.7±0.7 | 12.4±0.5 | 8.7±0.4 | 11.8±0.6 | 9.1±0.5 |

State-space mixing ~1-2 orders of magnitude faster; structural ~2-10x.

## 2D REX umbrella sampling, alanine dipeptide (Table IV)

$K=101$ replicas on a $10\times10$ toroidal grid + 1 unbiased. von Mises bias, $\kappa=(2\pi/30)^{-2}\beta^{-1}$ (neighbors separated by $3\sigma$). 300 K, 2 fs step, 5 ps between exchanges, 2000 iterations (first 100 discarded).

| scheme | $\tau_2$ | $\tau_{ac}$ | $\tau_{end}$ | $\tau_{\cos\phi}$ | $\tau_{\sin\phi}$ | $\tau_{\cos\psi}$ | $\tau_{\sin\psi}$ |
|---|---|---|---|---|---|---|---|
| Neighbor | 82±4 | 31.0±0.9 | 350±30 | 47±2 | 57±2 | 26.4±0.8 | 27.1±0.9 |
| Independence | 24.2±0.3 | 5.45±0.06 | 175±6 | 8.92±0.09 | 9.9±0.1 | 5.63±0.04 | 6.09±0.04 |

State relaxation reduced 2-6x; structural correlation reduced 4-5x.
