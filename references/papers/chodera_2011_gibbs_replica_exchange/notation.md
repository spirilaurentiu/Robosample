# Notation - Chodera & Shirts 2011

Reduced/dimensionless convention: all state energetics enter through the reduced potential $u_k(x)$, which is in units of $k_B T$. Acceptance ratios and weights are therefore dimensionless. Log-weights $g_k$ are in units of $k_B T$.

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| $x$ | microscopic configuration | point in $\Omega$ (continuous or discrete) | full coordinate block |
| $\Omega$ | accessible configuration space | set | may be continuous or discrete |
| $k, i, j$ | thermodynamic state index | integer $\in\{1,\dots,K\}$ | discrete states |
| $K$ | number of thermodynamic states | integer | |
| $\lambda$ | vector of thermodynamic parameters | $\{\beta,H,p,\mu,\dots\}$ | parameterizes a state |
| $\lambda_k$ | parameters of state $k$ | vector | |
| $u(x), u_k(x)$ | reduced potential | dimensionless ($k_BT$ units) | Eq. 1; $u_k$ for state $k$ |
| $\beta$ | inverse temperature | $1/(k_B T)$ | $\beta=(k_B T)^{-1}$ |
| $H(x)$ | Hamiltonian / potential energy | energy | |
| $p$ | pressure | pressure | conjugate to $V$ |
| $V(x)$ | volume | volume | constant-pressure ensemble |
| $\mu_i$ | chemical potential of component $i$ | energy | conjugate to $n_i$ |
| $n_i(x)$ | number of molecules of component $i$ | integer | (semi)grand ensemble |
| $q(x), q_k(x)$ | unnormalized density | dimensionless | $q=e^{-u}$ (Boltzmann) |
| $\pi(x), \pi_k(x)$ | normalized density | probability density | $\pi=Z^{-1}q$ |
| $\pi(x,k)$ | joint config-state density | | expanded ensemble |
| $\pi(k|x)$ | conditional state given config | probability over $K$ states | Eq. 7 - state-update target |
| $Z, Z_k$ | partition function | dimensionless | Eq. 2 |
| $g_k$ | state log-weight (bias) | $k_B T$ units | $g_k=-\ln Z_k$ for equal-probability states |
| $X$ | vector of replica configurations | $\{x_1,\dots,x_K\}$ | replica exchange |
| $S$ | permutation of state indices | $\{s_1,\dots,s_K\}\in\mathcal{S}_K$ | maps replica->state |
| $s_i$ | state assigned to replica $i$ | integer | |
| $\mathcal{S}_K$ | symmetric group ($K!$ permutations) | set | |
| $\mathbf{U}$ | reduced-potential matrix | $K\times K$, $u_{ij}=u_i(x_j)$ | precomputed for swap loop |
| $\alpha(j|x,i)$ | proposal probability $i\to j$ | probability | state move |
| $P_{\text{accept}}$ | Metropolis acceptance probability | $\in[0,1]$ | |
| $\mathcal{S}_i$ | restricted proposal set for state $i$ | subset of $\{1,\dots,K\}$ | symmetric: $i\in\mathcal{S}_j\iff j\in\mathcal{S}_i$ |
| $n$ | restricted-range half-width | integer, $n\ll K$ | $\mathcal{S}_i=\{i-n,\dots,i+n\}$ |
| $\mathbf{T}$ | empirical state transition matrix | $K\times K$ row-stochastic | symmetrized (Eq. 25) |
| $N_{ij}$ | count of $i\to j$ transitions | integer | pooled over replicas |
| $\mu_1\ge\dots\ge\mu_K$ | eigenvalues of $\mathbf{T}$ | real, $\mu_1=1$ | $\mu_2=1\Rightarrow$ decomposable |
| $\tau$ | time between exchange attempts | time (ps) | |
| $\tau_2$ | relaxation time from $\mu_2$ | time (ps) | $\tau/(1-\mu_2)$; lower bound |
| $\tau_{ac}$ | integrated autocorrelation time of state index | time (ps) | statistical inefficiency $=2\tau_{ac}+1$ |
| $\tau_{end}$ | avg end-to-end transit time of state index | time (ps) | $k=1\leftrightarrow k=K$ |
| $\tau_N$ | autocorr time of coordination number $N$ | time (ps) | structural probe |
| $U(x)$ | 1D model potential | energy | $10(x-1)^2(x+1)^2$ |
| $\lambda$ (alchemical) | coupling parameter | $\in[0,1]$ | $1$=coupled, $0$=decoupled |
| $\alpha$ (softcore) | softcore LJ parameter | dimensionless | Eq. 30 |
| $\phi,\psi$ | backbone torsion angles | radians | alanine dipeptide slow DOF |
| $\kappa$ | umbrella force constant | energy | von Mises bias (Eq. 31) |
