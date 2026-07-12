# Notation - Sugita & Okamoto 1999 (Replica-Exchange MD)

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $N$ | number of atoms in one replica | integer | |
| $M$ | number of replicas = number of temperatures | integer | exactly one replica per temperature |
| $q = \{q_1,\dots,q_N\}$ | atomic coordinates | length (Cartesian), $\mathbb{R}^{3N}$ | Cartesian coordinates |
| $p = \{p_1,\dots,p_N\}$ | atomic momenta | momentum, $\mathbb{R}^{3N}$ | |
| $m_k$ | mass of atom $k$ | mass | |
| $x = (q,p)$ | phase-space microstate | | |
| $H(q,p)$ | Hamiltonian | energy | $H = K + E$ |
| $K(p)$ | kinetic energy | energy | $\sum_k p_k^2 / 2m_k$ |
| $E(q)$ | potential energy | energy (kcal/mol in tests) | AMBER all-atom force field, dielectric = 1 |
| $T, T_m$ | temperature (of slot $m$) | K | $\beta_1<\beta_2<\dots<\beta_M$ i.e. $T_1>T_2>\dots>T_M$ |
| $\beta = 1/k_B T$ | inverse temperature | 1/energy | $\beta_m = 1/k_B T_m$ |
| $k_B$ | Boltzmann constant | energy/K | |
| $W_B(x;T)$ | canonical Boltzmann weight | dimensionless | $e^{-\beta H}$ |
| $i$ | replica label | integer $1..M$ | superscript $[i]$ in $x_m^{[i]}$ |
| $m$ | temperature label / slot | integer $1..M$ | subscript in $x_m^{[i]}$ |
| $f, f^{-1}$ | permutation functions replica<->temperature | | $i=f(m)$, $m=f^{-1}(i)$; updated on accepted swap |
| $f^{-1}(i;t)$, $f(m;t)$ | time-dependent permutation at measurement $t$ | | tracked per node in parallel implementation |
| $x_m^{[i]}$ | microstate of replica $i$ at temperature $T_m$ | | superscript=replica, subscript=temperature |
| $X$ | full generalized-ensemble state (all $M$ replicas) | | |
| $W_{\text{REM}}(X)$ | weight of full state | dimensionless | product of per-replica Boltzmann factors (a priori known) |
| $p^{[i]'}$ | rescaled momentum after swap | momentum | $\sqrt{T_n/T_m}\,p^{[i]}$ (velocity scaling) |
| $\Delta$ | reduced swap-acceptance cost | dimensionless | $(\beta_n-\beta_m)(E(q^{[i]})-E(q^{[j]}))$ |
| $w(X\to X')$ | swap acceptance probability | dimensionless | Metropolis criterion |
| $A$ | physical observable | | |
| $\langle A\rangle_{T_m}$ | canonical expectation at $T_m$ | | arithmetic mean over samples |
| $N_{\text{sim}}$ | measurements per replica | integer | $=10^5$ in the test |
| $N_m(E)$ | energy histogram of run $m$ | counts | |
| $n_m$ | total samples in run $m$ | integer | $n_m = N_{\text{sim}}$ in REM |
| $P(E;\beta)$ | reweighted density of states | | WHAM, self-consistent |
| $f_m$ | dimensionless free energy of run $m$ | dimensionless | self-consistent normalization |
| $g_m = 1+2\tau_m$ | statistical inefficiency of run $m$ | dimensionless | set const for biomolecules |
| $\tau_m$ | integrated autocorrelation time at $T_m$ | steps/time | |
| $R$ | number of independent reweighting runs | integer | $R = M$ in the test |
| $\delta_{k,l}$ | Kronecker delta | | |

## Conventions

- Temperatures are ordered so that $\beta_1 < \beta_2 < \dots < \beta_M$, i.e.
  $T_1 > \dots > T_M$ (index 1 is hottest). Only neighboring temperature pairs
  ($n = m\pm1$) are exchanged because acceptance decays exponentially with
  $|\beta_n - \beta_m|$.
- Momentum rescaling on a swap uses the square root of the temperature ratio; this
  is exactly a uniform velocity rescale and is what makes kinetic energy cancel in
  the acceptance ratio, reducing the MD-REM criterion to the MC-REM criterion
  (potential energy only).
- The swap criterion $\Delta$ (eq:15) contains no momentum/kinetic term: only
  potential energies and inverse temperatures enter. An implementation exchanges
  temperature labels (or, equivalently, configurations + rescaled velocities).
