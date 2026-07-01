# Notation - Pear & Weiner 1979

**Global conventions.** Reduced units: carbon-atom mass $m$, C-C bond length $l$, and
reference temperature $T_{\mathrm{ref}}=600$ K are all set to 1. $\beta = 1/(k_B T)$.
Chain has $N$ bonds, $N+1$ atoms $C_0,\dots,C_N$, $N-1$ rigid bodies. Generalized
coordinates $\phi_0,\dots,\phi_N$: $\phi_0,\phi_1,\phi_2$ are Bryant (Euler-like)
orientation angles of the first body; $\phi_3,\dots,\phi_N$ are internal dihedral angles.
$\phi_i = 0$ = trans conformation. Vectors bold; tensors rank 0/1/2 in underlined matrix
elements. Cross product $\times$; dyadic (outer product) written as adjacency $\mathbf{ab}$.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $m$ | atom mass (carbon) | reduced (=1) | interior atoms; ends scaled by $\alpha$ |
| $l$ | C-C bond length | reduced (=1) | fixed in rigid model |
| $\theta$ | equilibrium valence angle | rad | 90 deg for all reported results |
| $\theta_i$ | valence angle at $C_i$ | rad | angle between $\mathbf{l}_i,\mathbf{l}_{i-1}$; $i=2,\dots,N$ |
| $\phi_i$ | dihedral angle | rad | rotation of $\mathbf{l}_i$ about $\mathbf{l}_{i-1}$; $\phi=0$ trans; $i=3,\dots,N$ |
| $\phi_0,\phi_1,\phi_2$ | Bryant orientation angles | rad | x-, then new-y, then new-z rotations of first body |
| $\alpha$ | end/interior atom mass ratio | dimensionless | mass of $C_0$/$C_1$ (and $C_3$/$C_2$); $\alpha=1$ or $10$ used |
| $N$ | number of bonds | int | 3 for reported results |
| $M$ | total chain mass | reduced | $\sum m_r$ |
| $\mathbf{x}_i$ | Cartesian position of atom $C_i$ | R^3, reduced | |
| $\mathbf{l}_i$ | bond vector $\mathbf{x}_i-\mathbf{x}_{i-1}$ | R^3 | $l_i=|\mathbf{l}_i|$ |
| $\mathbf{R}$ | end-to-end vector | R^3 | $R(\phi)=l(3+2\cos\phi)^{1/2}$ for 3-bond 90-deg chain |
| $\eta$ | viscosity / damping coefficient | reduced | Langevin friction; reported as $\eta/m\omega$ |
| $\mathbf{L}(t)$ | Langevin fluctuating force | R^3 | $\langle L_i L_j\rangle = 2\eta k_B T\,\delta(t-t')\delta_{ij}$ |
| $\mathbf{B}(\delta t)$ | Langevin impulse over $\delta t$ | R^3 | variance $2\eta k_B T\,\delta t$ per component |
| $\gamma_i^{(j)}$ | random draw | standard normal | mean 0, std 1 |
| $\delta t$ | integration time step | reduced ($1/\omega$ units) | $\delta t=0.1/\omega,0.05/\omega,0.025/\omega$ |
| $\omega$ | reference frequency | $1.33\times10^{13}$ s$^{-1}$ | sets time-step scale |
| $T$ | temperature | K (or reduced via $T_{\mathrm{ref}}$) | |
| $k_B$ | Boltzmann constant | | $\beta=1/(k_B T)$ |
| $\underline{T}$ | connectivity matrix | $(N-1)\times(N-1)$ int | $T_{ij}=-1$ if $C_j$ on path (Eq. 2.1) |
| $\underline{A}$ | generalized mass matrix | $(N+1)\times(N+1)$ | symmetric positive-definite |
| $\underline{B}$ | generalized force vector | $N+1$ | RHS of eq. of motion |
| $\underline{\mathbf{p}}$ | hinge-axis matrix | $(N+1)\times(N-1)$ of R^3 | Eq. 2.9 |
| $\mathbf{p}_k$ | unit hinge/Bryant axis | R^3 unit | rotation axis for $\phi_k$ |
| $\underline{K}$, $\mathbf{K}_{ij}$ | inertia-tensor matrix | tensors R^{3x3} | Eq. 2.10 |
| $\mathbf{K}_i^*$ | augmented-body inertia tensor | R^{3x3}, about barycenter | Eqs. A5, A7 |
| $\mathbf{b}_{j0},\mathbf{b}_{jj},\mathbf{b}_{jN}$ | body vectors (lower hinge, CoM, upper hinge) | R^3 | rel. barycenter of $B_j$; $\mathbf{b}_{10}=\mathbf{b}_{N-1,N}=0$ |
| $\boldsymbol{\rho}_i$ | vector body-CoM to end atom | R^3 | terminal-body corrections |
| $\Omega_i$ | relative angular velocity | R^3 | body $C_i$ rel. $C_{i-1}$ |
| $\omega_i$ | absolute angular velocity | R^3 | $\sum_{j\le i}\Omega_j$ |
| $\tau_i$ | internal hinge torque | scalar along $\mathbf{p}_i$ | $\partial V_\phi/\partial\phi_i$; 0 for $i=0,1,2$ |
| $\underline{\mathbf{f}}$ | kinematic term | R^3 per body | Eq. 2.11 |
| $\underline{\mathbf{M}}'$ | gyroscopic moments | R^3 per body | Eq. 2.12 |
| $\underline{\mathbf{M}}_\eta$ | damping moments | R^3 per body | Eqs. 2.17-2.20 |
| $\mathbf{E}$ | 3D identity tensor | R^{3x3} | |
| $g(\phi_3)$ | metric determinant | scalar | $=|G_{ij}|$; Eq. 3.6 closed form |
| $G_{ij}$ | covariant metric tensor | R^{4x4} | Eq. 4.5; mass-weighted |
| $G^{ij}$ | contravariant metric tensor | $(G_{ij})^{-1}$ | $G^{33}$ used in rate formulas |
| $H_{ij}$ | Fixman constraint matrix | R^{5x5} | Eq. 3.5; $g\propto1/|H|$ |
| $c_i$ | constrained coordinate | | bond lengths + valence-angle proxies |
| $U$ | Fixman compensating potential | energy | $k_B T\ln\sqrt{g}$ |
| $V(\phi_3)$ | rotational barrier potential | energy | Eq. 4.1, single barrier |
| $E_b$ | barrier height | energy | $=k/4$ |
| $k$ | barrier stiffness | | Eq. 4.1 |
| $\phi_b$ | barrier-peak position | rad | $\phi_b=0$ or $\pi$ |
| $\sigma$, $\boldsymbol{\sigma}$ | applied end force (magnitude / vector) | force | $\pm\boldsymbol{\sigma}$ on $C_0,C_3$; reported as $\sigma l/E_b$ |
| $\Theta$ | angle between $\mathbf{R}$ and $\boldsymbol{\sigma}$ | rad | |
| $\Phi,\Psi$ | azimuth of $\mathbf{R}$ about $\sigma$ / spin about $\mathbf{R}$ | rad | orientation coords under stress |
| $f_{TS}$ | transition-state crossing rate | 1/time | Eqs. 4.10, 4.11, 5.9 |
| $I_c$ | configurational normalization | | Eqs. 4.8, 4.11, 5.10 |
| $k_l,k_\theta$ | bond / angle spring constants (flexible) | | large -> rigid limit |
| $T(\phi)$ (matrix) | inter-body transform | R^{3x3} | Eq. A2; distinct from kinetic energy $T$ |
| $a^{01}$ | first-body-to-reference transform | R^{3x3} | Eq. A3 |
