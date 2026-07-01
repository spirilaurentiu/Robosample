# Notation — Jain 1997, Compensating Mass Matrix Potential

Spatial operator algebra (SOA) for tree-topology multibody / internal-coordinate
molecular dynamics. Serial chain of $n$ clusters used for derivation; a
single-dof rotational hinge between clusters unless stated. Base cluster carries
6 dof, so a serial chain has $\mathcal{N}=n+5$ dof.

Convention: $*$ denotes transpose / adjoint. $\tilde v$ is the $3\times3$ skew
(cross-product) matrix of a 3-vector $v$; $v\times w = \tilde v\, w$. "outboard"
= toward the tips, away from the base. Operators $\phi,\psi,K,\mathcal{E}$ are
$6n\times6n$ block; $H$ maps spatial ($6n$) to hinge ($\mathcal{N}$) space.

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| $\boldsymbol{\theta}$ | internal (hinge) configuration coordinates | vector, $\mathbb{R}^{\mathcal{N}}$ | generalized coordinates |
| $\dot{\boldsymbol{\theta}}$ | generalized velocities | $\mathbb{R}^{\mathcal{N}}$ | |
| $\ddot{\boldsymbol{\theta}}$ | generalized accelerations | $\mathbb{R}^{\mathcal{N}}$ | solved by $O(\mathcal{N})$ algo |
| $\beta$ | internal velocity coordinates vector | $\mathbb{R}^{\mathcal{N}}$ | $=\dot\theta$ up to hinge map |
| $p$ | conjugate momenta | $\mathbb{R}^{\mathcal{N}}$ | $p=\mathcal{M}\beta$ |
| $q_i$ | Cartesian config coordinate (unconstrained model) | | |
| $n$ | number of clusters / hinges (serial chain) | int | |
| $\mathcal{N}$ | number of constrained degrees of freedom | int | $=n+5$ for serial single-dof chain |
| $\mathcal{M}(\boldsymbol{\theta})$ | system mass matrix / metric tensor | $\mathbb{R}^{\mathcal{N}\times\mathcal{N}}$ | symmetric PD; config-dependent (constrained), constant+diagonal (unconstrained) |
| $\mathcal{V}$ | standard potential energy | energy | |
| $\mathcal{V}_c$ | compensating (Fixman / metric-tensor) potential | energy | $=\tfrac12\ln\det\mathcal{M}=\tfrac12\sum_i\ln\det D(i)$ |
| $\mathcal{V}'$ | modified potential $=\mathcal{V}+\mathcal{V}_c$ | energy | used for forces in constrained MD |
| $T$ | standard hinge torque $=\nabla_\theta\mathcal{V}$ | $\mathbb{R}^{\mathcal{N}}$ | |
| $T_c$ | compensating hinge torque $=\nabla_\theta\mathcal{V}_c$ | $\mathbb{R}^{\mathcal{N}}$ | $T_c=0$ on the 6 base-cluster dof |
| $T'$ | total hinge torque $=T+T_c$ | $\mathbb{R}^{\mathcal{N}}$ | |
| $\mathcal{C}(\theta,\dot\theta)$ | Coriolis/centrifugal/gyroscopic + Cartesian force vector | $\mathbb{R}^{\mathcal{N}}$ | in $T=\mathcal{M}\ddot\theta+\mathcal{C}$ |
| $\mathscr{Z},\mathcal{Z}'$ | partition function (unconstrained / constrained) | scalar | |
| $\mathscr{T}$ | temperature | K | |
| $k$ | Boltzmann constant | — | appears as $k\mathscr{T}$ |
| $h$ | Planck constant (partition-fn prefactor $h^{-n}$) | — | NOT the hinge axis $h(i)$ |
| $H$ | joint-map / hinge axis operator | $\mathcal{N}\times 6n$ | block; $H(k)$ picks hinge dof |
| $\phi$ | rigid-body transform propagation operator | $6n\times6n$ | block lower triangular; $\phi(k,k-1)$ adjacent transform |
| $M$ | spatial (link) inertia operator | $6n\times6n$ block diagonal | blocks $M(k)\in\mathbb{R}^{6\times6}$ |
| $P(k)$ | articulated body inertia outboard of hinge $k$ | $\mathbb{R}^{6\times6}$ | from Riccati recursion Eq. (3.6) |
| $P^+(k)$ | articulated inertia after projecting out hinge $k$ | $\mathbb{R}^{6\times6}$ | $=\bar\tau(k)P(k)$ |
| $D(k)$ | hinge inertia | $\mathbb{R}^{d_k\times d_k}$ | $=H(k)P(k)H^*(k)$; $d_k$=hinge dof (1 for single-dof) |
| $G(k)$ | Kalman gain block | $6\times d_k$ | $=P(k)H^*(k)D^{-1}(k)$ |
| $K$ | Kalman gain operator | $6n\times\mathcal{N}$ | nonzero blocks $K(k,k-1)$ on first subdiagonal |
| $\tau(k),\bar\tau(k)$ | hinge projection / complementary projection | $6\times6$ | $\bar\tau=I-\tau$, $\tau=G(k)H(k)$ |
| $\psi$ | articulated-body transform propagation operator | $6n\times6n$ | $=(I-\mathcal{E}_\psi)^{-1}$; $\psi(k+1,k)=\phi(k+1,k)\bar\tau(k)$ |
| $\tilde\psi$ | $\psi-I$ | | strictly lower triangular |
| $\mathcal{E}_\phi,\mathcal{E}_\psi$ | shift operators | $6n\times6n$ | $\mathcal{E}_\psi=\mathcal{E}_\phi\bar\tau$ |
| $\Omega$ | articulated-body operator | $6n\times6n$ | $=\psi^*H^*D^{-1}H\psi$ |
| $Y(k)$ | block-diagonal element for $T_c$ recursion | $\mathbb{R}^{6\times6}$ | from Eq. (4.16) |
| $\mathbb{H}_\delta^i$ | single-block selector at hinge $i$ | $6n\times6n$ | zero except block $\mathbb{H}(i)$ at $i$th diagonal |
| $\mathbb{H}(i)$ | $6\times6$ block $\operatorname{diag}(\tilde h(i),\tilde h(i))$ | $\mathbb{R}^{6\times6}$ | built from hinge axis |
| $h(i)$ | hinge $i$ rotational axis unit vector | $\mathbb{R}^{3}$, unit | rotational hinge |
| $\tilde h(i)$ | skew matrix of $h(i)$ | $\mathbb{R}^{3\times3}$ | |
| $\mathcal{M}_{\theta_i}$ | mass-matrix sensitivity $\partial\mathcal{M}/\partial\theta_i$ | $\mathbb{R}^{\mathcal{N}\times\mathcal{N}}$ | shorthand $\mathcal{M}_{\theta(k)}$ |
| $\mathcal{F}[A]$ | axial-vector map | $\mathbb{R}^{3\times3}\!\to\!\mathbb{R}^3$ | $v=\mathcal{F}[A]\iff\tilde v=A-A^*$ |
| $Q_{jk}$ | $3\times3$ blocks of $P(i)Y(i)$ | $\mathbb{R}^{3\times3}$ | |
| $V(k),\alpha(k)$ | spatial velocity / acceleration of link $k$ | $\mathbb{R}^{6}$ | |
| $a(k),b(k)$ | Coriolis/velocity accel term / gyroscopic spatial force | $\mathbb{R}^{6}$ | |
| $\hat f_c(k)$ | applied Cartesian spatial force on link $k$ | $\mathbb{R}^{6}$ | |
| $z(k),z^+(k)$ | residual spatial forces ($O(\mathcal{N})$ solve) | $\mathbb{R}^{6}$ | |
| $\varepsilon(k),\nu(k)$ | innovation / $D^{-1}$-weighted innovation | $\mathbb{R}^{d_k}$ | |

## Sign / convention notes

- Compensating potential adds to the potential: $\mathcal{V}'=\mathcal{V}+\mathcal{V}_c$, and the extra force is the gradient $T_c=\nabla_\theta\mathcal{V}_c$ (Eq. 2.10). In an MD force = $-\nabla\mathcal{V}'$ convention, subtract $T_c$; here $T$ is defined as the gradient (torque), so signs follow Eq. (2.10) directly.
- $T_c$ on the base cluster's 6 rigid-body dof is identically zero.
- $\det\{\mathcal{M}\}=\det\{D\}=\prod_i\det\{D(i)\}$ (from innovations factorization with unit-diagonal triangular factor), so $\mathcal{V}_c$ is a sum of small log-dets — the key computational reduction.
- Hinge inertia $D(k)$ blocks are $d_k\times d_k$; for the single-dof rotational hinges used in the serial-chain derivation $D(k)$ is a scalar and $\ln\det D(k)=\ln D(k)$.
- Multi-dof hinges are modeled as sequences of single-dof hinges joined by zero-mass/zero-inertia pseudo-clusters.
