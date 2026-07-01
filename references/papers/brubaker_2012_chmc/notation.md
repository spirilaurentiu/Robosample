# Notation - brubaker_2012_chmc

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $q$ | configuration / state on manifold | $\mathbb{R}^n$ | must satisfy $c(q)=0$ |
| $p$ | momentum | $\mathbb{R}^n$ | conjugate to $q$; $p=\partial\mathcal{L}/\partial\dot q$ |
| $\dot q$ | velocity | $\mathbb{R}^n$ | $\dot q = M^{-1}p = \partial\mathcal{H}/\partial p$ |
| $\mathcal{M}$ | constraint manifold | $\{q\in\mathbb{R}^n \mid c(q)=0\}$ | connected, smooth, differentiable |
| $c(q)$ | constraint function | $\mathbb{R}^n\to\mathbb{R}^m$ | $c(q)=0$ defines $\mathcal{M}$ |
| $C(q)$ | constraint Jacobian $\partial c/\partial q$ | $m\times n$ | full rank everywhere |
| $\mathcal{T}\mathcal{M}$ | tangent bundle | $\{(q,\dot q)\mid c(q)=0,\ C(q)\dot q=0\}$ | |
| $\mathcal{T}_q\mathcal{M}$ | tangent space at $q$ | $\{\dot q\mid C(q)\dot q=0\}$ | |
| $\mathcal{T}^*\mathcal{M}$ | cotangent bundle | $\{(p,q)\mid c(q)=0,\ C(q)M^{-1}p=0\}$ | momentum phase space |
| $\mathcal{T}^*_q\mathcal{M}$ | cotangent space at $q$ | $\{p\mid C(q)M(q)^{-1}p=0\}$ | momentum lives here |
| $M(q)$ | mass matrix | $n\times n$, SPD | symmetric positive definite; may depend on $q$ |
| $T(p,q)$ | kinetic energy | scalar | $\frac12 p^T M(q)^{-1}p$ |
| $U(q)$ | potential energy (acceptance) | scalar | $\frac12\log|M(q)| - \log\pi(q)$ |
| $\hat U(q)$ | potential energy (guidance) | scalar | free choice for simulation Hamiltonian |
| $\mathcal{H}$ | acceptance Hamiltonian | scalar | used in Metropolis test; based on $\pi$ |
| $\hat{\mathcal{H}}$ | guidance Hamiltonian | scalar | used for simulation; must be $\mathcal{C}^2$ |
| $\mathcal{L}$ | continuous Lagrangian | scalar | $T-U-\lambda^T c(q)$ |
| $\mathcal{L}'_h$ | discrete Lagrangian of integrator | $\mathcal{M}\times\mathcal{M}\to\mathbb{R}$ | |
| $\lambda$ | Lagrange multiplier, state constraint | $\mathbb{R}^m$ | enforces $c(q_1)=0$ |
| $\mu$ | Lagrange multiplier, momentum constraint | $\mathbb{R}^m$ | enforces $C(q_1)\dot q_1=0$ |
| $\pi(q)$ | target (unnormalized) density | scalar $\ge 0$ | $\int_\mathcal{M}\pi\,dq=1$; strictly positive for convergence |
| $\bar\pi(p,q)$ | augmented density | scalar | $\exp(-\mathcal{H}(p,q))=\pi(q)\mathcal{N}(p|0,M)$ |
| $\Phi_h^{\mathcal{H}}$ | numerical integrator | $\mathcal{T}^*\mathcal{M}\to\mathcal{T}^*\mathcal{M}$ | symmetric, symplectic, consistent (RATTLE, order 2) |
| $h$ | integration step size | scalar | must be small enough for convergence proofs |
| $L$ | number of leapfrog/RATTLE steps per proposal | integer | $L=1$ -> Langevin; $L=1,h=1$ -> Metropolis |
| $r$ | integrator order | integer | consistent if $r\ge1$; RATTLE $r=2$ |
| $e_h$ | bounded error function | scalar | in consistency condition eq:2 |
| $\rho$ | momentum-flip map | $(p,q)\mapsto(-p,q)$ | dynamics are $\rho$-reversible |
| $J$ | canonical symplectic matrix | $2n\times2n$ | $[[0,I],[-I,0]]$ |
| $T(q\to q')$ | transition kernel | scalar | Markov chain |
| $\mathcal{N}(\cdot|\mu,\Sigma)$ | Gaussian density | | mean $\mu$, covariance $\Sigma$ |
| $d, A$ | BvMF location vector / spread matrix | $\mathbb{R}^n$ / $n\times n$ | $A$ defined up to $A+\alpha I$ on sphere |
| $\sigma_p, \sigma_m, \sigma_j$ | noise std devs | scalar | prediction / measurement / eigenvalue std |
| $\mathcal{B}_\ell(q)$ | geodesic ball radius $\ell$ | $\subset\mathcal{M}$ | used in irreducibility proof |

Conventions:
- Energies are in reduced/natural units: $\exp(-\mathcal{H})$ is the Boltzmann weight with $\beta=1$ implied (no explicit temperature).
- $\partial\mathcal{H}/\partial p = M(q)^{-1}p = \dot q$; the momentum constraint $C(q)\dot q=0$ equals $C(q)M(q)^{-1}p=0$.
- The $\frac12\log|M(q)|$ term in $U$ is ONLY in the acceptance Hamiltonian and only matters when $M$ depends on $q$ (Riemann-manifold / state-dependent mass); with constant $M$ it is a dropped constant.
- Constrained HMC reduces to standard HMC when $\mathcal{M}=\mathbb{R}^n$ (no constraint, $c\equiv0$, RATTLE reduces to leapfrog).
