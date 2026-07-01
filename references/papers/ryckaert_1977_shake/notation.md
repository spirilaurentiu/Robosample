# Notation - Ryckaert, Ciccotti, Berendsen (1977)

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $N$ | number of point particles (force centers) | int | 256 in butane test (64 mols x 4) |
| $l$ | number of holonomic constraints | int | 5 per butane (3 bonds + 2 angles); $2n-3$ per n-alkane |
| $n$ | chain length (number of CH2/CH3 groups) in Appendix | int | also used as Taylor/integrator order elsewhere |
| $\mathbf{r}_i$ | position of particle i | $\mathbb{R}^3$, Å | Cartesian, single global frame |
| $\dot{\mathbf{r}}_i$ | velocity of particle i | $\mathbb{R}^3$, Å/s | |
| $\ddot{\mathbf{r}}_i$ | acceleration | $\mathbb{R}^3$ | $=\mathbf{F}_i/m_i$ unconstrained |
| $\mathbf{r}_{ij}$ | bond vector $\mathbf{r}_i-\mathbf{r}_j$ | $\mathbb{R}^3$ | direction of constraint force |
| $m_i$ | mass of particle i | mass | uniform $m$ per group in Appendix |
| $V$ | potential energy | energy | $\mathbf{F}_i=-\nabla_i V$ |
| $\mathbf{F}_i$ | physical (potential) force on i | $\mathbb{R}^3$ | |
| $\mathbf{G}_i$ | constraint force on i | $\mathbb{R}^3$ | $=-\sum_k\lambda_k\nabla_i\sigma_k$ |
| $\sigma_k$ | kth constraint function | scalar | $\sigma_k=(\mathbf{r}_i-\mathbf{r}_j)^2-d_{ij}^2$ for bonds; note factor: $\nabla_i\sigma_k=2\mathbf{r}_{ij}$ |
| $\lambda_k$ | kth Lagrange multiplier | scalar(t) | true (analytic) multiplier |
| $\lambda_k^{(s)}$ | sth time derivative of $\lambda_k$ at $t_0$ | scalar | Taylor coefficient |
| $\gamma_k$ | undetermined parameter replacing $\lambda_k^{(n-2)}$ | scalar | solved so constraints hold at $t_0+\Delta t$; $\gamma_k\approx\lambda_k+O(h^2)$ |
| $g_{ij}$ | SHAKE scalar constraint magnitude for pair (i,j) | scalar | $g_{ij}=-2(\Delta t)^2\gamma_k$; symmetric $g_{ij}=g_{ji}$ |
| $d_{ij}$ | fixed bond length of constraint k | Å | |
| $a$ | C-C bond length | 1.53 Å | Appendix |
| $b$ | fixed 1-3 distance encoding bond angle | Å | $b=2a\sin(\theta/2)$ |
| $\theta$ | C-C-C bond angle | 109.28deg (109deg28') | Appendix |
| $h,\Delta t$ | integration time step | s | Verlet; $h=1.95\times10^{-15}$ s (butane) |
| $\mathbf{r}_i'$ | unconstrained (pre-correction) new position | $\mathbb{R}^3$ | from plain Verlet step |
| $\delta\mathbf{r}_i$ | constraint correction displacement | $\mathbb{R}^3$ | $\mathbf{r}_i=\mathbf{r}_i'+\delta\mathbf{r}_i$ |
| $A$ | banded constraint coupling matrix | $(2n-3)\times(2n-3)$ | Appendix eq A.7; elements A.8 |
| $B$ | RHS vector (nonlinear-in-$\gamma$ terms) | $\mathbb{R}^{2n-3}$ | Appendix eq A.7 |
| $m$ (order) | integrator coordinate-error order: error $O(\Delta t^{m+1})$ | int | Verlet: local error $O(h^4)$ so $m=3$ effectively; force derivs to order $n-2$ |

## Conventions

- Single global Cartesian frame; all particles share it (contrast with generalized
  coordinates / Eulerian angles of method 1).
- Constraint sign: $\mathbf{G}_i=-\sum_k\lambda_k\nabla_i\sigma_k$. For a bond
  constraint $\sigma_k=|\mathbf{r}_i-\mathbf{r}_j|^2-d_{ij}^2$, gradient
  $\nabla_i\sigma_k=+2(\mathbf{r}_i-\mathbf{r}_j)=2\mathbf{r}_{ij}$.
- Angle constraints are implemented as fixed 1-3 distances (eq A.2), not as
  explicit angle potentials - this is the RATTLE/SHAKE convention for rigid angles.
- Verlet is not self-starting: needs $\mathbf{r}(0)$ and $\mathbf{r}(h)$; bootstrap
  $\mathbf{r}(h)$ with a second-order Taylor step at $O(h^3)$.
- SHAKE iterates constraints sequentially; the first-order-in-$g$ solve of (5.6)
  (dropping $g^2$) is the standard cheap update, with outer iteration to tolerance.
- Reference bond vector $\mathbf{r}=\mathbf{r}_{ij}(t_0)$ (old positions) sets the
  correction direction; $\mathbf{r}'$ is the current pair vector being corrected.
