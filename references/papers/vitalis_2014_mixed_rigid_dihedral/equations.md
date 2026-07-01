# Equations - Vitalis & Pappu 2014, mixed rigid-body/dihedral integrator

Coordinate model: each molecule is a tree of rigid bodies connected by dihedral
"hinges" plus 6 rigid-body DOF (3 translation, 3 rotation) for the base. The full
generalized coordinate vector is `phi`, dimension K <= N_at. `omega` are the
generalized velocities (angular velocities for hinges/rotations).

The central algorithmic idea: **replace the true (non-diagonal) mass-metric tensor
G = J^T M J with its diagonal `I_D`**, decoupling every DOF, and integrate each DOF
independently by imposing per-DOF energy conservation. This yields "MMT-artifact-free"
configurational sampling by construction, at the cost of correct dynamics.

---

<!-- eq:1 -->
$$ \dot{\mathbf{p}} = -\nabla U(\mathbf{r}), \qquad \dot{\mathbf{r}} = \mathbf{M}^{-1}\mathbf{p} $$
- **what:** Cartesian Hamilton equations of motion (reference/baseline). M diagonal, KE = (1/2) p^T M^{-1} p.
- **symbols:** p - Cartesian momenta (R^{3N_at}); r - Cartesian positions; U - potential energy; M - diagonal mass matrix.

<!-- eq:2 -->
$$ \dot{\mathbf{p}}_\phi = \frac{d}{dt}\left[(\mathbf{J}^{\mathrm{T}}\mathbf{M}\mathbf{J})\boldsymbol{\omega}\right] = -\mathbf{J}\nabla U(\mathbf{r}) = \mathbf{F}_\phi $$
- **what:** exact generalized-coordinate EOM. The d/dt of (J^T M J) produces bias torques (via dJ/dt) not contained in F_phi; the paper's approximation is designed to hide these.
- **symbols:** p_phi - generalized momenta conjugate to omega; J - Jacobian dr/dphi (3N_at x K); G = J^T M J - mass-metric tensor (MMT); F_phi - generalized (projected) force; omega - generalized velocities.

<!-- eq:3 -->
$$ (J_{3i-2,k}\quad J_{3i-1,k}\quad J_{3i,k}) = \frac{\partial \vec{r}_i}{\partial \phi_k} = \vec{a}_k \times (\vec{r}_i - \vec{b}_k) $$
- **what:** Jacobian block for angular DOF k acting on atom i (rotation about an axis). Same form as eq:A2 for Y.
- **symbols:** a_k - unit rotation axis of DOF k; b_k - reference point on the axis; r_i - position of atom i.

<!-- eq:5 -->
$$ E'_k = \frac{1}{2}\sum_{i}^{K} I_{ii}\,\omega_i^2 = \frac{1}{2}\boldsymbol{\omega}^{\mathsf{T}}\mathbf{I}_D\boldsymbol{\omega} \;\neq\; \frac{1}{2}\mathbf{p}^{\mathsf{T}}\mathbf{M}^{-1}\mathbf{p} $$
- **what:** KEY APPROXIMATION - modified kinetic energy using a diagonal effective-mass matrix I_D (diagonal of the MMT), removing kinetic coupling. Not equal to true Cartesian KE.
- **symbols:** E'_k - modified kinetic energy; I_D - diagonal effective-mass matrix (elements I_kk = I_ii); omega - generalized velocities; K - number of free DOF.

<!-- eq:6 -->
$$ \operatorname{diag}(\mathbf{I}_D) = \begin{bmatrix} M^1 & M^1 & M^1 & I_x^1 & I_y^1 & I_z^1 & I_{\phi 1}^1 & \dots & I_{\phi k}^1 & M^2 & \dots \end{bmatrix} $$
- **what:** layout of the diagonal effective-mass vector, per molecule: 3 translational masses (=molecular mass M), 3 rotational inertias about lab axes through the COM, then one effective inertia per dihedral hinge; then next molecule.
- **symbols:** superscript = molecule index; M - molecular mass; I_{x/y/z} - rotational inertia about lab axes through COM; I_{phi i} - effective inertia of dihedral i.

<!-- eq:7 -->
$$ 0 = \boldsymbol{\omega}(t_2)^{\mathsf{T}}\mathbf{I}_D(t_2)\boldsymbol{\omega}(t_2) - \boldsymbol{\omega}(t_1)^{\mathsf{T}}\mathbf{I}_D(t_1)\boldsymbol{\omega}(t_1) - \delta t\,(\boldsymbol{\omega}(t_1)+\boldsymbol{\omega}(t_2))\mathbf{F}_\phi(t_{1.5}) $$
- **what:** discrete energy-conservation condition over a step t1 -> t2 (KE change balances work by the midpoint force). Defines the constant-energy (NVE) integrator.
- **symbols:** t_1, t_2 - step endpoints; t_{1.5} - midpoint at which force is evaluated; delta_t - time step; F_phi(t_{1.5}) - midpoint generalized force.

<!-- eq:8 -->
$$ \omega_k(t_2) = \delta t\,\frac{F_{\phi,k}(t_{1.5})}{2 I_{kk}(t_2)} \pm \frac{\sqrt{F_{\phi,k}(t_{1.5})^2\,\delta t^2 + 4 I_{kk}(t_2)\,\omega_k(t_1)\,F_{\phi,k}(t_{1.5})\,\delta t + 4 I_{kk}(t_2)\,I_{kk}(t_1)\,\omega_k(t_1)^2}}{2 I_{kk}(t_2)} $$
- **what:** CORE per-DOF velocity update. Solving eq:7 independently for each DOF gives a quadratic in omega_k(t2); this is its exact solution. Implicit (needs I_kk(t2)); time-reversible; may lack a real root (fixed by eq:9). Reduces to leapfrog if masses are constant.
- **symbols:** omega_k - generalized velocity of DOF k; F_{phi,k} - projected force on DOF k (eq:B1); I_kk(t1),I_kk(t2) - effective mass at step endpoints; +/- - branch chosen using eq:9 as a guide.

<!-- eq:9 -->
$$ \omega_k(t_2) \approx \delta t\,\frac{F_{\phi,k}(t_{1.5})}{I_{kk}(t_2)} + \sqrt{\frac{I_{kk}(t_1)}{I_{kk}(t_2)}}\;\omega_k(t_1) $$
- **what:** approximate, single-valued velocity update (approximates I_kk(t2) in the cross term by sqrt(I_kk(t1) I_kk(t2))). NOT time-reversible; used to pick the correct root/sign of eq:8 and eq:11. Explicit conservation of omega_k * I_kk^{1/2} is visible here.
- **symbols:** as eq:8.

<!-- eq:10 -->
$$ \phi_k(t_{2.5}) = \phi_k(t_{1.5}) + \delta t\,\omega_k(t_2) $$
- **what:** position (dihedral) update - simple leapfrog-style increment. Exception: rigid-body rotation uses the quaternion eq:12 instead of an explicit phi.
- **symbols:** phi_k - generalized coordinate; t_{2.5} - half-step ahead of t_2.

<!-- eq:11 -->
$$ \omega_k(t_i) = \tau_\Lambda\,\frac{F_{\phi,k}(t_{1.5})}{2 I_{kk}(t_i)} \pm \frac{\sqrt{F_{\phi,k}(t_{1.5})^2\,\tau_\Lambda^2 + 4 I_{kk}(t_i)\,\omega_k(t_{i-\Lambda^{-1}})\,F_{\phi,k}(t_{1.5})\,\tau_\Lambda + 4 I_{kk}(t_i)\,I_{kk}(t_{i-\Lambda^{-1}})\,\omega_k(t_{i-\Lambda^{-1}})^2}}{2 I_{kk}(t_i)} $$
with $\tau_\Lambda = \delta t/\Lambda$ and $i = 1+\Lambda^{-1},\,1+2\Lambda^{-1},\,\dots,\,2$.
- **what:** CORE PRODUCTION INTEGRATOR (all simulations use this). Sub-steps eq:8 into Lambda segments assuming linear evolution of I_kk between known values at t_1, t_{1.5}, t_2. Force held constant (force-explicit). Not time-reversible; eq:9-analog picks the root each sub-step.
- **symbols:** Lambda - number of velocity sub-steps (multiple of 2; benefits taper for large Lambda; typical Lambda=4); tau_Lambda - sub-step; I_kk(t_i) linearly interpolated between t_1, t_{1.5}, t_2.

<!-- eq:12 -->
$$ q_{rot}(t_{1.5}\to t_{2.5}) = \left[\,c\quad \sin\!\left(\tfrac{1}{2}\delta t\,\omega_x\right)\quad \sin\!\left(\tfrac{1}{2}\delta t\,\omega_y\right)\quad \sin\!\left(\tfrac{1}{2}\delta t\,\omega_z\right)\right] $$
- **what:** rigid-body rotation update as a quaternion built directly from the three rotational angular velocities at t_2 (about lab axes through the molecular COM). c set by unit-norm constraint.
- **symbols:** q_rot - unit rotation quaternion [w, x, y, z]; omega_{x/y/z} - rigid-body rotational velocities at t_2 about lab axes through COM; c = sqrt(1 - sum of the three sin^2 terms) (unit-length constraint). <!-- CHECK: paper states "c is determined by the constraint that the quaternion be of unit length"; explicit c=sqrt(1-...) inferred. -->

<!-- eq:13 -->
$$ \forall k:\quad \langle \omega_k^2\, I_{kk}\rangle = \beta^{-1} $$
- **what:** equipartition target enforced per DOF by the thermostats (each DOF carries kT of the modified KE). Must hold approximately for the artifact-free argument to be valid.
- **symbols:** beta = 1/(k_B T); <.> - ensemble average.

<!-- eq:17 -->
$$ I_{kk}(t)^{1/2}\,\dot{\Omega}_k(t) = -\frac{\partial U(\phi(t))}{\partial \phi_k}, \qquad I_{kk}(t)^{1/2}\,\dot{\phi}_k(t) = \Omega_k(t) $$
- **what:** underlying continuous EOM of the scheme in the mass-weighted velocity Omega_k = omega_k I_kk^{1/2}. The I_kk^{1/2} factors are what preserve volume in phi (hence no MMT artifacts). Generally artificial dynamics (Phi_k substitution fails for non-constant mass).
- **symbols:** Omega_k = omega_k * I_kk(t)^{1/2} - mass-weighted generalized velocity (dynamical variable); phi_k - generalized coordinate.

<!-- eq:20 -->
$$ Q_S = C(T)\int_{\phi_{\min}}^{\phi_{\max}} \exp\left[-\beta\,U(\phi(\mathbf{r}))\right]\,d\phi $$
- **what:** configurational partition function of the scheme - factorizes into a temperature-only prefactor C(T) times the pure Boltzmann configurational integral, i.e. NO det(G) / Fixman weighting. This is the "artifact-free" statement.
- **symbols:** Q_S - subsystem partition function; C(T) - thermal prefactor (independent of phi); phi_min/phi_max - integration bounds per coordinate.

<!-- eq:A1 -->
$$ \vec{r}_i = \vec{r}_j + |\vec{r}_i - \vec{r}_j|\cdot\left[(\vec{c}_1\times\vec{a}_k)\sin\alpha\cos\phi + \vec{c}_1\sin\alpha\sin\phi - \vec{a}_k\cos\alpha\right] $$
$$ \vec{c}_1 = \frac{(\vec{r}_l - \vec{r}_m)\times(\vec{r}_j - \vec{r}_l)}{|(\vec{r}_l - \vec{r}_m)\times(\vec{r}_j - \vec{r}_l)|}, \qquad \vec{a}_k = \frac{\vec{r}_j - \vec{r}_l}{|\vec{r}_j - \vec{r}_l|} $$
- **what:** Z-matrix -> Cartesian placement (NeRF-style): place atom i from bond length |r_i-r_j|, bond angle alpha (i-j-l), dihedral phi (i-j-l-m) using three previously-built reference atoms j,l,m. Applied hierarchically = operator A^{-1}. Requires j,l,m non-collinear.
- **symbols:** r_j,r_l,r_m - three reference atoms already placed; alpha - bond angle at r_j; phi - dihedral; a_k - unit bond axis (r_j->r_l direction); c_1 - unit normal defining the perpendicular direction.

<!-- eq:A2 -->
$$ (Y_{3i-2,k}\quad Y_{3i-1,k}\quad Y_{3i,k}) = \frac{\partial \vec{r}_i}{\partial \phi_k} = \vec{a}_k\times(\vec{r}_i - \vec{b}_k) $$
- **what:** reduced covariant-basis matrix Y (3N_at x K) for angular DOF; identical form to J block eq:3 but restricted to free DOF. Sign of a_k and which terms are zero change with choice of base.
- **symbols:** Y - reduced Jacobian (free DOF only); a_k,b_k as eq:3.

<!-- eq:A3 -->
$$ \frac{d\vec{r}_i}{dt} = \sum_k \frac{\partial \vec{r}_i}{\partial \phi_k}\,\omega_k $$
- **what:** instantaneous Cartesian velocity of atom i; sum over DOF toward the base in the same branch (incl. parent branches and rigid-body motion). Used to compute the TRUE KE (1/2)p^T M^{-1} p for diagnostics.
- **symbols:** as eq:A2; sum restricted to base-ward DOF of atom i.

<!-- eq:B1 -->
$$ F_{\phi,k} = \vec{a}_k\cdot\sum_i \vec{r}_i\times\vec{F}_{r,i} - \vec{a}_k\cdot\left(\vec{b}_k\times\sum_i \vec{F}_{r,i}\right) $$
- **what:** CORE recursion - projected generalized force on hinge k = axis-projected net torque about b_k. Inward (tip->base) recursion; sums accumulate over all atoms tip-ward of k (incl. sub-branches); branch merges combine partial sums. O(N_at) because sums have no k-specific terms. Rigid rotation is the last recursion step (sum over all atoms).
- **symbols:** F_{r,i} - Cartesian force on atom i (= -grad_i U); a_k - hinge axis; b_k - point on axis; sum_i over tip-ward atoms.

<!-- eq:B2 -->
$$ I_{kk} = \sum_i m_i \vec{r}_i^{\,2} + \left[\vec{b}_k^{\,2} - (\vec{a}_k\cdot\vec{b}_k)^2\right]\sum_i m_i + \left[2(\vec{a}_k\cdot\vec{b}_k)\vec{a}_k - 2\vec{b}_k\right]\cdot\sum_i m_i \vec{r}_i - \left((\vec{a}_k\otimes\vec{a}_k),\;\sum_i m_i(\vec{r}_i\otimes\vec{r}_i)\right) $$
- **what:** CORE recursion - effective inertia (diagonal MMT element) of hinge k. Last term is the Frobenius inner product of outer-product matrices. Same inward recursion / accumulated sums as eq:B1 (O(N_at)); sums have no k-specific terms.
- **symbols:** m_i - atomic mass; r_i - atomic position; a_k - unit axis; b_k - point on axis; (X,Y) - Frobenius inner product; (x ⊗ x) - outer product (3x3).

<!-- eq:B3 -->
$$ \frac{d\vec{r}_i}{dt} = -\vec{r}_i\times\sum_k \omega_k\vec{a}_k + \sum_k \vec{b}_k\times\omega_k\vec{a}_k $$
- **what:** CORE recursion - Cartesian atomic velocities via OUTWARD (base->tip) recursion; sums over DOF toward the base (see eq:A3). O(N_at); no atom-specific terms in the sums.
- **symbols:** as eq:A3/B1; sum_k over base-ward DOF of atom i.

<!-- eq:21 -->
$$ \frac{K}{2}k_B \approx \frac{\langle(\boldsymbol{\omega}^{\mathsf{T}}\mathbf{I}_D\boldsymbol{\omega})^2\rangle - \langle\boldsymbol{\omega}^{\mathsf{T}}\mathbf{I}_D\boldsymbol{\omega}\rangle^2}{4 k_B\langle T\rangle^2} < \frac{\langle(\mathbf{p}^{\mathsf{T}}\mathbf{M}^{-1}\mathbf{p})^2\rangle - \langle\mathbf{p}^{\mathsf{T}}\mathbf{M}^{-1}\mathbf{p}\rangle^2}{4 k_B\langle T\rangle^2} $$
- **what:** ideal (kinetic) heat-capacity relation. Only the modified KE (omega^T I_D omega) reproduces the ideal K/2 k_B; the true-Cartesian-KE fluctuations are strictly larger (inequality). Diagnostic/consistency check.
- **symbols:** K - number of constrained-system DOF; T - temperature; other terms as above.

<!-- eq:22 -->
$$ \left\langle \sum_i^K \sum_{j\neq i}^K G_{S,ij}\,\omega_i\,\omega_j \right\rangle \approx 0 $$
- **what:** condition under which <omega^T I_D omega> ~ <p^T M^{-1} p>: the ensemble-averaged off-diagonal MMT contribution vanishes (uncoupled branches, weak long-time velocity cross-correlations, or cancellation).
- **symbols:** G_S - MMT of the free subsystem; G_{S,ij} - off-diagonal elements.

---

## Derivations (not implemented)
Eqs. 14-16, 18-20 are the derivation of the underlying continuous EOM (eq:17) and
partition function (eq:20) from the discrete scheme (eq:8): substituting
Omega_k = omega_k I_kk^{1/2} into eq:8 (eq:14), rewriting the finite differences
(eq:15), taking delta_t -> 0 (eq:16), recognizing the Lagrangian
L = sum_k (1/2)Omega_k^2 - U (eq:18), and integrating the momenta in Q_S
(eq:19 -> eq:20). These justify the artifact-free property; they are not part of
the integrator loop. Implement eq:8/eq:11 (velocities), eq:10/eq:12 (positions),
eq:B1/eq:B2/eq:B3 (recursions), eq:A1 (Z-matrix build).
