# Equations - Pear & Weiner 1979, Brownian dynamics of linked rigid bodies

Reduced units throughout: carbon mass $m$, C-C bond length $l$, and $T_{\mathrm{ref}}=600$ K
all set to unity. $\beta = 1/(k_B T)$. Chain has $N$ bonds, $N+1$ atoms $C_0,\dots,C_N$,
$N-1$ rigid bodies, and generalized coordinates $\phi_0,\dots,\phi_N$ ($\phi_0,\phi_1,\phi_2$
= Bryant orientation angles; $\phi_3,\dots,\phi_N$ = internal dihedrals).

## Langevin dynamics

<!-- eq:1.1 -->
$$ m\ddot{\mathbf{x}} = -\nabla U - \eta \, \dot{\mathbf{x}} + \mathbf{L}(t) $$
- **what:** Langevin equation of motion for a particle of mass $m$ in potential $U$.
- **symbols:** $\mathbf{x}$ - Cartesian position (R^3); $U$ - potential energy; $\eta$ - viscosity (damping coefficient); $\mathbf{L}(t)$ - random fluctuating Langevin force; overdot - time derivative.

<!-- eq:1.2 -->
$$ \langle L_i(t) \rangle = 0 , \qquad \langle L_i(t)L_j(t') \rangle = 2\eta k_B T \, \delta(t - t')\, \delta_{ij} $$
- **what:** first and second moments of the Langevin force (fluctuation-dissipation).
- **symbols:** $L_i$ - Cartesian component $i$ of $\mathbf{L}$; $T$ - temperature; $k_B$ - Boltzmann constant; $\langle\cdot\rangle$ - ensemble average; $\delta$ - Dirac delta / Kronecker delta.

## Wittenburg linked-rigid-body formulation

<!-- eq:2.1 -->
$$ T_{ij} = \begin{cases} -1, & j \ge i \\ 0, & j < i \end{cases}, \quad i, j = 1, \dots, N-1 $$
- **what:** connectivity matrix for the linear chain (tree structure).
- **symbols:** $T_{ij}$ - connectivity matrix element ($-1$ if $C_j$ is on path from $C_i$ to the terminal body). <!-- CHECK: printed as "0, i<i"; corrected to j<i from the tree definition -->

<!-- eq:2.2 -->
$$ \Omega_i = \mathbf{p}_{i+1} \dot{\phi}_{i+1}, \quad i = 2, \dots, N-1 $$
- **what:** relative angular velocity of body $C_i$ w.r.t. body $C_{i-1}$ (single hinge DOF).
- **symbols:** $\Omega_i$ - relative angular velocity (R^3); $\mathbf{p}_{i+1}$ - unit vector along the hinge axis (bond $C_{i-1}C_i$); $\dot\phi_{i+1}$ - dihedral angular rate.

<!-- eq:2.3 -->
$$ \Omega_1 = \mathbf{p}_0 \dot{\phi}_0 + \mathbf{p}_1 \dot{\phi}_1 + \mathbf{p}_2 \dot{\phi}_2 $$
- **what:** angular velocity of the first body from the three Bryant angles.
- **symbols:** $\mathbf{p}_0,\mathbf{p}_1,\mathbf{p}_2$ - unit vectors along the Bryant-angle axes; $\dot\phi_0,\dot\phi_1,\dot\phi_2$ - Bryant angular rates.

<!-- eq:2.4 -->
$$ \omega_i = \sum_{j=1}^i \Omega_j , \quad i = 1, \dots, N-1 $$
- **what:** absolute angular velocity of body $i$ = sum of relative angular velocities up the chain.
- **symbols:** $\omega_i$ - absolute angular velocity of body $C_i$ (R^3).

<!-- eq:2.5 -->
$$ \tau_i = \partial V_{\phi} / \partial \phi_i , \quad i = 3, \dots, N ; \qquad \tau_i = 0,\ i=0,1,2 $$
- **what:** internal hinge torque magnitude driving dihedral $\phi_i$; zero for the orientation angles.
- **symbols:** $\tau_i$ - torque magnitude along $\mathbf{p}_i$ acting on body $C_{i-1}$; $V_\phi$ - rotational (dihedral) potential.

<!-- eq:2.6 -->
$$ \underline{A}\, \ddot{\underline{\phi}} = \underline{B} $$
- **what:** the equations of motion in generalized coordinates (solve for angular accelerations).
- **symbols:** $\underline{A}$ - symmetric positive-definite $(N+1)\times(N+1)$ mass matrix; $\ddot{\underline{\phi}}$ - column of $N+1$ angular accelerations; $\underline{B}$ - column of $N+1$ generalized forces.

<!-- eq:2.7 -->
$$ \underline{A} = (\underline{\mathbf{p}}\, \underline{T}) \cdot \underline{K} \cdot (\underline{\mathbf{p}}\, \underline{T})^T $$
- **what:** generalized mass matrix from hinge-axis projection of the inertia matrix.
- **symbols:** $\underline{\mathbf{p}}$ - $(N+1)\times(N-1)$ matrix of hinge-axis vectors (Eq. 2.9); $\underline{T}$ - connectivity matrix (Eq. 2.1); $\underline{K}$ - matrix of second-rank inertia tensors (Eq. 2.10).

<!-- eq:2.8 -->
$$ \underline{B} = -(\underline{\mathbf{p}}\, \underline{T}) \cdot \{ \underline{K} \cdot (\underline{T}^T \underline{\mathbf{f}}) + \underline{\mathbf{M}}' + \underline{\mathbf{M}}_\eta + \underline{\mathbf{L}}(t) \} - \underline{\tau} $$
- **what:** generalized force vector (kinematic, gyroscopic, damping, Langevin, and hinge-torque terms).
- **symbols:** $\underline{\mathbf{f}}$ - kinematic term (Eq. 2.11); $\underline{\mathbf{M}}'$ - gyroscopic moments (Eq. 2.12); $\underline{\mathbf{M}}_\eta$ - damping moments (Eqs. 2.17-2.20); $\underline{\mathbf{L}}$ - Langevin moments (Eqs. 2.13-2.15); $\underline{\tau}$ - hinge torques (Eq. 2.5).

<!-- eq:2.9 -->
$$ \underline{\mathbf{p}}^{T} = \begin{bmatrix} \mathbf{p}_{0} & \mathbf{p}_{1} & \mathbf{p}_{2} & 0 & \dots & 0 \\ 0 & 0 & 0 & \mathbf{p}_{3} & & \vdots \\ \vdots & & & & \ddots & \vdots \\ 0 & \dots & \dots & & & \mathbf{p}_{N} \end{bmatrix} $$
- **what:** the $(N+1)\times(N-1)$ matrix of hinge-axis vectors (its transpose shown).
- **symbols:** $\mathbf{p}_k$ - unit vector along hinge axis $k$ (R^3 entries).

<!-- eq:2.10 -->
$$ \mathbf{K}_{ij} = \begin{cases} \mathbf{K}_{i}^{*}, & j = i, \\ M\,[\mathbf{b}_{j0}\mathbf{b}_{iN} - (\mathbf{b}_{j0} \cdot \mathbf{b}_{iN})\mathbf{E}], & j > i, \end{cases} $$
- **what:** matrix of inertia tensors; diagonal = augmented-body inertia, off-diagonal = coupling. For $j<i$, $\mathbf{K}_{ij}=\mathbf{K}_{ji}^T$.
- **symbols:** $\mathbf{K}_i^*$ - augmented-body inertia tensor about barycenter; $M$ - total chain mass; $\mathbf{b}_{j0},\mathbf{b}_{iN}$ - body vectors (lower/upper hinge positions rel. barycenter); $\mathbf{E}$ - $3\times3$ identity tensor; $\mathbf{b}_{j0}\mathbf{b}_{iN}$ - dyadic (outer product).

<!-- eq:2.11 -->
$$ f_{1} = \sum_{k=0}^{2} \sum_{l=0}^{2} \frac{\partial \mathbf{p}_{l}}{\partial \phi_{k}} \dot{\phi}_{k} \dot{\phi}_{l} , \qquad f_{i} = \omega_{i-1} \times \Omega_{i} , \quad i = 2, \dots, N-1 $$
- **what:** kinematic (centripetal/Coriolis) term of the equations of motion.
- **symbols:** $\partial\mathbf{p}_l/\partial\phi_k$ - hinge-axis derivatives (Appendix A); $\times$ - cross product. <!-- CHECK: subscript t-1 read as i-1 -->

<!-- eq:2.12 -->
$$ \mathbf{M}_{i}' = -\boldsymbol{\omega}_{i} \times \mathbf{K}_{i}^{*} \cdot \boldsymbol{\omega}_{i} + M \left[ \mathbf{b}_{i0} \times \sum_{j=1}^{i-1} \boldsymbol{\omega}_{j} \times (\boldsymbol{\omega}_{j} \times \mathbf{b}_{jN}) + \mathbf{b}_{iN} \times \sum_{j=i+1}^{N-1} \boldsymbol{\omega}_{j} \times (\boldsymbol{\omega}_{j} \times \mathbf{b}_{j0}) \right] $$
- **what:** gyroscopic moment term $\underline{\mathbf{M}}'$ for body $i$ ($i=1,\dots,N-1$).
- **symbols:** as above; nested cross products give centrifugal contributions.

<!-- eq:2.13 -->
$$ \mathbf{L}_{i} = \mathbf{b}_{i0} \times \sum_{j=0}^{i-1} \mathbf{L}_{C_{j}} + \mathbf{b}_{ii} \times \mathbf{L}_{C_{i}} + \mathbf{b}_{iN} \times \sum_{j=i+1}^{N} \mathbf{L}_{C_{j}}, \quad i = 2, \dots, N-2 $$
- **what:** Langevin moment on interior body $i$ from per-atom Langevin forces.
- **symbols:** $\mathbf{L}_{C_j}$ - Langevin force on atom $C_j$; $\mathbf{b}_{ii}$ - center-of-mass body vector. <!-- CHECK: last term printed without explicit cross; supplied by analogy with (2.17) -->

<!-- eq:2.14 -->
$$ \mathbf{L}_{1} = \mathbf{b}_{11} \times (\mathbf{L}_{C_{0}} + \mathbf{L}_{C_{1}}) + \mathbf{b}_{1N} \times \sum_{j=2}^{N} \mathbf{L}_{C_{j}} + \boldsymbol{\rho}_{0} \times \mathbf{L}_{C_{0}} + \boldsymbol{\rho}_{1} \times \mathbf{L}_{C_{1}} $$
- **what:** Langevin moment on the first (terminal) body, which contains atoms $C_0$ and $C_1$.
- **symbols:** $\boldsymbol{\rho}_0,\boldsymbol{\rho}_1$ - vectors from first-body center of mass to $C_0, C_1$.

<!-- eq:2.15 -->
$$ \mathbf{L}_{N-1} = \mathbf{b}_{N-1, N-1} \times (\mathbf{L}_{C_{N-1}} + \mathbf{L}_{C_{N}}) + \mathbf{b}_{N-1, 0} \times \sum_{j=0}^{N-2} \mathbf{L}_{C_{j}} + \boldsymbol{\rho}_{N-1} \times \mathbf{L}_{C_{N-1}} + \boldsymbol{\rho}_{N} \times \mathbf{L}_{C_{N}} $$
- **what:** Langevin moment on the last (terminal) body, which contains atoms $C_{N-1}$ and $C_N$.
- **symbols:** $\boldsymbol{\rho}_{N-1},\boldsymbol{\rho}_N$ - vectors from last-body center of mass to $C_{N-1}, C_N$.

<!-- eq:2.16 -->
$$ \mathbf{F}_i = -\eta\, \dot{\mathbf{r}}_i $$
- **what:** viscous damping force on atom $C_i$.
- **symbols:** $\mathbf{r}_i$ - position of atom $C_i$ relative to chain center of mass; $\eta$ - viscosity.

<!-- eq:2.17 -->
$$ \mathbf{M}_{\eta_{i}} = \mathbf{b}_{i0} \times \sum_{j=0}^{i-1} \mathbf{F}_{j} + \mathbf{b}_{ii} \times \mathbf{F}_{i} + \mathbf{b}_{iN} \times \sum_{j=i+1}^{N} \mathbf{F}_{j} , \quad i = 2, \dots, N-2 $$
- **what:** damping moment on interior body $i$ from per-atom damping forces.
- **symbols:** $\mathbf{F}_j$ - damping force on atom $C_j$ (Eq. 2.19).

<!-- eq:2.18 -->
$$ \mathbf{M}_{\eta_{1}} = \mathbf{b}_{11} \times (\mathbf{F}_{0} + \mathbf{F}_{1}) + \mathbf{b}_{1N} \times \sum_{j=2}^{N} \mathbf{F}_{j} + \boldsymbol{\rho}_{0} \times \mathbf{F}_{0} + \boldsymbol{\rho}_{1} \times \mathbf{F}_{1} $$
$$ \mathbf{M}_{\eta_{N-1}} = \mathbf{b}_{N-1, N-1} \times (\mathbf{F}_{N-1} + \mathbf{F}_{N}) + \mathbf{b}_{N-1, 0} \times \sum_{j=0}^{N-2} \mathbf{F}_{j} + \boldsymbol{\rho}_{N-1} \times \mathbf{F}_{N-1} + \boldsymbol{\rho}_{N} \times \mathbf{F}_{N} $$
- **what:** damping moments on the two terminal bodies.
- **symbols:** as in Eqs. (2.14)-(2.15) with damping forces $\mathbf{F}$.

<!-- eq:2.19 -->
$$ \mathbf{F}_{i} = \begin{cases} -\eta \left[ \sum_{j=2}^{N-1} \boldsymbol{\omega}_{j} \times \mathbf{b}_{j0} + \boldsymbol{\omega}_{1} \times (\mathbf{b}_{11} + \boldsymbol{\rho}_{i}) \right], & i = 0, 1, \\[4pt] -\eta \left( \sum_{j=1}^{i-1} \boldsymbol{\omega}_{j} \times \mathbf{b}_{jN} + \boldsymbol{\omega}_{i} \times \mathbf{b}_{ii} + \sum_{j=i+1}^{N-1} \boldsymbol{\omega}_{j} \times \mathbf{b}_{j0} \right), & i = 2, \dots, N-2, \\[4pt] -\eta \left[ \sum_{j=1}^{N-2} \boldsymbol{\omega}_{j} \times \mathbf{b}_{jN} + \boldsymbol{\omega}_{N-1} \times (\mathbf{b}_{N-1, N-1} + \boldsymbol{\rho}_{i}) \right], & i = N-1, N. \end{cases} $$
- **what:** per-atom damping force evaluated from body angular velocities and body vectors.
- **symbols:** as above; branch depends on which body the atom belongs to.

<!-- eq:2.20 -->
$$ \mathbf{M}_{\eta_{i}} = -\eta\, \mathbf{K}_{i}^{*} \cdot \boldsymbol{\omega}_{i} + (N+1)\eta \left\{ \mathbf{b}_{i0} \times \sum_{j=1}^{i-1} \boldsymbol{\omega}_{j} \times \mathbf{b}_{jN} + \mathbf{b}_{iN} \times \sum_{j=i+1}^{N-1} \boldsymbol{\omega}_{j} \times \mathbf{b}_{j0} \right\}, \quad i = 1, \dots, N-1 $$
- **what:** damping moment $\underline{\mathbf{M}}_\eta$ reduced form when all masses are equal.
- **symbols:** as above. <!-- CHECK: last term printed without explicit cross; supplied by analogy -->

## Langevin impulse simulation

<!-- eq:2.21 -->
$$ \mathbf{B}(\delta t) = \int_{t}^{t+\delta t} \mathbf{L}(\tau)\, d\tau $$
- **what:** Langevin force impulse over a time step (independent of $t$).
- **symbols:** $\mathbf{B}(\delta t)$ - impulse (R^3); $\delta t$ - time step.

<!-- eq:2.22 -->
$$ \rho[B_j(\delta t)] = \frac{1}{(4\pi \eta k_B T \delta t)^{1/2}} \exp\left[-\frac{|B_j(\delta t)|^2}{4\eta k_B T \delta t}\right] $$
- **what:** Gaussian distribution of each Cartesian impulse component (variance $2\eta k_B T\,\delta t$).
- **symbols:** $B_j$ - Cartesian component $j$ of the impulse.

<!-- eq:2.23 -->
$$ v_i^{(j)} = (v_i^{(j)})^{RK} + \Delta v_i^{(j)}, \quad j = 0, \dots, N, \quad i = 1, 2, 3 $$
- **what:** velocity update (flexible model): Runge-Kutta velocity plus Langevin impulse.
- **symbols:** $v_i^{(j)}$ - component $i$ of velocity of atom $C_j$; superscript RK - Runge-Kutta (no-Langevin) result.

<!-- eq:2.24 -->
$$ \Delta v_i^{(j)} = \frac{1}{m_j} (2\eta k_B T \delta t)^{1/2} \gamma_i^{(j)} $$
- **what:** Langevin velocity increment for atom $C_j$ at temperature $T$.
- **symbols:** $m_j$ - mass of atom $C_j$; $\gamma_i^{(j)}$ - standard normal random variable (mean 0, std 1). <!-- CHECK: printed m_i; should be m_j (mass of atom C_j) -->

<!-- eq:2.25 -->
$$ \dot{\phi}_i = (\dot{\phi}_i)^{RK} + \Delta \dot{\phi}_i , \quad i = 0, \dots, N , \qquad \Delta \dot{\phi}_i = [\underline{A}^{-1}(-\underline{\mathbf{p}}\,\underline{T} \cdot \hat{\underline{\mathbf{L}}})]_i $$
- **what:** generalized-velocity update (rigid model): RK velocity plus transformed Langevin impulse.
- **symbols:** $\hat{\underline{\mathbf{L}}}$ - impulse matrix formed like Eqs. (2.13)-(2.15) with Langevin impulses in place of forces; $\underline{A}^{-1}$ - inverse generalized mass matrix.

## Equilibrium distributions and the metric determinant

<!-- eq:3.1 -->
$$ \rho(\phi_3, l_1, l_2, l_3, \theta_2, \theta_3) = C \exp \left\{ -\frac{\beta}{2} \left[ \sum_{i=1}^{3} k_i (l_i - l)^2 + \sum_{i=2}^{3} k_{\theta} (\theta_i - \theta)^2 \right] \right\} $$
- **what:** flexible-model internal-coordinate equilibrium distribution (three-bond chain).
- **symbols:** $C$ - normalization; $\beta=1/(k_B T)$; $k_i$ - bond spring constants; $k_\theta$ - valence-angle spring constant; $l,\theta$ - equilibrium bond length / valence angle.

<!-- eq:3.2 -->
$$ \rho_{F}(\phi_3) = \frac{1}{2\pi} $$
- **what:** flexible-model dihedral distribution is uniform (even as $k_i,k_\theta\to\infty$).
- **symbols:** $\rho_F$ - flexible dihedral probability density.

<!-- eq:3.3 -->
$$ \rho_{R}(\phi_3) = C_{R}\sqrt{g(\phi_3)} $$
- **what:** rigid-model dihedral distribution is proportional to the square root of the metric determinant.
- **symbols:** $\rho_R$ - rigid dihedral density; $C_R$ - normalization (Eq. 3.4); $g(\phi_3)$ - metric determinant.

<!-- eq:3.4 -->
$$ C_R^{-1} = \int_0^{2\pi} \sqrt{g(\phi_3)} \, d\phi_3 $$
- **what:** normalization constant for the rigid dihedral distribution.
- **symbols:** as above.

<!-- eq:3.5 -->
$$ H_{ij} = \sum_{r=0}^{3} \frac{2}{m_r} \left( \frac{\partial c_i}{\partial \mathbf{x}_r} \cdot \frac{\partial c_j}{\partial \mathbf{x}_r} \right), \quad i, j = 1, \dots, 5 $$
- **what:** Fixman's matrix relating the metric determinant to the constrained coordinates.
- **symbols:** $m_r$ - mass of atom $C_r$; $\mathbf{x}_r$ - Cartesian coordinate of $C_r$; $c_i$ - the 5 constrained coordinates: bond lengths $|\mathbf{l}_i|$ ($i=1,2,3$) and valence-angle proxies $|\mathbf{l}_{i-3}+\mathbf{l}_{i-2}|$ ($i=4,5$). Then $g \propto 1/|H|$ (Fixman theorem).

<!-- eq:3.6 -->
$$ g = C_{g}\,(a + b\cos\phi + c\cos^2\phi + \cos^4\phi) $$
- **what:** closed-form metric determinant for the three-bond 90-degree chain.
- **symbols:** $C_g$ - constant absorbed in normalization; $\phi$ = dihedral $\phi_3$; $a,b,c$ below.

<!-- eq:3.6a -->
$$ a = \frac{1 + 7\alpha + 15\alpha^2 + 10\alpha^3 + 2\alpha^4}{\alpha^4}, \quad b = \frac{2(1 + \alpha)}{\alpha^2}, \quad c = -\frac{3\alpha^2 + 7\alpha + 6}{\alpha^2} $$
- **what:** coefficients of the metric determinant Eq. (3.6).
- **symbols:** $\alpha$ - ratio of end-atom mass to interior-atom mass. <!-- CHECK: c denominator printed alpha^2; consistent with b -->

<!-- eq:3.7 -->
$$ U = k_B T \ln \sqrt{g(\phi_3)} $$
- **what:** Fixman compensating potential added to the internal rotational potential to cancel the metric determinant.
- **symbols:** $U$ - compensating potential energy.

<!-- eq:3.8 -->
$$ \rho_{RF}(\phi_3) = \frac{1}{2\pi} $$
- **what:** rigid model + Fixman potential gives the uniform dihedral distribution.
- **symbols:** $\rho_{RF}$ - compensated rigid dihedral density.

## Barrier crossing rates

<!-- eq:4.1 -->
$$ V(\phi_3) = \begin{cases} \dfrac{1}{2\pi^2} k [(\phi - \phi_b) + \pi]^2, & \phi - \phi_b < -\pi/2, \\[4pt] E_b - \dfrac{k}{2\pi^2} (\phi - \phi_b)^2, & -\pi/2 \le \phi - \phi_b \le \pi/2, \\[4pt] \dfrac{1}{2\pi^2} k [(\phi - \phi_b) - \pi]^2, & \phi - \phi_b > \pi/2, \end{cases} $$
- **what:** single-barrier piecewise-quadratic rotational potential on $[-\pi,\pi]$ (continuous with continuous derivative).
- **symbols:** $k$ - barrier stiffness; $E_b = k/4$ - barrier height; $\phi_b$ - barrier-peak position; $\phi=\phi_3$.

<!-- eq:4.2 -->
$$ f_{TS} = \int_0^\infty \dot{\phi}_3\, d\dot{\phi}_3 \iiint d\dot{\phi}_0\, d\dot{\phi}_1\, d\dot{\phi}_2\; \rho(\phi_b, \dot{\underline{\phi}}) $$
- **what:** exact transition-state barrier-crossing rate (flux of dihedral velocity through the barrier peak).
- **symbols:** $f_{TS}$ - transition-state rate; $\dot{\underline{\phi}}=[\dot\phi_0\,\dot\phi_1\,\dot\phi_2\,\dot\phi_3]$; $\rho$ - phase-space distribution (Eq. 4.3) evaluated at $\phi_3=\phi_b$.

<!-- eq:4.3 -->
$$ \rho(\phi_3, \dot{\underline{\phi}}) = C\, g(\phi_3) \exp\{-\beta[T(\phi_3, \dot{\underline{\phi}}) + V(\phi_3)]\} $$
- **what:** canonical phase-space distribution (orientation angles integrated into $C$).
- **symbols:** $T$ - kinetic energy (Eq. 4.4); $V$ - potential (Eq. 4.1); $C$ - normalization (Eq. 4.7).

<!-- eq:4.4 -->
$$ T(\phi_3, \dot{\underline{\phi}}) = \sum_{i,j=0}^{3} G_{ij}\, \dot{\phi}_i \dot{\phi}_j $$
- **what:** kinetic energy as a quadratic form in the generalized velocities.
- **symbols:** $G_{ij}$ - covariant metric tensor (Eq. 4.5).

<!-- eq:4.5 -->
$$ G_{ij} = \frac{1}{2} \sum_{r=0}^{3} m_r \frac{\partial \mathbf{x}_r}{\partial \phi_i} \cdot \frac{\partial \mathbf{x}_r}{\partial \phi_j} $$
- **what:** covariant metric tensor (mass-weighted Jacobian); $g(\phi_3)=|G_{ij}|$.
- **symbols:** $\partial\mathbf{x}_r/\partial\phi_i$ - derivative of atom position w.r.t. generalized coordinate. <!-- CHECK: second factor printed with phi_i; should be phi_j -->

<!-- eq:4.6 -->
$$ \iiint_{-\infty}^{\infty} d\dot{\underline{\phi}}\; \exp(-\beta G_{ij} \dot{\phi}_{i} \dot{\phi}_{j}) = \frac{\pi^{2}}{\beta^{2} \sqrt{g(\phi_{3})}} $$
- **what:** Gaussian integral over all four generalized velocities.
- **symbols:** $g(\phi_3)=|G_{ij}|$; summation over repeated $i,j$.

<!-- eq:4.7 -->
$$ C = \frac{I_c}{(\pi k_B T)^2} $$
- **what:** normalization constant of the distribution Eq. (4.3).
- **symbols:** $I_c$ - configurational normalization (Eq. 4.8).

<!-- eq:4.8 -->
$$ I_c = \left(\int_0^{2\pi} \sqrt{g(\phi)}\, e^{-\beta V(\phi)}\, d\phi\right)^{-1} $$
- **what:** configurational normalization integral (uncompensated rigid model).
- **symbols:** as above.

<!-- eq:4.9 -->
$$ \int_{0}^{\infty} \dot{\phi}_{3}\, d\dot{\phi}_{3} \iiint_{-\infty}^{\infty} d\dot{\phi}_{0}\, d\dot{\phi}_{1}\, d\dot{\phi}_{2}\; \exp(-\beta G_{ij} \dot{\phi}_{i} \dot{\phi}_{j}) = \frac{1}{2} \left[ \frac{\pi^{3} G^{33}}{\beta^{5}\, g(\phi_{3})} \right]^{1/2} $$
- **what:** flux integral over velocities (one-sided in $\dot\phi_3$).
- **symbols:** $G^{33}$ - (3,3) component of the contravariant metric tensor $G^{ij}=(G_{ij})^{-1}$.

<!-- eq:4.10 -->
$$ f_{TS} = \left(\frac{I_c}{2\sqrt{\pi}}\right) [g(\phi_b)\, G^{33} k_B T]^{1/2}\, e^{-E_b/k_B T} $$
- **what:** transition-state rate, uncompensated rigid model (splits with barrier position via $g(\phi_b)$).
- **symbols:** $g(\phi_b)$ - metric determinant at the barrier; $E_b$ - barrier height.

<!-- eq:4.11 -->
$$ f_{TS} = \frac{I_c}{2\sqrt{\pi}} (k_B T\, G^{33})^{1/2}\, e^{-E_b/k_B T}, \qquad I_c = \left(\int_0^{2\pi} e^{-V(\phi)/k_B T}\, d\phi\right)^{-1} $$
- **what:** transition-state rate with the Fixman potential ($g(\phi_b)$ dropped out; independent of barrier position).
- **symbols:** $I_c$ - configurational normalization without the metric determinant.

## Applied stress

<!-- eq:5.1 -->
$$ V_{\sigma} = V(\phi_3) - \boldsymbol{\sigma} \cdot \mathbf{R} $$
- **what:** total potential (excluding Fixman potential) with end forces $\pm\boldsymbol{\sigma}$ on atoms $C_0,C_3$.
- **symbols:** $\boldsymbol{\sigma}$ - applied force; $\mathbf{R}$ - end-to-end vector of the chain.

<!-- eq:5.2 -->
$$ V_{\sigma}(\phi_3, \Theta) = V(\phi_3) - \sigma R(\phi) \cos \Theta, \qquad R(\phi) = l\,(3 + 2\cos\phi)^{1/2} $$
- **what:** stress potential with $\boldsymbol{\sigma}$ as the reference orientation axis.
- **symbols:** $\Theta$ - angle between $\mathbf{R}$ and $\boldsymbol{\sigma}$; $R(\phi)$ - end-to-end distance of the three-bond 90-degree chain; $\sigma=|\boldsymbol{\sigma}|$.

<!-- eq:5.4 -->
$$ \rho(\Theta,\phi_3,\dot{\underline{\phi}}) = C\sqrt{g(\phi_3)}\, \exp\{-\beta[T(\phi_3,\dot{\underline{\phi}}) + V_{\sigma}(\Theta,\phi_3)]\} $$
- **what:** canonical distribution under applied stress, with Fixman compensation folded in ($\sqrt{g}$ instead of $g$).
- **symbols:** as in Eqs. (4.3)-(4.4) and (5.2).

<!-- eq:5.9 -->
$$ f_{TS}(\sigma) = \frac{I_c(\sigma)}{2\sqrt{\pi}} (k_B T\, G^{33})^{1/2}\, \sinh\!\frac{\sigma R(\phi_b)}{k_B T}\; \frac{e^{-E_b/k_B T}}{R(\phi_b)} $$
- **what:** exact two-dimensional transition-state rate under applied stress (predicts curved Arrhenius plots).
- **symbols:** $I_c(\sigma)$ - stress-dependent normalization (Eq. 5.10); reduces to Eq. (4.11) as $\sigma\to0$.

<!-- eq:5.10 -->
$$ I_c(\sigma) = \left[ \int_0^{2\pi} d\phi\, \frac{1}{R(\phi)} \sinh\!\frac{\sigma R(\phi)}{k_B T}\, e^{-V(\phi)/k_B T} \right]^{-1} $$
- **what:** stress-dependent configurational normalization (evaluated numerically).
- **symbols:** $R(\phi) = l(3+2\cos\phi)^{1/2}$.

<!-- eq:5.11 -->
$$ \mathbf{M}_{\sigma_1} = (\mathbf{b}_{11} - \mathbf{b}_{21}) \times \boldsymbol{\sigma} + \boldsymbol{\rho}_0 \times \boldsymbol{\sigma}, \qquad \mathbf{M}_{\sigma_2} = (\mathbf{b}_{12} - \mathbf{b}_{22}) \times \boldsymbol{\sigma} + \boldsymbol{\rho}_3 \times \boldsymbol{\sigma} $$
- **what:** moments added to $\underline{\mathbf{M}}'$ to include the applied end forces (three-bond chain).
- **symbols:** $\mathbf{b}_{ij}$ - body vectors; $\boldsymbol{\rho}_0,\boldsymbol{\rho}_3$ - vectors from body CoM to the loaded end atoms.

## Appendix A: coordinate transformations and body descriptions

<!-- eq:A1 -->
$$ \mathbf{V}^{(i-1)} = T(\phi_{i+1})\, \mathbf{V}^{(i)} $$
- **what:** transform a vector from body-$C_i$ frame to body-$C_{i-1}$ frame.
- **symbols:** $\mathbf{V}^{(i)}$ - components in frame $i$; $T$ - transformation matrix (Eq. A2).

<!-- eq:A2 -->
$$ T(\phi) = \begin{pmatrix} \cos\theta & -\sin\theta\cos\phi & \sin\theta\sin\phi \\ -\sin\theta & -\cos\theta\cos\phi & \cos\theta\sin\phi \\ 0 & -\sin\phi & -\cos\phi \end{pmatrix} $$
- **what:** inter-body transformation matrix (fixed valence angle $\theta$, variable dihedral $\phi$).
- **symbols:** $\theta$ - constant valence angle; $\phi$ - dihedral angle.

<!-- eq:A3 -->
$$ a^{01} = \begin{pmatrix} c_1 c_2 & -c_1 s_2 & s_1 \\ c_0 s_2 + s_0 s_1 c_2 & c_0 c_2 - s_0 s_1 s_2 & -s_0 c_1 \\ s_0 s_2 - c_0 s_1 c_2 & s_0 c_2 + c_0 s_1 s_2 & c_0 c_1 \end{pmatrix} $$
- **what:** transformation from first-body frame to the reference frame (Bryant angles).
- **symbols:** $c_i=\cos\phi_i$, $s_i=\sin\phi_i$ for $i=0,1,2$.

<!-- eq:A4 -->
$$ \mathbf{b}_{i0} = \frac{l}{2M} [-\beta_i \cos\theta - (N + \alpha - i),\ \ \beta_i \sin\theta,\ \ 0]^T $$
$$ \mathbf{b}_{ii} = \frac{l}{2M} [-\beta_i \cos\theta + (\alpha - 1 + i),\ \ \beta_i \sin\theta,\ \ 0]^T $$
$$ \mathbf{b}_{iN} = \frac{l}{2M} [(\alpha + i)\cos\theta + (\alpha - 1 + i),\ \ -(\alpha + i)\sin\theta,\ \ 0]^T $$
- **what:** augmented body vectors (lower hinge, CoM, upper hinge) for interior bodies, in the body frame ($i=2,\dots,N-2$).
- **symbols:** $M$ - total mass; $\beta_i = N+\alpha-1-i$; $\alpha$ - end/interior mass ratio; $l$ - bond length; $\theta$ - valence angle. <!-- CHECK: OCR row-order ambiguous; components grouped as [x,y,z] with z=0 -->

<!-- eq:A5 -->
$$ \mathbf{K}_{i}^{*} = (\alpha - 1 + i)(b_{i0}^{2}\mathbf{E} - \mathbf{b}_{i0}\mathbf{b}_{i0}^{T}) + (b_{ii}^{2}\mathbf{E} - \mathbf{b}_{ii}\mathbf{b}_{ii}^{T}) + (N + \alpha - i - 1)(b_{iN}^{2}\mathbf{E} - \mathbf{b}_{iN}\mathbf{b}_{iN}^{T}) $$
- **what:** augmented-body inertia tensor for interior bodies (parallel-axis sums over collapsed point masses).
- **symbols:** $b^2 = \mathbf{b}\cdot\mathbf{b}$; $\mathbf{b}\mathbf{b}^T$ - dyadic; $\mathbf{E}$ - identity.

<!-- eq:A6 -->
$$ \mathbf{b}_{10} = 0, \quad \mathbf{b}_{11} = \frac{-(N-2+\alpha)l}{M} \left(\frac{\cos\theta}{2} + \frac{\alpha}{1+\alpha},\ -\frac{\sin\theta}{2},\ 0 \right)^T, \quad \mathbf{b}_{12} = \frac{l}{M} \left( \frac{\alpha + 1}{2} \cos \theta + \alpha,\ -\frac{\alpha + 1}{2} \sin \theta,\ 0 \right)^{T} $$
- **what:** body vectors for the first (terminal) body.
- **symbols:** as above. <!-- CHECK: component grouping reconstructed from OCR; z=0 -->

<!-- eq:A7 -->
$$ \mathbf{K}_{1}^{*} = \begin{pmatrix} a & b & 0 \\ b & c & 0 \\ 0 & 0 & a+c \end{pmatrix} $$
- **what:** augmented inertia tensor of the first body.
- **symbols:** $a,b,c$ from Eq. (A8).

<!-- eq:A8 -->
$$ a = \frac{(N-2+\alpha)(\alpha+1)}{4M} \sin^2 \theta, \quad b = \frac{N-2+\alpha}{2M} \sin\theta \left(\frac{\alpha+1}{2}\cos\theta + \alpha\right), \quad c = \frac{1}{M} \left[ (N-2+\alpha) \cos\theta \left( \frac{\alpha+1}{4} \cos\theta + \alpha \right) + \alpha (N+\alpha-1) \right] $$
- **what:** entries of the first-body inertia tensor Eq. (A7).
- **symbols:** as above.

<!-- eq:A9 -->
$$ \mathbf{R} = \begin{pmatrix} \cos\theta & \sin\theta & 0\\ \sin\theta & -\cos\theta & 0\\ 0 & 0 & -1 \end{pmatrix} $$
- **what:** transformation relating the two terminal bodies: $\mathbf{b}_{N-1,N}=0$, $\mathbf{b}_{NN}=\mathbf{R}\mathbf{b}_{11}$, $\mathbf{b}_{N0}=\mathbf{R}\mathbf{b}_{12}$, $\mathbf{K}_N^*=\mathbf{R}\mathbf{K}_N^*\mathbf{R}^T$.
- **symbols:** $\theta$ - valence angle.

<!-- eq:A10 -->
$$ \mathbf{p}_0 = (\cos\phi_1 \cos\phi_2,\ -\cos\phi_1 \sin\phi_2,\ \sin\phi_1)^T, \quad \mathbf{p}_1 = (\sin\phi_2,\ \cos\phi_2,\ 0)^T, \quad \mathbf{p}_2 = (0,\ 0,\ 1)^T $$
- **what:** hinge (Bryant-axis) vectors in the frame of the first body $C_1$.
- **symbols:** $\phi_1,\phi_2$ - Bryant angles. <!-- CHECK: p0 middle component sign from OCR "-cos phi_1 sin phi_2" -->

## Appendix B: flexible chain formalism

<!-- eq:B1 -->
$$ l_i = |\mathbf{l}_i| = |\mathbf{x}_i - \mathbf{x}_{i-1}|, \quad \theta_i = \cos^{-1}\!\left(\frac{\mathbf{l}_i \cdot \mathbf{l}_{i-1}}{l_i\, l_{i-1}}\right), \quad \phi_{i} = \cos^{-1}\!\left[ -\frac{(\mathbf{l}_{i-1} \times \mathbf{l}_{i-2})}{|\mathbf{l}_{i-1} \times \mathbf{l}_{i-2}|} \cdot \frac{(\mathbf{l}_{i} \times \mathbf{l}_{i-1})}{|\mathbf{l}_{i} \times \mathbf{l}_{i-1}|} \right] $$
- **what:** generalized coordinates (bond length, valence angle, dihedral) from Cartesian atom positions.
- **symbols:** $\mathbf{l}_i = \mathbf{x}_i - \mathbf{x}_{i-1}$ - bond vector; ranges $i=1,\dots,N$ (length), $i=2,\dots,N$ (angle), $i=3,\dots,N$ (dihedral).

<!-- eq:B2 -->
$$ m_{i}\ddot{\mathbf{x}}_{i} = -\nabla_{i} \left[ \sum_{j=1}^{N} V_{l}(l_{j}) + V_{\theta}(\theta_{j}) \right] - \nabla_{i} V_{\phi}(\phi_{3}, \dots, \phi_{N}) - \eta \dot{\mathbf{x}}_{i} + \mathbf{L}_{i}(t) , \quad i = 0, \dots, N $$
- **what:** flexible-model Cartesian Langevin equations of motion.
- **symbols:** $\nabla_i = (\partial/\partial x_1^{(i)}, \partial/\partial x_2^{(i)}, \partial/\partial x_3^{(i)})$; $V_l,V_\theta,V_\phi$ - bond/angle/dihedral potentials.

<!-- eq:B3 -->
$$ V_{l}(l_{j}) = \frac{1}{2} k_{l}(l_{j} - l)^{2}, \quad j = 1, \dots, N $$
- **what:** harmonic bond-stretch potential.
- **symbols:** $k_l$ - bond spring constant; $l$ - equilibrium bond length.

<!-- eq:B4 -->
$$ V_{\theta}(\theta_{j}) = \begin{cases} 0, & j = 1, \\ \frac{1}{2}k_{\theta}(\theta_{j} - \theta)^{2}, & j = 2, \dots, N, \end{cases} $$
- **what:** harmonic valence-angle potential.
- **symbols:** $k_\theta$ - angle spring constant; $\theta$ - equilibrium valence angle.

<!-- eq:B5 -->
$$ V_{\phi}(\phi_3, \ldots, \phi_N) = \sum_{i=3}^{N} V(\phi_i) $$
- **what:** total rotational potential (independent per-dihedral contributions).
- **symbols:** $V(\phi_i)$ - single-dihedral potential (e.g., Eq. 4.1).

<!-- eq:B6 -->
$$ m_{i}\ddot{\mathbf{x}}_{i} = \mathbf{F}_{l_i}^{i} + \mathbf{F}_{l_{i+1}}^{i} + \sum_{k=0}^{2} \mathbf{F}_{\theta_{i+k}}^{i} + \sum_{k=0}^{3} \mathbf{F}_{\phi_{i+k}}^{i} - \eta \dot{\mathbf{x}}_{i} + \mathbf{L}_{i}(t) , \quad i = 0, \dots, N $$
- **what:** flexible equations of motion after differentiating potentials into Cartesian forces.
- **symbols:** $\mathbf{F}_{\alpha_k}^{j}$ - force on atom $C_j$ from generalized coordinate $\alpha_k$ (zero if $\alpha_k$ undefined for that $k$).

<!-- eq:B7 -->
$$ \mathbf{F}_{l_i}^i = -\,k_l\, \mathbf{l}_i\, (1 - l/l_i), \quad i=1,\ldots,N; \qquad \mathbf{F}_{l_{i+1}}^{i} = k_{l}\, \mathbf{l}_{i+1}\, (1 - l/l_{i+1}), \quad i = 1, \dots, N-1 $$
- **what:** bond-stretching forces on atom $C_i$ from its two bonds.
- **symbols:** $\mathbf{l}_i$ - bond vector; $l$ - equilibrium bond length.

<!-- eq:B8 -->
$$ \mathbf{F}_{\theta_{i}}^{i} = -K_{\theta_{i}}(\mathbf{l}_{i-1} - d_{i}^{i-1}\mathbf{e}_{i}), \quad \mathbf{F}_{\theta_{i+1}}^{i} = -K_{\theta_{i+1}}[\mathbf{l}_{i+1} - \mathbf{l}_{i} + d_{i}^{i+1}(\mathbf{e}_{i+1} - \mathbf{e}_{i})], \quad \mathbf{F}_{\theta_{i+2}}^{i} = K_{\theta_{i+2}}(\mathbf{l}_{i+2} - d_{i+1}^{i+2} \mathbf{e}_{i+2}) $$
$$ K_{\theta_j} = \frac{k_{\theta}}{|\mathbf{l}_j \times \mathbf{l}_{j-1}|} (\theta_j - \theta) $$
- **what:** valence-angle forces on atom $C_i$ from the three angles it participates in.
- **symbols:** $\mathbf{e}_j \equiv \mathbf{l}_j/l_j^2$ (Eq. B9); $d_j^k \equiv \mathbf{l}_j\cdot\mathbf{l}_k$ (Eq. B10).

<!-- eq:B9 -->
$$ \mathbf{e}_{j} \equiv \mathbf{l}_{j} / l_{j}^{2}, \qquad d_{j}^{k} = d_{k}^{j} \equiv \mathbf{l}_{j} \cdot \mathbf{l}_{k} $$
- **what:** auxiliary vectors/scalars for the valence-angle and dihedral forces.
- **symbols:** as above.

<!-- eq:B11 -->
$$ a_{j}^{k} = l_{j}^{2} l_{k}^{2} - (d_{j}^{k})^{2} $$
- **what:** auxiliary scalar (squared area of the bond parallelogram) for dihedral forces.
- **symbols:** $l_j^2 = \mathbf{l}_j\cdot\mathbf{l}_j$; $d_j^k = \mathbf{l}_j\cdot\mathbf{l}_k$.

<!-- eq:B13 -->
$$ \mathbf{F}_{\phi_{i}}^{i} = K_{\phi_{i}}\,^{i}\mathbf{A}_{i-2}^{i-1}, \quad \mathbf{F}_{\phi_{i+1}}^{i} = K_{\phi_{i+1}}(^{i}\mathbf{A}_{i+1}^{i-1} - ^{i+1}\mathbf{A}_{i-1}^{i} + ^{i}\mathbf{A}_{i-1}^{i+1}), $$
$$ \mathbf{F}_{\phi_{i+2}}^{i} = -K_{\phi_{i+2}}(^{i+1}\mathbf{A}_{i}^{i+2} - ^{i}\mathbf{A}_{i+2}^{i+1} + ^{i+1}\mathbf{A}_{i+2}^{i}), \quad \mathbf{F}_{\phi_{i+3}}^{i} = -K_{\phi_{i+3}}\,^{i+1}\mathbf{A}_{i+3}^{i+2} $$
- **what:** dihedral (rotational-potential) forces on atom $C_i$ from the (up to four) dihedrals it enters.
- **symbols:** $^i\mathbf{A}_b^j$ - auxiliary vector (Eq. B12); $K_{\phi_k}$ - dihedral force coefficient $\partial V/\partial\phi_k$ scaled by geometry. <!-- CHECK: leading superscripts/subscripts heavily OCR-garbled in B12/B13; indices reconstructed from the stated symmetry pattern -->
