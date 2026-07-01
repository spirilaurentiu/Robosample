# Equations — Mazur 1991, explicit internal-coordinate MD

Implementable equations. Intermediate derivation results (eqs. 6-12, 15) are kept
in `paper.md` under `Derivation (not implemented)`; the final assembled equations
of motion are eq:16 and eq:18.

<!-- eq:1 -->
$$\frac{d}{dt}\left(\frac{\partial L}{\partial \dot{\theta}_k}\right) - \frac{\partial L}{\partial \theta_k} = 0, \qquad k = 1, 2, \dots, n$$
- **what:** Lagrange's equations of motion for each generalized (internal) coordinate; L = T − U.
- **symbols:** L - Lagrangian (kinetic minus potential energy); θ_k - k-th generalized internal coordinate; θ̇_k - its time derivative (generalized velocity); n - number of free generalized coordinates.

<!-- eq:2 -->
$$\sum_{\alpha \in d_k} m_{\alpha} \left[ \frac{d}{dt} \left( \dot{\mathbf{r}}_{\alpha} \frac{\partial \dot{\mathbf{r}}_{\alpha}}{\partial \dot{\theta}_{k}} \right) - \dot{\mathbf{r}}_{\alpha} \left( \frac{\partial \dot{\mathbf{r}}_{\alpha}}{\partial \theta_{k}} \right) \right] = -\frac{\partial U}{\partial \theta_{k}}$$
- **what:** Lagrange equation after substituting L = Σ ½ m_α ṙ_α² − U(θ); LHS is the inertial term (mass matrix + velocity terms), RHS is the negative conformational-energy gradient.
- **symbols:** d_k - set of atoms whose position depends on θ_k; m_α - mass of atom α (amu); **r**_α - Cartesian position of atom α (Å); ṙ_α - its velocity; U - potential energy (kcal/mol); ∂U/∂θ_k - generalized force on coordinate k.

<!-- eq:3 -->
$$d\mathbf{r}_{\alpha} = \sum_{i=1}^{n_{\alpha}} \left[ S_{i}\, \mathbf{e}_{i}^{\alpha} \times (\mathbf{r}_{\alpha} - \mathbf{r}_{i}^{\theta})\, d\theta_{i}^{\alpha} + (1 - S_{i})\, \mathbf{e}_{i}^{\alpha}\, d\theta_{i}^{\alpha} \right]$$
- **what:** infinitesimal Cartesian displacement of atom α as a sum over the chain of internal coordinates that affect it; rotation term (S=1, angle) + translation term (S=0, bond length).
- **symbols:** n_α - number of variables in chain V_α determining r_α; S_i - variable-type indicator (1 for angle, 0 for bond length); **e**_i^α - unit vector of variable i (rotation axis for angles, translation direction for bonds); **r**_i^θ - radius vector of the node of variable i; θ_i^α - value of variable i in chain V_α.

<!-- eq:4 -->
$$\dot{\mathbf{r}}_{\alpha} = \sum_{i=1}^{n_{\alpha}} \left[ S_{i}\, \mathbf{e}_{i} \times (\mathbf{r}_{\alpha} - \mathbf{r}_{i}^{\theta}) \cdot \dot{\theta}_{i} + (1 - S_{i})\, \mathbf{e}_{i} \cdot \dot{\theta}_{i} \right]$$
- **what:** Cartesian velocity of atom α (time derivative of eq:3); upper index α on variables dropped for brevity.
- **symbols:** θ̇_i - generalized velocity of variable i; other symbols as in eq:3. Note **r**_α − **r**_i^θ = **r**_{α/i} (see eq:5).

<!-- eq:5 -->
$$\frac{\partial \dot{\mathbf{r}}_{\alpha}}{\partial \dot{\theta}_{k}} = S_{k}\, \mathbf{e}_{k} \times \mathbf{r}_{\alpha/k} + (1 - S_{k})\, \mathbf{e}_{k}$$
- **what:** partial derivative of atom velocity w.r.t. generalized velocity θ̇_k; this is the Jacobian column (per-atom) used to reduce Cartesian forces to generalized forces and to build the mass matrix.
- **symbols:** **r**_{α/k} = **r**_α − **r**_k^θ - position of atom α relative to node k (Å); S_k, **e**_k as before.

<!-- eq:13 -->
$$\frac{\partial \mathbf{e}_m}{\partial \theta_k} = \begin{cases} S_k(\mathbf{e}_k \times \mathbf{e}_m) & \text{for } m > k \\ 0 & \text{for } m \leq k \end{cases}$$
- **what:** derivative of a unit vector w.r.t. an earlier coordinate; only lower-index (earlier) angle variables rotate later unit vectors.
- **symbols:** **e**_m - unit vector of variable m; **e**_k - unit vector (rotation axis) of variable k; S_k - indicator.

<!-- eq:14 -->
$$\frac{\partial \mathbf{r}_{\alpha/m}}{\partial \theta_k} = \begin{cases} S_k(\mathbf{e}_k \times \mathbf{r}_{\alpha/m}) & \text{for } m > k \\ S_k(\mathbf{e}_k \times \mathbf{r}_{\alpha/k}) + (1 - S_k)\, \mathbf{e}_k & \text{for } m \leq k \end{cases}$$
- **what:** derivative of relative atom position **r**_{α/m} w.r.t. coordinate θ_k; drives the velocity-dependent (Coriolis/centrifugal) terms.
- **symbols:** **r**_{α/m} = **r**_α − **r**_m^θ; **r**_{α/k} = **r**_α − **r**_k^θ; **e**_k rotation axis of k; S_k indicator.

<!-- eq:16 -->
$$\begin{aligned}
\sum_{\alpha \in d_{k}} m_{\alpha} \Big\{ &\sum_{i=1}^{n_{\alpha}} \left[ S_{k} S_{i}(\mathbf{e}_{k} \times \mathbf{r}_{\alpha/k})(\mathbf{e}_{i} \times \mathbf{r}_{\alpha/i}) + S_{k} (1 - S_{i})(\mathbf{e}_{k} \times \mathbf{r}_{\alpha/k}) \mathbf{e}_{i} + (1 - S_{k}) S_{i}\, \mathbf{e}_{k}(\mathbf{e}_{i} \times \mathbf{r}_{\alpha/i}) + (1 - S_{k})(1 - S_{i})\, \mathbf{e}_{k} \mathbf{e}_{i} \right] \cdot \ddot{\theta}_{i} \\
&+ \sum_{i=1}^{n_{\alpha}} \left[ S_{k} S_{i}^{2}(\mathbf{e}_{k} \times \mathbf{r}_{\alpha/k})(\mathbf{e}_{i} \times \mathbf{e}_{i} \times \mathbf{r}_{\alpha/i}) + (1 - S_{k}) S_{i}^{2}\, \mathbf{e}_{k}(\mathbf{e}_{i} \times \mathbf{e}_{i} \times \mathbf{r}_{\alpha/i}) \right] \cdot \dot{\theta}_{i}^{2} \\
&+ 2 \sum_{i=2}^{n_{\alpha}} \sum_{m=1}^{i-1} \big[ S_{k} S_{i} S_{m}(\mathbf{e}_{k} \times \mathbf{r}_{\alpha/k})(\mathbf{e}_{m} \times \mathbf{e}_{i} \times \mathbf{r}_{\alpha/i}) + S_{k} (1 - S_{i}) S_{m}(\mathbf{e}_{k} \times \mathbf{r}_{\alpha/k})(\mathbf{e}_{m} \times \mathbf{e}_{i}) \\
&\qquad + (1 - S_{k}) S_{i} S_{m}\, \mathbf{e}_{k}(\mathbf{e}_{m} \times \mathbf{e}_{i} \times \mathbf{r}_{\alpha/i}) + (1 - S_{k})(1 - S_{i}) S_{m}\, \mathbf{e}_{k}(\mathbf{e}_{m} \times \mathbf{e}_{i}) \big] \cdot \dot{\theta}_{m} \dot{\theta}_{i} \Big\} = -\frac{\partial U}{\partial \theta_{k}}
\end{aligned}$$
- **what:** THE explicit equation of motion for internal coordinate θ_k (local per-α numbering inside braces). Coefficient of θ̈_i is the mass-matrix element a_ki; the θ̇_i² and θ̇_m θ̇_i terms are the velocity-dependent (b, c) coefficients; RHS is the generalized force. Products like (**a**)(**b**) are dot products; **a** × **b** × **c** means **a** × (**b** × **c**).
- **symbols:** m_α - atom mass; S_k,S_i,S_m - variable-type indicators (1 angle / 0 bond); **e**_k,**e**_i,**e**_m - unit vectors; **r**_{α/k},**r**_{α/i} - relative atom positions; θ̈_i - generalized acceleration; θ̇_i,θ̇_m - generalized velocities; d_k - atoms depending on θ_k; U - potential.

<!-- eq:17 -->
$$\sum_{i=1}^{n_\alpha-1} \sum_{m=i+1}^{n_\alpha} F_{im} = \sum_{i=2}^{n_\alpha} \sum_{m=1}^{i-1} F_{mi}$$
- **what:** summation-order swap identity used to simplify the double sums into eq:16 (upper vs lower triangular reindexing).
- **symbols:** F_{im} - arbitrary term indexed by ordered pair (i,m); n_α - chain length.

<!-- eq:18 -->
$$\sum_{i} a_{ki}\, \ddot{\theta}_{i} = -\frac{\partial U}{\partial \theta_{k}} - \sum_{i} b_{ki}\, \dot{\theta}_{i}^{2} - \sum_{i} \sum_{m} c_{kim}\, \dot{\theta}_{m} \dot{\theta}_{i}$$
- **what:** compact global-numbering form of the equations of motion; a linear system A θ̈ = force solved for the generalized accelerations each step (A is the internal-coordinate mass matrix). Coefficients a_ki, b_ki, c_kim assembled by summation over the tree topology from eq:16.
- **symbols:** a_ki - mass-matrix element (symmetric, coefficient of θ̈_i in equation k); b_ki - coefficient of θ̇_i² (centrifugal); c_kim - coefficient of the cross velocity product θ̇_m θ̇_i (Coriolis); indices run 1..n_var (global). Solve for θ̈ via Cholesky since A is symmetric positive-definite.

<!-- eq:temperature -->
$$T = \frac{2 \langle K \rangle}{N\, k_{\rm B}}$$
- **what:** instantaneous/average temperature from kinetic energy using N internal DOF (equipartition); note N is the number of internal degrees of freedom, NOT 3× atoms.
- **symbols:** K - instantaneous kinetic energy; ⟨K⟩ - its average; N - number of internal degrees of freedom; k_B - Boltzmann constant.

<!-- eq:deltaE -->
$$\delta_{E} = \frac{\sqrt{\langle \Delta E^{2} \rangle}}{\langle E \rangle}$$
- **what:** relative RMS energy-conservation error used to compare models/time steps (dimensionless).
- **symbols:** ΔE - deviation of total energy from its mean; ⟨ΔE²⟩ - mean-square deviation; ⟨E⟩ - average total energy over the trajectory.
