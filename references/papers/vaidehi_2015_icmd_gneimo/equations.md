# Equations - Vaidehi & Jain 2015, ICMD / GNEIMO

<!-- eq:1 -->
$$ \mathcal{M}(\theta)\ddot{\theta} + C(\theta, \dot{\theta}) = \mathcal{T}(\theta) $$
- **what:** Coupled internal-coordinate (torsional MD) equations of motion when bond lengths and bond angles are rigid; solve for $\ddot\theta$ then integrate.
- **symbols:** $\theta$ - generalized internal coordinates, e.g. torsion angles (vector, R^N); $\dot\theta,\ddot\theta$ - generalized velocity/acceleration; $\mathcal{M}(\theta)$ - configuration-dependent mass matrix / moment-of-inertia tensor (N x N, symmetric PD); $C(\theta,\dot\theta)$ - velocity-dependent Coriolis + centrifugal force vector (R^N); $\mathcal{T}(\theta)$ - generalized forces / torques (R^N).
<!-- Note: raw OCR had a stray extra "(" as "M(θ)(θ̈ + C"; corrected to M θ̈ + C = T per surrounding prose. -->

<!-- eq:2 -->
$$ \ddot{\theta} = [I - \mathcal{H}\psi\mathcal{K}]\, \mathcal{D}^{-1} [\mathcal{T} - \mathcal{H}\psi(\mathcal{K}\mathcal{T} + \mathcal{P}a + b)] - \mathcal{K}^{*}\psi^{*}a $$
- **what:** GNEIMO spatial-operator-algebra closed form for the joint accelerations that avoids inverting the dense mass matrix; RHS is evaluated by O(N) recursive articulated-body sweeps (base-to-tip and tip-to-base).
- **symbols:** $\mathcal{H}$ - joint/hinge map matrix (projects spatial velocities onto joint dofs); $\psi$ - articulated-body/rigid-body transformation propagation operator (from Newton-Euler factorization); $\mathcal{K}$ - Kalman-gain-like operator from the mass-matrix factorization; $\mathcal{D}$ - block-diagonal articulated hinge inertia (invertible per hinge); $\mathcal{P}$ - articulated-body inertia operator; $a$ - Coriolis/gyroscopic spatial acceleration bias vector; $b$ - spatial-force bias (gyroscopic force) vector; $I$ - identity; $\mathcal{K}^{*},\psi^{*}$ - adjoints/transposes. Full definitions in Jain 2010 (Robot and Multibody Dynamics) and Jain-Vaidehi-Rodriguez 1993.
<!-- CHECK: operator-algebra identities (H, psi, K, D, P, a, b) are only named in this review; exact per-symbol definitions must be taken from Jain 2010 / Jain 1993 before implementing. -->

<!-- eq:3 -->
$$ \mathcal{M}(\theta)\ddot{\theta} + C(\theta, \dot{\theta}) + \mathcal{F}(\theta, \dot{\theta}) = \mathcal{T}(\theta) $$
- **what:** Nosé-Hoover constant-temperature (NVT) torsional MD equations of motion; adds a thermostat frictional force to eq 1. Same spatial-operator factorization as eq 1/2 applies.
- **symbols:** $\mathcal{F}(\theta,\dot\theta)$ - additional frictional (thermostat) force term for the canonical ensemble, proportional to the thermostat variable $\eta$ times generalized momenta; other symbols as in eq 1.

<!-- eq:4 -->
$$ \dot{\eta} = \frac{1}{\tau^{2}} \left[ \frac{T}{T_{B}} - 1 \right] $$
- **what:** Nosé-Hoover thermostat evolution: the thermostat variable's rate is driven by the ratio of instantaneous to target temperature.
- **symbols:** $\eta$ - thermostat dynamic variable (its time derivative $\dot\eta$ enters $\mathcal{F}$); $\tau$ - thermostat mass parameter (units of time); $T$ - instantaneous temperature; $T_{B}$ - thermostat/target (bath) temperature. Recommended $\tau = 10 \times \Delta t$ (time step).
<!-- CHECK: raw printed "η = (1/τ²)[T/T_B − 1]" with no dot; written as η̇ (rate) since η is described as the dynamic thermostat variable, matching standard Nosé-Hoover form. -->

<!-- eq:5 -->
$$ \mathcal{U}_{f}(\theta) \triangleq \frac{1}{2} k T \ln \frac{\det\{\mathcal{M}(\theta)\}}{\det\{\mathcal{M}_{B}(\theta, q_{0})\}} $$
- **what:** Fixman compensating potential removing constraint-induced bias in the partition function / probability density when stiff bond angles are treated as rigid; its gradient gives the Fixman torque applied during constrained MD.
- **symbols:** $\mathcal{U}_f(\theta)$ - Fixman potential (energy); $k$ - Boltzmann constant; $T$ - temperature; $\mathcal{M}(\theta)$ - mass matrix in the reduced (constrained) internal coordinates (eq 1); $\mathcal{M}_{B}(\theta,q_0)$ - mass matrix in the full BAT coordinates evaluated with frozen dofs at fixed values $q_0$; $q_0$ - coordinates of the frozen (constrained) degrees of freedom; $\det\{\cdot\}$ - determinant.
