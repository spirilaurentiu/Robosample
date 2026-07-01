# Equations - Wagner 2013, Advanced Techniques for Constrained ICMD (GNEIMO)

Spatial Operator Algebra (SOA) notation follows Jain (Robot and Multibody Dynamics, 2010).
Operators act over the stacked per-cluster (per-body) state; `*` denotes transpose/adjoint.

<!-- eq:1 -->
$$ \Re_e = \frac{1}{2}(\mathcal{N} - 6)\, k T $$
- **what:** target total thermal (kinetic) energy for a system with N generalized DOF, at temperature T; the `-6` removes the six rigid-body DOF of the whole system.
- **symbols:** $\Re_e$ - system kinetic energy (target); $\mathcal{N}$ - number of generalized DOF; $k$ - Boltzmann constant; $T$ - temperature.

<!-- eq:2 -->
$$ \Re_e = \frac{1}{2}\dot{\theta}^{*}\mathcal{M}\dot{\theta} = \frac{1}{2}\dot{\theta}^{*}[I + H\phi\mathcal{K}]\,\mathcal{D}\,[I + H\phi\mathcal{K}]^{*}\dot{\theta} $$
- **what:** CICMD kinetic energy as a quadratic form in generalized velocities; the mass matrix M is factored via SOA into `[I + HφK] D [I + HφK]*` (Newton-Euler operator factorization).
- **symbols:** $\dot{\theta}$ - generalized (hinge) velocities ($\mathcal{R}^{\mathcal{N}}$); $\mathcal{M}$ - mass matrix (dense, configuration-dependent); $I$ - identity; $H$ - joint/hinge map operator; $\phi$ - rigid-body transformation (SOA shift) operator; $\mathcal{K}$ - SOA gain operator; $\mathcal{D}$ - block-diagonal articulated hinge inertia.

<!-- eq:3 -->
$$ v \triangleq \mathcal{D}^{\frac{1}{2}}[I + H\phi\mathcal{K}]^{*}\dot{\theta} $$
- **what:** definition of modal velocity coordinates v (the transform that diagonalizes the kinetic energy); assign these from a Boltzmann distribution because they are independent.
- **symbols:** $v$ - modal velocity coordinates ($\mathcal{R}^{\mathcal{N}}$); $\mathcal{D}^{1/2}$ - matrix square root of block-diagonal articulated inertia; others as in eq:2.

<!-- eq:4 -->
$$ \dot{\theta} = [I - H\psi \mathcal{K}]^{*}\, \mathcal{D}^{-\frac{1}{2}}\, v $$
- **what:** inverse transform recovering generalized velocities from modal velocities; computable by an O(N) recursive base-to-tips scatter sweep.
- **symbols:** $\psi$ - SOA articulated transformation operator ($\psi = \phi(I - \mathcal{K}H\phi)^{-1}$-type; distinct from $\phi$); other symbols as above.

<!-- eq:5 -->
$$ \Re_e = \frac{1}{2} v^{*} v = \frac{1}{2} \sum_{k=1}^{\mathcal{N}} v(k)^{*}\, v(k) $$
- **what:** kinetic energy in modal coordinates is a plain sum of squares - so equipartition holds mode-by-mode and each v(k) is drawn independently.
- **symbols:** $v(k)$ - k-th modal velocity component; $\Re_e$ - kinetic energy.

<!-- eq:6 -->
$$ M_{S}\, \mathcal{V}_{CM} = \mathfrak{h}_{S} $$
- **what:** linear system relating system spatial inertia, CM spatial velocity, and system spatial momentum; solve for V_CM to detect/reset CM drift.
- **symbols:** $M_S$ - system 6x6 spatial inertia (referenced to base cluster); $\mathcal{V}_{CM}$ - CM spatial velocity (6-vector: angular+linear); $\mathfrak{h}_S$ - system spatial momentum (6-vector).

<!-- eq:7 -->
$$ M_{S} = \sum_{k=1}^{n} \phi(n,k)\, M(k)\, \phi^{*}(n,k) = E\phi M \phi^{*} E^{*} = E\Big[ \mathcal{R} + \tilde{\phi}\mathcal{R} + \mathcal{R}\tilde{\phi}^{*} \Big] E^{*} = E\mathcal{R}E^{*} = \mathcal{R}(n) $$
- **what:** the whole-system 6x6 spatial inertia equals the composite-rigid-body inertia R(n) accumulated at the base cluster n; computed by an O(N) tip-to-base sweep.
- **symbols:** $n$ - number of clusters (base cluster index); $\phi(j,k)$ - 6x6 rigid-body transform between cluster j and k; $M(k)$ - spatial inertia of cluster k; $\mathcal{R}(k)$ - composite rigid-body spatial inertia of cluster k and all its children; $\tilde{\phi} = \phi - I$; $E \triangleq [0_6,\dots,0_6, I_6] \in \mathcal{R}^{6\times 6n}$ - base pick-off operator; $\mathcal{R}$ - stacked composite-inertia operator.

<!-- eq:7b -->
$$ E \triangleq [0_6, \dots, 0_6, I_6] \in \mathcal{R}^{6 \times 6n}, \qquad E\phi = [\phi(n,1), \dots, \phi(n,n)], \qquad \tilde{\phi}\mathcal{R}E^{*} = 0 $$
- **what:** identities used to derive eq:7 (definition of the base pick-off operator E and two operator facts).
- **symbols:** as in eq:7; $I_6$ - 6x6 identity; $0_6$ - 6x6 zero block.

<!-- eq:hS -->
$$ \mathfrak{h}_{S} = \sum_{k=1}^{n} \phi(n,k)\, M(k)\, \mathcal{V}(k) = E\phi M \mathcal{V} = E\phi M \phi^{*} H^{*}\dot{\theta} = E\phi \mathcal{R} H^{*}\dot{\theta} $$
- **what:** base-cluster-referenced system spatial momentum, expressed via spatial velocities V(k) and then via generalized velocities.
- **symbols:** $\mathcal{V}(k)$ - spatial velocity (6-vector) of cluster k; $H^{*}$ - transpose hinge map; other symbols as above.

<!-- eq:hS2 -->
$$ \mathfrak{h}_{S} = \mathcal{R}(n)\,\mathcal{V}(n) + \sum_{k=1}^{n-1} \phi(n,k)\, \mathcal{R}(k)\, H^{*}(k)\, \dot{\theta}(k) $$
- **what:** momentum split into the base-cluster contribution plus contributions from all other hinges; uses that for a 6-DOF base hinge $H^{*}(n)=I$ and $\dot{\theta}(n)=\mathcal{V}(n)$.
- **symbols:** $\mathcal{V}(n)$ - base cluster spatial velocity; $H^{*}(k)$ - hinge map of cluster k; others as above.

<!-- eq:8 -->
$$ \mathfrak{h}_{S} = M_{S}\, \mathcal{V}_{CM} = \mathcal{R}(n)\, \mathcal{V}_{CM} $$
- **what:** restatement of eq:6 using $M_S = \mathcal{R}(n)$ from eq:7.
- **symbols:** as above.

<!-- eq:9 -->
$$ \mathcal{V}_{CM} = \mathcal{R}^{-1}(n)\left[ \mathcal{R}(n)\mathcal{V}(n) + \sum_{k=1}^{n-1} \phi(n,k)\, \mathcal{R}(k)\, H^{*}(k)\, \dot{\theta}(k) \right] = \mathcal{V}(n) + \mathcal{R}^{-1}(n) \sum_{k=1}^{n-1} \phi(n,k)\, \mathcal{R}(k)\, H^{*}(k)\, \dot{\theta}(k) $$
- **what:** explicit CM spatial velocity for an isolated molecule; used to compute the correction. Nulling: apply an extra base-cluster spatial velocity $\delta_V = -\mathcal{V}_{CM}$ (see eq:deltaV) to zero total spatial momentum.
- **symbols:** $\mathcal{R}^{-1}(n)$ - inverse of composite base inertia; others as above.

<!-- eq:deltaV -->
$$ 0 = \mathfrak{h}_{S} + \mathcal{R}(n)\delta_{V} \;\Rightarrow\; \delta_{V} = -\mathcal{R}^{-1}(n)\,\mathfrak{h}_{S} = -\mathcal{V}_{CM} $$
- **what:** the base-cluster spatial-velocity correction that resets system spatial momentum to zero.
- **symbols:** $\delta_V$ - spatial velocity increment applied to base cluster; others as above.

<!-- eq:10 -->
$$ \frac{d\, m_i \mathbf{v}_i}{dt} = \mathbf{F}_i - \zeta\, m_i \mathbf{v}_i $$
- **what:** Nose-Hoover equation of motion (mass-weighted) for particle i in Cartesian coordinates.
- **symbols:** $m_i$ - mass of particle i; $\mathbf{v}_i$ - velocity of particle i; $\mathbf{F}_i$ - internal force on particle i; $\zeta$ - thermostat friction (bath) coefficient.

<!-- eq:11 -->
$$ \frac{d\, \mathcal{M}\mathcal{V}_{CM}}{dt} = -\zeta\, \mathcal{M}\, \mathcal{V}_{CM} $$
- **what:** summing eq:10 over all particles (using $\sum_i \mathbf{F}_i = 0$) gives the CM-momentum evolution under the thermostat.
- **symbols:** $\mathcal{M} = \sum_i m_i$ - total system mass; $\mathcal{V}_{CM}$ - CM velocity; $\zeta$ - thermostat friction.

<!-- eq:12 -->
$$ \mathcal{V}_{CM}(t) = \mathcal{V}_{CM}(t=0)\, \exp\!\left[-\int_{0}^{t} \zeta(t')\,dt'\right] = \mathcal{V}_{CM}(t=0)\, \exp\!\left[-\ln(s_t)\right] $$
- **what:** integrated CM velocity; negative $\ln(s_t)$ makes the exponential grow -> flying-ice-cube instability in Nose-Hoover.
- **symbols:** $s_t$ - Nose bath variable (proportional to bath potential energy); $\zeta(t')$ - time-dependent friction; $t$ - time.

<!-- eq:13 -->
$$ KE_{CM} = c\, \exp(-2 \ln s_t) $$
- **what:** CM kinetic energy vs bath variable; fit constant c to simulation to diagnose the flying-ice-cube growth.
- **symbols:** $KE_{CM}$ - CM kinetic energy; $c$ - constant proportional to initial CM kinetic energy; $s_t$ - Nose bath variable.

<!-- eq:bfactor -->
$$ B = \frac{8\pi^{2}}{3}\, \mathrm{RMSF}^{2} $$
- **what:** relation between crystallographic B-factor and per-residue RMSF; used to convert experimental B-factors to expected fluctuations for validation.
- **symbols:** $B$ - crystallographic B-factor of a residue; $\mathrm{RMSF}$ - root-mean-square fluctuation of the residue position.
