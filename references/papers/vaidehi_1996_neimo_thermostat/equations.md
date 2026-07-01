# Equations — Vaidehi, Jain, Goddard 1996 (NEIMO constant-temperature)

Reduced/MD units: length = Å, time unit = 0.0488 ps, energy in kcal/mol, $k$ = Boltzmann constant. Primes = virtual (Nosé) variables; unprimed = real variables (opposite the convention of Nosé 1984). Repeated indices $i,j$ sum over the $\mathcal{N}$ generalized coordinates; $\alpha$ sums over the $3N$ Cartesian coordinates.

<!-- eq:1 -->
$$ \mathcal{M}(\boldsymbol{\theta})\ddot{\boldsymbol{\theta}} + \mathbf{C}(\boldsymbol{\theta},\dot{\boldsymbol{\theta}}) = \mathcal{T}(\boldsymbol{\theta}) $$
- **what:** Constrained (internal-coordinate) Newton equation of motion; mass matrix times acceleration plus Coriolis equals generalized force.
- **symbols:** $\boldsymbol{\theta}$ - generalized coords (torsions), $\mathbb{R}^{\mathcal{N}}$; $\dot{\boldsymbol{\theta}},\ddot{\boldsymbol{\theta}}$ - generalized velocity, acceleration; $\mathcal{M}$ - mass matrix / moment-of-inertia tensor ($\mathcal{N}\times\mathcal{N}$); $\mathbf{C}$ - Coriolis (velocity-dependent) forces; $\mathcal{T}$ - generalized force (torques).

<!-- eq:2 -->
$$ \ddot{\boldsymbol{\theta}} = \mathcal{M}^{-1}[\mathcal{T}(\boldsymbol{\theta}) - \mathbf{C}(\boldsymbol{\theta},\dot{\boldsymbol{\theta}})] $$
- **what:** Solve eq:1 for the acceleration. Direct $\mathcal{M}^{-1}$ costs $O(\mathcal{N}^3)$; NEIMO computes this in $O(\mathcal{N})$ without forming/inverting $\mathcal{M}$.
- **symbols:** as eq:1; $\mathcal{M}^{-1}$ - inverse mass matrix. <!-- CHECK: raw printed "M^1"; superscript is -1 -->

<!-- eq:3 -->
$$ \frac{\mathrm{d}}{\mathrm{d}t} = s\,\frac{\mathrm{d}}{\mathrm{d}t_s} $$
- **what:** Nosé real-time / virtual-time relation ($\mathrm{d}t = \mathrm{d}t_s / s$).
- **symbols:** $t$ - real time; $t_s$ - virtual (Nosé) time; $s$ - time-scaling bath coordinate (dimensionless).

<!-- eq:4 -->
$$ \theta_i = \theta_i', \qquad \dot{\theta}_i = s\,\dot{\theta}_i' $$
- **what:** Definition of virtual Nosé variables (eq:4a, eq:4b): positions equal, real velocity = $s$ × virtual velocity.
- **symbols:** primed = virtual, unprimed = real.

<!-- eq:5a -->
$$ \mathcal{L} = \mathrm{KE} - \mathrm{PE} = \frac{1}{2}s^2\sum_{i,j=1}^{\mathcal{N}}\dot{\theta}_i' M_{ij}(\boldsymbol{\theta}')\dot{\theta}_j' - \Phi(\boldsymbol{\theta}') + \frac{1}{2}Q(\dot{s}')^2 - gkT_B\ln s $$
- **what:** Nosé extended-system Lagrangian in virtual variables: physical KE, physical PE, bath kinetic term, bath potential term.
- **symbols:** $M_{ij}$ - mass-matrix elements; $\Phi$ - physical potential energy; $Q$ - Nosé (thermostat) mass; $g$ - thermostat dof count; $k$ - Boltzmann constant; $T_B$ - bath temperature; $s$ - bath coordinate; $\dot{s}'$ - virtual bath velocity.

<!-- eq:5b -->
$$ g = \mathcal{N} + 1 $$
- **what:** Nosé choice for the extended-system dof so the marginal is exactly canonical (Nosé variables). For the real-variable/Hoover formulation use $g = \mathcal{N}$ instead (see eq:29).
- **symbols:** $\mathcal{N}$ - number of generalized (internal) dof.

<!-- eq:12 -->
$$ \sum_j M_{kj}\ddot{\theta}_j' + \mathbf{C}_N(\boldsymbol{\theta},\dot{\boldsymbol{\theta}}') + \mathbf{F}_N(\boldsymbol{\theta},\dot{\boldsymbol{\theta}}') = \mathcal{T}_N(\boldsymbol{\theta}) $$
- **what:** NEIMO–Nosé equation of motion for the generalized coordinates (row $k$). Same structure as eq:1 with an extra friction-like term $\mathbf{F}_N$.
- **symbols:** $\mathbf{C}_N$ - Coriolis (eq:13); $\mathbf{F}_N$ - velocity-dependent friction-like force (eq:14); $\mathcal{T}_N$ - scaled generalized force (eq:15).

<!-- eq:13 -->
$$ \mathbf{C}_N(\boldsymbol{\theta},\dot{\boldsymbol{\theta}}') = \sum_j \dot{M}_{kj}\dot{\theta}_j' - \frac{1}{2}\sum_{ij}\dot{\theta}_i'\frac{\partial M_{ij}}{\partial\theta_k}\dot{\theta}_j' + \frac{1}{s^2}\sum_\alpha \frac{\partial\Phi_{\rm NB}}{\partial X_\alpha}\frac{\partial X_\alpha}{\partial\theta_k} $$
- **what:** NEIMO–Nosé Coriolis term; last term carries nonbond/external forces mapped from Cartesian to internal coordinates.
- **symbols:** $\dot{M}_{kj}$ - time derivative of mass-matrix element; $\Phi_{\rm NB}$ - nonbond potential (Coulomb + vdW + external); $X_\alpha$ - Cartesian coordinate $\alpha$. <!-- CHECK: raw eq(13) printed the derivative as ∂M_ij/∂θ̇'_k; from eq(7)/eq(9)/eq(24) it is ∂M_ij/∂θ_k -->

<!-- eq:14 -->
$$ \mathbf{F}_N(\boldsymbol{\theta},\dot{\boldsymbol{\theta}}') = \frac{2}{s}\frac{\mathrm{d}s}{\mathrm{d}t_s}\sum_j M_{kj}\dot{\theta}_j' $$
- **what:** NEIMO–Nosé friction-like force (may be folded into the Coriolis term).
- **symbols:** $\mathrm{d}s/\mathrm{d}t_s$ - virtual-time derivative of bath coordinate.

<!-- eq:15 -->
$$ \mathcal{T}_N(\boldsymbol{\theta}) = -\frac{1}{s^2}\frac{\partial\Phi_{\rm int}}{\partial\theta_k'} = -\frac{1}{s^2}\frac{\partial\Phi_{\rm int}}{\partial\theta_k} $$
- **what:** NEIMO–Nosé generalized force (internal/torsional part only, scaled by $1/s^2$).
- **symbols:** $\Phi_{\rm int}$ - internal (torsional) potential energy. <!-- CHECK: raw printed "1/c^2"; c is OCR of s -->

<!-- eq:16 -->
$$ \ddot{s} = \frac{1}{Qs}\left[ s^2\sum_{ij}\dot{\theta}_i' M_{ij}\dot{\theta}_j' - gkT_B \right] $$
- **what:** Nosé equation of motion for the bath coordinate $s$ (from the $s$ Lagrange equation).
- **symbols:** $\ddot{s}$ - bath acceleration; other symbols as eq:5a.

<!-- eq:17 -->
$$ \mathrm{KE} = \frac{1}{2}s^2\sum_{ij}\dot{\theta}_i' M_{ij}(\boldsymbol{\theta}')\dot{\theta}_j' = \frac{1}{2}\mathcal{N}kT $$
- **what:** Physical kinetic energy; defines the instantaneous temperature $T$.
- **symbols:** $T$ - instantaneous temperature; $\mathcal{N}$ - number of internal dof.

<!-- eq:18 -->
$$ \ddot{s} = \frac{1}{Qs}\left[ \mathcal{N}kT - gkT_B \right] $$
- **what:** Bath equation of motion with KE expressed via instantaneous temperature (eq:16 + eq:17).
- **symbols:** as above.

<!-- eq:19 -->
$$ \tau_s^2 = \frac{Q}{\mathcal{N}kT_B} $$
- **what:** Thermostat mass expressed as a relaxation time constant $\tau_s$ of the bath variable $s$.
- **symbols:** $\tau_s$ - bath relaxation time (ps); $Q$ - thermostat mass.

<!-- eq:20 -->
$$ \ddot{s} = \frac{1}{s\,\tau_s^2}\left[ \frac{T}{T_B} - \frac{g}{\mathcal{N}} \right] $$
- **what:** Bath equation of motion in terms of $\tau_s$ (substitute eq:19 into eq:18).
- **symbols:** as above.

<!-- eq:21 -->
$$ \ddot{s} = \frac{1}{s\,\tau_s^2}\left[ \left(\frac{T}{T_B} - 1\right) - \frac{1}{\mathcal{N}} \right] $$
- **what:** NEIMO–Nosé bath equation of motion with $g = \mathcal{N}+1$ (eq:5b). Fundamental Nosé eq together with eq:12.
- **symbols:** as above.

<!-- eq:23 -->
$$ \sum_j M_{kj}\ddot{\theta}_j + \mathbf{C}_H(\boldsymbol{\theta},\dot{\boldsymbol{\theta}}) + \mathbf{F}_H(\boldsymbol{\theta},\dot{\boldsymbol{\theta}}) = \mathcal{T}_H $$
- **what:** NEIMO–Hoover equation of motion (real variables); same form as eq:1 with friction term $\mathbf{F}_H$.
- **symbols:** $\mathbf{C}_H$ - Coriolis (eq:24); $\mathbf{F}_H$ - friction (eq:25); $\mathcal{T}_H$ - generalized force (eq:26).

<!-- eq:24 -->
$$ \mathbf{C}_H(\boldsymbol{\theta},\dot{\boldsymbol{\theta}}) = \sum_j \dot{M}_{kj}\dot{\theta}_j - \frac{1}{2}\sum_{ij}\dot{\theta}_i\frac{\partial M_{ij}}{\partial\theta_k}\dot{\theta}_j + \sum_\alpha \frac{\partial\Phi_{\rm NB}}{\partial X_\alpha}\frac{\partial X_\alpha}{\partial\theta_k} $$
- **what:** NEIMO–Hoover Coriolis term (real variables, no $1/s^2$ scaling on nonbond term).
- **symbols:** as eq:13.

<!-- eq:25 -->
$$ \mathbf{F}_H(\boldsymbol{\theta},\dot{\boldsymbol{\theta}}) = \zeta \sum_j M_{kj}\dot{\theta}_j $$
- **what:** NEIMO–Hoover friction-like force; $\zeta$ is the thermostat friction coefficient. Folded into the Coriolis terms.
- **symbols:** $\zeta$ - Hoover friction coefficient (1/time).

<!-- eq:26 -->
$$ \mathcal{T}_H(\boldsymbol{\theta}) = -\frac{\partial\Phi_{\rm int}}{\partial\theta_k} $$
- **what:** NEIMO–Hoover generalized force (internal/torsional gradient only; nonbond forces enter via $\mathbf{C}_H$).
- **symbols:** as eq:15.

<!-- eq:27 -->
$$ \dot{\zeta} = \frac{1}{\tau_s^2}\left[ \frac{T}{T_B} - \frac{g}{\mathcal{N}} \right] $$
- **what:** Hoover friction equation of motion (real variables; transform of eq:20).
- **symbols:** $\dot{\zeta}$ - time derivative of friction coefficient.

<!-- eq:28 -->
$$ \zeta = \frac{1}{s}\frac{\mathrm{d}s}{\mathrm{d}t_s} $$
- **what:** Definition of the Hoover friction coefficient in terms of the bath coordinate.
- **symbols:** as above.

<!-- eq:29 -->
$$ \dot{\zeta} = \frac{1}{\tau_s^2}\left[ \frac{T}{T_B} - 1 \right] $$
- **what:** NEIMO–Hoover friction equation of motion with $g = \mathcal{N}$ (Liouville-conserved density). Fundamental Hoover eq with eq:23.
- **symbols:** as above.

<!-- eq:30 -->
$$ \mathcal{M}^{-1} = [\mathbf{I} - \mathbf{H}\psi\mathbf{K}]^{\mathrm{T}}\,\mathcal{D}^{-1}\,[\mathbf{I} - \mathbf{H}\psi\mathbf{K}] $$
- **what:** Innovations operator factorization: closed-form block $\mathbf{LDL}^T$ decomposition of the inverse mass matrix (basis of the $O(\mathcal{N})$ NEIMO recursion).
- **symbols:** $\mathbf{I}$ - identity; $\mathbf{H}$ - hinge matrix ($m$ dof/hinge, $m=1$ torsion); $\psi$ - lower-diagonal spatial transformation matrix; $\mathbf{K}$ - spatial operator (body inertia + hinge characteristics); $\mathcal{D}$ - block-diagonal, $\mathcal{N}$ subblocks of $m\times m$.

<!-- eq:32 -->
$$ \ddot{\boldsymbol{\theta}} = \mathcal{M}^{-1}[\mathcal{T} - \mathbf{C}_H - \zeta\mathcal{M}\dot{\boldsymbol{\theta}}] $$
- **what:** NEIMO–Hoover acceleration from eq:23.
- **symbols:** as above.

<!-- eq:33 -->
$$ \ddot{\boldsymbol{\theta}} = \mathcal{M}^{-1}[\mathcal{T} - \mathbf{C}_H] - \zeta\dot{\boldsymbol{\theta}} $$
- **what:** NEIMO–Hoover acceleration = microcanonical acceleration minus friction correction $\zeta\dot{\boldsymbol{\theta}}$. The single change from microcanonical NEIMO.
- **symbols:** as above.

<!-- eq:34 -->
$$ \ddot{\boldsymbol{\theta}} = [\mathbf{I} - \mathbf{H}\psi\mathbf{K}]^{\mathrm{T}}\mathcal{D}^{-1}[\mathcal{T} - \mathbf{H}\psi(\mathbf{K}\mathcal{T} + \mathbf{P}\mathbf{a} + \mathbf{b} + \hat{\mathbf{f}}_c)] - \mathbf{K}^{\mathrm{T}}\psi^{\mathrm{T}}\mathbf{a} - \zeta\dot{\boldsymbol{\theta}}$$
- **what:** Operator form of the NEIMO–Hoover acceleration (eq:30 into eq:33).
- **symbols:** $\mathbf{P}$ - articulated-body inertia; $\mathbf{a}$ - Coriolis force; $\mathbf{b}$ - gyroscopic force; $\hat{\mathbf{f}}_c$ - Cartesian/external force contribution.

<!-- eq:34seq -->
$$ \mathbf{z} = \psi[\mathbf{K}\mathcal{T} + \mathbf{P}\mathbf{a} + \mathbf{b} + \hat{\mathbf{f}}_c],\quad \epsilon = \mathcal{T} - \mathbf{H}\mathbf{z},\quad \mathbf{v} = \mathcal{D}^{-1}\epsilon,\quad \alpha = \psi[\mathbf{H}^{t}\mathbf{v} + \mathbf{a}],\quad \ddot{\boldsymbol{\theta}} = \mathbf{v} - \mathbf{K}^{t}\alpha - \zeta\dot{\boldsymbol{\theta}} $$
- **what:** Recursive $O(\mathcal{N})$ evaluation sequence building eq:34. Each step is a base-to-tips or tips-to-base sweep. Final $-\zeta\dot{\boldsymbol{\theta}}$ is the only thermostat addition to microcanonical NEIMO.
- **symbols:** $\mathbf{z},\epsilon,\mathbf{v},\alpha$ - intermediate spatial vectors; $\mathbf{H}^t,\mathbf{K}^t$ - transposes.

<!-- eq:verlet-est -->
$$ \dot{\boldsymbol{\theta}}_n = 1.5\,\dot{\boldsymbol{\theta}}_{n-1/2} - 0.5\,\dot{\boldsymbol{\theta}}_{n-3/2} $$
- **what:** Predictor for current-step velocity (needed by NEIMO to form Coriolis forces before computing accelerations).
- **symbols:** subscripts denote leapfrog half/full steps.

<!-- eq:verlet-recorrect -->
$$ \dot{\boldsymbol{\theta}}_n = 0.5\,\dot{\boldsymbol{\theta}}_{n-1/2} + 0.5\,\dot{\boldsymbol{\theta}}_{n+1/2} $$
- **what:** Corrector re-estimate of current-step velocity after computing $\dot{\boldsymbol{\theta}}_{n+1/2}$; iterate predictor/corrector until velocity change < 0.001 MD units (1–2 iterations typical).
- **symbols:** as above.

<!-- eq:35 -->
$$ \zeta_{n+1} = \zeta_n + \delta\,D_{n+1/2} $$
- **what:** Verlet update of the Hoover friction coefficient.
- **symbols:** $\delta$ - time step; $D_{n+1/2}$ - half-step $\dot{\zeta}$ (eq:36).

<!-- eq:36 -->
$$ D_{n+1/2} = \frac{1}{\tau_s^2}\left[ \frac{T_{n+1/2}}{T_B} - 1 \right] $$
- **what:** Half-step $\dot{\zeta}$ value (eq:29 evaluated at the half-step temperature).
- **symbols:** $T_{n+1/2}$ - instantaneous temperature at the $n+1/2$ step.

<!-- eq:38 -->
$$ \tau_s \ge \frac{10}{2\pi}\delta = 1.6\,\delta $$
- **what:** Lower bound on $\tau_s$: at least 10 integration steps per $\tau_s$ period (from eq:37 $10\delta \le 2\pi\tau_s$). Below this the dynamics blows up.
- **symbols:** $\delta$ - integration time step; $\tau_s$ - bath relaxation time.

<!-- eq:40 -->
$$ \tau_{\rm slong} = \sqrt{\mathcal{N}}\,\tau_s $$
- **what:** Characteristic long-time relaxation time of $\langle s\rangle$ (harmonic behavior from eq:39).
- **symbols:** $\tau_{\rm slong}$ - long relaxation time.

<!-- eq:41 -->
$$ t_{\rm total} = 20(2\pi\tau_{\rm slong}) = 40\pi\sqrt{\mathcal{N}}\,\tau_s $$
- **what:** Minimum total simulation time for good averaging (20 periods of the long relaxation).
- **symbols:** $t_{\rm total}$ - required simulation length.

<!-- eq:KE-avg -->
$$ \langle\mathrm{KE}\rangle = \frac{\mathcal{N}}{2}kT_B $$
- **what:** Equipartition: mean kinetic energy for $\mathcal{N}$ dof.
- **symbols:** as above.

<!-- eq:KE-fluc -->
$$ \langle(\delta\mathrm{KE})^2\rangle = \langle\mathrm{KE}^2\rangle - \langle\mathrm{KE}\rangle^2 = \frac{\mathcal{N}}{2}(kT_B)^2 $$
- **what:** Canonical kinetic-energy fluctuation.
- **symbols:** as above.

<!-- eq:43 -->
$$ \langle T_{\rm calc}\rangle = \frac{2}{\mathcal{N}k}\langle\mathrm{KE}\rangle = T_B $$
- **what:** Calculated instantaneous temperature averages to the bath temperature.
- **symbols:** $\langle T_{\rm calc}\rangle$ - time-averaged temperature.

<!-- eq:dT -->
$$ \langle\delta T_{\rm calc}^2\rangle = \frac{2}{\mathcal{N}}T_B^2 $$
- **what:** Mean-square deviation of the calculated temperature.
- **symbols:** as above.

<!-- eq:44 -->
$$ \sqrt{\frac{\mathcal{N}}{2}\langle\delta T_{\rm calc}^2\rangle} = T_B $$
- **what:** Normalized temperature-fluctuation observable that should equal $T_B$ at equilibrium (test statistic used in Tables 1–3).
- **symbols:** as above.
