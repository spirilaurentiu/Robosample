# Equations - Kandel et al. 2016, Fixman potential + hybrid ICMD

Mass matrix is denoted $\mathcal{M}$ (constrained, $N\times N$) and $\mathcal{M}_B$
(full $3n$-dimensional / unconstrained BAT mass matrix). `*` denotes transpose.
`kT` is Boltzmann's constant times temperature.

<!-- eq:1 -->
$$ \mathcal{M}(\alpha)\ddot{\alpha} + \mathcal{C}(\alpha,\dot{\alpha}) = \mathcal{T}(\alpha) $$
- **what:** Coupled ICMD equations of motion in generalized (internal) coordinates.
- **symbols:** $\alpha$ - vector of generalized coords (torsions + any open bond angles); $\ddot\alpha$ - generalized acceleration; $\dot\alpha$ - generalized velocity; $\mathcal{M}(\alpha)$ - mass matrix / articulated moment-of-inertia tensor ($R^{N\times N}$); $\mathcal{C}(\alpha,\dot\alpha)$ - velocity-dependent Coriolis forces; $\mathcal{T}(\alpha)$ - generalized forces (torques).

<!-- eq:2 -->
$$ \ddot{\alpha} = \left[I - H\psi \mathcal{K}\right]^{*} \mathcal{D}^{-1} \left[\mathcal{T} - H\psi (\mathcal{K}\mathcal{T} + \mathcal{P}\mathfrak{a} + \mathfrak{b})\right] - \mathcal{K}^{*} \psi^{*} \mathfrak{a} $$
- **what:** GNEIMO spatial-operator-algebra closed form for the acceleration; evaluated by recursive O(N) sweeps (same recursion whether or not bond angles are open).
- **symbols:** $H$ - joint map matrix; $\psi$ - spatial propagation operator; $\mathcal{K}$ - Kalman-gain-like operator; $\mathcal{D}$ - articulated hinge inertia (block diagonal); $\mathcal{P}$ - articulated body inertia operator; $\mathfrak{a}$ - Coriolis/gyroscopic spatial acceleration; $\mathfrak{b}$ - spatial gyroscopic force; $I$ - identity. Factorization details in Jain refs 21, 26.
<!-- CHECK: operator terms transcribed from OCR; exact operator identities defined in Jain 1993/2010, not fully self-contained here -->

<!-- eq:3 -->
$$ \mathcal{H}_u(\alpha,q,\mathfrak{p}) = \tfrac{1}{2}\,\mathfrak{p}^{*}\mathcal{M}_B^{-1}(\alpha,q)\,\mathfrak{p} + \mathcal{U}(\alpha,q) $$
- **what:** Unconstrained Hamiltonian in BAT coordinates.
- **symbols:** $\alpha$ - unconstrained coords; $q$ - (to-be-)constrained coords; $\mathfrak{p}$ - canonical momenta conjugate to $(\alpha,q)$; $\mathcal{M}_B$ - full $3n$-dim mass matrix; $\mathcal{U}$ - forcefield potential.

<!-- eq:4 -->
$$ \mathcal{Z}_u(T) = c_1 \int dp\; d\alpha\; dq\; e^{-\mathcal{H}_u(\alpha,q,p)/kT} $$
- **what:** Unconstrained canonical partition function.
- **symbols:** $\mathcal{Z}_u$ - partition function; $c_1$ - scaling constant; $k$ - Boltzmann constant; $T$ - temperature.

<!-- eq:5 -->
$$ \mathcal{Z}_u(T) = c_2 \int d\alpha\; dq\; \det\!\left\{\mathcal{M}_B^{\frac{1}{2}}(\alpha,q)\right\} e^{-\mathcal{U}(\alpha,q)/kT} $$
- **what:** Partition function after Gaussian integration over momenta; $\det\{\mathcal{M}_B^{1/2}\}=(\det\{\mathcal{M}_B\})^{1/2}$.
- **symbols:** $c_2$ - constant absorbing momentum-integral factors.

<!-- eq:6 -->
$$ \rho_u(\alpha,q) \;\propto\; \det\!\left\{\mathcal{M}_B^{\frac{1}{2}}(\alpha,q)\right\} e^{-\mathcal{U}(\alpha,q)/kT} $$
- **what:** Unconstrained configuration pdf (reference distribution to recover).
- **symbols:** $\rho_u$ - unconstrained configuration probability density.

<!-- eq:7 -->
$$ \det\{\mathcal{M}_B\} = \sin^2\gamma_2\; d_2^4 \prod_{i=3}^{n} d_i^4 \sin^2\theta_i \prod_{i=1}^{n} m_i^3 $$
- **what:** Closed-form determinant of the full BAT mass matrix; note it is independent of torsion angles.
- **symbols:** $d_i$ - the $(n-1)$ bond lengths; $\theta_i$ - the $(n-2)$ bond angles; $m_i$ - atomic masses; $(\gamma_1,\gamma_2,\gamma_3)$ - ZXZ Euler angles for overall molecule orientation; $n$ - number of atoms.

<!-- eq:8 -->
$$ \det\{\mathcal{M}_B(\alpha,q)\} = f_1(\alpha)\,f_2(q) $$
- **what:** The full BAT determinant factorizes into an $\alpha$-only and a $q$-only part (consequence of eq:7).
- **symbols:** $f_1,f_2$ - factor functions of the unconstrained and constrained coords respectively.

<!-- eq:9 -->
$$ \mathcal{Z}_c(T) = c_3 \int d\alpha\; \det\!\left\{\mathcal{M}^{\frac{1}{2}}(\alpha,q_0)\right\} e^{-\mathcal{U}(\alpha,q_0)/kT} $$
- **what:** Constrained partition function ($q$ frozen at $q_0$); uses the reduced constrained mass matrix.
- **symbols:** $\mathcal{M}(\alpha,q_0)\in R^{N\times N}$ - constrained-model mass matrix; $q_0$ - frozen constrained-coord values; $c_3$ - constant.

<!-- eq:10 -->
$$ \rho_c(\alpha,q_0) \;\propto\; \det\!\left\{\mathcal{M}^{\frac{1}{2}}(\alpha,q_0)\right\} e^{-\mathcal{U}(\alpha,q_0)/kT} $$
- **what:** Constrained configuration pdf; the $\det\{\mathcal{M}\}$ factor does NOT factorize -> intrinsic distortion (biased vs $\rho_u$).
- **symbols:** $\rho_c$ - constrained configuration pdf.

<!-- eq:11 -->
$$ \mathcal{U}_f(\alpha) = \tfrac{1}{2}\, kT \, \ln \frac{\det\{\mathcal{M}(\alpha,q_0)\}}{\det\{\mathcal{M}_B(\alpha,q_0)\}} $$
- **what:** The Fixman compensating potential; add to the forcefield ($U'(\alpha)=U(\alpha,q_0)+U_f(\alpha)$) in constrained dynamics to cancel the intrinsic mass-matrix bias.
- **symbols:** $\mathcal{U}_f$ - Fixman potential (energy units); $\mathcal{M}$ - constrained mass matrix; $\mathcal{M}_B$ - full BAT mass matrix; ratio evaluated at frozen $q_0$.

<!-- eq:12 -->
$$ \rho_f(\alpha,q_0) \;\propto\; \det\!\left\{\mathcal{M}_B^{\frac{1}{2}}(\alpha,q_0)\right\} e^{-\mathcal{U}(\alpha,q_0)/kT} $$
- **what:** Fixman-compensated constrained pdf; matches $\rho_u$ (eq:6) when $q\approx q_0$ / separable / no-force conditions hold.
- **symbols:** $\rho_f$ - Fixman-corrected constrained configuration pdf.

<!-- eq:13 -->
$$ \mathcal{U}_f(\alpha) = c_f + \tfrac{1}{2}\, kT \, \ln \det\{\mathcal{M}(\alpha)\}, \qquad \rho_f(\alpha) \;\propto\; e^{-\mathcal{U}(\alpha,q_0)/kT} $$
- **what:** Simplified Fixman potential when unconstrained coords $\alpha$ are pure torsions (TMD): since $\det\{\mathcal{M}_B\}$ is torsion-independent, only $\det\{\mathcal{M}(\alpha)\}$ (constrained) varies; $c_f$ is a constant bond/angle contribution.
- **symbols:** $c_f$ - constant (bond-length + bond-angle contributions of $\det\{\mathcal{M}_B\}$).

<!-- eq:14 -->
$$ \rho_u(\alpha_i) = \frac{1}{2\pi} $$
- **what:** With zero torsional force ($U(\alpha)=0$), the unconstrained torsion-angle pdf is uniform; regression target for Fixman correctness.
- **symbols:** $\alpha_i$ - a single torsion angle in radians.

<!-- eq:15 -->
$$ U(\alpha) = k_\alpha\bigl(1 + \cos(\alpha - \alpha_0)\bigr) $$
- **what:** Single-barrier harmonic torsional test potential (separable) used for the C4 chain.
- **symbols:** $k_\alpha$ - barrier peak amplitude (kcal/mol); $\alpha_0$ - location of barrier peak; $\alpha$ - torsion angle.

<!-- eq:16 -->
$$ f_{TS}(x_0) = \int_0^\infty \dot{x}\, \rho(x=x_0,\dot{x})\, d\dot{x} $$
- **what:** Transition-state barrier-crossing rate for a generic 1D well (flux of positive-velocity crossings at barrier center).
- **symbols:** $x$ - reaction coordinate; $x_0$ - barrier center; $\dot x$ - its velocity; $\rho(x,\dot x)$ - phase-space pdf.

<!-- eq:17 -->
$$ f_{TS}(\alpha_0) = \int_0^\infty \dot{\alpha}\, d\dot{\alpha} \iiint_{-\infty}^{\infty} d\dot{\gamma}_0\, d\dot{\gamma}_1\, d\dot{\gamma}_2\; \rho(\alpha=\alpha_0,\dot{\gamma},\dot{\alpha}) $$
- **what:** Barrier-crossing rate for the C4 torsion, integrating over Euler-angle velocities.
- **symbols:** $\alpha$ - torsion; $\alpha_0$ - barrier peak; $\gamma=(\gamma_0,\gamma_1,\gamma_2)$ - Euler angles; $\dot\gamma_i$ - Euler-angle velocities; $\rho$ - pdf.

<!-- eq:18 -->
$$ f_{TS,\text{cons}}(\alpha_0) = C\, e^{-2k_\alpha/kT} \left[(2\pi)^3 (kT)^5 \det\{\mathcal{M}(\alpha_0)\}\, S^{-1}(\alpha_0)\right]^{1/2} $$
- **what:** Constrained-model barrier-crossing rate; depends on mass-matrix determinant (intrinsic bias) AND on $S^{-1}(\alpha_0)$.
- **symbols:** $C$ - normalization constant; $k_\alpha$ - barrier amplitude; $\det\{\mathcal{M}(\alpha_0)\}$ - constrained mass-matrix determinant at peak; $S^{-1}$ - see eq below.

<!-- eq:18b -->
$$ S^{-1}(\alpha_0) = \left[\mathcal{M}^{-1}(\alpha_0)\right]_{\alpha} $$
- **what:** $S^{-1}$ is the $(\alpha,\alpha)$ sub-block of the inverse constrained mass matrix (Schur complement inverse).
- **symbols:** $\left[\mathcal{M}^{-1}\right]_\alpha$ - torsion-block of $\mathcal{M}^{-1}$.

<!-- eq:19 -->
$$ f_{TS,\text{fix}}(\alpha_0) = C\, e^{-2k_\alpha/kT} \left[(2\pi)^3 (kT)^5\, S^{-1}(\alpha_0)\right]^{1/2} $$
- **what:** Fixman-corrected barrier-crossing rate: $\det\{\mathcal{M}\}$ cancels, but the residual $S^{-1}(\alpha_0)$ dependence on barrier location remains -> Fixman only PARTIALLY corrects velocity-dependent quantities.
- **symbols:** as eq:18.

<!-- eq:20 -->
$$ U_{\text{coul}} = \tfrac{1}{2}\, k_{\text{coul}}\, \frac{q_1 q_2}{r} $$
- **what:** Coulombic coupling test potential added to C4 terminal atoms (makes potential non-separable in torsion vs bond angle -> extrinsic distortion).
- **symbols:** $k_{\text{coul}}$ - Coulomb constant (332.06 kcal Å / e^2); $q_1,q_2$ - terminal charges (e); $r$ - distance between the two terminal beads (function of torsion AND bond angles).
<!-- CHECK: paper writes the 1/2 prefactor explicitly; standard Coulomb has no 1/2 -->

<!-- eq:21 -->
$$ \mathcal{U}_\theta = k_\theta (\theta - \theta_0)^2 $$
- **what:** Harmonic bond-angle spring potential (note: no 1/2 prefactor as written by this paper's convention).
- **symbols:** $\theta$ - bond angle; $\theta_0$ - equilibrium angle; $k_\theta$ - angle spring constant (kcal, tunes stiffness).

<!-- eq:22 -->
$$ \mathcal{U}_{ff} = \tfrac{1}{2}\sum_{i=1}^{N_{bonds}} K_{r_i}(r_i - r_i^0)^2 + \tfrac{1}{2}\sum_{i=1}^{N_{angles}} K_{\theta_i}(\theta_i - \theta_i^0)^2 + \mathcal{U}_{ff}^{tors}(\alpha_i) + \mathcal{U}_{ff}^{long\text{-}range} $$
- **what:** General all-atom forcefield form used; long-range term couples all BAT degrees of freedom.
- **symbols:** $K_{r_i}$ - bond spring constants; $r_i,r_i^0$ - bond lengths and equilibria; $K_{\theta_i}$ - angle spring constants; $\theta_i,\theta_i^0$ - bond angles and equilibria; $\mathcal{U}_{ff}^{tors}$ - dihedral potential; $\mathcal{U}_{ff}^{long\text{-}range}$ - Coulomb + van der Waals + implicit solvent.

## Derivations (not implemented) - Appendix identities

The Appendix derives eq:18/eq:19 from the constrained-model pdf. Implementable pieces:

<!-- eq:A2 -->
$$ \rho(\alpha,\dot{\alpha},\dot{\gamma}) = C\, \det\{\mathcal{M}(\alpha)\}\, e^{-[E_K(\alpha,\dot{\alpha},\dot{\gamma}) + \mathcal{U}(\alpha)]/kT} $$
- **what:** Velocity-space pdf of the constrained model; the $\det\{\mathcal{M}\}$ prefactor comes from the momentum->velocity Jacobian.
- **symbols:** $E_K$ - kinetic energy (eq:A3); $C$ - normalization.

<!-- eq:A3 -->
$$ E_k(\alpha,\dot{\alpha},\dot{\gamma}) = \tfrac{1}{2}[\dot{\alpha}^{*},\ \dot{\gamma}^{*}]\, \mathcal{M}(\alpha) \begin{bmatrix} \dot{\alpha} \\ \dot{\gamma} \end{bmatrix} $$
- **what:** Quadratic kinetic energy in internal + Euler-angle velocities.
- **symbols:** $\dot\alpha$ - torsion velocity; $\dot\gamma$ - Euler-angle velocity vector.

<!-- eq:A4 -->
$$ \mathcal{M}(\alpha) = \begin{bmatrix} S_0 & V \\ V^{*} & W_0 \end{bmatrix} $$
- **what:** Block partition of mass matrix ($\alpha$ dimension $m$, $\gamma$ dimension $n$).
- **symbols:** $S_0$ ($m\times m$) - torsion block; $W_0$ ($n\times n$) - Euler block; $V$ - coupling block.

<!-- eq:A5 -->
$$ \det\{\mathcal{M}\} = \det\{W_0\}\,\det\{S\}, \qquad S = S_0 - V W_0^{-1} V^{*} $$
- **what:** Schur-complement determinant factorization ($S$ = Schur complement of $W_0$).
- **symbols:** $S,W$ - Schur complements; $S^{-1}=[\mathcal{M}^{-1}]_\alpha$ (eq:A6).

<!-- eq:A9a -->
$$ \int_{-\infty}^{\infty} e^{-x^{*}Ax/2}\, dx = \left[\frac{(2\pi)^p}{\det\{A\}}\right]^{1/2} $$
- **what:** Multivariate Gaussian integral used to integrate out Euler-angle velocities.
- **symbols:** $x$ - $p$-vector; $A$ - $p\times p$ SPD matrix; $p$ - dimension.

<!-- eq:A9b -->
$$ \int_{0}^{\infty} y\, e^{-s y^{2}/2}\, dy = 1/s $$
- **what:** Half-line first-moment Gaussian integral used for the positive-velocity torsion flux.
- **symbols:** $y$ - scalar velocity; $s$ - positive scalar ($=S(\alpha_0)/kT$).
