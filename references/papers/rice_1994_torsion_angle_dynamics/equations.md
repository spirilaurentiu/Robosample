# Equations — Rice & Brünger 1994, Torsion Angle Dynamics

Implementable equations with co-located symbols. Pure proof-algebra steps
(Eqs. 23, 25, 26 — the D'Alembert expansion) are kept in `paper.md` under
"Derivation (not implemented)"; the reusable results are Eqs. 21, 22, 27 and the
recursion (9-18).

---

## Refinement target and energy

<!-- eq:1 -->
$$ E = E_{\text{chem}} + w_{\text{X-ray}}\, E_{\text{X-ray}} $$
- **what:** total refinement objective: empirical chemistry energy plus weighted X-ray agreement term.
- **symbols:** $E$ - total target energy; $E_{\text{chem}}$ - empirical chemical energy; $E_{\text{X-ray}}$ - X-ray data mismatch; $w_{\text{X-ray}}$ - scalar weight balancing the two force sources.

<!-- eq:2 -->
$$ E_{\text{chem}} = \sum_{\text{bonds}} k_{\text{b}}(r - r_0)^2
+ \sum_{\text{angles}} k_{\theta}(\theta - \theta_0)^2
+ \sum_{\text{dihedrals}} k_{\phi} \cos(n\phi + d)
+ \sum_{\text{chiral,planar}} k_{\omega}(\omega - \omega_0)^2
+ \sum_{\text{atom-pairs}} \left( a\,r^{-12} + b\,r^{-6} + c\,r^{-1} \right) $$
- **what:** empirical covalent + nonbonded energy (harmonic bonds/angles, cosine dihedrals, harmonic chiral/planar improper, and 12-6-1 nonbonded). In geometric refinement the nonbonded sum is replaced by a purely repulsive quartic term.
- **symbols:** $k_{\text{b}},k_{\theta},k_{\phi},k_{\omega}$ - force constants; $r,\theta,\phi,\omega$ - bond length, bond angle, dihedral, improper (chiral/planar) coordinate; $r_0,\theta_0,\omega_0$ - equilibrium values; $n$ - dihedral multiplicity; $d$ - dihedral phase; $a,b,c$ - nonbonded coefficients (repulsion, dispersion, electrostatic); $r$ (last sum) - interatomic distance.

<!-- eq:3 -->
$$ E_{\text{X-ray}} = \sum_{\mathbf{h}} \left[ |F_{\text{obs}}(\mathbf{h})| - k\,|F_{\text{calc}}(\mathbf{h})| \right]^2 $$
- **what:** crystallographic residual — sum of squared amplitude differences over reflections; used for high-resolution refinements.
- **symbols:** $\mathbf{h}$ - reciprocal-space Miller index (reflection); $|F_{\text{obs}}|,|F_{\text{calc}}|$ - observed / model structure-factor amplitudes; $k$ - overall scale factor.

<!-- eq:4 -->
$$ E_{\text{X-ray}} = \sum_{\mathbf{h}} \left[ |F_{\text{obs}}(\mathbf{h})| - k\,|F_{\text{calc}}(\mathbf{h})| \right]^2
+ w_{\text{p}} \sum_{\mathbf{h}} f\!\left[\phi_{\text{obs}}(\mathbf{h}) - \phi_{\text{calc}}(\mathbf{h})\right] $$
- **what:** amplitude residual plus a phase-restraint penalty.
- **symbols:** $w_{\text{p}}$ - phase-restraint weight; $f[\cdot]$ - square-well function of width $\arccos(fom(\mathbf{h}))$; $\phi_{\text{obs}},\phi_{\text{calc}}$ - observed / model phases; $fom(\mathbf{h})$ - figure of merit.

<!-- eq:5 -->
$$ E_{\text{X-ray}} = \sum_{\mathbf{h}} fom(\mathbf{h}) \left\{
\left[A_{\text{obs}}(\mathbf{h}) - k A_{\text{calc}}(\mathbf{h})\right]^2
+ \left[B_{\text{obs}}(\mathbf{h}) - k B_{\text{calc}}(\mathbf{h})\right]^2 \right\} $$
- **what:** vector residual — restrains real and imaginary structure-factor parts, weighted by figure of merit; used for medium-resolution refinements.
- **symbols:** $A_{\text{obs}},A_{\text{calc}}$ - real parts of observed / model structure factor; $B_{\text{obs}},B_{\text{calc}}$ - imaginary parts; $fom(\mathbf{h})$ - figure of merit; $k$ - scale.

---

## Molecular dynamics and thermostat

<!-- eq:6 -->
$$ m_i \frac{\partial^2 \mathbf{r}_{i}}{\partial t^2} = -\frac{\partial E}{\partial \mathbf{r}_{i}} $$
- **what:** Newton's equation of motion for atom $i$ (unconstrained Cartesian MD).
- **symbols:** $m_i$ - atom mass; $\mathbf{r}_i$ - atom position ($\mathbb{R}^3$); $E$ - potential/target energy; $t$ - time.

<!-- eq:7 -->
$$ \mathbf{f}_{i}^{\text{frict}} = \beta_{i}\left(\frac{T_{0}}{T} - 1\right)\mathbf{v}_{i} $$
- **what:** Berendsen-style temperature-coupling pseudo-friction force; adds/removes heat by driving velocities toward the bath temperature.
- **symbols:** $\mathbf{f}_i^{\text{frict}}$ - friction force on atom $i$; $\beta_i$ - coupling force constant; $T_0$ - bath (target) temperature; $T$ - instantaneous kinetic temperature; $\mathbf{v}_i$ - atom velocity.

<!-- eq:8 -->
$$ \mathbf{a}(\mathbf{r},\dot{\mathbf{r}}) = \mathbf{M}^{-1}(\mathbf{r})\,\mathbf{Q}(\mathbf{r},\dot{\mathbf{r}}) $$
- **what:** constrained (torsion-angle) system acceleration; unlike unconstrained MD the acceleration depends on velocity, so a 4th-order Runge-Kutta integrator is used (not Verlet).
- **symbols:** $\mathbf{a}$ - system acceleration vector; $\mathbf{M}$ - system inertia (mass) matrix (function of positions); $\mathbf{Q}$ - generalized force vector (function of positions and velocities); $\mathbf{r},\dot{\mathbf{r}}$ - positions, velocities.

---

## Two-body torsion kinematics (recursion primitives)

<!-- eq:9 -->
$$ \boldsymbol{\omega}_j = \boldsymbol{\omega}_i + \hat{\mathbf{h}}_{ij}\,\dot{q}_{ij} $$
- **what:** angular-velocity constraint — child body $j$ shares parent $i$'s angular velocity plus rotation about the shared bond (single torsional DOF).
- **symbols:** $\boldsymbol{\omega}_i,\boldsymbol{\omega}_j$ - lab-frame angular velocities of bodies $i$ (parent), $j$ (child); $\hat{\mathbf{h}}_{ij}=\mathbf{h}_{ij}/|\mathbf{h}_{ij}|$ - unit bond axis; $\dot{q}_{ij}$ - torsion angular velocity.

<!-- eq:10 -->
$$ \mathbf{r}_{j} = \mathbf{r}_{i} + \mathbf{s}_{ij} + |\mathbf{h}_{ij}|\,\hat{\mathbf{h}}_{ij} - \mathbf{s}_{ji} $$
- **what:** child center-of-mass position from parent position through the two attachment points and the bond.
- **symbols:** $\mathbf{r}_i,\mathbf{r}_j$ - COM positions; $\mathbf{s}_{ij}$ - attachment point in body $i$ (from $i$'s COM); $\mathbf{s}_{ji}$ - attachment point in body $j$ (from $j$'s COM); $\mathbf{h}_{ij}$ - bond vector; $|\mathbf{h}_{ij}|$ - fixed bond length.

<!-- eq:11 -->
$$ \dot{\mathbf{r}}_{j} = \dot{\mathbf{r}}_{i}
- \mathbf{r}_{ij} \times \boldsymbol{\omega}_{i}
- (\hat{\mathbf{h}}_{ij} \times \mathbf{s}_{ji})\,\dot{q}_{ij} $$
- **what:** child COM velocity from parent COM velocity, parent angular velocity, and the torsion rate.
- **symbols:** $\dot{\mathbf{r}}_i,\dot{\mathbf{r}}_j$ - COM velocities; $\mathbf{r}_{ij}=\mathbf{r}_j-\mathbf{r}_i$; $\boldsymbol{\omega}_i$ - parent angular velocity; $\mathbf{s}_{ji}$ - child attachment point; $\dot{q}_{ij}$ - torsion rate. <!-- CHECK: sign of the r_ij x omega_i term — paper writes -r_ij x omega_i, equivalent to +omega_i x r_ij -->

<!-- eq:12 -->
$$ \mathbf{Y}_i = \begin{bmatrix} \dot{\mathbf{r}}_i \\ \boldsymbol{\omega}_i \end{bmatrix} $$
- **what:** spatial velocity of body $i$ (stacked linear + angular).
- **symbols:** $\mathbf{Y}_i$ - $6\times1$ spatial velocity; $\dot{\mathbf{r}}_i$ - COM linear velocity ($\mathbb{R}^3$); $\boldsymbol{\omega}_i$ - angular velocity ($\mathbb{R}^3$).

<!-- eq:14 -->
$$ \delta \mathbf{Z}_i = \begin{bmatrix} \delta \mathbf{r}_i \\ \delta \boldsymbol{\pi}_i \end{bmatrix} $$
- **what:** spatial virtual displacement of body $i$ (translational + rotational).
- **symbols:** $\delta \mathbf{Z}_i$ - $6\times1$ virtual displacement; $\delta \mathbf{r}_i$ - translational virtual displacement; $\delta \boldsymbol{\pi}_i$ - rotational virtual displacement.

<!-- eq:15 -->
$$ \mathbf{Y}_{j} = \begin{bmatrix} \mathbf{I} & -\tilde{\mathbf{r}}_{ij} \\ \mathbf{O} & \mathbf{I} \end{bmatrix} \mathbf{Y}_{i}
+ \begin{bmatrix} -\hat{\mathbf{h}}_{ij} \times \mathbf{s}_{ji} \\ \hat{\mathbf{h}}_{ij} \end{bmatrix} \dot{\mathbf{q}}_{ij}
= \mathbf{B}_{ij}^{(1)} \mathbf{Y}_{i} + \mathbf{B}_{ij}^{(2)} \dot{\mathbf{q}}_{ij} $$
- **what:** spatial-velocity propagation parent -> child; $\mathbf{B}^{(1)}$ is the rigid shift (spatial transform), $\mathbf{B}^{(2)}$ is the joint (hinge) map.
- **symbols:** $\mathbf{B}_{ij}^{(1)}$ - $6\times6$ spatial transform (top-right block is $-\tilde{\mathbf{r}}_{ij}$, the skew-symmetric cross-product matrix of $\mathbf{r}_{ij}$); $\mathbf{B}_{ij}^{(2)}$ - $6\times1$ joint axis map; $\mathbf{I}$ - $3\times3$ identity; $\mathbf{O}$ - $3\times3$ zero; $\dot{\mathbf{q}}_{ij}$ - torsion rate ($1\times1$).

<!-- eq:16 -->
$$ \delta \mathbf{Z}_{j} = \mathbf{B}_{ij}^{(1)}\, \delta \mathbf{Z}_{i} + \mathbf{B}_{ij}^{(2)}\, \delta \mathbf{q}_{ij} $$
- **what:** same map applied to virtual displacements (used in the D'Alembert reduction).
- **symbols:** as Eq. 15; $\delta \mathbf{q}_{ij}$ - virtual torsion displacement.

<!-- eq:17 -->
$$ \dot{\mathbf{Y}}_{j} = \mathbf{B}_{ij}^{(1)} \dot{\mathbf{Y}}_{i} + \mathbf{B}_{ij}^{(2)} \ddot{\mathbf{q}}_{ij} + \mathbf{D}_{ij},
\qquad
\mathbf{D}_{ij} = \dot{\mathbf{B}}_{ij}^{(1)} \mathbf{Y}_{i} + \dot{\mathbf{B}}_{ij}^{(2)} \dot{\mathbf{q}}_{ij} $$
- **what:** spatial-acceleration propagation parent -> child; $\mathbf{D}_{ij}$ is the velocity-dependent (Coriolis/gyroscopic) bias term.
- **symbols:** $\dot{\mathbf{Y}}_i,\dot{\mathbf{Y}}_j$ - spatial accelerations; $\ddot{\mathbf{q}}_{ij}$ - torsion angular acceleration; $\mathbf{D}_{ij}$ - bias acceleration; $\dot{\mathbf{B}}^{(1)},\dot{\mathbf{B}}^{(2)}$ - time derivatives of the transform/joint maps. <!-- CHECK: OCR gave D_ij = B^(1)_dot Y_i + B^(2)_dot q_dot; reconstructed as time-derivatives of B^(1),B^(2) -->

---

## D'Alembert dynamics and articulated-body solve

<!-- eq:18 -->
$$ \sum_{\text{bodies } k} \delta \mathbf{Z}_k^{T}\left(\mathbf{M}_k \dot{\mathbf{Y}}_k - \mathbf{Q}_k\right) = 0 $$
- **what:** D'Alembert principle — virtual work of applied forces vanishes over all constraint-consistent displacements; the governing equation of constrained motion.
- **symbols:** $\delta \mathbf{Z}_k$ - constraint-consistent virtual displacement of body $k$; $\mathbf{M}_k$ - $6\times6$ spatial inertia; $\dot{\mathbf{Y}}_k$ - spatial acceleration; $\mathbf{Q}_k$ - generalized force.

<!-- eq:21 -->
$$ \mathbf{M}_{i} = \begin{bmatrix}
m_{i} & 0 & 0 & 0 & 0 & 0 \\
0 & m_{i} & 0 & 0 & 0 & 0 \\
0 & 0 & m_{i} & 0 & 0 & 0 \\
0 & 0 & 0 & I_{xx} & I_{xy} & I_{xz} \\
0 & 0 & 0 & I_{xy} & I_{yy} & I_{yz} \\
0 & 0 & 0 & I_{xz} & I_{yz} & I_{zz}
\end{bmatrix} $$
- **what:** spatial (6x6) inertia matrix of a rigid body about its center of mass — block-diagonal mass and inertia tensor (no coupling because it is taken about the COM).
- **symbols:** $m_i$ - body mass; $I_{xx},\dots,I_{zz}$ - components of the $3\times3$ inertia tensor $\mathbf{I}$ about the COM.

<!-- eq:22 -->
$$ \mathbf{Q}_{i} = \begin{bmatrix}
F_{i}^{x} \\ F_{i}^{y} \\ F_{i}^{z} \\
N_{i}^{x} - (\boldsymbol{\omega}_{i} \times \mathbf{I}\boldsymbol{\omega}_{i})_{x} \\
N_{i}^{y} - (\boldsymbol{\omega}_{i} \times \mathbf{I}\boldsymbol{\omega}_{i})_{y} \\
N_{i}^{z} - (\boldsymbol{\omega}_{i} \times \mathbf{I}\boldsymbol{\omega}_{i})_{z}
\end{bmatrix} $$
- **what:** spatial generalized force on body $i$ — net force stacked with net torque minus the gyroscopic (Euler) term, so that $\mathbf{M}_i\dot{\mathbf{Y}}_i = \mathbf{Q}_i$ holds with only $\mathbf{I}\dot{\boldsymbol{\omega}}_i$ on the LHS.
- **symbols:** $F_i^{x,y,z}$ - net force components on body $i$; $N_i^{x,y,z}$ - net lab-frame torque about the COM; $\boldsymbol{\omega}_i \times \mathbf{I}\boldsymbol{\omega}_i$ - gyroscopic/Euler correction; $\mathbf{I}$ - $3\times3$ inertia tensor.

<!-- eq:27 -->
$$ \ddot{\mathbf{q}}_{ij} = -\left(\mathbf{B}_{ij}^{(2)T}\mathbf{M}_{j}\mathbf{B}_{ij}^{(2)}\right)^{-1}
\left\{\mathbf{B}_{ij}^{(2)T}\mathbf{M}_{j}\mathbf{B}_{ij}^{(1)}\dot{\mathbf{Y}}_{i}
+ \mathbf{B}_{ij}^{(2)T}\left[\mathbf{M}_{j}\mathbf{D}_{ij} - \mathbf{Q}_{j}\right]\right\} $$
- **what:** relative (joint) acceleration of the tip body in terms of the inboard body's acceleration — the inward-reduction step of the recursive articulated-body solve. The bracketed inverse is the (here scalar) articulated joint inertia; setting the $\delta\mathbf{q}$ coefficient to zero gives this.
- **symbols:** $\ddot{\mathbf{q}}_{ij}$ - torsion angular acceleration; $\mathbf{B}_{ij}^{(1)},\mathbf{B}_{ij}^{(2)}$ - transform/joint maps (Eq. 15); $\mathbf{M}_j$ - child spatial inertia; $\mathbf{D}_{ij}$ - bias acceleration (Eq. 17); $\mathbf{Q}_j$ - child generalized force; $\dot{\mathbf{Y}}_i$ - inboard spatial acceleration (solved after full inward reduction).

---

## Derivations note

Eqs. 23, 25, 26 are the intermediate D'Alembert expansion for the two-body case
(substituting Eq. 24's constraint relations and factoring the $\delta\mathbf{Z}_i$
and $\delta\mathbf{q}_{ij}$ terms). They are algebra, not standalone
implementables; the extractable results are the inertia/force blocks (Eqs. 21,
22) and the joint-acceleration solve (Eq. 27), together with the velocity/
acceleration propagation maps (Eqs. 15-17). See `paper.md` for the full
derivation.
