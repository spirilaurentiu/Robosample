# Notation — Rice & Brünger 1994, Torsion Angle Dynamics

Conventions: rigid bodies are indexed with $i$ = inboard (parent, toward root),
$j$ = outboard (child, toward tip). Spatial (6-vector) quantities stack
translation on top of rotation. Inertia is taken about each body's center of
mass, so the spatial inertia $\mathbf{M}_i$ is block-diagonal. Angular
quantities are in the lab (inertial) frame unless noted. Integration is 4th-order
Runge-Kutta (not Verlet) because acceleration depends on velocity.

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| $\mathbf{r}_i$ | center-of-mass position of body $i$ | Angstrom, $\mathbb{R}^3$ | inertial frame |
| $\dot{\mathbf{r}}_i$ | COM linear velocity of body $i$ | Angstrom/ps, $\mathbb{R}^3$ | |
| $\boldsymbol{\omega}_i$ | angular velocity of body $i$ | rad/ps, $\mathbb{R}^3$ | lab frame |
| $\mathbf{r}_{ij}$ | $\mathbf{r}_j - \mathbf{r}_i$ | Angstrom, $\mathbb{R}^3$ | child minus parent |
| $\tilde{\mathbf{r}}_{ij}$ | skew-symmetric cross-product matrix of $\mathbf{r}_{ij}$ | $3\times3$ | $\tilde{\mathbf{r}}\,\mathbf{x}=\mathbf{r}\times\mathbf{x}$ |
| $\mathbf{s}_{ij}$ | attachment point of bond in body $i$ (from $i$'s COM) | Angstrom, $\mathbb{R}^3$ | body-fixed in $i$ |
| $\mathbf{s}_{ji}$ | attachment point of bond in body $j$ (from $j$'s COM) | Angstrom, $\mathbb{R}^3$ | body-fixed in $j$ |
| $\mathbf{h}_{ij}$ | bond vector connecting bodies $i,j$ | Angstrom, $\mathbb{R}^3$ | fixed length $|\mathbf{h}_{ij}|$; body-fixed in $i$ |
| $\hat{\mathbf{h}}_{ij}$ | unit bond axis $\mathbf{h}_{ij}/|\mathbf{h}_{ij}|$ | dimensionless, $\mathbb{R}^3$ | the single torsional hinge axis |
| $q_{ij}$ | relative torsion angle about the bond | rad, scalar | the only inter-body DOF |
| $\dot{q}_{ij},\ddot{q}_{ij}$ | torsion angular velocity / acceleration | rad/ps, rad/ps^2 | |
| $\mathbf{q}_{ij}$ | $1\times1$ vector holding $q_{ij}$ | $1\times1$ | vectorized for matrix ops |
| $\mathbf{Y}_i$ | spatial velocity $[\dot{\mathbf{r}}_i;\boldsymbol{\omega}_i]$ | $6\times1$ | translation over rotation |
| $\dot{\mathbf{Y}}_i$ | spatial acceleration | $6\times1$ | |
| $\delta\mathbf{Z}_i$ | spatial virtual displacement $[\delta\mathbf{r}_i;\delta\boldsymbol{\pi}_i]$ | $6\times1$ | translation over rotation |
| $\mathbf{B}_{ij}^{(1)}$ | spatial transform (rigid shift) parent->child | $6\times6$ | top-right block $-\tilde{\mathbf{r}}_{ij}$ |
| $\mathbf{B}_{ij}^{(2)}$ | joint (hinge) map for the torsion DOF | $6\times1$ | $[-\hat{\mathbf{h}}_{ij}\times\mathbf{s}_{ji};\ \hat{\mathbf{h}}_{ij}]$ |
| $\mathbf{D}_{ij}$ | velocity-dependent bias (Coriolis/gyroscopic) acceleration | $6\times1$ | $\dot{\mathbf{B}}^{(1)}\mathbf{Y}_i+\dot{\mathbf{B}}^{(2)}\dot{q}_{ij}$ |
| $\mathbf{M}_i$ | spatial inertia of body $i$ | $6\times6$ | about COM; block-diag mass, inertia |
| $m_i$ | mass of body $i$ (or atom $i$ in Cartesian MD) | amu | |
| $\mathbf{I}$ | inertia tensor $I_{xx}\dots I_{zz}$ | amu·Angstrom^2, $3\times3$ | symmetric, about COM |
| $\mathbf{Q}_i$ | spatial generalized force | $6\times1$ | $[\,\mathbf{F}_i;\ \mathbf{N}_i-\boldsymbol{\omega}_i\times\mathbf{I}\boldsymbol{\omega}_i\,]$ |
| $\mathbf{F}_i$ | net force on body $i$ | force units | reduced from per-atom forces |
| $\mathbf{N}_i$ | net torque on body $i$ about COM | torque units | lab frame |
| $\mathbf{a}$ | system acceleration vector | $\mathbf{M}^{-1}\mathbf{Q}$ | depends on $\mathbf{r},\dot{\mathbf{r}}$ |
| $E$ | total refinement target energy | energy units | Eq. 1 |
| $E_{\text{chem}}$ | empirical chemical energy | energy units | $E^{*}_{\text{chem}}$ = repulsive-only variant |
| $E_{\text{X-ray}}$ | X-ray data mismatch term | energy units | residual (Eq. 3) or vector residual (Eq. 5) |
| $w_{\text{X-ray}}$ | X-ray weight | scalar | balances forces |
| $w_{\text{p}}$ | phase-restraint weight | scalar | Eq. 4 |
| $T_0$ | bath (target) temperature | K | thermostat setpoint |
| $T$ | instantaneous kinetic temperature | K | from KE |
| $\beta_i$ | temperature-coupling force constant | | Berendsen coupling |
| $\mathbf{h}$ | reflection (Miller index) | reciprocal space | X-ray sums run over $\mathbf{h}$ |
| $F_{\text{obs}},F_{\text{calc}}$ | observed / calculated structure factor amplitude | | $k$ = overall scale |
| $A,B$ | real / imaginary parts of structure factor | | Eq. 5 |
| $fom(\mathbf{h})$ | figure of merit | $[0,1]$ | phase reliability weight |
| $\phi_{\text{obs}},\phi_{\text{calc}}$ | observed / model phase | rad | |
| $R_{\text{free}}$ | free R value (cross-validation) | dimensionless | model-quality selector |
| $\Delta\varphi_{\text{calc}}$ | fom-weighted mean phase error vs. crystal structure | degrees | convergence measure |

## Sign / frame notes

- Eq. 22 subtracts the gyroscopic term $\boldsymbol{\omega}_i\times\mathbf{I}\boldsymbol{\omega}_i$
  from the torque so that the mass-matrix form $\mathbf{M}_i\dot{\mathbf{Y}}_i=\mathbf{Q}_i$
  uses only $\mathbf{I}\dot{\boldsymbol{\omega}}_i$ (Euler's rigid-body equation:
  $\mathbf{N}_i=\mathbf{I}\dot{\boldsymbol{\omega}}_i+\boldsymbol{\omega}_i\times\mathbf{I}\boldsymbol{\omega}_i$).
- Cross-product with a vector on the right, $\mathbf{v}\times\boldsymbol{\omega}$,
  equals $-\tilde{\mathbf{v}}\,\boldsymbol{\omega}$; watch this sign against
  Robosample's spatial-operator conventions (Jain 1993 / Simbody).
- The $6\times6$ inertia is block-diagonal only because it is expressed about the
  center of mass; a body-frame-at-joint convention would introduce off-diagonal
  coupling.
