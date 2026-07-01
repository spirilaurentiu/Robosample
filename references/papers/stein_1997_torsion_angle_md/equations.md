# Equations - Stein, Rice & Brünger 1997 (Torsion-Angle MD)

<!-- eq:1 -->
$$E = E_{\text{chem}} + E_{\text{nmr}}$$
- **what:** Total hybrid energy = chemical (force-field) energy plus NMR-restraint energy.
- **symbols:** E - total potential energy (kcal/mol); E_chem - stereochemistry + nonbonded energy; E_nmr - NMR-restraint energy.

<!-- eq:2 -->
$$E_{\text{nmr}} = w_{\text{NOE}} E_{\text{NOE}} + w_{\text{dihedral}} E_{\text{dihedral}}$$
- **what:** NMR-restraint energy is a weighted sum of NOE-distance and dihedral-angle restraint energies.
- **symbols:** w_NOE - weight on NOE term (dimensionless, e.g. 150); E_NOE - NOE distance-restraint energy (Eq. 6); w_dihedral - weight on dihedral term (e.g. 100); E_dihedral - dihedral-restraint energy (Eq. 8).

<!-- eq:3 -->
$$E_{\text{chem}} = E_{\text{geom}} + w_{\text{vdw}} E_{\text{vdw}}$$
- **what:** Chemical energy = ideal-geometry term (bonds, angles, planarity, chirality) plus weighted van der Waals term.
- **symbols:** E_geom - ideal-geometry energy (bonds/angles/planarity/chirality); w_vdw - weight on vdW term (0.1 -> 1.0 during protocol); E_vdw - van der Waals energy (Eq. 4 or 5).

<!-- eq:4 -->
$$E_{\text{vdw}} = 4\epsilon \left[ \left( \frac{\sigma}{R} \right)^{12} - \left( \frac{\sigma}{R} \right)^{6} \right]$$
- **what:** Standard Lennard-Jones 12-6 potential; used only for final analysis of refined structures, not during the protocol.
- **symbols:** epsilon - LJ well depth for the atom pair (kcal/mol); sigma - LJ distance parameter (Angstrom); R - interatomic distance (Angstrom).

<!-- eq:5 -->
$$E_{\text{vdw}} = \left[ \left(0.8\,\sigma\,\sqrt[6]{2}\right)^2 - R^2 \right]^2$$
- **what:** Purely repulsive quartic ("repel") potential used during structure calculation in place of Lennard-Jones; nonzero only when R is below the reduced contact distance 0.8*sigma*2^(1/6), else 0.
- **symbols:** sigma - LJ distance parameter of the atom pair (Angstrom); R - interatomic distance (Angstrom); 2^(1/6)*sigma - LJ minimum distance; 0.8 - van der Waals radius scale factor.
<!-- CHECK: paper gives the squared form as written; effective potential is E=0 for R > 0.8*sigma*2^(1/6). A vdW force-constant prefactor (w_vdw of Eq.3) multiplies this. -->

<!-- eq:6 -->
$$E_{\text{NOE}} = \begin{cases} \Delta^2 & R < d_{\text{upper}} + 0.5 \\[4pt] a + \dfrac{b}{\Delta^{\text{softexp}}} + \Delta & R > d_{\text{upper}} + 0.5 \end{cases}$$
- **what:** Flat-bottomed parabolic (square-well) NOE distance-restraint energy with a soft asymptote beyond d_upper+0.5 Angstrom; summed over all NOEs.
- **symbols:** Delta - distance violation (Eq. 7, Angstrom); R - model distance between the spin pair (Angstrom); d_upper - upper distance bound (Angstrom); a, b - constants chosen so E_NOE is differentiable at R = d_upper + 0.5; softexp - soft-asymptote exponent.
<!-- CHECK: original OCR wrapped this in "min{...}"; it is a piecewise (flat-bottom harmonic switching to a soft asymptote), reconstructed from the prose "differentiable at R = d_upper + 0.5". -->

<!-- eq:7 -->
$$\Delta = \begin{cases} R - d_{\text{upper}} & d_{\text{upper}} < R \\ 0 & d_{\text{lower}} < R < d_{\text{upper}} \\ d_{\text{lower}} - R & R < d_{\text{lower}} \end{cases}$$
- **what:** Distance violation Delta: zero inside the [d_lower, d_upper] flat well, otherwise the signed distance outside the nearest bound.
- **symbols:** R - model spin-pair distance (Angstrom); d_lower - lower distance bound (Angstrom); d_upper - upper distance bound (Angstrom).

<!-- eq:8 -->
$$E_{\text{dihedral}} = \sum_{\text{dihedrals}} \begin{cases} (\phi - \phi_{\text{upper}})^2 & \phi_{\text{upper}} < \phi \\ 0 & \phi_{\text{lower}} < \phi < \phi_{\text{upper}} \\ (\phi_{\text{lower}} - \phi)^2 & \phi < \phi_{\text{lower}} \end{cases}$$
- **what:** Flat-bottomed harmonic dihedral-angle restraint energy from J-coupling-derived bounds; summed over all restrained dihedrals.
- **symbols:** phi - model dihedral angle (rad or deg); phi_lower - lower dihedral bound; phi_upper - upper dihedral bound.

<!-- eq:9 -->
$$m_i \frac{\partial^2 r_{i,u}}{\partial t^2} = -\frac{\partial E}{\partial r_{i,u}}$$
- **what:** Newton's equation of motion for atom i (Cartesian MD); force is the negative gradient of the hybrid energy.
- **symbols:** m_i - mass of atom i; r_{i,u} - coordinate u (x/y/z) of atom i (Angstrom); t - time (ps); E - hybrid energy (Eq. 1).

<!-- eq:10 -->
$$m_{i} \frac{\partial^{2} r_{i,u}}{\partial t^{2}} = -\frac{\partial E}{\partial r_{i,u}} + \beta_{i} \left(\frac{T_{0}}{T} - 1\right) \vec{v}_{i}$$
- **what:** Berendsen weak-coupling thermostat: adds a velocity-proportional force that drives the system temperature T toward bath temperature T_0 (used for simulated annealing).
- **symbols:** beta_i - coupling force constant for atom i; T_0 - bath (target) temperature (K); T - instantaneous system temperature (K); v_i - velocity vector of atom i (Angstrom/ps).

<!-- eq:11 -->
$$\boldsymbol{\omega}_{j} = \boldsymbol{\omega}_{i} + \hat{\mathbf{h}}_{ij}\,\dot{q}_{ij}$$
- **what:** Angular-velocity constraint for two rigid bodies connected by a torsion bond: body j's angular velocity equals body i's plus the relative rotation rate about the bond axis (single rotational DOF per joint).
- **symbols:** omega_i, omega_j - angular velocities of bodies i and j (rad/ps); h_hat_ij - unit vector along the connecting bond h_ij; q̇_ij - relative torsion-angle rate (rad/ps).

<!-- eq:12 -->
$$\mathbf{r}_{j} = \mathbf{r}_{i} + \mathbf{r}_{ij} = \mathbf{r}_{i} + \mathbf{s}_{ij} + |\mathbf{h}_{ij}|\, \hat{\mathbf{h}}_{ij} - \mathbf{s}_{ji}$$
- **what:** Position of body j's center of mass in terms of body i's, walking i-COM -> bond start (s_ij) -> along bond (|h_ij| h_hat_ij) -> back to j-COM (-s_ji).
- **symbols:** r_i, r_j - centers of mass of bodies i, j (Angstrom, inertial frame); r_ij = r_j - r_i; s_ij - vector from i-COM to the bond endpoint on body i; s_ji - vector from j-COM to the bond endpoint on body j; |h_ij| - fixed bond length; h_hat_ij - unit bond vector.

<!-- eq:13 -->
$$\dot{\mathbf{r}}_{j} = \dot{\mathbf{r}}_{i} - \mathbf{r}_{ij} \times \boldsymbol{\omega}_{i} - (\hat{\mathbf{h}}_{ij} \times \mathbf{s}_{ji})\, \dot{q}_{ij}$$
- **what:** Center-of-mass velocity of body j in terms of body i's linear velocity, body i's angular velocity, and the single joint rate q̇_ij; the recurrence propagated inward/outward over the kinematic tree (differentiate again for accelerations, integrate for positions).
- **symbols:** r_dot_i, r_dot_j - COM velocities of bodies i, j (Angstrom/ps); r_ij - COM offset r_j - r_i; omega_i - angular velocity of body i; h_hat_ij - unit bond vector; s_ji - vector from j-COM to bond endpoint on body j; q̇_ij - relative torsion-angle rate.
<!-- CHECK: compact final form recovered from a corrupted multi-line derivation (see paper.md Derivation block); signs follow r_ij x omega_i and (h_hat x s_ji) q̇. -->
