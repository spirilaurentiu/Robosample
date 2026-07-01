# Equations - Gō & Scheraga 1976 (rigid vs flexible partition functions)

Reduced conventions: $\beta = 1/kT$, $k$ = Boltzmann constant, $T$ = absolute temperature,
$\hbar$ = reduced Planck constant. "Hard variables" = bond lengths + bond angles; "soft
variables" = dihedral angles + 6 external (rigid-body) variables.

<!-- eq:1 -->
$$ Z = (\text{constant}) \int \exp[-\beta F(Q)]\, dQ $$
- **what:** Classical *flexible*-model partition function: a flat integral of the Boltzmann factor over dihedral angles only, with NO Jacobian/metric weight. This is the form the paper concludes is the better approximation.
- **symbols:** $Z$ - configurational partition function; $Q=(q_1,\dots,q_m)$ - soft variables (dihedral + external), integration domain; $F(Q)$ - conformational (free) energy / potential of mean force; $\beta=1/kT$.

<!-- eq:2 -->
$$ Z = (\text{constant}) \int \left[ \frac{1}{\det \mathbf{G}} \right]^{1/2} \exp[-\beta F(Q)]\, dQ $$
- **what:** Classical *rigid*-model partition function: same integral but weighted by $(\det \mathbf{G})^{-1/2}$, the metric-tensor (Fixman-type) factor. Same as eq 13. The det G factor is Q-dependent and cannot be pulled outside the integral.
- **symbols:** $\mathbf{G}$ - $m\times m$ kinetic-energy matrix in soft-variable momenta, $\mathbf{G}=(\mathbf{H}^0)^{-1}$; other symbols as in eq 1.

<!-- eq:3 -->
$$ F(Q,Q') = F_0(Q) + \frac{1}{2} \sum_{i,j=1}^{l} f_{ij}''\,(q_i' - q_{i0}')(q_j' - q_{j0}') $$
- **what:** Conformational energy of the flexible model: minimum energy $F_0(Q)$ over hard variables plus a harmonic (quadratic) restoring term in the hard-variable displacements.
- **symbols:** $Q'=(q_1',\dots,q_l')$ - hard variables; $q_{i0}'$ - strain-free (minimum-energy) hard-variable values, approximated as Q-independent; $f_{ij}''$ - hard-variable force constants (Q-independent); $F_0(Q)$ - min energy for fixed $Q$; $l$ - number of hard variables.

<!-- eq:4 -->
$$ K_f = \frac{1}{2} \sum_{k=1}^{n} \sum_{\alpha=1}^{3} m_k\, \dot{x}_{k\alpha}^2 $$
- **what:** Kinetic energy of the flexible model in Cartesian velocities.
- **symbols:** $K_f$ - kinetic energy (flexible model); $n$ - number of atoms; $m_k$ - mass of atom $k$; $\dot{x}_{k\alpha}$ - Cartesian velocity component $\alpha\in\{1,2,3\}$ of atom $k$; note $3n=l+m$.

<!-- eq:5 -->
$$ K_f = \frac{1}{2}\dot{Q}^{\dagger} \mathcal{H}\dot{Q} = \frac{1}{2}(\dot{Q}^{+}, \dot{Q}'^{+}) \begin{pmatrix} \mathbf{H}^0 & \mathbf{H}' \\ \mathbf{H}'^{+} & \mathbf{H}'' \end{pmatrix} \begin{pmatrix} \dot{Q} \\ \dot{Q}' \end{pmatrix} $$
- **what:** Kinetic energy of the flexible model in internal-coordinate velocities; $\mathcal{H}$ is the full $3n\times3n$ mass-metric block matrix partitioned into soft/hard blocks.
- **symbols:** $\dot{Q}$ - $m$-vector of soft velocities; $\dot{Q}'$ - $l$-vector of hard velocities; $\mathbf{H}^0$ ($m\times m$), $\mathbf{H}'$ ($m\times l$), $\mathbf{H}''$ ($l\times l$) - blocks of the mass-metric matrix $\mathcal{H}$; superscript $+$ = transpose.

<!-- eq:6 -->
$$ K_r = \frac{1}{2}\dot{Q}^{+} \mathbf{H}^0 \dot{Q} $$
- **what:** Kinetic energy of the rigid model: set hard-variable velocities $\dot{Q}'=0$ in eq 5, leaving only the soft-soft block $\mathbf{H}^0$.
- **symbols:** $K_r$ - kinetic energy (rigid model); $\mathbf{H}^0$ - soft-soft block of the mass metric.

<!-- eq:9 -->
$$ K_r = \frac{1}{2}P_r^{+}\, \mathbf{G}\, P_r $$
- **what:** Rigid-model kinetic energy in terms of soft momenta; $\mathbf{G}=(\mathbf{H}^0)^{-1}$ is the inverse mass metric.
- **symbols:** $P_r = \partial K_r/\partial \dot{Q} = \mathbf{H}^0 \dot{Q}$ - generalized momentum conjugate to soft variable $Q$; $\mathbf{G}$ - $m\times m$ matrix $=(\mathbf{H}^0)^{-1}$.

<!-- eq:10 -->
$$ \mathbf{G} \left[ \equiv (\mathbf{H}^0)^{-1} \right] = \mathbf{G}^0 - \mathbf{G}'\, \mathbf{G}''^{-1}\, \mathbf{G}'^{+} $$
- **what:** Schur-complement expression for $\mathbf{G}$: the soft-block of the inverse of the full mass metric $\mathcal{H}$ (NOT simply the inverse of $\mathbf{H}^0$'s own block unless coupling vanishes).
- **symbols:** $\mathbf{G}^0$ ($m\times m$), $\mathbf{G}'$ ($m\times l$), $\mathbf{G}''$ ($l\times l$) - blocks of $\mathcal{G}=\mathcal{H}^{-1}$ (eq 11).

<!-- eq:11 -->
$$ \mathcal{G} = \begin{pmatrix} \mathbf{G}^0 & \mathbf{G}' \\ \mathbf{G}'^{+} & \mathbf{G}'' \end{pmatrix} = \mathcal{H}^{-1} $$
- **what:** Full $3n\times3n$ inverse mass-metric matrix, block-partitioned into soft/hard.
- **symbols:** $\mathcal{G}$ - inverse of $\mathcal{H}$ from eq 5.

<!-- eq:12 -->
$$ Z_r = \left(\frac{1}{2\pi\hbar}\right)^m \int \exp\left[-\beta \{K_r + F_0(Q)\}\right]\, dP_r\, dQ $$
- **what:** Phase-space definition of the rigid-model partition function (integrate over soft momenta and soft coordinates).
- **symbols:** $m$ - number of soft variables; $dP_r$ - momentum measure; $F_0(Q)$ - conformational energy of rigid model.

<!-- eq:13 -->
$$ Z_r = \left(\frac{kT}{2\pi\hbar^2}\right)^{m/2} (8\pi^2 V) \int \left[\frac{1}{\det \mathbf{G}}\right]^{1/2} \exp[-\beta F_0(Q)]\, dQ $$
- **what:** Rigid-model partition function after Gaussian momentum integration and integrating out the 6 external variables (giving $8\pi^2 V$). Concrete form of eq 2; carries the $(\det\mathbf{G})^{-1/2}$ weight.
- **symbols:** $V$ - system volume (from translational integration); $8\pi^2$ - from rotational (Euler-angle) integration; $\det\mathbf{G}$ - determinant of the $m\times m$ inverse mass metric.

<!-- eq:8 -->
$$ Z_f = \left[ \prod_{k=1}^{n}\left( \frac{m_k kT}{2\pi\hbar^2} \right)^{3/2} \right]\left[ \frac{(2\pi kT)^l}{\det \mathbf{F}''} \right]^{1/2}\left[ 8\pi^2 V\, D(Q_0') \right]\int \exp[-\beta F_0(Q)]\, dQ $$
- **what:** Flexible-model partition function: momentum + hard-variable Gaussian integrations produce Q-independent prefactors, leaving the flat integral of eq 1. The det-F'' and Jacobian prefactors are constants pulled outside.
- **symbols:** $\mathbf{F}''$ - $l\times l$ matrix, $(i,j)$ element $f_{ij}''$; $D(Q_0')$ - Jacobian evaluated at strain-free hard values; $l$ - number of hard variables; product over atom masses from Cartesian momentum integration.

<!-- eq:14 -->
$$ D = \left[ \det \mathbf{G}\; \det \mathbf{G}'' \left( \prod_{k=1}^{n} m_k \right)^3 \right]^{-1/2} $$
- **what:** Jacobian of the Cartesian -> internal-coordinate transformation, expressed via the metric determinants. Key link between flexible and rigid forms.
- **symbols:** $D$ - Jacobian $|\partial x / \partial (Q,Q',\text{ext})|$; $\det\mathbf{G}$, $\det\mathbf{G}''$ - metric determinants; $\prod_k m_k$ - product of atomic masses (cubed for 3 Cartesian components).

<!-- eq:15 -->
$$ \det\!\left( \mathbf{F}''\mathbf{G}'' \right) = \prod_{i=1}^{l} (2\pi\nu_i)^2 $$
- **what:** Wilson GF-matrix relation: product of squared angular vibrational frequencies of the hard variables equals det of the force-constant times inverse-mass metric.
- **symbols:** $\nu_i$ - vibrational frequency (Hz) of $i$th hard-variable normal mode; $\mathbf{F}''$ - hard-variable force constants; $\mathbf{G}''$ - hard-block of inverse mass metric.

<!-- eq:16 -->
$$ Z_f = \left(\frac{kT}{2\pi\hbar^2}\right)^{m/2} (8\pi^2 V) \int \left[\prod_{i=1}^{l} \frac{kT}{2\pi\hbar\nu_i}\right]\left[\frac{1}{\det \mathbf{G}}\right]^{1/2} \exp[-\beta F_0(Q)]\, dQ $$
- **what:** Alternate exact form of $Z_f$ obtained by substituting eq 14 and eq 15 into eq 8. Reveals: flexible model = rigid model $(\det\mathbf{G})^{-1/2}$ TIMES a product of classical harmonic-oscillator partition functions with conformation-dependent frequencies $\nu_i$. If those frequencies were Q-independent the two models coincide.
- **symbols:** $kT/2\pi\hbar\nu_i$ - classical harmonic-oscillator partition function of mode $i$; other symbols as in eq 13.

<!-- eq:17 -->
$$ Z_{\rm QM} = \left(\frac{kT}{2\pi\hbar^2}\right)^{m/2} (8\pi^2 V) \int \left[ \prod_{i=1}^{l}\left( 2\sinh\frac{\pi\hbar\nu_i}{kT} \right)^{-1} \right]\left[ \frac{1}{\det \mathbf{G}} \right]^{1/2}\exp[-\beta F_0(Q)]\, dQ $$
- **what:** Quantum-mechanically-correct partition function: replace each classical oscillator factor $kT/2\pi\hbar\nu_i = (2\pi\hbar\nu_i/kT)^{-1}$ in eq 16 by the quantum factor $[2\sinh(\pi\hbar\nu_i/kT)]^{-1}$. Both $Z_r$ and $Z_f$ are approximations to this.
- **symbols:** $[2\sinh(\pi\hbar\nu_i/kT)]^{-1}$ - QM harmonic-oscillator partition function (including zero-point energy $\pi\hbar\nu_i = \tfrac12 h\nu_i$).

<!-- eq:18 -->
$$ \Gamma_r = \prod_{i=1}^{l}\left( \frac{1}{2\sinh (x_i/2)} \right) $$
- **what:** Ratio of the QM vibrational partition function to the corresponding factor (unity) in the rigid model $Z_r$. Its conformation dependence measures the rigid-model error.
- **symbols:** $x_i = 2\pi\hbar\nu_i/kT = h\nu_i/kT$ - dimensionless mode variable; $l$ - number of hard modes.

<!-- eq:19 -->
$$ \Gamma_f = \prod_{i=1}^{l}\left( \frac{x_i}{2\sinh (x_i/2)} \right) $$
- **what:** Ratio of the QM vibrational partition function to the corresponding classical factor in the flexible model $Z_f$. Its conformation dependence measures the flexible-model error; equals Flory's $\Gamma$.
- **symbols:** $x_i = 2\pi\hbar\nu_i/kT$ as above.

<!-- eq:20 -->
$$ \ln\!\left( \Gamma_r / \Gamma_{r_0} \right) = \sum_{i=1}^{l} g_r(x_{i0})\, \Delta x_i $$
- **what:** First-order (log) change of the rigid-model correction factor between a conformation and a reference conformation.
- **symbols:** $\Gamma_{r_0}$ - reference-conformation value; $\Delta x_i = x_i - x_{i0}$ - shift of mode $i$; $g_r$ - log-derivative from eq 22.

<!-- eq:21 -->
$$ \ln\!\left( \Gamma_f / \Gamma_{f_0} \right) = \sum_{i=1}^{l} g_f(x_{i0})\, \Delta x_i $$
- **what:** First-order (log) change of the flexible-model correction factor between a conformation and reference.
- **symbols:** $g_f$ - log-derivative from eq 23; other symbols as eq 20.

<!-- eq:22 -->
$$ g_r(x_{i0}) = \frac{d}{dx}\ln\!\left( \frac{1}{2\sinh (x/2)} \right)\Bigg|_{x=x_{i0}} = -\frac{1}{2}\coth\frac{x_{i0}}{2} $$
- **what:** Sensitivity of the rigid-model log-weight per mode. Magnitude $\to \tfrac12$ as $x\to\infty$; diverges as $x\to0$ (rigid model fails for softened low-frequency modes).
- **symbols:** $x_{i0}$ - reference value of $x_i$; $\coth$ - hyperbolic cotangent.

<!-- eq:23 -->
$$ g_f(x_{i0}) = \frac{d}{dx}\ln\!\left( \frac{x}{2\sinh (x/2)} \right)\Bigg|_{x=x_{i0}} = \frac{1}{x_{i0}} - \frac{1}{2}\coth\frac{x_{i0}}{2} $$
- **what:** Sensitivity of the flexible-model log-weight per mode. Magnitude $\to\tfrac12$ as $x\to\infty$ but $\to0$ as $x\to0$, so the flexible model is well-behaved for softened modes. Always $|g_f|<|g_r|$, hence flexible model is the better approximation.
- **symbols:** as eq 22.

## Derivations (not implemented)

The following are proof/algebra steps, not standalone implementable formulas:

- **eq:C-3** (Kirkwood diffusion form): with all masses equal $m_0$, $\mathbf{H}^0 = m_0\, g$ where $g$ is the pure geometric metric tensor of the constrained space, giving $Z_r \propto \int (\det g)^{1/2}\exp[-\beta F_0]\,dQ$. Note the sign flip: $(\det g)^{+1/2}$ because $g=\mathbf{G}^{-1}/m_0$, i.e. $\det\mathbf{G}\propto 1/\det g$.
- **Appendix A** (eqs A-1..A-8): explicit mass-weighted expressions for the metric blocks $\mathbf{H}^0,\mathbf{H}',\mathbf{H}''$ and $\mathbf{G}^0,\mathbf{G}',\mathbf{G}''$ from Cartesian-to-internal Jacobians; heavily OCR-corrupted in the raw file.
- **Appendix B** (eqs B-2..B-15): proof that the Cartesian->internal Jacobian $D = \prod_{k}d_k^{-2} \cdot \prod_k \sin\tau_k$ (bond lengths $d_k$, bond angles $\tau_k$) is independent of dihedral angles; used to derive eq 8. OCR-corrupted.
- **Appendix D** (eqs D-1..D-5): rewriting the flexible Hamiltonian as decoupled hard-oscillator + soft-variable parts, justifying replacement of classical by QM oscillator factors (eq 16 -> eq 17). OCR-corrupted.

<!-- CHECK: Appendix A/B/D display equations in raw.md are severely OCR-garbled (Greek/subscripts lost); the clean forms above are reconstructed from the main-text prose and standard Wilson GF-matrix theory. -->
