# Equations - Forrest & Suter 1994, Generalized coordinate hybrid Monte Carlo

Reduced units: bead mass `mb = 1`; `β = 1/(kB T)`. Tildes denote fictitious
(computational) variables. Angles `φ_k`: `k=1,2,3` are Euler angles, `4 ≤ k ≤ Nb`
are torsional angles.

<!-- eq:1 -->
$$\mathscr{E}(\{\boldsymbol{r}_{1}^{(c)}, \{\phi_{k}^{(c)}\}\}) = \sum_{c=1}^{N_{c}} \sum_{k=4}^{N_{b}} V_{\phi}(\phi_{k}^{(c)}) + \sum_{c=1}^{N_{c}} \sum_{i=5}^{N_{b}} \sum_{j=1}^{i-4} V_{LJ}(|\boldsymbol{r}_{i}^{(c)} - \boldsymbol{r}_{j}^{(c)}|) + \sum_{c=2}^{N_{c}} \sum_{i=1}^{N_{b}} \sum_{c'=1}^{c-1} \sum_{j=1}^{N_{b}} V_{LJ}(|\boldsymbol{r}_{i}^{(c)} - \boldsymbol{r}_{j}^{(c')}|)$$
- **what:** Total potential energy of the polybead melt: torsional + intramolecular LJ (beads ≥4 apart) + intermolecular LJ.
- **symbols:** E - total potential energy; Nc - number of chains; Nb - beads per chain; r_i^(c) - position of bead i of chain c; φ_k^(c) - k-th generalized angle of chain c; V_φ - torsional potential; V_LJ - Lennard-Jones potential.

<!-- eq:2 -->
$$V_{\phi}(\phi) = C \sum_{n=0}^{5} a_n \cos^n(\phi)$$
- **what:** Ryckaert-Bellemans torsional potential (5th-order cosine polynomial).
- **symbols:** C = 9.0 kJ/mol; a0=1, a1=1.31, a2=-1.414, a3=-0.3297, a4=2.828, a5=-3.3943 (dimensionless); φ - dihedral angle (rad).

<!-- eq:3 -->
$$V_{\rm LJ}(r_{ij}) = 4\varepsilon \left[ \left( \frac{\sigma}{r_{ij}} \right)^{12} - \left( \frac{\sigma}{r_{ij}} \right)^{6} \right]$$
- **what:** 12-6 Lennard-Jones potential (same params for all beads).
- **symbols:** ε = 410 J/mol (well depth); σ = 3.94 Å (diameter); r_ij - inter-bead distance (Å).

<!-- eq:4 -->
$$\langle O \rangle \equiv \langle O(\{x\}) \rangle = \frac{1}{C} \int dp \exp\left[-\frac{p^2}{2m_b k_B T}\right] \int dx \, O(\{x\}) \exp\left[-\frac{\mathscr{E}(\{x\})}{k_B T}\right]$$
- **what:** Full phase-space canonical average in Cartesian coordinates; momentum Gaussian factors out.
- **symbols:** O - observable; {x} - 3 Nc Nb Cartesian coords; {p} - conjugate momenta; mb - bead mass; kB - Boltzmann constant; T - temperature; C - normalization constant.

<!-- eq:5 -->
$$\langle O(\{x\}) \rangle = \frac{1}{Z} \int dx \ O(\{x\}) \exp \left[ -\frac{\mathscr{E}(\{x\})}{k_{\rm B} T} \right]$$
- **what:** Configurational canonical average after integrating out Cartesian momenta.
- **symbols:** Z = ∫ dx exp[-E({x})/kB T] - configurational partition function.

<!-- eq:6 -->
$$\langle O \rangle = \frac{1}{Z_q} \int dq \ D(\{q\})O(\{q\}) \exp \left[ -\frac{\mathscr{E}(\{q\})}{k_B T} \right]$$
- **what:** Average in generalized (rigid-model) coordinates; D({q}) is the coordinate-transform Jacobian plus Fixman correction.
- **symbols:** {q} = {r_1^(c), φ_k^(c)} generalized coords; D({q}) - Jacobian × Fixman weight (dimensionless); Zq = ∫ dq D({q}) exp[-E({q})/kB T].

<!-- eq:7 -->
$$P_{\text{acc}} = \min \left\{ 1, \frac{P_{\text{S}}[\{q_n\} \to \{q_0\}] D(\{q_n\}) \exp(-\mathscr{E}(\{q_n\})/k_{\text{B}}T)}{P_{\text{S}}[\{q_0\} \to \{q_n\}] D(\{q_0\}) \exp(-\mathscr{E}(\{q_0\})/k_{\text{B}}T)} \right\}$$
- **what:** General Metropolis acceptance for a move in {q} space including the D-factor; satisfies detailed balance.
- **symbols:** PS[a→b] - proposal probability a to b; {q0} - old config; {qn} - new candidate config.

<!-- eq:8 -->
$$\tilde{K}(\{p_1^{(c)}\}, \{\tilde{\pi}_k^{(c)}\}) = \sum_{c=1}^{N_c} \frac{(p_1^{(c)})^2}{2m_c} + \sum_{c=1}^{N_c} \sum_{k=1}^{N_b} \frac{(\tilde{\pi}_k^{(c)})^2}{2\tilde{I}_k^{(c)}}$$
- **what:** Fictitious Cartesian-like kinetic energy: translational (chain origins) + diagonal angular terms with effective moments of inertia.
- **symbols:** K̃ - fictitious kinetic energy; p_1^(c) - Cartesian momentum of chain-c origin; mc - mass of chain c; π̃_k^(c) - fictitious momentum conjugate to angle φ_k^(c); Ĩ_k^(c) - effective moment of inertia (freely chosen constant).

<!-- eq:9 -->
$$\frac{\mathrm{d}\tilde{\pi}_{k}^{(c)}}{\mathrm{d}\tilde{t}} = -\frac{\partial \mathscr{E}}{\partial \phi_{k}^{(c)}}, \qquad \frac{\mathrm{d}\phi_{k}^{(c)}}{\mathrm{d}\tilde{t}} = \frac{\tilde{\pi}_{k}^{(c)}}{\tilde{I}_{k}^{(c)}}, \qquad \frac{\mathrm{d}\boldsymbol{p}_{1}^{(c)}}{\mathrm{d}\tilde{t}} = -\frac{\partial \mathscr{E}}{\partial \boldsymbol{r}_{1}^{(c)}}, \qquad \frac{\mathrm{d}\boldsymbol{r}_{1}^{(c)}}{\mathrm{d}\tilde{t}} = \frac{\boldsymbol{p}_{1}^{(c)}}{m_{c}}$$
- **what:** Decoupled Hamiltonian equations of motion from K̃ (eq 8); integrate with leap-frog. Torque = -∂E/∂φ_k, angular velocity = π̃/Ĩ.
- **symbols:** t̃ - fictitious (computer) time; ∂E/∂φ_k^(c) - generalized force/torque on angle k; ∂E/∂r_1^(c) - Cartesian force on chain origin.

<!-- eq:10 -->
$$\Pr\left(\left\{\boldsymbol{p}_{1}, \tilde{\boldsymbol{\pi}}_{k}\right\}\right) \propto \prod_{c=1}^{N_{c}} \left[\exp\left(-\frac{\beta(\boldsymbol{p}_{1}^{(c)})^{2}}{2m_{c}}\right) \prod_{k=1}^{N_{b}} \exp\left(-\frac{\beta(\tilde{\boldsymbol{\pi}}_{k}^{(c)})^{2}}{2\tilde{I}_{k}}\right)\right]$$
- **what:** Maxwell-Boltzmann sampling of translational and fictitious angular momenta at start of each MC step (momentum refresh). Draw p_1 ~ N(0, mc/β), π̃_k ~ N(0, Ĩ_k/β).
- **symbols:** β = 1/(kB T); other symbols as in eq 8.

<!-- eq:11 -->
$$P_{\text{acc}} = \Pr\left(\boldsymbol{q}(\tilde{t}) \to \boldsymbol{q}(\tilde{t} + \Delta \tilde{t})\right) = \min\left\{1, \left[\prod_{c=1}^{N_c} \frac{\sin\left(\phi_2^{(c)}(\tilde{t} + \Delta \tilde{t})\right)}{\sin(\phi_2^{(c)}(\tilde{t}))}\right] \exp\left(-\beta \Delta \tilde{\mathcal{H}}\right)\right\}$$
- **what:** GC-HMC acceptance probability: Boltzmann factor of the total-Hamiltonian discretization error times the Euler-angle Jacobian ratio (product of sin φ2 over chains).
- **symbols:** φ2^(c) - second Euler angle of chain c; Δt̃ = N_MD δt̃_MD; ΔH̃ = [E(t̃+Δt̃)-E(t̃)] + [K̃(t̃+Δt̃)-K̃(t̃)] - discretization error of total Hamiltonian.

<!-- eq:12 -->
$$K_{\rm ang}(\{\phi_k^{(c)}\}, \{\omega_k^{(c)}\}) = \frac{1}{2} \sum_{i,j} \sum_{b=1}^{N_b} m_b \left(\frac{\partial \mathbf{r}_b}{\partial \phi_i}\right) \cdot \left(\frac{\partial \mathbf{r}_b}{\partial \phi_j}\right) \omega_i \omega_j \equiv \frac{1}{2} \sum_{i,j} I_{ij}(\{\phi\}) \omega_i \omega_j$$
- **what:** (Contrast, not used in the sampler) true angular kinetic energy in generalized coords; defines the configuration-dependent inertia tensor.
- **symbols:** ω_i - actual angular velocity of angle i; I_ij({φ}) = Σ_b mb (∂r_b/∂φ_i)·(∂r_b/∂φ_j) - inertia tensor; r_b - position of bead b.

<!-- eq:13 -->
$$\sum_{j} I_{ij} \frac{\mathrm{d}\omega_{j}}{\mathrm{d}t} + \sum_{j,k} \frac{\partial I_{ik}}{\partial \phi_{j}} \omega_{k} \frac{\mathrm{d}\phi_{j}}{\mathrm{d}t} = -\frac{\partial \mathscr{E}}{\partial \phi_{i}} + \frac{1}{2} \sum_{j,k} \frac{\partial I_{jk}}{\partial \phi_{i}} \omega_{j} \omega_{k}$$
- **what:** (Contrast, not used) coupled true equations of motion requiring matrix inversion; no simple time-reversible/area-preserving discretization.
- **symbols:** as eq 12; ∂I_ik/∂φ_j - derivative of inertia tensor; t - real time.

<!-- eq:14 -->
$$\mathscr{I}_k \equiv \sum_{b \in M_k} m_b [e_k \times (r_b - P_k)]^2$$
- **what:** Instantaneous moment of inertia of angle φ_k: rotate only φ_k, freeze rest. Sum of mb × (perpendicular distance to axis)^2 over moved beads. Its equilibrium mean ⟨I_k⟩ is the recommended choice for Ĩ_k.
- **symbols:** M_k - set of beads moved by angle φ_k; e_k - unit rotation axis of φ_k; P_k - point the axis passes through; r_b - bead position. For Euler angles: e_k = axis through origin r1, P_k = r1. For torsions: e_k = (r_{k-1}-r_{k-2})/||r_{k-1}-r_{k-2}||, P_k = r_{k-1}. Outermost torsion: I = mb lb^2 sin^2 θ.

<!-- eq:15 -->
$$\Delta t \equiv N_{\rm MD} \delta t_{\rm MD} \approx \left( \frac{\langle \mathscr{I}_k \rangle \langle (\Delta \phi_k)^2 \rangle}{k_{\rm B} T} \right)^{1/2}$$
- **what:** Effective real-time step per MC step from equipartition (½⟨I_k⟩⟨φ̇^2⟩ ≈ ½ kB T); when Ĩ_k = ⟨I_k⟩ this is ~uniform across DOF, mapping computer time to real time.
- **symbols:** Δt - real time per MC step; N_MD - MD steps per MC step; δt_MD - MD step (fictitious); ⟨I_k⟩ - mean instantaneous moment of inertia; ⟨(Δφ_k)^2⟩ - mean-square angle change per MC step.

<!-- eq:16 -->
$$f_{\text{EEV}}(t) = \left\langle \frac{(\mathbf{r}_{N_b}(t) - \mathbf{r}_1(t)) \cdot (\mathbf{r}_{N_b}(0) - \mathbf{r}_1(0))}{\|\mathbf{r}_{N_b}(t) - \mathbf{r}_1(t)\| \|\mathbf{r}_{N_b}(0) - \mathbf{r}_1(0)\|} \right\rangle$$
- **what:** Normalized end-to-end-vector autocorrelation function (sampling-efficiency diagnostic).
- **symbols:** r_Nb - last bead position; r_1 - first bead position; t - time.

<!-- eq:17 -->
$$f_{\text{BCF; 1}}(t) = \langle \mathbf{u}_i(t) \cdot \mathbf{u}_i(0) \rangle, \qquad f_{\text{BCF; 2}}(t) = \frac{3}{2} \langle (\mathbf{u}_i(t) \cdot \mathbf{u}_i(0))^2 \rangle - \frac{1}{2}$$
- **what:** First- and second-degree bond orientational autocorrelation functions (2nd is the P2 Legendre form), averaged over all bonds.
- **symbols:** u_i - unit vector along i-th bond of the chain; average over all bonds implied.

<!-- eq:18 -->
$$\left\langle \left[ \prod_{c=1}^{N_c} \frac{\sin\left(\phi_2^{(c)}(\tilde{t} + \Delta \tilde{t})\right)}{\sin\left(\phi_2^{(c)}(\tilde{t})\right)} \right] \exp\left(-\frac{\Delta \tilde{\mathcal{H}}}{k_B T}\right) \right\rangle = 1$$
- **what:** Modified fluctuation identity for GC-HMC (analog of ⟨exp(-ΔH̃/kBT)⟩=1) including the Euler-angle Jacobian; a code correctness check. Defining ΔH̃* = ΔH̃ - kB T ln(Πc) recovers ⟨exp(-ΔH̃*/kBT)⟩ = 1.
- **symbols:** Πc = Π_c sin(φ2(t̃+Δt̃))/sin(φ2(t̃)); other symbols as eq 11.

<!-- eq:19 -->
$$\langle P_{\rm acc} \rangle \approx \operatorname{erfc} \left( \frac{1}{2} \langle \Delta \widetilde{\mathcal{H}} / k_{\rm B} T \rangle^{1/2} \right)$$
- **what:** Gupta et al. relation for mean acceptance vs mean energy error; here use ΔH̃* in place of ΔH̃.
- **symbols:** erfc - complementary error function; ⟨ΔH̃/kBT⟩ - mean discretization error.

## Appendix (generic HMC restatement)

<!-- eq:A1 -->
$$\tilde{\mathscr{H}}(\{\tilde{\pmb{\pi}}\},\{q\}) = \mathscr{E}(\{q\}) + \sum_{i} \tilde{\pmb{\pi}}_{i}^{2}/2m_{i}$$
- **what:** Extended-phase-space Hamiltonian with diagonal fictitious kinetic energy.
- **symbols:** π̃_i - fictitious momentum of DOF i; m_i - constant coefficient (unity in Duane et al.).

<!-- eq:A2 -->
$$P_{\rm G}(\{\tilde{\pi}_0\})[\mathrm{d}\tilde{\pi}_0] \propto \exp\left(-\sum_i (\tilde{\pi}_0)_i^2/2m_ik_{\rm B}T\right)[\mathrm{d}\tilde{\pi}_0]$$
- **what:** Gaussian generation of initial fictitious momenta.
- **symbols:** π̃0 - initial momenta; m_i, kB, T as above.

<!-- eq:A3 -->
$$P_{\mathbf{A}} = \min\{1, \exp(-\beta \Delta \widetilde{\mathscr{H}})\}$$
- **what:** Plain HMC acceptance (no D-factor).
- **symbols:** β = 1/kB T; ΔH̃ = H̃({π̃n},{qn}) - H̃({π̃0},{q0}).

<!-- eq:A5 -->
$$P_{\mathbf{A}} = \min \left\{ 1, \exp\left(-\beta \Delta \widetilde{\mathcal{H}}\right) D(\left\{q_{n}\right\}) / D(\left\{q_{0}\right\}) \right\}$$
- **what:** HMC acceptance corrected by the coordinate-transform factor to sample Cartesian Boltzmann via {q} space.
- **symbols:** D({q}) - Jacobian (+Fixman) weight from eq 6.
