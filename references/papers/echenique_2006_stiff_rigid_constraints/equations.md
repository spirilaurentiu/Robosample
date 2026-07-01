# Equations - Echenique, Calvo, Alonso 2006 (stiff vs rigid constraints, mass-metric tensors, Fixman)

Per-mole energy units throughout: `β := 1/RT` (R = gas constant), NOT `1/k_B T`.
Einstein summation on repeated indices. Index ranges: see `notation.md`.

---

<!-- eq:1 -->
$$ q^{I} = f^{I}(q^{i}) \qquad I = M+7, \dots, N $$
- **what:** the L = N-M-6 holonomic constraints defining the hypersurface Σ; each hard coordinate is a function of the soft internal coordinates.
- **symbols:** q^I - hard internal coordinate; q^i - soft internal coordinate (i = 7..M+6); f^I - constraint function; M - number of soft internal DOF; N=3n; L=N-M-6 number of hard coords.

<!-- eq:2 -->
$$ V(q^{i}, q^{I}) = \underbrace{V(q^{i}, f^{I}(q^{i}))}_{V_{\Sigma}(q^{i})} + \underbrace{\left[V(q^{i}, q^{I}) - V(q^{i}, f^{I}(q^{i}))\right]}_{V_{c}(q^{i}, q^{I})} $$
- **what:** split of total potential into the value on the constraint surface (PES) plus the constraining potential V_c (zero on Σ by construction).
- **symbols:** V_Σ(q^i) - potential energy on Σ (the PES from constrained geometry optimization); V_c - constraining potential, V_c(q^i, f^I(q^i))=0.

<!-- eq:3 -->
$$ V_{c}(q^{i}, q^{I}) \simeq V_{c}(q^{i}, f^{I}(q^{i})) + \left[\frac{\partial V_{c}}{\partial q^{J}}\right]_{\Sigma} (q^{J} - f^{J}(q^{i})) + \frac{1}{2} \underbrace{\left[\frac{\partial^{2} V_{c}}{\partial q^{J} \partial q^{K}}\right]_{\Sigma}}_{\mathcal{H}_{JK}(q^{i})} (q^{J} - f^{J}(q^{i})) (q^{K} - f^{K}(q^{i})) $$
- **what:** second-order Taylor expansion of V_c about Σ; 0th and 1st order terms vanish (conditions i, ii), defining the partial Hessian H_JK over hard coords only.
- **symbols:** [·]_Σ - evaluated on Σ (q^I=f^I(q^i)); H_JK(q^i) - Hessian of V_c wrt hard coords only, positive definite -> det H > 0.

<!-- eq:4 -->
$$ H_{s}(q^{\mu}, p_{\mu}) := \frac{1}{2} p_{\nu} G^{\nu\rho}(q^{u}, q^{I}) p_{\rho} + V_{\Sigma}(q^{i}) + \frac{1}{2} \mathcal{H}_{JK}(q^{i}) (q^{J} - f^{J}(q^{i})) (q^{K} - f^{K}(q^{i})) $$
- **what:** stiff-model Hamiltonian: kinetic (full inverse mass-metric tensor) + PES + harmonic constraining term.
- **symbols:** H_s - stiff Hamiltonian; p_μ - momentum conjugate to q^μ; G^{νρ} - inverse mass-metric tensor. <!-- CHECK: raw eq(4) shows a spurious factor k/2 on the Hessian term; from eq(3) the coefficient is 1/2, used consistently in eqs (8),(9). Corrected to 1/2. -->

<!-- eq:5 -->
$$ G_{\nu\rho}(q^{u}, q^{I}) := \sum_{\sigma=1}^{N} \frac{\partial x^{\sigma}(q^{\mu})}{\partial q^{\nu}} \, m_{\sigma} \, \frac{\partial x^{\sigma}(q^{\mu})}{\partial q^{\rho}} $$
- **what:** covariant mass-metric tensor: mass-weighted inner product of Jacobians of Euclidean coords wrt generalized coords. Core object for constraint corrections.
- **symbols:** G_{νρ} - mass-metric tensor (N×N); x^σ - Euclidean coordinate σ (σ=1..N=3n); m_σ - mass associated with Euclidean coord σ (atom mass of atom owning that coord).

<!-- eq:6 -->
$$ G^{\nu\sigma}(q^u, q^I)\, G_{\sigma\rho}(q^u, q^I) = \delta^{\nu}_{\rho} $$
- **what:** definition of inverse mass-metric tensor G^{νσ} via Kronecker delta.
- **symbols:** δ^ν_ρ - Kronecker delta.

<!-- eq:7 -->
$$ Z_{\rm s} = \frac{\alpha_{QM}}{h^N} \int dq^{\mu}\, dp_{\mu} \, \exp\left[-\beta H_{\rm s}(q^{\mu}, p_{\mu})\right] $$
- **what:** stiff partition function over full phase space.
- **symbols:** Z_s - stiff partition function; h - Planck constant; α_QM - quantum-indistinguishability combinatorial factor (e.g. 1/N! for N identical particles); β=1/RT.

<!-- eq:10 -->
$$ Z_{s} = \chi_{s}(T) \int dq^{u} \exp\left[-\beta \left(V_{\Sigma}(q^{i}) + T\frac{R}{2}\ln[\det \mathcal{H}(q^{i})] - T\frac{R}{2}\ln[\det G(q^{u}, f^{I}(q^{i}))]\right)\right] $$
- **what:** stiff partition function reduced to a coordinate-only integral over soft coords q^u, after integrating hard coords and momenta; two conformation-dependent entropy corrections appear.
- **symbols:** χ_s(T) - temperature prefactor (eq 11); det H - determinant of partial Hessian; det G - determinant of mass-metric tensor evaluated on Σ. <!-- CHECK: sign of the det G term: it enters F_s as -T S_s^k with S_s^k=+(R/2)ln det G (eq 12c), so in the exponent it is -T(R/2)ln det G; reconstructed accordingly. -->

<!-- eq:11 -->
$$ \chi_{s}(T) := \left(\frac{2\pi}{\beta}\right)^{\frac{N+L}{2}} \frac{\alpha_{QM}}{h^{N}} $$
- **what:** temperature-dependent multiplicative prefactor of the stiff coordinate partition function.
- **symbols:** N=3n; L=N-M-6.

<!-- eq:12a -->
$$ F_{s}(q^{u}) := V_{\Sigma}(q^{i}) - T\left(S_{s}^{c}(q^{i}) + S_{s}^{k}(q^{u})\right) $$
- **what:** stiff effective free energy = PES minus temperature times (conformational + kinetic entropies). This is the effective potential to sample for the stiff model.
- **symbols:** F_s - stiff effective free energy; S_s^c - conformational entropy (Hessian); S_s^k - kinetic entropy (mass-metric G).

<!-- eq:12b -->
$$ S_{s}^{c}(q^{i}) := -\frac{R}{2} \ln[\det \mathcal{H}(q^{i})] $$
- **what:** stiff conformational entropy from the determinant of the constraining Hessian (found to be the LARGEST correction).
- **symbols:** R - gas constant; H(q^i) - partial Hessian over hard coords.

<!-- eq:12c -->
$$ S_{s}^{k}(q^{u}) := \frac{R}{2} \ln\left[\det G(q^{u}, f^{I}(q^{i}))\right] $$
- **what:** stiff kinetic entropy from the determinant of the full mass-metric tensor G on Σ (found to be the SMALLEST correction).
- **symbols:** G - mass-metric tensor (eq 5) evaluated at hard coords = f^I(q^i).

<!-- eq:13 -->
$$ P_{s}(q^{u}) = \frac{\exp[-\beta F_{s}(q^{u})]}{Z'_{s}}, \qquad Z'_{s} := \int dq^{u} \exp[-\beta F_{s}(q^{u})] $$
- **what:** stiff equilibrium probability density over soft subspace E×Σ.
- **symbols:** P_s - stiff equilibrium density; Z'_s - configurational normalizer.

<!-- eq:14 -->
$$ H_{\rm r}(q^u, \eta_u) := \frac{1}{2} \eta_{\nu} \, g^{\nu w}(q^u) \, \eta_w + V_{\Sigma}(q^i) $$
- **what:** rigid-model Hamiltonian: reduced-kinetic + PES only (no Hessian term; constraints exact/holonomic).
- **symbols:** H_r - rigid Hamiltonian; η_u - momentum conjugate to soft coord q^u in the reduced space; g^{νw} - inverse reduced mass-metric tensor.

<!-- eq:15 -->
$$ g_{vw}(q^{u}) = G_{vw} + \frac{\partial f^{J}}{\partial q^{v}} G_{JK} \frac{\partial f^{K}}{\partial q^{w}} + G_{vK} \frac{\partial f^{K}}{\partial q^{w}} + \frac{\partial f^{J}}{\partial q^{v}} G_{Jw} = \frac{\partial \tilde{f}^{\mu}}{\partial q^{\nu}} G_{\mu\nu} \frac{\partial \tilde{f}^{\nu}}{\partial q^{w}} $$
- **what:** reduced mass-metric tensor on E×Σ: pull-back of full G through the constraint map. All G blocks evaluated at (q^u, f^I(q^i)).
- **symbols:** g_{vw} - reduced mass-metric tensor ((M+6)×(M+6)); indices v,w soft (1..M+6); J,K hard; ∂f^J/∂q^v - constraint Jacobian; f̃^μ - eq 16.

<!-- eq:16 -->
$$ \tilde{f}^{\mu} := \begin{cases} q^{u} & u := \mu = 1, \dots, M+6 \\ f^{I}(q^{i}) & I := \mu = M+7, \dots, N \end{cases} $$
- **what:** the combined map f̃ giving all Euclidean coords as functions of soft coords (identity on soft, constraint on hard).
- **symbols:** f̃^μ - full parametrization map.

<!-- eq:17 -->
$$ H(q^{\mu}, p_{\mu}) := \frac{1}{2} p_{\nu} G^{\nu \rho}(q^{\mu}) p_{\rho} + V(q^{a}) $$
- **what:** unconstrained (fully flexible) Hamiltonian in E×I, from which the rigid one is derived.
- **symbols:** H - unconstrained Hamiltonian; V(q^a) - full potential over all internal coords.

<!-- eq:18 -->
$$ \dot{q}^I := \frac{\partial f^I(q^i)}{\partial q^{j}} \, \dot{q}^{j} $$
- **what:** time derivative of the constraint relation (velocities of hard coords slaved to soft velocities).
- **symbols:** overdot - time derivative.

<!-- eq:19 -->
$$ \eta_{\nu} := g_{\nu w}(q^{u}) \, \dot{q}^{w} = g_{\nu w}(q^{u}) \, G^{w \mu}(q^{u}, f^{I}(q^{i})) \, p_{\mu} $$
- **what:** definition of reduced momenta conjugate to soft coordinates.
- **symbols:** η_ν - reduced momentum.

<!-- eq:20 -->
$$ Z_{\rm r} = \frac{\alpha_{QM}}{h^{M+6}} \int dq^u \, d\eta_u \exp\left[-\beta \left(\frac{1}{2}\eta_{\nu}\,g^{\nu w}(q^u)\,\eta_{w} + V_{\Sigma}(q^i)\right)\right] $$
- **what:** rigid partition function over reduced phase space (dimension M+6 coords).
- **symbols:** Z_r - rigid partition function.

<!-- eq:21 -->
$$ Z_{\rm r} = \chi_{\rm r}(T) \int dq^u \exp\left[-\beta \left(V_{\Sigma}(q^i) - T\frac{R}{2}\ln[\det g(q^u)]\right)\right] $$
- **what:** rigid partition function reduced to coordinate-only integral; single kinetic-entropy correction from det g.
- **symbols:** det g - determinant of reduced mass-metric tensor.

<!-- eq:22 -->
$$ \chi_{\rm r}(T) := \left(\frac{2\pi}{\beta}\right)^{\frac{M+6}{2}} \frac{\alpha_{QM}}{h^{M+6}} $$
- **what:** temperature prefactor of the rigid coordinate partition function.
- **symbols:** M+6 - number of soft (external+internal) coordinates. <!-- CHECK: raw shows h^{(M+6)/2}; by dimensional analogy with eq(11) χ_s ~ α/h^N it should be h^{M+6}. Flagged. -->

<!-- eq:23a -->
$$ F_{\rm r}(q^u) := V_{\Sigma}(q^i) - T\,S_{\rm r}^{\rm k}(q^u) $$
- **what:** rigid effective free energy = PES minus temperature times kinetic entropy. Effective potential to sample for the rigid model.
- **symbols:** F_r - rigid effective free energy; S_r^k - rigid kinetic entropy.

<!-- eq:23b -->
$$ S_{\rm r}^{\rm k}(q^{u}) := \frac{R}{2} \ln[\det g(q^{u})] $$
- **what:** rigid kinetic entropy from determinant of reduced mass-metric tensor g.
- **symbols:** g - reduced mass-metric tensor (eq 15).

<!-- eq:24 -->
$$ P_{r}(q^{u}) = \frac{\exp[-\beta F_{r}(q^{u})]}{Z'_{r}}, \qquad Z'_{r} := \int dq^{u} \exp[-\beta F_{r}(q^{u})] $$
- **what:** rigid equilibrium probability density over soft subspace.
- **symbols:** P_r - rigid equilibrium density.

<!-- eq:25 -->
$$ V_{\mathrm{F}}(q^{u}) := T S_{\rm r}^{\rm k}(q^{u}) - T S_{\rm s}^{c}(q^{i}) - T S_{\rm s}^{\rm k}(q^{u}) = \frac{RT}{2} \ln \left[ \frac{\det G(q^{u})}{\det \mathcal{H}(q^{i}) \, \det g(q^{u})} \right] $$
- **what:** FIXMAN COMPENSATING POTENTIAL = F_s - F_r. Add this to V_Σ in a rigid MD/MC simulation to recover the stiff equilibrium distribution. Central implementable quantity.
- **symbols:** V_F - Fixman potential; det G, det H, det g as above; RT=1/β.

<!-- eq:26 -->
$$ \det G = \left(\prod_{\alpha=1}^{n} m_{\alpha}^{3}\right) \sin^{2} \theta \left(\prod_{\alpha=2}^{n} r_{\alpha}^{4}\right) \left(\prod_{\alpha=3}^{n} \sin^{2} \theta_{\alpha}\right) $$
- **what:** closed form for det G in SASMIC internal (Z-matrix) coordinates: product of mass factors, one external-orientation sin^2θ, bond-length^4, and bond-angle sin^2 factors. Independent of dihedral angles explicitly.
- **symbols:** m_α - mass of atom α; n - number of atoms; r_α - bond length of atom α (α=2..n); θ_α - bond angle of atom α (α=3..n); θ - one external orientation angle. Go-Scheraga / Volkenstein serial-polymer result generalized.

<!-- eq:27 -->
$$ S_{\rm s}^{\rm k}(q^i) = \frac{R}{2} \left[ \sum_{\alpha=2}^n \ln(r_{\alpha}^4) + \sum_{\alpha=3}^n \ln(\sin^2 \theta_{\alpha}) \right] $$
- **what:** stiff kinetic entropy from det G, external factors and mass constants dropped (conformation-independent). Directly computable from bond lengths and bond angles.
- **symbols:** r_α - bond lengths; θ_α - bond angles. Additive constants omitted.

<!-- eq:28 -->
$$ \det g = \sin^2 \theta \; \det g_2(q^i) $$
- **what:** det of reduced mass-metric tensor factorizes into external sin^2θ times an internal-only determinant det g_2.
- **symbols:** g_2 - internal reduced mass-metric matrix (eq 29); θ - external orientation angle.

<!-- eq:29 -->
$$ g_{2} = \begin{pmatrix} m_{\text{tot}} I^{(3)} & m_{\text{tot}} v(\vec{R}) & \cdots & m_{\text{tot}} \dfrac{\partial \vec{R}}{\partial q^{j}} & \cdots \\[2mm] m_{\text{tot}} v^{T}(\vec{R}) & \mathcal{J} & \cdots & \displaystyle\sum_{\alpha} m_{\alpha} \frac{\partial \vec{x}_{\alpha}'}{\partial q^{j}} \times \vec{x}_{\alpha}' & \cdots \\[2mm] \vdots & \vdots & \ddots & \vdots & \\[1mm] m_{\text{tot}} \dfrac{\partial \vec{R}}{\partial q^{i}} & \displaystyle\sum_{\alpha} m_{\alpha} \left( \frac{\partial \vec{x}_{\alpha}'}{\partial q^{i}} \times \vec{x}_{\alpha}' \right)^{\mathsf{T}} & \cdots & \displaystyle\sum_{\alpha} m_{\alpha} \frac{\partial \vec{x}_{\alpha}'^{T}}{\partial q^{i}} \frac{\partial \vec{x}_{\alpha}'}{\partial q^{j}} & \cdots \\[1mm] \vdots & \vdots & & \vdots & \ddots \end{pmatrix} $$
- **what:** internal reduced mass-metric matrix g_2 in block form: (translation-translation) m_tot I; (translation-rotation) m_tot v(R); (translation-internal) m_tot ∂R/∂q^j; (rotation-rotation) inertia tensor J; (rotation-internal) Σ m_α (∂x'_α/∂q^j)×x'_α; (internal-internal) Σ m_α (∂x'_α/∂q^i)·(∂x'_α/∂q^j).
- **symbols:** m_tot=Σm_α total mass; I^(3) - 3×3 identity; v(R) - skew matrix (eq 31); R - center of mass in primed frame; J - inertia tensor (eq 30); x'_α - position of atom α in body-fixed (primed) frame; i,j soft internal indices. <!-- CHECK: raw eq(29) is heavily OCR-garbled with duplicated rows; block structure reconstructed from surrounding prose (m_tot, v(R), J, cross-product and dot-product blocks). Verify exact index placement and transpose conventions against ref 58 before coding. -->

<!-- eq:30 -->
$$ \mathcal{J} := \begin{pmatrix} \sum_{\alpha} m_{\alpha} ((x_{\alpha}'^{2})^{2} + (x_{\alpha}'^{3})^{2}) & -\sum_{\alpha} m_{\alpha} x_{\alpha}'^{1} x_{\alpha}'^{2} & -\sum_{\alpha} m_{\alpha} x_{\alpha}'^{1} x_{\alpha}'^{3} \\ -\sum_{\alpha} m_{\alpha} x_{\alpha}'^{1} x_{\alpha}'^{2} & \sum_{\alpha} m_{\alpha} ((x_{\alpha}'^{1})^{2} + (x_{\alpha}'^{3})^{2}) & -\sum_{\alpha} m_{\alpha} x_{\alpha}'^{2} x_{\alpha}'^{3} \\ -\sum_{\alpha} m_{\alpha} x_{\alpha}'^{1} x_{\alpha}'^{3} & -\sum_{\alpha} m_{\alpha} x_{\alpha}'^{2} x_{\alpha}'^{3} & \sum_{\alpha} m_{\alpha} ((x_{\alpha}'^{1})^{2} + (x_{\alpha}'^{2})^{2}) \end{pmatrix} $$
- **what:** inertia tensor in the body-fixed (primed) reference frame.
- **symbols:** x'^k_α - k-th Cartesian component (k=1,2,3) of atom α in primed frame; m_α - atom mass.

<!-- eq:31 -->
$$ v(\vec{R}) := \begin{pmatrix} 0 & -R^3 & R^2 \\ R^3 & 0 & -R^1 \\ -R^2 & R^1 & 0 \end{pmatrix} $$
- **what:** skew-symmetric cross-product matrix of the center-of-mass position R (so v(R) w = R × w).
- **symbols:** R^k - k-th component of center of mass R in primed frame.

<!-- eq:32 -->
$$ S_{\rm r}^{\rm k}(q^i) = \frac{R}{2} \ln[\det g_2(q^i)] $$
- **what:** rigid kinetic entropy depending only on soft internals, after integrating out external sin^2θ.
- **symbols:** g_2 - internal reduced mass-metric matrix (eq 29); R - gas constant.

<!-- eq:33 -->
$$ N_{\rm res} = \left(\frac{RT}{d_{12}}\right)^2 $$
- **what:** maximum number of residues in an additive per-residue polypeptide potential before the accumulated statistical distance from dropping a correction term exceeds thermal energy RT.
- **symbols:** N_res - residue limit; d_12 - statistical energy distance between reference V_1 and approximation V_2; RT≈0.6 kcal/mol at 300 K.

<!-- eq:34 -->
$$ d_{12} = \sqrt{2}\,\sigma_2 \left(1 - r_{12}^2\right)^{1/2} $$
- **what:** statistical distance between two energy functions (typical error in energy differences when using V_2 for V_1, allowing linear rescaling).
- **symbols:** σ_2 - standard deviation of V_2 over the working set; r_12 - Pearson correlation coefficient between V_1 and V_2.

<!-- eq:A1 -->
$$ h^{IJ}(q^{\mu}) := \sum_{\sigma=1}^{N} \frac{\partial q^{I}}{\partial x^{\sigma}} \frac{1}{m_{\sigma}} \frac{\partial q^{J}}{\partial x^{\sigma}} $$
- **what:** Fixman's sparse matrix over hard coords; under approximation (iii) det G/det g = 1/det h, giving V_F = (RT/2) ln det h. Sparse because each internal coord involves few atoms.
- **symbols:** h^{IJ} - Fixman hard-coord matrix; q^I,q^J - hard coords; x^σ - Euclidean coords; m_σ - mass.

<!-- eq:A2 -->
$$ V_{\mathrm{ff}}(q^{a}) := \frac{1}{2} \sum_{\alpha=1}^{N_{r}} K_{r_{\alpha}} (r_{\alpha} - r_{\alpha}^{0})^{2} + \frac{1}{2} \sum_{\alpha=1}^{N_{\theta}} K_{\theta_{\alpha}} (\theta_{\alpha} - \theta_{\alpha}^{0})^{2} + V_{\mathrm{ff}}^{\mathrm{tors}}(\phi_{\alpha}) + V_{\mathrm{ff}}^{\mathrm{long\text{-}range}}(q^{a}) $$
- **what:** generic classical force-field potential form (harmonic bonds + angles + torsions + long-range), used to argue det H is conformation-dependent even for force fields.
- **symbols:** K_{r_α}, K_{θ_α} - force constants; r_α^0, θ_α^0 - equilibrium bond length / angle; N_r, N_θ - counts; V_ff^tors torsional term; V_ff^long-range Coulomb/vdW.

<!-- eq:A3 -->
$$ Q^{i} := q^{i} \quad (i = 7, \dots, M+6), \qquad Q^{I} := q^{I} - f^{I}(q^{i}) + C^{I} \quad (I = M+7, \dots, N) $$
- **what:** definition of "exactly separable" hard/soft coordinates Q^a such that the hard coords equal constants C^I on Σ.
- **symbols:** Q^a - exactly-separable coords; C^I - arbitrary constants.
