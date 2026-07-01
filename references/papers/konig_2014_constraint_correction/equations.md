# Equations - Correcting for the free energy costs of bond or angle constraints

Convention: reduced thermodynamic units, `β = 1/(k_B T)`. Energies in kcal/mol.
"Constrained" = a hard degree of freedom (bond/angle) frozen to a fixed value
(Dirac delta in the integral). "Releasing" a constraint restores that DOF.

<!-- eq:1 -->
$$ U(\Delta q) = U_0 + \frac{K}{2}\Delta q^2 $$
- **what:** Harmonic potential energy of a single internal coordinate (bond/angle) as a function of the deviation from the energy minimum.
- **symbols:** U - potential energy (kcal/mol); U_0 - zero-point/reference energy (kcal/mol); K - force constant (kcal/mol/A^2 for bonds, kcal/mol/rad^2 for angles); Delta q - deviation of the internal coordinate q from its energy minimum.

<!-- eq:2 -->
$$ Z^{h.o.} = \int e^{-\beta U(\Delta q)}\, d\Delta q = e^{-\beta U_0} \int e^{-\beta \frac{K}{2} \Delta q^2}\, d\Delta q $$
- **what:** Partition function of the harmonic oscillator, integrating over the single DOF.
- **symbols:** Z^{h.o.} - partition function of harmonic oscillator; beta = 1/(k_B T) (mol/kcal); other symbols as in eq:1.

<!-- eq:3 -->
$$ G^{h.o.} = -\beta^{-1} \ln Z^{h.o.} = U_0 - \beta^{-1} \ln \int e^{-\beta \frac{K}{2} \Delta q^2}\, d\Delta q $$
- **what:** Absolute free energy of the harmonic oscillator (before evaluating the Gaussian integral).
- **symbols:** G^{h.o.} - absolute free energy (kcal/mol); rest as above. <!-- CHECK: raw printed force constant as x in the exponent; corrected to K per eq:1/eq:2 -->

<!-- eq:4 -->
$$ G^{h.o.} = U_{0} - \beta^{-1}\ln\sqrt{\frac{2\pi}{\beta K}} $$
- **what:** Closed form of eq:3 using the Gaussian integral result.
- **symbols:** as eq:3. Gaussian integral: `∫ exp(-βK Δq²/2) dΔq = sqrt(2π/(βK))`.

<!-- eq:5 -->
$$ Z_{\text{cons}}^{\text{h.o.}} = e^{-\beta U(\Delta q_{\text{cons}})} $$
- **what:** Partition function of the constrained system: the integral collapses to a single point (Dirac delta), value = Boltzmann factor at the constrained coordinate.
- **symbols:** Z_cons^{h.o.} - constrained partition function; Delta q_cons - value of the internal coordinate at which the constraint is imposed (deviation from minimum).

<!-- eq:6 -->
$$ G_{cons}^{h.o.} = -\beta^{-1} \ln Z_{cons}^{h.o.} = U(\Delta q_{cons}) $$
- **what:** Free energy of the constrained system equals the potential energy at the constrained coordinate.
- **symbols:** G_cons^{h.o.} - constrained free energy (kcal/mol).

<!-- eq:7 -->
$$ \Delta G_{cons}^{h.o.} = G_{cons}^{h.o.} - G^{h.o.} = U(\Delta q_{cons}) - U_0 + \beta^{-1} \ln \sqrt{\frac{2\pi}{\beta K}} $$
- **what:** Free energy of IMPOSING a constraint = constrained minus unconstrained free energy.
- **symbols:** Delta G_cons^{h.o.} - free energy cost of imposing the constraint (kcal/mol).

<!-- eq:8 -->
$$ \Delta H^{h.o.} = U(\Delta q_{cons}) - U_0 $$
- **what:** Enthalpic contribution to the constraint free energy (temperature-independent). Note U(Δq_cons) inherently contains U_0.
- **symbols:** Delta H^{h.o.} - enthalpic contribution (kcal/mol); >= 0 always.

<!-- eq:9 -->
$$ \Delta G_{harm}^{h.o.} = \beta^{-1} \ln \sqrt{\frac{2\pi}{\beta K}} $$
- **what:** Entropic (vibrational) contribution from the harmonic potential.
- **symbols:** Delta G_harm^{h.o.} - vibrational entropy contribution (kcal/mol); always <= 0 (no entropy gained by adding a bond).

<!-- eq:10 -->
$$ \Delta G_{cons}^{h.o.} = \Delta H^{h.o.} + \Delta G_{harm}^{h.o.} $$
- **what:** Total constraint free energy = enthalpy + vibrational entropy contributions.
- **symbols:** as eqs 7-9.

<!-- eq:11 -->
$$ U(\Delta q) = U(\Delta q_{cons}) + U'(\Delta q_{cons})(\Delta q - \Delta q_{cons}) + \frac{U''(\Delta q_{cons})}{2}(\Delta q - \Delta q_{cons})^2 + \dots $$
- **what:** Taylor expansion of the anharmonic potential about the constrained conformation.
- **symbols:** U'(Delta q_cons) - first derivative (gradient) at constrained conf; U''(Delta q_cons) - second derivative (Hessian) at constrained conf.

<!-- eq:12 -->
$$ U(\Delta q) \approx U(\Delta q_{\rm cons}) + U'(\Delta q_{\rm cons})(\Delta q - \Delta q_{\rm cons}) + \frac{U''(\Delta q_{\rm cons})}{2}(\Delta q - \Delta q_{\rm cons})^{2} $$
- **what:** Harmonic (2nd-order) truncation of eq:11; anharmonic terms dropped.
- **symbols:** as eq:11. <!-- CHECK: raw eq (12) printed U'(Δq_cons)/2 in the quadratic coeff; corrected to U''(Δq_cons)/2 per eq:11 -->

<!-- eq:13 -->
$$ \Delta q_{\rm cons} \approx \frac{U'(\Delta q_{\rm cons})}{U''(\Delta q_{\rm cons})} $$
- **what:** One-step Newton-Raphson: displacement from the constrained coordinate to the true (unconstrained) energy minimum.
- **symbols:** as eq:11. This is the distance the DOF relaxes when the constraint is released.

<!-- eq:14 -->
$$ \Delta H^{a.o.} \approx \frac{U'(\Delta q_{cons})^2}{2\,U''(\Delta q_{cons})} $$
- **what:** Enthalpic contribution for the anharmonic oscillator, from the Newton-Raphson step (via eq:1).
- **symbols:** Delta H^{a.o.} - anharmonic-oscillator enthalpy contribution (kcal/mol). <!-- CHECK: raw printed U''(Δq_cons)^2 in denominator; corrected to U''(Δq_cons) so that ΔH = (1/2)K Δq² = U'^2/(2U'') is consistent with eqs 1 and 13 -->

<!-- eq:15 -->
$$ \Delta G_{harm}^{a.o.} \approx \beta^{-1} \ln \sqrt{\frac{2\pi}{\beta\, U''(\Delta q_{cons})}} $$
- **what:** Harmonic entropy contribution for the anharmonic oscillator, using U'' as the effective force constant.
- **symbols:** Delta G_harm^{a.o.} - anharmonic-oscillator vibrational entropy (kcal/mol).

<!-- eq:16 -->
$$ H = C^{\dagger} F C $$
- **what:** Reduced (m x m) Hessian for the constrained DOFs, obtained by projecting the full mass-weighted second-derivative matrix onto the reduced basis.
- **symbols:** H - reduced Hessian (m x m); C - orthonormal reduced basis set (3n x m), columns are mass-weighted Cartesian-displacement vectors c_i for each constrained DOF; F - full mass-weighted energy second-derivative (Hessian) matrix (3n x 3n); n - number of atoms; m - number of constrained DOFs; dagger - transpose/adjoint.

<!-- eq:17 -->
$$ \mathbf{g} = C^{\dagger} \mathbf{f} $$
- **what:** Reduced-basis force (gradient) vector, projecting the full-system force into the basis.
- **symbols:** g - reduced gradient vector (length m); f - full-system force vector (length 3n); C as eq:16.

<!-- eq:18 -->
$$ \Lambda = U^{\dagger} H U $$
- **what:** Eigendecomposition of the reduced Hessian (U is unitary here, an eigenvector matrix - NOT potential energy); diagonal Lambda holds the force constants/eigenvalues of the constrained DOFs.
- **symbols:** Lambda - diagonal eigenvalue matrix (force constants of decoupled constrained modes); U - unitary eigenvector matrix; H as eq:16.

<!-- eq:19 -->
$$ \Delta H \approx \frac{1}{2} \left\| \mathbf{g}^{\dagger} \Lambda^{-1} \mathbf{g} \right\|_{1} $$
- **what:** Multi-constraint enthalpy correction (matrix form of eq:14). Note: g must be expressed in the eigenbasis of Lambda (i.e. rotate g by U) before this contraction.
- **symbols:** ||...||_1 - sum over all elements; Lambda^{-1} - inverse of diagonal eigenvalue matrix; g as eq:17.

<!-- eq:20 -->
$$ \Delta G_{harm} \approx \left\| \frac{1}{2\beta} \ln \left( \frac{2\pi}{\beta} \Lambda^{-1} \right) \right\|_{1} $$
- **what:** Multi-constraint vibrational-entropy correction (matrix form of eq:15). Log applied element-wise since Lambda is diagonal.
- **symbols:** as eqs 15,19.

<!-- eq:21 -->
$$ \Delta G_{rls} = -\Delta G_{cons} = -\Delta H - \Delta G_{harm} - \Delta G_{Jacobian} $$
- **what:** Total free energy of RELEASING constraints in Cartesian space = negative of imposing; sums enthalpy, vibrational entropy, and Jacobian contributions.
- **symbols:** Delta G_rls - free energy of releasing constraints (kcal/mol); Delta G_Jacobian - Jacobian (internal->Cartesian) contribution.

<!-- eq:22 -->
$$ \Delta G_{Jacobian} = -\beta^{-1} \ln \left( \prod_{i} \frac{J_{i}^{after}}{J_{i}^{before}} \right) $$
- **what:** Jacobian free energy change from internal->Cartesian coordinate transformation, over all constrained DOFs.
- **symbols:** J_i - Jacobian factor for the i-th constraint's DOF; "before"/"after" - coordinates before/after applying the correction (constraint release).

<!-- eq:23 -->
$$ J_r = r_{ik}^2 $$
- **what:** Jacobian factor for a bond between atoms i and k (rigid-rotor analysis, Herschbach et al.).
- **symbols:** J_r - bond Jacobian; r_ik - bond length between atoms i and k (A).

<!-- eq:24 -->
$$ J_{\theta} = \sin \theta_{jkl} $$
- **what:** Jacobian factor for an angle theta between atoms j, k, l forming a linear chain.
- **symbols:** J_theta - angle Jacobian; theta_jkl - bond angle at vertex k (rad).

<!-- eq:25 -->
$$ J_{\theta'} = \sin^{-1}\theta'_{jk\ell} $$
- **what:** Jacobian factor for a branching angle theta' (terminal atom k branching off a chain), between bond r_jk and the plane defined by bonds r_lj and r_je. Code only handles linear-chain angles (eq:24); branched case is not implemented.
- **symbols:** J_theta' - branching-angle Jacobian; theta' - angle between the branch bond and the reference plane (rad). <!-- CHECK: raw index printed as θ'_{1kℓ'}; read as θ'_{jkℓ}. sin^{-1} here denotes reciprocal sine (1/sin), not arcsin, per the rigid-rotor Jacobian context -->

<!-- eq:26 -->
$$ Z = Z_{cons} \cdot Z_{add} = \int_{cons} \int_{add} e^{-\beta U(\mathbf{X})} = \int_{cons} e^{-\beta U_0(\Delta \mathbf{X}_{cons})} \int_{add} e^{-\beta U(\Delta \mathbf{X})} $$
- **what:** Factorization of the unconstrained partition function into constrained DOFs and the additional DOFs freed by releasing constraints.
- **symbols:** Z - unconstrained partition function; Z_cons - constrained-DOF part; Z_add - part from released (added) DOFs; X - full coordinate set; U_0(ΔX_cons) - potential depending only on the constrained coordinates.

<!-- eq:27 -->
$$ Z = \int_{cons} e^{-\beta U_0(\Delta \mathbf{X}_{cons})} e^{-\beta \Delta G_{cons}(\mathbf{X}_{cons})} $$
- **what:** Partition function with the inner (added-DOF) integral replaced by the per-configuration constraint free energy.
- **symbols:** Delta G_cons(X_cons) - free energy of releasing the constraint evaluated at each point of the constrained-DOF integral.

<!-- eq:28 -->
$$ \Delta G_{rls} = -\Delta G_{cons} = \beta^{-1} \ln \frac{Z}{Z_{cons}} = \beta^{-1} \ln \frac{\int_{cons} e^{-\beta U_0(\Delta \mathbf{x}_{cons})} e^{-\beta \Delta G_{cons}(\mathbf{x}_{cons})}}{\int_{cons} e^{-\beta U_0(\Delta \mathbf{x}_{cons})}} $$
- **what:** Free energy of releasing all constraints as the ratio of unconstrained to constrained partition functions.
- **symbols:** as eqs 26-27.

<!-- eq:29 -->
$$ \Delta G_{rls} = \beta^{-1} \ln \left\langle e^{-\beta \Delta G_{\text{cons}}(\mathbf{x}_{\text{cons}})} \right\rangle_{\text{cons}} $$
- **what:** Key working formula: releasing free energy is an exponential (Zwanzig / thermodynamic perturbation) average of the per-frame constraint free energy over the constrained ensemble. This is what is evaluated during trajectory post-processing.
- **symbols:** <...>_cons - ensemble average using constraints (over trajectory frames); Delta G_cons(x_cons) - per-frame constraint free energy from eqs 19-22.

## Auxiliary geometric relations (from prose, section 2.1)

<!-- eq:cbasis-bond -->
$$ \mathbf{c}_i = \frac{\partial \mathbf{x}}{\partial r_i} $$
- **what:** Basis vector for a constrained bond r_i: unit vector (mass-weighted, normalized) from atom 1 to atom 2 of the bond (and vice versa).
- **symbols:** c_i - basis vector (3n); x - 3n Cartesian coords; r_i - i-th constrained bond length.

<!-- eq:cbasis-angle -->
$$ \mathbf{c}_i = \frac{\partial \mathbf{x}}{\partial \theta_i} $$
- **what:** Basis vector for a constrained angle theta_i (atoms j,k,l, constraint on atom l): tangential vector normal to bond r_kl in the j-k-l plane.
- **symbols:** c_i - basis vector; theta_i - i-th constrained angle.

<!-- eq:angle-curvature -->
$$ \Delta \mathbf{x} = \sin(\Delta\theta_i)\,\mathbf{c}_i\, r_{kl} + (1-\cos(\Delta\theta_i))\,\mathbf{r}_{lk} $$
- **what:** Cartesian displacement for an angle change Delta theta_i that PRESERVES the bond length r_kl (adds the curvature term). Used by the CANG modes. Bonds (CBND) must be corrected before angles (CANG) so r_kl is current.
- **symbols:** Delta x - Cartesian displacement; Delta theta_i - angle change (rad); c_i - angle tangential basis vector; r_kl - bond length between atoms k and l (A); r_lk - Cartesian vector from atom l to atom k.
