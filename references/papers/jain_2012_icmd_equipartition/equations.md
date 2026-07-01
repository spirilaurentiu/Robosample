# Equations — Jain 2012, Equipartition Principle for ICMD

Transpose is denoted `x*`. Reduced/statistical-mechanics conventions:
`k` = Boltzmann constant, `T` = temperature, `kT` has units of energy.

<!-- eq:1 -->
$$ p_i = \frac{\partial \mathcal{L}(q, \dot{q})}{\partial \dot{q}_i} $$
- **what:** Conjugate momentum is the derivative of the Lagrangian w.r.t. generalized velocity.
- **symbols:** p_i - ith conjugate momentum; L(q,q̇) - Lagrangian; q̇_i - ith generalized velocity; q ∈ R^n generalized coords.

<!-- eq:2 -->
$$ \mathcal{Z}(T) = \frac{1}{h^n} \int_{-\infty}^{\infty} \mathrm{d}p_1 \dots \int_{-\infty}^{\infty} \mathrm{d}p_n \int_{-\alpha_1}^{\gamma_1} \mathrm{d}q_1 \dots \int_{-\alpha_n}^{\gamma_n} \mathrm{d}q_n \; e^{-H(q,p)/kT} $$
- **what:** Canonical partition function for an n-DOF system.
- **symbols:** Z(T) - partition function; h - Planck constant; H(q,p) - Hamiltonian; α_i, γ_i - geometry-determined integration limits; n - number of DOF.

<!-- eq:3 -->
$$ \langle f(q, p) \rangle = \frac{1}{h^n \mathcal{Z}(T)} \int_{-\infty}^{\infty} dp_1 \dots \int_{-\infty}^{\infty} dp_n \int_{-\alpha_1}^{\gamma_1} dq_1 \dots \int_{-\alpha_n}^{\gamma_n} dq_n \; f(q, p)\, e^{-H(q, p)/kT} $$
- **what:** Canonical ensemble average of a phase-space function f(q,p).
- **symbols:** ⟨·⟩ - ensemble average; f(q,p) - observable.

<!-- eq:4 -->
$$ \left\langle y_i \frac{\partial H(q, p)}{\partial y_j} \right\rangle = kT\,\delta_{ij} $$
- **what:** Generalized equipartition theorem (valid only for CANONICAL coordinates y ∈ {q,p}).
- **symbols:** y_i, y_j - any pair of phase-space coordinates (elements of q or p); δ_{ij} - Kronecker delta (1 if i=j else 0). <!-- CHECK: paper's "δ_{i=j}" mislabeled a "Dirac delta"; it is the Kronecker delta here. -->

<!-- eq:5 -->
$$ \langle y\, \nabla_{y} H(q, p) \rangle = kT\, I_{2n} $$
- **what:** Matrix form of the equipartition theorem gathering all i,j; the outer-product matrix of ⟨y_i ∂H/∂y_j⟩ equals kT times identity.
- **symbols:** y - full phase-space vector (2n components); ∇_y H - gradient of H; I_{2n} - 2n×2n identity.

<!-- eq:6 -->
$$ H(x, p) = \sum_{i=1}^{3n} p_i^2 / 2m_i + \mathcal{U}(x) $$
- **what:** Separable Cartesian MD Hamiltonian for n atoms (3n DOF).
- **symbols:** x ∈ R^{3n} - atom positions; p_i - linear momentum component; m_i - atom mass; U(x) - potential energy.

<!-- eq:7 -->
$$ \langle p_i^2 / 2m_i \rangle = kT/2 $$
- **what:** Cartesian equipartition of kinetic energy — kT/2 per momentum coordinate.
- **symbols:** as eq:6.

<!-- eq:8 -->
$$ \mathcal{K}_e(\theta, p) = \frac{1}{2}\dot{\theta}^{*}\mathcal{M}(\theta)\dot{\theta} = \frac{1}{2}p^{*}\mathcal{M}^{-1}(\theta)p $$
- **what:** ICMD kinetic energy: quadratic in θ̇ via mass matrix, or in p via its inverse.
- **symbols:** K_e - kinetic energy; θ ∈ R^N - generalized coords; θ̇ - generalized velocities; M(θ) ∈ R^{N×N} - configuration-dependent, symmetric positive-definite mass matrix; p - conjugate momenta; N - number of ICMD DOF.

<!-- eq:9 -->
$$ p = \mathcal{M}(\theta)\dot{\theta} $$
- **what:** ICMD conjugate momenta from generalized velocities.
- **symbols:** as eq:8.

<!-- eq:10 -->
$$ H(\theta, p) = \frac{1}{2} p^* \mathcal{M}^{-1}(\theta) p + \mathcal{U}(\theta) $$
- **what:** ICMD Hamiltonian — NON-separable because K_e depends on θ via M^{-1}(θ).
- **symbols:** U(θ) - potential energy; others as eq:8.

<!-- eq:11 -->
$$ \mathcal{Z}(T) = \frac{1}{h^{N}} \int_{-\infty}^{\infty} \mathrm{d}p_1 \dots \int_{-\infty}^{\infty} \mathrm{d}p_{N} \int_{-\alpha_1}^{\gamma_1} \mathrm{d}\theta_1 \dots \int_{-\alpha_N}^{\gamma_N} \mathrm{d}\theta_{N} \; e^{-\left[\frac{1}{2}p^*\mathcal{M}^{-1}(\theta)p + \mathcal{U}(\theta)\right]/kT} $$
- **what:** ICMD canonical partition function.
- **symbols:** as eqs:8,10.

<!-- eq:12 -->
$$ \frac{\partial H(\theta, p)}{\partial p_j} = \frac{\partial p^*}{\partial p_j} \mathcal{M}^{-1}(\theta)\, p = e_j^* \dot{\theta} = \dot{\theta}_j $$
- **what:** Hamilton's equation: ∂H/∂p_j equals the jth generalized velocity.
- **symbols:** e_j - unit vector (1 in jth slot); θ̇_j - jth generalized velocity.

<!-- eq:13 -->
$$ \left\langle p_i \frac{\partial H(\theta, p)}{\partial p_j} \right\rangle = \left\langle p_i \dot{\theta}_j \right\rangle = kT\, \delta_{ij} $$
- **what:** Equipartition applied to ICMD momenta — has energy units but is NOT interpretable as kinetic-energy equipartition (p_i couples multiple velocities).
- **symbols:** as eqs:8,12.

<!-- eq:14 -->
$$ \langle \dot{\theta}\, p^* \rangle = \langle \mathcal{M}^{-1}(\theta) p\, p^* \rangle = kT\, I_{N} $$
- **what:** Matrix form of eq:13; no clean way to distribute thermal energy across DOF.
- **symbols:** I_N - N×N identity.

<!-- eq:15 -->
$$ \mathcal{M}(\theta) = m(\theta)\, m^*(\theta) $$
- **what:** Factorization of the SPD mass matrix (m is an invertible "square root"-type factor, not necessarily Cholesky).
- **symbols:** m(θ) ∈ R^{N×N} - invertible factor.

<!-- eq:16 -->
$$ m^{-1}(\theta) = l(\theta) \quad\text{and}\quad \mathcal{M}^{-1}(\theta) = l^*(\theta)\, l(\theta) $$
- **what:** l is the inverse of m; the inverse mass matrix factorizes as l* l.
- **symbols:** l(θ) ∈ R^{N×N} - inverse factor.

<!-- eq:17 -->
$$ \nu(\theta, p) \stackrel{\Delta}{=} m^*(\theta)\dot{\theta} = l(\theta)\, p $$
- **what:** DEFINITION of modal coordinates ν (noncanonical; decouple the kinetic energy).
- **symbols:** ν ∈ R^N - modal velocity coordinates; others as eqs:8,15,16.

<!-- eq:18 -->
$$ H(\theta, \nu) = \frac{1}{2} \nu^* \nu + \mathcal{U}(\theta) $$
- **what:** Hamiltonian in modal coordinates — kinetic energy decoupled, no cross-terms, no explicit θ dependence.
- **symbols:** as eq:17. Kinetic energy of modal component i is ν_i²/2.

<!-- eq:19 -->
$$ \det\{\mathcal{M}(\theta)\}^{1/2}\, d\nu_1 \dots d\nu_N = dp_1 \dots dp_N $$
- **what:** Volume-element (Jacobian) relation between modal and momentum coordinate spaces.
- **symbols:** det{M(θ)}^{1/2} - Jacobian of the p→ν change of variables.

<!-- eq:20 -->
$$ \mathcal{Z}(T) = \frac{1}{h^{N}} \int_{-\infty}^{\infty} d\nu_1 \dots \int_{-\infty}^{\infty} d\nu_N \int_{-\alpha_1}^{\gamma_1} d\theta_1 \dots \int_{-\alpha_N}^{\gamma_N} d\theta_N \; \det\{\mathcal{M}(\theta)\}^{1/2}\, e^{-\left[\frac{1}{2}\nu^*\nu + \mathcal{U}(\theta)\right]/kT} $$
- **what:** Partition function in modal coordinates (with det{M}^{1/2} weight).
- **symbols:** as eqs:18,19.

<!-- eq:21 -->
$$ \langle f(\theta, \nu) \rangle = \frac{1}{h^{N} \mathcal{Z}(T)} \int_{-\infty}^{\infty} d\nu_1 \dots \int_{-\infty}^{\infty} d\nu_N \int_{-\alpha_1}^{\gamma_1} d\theta_1 \dots \int_{-\alpha_N}^{\gamma_N} d\theta_N \; f(\theta, \nu)\, \det\{\mathcal{M}^{1/2}(\theta)\}\, e^{-\left[\frac{1}{2}\nu^*\nu + \mathcal{U}(\theta)\right]/kT} $$
- **what:** Ensemble average in modal coordinates.
- **symbols:** as eq:20.

<!-- eq:22 -->
$$ \mathcal{U}'(\theta) \stackrel{\Delta}{=} \mathcal{U}(\theta) + \mathcal{U}_c(\theta), \quad \mathcal{U}_c(\theta) \stackrel{\Delta}{=} \frac{1}{2} \ln \det \{ \mathcal{M}(\theta) \} $$
- **what:** Corrected (effective) potential absorbing the mass-matrix determinant — the Fixman-type compensating term.
- **symbols:** U'(θ) - effective potential; U_c(θ) - configurational (Fixman) correction.

<!-- eq:23 -->
$$ \langle f(\theta, \nu) \rangle = \frac{1}{h^{N} \mathcal{Z}(T)} \int_{-\infty}^{\infty} d\nu_1 \dots \int_{-\infty}^{\infty} d\nu_N \int_{-\alpha_1}^{\gamma_1} d\theta_1 \dots \int_{-\alpha_N}^{\gamma_N} d\theta_N \; f(\theta, \nu)\, e^{-\left[\frac{1}{2}\nu^*\nu + \mathcal{U}'(\theta)\right]/kT} $$
- **what:** Ensemble average with the corrected potential U' — factorizes cleanly into modal (Gaussian) and configurational parts.
- **symbols:** as eq:22.

<!-- eq:25 -->
$$ \langle \nu_i \nu_j \rangle = kT\, \delta_{ij} $$
- **what:** EQUIPARTITION PRINCIPLE FOR ICMD (modal coordinates): modal components uncorrelated, each carries kT/2 kinetic energy on average. Analog of Cartesian eq:7.
- **symbols:** as eq:17; δ_{ij} - Kronecker delta.

<!-- eq:26 -->
$$ \langle \nu\, \nu^* \rangle = kT\, I_{N} $$
- **what:** Matrix form of the ICMD equipartition principle (used for velocity initialization).
- **symbols:** I_N - N×N identity.

<!-- eq:27 -->
$$ \nu = m^*(\theta)\dot{\theta} \quad\text{and}\quad \dot{\theta} = l^*(\theta)\nu $$
- **what:** Physical ↔ modal velocity transformations.
- **symbols:** as eqs:15,16,17.

<!-- eq:28 -->
$$ \dot{x} = J(\theta)\dot{\theta}, \quad J(\theta) \stackrel{\Delta}{=} \left[ \frac{\partial x_a}{\partial \theta_b} \right]_{a=1..3n,\, b=1..N} $$
- **what:** Cartesian-to-ICMD velocity Jacobian relation and its definition.
- **symbols:** ẋ ∈ R^{3n} - Cartesian velocities; J(θ) ∈ R^{3n×N} - Jacobian; ∂x_a/∂θ_b - element.

<!-- eq:29 -->
$$ \mathcal{M}(\theta) = J^*(\theta)\, M\, J(\theta) $$
- **what:** ICMD mass matrix built from the Jacobian and the constant Cartesian atom-mass matrix.
- **symbols:** M ∈ R^{3n×3n} - constant DIAGONAL matrix of atom masses (Cartesian mass matrix); J(θ) as eq:28.

<!-- eq:30 -->
$$ m(\theta) = J^*(\theta)\, M^{1/2} $$
- **what:** For BAT models (N=3n, J square invertible), the factor m(θ) is directly this — square and invertible.
- **symbols:** M^{1/2} - diagonal matrix of sqrt atom masses.

<!-- eq:31 -->
$$ \mathcal{M} = \mathcal{H}\phi \mathcal{M}\phi^* \mathcal{H}^* $$
$$ \mathcal{M} = [I + \mathcal{H}\phi \mathcal{K}]\, \mathcal{D}\, [I + \mathcal{H}\phi \mathcal{K}]^* $$
$$ [I + \mathcal{H}\phi \mathcal{K}]^{-1} = [I - \mathcal{H}\psi \mathcal{K}] $$
$$ \mathcal{M}^{-1} = [I - \mathcal{H}\psi \mathcal{K}]^*\, \mathcal{D}^{-1}\, [I - \mathcal{H}\psi \mathcal{K}] $$
- **what:** Spatial-operator factorizations. (1) Newton–Euler factorization (nonsquare factors). (2) Square factorization: block-diagonal D, block-lower-triangular [I+HφK] (AB algorithm). (3) Analytical inverse of [I+HφK]. (4) Analytical inverse mass matrix.
- **symbols:** H - hinge-articulation operator; φ - rigid-body propagation operator; M (script) - link spatial-inertia operator; K - gain operator; ψ - articulated-body propagation operator; D - block-diagonal articulated-body inertia operator; I - identity. (See ref 19, Jain 2011, for operator definitions.)

<!-- eq:32 -->
$$ m(\theta) = [I + \mathcal{H}\phi\mathcal{K}]\,\mathcal{D}^{1/2}, \qquad l(\theta) = \mathcal{D}^{-1/2}[I - \mathcal{H}\psi\mathcal{K}] $$
- **what:** Analytical spatial-operator transformation matrices (valid for general tree-topology ICMD, constrained or not).
- **symbols:** as eq:31; D^{1/2}, D^{-1/2} - block-wise sqrt / inverse-sqrt of block-diagonal D.

<!-- eq:33 -->
$$ \nu(\theta) = \mathcal{D}^{1/2}[I + \mathcal{H}\phi\mathcal{K}]^*\, \dot{\theta}, \qquad \dot{\theta} = [I - \mathcal{H}\psi\mathcal{K}]^*\, \mathcal{D}^{-1/2}\nu $$
- **what:** Modal↔physical velocity transformations in spatial-operator form (O(N²) naive; O(N) via Table 1 recursion).
- **symbols:** as eqs:31,32.

## Table 1 — Recursive base-to-tip O(N) algorithms (serial topology)

Intermediate `V(k) ∈ R^6` is the combined angular+linear velocity of the kth
coordinate frame; loop runs k = N … 1 with `V(N+1) = 0`.

**Left: physical → modal**, $\nu = \mathcal{D}^{1/2}[I + \mathcal{H}\phi\mathcal{K}]^*\dot{\theta}$
<!-- eq:t1-left -->
$$ \mathcal{V}^+(k) = \phi^*(k+1,k)\,\mathcal{V}(k+1) $$
$$ \nu(k) = \mathcal{D}^{1/2}(k)\left[\dot{\theta}(k) + \mathcal{G}^*(k)\,\mathcal{V}^+(k)\right] $$
$$ \mathcal{V}(k) = \mathcal{V}^+(k) + \mathcal{H}^*(k)\,\dot{\theta}(k) $$

**Right: modal → physical**, $\dot{\theta} = [I - \mathcal{H}\psi\mathcal{K}]^*\mathcal{D}^{-1/2}\nu$
<!-- eq:t1-right -->
$$ \mathcal{V}^+(k) = \phi^*(k+1,k)\,\mathcal{V}(k+1) $$
$$ \dot{\theta}(k) = \mathcal{D}^{-1/2}(k)\,\nu(k) - \mathcal{G}^*(k)\,\mathcal{V}^+(k) $$
$$ \mathcal{V}(k) = \mathcal{V}^+(k) + \mathcal{H}^*(k)\,\dot{\theta}(k) $$
- **what:** O(N) base-to-tip recursions implementing eq:33 without forming m, l, or M.
- **symbols:** φ*(k+1,k) - transpose rigid-body transform from body k+1 to k; H*(k) - hinge map of body k; D^{±1/2}(k) - per-body (block) sqrt / inverse-sqrt of articulated inertia; G*(k) - Kalman-gain-type operator (from K/ψ). <!-- CHECK: G(k) appears only in Table 1; identified from spatial-operator AB algorithm (ref 19/25) as the articulated-body gain; not defined in-text. -->

<!-- eq:34 -->
$$ \langle M^{1/2} \dot{x}\, \dot{x}^* M^{1/2} \rangle = kT\, I_{3n} $$
- **what:** Cartesian velocity-initialization identity (Boltzmann assignment, mean kT/2 per coordinate).
- **symbols:** as eqs:6,29,30.

<!-- eq:35 -->
$$ \dot{\theta} = J^{-1}(\theta)\dot{x} $$
- **what:** BAT velocity initialization by mapping Cartesian velocities (only valid when J is square/invertible, N=3n).
- **symbols:** as eq:28.

<!-- eq:36 -->
$$ \langle \nu\, \nu^* \rangle = \langle M^{1/2} J \dot{\theta}\, \dot{\theta}^* J^* M^{1/2} \rangle = \langle M^{1/2} \dot{x}\, \dot{x}^* M^{1/2} \rangle = kT\, I_{3n} $$
- **what:** Verification that Cartesian-then-BAT initialization satisfies the ICMD equipartition principle eq:26.
- **symbols:** as eqs:30,34,35.

<!-- eq:37 -->
$$ \dot{\theta} = P\dot{x} $$
- **what:** Constrained-ICMD mapping of Cartesian velocities to lower-dimensional ICMD velocity space.
- **symbols:** P ∈ R^{N×3n} - mapping matrix (N < 3n).

<!-- eq:38 -->
$$ \langle \nu\, \nu^* \rangle = m^* P\, P^* m\; kT $$
- **what:** ⟨νν*⟩ under the constrained mapping P; equals kT·I_N iff m*PP*m = I.
- **symbols:** as eqs:15,37.

<!-- eq:39 -->
$$ P = \mathcal{M}^{-1} \mathcal{H} \phi\, M^{1/2} $$
- **what:** A valid mapping matrix giving thermodynamically-rigorous constrained velocity initialization (satisfies m*PP*m=I). Downside: needs the mass-matrix inverse.
- **symbols:** M^{-1} - ICMD inverse mass matrix; H, φ - spatial operators (eq:31); M^{1/2} - sqrt Cartesian mass matrix.

## Derivations (not implemented)

- Full derivation of eq:4 (generalized equipartition via integration by parts,
  Tolman) is in the paper's Supporting Information — not reproduced here.
- The identity `m*PP*m = m*M^{-1}m = I` (verifying eq:39) is pure algebra using
  eq:29, eq:31, eq:15, eq:16 and is not a separate implementable quantity.
