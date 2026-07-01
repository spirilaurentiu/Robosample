# Simbody - implementable equations

Sign/unit convention: SI throughout (lengths m, forces N, moduli Pa). Generalized
coordinates `q`, generalized speeds `u`, related by a block-diagonal kinematic map `N`.
Constraint count `m`, generalized-coordinate count `n`; the multibody solver is O(n)
recursive, with an O(m^3) constraint term.

## Time-stepper view (DAE / DEM)

<!-- eq:1 -->
$$ \dot{y} = f(d; t, y) $$
- **what:** ODE for the continuous state vector y (positions, speeds, aux vars), with fixed discrete state d.
- **symbols:** y - continuous state vector = {q, u, z}; d - discrete state variables; t - time; f - drift function.

<!-- eq:2 -->
$$ 0 = c(d; t, y) $$
- **what:** algebraic constraint equations (loop closure, couplers, prescribed motion, contact).
- **symbols:** c - constraint residual vector; comprises p, pdot, v (position + velocity constraints).

<!-- eq:3 -->
$$ 0 = e(d; t, y) $$
- **what:** event trigger functions; each changes sign (crosses zero) when an event occurs.
- **symbols:** e - vector of event witness functions.

<!-- eq:4 -->
$$ c = 0 \Rightarrow \dot{c} = 0 $$
- **what:** differential-equation-on-a-manifold (DEM) condition; ODE (1) must preserve constraint derivatives so a consistent start stays on the manifold.
- **symbols:** c - constraint residual; cdot - its time derivative.

## Accuracy control

<!-- eq:5 -->
$$ \alpha = 10^{-n} $$
- **what:** map a requested number of accuracy digits n to the single scalar accuracy knob alpha (≈ relative error).
- **symbols:** alpha - accuracy request (dimensionless, ~relative tolerance); n - desired significant digits.

<!-- eq:6 -->
$$ \left\| \mathbf{W}\, \varepsilon_{y} \right\|_{RMS} \le \alpha $$
- **what:** step-acceptance test on weighted state error (RMS norm ≤ accuracy).
- **symbols:** W - diagonal weighting matrix mapping each error to unit error; epsilon_y - integrator estimate of absolute error in each element of y; alpha - accuracy request.

<!-- eq:7 -->
$$ \left\| \mathbf{T}\, c(d; t, y) \right\|_{RMS} \le \alpha $$
- **what:** step-acceptance test on weighted constraint violation.
- **symbols:** T - diagonal weighting matrix for constraint errors; c - constraint error function (Eq. 2).

## Multibody kinematics

<!-- eq:8 -->
$$ \dot{q}_i = \mathbf{N}_i(q_i)\, u_i $$
- **what:** per-mobilizer kinematic differential equation relating coordinate rates to generalized speeds.
- **symbols:** q_i - generalized coordinates of mobilizer i (nq_i of them); u_i - generalized speeds of mobilizer i (nu_i of them); N_i - nq_i x nu_i map (invertible when nq_i = nu_i, e.g. non-quaternion).

<!-- eq:9 -->
$$ \dot{q} = \mathbf{N}(q)\, u $$
- **what:** whole-system kinematic differential equation; N is block-diagonal over mobilizers.
- **symbols:** q - full generalized-coordinate vector; u - full generalized-speed vector; N - block-diagonal kinematic map.

## Constraints (position / velocity / acceleration)

<!-- eq:10 -->
$$ p(t; q) = 0 $$
- **what:** holonomic (position) constraints defining the position manifold restricting q.
- **symbols:** p - holonomic constraint residual.

<!-- eq:11 -->
$$ v(t, q; u) = 0 $$
- **what:** nonholonomic (velocity) constraints.
- **symbols:** v - nonholonomic constraint residual.

<!-- eq:12 -->
$$ a(t, q, u; \dot{u}) = \mathbf{A}\dot{u} - \mathbf{b}_{a} = 0 $$
- **what:** acceleration-only constraints, linear in the generalized accelerations.
- **symbols:** A - acceleration constraint Jacobian block; b_a - residual not depending on udot; udot - generalized acceleration.

<!-- eq:13 -->
$$ \dot{p} = \mathbf{P}u + \frac{\partial p}{\partial t} = 0 $$
- **what:** time derivative of holonomic constraints (velocity-level); P = (dp/dq) N.
- **symbols:** P - holonomic velocity Jacobian = (∂p/∂q) N; ∂p/∂t - explicit time dependence.

<!-- eq:14 -->
$$ \ddot{p} = \mathbf{P}\dot{u} - \mathbf{b}_{p} = 0 $$
- **what:** second time derivative of holonomic constraints (acceleration-level).
- **symbols:** b_p - holonomic acceleration residual (terms independent of udot).

<!-- eq:15 -->
$$ \dot{v} = \mathbf{V}\dot{u} - \mathbf{b}_{v} = 0 $$
- **what:** time derivative of nonholonomic constraints; V = ∂v/∂u.
- **symbols:** V - nonholonomic Jacobian = ∂v/∂u; b_v - residual independent of udot.

<!-- eq:16 -->
$$ g(t, q, u; \dot{u}) = \mathbf{G}\dot{u} - \mathbf{b} = 0 $$
- **what:** all acceleration constraints stacked; G is the acceleration constraint Jacobian (may be rank-deficient).
- **symbols:** G (m x n) = [P V A]^T; b (m x 1) = [b_p b_v b_a]^T; udot - generalized acceleration.

## Dynamics (equations of motion)

<!-- eq:17 -->
$$ \mathbf{M}\dot{u} + \mathbf{G}^{\mathrm{T}}\lambda = \mathbf{f}_{applied} - \mathbf{f}_{inertial} $$
- **what:** constrained equations of motion in internal (generalized) coordinates.
- **symbols:** M(q) - n x n symmetric positive-definite mass matrix; lambda - m Lagrange multipliers (constraint forces); f_applied(t,q,u,z) - applied generalized forces (incl. gravity, motor torques); f_inertial(q,u) - Coriolis/gyroscopic generalized forces.

<!-- eq:18 -->
$$ \dot{z} = \dot{z}(t, q, u, z) $$
- **what:** auxiliary first-order ODEs (e.g. muscle dynamics, controller states).
- **symbols:** z - auxiliary continuous state variables.

## Solving for accelerations (O(n) recursive + O(m^3) constraint)

<!-- eq:19 -->
$$ \mathbf{G}\mathbf{M}^{-1}\mathbf{G}^{T}\lambda = \mathbf{G}\dot{u}_{0} - \mathbf{b} $$
- **what:** eliminate udot from (16)+(17) to get a linear system for the multipliers lambda.
- **symbols:** udot_0 = M^{-1}(f_applied - f_inertial) unconstrained acceleration; RHS = g_0 = acceleration constraint errors of the unconstrained system.

<!-- eq:udot0 -->
$$ \dot{u}_0 = \mathbf{M}^{-1}(\mathbf{f}_{applied} - \mathbf{f}_{inertial}) $$
- **what:** unconstrained generalized acceleration; the M^{-1} v operator is available in O(n) time.
- **symbols:** udot_0 - unconstrained acceleration vector.

<!-- eq:20 -->
$$ \mathbf{Y}\lambda = g_0 $$
- **what:** condensed constraint system; Y = G M^{-1} G^T built column-by-column, total cost O(mn + m^2).
- **symbols:** Y (m x m) = G M^{-1} G^T; g_0 - unconstrained acceleration constraint errors.

<!-- eq:21 -->
$$ \lambda = \mathbf{Y}^{+} g_0 $$
- **what:** least-squares multiplier solution for redundant/singular Y (pseudoinverse via complete orthogonal factorization QTZ, ~5x faster than SVD).
- **symbols:** Y^+ - Moore-Penrose pseudoinverse of Y; g_0 - RHS from Eq. 19.

<!-- eq:udot-final -->
$$ \dot{u} = \mathbf{M}^{-1}(\mathbf{f}_{applied} - \mathbf{f}_{inertial} - \mathbf{f}_{constraint}), \quad \mathbf{f}_{constraint} = \mathbf{G}^{T}\lambda $$
- **what:** final constrained acceleration after multipliers known; one more O(n) M^{-1} application.
- **symbols:** f_constraint = G^T lambda - generalized constraint force.

## Compliant contact (Hertz / Hunt-Crossley / Stribeck)

<!-- eq:22 -->
$$ \mathbf{f}_{contact} = \mathbf{f}_{stiffness} + \mathbf{f}_{dissipation} + \mathbf{f}_{friction} $$
- **what:** contact force = stiffness + dissipation + friction contributions.
- **symbols:** each term a spatial force at contact point P.

<!-- eq:23 -->
$$ f_{stiffness} = f_{Hz} = \left(\tfrac{4}{3}\,\sigma\, R^{1/2} E^*\right) x^{3/2} $$
- **what:** Hertz normal stiffness force; nonlinear in deformation x even for linear-elastic materials.
- **symbols:** x - normal deformation (m, x>0 when contacting); R - composite relative radius of curvature; E* - composite plane-strain modulus (Pa); sigma - eccentricity factor (=1 for circular contact, else Eq. 28).

<!-- eq:24 -->
$$ f_{HC} = \tfrac{3}{2}\, f_{Hz}\, c^*\, \dot{x} $$
- **what:** Hunt-Crossley dissipation force (signed; gives hysteresis).
- **symbols:** c* - effective dissipation coefficient (s/m); xdot - normal approach velocity; f_Hz - Hertz stiffness force.

<!-- eq:25 -->
$$ f_{dissipation} = \max(f_{HC}, -f_{Hz}) $$
- **what:** clamp so total normal force f_Hz + f_dissipation stays >= 0 (no pulling).
- **symbols:** f_HC - Hunt-Crossley force (Eq. 24); f_Hz - Hertz force (Eq. 23).

<!-- eq:26 -->
$$ f_{normal} = f_{stiffness} + f_{dissipation} $$
- **what:** total normal contact force.
- **symbols:** as above.

<!-- eq:27 -->
$$ f_{friction} = \mu(v)\, f_{normal} $$
- **what:** Stribeck friction magnitude; mu depends only on slip speed, parameterized by static/dynamic/viscous coefficients and a transition speed.
- **symbols:** v - slip speed = |v| in contact plane; mu(v) - effective friction coefficient (C2 spline: two quintic segments + sliding); f_normal - Eq. 26.

## Hertz combining rules and contact-point split (Appendix A)

<!-- eq:Estar-material -->
$$ E_i^* = E_i / (1 - \nu_i^2) $$
- **what:** plane-strain modulus of material i.
- **symbols:** E_i - Young's modulus (Pa); nu_i - Poisson ratio.

<!-- eq:combining -->
$$ R = \frac{R_1 R_2}{R_1 + R_2}, \quad E^* = \left(\frac{E_1^{*2/3} E_2^{*2/3}}{E_1^{*2/3} + E_2^{*2/3}}\right)^{3/2}, \quad s_1 = \frac{E_2^{*2/3}}{E_1^{*2/3} + E_2^{*2/3}}, \quad s_2 = 1 - s_1 $$
- **what:** Hertz combining rules for relative curvature R, composite modulus E*, and deformation split fractions s1,s2. <!-- CHECK: paper writes s1 with E_2^{2/3} (star omitted in OCR); the derivation uses starred plane-strain moduli E_i^* throughout, so s1 = E_2^{*2/3}/(E_1^{*2/3}+E_2^{*2/3}). -->
- **symbols:** R_i - principal radius of body i; E_i* - plane-strain modulus; s_i - fraction of total deformation borne by body i.

<!-- eq:x-split -->
$$ x_1 = s_1 x, \quad x_2 = s_2 x $$
- **what:** split total deformation x between the two bodies to locate contact point P.
- **symbols:** x - total normal deformation; x_1,x_2 - per-body deformations.

<!-- eq:cstar -->
$$ c^* \dot{x} = c_1 \dot{x}_1 + c_2 \dot{x}_2 \Rightarrow c^* = c_1 s_1 + c_2 s_2 $$
- **what:** effective dissipation coefficient from per-material coefficients weighted by deformation split.
- **symbols:** c_i - dissipation coefficient of material i; s_i - split fractions.

Note: for the **Elastic Foundation Model** (linear elements) the standard combining
rule `E* = E1* E2* / (E1* + E2*)` (and corresponding s1, s2) is used instead of the
2/3-power Hertz rule above.

## Elliptical Hertz contact correction (Appendix A.2)

<!-- eq:28 -->
$$ \sigma = \frac{\pi\, E(m)^{1/2}\, k}{2\, K(m)^{3/2}} $$
- **what:** eccentricity correction factor for elliptical contact patch.
- **symbols:** k = a/b >= 1 (semi-major/semi-minor axis ratio); m = 1 - (1/k)^2; E(m), K(m) - complete elliptic integrals of the second and first kind respectively. <!-- CHECK: prose says "E and K are complete elliptic integrals of the first and second kinds, resp." but standard notation has K = first kind, E = second kind; the formula uses E and K as second/first kind respectively. -->

<!-- eq:29 -->
$$ \frac{B}{A} = \frac{k^2 E(m) - K(m)}{K(m) - E(m)} $$
- **what:** relates ellipse axis ratio k to the principal semi-curvatures A,B of the separation paraboloid; solve numerically for k (or use [70] approximations for 5-decimal accuracy).
- **symbols:** A, B - principal semi-curvatures of separation paraboloid, B >= A; k = a/b; E(m), K(m) - elliptic integrals.
