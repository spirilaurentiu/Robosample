# Equations - Fixman compensating potential for general branched molecules (Jain et al. 2013)

<!-- eq:1 -->
$$ \mathcal{H}(\alpha, \mathbf{q}, \mathbf{p}) = \frac{1}{2} p^* \mathcal{M}_B^{-1}(\alpha, \mathbf{q})\, p + \mathcal{U}(\alpha, \mathbf{q}) $$
- **what:** Hamiltonian of the unconstrained (flexible) BAT-coordinate model: kinetic term with inverse mass matrix plus force-field potential.
- **symbols:** alpha - unconstrained BAT coordinates (torsions), N of them; q - the (3n-N) coordinates to be constrained (bond lengths/angles); p - canonical momenta conjugate to full BAT set; M_B - BAT-coordinate mass matrix (R^{3n x 3n}); U - force-field potential; `*` denotes transpose.

<!-- eq:2 -->
$$ \mathcal{Z}(T) = c_1 \int dp\; d\alpha\; dq\; e^{-\mathcal{H}(\alpha,q,p)/kT} $$
- **what:** full phase-space partition function at temperature T.
- **symbols:** Z(T) - partition function; c_1 - scaling constant; k - Boltzmann constant; T - temperature.

<!-- eq:3 -->
$$ \mathcal{Z}(T) = c_2 \int d\alpha\; dq\; \det\!\big\{ \mathcal{M}_B^{\frac{1}{2}}(\alpha,q) \big\}\, e^{-\mathcal{U}(\alpha,q)/kT} $$
- **what:** configuration-space partition function after Gaussian integration over momenta; introduces the det(M_B^{1/2}) Jacobian factor. This det factor is the origin of the constraint bias.
- **symbols:** c_2 - constant absorbing momentum integral; det{M_B^{1/2}} = det{M_B}^{1/2}.

<!-- eq:4 -->
$$ \rho(\alpha, q) \propto \det\!\left\{ \mathcal{M}_B^{\frac{1}{2}}(\alpha, q) \right\} e^{-\mathcal{U}(\alpha, q)/kT} $$
- **what:** configuration pdf for the unconstrained (flexible) model.
- **symbols:** rho - probability density over configuration coordinates.

<!-- eq:5 -->
$$ \mathcal{M}_B = \mathcal{J}_B^* \, \mathcal{M}_c \, \mathcal{J}_B $$
- **what:** BAT mass matrix built from Cartesian mass matrix pushed through the BAT->Cartesian Jacobian.
- **symbols:** J_B - Jacobian of BAT->Cartesian transform (R^{3n x 3n}); M_c - constant diagonal Cartesian mass matrix (atom masses on diagonal, R^{3n x 3n}). <!-- CHECK: raw shows J_R^* M_c J_B; per Eq.5 context both factors are J_B, so J_R is an OCR typo for J_B -->

<!-- eq:6 -->
$$ \det\{\mathcal{J}_B\} = \sin\theta_{ex}\; d_2^2 \prod_{i=3}^{n} d_i^2 \sin\theta_i $$
- **what:** Go-Scheraga closed form for the BAT->Cartesian Jacobian determinant; note it is independent of torsion angles - depends only on bond lengths and bond angles.
- **symbols:** d_i - the (n-1) bond lengths; theta_i - the (n-2) bond angles; theta_ex - overall molecule orientation coordinate; n - number of atoms.

<!-- eq:7 -->
$$ \det\{\mathcal{M}_B\} = \det\{\mathcal{J}_B\}^2 \prod_{i=1}^{n} \mathfrak{m}_i^3 = \sin^2\theta_{ex}\; d_2^4 \prod_{i=3}^{n} d_i^4 \sin^2\theta_i \prod_{i=1}^{n} \mathfrak{m}_i^3 $$
- **what:** closed-form determinant of the full BAT mass matrix; torsion-independent, cheap to evaluate.
- **symbols:** m_i - mass of the ith atom (raw uses fraktur m); other symbols as in eq:6.

<!-- eq:8 -->
$$ \rho(\beta_i) = \frac{1}{2\pi} $$
- **what:** for the flexible model with only bond/angle potentials, every torsion angle has a uniform pdf on [0, 2pi). Key ground-truth target for validation.
- **symbols:** beta_i - ith torsion angle.

<!-- eq:9 -->
$$ \mathcal{Z}(T) = c_3 \int d\alpha\; \det\!\left\{ \mathcal{M}^{\frac{1}{2}}(\alpha) \right\} e^{-\mathcal{U}(\alpha, q_0)/kT} $$
- **what:** configuration partition function of the constrained model (q frozen at q_0).
- **symbols:** M(alpha) - constrained-model mass matrix (R^{N x N}), the alpha-alpha sub-block of M_B; q_0 - fixed constrained-coordinate values; c_3 - constant.

<!-- eq:10 -->
$$ \rho(\alpha) \propto \det\!\big\{ \mathcal{M}^{\frac{1}{2}}(\alpha) \big\} e^{-\mathcal{U}(\alpha, q_0)/kT} $$
- **what:** pdf of the constrained model; det{M(alpha)} DOES depend on torsions, causing the bias vs eq:4.
- **symbols:** as above.

<!-- eq:11 -->
$$ \mathcal{U}'(\alpha) \triangleq \mathcal{U}(\alpha, q_0) + \mathcal{U}_f(\alpha), \qquad \mathcal{U}_f(\alpha) \triangleq \frac{1}{2} kT \ln \frac{\det\{\mathcal{M}(\alpha)\}}{\det\{\mathcal{M}_B(\alpha, q_0)\}} $$
- **what:** Fixman-corrected potential: add the Fixman compensating potential U_f to the constrained force-field potential so the constrained pdf matches the flexible one.
- **symbols:** U' - modified (corrected) potential; U_f - Fixman compensating potential; "triangle-equals" = "defined as".

<!-- eq:12 -->
$$ \rho(\alpha) \propto \det\!\left\{ \mathcal{M}_B^{\frac{1}{2}}(\alpha, q_0) \right\} e^{-\mathcal{U}(\alpha, q_0)/kT} $$
- **what:** compensated pdf after substituting U' into eq:10; now agrees with flexible pdf eq:4 at q=q_0.
- **symbols:** as above.

<!-- eq:13 -->
$$ \mathcal{U}_f(\alpha) = c_f + \frac{1}{2} \ln \det\{ \mathcal{M}(\alpha) \}, \qquad \rho(\alpha) \propto e^{-\mathcal{U}(\alpha, q_0)/kT} $$
- **what:** simplification when the unconstrained coordinates alpha are only torsion angles: the M_B term becomes a constant, so U_f reduces to a constant plus half the log-det of the constrained mass matrix. This is the working form implemented.
- **symbols:** c_f - constant (absorbs torsion-independent M_B contribution and kT prefactor convention). Note: this drops the kT prefactor of eq:11 into c_f / absorbs it; U_f here is in units where the ln det term has coefficient 1/2. <!-- CHECK: eq:11 has (1/2)kT ln(...) but eq:13 writes (1/2) ln det{M}; the kT factor is implicit/absorbed - verify unit convention when porting -->

<!-- eq:14 -->
$$ \mathcal{T}(i) = -\frac{\partial \mathcal{U}_f(\alpha)}{\partial \alpha_i} $$
- **what:** Fixman torque on coordinate i = negative gradient of the Fixman potential w.r.t. torsion alpha_i.
- **symbols:** T(i) - Fixman torque for ith unconstrained coordinate.

<!-- eq:15 -->
$$ \mathcal{M} = H\phi M\phi^* H^* $$
$$ \mathcal{M} = [I + H\phi \mathcal{K}]\, \mathcal{D}\, [I + H\phi \mathcal{K}]^* $$
$$ [I + H\phi \mathcal{K}]^{-1} = [I - H\psi \mathcal{K}] $$
$$ \mathcal{M}^{-1} = [I - H\psi \mathcal{K}]^*\, \mathcal{D}^{-1}\, [I - H\psi \mathcal{K}] $$
- **what:** Spatial Operator Algebra (SOA) factorizations of the constrained mass matrix: (1) Newton-Euler operator factorization; (2) alternative square-factor factorization via articulated-body (AB) recursion; (3) closed-form inverse of the [I+HphiK] operator; (4) resulting operator factorization of M^{-1} (basis of linear-cost GNEIMO dynamics).
- **symbols:** H - block-diagonal hinge articulation operator (torsional axes per DOF); phi - lower-triangular rigid-body force-propagation operator (phi(i,j) propagates spatial force from cluster j to cluster i); M - block-diagonal link spatial inertia operator; D - block-diagonal AB articulated inertia; K - AB Kalman-gain-like operator; psi = phi(I - K H) propagation with AB feedback; I - identity; `*` transpose.

<!-- eq:16 -->
$$ \det\{\mathcal{M}\} = \det\{I + H\phi \mathcal{K}\}^2 \det\{\mathcal{D}\} = \prod_{i=0}^{\mathcal{N}-6} \det\{\mathcal{D}(i)\} = \det\{\mathcal{D}(0)\} \prod_{i=1}^{\mathcal{N}-6} \mathcal{D}(i) $$
- **what:** determinant of constrained mass matrix reduces to product of AB articulated-inertia block determinants, since det{I+HphiK}=1 and D is block-diagonal. All D(i) except the base D(0) are scalars.
- **symbols:** D(i) - articulated-inertia block for cluster i (scalar for i>=1, a 6x6 matrix for base cluster i=0); N - number of unconstrained coordinates (index runs to N-6 accounting for 6 base DOF). <!-- CHECK: index upper limit N-6 as printed -->

<!-- eq:17 -->
$$ \mathcal{U}_f(\alpha) = c_f + \frac{1}{2} \ln \det\{ \mathcal{D}(0) \} + \frac{1}{2} \sum_{i=1}^{N-6} \ln \mathcal{D}(i) $$
- **what:** GNEIMO-Fixman working formula for the Fixman potential: sum of logs of the scalar articulated inertias D(i) plus the base 6x6 log-det, all available as by-products of the constrained-dynamics AB recursion. Linear cost.
- **symbols:** as eq:16; c_f constant.

<!-- eq:18 -->
$$ \mathcal{T}(i) = -\frac{1}{2} \frac{\partial \ln \det\{\mathcal{M}(\alpha)\}}{\partial \alpha_i} $$
- **what:** Fixman torque as half the negative derivative of ln det of the constrained mass matrix (combine eq:13 and eq:14).
- **symbols:** as above.

<!-- eq:19 -->
$$ \frac{\partial g(X(y))}{\partial y} = \operatorname{Trace}\!\left\{ \left[ \frac{\partial g}{\partial X} \right]^* \frac{\partial X(y)}{\partial y} \right\} $$
- **what:** matrix-calculus chain rule for a scalar function g of a matrix X(y).
- **symbols:** g - smooth scalar function; X - matrix (R^{m x n}); y - scalar parameter; dg/dX and dX/dy are m x n matrices (elementwise, see eq:20).

<!-- eq:20 -->
$$ \frac{\partial g}{\partial X}(i,j) \triangleq \frac{\partial g}{\partial X(i,j)}, \qquad \frac{\partial X}{\partial y}(i,j) \triangleq \frac{\partial X(i,j)}{\partial y} $$
- **what:** elementwise definitions of the matrix derivatives used in eq:19.
- **symbols:** (i,j) matrix element indices.

<!-- eq:21 -->
$$ \frac{\partial \ln \det\{X\}}{\partial X} = \{X^*\}^{-1} $$
- **what:** standard identity: gradient of log-det is inverse-transpose. Used to turn eq:18 into explicit SOA form.
- **symbols:** X - invertible matrix.

<!-- eq:22 -->
$$ \mathcal{T}(i) = -\operatorname{Trace}\{\mathcal{P}(i)\,\Upsilon(i)\,\widetilde{H}_{\omega}^{*}(i)\} $$
- **what:** SOA explicit Fixman torque for the ith torsion; all three factors are 6x6 matrices per cluster.
- **symbols:** P(i) - articulated-body quantity for cluster i (by-product of GNEIMO), 6x6; Upsilon(i) - 6x6 matrix computed by an extra recursive scatter algorithm from the base cluster; H-tilde_omega(i) - rotational part of the hinge map operator, 6x6.

<!-- eq:23 -->
$$ \mathcal{T}(i) = h^*(i)\, \mathcal{F}[Q_{11} + Q_{22}], \qquad \mathcal{P}(i)\Upsilon(i) = \begin{bmatrix} Q_{11} & Q_{12} \\ Q_{21} & Q_{22} \end{bmatrix},\; Q_{ij} \in \mathcal{R}^{3\times3} $$
- **what:** simpler equivalent Fixman-torque form; partition the 6x6 product P(i)Upsilon(i) into 3x3 blocks, take the "uncross" of the sum of diagonal blocks, dot with the hinge axis. Holds for branched systems of any size.
- **symbols:** h(i) - torsion hinge axis (3-vector); F(A) - operator mapping a 3x3 matrix A to the 3-vector v whose cross-product (skew) matrix equals (A - A^*); Q_{jk} - 3x3 sub-blocks of P(i)Upsilon(i).

<!-- eq:24 -->
$$ \det\{\mathcal{M}(\alpha)\} = c_5\big(35 + 4\cos(\alpha) - 16\cos^2(\alpha) + \cos^4(\alpha)\big) $$
- **what:** Pear-Weiner closed-form mass-matrix determinant for the idealized C4 3-bond serial chain with 90-degree bond angles; single torsion alpha. Used as an analytical ground truth to validate eq:17.
- **symbols:** alpha - single torsion angle; c_5 - constant depending on bond lengths and masses. <!-- CHECK: raw is missing a closing paren after cos^4(alpha); restored -->

<!-- eq:25 -->
$$ \rho(\alpha) \propto \det\!\left\{\mathcal{M}^{\frac{1}{2}}(\alpha)\right\} = \left[\det\{\mathcal{D}(0)\} \prod_{i=1}^{N-6} \mathcal{D}(i)\right]^{\frac{1}{2}} $$
- **what:** predicted torsion pdf for the TORSIONAL (constrained, no Fixman) simulations - the biased distribution to be corrected.
- **symbols:** as eq:16/eq:17.

<!-- eq:26 -->
$$ M\ddot{x} = -\nabla_x U - \gamma M\dot{x} + \sqrt{2\gamma kT}\, M^{\frac{1}{2}}\, dW $$
- **what:** Langevin dynamics equation of motion for the flexible Cartesian model (integrated with BBK).
- **symbols:** x - Cartesian coordinates; M - Cartesian mass matrix; gamma - damping coefficient (0.01/fs); U - bond/angle potential; dW - vector of independent Wiener processes; k Boltzmann, T temperature.

## Derivations note

The momentum-integration step turning eq:2 into eq:3 (Gaussian integral yielding det{M_B^{1/2}}), the substitution of U' into eq:10 to recover eq:12, and the algebraic manipulations of eqs:18-21 into the explicit SOA torque eqs:22-23 are derivation steps; full derivations are in Refs. 27 and 35 (Jain 1997; Jain, Robot and Multibody Dynamics 2011) and are not reproduced here.
