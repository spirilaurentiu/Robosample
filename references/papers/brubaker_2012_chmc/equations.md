# Equations - brubaker_2012_chmc (Constrained HMC on manifolds)

Note: appendix proof-algebra steps (marginalization identity, discrete
Euler-Lagrange derivations, partial-integration steps for Theorem 2) are kept in
`paper.md` under "Derivation (not implemented)" and are not reproduced here.

<!-- eq:1 -->
$$\mathcal{H}(p,q) = T(p,q) + U(q) + \lambda^T c(q)$$
- **what:** constrained Hamiltonian = kinetic + potential + Lagrange-multiplier constraint term.
- **symbols:** $p$ - momentum ($\mathbb{R}^n$); $q$ - configuration on manifold ($\mathbb{R}^n$, $c(q)=0$); $T(p,q)$ - kinetic energy; $U(q)$ - potential energy; $c(q)$ - constraint function ($\mathbb{R}^n\to\mathbb{R}^m$); $\lambda$ - Lagrange multipliers ($\mathbb{R}^m$).

<!-- eq:kinetic -->
$$T(p,q) = \tfrac{1}{2}\,p^T M(q)^{-1} p$$
- **what:** kinetic energy in momentum form (used by both acceptance and guidance Hamiltonians).
- **symbols:** $M(q)$ - symmetric positive-definite mass matrix ($n\times n$), may depend on state $q$.

<!-- eq:potential -->
$$U(q) = \tfrac{1}{2}\log|M(q)| - \log\pi(q)$$
- **what:** potential energy of the ACCEPTANCE Hamiltonian; the $\log|M(q)|$ term corrects the state-dependent Gaussian normalization.
- **symbols:** $|M(q)|$ - determinant of mass matrix; $\pi(q)$ - target (unnormalized) density on manifold. For constant $M$ the log-det term is a constant and drops out.

<!-- eq:augmented -->
$$\exp(-\mathcal{H}(p,q)) = \pi(q)\,\mathcal{N}(p\,|\,0,M(q))$$
- **what:** the augmented (joint) distribution factorizes into target times Gaussian momentum; marginal over $p$ recovers $\pi(q)$.
- **symbols:** $\mathcal{N}(\cdot|\mu,\Sigma)$ - multivariate Gaussian density, mean $\mu$, covariance $\Sigma$.

<!-- eq:ham-dynamics -->
$$\dot{p} = -\frac{\partial \mathcal{H}}{\partial q} , \qquad \dot{q} = \frac{\partial \mathcal{H}}{\partial p} , \qquad c(q) = 0$$
- **what:** constrained Hamilton equations of motion (continuous dynamics being integrated).
- **symbols:** $\dot{p},\dot{q}$ - time derivatives; constraint $c(q)=0$ enforced throughout.

<!-- eq:symplectic -->
$$F(x)^T J F(x) = J , \qquad J = \begin{bmatrix} 0 & \mathbf{I}_{n\times n} \\ -\mathbf{I}_{n\times n} & 0 \end{bmatrix}$$
- **what:** symplecticity condition for a map $f:\mathbb{R}^{2n}\to\mathbb{R}^{2n}$; implies volume preservation ($\det F(x)^2 = 1$). Required of the integrator.
- **symbols:** $x=(p,q)$; $F = \partial f/\partial x$ - Jacobian ($2n\times 2n$); $J$ - canonical symplectic matrix.

<!-- eq:momentum-resample -->
$$p_0 \sim \mathcal{N}(0, M(q_0)) \ \ \text{subject to} \ \ C(q_0)\,M(q_0)^{-1} p_0 = 0$$
- **what:** momentum refresh step (Gibbs on momentum): sample unconstrained Gaussian then project onto the cotangent space $\mathcal{T}^*_{q_0}\mathcal{M}$.
- **symbols:** $C(q)=\partial c/\partial q$ - constraint Jacobian ($m\times n$, full rank); note $M^{-1}p = \partial\mathcal{H}/\partial p = \dot q$, so this enforces $C(q)\dot q = 0$.

<!-- eq:accept -->
$$a = \min\!\left(1,\ \exp\{\mathcal{H}(p_0,q_0) - \mathcal{H}(p_L,q_L)\}\right)$$
- **what:** Metropolis acceptance probability for the CHMC proposal $q_L$ after $L$ integration steps.
- **symbols:** $(p_0,q_0)$ - start of trajectory; $(p_L,q_L)$ - end after $L$ steps of size $h$; accept $q_L$ if $u\sim U(0,1)\le a$, else keep $q_0$.

<!-- eq:rattle -->
$$\begin{split}
p_{1/2} &= p_0 - \frac{h}{2} \left( \frac{\partial \hat{\mathcal{H}}(p_{1/2}, q_0)}{\partial q} + C(q_0)^T \lambda \right) , \\
q_1 &= q_0 + \frac{h}{2} \left( \frac{\partial \hat{\mathcal{H}}(p_{1/2}, q_0)}{\partial p} + \frac{\partial \hat{\mathcal{H}}(p_{1/2}, q_1)}{\partial p} \right) , \\
0 &= c(q_1) , \\
p_1 &= p_{1/2} - \frac{h}{2} \left( \frac{\partial \hat{\mathcal{H}}(p_{1/2}, q_1)}{\partial q} + C(q_1)^T \mu \right) , \\
0 &= C(q_1) \frac{\partial \hat{\mathcal{H}}(p_1, q_1)}{\partial p} .
\end{split}$$
- **what:** one step of the generalized RATTLE integrator $\Phi_h^{\hat{\mathcal{H}}}$ (constrained leapfrog); solve for unknowns $p_{1/2}, q_1, p_1, \lambda, \mu$ with Newton's method. Symplectic, symmetric, order 2. First 3 equations are independent of $(p_1,\mu)$; last 2 are linear in $(p_1,\mu)$ for quadratic kinetic energy.
- **symbols:** $h$ - step size; $\hat{\mathcal{H}}$ - guidance/simulation Hamiltonian; $\lambda$ - state-constraint multiplier at end of step; $\mu$ - momentum-constraint multiplier; $p_{1/2}$ - half-step momentum; $\partial\hat{\mathcal{H}}/\partial p = M^{-1}p$ for quadratic kinetic energy.

<!-- eq:2 -->
$$\mathcal{L}'_h(q_0, q_1) = \int_0^h \mathcal{L}(q(t), \dot{q}(t))\, dt + h^r e_h(q_0, q_1)$$
- **what:** consistency/order condition: discrete Lagrangian equals the exact action integral plus $O(h^r)$ error; integrator is consistent if $r\ge 1$ (RATTLE has $r=2$).
- **symbols:** $\mathcal{L}'_h:\mathcal{M}\times\mathcal{M}\to\mathbb{R}$ - discrete Lagrangian; $\mathcal{L}=T-U-\lambda^Tc$ - continuous Lagrangian; $q(t)$ - Euler-Lagrange solution with $q(0)=q_0,q(h)=q_1$; $r$ - integrator order; $e_h$ - bounded error function.

<!-- eq:bvmf -->
$$\pi(q) \propto \exp\!\left(d^T q + q^T A q\right), \qquad \mathbb{S}^{n-1} = \{q \in \mathbb{R}^n \mid q^T q = 1\}$$
- **what:** Bingham-von Mises-Fisher target density on the unit sphere (test distribution). $d=0$ gives Bingham; $A=0$ gives von Mises-Fisher.
- **symbols:** $d$ - location/asymmetry vector ($\mathbb{R}^n$); $A$ - symmetric spread/concentration matrix ($n\times n$), defined up to $A\mapsto A+\alpha I$ (same distribution on sphere). Constraint here $c(q) = q^Tq - 1$.

<!-- eq:collab -->
$$\pi(q) \propto \prod_{(i,j) \in \mathcal{E}} \exp\!\left(-\left(f\!\left(\mathbf{U}_i^T \mathbf{S} \mathbf{V}^j\right) - \mathbf{Y}_{i,j}\right)^2 / \sigma_p^2\right)$$
- **what:** collaborative-filtering target: Gaussian likelihood over observed matrix entries under an orthonormal low-rank factorization $\mathbf{Y}=f(\mathbf{U}^T\mathbf{S}\mathbf{V})$.
- **symbols:** $\mathcal{E}$ - set of observed $(i,j)$ index pairs; $\mathbf{U}_i$ - $i$-th column of $\mathbf{U}$ ($r\times N$); $\mathbf{V}^j$ - $j$-th row of $\mathbf{V}$ ($r\times M$); $\mathbf{S}$ - diagonal ($r\times r$); $f$ - identity or logistic link; $\sigma_p$ - expected prediction error; constraints $\mathbf{U}\mathbf{U}^T=\mathbf{V}\mathbf{V}^T=\mathbf{I}_{r\times r}$.

<!-- eq:5 -->
$$\|q^i - q^j\|_2^2 = l_{i,j}^2 , \qquad \forall (i,j) \in \mathcal{J}$$
- **what:** fixed-limb-length constraints defining the human-pose manifold.
- **symbols:** $q^i\in\mathbb{R}^3$ - 3D position of joint $i$; $l_{i,j}$ - known limb length between joints $i,j$; $\mathcal{J}$ - set of limbs (index pairs).

<!-- eq:6 -->
$$\pi(q) \propto \prod_{i=1}^{N} \exp\!\left(-\|\hat{x}^{i}(q^{i}) - x^{i}\|^{2} / \sigma_{m}^{2}\right) \cdot \prod_{j=1}^{3N} \exp\!\left(-\left(\mathbf{P}_{j}^{T}(q - q_{0})\right)^{2} / \sigma_{j}^{2}\right)$$
- **what:** pose-estimation target: reprojection-error likelihood times PCA-based Gaussian prior over pose.
- **symbols:** $N$ - number of joints; $\hat{x}^i(q^i)$ - projected 2D position of joint $i$ (eq:7); $x^i\in\mathbb{R}^2$ - observed 2D joint location; $\sigma_m^2$ - image-measurement variance; $\mathbf{P}_j$ - $j$-th eigenpose (PCA column vector); $\sigma_j^2$ - corresponding eigenvalue; $q_0$ - mean training pose.

<!-- eq:7 -->
$$\hat{x}^{i}(q^{i}) = \begin{pmatrix} (\mathbf{A}^{1}q^{i})/(\mathbf{A}^{3}q^{i}) \\ (\mathbf{A}^{2}q^{i})/(\mathbf{A}^{3}q^{i}) \end{pmatrix}$$
- **what:** perspective projection of 3D joint $i$ to image coordinates (camera assumed calibrated, pose in camera frame).
- **symbols:** $\mathbf{A}$ - internal camera parameter matrix; $\mathbf{A}^k$ - $k$-th row of $\mathbf{A}$; $q^i\in\mathbb{R}^3$ - 3D joint position (homogeneous division by the third-row projection).

<!-- eq:detailed-balance -->
$$\int_{Q'} \int_{Q} \pi(q) T(q \to q')\, dq\, dq' = \int_{Q} \int_{Q'} \pi(q') T(q' \to q)\, dq'\, dq$$
- **what:** detailed-balance condition CHMC satisfies (Theorem 1) - the correctness guarantee to test a sampler against.
- **symbols:** $T(q\to q')$ - transition kernel; $Q,Q'\subset\mathcal{M}$ - measurable regions.
