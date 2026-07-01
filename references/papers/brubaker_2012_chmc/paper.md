# A Family of MCMC Methods on Implicitly Defined Manifolds

Marcus A Brubaker, Mathieu Salzmann, Raquel Urtasun (AISTATS 2012)

## Abstract

Traditional MCMC methods are only applicable to distributions defined on
$\mathbb{R}^n$. However, there exist many application domains where the
distributions cannot easily be defined on a Euclidean space. To address this
limitation, we propose a general constrained version of Hamiltonian Monte Carlo
(CHMC), and give conditions under which the Markov chain is convergent. Based on
this general framework we define a family of MCMC methods which can be applied to
sample from distributions on non-linear manifolds. We demonstrate the
effectiveness of our approach on sampling from the Bingham-von Mises-Fisher
distribution, collaborative filtering and pose estimation.

## 1 Introduction

Markov Chain Monte Carlo (MCMC) samples complex distributions requiring only the
ability to evaluate the unnormalized probability density. Traditional MCMC targets
distributions defined on $\mathbb{R}^n$, but distributions over non-Euclidean
spaces arise in many domains (protein conformation modelling with the
Fisher-Bingham distribution, texture analysis using distributions over rotations,
fixed-rank matrix factorization for collaborative filtering).

Prior rejection/Gibbs approaches exist for specific manifolds (Gaussian-like
distributions over unit-length vectors and orthogonal matrices) but are highly
specific. In computational physics and molecular dynamics, schemes to sample from
constrained systems were proposed. This paper constructs a family of MCMC methods
for distributions defined on manifolds: it derives the Constrained Hamiltonian
Monte Carlo (CHMC) algorithm which generalizes several HMC methods, and defines a
Constrained Metropolis Monte Carlo sampler for sampling on constrained spaces
without gradients of the target posterior.

## 2 Hamiltonian Dynamics

Let $\mathcal{M} = \{q \in \mathbb{R}^n | c(q) = 0\}$ be a connected,
differentiable submanifold of $\mathbb{R}^n$, where $C(q) = \frac{\partial c}{\partial q}$
is the Jacobian of the constraints, assumed to have full rank everywhere. The
tangent bundle of $\mathcal{M}$ is
$\mathcal{T}\mathcal{M} = \{(q,\dot{q}) | c(q) = 0 \text{ and } C(q)\dot{q} = 0\}$.
The tangent space of $\mathcal{M}$ at a point $q \in \mathcal{M}$ is
$\mathcal{T}_{q}\mathcal{M} = \{\dot{q} | C(q)\dot{q} = 0\}$, with the constraint
$C(q)\dot{q} = 0$ found by differentiating the constraint $c(q)$.

The Lagrangian $\mathcal{L}:\mathcal{TM}\to\mathbb{R}$ of a constrained mechanical
system is the difference between kinetic and potential energies,
$\mathcal{L} = T - U - \lambda^T c(q)$, where $\lambda$ is a vector of Lagrange
multipliers. We assume the potential energy $U(q)$ does not depend on velocity,
and the kinetic energy has the form
$T(q,\dot{q}) = \frac{1}{2}\dot{q}^T M(q)\dot{q}$, where $M(q)$ is a symmetric
positive definite mass matrix. The dynamics are the solution to the Euler-Lagrange
equation $\frac{\partial \mathcal{L}}{\partial q} = \frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{q}}$
coupled with the constraint $c(q) = 0$.

The Hamiltonian $\mathcal{H}: \mathcal{T}^*\mathcal{M} \to \mathbb{R}$ is
constructed by defining the momentum $p = \frac{\partial \mathcal{L}}{\partial \dot{q}}$
and taking the Legendre transformation, $\mathcal{H} = p^T \dot{q} - \mathcal{L}$,
giving

<!-- eq:1 -->
$$\mathcal{H}(p,q) = T(p,q) + U(q) + \lambda^T c(q) , \qquad (1)$$

where the kinetic energy $T(p,q) = \frac{1}{2}p^T M(q)^{-1} p$ is now defined in
terms of the momentum. The cotangent bundle is
$\mathcal{T}^*\mathcal{M} = \{(p,q) | c(q)=0 \text{ and } C(q)\frac{\partial \mathcal{H}}{\partial p}(p,q)=0\}$,
and the cotangent space at $q$ is
$\mathcal{T}_q^*\mathcal{M} = \{p | C(q)\frac{\partial \mathcal{H}}{\partial p}(p,q)=0\}$.
The dynamics in terms of the Hamiltonian are

<!-- eq:ham-dynamics -->
$$\dot{p} = -\frac{\partial \mathcal{H}}{\partial q} , \qquad
\dot{q} = \frac{\partial \mathcal{H}}{\partial p} , \qquad
c(q) = 0 .$$

Properties of Hamiltonian dynamics: (1) symmetry - forward simulation can be
inverted by reversing the direction of time; (2) $\rho$-reversibility with respect
to the map $\rho(p,q)=(-p,q)$; (3) conservation of the Hamiltonian,
$d\mathcal{H}/dt=0$; (4) symplecticity. A mapping
$f:\mathbb{R}^{2n}\to\mathbb{R}^{2n}$ is symplectic if $F(x)^T J F(x)=J$ where
$F=\frac{\partial f}{\partial x}$ is the Jacobian of $f$ and
$J=\begin{bmatrix} 0 & \mathbf{I}_{n\times n} \\ -\mathbf{I}_{n\times n} & 0 \end{bmatrix}$.

Symplectic dynamics imply volume preservation of the cotangent bundle: since
$F(x)^T J F(x) = J$ implies $\det(F(x))^2 = 1$.

## 3 Constrained HMC

Let $\pi: \mathbb{R}^n \to \mathbb{R}$ be a continuous probability measure on
$\mathcal{M}$ such that $\int_{\mathcal{M}} \pi(q) dq = 1$ and $\pi(q) \geq 0$ for
all $q \in \mathcal{M}$. We define two Hamiltonians: the **acceptance Hamiltonian**
$\mathcal{H}$ (used to compute the acceptance probability, based on the target
distribution $\pi(q)$) and the **guidance Hamiltonian** $\hat{\mathcal{H}}$ (used
for simulation).

For the acceptance Hamiltonian, the kinetic energy is
$T(p,q) = \frac{1}{2}p^T M(q)^{-1} p$ and the potential energy is
$U(q) = \frac{1}{2}\log|M(q)| - \log\pi(q)$. This choice means

<!-- eq:augmented -->
$$\exp(-\mathcal{H}(p,q)) = \pi(q)\,\mathcal{N}(p|0,M(q)) ,$$

where $\mathcal{N}(\cdot|\mu,\Sigma)$ is a multivariate Gaussian density. The
guidance Hamiltonian $\hat{\mathcal{H}}$ uses the same kinetic energy but can vary
in its potential energy function $\hat{U}(q)$.

A step of CHMC: (1) draw a new momentum from $\mathcal{N}(p_0|0,M(q_0))$ subject to
$C(q)\frac{\partial \mathcal{H}}{\partial p}(p_0,q_0)=0$ (sample the unconstrained
Gaussian, then project onto $\mathcal{T}_{q_0}^*\mathcal{M}$ - a Gibbs sampler of
the momentum for the augmented distribution $\exp(-\mathcal{H}(p,q))$);
(2) starting at $(p_0,q_0)$, simulate the guidance Hamiltonian for $L$ steps with
step size $h$ ending at $(p_L,q_L)$; (3) accept $q_L$ with probability
$\min(1, \exp\{\mathcal{H}(p_0,q_0) - \mathcal{H}(p_L,q_L)\})$, otherwise keep $q_0$.

### Algorithm 1 Constrained Hamiltonian Monte Carlo

```
Input: q_0, M(q), h, L, H(p,q), Hhat(p,q)
p_0 ~ N(0, M(q_0) | C(q_0) M(q_0)^{-1} p_0 = 0)
for i = 1, ..., L do
    (p_i, q_i) <- Phi_h^Hhat(p_{i-1}, q_{i-1})
end for
u ~ U(0, 1)
if u <= min(1, exp{ H(p_0, q_0) - H(p_L, q_L) }) then
    return q_L   {accept the proposal}
else
    return q_0   {reject the proposal}
end if
```

### 3.1 Numerical Simulation

We denote by $\Phi_h^{\mathcal{H}}: \mathcal{T}^*\mathcal{M} \to \mathcal{T}^*\mathcal{M}$
the numerical integrator which approximates the dynamics of $\mathcal{H}$ a time
$h$ into the future. It must satisfy the state constraints $c(q) = 0$ and the
momentum constraints $C(q) \frac{\partial \mathcal{H}}{\partial p}(p_0, q_0) = 0$.
We further require the integrator to be **symmetric**
($(p,q) = \Phi_{-h}^{\mathcal{H}}(\Phi_h^{\mathcal{H}}(p,q))$) and **symplectic**.

A symplectic integrator implies there exists a discrete Lagrangian
$\mathcal{L}'_h$ which the numerical method integrates: starting at $(p_0, q_0)$,
$q_1$ is the solution to
$p_0 + \frac{\partial}{\partial q_0} \mathcal{L}'_h(q_0, q_1) = \lambda_0^T C(q_0)$.

The integrator must also be **consistent**. A symplectic integrator of order $r$
satisfies

<!-- eq:2 -->
$$\mathcal{L}'_h(q_0, q_1) = \int_0^h \mathcal{L}(q(t), \dot{q}(t)) dt + h^r e_h(q_0, q_1) , \quad (2)$$

where $\mathcal{L}'_h: \mathcal{M} \times \mathcal{M} \to \mathbb{R}$ is the
discrete Lagrangian, $q(t)$ is the solution to the Euler-Lagrange equation with
boundary conditions $q(0) = q_0$, $q(h) = q_1$, and $e_h$ is a bounded error
function. An integrator is consistent if $r \geq 1$. Numerical integration does
not exactly conserve the Hamiltonian; this error is corrected by the Metropolis
acceptance test.

One integration method satisfying these conditions is **RATTLE**, a generalization
of the Leapfrog integrator to handle manifold constraints and the more general
Hamiltonian form. A step of the generalized RATTLE algorithm solves the system of
non-linear equations

<!-- eq:rattle -->
$$\begin{split}
p_{1/2} &= p_0 - \frac{h}{2} \left( \frac{\partial \hat{\mathcal{H}}(p_{1/2}, q_0)}{\partial q} + C(q_0)^T \lambda \right) , \\
q_1 &= q_0 + \frac{h}{2} \left( \frac{\partial \hat{\mathcal{H}}(p_{1/2}, q_0)}{\partial p} + \frac{\partial \hat{\mathcal{H}}(p_{1/2}, q_1)}{\partial p} \right) , \\
0 &= c(q_1) , \\
p_1 &= p_{1/2} - \frac{h}{2} \left( \frac{\partial \hat{\mathcal{H}}(p_{1/2}, q_1)}{\partial q} + C(q_1)^T \mu \right) , \\
0 &= C(q_1) \frac{\partial \hat{\mathcal{H}}(p_1, q_1)}{\partial p} ,
\end{split}$$

for the unknowns $p_{1/2}, q_1, p_1, \lambda, \mu$ where $\hat{\mathcal{H}}$ is the
simulation Hamiltonian, and $\lambda$ and $\mu$ are the Lagrange multipliers
associated with the state and momentum constraints at the end of the step. A
solution can be obtained using Newton's method. This method is symplectic,
symmetric, of order 2 (consistent), and respects the manifold constraints,
ensuring the solution lies in $\mathcal{T}^*\mathcal{M}$. It naturally handles a
state-dependent mass matrix. Solving can be done more efficiently by noting that
the first three equations are independent of $p_1, \mu$ and the last two equations
are linear in $p_1, \mu$ for quadratic kinetic energies. Higher-order methods on a
manifold can be obtained e.g. via a partitioned Runge-Kutta method using the
Lobatto IIIA-IIIB pair.

### 3.2 Convergence of CHMC

We prove CHMC converges to the target posterior $\pi(q)$ from any starting point
$q \in \mathcal{M}$.

**Theorem 1 (Detailed Balance).** Let $\hat{\mathcal{H}}$ be $\mathcal{C}^2$
continuous, $\mathcal{M} = \{q \in \mathbb{R}^n | c(q) = 0\}$ be a connected,
smooth and differentiable manifold with $\frac{\partial c}{\partial q}$ full-rank
everywhere, $M(q)$ be positive definite on $\mathcal{M}$, and $\pi(q)$ be smooth.
If $\Phi_h^{\hat{\mathcal{H}}}$ is symmetric and symplectic, then

<!-- eq:3 -->
$$\int_{Q'} \int_{Q} \pi(q) T(q \to q') dq dq' = \int_{Q} \int_{Q'} \pi(q') T(q' \to q) dq' dq , \qquad (3)$$

where $Q, Q' \subset \mathcal{M}$ and $T$ is the transition kernel.

**Theorem 2 (Accessibility).** Under the same assumptions plus consistency of
$\Phi_h^{\hat{\mathcal{H}}}$, for any $q_0, q_1 \in \mathcal{M}$ and $h$
sufficiently small, there exist finite $p_0 \in \mathcal{T}_{q_0}^* \mathcal{M}$,
$p_1 \in \mathcal{T}_{q_1}^* \mathcal{M}$ and Lagrange multipliers
$\lambda_0, \lambda_1$ such that $(p_1, q_1) = \Phi_h^{\hat{\mathcal{H}}}(p_0, q_0)$.

For a single step of RATTLE on a sphere there will generally be two choices of
$\lambda$ corresponding to points on opposite hemispheres. To prove irreducibility
we assume that for small regions around every point on the manifold the integrator
can uniquely choose the Lagrange multiplier which moves between points within the
region.

**Theorem 3 (Irreducibility).** Under the above assumptions plus $\pi(q)$ strictly
positive and the uniqueness-of-Lagrange-multiplier assumption within balls
$\mathcal{B}_{\ell}(q) = \{q' \in \mathcal{M} | d(q', q) \leq \ell\}$, for $h$
sufficiently small and any $q_0, q_1 \in \mathcal{M}$ with $\pi(q_1) > 0$, there
exists $n \in \mathbb{N}$ such that

<!-- eq:4 -->
$$T^n(q_0 \to q_1) > 0. \qquad (4)$$

**Theorem 4 (Convergence).** Under the above assumptions with $\pi(q)$ smooth and
strictly positive on $\mathcal{M}$, for all $q_0 \in \mathcal{M}$

<!-- eq:convergence -->
$$\lim_{n\to\infty} ||T^n(q_0\to\cdot)-\pi(\cdot)||=0.$$

### 3.3 Variations of CHMC

Unconstrained HMC and Riemann Manifold HMC are both instances of CHMC when
$\mathcal{M} = \mathbb{R}^n$.

**Constrained Langevin MC:** a special case of CHMC with $\hat{\mathcal{H}}=\mathcal{H}$
and a single simulation step, $L=1$. Requires fewer gradient evaluations, but can
exhibit random-walk-like behavior.

**Constrained Metropolis Monte Carlo:** the unconstrained Metropolis algorithm
with proposal covariance $\Sigma$ is a special case of HMC with simulation
Hamiltonian $\hat{\mathcal{H}}(p,q) = \frac{1}{2}p^T M^{-1} p$ and parameters $h=1$,
$L=1$, $M=\Sigma^{-1}$. This yields an MCMC method on a constrained space which
does not require the gradient of the target posterior.

**Constrained Riemann Manifold HMC:** the mass matrix $M(q)$ becomes dependent on
state $q$, exploiting geometric information; naturally handled in this framework.

## 4 Experimental Evaluation

### 4.1 Linearly Constrained Gaussian

Sample from a multivariate Gaussian subject to linear equality constraints:
$\pi(q) \propto \mathcal{N}(q|\mu,\Sigma)$ subject to $c(q) = Aq - b$, where
$A \in \mathbb{R}^{D \times n}$ and $b \in \mathbb{R}^D$. Parameters:
$\mu = (0,0,0,0)^T$, $\Sigma = \text{diag}(1,1,0.01,0.01)$, two constraints
$A_1 = (1,1,1,1)$, $A_2 = (1,1,-1,1)$, $b = (0,0)^T$. Chains initialized at
$q = (9, -9, 11, -11)^T$; CHMC and CLangevin converge quickly to the correct value
0. When initialized at the mode $q=(0,0,0,0)^T$, CHMC is more efficient than
CMetropolis and CLangevin (autocorrelation of $-\log\pi(q)$).

### 4.2 Bingham-von Mises-Fisher Distribution

The Bingham-von Mises-Fisher distribution is a Gaussian in $\mathbb{R}^n$
restricted to the unit sphere $\mathbb{S}^{n-1}$. Density:

<!-- eq:bvmf -->
$$\pi(q) \propto \exp\left(d^T q + q^T A q\right) \quad \text{restricted to} \quad \mathbb{S}^{n-1} = \{q \in \mathbb{R}^n | q^T q = 1\}.$$

If $d$ is the zero-vector it reduces to the Bingham distribution; if $A$ is the
zero-matrix it is the von Mises-Fisher distribution. If $d$ is non-zero the
distribution is antipodally asymmetric ($\pi(q) \neq \pi(-q)$) with a bias towards
values pointing in the same direction as $d$. $A$ is not uniquely defined: $A$ and
$A + \alpha I$ describe the same distribution for any scalar $\alpha$, since
$\exp(q^T(A + \alpha I)q) = \exp(q^T A q + \alpha q^T q) \propto \exp(q^T A q)$ on
the sphere.

CLangevin outperformed CHMC here, contrasting with the general belief in
unconstrained MCMC that multiple steps perform better; on closed compact spaces
such as $\mathbb{S}^{n-1}$ the ability to reach distant points is less valuable,
and longer simulations begin to oscillate ($L>2$).

### 4.3 Collaborative Filtering

Given a matrix $\mathbf{Y} \in \mathbb{R}^{N \times M}$ of observed ratings, seek a
low-rank decomposition $\mathbf{Y} = f(\mathbf{U}^T \mathbf{S} \mathbf{V})$ with
$\mathbf{U} \in \mathbb{R}^{r \times N}$, $\mathbf{V} \in \mathbb{R}^{r \times M}$,
and $\mathbf{S}$ a diagonal matrix. $f()$ is often the identity, but the logistic
function is also possible. We sample $\mathbf{U}$, $\mathbf{V}$, and the diagonal
of $\mathbf{S}$ for fixed rank $r$ under orthonormality constraints
$\mathbf{U}\mathbf{U}^T = \mathbf{I}_{r \times r}$ and
$\mathbf{V}\mathbf{V}^T = \mathbf{I}_{r \times r}$. The vector $q$ contains the
concatenated vectorized forms of $\mathbf{U}$, $\mathbf{V}$, and $\text{diag}(\mathbf{S})$.
Given a set $\mathcal{E}$ of index pairs of known entries, the density is

<!-- eq:collab -->
$$\pi(q) \propto \prod_{(i,j) \in \mathcal{E}} \exp\left(-\left(f\left(\mathbf{U}_i^T \mathbf{S} \mathbf{V}^j\right) - \mathbf{Y}_{i,j}\right)^2 / \sigma_p^2\right) ,$$

where $\mathbf{U}_i$ is the $i$-th column of $\mathbf{U}$, $\mathbf{V}^j$ is the
$j$-th row of $\mathbf{V}$, and $\sigma_p$ is the expected prediction error.

Tested on 1M MovieLens and EachMovie under the weak generalization setting (a
single rating per user withheld). NMAE normalization constants: 1.6 for MovieLens,
1.944 for EachMovie. Results are the mean predictions over 2000 samples.

### 4.4 Human Pose Estimation

3D human pose estimation from monocular 2D observations. $q$ is the vector of 3D
coordinates of $N$ joints. The constraints encode fixed limb lengths:

<!-- eq:5 -->
$$\|q^i - q^j\|_2^2 = l_{i,j}^2 , \ \forall (i,j) \in \mathcal{J} , \qquad (5)$$

where $q^i$ encodes the 3D position of joint $i$, $l_{i,j}$ is the known limb
length, and $\mathcal{J}$ is the set of limbs. Given noisy image locations
$x^i \in \mathbb{R}^2$ of the joints and camera parameter matrix $\mathbf{A}$, with
a linear pose model from PCA, the density is

<!-- eq:6 -->
$$\pi(q) \propto \prod_{i=1}^{N} \exp\left(-\|\hat{x}^{i}(q^{i}) - x^{i}\|^{2} / \sigma_{m}^{2}\right) \cdot \prod_{j=1}^{3N} \exp\left(-\left(\mathbf{P}_{j}^{T}(q - q_{0})\right)^{2} / \sigma_{j}^{2}\right) , \qquad (6)$$

where $\mathbf{P}_j$ is the column vector containing the $j$-th eigenpose from PCA,
$\sigma_j^2$ is the corresponding eigenvalue, $q_0$ is the mean pose of the
training data, and

<!-- eq:7 -->
$$\hat{x}^{i}(q^{i}) = \begin{pmatrix} (\mathbf{A}^{1}q^{i})/(\mathbf{A}^{3}q^{i}) \\ (\mathbf{A}^{2}q^{i})/(\mathbf{A}^{3}q^{i}) \end{pmatrix} \qquad (7)$$

is the projection of joint $i$, with $\mathbf{A}^k$ the $k$-th row of $\mathbf{A}$,
and $\sigma_m^2$ the expected variance of the image measurements. The first factor
encodes the reprojection error; the second is the PCA-based regularizer.

Experiments on the walking sequence of subject 1 in HumanEva. One circle used to
learn the linear pose model; noise std 0 to 10 pixels added to projections.
$q$ initialized with a random training pose; $L = 200$ steps; the Hessian of the
negative log of the observations used as a (state-dependent) mass matrix. CHMC
clearly outperforms a projected-gradient-descent constrained-optimization baseline
which tends to get stuck in local optima.

## 5 Discussion

A general framework for constructing Markov chains on manifolds defined by implicit
constraints, yielding a family of samplers: Constrained Hamiltonian Monte Carlo,
Constrained Metropolis and Constrained Langevin. Conditions necessary for
convergence are explicitly stated. Traditional HMC is a special case where
$\mathcal{M} = \mathbb{R}^n$. The framework allows state-dependent mass matrices,
extending Riemann Manifold HMC to constrained configuration spaces, and the method
and proofs are applicable to any appropriate integration method applied to any
guidance Hamiltonian. Matlab code implementing CHMC is available from the authors.

## Derivation (not implemented) - Appendix proofs

**Theorem 1 (Detailed Balance), proof.** CHMC first satisfies detailed balance
with respect to $\bar{\pi}(p,q) = \exp(-\mathcal{H}(p,q))$ on the augmented state
space $\mathcal{T}^*\mathcal{M}$. Let $R,R'\subset\mathcal{T}^*\mathcal{M}$ be
sufficiently small regions such that $\bar{\pi}$ is constant over them, with values
$\bar{\pi}(R)$ and $\bar{\pi}(R')$, and $R'$ the image of $R$ under $T$. Since
$\Phi_h^{\hat{\mathcal{H}}}$ is symmetric, the image of $R'$ under $T$ is $R$ with
the sign of the momentum reversed. Since integration is symplectic, phase-space
volume is preserved, $\int dR = \int dR' = \delta V$. Thus

$$\begin{split}
& \int_{R'} \int_{R} \tfrac{1}{Z_{\bar{\pi}}} \bar{\pi}(r) T(r \to r') dr dr' \\
& = \tfrac{1}{Z_{\bar{\pi}}} \bar{\pi}(R) \delta V \min\left(1, \tfrac{\bar{\pi}(R')}{\bar{\pi}(R)}\right) \\
& = \tfrac{1}{Z_{\bar{\pi}}} \bar{\pi}(R') \delta V \min\left(1, \tfrac{\bar{\pi}(R)}{\bar{\pi}(R')}\right) \\
& = \int_{R} \int_{R'} \tfrac{1}{Z_{\bar{\pi}}} \bar{\pi}(r') T(r' \to r) dr' dr ,
\end{split}$$

where $Z_{\bar{\pi}} = \int_{\mathcal{T}^*\mathcal{M}} \bar{\pi}(r) dr$. Then
$\pi(q)$ is the marginal of $\bar{\pi}(p,q)$:

$$\int \bar{\pi}(p,q)dp = \int_{\mathcal{T}_q^* \mathcal{M}} \pi(q) \frac{1}{Z_{\mathcal{N}}} \mathcal{N}(p|0, M(q)) dp = \pi(q) ,$$

where $Z_{\mathcal{N}} = \int_{\mathcal{T}_q^*\mathcal{M}} \mathcal{N}(p|0,M(q)) dp$
is the normalizing constant of the Gaussian restricted to the tangent space.
Detailed balance with respect to $\pi(q)$ follows by ignoring the momentum.

**Theorem 2 (Accessibility), proof.** If the integrator is symplectic there exists
a discrete Lagrangian $\mathcal{L}'_h(q_0,q_1)$ for which
$p_0 + \frac{\partial \mathcal{L}'_h(q_0,q_1)}{\partial q_0} = \lambda_0^T C(q_0)$.
For given $q_0, q_1$, $p_0 = -\frac{\partial \mathcal{L}'_h(q_0,q_1)}{\partial q_0} + \lambda_0^T C(q_0)$
with $\lambda_0$ chosen so $p_0 \in \mathcal{T}^*_{q_0}\mathcal{M}$; by symmetry
$p_1 = -\frac{\partial \mathcal{L}'_{-h}(q_1,q_0)}{\partial q_1} + \lambda_1^T C(q_1)$.
Such choices exist so long as $C(q_0), M(q_0), C(q_1), M(q_1)$ all have full rank.
Using consistency (Eq. 2) and partial integration one derives
$\frac{\partial \mathcal{L}_h'}{\partial q_0}(q_0,q_1) = -\frac{\partial \mathcal{L}}{\partial \dot{q}}(q_0,\dot{q}_0) + h^r \frac{\partial e_h}{\partial q_0}(q_0,q_1)$,
so the momenta $p_0, p_1$ exist and are finite.

**Theorem 3 (Irreducibility), proof.** Because $\mathcal{M}$ is connected there is
a geodesic between any $q, q'$. Divide it into $n = \lceil d(q,q')/\ell \rceil$
chunks with $d(q_{i-1},q_i) \leq \ell$; by Theorem 2 each transition has non-zero
probability, and $T^n(q \to q') = \prod_{i=1}^{n} T(q_{i-1} \to q_i) > 0$ since
$\pi$ is strictly positive.

**Lemma 1 (Aperiodicity).** Proven by contradiction as an almost direct consequence
of irreducibility.

**Theorem 4 (Convergence), proof.** Since CHMC satisfies detailed balance (Thm 1),
is $\pi$-irreducible (Thm 3), and is aperiodic (Lemma 1), convergence follows by
Theorem 1 of Tierney (1994).
