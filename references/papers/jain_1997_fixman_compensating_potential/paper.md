# Compensating Mass Matrix Potential for Constrained Molecular Dynamics

**Author:** Abhinandan Jain (Jet Propulsion Laboratory / Caltech)
**Venue:** Journal of Computational Physics, 1997. DOI: 10.1006/jcph.1997.5731

## Abstract

Rigid internal constraints are used in molecular models to speed up molecular
dynamics (MD) simulations. Statistical averages from constrained MD simulations
differ by a metric-tensor-dependent term from averages computed using
conventional unconstrained MD. Fixman proposed augmenting the standard potential
with a compensating term depending on the metric tensor (mass matrix) to nullify
this bias. In the absence of tractable algorithms to compute this compensating
potential and its gradient, its use has been impractical. This paper derives a
new algorithm for computing the compensating potential and its gradient for
tree-topology molecular systems. The algorithm is an extension of the spatial
operator based $O(\mathcal{N})$ algorithm for constrained dynamics; the
compensating potential is closely related to and computed from the articulated
body inertia quantities already available in that $O(\mathcal{N})$ algorithm.

## 1. Introduction

Rigid internal constraints are used in molecular models to eliminate high
frequency modes and enable larger integration time steps. Fixman pointed out
that ensemble averages from constrained vs. unconstrained models differ due to a
metric-tensor-dependent term in the partition function for constrained models.
Fixman proposed augmenting the standard potential by a compensating metric-tensor
potential in constrained MD to offset this bias. The prohibitive complexity of
computing the metric tensor potential and its gradient has been the major hurdle
to broad use.

This paper derives substantially simpler closed-form expressions for the
compensating potential and its gradient for tree-topology models, and derives
$O(\mathcal{N})$ computational algorithms. The compensating potential and its
gradient are readily computable from the articulated body inertia quantities
available from the spatial-operator $O(\mathcal{N})$ constrained-dynamics
algorithm.

## 2. Ensemble Averages in Constrained Dynamics Simulations

For an $n$-degree-of-freedom **unconstrained** Cartesian molecular model, the
partition function is given by Eq. (2.1). The system kinetic energy is
$\tfrac{1}{2}p^*\mathcal{M}^{-1}p$; $\mathcal{M}\in\mathbb{R}^{n\times n}$ is the
system **mass matrix** (metric tensor). For unconstrained Cartesian dynamics
$\mathcal{M}$ is constant and diagonal (configuration independent), so the
ensemble average of $f(q)$ takes the form of Eq. (2.2), where the momentum
variables have been integrated out.

In **constrained** dynamics models the molecular system is a collection of rigid
**clusters** coupled by articulable **hinges**. The mass matrix
$\mathcal{M}(\boldsymbol{\theta})\in\mathbb{R}^{\mathcal{N}\times\mathcal{N}}$ is
now a function of the internal configuration coordinates $\boldsymbol{\theta}$,
with $\mathcal{N}$ the number of constrained degrees of freedom. Kinetic energy
is given by Eq. (2.3), with conjugate momenta Eq. (2.4). The constrained
partition function Eq. (2.5) is obtained; a diagonalizing momentum transformation
Eq. (2.6) lets one integrate over momenta to yield Eq. (2.7). The resulting
constrained ensemble average Eq. (2.8) contains an extra determinant factor
$\det\{\mathcal{M}^{1/2}(\boldsymbol{\theta})\}$ absent from the unconstrained
case; this introduces a bias in statistical averages.

Fixman's remedy is to replace the potential $\mathcal{V}(\boldsymbol{\theta})$ by
the modified potential $\mathcal{V}'(\boldsymbol{\theta})$ of Eq. (2.9), adding
the **compensating mass matrix potential** (metric-tensor potential)
$\mathcal{V}_c(\boldsymbol{\theta})=\tfrac12\ln\det\{\mathcal{M}(\boldsymbol{\theta})\}$.
Substituting $\mathcal{V}'$ eliminates the metric tensor from the partition
function and ensemble average. Remaining differences between constrained and
unconstrained averages are then due only to coarser sampling of conformational
space.

Using $\mathcal{V}'$ means its gradient is used for the forces. The overall hinge
torque vector $T'=\nabla_\theta\mathcal{V}'=T+T_c$ (Eq. 2.10). The compensating
hinge torque $T_c$ arises from $\mathcal{V}_c$; its $k$th element is Eq. (2.11).

Consensus from prior simulation work: the compensating potential is less
important for rigid bond-stretch constraints but a significant factor for rigid
bond-angle constraints. For $n$-butane its inclusion bridged the gap in the
number of dihedral transitions between constrained and unconstrained MD. Despite
its accepted importance, it is rarely used in practice because of the lack of a
tractable method for computing $T_c$ for all but simple systems.

### Derivation (not implemented): expression for the compensating torque $T_c(k)$

For a scalar function $g(X)$ of a matrix $X\in\mathbb{R}^{m\times n}$, the
derivative w.r.t. a scalar $y$ is Eq. (2.12) (a trace). For
$g(X)=\ln\det\{X\}$, the standard matrix-calculus identity Eq. (2.13),
$\partial\ln\det\{X\}/\partial X=\{X^*\}^{-1}$, holds. Substituting (2.13) and
(2.12) into (2.11) gives the trace form Eq. (2.14),
$T_c(k)=\tfrac12\operatorname{Trace}\{\mathcal{M}^{-1}\mathcal{M}_{\theta(k)}\}$,
with the shorthand $\mathcal{M}_{\theta(k)}\equiv\partial\mathcal{M}/\partial\theta(k)$.

## 3. Spatial Operator Form of Mass Matrix

For notational simplicity the initial discussion is an $n$-cluster serial chain
with single-degree-of-freedom rotational hinges between clusters; the number of
degrees of freedom is $\mathcal{N}=n+5$ (the base cluster contributes 6 dof, and
each of the remaining serial hinges contributes 1). Extension to general
tree-topology systems and multi-dof hinges is discussed later.

The recursive Newton–Euler equations of motion (Eq. 3.1) have a base-to-tip
velocity/acceleration sweep and a tip-to-base force sweep. Introducing **spatial
operators** yields the concise operator form Eq. (3.2), so that the equations of
motion have the form Eq. (3.3), $T=\mathcal{M}(\boldsymbol{\theta})\ddot{\boldsymbol{\theta}}+\mathcal{C}(\boldsymbol{\theta},\dot{\boldsymbol{\theta}})$,
with the **Newton–Euler operator factorization** of the mass matrix Eq. (3.4a),
$\mathcal{M}=H\phi M\phi^*H^*$, and the Coriolis/gyroscopic/Cartesian force
vector Eq. (3.4b). $\mathcal{M}$ and $\mathcal{C}$ are nonlinear in
$\boldsymbol{\theta},\dot{\boldsymbol{\theta}}$.

Directly forming $\mathcal{M}$ and solving Eq. (3.3) is $O(\mathcal{N}^3)$; the
alternative recursive algorithm is $O(\mathcal{N})$.

### Innovations factorization of the mass matrix

**Lemma 3.1.** The Innovations Operator Factorization Eq. (3.5) gives closed-form
operator expressions for the square factorization (block $LDL^*$) and inversion
of the mass matrix. The factor $[I+H\phi K]$ is square, block lower triangular,
nonsingular, with identity blocks on the diagonal; $D$ is block diagonal. The
inverse factor is $[I+H\phi K]^{-1}=[I-H\psi K]$. Proof: see Ref. [2].

**Algorithm 3.1 (articulated body inertia recursion).** The articulated body
inertia quantities $P(\cdot),D(\cdot),G(\cdot),K(\cdot),\tau(\cdot),\bar\tau(\cdot),P^+(\cdot),\psi(\cdot,\cdot)$
are computed by the tip-to-base recursion Eq. (3.6). This is the classical
Riccati equation of Kalman filtering; $P(k)$ is the articulated body inertia of
the part of the system outboard of hinge $k$. The operator $P$ is block-diagonal
$6n\times6n$ with $P(k)\in\mathbb{R}^{6\times6}$. Auxiliary operator definitions
are Eq. (3.7). $D,G,\bar\tau$ are block diagonal; $K,\mathcal{E}_\psi$ have
nonzero blocks only on the first subdiagonal.

**Lemma 3.2.** The generalized accelerations $\ddot{\boldsymbol{\theta}}$ in terms
of hinge forces $T$ and Cartesian spatial forces $\hat f_c$ are given by
Eq. (3.8). Proof: see Ref. [2].

The recursive implementation of Eq. (3.8) is the $O(\mathcal{N})$ algorithm for
the accelerations: a tip-to-base residual-force sweep Eq. (3.9) followed by a
base-to-tip acceleration sweep Eq. (3.10). Summary of steps:

1. Base-to-tip recursion (Eq. 3.1) computes orientation, location, spatial
   velocities $V(k)$, and Coriolis/gyroscopic terms $a(k),b(k)$.
2. Tip-to-base recursion (Eq. 3.6) computes $P(k)$ etc.
3. Tip-to-base recursion (Eq. 3.9) computes residual forces $z(k)$ (can be fused
   with step 2 into one tip-to-base sweep).
4. Base-to-tip recursion (Eq. 3.10) computes $\ddot{\theta}(k)$.

Cost is linear in the number of clusters; the structure resembles Kalman
filtering and smoothing.

## 4. Compensating Mass Matrix Torque $T_c(i)$

Starting from Eq. (2.14), a simple-to-compute expression for $T_c(i)$ is
developed. An expression is needed for the derivative of the mass matrix with
respect to hinge coordinates.

### Spatial operator expression for $\mathcal{M}_\theta$

**Lemma 4.1.** The mass matrix sensitivity is the closed-form operator
expression Eq. (4.11),
$\mathcal{M}_{\theta_i}=H\phi[\mathbb{H}_\delta^i\phi M-M\phi^*\mathbb{H}_\delta^i]\phi^*H^*$.
Proof: see Ref. [26]. Here $\mathbb{H}_\delta^i$ is the $6n\times6n$ matrix that
is all zero except a single $6\times6$ block $\mathbb{H}(i)$ at the $i$th diagonal
location, where $i$ is the joint index for $\theta_i$. The nonzero block is
Eq. (4.12), built from the hinge rotational axis unit vector $h(i)$ via the
cross-product tensor $\tilde h(i)$. The tilde $\tilde v$ denotes the $3\times3$
cross-product (skew) matrix of a 3-vector $v=(x,y,z)^\top$.

That the formula is closed-form is important: the mass matrix derivatives can be
computed with spatially recursive algorithms similar to those computing the mass
matrix itself.

### Spatial operator expressions for $T_c(i)$

**Lemma 4.2.** $T_c(i)=\operatorname{Trace}\{P\Omega\mathbb{H}_\delta^i\}$ with
$\Omega\equiv\psi^*H^*D^{-1}H\psi\in\mathbb{R}^{6n\times6n}$ (Eq. 4.13).

#### Derivation (not implemented): proof of Lemma 4.2

Two spatial operator identities (derived in the appendix of Ref. [27]) are used:
Eq. (4.14a) $[I-H\psi K]H\phi=H\psi$ and Eq. (4.14b)
$\phi M\Omega=(\phi-\psi)+P\Omega$. Substituting Eq. (2.11), the inverse
factorization Eq. (3.5c) and the sensitivity Eq. (4.11), and repeatedly using
$\operatorname{Trace}\{AB\}=\operatorname{Trace}\{BA\}$, identity (4.14a), and
identity (4.14b), reduces the trace to
$\operatorname{Trace}\{P\Omega\mathbb{H}_\delta^i\}$. The step
$\operatorname{Trace}\{(\phi-\psi)\mathbb{H}_\delta^i\}=0$ holds because
$(\phi-\psi)$ is strictly lower triangular and $\mathbb{H}_\delta^i$ is block
diagonal.

It has been shown (Ref. [19]) that $\Omega$ decomposes as Eq. (4.15),
$\Omega=Y+\tilde\psi^*Y+Y\tilde\psi$ with $\tilde\psi\equiv\psi-I$, where the
block-diagonal $Y(k,k)\in\mathbb{R}^{6\times6}$ obeys the base-to-tip recursion
Eq. (4.16).

**Lemma 4.3.** $T_c(i)=\operatorname{Trace}\{P(i)Y(i)\mathbb{H}(i)\}$ (Eq. 4.17).

#### Derivation (not implemented): proof of Lemma 4.3

Using Eq. (4.15) in Eq. (4.13):
$T_c(i)=\operatorname{Trace}\{P(Y+\tilde\psi^*Y+Y\tilde\psi)\mathbb{H}_\delta^i\}
=\operatorname{Trace}\{PY\mathbb{H}_\delta^i\}=\operatorname{Trace}\{P(i)Y(i)\mathbb{H}(i)\}$,
because $\operatorname{Trace}\{PY\tilde\psi\mathbb{H}_\delta^i\}=0$ ($\tilde\psi$
strictly lower triangular).

**Lemma 4.4 (final simplification).** Partition the $6\times6$ matrix $P(i)Y(i)$
into $3\times3$ blocks $Q_{jk}$. Then Eq. (4.18),
$T_c(i)=-h^*(i)\,\mathcal{F}[Q_{11}+Q_{22}]$, where the map
$\mathcal{F}[\cdot]:\mathbb{R}^{3\times3}\to\mathbb{R}^3$ is defined by
$v=\mathcal{F}[A]$ iff $\tilde v=A-A^*$.

#### Derivation (not implemented): proof of Lemma 4.4

$P(i)Y(i)\mathbb{H}(i)$ has the block form Eq. (4.19), so
$\operatorname{Trace}\{P(i)Y(i)\mathbb{H}(i)\}=\operatorname{Trace}\{(Q_{11}+Q_{22})\tilde h(i)\}$.
Using the identity $\operatorname{Trace}\{A\tilde v\}=-v^*\mathcal{F}[A]$ gives the
result.

This expression for $T_c(i)$ is vastly simpler than the original mass-matrix
inverse/sensitivity form and reuses articulated body inertia quantities already
available. The compensating torque for the six base-cluster degrees of freedom
is zero.

## O(N) Constrained Dynamics Algorithm with Compensating Potential

1. Base-to-tip recursion (Eq. 3.1): compute $V(k),a(k),b(k)$ for all clusters.
2. Tip-to-base recursion (Eq. 3.6): compute all articulated body inertia
   quantities $P,D$, etc.
3. Base-to-tip recursion (Eq. 4.16): compute $Y(k)$ for all links; compute the
   compensating torque $T_c(k)$ (Eq. 4.18) and $T'(k)$ (Eq. 2.10) simultaneously.
4. Recursions (Eq. 3.9) and (Eq. 3.10): solve for $\ddot\theta(k)$ with $T(k)$
   replaced by $T'(k)$.

The only significant change from the constrained-dynamics algorithm of Ref. [2]
is the additional Step 3, which is also $O(\mathcal{N})$; overall cost remains
$O(\mathcal{N})$ and marginal since it reuses articulated body inertia
quantities. Implemented as part of the **NEIMO** (Newton–Euler Inverse Mass
Operator) software package.

### Extensions

Ref. [28] gives an alternative for computing $Y(k)$ using dual articulated body
inertias instead of Eq. (4.16), advantageous for parallel implementation
(concurrent with Step 2 rather than sequential).

Extension to tree-topology molecules is identical to Ref. [2]: the system has
multiple tips and a single designated base cluster; all tip-to-base and
base-to-tip recursions become tips-to-base and base-to-tips recursions, with
"scatter"/"gather" steps at branching hinges. Multi-dof hinges are modeled as a
sequence of single-dof hinges interconnected by pseudo-clusters of zero mass and
inertia.

### Spatial operator expression for $\mathcal{V}_c$

**Lemma 4.5.** $\mathcal{V}_c(\boldsymbol{\theta})=\tfrac12\sum_{i=1}^n\ln\det\{D(i)\}$
(Eq. 4.20).

#### Derivation (not implemented): proof of Lemma 4.5

From the innovations factorization Eq. (3.5a),
$\det\{\mathcal{M}\}=\det\{I+H\phi K\}^2\det\{D\}$. Since $[I+H\phi K]$ is lower
triangular with identity blocks on the diagonal, its determinant is 1, hence
$\det\{\mathcal{M}\}=\det\{D\}$. Taking the logarithm gives the result.

This allows easy computation of $\mathcal{V}_c$ using Algorithm 3.1 for the
required $D(k)$ quantities.

## 5. Conclusions

Internal rigid constraints with internal coordinates speed up MD but introduce
systematic biases into statistical averages; a compensating mass matrix potential
offsets these biases. The lack of a tractable method for computing the
compensating terms had rendered them impractical. This paper derives analytical
expressions and $O(\mathcal{N})$ algorithms for the compensating mass matrix
potential and its gradient for general tree-topology molecular systems, as an
extension of the $O(\mathcal{N})$ internal-coordinate MD algorithm. Extension to
closed-topology molecular systems is ongoing research.
