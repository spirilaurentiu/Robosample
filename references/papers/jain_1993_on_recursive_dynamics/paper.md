# A Fast Recursive Algorithm for Molecular Dynamics Simulation

**Authors:** A. Jain, N. Vaidehi, G. Rodriguez (Jet Propulsion Laboratory / Caltech).
J. Comput. Phys. 106(2), 258 (1993).

## Abstract

A recursive $O(\mathcal{N})$ algorithm for solving the dynamical equations of
motion for molecular systems using internal-variable models. Internal-variable
models already reduce MD cost by an order of magnitude vs. Cartesian models; the
$O(\mathcal{N})$ algorithm gives an additional large speedup. The method uses
spatial operator algebra: spatial operators derive the equations of motion and
an operator expression for the system mass matrix. An alternative square
factorization of the mass matrix (the *innovations* factorization) yields a
closed-form expression for its inverse, from which follows a recursive algorithm
for the generalized accelerations whose cost grows **linearly** with the number
of degrees of freedom (vs. cubic for conventional constrained dynamics). For a
polypeptide of 400 residues the $O(\mathcal{N})$ algorithm is ~450x faster than
the $O(\mathcal{N}^3)$ algorithm. A simplified handling of potential-function
gradients avoids computing gradients w.r.t. internal coordinates.

## 1. Introduction

Large computation times limit the time-scales of MD simulations. Hard
constraints on bonds eliminate lightly excited high-frequency degrees of freedom
(inter-atomic oscillations, rotations about double bonds) that force small
integration step-sizes.

With constraints, the Cartesian equations of motion become
differential-algebraic equations (DAEs) rather than ODEs. SHAKE handles
inter-atomic constraints iteratively. An alternative is to use *internal
coordinates* to directly incorporate the constraints, giving a smaller set of
generalized coordinates and equations of motion that remain ODEs.

The internal-variable model here consists of rigid subunits of atoms — called
**clusters** — whose relative motion is described by internal coordinates.
Attention is restricted to **tree-topology** (branched, no closed loops among
clusters) models. For tree topologies the generalized coordinates are of minimum
dimension and the EOM are ODEs, so integration is simpler.

The internal-variable EOM for a tree-topology system with $\mathcal{N}$ degrees
of freedom is $\mathcal{M}(\theta)\ddot{\theta} + \mathcal{C}(\theta,\dot\theta)
= T(\theta)$ (eq:1.1). $\mathcal{M}$ is $\mathcal{N}\times\mathcal{N}$, symmetric,
positive definite, but non-diagonal and depends nonlinearly on $\theta$.
Conventional dynamics algorithms compute $\mathcal{M}$, $\mathcal{C}$, $T$
explicitly and solve the linear system — cost is $O(\mathcal{N}^3)$.

The algorithm here avoids explicit computation of $\mathcal{M}$ and the explicit
linear solve, and provides a simplified method for potential-function forces. It
is based on spatial operator algebra developed for multibody dynamics of robots,
spacecraft, and vehicles; the operators are closely related to those used in
Kalman filtering and smoothing. Extensions to closed-topology models are
straightforward.

## 2. Dynamical Models for Molecular Systems

All MD requires: (i) computing forces for the current conformation; (ii) solving
the EOM and integrating to obtain trajectories.

Potentials split into bonding (quadratic/harmonic/Morse) and non-bonding (van der
Waals, electrostatic, dipole-dipole, dispersion). The potential energy splits
into an internal-coordinate part and a Cartesian part:
$PE = \mathcal{P}[\theta] + \mathcal{P}[x]$ (eq:PE).

**Free-atom Cartesian model.** With $n$ atoms and $3n$-vector of inter-atomic
forces $f$, Newton's law gives $\mathcal{M}_c \ddot{x} = f$ (eq:2.1), where
$\mathcal{M}_c$ is $3n\times 3n$ diagonal (atomic masses) and $x\in\Re^{3n}$.

**Constraints.** $m$ hard constraints are expressed as instantaneous linear
constraints on atomic velocities: $A\dot{x} = B$ (eq:2.2), with
$A\in\Re^{m\times 3n}$ and $B\in\Re^m$ both configuration-dependent. Solving for
$\ddot x$ subject to constraints makes the system a DAE (handled by SHAKE or
Lagrange multipliers).

**Reduced form.** Eliminating components via eq:2.2 gives a generalized velocity
of minimal dimension $\mathcal{N} \triangleq 3n - m$ and a reduced ODE
$\mathcal{M}_r \ddot{x}_r = f_r$ (eq:2.3). $\mathcal{M}_r$ is
$\mathcal{N}\times\mathcal{N}$ but no longer diagonal.

Internal-variable models give the same low-dimensional ODE form more naturally
while preserving structural properties of the mass matrix.

## 3. Internal Variable Models

A **cluster** is a rigid body of atoms (or a section with frozen covalent
structure, e.g. benzene rings, amino-acid rings). A single atom is a cluster with
zero extent and zero rotational inertia.

One cluster is the **base cluster** (chosen to minimize the maximum branch length
in number of clusters — roughly mid-centered). Adjoining clusters have a
parent/child designation: the one on the path to the base is the parent. Each
cluster has a unique parent but zero or more children.

A **hinge** characterizes permissible relative motion between adjoining clusters,
with 0 to 6 degrees of freedom. 0-dof = rigid coupling; 6-dof = no coupling
(independent molecules). A torsional bond = 1-dof rotational hinge; a variable
bond length = 1-dof translational (sliding) hinge.

Each hinge assigns coordinate frames $\mathcal{O}^+$ and $\mathcal{O}^-$ to the
two clusters. For an $m$-dof hinge, $\theta$ is the $m$-vector of generalized
coordinates. The generalized hinge velocity $\beta$ is usually $\dot\theta$; for
multi-dof hinges a "quasi-coordinate" velocity may be preferable (e.g. relative
angular velocity for a ball-and-socket). Full base-cluster mobility is a 6-dof
hinge to the inertial frame. There is always an invertible kinematic map
$\dot\theta = \Pi(\theta)\beta$; the paper assumes $\beta = \dot\theta$ throughout
with no loss of generality.

The 6-dimensional relative spatial velocity across a hinge is $H^*\dot\theta$,
where $H^*$ is the $(6\times m)$ hinge matrix. Example $H^*$ columns select the
torsional row, the stretching row, or both (see notation.md). When axes of a
multi-dof hinge depend on each other, the hinge is decomposed into component
hinges separated by **pseudoclusters** (zero mass, zero extent) so each component
has a constant $H$ matrix. This decomposition is assumed throughout.

## 4. Equations of Motion

Spatial operator algebra with a Newton-Euler approach (reveals recursion). First
a **serial chain** (branchless); tree extensions later. Clusters are numbered
1..n from tip to base; the $k$th hinge couples clusters $k{+}1$ and $k$ with
frames $\mathcal{O}_k^+$, $\mathcal{O}_k^-$. Parent of cluster $k$ is $k{+}1$,
child is $k{-}1$. The base cluster $n$ has a 6-dof hinge to the inertial frame,
which carries index $n{+}1$. Shorthand: $x(k) \equiv x(\mathcal{O}_k^-)$.

Recursive spatial velocity (eq:4.2), spatial acceleration (eq:4.3), Coriolis
acceleration $a(k)$ (eq:4.4). Since $H(k)$ is constant in frames
$\mathcal{O}_k^\pm$, $\omega(k{+}1)$ may replace $\omega(k)$ in the second term of
eq:4.4.

The generalized force for the $k$th hinge is $T(k) \triangleq
\nabla_{\theta(k)}\mathcal{P}[\theta]$ (eq:4.5), easy to compute. The complex
gradient $\nabla_{\theta(k)}\mathcal{P}[x]$ is **avoided**. Instead the algorithm
uses the Cartesian atomic-force gradient $\hat f_i(k) =
\nabla_{x_i(k)}\mathcal{P}[x]$ and combines the atomic forces of a cluster into a
single effective 6-dimensional spatial force $\hat f_c(k)$ (eq:4.6), where
$l[\mathcal{O}_k^-, x_i(k)]$ is the vector from the hinge frame to atom $i$.

Force recursion and hinge force (eq:4.7): $f(k)$ is the spatial interaction force
at $\mathcal{O}_k^-$ between clusters $k{+}1$ and $k$; $M(k)$ is the spatial
inertia about $\mathcal{O}_k^-$; $b(k)$ is the gyroscopic spatial force (eq:A.6).
Incorporating $T(k)$ and $\hat f_c(k)$ separately avoids computing
$\nabla_{\theta(k)}\mathcal{P}[x]$.

The full Newton-Euler recursion for the whole serial chain is a base-to-tip sweep
for $V,\alpha$ then a tip-to-base sweep for $f,T$ (eq:4.8).

**Spatial operators.** $M = \mathrm{diag}\{M(1)\cdots M(n)\}$ ($6n\times 6n$),
$H = \mathrm{diag}\{H(1)\cdots H(n)\}$ ($\mathcal{N}\times 6n$), and the shift
operator $\mathcal{E}_\phi$ (strictly lower block-bidiagonal, eq:4.9). Stacking
$V,\alpha,\dots$ into $6n$-vectors, eq:4.2 becomes $V = \mathcal{E}_\phi^* V +
H^*\dot\theta$ (eq:4.10). The inverse $[I - \mathcal{E}_\phi]^{-1} = \phi$ is the
lower-triangular operator (eq:4.11) with $\phi(i,j) = \phi(i,i{-}1)\cdots
\phi(j{+}1,j)$ for $i>j$ (semigroup property). Hence $V = \phi^* H^*\dot\theta$
(eq:4.12).

All component equations become operator expressions (eq:4.13), giving the EOM
$T = \mathcal{M}(\theta)\ddot\theta + \mathcal{C}(\theta,\dot\theta)$ (eq:4.14)
with the **Newton-Euler operator factorization** of the mass matrix
$\mathcal{M} = H\phi M\phi^* H^*$ (eq:4.15a) and
$\mathcal{C} = H\phi(M\phi^* a + b + \hat f_c)$ (eq:4.15b). This factorization is
equivalent to the Newton-Euler recursion. Solving eq:4.14 for $\ddot\theta$ by
forming $\mathcal{M},\mathcal{C}$ and a linear solve is $O(\mathcal{N}^3)$.

## 5. Recursive Solution of Equations of Motion

An alternative **innovations operator factorization** of $\mathcal{M}$ whose
factors are square and invertible, giving a closed-form mass-matrix inverse and
an $O(\mathcal{N})$ algorithm for $\ddot\theta$.

### 5.1 Operator Expression for the Mass Matrix Inverse

A recursive (discrete Riccati) algorithm defines $P(k), D(k), G(k), K(k{+}1,k),
\bar\tau(k), P^+(k), \psi(k{+}1,k)$ per cluster (eq:5.1). These form the operators
$D, G, K, \bar\tau$ (eq:5.2), the shift operator $\mathcal{E}_\psi =
\mathcal{E}_\phi\bar\tau$ (eq:5.3), and $\psi = [I-\mathcal{E}_\psi]^{-1}$
(eq:5.4) with the semigroup property $\psi(i,j) = \psi(i,i{-}1)\cdots
\psi(j{+}1,j)$ (eq:5.5). $\psi$ has the same structure as $\phi$, so operator
expressions map directly to recursions and $\psi$ need not be formed explicitly.

Innovations factorization (Lemma 5.1): $\mathcal{M} = [I + H\phi K] D [I + H\phi
K]^*$ (eq:5.6) — a closed-form block $LDL^*$ decomposition. Inverse of the factor
(Lemma 5.2): $[I + H\phi K]^{-1} = [I - H\psi K]$ (eq:5.7). Hence the mass-matrix
inverse (Lemma 5.3): $\mathcal{M}^{-1} = [I - H\psi K]^* D^{-1} [I - H\psi K]$
(eq:5.8).

### 5.2 Recursive Computational Algorithm

Generalized accelerations (Lemma 5.4, eq:5.9) decompose into the operator
sequence eq:5.10, whose recursive implementation is the $O(\mathcal{N})$
algorithm eq:5.11a (tip-to-base) and eq:5.11b (base-to-tip). Steps:

1. Base-to-tip recursion: orientation, location, spatial velocities $V(k)$, and
   Coriolis/gyroscopic terms $a(k), b(k)$ (eq:4.2, eq:4.4, eq:A.6).
2. Tip-to-base recursion (eq:5.1): compute $P(k)$.
3. Tip-to-base recursion (eq:5.11a): residual forces $z(k)$ — can be merged with
   step 2 into a single tip-to-base sweep.
4. Base-to-tip recursion (eq:5.11b): accelerations $\ddot\theta(k)$.

No explicit mass matrix and no linear solve are needed. Cost is linear in the
number of clusters. Because it solves the EOM exactly, it is numerically more
stable than iterative SHAKE. The structure mirrors Kalman filtering/smoothing.

### 5.3 Extensions to Branched Molecular Structures

For tree topology each cluster may have several children: a tips-to-base
recursion (summing results from all children at each cluster) followed by a
base-to-tips recursion (continuing separately along each outgoing child branch).
Closed-topology extensions reduce to solving a tree-topology subsystem with
exactly this algorithm.

### 5.4 Computational Costs

$O(\mathcal{N})$ cost is maximal for a serial chain with no point-mass clusters
and only single-dof hinges: ~$500\mathcal{N}$ flops. Point masses, multi-dof
hinges, or branches reduce the cost. The $O(\mathcal{N}^3)$ method costs roughly
$\mathcal{N}^3/3 + 19\mathcal{N}^2 + 350\mathcal{N}$ flops (serial chain,
single-dof hinges, no point masses). For a polypeptide with each residue a rigid
cluster and two bending dof between neighbors ($\mathcal{N}\approx 2\times$
residues), a 400-residue molecule gives ~450x speedup.

## 6. Conclusions

Internal-variable models allow larger integration time steps than Cartesian
models. The $O(\mathcal{N})$ algorithm from spatial operator algebra factorizes
and inverts the mass matrix in closed form, giving closed-form generalized
accelerations. It is recursive (tips-to-base and base-to-tips sweeps), needs
neither the mass matrix nor a linear solve, and solves the EOM exactly (better
numerical stability than SHAKE). Only gradients w.r.t. natural coordinates are
needed (simpler than internal-coordinate gradients). It extends naturally to
closed topologies and allows changing constraints during simulation (bond
making/breaking).

## Appendix A: Spatial Notation

Spatial velocity $V(\mathcal{O})$ combines angular and linear velocity
$\omega,v$; spatial force $f(\mathcal{O})$ combines moment and force $N,F$
(eq:A.1). Spatial acceleration $\alpha(\mathcal{O}) = \dot V(\mathcal{O})$.
$l(\mathcal{O}_x,\mathcal{O}_y)$ is the vector between frame origins. The
$6\times 6$ spatial transformation $\phi(\mathcal{O}_x,\mathcal{O}_y)$ (eq:A.2)
transforms forces and velocities between frames (eq:A.3). Spatial inertia
$M(\mathcal{O})$ (eq:A.4); rigid-body EOM about the frame $f = M\alpha + b$
(eq:A.5) with gyroscopic force $b(\mathcal{O})$ (eq:A.6).

### Derivation (not implemented) — Appendix B: Proofs of the Lemmas

The proofs parallel those in Refs. [13,21] for rigid multibody systems and are
algebraic manipulations of the spatial operators; they are not needed for
implementation. Key intermediate identities:

- $\bar\tau P\bar\tau^* = \bar\tau P$, giving the Riccati rewrite $M = P -
  \mathcal{E}_\psi P \mathcal{E}_\psi^* = P - \mathcal{E}_\phi P\mathcal{E}_\phi^*
  + KDK^*$ (eq:B.1); pre/post-multiplying by $\phi,\phi^*$ yields eq:5.6 (Lemma
  5.1).
- Matrix identity $[I + H\phi K]^{-1} = I - H\phi[I + KH\phi]^{-1}K$ (eq:B.2)
  together with $\psi^{-1} = \phi^{-1} + KH$ (eq:B.3) and $\psi^{-1}\phi = I +
  KH\phi$ gives Lemma 5.2 (eq:5.7).
- $[I - H\psi K]H\phi = H\psi$ (eq:B.5) and $\psi M\phi^* = \psi P + P\tilde\phi^*$
  (eq:B.7) reduce eq:B.4 to Lemma 5.4 (eq:5.9).
