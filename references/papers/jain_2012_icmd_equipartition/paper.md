# Equipartition Principle for Internal Coordinate Molecular Dynamics

**Authors:** Abhinandan Jain, In-Hee Park, Nagarajan Vaidehi
**Venue:** Journal of Chemical Theory and Computation (JCTC), 2012
**DOI:** 10.1021/ct3002046

## Abstract

The principle of equipartition of (kinetic) energy for all-atom Cartesian
molecular dynamics states that each momentum phase space coordinate on the
average has kT/2 of kinetic energy in a canonical ensemble. This principle is
used in molecular dynamics simulations to initialize velocities, and to
calculate statistical properties such as entropy. Internal coordinate molecular
dynamics (ICMD) models differ from Cartesian models in that the overall kinetic
energy depends on the generalized coordinates and includes cross-terms. Due to
this coupled structure, no such equipartition principle holds for ICMD models.
In this paper, we introduce noncanonical **modal coordinates** to recover some
of the structural simplicity of Cartesian models and develop a new equipartition
principle for ICMD models. We derive low-order recursive computational
algorithms for transforming between the modal and physical coordinates. The
equipartition principle in modal coordinates provides a rigorous method for
initializing velocities in ICMD simulations, thus replacing the ad hoc methods
used until now. It also sets the basis for calculating conformational entropy
using internal coordinates.

## 1. Introduction

The equipartition theorem for canonical ensembles is a fundamental principle of
statistical mechanics. For Cartesian molecular models, each momentum coordinate
in the canonical phase space has kT/2 of thermal energy on average (k is the
Boltzmann constant, T the thermodynamic temperature). This is used to initialize
atomic velocities with a Boltzmann distribution and to connect the classical
variance to the quantum harmonic oscillator for absolute-entropy calculations.

In contrast, internal-coordinate molecular dynamics (ICMD) models use relative
(rather than absolute) generalized coordinates. Examples:

- **BAT coordinates** (bond length / angle / torsion): an alternative
  representation that retains the full system degrees of freedom (N = 3n).
- **TAMD** (torsion angle molecular dynamics): BAT models with additional
  holonomic constraints that freeze bond and angle coordinates, giving fewer
  coordinates (N < 3n). Attractive for tracking low-frequency biomolecular
  motion, but with qualitatively different dynamics from Cartesian models.

Unlike the Cartesian case, the ICMD Hamiltonian is **not separable**: the
kinetic energy depends not only on the momentum but also on the generalized
coordinates, giving dynamical cross-coupling among coordinates. This complicates
energy-conserving integration schemes analogous to Velocity Verlet.

The intuitive presumption that constrained ICMD models behave as the limiting
case of increasingly stiff Cartesian models is incorrect: statistical-mechanics
ensemble averages of conformation-dependent quantities from constrained ICMD
models differ systematically from those obtained with arbitrarily stiff
Cartesian models. This led Fixman to propose a mass-matrix-tensor-based
compensating potential to correct such discrepancies.

This work shows that applying the equipartition theorem to the ICMD model does
NOT yield an equipartition principle analogous to that for Cartesian models.
Instead the ensemble averages involve configuration-dependent, coupled
coordinates that are hard to interpret. We introduce a coordinate transformation
defining new ICMD **modal coordinates** that decouple the kinetic energy (like
Cartesian models) and for which an equipartition principle holds — even though
these coordinates are **not canonical** in the Hamiltonian sense.

The transformations between physical ICMD velocity/momentum coordinates and
modal coordinates involve the square root of the configuration-dependent ICMD
mass matrix and its inverse. Using spatial operators we develop analytical
closed-form factorizations of the mass matrix, leading to recursive, low-order
O(N) algorithms that avoid explicit evaluation of the transformation matrices.
The modal equipartition principle is applied to assigning initial velocities in
TAMD simulations of protein folding and refinement.

### 1.1. Generalized Equipartition Theorem

Let H(q, p) and L(q, q̇) denote the Hamiltonian and the Lagrangian for an
n-dimensional dynamical system, with q ∈ R^n the generalized coordinates and
p ∈ R^n the conjugate momenta. The momenta satisfy eq (1). The canonical
ensemble partition function Z(T) is defined in eq (2), and the ensemble average
of a function f(q, p) in eq (3), where h is Planck's constant and the α_i, γ_i
integration limits are set by the problem geometry.

With y_i, y_j a pair of phase space coordinates, the equipartition theorem
states the ensemble average of f(q, p) = y_i ∂H/∂y_j is given by eq (4). Here
the y_i, y_j can be any of the elements of q or p, but the theorem **requires
that q and p be canonical phase space coordinates**. Gathering the averages for
all i, j gives the matrix form eq (5). A derivation via integration by parts
appears in the Supporting Information (Tolman's derivation).

### 1.2. Cartesian Molecular Dynamics Model

The Cartesian MD Hamiltonian for n atoms is eq (6), with x ∈ R^{3n} the atom
positions, U(x) the potential energy, and p_i a linear momentum component for an
atom of mass m_i. Applying eq (4) with y_i = y_j = p_i gives eq (7): the average
kinetic energy in each conjugate momentum coordinate is kT/2 (the *principle of
equipartition of kinetic energy*). This follows from three special properties of
the Cartesian Hamiltonian:

- the quadratic (harmonic) form of the kinetic energy,
- the absence of cross-terms in the kinetic energy,
- the separable nature of the Hamiltonian (kinetic energy independent of the
  generalized coordinates).

## 2. ICMD Models

Cartesian models use the absolute atom locations as generalized coordinates;
ICMD models use atom–atom relative coordinates. BAT coordinates retain the
dimensionality of the Cartesian coordinate space (N = 3n). TAMD models are BAT
models with additional holonomic constraints freezing bond and angle DOF, giving
N < 3n.

For an ICMD model with N configuration DOF, θ ∈ R^N is the set of generalized
coordinates. The kinetic energy depends on a configuration-dependent mass matrix
M(θ) ∈ R^{N×N} and takes the form eq (8); the conjugate momenta are eq (9). The
notation x* denotes transpose. M(θ) is symmetric and positive definite, but
unlike the Cartesian case it is configuration dependent and dense, so the
kinetic energy includes cross-terms. The ICMD Hamiltonian eq (10) is **not
separable**. The partition function is eq (11).

From eq (10), ∂H/∂p_j = θ̇_j (eq 12), where e_j* is a unit vector with 1 in the
jth element. Using the equipartition theorem eq (4) gives eq (13). While the LHS
ensemble average has units of energy, it does **not** admit an interpretation as
equipartition of kinetic energy, because p_i depends on multiple velocity
coordinates, not just θ̇_i. Combining over all i, j gives the matrix form
eq (14). Unlike the Cartesian eq (7), there is no clear way to use this relation
to assign thermal energy across the DOF with a Boltzmann distribution.

### 2.1. Modal Coordinates

Since M(θ) is symmetric positive definite, there is an invertible
m(θ) ∈ R^{N×N} such that M = m m* (eq 15). Let l(θ) be the inverse of m(θ),
so M^{-1} = l* l (eq 16). Define the modal coordinates v ∈ R^N by eq (17):
v = m*(θ) θ̇ = l(θ) p. Using v, the Hamiltonian eq (10) becomes eq (18):
H = ½ v* v + U(θ). The kinetic energy is now **decoupled** (no cross-terms) and
no longer explicitly depends on θ — hence "modal" coordinates.

Unlike p, the v coordinates are **not** conjugate momenta of θ, so (θ, v) is
**not** a canonical phase space pair. Changing integration variables from p to v,
eq (17) gives the volume-element relation eq (19). Substituting into the
partition function gives eq (20) and the ensemble average eq (21). Defining the
corrected potential U'(θ) via eq (22) — where U_c(θ) = ½ ln det{M(θ)} is the
Fixman-like potential term — the ensemble average simplifies to eq (23).

For f(v) = v_i v_j the ensemble average is eq (24), which integrates to eq (25):
⟨v_i v_j⟩ = kT δ_{i=j}. Since v_i²/2 is the kinetic energy of modal coordinate
v_i, this means the modal components are uncorrelated and each carries kT/2 of
kinetic energy on average — the **equipartition principle for ICMD models**,
the analog of the Cartesian eq (7). Remarkably this holds even though (θ, v) are
noncanonical (no such principle holds for the canonical (θ, p)). The matrix form
is eq (26): ⟨v v*⟩ = kT I_N.

## 3. Physical to Modal Coordinate Transformations

During simulations one must transform between the v modal and θ̇ physical
velocities (eq 27): v = m*(θ) θ̇ and θ̇ = l*(θ) v. The configuration-dependent
Jacobian J(θ) ∈ R^{3n×N} relates Cartesian velocities ẋ to ICMD velocities θ̇
via eq (28): ẋ = J(θ) θ̇. This gives the ICMD mass matrix eq (29):
M(θ) = J*(θ) M J(θ), where M ∈ R^{3n×3n} is the constant diagonal Cartesian atom
mass matrix.

For **BAT models** N = 3n and J(θ) is square and invertible, so m(θ) is simply
eq (30): m(θ) = J*(θ) M^{1/2}, square and invertible.

For **constrained (TAMD) models** with holonomic constraints, N < 3n; eqs (28)
and (29) still hold but J(θ) is neither square nor invertible, so eq (30) yields
a nonsquare, noninvertible m(θ). One option is to evaluate M explicitly via
eq (29) then factorize numerically to obtain a square invertible m(θ).

### 3.1. Analytical Expressions for m(θ) and l(θ)

Spatial operators give analytical expressions leading to recursive, low-order
algorithms for transforming between v and θ̇ **without explicitly evaluating**
m(θ) and l(θ). These are valid with or without holonomic constraints (essential
when constraints are present) and assume a tree-topology ICMD model (BAT, TAMD);
closed-chain topologies are handled via constraint-embedding techniques (not
discussed).

Key analytical spatial-operator expressions for square factorization and
inversion of the mass matrix are eq (31). The first line is the Newton–Euler
factorization in terms of the H hinge-articulation operator, the φ rigid-body
propagation operator, and the M link spatial-inertia operator (nonsquare
factors). The second line is an alternative factorization with **only square
factors**: block-diagonal D and block-lower-triangular [I + H φ K], associated
with the articulated-body (AB) forward-dynamics algorithm. The third and fourth
lines give the analytical inverse of [I + H φ K] and hence of M. These hold
generally for tree-topology systems regardless of number of bodies, hinge types,
topology, or even nonrigid links.

From the second line of eq (31) one identifies a square m(θ) satisfying
M = m m* (eq 15), giving the analytical transformation matrices eq (32):
m(θ) = [I + H φ K] D^{1/2}, l(θ) = D^{-1/2} [I − H ψ K]. These lead to the
transformation expressions eq (33). Direct evaluation would cost O(N²), but the
special internal structure of the spatial operators allows the matrix/vector
products to be carried out recursively at **O(N) cost**. The recursive
base-to-tip algorithms are in Table 1. The intermediate quantity V(k) ∈ R^6 is
the combined angular and linear velocity of the kth coordinate frame. The table
is for serial topology; general tree topology is in ref 19.

## 4. Theoretical Implications for Velocity Initialization

At the start of an MD simulation, initial velocities must be assigned consistent
with the Boltzmann distribution of thermal energy; the simulation temperature
sets the overall kinetic energy.

### 4.1. Cartesian Models

The Cartesian equipartition principle eq (7) lets kinetic energy be assigned
independently to each velocity coordinate ẋ(i) with a Boltzmann distribution of
mean kT/2, giving eq (34): ⟨M^{1/2} ẋ ẋ* M^{1/2}⟩ = kT I_{3n}.

### 4.2. Unconstrained ICMD Models

Because equipartition does not hold for physical (θ, p) coordinates, θ̇ cannot be
directly initialized. A common option: assign ẋ Cartesian atom velocities via
the Cartesian principle, then transition to θ̇. For full-dimensional BAT models
this uses eq (35): θ̇ = J^{-1}(θ) ẋ. Correctness is checked by verifying the ICMD
equipartition principle eq (26) holds — see eq (36), which confirms that for BAT
coordinates initializing Cartesian velocities then transforming via eq (35) is
statistically correct.

### 4.3. Constrained ICMD Models

For constrained models N < 3n and J(θ) is nonsquare/noninvertible, so eq (35)
fails. Cartesian velocities must be mapped to the lower-dimensional ICMD velocity
space via a mapping P ∈ R^{N×3n} (eq 37): θ̇ = P ẋ. Evaluating ⟨v v*⟩ gives
eq (38): ⟨v v*⟩ = m* P P* m kT. For consistency with the ICMD equipartition
principle eq (26), P must satisfy m* P P* m = I. One valid choice is eq (39):
P = M^{-1} H φ M^{1/2}, verified because m* P P* m = m* M^{-1} m = I. Unlike ad
hoc projection techniques, this gives a thermodynamically rigorous
initialization. A practical disadvantage is that eq (39) requires the expensive
mass-matrix inverse.

**Preferred direct procedure:** exploit the modal equipartition principle
eq (26). Initialize each modal velocity v_i with mean kT/2 energy distributed by
the Boltzmann distribution, then convert to θ̇ physical velocities via the O(N)
recursive algorithm on the right of Table 1 (section 3.1). Unlike eq (39), this
avoids the explicit mass-matrix inverse.

## 5. Conclusions

We derived the equipartition principle for ICMD models by introducing a set of
noncanonical modal coordinates, applicable to both unconstrained and constrained
ICMD models. This enables rigorous velocity initialization and sets the stage for
calculating torsional entropy from ICMD trajectories. Using spatial operators we
obtained analytical transformations between modal and physical velocity
coordinates carried out by low-order O(N) recursive algorithms that avoid
computing the transformation matrices.

Modal coordinates were originally identified as diagonalizing coordinates
decoupling the ICMD (non-Hamiltonian) equations of motion. A remarkable property
is their orthogonality to the Coriolis-forces term in the equations of motion — a
property that does not hold for physical velocity coordinates — which simplifies
control and stability analysis for constrained ICMD-like systems.
