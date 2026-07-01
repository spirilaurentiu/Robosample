# Derivation and Testing of Explicit Equations of Motion for Polymers Described by Internal Coordinates

Mazur, Dorofeev, Abagyan — J. Comput. Phys. 92 (1991).

## Abstract

General Lagrange's equations of motion for a system of polymeric molecules are
obtained in explicit form. They can be used for simulating molecular dynamics of
large molecules. The molecular conformations are described by internal
coordinates, i.e., bond lengths, valence angles, and torsion angles. The
equations derived permit any internal degrees of freedom to be frozen. The method
is applied to an oligopeptide in an α-helical conformation. Three models of the
molecule with different degrees of fixation are compared. It is shown that the
method permits one to increase significantly the time step in molecular dynamics
calculations.

## 1. Introduction

Standard MD for polymers uses Cartesian coordinates of atoms and Newton's
equations. This is effectively "atomic dynamics" rather than molecular dynamics.
Cartesian coordinates introduce limitations: many of the 3N degrees of freedom of
an N-atomic molecule are weakly excited at room temperature, yet they are treated
as classical oscillators forced to oscillate at very high frequencies. This limits
the integration time step and restricts the accessible time scale to the
subnanosecond range, while many interesting events (intramolecular rotations
around single bonds) take place over longer intervals.

One remedy is to freeze the covalent structure (bond lengths, valence angles,
aromatic rings). In the atomic representation this is done by imposing constraints
on Cartesian coordinates (e.g. fixing atom-atom distances, SHAKE-type methods),
but the computational cost grows rapidly as constraints are added; in practice
only bond-length fixation is efficient for large molecules.

A preferable approach freezes all fast internal vibrations by introducing rigid
bonds and angles directly into the molecular description, using generalized
coordinates and the Lagrange-Hamilton formalism. Historically this was applied
only to small molecules (e.g. n-butane) because deriving the Lagrange equations
for a general polymer was considered prohibitively complicated. This paper
presents a complete derivation of universal explicit equations of motion in
internal coordinates that allow any set of internal coordinates (bond lengths,
bond angles, dihedral angles) to be frozen, together with tests on a small protein
fragment.

## 2. The equations of motion for a system of branched polymer molecules

### 2.1. Formal description and original equations

Any number of molecules with any set of frozen internal coordinates (bond lengths
inclusive) can be described as a unified tree (a **BKS-tree**). The tree is
composed of atomic groups called **rigid bodies** whose internal structure is kept
fixed. The construction disconnects cycles formed by rigid bodies, applies a
special numeration of rigid bodies, and introduces virtual atoms and virtual bonds
to connect molecules within the system and impose a tree topology on it.

**Fig. 1 (semantic note):** A molecular system is represented as a BKS-tree. Rigid
bodies are groups of atoms depending on the same set of internal variables. The
tree is constructed from its origin, coinciding with the origin of the global
coordinate frame. Its structure is determined by internal coordinates (bond
lengths, planar and dihedral angles) that can be arbitrarily frozen. Rigid bodies
are connected at hinges/nodes to which the internal coordinates are attributed.

The non-fixed internal coordinates of the BKS-tree form the set of generalized
coordinates {θ_k} of the system. Each variable is attributed to a particular node
where a unit vector determining the infinitesimal displacements due to that
variable is defined.

**Ordering rule (critical).** There exists a natural order of variables owing to
the relations between variables and their unit vectors: if one variable influences
another's unit vector (only one direction of influence can exist), the influencing
variable always has the smaller index. The order of variables attributed to the
**same node** is (Fig. 2): first the torsion and phase angles (φ and Φ), then the
planar angle (ω), and finally the bond length (b). The phase angle Φ is a dihedral
angle added to the torsion angle φ to construct various branches coming from the
same node.

**Per-atom / per-variable sets.** The position of a particular atom α (radius
vector **r**_α) is determined by the chain V_α of generalized coordinates θ_i^α,
ordered as above at each node, with nodes ordered according to the rigid bodies
they define. The lower index i runs from 1 to n_α. For each variable θ_k, define
the set d_k of atoms whose positions depend on θ_k.

The derivation starts from the Lagrangian equations of motion (see eq:1) with
Lagrangian L = T(θ, θ̇) − U(θ), the difference of kinetic and potential energy.
Substituting L into eq:1 gives eq:2, whose right-hand side (the conformational
energy gradient) can be computed rapidly by previously described algorithms. The
aim of the paper is the explicit left-hand side (the mass-matrix / inertial
terms).

The infinitesimal displacement of atom α is given by eq:3. There, **e**_i^α are
the unit vectors of variable θ_i^α, **r**_i^θ is the radius vector of its node, and
S_i is an indicator: S_i = 1 for angle variables, S_i = 0 for bond-length
variables. The first term describes rotation of atom α about the axis **e**_i; the
second describes translation along the unit vector of a variable bond.

### 2.2. Computation of (d/dt)[ṙ_α (∂ṙ_α/∂θ̇_k)]

### Derivation (not implemented)

Differentiating eq:3 gives the velocity eq:4 and the partial derivative eq:5. Time
derivatives of the unit vectors and node positions are needed:

$$d\mathbf{e}_{i} = \sum_{m=1}^{i-1} S_{m}\, \mathbf{e}_{m} \times \mathbf{e}_{i} \cdot d\theta_{m}$$
(eq. 6)

$$d\mathbf{r}_{i}^{\theta} = \sum_{m=1}^{i-1} \left[ S_{m} \mathbf{e}_{m} \times (\mathbf{r}_{i}^{\theta} - \mathbf{r}_{m}^{\theta}) \cdot d\theta_{m} + (1 - S_{m}) \mathbf{e}_{m} \cdot d\theta_{m} \right]$$
(eq. 7)

$$d\mathbf{r}_{\alpha/i} = \sum_{m=1}^{i-1} S_m \mathbf{e}_m \times \mathbf{r}_{\alpha/i} \cdot d\theta_m + \sum_{m=i}^{n_\alpha} \left[ S_m \mathbf{e}_m \times \mathbf{r}_{\alpha/m} \cdot d\theta_m + (1 - S_m) \mathbf{e}_m \cdot d\theta_m \right]$$
(eq. 8)

with corresponding time derivatives (eqs. 9-11). Using these, the differentiated
inertial term (d/dt)[ṙ_α (∂ṙ_α/∂θ̇_k)] is expressed as 4 single sums and 16 double
sums (eq. 12; intervening computations omitted in the original). Equation (12) is a
long intermediate result superseded by the simplified final form eq:16 and is not
reproduced term-by-term here.

### 2.3. Computation of ½(∂ṙ_α²/∂θ_k)

### Derivation (not implemented)

The derivatives of **e**_m and **r**_{α/m} with respect to θ_k (eqs. 13, 14) follow
from eqs. 6 and 8. Squaring eq:4 and differentiating using eqs. 13-14 yields the
second (potential-of-velocity) term ½(∂ṙ_α²/∂θ_k) as eq. 15, another long
intermediate expression superseded by eq:16.

### 2.4. The ultimate form of the equations

Substituting eqs. 12 and 15 into eq:2 and simplifying gives the explicit equations
of motion eq:16. The only algebraic trick used is the change of summation order in
double sums (eq:17). The brace symbol V_α indicates that within the braces the
variables are numbered in succession according to their order in the chain V_α; the
indices i, k, m within the braces are functions of α (local numbering).

When **global numbering** is used (each variable has a unique index and i, j, k run
from 1 to the total number of variables n_var), the equations take the compact
linear-system form eq:18, with coefficient arrays a_ki (mass matrix), b_ki, and
c_kim assembled by summation over the tree topology. Using a procedure that
exploits the tree topology, assembling these coefficients takes only a minor part
of the total computation time.

## 3. A test example: molecular dynamics of an α-helix

The method was applied to an oligopeptide (Ala)_n in an α-helical starting
conformation, one molecule in vacuum (solvent neglected). Several probe
integrations were run in forward and backward time directions for each model, and
conservation of energy, momentum, and angular momentum was checked. For a
completely unfixed molecule, trajectories from eqs. 16 and 18 were compared with a
traditional Cartesian-coordinate MD simulation. All tests confirmed the
correctness of the derived equations.

Trajectories of 110 time steps were obtained by integrating eqs. 16 and 18. The
molecule was slowly heated to ~300 K and equilibrated during 4 ps (sufficient for
the mean fluctuation of kinetic energy ⟨δK⟩ to stabilize). Temperature was computed
from eq:temperature.

Since eq:18 is not resolved with respect to θ̈, the corresponding linear algebraic
system is solved at each time step by the Cholesky ("Kholezki") algorithm. The
equations are integrated by Beeman's method with fourth-order prediction formulas
for generalized velocities. Empirical potentials were compiled from ECEPP and
CHARMM.

**Three models (see checks.md for numbers):**
- Model 1: 93 atoms, all hydrogens explicit; all bond lengths fixed; phase angles
  fixed (no proper empirical potentials for them). 133 degrees of freedom.
- Model 2: as model 1 but all valence angles also fixed. 42 variables.
- Model 3: united atoms for alanine side chains (excluding methyl torsions);
  C-terminus hydroxyl torsion fixed. 32 degrees of freedom.

Accuracy is measured by conservation of total energy via the relative RMS
deviation eq:deltaE (used instead of ⟨E⟩ because ⟨E⟩ varies strongly between
models; only the varying parts of the potential energy are computed, since
intra-rigid-body interactions are always discarded).

Freezing the fastest degrees of freedom in each successive model both improves
energy conservation at a given step and makes the computational task cheaper
(unlike Newtonian equations). One integration step for the three models took 6.0,
1.1, and 3.5 s respectively (EC1051 computer): for model 1 solving the 133-equation
linear system dominates; for the smaller models the energy-gradient calculation
dominates. The main open problem is an effective method for solving the linear
system eq:18.
