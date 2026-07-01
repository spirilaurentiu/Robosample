# Fast Flexible Modeling of RNA Structure Using Internal Coordinates

Samuel Coulbourn Flores, Michael A. Sherman, Christopher M. Bruns, Peter Eastman, Russ Biagio Altman.
*IEEE/ACM Trans. Comput. Biol. Bioinform.* 2011; 8(5): 1247-1257. doi:10.1109/TCBB.2010.104.

## Abstract

Modeling the structure and dynamics of large macromolecules remains a critical
challenge. Molecular dynamics (MD) simulations are expensive because they model
every atom independently, and are difficult to combine with experimentally
derived knowledge. Assembly of molecules using fragments from libraries relies
on the database of known structures and thus may not work for novel motifs.
Coarse-grained modeling methods have yielded good results on large molecules but
can suffer from difficulties in creating more detailed full atomic realizations.
RNABuilder works in the internal coordinate space of dihedral angles and thus
has time requirements proportional to the number of moving parts rather than the
number of atoms. It provides accurate physics-based response to applied forces,
but also allows user-specified forces for incorporating experimental
information. A particular strength of RNABuilder is that all Leontis-Westhof
basepairs can be specified as primitives by the user to be satisfied during
model construction. RNABuilder predicts the structure of a 160-base RNA molecule
from its secondary structure plus experimental information, matching the known
structure to 10.2 Angstroms RMSD at low computational expense.

## 1. Introduction

RNA is central to gene regulation and expression, but understanding RNA
structure is limited because the molecules are difficult to crystallize and are
prone to misfolding and long-lived kinetic traps. Theoretical approaches are
challenged by counterions, molecule size, the delicate energetic balance between
alternative conformations, and long equilibration times during folding.

Several existing methods predict smaller RNA structures using knowledge-based
methods (FARNA, MC-Sym, DMD, NAST). NAST folds the ~160-nt P4/P6 domain of the
*Tetrahymena* group I intron to 16.3 Å RMSD (best computational prediction prior
to this work) but required 300 cpu-hours; the present method converges in a few
hours on a single processor with greater accuracy.

RNABuilder is an internal coordinate (IC) dynamics code. Because of its IC
framework, its computational requirements for the most part grow linearly with
system size (O(n)). It does not use fragment libraries. A key design feature is
to let the user make all modeling decisions (which forces are applied, which
regions are rigid vs flexible, how sterics are treated).

The software stack is: the open-source **Simbody** IC mechanics library at the
base, the **Molmodel** API layer providing a molecular interface to Simbody, and
**RNABuilder** built on Molmodel incorporating knowledge of RNA bases, geometry,
physics, force field, steric exclusion, and supported polymers.

### Internal coordinate scheme (method prose)

Dynamics can be computed in any coordinate set *q* in which Newton's second law
can be satisfied, with Cartesian coordinates obtained as needed via known
functions *x*(*q*), *y*(*q*), *z*(*q*). Under an IC scheme, atoms are partitioned
onto rigid bodies, and the bodies are interconnected via joints into an open tree
structure, with the coordinates *q* representing nonlinear joint coordinates
relating each body to its "parent" body within the tree. Algebraic constraints
*g*(*q*) = 0 are adjoined if there are interbody loops. Bond lengths and bond
angles are typically fixed, leaving the bodies free to rotate only about
rotatable bond axes; a major strength of the scheme is the freedom to choose
these mobilities explicitly. Recursive O(n) methods (linear in the number of
independent bodies) make IC dynamics practical; earlier Mazur et al. work
claimed O(n^3) scaling, which caused these methods to be abandoned in favor of
atomistic MD.

## 2. Methods

### 2.1 Simbody Internal Coordinate Mechanics Code

Simbody describes a system as a set of rigid bodies interconnected by
*mobilizers*. A rigid body has no inherent degrees of freedom and cannot move;
the only degrees of freedom present are those explicitly *granted* by a
mobilizer. A pin mobilizer (used to represent a torsional bond) defines a single
degree of freedom, the rotation angle around its axis. The more rigid the system,
the fewer degrees of freedom, and the more efficient the simulation. The *i*th
mobilizer defines a small number (1-6) of *generalized coordinates q_i* and
*generalized speeds u_i*, depending on the number of degrees of freedom it
introduces. These are aggregated into *q* = {*q_i*} and *u* = {*u_i*}, the
complete set of *nq* internal coordinates and *n* internal velocities. The
multibody tree system's equations of motion are:

<!-- eq:1 -->
$$ \dot{q}=N(q)\,u, \qquad M(q)\,\dot{u}=f(t,q,u), $$

where *M* is the *n* x *n* composite system mass matrix, *f* is *n* generalized
forces including applied and Coriolis forces, *N* is an *nq* x *n* block-diagonal
kinematic coupling matrix, and *t* is time. If the system is also subject to
constraints, its equations of motion become

<!-- eq:2 -->
$$ \dot{q} = N(q)\,u, \qquad M(q)\,\dot{u} = f(t,q,u) - G^{T}\lambda, \qquad g(t,q) = 0 $$

where *g*(*t,q*) is a set of *m* constraint equations, *G* = ∂*g*/∂*q*, and the
Lagrange multipliers λ represent the unknown constraint forces. Equation (1) is a
set of ODEs while (2) is a numerically challenging set of mixed differential and
algebraic equations (DAEs) of index 3.

Simbody uses variable-step, error-controlled integrators. This work uses the
**fourth-order Runge-Kutta-Merson** integrator: a variable step-size integrator
using five force evaluations per time step to produce a fourth-order accurate
trajectory and a third-order accurate error estimate. RNABuilder also offers the
**Velocity Verlet** integrator (conserves energy in optional fixed time-step
mode). Constraints are stabilized using the method of coordinate projection.
Monte Carlo simulation is also possible, but in IC mechanics a small torsion
change near the root of the biopolymer amplifies into a large displacement at the
other end, so MC has a high rejection rate.

Steric interactions are represented by elastic spheres placed at points on the
molecule that apply a repulsive force on contact (very short range, going to zero
as soon as spheres are no longer in contact), computed with Simbody's
Hunt-Crossley contact model.

Total forces *f*(*t,q,u*) combine modular component forces. RNABuilder uses three
main force subsystems: (1) contact forces for collision-detecting spheres to
prevent steric clashes, (2) a base-pairing module that uses a force-torque pair
to bring bases into the desired interaction geometry, and (3) a Tinker-style
force field with Amber99 parameters (by default only the bond-stretching term is
active). Temperature is maintained with either a velocity-rescaling thermostat or
a Nosé-Hoover thermostat.

### 2.2 Molmodel Extension to Simbody

Molmodel is a molecular modeling API layered on Simbody for modeling and
simulating molecules with customizable flexibility (all-atom Cartesian, internal
coordinate, fully rigid, and hybrid models). RNABuilder generally assumes bases
are fully rigid. There is a single ring-closing bond in the ribose ring that is
free to change length, angle, and dihedral to accommodate puckering motions; the
force field's bond-stretching term keeps this bond length within the normal
range. All remaining bonds have fixed lengths and angles but are free to rotate.

### 2.3 Enforcing Base-Base Interactions

RNABuilder can apply interactions between bases which at equilibrium reproduce
any of the base pair types classified in the Leontis-Stombaugh-Westhof catalog.
These consist of a force and torque that tend to align an *attachment frame* on
the first residue's base with a *body frame* on the second residue's base
(centered on the glycosidic nitrogen). The parameterization task is primarily
choosing the position and orientation of the *attachment frame* with respect to a
frame of reference fixed on the first residue's glycosidic nitrogen, to reproduce
a desired base-pairing geometry. The RNABuilder parameter file contains the
X,Y,Z distances and rotation angles of *attachment frames* A1 needed to generate
any of the Leontis and Westhof base pairs, as well as stacking and other
interactions.

The rotation that must be applied to align frame B2 with A1 is computed as:

<!-- eq:3 -->
$$ {}^{A1}R^{B2} = {}^{A1}R^{G} \cdot {}^{G}R^{B2} = \left( {}^{G}R^{A1} \right)^{-1} \cdot {}^{G}R^{B2}, $$

<!-- eq:4 -->
$$ {}^{G}R^{A1} = {}^{G}R^{B1} \cdot {}^{B1}R^{A1}. $$

Here we are given ${}^{G}R^{B1}$ and ${}^{G}R^{B2}$, the body frame orientations
with respect to ground, which are known as functions of the generalized
coordinates *q*, and ${}^{B1}R^{A1}$ the constant orientation of attachment frame
A1 in residue 1's body frame, known from the model of the particular base-pairing
interaction being enforced.

From Euler's rotation theorem, ${}^{A1}R^{B2}$ can be expressed as a rotation of
scalar angle θ about a unit axis. The following potential minimizes θ as well as
*r*, the translational distance between A1 and B2:

<!-- eq:5 -->
$$ U(r,\theta) = \left[ \frac{\theta^{2} \cdot \kappa}{2 \cdot k} + 1 \right] \cdot g(r, k, c) \cdot m, \qquad -\pi < \theta < \pi $$

where

<!-- eq:6 -->
$$ g(r,k,c) = \begin{cases} -\dfrac{k \cdot r^{2}}{2 \cdot c^{2}} + \dfrac{3 \cdot k}{2}, & 0 \leq r < c \\[2mm] \dfrac{k \cdot c}{r}, & r \geq c \end{cases} $$

<!-- eq:7 -->
$$ \vec{r} = r \cdot \hat{r} = \vec{x}_{B2} - \vec{x}_{A1}. $$

The constants κ and *k* are set separately for each interaction type in the
parameter file, where κ is typically positive and *k* is typically negative for
this potential type. The radial range *c* is set globally by the user, as is the
scaling factor *m*. The function *g* is harmonic at short range, decays with
inverse radius at long range (like electrostatic forces), and has a maximum
derivative (hence maximum force) at its inflection point at *c*. The force is

<!-- eq:8 -->
$$ \vec{F} = -\vec{\nabla}U = -\frac{\partial}{\partial\theta}U \cdot \hat{\theta} - \frac{\partial}{\partial r}U \cdot \hat{r} = -\frac{\theta \cdot \kappa}{k} \cdot g(r, k, c) \cdot m \cdot \hat{\theta} - \left[\frac{\theta^{2} \cdot \kappa}{2 \cdot k} + 1\right] \cdot g'(r, k, c) \cdot m \cdot \hat{r}, $$

where

<!-- eq:9 -->
$$ g'(r,k,c) = \begin{cases} -\dfrac{k \cdot r}{c^{2}}, & r < c, \\[2mm] -\dfrac{k \cdot c}{r^{2}}, & r \geq c. \end{cases} $$

The *translational* force is therefore

<!-- eq:10 -->
$$ \vec{f}_{A1} = \left[ \frac{\theta^{2} \cdot \kappa}{2 \cdot k} + 2 \right] \cdot g'(r, k, c) \cdot m \cdot \hat{r} = -\vec{f}_{B2}. $$

The angular dependence of this expression did not appear in the earlier work.

As a technical matter, the force is not applied to A1 directly but onto the
*body origin* of the first base, O1 (and forces on B2 are applied to O2). Moving
the point of application results in a torque that must be removed. The adjusted
torques are

<!-- eq:11 -->
$$ \vec{\tau}_{A1}^{*} = \frac{\theta \cdot \kappa \cdot m}{k} \cdot g(r, k, c) \cdot m \cdot \hat{\theta} + (\vec{x}_{A1} - \vec{x}_{O1}) \times \vec{f}_{A1}, $$

$$ \vec{\tau}_{B2}^{*} = \frac{\theta \cdot \kappa \cdot m}{k} \cdot g(r, k, c) \cdot m \cdot \hat{\theta} - (\vec{x}_{B2} - \vec{x}_{O2}) \times \vec{f}_{A1}. $$

The first term on the right-hand side in each expression is recognizable from
(8) and the second constitutes the described adjustment.

**Units.** RNABuilder inherits its system of units (picoseconds, nanometers,
kJ/mol, and Daltons) from Molmodel. The depth *k* and range *c* of the potential
are motivated by physicochemical experiments and statistical studies. Time,
energy, and temperature are not physically meaningful since the interactions are
imposed by the user and the dimensionality of the kinematics is reduced.

### 2.4 Steric Exclusion

Two collision-detecting contact-sphere schemes prevent atomic nuclei from
approaching too closely:

- **Reduced (*SelectedAtoms*) scheme:** contact spheres on the phosphorus, C4*,
  and glycosidic nitrogen atoms (up to four atoms, user-modifiable identity,
  radius, and stiffness per residue). Default radii determined by iteratively
  forming 10-base-pair helices and minimizing RMSD vs idealized helices from the
  make_na server. Designed for economy at some cost in accuracy.
- **Full (*AllHeavyAtomSterics*) scheme:** every atom except hydrogens gets a
  contact sphere, all with the same (user-adjustable) radius and stiffness.
  Default radius optimized similarly.

### 2.5 Protein Modeling

RNABuilder can also create protein chains (20 canonical amino acids, single-letter
code). The base-pairing force field and reduced sterics scheme are RNA-specific,
but all other features apply to proteins. The objective is a basic treatment of
the protein component of protein-RNA complexes such as the ribosome.

### 2.6 Modeling in Stages

Modeling parameters (temperature, size, reporting intervals, Amber99 force-field
term weights, etc.) can be grouped into "stages" to apply a multistep strategy in
a single run: change temperature across stages, apply forces/mobilities/sterics
in sequence, rigidify regions after convergence, apply time-ordered temperature
profiles for simulated annealing.

### 2.7 Computational Complexity and Memory Requirements

The computer time to integrate the equations of motion is O(*n*) for *n*
mobilities (degrees of freedom) in Simbody's formulation. Time to compute 1 ps of
dynamics for extended RNA chains of varying length (single core, Intel Nehalem)
was approximately O(*n*) for rigid chains and flexible chains with and without
AllHeavyAtomSterics; rigidifying bonds gives considerable savings. Memory: the
RNA Biopolymer required about 226 KBytes of real memory per residue; chains of at
least 13,000 residues were instantiated before running out of free memory
(32-bit mode, 4 GB RAM).

## 3. Application: Folding P4/P6

RNABuilder folded the P4/P6 domain of the *Tetrahymena* group I intron by
applying base-pairing contacts obtained from experiments and calculations that
did not make explicit reference to the crystallographic structure (UV-crosslinking,
dimethyl-sulfate and other protection assays, NMR, structural bioinformatics,
phylogenetics). The run (single core, Intel Nehalem) was repeated four times with
randomized initial velocities, obtaining 9.3, 10.1, 11.2, and 10.3 Å RMSD
averaged over the last nanosecond, some 6 Å lower than the best previously
published computational prediction. Each run required about 10.5 hours. The main
discrepancy was in the geometry and topology around P5c (whose reference-structure
position is affected by a crystal contact). No particular rigidification was
applied; rigidification beyond the default configuration diminished accuracy.

## 4. Conclusions

RNABuilder generates accurate structures of tRNA, P4/P6, and the entire
*Azoarcus* group I intron. Reducing the degrees of freedom and using
collision-detecting spheres rather than physical pairwise interactions achieved
convergence in a few hours, with RMSD 6 Å lower than previously published methods
and an order-of-magnitude lower computational expense. Empirical results confirm
O(*n*) computational complexity and that rigidification greatly reduces cost.
RNABuilder should be useful up to around 13,000 residues on a laptop in 32-bit
mode.

## 5. Availability

Binary distributions (Linux, Mac, Windows) and source code of RNABuilder are
available from the RNA-Toolbox project at https://simtk.org/home/rnatoolbox.
This work used RNABuilder revision 284, Molmodel revision 650, and Simbody
revision 1030.

## Figure notes (semantics retained)

- **Fig. 1:** An internal coordinate multibody system. Mobilizers (sticks) define
  relative motion between bodies (filled ellipses).
- **Fig. 3:** Default bond mobility for Molmodel/RNABuilder RNA residues. Backbone
  bonds, most ribose ring bonds, the C2'-O2' bond, and the glycosidic nitrogen
  bond are set to **Torsion** (fixed bond lengths and angles). The O4'-C1' bond is
  set to **Free** (no restriction). Base bonds and bonds of all single-coordinated
  atoms (hydrogens and some phosphate oxygens) are set to **Rigid** (no freedom).
- **Fig. 4:** The base-pairing interaction includes a force and a torque that act
  to align the *attachment frame* A1 of residue 1 with the *body frame* B2 of
  residue 2.
- **Fig. 5:** The base-pairing potential. At short range (*r* < *c*) the potential
  is quadratic in *r*; at long range it goes like the inverse, approaching zero at
  infinite *r* for all θ. Inflection point at *r* = *c* and θ = 0, where *U* = *k*.
  The dependence on θ is quadratic everywhere.
