# Internal Coordinate Molecular Dynamics: A Foundation for Multiscale Dynamics

Nagarajan Vaidehi, Abhinandan Jain. Feature Article, J. Phys. Chem. B (2015).

## Abstract

Internal coordinates such as bond lengths, bond angles, and torsion angles (BAT)
are natural coordinates for describing a bonded molecular system. The molecular
dynamics (MD) methods widely used for proteins, DNA, and polymers are based on
Cartesian coordinates owing to the mathematical simplicity of the equations of
motion, but constraints are often needed with Cartesian MD to enhance
conformational sampling, which makes the equations of motion
differential-algebraic and adversely impacts complexity and robustness. In BAT
coordinates constraints are easily placed by removing degrees of freedom.
Internal coordinate MD (ICMD) is an attractive alternative for developing
multiscale MD. Torsional MD is a special adaptation of ICMD where all bond
lengths and bond angles are kept rigid. Advantages: longer time step from
freezing high-frequency degrees of freedom, and conformational search in the
important low-frequency torsional degrees of freedom. This review summarizes the
authors' mathematical advances (spatial operator algebra / GNEIMO) toward a
robust long-time-scale ICMD toolkit, and applications to protein conformational
change and structure refinement.

## 1. Introduction

MD simulations are used for (a) dynamics of protein structures, (b) protein
structure prediction, and (c) thermodynamic properties (free energies, enthalpy,
entropy) of conformational states. All-atom Cartesian MD is a classical-mechanics
toolkit using Cartesian coordinates as degrees of freedom; its attraction is
mathematical simplicity. Constraints and/or bias potentials are used in Cartesian
MD to increase step size, enhance sampling, and simulate large conformational
changes, but adding constraints makes the equations of motion differential
algebraic, requiring DAE solvers that hurt robustness and complexity.

Bond length, bond angle, and torsion angle (BAT) relative coordinates are more
natural than Cartesian absolute coordinates for bonded proteins. MD in BAT
coordinates is called internal coordinate MD (ICMD). Nonessential high-frequency
degrees of freedom (e.g. bond lengths) can be constrained by simply excluding
them; the resulting ICMD models have fewer degrees of freedom and retain the
simpler ordinary-differential structure instead of differential-algebraic.

Other advantages of ICMD:

1. Low-frequency torsional coordinates allow larger integration time steps.
2. Conformational search is more effective in the low-frequency torsional
   degrees of freedom and leads to significant conformational changes.
3. Enhanced sampling methods are more effective in torsional space.
4. The six translation/orientation degrees of freedom are explicit coordinates
   (useful for conformational entropy via quasiharmonic analysis).
5. They provide a large range of options for selecting/controlling model
   granularity.
6. Fixman potential corrections for constraint-induced biases in the partition
   function are easier to apply.

### Challenges and solutions to the bottlenecks

Torsional MD constrains all bond lengths and bond angles. Two long-standing
bottlenecks: (a) mathematical/computational complexity of the equations of
motion; (b) increased rigidity from freezing bond lengths and angles, affecting
transition barriers and probability density functions.

When all bond lengths and bond angles are constrained, naive solution of the
dihedral-space equations of motion scales as the **cubic** power of the number
of torsion degrees of freedom. The authors developed the generalized
Newton-Euler inverse mass operator (**GNEIMO**) method, based on spatial operator
algebra originally developed for spacecraft and robot dynamics. Key insight:
analytical factorization and inversion of the mass matrix for tree-topology
systems, reducing cost to **linear** instead of cubic in the number of degrees
of freedom. GNEIMO has been used in CYANA (NMR refinement) and NIH-XPLOR (X-ray
refinement).

The second bottleneck (rigidity) leads to fewer dihedral transitions and
systematic errors in probability density functions. Rigid bond-length and
bond-angle constraints alter the potential and free energy surface versus
unconstrained MD. Fixman (1970s) proposed a compensating potential that
rigorously corrects treating stiff, uncoupled bond angles as rigid, generating a
partition function whose probability density functions approach all-atom
Cartesian simulations. The authors developed a spatial-operator-algebra,
general-purpose, low-cost method to compute the Fixman potential for all linear
and branched molecules, solving the long-standing tractability problem.

## 2. Methods

ICMD is MD in BAT (internal) coordinates. GNEIMO is a generalized ICMD method
based on spatial operator algebra for multibody dynamics of macromolecules.
Constraints on high-frequency bond lengths can be placed to run ICMD with bond
angles and torsions as degrees of freedom. If both bond lengths and bond angles
are kept rigid, the resulting **torsional MD** method is supported by GNEIMO. The
macromolecule of tree topology is modeled as a collection of rigid bodies
("clusters", varied sizes: single atom, methyl group, helix, or entire domain)
connected by flexible hinges. Hinges have one to six degrees of freedom; one dof
= just the torsion angle; six dof allows bond stretch, bond-angle bending, and
torsion about the connecting bond.

**Note (Fig 1):** Standard clustering scheme for torsional MD in GNEIMO,
illustrated for tripeptide Ala-Tyr-Ala. Rigid clusters are single-colored and
move as one unit; rods/arrows are hinges connecting two clusters.

### Equations of motion (rigid bond lengths and angles)

When bond lengths and bond angles are treated as rigid, the equations of motion
in ICMD become coupled:

<!-- eq:1 -->
$$ \mathcal{M}(\theta)\ddot{\theta} + C(\theta, \dot{\theta}) = \mathcal{T}(\theta) $$

where $\theta$ is the vector of generalized coordinates (e.g. torsional angles),
$\mathcal{T}$ is the vector of generalized forces (e.g. torques),
$\mathcal{M}(\theta)$ is the mass matrix (moment of inertia tensor), and
$C(\theta,\dot{\theta})$ includes the velocity-dependent Coriolis forces. The
dynamics is obtained by solving for $\ddot{\theta}$ and integrating to obtain new
velocities and coordinates. Conventional algorithms invert the dense mass matrix
at cubic cost in the number of degrees of freedom.

GNEIMO uses spatial operator algebra to derive an analytical expression for the
inverse of the mass matrix, giving the following expression for $\ddot{\theta}$:

<!-- eq:2 -->
$$ \ddot{\theta} = [I - \mathcal{H}\psi\mathcal{K}]\, \mathcal{D}^{-1} [\mathcal{T} - \mathcal{H}\psi(\mathcal{K}\mathcal{T} + \mathcal{P}a + b)] - \mathcal{K}^{*}\psi^{*}a $$

The $\mathcal{H}$, $\psi$, $\mathcal{K}$, etc., terms are the mass-matrix-related
spatial-operator factorizations (see Jain 2010, Jain-Vaidehi-Rodriguez 1993).
This right-hand side is evaluated by recursive algorithms whose cost scales
**linearly** with the number of degrees of freedom, avoiding the dense mass
matrix inversion of eq 1.

### Nosé-Hoover constant-temperature ICMD

The method was extended to the canonical (N, V, T) ensemble with the Nosé-Hoover
thermostat. The torsional MD equations of motion for the Nosé-Hoover ICMD method:

<!-- eq:3 -->
$$ \mathcal{M}(\theta)\ddot{\theta} + C(\theta, \dot{\theta}) + \mathcal{F}(\theta, \dot{\theta}) = \mathcal{T}(\theta) $$

<!-- eq:4 -->
$$ \dot{\eta} = \frac{1}{\tau^{2}} \left[ \frac{T}{T_{B}} - 1 \right] $$

Here $\mathcal{F}$ is the additional frictional force term due to the canonical
ensemble, dependent on $\eta$, the dynamic variable representing the thermostat;
$\tau$ is the mass parameter of the thermostat; $T$ is the instantaneous
temperature; $T_{B}$ is the thermostat (target) temperature. Because these
equations are similar in form to eq 1, all the spatial-operator equations and
factorizations hold. The thermostat mass parameter $\tau$ was optimized to be
**10 times the time step size** for torsional MD. Accuracy/stability were
measured by conservation of the total Hamiltonian and temperature fluctuations.

### Fixman compensating potential

To correct systematic biases in the probability density function caused by
treating stiff bond angles as rigid, Fixman proposed a compensating potential:

<!-- eq:5 -->
$$ \mathcal{U}_{f}(\theta) \triangleq \frac{1}{2} k T \ln \frac{\det\{\mathcal{M}(\theta)\}}{\det\{\mathcal{M}_{B}(\theta, q_{0})\}} $$

where $\mathcal{M}_{B}$ is the mass matrix in the full BAT coordinates, $q_0$ the
coordinates for the frozen degrees of freedom, $k$ the Boltzmann constant, and
$T$ the temperature. The Fixman potential removes such biases but was
computationally intractable for generalized branched molecules (used only for
small models such as C4/C5). The authors derived a spatial-operator-algebra
algorithm (**GNEIMO-Fixman**) that computes the Fixman potential with only ~24%
additional computation time, applicable to general linear and branched systems.
GNEIMO-Fixman also computes partial derivatives of the Fixman potential, defining
additional forces (the **Fixman torque**) applied within constrained MD.

**Note (Fig 2):** Joint probability distribution of the two backbone torsion
angles in alanine dipeptide: (a) Cartesian, (b) ICMD, (c) ICMD + Fixman. Torsional
MD introduces biases in the joint PDF that are removed with the Fixman potential
and torques. The magnitude of the Fixman potential is much smaller than the
all-atom force-field potential energy. It corrects errors from stiff, uncoupled
bond angles, but not soft bond angles coupled to torsions or nonbond interactions.

### Handling coupled (soft) bond angles

Two approaches to eliminate bias from treating torsion-coupled bond angles as
rigid: (1) open up some bond angles as movable degrees of freedom ("hybrid
internal coordinate molecular dynamics"), reducing rigidity and PDF/transition
error with little impact on time step; (2) correction torsional potentials
(e.g. ECEPP force fields, or ICMFF) that refit torsional force constants from the
rigid model to reproduce the flexible-model torsion energy barriers. ICMFF is not
rigorous and is system specific. Thermodynamic accuracy in the PDF requires
treating strongly-coupled bond angles as flexible; this is not required for
structure prediction (where the goal is enrichment of the native ensemble).

### Modal-coordinate equipartition

Constrained ICMD models are not the limiting case of stiff Cartesian models.
Applying the conventional Cartesian equipartition theorem to the ICMD model does
not yield an analogous equipartition principle; ensemble averages involve
configuration-dependent coupled coordinates that are hard to interpret. The
authors introduced a coordinate transformation to **modal coordinates** that
decouples the kinetic energy into Cartesian-like form. Using "modal velocity"
coordinates, an equipartition principle for ICMD analogous to the Cartesian one
is derived (holds even though modal coordinates are not canonical). This provides
thermodynamically correct velocity initialization for ICMD simulations.

### Software and scalability (GneimoSim)

GneimoSim is a modular software package for ICMD. It includes the ICMD
equipartition principle and Fixman-potential methods, and interfaces to LAMMPS,
OpenMM, and Rosetta force fields via a Python interface to underlying C++
classes.

**Note (Fig 3):** GNEIMO equations-of-motion solver cost scales linearly with
number of clusters. For AMBER + generalized Born on GPU (OpenMM), the GNEIMO EOM
solver cost is comparable to the force-field cost; for larger systems the GNEIMO
solver is linear while force calculation grows superlinearly. Extra ICMD cost is
partly compensated by larger time steps than Cartesian.

Implemented GneimoSim capabilities:

1. Advanced integrators: Runge-Kutta, Lobatto, adaptive CVODE, and Verlet, for
   long-time-scale (microsecond) torsional MD.
2. Generalized Born solvation (GBSA) for implicit solvation.
3. Multiple molecules of any type, including explicit solvent.
4. Temperature-based replica-exchange (REMD), with random or
   Metropolis-probabilistic temperature switching.
5. Periodic boundary conditions for explicit water.
6. Standard Cartesian simulations for comparison.
7. User-defined harmonic distance restraints between atom pairs.
8. Langevin dynamics and accelerated MD (aMD).

## 3. Results and Discussion

### Multiscale simulation with ICMD (dynamic clustering)

BAT coordinates readily allow freezing and thawing any internal coordinate degree
of freedom (e.g. freezing a whole helix). This freeze/thaw can be done at the
start or on the fly ("dynamic clustering"). Default clustering treats all
backbone/main-chain torsions as degrees of freedom. Adaptive clustering that
changes the cluster model during simulation is more suitable for large
conformational changes (e.g. Poursina et al. for RNA; Wagner et al. freeze
secondary structure at high-T replicas, thaw at low-T). Because freeze/thaw can
alter the pathway, it is advisable for conformational search only. Other
freeze/thaw criteria: an upper velocity threshold (hinge treated flexible if
velocity stays near threshold), or monitoring stress-force accumulation at a
rigid hinge (not yet implemented).

### Freeze-and-thaw folding

Four proteins were folded from extended structures with predicted secondary
structure elements, using GNEIMO torsional MD + REMD with the adaptive-time-step
CVODE (Adams-Moulton) integrator. For the B domain of staphylococcal protein A
(1BDD): 12 temperature replicas (300-1050 K). Folding begins from extended
structure with predicted helices; a few interhelical contacts collapse the
structure to 12-16 Å RMSD; at 8-11 Å slightly incorrect topologies are explored;
finally folds to correct topology below 7 Å. Most conformations fall between
5 and 7 Å.

### Domain motion (calmodulin)

GNEIMO-REMD (no bias potential) was applied to the whole calmodulin protein.
Calmodulin (EF-hand family) has N-terminal and C-terminal domains connected by a
long helical stretch. Upon calcium removal: (1) the central connecting helix
collapses, (2) the N-terminal domain moves relative to the static C-terminal
domain. Cartesian MD in explicit solvent showed the helix collapse but not the
full ensemble; GNEIMO-REMD captured both the collapse and the N-terminal domain
flexibility. GNEIMO conformations cover most of the NMR ensemble (PDB 1DMO). About
half of the average hydrogen-bond distances between residues fall within one
standard deviation of the NMR structures.

For BPTI (NVT torsional MD at 310 K, 100 ns), correlations in backbone torsion
angles between residues 9-18 and residues 35-40 (two loops connected by a
disulfide bond) were captured, matching millisecond-scale Shaw et al. results.

### Structure refinement

GNEIMO is implemented in CYANA (NMR refinement) and NIH-XPLOR (X-ray refinement).
Homology/template modeling requires >60% sequence similarity to the template.
GNEIMO-REMD with generalized Born solvation was evaluated on 30 CASP target
proteins and other small proteins: refinement of up to 1.3 Å for 28 of 30 CASP
targets, without experimental restraints. For target T0453, long-range contacts
between loop residues 30-40 and residues 40-60 improved from 14-16 Å to 2-4 Å.
Unlike unconstrained Cartesian REMD (which unravels/reforms secondary structure),
GNEIMO refines loop packing without unraveling structures. Overall refinement was
modest and needs improvement; selecting the best refined structure remains open.

## 4. Conclusions

ICMD is not a replacement for Cartesian MD. The vision is an ICMD method in BAT
coordinates suitable for selecting/controlling model granularity, as a foundation
for multiscale simulation of protein macromolecular complexes and polymers, and
for real-space refinement against low-resolution crystallography/EM data. GNEIMO
(and GneimoSim) is a robust long-time-scale ICMD toolkit; demonstrated on protein
domain motion and homology-model refinement.

Possible extensions:

1. Accurate conformational entropy from GNEIMO torsional MD trajectories via
   quasi-harmonic analysis (QHA), including metric tensor correction terms.
2. Allow movement of bond angles strongly coupled to dihedrals (e.g. the bond
   angle hinged on Cα atoms), yielding hybrid ICMD models bridging coarse-grain
   torsional MD and fine-grain all-atom MD; hybrid models with all bond angles
   open can also compute conformational entropy accurately.
3. Combine GNEIMO with other enhanced sampling (steered MD, umbrella sampling,
   accelerated MD) and with torsional Monte Carlo for structure prediction.
