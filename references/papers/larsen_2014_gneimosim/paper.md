# GneimoSim: A Modular Internal Coordinates Molecular Dynamics Simulation Package

Adrien B. Larsen, Jeffrey R. Wagner, Saugat Kandel, Romelia Salomon-Ferrer, Nagarajan Vaidehi, Abhinandan Jain. *J. Comput. Chem.* 2014, 35, 2245-2255. DOI: 10.1002/jcc.23743

> Routing-only note: this is a software-description paper. It contains no
> implementable equations; its value is pointing the implementer at the
> underlying GNEIMO method papers (SOA linear-cost ICMD solver, the ICMD
> equipartition principle, and the Fixman correction). See `depends_on` in
> `meta.yaml`.

## Abstract

The generalized Newton-Euler inverse mass operator (GNEIMO) method is an
advanced method for internal coordinates molecular dynamics (ICMD). GNEIMO
includes several theoretical and algorithmic advancements that address
longstanding challenges with ICMD simulations. This article describes the
GneimoSim ICMD software package that implements the GNEIMO method. GneimoSim is
described as the first software package to include advanced features such as the
equipartition principle derived for internal coordinates, and a method for
including the Fixman potential to eliminate systematic statistical biases
introduced by the use of hard constraints. GneimoSim is extensible and can be
interfaced with third party force field packages for ICMD simulations. It
includes interfaces to LAMMPS, OpenMM, and Rosetta force field calculation
packages. A comprehensive Python interface to the underlying C++ classes and
methods lets users write simulation scripts to configure and control the
simulation flow.

## Introduction (context and motivation)

All-atom (Cartesian) MD uses absolute coordinates with a simple dynamics model,
but this simplicity suffers when constraints and/or bias potentials are added.
Adding hard constraints requires differential-algebraic equation (DAE) solvers
that adversely impact robustness and complexity.

Bond, Angle, Torsion (BAT) relative coordinates are more natural for the bonded
structure of a protein. MD in BAT coordinates is called internal coordinate
molecular dynamics (ICMD). In the ICMD model, high-frequency bond-length degrees
of freedom can be constrained by simply excluding them. The resulting models are
of smaller dimension and retain ordinary-differential-equation (ODE) structure
instead of the DAE structure required for constrained Cartesian models. Torsional
molecular dynamics (TMD) is the ICMD model where bond length and bond angle
coordinates are frozen and only torsional degrees of freedom are free.

In ICMD models the six translation/orientation degrees of freedom for a molecule
are explicit coordinates rather than implicit as in Cartesian models. TMD
advantages: (1) low-frequency torsional coordinates allow larger time steps;
(2) conformational search in low-frequency torsions leads to significant
conformational changes; (3) enhanced sampling is more effective in torsional
space.

Longstanding ICMD challenges the GNEIMO method addresses:

1. Strong coupling among BAT coordinates increases analytical complexity of the
   dynamics model.
2. Computational cost of the dynamics solution grows cubically with the number
   of degrees of freedom, limiting scalability.
3. Rigidity from freezing degrees of freedom affects transition barriers and the
   probability density function of conformational states. Fixman proposed a
   compensating potential that rigorously corrects this bias, but it has been too
   complex to compute for even moderate sized molecules.
4. Availability of ICMD software and application examples is very limited.

## The GNEIMO ICMD method

In the GNEIMO ICMD model, each molecule is modeled as a collection of rigid
bodies (termed **clusters**) connected by one-to-six degree-of-freedom **hinges**.
A cluster is a group of atoms that move as a rigid unit; it can be a single atom,
a methyl group, a phenyl ring, an alpha helix, or an entire protein domain.
Different cluster choices control the granularity of the ICMD model. A hinge can
have one (TMD) to six degrees of freedom, and hinges can be frozen or thawed
during a simulation to decrease or increase the number of degrees of freedom.
The default clustering removes high-frequency modes of motion, leaving the
relatively low-force torsion terms of the force field to act on the model,
allowing longer time steps.

GNEIMO advancements:

1. The spatial operator algebra (SOA)-based GNEIMO ICMD algorithm reduces the
   cost of solving the ICMD equations of motion to **linear** (O(N)) instead of
   cubic in the number of degrees of freedom. Graph-theory extensions generalize
   the mass-matrix factorization for the low-cost recursive solution.
2. GNEIMO uses a new **equipartition principle** that generalizes the classical
   equipartition principle to ICMD models. This forms the basis for "modal
   velocity coordinates" that give a rigorous method for thermodynamically
   correct initialization of velocities in ICMD simulations.
3. GNEIMO includes a low-cost, general-purpose SOA-based algorithm for including
   the **Fixman correction potential** for bias-free ICMD simulations. Inclusion
   of the Fixman potential recovers the equilibrium probability density function
   of conformational states, transition barrier crossing rates, and the free
   energy surface for general serial and branched polymers (verified without
   forcefields).
4. GNEIMO ICMD has been expanded beyond TMD to allow freeing up bond angle
   degrees of freedom. These hybrid ICMD models reduce the rigidity encountered
   with TMD models.

The GNEIMO method has been tested on proteins with 40-300 residues for long time
scale dynamics (500 ns-1 ms each), used for protein homology model refinement
(refines up to 1.5 A without additional experimental restraints), and combined
with temperature-based replica exchange MD (REMD) for enhanced conformational
sampling (calmodulin, fasciculin) and ab initio folding of simple proteins.

## The GneimoSim ICMD software

Design philosophy: implement high-quality ICMD algorithms and leverage (rather
than reimplement) established forcefield, solvation, and enhanced-sampling
modules from the MD community. Additional included algorithms:

1. Extension of the Nose-Hoover thermostat (NVT) for ICMD.
2. Several integrators: Runge-Kutta, Lobatto, adaptive CVODE, and Verlet.
   Stability and Hamiltonian conservation verified for 40 proteins, 30-300
   residues.
3. Replica Exchange MD (REMD) for enhanced sampling.
4. Accelerated molecular dynamics (aMD) enhanced sampling.
5. Langevin dynamics.
6. Generalized Born Solvation (GBSA) module (from Simbios); periodic boundary
   conditions for explicit solvent.

Core is implemented in C++ for speed; a comprehensive Python interface is
auto-generated with SWIG.

### GneimoSim architecture

- **Modularity.** Object-oriented, modular architecture based on a functional
  decomposition of the ICMD problem, allowing mixing/matching of component
  implementations.
- **Extensibility.** Interface third-party software without reimplementing it
  inside GneimoSim (e.g. OpenMM and Rosetta forcefield modules).
- **Configurability.** C++ core with a comprehensive Python interface for
  selecting ensembles, enhanced sampling methods, integrators and forcefield
  variants, and building simulation loops with conditionals and staging.
- **Computational speed.** Core dynamics uses the DARTS (dynamics algorithms for
  real time simulation) module, an SOA-based high-performance multibody dynamics
  solver implementing SOA's linear-cost algorithm for the ICMD equations of
  motion.

### Functional decomposition and base classes

Base classes establishing GneimoSim's functional interfaces:

- **Gneimo**: top-level simulation manager class and entry point. Instantiated
  with data defining the ICMD molecular model plus configuration options; has
  methods to run, load and save a simulation. Constructor requires files
  describing masses/number of atoms (system), starting coordinates (coords),
  molecule topology, and ICMD clusters. Options select model type, integrator,
  and ensemble type. Methods add forcefield(s), toggle features, initialize
  temperature, and run for a number of steps.
- **GneimoModel**: defines the clustering ICMD model. Options are the ICMD
  (internal coordinates) model or a Cartesian model (each atom is its own
  cluster, enabling all-atom MD and comparison with Cartesian MD packages).
- **GneimoIntegrator**: base class for numerical integrators (native
  implementations plus third-party CVODE from SUNDIALS).
- **GneimoForceField**: base class for forcefields; computes the net force on
  each atom given current model coordinates. Integrator queries forces to
  determine accelerations. Multiple forcefields can be instanced and combined
  with user-specified weights.
- **GneimoReporter**: records data (Cartesian coordinates, dihedral velocities,
  energies, diagnostics). Any number of reporters can be used.

### DARTS ICMD dynamics solver

DARTS contains the ICMD model and SOA-based algorithms: solving the ICMD
equations of motion, kinetic energy computation, NVT bath dynamics terms, and
Fixman potential and torque computations. The equations-of-motion solver is a
recursive algorithm whose cost grows linearly with the number of degrees of
freedom, and it allows freezing/thawing of individual degrees of freedom during
a simulation. DARTS also contains velocity-initialization algorithms based on the
constrained-system equipartition principle.

Fixman's compensating potential depends on the mass matrix of the constrained
system and also induces additional torques at the hinges. GneimoSim implements a
computationally efficient SOA-based algorithm for the Fixman potential and
associated torque for general branched molecules; it reuses several terms from
the DARTS solver, so including the Fixman correction has only a modest impact on
overall cost.

## GneimoSim configuration options

### ICMD models and clustering

A GneimoModel instance defines an ICMD model of clusters connected by hinges.
Internal atom positions within a cluster are set after conjugate-gradient
minimization during structure preparation (sets harmonic bond lengths and bond
angles to equilibrium values). Proline rings are handled as open tree structures
(the ring is broken into clusters, not kept rigid as a whole). In the Cartesian
model each atom is its own cluster.

The default ICMD cluster model file (generated by a script from an input PDB
file) defines a torsional dynamics model with all bond lengths and bond angles
frozen. Users can edit the cluster file to change cluster definitions (freeze
whole domains, sample torsions connecting domains). Hinges can be frozen and
thawed at any time during a simulation ("dynamic clustering"), enabling adaptive
coarse-graining, e.g. freezing backbone torsions of a formed helical turn at high
temperature replicas and thawing at lower temperature replicas. GneimoSim also
supports freeing up any or all bond angle degrees of freedom. Multiple chains
(multi-subunit proteins, explicit water) are supported. All-atom GneimoSim
simulations were cross-validated against LAMMPS and OpenMM.

### Force field modules

A forcefield is integrated by deriving an interface class from GneimoForceField
that transfers coordinates and forces between the Gneimo manager and the
forcefield module. Available interfaces:

- **GneimoForceFieldLammps**: LAMMPS forcefields (fastest CPU implementation
  used); select between standard Amber and CHARMM parameters.
- **GneimoForceFieldLammpsGBSA**: CPU Generalized Born solvent model, works with
  GneimoForceFieldLammps.
- **GneimoForceFieldOpenMM**: OpenMM forcefields; recommended for CUDA GPUs
  (significant force-calculation speedup); includes GBSA.
- **GneimoForceFieldRosetta**: Rosetta forcefield (Rosetta energy function), used
  for structure prediction and homology-model refinement.
- **GneimoForceFieldSpring**: custom harmonic potential for user-defined pairwise
  constraints (NOE restraints, steered dynamics harmonic restraint).
- **GneimoForceFieldLangevin**: temperature-dependent random forces for Langevin
  dynamics.

Usable forcefields include AMBER99SB, CHARMM, and Rosetta. Multiple forcefields
can be registered; forces are combined with user-specified weights (e.g. LAMMPS +
LAMMPS-GBSA for implicit solvent added to an all-atom forcefield).

### Solvent model

Generalized Born implicit solvation (via GneimoForceFieldLammpsGBSA) and explicit
solvent models are supported. Explicit solvent works with both LAMMPS and OpenMM
modules. Explicit water molecules are treated as individual three-atom rigid-body
cluster bodies in ICMD mode, with parameters from a common forcefield (e.g.
Amber), input via the cluster file. GneimoSim uses the forcefield module's Ewald
summation and periodic boundary conditions.

### Ensembles

Microcanonical (NVE) and canonical (NVT) ensembles are included. NVT via the
Nose-Hoover thermostat, plus a Berendsen thermostat and Langevin dynamics.
Selected via Gneimo methods such as `runNVE` or `runNoseHoover`.

### Integrators

Fixed-time-step Lobatto and Runge-Kutta 4 (RK4); a Brunger-Brooks-Karplus
integrator for Langevin dynamics (with GneimoForceFieldLangevin); and the CVODE
adaptive-time-step integrator (most suitable for stable long-time-scale TMD
simulations, tested on folding small proteins).

### Enhanced sampling methods

REMD (adds replicas at higher temperatures to escape barriers) and aMD
(especially successful applied specifically to dihedral torsions; GneimoSim TMD
supports aMD applied to torsional modes).

### Input/output

Model initialization requires starting Cartesian coordinates, atom masses, and
cluster topology. Starting coordinates from PDB, Amber crd, or LAMMPS data file.
Masses from a LAMMPS data file, an OpenMM serialized XML system, or a simple
gsystem mass set (generated from a PDB file with an included script). Output logs
are ASCII; trajectories in uncompressed crd or binary DCD (readable by VMD,
PyMOL, MDAnalysis). All output goes through GneimoReporter instances.

## GneimoSim usage (example script, semantics)

A simulation is a Python script that creates and configures GneimoSim class
instances during setup, then runs an execution loop and logs data. Example: an
NVT simulation of a 20-alanine polypeptide with the Lobatto integrator.

1. Create the Gneimo simulation manager, selecting system source
   (`SYSTEM_LAMMPS`), coordinate source (`COORDS_LAMMPS`), the cluster file, the
   model type (`MODEL_ICMD` selects the standard ICMD cluster model), and the
   integrator (`INTEGRATOR_LOBATTO`).
2. Create and register forcefield instances (e.g. a LAMMPS forcefield plus a
   LAMMPS-GBSA forcefield) via `g.addForceField(...)`. Registered forcefields are
   all processed during the simulation.
3. Create reporters (log reporter, DCD trajectory reporter) with output filenames
   and frequencies.
4. Initialize velocities for the target temperature per the ICMD equipartition
   principle: `g.initTemperature(temperature=300, random_seed=111)`. The random
   seed makes velocity initialization repeatable or intentionally distinct across
   parallel runs.
5. Run with the ensemble-specific run command, e.g. a Nose-Hoover NVT run:
   `g.runNoseHoover(step_size=5.0, num_step=200000, cm_reset_frequency=100,
   bath_temperature=300.0, temperature_relax_scale=500)` — a 200,000-step run at
   5.0 fs steps, resetting overall linear/angular momentum every 100 steps to
   remove numerical drift, maintaining 300 K with a bath relaxation coefficient
   of 500 fs.

Variants: swap in OpenMM GPU-accelerated AMBERff (`GneimoForceFieldOpenMM`) or
the Rosetta forcefield (`GneimoForceFieldRosetta`) for torsional-MD homology
refinement. Bond angle degrees of freedom can be freed at runtime, e.g.
`cluster.addBondAngleDof(atom_number=3)` on a cluster obtained from the
GneimoModel.

## GneimoSim performance

Performance is analyzed via the cost of the ICMD equations of motion and the cost
of the force calculation. A series of PDB proteins of increasing size was used;
the largest system (7165 clusters) is human alpha-2-macroglobulin, 20,426 atoms.
For each protein, five independent 300 ns NVE ICMD simulations were run with the
Lobatto integrator, 1 fs time step, 17 A cutoff for long-range interactions, no
explicit solvent. Reported run time is the average time per step over the last
100 ns of each of the 5 simulations. Hardware: Intel Xeon E5-2670 CPU and one
Nvidia Tesla K20m GPU.

The dynamics-only run time (no forcefield) shows effective linear (O(N)) scaling
of the GNEIMO dynamics solver. Force calculation (OpenMM+GBSA on GPU; LAMMPS and
Rosetta on CPU) scales at a higher order; for LAMMPS Amber/Rosetta on CPU the
force cost dominates from a small number of clusters. The relative cost of the
dynamics decreases with increasing system size while force computation cost
increases. The higher up-front ICMD cost is compensated by larger time steps:
GneimoSim has run stably with integration steps as large as 10 fs, versus ~2 fs
for all-atom Cartesian MD with SHAKE constraints on bond lengths. For CPU
implementations, larger time steps give similar simulation speed for ICMD and
Cartesian MD. On GPUs, Cartesian simulations are significantly faster because no
ICMD implementations on the GPU are available to date.

## GneimoSim applications

- **Folding / torsional MD.** Folding of small proteins from extended structures
  with NVT Nose-Hoover dynamics and REMD (12 temperature replicas), Lobatto
  integrator. Dynamic clustering demonstrated on PDB IDs 1BDD (res. 11-56), 1EON
  (res. 7-31), 1PRB (res. 11-53), and 1UB — helical regions kept rigid, the rest
  flexible. Trp-cage folding with partially formed helices frozen as the
  simulation proceeded (CVODE integrator) showed better sampling of near-native
  structures than all-torsion constrained MD. Four proteins folded from extended
  to molten-globule-like native structures within 4-5 A of the crystal.
- **Protein structure refinement.** GneimoSim + temperature REMD + NVT Hoover +
  Lobatto refined homology models for various CASP targets without knowledge of
  the predicted structures; GNEIMO TMD led to refinement in most cases with REMD,
  in contrast to all-atom MD (which does not refine without known-structure
  restraints).
- **Large-scale conformational changes.** GNEIMO TMD + REMD sampled two
  experimentally established conformational substates of fasciculin, and the
  calmodulin Ca2+-bound to Ca2+-free transition occurred readily. Unconstrained
  all-atom Cartesian simulations failed to sample these transitions.

## Conclusions

GneimoSim is a computationally efficient implementation of the GNEIMO method for
ICMD simulations, described as the first package to include the ICMD
equipartition principle and methods for including the Fixman potential.
Its extensible architecture interfaces with third-party forcefield packages
(LAMMPS, OpenMM, Rosetta) and exposes a comprehensive Python interface to the C++
classes. Planned work: support the full range of ICMD models (including free
nontorsional degrees of freedom), an analysis toolkit operating directly on BAT
coordinate trajectories (e.g. entropy computations), and integration with
third-party specialized forcefields and enhanced-sampling modules.
