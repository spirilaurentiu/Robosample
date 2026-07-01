# The Molecular Modeling Toolkit: A New Approach to Molecular Simulations

**Author:** Konrad Hinsen (Centre de Biophysique Moléculaire, CNRS, Orléans, France)
**Venue:** J. Comput. Chem. 21(2): 79–85 (2000)

> Routing-only note: this is a software/overview paper describing the MMTK
> library. It contains no implementable equations of its own; its value is in
> routing to the underlying methods (Amber 94 force field, velocity Verlet with
> constraints, Ewald summation, normal modes, etc.) and in documenting library
> design and capability. No `equations.md`, `notation.md`, or `checks.md` is
> produced.

## Abstract

The Molecular Modeling Toolkit (MMTK) is a library that implements common
molecular simulation techniques, with an emphasis on biomolecular simulations.
It uses object-oriented design and a high-level language (Python) to overcome
limitations of large monolithic simulation programs. Principal advantages:
(1) easy extension and combination with other libraries via modular design;
(2) a single high-level general-purpose language (Python) for both library
implementation and application scripts; (3) documented, machine-independent
formats for all data files; and (4) interfaces to other simulation and
visualization programs. Description is based on release 2.0.

## Design philosophy

MMTK is built on two software-engineering techniques: object-oriented design
(software organized around data structures, with procedures attached to the
data they operate on) and high-level languages (data structures matched to the
problem rather than the machine). MMTK is written in Python; time-critical
sections such as force-field evaluation are written in C, giving a
mixed-language design where only short time-critical parts are low-level. The
highly modular design allows new techniques and new force fields to be added
without modifying existing code.

## Object-oriented representation of chemical systems

The core is a set of classes describing chemical systems: atoms connected by
bonds form molecules, which combine into complexes. A `Group` class represents
functional groups. Specializations exist, e.g. proteins are complexes of
peptide chains, which are molecules made of amino-acid residue sequences; a
peptide chain can be asked for its fifth residue.

Example (create two water molecules, compute a distance):

```python
from MMTK import *
molecule_1 = Molecule('water', position=Vector(0., 0., 0.))
molecule_2 = Molecule('water', position=Vector(1., 0., 0.))
molecule_2.rotateAroundCenter(Vector(0., 0., 1.), 45*Units.deg)
distance = molecule_1.O.position() - molecule_2.H1.position()
print distance.length()/Units.Ang
```

Distances are in nm internally (the shift of 1 nm along x); `Units.Ang`
converts to Ångströms. A molecule definition is a small Python program in a
database, e.g. a minimal water:

```python
O = Atom('O')
H1 = Atom('H')
H2 = Atom('H')
bonds = [Bond(O, H1), Bond(O, H2)]
```

Biomolecules are usually built from PDB files. A protein can be built
automatically, creating a peptide chain per chain, placing missing hydrogens,
and detecting sulfur bridges. Models offered: all-atom, polar-hydrogen-only,
no-hydrogen, and (for proteins) C-alpha-only.

```python
from MMTK.Proteins import Protein
lysozyme = Protein('135l.pdb')
from MMTK.ForceFields import Amber94ForceField
universe = InfiniteUniverse(Amber94ForceField())
universe.addObject(lysozyme)
print universe.energy()
```

A `universe` represents a complete system: molecules, a geometry (infinite or
periodic), an optional force field, and environment objects (thermostats,
barostats). Method dispatch (`universe.energy()`) chooses the algorithm based
on object type (infinite vs. periodic).

Hierarchical access uses indexing (indices start at zero), e.g.
`lysozyme[0][4].sidechain.C_beta` for the C-beta of residue 5 of chain 1, or
`lysozyme[0][0:6].mass()` for the mass of residues 1–6.

## Data storage

Any object (or combination) can be stored in files using the documented,
machine-independent format from the Python standard library (pickle).
Trajectories, produced by MD integrators, are generated directly as disk files
(too large for memory) and use the netCDF format via the netCDF library.
netCDF files are machine-independent and self-describing; an MMTK trajectory
contains a complete description of the chemical system. Trajectories support
access by time step and access by quantity (e.g. one atom's position across all
steps). MMTK also reads/writes PDB and the DCD trajectory format used by
CHARMM and X-PLOR.

## Force fields

Force fields are objects defined by classes; adding a new force field means
writing classes without recompiling existing code. Force fields can be combined
and defined as modifications of existing ones. Force-field evaluation is in C
for speed, but terms can be implemented in Python, where energy gradients and
second derivatives can be obtained by automatic differentiation.

Implemented force fields:
- Amber 94 force field (most extensive).
- A deformation force field for normal-mode calculations in large proteins.
- A simple Lennard–Jones force field for noble gases (illustration).

For all force fields, first and second derivatives of the potential energy are
available.

Electrostatics options:
- direct evaluation of all N^2 pair interactions;
- direct evaluation with a cutoff and charge neutralization (Wolf method);
- Ewald summation for periodic systems;
- fast-multipole method via an interface to the DPMTA library.

## Algorithms

Provided algorithms include:
- Energy minimization: steepest descent and conjugate gradients.
- Molecular dynamics: velocity Verlet adapted to systems with optional
  distance constraints and optional thermostat/barostat, based on the equations
  of motion of Kneller and Mülders.
- Low-frequency normal modes using Fourier bases; standard normal modes; normal
  modes in arbitrary subspaces; sparse force-constant matrices; force constants
  by numerical differentiation.
- Quaternion-based structure superposition fits.
- Stable SVD-based partial-charge fits.
- Molecular surface calculations via the Analytic Surface Calculation Package.

The emphasis is on facilitating implementation and testing of new algorithms.
Many analysis algorithms are a few lines of Python using MMTK plus Numerical
Python (arrays, linear algebra, FFT) and Scientific Python (visualization,
statistics, geometry, interpolation, automatic differentiation).

## Visualization

MMTK uses external programs for visualization (direct OpenGL planned). Any
program reading PDB works; special support exists for VMD (animations from
trajectories, normal modes). VRML output (VRML 1 and VRML 97) can be written or
fed to a VRML browser, generated via an intermediate graphics module that maps
basic graphics objects (lines, spheres) to VRML; other formats can be added by
writing an equivalent module.

## Applications

Application programs built on MMTK include DomainFinder (identifying dynamical
domains in proteins) and a new version of nMOLDYN (neutron-scattering quantities
from MD trajectories). Typical usage is a Python script combining MMTK
algorithms with newly developed ones. Standard application areas: MD in the NVE,
NVT, and NPT ensembles (with optional distance constraints); energy
minimization; normal-mode calculations; and analysis of macromolecular systems.

## Computational methods implemented in MMTK 2.0 (summary)

- Representation: atoms, functional groups, molecules, complexes, peptide
  chains, proteins, nucleotide chains, collections, universes (infinite and
  periodic), thermostats, barostats, mechanical constraints (distance
  constraints and fixed atoms), scalar/vectorial atom properties,
  configurations, force fields.
- Inquiry functions (any subset): number of atoms, degrees of freedom, total
  mass, total charge, dipole moment, distances, angles, dihedral angles, center
  of mass, tensor of inertia, bounding box/sphere, molecular surface and volume,
  RMS distance, potential energy, forces, force constants, kinetic energy,
  temperature, momentum, angular momentum, angular velocity.
- Universe inquiry: elementary cell shape and volume, reciprocal basis vectors,
  maximal Cartesian distance.
- Coordinate manipulation: translation, rotation, general linear coordinate
  transformations, rigid-body superposition.
- Object selection: geometrical (boxes, spherical shells) or via any inquiry
  function.
- Random numbers: random points in a universe, random directions, random
  rotations, random atom velocities.
- Energy minimization: steepest descent, conjugate gradients.
- Molecular dynamics: velocity Verlet with optional pair distance constraints,
  Nosé thermostat, Andersen barostat, velocity scaling, heating, removal of
  global motions.
- Normal modes: standard, in arbitrary subspaces, Fourier bases for
  low-frequency modes, sparse force-constant matrices and force-constant
  calculation by numerical differentiation.
- Trajectory operations: generation for complete/partial systems, writing
  individual configurations, output during minimization and MD, reading by step,
  reading by atom, extraction of rigid-body motions.
- Interfacing: I/O to/from PDB and CHARMM/X-PLOR (DCD) trajectory formats;
  output to VRML; visualization with PDB/VRML viewers; animations with XMol and
  VMD.
- Miscellaneous: charge fits to electrostatic potential surfaces, deformation
  analysis, motion analysis by subspace projection, tensor fields from atomic
  quantities, solvation of macromolecules.
