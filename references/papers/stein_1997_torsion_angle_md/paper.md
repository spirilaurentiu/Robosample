# Torsion-Angle Molecular Dynamics as a New Efficient Tool for NMR Structure Calculation

Evan G. Stein, Luke M. Rice, Axel T. Brünger (1997), *Journal of Magnetic Resonance* series B.

## Abstract

Molecular dynamics in torsion-angle space was applied to NMR structure
calculation using NOE-derived distances and J-coupling-constant-derived dihedral
angle restraints. Compared to molecular dynamics in Cartesian space (SA) and
metric-matrix distance geometry combined with Cartesian molecular dynamics
(DGSA), the method shows increased computational efficiency and success rate for
large proteins, and a dramatically increased radius of convergence for DNA. The
algorithm starts from an extended strand conformation and proceeds in four
stages: high-temperature torsion-angle MD, slow-cooling torsion-angle MD,
Cartesian MD, and minimization. For villin 14T (126 residues) the success rate
is 85%, more than twofold over other methods. For a 12 base-pair DNA duplex the
torsion-angle success rate was 52% while Cartesian MD and metric-matrix distance
geometry always failed.

## Introduction

NMR structure calculation aims to simultaneously satisfy experimentally observed
NMR data (NOEs, J coupling constants, chemical shifts) and chemical information
(stereochemistry, nonbonded interactions). Methods based on metric-matrix
distance geometry or Cartesian-space MD can show low success rates for large or
poorly determined systems, and are often combined to improve convergence.

Fixed-length and fixed-angle constraints reduce the number of adjustable
parameters roughly 10-fold, improving the observable-to-parameter ratio. These
constraints are most useful at high simulated-annealing temperatures where
conventional MD allows significant deviations from ideal geometry. Efficient and
robust algorithms for torsion-angle-space MD have become available (Bae & Haug;
Jain, Vaidehi & Rodriguez; Rice & Brünger; Mathiowetz et al.).

This paper describes an NMR structure-calculation method using MD constrained to
torsion-angle space, converging from an extended strand. Stage 1: an initial
high-temperature search of torsion-angle space with decreased weight on the
repulsive energy term. Stage 2: torsion-angle dynamics with temperature
gradually reduced while the repulsive-term weight is linearly increased to unity.
Stages 3 and 4: bond lengths and bond angles are allowed to relax (Cartesian MD
then minimization).

## Methods

### Energy function

NMR structure calculation is formulated as a hybrid-energy-function optimization
problem. The total energy is the sum of a chemical term and an NMR-restraint
term (Eqs. [1]-[3]). `E_chem` describes agreement with expected values for bond
lengths, bond angles, planarity, chirality, and nonbonded interactions (van der
Waals, hydrogen bonding, electrostatics). Because solvent is neglected,
electrostatic interactions are excluded and hydrogen bonds are modeled as
pseudo-NOEs.

Rather than the Lennard-Jones potential (Eq. [4]), the van der Waals
interactions were described by a purely repulsive quartic potential (Eq. [5]),
where `R` is the distance between two atoms and `epsilon`, `sigma` are the
Lennard-Jones parameters for a particular atom pair. For the final analysis of
refined structures, Eq. [4] (Lennard-Jones) was used. Force field parameters
were taken from parameter sets designed for NMR refinement of proteins
(parallhdg.pro) and nucleic acids (parallhdg.dna).

NOE-derived distance restraints were described by a flat-bottomed parabolic
(square-well) function with a soft asymptote (Eq. [6]), where `Delta` is defined
by Eq. [7]. `R` is the distance between a particular pair of spins in the model,
`d_lower` and `d_upper` are the lower and upper distance bounds, and `a`, `b`
are determined such that `E_NOE` is differentiable at `R = d_upper + 0.5`. The
sum runs over all NOEs.

Dihedral angle restraints derived from J coupling constant measurements were
described by Eq. [8]. The method in principle allows other functional forms:
direct refinement against NOEs, J-coupling values via the Karplus equation, and
restraints derived from chemical shifts.

### Molecular Dynamics

Structure calculation based on MD consists of the numerical integration of
Newton's equations of motion (Eq. [9]), where `r_{i,u}` and `m_i` are the
coordinates and mass of atom `i`, and `E` is the hybrid energy function (Eq. [1]).
Temperature control, required for simulated annealing, was performed by
temperature (Berendsen) coupling (Eq. [10]), where `T_0` is the bath temperature,
`beta_i` is a force constant, `T` is the system temperature, and `v_i` is the
velocity of atom `i`. Temperature coupling causes "heat" to be added or removed
(as kinetic energy) as needed to maintain the temperature.

### Torsion-Angle Molecular Dynamics

What follows is a simplified sketch of one implementation of torsion-angle
constrained molecular dynamics, following the algorithm of Bae and Haug. Consider
two bodies `i` and `j` connected by a bond of fixed length `|h_ij|` (see Fig. 1).
Let `r_i` and `r_j` locate (with respect to an arbitrary inertial "lab" frame)
the centers of mass of bodies `i` and `j`. Let `s_ij` (`s_ji`) locate the
endpoint of `h_ij` on body `i` (`j`) with respect to its center of mass; thus
`s_ij` is a vector from the center of mass of body `i` to the end of `h_ij`. The
position of the center of mass of body `j` with respect to body `i` is
`r_ij = r_j - r_i`. The scalar `q_ij` measures the relative angle of rotation
about the bond `h_ij`.

The assumption that the only allowable relative motion between the two bodies is
a rotation about the connecting bond implies a relationship between the angular
velocities of their centers of mass measured in the inertial frame (Eq. [11]),
where `q̇_ij` (denoted `q_ij` with a dot) is the time derivative of the relative
angle and `ĥ_ij = h_ij / |h_ij|` is the unit vector along the bond.

The expression for `r_j` (Eq. [12]) can be differentiated and rearranged,
resulting in an expression for the center-of-mass velocity of body `j` in terms
of that of body `i` (Eq. [13]). Thus, assuming certain constraints act between
atoms or groups of atoms, one can obtain an expression for the velocity of one
group in terms of the velocity of another. This relationship can be
differentiated to give a relationship between accelerations, and integrated to
give a relationship between positions.

The current implementation cannot treat nonrigid closed bonding networks
exactly; it introduces an approximation whereby one bond in the closed network is
allowed to vibrate. This could cause numerical instabilities at high temperatures
for nucleotide ribose rings, so lower simulation temperatures are required for
nucleic acids.

### Derivation (not implemented): velocity relation between connected bodies

Starting from the position relation

$$\mathbf{r}_{j} = \mathbf{r}_{i} + \mathbf{r}_{ij} = \mathbf{r}_{i} + \mathbf{s}_{ij} + |\mathbf{h}_{ij}| \hat{\mathbf{h}}_{ij} - \mathbf{s}_{ji},$$

differentiating and using $\dot{\mathbf{s}}_{ij} = \boldsymbol{\omega}_i \times \mathbf{s}_{ij}$,
$\dot{\hat{\mathbf{h}}}_{ij} = \boldsymbol{\omega}_i \times \hat{\mathbf{h}}_{ij}$,
and $\dot{\mathbf{s}}_{ji} = \boldsymbol{\omega}_j \times \mathbf{s}_{ji}$ with
$\boldsymbol{\omega}_j = \boldsymbol{\omega}_i + \hat{\mathbf{h}}_{ij}\dot{q}_{ij}$:

$$\dot{\mathbf{r}}_{j} = \dot{\mathbf{r}}_{i} + \dot{\mathbf{s}}_{ij} + |\mathbf{h}_{ij}| \dot{\hat{\mathbf{h}}}_{ij} - \dot{\mathbf{s}}_{ji}$$
$$= \dot{\mathbf{r}}_{i} + \boldsymbol{\omega}_{i} \times \mathbf{s}_{ij} + |\mathbf{h}_{ij}| \, \boldsymbol{\omega}_{i} \times \hat{\mathbf{h}}_{ij} - \boldsymbol{\omega}_{j} \times \mathbf{s}_{ji}$$
$$= \dot{\mathbf{r}}_{i} - \mathbf{s}_{ij} \times \boldsymbol{\omega}_{i} - |\mathbf{h}_{ij}| \, \hat{\mathbf{h}}_{ij} \times \boldsymbol{\omega}_{i} + \mathbf{s}_{ji} \times \boldsymbol{\omega}_{i} - \dot{q}_{ij}\, \hat{\mathbf{h}}_{ij} \times \mathbf{s}_{ji}$$
$$= \dot{\mathbf{r}}_{i} - \mathbf{r}_{ij} \times \boldsymbol{\omega}_{i} - (\hat{\mathbf{h}}_{ij} \times \mathbf{s}_{ji})\, \dot{q}_{ij}.$$

<!-- CHECK: the sign/collection of the intermediate cross-product terms is reconstructed from corrupted OCR; the final compact form (Eq. [13]) is the implementable result. -->

## Test Cases

Structure calculations were carried out on protein G, BPTI, interleukin-8 (IL8),
villin 14T, and a 12 base-pair duplex of DNA (CGCGAATTCGCG) (Table 1). Nearly all
phi, psi, and chi1 dihedrals are restrained for both monomers of IL8; most phi,
chi1, chi2 dihedrals are restrained for protein G. Villin 14T has ~1 dihedral
restraint per residue. BPTI has no dihedral restraints (simulated NMR data). The
DNA dodecamer includes 136 dihedral and 30 hydrogen-bond restraints in addition
to NOE-derived distance restraints.

## Torsion-Angle Molecular Dynamics Protocol

Initial structures were extended strands generated by sequentially placing all
atoms along the x axis at 0.1 Angstrom intervals, with y and z coordinates set to
random numbers between zero and one. Initial coordinates were regularized using
simulated annealing and conjugate-gradient minimization against `E_chem` (Eq. [3])
to obtain good local geometry.

The protocol (Table 2) has four stages:

1. **High-temperature torsion-angle MD:** 15 ps at 50,000 K (20,000 K for nucleic
   acids) using the hybrid energy function E (Eq. [1]); `w_vdw` set to 0.1 to
   facilitate rotational barrier crossings. Time step 0.015 ps.
2. **Slow-cooling torsion-angle MD:** temperature reduced from 50,000 to 1,000 K
   over 15 ps while `w_vdw` linearly increased from 0.1 to 1.0. Time step 0.015 ps.
3. **Slow-cooling Cartesian MD:** 1,000 to 300 K for 6 ps. Time step 0.003 ps.
4. **Minimization:** 1000 steps conjugate gradient.

The protocol was repeated with different initial velocities drawn from a random
Maxwellian distribution to obtain an ensemble. The parameters that most affected
the protocol are the temperature and the duration of the torsion-angle stages.
Temperatures above 50,000 K accelerated convergence for large molecules but
lowered success for smaller molecules; 50,000 K was a good compromise for protein
structures with 10-15 NOE restraints per residue. The protocol suffices to refine
proteins of 1000-3000 atoms with ~10-15 NOE restraints per residue; larger
structures may need longer torsion-angle stages.

For nucleic acids, due to vibrations in nonrigid ribose rings: simulation
temperature reduced to 20,000 K for torsion-angle stages, the dihedral-restraint
weight (`E_cdih`) reduced from 100 to 5 during both torsion-angle stages, and the
length of both torsion-angle stages tripled.

A major advantage of the method is its simplicity: only four stages with two
changing parameters (`w_vdw` and `w_dihedral`). By contrast, SA has five stages
with many changing parameters, and DGSA has eight stages.

## Comparisons

The torsion-angle MD algorithm was compared to a Cartesian-MD-based simulated
annealing method (SA) and a protocol using metric-matrix distance geometry
combined with Cartesian MD (DGSA). Both comparison algorithms are implemented in
X-PLOR version 3.1 (files SA.INP and DGSA.INP) and were unmodified except to
extend the Cartesian MD stage for interleukin-8 and villin 14T by a factor of 4
to obtain a reasonable acceptance rate.

## Acceptance Criterion

The three algorithms were repeated with different initial velocities until each
produced 50 acceptable structures. An acceptable structure has no NOE-restraint
violations greater than 0.5 Angstrom and no dihedral angle violations greater
than 5 degrees. Structures were rejected if RMS deviation of bonds from ideal
exceeded 0.02 Angstrom, or RMS deviation of angles exceeded 2.0 degrees. Success
rate is the ratio of accepted structures to total trials. Computational
efficiency is the average computing time to obtain one acceptable structure.

## Results and Discussion

Dial plots were used to compare sampling of conformational space. The torsion-angle
dial plots are significantly better sampled than the other protocols. The RMS
difference from the average structure is approximately the same (within a standard
deviation) for all three methods (Table 3), and pairwise RMS differences between
average structures are similar (Table 4). The ensembles satisfy experimental and
chemical restraints to the same degree except that van der Waals energies are
lower for the torsion-angle MD structures (Table 5).

The success rate and computational efficiency of torsion-angle MD is higher than
the other two methods for larger proteins (Table 6). For interleukin-8 and villin
14T, SA takes 2-4x longer to generate an acceptable structure. DGSA is faster at
generating a single structure but takes about twice as long to generate an
acceptable one.

For the DNA dodecamer (Table 7), neither DGSA nor SA generated any acceptable
structures; torsion-angle MD produced acceptable structures in 52% of trials. The
torsion-angle ensemble agrees most closely with the original structure (RMS
deviation 2.67 Angstrom), while SA and DGSA helices deviate significantly and have
large NOE violations. When the original structure is subjected to a brief (6 ps)
Cartesian MD refinement and 1000-step minimization against the same energy
function used here, the ensemble moves away from the original (RMS 2.26 Angstrom)
toward the torsion-angle MD structure (RMS 1.25 Angstrom), showing the differences
are artifacts of a different energy function.

## Conclusions

Molecular dynamics constrained to torsion angles provides a powerful tool for NMR
structure calculation. It has a higher success rate and efficiency than
conventional SA or DGSA. A significant difference in computing time appears for
proteins larger than 100 residues. Furthermore, torsion-angle MD folds extended
DNA strands into B-form DNA with correct helicity without additional restraints
and without starting from A- or B-form DNA. As NMR-analyzed structures increase
in size, the advantage of torsion-angle MD is expected to become increasingly
important. All calculations were carried out with X-PLOR.
