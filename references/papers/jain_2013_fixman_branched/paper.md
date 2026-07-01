# Fixman compensating potential for general branched molecules

Abhinandan Jain, Saugat Kandel, Jeffrey Wagner, Adrien Larsen, Nagarajan Vaidehi
J. Chem. Phys. 139, 244103 (2013). DOI: 10.1063/1.4851315

## Abstract

Constraining high-frequency modes of molecular motion increases simulation time
scale and improves conformational sampling in MD. However, constraints on
high-frequency modes such as bond lengths and bond angles stiffen the molecular
model and introduce systematic biases in the statistical behavior of the
simulations. Fixman proposed a compensating potential to remove such biases in
thermodynamic and kinetic properties. Previous implementations were limited to
short serial-chain systems. This paper presents a spatial-operator-algebra (SOA)
based algorithm (GNEIMO-Fixman) to compute the Fixman potential and its gradient
(the Fixman torque) within constrained dynamics for branched-topology molecules
of any size. Numerical studies validate the algorithm by recovering the dihedral
angle probability distribution function for systems from serial chains to protein
molecules. The Fixman potential recovers the free energy surface of a serial-chain
polymer, annulling the biases from constraining bond lengths and bond angles, at
only a modest increase in computational cost. This is the first use of the Fixman
potential for general branched systems.

## I. Introduction

Rigid constraints on higher-frequency DOF (bond lengths, bond angles) are common
in MD; SHAKE, RATTLE, and torsional MD are examples. Torsional MD speeds up
simulations and enhances sampling of low-frequency torsional DOF, and provides a
way to vary model coarseness.

The use of rigid constraints leads to systematic biases in the statistical behavior
of MD (Scheraga, Fixman and others). Fixman proposed adding a compensating
potential - the *Fixman potential* - to remove these biases. The Fixman potential
depends on the determinant of the mass matrix for the constrained molecular model.

The mass-matrix determinant is computationally difficult: large, structurally
complex, and configuration dependent. Prior numerical studies were therefore mostly
limited to small idealized serial chains, where the Fixman potential recovers the
uniform torsion-angle pdf. Patriciu et al. showed the maximum variation of the
Fixman potential grows with chain length. Echenique et al. found the Fixman
potential becomes significant for peptides with more than 2 residues. Brooks and
Abagyan proposed alternate ad hoc corrective torsional potentials that are system
dependent and not generalizable.

Previously the authors developed SOA techniques for the Generalized Newton-Euler
Inverse Mass Operator (GNEIMO) constrained MD method, whose cost scales linearly
with the number of DOF (vs cubic for prior methods). SOA techniques were then used
to extend GNEIMO to the GNEIMO-Fixman method for computing the Fixman potential and
its partial derivatives (the Fixman torque, the additional forces applied within
constrained MD). GNEIMO-Fixman works for branched molecules of arbitrary size at
modest additional cost. This work validates GNEIMO-Fixman on serial and branched
systems via Langevin dynamics.

## II. Computational methods

### A. Unconstrained (flexible) models

Cartesian (absolute) coordinates are common in all-atom MD. Curvilinear
bond/angle/torsion (BAT) coordinates are relative/internal coordinates; conformational
motion is dominated by torsions, so BAT coordinates focus models on dominant motion.
Rigid constraints can freeze/eliminate bond-length and bond-angle BAT coordinates
(torsional MD freezes all bond and angle coordinates, leaving only torsions).

Partition the 3n BAT coordinates into N unconstrained coordinates `alpha` and
`(3n-N)` coordinates `q` to be constrained. In the flexible model both q and alpha
vary; in the constrained model q is fixed at q_0 and only alpha varies.

With momenta p, the Hamiltonian has the form (see eq:1). The partition function
(eq:2), after integrating over momenta, yields the configuration-space partition
function (eq:3) with the det{M_B^{1/2}} factor, giving the flexible pdf (eq:4).

The BAT mass matrix factors as `M_B = J_B^* M_c J_B` (eq:5). Go and Scheraga
derived a simple closed form for det{J_B} (eq:6) that is independent of torsion
angles - it depends only on bond lengths and bond angles. Combining gives a simple
closed form for det{M_B} (eq:7). When U is torsion-independent (bond and angle
potentials only), the flexible torsion pdf is uniform, rho(beta_i) = 1/2pi (eq:8).

### B. Constrained models

Hard constraints freeze the (3n-N) coordinates q at q_0. The constrained partition
function (eq:9) and pdf (eq:10) involve det{M(alpha)} of the constrained mass
matrix M (the alpha-alpha sub-block of M_B). Unlike det{M_B}, det{M} DOES depend on
the torsion angles, so the flexible and constrained pdfs differ - the constraint
bias. Torsions that are uniform in the flexible model become non-uniform under
constraints.

Fixman proposed the modified potential U' = U + U_f (eq:11), with the Fixman
compensating potential U_f = (1/2)kT ln( det{M} / det{M_B} ). Substituting U' into
eq:10 recovers agreement with the flexible pdf (eq:12) at q=q_0.

When the unconstrained coordinates are only torsions, the M_B contribution is a
constant and U_f simplifies to `c_f + (1/2) ln det{M(alpha)}` (eq:13). The rest of
the paper assumes alpha are just torsion angles.

The key computational challenge is det{M}: M is dense and configuration dependent
with non-square Jacobians. Symbolic determinants are only feasible for very small
systems. Fixman's own serial-chain method exploits sparsity of the constrained
sub-block of M_B^{-1} but does not generalize to branched systems.

The Fixman torque for coordinate alpha_i is T(i) = -dU_f/dalpha_i (eq:14). No
general method existed for it before GNEIMO-Fixman.

### C. GNEIMO-Fixman method for calculating the Fixman potential

GNEIMO-Fixman (Jain, Ref. 27) evaluates the Fixman potential and torque for
arbitrary serial and branched constrained models using SOA factorization/inversion
of the constrained mass matrix M. The key SOA operator expressions (eq:15) are: the
Newton-Euler operator factorization `M = H phi M phi^* H^*`; the square-factor
factorization `M = [I+HphiK] D [I+HphiK]^*`; the closed-form inverse
`[I+HphiK]^{-1} = [I-HpsiK]`; and the resulting `M^{-1} = [I-HpsiK]^* D^{-1}
[I-HpsiK]`. Here H is the block-diagonal hinge articulation operator (torsional
axes), phi the lower-triangular rigid-body force-propagation operator (phi(i,j)
propagates spatial force from cluster j rigidly to cluster i), M the block-diagonal
link spatial-inertia operator, and D the block-diagonal articulated-body (AB)
inertia. The last expression is the basis of the linear-cost GNEIMO constrained
dynamics recursion.

From the square factorization, det{M} = det{I+HphiK}^2 det{D}, and since
det{I+HphiK}=1 and D is block-diagonal, det{M} = product of det{D(i)} (eq:16). All
D(i) are scalars except the 6x6 base cluster D(0). The D(i) are by-products of the
GNEIMO dynamics algorithm. Substituting into eq:13 gives the GNEIMO-Fixman working
formula for the Fixman potential (eq:17), computable at negligible extra cost for
branched molecules of arbitrary size.

### D. GNEIMO-Fixman method for calculating the Fixman torque

Combining eq:13 and eq:14 gives the torque as T(i) = -(1/2) d ln det{M}/dalpha_i
(eq:18). Using matrix-calculus identities (eq:19-21, including d ln det{X}/dX =
{X^*}^{-1}), Refs. 27 and 35 derive explicit SOA expressions for the Fixman torque:
`T(i) = -Trace{ P(i) Upsilon(i) H-tilde_omega^*(i) }` (eq:22), with all factors 6x6
matrices. A simpler equivalent form partitions the 6x6 product P(i)Upsilon(i) into
3x3 blocks: `T(i) = h^*(i) F[Q_11 + Q_22]` (eq:23), where F maps a 3x3 matrix A to
the 3-vector v whose skew (cross-product) matrix equals A - A^*. These hold for
branched systems of any size. Reference 35 generalizes eq:23 to the case where
alpha also includes bond-length coordinates.

The P(i) matrix is a by-product of the GNEIMO algorithm. Computing the 6x6
Upsilon(i) matrices requires an additional recursive scatter algorithm starting at
the base cluster, also using GNEIMO quantities, at linear cost. The overall cost of
the Fixman torque is linear and adds only marginally to GNEIMO cost.

## III. Results and discussion

### A. Numerical validation of the GNEIMO-Fixman potential and torque

Pear and Weiner derived a closed-form det{M(alpha)} for the C4 idealized 3-bond
serial chain with fixed bond lengths and 90-degree bond angles (eq:24). Comparing
the Fixman potential from eq:24 with GNEIMO-Fixman eq:17 shows excellent agreement
(Fig 1a). The Fixman torque from eq:23 matches the numerical derivative of the
eq:17 potential exactly for C4 (Fig 1b).

Patriciu et al. computed the Fixman potential for idealized serial chains. The
authors' Fixman-potential contour for the two torsions of the C5 four-bond chain
agrees excellently with Patriciu et al. (Fig 2), validating GNEIMO-Fixman for serial
systems. Since there is no prior Fixman-torque data, the torque expression (eq:23)
was validated against numerical differentiation of eq:17, with similar agreement for
larger systems (e.g. C5).

### B. Descriptions of the MD simulations

MD simulations study how effectively the Fixman potential recovers the thermodynamic
pdf of torsion angles for serial and branched systems. Langevin dynamics enhances
sampling; all simulations were isolated single molecules at 300 K with damping
constant 0.01/fs. Three simulation sets per system:

1. FLEXIBLE: flexible Cartesian Langevin with only bond-length and bond-angle
   potentials (nonbonded Coulomb and vdW off). U(alpha,q) is torsion-independent, so
   every torsion has the uniform pdf rho(alpha_i) = 1/2pi (eq:8).
2. TORSIONAL: constrained Langevin WITHOUT the Fixman potential; bond and angle
   coordinates are hard constraints (their potentials unused), nonbonded off.
   Expected torsion pdf given by eqs:10 and 16 -> eq:25 (biased).
3. FIXMAN: same as TORSIONAL but WITH the Fixman potential included; expected torsion
   pdf simplifies back to uniform, rho(alpha_i) = 1/2pi.

FLEXIBLE establishes the expected uniform pdf; TORSIONAL measures the bias; FIXMAN
assesses recovery. The flexible Langevin dynamics has the form of eq:26 and was
integrated with the BBK algorithm.

### C. MD simulation results for serial chains

FLEXIBLE, TORSIONAL, FIXMAN simulations were run for C4, C5, C11, C15 serial chains.
Each: timestep 1 fs, 50 ns total, coordinates/energy every 100 steps. Beads: equal
mass 14 amu, bond length 1.54 Angstrom, 90-degree bond angle. FLEXIBLE spring
constants: bond 83.66 kcal/A^2, angle 43.46 kcal.

For C4 (single torsion): FLEXIBLE shows uniform pdf; TORSIONAL shows a bimodal
distribution with maxima at approx +/-83 deg and minima at 0 and +/-180 deg; FIXMAN
recovers the uniform pdf. The bimodal TORSIONAL structure illustrates the constraint
bias, which affects torsion pdfs, conformational transition rates, and the free
energy surface.

Free energy G(alpha_i)/kT = -ln(rho(alpha_i)); the uniform level equals -ln(1/360).
The negative Fixman potential -U_f closely matches the TORSIONAL free-energy profile
(equal magnitude, opposite sign), and the FIXMAN net free energy agrees with
FLEXIBLE (Fig 4).

For C5, C11, C15 the Fixman potential likewise recovered the uniform distribution
for every torsion. For C5 and C11, the Fixman-potential global minimum corresponds
to a self-intersecting planar conformation (each torsion 0 deg) and the maximum to a
self-intersecting spatial conformation (each torsion approx +/-83 deg), consistent
with Patriciu et al.

### D. MD simulations for branched systems

Langevin FLEXIBLE, TORSIONAL, FIXMAN simulations were run on three realistic
branched peptides: alanine dipeptide, valine dipeptide, and the ten-residue peptide
chignolin.

#### 1. Small branched peptides

Peptides are modeled as rigid *clusters* (atoms with frozen bond lengths/angles)
connected by flexible torsion hinges. Terminal bonds and aromatic rings are rigid.
Alanine dipeptide -> 8 clusters; valine dipeptide -> 10 clusters. Bond/angle
parameters from AMBER99SB. Per set: 3 runs, timestep 1 fs, 100 ns each (300 ns
total), sampled every 200 steps.

For the C1-N2-Ca2-C2 torsion in both peptides, TORSIONAL gives bimodal
distributions while FIXMAN recovers the uniform FLEXIBLE pdf. RMS deviations of
torsion pdfs from uniform (R_flex, R_tor, R_Fix): alanine dipeptide R_tor=1.0e-3,
R_Fix=8.7e-5, R_flex=6.5e-5; valine dipeptide R_tor=8.4e-4, R_Fix=7.9e-5,
R_flex=6.6e-5. Joint pdfs of torsion pairs (Ca1-C1-N2-Ca2 vs C1-N2-Ca2-C2) are flat
for FLEXIBLE and FIXMAN and biased for TORSIONAL, so the Fixman potential recovers
joint as well as single torsion pdfs.

#### 2. Moderate sized branched peptide

Chignolin (ten-residue beta-hairpin), starting from NMR structure PDB 1UAO.
TORSIONAL, FLEXIBLE, FIXMAN Langevin simulations; terminal bonds and rings
clustered. Per set: 4 runs, timestep 1 fs, 100 ns each (300 ns total), sampled every
200 steps. Torsions C1-N2-Ca2-C2, C4-N5-Ca5-C5, C7-N8-Ca8-C8 (at different parts of
the peptide) are all bimodal under TORSIONAL and become nearly flat under FIXMAN.
pdf deviation measures: R_tor=1.2e-3, R_Fix=8.6e-5.

### E. Computational cost

GNEIMO dynamics scales linearly with system size. Adding GNEIMO-Fixman increases
cost by 24% (Fig 10), modest because it reuses GNEIMO by-products; most extra cost is
the recursive Upsilon(i) evaluation per cluster. For all-atom force fields, a
TORSIONAL step costs about 2x a FLEXIBLE step (added constrained-EOM complexity) but
enables 5-20 fs timesteps vs 0.5-2 fs for FLEXIBLE. With the Fixman potential a
GNEIMO-Fixman step costs about 2.24x a FLEXIBLE step - a modest increase overcome by
the larger timestep.

## IV. Conclusions

The Fixman compensating potential is critical for accuracy of constrained MD, but
lack of practical algorithms limited its use. The SOA-based GNEIMO-Fixman algorithm
computes the Fixman potential and torque at linear cost. It was cross-validated
against prior published data for serial systems, and demonstrated to recover the
single and joint torsion-angle pdfs for C4, C5, C11, C15 serial chains and for
realistic branched molecules (alanine dipeptide, valine dipeptide, chignolin). This
is the first use of the Fixman potential and torque for realistic branched molecules
of arbitrary size, establishing viability of routine inclusion of the Fixman
correction in constrained MD to remove statistical biases. Future work: interplay of
the Fixman potential with all-atom force fields, and improved entropy/free-energy
computations using constrained dynamics.
