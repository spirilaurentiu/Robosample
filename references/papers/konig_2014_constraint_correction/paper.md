# Correcting for the free energy costs of bond or angle constraints in molecular dynamics simulations

Gerhard Konig, Bernard R. Brooks. Laboratory of Computational Biology, NHLBI, NIH.
Biochim. Biophys. Acta (2014). DOI: 10.1016/j.bbagen.2014.09.001.

## Abstract

Free energy simulations allow calculation of thermodynamic properties of binding or
enzymatic reactions. This paper introduces methods to increase the accuracy and
precision of free energy calculations by calculating the free energy **costs of
constraints during post-processing**. The primary purpose of employing constraints
is to increase phase-space overlap between ensembles (required for accuracy and
convergence).

Method: the free energy costs of applying/removing constraints are computed as
additional explicit steps in the free energy cycle. The techniques focus on **hard
degrees of freedom** and use both gradients and Hessian estimation. Enthalpy,
vibrational entropy, and Jacobian free energy terms are considered.

Results: demonstrated on harmonic/anharmonic oscillators, four-atomic benchmark
systems, ethane->methanol alchemical mutation, and alanine<->serine free energy
simulations. Analytical-case errors all below 0.0007 kcal/mol; ethane->methanol
accuracy improved from 0.15 to 0.04 kcal/mol. Unconstrained alanine->serine overlaps
0.15-0.9%; constraints raise overlap to 2.05% (+94% average), doubling precision. The
approach reduces constraint-induced errors by ~an order of magnitude.

General significance: primary utility is free energies for systems with disparate
energy surfaces and bonded terms, especially multiscale MM/QM simulations.

## 1. Introduction

Constraints are used in most MD simulations because the maximum time step is
restricted by the fastest motions; freezing rapid vibrational modes lets one use
longer time steps without losing energy conservation. Bond constraints can reduce
required computer time by a factor of three. Constraints also improve convergence and
efficiency of free energy simulations (e.g. Simplified Confinement, Confine-and-Release).

Available constraint algorithms include SHAKE (two-body, iterative Gauss-Seidel;
recursive nature makes it impractical for many coupled bonds/rings), RATTLE (adds
velocities), LINCS (linear solver), and the Shape rigid-body integrator (SHAKE
accuracy with more/higher-order constraints - rigid bodies of >=3 centers, i.e. bonds,
angles, dihedrals).

Imposing constraints restricts phase space: constraining both bonds and angles can
quench dihedral-angle transitions (van Gunsteren), lower efficiency, and change
trans-gauche transition rates (Toxvaerd on decane). Constraints shift normal-mode
frequencies of biomolecules (between 100 and 1400 cm^-1). Deviations of 0.2-0.5
kcal/mol appear in solvation free energies when bond-length constraints are employed -
motivating explicit accounting for constraint free energy changes.

Existing constraint corrections were mostly developed for Thermodynamic Integration
and cannot be used straightforwardly with the increasingly popular Bennett's
acceptance ratio (BAR), multistate-BAR (MBAR), or Non-Boltzmann Bennett (NBB) methods.
For BAR, the utility of constraints is to **increase phase-space overlap** between two
distinct energy surfaces; since the variance of the free energy estimate is directly
linked to overlap, constraints yield more accurate, faster-converging results. This is
particularly useful for coupling disparate surfaces (QM/MM with MM).

Approach: correct for constraints by post-processing the trajectories - explicitly
calculating the free energy costs of adding/removing constraints as additional steps
in the free energy cycle, via a gradient calculation combined with normal mode
analysis to approximate the constrained-DOF partition-function contributions.
Compatible with most quantum packages and QM/MM.

Inspiration (Go and Scheraga): bond lengths and angles can be treated as functions of
the dihedral angles. **Hard** normal modes (high force constant) are insensitive to
conformation and can be regarded as functions of the **soft** variables; they can be
frozen (constrained) during dynamics, while soft modes are treated classically. A
post-processing step then changes bond lengths/angles according to their soft
environment and accounts for the associated free energy decrease. Valid only near the
equilibrium conformation.

## 2. Theory

Start with a classical harmonic oscillator with internal coordinate q and potential
energy eq:1. The partition function (eq:2) integrates over the single DOF; absolute
free energy follows via the Gaussian integral (eqs 3-4).

Introducing a constraint removes DOF q entirely: the integral collapses to a Dirac
delta point (using ∫δ(x)dx = 1), giving the constrained partition function (eq:5) and
free energy (eq:6). The free energy of imposing a constraint (eq:7) splits into an
enthalpic contribution ΔH (eq:8, temperature-independent, note U(Δq_cons) inherently
contains U_0) and an entropic vibrational contribution ΔG_harm (eq:9). Total: eq:10.

### 2.1. Implementation using normal mode analysis

The isolated 1-D analysis fails once non-bonded interactions are present (Coulomb, LJ,
water box) - the system becomes an anharmonic oscillator with unknown potential.
Starting from an arbitrary constrained structure Δq_cons, expand the potential via
Taylor series (eq:11). At room temperature the harmonic term dominates (bond will not
break; non-bonded interactions mostly shift the energy minimum, not the frequency), so
truncate at second order (eq:12), ignoring anharmonicity.

To find the unconstrained energy minimum, use a one-step Newton-Raphson with the
gradient and Hessian at the constrained structure; since the gradient is zero at the
minimum, the displacement is eq:13. This yields the enthalpic contribution eq:14 and
harmonic entropy eq:15. If Δq_cons is large, re-evaluate U'' at the minimum.

For **multiple simultaneous constraints**, use a reduced-basis harmonic analysis:
partition the full 3n x 3n Hessian into relevant (constrained) and irrelevant parts,
approximately block-diagonalizing it (valid because soft->hard coupling is small and
most hard terms are uncoupled, except the Urey-Bradley term in bond angles). The m
constrained DOFs form an orthonormal sub-basis C of m mass-weighted, normalized
Cartesian-displacement vectors c.

- Constrained bond r_i: basis vector eq:cbasis-bond (vector from atom 1 to atom 2 and
  vice versa).
- Constrained angle theta_i (atoms j,k,l; constraint on atom l): basis vector
  eq:cbasis-angle (tangential, normal to bond r_kl in the j-k-l plane). A displacement
  along it changes bond length r_kl, corrected by including angle curvature via
  eq:angle-curvature (implemented in the new CANG modes).

The reduced Hessian is eq:16 and reduced gradient eq:17. Because of coupling
(off-diagonal H elements), an eigendecomposition eq:18 diagonalizes to Λ (force
constants/frequencies of the constrained DOFs). CHARMM offers two routes: **REDUce**
(VIBRAN) uses the full 3N x m matrix C (keeps coupling); **RAYLeigh** analyzes each
basis vector separately (neglects intra-constraint coupling). Enthalpy and entropy
corrections in matrix form: eq:19 and eq:20 (‖...‖_1 = sum over all elements; log
applied element-wise since Λ is diagonal).

### 2.2. Jacobian factors

The correction so far applies to internal coordinates, but MD runs in Cartesian
coordinates. The internal->Cartesian conversion adds an entropic Jacobian term (longer
bonds and obtuse angles have more available phase space). Only DOFs that change
(affected by constraints) need it. Cartesian-space release free energy: eq:21; general
Jacobian change eq:22.

Analytic Jacobian factors (rigid-rotor analysis, Herschbach et al.): bond eq:23
(J_r = r_ik^2); linear-chain angle eq:24 (J_theta = sin theta); branching angle eq:25
(J_theta' = 1/sin theta'). The current code only handles angle constraints of linear
chains - not a major restriction since SHAKE is usually not applied to angles of
branched structures (poor convergence).

### 2.3. Application to free energy simulations

Running a free energy simulation with constraints explores only an infinitesimally
slim slice of the constrained-DOFs' phase space. Accounting for the cost of releasing
constraints adds those extra dimensions back (analogous to Straatsma-McCammon
treatment of rotational isomers).

Factor the unconstrained partition function into constrained-DOF and released-DOF parts
(eq:26), replace the inner integral with the per-configuration constraint free energy
(eq:27), and express the release free energy as a partition-function ratio (eq:28).
This reduces to the key working formula eq:29: an exponential (Zwanzig / thermodynamic
perturbation) average of the per-frame constraint free energy over the constrained
ensemble - implemented in most major simulation packages.

Thus the harmonic analysis performs an a-posteriori sampling of the constrained DOFs.
If the target potential is harmonic, it performs *perfect* sampling of that DOF given a
set of soft coordinates; the only errors arise from anharmonicity and incomplete
sampling of soft DOFs.

## 3. Methods

### 3.1. (An-)harmonic oscillators
H atom bonded to a fixed non-interacting atom (equilibrium bond length). Three cases:
gas phase; with a fixed Na+ at 2.5 A; with a fixed Cl- at 2.5 A (b,c are anharmonic).
Five CGenFF hydrogen types: HGA3 (ethane), HGA5 (ethene), HGR61 (benzene), HGP1
(methanol OH), HGP3 (ethanethiol SH). Absolute free energies via numerical integration
of the partition function (Mathematica NIntegrate, AccuracyGoal=Infinity), integrating
within +/-1.9 A of equilibrium. Gradients/Hessians by numerical differentiation (ND,
scale 0.0001 A) at the minimum (FindMinimum).

### 3.2. Four-atomic benchmark systems
Linear four-atom system, no non-bonded interactions -> analytic reference via
rigid-rotor. Init: all bonds 2 A (force const 200 kcal/(mol A^2)); angles 110 deg
(force const 50 kcal/(mol rad^2)); dihedral force const 1 kcal/mol, multiplicity 3.
Corrections via the new OPTI command in CHARMM VIBRAN.

### 3.3. Water boxes
TIP3P water: gas-phase pentamer (global minimum structure) plus boxes of 395, 787,
1636, 3290 molecules (cubic, sides 24.01/30.25/38.11/48.02 A, CHARMM-GUI). No Ewald
(VIBRAN limitation); cutoff 10-12 A, force-shifting electrostatics, VdW shifting. SHAKE
imposed at parameter equilibrium bond length; minimized (100 steepest descent + 2000
ABNR); SHAKE released; correction via CBND and CANG modes; re-minimized to true
minimum.

### 3.3.1. Ethane-Methanol
Solvation free energy difference. CHARMM22, dual-topology hybrid via MSCALE (bonds/
angles in main process; dihedral + non-bonded in ethane and methanol subprocesses,
mixed by lambda). Gas: Langevin (friction 5 ps^-1, 300 K, cutoff 998 A). Solution: 862
TIP3P, octahedral box 32.168 A, Nose-Hoover 300 K, LJ switched 10-12 A, PME, 1 fs step;
runs with and without SHAKE on all hydrogens. BAR from 5 ns gas / 1 ns solution; 5
lambda (gas), 11 lambda (solution); 4 repeats.

### 3.3.2. Alanine-Serine
N-acetyl-methylamide alanine and serine in water. AMBER Cornell and CHARMM22, with and
without SHAKE. 243 TIP3P, truncated octahedron (cube side 21.4 A), Nose-Hoover 300 K,
LJ switched 9-10 A, PME. Trajectories every 10 steps; 4 repeats. BAR energies via
BLOCK EAVG; constraint corrections via VIBRAN; analyzed with the Zwanzig equation.

## 4. Results and discussion

### 4.1. Analytical results for (an-)harmonic oscillators
See checks.md Table 1. In gas phase (constraint at equilibrium) the only contribution
is vibrational entropy ΔG_harm (directly related to force constant): weakest for HGP3
(thiol, -1.481), strongest confinement for HGP1 (hydroxyl, -1.684). Adding Na+ shifts
the minimum -> nonzero ΔH (up to 0.211) and Jacobian (below -0.005); ΔG_harm dominates.

Neglecting ΔH and ΔG_Jacobian is tempting (few percent of total) but wrong in a
thermodynamic cycle: ΔG_harm largely cancels between states (differences <0.003
kcal/mol), leaving ΔH (0.017-0.211) and Jacobian (-0.002 to -0.005) as the relevant
non-canceling terms. Cl- causes steric clash -> larger ΔH (0.052-0.383); polar HGP1/
HGP3 (small VdW radii, large positive charge) are attracted, stretching the bond ->
positive Jacobian. Signs: ΔH >= 0 always; ΔG_harm always negative; ΔG_Jacobian either
sign.

Anharmonicity is the main approximation error. Gas phase has none. With non-bonded
forces the worst case (HGR61 + Cl-) gives 0.0007 kcal/mol - three orders of magnitude
below ΔH, one below the Jacobian. Anharmonic error can be either sign (HGP3+Cl-
negative), allowing cancellation. Not a limiting factor for SHAKE-on-hydrogens at room
temperature.

### 4.1.1. Four-atomic benchmark systems
See checks.md Table 2. Uncorrected errors 2-476 kcal/mol. Standard normal modes (NMD,
no Jacobian/curvature) reduce to 6.39E-05 to 85 kcal/mol; residual from missing
Jacobian (bonds) and bond-length distortion from using only the tangential vector for
large angle changes (up to 24 kcal/mol remaining). Dedicated CBND+CANG modes (C-NMD)
reach machine precision (3.4E-16 to 4.1E-19) because there is no anharmonicity. **Must
correct bonds (CBND) before angles (CANG)** so bond vectors are current; angle
correction uses trigonometric functions (quadrant care needed).

### 4.2. Impact of soft degrees of freedom
The correction applies only to hard DOFs; reaching the local minimum usually also needs
soft-DOF (dihedral/translation/rotation) rearrangement, which cannot be corrected
(non-harmonic, external-force-dependent; changing them would bias Boltzmann weights).
Residual error is therefore unavoidable. Tested on water boxes (perturb bonds 0.05 A,
angles 1 deg; see checks.md Table 3): uncorrected ~0.6 kcal/mol per water; ~90% error
reduction consistently; residual ~0.05 kcal/mol per water. Hessian cost scales
cubically; fast path (Hessian = force-field force constant) does 3290 waters in 0.56 s
vs 16494.8 s, only ~30% over a single energy call. Embarrassingly parallel (each frame
independent; RAYLeigh makes modes independent too).

### 4.3. Application to free energy calculations (ethane->methanol)
See checks.md Table 4. Unconstrained deviation from experiment 0.04 kcal/mol; SHAKE
raises error to ~0.19 (most, 0.18, on methanol side - missing hydroxyl-bond stretching
to accommodate H-bonds with water). Correction restores deviation to 0.04 kcal/mol
(difference from unconstrained not statistically significant).

### 4.4. Free energy difference between alanine and serine
See checks.md. Alanine->serine is among the largest single-step alchemical mutations;
unconstrained overlap only 0.15% (critically low). Constraints raise overlap to 2.05%
(+94% avg). Constraint-correction std devs <0.003 kcal/mol (fast convergence despite
exponential averaging). Cube cycle-closing error 0.22 kcal/mol total (<0.32 propagated
-> not significant); constrained trajectories close to 0.02 vs 0.22 unconstrained.
Constrained std devs ~33% lower on average; precision effectively doubled. Cost: full
Hessian ~14x a plain energy eval; force-constant approximation ~40% slower. Constraints
enable ~3x larger time step; with 2 fs, ~8x cheaper than equivalent unconstrained sim.

## 5. Conclusions

Calculates the free energy of releasing constraints during post-processing via harmonic
analysis. Three assumptions: (a) constrained DOFs approximable by a harmonic potential
(true for bonds/angles in MM); (b) constrained structure close to the energy minimum
(select constraints near equilibrium bond length; polar groups may need bond-stretch
accounting); (c) coupling of constrained to unconstrained DOFs is small (angle
constraints must be obtuse enough and frequency-separated from dihedrals to avoid
altering dihedral transition rates).

Validated to 0.0007 kcal/mol (anharmonic oscillators) and machine precision (harmonic
four-atomic). Reduces water-box deviation by >90%; Hessian cost is the main limitation
(cubic scaling) but is mitigable via force-constant approximation or subsampling
frames, and is embarrassingly parallel. Improves ethane->methanol accuracy from 0.15 to
0.04 kcal/mol and doubles alanine<->serine phase-space overlap and precision.
Particularly useful for multiscale (QM/MM) free energy calculations, allowing larger QM
regions and larger MM/QM energy-surface differences, fewer lambda steps, and shorter
simulations.

### Relevance to Robosample (implementation note)
This paper's machinery maps directly onto rigid-body / constrained internal-coordinate
sampling: when a molecule is represented as rigid bodies with frozen hard DOFs
(bonds/angles), the Boltzmann-correct free energy relative to a fully flexible model
requires exactly the ΔH + ΔG_harm + ΔG_Jacobian correction (eqs 19-25) evaluated
per-frame and combined via the Zwanzig average eq:29. The Jacobian factors (eqs 23-25)
are the internal<->Cartesian coordinate-transformation terms that also arise in
constrained-ensemble sampling.
