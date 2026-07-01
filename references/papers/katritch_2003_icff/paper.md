# ICFF: A New Method to Incorporate Implicit Flexibility into an Internal Coordinate Force Field

Vsevolod Katritch, Maxim Totrov, Ruben Abagyan. J Comput Chem 24: 254-265, 2003.

## Abstract

A method to accurately "project" a Cartesian force field onto an internal
coordinate molecular model with fixed-bond geometry. The algorithm generates
the Internal Coordinate Force Field (ICFF), a close approximation of the
"source" Cartesian force field. ICFF reduces the number of free variables by
at least 10-fold and facilitates fast convergence of geometry optimizations
(useful for flexible-ligand docking and macromolecule conformational modeling).
Although covalent geometry is fixed, implicit flexibility is incorporated into
the force field parameters in two ways. First, an empirical torsion energy term
is formulated as a sixfold Fourier series, with Fourier coefficients computed
from the conformational energy profiles of the fully flexible Cartesian model;
the ICFF torsion parameters represent not only the torsion component of the
source force field but also bond bending, bond stretching, and "1-4" van der
Waals interactions. Second, a soft polynomial repulsion function is used for
"1-5" and "1-6" interactions to mimic the flexibility of the connecting bonds.
A local part of the Cartesian force field is used to automatically generate
fixed covalent geometries compatible with the ICFF energy function. The
implementation uses MMFF94s as the source. Benchmarking reproduces MMFF94s
equilibrium conformational energy differences (RMSD ~ 0.64 kcal) and torsion
energy profiles (RMSD ~ 0.37 kcal). All ICFF parameters (except one scaling
factor in the "1-5,1-6" repulsion term) are derived directly from the source
Cartesian force field and do not depend on any particular molecular set. In
contrast, the rigid geometry model with the MMFF94s energy function yields
biased estimations: RMSD > 1.2 kcal for equilibrium energies and ~3.4 kcal for
torsion profiles.

## Introduction / Motivation

Force field-based molecular mechanics (MM) and molecular dynamics (MD) rely on
the accuracy of the force field to discriminate low- from high-energy
conformations, and on the ability to reach convergence of global optimization
or MD equilibrium in reasonable time, which requires ample conformational-space
sampling.

### Modeling in internal coordinates

The size of conformational space and complexity of energy landscapes make
extensive sampling prohibitive for large systems. Constraining high-frequency
modes reduces the space: switching from Cartesian to an internal coordinate
representation with fixed bond lengths and bond angles cuts free variables by
about sixfold; enforcing planar/tetrahedral covalent geometry, planarity of
double bonds and aromatic rings adds up to a ~10-fold decrease. Activating only
a subset of torsion angles (e.g. side chains + loop backbones for homology
modeling, or binding-site residues + ligand for docking) reduces the space
further. Torsion-angle representation gives smaller dimensionality, faster energy
evaluation, and more reliable local minimizations (larger radii of convergence
than Cartesian). Combined with Monte Carlo Minimization (MCM) or Optimal Bias
MCM, this yields reproducible global convergence for demanding problems (ab
initio folding of a 23-residue beta-beta-alpha peptide in a few hours; flexible
docking in minutes).

### The problem: applying Cartesian force fields to rigid-geometry torsion models

Developing accurate transferable force fields compatible with rigid-geometry is
the major bottleneck. Historically ECEPP/ECEPP-3 were derived in torsion space
and adjusted empirically, but are hard to extend beyond peptides. Using the
torsion component of a Cartesian force field directly for rigid models is
questionable for three reasons:

1. **"1-4" contacts.** Van der Waals repulsion between atoms three bonds apart is
   partially relaxed in the flexible Cartesian model by bond/angle flexibility;
   fixing geometry creates inappropriate high-energy clashes. The torsion energy
   function must compensate.
2. **"1-5"/"1-6" clashes.** Steric clashes between atoms four or five bonds apart
   depend on more than one torsion angle, so they cannot be compensated by the
   (single-variable) torsion function.
3. **Covalent geometry definition.** Rigid-model energy is hypersensitive to bond
   lengths/angles; small variations cause severe vdW clashes, so geometry must be
   derived individually to comply with the force field.

ICFF addresses these by accurately *projecting* a source Cartesian force field
into the torsion-angle energy function: "1-4" vdW is folded into the composite
torsion term (with bond bending, stretching, out-of-plane); "1-5,1-6" vdW is a
soft empirical bonded-flexibility term; covalent geometry comes from local
Cartesian optimization. Parameterization is derived entirely from the source
force field with a single adjustable scaling factor, suggesting high
transferability. Parameters are generated at run time ("on the fly") per torsion
type, avoiding the combinatorial explosion of storing parameters for all
four-atom-type combinations.

## ICFF Formulation

### Torsion energy term

The composite torsion term implicitly accounts for bond bending, bond stretching
and "1-4" van der Waals interactions. It is defined as a sixfold Fourier series
(see equations.md eq:1), where $\theta$ is the torsion angle about a specific
bond. Fourier coefficients per bond are computed by:

1. Copy the two atoms flanking the bond and all their immediate neighbors (4 to 8
   atoms) into a separate object, a *torsion fragment*, preserving original atom
   types and bond topology.
2. Grid-search the torsion angle from 0° to 360° in 12° increments -> 30
   conformations of the torsion fragment.
3. For each constrained torsion value, perform a Cartesian *minimization* of the
   fragment with the *local* part of the source force field (all bonded terms +
   "1-4" vdW; no other vdW, no electrostatics).
4. The best (minimized) energies over the grid form the *energy profile*.
5. Approximate this profile by the Fourier series (eq:1) to yield $\{A_k, B_k\}$.
   A direct Fourier transform usually suffices, but a weighted least-squares fit
   (weight = inverse relative energy, see eq:weight) is more accurate in the
   low-energy part.

Torsion parameters do not depend on the model's initial covalent geometry because
fragment geometries are reoptimized per torsion value; only connection topology
must be correct. Parameters are derived per ICFF *bond type* (atom types +
connection topology of the torsion fragment).

### Empirical term for "1-5" and "1-6" interactions

Van der Waals repulsion between atoms four or five bonds apart is treated with a
special soft empirical term (see equations.md eq:2), applied only for
$R_{IJ} < R_{IJ}^*$. Here $\varepsilon_{IJ}$ is the vdW well depth for atom types
$I,J$, $R_{IJ}^*$ the minimum-distance, $R_{IJ}^0$ the zero-crossing distance, and
$C \in [0,1]$ a scaling factor (the only adjustable parameter) balancing harmonic
and cubic contributions. The term is built to satisfy: (1) $E = \varepsilon_{IJ}$
with zero derivative at $R_{IJ}=R_{IJ}^*$; (2) $E$ crosses zero at $R_{IJ}=R_{IJ}^0$,
the same point as the original vdW function. It forms a continuous curve with the
unchanged attractive branch ($R_{IJ}>R_{IJ}^*$), behaving as a soft cubic
polynomial to mimic bond flexibility. Hydrogen-bond atom pairs keep the original
"hard" vdW term (to preserve donor-acceptor attraction). The soft term is finite
at $R_{IJ}=0$, so an appropriately capped electrostatic term is required to
prevent highly charged atoms from collapsing in MC/MD.

### Covalent geometry

Fixed covalent geometry for an ICFF model is obtained by Cartesian optimization
with the local part of the source potential (bonded + "1-4" vdW). Without
nonbonded interactions the whole-molecule optimization decouples into many nearly
independent small problems. A unique initial guess is built by setting torsion
angles to the best-energy values from the torsion profiles; this plus the
simplified energy function guarantees fast convergence and a unique covalent
geometry. Such local-energy geometries also give better average ICFF accuracy
than geometries from the full Cartesian force field or MP2 QM. The procedure can
build an ICFF "residue library" for peptides/nucleic acids/biopolymers by storing
optimal residue covalent geometries.

### MMFF94s as source Cartesian force field

The full MMFF94 potential (see equations.md eq:3) sums bond stretching, bond
bending, stretch-bend, out-of-plane bending, torsion, van der Waals, and
electrostatic terms; detailed term definitions are in eqs. (1)-(14) of ref. 28
(Halgren 1996). MMFF94s (a variant for MM applications, ref. 34) is used.
Reasons for choosing MMFF94: general-purpose self-consistent parameterization for
diverse organics; fully documented functional form/parameters/atom-typing enabling
exact implementation; a published diverse-molecule benchmark (ref. 33) available
electronically as a reference.

## Computational Methods

Development used the ICM molecular modeling package (Molsoft). MMFF94s energy,
atom typing, and parameterization were implemented per published guidelines and
verified against the CCL validation suite. Cartesian optimizations use ICM's
quasi-Newton method with analytic first derivatives. To constrain a torsion angle
in profile calculations, a harmonic restraint $E_r = C_r(\theta-\theta^0)^2$ with
$C_r = 10000$ kcal was applied (deviation typically < 0.1°).

The ICFF module in ICM has four parts: (1) assign an ICFF type to a torsion
fragment from MMFF94 atom types + topology, reusing precomputed params if
available; (2) compute torsion params $\{A_k,B_k\}$ for new types; (3) generate
covalent geometry by local-MMFF94s Cartesian optimization, then build a molecular
tree fixing bond lengths, bond angles, phase angles, leaving torsions free; (4)
compute ICFF torsion, "1-5,1-6", vdW, and electrostatic terms with first
derivatives in internal coordinates, usable in ICM local or MC global
minimization.

## Results and Discussion

Accuracy of the Cartesian-to-torsion projection was assessed by comparing ICFF
energies against MMFF94s (original Cartesian and an artificial fixed-geometry
version). Two benchmarks: equilibrium conformational energy differences, and
detailed torsion energy profiles. The test set is Halgren's diverse organic set
(refs 33, 66); ring-conformation pairs are excluded (rings need a different
algorithm). Numeric results are in checks.md (Tables 1-4).

### Equilibrium energy comparisons

Fixed covalent geometries were defined by local MMFF94s Cartesian optimization.
Different equilibrium conformers were produced by setting torsions to published
MP2/6-31G* geometry values then optimizing torsions with ICFF (MP2 conformations
were not used directly as starting points because of slightly different covalent
geometry). Optimizations were done first with mild angular constraints (to keep
the molecule in a given minimum) then without.

ICFF reproduces equilibrium MMFF94s geometries with mean distance RMSD ~0.04 Å.
Maximum RMSD 0.27 Å occurs in molecules with tri-coordinate nitrogens (ICFF forces
planar configuration). ICFF reproduces relative equilibrium energies much better
than rigid MMFF94s (RMSD ~0.64 vs ~1.21 kcal). Excluding eight high-energy points
($\Delta E_{\text{MMFF94s}} > 4$ kcal) gives ICFF RMSD 0.46 kcal.
RMSD(ICFF vs QM) = 0.74 kcal is close to the quadrature of RMSD(ICFF vs MMFF94s)
= 0.64 and RMSD(MMFF94s vs QM) = 0.31, suggesting independent error contributions.

### Derivation (not implemented): added-error estimate

The extremely low MMFF94s RMSD = 0.31 kcal vs QM is specific to this set (used to
parameterize MMFF94). A fairer unbiased measure of Cartesian-FF accuracy for this
set is the second-best RMSD, 1.03 kcal (MM3). Assuming ICFF projection adds an
independent RMSD contribution of 0.64 kcal on top of a 1.03-kcal Cartesian FF, the
combined RMSD is $\sqrt{1.03^2 + 0.64^2} \approx 1.21$, i.e. ICFF projection adds
only about 18% ($1.21/1.03 \approx 1.18$) to the overall energy RMSD.

### Torsion profile comparisons

All 110 dihedral angles responsible for conformational changes were profiled.

**Accuracy of the ICFF torsion term.** "Local energy" excludes electrostatics, vdW
beyond "1-4", and interactions of atoms four+ bonds apart. ICFF local energy is the
torsion term alone; MMFF94s local energy = torsion + bend + stretch + OOP + "1-4"
vdW. In the moderate-energy range ($E_{\text{MMFF94s}} < 5$ kcal) ICFF RMSD is
0.26 kcal vs 0.58 kcal for the rigid MMFF94s model. Even in simple cases the rigid
Cartesian model can fail (e.g. loses the C5 minimum of the glycine dipeptide analog
that bond flexibility - explicit in Cartesian, implicit in ICFF - preserves).
Residual ICFF errors arise because the torsion function is one-dimensional while
interactions of fragment atoms with outside atoms depend on two+ torsions (up to
~1 kcal in rare highly-flexible clashed cases). RMSD is computed only where
$E_{\text{MMFF94s}} < E_{\text{cut}}$ (e.g. 5 kcal), since both ICFF and the
reference lose accuracy at high strain.

**Fourier series approximation.** A sixfold Fourier series reproduces torsion
profiles to ~0.01 kcal (fivefold ~0.08, fourfold ~0.2). Small coefficients (< 0.01)
can be nullified. Weighted least-squares fitting (weight ∝ inverse relative energy,
eq:weight) improves fit RMSD to ~0.001 kcal below the 5-kcal cutoff and reduces the
effect of high-energy deviations.

**Accuracy of the full ICFF energy function.** "1-5,1-6" interactions cannot be
folded into the torsion term and must be computed explicitly, but the original
MMFF94 vdW term is inappropriate in fixed geometry: e.g. a 1-6 O-H clash in the
C=C-O-C trans vinyl formate reaches 24 kcal with hard vdW but stays under 1 kcal
when bonded geometry is relaxed. The soft polynomial term (eq:2) brings the energy
close to the flexible-model value. Coulomb electrostatics between 1-5/1-6 neighbors
does not introduce substantial rigid-model error and needs no special treatment.
The blend factor $C$ was fit to minimize RMSD at $E_{\text{cut}} = 3, 5, 10$ kcal
and for the equilibrium set; optimal RMSD occurs for $0.50 < C < 0.60$ with weak
dependence, so consensus $C = 0.55$ was chosen (RMSD = 0.37 kcal at 5-kcal cutoff).

Interactions of atoms six+ bonds apart keep the original MMFF94 vdW term, assuming
enough torsional freedom to avoid steric hindrance (valid for the test set; in rare
stiff conjugated ring systems the soft potential would be more appropriate).

**Effect of covalent geometry.** ICFF accuracy is relatively insensitive to how
covalent geometry is generated (local MMFF94s, full MMFF94s, or MP2 QM): energy-RMSD
differences stay under 25% across cutoffs. Results favor "local energy" Cartesian
optimization (bonded + "1-4" vdW only): omitting electrostatic/long-range vdW makes
geometry more "generic" (less conformationally biased) and gives fast reliable
convergence with a deterministic unique result.

### ICFF performance

The projection is inexpensive and automated; parameters generate "on the fly",
precalculating fixed geometry + torsion params for an average drug-like compound on
a ~0.1 s scale. The same procedure can build geometry/torsion libraries for proteins
and biopolymers. ICFF is ~50% faster than ECEPP, due to implicit "1-4" vdW treatment
and the simplified "1-5,1-6" repulsion (fewer pairwise interactions).

## Conclusions

ICFF derives a conformational energy function for a fixed-covalent-geometry model
directly from a source Cartesian force field. Using MMFF94s as source, the ICFF
projection into torsion coordinates mimics the Cartesian force field with RMSD
~0.3-0.5 kcal, limited by the source quality. "Implicit flexibility" built into the
internal-coordinate parameters is critical: torsion-profile energy RMSD drops from
3.5 kcal (rigid) to 0.37 kcal (ICFF). Flexibility enters two ways: (1) the composite
torsion term absorbs bond stretching, bending, out-of-plane, and "1-4" vdW; (2) a
soft "1-5,1-6" repulsion mimics connecting-bond flexibility. Fixed covalent geometry
is obtained by fast deterministic minimization of the "local" part of the source
Cartesian force field.
