# MMB-GUI: a fast morphing method demonstrates a possible ribosomal tRNA translocation trajectory

Alex Tek, Andrei A. Korostelev, Samuel Coulbourn Flores. Nucleic Acids Research (2016), DOI 10.1093/nar/gkv1457.

> Note: This is primarily a software/application paper describing a Chimera-based GUI for MMB
> (MacroMoleculeBuilder) and its internal-coordinate multiscale morphing method. Its implementable
> content is limited to one benchmark metric (the morphing *improvement* score) plus concrete
> benchmark numbers and simulation parameters (see `equations.md` and `checks.md`). Its main value to
> an implementer is routing: it points at the MMB internal-coordinate engine (Flores 2011, built on
> Simbody), the flexibility/physics-zone multiscale scheme, and the Weiss & Levitt morphing benchmark.

## Abstract

Easy-to-use macromolecular viewers, such as UCSF Chimera, are a standard tool in structural biology.
They allow rendering and geometric operations on large complexes (viruses, ribosomes). Dynamical
simulation codes enable modeling of conformational changes but may require considerable time and many
CPUs. MMB-GUI targets the middle ground: visualization combined with quick, interactive modeling of
conformational changes, even of large complexes. MMB uses an internal-coordinate, multiscale approach,
yielding as much as a 2000-fold speedup over conventional simulation methods. Chimera is used as an
interactive graphical interface to control MMB. The authors demonstrate morphing of heterogeneous
macromolecules (varying biopolymer type, sequence, chain count), recapitulating structural
intermediates, and build a possible trajectory of EF-G mediated tRNA translocation in the ribosome
(~150 000 atoms) with all-atom structures.

## Materials and Methods

### The GUI (Chimera platform) controls MMB features

MMB is a multiscale, internal-coordinate macromolecular modeling code. It lets the user fully control
system flexibility at the level of domains, residues or individual bonds. It supports a wide array of
constraints and forces: base-pairing forces (RNA folding), springs (homology modeling), density-based
forces (fitting to density maps), and a conventional MD force field. Multiscale processing is key to
speed on large complexes: rigid regions become kinematically a single body, giving up to 2000-fold
savings in compute time. *Flexibility zones* (hinges, interfaces, mutation/active sites) can be created,
while the MD force field is limited to *physics zones* surrounding these regions.

Three steps to run a simulation: (1) initialize MMB with biopolymers (RNA, DNA, protein chains), each
given a type, chain ID, sequence and structure (from a form, command line, input file, PDB, or Chimera
model); (2) set flexibility, constraints, forces, and physics zones via widgets; (3) run the simulation
(start/stop/restart, review trajectory).

### Morphing

A morph is an interpolated trajectory connecting two or more conformations of a macromolecule. Morphing
generates a single, clash-free, physically rational, directed trajectory connecting an initial (e.g.
open) to a final (e.g. closed) conformation. It does not compute thermodynamic quantities (unlike biased
MD: Steered MD, Umbrella Sampling, PMF). RMSD of the trajectory from an experimentally observed
intermediate structure, at closest approach, is the objective quality measure.

Procedure:
1. Rigidly align the initial and final structures.
2. MMB computes a gapped alignment (SeqAn) between user-specified biopolymers in the initial and final
   structures and connects corresponding atoms with springs (same springs used for homology modeling and
   rigid alignment).
3. The user defines flexible segments in the initial structure while keeping the final structure rigid.
   Initial and final structures may comprise protein, DNA, RNA, water, ions, small molecules, and may
   have different chain counts. Collision-detecting spheres can be placed at hinges/interfaces; physics
   zones with an active MD force field can be specified.
4. Driven by the springs, the semiflexible initial structure is aligned with the final structure. The
   process stops after a specified simulation time, upon a convergence criterion, or by user command.

Hinges are often annotated in the literature or predicted using web servers such as HingeMaster.

### Recapitulating intermediate structures with morphing (benchmark)

Tested against the Weiss & Levitt benchmark of morphing methods. Five initial structures were submitted
to HingeMaster to predict hinge residues; results were manually filtered to define up to three flexible
groups of residues. Morph simulations applied a physics zone of 10 Å around all flexible residues. Each
pair of starting/target structures was rigidly aligned using the Chimera *match* command.

The *improvement* score (defined in Weiss & Levitt) measures how much better an interpolated structure is
as an approximation of the known intermediate, relative to the initial and final structures — see
`equations.md`. All RMSD values are based on Cα atoms and computed with Chimera. sRMSD (aligning on one
of two domains, then computing RMSD on the other) is also reported to emphasize large-scale rearrangement.

MMB gave the highest improvement score in three of five cases and was close to the best method in the
others. It produces smooth trajectories and reaches convergence in a few minutes on a single laptop CPU
core (Macbook Air, Intel Core i5 at 2.8 GHz, OS X 10.8.5).

### Trajectory of tRNA translocation induced by EF-G

To demonstrate applicability on a large system, an atomistic trajectory of tRNA translocation in the 70S
bacterial ribosome, catalyzed by elongation factor G (EF-G), was generated. The system comprises more
than 150 000 non-hydrogen atoms. Conformational rearrangements can be modeled with surprisingly few
flexible residues; rearrangements are largely describable as domain motions.

During protein synthesis, tRNAs translocate from A (aminoacyl) to P (peptidyl) to E (exit) site,
traversing > 100 Å within the ribosome. Translocation proceeds in two global steps: (1) displacement of
the acceptor ends of the A- and P-site tRNAs to the P and E sites of the large 50S subunit, while the
anticodon stem loops (ASL) remain in the A and P sites on the small 30S subunit (A/P and P/E hybrid
states), coupled with intersubunit rotation (small subunit rotates clockwise up to 12° relative to the
large subunit); (2) during the reverse rotation, EF-G triggers concerted movement of mRNA and tRNAs from
A/P and P/E into the classical P and E sites.

Morphing was calculated between six experimentally determined structures (states 1–6) representing
successive translocation states. The 3.3 Å crystal structure from *T. thermophilus* (PDB 2WDG, 2WDI) was
used as the base structure. Each subsequent morph used the resulting structure of the previous morph as
its starting point.

In a first morph all ribosomal proteins other than EF-G were excluded. Flexibility was allowed at the
base of the neck and beak of 16S rRNA, the base of the L1 stalk and the A-site finger (H38) on 23S rRNA,
and at the base of the anticodon stem-loops of the tRNAs. Collision-detecting spheres around these zones
prevent clashes. mRNA was made entirely flexible with imposed Watson–Crick interactions between codon
bases and their anticodon bases on the tRNAs. EF-G was made flexible between domains II–III and III–IV.
EF-G was placed at a distance from the ribosome after steps 2 and 5 to mimic binding (step 2→3) and
release (step 5→6).

For each morph, the final structure was rigidly aligned to the initial based on the 16S RNA
(Chimera MatchMaker). Each chain of the initial structure was morphed to its correspondent on the final
via MMB's *gappedThreading* command (auto springs from a gapped alignment); only mRNA used the *threading*
command (manual alignment) to ensure correct translocation of codon bases. Simulations stopped when the
energy difference between two consecutive frames was below 50 kJ/mol during five frames. Each morph
included between 144 791 and 155 378 moving atoms and took 20 to 43 min to converge on a laptop.

A full-ribosome trajectory (with most proteins available in crystal structures of states 4 and 5) was
also generated (249 313 atoms; ~45 min). Most proteins were kept rigid and welded to the rRNA of their
subunit; only protein S7 (near the E site) was made partly flexible. The threading force constant *F* for
the P-site tRNA was increased from 30 to 60 to overcome transient contacts with S7; with *F* = 30 the
tRNA became stuck in the gate.

Finally, to test predictive potential, a morph from state 1 directly to state 3 (omitting state 2)
partially recapitulated state 2, despite state 2 lying off the linear interpolation pathway.

## Results and Discussion

MMB-GUI adds MMB's modeling capabilities (flexible fitting, homology modeling, ΔΔG prediction of
protein–protein binding, macromolecular folding) to Chimera through a visual, interactive interface. It
can economically morph large macromolecular complexes heterogeneous in sequence and chain count. The
morph method accurately recapitulates known intermediates for structurally distinct proteins (RNase III,
ribose-binding protein, myosin, Ca²⁺-ATPase, 5′-nucleotidase). The best-performing cases (5′-nucleotidase,
myosin, RBP) are domain-hinge-bending proteins with clear domain boundaries and few intradomain
rearrangements; RBP appears in the Hinge Atlas Gold dataset.

Ca²⁺-ATPase is a harder case: under a uniform domain-hinge-bending protocol (HingeMaster hinges) the
minimum RMSD vs final was 5.13 Å. Ca²⁺-ATPase has multiple domains with unclear boundaries, intradomain
rearrangements, and secondary-structure changes; it does not meet the definition of domain hinge bending
and is less amenable to multiscale treatment. Flexibilizing residues 42–47, 57–59, 80–84, 112–114, and
122–126 decreased RMSD vs final to 4.15 Å, but RMSD vs intermediate increased to 9.09 Å.

For tRNA translocation, superposition of tRNA intermediates suggests A- and P-tRNAs undergo similar
motions during translocation to the P and E sites. As the A-tRNA moves to the P site, its CCA 3′-end and
elbow travel up to 10 and 40 Å respectively; the P-tRNA into the E site requires the CCA end and elbow to
traverse ~40 and 60 Å. Residue 34 at the tip of the anticodon loop sweeps ~14 Å during A→P translocation
and ~17 Å during P→E translocation. The ASL passes through a narrow channel (gate) between the body and
head of the small subunit (16S rRNA residues 789–791 and 1338–1342). The channel (~13 Å wide in
non-rotated/partially rotated ribosomes) widens to more than 20 Å due to head swiveling. Morphing
suggests movement through the gate is coupled with 'untwisting' of the ASL relative to the D stem; the
ASL is rotated by nearly 30° relative to the pE/E state.

With all proteins included, complete P→E translocation initially failed (ASL stuck at the gate) because
S7 (near the E site) transiently repels the ASL. Moderately increasing the threading force (F: 30→60)
achieved translocation; a β-hairpin of S7 (aa 77–84) engages the translocating tRNA after gate passage.

Morphing is unlikely to have the predictive power of MD (flexibility and physics are limited), but the
economical GUI-controlled morphing may yield functionally relevant hypotheses meriting validation with
MD and experiment.

## Availability

- Input files, final trajectory, and video: https://simtk.org/home/efgtranslocat
- MMB 2.15 documentation, source, and binaries (OSX, Linux): https://simtk.org/home/rnatoolbox
