# Protein Structure Refinement of CASP Target Proteins Using GNEIMO Torsional Dynamics Method

**Authors:** Adrien B. Larsen, Jeffrey R. Wagner, Abhinandan Jain, Nagarajan Vaidehi
**Venue:** J. Chem. Inf. Model. 2014, DOI 10.1021/ci400484c

> Routing / application paper. No implementable equations of its own; it is an
> application of the GNEIMO internal-coordinate (torsional) MD method combined
> with temperature replica exchange (REXMD) to protein-structure refinement.
> Method-defining equations live in the GNEIMO source papers (see `depends_on`).
> This file keeps the protocol, parameters and qualitative findings; concrete
> benchmark numbers are in `checks.md`.

## Abstract

A longstanding challenge in computational protein-structure prediction is refining
low-resolution comparative-modeling structures into accurate atomistic models.
The internal-coordinate MD technique GNEIMO (generalized Newton-Euler inverse mass
operator) freezes high-frequency degrees of freedom and models the protein as a
collection of rigid clusters connected by torsional hinges. This allows larger
integration time steps and focuses the conformational search on the low-frequency
torsional degrees of freedom. GNEIMO with temperature replica exchange was applied
to refine low-resolution models of 30 CASP target proteins. GNEIMO torsional MD
gave refinement of up to 1.3 Angstrom RMSD in coordinates for the 30 CASP targets
**without** using any experimental data as restraints, in contrast to unconstrained
all-atom Cartesian MD which required restraints for refinement under the same
conditions.

## Introduction (method context)

Comparative (homology / template-based) modeling produces protein models that can
deviate significantly from crystal structures, especially in local regions, and
must be refined for functional analysis and drug design. Torsional Monte Carlo
methods are successful for refinement but limited by their energy-driven
conformational search. Force-driven MD can cross energy barriers. All-atom
(Cartesian) MD has shown limited refinement success without knowledge-based
potentials or experimental restraints; restraints to the starting structure give
better refinement than unrestrained Cartesian MD. This work focuses on MD-based
refinement **without** restraints.

GNEIMO is an internal-coordinate MD method. Torsional dynamics is one application:
high-frequency degrees of freedom are held rigid with hard holonomic constraints,
and the protein is a collection of rigid clusters (single atoms up to whole domains,
user defined) connected by flexible hinges of one to six degrees of freedom.
Coupling GNEIMO torsional MD with temperature replica-exchange MD (REXMD) enables
efficient sampling in the low-frequency torsional space.

## Computational Methods

### GNEIMO constrained dynamics (summary)

GNEIMO is a constrained internal-coordinate MD method. The equations of motion in
internal coordinates are coupled; naively solving them scales as the cube of the
number of degrees of freedom, but the GNEIMO recursive algorithm reduces this to
**linear** scaling in the number of degrees of freedom, making torsional MD feasible
for proteins. Hinges can carry 1-6 DOF; clusters range from single atoms to helices
to whole domains. (The recursion and mass-operator algebra are in the GNEIMO source
papers, `jain_1993`, `vaidehi_1996`, `wagner_2013`, `jain_2010`.)

### All-torsion GNEIMO-REXMD refinement protocol

- Force field: **AMBER99SB**.
- Solvation: **Generalized Born / surface area (GB/SA), OBC** implicit model.
  - Solute interior dielectric = **1.5**; solvent exterior dielectric = **78.3**.
  - Solvent probe radius = **1.4 Angstrom** for the nonpolar GB/SA component.
- Nonbonded forces switched off at a cutoff radius of **20 Angstrom**.
- All torsional degrees of freedom active.
- Constant temperature via **Nose-Hoover** thermostat.
- Integrator: **Lobatto**, time step **5 fs** (text also cites stable steps up to 10 fs).
- **REXMD:** 32 replicas spanning **310-415 K**; replica temperature swaps attempted
  via the Metropolis criterion every **5 ps**.
- Total simulation time **15-100 ns per replica** (32 replicas) per target.
- Starting decoys minimized first with all-atom conjugate gradient ("sander",
  AMBER99SB) before GNEIMO-REXMD.

### Targets

30 CASP8 and CASP9 targets: 23 from the structure-refinement category (names begin
"TR", sizes 63-192 residues) and 7 from the structure-prediction category (names
begin "T0"). Refinement-category decoys taken from the CASP website. For prediction
targets, homology models were built with **MODELLER**: templates chosen by PDB
sequence query at 30-80% identity (target and close homologues published after the
CASP competition removed to avoid bias); 100 models per target clustered into 5
groups; best representative by PROCHECK G-factor chosen as the starting decoy, then
minimized (sander, AMBER FF99SB) and refined.

### Assessment metrics

- **RMSD**: root-mean-square deviation of backbone (or C-alpha) coordinates to the
  X-ray/NMR reference, over the combined REXMD trajectory of all replicas. Computed
  with MDAnalysis. Reported both for secondary-structure regions and whole structure.
- **Percent native contacts**: build the N x N C-alpha(i)-C-alpha(j) distance matrix
  (N = number of residues). A native contact = pair with distance < 8 Angstrom and
  more than 4 residues apart in sequence (index 1; else 0). A simulation snapshot pair
  is counted as the same contact if its distance is within 0.5 Angstrom of the native
  contact distance. Percent native contacts = identical contacts / total native contacts.
- **GDT_TS**: average number of aligned C-alpha atoms fitting under distance-to-native
  cutoffs, using cutoff set {8, 4, 2, 1} Angstrom.
- **TM-score**: designed to correlate with expert assessment; computed (with GDT_TS)
  using MaxCluster.
- Reference structures downloaded from the PDB; for NMR ensembles the top-ranked model
  was used; missing residues excluded from scoring.

## Results and Discussion

### Refinement category (TR targets)

GDT and TM scores of the best GNEIMO structure improved over the starting decoy for
**19 of 23** proteins. GDT increased by up to **14.0** points; TM increased by up to
**0.13**. Average improvement over 23 targets: **+4.9 GDT_TS**, **+0.04 TM**,
**0.52 Angstrom RMSD**. Targets with > 5.0 GDT improvement (TR429, TR435, TR453,
TR454, TR464, TR476, TR530, TR557, TR574, TR624) all had at least 55% secondary
structure content. Targets with little/no refinement had less than 40% secondary
structure content. TR462 (two-domain, linker) improved 3.3 GDT overall but 5.7 and
6.2 in the individual domains, the linker limiting the global score. TR576 (parallel
beta-sheet native, decoy misfolded to antiparallel; crystal waters implicated) was
not refined; explicit water was not used.

### Structure-prediction category (T0 targets)

Homology models refined by GNEIMO-REXMD. All predicted structures ended within 4
Angstrom RMSD of the crystal structure except T0488. Average improvements:
**+4.5 GDT**, **+0.04 TM**, **0.7 Angstrom RMSD** (MODELLER starting RMSD 1-9 Angstrom).

### Enrichment and comparison to Cartesian MD

For TR429 and TR454 (both > 50% secondary structure), > 50% of the sampled
population shifted toward the native structure in both GDT and TM. TR568 and TR624
showed ~10-20% refined population. TR576/TR606/TR614/TR622 (< 40% secondary
structure, large loops) showed little shift. Using the identical force field and GBSA
solvation, the GNEIMO ensemble had a larger refined population than Cartesian MD;
torsional sampling was more effective at approaching native structures than all-atom
Cartesian sampling. (Shaw and co-workers observed multi-microsecond Cartesian MD
unraveled homology models on 25 CASP targets, 21 shared with this study.)

### Energy functions

All-atom AMBER99SB energies and Rosetta knowledge-based energies were computed for
conformations of three targets. Rosetta (CHARMM-based plus knowledge-based H-bond
terms) showed a more funnel-like character for some targets (near-native = lowest
energy), suggesting better best-structure selection, though neither AMBER nor Rosetta
was funnel-like for many targets. Best structures reported here were selected by RMSD
to native, not by an energy/scoring function.

### Assessment as a refinement tool

GNEIMO's larger stable time steps (up to 10 fs) let it explore more conformational
space per CPU cycle; teams replacing Cartesian MD with GNEIMO could simulate roughly
one order of magnitude longer in the same wall-clock time. Forces (not random Monte
Carlo moves) govern GNEIMO torsional moves; the search focuses on low-frequency DOF,
and REXMD supplies thermal energy to cross barriers arising from the stiff frozen
high-frequency DOF. GNEIMO is force-field agnostic (modular interface), integrator
modular, and can rigidify or free any DOF via the generalized coordinate system.

## Conclusions

GNEIMO-REXMD refined 30 CASP targets by up to 1.3 Angstrom starting from variable-
resolution homology models, refining 21 of 23 refinement targets; average refinement
over 23 targets was 4.0 GDT, 0.04 TM, 0.5 Angstrom RMSD, all without experimental
restraints. The extent of refinement was independent of starting decoy resolution
(> 5 Angstrom low-res vs < 3 Angstrom high-res). Refinement was effective for
secondary-structure regions and their packing; loop regions were the main failure
mode. Future directions: distance restraints from experiment, better scoring/energy
functions (including Rosetta-derived force fields), side-chain rotamer reassignment,
ensemble averaging, and combination with torsional Monte Carlo.
