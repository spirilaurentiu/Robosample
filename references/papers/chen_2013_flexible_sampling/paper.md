# Tackling the conformational sampling of larger flexible compounds and macrocycles in pharmacology and drug discovery

I-Jen Chen, Nicolas Foloppe (Vernalis (R&D) Ltd, Cambridge, UK). Bioorg. Med. Chem. 2013, 21, 7898-7920. DOI: 10.1016/j.bmc.2013.10.003

> Routing / benchmark paper. Value is (i) a comparison of mainstream low-mode-based conformational sampling methods and (ii) concrete best-practice hyperparameters and benchmark numbers (see `checks.md`). No implementable equations of its own; the underlying algorithms live in the cited method papers (`depends_on`).

## Abstract

Computational conformational sampling underpins molecular modeling and design in pharmaceutical work. Sampling of smaller drug-like compounds is well studied, but few studies test sampling of larger, more flexible compounds relevant to drug discovery (therapeutic peptides, macrocycles, protein-protein-interaction inhibitors). This work tests mainstream conformational sampling methods on three curated compound sets: 'Drug-like' (253 compounds), 'Flexible' (50 compounds, >=12 rotatable bonds, non-macrocyclic), and 'Macrocycle' (30 compounds). Compared methods: Stochastic Search and LowModeMD from MOE; the low-mode variants LMOD, LLMOD, MT/LMOD, MT/LLMOD from MacroModel; and MD/LLMOD (Schrodinger) for macrocycles. Performance is assessed by (i) reproduction of X-ray bioactive structures, (ii) size/coverage/diversity of output ensembles, (iii) compactness/extendedness (radius of gyration), and (iv) ability to locate the global energy minimum. Search parameters enhanced above defaults give much better results while staying tractable. LowModeMD emerged as the method of choice in MOE; mixed torsional/low-mode (MT/LMOD) matched it in MacroModel; MD/LLMOD performed well for macrocycles.

## 1. Introduction (method landscape)

To generate conformers, methods rely on an energy model plus a sampling algorithm acting on conformational degrees of freedom. Search methods typically have a stochastic element:

- **Torsional Monte Carlo moves** (Chang, Guida, Still 1989).
- **Random pulses in Cartesian coordinates** (Ferguson & Raber 1989).
- **Seed inter-atomic distances** in distance-geometry calculations.
- **Initial velocities** of a molecular dynamics simulation.
- **Chosen low-modes and seed conformer** in low-mode searches.

### Low-mode search principle

The low-frequency vibrational-mode eigenvectors point along low-energy paths connecting energy minima via saddle points on the conformational energy surface. Moving the coordinates along the low-frequency modes is an efficient way to cross energy barriers between minima; once a low-mode move locates a new energy well, energy minimization is performed and another search cycle begins. The search works in the reduced-dimensionality space of the low-frequency modes and is well-adapted to cyclic topologies. Searches can be tuned by controlling how frequently eigenvectors are re-computed.

Two ways to apply low-mode search to large systems:

- **LLMOD (Large-scale low-mode)** generates eigenvectors *without* explicitly diagonalizing the entire Hessian (Kolossvary & Keseru 2001).
- **LowModeMD** (Labute 2010, in MOE) does *not* compute low-frequency modes explicitly; it efficiently channels and amplifies atomic motions along directions of low curvature of the potential energy surface, via a short MD run at the beginning of every iteration.

MacroModel variants:
- **LMOD** - plain low-mode search.
- **LLMOD** - large-scale low-mode (approximate eigenvectors, cheaper for large systems).
- **MT/LMOD** - Mixed torsional / low-mode (adds random torsional moves). Default in MacroModel/Maestro.
- **MT/LLMOD** - Mixed torsional / large-scale low-mode.
- **MD/LLMOD** - Schrodinger's specialized macrocycle protocol: high-temperature MD-based simulated annealing, then LLMOD.

## 2. Methods

### 2.1 Compound test sets

- **Drug-like set:** 253 diverse compounds, MW and rotatable-bond count in conventional drug-like range (1-13 rotatable bonds).
- **Flexible set:** 50 diverse non-macrocyclic ligands with 12 <= opr_nrot <= 20 (in practice 12-17 rotatable bonds), bound to 32 diverse protein families. Filtered from high-quality X-ray structures (resolution 2 A or better, structure factors deposited). Removed metal-containing, >2 charged centers at pH 7, alkyl chains >=4 carbons, >10 chiral centers, and >20 noncyclic nonterminal rotatable bonds.
- **Macrocycle set:** 30 compounds (macrocycle = ring of at least 9 atoms), 9-30 rotatable bonds. Filtered by b_1rotN < 10, clustered at 70% MACCS Tanimoto similarity.

All compounds prepared in MOE (bond orders, pH-7 protonation, tautomers consistent with X-ray binding). Converted to 2D before sampling to erase memory of the X-ray conformation; R/S stereochemistry of chiral centers preserved. Every compound except one (PDB 3OMJ, bound to DNA) is protein-bound, noncovalently.

### 2.2-2.4 Search algorithms

- **MOE** version 2011.10; LowModeMD and Stochastic Search with MMFF94x. Three independent runs per protocol.
- **MacroModel** BatchMin V9.9 (Maestro 9.3.5). Generic methods LMOD, LLMOD, MT/LMOD, MT/LLMOD. Energy-minimize with Polak-Ribiere conjugate gradient; convergence criterion energy gradient <= 0.05 kJ/mol/A; number of minimization iterations set to 3000 (default 500 insufficient). 3D input via LigPrep. Three independent runs per protocol.
- **MD/LLMOD** for macrocycles: stage 1 = high-T MD-based simulated annealing (default 5000 cycles, 1000 K -> 300 K, then minimization); stage 2 = LLMOD (5000 search steps default), eigenvectors recomputed on each new global energy minimum by default. One run per protocol.

### 2.5-2.6 Search parameters and energy models

Tunable parameters investigated (see Table 1 in `checks.md`): energy window DE (kcal/mol), solvation (distance-dependent dielectric "Diel" vs generalized Born "GB"), Duplicate RMS cutoff (A) for removing similar conformers, Max-Iterations, and RotSteps (max moves per rotatable bond, MacroModel only). Force fields: MMFF94x/MMFFs (Merck), OPLS2005, OPLS2.0.

### 2.7 Reproduction of bioactive structures

For each compound, find the lowest heavy-atom RMSD between any ensemble member and the X-ray bioactive reference after best fit. %BioConf_Rep = percentage of bioactive structures reproduced within a given RMSD threshold (0.5, 1.0, 1.5, 2.0 A). The 1.0 A threshold is emphasized. Significance of a given RMSD depends on molecule size (a fixed RMSD is better performance for a larger compound).

### 2.8 Conformational coverage via 3D descriptors

- **NbConfs:** average number of generated conformers per compound (higher = broader coverage).
- **Radius of gyration Rgyr:** quantifies compactness/extendedness range. Rgyr_X-ray compared to computed Rgyr_min / Rgyr_max.
- **3D three-point pharmacophores:** counted with Schrodinger Canvas fingerprints, features = H-bond acceptor, H-bond donor, hydrophobic, negative charge, positive charge, aromatic ring. Distances binned in 2 A bins [0,2), [2,4), ..., [20, inf). Total number of distinct nonredundant pharmacophores across three merged runs = coverage estimate.

### 2.9 Global energy minima and convergence

The plausible global energy minimum per compound and energy model is identified empirically as the lowest-energy conformer across all aggregated runs (MOE: 12/15/15 runs for Drug-like/Flexible/Macrocycle; MacroModel: 36/39/24 runs). A run is deemed to have "found" the global minimum if it produces a conformer within 0.5 kcal/mol and 0.5 A of the reference. %GlobMin_found = percent of compounds for which a run located the global minimum. Frequency of finding the global minimum across independent runs is a (simplified) convergence test.

## 3. Results and discussion (key findings)

- **Solvation dominates.** GB consistently yields more conformers and reproduces more bioactive structures than distance-dependent dielectric (Diel). In MOE (Diel default) switching to GB raised Drug-like %BioConf_Rep (1 A) from 77% to 91% for LowModeMD; Flexible from 9% to 47%. Recommendation: use GB.
- **Energy window.** Widening DE from default (7 kcal/mol MOE, 5 kcal/mol MacroModel) to 15 kcal/mol markedly improves %BioConf_Rep for Flexible and Macrocycle compounds, then plateaus (little gain 15 -> 20). For Drug-like the effect is marginal.
- **RotSteps** (MacroModel). Increasing moves per rotatable bond from 100 (default) to 400 dramatically raises %BioConf_Rep for Flexible (45% -> 65%) and Macrocycle (64% -> 79%) sets; more impactful than choosing MT/LMOD vs MT/LLMOD. Doubling RotSteps roughly doubles NbConfs and compute time.
- **Method ranking.** For Drug-like compounds LowModeMD and Stochastic Search are equivalent. For Flexible/Macrocycle sets LowModeMD beats Stochastic Search (larger ensembles, more often finds global min). MT/LMOD approx MT/LLMOD; both beat pure LMOD/LLMOD (torsional moves help). MD/LLMOD default settings match enhanced LowModeMD for macrocycles; its wide default DE (10 kcal/mol) and the initial MD simulated-annealing phase are essential (turning off the MD stage drops %BioConf_Rep from 77% to 60%).
- **Force field.** MMFF, OPLS2005, OPLS2.0 give similar %BioConf_Rep and NbConfs; sampling extent and electrostatics treatment dominate, not the force field.
- **Coverage.** Enhanced protocols broaden Rgyr_max (more extended) and lower Rgyr_min (more compact), covering the Rgyr of the bioactive X-ray structures. They also visit more 3D pharmacophores (e.g. Flexible: enhanced LowModeMD 28,185 vs default 19,021).
- **Convergence.** %GlobMin_found >= 97% for enhanced protocols on Drug-like (excellent convergence). For Flexible/Macrocycle it is lower and more variable across runs; collate several independent runs in practice.
- **Cost.** Enhanced MT/LMOD and enhanced LowModeMD run about 6x and 7x longer than their default counterparts.

## 4. Conclusions

With suitable protocols, conformational sampling is hardly a bottleneck for Drug-like compounds. For larger, more flexible compounds and macrocycles, mainstream low-mode-based methods give robust and encouraging results provided parameters are enhanced above defaults: use GB solvation, DE = 15 kcal/mol, Duplicate RMS = 0.25 A, Max-Iterations = 10,000, and (MacroModel) RotSteps = 400. LowModeMD (MOE) and MT/LMOD or MT/LLMOD (MacroModel) are the recommended methods; MD/LLMOD default settings are well-balanced for macrocycles.
