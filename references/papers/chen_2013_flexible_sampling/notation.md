# Notation and conventions

| symbol / term | meaning | units / convention |
|---|---|---|
| DE (Delta E) | allowed relative conformational energy window above the lowest-energy conformer; conformers outside are discarded | kcal/mol |
| Duplicate RMS | heavy-atom RMSD cutoff for removing duplicate conformers | Angstrom (A) |
| RMSD (vs bioactive) | lowest heavy-atom RMSD of any ensemble member to the X-ray bioactive structure after best fit | Angstrom (A); thresholds 0.5, 1.0, 1.5, 2.0 |
| %BioConf_Rep | percent of compounds whose bioactive X-ray structure is reproduced within an RMSD threshold | percent |
| NbConfs | average number of conformers generated per compound per run | count |
| %GlobMin_found | percent of compounds for which a run located the global energy minimum (within 0.5 kcal/mol AND 0.5 A of reference) | percent |
| Max-Iterations | maximum total number of search iterations per compound | count |
| RotSteps | maximum number of search moves per rotatable bond (MacroModel generic methods only) | count |
| opr_nrot | Oprea number of rotatable bonds (MOE descriptor); counts exocyclic nonterminal single bonds and assigns some flexibility to aliphatic rings (1 per 5-ring, 2 per 6-ring, ... up to 8-membered; for >=9-membered rings all non-shared single bonds count) | count |
| b_1rotN | MOE descriptor counting rotatable bonds outside rings | count |
| Rgyr | radius of gyration of a conformer (compactness/extendedness measure) | Angstrom |
| Rgyr_min / Rgyr_max | min/max Rgyr across a computed ensemble for a compound | Angstrom |
| Rgyr_X-ray | Rgyr of the X-ray bioactive structure | Angstrom |
| Diel | distance-dependent dielectric solvation model | - |
| GB | generalized Born solvation model (recommended default) | - |
| macrocycle | ring of at least 9 atoms | - |
| Flexible compound | non-macrocyclic molecule with >= 12 rotatable bonds | - |

## Force fields

| name | note |
|---|---|
| MMFF94x | Merck Molecular Force Field, MOE default |
| MMFFs | Merck Molecular Force Field variant, MacroModel option (largely equivalent to MMFF94x) |
| OPLS2005 | MacroModel/Schrodinger default |
| OPLS2.0 | recent OPLS reparameterization |

## Software / method acronyms

| acronym | expansion |
|---|---|
| MOE | Molecular Operating Environment (Chemical Computing Group) |
| LowModeMD | low-mode + short-MD search method (MOE); channels motion along low-curvature directions without explicit mode calculation |
| Stochastic Search | random torsional-move search (MOE); no low-mode moves |
| LMOD | plain low-mode search (MacroModel) |
| LLMOD | large-scale low-mode search; eigenvectors without full Hessian diagonalization |
| MT/LMOD | mixed torsional / low-mode (MacroModel default) |
| MT/LLMOD | mixed torsional / large-scale low-mode |
| MD/LLMOD | MD-based simulated annealing then LLMOD; Schrodinger macrocycle protocol |
