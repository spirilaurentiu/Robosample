# Checks / Fixtures - Chen, Im, Brooks 2005 (TAMD)

All simulations: CHARMM with PARAM22 all-atom force field, vacuum (dielectric 1.0) unless noted, no nonbonded cutoff unless noted.

## ICFF construction parameters (charge neutralization when deleting bonded atoms)
- Valine $C_\gamma$ atoms: set charge to **0.0** (from **-0.27e**) after deleting 3 bonded hydrogens in the $\phi/\chi_1$,$\psi/\chi_1$ fragment.
- Serine/Threonine $O_\gamma$: set charge to **-0.23e** (from **-0.66e**) after deleting bonded H (H charge **0.43e**).
- Backbone $\phi/\psi$ crossterm fragment: `-CO-NH-C_\alpha H-NH-CO-`. Single shared $\phi/\psi$ correction map for all residues except proline.
- Side-chain fragments: `-CO-NH-C_\alpha-C_\beta H(X-)(Y-)` and `C_\beta H(X-)(Y-)-C_\alpha-NH-CO-` for $\phi/\chi_1$ and $\psi/\chi_1$. Val and Ile share the same map (grouped by $C_\beta$ connectivity).
- Correction maps applied via CHARMM CMAP facility; require fixed peptide planes.

## ICFF accuracy (valine dipeptide, and all standard residues except proline)
- RMSD between final ICFF $\phi/\psi$ map and original Cartesian map:
  - **0.67 kcal/mol** over low-energy region (within 5 kcal/mol of global minimum)
  - **1.6 kcal/mol** over the whole surface
- Backbone barrier increases (fixed-geometry alanine dipeptide, source of distortion):
  - near $(\phi,\psi)=(0^\circ,0^\circ)$: clash N-term carbonyl O vs C-term amide H (separated by 6 bonds)
  - near $\phi=120^\circ$: clash N-term carbonyl O vs side-chain methyl (separated by 4-5 bonds)

## Energy fluctuation / drift metric (eq:11), NVE, computed over >1.0 ps
- TAMD energy fluctuations are **2-3 orders of magnitude smaller** than CMD in all cases.
- CMD (SHAKE on H-bonds): conserves energy well up to **3 fs**; diverges quickly for time steps **> 5 fs**.
- TAMD: conserves energy well up to **10 fs** (small peptides) and up to **5 fs** (larger proteins); does **not diverge even at 20 fs**.
- Long-term drift (linear fit over 100 ps NVE, GB1): TAMD drift **< 0.1 kcal/ps** at 5 fs; lower than CMD at >=2 fs.
- CMD 5-fs run of GB1 diverged shortly after ~3 ps.
- Hydrogen mass set to **6.0 amu** in all TAMD runs (enables 10 fs even for compact GB1; less effective for large compact DHFR).

## Test systems (NVE energy-conservation study)
- Met-enkephalin: Ace-Tyr-Gly-Gly-Phe-Met-NMe (pentapeptide, collapsed coil)
- (Val)$_{10}$ polyvaline, canonical $\alpha$-helix
- GB1: 56-residue B1 domain of protein G, native $\alpha/\beta$ fold
- DHFR: 159-residue dihydrofolate reductase, native $\alpha/\beta$ with long loops
- GB1/DHFR native unstable in vacuum -> weak harmonic restraints on backbone heavy atoms **0.1 kcal/mol/Å$^2$**.

## REX folding simulations
- (Ala)$_{10}$: 8 replicas, 300-800 K, vacuum, no cutoff, TAMD time steps 2-20 fs. Time steps up to **20 fs** reliable -> at least **10-fold** speedup vs typical CMD. Average "folding" time = **29 ± 5 ps** (control repeats at 5 and 20 fs).
- (Val)$_{10}$: same setup as (Ala)$_{10}$. Average "folding" time **~250-500 ps** (backbone RMSD 1.0 Å criterion for folded helix), CMD and TAMD 2 fs.
- WALP16: sequence Ace-GWW(LA)$_5$WWA-NMe, GBSW implicit membrane (thickness **28.0 Å**), nonbonded switched off at **20 Å**. 16 replicas; CMD 300-600 K, TAMD 300-1000 K. TAMD REX at 5 fs fails to sample fully folded helix.
- Optimal number of replicas $\propto \sqrt{N_{\mathrm{dof}}}$; internal-coord $N_{\mathrm{dof}}$ ~ 1/10 of Cartesian.

## NMR refinement (HSP redox-switch domain, REX/GB, GBSW implicit solvent)
Setup: 16 replicas, 300-600 K, MD length 1.0 ps per REX step, 1000 REX steps, last 200 as production. Structured region residues 7-60 for RMSD.

Table 1 (RMSD in Å = backbone RMSD from final NMR structure ± RMS fluctuation; NOE = avg # restraints violated >0.2 Å / RMSD of NOE restraints in Å):

| Ensemble | RMSD (Å) | NOE |
|---|---|---|
| Initial | 8.8 ± 5.8 | 2.1 / 0.021 |
| CMD/2 fs | 2.2 ± 2.7 | 4.4 / 0.020 |
| TAMD/2 fs | 3.0 ± 2.4 | 5.0 / 0.024 |
| TAMD/3 fs | 3.1 ± 0.9 | 4.5 / 0.021 |
| TAMD/4 fs | 3.3 ± 2.7 | 5.7 / 0.026 |
| TAMD/5 fs | 3.1 ± 4.3 | 5.1 / 0.025 |

No NOE restraint violated by >0.5 Å in any ensemble. Time steps >=5 fs show limitations in sampling native states and ranking conformers.

## Algorithmic complexity
- Direct solve of eq:1 (dense mass-matrix inversion): $O(n^3)$.
- Recursive NEIMO / innovations factorization (eq:9-eq:10): $O(n)$ per step, $\Omega$ never formed explicitly.
- CMD time-step ladder rationale: H-bond constraints (SHAKE) -> up to 3 fs; all bond+angle constraints -> ~5 fs; beyond that, tip-group rotations (hydroxyl/methyl) and nonbonded heavy-atom collisions dominate.
