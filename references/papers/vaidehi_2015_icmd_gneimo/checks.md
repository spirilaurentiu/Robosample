# Checks / Fixtures - Vaidehi & Jain 2015, ICMD / GNEIMO

This is a review; most numbers are application benchmarks, not tight unit-test
fixtures. Still useful as sanity anchors for an ICMD/GNEIMO/Fixman implementation.

## Algorithm / parameter fixtures

- **GNEIMO solver scaling:** equations-of-motion solve is **O(N) (linear)** in
  the number of degrees of freedom (clusters); the naive dense mass-matrix
  inversion of eq 1 is **O(N^3) (cubic)**. Given a system where dof are doubled,
  expect ~2x GNEIMO EOM cost (not ~8x).
- **Nosé-Hoover thermostat mass:** given time step $\Delta t$, expect optimized
  $\tau = 10 \times \Delta t$ for torsional MD.
- **GNEIMO-Fixman overhead:** computing the Fixman potential (eq 5) + torque adds
  **~24%** to computation time versus the base GNEIMO solver.
- **Fixman validity:** including the Fixman potential should recover the
  unconstrained-Cartesian joint probability density of the two backbone torsion
  angles of **alanine dipeptide** (Fig 2: ICMD biased -> ICMD+Fixman matches
  Cartesian). PDF histogram bin size used: $d\theta = 18^\circ$ per axis.
- **System size sanity:** 7165 clusters corresponds to human
  alpha-2-macroglobulin with **20,426 atoms** (standard clustering).

## Application benchmarks (qualitative regression targets)

- **1BDD folding** (B domain of staphylococcal protein A), GNEIMO torsional MD +
  REMD, adaptive CVODE (Adams-Moulton): 12 temperature replicas from **300 to
  1050 K**. Trajectory: extended -> collapse to 12-16 Å backbone RMSD -> 8-11 Å
  (incorrect topologies) -> below 7 Å (correct topology); most conformations fall
  **5-7 Å** RMSD from crystal.
- **Calmodulin** (apo, no calcium) GNEIMO-REMD, compared to NMR ensemble PDB
  **1DMO**: captures central-helix collapse AND N-terminal domain flexibility;
  ~half of average residue hydrogen-bond distances fall within one standard
  deviation of NMR values.
- **BPTI** NVT torsional MD at **310 K, 100 ns**: reproduces backbone torsion-angle
  correlations between residues **9-18** and **35-40** (loops joined by a
  disulfide bond), matching millisecond-scale reference simulations.
- **CASP structure refinement** (GNEIMO-REMD + generalized Born, no experimental
  restraints): refinement of **up to 1.3 Å for 28 of 30** CASP target proteins.
  For target **T0453**, long-range contacts between loop residues 30-40 and
  residues 40-60 improve from **14-16 Å to 2-4 Å**.
- Homology modeling requires **>60%** sequence similarity to the template.
