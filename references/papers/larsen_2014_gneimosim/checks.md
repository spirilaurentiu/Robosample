# Checks / fixtures - Larsen 2014 (GneimoSim)

This is a software-description paper; it states no benchmark energies or
per-equation test values. The concrete numbers below are configuration/system
facts and reported qualitative scaling, useful as sanity references rather than
numerical regression fixtures.

## System / clustering facts

- Human alpha-2-macroglobulin: given the standard (default TMD) clustering
  scheme, expect **7165 clusters** from **20,426 atoms** (largest performance
  system). Rough ratio ~2.85 atoms/cluster for the default backbone-torsion
  clustering.

## Reported time-step regime

- ICMD/TMD stable integration time step used in GneimoSim: up to **10 fs**
  (default performance runs used **1 fs**; the usage example used **5.0 fs**).
- All-atom Cartesian MD with SHAKE on bond lengths: limited to ~**2 fs**.

## Performance-benchmark run parameters (Fig. 3, not per-value fixtures)

- 5 independent NVE ICMD simulations per protein, **300 ns** each.
- Integrator: Lobatto. Time step: **1 fs**. Long-range cutoff: **17 A**.
- No explicit solvent. Average time-per-step measured over the last **100 ns**.
- Hardware: Intel Xeon E5-2670 CPU + one Nvidia Tesla K20m GPU.
- Reported scaling: GNEIMO dynamics-solver cost is **linear (O(N))** in number of
  clusters/degrees of freedom; forcefield cost scales at higher order.

## Usage-example NVT run parameters (Nose-Hoover)

- step_size = 5.0 fs; num_step = 200000; cm_reset_frequency = 100 steps;
  bath_temperature = 300.0 K; temperature_relax_scale = 500 fs.
- Velocity initialization: `initTemperature(temperature=300, random_seed=111)`
  per the ICMD equipartition principle.

## Validation claims (qualitative, for cross-checking a port)

- All-atom GneimoSim simulations were cross-validated against LAMMPS and OpenMM
  all-atom simulations (expect agreement of a Cartesian-mode GneimoSim model with
  standard MD engines).
- Homology-model refinement: GNEIMO TMD + REMD refines models up to **1.5 A**
  without additional experimental restraints.
- Dynamic-clustering folding demonstrated on PDB IDs 1BDD (res. 11-56), 1EON
  (res. 7-31), 1PRB (res. 11-53), 1UB; folded to within **4-5 A** of crystal.
