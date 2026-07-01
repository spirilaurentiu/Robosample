# Checks / fixtures - Vitalis & Pappu 2014

Regression targets for a port of the mixed rigid-body/dihedral integrator.
Reduced-unit reminders: energies kcal/mol, T in K, dt in fs, mass in Da, lengths A.

## Flatness / MMT-artifact test (zero potential)
- System: linear 18-atom polymer (PEG-like heavy-atom geometry), **15 dihedral DOF + 6 rigid-body DOF**.
- Integrator: eq:11 with **Lambda = 4, dt = 5 fs**; Andersen thermostat, tau_T = 1 ps.
- Mass distributions tested: (i) equal 10 Da; (ii) 4->38 Da along chain in +2 Da steps; (iii) 6 triplets of (10, 5, 20) Da.
- Base of motion: atoms 8-10 (middle) or atoms 1-3 (terminus).
- EXPECT: with U = 0, all 15 dihedral-angle histograms are **flat (uniform)** for all mass distributions and both bases (no MMT artifacts). A Cartesian Langevin + SHAKE (bond+angle constraints) reference shows the artifact (non-flat). Bond-length-only constraints show negligible artifact.
- Sampling per line: 50 runs x 10 ns = 0.5 us cumulative; SHAKE converged to abs error < 1e-6 A (>100 iterations avg for full dihedral constraints).

## Integrator stability - two capped (GS)_50 polypeptides (excluded-volume)
- Box: cubic 200 A, PBC; U = 12th-power repulsion (cutoff 10 A) + amide-planarity dihedral potentials.
- Integrator eq:11, Lambda = 4. Metric: time to leave 100.5% of initial total energy (threshold drift ~6.5 kcal/mol/ns), 20 runs x 1 ns.
- EXPECT (correct atomic masses, ~340 K): stable up to **dt ~ 6 fs**; beyond, majority diverge within 1 ns. N-base vs M-base not significantly different.
- EXPECT (serine hydroxyl O and H masses both set to 8.5 Da): stable up to **dt ~ 16 fs** (serine chi_2 is the dominant error source).

## Liquid TIP4P water - Table I (thermal variables)
- System: 1095 TIP4P molecules, cubic box 32 A, PBC; 12 A nonbonded cutoff; reaction-field electrostatics; velocity-rescaling (Bussi) thermostat tau_T = 1 ps.
- Ideal temperature fluctuation for this size at 300 K: **sigma(T) = 5.23 K**.
- Reference: SETTLE-constrained Cartesian leapfrog.

dt = 2 fs (from 5e6 samples):
| variant | <T> (K) | <T_c> (K) | sigma(T) | sigma(T_c) | <U> per mol (kcal/mol) | C_v | <Delta T> (K) |
|---|---|---|---|---|---|---|---|
| Lambda=1 | 300.05 | 299.97 | 5.27 | 5.40 | -9.900 | 20.4 | 1.26 |
| Lambda=2 | 300.16 | 300.08 | 5.23 | 5.36 | -9.898 | 19.8 | 1.39 |
| Lambda=4 | 300.10 | 300.01 | 5.25 | 5.38 | -9.898 | 20.0 | 1.28 |
| SETTLE   | N/A    | 300.13 | N/A  | 5.26 | -9.899 | 20.5 | 1.40 |

dt = 5 fs (from 2e6 samples):
| variant | <T> (K) | <T_c> (K) | sigma(T) | sigma(T_c) | <U> per mol | C_v | <Delta T> (K) |
|---|---|---|---|---|---|---|---|
| Lambda=1 | 304.40 | 304.79 | 5.34 | 5.49 | -9.854 | 20.0 | 8.54 |
| Lambda=2 | 300.95 | 301.28 | 5.26 | 5.41 | -9.900 | 19.9 | 7.15 |
| Lambda=4 | 300.58 | 300.89 | 5.25 | 5.40 | -9.906 | 19.9 | 6.96 |
| SETTLE   | N/A    | 300.19 | N/A  | 5.26 | -9.919 | 19.8 | 8.14 |
- EXPECT: at dt=2 fs Lambda in {1,2,4} all match SETTLE; <T_c> fluctuation ~0.15 K larger than <T>; means within ~0.025%. At dt=5 fs larger Lambda gives better stability (temperature closest to 300 K). C_v ~ 20.0 cal/mol/K/molecule (MC reference Jorgensen & Madura: 20.0). C_v noisy to ~0.7 kcal/mol.
- SETTLE K_rot/K_trans split: rotational velocity from arcsin of quaternion vector part; I_xyz is full (non-diagonal) molecular inertia tensor.

## Liquid TIP4P water - Table II (dynamics)
dt = 2 fs:
| variant | D (1e-5 cm^2/s) | eps_r | tau_rot (ps) |
|---|---|---|---|
| Lambda=1 | 3.47 | 52.24 | 2.23 |
| Lambda=2 | 3.45 | 53.78 | 2.22 |
| Lambda=4 | 3.46 | 52.93 | 2.22 |
| SETTLE   | 3.48 | 52.48 | 2.20 |

dt = 5 fs:
| variant | D | eps_r | tau_rot (ps) |
|---|---|---|---|
| Lambda=1 | 3.53 | 51.65 | 2.16 |
| Lambda=2 | 3.37 | 53.56 | 2.30 |
| Lambda=4 | 3.34 | 53.52 | 2.33 |
| SETTLE   | 3.27 | 53.19 | 2.34 |
- D from linear MSD fits over 100 ps, restarts every 10 ps; error ~0.05e-5 cm^2/s. tau_rot error ~0.02 ps. ~2e4 samples via GROMACS utilities.
- EXPECT: at dt=2 fs all rigid-body integrators match SETTLE; at dt=5 fs discretization slows dynamics (D drops), clearest for SETTLE.

## Constant-energy water drift (Fig. 3c, last 100 ps)
- Total-energy drift for eq:11 with Lambda = 1, 2, 4, 8: **13.3, 2.5, 1.6, 1.6 kcal/mol/ps** respectively. Drift decreases sharply then plateaus by Lambda=4.

## TIP4P constant-energy stability limit (literature reference)
- SETTLE + correct masses: max stable dt ~ 7 fs for TIP4P.

## FS peptide (helix-coil), ABSINTH implicit solvent
- Sequence: N-Acetyl-A5(AAARA)3A-N'-methylamide (21 residues). Spherical droplet radius 40 A, explicit counterions + ~0.15 M NaCl; half-harmonic boundary spring 0.05 kcal/mol/A^2; 12 A cutoff; Andersen thermostat tau_T = 10 ps; eq:11 Lambda = 4.
- Blocking potential on backbone phi of every residue (to speed convergence away from left-handed basin).
- dt vs target T: **10.1 fs at 220 K -> 7.8 fs at 374 K** (decreases with T).
- Run length: 1.08e8 steps, first 1.8e7 discarded (0.7-0.9 us production).
- Helical segment defs: >=2 consecutive residues in alpha-basin = a segment (contributes to N_s); a segment of length N_alpha contributes N_alpha - 2 H-bonds to N_h; length-1 alpha runs contribute only to N_1.
- EXPECT: <N_h>, <N_s>, <N_1> vs <T> overlap the REMC reference within error; same melting temperature; M-base has smaller errors than N-/C-base. Ensembles thermodynamically identical, kinetically distinct (base-dependent transition rates).

## Timing (single core Intel Xeon E5410)
- TIP4P water, dt=5 fs: reference (Cartesian+SETTLE) and this integrator (Lambda=4) both ~1.6 ns/day.
- FS peptide, dt=2 fs: internal-coordinate (Lambda=4) ~40 ns/day vs Cartesian leapfrog+SHAKE ~41 ns/day.
- Zero-force system: 11.4 us/day (Cartesian Langevin+SHAKE) vs 9.5 us/day (eq:11).
- EXPECT: auxiliary O(N_at) recursions add negligible cost vs force evaluation (unlike GNEIMO-Fixman which costs ~2.24x FLEXIBLE vs 2x for plain torsional).
