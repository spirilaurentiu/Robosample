# Checks / fixtures - Constraint free-energy correction (Konig & Brooks 2014)

All energies in kcal/mol. T = 300 K. Force fields: CHARMM CGenFF / CHARMM22 / AMBER Cornell.

## Table 1 - (An)harmonic oscillator single-bond constraint free energies

Setup: one H atom bonded to a fixed non-interacting atom; constraint set at the
equilibrium bond length. `ΔG_cons = ΔH + ΔG_harm + ΔG_Jacobian`. Ion (Na+/Cl-) at
2.5 A from the H's initial position, projected along the bond. "Anharm. error" =
deviation of the harmonic approximation from the exact analytical result.

Gas phase (at equilibrium bond length -> ΔH = 0, ΔG_Jacobian = 0, anharm = 0):
| H type | ΔG_cons | ΔH | ΔG_harm | ΔG_Jacobian | Anharm err |
|---|---|---|---|---|---|
| HGA3 (ethane)     | -1.528 | 0 | -1.528 | 0 | 0 |
| HGA5 (ethene)     | -1.565 | 0 | -1.565 | 0 | 0 |
| HGR61 (benzene)   | -1.544 | 0 | -1.544 | 0 | 0 |
| HGP1 (methanol)   | -1.684 | 0 | -1.684 | 0 | 0 |
| HGP3 (ethanethiol)| -1.481 | 0 | -1.481 | 0 | 0 |

With Na+ at 2.5 A:
| H type | ΔG_cons | ΔH | ΔG_harm | ΔG_Jacobian | Anharm err |
|---|---|---|---|---|---|
| HGA3  | -1.512 | 0.017 | -1.530 | -0.002 | 0.0001 |
| HGA5  | -1.488 | 0.080 | -1.568 | -0.004 | 0.0001 |
| HGR61 | -1.519 | 0.027 | -1.546 | -0.002 | 0.0001 |
| HGP1  | -1.477 | 0.211 | -1.687 | -0.005 | 0.0003 |
| HGP3  | -1.423 | 0.060 | -1.483 | -0.004 | 0.0005 |

With Cl- at 2.5 A (steric clash):
| H type | ΔG_cons | ΔH | ΔG_harm | ΔG_Jacobian | Anharm err |
|---|---|---|---|---|---|
| HGA3  | -1.390 | 0.179 | -1.570 | -0.005 | 0.0000 |
| HGA5  | -1.576 | 0.019 | -1.595 | -0.001 | 0.0002 |
| HGR61 | -1.362 | 0.227 | -1.589 | -0.006 | 0.0007 |
| HGP1  | -1.461 | 0.218 | -1.678 | -0.005 | 0.0004 |
| HGP3  | -1.428 | 0.052 | -1.480 | +0.003 | -0.0001 |

- Max anharmonic error across all cases: 0.0007 kcal/mol (HGR61 + Cl-).
- Integration (Mathematica NIntegrate, AccuracyGoal=Infinity) done within +/-1.9 A of equilibrium bond length; result invariant to cutoff choice in [1,2] A.
- Numerical differentiation scale of variation: 0.0001 A.

## Table 2 - Four-atomic benchmark (no non-bonded terms; analytical reference)

Init state: all equilibrium bonds 2 A, bond force const 200 kcal/(mol A^2); angles
110 deg, angle force const 50 kcal/(mol rad^2); dihedral force const 1 kcal/mol,
multiplicity 3. Columns = deviation from analytical free energy (kcal/mol) for:
Uncorrected / NMD (standard normal modes, no Jacobian/curvature) / C-NMD (CBND+CANG,
full correction).

Bond-length changes Δr12:
| Δr12 (A) | Uncorrected | NMD | C-NMD |
|---|---|---|---|
| 0.10 | 2.00E+00 | 6.39E-05 | 6.91E-19 |
| 0.25 | 1.25E+01 | 2.38E-03 | 1.66E-18 |
| 0.50 | 5.00E+01 | 3.50E-02 | 5.38E-19 |
| 1.00 | 2.00E+02 | 4.71E-01 | 6.04E-18 |

Bond-angle changes Δφ123:
| Δφ123 (deg) | Uncorrected | NMD | C-NMD |
|---|---|---|---|
| 1  | 1.52E-02 | 6.74E-06 | 9.86E-19 |
| 5  | 3.81E-01 | 3.98E-03 | 3.34E-16 |
| 10 | 1.52E+00 | 5.95E-02 | 3.34E-16 |
| 25 | 9.52E+00 | 1.93E+00 | 3.37E-16 |
| 50 | 3.81E+01 | 2.40E+01 | 4.05E-19 |

Both terminal bonds+angles changed (Δr12=+1.0, Δr34=-1.0 A, Δφ123=+50, Δφ234=-50 deg):
| Uncorrected | NMD | C-NMD |
|---|---|---|
| 4.76E+02 | 8.53E+01 | 7.37E-18 |

- Full correction (C-NMD) reaches machine precision (errors 3.4E-16 to 4.1E-19 kcal/mol) because there is no anharmonicity.
- MUST correct bonds (CBND) before angles (CANG) for the angle-curvature result to hold.

## Table 3 - Water boxes (TIP3P), enthalpy-only correction (ΔH)

Perturbation: stretch all bonds by 0.05 A and angles by 1 deg; then release SHAKE and
apply ΔH correction. Deviation (kcal/mol) from the local potential-energy minimum.
CPU: single Q9300 2.50 GHz core.

| N_H2O | Uncorrected | Corrected | % Reduction | Time (s) |
|---|---|---|---|---|
| 5    | 11.2   | 0.2   | 98.2 | <0.01 |
| 395  | 263.6  | 20.0  | 92.4 | 29.0 |
| 787  | 538.7  | 53.9  | 90.0 | 227.4 |
| 1636 | 1111.0 | 105.0 | 90.5 | 2038.7 |
| 3290 | 2256.4 | 205.9 | 90.9 | 16494.8 |

- Uncorrected error ~0.6 kcal/mol per water; residual after correction ~0.05 kcal/mol per water.
- Error reduction consistently ~90%.
- Hessian cost scales cubically with atom count.
- Fast path (Hessian = force-field force constant): 3290 waters in 0.56 s (vs 16494.8 s full Hessian). A single energy call for that system = 1.87 s.

## Table 4 - Ethane -> methanol solvation free energy (kcal/mol)

Columns: 1.0 fs (unconstrained), 1.0 fs/SHAKE, 1.0 fs/corr (SHAKE + correction).
BAR; 5 ns gas / 1 ns solution. Experimental reference ΔΔA_solv = -6.93 kcal/mol.

| quantity | 1.0 fs | 1.0 fs/SHAKE | 1.0 fs/corr |
|---|---|---|---|
| ΔA_H2O | -0.86 +/- 0.02 | -0.70 +/- 0.03 | -0.90 +/- 0.03 |
| ΔA_gas | 6.024 +/- 0.006 | 6.039 +/- 0.002 | 6.027 +/- 0.002 |
| ΔΔA_solv | -6.89 +/- 0.02 | -6.74 +/- 0.03 | -6.93 +/- 0.03 |
| Deviation from exp | - | 0.15 | 0.04 |

Per-endpoint constraint free energies (ΔA_cons): ethane side 0.008; methanol side -0.180 kcal/mol.
- Unconstrained deviation from experiment: 0.04 kcal/mol.
- SHAKE without correction: deviation 0.15-0.19 kcal/mol (most error 0.18 on methanol side).
- With correction: deviation reduced back to 0.04 kcal/mol.

## Alanine <-> Serine free energy (Fig. 2 numbers)

- Unconstrained horizontal ΔG (ala->ser): 4.94 and -11.62 kcal/mol (agree with prior literature).
- Phase space overlaps unconstrained: 0.15% (CHARMM ala->ser), 0.55% (AMBER), 0.90% (CHARMM<->AMBER alanine), 0.55% (serine). All below the 1% rule-of-thumb threshold.
- Constraints raise overlap up to 2.05%; on average +94% (roughly double).
- Constraint-correction std devs < 0.003 kcal/mol (very low fluctuations of ΔG_cons).
- Total cube cycle-closing error: 0.22 kcal/mol (< 0.32 propagated std dev -> not significant).
  - constrained trajectories cycle error: 0.02 kcal/mol; unconstrained: 0.22 kcal/mol.
  - alanine trajectories: 0.04 kcal/mol; serine: 0.24 kcal/mol.
- Constrained std devs ~33% lower than unconstrained (avg); precision effectively doubled.
- Cost: full-Hessian correction ~14x a plain energy eval (3 h vs 12.8 min for 60,000 serine-in-water frames). Force-constant approximation: 17.7 min (~40% slower than plain eval).
- SHAKE lets simulation use ~3x larger time step; combined with 2 fs vs 1 fs -> ~8x cheaper than equivalent unconstrained free energy sim.

## Simulation parameters (for reproduction)

- Ethane->methanol: CHARMM22, dual-topology via MSCALE; gas Langevin (friction 5 ps^-1, 300 K, cutoff 998 A); solution 862 TIP3P waters, octahedral box 32.168 A, Nose-Hoover 300 K, LJ switched 10-12 A, PME, 1 fs step. Gas 5 lambda pts (0,0.25,0.5,0.75,1.0); solution 11 lambda pts (0.0..1.0 by 0.1). 4 repeats.
- Alanine->serine: N-acetyl-methylamide; AMBER Cornell + CHARMM22; 243 TIP3P waters, truncated octahedron (cube side 21.4 A), Nose-Hoover 300 K, LJ switched 9-10 A, PME. Trajectories every 10 steps. 4 repeats.
- Four-atomic and water box: VIBRAN, cutoff 10-12 A (no Ewald), force-shifting electrostatics.
