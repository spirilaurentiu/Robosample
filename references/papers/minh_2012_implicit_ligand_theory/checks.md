# Checks / Fixtures - Implicit Ligand Theory (Minh 2012)

System: Cucurbit[7]uril (CB[7]) receptor with 12 ligands (adamantanes AD1-AD5,
bicyclooctanes B02/B05/B11, ferrocenes F01/F02/F03/F06) in water.
T = 300 K. Free energies in kcal/mol. Values are mean (std dev).
Solvent: GBSA implicit (surface tension 0.006 kcal/mol/Å², receptor dielectric 1.0,
solvent dielectric 78.5). Alchemical coupling run in vacuum.

## Demonstration hyperparameters (Methods)
- MD engine: modified NAMD 2.9, Langevin dynamics, 1 fs time step.
- 1-4 electrostatics scaled by 0.5; nonbonded cutoff 999 Å (effectively no cutoff).
- Receptor: minimized 2500 steps; thermalized 0→300 K (+10 K, reinit velocities every 100 steps); snapshots every 0.1 ns over 10 ns → 100 receptor snapshots.
- λ ladder: {0, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0}. Soft-core vdW shift coefficient 5. Electrostatics turned on at λ=0.5.
- Restraint: flat-bottom harmonic, spring constant 10 kcal·mol⁻¹·Å⁻¹, flat region to 0.75 Å, on COM distance (ligand core heavy atoms ↔ receptor heavy atoms).
- Binding-site volume Ω ≈ (4/3)π(0.75³)(8π²).
- HREX: 1000 exchange attempts per 5 ps per λ; 25 cycles; snapshots every 0.5 ps → 2 ns total per binding PMF. Equilibration cutoff: energy of fully coupled state within 20 k_BT of final snapshot.
- Ligand reservoir: 10 ns simulation, snapshots every 10 ps; reservoir swap at λ=0 satisfies detailed balance.
- PBSA outlier filter: drop ligand snapshots whose total PBSA energy exceeds final-snapshot PBSA by ≥ 100 k_BT.
- Poisson-Boltzmann (UHBD): grid spacing 0.18 Å, molecule ≤ 0.7 of grid.
- Binding PMF std dev across ligands: ranges 0.12 to 1.63 kcal/mol.
- Force-field sensitivity: ligand B11 binding PMF shifts ~40 kcal/mol between force fields.
- Convergence: mean B(r_R) stops shifting after ~0.75 ns; ΔG° average stabilizes after ~15 receptor snapshots.

## Table I - Binding PMF B(r_R) for minimized CB[7], Eq (23), mean(std)
| Ligand | NAMD | M2 | PB | PBSA | min{Ψ(r_RL)} |
|---|---|---|---|---|---|
| AD1 | -14.1 (0.79) | -22.0 (0.51) | -23.0 (0.82) | -25.5 (0.83) | -31.3 (0.55) |
| AD2 | -32.5 (0.15) | -29.0 (0.13) | -26.8 (0.12) | -29.4 (0.12) | -36.9 (0.30) |
| AD3 | -31.0 (0.16) | -30.7 (0.18) | -28.9 (0.23) | -31.6 (0.23) | -40.3 (0.28) |
| AD4 | -44.0 (0.94) | -36.7 (1.11) | -24.0 (1.12) | -26.9 (1.12) | -36.1 (0.45) |
| AD5 | -32.2 (0.68) | -29.0 (0.25) | -26.0 (0.14) | -28.5 (0.14) | -36.2 (0.29) |
| B02 | -12.8 (0.41) | -18.8 (0.38) | -19.8 (0.53) | -22.6 (0.53) | -30.6 (0.54) |
| B05 | -40.4 (0.29) | -30.6 (0.40) | -19.5 (0.50) | -22.3 (0.50) | -34.3 (0.63) |
| B11 | -52.4 (1.50) | -38.5 (1.81) | -14.1 (1.72) | -17.5 (1.63) | -39.3 (1.90) |
| F01 | -1.8 (1.57) | -5.5 (0.81) | -10.9 (0.53) | -13.6 (0.53) | -24.7 (0.33) |
| F02 | -14.8 (0.96) | -14.1 (0.67) | -16.2 (0.42) | -19.2 (0.42) | -31.7 (0.38) |
| F03 | -16.3 (1.61) | -13.1 (1.03) | -16.4 (0.95) | -19.5 (0.94) | -31.1 (0.71) |
| F06 | -30.6 (0.18) | -20.0 (0.21) | -21.9 (0.18) | -25.4 (0.18) | -37.2 (0.24) |
| R²_ITC | 0.884 | 0.750 | 0.454 | 0.490 | 0.883 |
| RMSE_ITC | 12.8 | 7.9 | 4.9 | 4.7 | 12.2 |
| R²_Gilson | 0.827 | 0.907 | 0.705 | 0.712 | 0.792 |
| RMSE_Gilson | 10.4 | 4.8 | 5.4 | 4.5 | 11.3 |
(Dominant-state approximation: ΔĜ° from a single B̂(r_R) as estimate of −β⁻¹ ln⟨e^{−βB}⟩.)

## Table II - ΔG° via Eq (22), 100 receptor snapshots (PBSA-based B), mean(bootstrap std)
| Ligand | ITC | Gilson(M2) | NAMD | M2 | PB | PBSA |
|---|---|---|---|---|---|---|
| AD1 | -14.1 | -18.2 | -9.4 (0.23) | -16.3 (0.15) | -17.6 (0.25) | -20.1 (0.25) |
| AD2 | -19.4 | -25.9 | -27.9 (0.19) | -24.3 (0.22) | -22.9 (0.27) | -25.4 (0.26) |
| AD3 | -20.4 | -25.6 | -35.7 (5.03) | -28.6 (1.87) | -23.5 (0.23) | -26.2 (0.23) |
| AD4 | -21.5 | -29.7 | -40.5 (0.21) | -33.7 (0.32) | -24.3 (1.11) | -27.1 (1.06) |
| AD5 | -19.1 | -24.1 | -29.5 (1.24) | -24.0 (0.20) | -22.0 (0.35) | -24.4 (0.34) |
| B02 | -13.4 | -12.0 | -9.0 (0.38) | -13.7 (0.16) | -15.4 (0.26) | -18.1 (0.25) |
| B05 | -19.5 | -23.1 | -38.0 (0.40) | -27.7 (0.27) | -18.6 (0.27) | -21.4 (0.27) |
| B11 | -20.6 | -22.4 | -51.2 (0.34) | -37.3 (0.24) | -17.2 (0.53) | -20.5 (0.51) |
| F01 | -12.9 | -10.2 | 0.3 (0.82) | -0.6 (0.34) | -4.9 (0.26) | -7.6 (0.25) |
| F02 | -16.8 | -12.4 | -12.0 (0.70) | -9.6 (0.75) | -11.7 (0.70) | -14.6 (0.71) |
| F03 | -17.2 | -12.2 | -10.2 (0.16) | -7.3 (0.24) | -10.2 (0.22) | -13.2 (0.22) |
| F06 | -21.0 | -17.8 | -24.1 (0.34) | -14.1 (0.46) | -16.2 (0.51) | -19.7 (0.52) |
| R²_ITC | | 0.782 | 0.870 | 0.745 | 0.671 | 0.704 |
| RMSE_ITC | | 4.6 | 14.0 | 9.0 | 4.4 | 4.5 |
| R²_Gilson | | | 0.841 | 0.892 | 0.923 | 0.925 |
| RMSE_Gilson | | | 11.3 | 5.9 | 3.4 | 2.4 |
(Bootstrap: 1000 random selections of 100 binding PMFs.)

## Table III - Mean potential-energy change on binding (PBSA), kcal/mol, mean(bootstrap std)
Columns: VDW, Coul, PB(electrostatic solvation), Val(bond+angle+dihedral), NP(nonpolar solvation), Total.
| Ligand | VDW | Coul | PB | Val | NP | Total |
|---|---|---|---|---|---|---|
| AD1 | -32.5 (0.471) | 0.1 (1.509) | 4.8 (1.547) | -5.0 (2.785) | -2.5 (0.011) | -35.2 (2.547) |
| AD2 | -33.6 (0.931) | -65.8 (1.032) | 64.9 (0.783) | -5.9 (1.910) | -2.5 (0.017) | -42.9 (1.693) |
| AD3 | -32.8 (0.718) | -64.4 (0.693) | 62.2 (0.855) | -5.7 (2.128) | -2.6 (0.009) | -43.4 (2.388) |
| AD4 | -38.1 (1.400) | -125.2 (3.283) | 124.4 (1.003) | 1.9 (4.475) | -2.7 (0.070) | -39.9 (2.817) |
| AD5 | -33.3 (1.374) | -65.1 (1.549) | 64.8 (1.415) | -4.9 (1.782) | -2.5 (0.025) | -40.9 (1.834) |
| B02 | -33.3 (0.622) | -5.8 (1.067) | 9.7 (0.770) | 1.4 (3.379) | -2.7 (0.022) | -30.6 (2.187) |
| B05 | -32.9 (0.896) | -138.2 (1.231) | 138.0 (1.151) | -2.3 (1.236) | -2.8 (0.013) | -38.1 (1.673) |
| B11 | -39.9 (1.192) | -199.3 (2.280) | 212.0 (1.066) | -5.6 (5.431) | -3.4 (0.075) | -36.2 (4.475) |
| F01 | -26.2 (0.497) | -8.2 (1.824) | 14.2 (1.119) | 8.7 (4.842) | -2.7 (0.017) | -14.3 (4.079) |
| F02 | -26.9 (1.518) | -65.7 (2.078) | 65.9 (0.923) | -0.9 (2.237) | -3.0 (0.012) | -30.6 (2.253) |
| F03 | -28.7 (0.832) | -58.0 (0.987) | 64.2 (0.649) | -0.4 (3.552) | -3.0 (0.015) | -26.1 (3.428) |
| F06 | -35.1 (1.154) | -116.1 (0.810) | 120.9 (0.651) | -8.8 (4.862) | -3.5 (0.013) | -42.6 (4.132) |

## Table IV - ΔG° (PBSA): dominant-state vs HREX
Columns: [min{Ψ}→min{B̂}], [min{Ψ}→EXP], [HREX→min{B̂}], [HREX→EXP(Eq 22)].
| Ligand | DS,min B̂ | DS,EXP | HREX,min B̂ | HREX,EXP |
|---|---|---|---|---|
| AD1 | -28.6 | -27.2 | -22.0 | -20.1 |
| AD2 | -36.4 | -34.6 | -27.6 | -25.4 |
| AD3 | -38.1 | -36.8 | -27.6 | -26.2 |
| AD4 | -43.1 | -40.4 | -29.8 | -27.1 |
| AD5 | -35.8 | -33.6 | -26.8 | -24.4 |
| B02 | -29.8 | -27.9 | -21.0 | -18.1 |
| B05 | -37.9 | -35.6 | -23.7 | -21.4 |
| B11 | -48.5 | -45.7 | -23.1 | -20.5 |
| F01 | -22.7 | -21.3 | -10.2 | -7.6 |
| F02 | -30.9 | -28.8 | -17.0 | -14.6 |
| F03 | -28.7 | -27.0 | -14.5 | -13.2 |
| F06 | -35.6 | -33.8 | -21.3 | -19.7 |
| R²_ITC | 0.849 | 0.855 | 0.684 | 0.704 |
| RMSE_ITC | 17.3 | 15.3 | 5.8 | 4.5 |
| R²_Gilson | 0.787 | 0.795 | 0.926 | 0.925 |
| RMSE_Gilson | 15.8 | 13.9 | 3.5 | 2.4 |
| R²_Exp | 0.723 | 0.736 | 0.996 | |
| RMSE_Exp | 15.5 | 13.6 | 2.3 | |

## Table SI - Binding PMF components (Eq 23), minimized CB[7], mean(std)
Charge = net formal charge. Columns: B_cpl, then B_RL/B_L per force field (NAMD, M2, PB, PBSA).
| Ligand | Charge | B_cpl | B_RL,NAMD | B_L,NAMD | B_RL,M2 | B_L,M2 | B_RL,PB | B_L,PB | B_RL,PBSA | B_L,PBSA |
|---|---|---|---|---|---|---|---|---|---|---|
| AD1 | 0 | -31.6 (0.09) | -115.9 (0.78) | -2.8 (0.01) | -111.5 (0.50) | -1.1 (0.01) | -128.1 (0.81) | -4.5 (0.01) | -122.3 (0.81) | -2.6 (0.01) |
| AD2 | 1 | -91.3 (0.12) | -123.7 (0.10) | -51.9 (0.04) | -112.8 (0.04) | -55.2 (0.04) | -123.4 (0.02) | -55.8 (0.03) | -117.6 (0.02) | -53.8 (0.03) |
| AD3 | 1 | -93.2 (0.17) | -118.9 (0.06) | -50.5 (0.04) | -108.9 (0.05) | -51.4 (0.04) | -122.3 (0.10) | -54.5 (0.04) | -116.4 (0.10) | -52.3 (0.04) |
| AD4 | 2 | -148.0 (0.41) | -202.5 (0.62) | -175.9 (0.63) | -192.2 (0.75) | -183.5 (0.58) | -188.8 (0.87) | -180.6 (0.54) | -182.8 (0.86) | -178.1 (0.54) |
| AD5 | 1 | -90.5 (0.13) | -124.3 (0.76) | -52.1 (0.05) | -111.9 (0.29) | -53.5 (0.02) | -123.5 (0.03) | -55.9 (0.03) | -117.6 (0.03) | -53.9 (0.03) |
| B02 | 0 | -31.9 (0.09) | -118.4 (0.37) | -6.9 (0.05) | -111.4 (0.33) | -4.5 (0.03) | -129.1 (0.49) | -9.0 (0.02) | -123.3 (0.49) | -6.8 (0.02) |
| B05 | 2 | -156.8 (0.13) | -187.3 (0.27) | -173.2 (0.11) | -176.4 (0.35) | -182.7 (0.18) | -173.4 (0.45) | -178.6 (0.15) | -167.5 (0.45) | -176.4 (0.15) |
| B11 | 4 | -220.7 (0.83) | -437.9 (1.24) | -475.7 (0.65) | -422.5 (1.20) | -484.7 (0.97) | -404.2 (1.93) | -478.7 (0.74) | -396.9 (1.80) | -474.3 (0.75) |
| F01 | 0 | -25.2 (0.27) | -117.2 (1.50) | -10.1 (0.10) | -65.0 (0.75) | 35.3 (0.07) | -92.5 (0.45) | 25.4 (0.07) | -86.6 (0.46) | 27.6 (0.08) |
| F02 | 1 | -84.2 (0.24) | -112.1 (0.87) | -50.9 (0.07) | -77.9 (0.57) | -28.0 (0.11) | -96.6 (0.32) | -32.5 (0.08) | -90.7 (0.32) | -29.9 (0.08) |
| F03 | 1 | -82.1 (0.51) | -111.8 (2.04) | -47.0 (0.04) | -76.2 (1.45) | -25.2 (0.07) | -97.0 (1.37) | -30.5 (0.05) | -91.1 (1.36) | -27.9 (0.05) |
| F06 | 2 | -144.6 (0.15) | -152.3 (0.09) | -135.7 (0.10) | -123.3 (0.12) | -127.9 (0.18) | -139.0 (0.09) | -129.5 (0.14) | -132.9 (0.09) | -126.3 (0.14) |

## Consistency identity (estimator sanity check)
For the observable O = 1, the interaction-weighted rigid-receptor expectation collapses to
⟨Θ̂(r_R)⟩_R = ⟨e^{−βB}⟩_R (an estimate of a constant returns that constant); this ⟨e^{−βB}⟩_R
is reused as the denominator of Eq (13).
