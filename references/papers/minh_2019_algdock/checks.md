# Checks - Minh 2019 (AlGDock)

Regression fixtures for an implementation. Free energies in $k_B T$; BPMF total is $f_{AE}$.

## Constants and hyperparameters (fixtures for the protocol code)

- Binding-site restraint: `k = 10000` kJ/(mol nm^2), `d_0 = 6.0` Å; given `d <= 6.0 Å` expect `u_I = 0`; given `d > 6.0 Å` expect `u_I = 0.5 * beta * k * (d - 6.0Å)^2`.
- Milestone temperatures: `T_T = 300 K`, `T_H = 600 K`.
- Temperature ramp: `T(alpha) = (300 - 600)*alpha + 600`; given `alpha=0` expect `T=600 K`; `alpha=1` expect `T=300 K`; `alpha=0.5` expect `T=450 K`.
- Soft-grid scaling: `alpha_sg(alpha) = -(2*alpha-1)^2 + 1`; given `alpha=0` expect `0`; `alpha=0.5` expect `1`; `alpha=1` expect `0`.
- Unperturbed-grid scaling: `alpha_g(alpha) = (2*alpha-1)^2 / (1 + exp(-1000*(alpha-0.5)))`; given `alpha=0` expect `~0` (numerator 1, denominator ~exp(500) -> ~0); given `alpha=1` expect `~1` (numerator 1, denominator ~1); given `alpha=0.5` expect `0` (numerator 0).
- Soft LJ repulsive cap: `v_max = 10.0` kJ mol^{-1/2}; soft transform `v = v_max * tanh(v_o / v_max)`.
- Thermodynamic speeds: `s_bc = 20.0`, `s_cd = 0.2`.
- vdW grid inverse-transform power = 4 (repulsive), none (attractive).
- PB grid: protein dielectric 2.0, solvent dielectric 80.0, solvent radius 1.4 Å, fine spacing 0.5 Å, T = 300 K.
- HMC: 50 velocity-Verlet steps, adaptive dt in [0.1, 5.0] fs; converged dt typically 2.75-3.75 fs.
- External-coordinate MCMC: translation std 0.6 Å per dimension; attempted only when `alpha < 0.01`.
- Time-step adaptation: acc > 0.8 -> dt += 0.125 fs; acc < 0.4 -> dt -= 0.25 fs; acc < 0.1 -> dt -= 0.5 fs; target acc in [0.4, 0.8].
- Replica exchange: 25 sweeps/cycle; pairs 1..min(5,K) apart; insert new state if neighbor exchange rate < 0.4; remove state k if <p_acc> > 0.99 during init; reselect (speed x 4/5) if <p_acc> < 0.4.
- Cycles (demonstrative): 8 for states BC, 15 for states CD; 1000 iterations/cycle; 50 snapshots per RE cycle.
- Clustering threshold: 1.0 Å (Hungarian symmetry-corrected heavy-atom RMSD).
- Native pose criterion: RMSD < 2 Å from crystal; DOCK 6 poses within 6 Å COM of crystal, min anchor size 5.

## Number of thermodynamic states (Astex set, 11 sims)

- States BC: Desolvated 67-182; Full 51-111. `sigma[N_states] < 2` for all systems (BC).

## BPMF precision (aggregate, demonstrative calculations)

- Within 1 kcal/mol (= 1.68 RT): 28.2% of systems.
- Within 4 kBT: 75.3% (Desolvated), 74.1% (Full).
- Within 8 kBT: 87.1% (Desolvated), 94.1% (Full).
- Desolvated/Full agree within error: 63 systems (74.1%).

## Convergence (RMSE vs final estimate)

- f_AB: < 1 kBT after 1 cycle (all); < 0.5 kBT after 8 cycles except 1p62 Desolvated (0.544 kBT).
- f_BC Desolvated: < 2 kBT after 1 cycle except 1jje (3.2 kBT); < 1 kBT after 8 cycles. Full: < 1 kBT after 1, < 0.15 kBT after 8.
- f_CD Desolvated: < 2.5 kBT after 15 cycles except 1t40 (2.65). Full: < 2.5 kBT after 15 except 1l7f (3.41).

## Native pose identification success (Astex; Table 1)

Fraction of calculations where a native pose is the minimum-energy / within cutoff. Format: value (std err of proportion). "min u" = min total energy, "min Ψ" = min interaction energy, "fe" = free-energy reweighted. Cutoffs 8, 4, 2, 0 kBT.

xtal + DOCK 6, DOCK6 grid score, min Ψ: 0.941(0.026), 0.906(0.032), 0.847(0.039), 0.812(0.042).
xtal + DOCK 6, milestone E, min u: 0.953, 0.953, 0.929(0.028), 0.929.
xtal + DOCK 6, milestone E, min Ψ: 0.953, 0.918(0.030), 0.882(0.035), 0.871(0.036).
Full, milestone E, min u (cutoff 0): 0.820(0.013).
Full, milestone E, min Ψ (cutoff 0): 0.846(0.012).
Desolvated, milestone E, min u (cutoff 0): 0.771(0.014).
All, milestone E, min Ψ (cutoff 0): 0.871(0.036).

External benchmarks on same Astex set: GOLD 80.5%, GLIDE 82%, ICM 91%. This work (milestone E free-energy ranking): Desolvated 75.8%, Full 83.9%.

Native pose observed during milestone D production: Desolvated 90.4%, Full 94.9%.

## Representative per-system free energies (Table S2, kBT; std dev in parens)

Decomposition: total `f_AE ≈ f_AB + f_BC + f_CD + f_DE` (reduced). Sample fixtures:

Desolvated:
- 1jla: fAB=62.88(0.04), fBC=-154.13(0.06), fCD=-137.10(0.15), fDE=-3.46(0.21), fAE=-49.31(0.26).
- 1tz8: fAB=23.80(0.02), fBC=3.33(0.05), fCD=39.84(0.17), fDE=-7.12(0.32), fAE=5.58(0.38).
- 1sqn: fAB=28.46(0.06), fBC=20.18(0.05), fCD=2.79(0.08), fDE=-17.37(0.38), fAE=-63.21(0.43).
- 1gpk: fAB=36.57(0.04), fBC=-142.54(0.05), fCD=-90.21(0.12), fDE=-54.08(1.18), fAE=-38.32(1.21).

Full:
- 1sqn: fAB=28.47(0.04), fBC=40.65(0.07), fCD=-17.42(0.11), fDE=23.67(0.25), fAE=-62.88(0.22).
- 1jla: fAB=62.90(0.08), fBC=-111.67(0.05), fCD=-177.93(0.17), fDE=79.81(0.14), fAE=-49.36(0.27).
- 1gpk: fAB=36.57(0.05), fBC=-77.85(0.05), fCD=-155.19(0.09), fDE=75.55(1.20), fAE=-38.36(1.16).

Consistency check: for the same system, Desolvated and Full fAE should agree within error (fAB nearly identical: 1jla 62.88 vs 62.90; 1sqn 28.46 vs 28.47; 1gpk 36.57 vs 36.57). fDE has opposite signs between pathways (e.g. 1gpk: -54.08 Desolvated vs +75.55 Full).

Note: the sum fAB+fBC+fCD+fDE does not exactly equal fAE in the table because BPMF uses fBC,L + f'_CD (eqs above) rather than fBC + fCD, and fAE is estimated separately. Use fAE as the reported total BPMF.
