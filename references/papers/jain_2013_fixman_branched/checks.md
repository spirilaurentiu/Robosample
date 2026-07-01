# Checks / fixtures - Jain et al. 2013, Fixman potential for branched molecules

Regression fixtures an implementation of the GNEIMO-Fixman potential/torque can be tested against.

## Analytical ground truth: C4 mass-matrix determinant (Pear-Weiner)

- System: idealized 3-bond serial chain (C4), fixed bond lengths, all bond angles = 90 deg, single torsion angle alpha.
- Expect (eq:24): `det{M(alpha)} = c_5 * (35 + 4 cos(alpha) - 16 cos^2(alpha) + cos^4(alpha))`, c_5 = const(bond lengths, masses).
- Fixman potential from eq:17 (GNEIMO-Fixman, U_f = c_f + 1/2 ln det{M}) must match, up to additive constant, U_f = c_f + 1/2 ln(det{M(alpha)}) using the eq:24 polynomial. Reported: "excellent agreement" (Fig 1a).
- Fixman torque from eq:23 must equal the numerical derivative -dU_f/dalpha of eq:17 potential: reported "identical" for C4 (Fig 1b).

## C4 torsion-angle distribution (serial chain)

- Bead mass 14 amu, bond length 1.54 A, bond angle 90 deg. Bond/angle spring constants: 83.66 kcal/A^2 (bond), 43.46 kcal (angle). Timestep 1 fs, 50 ns total, sample every 100 steps. Langevin, 300 K, gamma = 0.01/fs.
- FLEXIBLE: torsion pdf uniform, expect rho(alpha_i) = 1/(2pi) (eq:8).
- TORSIONAL (constrained, no Fixman): bimodal biased pdf; maxima at approx +/-83 deg, minima at 0 deg and +/-180 deg. Shape given by eq:25.
- FIXMAN (constrained + Fixman): recovers uniform pdf 1/(2pi).
- Free energy: G(alpha_i)/kT = -ln(rho(alpha_i)); uniform free-energy level = -ln(1/360) (using 7.2 deg bins => 360/7.2 = 50 bins... note: printed as -ln(1/360)). The negative Fixman potential -U_f nearly equals the TORSIONAL free-energy profile (equal magnitude, opposite sign).

## C5, C11, C15 serial chains

- C5: four-bond chain, 109 deg bond angle. Fixman potential minimum -> self-intersecting planar conformation (each torsion = 0 deg); maximum -> self-intersecting spatial conformation (each torsion approx +/-83 deg). Contour matches Patriciu et al.
- C11: ten-bond chain; same min/max conformation characterization as C5.
- C15: fourteen-bond chain; FIXMAN recovers uniform pdf for all torsions where TORSIONAL is bimodal.
- All: histogram bin dalpha = 7.2 deg.

## Branched peptides (AMBER99SB force field, Langevin, 300 K, gamma 0.01/fs)

Clustering: rigid clusters = groups with frozen bonds/angles; hinges = torsions. Terminal bonds rigid; aromatic rings rigid. Alanine dipeptide -> 8 clusters; valine dipeptide -> 10 clusters. Timestep 1 fs; 3 runs x 100 ns = 300 ns; sample every 200 steps.

RMS deviation of torsion pdf from uniform (smaller = closer to uniform):

| system | R_tor (TORSIONAL) | R_Fix (FIXMAN) | R_flex (FLEXIBLE) |
|---|---|---|---|
| alanine dipeptide | 1.0e-3 | 8.7e-5 | 6.5e-5 |
| valine dipeptide | 8.4e-4 | 7.9e-5 | 6.6e-5 |
| chignolin (10-aa) | 1.2e-3 | 8.6e-5 | (n/a) |

- Interpretation fixture: R_Fix should drop by ~order of magnitude vs R_tor and land close to R_flex (Fixman correction recovers flexible-model uniform pdf).
- Chignolin: PDB ID 1UAO (NMR structure); 4 runs x 100 ns; torsions validated: C1-N2-Ca2-C2, C4-N5-Ca5-C5, C7-N8-Ca8-C8. TORSIONAL bimodal -> FIXMAN nearly flat.

## Computational cost

- Adding GNEIMO-Fixman to GNEIMO constrained dynamics: +24% cost; scales linearly with number of clusters.
- With full all-atom force field: a GNEIMO-Fixman timestep costs approx 2.24x a FLEXIBLE timestep (vs 2.0x for plain TORSIONAL). Cost dominated by the recursive Upsilon(i) evaluation per cluster.
