# Checks / fixtures - Kandel et al. 2016

## Physical / mathematical invariants

- **Uniform torsion pdf (no torsional force):** given $U(\alpha)=0$ and Fixman
  potential applied to a constrained (TMD) chain, expect torsion-angle pdf
  $\rho_u(\alpha_i)=1/(2\pi)$ (uniform). (eq:14)
- **Fixman cancels intrinsic bias in config pdf:** for a separable potential,
  the Fixman-corrected constrained config pdf (eq:12) equals the unconstrained
  pdf (eq:6). RMS deviation of the Fixman torsion pdf from FLEXIBLE -> 0 and its
  dependence on barrier-peak location -> removed.
- **Fixman only PARTIALLY corrects barrier-crossing rates:** even with Fixman,
  $f_{TS,\text{fix}}$ retains $S^{-1}(\alpha_0)$ (eq:19), so the corrected rate
  still depends on barrier-peak location, unlike the unconstrained model.
- **Gaussian integral fixtures (Appendix):** $\int_{-\infty}^{\infty}
  e^{-x^*Ax/2}dx = [(2\pi)^p/\det\{A\}]^{1/2}$; $\int_0^\infty y\,e^{-sy^2/2}dy = 1/s$.
- **Schur determinant:** $\det\{\mathcal{M}\}=\det\{W_0\}\det\{S\}$ with
  $S=S_0-VW_0^{-1}V^*$ (eq:A5).

## C4 idealized serial chain (test system)

Fixed geometry / masses:
- bond angles = 90°, bond lengths = 1.54 Å, atom masses = 14.01 amu.
- barrier torsional potential eq:15 with $k_\alpha = 0.30$ kcal/mol.
- barrier peak locations tested: $\alpha_0 \in \{0°,45°,90°,135°,180°\}$.

Langevin dynamics settings (separable-potential runs):
- three 20 ns simulations, time step = 1 fs, damping constant = 0.1/fs.
- FLEXIBLE bond spring = 303.1 kcal/Å², angle spring = 63.21 kcal.
- example analysis point: $\alpha_0=90°$, $T=800$ K (Fig 2a torsional pdf).

Non-separable (Coulomb-coupled) C4 runs:
- Coulomb constant $k_{\text{coul}} = 332.06$ kcal·Å/e².
- terminal charges $= +0.2e$ and $-0.2e$.
- bond lengths 1.54 Å, bond angles 90°, masses 14.01 amu.
- angle spring $k_\theta$ swept; AMBER99SB realistic range = 30–100 kcal.
  Fig 3 uses $k_\theta = 95$ kcal.
- Expect: FIXMAN torsion pdf matches FLEXIBLE only for very large $k_\theta$
  (stiff angle); for realistic $k_\theta$ the extrinsic distortion persists.

## Alanine dipeptide (hybrid ICMD)

- Forcefield: AMBER99SB + GBSA implicit solvent, GneimoSim software.
- Backbone dihedrals: $\phi$ = C–N–Cα–C, $\psi$ = N–Cα–C–N; initial conformation
  alpha-helical $(\phi,\psi)=(-60°,-40°)$.
- 20 FLEXIBLE + 20 TMD sims each at 300 K and 800 K, 20 ns each.
- FIXMAN sims at 300 K; Langevin damping 0.1/fs, time step 2 fs, 20 sims × 20 ns.
- Result: opening backbone angles **C–N–Cα and N–Cα–C together** recovers the
  FLEXIBLE free-energy surface (samples 1st & 4th quadrants). Opening only one of
  the two, or additional angles, is insufficient / no extra benefit.

## Dipeptide series - required open angles (min Hellinger distance)

Opening C–N–Cα + N–Cα–C suffices for PHE, ALA, LEU, TRP.
For VAL, ILE, MET, TYR additionally open Cα–C–N.

| dipeptide | open backbone angles | Hellinger distance to FLEXIBLE |
|---|---|---|
| VAL | C–N–Cα, N–Cα–C, Cα–C–N | 0.14 |
| MET | C–N–Cα, N–Cα–C, Cα–C–N | 0.13 |
| TYR | C–N–Cα, N–Cα–C, Cα–C–N | 0.17 |
| ILE | C–N–Cα, N–Cα–C, Cα–C–N | 0.16 |
| PRO (FIXMAN, all angles fixed) | none | 0.45 |
| PRO (hybrid) | N–Cα–C + sidechain Cα–Cβ–Cγ | 0.26 |

- Opening only N–Cα–C alone gives pdfs far from FLEXIBLE (insufficient).

## CLN025 10-residue peptide

- NMR reference structure PDB ID **2RVD** (Honda et al. 2008).
- 8 sims × 10 ns each for FIXMAN / hybrid ICMD / FLEXIBLE.
- 300 K via Nose–Hoover thermostat, time step 1 fs.
- GBSA solvation: internal dielectric = 4.0, external dielectric = 78.0.
- hybrid ICMD open angles: N–Cα–C for proline; C–N–Cα–C and N–Cα–C for others.

### Table I - % of confs within 1 std of NMR mean H-bond distance

| Hydrogen bond | NMR mean (Å) | NMR std (Å) | % FIXMAN | % hybrid ICMD | % FLEXIBLE |
|---|---|---|---|---|---|
| Thr8.O–Tyr10.N  | 3.62 | 0.25 | 58.46 | 22.46 | 8.32  |
| Asp3.OD1–Thr6.N | 3.38 | 0.14 | 0.55  | 7.27  | 10.25 |
| Asp3.O–Thr8.N   | 2.89 | 0.16 | 9.78  | 65.51 | 57.60 |
| Asp3.O–Gly7.N   | 3.56 | 0.20 | 0.18  | 31.20 | 25.29 |
| Asp3.O–Thr6.N   | 5.58 | 1.14 | 0.19  | 45.98 | 43.62 |
| Asp3.O–Thr8.O   | 3.51 | 0.42 | 1.16  | 0.41  | 0.49  |

- Trend: hybrid ICMD % much closer to FLEXIBLE than FIXMAN for most H-bonds.

## Time step / stability

- Hybrid ICMD with key backbone angles open remains stable at **5 fs** time step
  (20 sims × 50 ns each) - much larger than all-atom Cartesian MD.
