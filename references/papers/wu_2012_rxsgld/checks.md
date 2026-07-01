# Checks / fixtures - Wu 2012 RXSGLD

## Special-case reductions (unit tests for the factor machinery)

- Given `λ_i = 0` (no guiding), eq:1 reduces to Langevin dynamics; guiding force
  `g_i = 0`.
- Given Langevin dynamics: `λ_lf = 1`, `λ_hf = 1`, `χ_lf = 1`, `χ_hf = 1`. Then
  eq:7 -> eq:13 (`Θ_SGLD = Θ_LD`), the reweighting factor `w_SGLD = 1`, and the
  exchange coefficients `μ_tilde_m = 0`, `μ_m = β_m = 1/(kT_m)`. Then the RXSGLD
  exchange probability eq:20 must collapse to the standard TRXLD form eq:21.

## Skewed double-well test system (Sec. III.A, eq:27-30)

Fixed parameters: `a = 20000 kT_0`, `b = 160 kT_0`, `w = 2 Å`, base `T_0 = 50 K`.
Energy barrier between the two wells ~ `10 kT_0`. Single argon atom. Simulation:
8 stages, 8 replicas, collision frequency `100/ps`, timestep `1 fs`, length
`100 ns`, `t_L = 0.2 ps`. TRXLD at `T = 50/100 K`; RXSGLD at `T = 50 K`,
`T_SG = 50/100 K`. Analytic references from eq:28-30.

Analytic identity: `E_xz = kT` exactly (eq:29a); at base `T = 50 K` this is the
x,z contribution. Reported `Epot` in Table I is in units of `kT` (total = x,z + y).

Table I - base-stage (`T = 50 K`) properties. `x_1` = fraction in well near
`y = 0 Å`.

| skew s | method | Epot (kT) | x_1 |
|---|---|---|---|
| 0 | analytic solution | 1.525 | 0.5 |
| 0 | TRXLD | 1.513 ± 0.002 | 0.477 ± 0.006 |
| 0 | RXSGLD | 1.520 ± 0.002 | 0.491 ± 0.006 |
| kT | analytic solution | 1.796 | 0.738 |
| kT | TRXLD | 1.760 ± 0.010 | 0.763 ± 0.005 |
| kT | RXSGLD | 1.827 ± 0.009 | 0.705 ± 0.006 |
| 2kT | analytic solution | 1.762 | 0.878 |
| 2kT | TRXLD | 1.754 ± 0.012 | 0.881 ± 0.003 |
| 2kT | RXSGLD | 1.772 ± 0.011 | 0.876 ± 0.003 |

Use these to validate a reweighting/RXSGLD implementation: base-stage `Epot` and
`x_1` must match the analytic solution within the quoted error bars.

Note: `s = 0` gives symmetric wells so `x_1 = 0.5` exactly (analytic). The wells
are at `y = 0 Å` and `y = w = 2 Å`.

## β-hairpin peptide, implicit solvent (Sec. III.B)

9-residue peptide Tyr-Gln-Asn-Pro-Asp-Gly-Ser-Gln-Ala; SCPISM implicit solvent.
8-stage TRXLD (`T = 274/400 K`) and 8-stage RXSGLD (`T = 274 K`,
`T_SG = 274/400 K`), 200 ns each, from extended conformation, collision
frequency `1/ps`.

Subset-indexing clustering (SIC): 16 backbone dihedral subsets (Tyr(1) has no φ,
Ala(9) has no ψ), region counts `k_i = 1,2,2,2,2,1,2,2,2,2,2,2,3,2,2,1`.
- Total possible clusters `N_c = prod k_i = 12288`. (Check: product of the 16
  values = 12288.)
- Clusters visited by all replicas: total 1145; TRXLD 1056, RXSGLD 730.
- Base-stage clusters: TRXLD 244, RXSGLD 283.
- Conformational searching relevancy CSR = base/all:
  TRXLD 244/1056 = 23.1%; RXSGLD 283/730 = 38.8%.

## β-hairpin peptide, explicit water (Sec. III.C)

829 TIP3P waters + 1 Na+, box `30 × 30 × 30 Å`, CHARMM 22 force field, 3D IPS
(local region radius 10 Å), collision frequency `1/ps`. 8 stages, 20 ns each.
Exchange attempts every 1000 steps.

Average replica-exchange acceptance ratios:

| temperature range | TRXLD accept | RXSGLD accept |
|---|---|---|
| 274/310 K | 31.1% | 65.3% |
| 274/350 K | 6.4% | 63.5% |
| 274/400 K | 5.2% | 70.2% |

(RXSGLD `T_SG` ranges 274/310, 274/350, 274/400 K at fixed `T = 274 K`.)

Replica-0 time to first reach top stage (stage 7):
- RXSGLD: within 0.1 ns for all three ranges.
- TRXLD: 0.6 ns (274/310 K), 2.35 ns (274/350 K), 3.28 ns (274/400 K).
