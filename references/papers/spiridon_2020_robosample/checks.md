# Checks / fixtures - Robosample (Spiridon et al. 2020)

Concrete numbers an implementation can be regression-tested against.

## Constants / targets
- Optimal HMC acceptance rate: **0.651** (target for time-step tuning; refs [46,48]).
- Generalized-coordinate count: `N_f = 3N - m` (N atoms, m constraints).
- Direct mass-matrix-tensor inversion complexity: O(N^3); SOA recursive: O(N).
- Cluster analysis: Ward's method hierarchical agglomerative clustering, RMSD flattened at threshold **2 Å**.
- Free-energy histograms: **15-degree** (φ,ψ) spacing; alanine runs of **120,000** Monte Carlo moves.

## Alanine dipeptide - optimized HMC parameters (all-bond dynamics)
given move type, expect (L steps × ε fs = T fs):
- fully flexible: 108 × 1.87 fs = 202.0 fs
- torsional (TD): 25 × 17.1 fs = 427.5 fs
- ball (spherical): 158 × 2.1 fs = 331.8 fs
- cylindrical: 47 × 6.1 fs = 286.7 fs

## Alanine dipeptide - optimized HMC parameters (Ramachandran dynamics)
- torsional: 11 × 44.73 fs = 492.0 fs
- ball: 51 × 8.957 fs = 456.8 fs
- cylindrical: 49 × 12.386 fs = 606.9 fs
(Ball and cylinder use larger time steps than flexible but smaller than torsional; larger rigid bodies of Rama dynamics allow larger steps.)

## Alanine dipeptide conformational region boundaries (Table 3, degrees)
given region, expect (φ_min, φ_max, ψ_min, ψ_max):
- C5:   (-180, -95, 105, 180)
- PPII: (-96, -45, 105, 180)
- C7eq: (-96, -45, -25, 104)
- αL:   (35, 85, -180, 25)
<!-- CHECK: Table 3 header prints "ψmax" twice; interpreted 3rd column as ψ_min. -->
- Simulations at **625 K** reproduce basin positions of prior work at 625 K (ref [25]) and 800 K (ref [24]).

## Alanine dipeptide MFPT (Table 4, in units of MD steps / 200; mean ± std over 3 runs)
Row = from state, column = to state (C5, PP2, C7eq, αL).

Fully flexible:
- C5:   3.2±0.2,  9.1±0.5,  14.3±1.7, 5380.0±285.0
- PP2:  9.6±0.9,  5.7±0.4,  10.9±1.8, 5380.0±285.0
- C7eq: 16.5±1.1, 12.6±0.6, 2.8±0.2,  5370.0±285.0
- αL:   502.0±195.0, 499.0±194.0, 491.0±194.0, 14.1±4.4

Mixed-RamaTD:
- C5:   1.8±0.03, 4.0±0.1, 2.8±0.1, 631.0±92.2
- PP2:  2.5±0.1,  3.2±0.1, 2.6±0.1, 631.0±92.2
- C7eq: 3.2±0.2,  4.7±0.2, 1.5±0.02, 630.0±92.2
- αL:   50.9±1.1, 52.4±1.2, 49.4±1.1, 8.3±1.3

Mixed-RamaCyl:
- C5:   2.3±0.03, 5.2±0.1, 3.5±0.02, 615.0±33.6
- PP2:  3.2±0.1,  4.1±0.1, 3.3±0.005, 615.0±33.4
- C7eq: 4.1±0.1,  5.9±0.1, 1.9±0.01, 614.0±33.6
- αL:   46.4±4.5, 48.4±4.6, 44.7±4.7, 12.0±1.4

Mixed-RamaBall:
- C5:   2.4±0.01, 5.2±0.01, 3.1±0.1, 518.0±18.6
- PP2:  3.1±0.01, 4.1±0.04, 3.0±0.1, 518.0±18.5
- C7eq: 3.7±0.01, 5.6±0.02, 1.9±0.02, 518.0±18.4
- αL:   37.0±1.4, 39.1±1.5, 35.8±1.6, 12.8±1.0

Interpretation: αL is the rare basin; crossing φ=0 barrier accelerated ~10x by Rama dynamics vs fully flexible (5380 -> ~518-631). Joint efficiency order: Torsion < Cylinder < Ball (Ball shortest MFPT for the rarest transition).

## E2N1 glycan - simulation setup
- Temperature: **300 K**; minimized to 1 kJ/mol (OpenMM v7.4); ff14SB (peptide) + GLYCAM06 (glycan).
- 4 worlds, HMC (trajectory length × step):
  - world 1 (fully flexible): 108 × 1 fs
  - world 2 (rigid monosaccharide rings): 25 × 8 fs
  - world 3 (mixed spherical+cylindrical on polysaccharide): 158 × 0.5 fs
  - world 4 (torsional at Asn anchor): 11 × 20 fs
- Each simulation repeated 3×.

## E2N1 glycan - timing benchmarks (Intel i9-9900K 3.60GHz, 2× RTX 2080Ti, 64 GB, Ubuntu 18.04)
given regime (DOFs), expect seconds/step:
- flexible MD (747 DOF): ~0.022 s/step
- RB-MD (105 DOF): ~0.008 s/step
- CDHMC world 2 (74 DOF): ~0.007 s/step
- CDHMC world 3: identical to RB-MD (~0.008 s/step)
- CDHMC world 4 (4 DOF): ~0.001 s/step
- Multi-world HMC reaches RMSD up to ~30 Å vs initial; fully-flexible MD/HMC stays ≤ ~15 Å.

## Alanine dipeptide model definitions
- Input: AmberTools16, ff14SB force field.
- Model (a) "all-bond": 7 mobile joints (all bonds except 2 terminal).
- Model (b) "Ramachandran": 2 joints (N-Cα and Cα-C bonds).
- Cycle structure: 1 fully-flexible world followed by 10 CDHMC worlds (TD/Cyl/Ball).
- Trial tuning: fixed step 2 fs (1 fs for flexible world), increase T then increase ε to hit ~0.651 acceptance.
