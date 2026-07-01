# Replica exchanging self-guided Langevin dynamics for efficient and accurate conformational sampling

Xiongwu Wu, Milan Hodoscek, Bernard R. Brooks. J. Chem. Phys. 137, 044106 (2012). DOI: 10.1063/1.4737094

## Abstract

This work presents a replica exchanging self-guided Langevin dynamics (RXSGLD)
simulation method for efficient conformational searching and sampling. Unlike
temperature-based replica exchanging simulations, which use high temperatures to
accelerate conformational motion, this method uses self-guided Langevin dynamics
(SGLD) to enhance conformational searching without the need to elevate
temperatures. A RXSGLD simulation includes a series of SGLD simulations, with
simulation conditions differing in the guiding effect and/or temperature. These
simulation conditions are called stages and the base stage is one with no guiding
effect. Replicas of a simulation system are simulated at the stages and are
exchanged according to the replica exchanging probability derived from the SGLD
partition function. Because SGLD causes less perturbation on conformational
distribution than high temperatures, exchanges between SGLD stages have much
higher probabilities than those between different temperatures. Therefore, RXSGLD
simulations have higher conformational searching ability than temperature based
replica exchange simulations.

## I. Introduction

Conformational searching and sampling is a fundamental process of molecular
systems. Raising temperatures can accelerate thermal motions (e.g. simulated
annealing, temperature-based replica exchange), but causes changes in
conformational distribution, often leading to protein unfolding and phase
transition.

Unlike high temperature simulations that accelerate all thermal motions, the
self-guided Langevin dynamics (SGLD) enhances only the low frequency motion that
is the most important for conformational searching and sampling. With a simple
local averaging scheme, SGLD selectively enhances molecular motions based on their
frequencies without modifying energy surfaces or raising temperatures.

The concept of SGLD: in normal (Langevin) dynamics, kinetic energy distributes
evenly among all degrees of freedom, i.e. `kT/2` per degree of freedom. A molecule
has high frequency motions (bond vibration, bond bending) and low frequency
motions (bond rotations about the `φ`, `ψ` dihedral angles). The low frequency
motions are the limiting steps for conformational searching. In SGLD, the guiding
forces enhance the low frequency motions (as defined by the local averaging time
`t_L`) while suppressing high frequency motions to maintain the overall
temperature. The low frequency motions gain more kinetic energy and become hotter;
the high frequency motions lose kinetic energy and become cooler; the overall
temperature remains unchanged.

Note (Fig. 1 semantics): In LD, all motions have kinetic energies equivalent to a
temperature `T`; in SGLD, low frequency motion gains more kinetic energy (hotter)
and high frequency motion loses kinetic energy (cooler), enhancing conformational
searching without raising temperature.

Under the effect of the guiding forces, SGLD has its own conformational
distribution, defined as the SGLD ensemble. The partition function of the SGLD
ensemble allows a conversion between the SGLD ensemble and the canonical ensemble
so that canonical ensemble properties can be calculated from SGLD simulations
through reweighting. Based on the SGLD partition function, the force-momentum based
SGLD (SGLDfp) method directly samples the canonical ensemble without reweighting.

Goals of this work: (1) describe the RXSGLD method and derive the exchange
probability from the SGLD partition function, (2) examine how accurately RXSGLD
samples the canonical ensemble, (3) evaluate the conformational searching ability
of RXSGLD vs a temperature-based replica exchange simulation method (TRXLD).

## II. Theory and methods

The RXSGLD simulation method is a combination of the replica exchange approach and
the SGLD method. By incorporating the SGLD partition function into the exchange
probability, RXSGLD maintains a canonical ensemble at its base stage while
achieving enhanced conformational searching and sampling.

### A. SGLD simulation method

For any particle `i`, the equation of the self-guided motion has the general form
(eq:1). Here `p_dot_i` is the time derivative of momentum, `f_i` is the interaction
force, and `R_i` is a random force related to mass `m_i`, collision frequency
`γ_i`, and temperature `T` by eq:2.

Equation (1) contains a guiding force `g_i` (eq:3) computed from the momentum `p_i`
and the low frequency momentum `p_tilde_i`. Here `λ_i` is the guiding factor
defining the strength of the guiding force; when `λ_i = 0`, eq:1 reduces to Langevin
dynamics. The parameter `ξ` is an energy conservation factor to eliminate any net
energy input from the guiding force (eq:4, eq:5). The guiding force does not cause
energy flow between the system and its environment; instead it causes energy flow
between different motion modes within the system.

The low frequency portion of any property `P` (denoted with a "~" cap, `P_tilde`)
is calculated as a progressive local average (eq:6). This is memory efficient: an
update with current instantaneous values. The local averaging acts as a low
frequency filter. `P - P_tilde` is the high frequency property. Every motion mode
contributes to both the low and high frequency properties, but the proportion
depends on its frequency.

The guiding-force effects are summarized by `λ_lf`, `λ_hf` (bias effects on the low
and high frequency energy surfaces) and `χ_lf`, `χ_hf` (effects on the low and high
frequency motions). Together, the SGLD ensemble has a configurational partition
function of the form eq:7. The energy factors `λ_lf`, `λ_hf` are average projections
of total forces along the interaction forces (eq:8, eq:9). The collision factors
`χ_lf`, `χ_hf` are computed from projections of the guiding forces along the friction
forces (eq:10, eq:11). `T_tilde` is the low frequency temperature (eq:12), and
`T_tilde_0` is the reference low frequency temperature (the value when all guiding
factors are zero).

In a Langevin dynamics (LD) simulation, `λ_lf = 1`, `λ_hf = 1`, `χ_lf = 1`,
`χ_hf = 1`, so eq:7 reduces to the canonical partition function (eq:13). The
canonical partition function `Θ_LD` relates to the SGLD partition function
`Θ_SGLD` via the SGLD reweighting factor `w_SGLD` (eq:14). Any canonical ensemble
average `<P>_LD` can be recovered from an SGLD simulation by reweighting (eq:15).

The self-guiding temperature `T_SG` (eq:16) provides a rough measure of the
conformational searching ability in units of temperature: an SGLD simulation with
self-guiding temperature `T_SG` has conformational search ability comparable to a
high temperature simulation at `T = T_SG`. One can either set the guiding factors
`{λ_i}` directly, or set a target self-guiding temperature `T_SG^0` and adjust
`{λ_i}` so `T_SG` approaches `T_SG^0`.

### B. Replica exchanging self-guiding Langevin dynamics simulation

A replica exchange simulation consists of parallel simulations of identical systems
(replicas) at different simulation conditions (stages). There are `k+1` stages. The
stage with the conditions of interest is stage 0 (base stage), with
`T_SG^(0) = T` and `T^(0) = T`. The other `k` stages have different guiding
temperatures `T_SG^(i)` and the same or different temperatures `T^(i)`. The stage
with the maximum condition is the top stage (`T^(k)`, `T_SG^(k)`).

To give all neighboring-stage exchanges similar acceptance ratios, temperatures are
exponentially distributed (eq:17a, eq:17b). A TRXLD simulation has different stage
temperatures `T^(i) ≥ T^(0)` but no guiding force (`T_SG^(i) = T^(0)`); a RXSGLD
simulation has different self-guiding temperatures `T_SG^(i) ≥ T^(0)`. In this work,
RXSGLD keeps temperatures constant, `T^(i) = T^(0)`.

The key quantity is the exchange probability. According to the SGLD partition
function (eq:7), the distribution probability of conformation `i` at stage `m`,
`ρ_SGLD(X_m^(i))`, is eq:18, with parameters `μ_tilde_m` (eq:19a) and `μ_m`
(eq:19b). When two replicas exchange, the exchange probability `π_RX` takes the
form eq:20. Here the approximation is that the low frequency energies at different
stages are the same for the same conformation:
`E_tilde_p(X_m^[j]) ≈ E_tilde_p(X_n^[j])` and
`E_tilde_p(X_n^[i]) ≈ E_tilde_p(X_m^[i])`; this is accurate if `T_m = T_n` (the
recommended condition for RXSGLD).

For TRXLD (`λ_lf = λ_hf = χ_lf = χ_hf = 1`), `μ_tilde_m = μ_tilde_n = 0` and
`μ_m = β_m = 1/(kT_m)`, so the exchange probability takes the well-known form
(eq:21). To evaluate the RXSGLD exchange probability, one needs the low frequency
exchange coefficient `μ_tilde_m = β_m(λ_lf^(m) χ_lf^(m) - λ_hf^(m) χ_hf^(m))` and
the high frequency exchange coefficient `μ_m = β_m λ_hf^(m) χ_hf^(m)`, which need
`λ_lf`, `λ_hf`, `χ_lf`, `χ_hf` at each stage. These can be computed from individual
SGLD presimulations, or (like in SGLDfp) estimated during the simulation as
evolving averages (eq:22), with estimation time `t_est` typically `10 t_L`.

At each exchange interval, the exchange probability between a pair of neighboring
stages `π_RX({X_m^[i], X_{m+1}^[j]} -> {X_m^[j], X_{m+1}^[i]})` is computed via
eq:20 and accepted by the Metropolis criterion
`min{1, π_RX(...)}`. Here `m` is alternately odd and even stage numbers. Once an
exchange is accepted, the momentum `p_i` is scaled by a temperature-scaling factor
(eq:23a, eq:23b, eq:24), and the low frequency momentum `p_tilde_i` is scaled by a
low frequency temperature-scaling factor (eq:25a, eq:25b, eq:26). Between exchanges,
standard SGLD simulations are performed at all stages.

### C. Simulation details

RXSGLD is implemented in CHARMM version c36 and AMBER version 12; results here are
from CHARMM. Eight stages are used for all replica exchange simulations. Exchanges
are attempted every 1000 time steps. All RXSGLD simulations use a local average time
of 0.2 ps. All RXSGLD simulations have `T_SG^(i) > T_SG^(0) = T^(0)` and
`T^(i) = T^(0)`; all TRXLD simulations have `T_SG^(i) = T^(0)` and
`T^(i) > T^(0)`. `T_SG^(i)` and `T^(i)` at stage `i` are set by eq:17a and eq:17b.

## III. Results and discussions

Only RXSGLD vs TRXLD is compared.

### A. The skewed double well system

The simplest system with an energy barrier. One particle moving on a skewed double
well energy surface (eq:27). Parameter `a` defines the energy surface in `x` and `z`;
`b` and `w` define a double well in `y` (wells at `y=0` and `y=w`); skew `s` defines
the energy difference between the wells. With three degrees of freedom, the partition
function separates per degree of freedom (eq:28a, eq:28b). Ensemble averages
(eq:29a, eq:29b) and distributions (eq:30a, eq:30b) follow, with
`r_xz = sqrt(x^2 + z^2)`.

Parameters: `a = 20000 kT_0`, `b = 160 kT_0`, `w = 2 Å`. Three skew parameters
`s = 0, kT_0, 2kT_0`. Base temperature `T_0 = 50 K`. Energy barrier between the two
wells is about `10 kT_0`. An argon atom was simulated. 8 stages and 8 replicas per
replica exchange simulation. TRXLD at `T = 50/100 K`; RXSGLD at `T = 50 K`,
`T_SG = 50/100 K`. Collision frequency `100/ps`. Time step 1 fs, length 100 ns.
Local averaging time `t_L = 0.2 ps`.

At the base stages both methods produce correct ensemble average energies (Table I).
RXSGLD has smaller energy differences across stages than TRXLD because SGLD only
enhances low frequency motions. In the `y` dimension (low frequency), guiding forces
and raising temperatures have similar effects; in the `x-z` dimensions (high
frequency), RXSGLD distributions are nearly unchanged across stages while TRXLD
distributions deviate. RXSGLD approaches the analytic solutions faster than TRXLD.

### B. The β-hairpin folding peptide with implicit solvent

A 9-residue β-hairpin folding peptide (Blanco et al.), modified from the β-hairpin
of α-amylase inhibitor tendamistat (residues 15-23). Sequence:
Tyr(1)-Gln(2)-Asn(3)-Pro(4)-Asp(5)-Gly(6)-Ser(7)-Gln(8)-Ala(9). Screened Coulomb
potential implicit solvent model (SCPISM). An 8-stage TRXLD (`T = 274/400 K`) and an
8-stage RXSGLD (`T = 274 K`, `T_SG = 274/400 K`), started from a fully extended
conformation, 200 ns each. Collision frequency `1/ps`.

#### Subset indexing clustering (SIC) method

(1) Subsets: separate the conformational variables into subsets (eq:32-subsets).
The variables are chosen by research interest (dihedral angles, hydrogen bonds,
secondary structures).

(2) Regions: define regions in the distribution of each subset variable. For subset
`s_i` with `k_i` regions `{R_i(1), ..., R_i(k_i)}`, the region index is eq:31.

(3) Clusters: a cluster is a unique list of region indexes of all subsets
`{I_1, I_2, ..., I_m}`. Total possible cluster count `N_c = prod_{i=1}^m k_i`. A
conformation's SIC is `SIC(x_1,...,x_n) = {I_1, I_2, ..., I_m}`. The method does not
evaluate pairwise properties, so its cost is order `N`.

For this peptide, the 16 `φ`, `ψ` dihedral angles (Tyr(1) has no `φ`, Ala(9) has no
`ψ`) are used as subsets, with region counts
`k_i = 1,2,2,2,2,1,2,2,2,2,2,2,3,2,2,1`, giving `N_c = 12288` possible clusters.
Total clusters visited by all replicas: 1145 (TRXLD 1056, RXSGLD 730). At the base
stage: TRXLD 244, RXSGLD 283. Conformational searching relevancy (CSR) = base-stage
clusters / all-stage clusters: TRXLD 244/1056 = 23.1%, RXSGLD 283/730 = 38.8%.
RXSGLD concentrates the search on high-population regions.

### C. The β-hairpin folding peptide in aqueous solutions

Temperature-based replica exchange is not size extensive: a large system gives a
large energy change for the same temperature difference, reducing exchange
probability exponentially (eq:21). The peptide was dissolved in a box of 829 TIP3P
water plus one sodium ion to neutralize; box size `30 × 30 × 30 Å`. Collision
frequency `1/ps`. CHARMM 22 force field; 3D IPS method with local region radius 10 Å
for electrostatic and Lennard-Jones energy. Three 8-stage TRXLD simulations
(`T = 274/310 K, 274/350 K, 274/400 K`) and three 8-stage RXSGLD simulations
(`T = 274 K`, `T_SG = 274/310 K, 274/350 K, 274/400 K`). All started from a fully
extended conformation, 20 ns each.

Acceptance ratios: TRXLD averages 31.1% (`274/310 K`), 6.4% (`274/350 K`), 5.2%
(`274/400 K`). RXSGLD averages 65.3%, 63.5%, 70.2% for the three ranges. Replica
diffusion: for RXSGLD, replica 0 reached stage 7 within 0.1 ns in all three cases;
for TRXLD, 0.6 ns (`274/310 K`), 2.35 ns (`274/350 K`), 3.28 ns (`274/400 K`). In
TRXLD the stage energy distributions barely overlap (low exchange probability); in
RXSGLD they overlap significantly (high acceptance). RXSGLD searches more clusters at
the base stage faster than TRXLD.

## IV. Conclusions

RXSGLD uses SGLD to enhance conformational searching and has high replica exchange
efficiency. By avoiding temperature elevation, it applies to large systems with high
replica exchange efficiency using relatively few replicas. By incorporating the SGLD
partition function into the exchange probability, it samples the canonical ensemble
distribution at the base stage. Thus RXSGLD is an alternative to SGLDfp for directly
sampling the canonical ensemble without reweighting. RXSGLD has better size
extensiveness than TRXLD.
