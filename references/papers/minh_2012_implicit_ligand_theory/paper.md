# Implicit Ligand Theory: Rigorous Binding Free Energies and Thermodynamic Expectations from Molecular Docking

David D. L. Minh (Department of Chemistry, Duke University), 2012 (arXiv:1208.4885).

## Abstract

A rigorous formalism is derived for estimating noncovalent binding free energies and
thermodynamic expectations from calculations in which receptor configurations are sampled
independently from the ligand. Because of this separation, receptor configurations only
need to be sampled once, facilitating binding free energy calculations in virtual screening.
Demonstrative calculations on a host-guest system (Cucurbit[7]uril) yield good agreement with
prior free energy calculations and isothermal titration calorimetry. Implicit ligand theory
guides improvements to molecular docking algorithms and clarifies induced fit vs conformational
selection in noncovalent macromolecular recognition.

## Introduction

Molecular docking predicts the most stable configuration of a noncovalent complex and assigns
a score to rank binding affinity. Current scoring functions predict binding free energies poorly,
so docking is used to filter large libraries; false positives/negatives are common. Improvement
has been hindered by the lack of a rigorous formalism for obtaining binding free energies from
docking to *rigid* receptor structures (existing formalisms require a flexible receptor). This
paper derives implicit ligand theory for estimating binding free energies and thermodynamic
expectations by docking ligands to rigid receptor structures, describes statistical estimation,
presents example calculations, and discusses how physics-based docking algorithms can exploit it.

## Theory

The standard binding free energy for R + L ⇋ RL is given by eq:1 (from concentrations) and,
via statistical thermodynamics, by a ratio of configurational partition functions, eq:2, with
the partition functions defined by eqs 3-5. The potential energy `U(r_X, r_S)` depends on the
internal coordinates `r_X` of receptor/ligand/complex (external DOF analytically integrated)
and on solvent coordinates `r_S`. The complex coordinates `r_{RL}` decompose into receptor
internal coords `r_R`, ligand internal coords `r_L`, and six external DOF `ξ_L` (relative
translation + rotation). The indicator function `I_ξ ≡ I(ξ_L) ∈ [0,1]` determines whether R and
L are complexed; for tight-binding complexes the binding free energy is insensitive to the
precise definition of `I_ξ`. Jacobians for the Cartesian→internal/external transformation are
not shown.

### Implicit Solvent Theory

Integrating the partition functions over solvent gives the implicit-solvent integrals eqs 6-7,
with the solvation PMF `W(r_X)` defined by eq:8. `W(r_X)` is the constant-pressure reversible
work of transferring species X from gas phase into solvent, typically approximated as a
Poisson-Boltzmann (or Generalized Born) electrostatic term plus a nonpolar term ∝ molecular
surface area. In terms of implicit-solvent integrals the binding free energy is eq:9. Implicit
solvent models fail to capture specific interactions (e.g. hydrogen bonding), so are generally
less accurate than explicit solvent, but have yielded promising agreement with experiment.

### Implicit Ligand Theory

Define the effective potential `𝒰(r_X) = U(r_X) + W(r_X)`, the effective interaction energy
`Ψ(r_{RL}) = 𝒰(r_{RL}) − 𝒰(r_R) − 𝒰(r_L)`, and the **binding PMF** `B(r_R)` by eq:10 - an
exponential average of `Ψ` over ligand internal + external coordinates at a fixed (rigid)
receptor. The angled-bracket notation `⟨...⟩_{X,...}^r` denotes an ensemble average over the
superscript coordinates with respect to density ∝ `q_{X,...}`; here
`q_{L,I}(r_L, ξ_L) = I_ξ e^{−β𝒰(r_L)}`. In terms of the binding PMF, the binding free energy
becomes eq:11, where `Ω = ∫ I_ξ dξ_L` is the binding-site volume,
`ΔG_ξ = −β⁻¹ ln(ΩC°/8π²)` is the free energy of confining the ligand external DOF to the
binding site, and `q_R(r_R) = e^{−β𝒰(r_R)}`. **Eqs 10 and 11 are the central theoretical results.**

Implicit ligand theory separates the sampling of receptor and ligand configurations: the receptor
probability density is independent of any ligand configuration, and the ligand-internal density is
independent of the receptor. The primary benefit is that the expensive step of sampling receptor
configurations only needs to be done once; predicting binding free energies for a chemical library
is then limited by the faster process of sampling ligand conformations. (In practice, sampling from
the noninteracting ligand distribution may converge slowly.)

### Thermodynamic Expectations

Define the interaction-weighted rigid-receptor expectation `Θ(r_R)` of an observable `O(r_{RL})`
by eq:12. If `O` depends only on the receptor, `Θ(r_R)` reduces to `O(r_R) e^{−βB(r_R)}`. The
bound-ensemble expectation of `O` (w.r.t. `q_{RL,I}(r_{RL}) = I_ξ e^{−β𝒰(r_{RL})}`) is eq:13.
Eqs 12-13 generalize implicit ligand sampling (Cohen et al.), which is recovered by choosing `O`
as a Dirac delta on the ligand center of mass, taking a log, and multiplying by `β⁻¹`.

## Estimation

Applying implicit ligand theory to predict binding free energies has three steps:
1. Sample receptor configurations.
2. Estimate the binding PMF `B(r_R)` for each receptor configuration.
3. Estimate `ΔG°` from the `B(r_R)` estimates.

### Receptor Configurations

Receptor configs may be drawn from `q_R(r_R)`, from any (possibly unnormalized) distribution
`q_{R,w}(r_R)` on the same support with computable `w(r_R) = q_R(r_R)/q_{R,w}(r_R)`, or from
multiple such distributions. Convergence requires representative sampling of both bound and
unbound receptor space. A simple protocol: run MD in the implicit solvent used for `W(r_R)`,
collecting snapshots at intervals longer than the correlation time. For larger fluctuations,
apply biasing potentials on order parameters; if the ligand significantly perturbs the receptor
ensemble, introduce alchemical intermediates (coupling parameter `λ`: noninteracting at λ=0,
fully interacting at λ=1) and enhance sampling with Hamiltonian replica exchange. Receptor
configs from a flexible-receptor HREX with one ligand may be reused for other ligands.

Caveat: implicit ligand theory does not justify docking to multiple experimentally-determined
structures or any set for which `w(r_R)` is unknown (e.g. homology modeling, flexible docking).

### Estimating a Binding PMF

The binding PMF is a rigid-receptor free energy difference (eq:14) and can be computed by FEP,
TI, or BAR. The most direct estimator is forward FEP, eq:15, sampling ligand coords from
`q_L(r_L, ξ_L) = e^{−β𝒰(r_L)}` and resampling `ξ_L` from `q_{ξ,I} = I_ξ`. In exponential
averages a small subset of samples may dominate; the limiting case is the **dominant state
approximation**, using a single value of `Ψ(r_{RL})` to estimate `B(r_R)`. A fourth-order
cumulant expansion (eq:16) is an alternative, with `δΨ = Ψ − ⟨Ψ⟩`.

If most ligand poses overlap atoms and have high `Ψ`, FEP converges slowly. One remedy is to
bias the external DOF with a confining potential `U_c(ξ_L)`, giving eq:17 with
`Ω_c = ∫ I_ξ e^{−βU_c(ξ_L)} dξ_L`. Good `U_c` (ascertainable from existing docking algorithms)
favor low-`Ψ` poses. Alternatively use the inverse (reverse-FEP) form eq:18, sampling from the
fully-interacting rigid-receptor complex; this is problematic due to the rarity of
separated/overlapping configs, so a confined reference state (eq:19) is preferred. Phase-space
overlap is generally resolved via multiple alchemical stages and HREX, summing adjacent-stage
FEP/TI/BAR estimates or using MBAR.

### Estimating the Binding Free Energy

If receptor configs are drawn from `q_R`, `ΔG°` is estimated by the sample mean eq:20 (the
dominant state approximation and cumulant expansion apply as for eq:10). If drawn from a biased
distribution, use the importance sampling identity eq:21 with `w(r) = q_T(r)/q_S(r)`, giving the
weighted estimator eq:22. Multiple biased distributions → MBAR.

### Thermodynamic Expectations (estimation)

Estimated from the same data. The estimator for `Θ(r_R)` depends on how ligand configs were
sampled (sample mean if from `q_{ξ,I}`, otherwise importance sampling / MBAR); then the estimator
for eq:13 depends on how receptor configs were sampled.

## Demonstration

Implicit ligand theory was applied to estimate `ΔG°` of various ligands (adamantanes,
bicyclooctanes, ferrocenes) to Cucurbit[7]uril (CB[7]) in water, benchmarked against isothermal
titration calorimetry (ITC) and second-generation mining minima (M2) calculations. Receptor
configs sampled by MD; binding PMFs by multi-stage alchemical calculation + MBAR; binding free
energy by eq:20 or the dominant state approximation. See `checks.md` for full hyperparameters.

Because alchemical coupling was performed in vacuum (NAMD cannot combine alchemical decoupling
with implicit solvent), the binding PMF was estimated via the decomposition eq:23:
`B(r_R) = B_cpl + B_RL − B_L − ΔU(r_R)`, where `B_cpl` (vacuum coupling FE) was estimated by MBAR
and `B_RL`, `B_L` (vacuum→target transfer FEs) by single-step FEP. This decomposition lets one
evaluate `B(r_R)` for many force fields from the same samples. Four force fields compared: NAMD
(GBSA), M2 (GBSA), PB (Poisson-Boltzmann from UHBD + M2 valence/coulomb/vdW), and PBSA (PB + M2
nonpolar surface-area). Rigid-receptor expectations were estimated by MBAR, eq:24. Accuracy
assessed by correlation coefficient R² and RMSE, eq:25.

### Results

- Binding PMF estimates are strongly force-field dependent (B11 changes ~40 kcal/mol across
  force fields); eq:23 entails computing a small difference between large values.
- With 2 ns total simulation, binding PMF std dev ranges 0.12-1.63 kcal/mol; means stabilize
  after ~0.75 ns. The slowest-converging component of eq:23 varies by ligand; `B(r_R)` and
  `min{Ψ(r_{RL})}` converge at similar rates, suggesting the limiting factor is finding the
  lowest-interaction-energy configuration.
- A single minimized receptor gives high correlation with experiment (R²=0.884 for NAMD) and M2
  (R²=0.827 for NAMD) but large RMSE (>10 kcal/mol). Using eq:23 with PBSA gives lower R² but
  lower RMSE.
- Using 100 receptor structures substantially improves estimates: R²_Exp=0.704/RMSE_Exp=4.5 and
  R²_Gilson=0.925/RMSE_Gilson=2.4 (PBSA). The average `ΔG°` stabilizes after ~15 receptor
  snapshots; further snapshots reduce variance until limited by binding-PMF variance.
- Mean potential-energy changes on complexation are consistent with M2 (Table III).

## Discussion

Good agreement with M2 is a proof of principle; convergence/accuracy will differ for
protein-ligand systems (representative receptor sampling and low-energy pose finding will need
more MD time, but many protein-ligand systems are less charged). The `B(r_R)` decomposition
(eq:23) makes it easy to integrate alternate/expensive potentials (e.g. QM, better nonpolar
solvation). Computations may be accelerated by MD packages that skip pairwise interactions
between rigid atoms and by implicit solvent models designed for rigid receptors.

Implicit ligand theory guides docking improvements: `Ψ(r_{RL})` gives a functional form
accounting for solvation and ligand strain frequently ignored by scoring functions; it shows how
to combine information from docking to multiple receptor snapshots (exponential average or
cumulant expansion of `B(r_R)` across snapshots, rather than taking the minimum). Current docking
packages rank by a single low-energy configuration - the crudest form of implicit ligand theory
(the dominant state approximation). The relaxed-complex cost may be reduced by clustering
snapshots and weighting cluster representatives by cluster size (assuming `B(r_R)` roughly
constant within a cluster).

Estimating `B(r_R)` by docking requires a paradigm shift from searching for a minimum to sampling
from a distribution. Matching algorithms (e.g. DOCK) can estimate `B(r_R)` by post-processing:
bias receptor-independent random ligand-orientation sampling with a confining potential `U_c(ξ_L)`
(eq:17; harmonic `U_c` ⇒ Gaussian orientation), or start rigid-receptor MD from the lowest-energy
match and use eq:18/eq:19. Docking-simulation methods (AutoDOCK, MCDOCK) must sample from a known
distribution; being Monte Carlo-based (often simulated annealing), they can be modified to compute
importance-sampling weights (Neal).

Implicit ligand theory also quantifies induced fit vs conformational selection. If the complex is
dominated by a single receptor structure `r_R*` with `B(r_R) = ∞` for all others, eq:11 simplifies
to `ΔG° = 𝒰(r_R*) + B(r_R*) + β⁻¹ ln Z_R + ΔG_ξ`. For conformational selection,
`p(r_R*) = e^{−β𝒰(r_R*)}/Z_R` is reasonably high; for induced fit, `𝒰(r_R*)` is less favorable and
`B(r_R*)` compensates to achieve the same `ΔG°`. Since all receptor configs have finite Boltzmann
probability, the distinction is a matter of degree.

Source code and data: https://simtk.org/home/implicit_ligand

## Supplemental: Hybrid Implicit-Explicit Solvent

A small number of explicit solvent molecules can be treated as part of the receptor during binding
PMF calculations. Separating solvent coords `r_S` into explicit `r_E` and implicit `r_I`, one
defines partition functions eqs 26-27, the effective interaction energy
`Ψ'(r_{RL}, r_E) = 𝒰(r_{RL}, r_E) − 𝒰(r_R, r_E) − 𝒰(r_L)`, and a binding PMF `B'(r_R, r_E)`, giving
the binding free energy eq:28 with `q_{R,E} = e^{−β𝒰(r_R,r_E)}`.

### Derivation (not implemented)

Eqs 26-28 are obtained by the same manipulation as the implicit-ligand derivation of eqs 6-11,
now carrying the explicit-solvent coordinates `r_E` alongside the receptor coordinates. The main
text focuses on implicit solvent with the understanding that explicit solvent may be readily
included.
