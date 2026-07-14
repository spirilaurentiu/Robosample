# Spec: Transition-aware Gibbs block selection for rare conformational transitions

Status: **PARTIAL / draft** — researcher artifact, ready for spec-review. Every grounded Claim (C1-C3,
C5-C7) rests on authority (`references/` + codebase); the adaptive-selection algebra (C4a) and the
containment conditions (C4b) await a CAS proof and ingestion before Objective 2 may ship. Derived from a
full read of the authority library (`references/papers/spiridon_2017_cdhmc_gibbs`,
`spiridon_2020_robosample`, `chodera_2011_gibbs_replica_exchange`, `duane_1987_hmc`,
`shirts_2012_ensemble_validation`, `wu_2011_sgld`, `wu_2012_rxsgld`, `wu_2016_sgld_gle`), the current
blocking pipeline (`python/robosample/autoblock.py`, `batstat.py`, `context.py`), the world/factorization
interface (`include/Context.hpp`, `include/World.hpp`, `include/GibbsSweep.hpp`), the refactor design docs
(`docs/specs/refactor/DOC-GibbsSweep.md`, `DOC-NMA.md`, `DOC-VelocityDistortion.md`), the Stage-0
orientation (`docs/RESEARCH-ROADMAP.md` Thread 1), and the sibling catalog
(`docs/specs/torsional-dynamics-efficiency-and-crossing.md`).

Requirement keywords follow RFC 2119 (`styles/spec.md`): SHALL / SHOULD / MAY / NOTE.

Tooling caveat (reduced resolution): the detailed-balance algebra was NOT machine-checked (no CAS
available this pass). Every symbolic invariance step is marked **LEMMA** and SHALL be re-derived with a CAS
and closed as a proof-grade oracle before the adaptive route (Objective 2) is implemented. The recommended
route (Objective 1) depends on no un-verified algebra.

Revision (2026-07-14, post spec-review): no Blocking findings. SF-1 applied — V2 now names the
absolute-fixture arm as load-bearing and routes shared-path measure perturbations to V1 (INV-1), not the
lag-0-vs-time-lagged diff. SF-2 — Paper 2 non-transfer control now requires per-pair, same-modality/
sample-size/lag/resolution noise floors with the test `VI(f,g) > max(floor_f, floor_g)`. The lag-0 code
diagnosis, the C1/C2/C3 grounding, and the OQ-3 Fixman determinant-sign contradiction were all
independently confirmed.

Scope: this spec decides **which objective defines a Gibbs block** so blocks help sample rare
conformational transitions rather than only within-basin mixing. It changes the block-**construction**
criterion and the characterization experiment; it does **not** change the HMC kernel, Fixman correction,
acceptance rule, or momentum draw on its recommended path. Thread-1 companion to
`docs/specs/torsional-dynamics-efficiency-and-crossing.md` (the *dynamics* levers); the two share one
invariant — every world SHALL be exactly `π`-invariant — and this spec reuses, not restates, that
catalog's per-world oracles.

## Motivation (Problem restatement)

### User's original phrasing

> Produce an implementation spec for Robosample Thread 1: selecting Gibbs blocks (subsets of
> generalized/BAT coordinates sampled together, conditioning on the rest) that help sample RARE
> CONFORMATIONAL TRANSITIONS, not just within-basin mixing. Diagnosis to encode: the current approach
> builds an equal-time correlation/covariance matrix from a within-basin trajectory and partitions
> coordinates by it; that objective encodes only within-basin co-fluctuation and cannot see
> transition-relevant coupling. Evaluate and choose among three objectives — time-lagged (tICA/VAMP),
> adaptive Gibbs blocking, transition-relevance (SGOOP/committor) — with correctness arguments. The key
> correctness question: does making the block partition a function of the current configuration break
> detailed balance, and what condition restores it? Also scope the companion characterization ("optimal
> within-basin blocking does not transfer across folds") on a CATH NMR-ensemble subset. Correct stationary
> distribution / detailed balance is the top priority.

### Restatement in the codebase's vocabulary

Robosample composes one Markov chain from a fixed scan of **worlds**. A world is a robot factorization
exposing a subset `S` of generalized (BAT) coordinates as free and welding the rest (`spiridon_2020` §2.1.1;
`include/World.hpp` `Selection.bondMobility`). A sampled torsional world runs generalized-coordinate HMC on
`S`, conditioning on the welded complement; it leaves invariant the conditional Boltzmann-with-Fixman
measure of `S` (Claim C1). A **Gibbs block** is that subset `S`. The **block-selection objective** is the
criterion deciding which coordinates share a block. The task is to replace the objective deciding `S` so
the blocks accelerate **barrier-crossing** (inter-basin) sampling, not only **within-basin decorrelation**.

### The defect, verified (not transcribed)

The current pipeline is fully specified in code, independent of any attribution:

1. `autoblock.py::run_pipeline` uses TICA + KMeans + an MSM + PCCA+ only to **detect** transition windows
   and metastable states.
2. It computes the blocking objective in `autoblock.py::transition_correlation_matrix` (lines 323-357): a
   **lag-0** Pearson correlation `np.corrcoef` of the `sin` of each dihedral, over the union of
   transition-window frames.
3. `batstat.py::chose_correlated_bonds` (404-423) runs Louvain community detection on the graph whose edge
   weights are `|corr[i,j]|` (`correlation_to_graph`, 395-402) and returns the communities as blocks.

So the operative objective is an **equal-time (lag-0) correlation matrix**, partitioned by modularity —
even when frames are restricted to transition windows.

**Why that objective is structurally blind to transitions (verified argument).** The equal-time covariance
`C(0) = ⟨δx δxᵀ⟩` is invariant under time reversal and carries no relaxation-timescale information: it ranks
coordinate pairs by instantaneous co-fluctuation *amplitude*, not by how slowly a collective coordinate
*relaxes*. The slowest relaxations — the coordinates whose motion **is** the barrier crossing — are the
leading eigenfunctions of the transfer/Koopman operator, recovered from the generalized eigenproblem
`C(τ) v = λ C(0) v` (tICA). At `τ = 0` this degenerates: every eigenvalue is 1 and eigenvectors are
arbitrary, so there is no slow-mode ranking. Restricting to transition windows re-weights `C(0)` toward
transition geometries but does not change its `τ=0` character. Hence "the matrix is non-stationary across
basins" is the **expected behavior of a lag-0 objective**, not a defect to repair by re-estimating the same
matrix.

NOTE — tier: the *weak* form (the code's objective is lag-0) is grounded in the codebase (authority). The
*strong* form (lag-0 provably cannot rank barrier-crossing coordinates, `C(τ)` provably can) rests on
tICA/VAMP theory **not** in `references/`; stated as proposed derivation D1 with an ingestion request; it
grounds no normative Claim.

### Provenance gap (do not fabricate)

The method has been attributed internally to "Islam & Venugopal." That reference could not be located in
scite or on the web under any query variant (2026-07-14). The closest published method of the class is
*Correlation-Based Feature Selection to Identify Functional Dynamics in Proteins*, JCTC 2022,
`10.1021/acs.jctc.2c00337`. This spec SHALL NOT assert the "Islam & Venugopal" citation. Confirming/
replacing it is OQ-1. NOTE: because the coded objective is fully specified, the *characterization* (Paper 2)
is blocked only on **attribution**, not method definition.

### Binding constraints and assumptions

- **C-STAT (top priority).** The composed chain SHALL leave the correct Boltzmann `π` invariant. A method
  raising crossing rate by perturbing `π` is out-of-constraint (Critical, `CLAUDE.md`). Rate is subordinate
  to correctness.
- **C-FIXED-KERNEL.** The recommended path SHALL NOT alter the HMC guidance dynamics, the Fixman term, the
  acceptance Hamiltonian, or the momentum draw. It changes only which `Selection` a world is built from.
- **C-SCALE.** Targets reach 1M atoms / 100k rigid bodies / 10k robots. The estimator SHALL be computable
  from short trajectories already generated, without a committor solve or biased sampling on the
  recommended path.

## Behavior

### Claims (normative)

- **C1 (per-block `π`-invariance).** Each sampled torsional world — generalized-coordinate HMC with the
  Fixman potential in the acceptance Hamiltonian and MH accept/reject — leaves invariant the conditional
  Boltzmann measure of its active coordinates `S`, `ρ(φ_S|φ_S̄) ∝ √det M_φ(S)·e^{−βU}`, conditioned on the
  welded complement. Grounding: **authority** — `spiridon_2017` eq:2, eq:3, eq:modified-hamiltonian, eq:4;
  `spiridon_2020` eq:6-8; codebase `tests/TestFixmanBoltzmann.cpp`/`World::calcFixman`.
- **C2 (fixed-scan composition preserves `π`).** A Gibbs sweep applying a **fixed** (state-independent) set
  of `π`-invariant block-worlds in fixed or randomized order leaves `π` invariant, regardless of the
  objective used to *construct* the blocks and whether they overlap or tile. Systematic-scan invariance
  needs only each factor `π`-invariant, not reversibility. Grounding: **authority** — `chodera_2011` §II.B,
  §III.A, eq:gibbs-update.
- **C3 (the recommended change is stationary-neutral).** Replacing the objective (lag-0 → time-lagged
  `C(τ)` slow modes, Objective 1) changes only which fixed `Selection` each world is built from. By C1
  every world stays `π`-invariant; by C2 the fixed scan stays `π`-invariant. So Objective 1 leaves the
  sampled distribution **exactly unchanged** with **no** new detailed-balance obligation — the
  dual-Hamiltonian freedom, with even the guidance untouched. Grounding: **authority** — `duane_1987_hmc`
  checks.md (any guidance preserves the target) + C1 + C2.
- **C4 (state-dependent selection breaks `π` unless corrected).** If block choice is a function of the
  current configuration `x` with weights `w_S(x)`, the mixture `Q(x,·) = Σ_S w_S(x) P_S(x,·)` does **not**
  leave `π` invariant in general, even with each `P_S` `π`-invariant, because `w_S(x)` couples to `P_S`
  inside the stationarity integral. Two restorers:
  - **C4a (exact).** `Q` is exactly `π`-invariant if each `w_S(x)` depends on `x` **only through
    coordinates block `S` welds** (`w_S(x)=w_S(φ_{S̄(S)})`) and `Σ_S w_S(x)=1` everywhere; equivalently,
    block selection MAY be a `π`-invariant MH proposal whose acceptance ratio includes the
    reverse-selection ratio `w_S(x')/w_S(x)` and the reverse-block-proposal probability. If `w_S` depends
    on **moved** coordinates and the correction is omitted, `π` is biased (Critical). Grounding:
    **verification pending (LEMMA — no CAS)**; authority anchors for the correction form: `duane_1987_hmc`
    eq:6, `spiridon_2017` eq:4, `chodera_2011` eqs:18-23 (proposal-set symmetry `i∈S_j ⟺ j∈S_i`).
  - **C4b (asymptotic).** If selection is a history-driven **adaptation parameter** (not an instantaneous
    function of moved coordinates), the adaptive chain converges to the common `π` **iff** (i) diminishing
    adaptation and (ii) containment both hold. Diminishing adaptation alone is **insufficient** (containment
    can fail, the chain converges wrong). Grounding: **authority for the warning** — `chodera_2011` §III.A
    ("History-dependent proposals (adaptive weights) break equilibrium unless handled carefully"),
    `wu_2016_sgld_gle` (running-average force breaks detailed balance; SGLD stationary distribution "not
    known explicitly"); **discovery for the conditions** — `10.1214/11-aap806`, arXiv:1801.09299.
- **C5 (transition-relevance objective: stationary-neutral but estimator-expensive).** SGOOP-class blocks
  (maximize MSM spectral gap / committor discrimination) are, like Objective 1, offline construction
  criteria; a fixed set is `π`-invariant by C2. Its estimator needs a committor/reaction coordinate or MSM
  spectral-gap optimization, generally requiring biased sampling or an observed crossing — violating
  C-SCALE as a *construction* driver. Grounding: **authority** (C2) + **discovery**
  (`10.1073/pnas.1600917113`, `10.1063/1.5064856`).
- **C6 (REX is the between-state layer; transition-aware blocks the within-state layer).** REX is a Gibbs
  sampler on `(X,S)` (configurations, state-index permutation): it alternates a per-replica world scan with
  a permutation swap. Replacing the blocks with transition-aware ones leaves the extended-ensemble
  stationary distribution unchanged (C1 + C2 + swap update). Grounding: **authority** — `chodera_2011`
  eqs:8-10, eq:24; `include/Context.hpp` `RunREX`.
- **C7 (transition-relevant blocking already accelerates the slow transition).** A hand-picked reduced
  block of exactly the barrier-crossing torsions (Ramachandran dynamics) accelerates the slowest
  alanine-dipeptide transition ~10× vs fully-flexible MD, at the same `π`. Authority evidence that *which
  coordinates share a block* controls crossing rate. Grounding: **authority** — `spiridon_2020`
  abstract/results, checks.md Table 4 (α_L MFPT).

### Derivation sketches (codebase notation)

- **C1/C2/C3.** `x=(φ_S,φ_S̄)`. `P_S` resamples `φ_S ~ ρ(·|φ_S̄)` via Fixman-corrected HMC (`spiridon_2017`
  eq:momenta-draw / eq:vv-* / eq:accept-pseudocode on `H'=H+U_F`), fixing `φ_S̄`, preserving the conditional
  (C1). For fixed `{S_i}`, `∫ π(dx)(P_{S_1}⋯P_{S_m})(x,A)=π(A)` factor-by-factor (C2). The construction
  objective enters only in choosing `{S_i}`, never in `P_{S_i}` (C3).
- **C4a (LEMMA — CAS pending).** With `π(dx)=π(φ_S̄)π(φ_S|φ_S̄)dφ`, `w_S=w_S(φ_S̄)`, and `P_S` fixing `φ_S̄`
  and preserving the conditional: `∫ π(dx) w_S(φ_S̄) P_S(x,A)=∫_A w_S(φ_S̄)π(dx)`; summing,
  `∫_A[Σ_S w_S]π(dx)=π(A)` iff `Σ_S w_S=1`. Moved-coordinate dependence breaks the equality; restore via
  the MH reverse-selection ratio (`duane_1987_hmc` eq:6). Close with a CAS (`FullSimplify[lhs−rhs]→0`,
  second-CAS numeric spot-check) before Objective 2 ships (V5, OQ-2).
- **C4b.** History-driven `γ_n` yields kernels `P_γ` all sharing `π` (C2); Roberts-Rosenthal ergodicity
  needs diminishing adaptation + containment. Safest realization is **between-run** adaptation: refit the
  objective from pooled history at run end, start a fresh fixed-scan chain — each run is a separate exactly-
  `π`-invariant chain (C2), no within-run condition needed; runs SHALL NOT be concatenated as one Markov
  chain across a selection change. Within-run adaptation is admissible only with a containment argument
  (OQ-2).

### Recommendation (the decision this spec makes)

**Adopt Objective 1 (time-lagged / tICA-VAMP slow-mode objective) as the block-construction method.**
Position Objective 3 (transition-relevance: MSM spectral gap / committor) as the **evaluation metric**, not
the construction driver. Treat Objective 2 (adaptive Gibbs blocking) as an **optional, guarded second
phase** admissible only under C4a or C4b with a verified condition. Rationale, correctness first:

1. **Correctness (decisive).** Objective 1 is stationary-neutral by C3 — it changes only an offline
   criterion, so C-STAT holds **by construction**, with no new detailed-balance obligation and no
   un-verified algebra. Objective 3 is equally stationary-neutral (C5). Objective 2 is the **only** one that
   puts `π` at risk (C4); not the first move.
2. **Repairs the diagnosed defect.** Objective 1 replaces `C(0)` with `C(τ) v = λ C(0) v`; its leading
   generalized eigenvectors approximate the slowest relaxation (barrier-crossing) modes the equal-time
   objective cannot see. Blocks become the coordinate groups loading on the slow tICA components.
3. **Reuses existing machinery.** TICA already runs in `autoblock.py::run_tica`; the change routes its
   slow-mode loadings (not a lag-0 correlation) into `chose_correlated_bonds`. Extends the NMA slow-mode
   plumbing (`DOC-NMA`, `DOC-VelocityDistortion`) and SGLD low-frequency separation (`wu_2011/2012/2016`).
4. **Estimable under C-SCALE.** `C(τ)` is estimable from short trajectories without a committor solve or
   biased sampling — unlike Objective 3 as a construction driver (C5).

Objective 3 earns its place as the **success criterion**: whether a block set is transition-relevant, and
whether blocks transfer across folds, SHALL be measured by the MSM spectral gap / crossing rate the blocks
achieve (Validation; Paper 2). Construct with Objective 1, judge with Objective 3, adapt (if at all) under
Objective 2's proven conditions.

NOTE (tempting but out-of-scope). Feeding tICA slow modes into the **momentum distortion**
(`DOC-VelocityDistortion`, `DistortOption::NMA`) is a *different* change: momentum distortion enters the
acceptance ratio (INV-5), and `wu_2016_sgld_gle` shows biasing dynamics toward slow modes breaks detailed
balance unless explicitly corrected. Not part of this recommendation; if pursued it SHALL carry its own
acceptance-accounting proof (sibling catalog OQ-1).

## Invariants (Correctness conditions)

- **INV-1 (per-world `π`-invariance).** Every world built from any objective SHALL leave its
  Boltzmann-with-Fixman conditional invariant (C1). Reused oracle: sibling O-INV / O-Shirts / O-C4. Critical
  if violated.
- **INV-2 (objective-neutrality of `π`).** Changing the objective SHALL leave sampled `π` unchanged: two
  block sets from two objectives SHALL produce the same basin populations and FES within CI (C3, C5).
  Load-bearing separation of a *rate* change (allowed) from a `π` change (Critical).
- **INV-3 (order-invariance).** For a fixed block set, two scan orders SHALL yield equal populations /
  `ΔF_AB` within CI (C2). Reused from sibling O-INV-2.
- **INV-4 (adaptive-selection gate).** Any state/history-dependent selection (Objective 2) SHALL satisfy
  C4a or C4b with a verified condition; absent that, the partition SHALL stay fixed per run. Critical if an
  unverified within-run adaptive selection ships.
- **INV-5 (no guidance/acceptance leakage).** On the recommended path the objective SHALL influence only
  the `Selection` a world is built from; it SHALL NOT alter the momentum draw, Fixman term, guidance
  integrator, or acceptance Hamiltonian (C-FIXED-KERNEL).

## Interface (Touch list)

Recommended path (Objective 1), the only functional change:

- `autoblock.py::transition_correlation_matrix` (323-357): replace the lag-0 `np.corrcoef` objective with a
  **time-lagged** objective — coordinate-space loadings of the leading tICA/VAMP components of `C(τ)` (the
  `tica` object already fit in `run_tica`, 108-116), reduced to a per-coordinate-pair slow-mode co-loading
  matrix that feeds the existing graph builder. Frame detection (`detect_transitions`) MAY be retained as
  the `C(τ)` sampling region.
- `batstat.py::chose_correlated_bonds` (404-423) and `correlation_to_graph` (395-402): unchanged in
  structure; they consume the new matrix (Louvain now groups by slow-mode co-loading).
- `context.py`: `build_flexibilities` (68-99) → `add_torsional_world`/`add_robotic_world` (783-836) consume
  the per-block bond lists unchanged.
- **Unchanged:** `include/GibbsSweep.hpp`, the world scan in `Context::RunREX`/`runREX`, `World::buildModel`,
  the Fixman/acceptance path, the momentum draw. The world set is still fixed before `initialize`
  (`Context.hpp:28-33`); blocks are one-`Selection`-per-world built at `add*World`.

Objective-2 feasibility (guarded second phase only), from the headers:

- The world set and each world's `Selection` (`sel_`) are **fixed at build time** (`Context.hpp:28-33`;
  `World.hpp:967-973` retains `sys_`/`sel_`/`rootMobilities_`). A public rebuild exists for **root**
  mobility (`World::setRootMobility(ies)`, `World.hpp:309-310`) but **no** public mutator for `sel_`, and
  rebuild "resets per-body sampler state" and must precede `add_sampler`/mass-scaling (`World.hpp:306,324`).
- A runtime *flexibility* change would need a new `setFlexibilities(sel)` rebuild (small, but resets sampler
  state and re-indexes bodies) **or**, preferred, be realized as **runtime selection over a pre-built menu
  of worlds**: instantiate candidate blocks as fixed worlds, let the sweep driver choose which to run each
  sweep. Because Cartesian coordinates are the sole stateless inter-world currency (INV-3;
  `GibbsSweep.hpp`; `World.hpp:391-402`), menu selection needs **no** rebuild — only a driver-level
  selection hook. This reduces "state-dependent partition" to a state-dependent **schedule** over a fixed
  menu, exactly what C4 governs. The menu route SHALL be preferred over a runtime `sel_` mutator.

### Guidance/acceptance split (the canonical bug)

On the recommended path the split is trivially clean: the time-lagged objective enters **neither** side —
consumed offline in block construction; the acceptance ratio stays `H'=H+U_F` (`spiridon_2017`
eq:modified-hamiltonian). For Objective 2, the selection weight `w_S(x)` is a **proposal** element: if it
depends on moved coordinates it belongs on the **acceptance** side as the reverse-selection ratio (C4a;
`duane_1987_hmc` eq:6). Putting a moved-coordinate-dependent selection weight only on the proposal side
without the acceptance correction is the canonical bias (Critical).

### Conventions at risk

- **Fixman determinant-exponent sign (contradiction within `references/`, surfaced not resolved).**
  `spiridon_2017` eq:2 writes marginal `ρ ∝ |M_φ|^{+1/2}e^{−βU}` with `M_φ=JᵀMJ` and eq:3
  `U_F=½β⁻¹ln(|M_{N_f}|/|M_{3N}|)`; `spiridon_2020` `notation.md` defines `M_tot=JᵀMJ` **identically** yet
  `equations.md` eq:4 writes `ρ ∝ |M_tot|^{−1/2}` and eq:7 the inverted ratio; `spiridon_2017`
  eq:momenta-draw gives `Σ=kT·M⁻¹`, inconsistent with `K=½pᵀM⁻¹p` (which implies `Σ=kT·M`). The codebase is
  the operative tie-breaker: `TestFixmanBoltzmann.cpp`/`World::calcFixman` use `p~N(0,M)`, marginal
  `∝ √det M`, `U_F=+½RT ln det M` (the `+1/2` convention). **Orthogonal** to this spec's recommendation;
  does **not** block it; blocks only the momentum-distortion extension. Surfaced as OQ-3.
- **Angular-vs-linear / circular statistics.** The current objective correlates `sin(dihedral)` (single
  circular projection); a time-lagged objective SHALL use a circular-consistent embedding (`sin`/`cos`
  pair, as `extract_dihedrals` already builds) so `C(τ)` respects the torus. A `sin`-only lag-0 matrix
  silently loses half the circular information; fix when the objective changes (feeds INV-2 fairness).

## Validation strategy (Verification plan)

Oracles tagged PRECONDITION / INVARIANT / LEMMA, graded proof (symbolic identity) or falsification
(numeric). Correctness oracles reuse `tests/TestEnsembleValidation.cpp`, `HmcDriver.hpp`, `StatTest.hpp`,
`TestFixmanBoltzmann.cpp` and the sibling catalog's oracles; only objective-neutrality and non-transfer are
new.

- **V1 — PRECONDITION (per-world `π`-invariance).** [falsification, model-independent] C1. Exercises: every
  world built by the new objective passes the Shirts slope test and the C4 uniform-torsion test. Fails on: a
  Selection/rebuild path dropping or mis-signing Fixman. Expected: slope `= −(β₂−β₁)` within one SE (butane
  fixture `0.13363±~4e−3`); uniform torsion histogram (chi-square α=1e−4).
- **V2 — INVARIANT (objective-neutrality of `π`).** [falsification, model-independent] C3/C5, INV-2.
  Exercises: blocks by (a) lag-0 and (b) time-lagged objective on the same system; run each fixed scan to
  convergence; compare distributions. Fails on: any construction leak into guidance/acceptance/momentum
  (INV-5) or a measure-perturbing rebuild. Expected: identical alanine `(φ,ψ)` FES within CI — basin free
  energies vs C7eq agree within 1 std across objectives (fixture: C5 `0.18`, PPII `2.15`, α_L `8.39`
  kJ/mol). NOTE on discriminating power: the **absolute-fixture arm** (basin free energies vs the
  `spiridon_2017` fixture) is the load-bearing check. The (a)-vs-(b) differential arm does NOT by itself
  catch a bias living in the **shared** construction/rebuild path (e.g. a rebuild that mis-signs Fixman,
  hit identically by both objectives) — both objectives then agree with each other while both are wrong.
  That shared-path measure perturbation is caught by **V1 (per-world Shirts test, INV-1) run on every world
  the new objective builds**, and by the absolute-fixture comparison — not by the a-vs-b diff.
- **V3 — INVARIANT (order-invariance).** [falsification, model-independent] C2, INV-3. Exercises: fixed
  transition-aware scan in two orders. Fails on: a per-world defect breaking C1. Expected: equal
  populations / `ΔF_AB` within CI. Reused O-INV-2.
- **V4 — LEMMA (crossing-rate gain).** [falsification] C7. Exercises: time-lagged blocks co-blocking the
  slow coordinates vs lag-0 blocks, at equal compute. Fails on: a time-lagged objective that does not
  co-block slow coordinates. Expected: more independent A↔B crossings per compute / lower MFPT to the slow
  basin (α_L MFPT fixture) **while V2 holds**. Worth is the fails-on only (a rate metric proves nothing
  about `π`).
- **V5 — LEMMA (adaptive detailed-balance algebra, CAS-pending).** [proof] C4a. Exercises: the
  exact-invariance condition and the MH reverse-selection correction. Oracle: CAS reduction
  `FullSimplify[lhs−rhs]→0` (Wolfram) + independent numeric spot-check (SymPy), re-anchored to codebase
  notation, emitted as a re-runnable string. Fails on: selection depending on moved coordinates without the
  acceptance correction (bias detectable by V2 on a two-block toy). **Open until the CAS runs (OQ-2);
  Objective 2 SHALL NOT ship before it closes.** Authority anchors: `duane_1987_hmc` eq:6, `chodera_2011`
  eqs:18-23.

### Paper 2 — companion characterization: block non-transfer across folds (CATH NMR subset)

Depends on the objective being characterized (OQ-1 attribution); because the coded objective is fully
specified, the experiment characterizes the **concrete coded objective** (lag-0 circular-correlation +
Louvain) and, as contrast, the time-lagged objective — not blocked on attribution.

- **Non-transfer metric (partition distance).** For each fold/conformer `f`, compute the optimal partition
  `P_f` of the shared coordinate set. Non-transfer = mean over fold pairs of normalized **Variation of
  Information** `VI(P_f,P_g)=H(P_f)+H(P_g)−2 I(P_f,P_g)` (a true metric, bounded by `log N`), normalized
  `VI/log N`. Complement with a **performance-transfer** matrix `T(f→g)=eff(blocks_f on g)/eff(blocks_g on
  g)`, `eff` = achieved MSM spectral gap or crossings-per-compute (Objective 3 as evaluation); `T ≪ 1`
  means functional non-transfer.
- **Control (mandatory — noise floor).** Within-fold reproducibility `VI(P_f^{(1)},P_f^{(2)})` between two
  independent trajectories of the **same** fold. The control SHALL be estimated under the **identical data
  modality, sample size (frame count), tICA lag, and Louvain resolution** as each cross-fold pair it gates;
  a control estimated under different settings does not bound the cross-fold statistic's noise. Because the
  design permits building the objective matrix from NMR models and/or short trajectories, a cross-fold `VI`
  computed across mismatched modalities (fold `f` from NMR models, `g` from a trajectory) can exceed a
  same-modality within-fold floor purely from data-type / sample-size artifacts. The test is therefore
  **per-pair**: claim non-transfer for `(f,g)` only when `VI(f,g) > max(floor_f, floor_g)` with CI, both
  floors estimated under that pair's settings — NOT a mean-cross-fold vs pooled-baseline comparison, which
  hides non-exceeding pairs. This is the difference between a null with error bars and a shrug.
- **Design.** CATH subset spanning distinct folds (reps across Class/Architecture/Topology), each with a
  multi-MODEL NMR ensemble. Per protein: build the objective matrix (NMR models and/or short trajectories),
  compute `P_f`, then `VI` across folds and `T(f→g)`. Report both objectives. Held-out discipline: the fold
  set used to tune hyperparameters (tICA lag, `n_components`, Louvain resolution) SHALL be disjoint from the
  fold set used to report `VI`/`T`.
- **Tags.** `VI`/`T` are [LEMMA | falsification] characterization oracles (measure *transfer*, not `π`).
  Correctness of every block set used is still gated by V1/V2/V3.

## Consequences and trade-offs

- Choosing Objective 1 (fixed time-lagged blocks) forecloses adaptive selection (Objective 2) as the *base*
  method; adaptive selection is admitted only behind INV-4. Trades a modest efficiency ceiling for
  correctness safety (consistent with the sibling catalog).
- Using Objective 3 as the **evaluation** metric rather than **construction** driver trades a
  sharper-but-expensive objective for one estimable under C-SCALE; the committor/SGOOP path stays available
  as a later upgrade if crossing rate plateaus.
- The time-lagged objective adds two hyperparameters (lag `τ`, retained components) whose mis-setting
  degrades *rate* but not *correctness* (INV-2 holds regardless) — the safe kind of knob.
- Preferring the pre-built-world-menu route for future adaptivity trades memory for a localized,
  sampler-state-preserving change and a clean mapping onto the C4 analysis.

## Open questions

- **OQ-1 (block-objective attribution).** "Islam & Venugopal" unlocatable (scite + web, 2026-07-14);
  closest is JCTC 2022 `10.1021/acs.jctc.2c00337`. Blocks Paper 2's *framing/attribution* only. Answer
  needed: the real citation or acceptance of the JCTC-2022 surrogate. Return to the human.
- **OQ-2 (adaptive detailed-balance algebra — CAS-pending, blocks Objective 2).** C4a/C4b un-verified this
  pass (no CAS). Blocks INV-4/V5 and any within-run adaptive selection. Answer needed: (i) a CAS-closed
  proof of C4a cross-checked in a second CAS, re-anchored to codebase notation; (ii) ingestion of
  `10.1214/11-aap806` and arXiv:1801.09299 so a Robosample-specific containment argument becomes
  authority-grade.
- **OQ-3 (Fixman determinant-exponent convention).** Contradiction within `references/` (`spiridon_2017`
  `|M|^{+1/2}`/`Σ=kT M⁻¹` vs `spiridon_2020` `|M|^{−1/2}` with `M_tot=JᵀMJ`), codebase operatively `+1/2`.
  Does **not** block this spec's recommendation; blocks any momentum-distortion extension (shared with
  sibling OQ-1). Answer needed: a reconciled statement of the marginal exponent, momentum covariance, and
  Fixman sign matching the code.

### Ingestion requests (external math needed before it can originate a Claim)

Discovery-tier here; human-verify before grounding a normative Claim. Peer-reviewed venues, scite-confirmed
(2026-07-14).

- **D1** (grounds C3 diagnosis strong-form + Objective-1 construction): VAMPnets `10.1038/s41467-017-02388-1`
  (Nat. Commun.); Variational Approach `10.1007/s00332-019-09567-y` (J. Nonlinear Sci.). Verify: leading
  generalized eigenvectors of `C(τ)v=λC(0)v` approximate transfer-operator slow eigenfunctions; `τ→0`
  degeneracy.
- **D2/D3** (grounds C4b; blocks Objective 2): Latuszynski-Roberts-Rosenthal `10.1214/11-aap806` (Ann.
  Appl. Probab.; diminishing adaptation + containment, AdapRSG counterexample); Chimisov-Latuszynski-Roberts
  arXiv:1801.09299 (adapts random-scan selection probabilities; containment nontrivial). Verify: the two
  conditions and the counterexample.
- **D4** (grounds Objective-3 evaluation metric): SGOOP `10.1073/pnas.1600917113` (PNAS); multi-CV SGOOP
  `10.1063/1.5064856` (J. Chem. Phys.). Verify: spectral-gap objective and estimator requirements.
- **D5** (surrogate for coded objective's provenance): `10.1021/acs.jctc.2c00337` (JCTC 2022). Verify: match
  to the coded lag-0 + community-detection pipeline.

## References

Authority (`references/`): `spiridon_2017_cdhmc_gibbs`, `spiridon_2020_robosample`,
`chodera_2011_gibbs_replica_exchange`, `duane_1987_hmc`, `shirts_2012_ensemble_validation`, `wu_2011_sgld`,
`wu_2012_rxsgld`, `wu_2016_sgld_gle`.

Codebase: `python/robosample/autoblock.py`, `batstat.py`, `context.py`; `include/Context.hpp`,
`include/World.hpp`, `include/GibbsSweep.hpp`; `tests/TestEnsembleValidation.cpp`,
`tests/TestFixmanBoltzmann.cpp`, `tests/StatTest.hpp`; `docs/specs/refactor/DOC-GibbsSweep.md`, `DOC-NMA.md`,
`DOC-VelocityDistortion.md`; `docs/specs/torsional-dynamics-efficiency-and-crossing.md`;
`docs/RESEARCH-ROADMAP.md`.

Discovery (to ingest, not yet authority): `10.1038/s41467-017-02388-1`, `10.1007/s00332-019-09567-y`,
`10.1214/11-aap806`, arXiv:1801.09299, `10.1073/pnas.1600917113`, `10.1063/1.5064856`,
`10.1021/acs.jctc.2c00337`.
