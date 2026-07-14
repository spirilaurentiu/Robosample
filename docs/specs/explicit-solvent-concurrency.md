# Spec: Native articulated-body HMC concurrent with unmodified explicit solvent

Status: **PARTIAL** — researcher artifact, ready for review. The discriminating diagnostic and its
authority-grounded acceptance backbone (C1–C3) are complete; the winning solution branch's full
derivation depends on discovery-grade external papers not yet in `references/` (D1–D3) and on symbolic
checks that require a CAS not available in the researcher environment. Detailed balance is the top
correctness priority: this spec commits to *diagnosis before remedy* and does not adopt any fix that
could bias the stationary distribution before the mechanism is known.

Revision (2026-07-14, post spec-review): B1 applied — D1 configurational-freezing no longer claims to
preserve C1; gated on OQ-2. S1 — decision table gates every M1 verdict on L-N extensivity. S2 — added
P-REVERSIBLE-B oracle so an M0 defect confined to the OpenMM path cannot be misread as M1. S3 — C3's
closed form relabeled as the Fermi/Barker acceptance (verification-tier heuristic); only the monotone
work-variance content is authority-grounded. S4 — INV-KERNEL regraded proof → falsification. The
three-mechanism reframe was independently confirmed against the analytic fixture.

Provenance: derived from a full read of the NCMC machinery (`include/World.hpp`
`ncmcMove`/`ncmcInnerGhmcStep`; `src/world/sampler/NcmcMove.cpp`; `include/NCMCProtocol.hpp`;
`python/robosample/context.py::add_ncmc_world`, `configure_ncmc_region`, `set_cartesian_solvent`,
`use_metropolized_inner`), the existing OpenMM-free oracle `tests/TestNcmcExplicitSolvent.cpp`, and
authority-tier equations in `references/`: `nilmeier_2011_ncmc` (eqs 11,12,14,15,16,18,19,20),
`crooks_1998` (eqs 5,6,9,10), `spiridon_2017_cdhmc_gibbs` (eqs 2,3, modified-hamiltonian, momenta-draw),
`wyczalkowski_pappu_2008` (eq 18), `ballard_2009_rens` (eq 2). Requirement keywords follow RFC 2119 per
`styles/spec.md`.

---

## Motivation

### Problem restatement (user's phrasing)

> Make native articulated-body (BAT / internal-coordinate) HMC sampling work concurrently with UNMODIFIED
> explicit solvent, at acceptable acceptance — attacked directly, not via an implicit-search /
> explicit-rescore posture. Prior local diagnosis may have conflated two separable acceptance-killing
> mechanisms: (1) integrator shadow work ~ dt²·nDOF, and (2) steric/excluded-volume clash of a concerted
> torsional move against frozen solvent. Design a decisive, cheap diagnostic that separates them first.

### The problem in the codebase's vocabulary, quantified

Robosample proposes concerted internal-coordinate (BAT) moves for a solute embedded in explicit
TIP3P/PME. The per-molecule NCMC switch (`World::ncmcMove`, a λ:1→0→1 palindrome, `NCMCProtocol.hpp`) is
the escape route: decouple the solute's nonbonded interaction with the bath, make a large uncaged
torsional move near λ=0, recouple, and accept exactly. In practice the move accepts at **≈0** in explicit
solvent (campaign `ncmc-explicit-solvent-acceptance`). Two documented remedies — soft-core decoupling and
`relax_solvent` — do not recover acceptance; `relax_solvent` under the endpoint-ΔH construction
(Construction I) is self-defeating because it *adds* integrated DOF.

The prior campaign labelled the killer "integrator shadow work" and built Construction II
(Metropolized-dynamics NCMC, `use_metropolized_inner`). Built and run at exhaustive tier on 2026-07-13, it
**still fails**: the primary OpenMM-free oracle `INV0_ConstructionIAndIIAgree` reports chi² ≈ 2.0×10⁵
against a critical value of 57; the φ₂ marginal under-mixes (moves accept but φ₂ barely moves).

Quantified observation that reframes the problem: `INV0` fails in `twoTorsionChainWithAtoms` — a
**two-torsion, five-atom analytic fixture with no explicit bath and no excluded-volume interaction**
(`tests/TestNcmcExplicitSolvent.cpp:84-108`, `HarmonicBridge`). In that fixture, bath-extensive shadow
work (mechanism 1) and steric clash (mechanism 2) are both **excluded by construction**: N_int ≈ 2 and
there is no LJ core. A failure there cannot be either mechanism the user named. The prior single "shadow
work" label therefore conflates **three** distinct mechanisms, not two:

- **M0 — kernel/algorithmic defect.** The composed move fails to sample correctly even with no bath and no
  cavity. Evidence: `INV0` above; and the campaign's own history of two successive frozen-coordinate bugs
  (a `restorePreStep` return-true bug, then a stale-acceleration-seed bug in `ncmcInnerGhmcStep`).
- **M1 — bath shadow work.** Integrator discretization heat, extensive in integrated solvent DOF,
  vanishing as dt→0. This is the user's mechanism (1).
- **M2 — steric clash / cavity mismatch.** A configurational-energy penalty from displacing the solute
  into frozen (or slowly-relaxing) solvent; dt-independent; vanishes without a dense excluded-volume
  neighbourhood. This is the user's mechanism (2).

The pivotal deliverable is a diagnostic that assigns the acceptance deficit to M0, M1, or M2
unambiguously and cheaply, because each implies a different remedy with a different correctness surface.
Target scale for the eventual remedy: FFAR1 (flexible GPCR) in explicit membrane/solvent, up to the
1M-atom regime, where a bath-extensive term is fatal and a solute-local term is not.

Binding constraints and assumptions:

- The explicit-solvent force field, PME, and periodic box SHALL remain unmodified (the concurrency
  requirement). Any remedy operates inside the NCMC/inner-kernel framework or by freezing *dynamics*,
  never the force field.
- Correct stationary distribution (detailed balance) is the top priority; a scheme that biases the
  distribution is worse than one that accepts rarely.
- The Robosample target is the internal-coordinate marginal ρ(φ) ∝ |M_φ|^{1/2} e^{−βU(φ)}
  (`spiridon_2017` eq:2), realized via the acceptance Hamiltonian H' = H + U_F
  (`spiridon_2017` eq:modified-hamiltonian, U_F from eq:3).

---

## Behavior

### Claims

**C1 (grounding: authority).** For a cyclic palindromic λ-protocol with a reversible (π_λ-invariant) inner
propagation kernel, the outer Metropolis test on protocol work alone, min(1, e^{−w}), leaves exp(−βH′)
invariant, where w is the perturbation-step work and H′ = H + U_F. Source: `nilmeier_2011_ncmc` eq:15
(per-step detailed balance) + eq:16 (work) + eq:19 (work replaces energy difference under reversible
propagation); H′ assembly from `spiridon_2017` eq:modified-hamiltonian. convention_scope: universal NCMC
result instantiated with the project-specific Fixman H′.

**C2 (grounding: authority).** The composed move "torsional proposal ∘ solvent-only relaxation ∘
acceptance" is the special case of C1 in which the inner kernel K_t acts only on solvent DOF at fixed λ.
It preserves the joint stationary distribution iff (i) the solvent-relaxation kernel is π_λ-invariant
(`nilmeier` eq:15) and (ii) the active/frozen region partition and protocol-selection probability are
symmetric under path reversal, i.e. the factor P(Λ̃|x̃_T,λ_T)/P(Λ|x_0,λ_0) in `nilmeier` eq:12 equals 1
(or is carried explicitly). The relaxation is *propagation* (K_t), not *perturbation* (α_t); its heat q
SHALL NOT enter the work accumulation. Source: `nilmeier_2011_ncmc` eqs 11,12,16,19.

**C3 (grounding: authority for the monotone content; heuristic for the closed form).** In the
near-reversible regime, any mechanism that inflates Var(W) — M1 bath shadow work or M2 clash — depresses
acceptance; the diagnostic attributes the dominant contributor. This monotone statement (Var(W)↑ ⇒
acceptance↓) is the grounded backbone the diagnostic needs. NOTE: the specific closed form
⟨A⟩ ≈ 1/2 − (β²/4)·Var(W) cited from `wyczalkowski_pappu_2008` eq:18 is the **Fermi/Barker (logistic)**
swap-acceptance `p_swap = f(βΔU_swap)`, whose intercept is 1/2 at zero work variance. Both constructions in
this spec use the **Metropolis** rule `min(1, e^{−w})`, whose mean acceptance near reversibility → 1, not
1/2. Only the monotone content transfers; the `1/2 − (β²/4)Var(W)` form is a verification-tier
near-equilibrium heuristic for a different acceptance function, not the Metropolis-rule acceptance the code
uses. Source (Fermi form): `wyczalkowski_pappu_2008` eq:18 (ΔU_swap ≡ W; C_λ = Var(∂U/∂λ)); the monotone
work-variance dependence for the Metropolis rule follows from `crooks_1998` (work fluctuation).

**C4 (grounding: verification — LEMMA; CAS unavailable in the researcher environment).** The two solvent
mechanisms carry orthogonal, measurable signatures:

- (i) Integrator shadow heat obeys ⟨Q_shadow⟩ ≈ A·N_int·ε^p with p ≥ 2, so ⟨Q_shadow⟩ → 0 as the
  timestep ε → 0 and grows ~linearly in the number of integrated DOF N_int. For equilibrium HMC the
  exponent is p = 4 with variance ∝ N·ε⁴ (discovery support: Beskos et al. 2013, 10.3150/12-bej414;
  Mangoubi & Smith 2018, arXiv:1802.08898); for the driven, un-Metropolized NCMC protocol the local
  O(ε²) drift does not average to zero and p ∈ [2,4]. The falsifiable content is the sign of the exponent
  (p > 0) and bath-extensivity, not its exact value — the diagnostic **measures** p (Axis-DT).
- (ii) Clash energy ΔU_clash is an instantaneous potential difference: dt-independent, growing
  super-linearly with move size Δθ once atoms enter the repulsive LJ core (∝ r^{−12}), and vanishing when
  the solute has no dense excluded-volume neighbourhood (vacuum, or implicit GB with no explicit waters).

Marked LEMMA because `mcp__wolfram__*` and the `verifier` sub-agent are absent in the researcher
environment; the driven-NCMC exponent and the clash tail SHALL be CAS-verified before the winning branch
is implemented (OQ-1).

### Derivation sketch (codebase notation)

The NCMC acceptance ratio (`nilmeier` eq:12) for a single expanded-ensemble state (ω_T = ω_0), symmetric
protocol selection, and torsion/λ perturbations with unit Jacobian ratio (`nilmeier` eq:14, α-ratio = 1)
reduces to A = min{1, e^{−ΔS − Δu}}. Two propagation regimes bracket the codebase's two constructions:

- **Construction I (symplectic Verlet):** ΔS = 0 (`nilmeier` eq:20) ⇒ A = min{1, e^{−Δu}}. With a cyclic
  palindrome (λ_T = λ_0 = 1), Δu = βΔH′ — the endpoint Hamiltonian change **including** all discretization
  error. Decompose βΔH′ = β(ΔU_config + ΔK + ΔU_F) + Q_shadow. Clash is the excluded-volume component of
  ΔU_config; shadow work is Q_shadow. Both are inside the accepted exponent.
- **Construction II (reversible/Metropolized inner kernel):** ΔS = −q (`nilmeier` eq:18) ⇒
  A = min{1, e^{−w}} (`nilmeier` eq:19). The propagation heat q — including Q_shadow — drops out; only
  perturbation work w remains. This is C1, and it is why Construction II is the exact remedy **for M1
  specifically**.

The three-mechanism separation follows from where each term lives and how it scales: Q_shadow scales as
N_int·ε^p (C4-i, in q, removed by Construction II); ΔU_clash scales with Δθ and vanishes without a cavity
(C4-ii, in w — NOT removed by Construction II, because relaxation must physically make room); and M0 is
any failure that survives when both N_int and the cavity are removed.

### Conditional solution architecture (Behavior of the remedy, branched on the diagnostic readout)

The spec commits to the diagnostic and to the branch-selection; it specifies each branch's correctness
conditions and leaves the numeric implementation of the winning branch to a follow-up once the diagnostic
reads out.

- **If M0 dominates:** repair the inner GHMC kernel to be per-step π_λ-invariant — correct reversible
  leapfrog, energy-assembly consistency between the Hbefore/Hafter evaluations (`ncmcInnerGhmcStep`), and
  reject-not-drift on corrector non-convergence. No new sampling architecture; the target is already C1.
  This SHALL precede any M1/M2 work.
- **If M1 (bath shadow work) dominates:** (a) per-step Metropolization (Construction II) removes Q_shadow
  exactly (C1); (b) shrink N_int by decoupling/relaxing only the moving substructure's neighbourhood —
  configurational freezing (D1); (c) protocol-length or mass-scaling to lower per-step ε-error. Branches
  (a) and (c) preserve C1. Branch (b) is **conditional on OQ-2**: a configuration-dependent freezing region
  generically breaks the `nilmeier` eq:12 reverse-selection symmetry `P(Λ̃|x̃_T,λ_T)/P(Λ|x_0,λ_0)`, so its
  detailed balance is UNVERIFIED until the region-selection probability is shown symmetric or the
  correction is carried explicitly. D1 SHALL NOT be built as if it preserves C1 before OQ-2 closes.
- **If M2 (steric clash) dominates:** insert a **new solvent-only relaxation Gibbs block** between
  torsional proposal and acceptance. Correctness is C2. Reference construction: Ribeiro et al. mixed MC/MD
  (D2), discovery-grade.
- **Cross-cutting alternative (M1 or M2):** replace the concerted-rotation proposal with a
  configurational-bias / concerted-rotation move (Dodd–Boone–Theodorou lineage, D3) that builds moves
  surviving dense excluded volume; requires the CBMC/Rosenbluth Jacobian in acceptance. Discovery-grade.

---

## Invariants

- **INV-1 (detailed balance).** Every accepted branch SHALL leave exp(−βH′) invariant (C1/C2). Highest
  priority.
- **INV-2 (unmodified solvent).** The explicit-solvent force field, PME, and box SHALL be unchanged by the
  diagnostic and by any remedy that stays inside the NCMC/inner-kernel framework. Region-freezing (D1)
  freezes dynamics, not the potential.
- **INV-3 (construction equivalence).** Constructions I and II SHALL sample the same marginal (existing
  `INV0` agreement arm).
- **INV-4 (identity-move nullity).** A closed palindrome with Δθ = 0 has w = 0 exactly (existing `L1`).
- **INV-5 (Fixman completeness).** The inner accept SHALL use the same H′ assembly (`currentTotalEnergy`)
  as the outer accept, structurally foreclosing the F1 bug (inner accept silently dropping U_F/pitch).

---

## Interface

Touch list (the diagnostic is primarily new tests + instrumentation; the only engine behavioral change is
in the M0 branch):

- **New test TU** `tests/TestNcmcMechanismDiagnostic.cpp`: the three system levels (S-A/S-B/S-C) × three
  axes (DT/MOVE/N) with the seven oracles below. Ungate the three `GTEST_SKIP`'d oracles in
  `tests/TestNcmcExplicitSolvent.cpp` once INV-KERNEL passes.
- **Instrumentation:** extend the per-move NCMC log (already emits `work`, `dH`, `gap`, `dPE/dKE/dKEsolv`;
  `ROBO_NCMC_DEBUG` per-substep dH) with ΔU_config (integrator-off potential delta), N_int, ε, and Δθ so
  the sweeps are extractable from `run.py` CSV.
- **Engine entry points** (read-mostly; behavioral change only in the M0 branch): `World::ncmcMove` /
  `ncmcInnerGhmcStep` (`include/World.hpp`; `src/world/sampler/NcmcMove.cpp`); the Construction-II toggle
  (`set_ncmc_construction_ii` / `SamplerConfig`); `set_cartesian_solvent` (N_int / relax path);
  `configure_ncmc_region` (Region-A / D1 freezing).
- **Python:** add a deterministic fixed-Δθ / instantaneous-proposal (T=1) / ε-sweep diagnostic mode to
  `add_ncmc_world`'s driver; the default behavior of `add_ncmc_world` is unchanged.
- **Guidance/acceptance split (the canonical bug this codebase is prone to):** shadow work and clash both
  live in **acceptance** (Δu / w / Var(W)); the Fixman U_F lives in acceptance only (H′, not guidance). No
  new term belongs in guidance. The M2 solvent-relaxation block adds *propagation* (guidance) steps whose
  heat q is removed from the acceptance exponent by C2 — the split SHALL keep the relaxation OUT of the
  work accumulation (it is K_t, not α_t).

---

## Validation strategy

The centerpiece is the discriminating diagnostic. Three nested system levels remove one candidate
mechanism at a time; three sweep axes expose the orthogonal signatures of C4; seven oracles read out
M0/M1/M2.

System levels (nested controls):

- **S-A** — analytic OpenMM-free fixture (`twoTorsionChainWithAtoms` + `HarmonicBridge`). No bath, no
  cavity. Isolates M0.
- **S-B** — solute (alanine dipeptide) in vacuum / implicit GB via the real OpenMM `ForceBridge`. Real
  integrator and intramolecular strain; no explicit cavity. Isolates M1-of-solute-DOF and is the clash
  control.
- **S-C** — the same solute in explicit TIP3P/PME. Full bath + cavity. Production target.

Sweep axes (fixed protocol shape unless noted):

- **Axis-DT** — vary ε (×{1, 1/2, 1/4, 1/8}) at fixed move geometry. Fit slope p of log⟨W⟩ (Construction
  II) or log⟨Q_shadow⟩ (Construction I; Q_shadow = ΔH′ − ΔU_config) vs log ε.
- **Axis-MOVE** — vary the deterministic active-torsion displacement Δθ at fixed small ε (or T=1
  instantaneous proposal, integrator off). Measure ΔU_config(Δθ) and A(Δθ).
- **Axis-N** — vary N_int (relax radius / number of relaxed waters via `set_cartesian_solvent`, or box
  padding) at fixed ε, Δθ. Fit ⟨W⟩ vs N_int.

Oracles (tag; what it exercises; fails-on; grade):

- **P-MOVE (PRECONDITION).** The proposal SHALL displace the active coordinate, |Δξ| ≥ Δξ_min. Fails-on:
  acc ≈ 100% with configuration bit-identical (the exact frozen-no-op illusion the campaign hit). Without
  it every downstream statistic measures nothing. Grade: falsification; model-independent.
- **P-ENERGY (PRECONDITION).** Identity move (Δθ = 0, closed palindrome) at the smallest ε SHALL give
  |w| → 0 and |ΔH′| → 0 within tol. Fails-on: a **dt-invariant** nonzero w/ΔH′ = energy-assembly /
  stale-accel-seed defect (M0). Discriminator: rerun across Axis-DT — dt-invariant nonzero ⇒ defect;
  ∝ ε^p ⇒ genuine shadow. Grade: falsification.
- **P-CONTROL (PRECONDITION).** On S-B, acceptance at the production ε SHALL exceed A_ctrl (high, e.g.
  ≥ 0.5). Fails-on: if S-B already collapses, the deficit is not solvent-specific and the explicit-solvent
  framing is premature (blocker is M0 or the solute integrator). Grade: falsification.
- **P-REVERSIBLE-B (PRECONDITION; M0-in-OpenMM-path gate).** On S-B (real OpenMM `ForceBridge` + PME code
  path, no cavity), a per-step palindrome round-trip SHALL close: propagate forward then reverse and recover
  the initial `(q,u)` and H′ within tol, and a closed identity-move palindrome SHALL give |w| → 0 across
  Axis-DT. Fails-on: a reversibility / energy-assembly M0 defect that lives ONLY in the OpenMM integration
  path — it passes INV-KERNEL (which runs on the OpenMM-free S-A fixture) and passes P-ENERGY's dt-invariance
  test if it scales as ε^p, then inflates the L-DT slope on S-C and is misread as M1. L-DT and L-N on S-C are
  INADMISSIBLE until P-REVERSIBLE-B passes. Grade: falsification. NOTE: INV-KERNEL (S-A) and P-REVERSIBLE-B
  (S-B) together bracket M0 across both the analytic and the OpenMM force paths; neither alone suffices.
- **INV-KERNEL (INVARIANT; M0 gate; repurposes existing INV0).** On S-A the composed move SHALL recover
  the Fixman-Boltzmann marginal (flat φ₂) within chi² at α = 1e-4, AND the Fixman-omitting arm SHALL
  diverge to the sqrt(det M) shape. This is π-invariance with no bath and no clash. If INV-KERNEL fails,
  the pipeline STOPS at M0; M1/M2 attribution on S-C is inadmissible until it passes. Grade: falsification
  (closed-form reference weights, finite-sample chi² at α=1e-4 — the reference weights are deductive but the
  check itself is a statistical falsification test with finite false-positive rate and power, not a proof).
- **L-DT (LEMMA; M1 signature).** On S-C, log⟨W⟩ (Construction II) or log⟨Q_shadow⟩ (Construction I) vs
  log ε has slope p ≥ 2 and →0 as ε→0. Fails-on / falsification: flat in ε ⇒ shadow work is NOT the
  dominant killer. Grade: falsification; model-independent for the flat-vs-sloped verdict.
- **L-N (LEMMA; M1 bath-extensivity).** On S-C, ⟨W⟩ grows ~linearly with N_int. Falsification: flat in
  N_int ⇒ not bath shadow work. Grade: falsification.
- **L-CLASH (LEMMA; M2 cavity mismatch).** The gap ΔU_config^{S-C}(Δθ) − ΔU_config^{S-B}(Δθ) is positive,
  grows with Δθ, and is dt-independent. Falsification: gap ≈ 0 ⇒ no clash. Grade: falsification.

Decision table (applied only after INV-KERNEL passes):

| | L-CLASH gap grows with Δθ | L-CLASH gap ≈ 0 |
|---|---|---|
| **L-DT slope p > 0 AND L-N extensive** | M1 + M2 both present | M1 shadow-dominant |
| **L-DT slope p > 0 AND L-N flat** | M2 + residual/solute shadow — re-examine | residual M0 (OpenMM path) or solute shadow — re-examine, NOT bath M1 |
| **L-DT slope p ≈ 0** | M2 clash-dominant | neither → re-examine residual M0 / proposal too small |

Every "M1" verdict is gated on **L-N extensivity** (⟨W⟩ ~ linear in N_int). A `p > 0` slope alone does not
imply bath shadow work: a residual M0 defect living in the OpenMM force/integration path but scaling as ε^p,
or solute-only shadow work, produces `p > 0, clash ≈ 0` without being the bath-extensive mechanism the M1
remedy targets. Only `p > 0 AND L-N extensive` selects the M1 (Construction II / D1) branch. Each remaining
cell selects the corresponding solution branch under Behavior. The diagnostic is cheap: S-A is
OpenMM-free and already exists; S-B is implicit GB; only S-C requires the explicit-solvent harness, and it
reuses the existing per-move logging.

---

## Consequences and trade-offs

- Diagnosis-first is deliberate: each remedy has a different correctness surface, and committing before the
  mechanism is known risks biasing the distribution (violating INV-1).
- Construction II (M1 remedy) costs per-step energy evaluations; it trades acceptance for force-eval cost
  (efficiency metric `nilmeier` eq:23).
- The M2 solvent-relaxation block adds a detailed-balance obligation (C2) and couples move quality to
  bath-dynamics quality.
- Region-freezing (D1) shrinks N_int but introduces a region-selection-probability correction that MUST be
  symmetric or explicitly carried (`nilmeier` eq:12).

---

## Open questions

- **OQ-1 (blocks M1 attribution).** Measured exponent p in ⟨Q_shadow⟩ ∝ N_int·ε^p on S-C. No CAS/verifier
  in the researcher environment; the driven-NCMC shadow-work exponent is a LEMMA measured by Axis-DT and
  SHALL be CAS-cross-checked before the M1 remedy is implemented.
- **OQ-2 (blocks D1).** Does Sindhikara configurational-freezing satisfy `nilmeier` eq:12 symmetry for a
  configuration-dependent freezing region? **Ingestion request:** 10.1021/ct500340b (JCTC, peer-reviewed)
  — verify the region-selection probability/Jacobian preserves detailed balance. Place the freezing math
  in *Proposed derivation (UNVERIFIED)*.
- **OQ-3 (blocks D2).** Ribeiro composed-move (MC/MD relaxation) detailed-balance argument. **Ingestion
  request:** 10.1002/jcc.22925 (J Comput Chem, peer-reviewed) — verify the composed-move balance for the
  new solvent-only relaxation Gibbs block.
- **OQ-4 (blocks D3).** Dodd–Boone–Theodorou concerted-rotation Jacobian for the CBMC alternative.
  **Ingestion request:** DBT concerted rotation (Mol. Phys. 1993) + parallel-rotation 10.1063/1.1371496.
- **OQ-5.** Is M0 fully resolved? The existing `INV0` still fails post stale-accel-seed fix; residual
  kernel bug vs mistuned knobs (`kNcmcSteps=8`, `kDt=4e-4`, extrapolated-not-run)? Answer: run INV-KERNEL
  with a knob sweep on S-A.
- **OQ-6.** `scite` returned no Smart-Citation tallies for D1/D2 (paywalled, not OA-indexed), so
  supported/contested status is UNVERIFIED; both are peer-reviewed. Obtain citation status on ingestion.
- **Provenance gap.** `docs/specs/ncmc_solvent_relax.md` is referenced by `context.py` and the `.pyi`
  stubs but is absent from the tree; the prior 5-file `docs/specs/ncmc-explicit-solvent/` spec cited in
  memory is also absent. Neither should be cited as authority until restored. `nilmeier_2011_ncmc` is now
  present in `references/` (contrary to the older memory note that it was missing).
- **Tooling gap.** No `mcp__wolfram__*` and no `verifier` sub-agent in the researcher environment; all
  symbolic scaling is hand-derived LEMMA. The p-exponent (OQ-1) and the C2 composed-move Jacobian (OQ-3)
  SHALL be CAS-verified before the winning branch is implemented.

---

## References

- Nilmeier, J. P., Crooks, G. E., Minh, D. D. L., & Chodera, J. D. (2011). Nonequilibrium candidate Monte
  Carlo is an efficient tool for equilibrium simulation. *PNAS*. https://doi.org/10.1073/pnas.1106094108
  (authority: `references/papers/nilmeier_2011_ncmc`)
- Crooks, G. E. (1998). Nonequilibrium measurements of free energy differences for microscopically
  reversible Markovian systems. *J. Stat. Phys.* https://doi.org/10.1023/A:1023208217925 (authority:
  `references/papers/crooks_1998`)
- Spiridon, L., & Minh, D. D. L. (2017). Constrained-dynamics HMC as Gibbs sampling. (authority:
  `references/papers/spiridon_2017_cdhmc_gibbs`)
- Wyczalkowski, M. A., & Pappu, R. V. (2008). Satisfying the fluctuation theorem in free-energy
  calculations with Hamiltonian replica exchange. *Phys. Rev. E*.
  https://doi.org/10.1103/PhysRevE.77.026104 (authority: `references/papers/wyczalkowski_pappu_2008`)
- Ballard, A. J., & Jarzynski, C. (2009). Replica exchange with nonequilibrium switches. *PNAS*.
  https://doi.org/10.1073/pnas.0900406106 (authority: `references/papers/ballard_2009_rens`)
- Sindhikara, D. J., et al. (2014). Nonequilibrium candidate Monte Carlo simulations with configurational
  freezing schemes. *JCTC*. https://doi.org/10.1021/ct500340b (discovery; citation status unavailable via
  scite; ingestion requested)
- Ribeiro, A. A. S. T., et al. (2011). Mixed Monte Carlo/molecular dynamics simulations in explicit
  solvent. *J. Comput. Chem.* https://doi.org/10.1002/jcc.22925 (discovery; citation status unavailable
  via scite; ingestion requested)
- Beskos, A., et al. (2013). Optimal tuning of the hybrid Monte Carlo algorithm. *Bernoulli*.
  https://doi.org/10.3150/12-BEJ414 (discovery, dt/dimension scaling)
- Mangoubi, O., & Smith, A. (2018). Dimensionally tight bounds for second-order Hamiltonian Monte Carlo.
  https://doi.org/10.48550/arXiv.1802.08898 (discovery, dt/dimension scaling)
