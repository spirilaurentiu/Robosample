# Spec: Fixman-potential validation on idealized serial chains + proline realism

Status: revised (post-empirical-validation round) — ready for review
Owner: coder (after review)
Scope: read-only research artifact. Reproduces the Fixman-potential validation of
Jain et al. 2013, Patriciu et al. 2004, and Spiridon & Minh 2017 on idealized
**serial (acyclic) bead chains** (C4, C5, C11, C15), then adds one **realistic
cyclic** check on the proline dipeptide. All systems run in **vacuum**. This spec
defines *what is correct*; it does not prescribe coding order.

Gate: every Tier-0 deterministic identity passes always-on; the Tier-1/Tier-2
statistical checks pass under `ROBOSAMPLE_SLOW_TESTS`; all wired into `nox -s tests`.

## 0. Relationship to existing Fixman tests (no duplication)

- `tests/TestFixmanBoltzmann.cpp` proves the *general* statement (OFF marginal ∝
  √det M, ON marginal flat) on a **synthetic** two-torsion chain, self-referencing
  `calcLogDetM`. That harness runs the C++ `HmcDriver` at a very large sample count,
  so its integrated-autocorrelation inflation is negligible; the Python full-stack
  tiers here do NOT have that luxury and SHALL apply §5's independence protocol.
- `docs/specs/ensemble-validation/30-tier2-*` validates torsion marginals of
  real small molecules against OpenMM MD. This spec is complementary: idealized
  beads with a closed-form metric, plus the cyclic (loop-closure) Fixman path.
- `docs/specs/singular-dof-fixman.md` governs the null-`D_b` edge case. The chains
  here are chosen to have full-rank M (assert `minEigD` bounded away from the lock,
  §6); they are not a phantom-DOF test.

## 1. Problem

The Fixman compensating potential removes the mass-metric bias of constrained
(torsional) sampling. The current suite validates its *behaviour* self-referentially
but pins it against **no external number**. Jain 2013 (eq:24) gives a closed-form
det M(α) for a specific C4 parametrization — an analytical oracle independent of the
engine — and a family of serial chains (C5/C11/C15) whose interior torsions carry a
genuine, conformation-dependent metric bias. Reproducing these is the strongest,
cheapest correctness anchor available for the Fixman metric.

Two premises in the pre-revision spec were empirically falsified and are corrected
here:

1. **det M is parametrization-dependent (§2), so the eq:24 oracle does NOT transfer
   through the full-stack pipeline.** The full-stack C4 uses a different decomposition
   whose det M(α) is flat; the eq:24 shape exists only in the Tier-0 hand-built 2+2
   parametrization. A Tier-1 test asserting the full-stack C4 TORSIONAL histogram
   matches eq:24 is physically impossible (§3.1 NOTE, §3.3, §5 Tier-1).
2. **HMC output is autocorrelated**, so the χ² goodness-of-fit gate inherited from
   `TestFixmanBoltzmann` is invalid when applied to raw per-round frames: the
   measured integrated autocorrelation time (IACT) at mdSteps=20/ts=0.001 is
   τ_int≈362 for the C4 torsion (~30 independent samples out of 12000 frames),
   inflating χ² into the thousands and rejecting every hypothesis including the true
   one (§5 statistical protocol).

## 2. Binding definitions

Vocabulary follows `ensemble-validation/00-foundations.md` §1–3 (worlds, the true
metric `M(φ) = Jᵀ M_cart J`, `calcLogDetM = ln det M(φ)`, the `+1/2`-power
convention). Terms specific to this spec:

- **Generalized-coordinate marginal** (Spiridon 2017 eq:2): `ρ(φ) ∝ |M(φ)|^{1/2}
  e^{-βU}`. Code binding: `World::reinitialize` (src/World.cpp ~:1596) seeds
  `u = √(RT)·M(q)^{-1/2}·g`, `g~N(0,I)`, so the sampled momentum is `p = M u ~
  N(0, RT·M(q))`; marginalizing `p` yields exactly this `√det M` configurational
  weight. This is the metric-aware draw that *produces* the bias Fixman removes.
- **Fixman compensating potential** (Spiridon eq:3, Jain 2013 eq:11/13):
  `U_F(φ_f) = ½ β⁻¹ ln( |M_{N_f}| / |M_{3N}| )`. Code binding: `World::calcFixman`
  (src/World.cpp:835, return at :864) computes
  `U_F = ½RT( ln|M_tree| − ln det(GM⁻¹Gᵀ) − ln|M_3N| )`. The sign of the tree term
  is **+** (energy is raised where det M is large), so `e^{-βU_F} ∝ 1/√det M_tree`
  cancels the marginal's `√det M`. For torsion-only `φ_f` with fixed bonds/angles
  `|M_{3N}|` is torsion-independent (Go-Scheraga; Jain eq:7), so the corrected
  marginal is `∝ e^{-βU}`.
- **Parametrization dependence of det M (load-bearing).** `det M(φ) = det(Jᵀ M_cart J)`
  depends on the *choice of generalized coordinates* (which tree decomposition and
  joint frames realize the same reachable Cartesian manifold), not only on the
  reachable configurations. The same physical chain decomposed two ways has two
  different `det M(φ)`, hence two different biased marginals. Consequences: (i) an
  external `det M` oracle such as eq:24 is a statement about a **specific**
  parametrization and does not carry across decompositions; (ii) whether a given
  torsion's marginal is biased is a property of the decomposition, and SHALL be
  established empirically for the full-stack pipeline rather than assumed (§5).
- **The three models** (each a Robosample world configuration):
  - **FLEXIBLE** — Cartesian/OpenMM world, harmonic bonds+angles, no constraints,
    `use_fixman` off (flat metric). The reference ensemble.
  - **TORSIONAL** — torsional world, `use_fixman=False`. Marginal `∝ √det M(φ)`; the
    *biased* distribution to be corrected. Non-flat only where `det M` is
    conformation-dependent (§2 parametrization NOTE, §5).
  - **FIXMAN** — torsional world, `use_fixman=True`. Marginal recovers FLEXIBLE.
- **Loop-closure Fixman term** (Spiridon 2017 via `Constraints::calcConstraintLogDet`):
  for a cyclic molecule with `m` loop constraints, `|M_{N_f}| = |M_tree| /
  det(G M⁻¹ Gᵀ)`, so `U_F` acquires `−½ β⁻¹ ln det(G M⁻¹ Gᵀ)` on top of the tree
  term. `m = 0` (acyclic) ⇒ this term is exactly 0. Active only in Tier 2.

## 3. Systems (all vacuum)

### 3.1 Idealized serial chains (acyclic; Jain 2013 §III.A, Spiridon 2017 §3.1)

United-atom CH₂ beads, single bead type, no nonbonded, no torsional force-field term.

| System | Beads | Bonds | Torsions | Bond angle |
|---|---|---|---|---|
| C4 | 4 | 3 | 1 | **90°** |
| C5 | 5 | 4 | 2 | **109°** |
| C11 | 11 | 10 | 8 | 109° |
| C15 | 15 | 14 | 12 | 109° |

Common parameters (Jain 2013 checks; Spiridon 2017 §3.1):
- Bead mass **14 amu**. Bond length **1.54 Å**. Bond spring **83.66 kcal/mol/Å²**,
  angle spring **43.46 kcal/mol**.
- **Force field (SHALL match the papers exactly).** The potential is harmonic
  **bonds + angles only**. No nonbonded (no vdW, no Coulomb), no torsion term, no
  1-4. The OpenMM `System` `load_amber` builds SHALL carry non-zero terms **only**
  in `HarmonicBondForce` + `HarmonicAngleForce`; `PeriodicTorsionForce` and
  `NonbondedForce` terms SHALL be identically zero. Consequence — the free torsional
  DOF see `U = 0`, so the TORSIONAL marginal is the **pure metric** `√det M(φ)` and
  FIXMAN targets it. Whether that marginal is actually non-flat depends on the
  decomposition (see NOTE below and §5): it is non-flat only for torsions whose
  `det M` is conformation-dependent.
- **NOTE (C4 metric is flat in the full-stack decomposition — verified, load-bearing).**
  In the shipped-prmtop full-stack decomposition, C4's single central torsion rotates
  a rigid outboard group whose moment of inertia about the bond axis is invariant
  under rotation about that same axis; well-mixed runs confirm `det M(α)` is
  α-independent there (TORSIONAL is flat within sampling noise and does NOT match
  eq:24). The α-dependent `det M(α) = c₅·(35 + 4cos α − 16cos²α + cos⁴α)` (eq:24)
  arises only in the deliberately different **2+2** hand-built decomposition of Tier-0
  (root owns 2 beads, leaf owns 2 — see `TestFixmanIdealizedChains.cpp`
  `buildC4RootPairLeafPair`). C4 therefore does **not** exhibit the metric bias
  through the full stack, and FIXMAN is a genuine no-op for it. The bias Fixman
  removes lives in the **interior torsions of C5+** (§5 Tier-1). This is exactly §2's
  parametrization-dependence: eq:24 is a property of the Tier-0 parametrization, not
  of the molecule.
- NOTE (which paper). We follow **Jain 2013 + Spiridon 2017** for the chains
  (springs 83.66 / 43.46, no nonbonded); not Kandel 2016's C4 variant.
- NOTE (unit conventions). Bond/angle spring constants matter only for the FLEXIBLE
  reference (Tier 1); the deterministic det M (Tier 0) depends only on masses, bond
  lengths, and bond angles at the frozen geometry. The ½-prefactor convention on the
  angle spring differs between sources; reconcile at implementation.
- NOTE (mass correction). The shipped C4 SHALL use bead mass 14 (det M's α-shape is
  mass-independent; mass enters only the constant `c₅` of eq:24).

### 3.2 Proline dipeptide (realistic cyclic; Kandel 2016 §III)

- **Ace-Pro-Nme** (one proline) is the smallest reproducible case; the shipped Tier-2
  example is the richer **Ace-Pro-Pro-Nme** (di-proline, two pyrrolidine rings → two
  loop constraints), which §3.2 already permits as a valid richer variant. Either
  satisfies the binding cyclic-path invariants (T2.0). Force field **AMBER ff99SB**
  applied under the **Spiridon 2017 §3.3 vacuum peptide protocol**: vacuum, `NoCutoff`,
  no periodic box, 300 K.
- NOTE (force-field decision). No paper ran proline in pure vacuum; Kandel 2016
  solvated it with GBSA. We take Kandel's molecule + ff (ff99SB) and Spiridon's
  vacuum protocol, dropping GBSA (an extrinsic solvent coupling Fixman cannot correct,
  Kandel §V). What is binding is that FLEXIBLE and FIXMAN run the **identical** force
  field; ff14SB/ff12SB are acceptable substitutes.
- Each pyrrolidine ring closes via the CD–N bond, classified `PROTEIN_RING_DIHEDRAL`
  and cut into a loop-closure `DistanceConstraint` (`World.cpp:679`). So FIXMAN here
  exercises **both** the tree term (`calcLogDetM`) and the loop term
  (`calcConstraintLogDet`) on a real molecule.

## 3.3 Test surface (binding — the C++ harness cannot load examples raw)

Building a model from an AMBER file goes through the Python `Context.load_amber`
pipeline (needs OpenMM); the Python API exposes no per-config det M / Fixman getter
(only `calc_openmm_potential_energy` and `atoms_x`). This splits the surface:

- **Tier 0 (deterministic det M) SHALL be C++**, on a **hand-built `RobotModel`**
  parameterized to §3.1 — precedent `TestFixmanBoltzmann.cpp` — because only C++
  reaches `RobotEngine::calcLogDetM`. It does **not** load `examples/fixman/*`.
- **Tier 1 / Tier 2 (statistical) SHALL be Python full-stack**, driving the shipped
  `examples/fixman/*` through the `run.py` path, histogramming torsions from
  `atoms_x`/DCD.
- **NOTE (Tier-0 ↔ Tier-1 relationship — corrected).** Tier-0 (C++, 2+2 hand-build)
  is the **only** place the eq:24 `det M(α)` oracle lives; it validates `calcLogDetM`
  against the closed form at machine precision. The full-stack C4 (Tier-1) uses a
  **different** decomposition whose `det M(α)` is flat (§3.1 NOTE), so the Tier-1 C4
  marginal **cannot** be compared against eq:24 — that comparison is physically
  impossible and any test asserting it can never pass. The tiers are tied not through
  eq:24 but through (i) identical C4 masses/geometry (§3.1) and (ii) the shared engine
  routine `calcLogDetM`. Drift between the hand-built and shipped C4 geometries
  surfaces as a Tier-0 `minEigD`/vacuity failure or a Tier-1 flatness failure, not as
  an eq:24 mismatch.

## 4. Example files and naming

- `examples/fixman/c4/`, `c5/`, `c11/`, `c15/` — each ships `<name>.mol2`,
  `<name>.frcmod`, `<name>.prmtop`, `<name>.rst7`, `README.md`, and the `tleap` input
  so the prmtop/rst7 regenerate deterministically. `c4/` ships with bead mass 14.
- `examples/fixman/` proline: `ace_pro_pro_nme.prmtop`/`.rst7` (shipped di-proline;
  Ace-Pro-Nme acceptable), `README.md`, built from the AMBER sequence.

## 5. Test plan — simplest to most complex

### Tier 0 — deterministic metric identities (always-on, no sampling) — C++

Unchanged and passing (`tests/TestFixmanIdealizedChains.cpp`). The keystone tier:
closed-form ground truth at machine precision, no simulator. Surface: hand-built
`RobotModel` (§3.3). This is the **only** tier that exercises eq:24.

- **T0.1 — C4 analytical det M (Jain 2013 eq:24 / Fig 1a).** SHALL. On the 2+2
  hand-built parametrization (`buildC4RootPairLeafPair`), over α ∈ [−π, π] assert
  `det M(α) = exp(calcLogDetM) ≈ c₅·(35 + 4cos α − 16cos²α + cos⁴α)` (relative
  residual below §6), plus the predicted extrema (√det M maximal at α ≈ ±82.8°, local
  minima at α = 0, ±180°, ratio f(0):f(180) = 24:16). NOTE the α convention is
  `α = π − dihedral(A,B,C,D)` (resolved in the C++ banner). This oracle is
  parametrization-specific (§2) and stays in C++.
- **T0.3 — multi-torsion determinant identity, C5/C11/C15.** SHALL. Assert the O(n)
  factorization `ln det M = Σ_b ln D_b` agrees with a dense `ln det(Jᵀ M_cart J)` to
  float/FD tolerance at several random torsion configurations. Guard:
  `calcConstraintLogDet == 0` exactly for every acyclic chain.
- **T0.2 / T0.4** (SHOULD, secondary): excluded with justification in the C++ file
  banner (unreachable / non-existent Fixman-torque path). Not required.
- Vacuity + `minEigD` guards as coded.

### Tier-1/Tier-2 statistical protocol (SHALL — applies to every statistical assertion below)

HMC output is autocorrelated; the goodness-of-fit gates are valid only on
(approximately) independent samples. The pre-revision recipe (χ² on raw per-round
frames at mdSteps=20/ts=0.001) is **invalid** and SHALL be replaced.

- **P1 (IACT).** SHALL estimate the integrated autocorrelation time `τ_int` of each
  analyzed torsion time series from the production trace (automated-windowing / Sokal
  or an emcee-style estimator). Effective sample size `N_eff = N_frames / τ_int`.
- **P2 (mixing).** SHALL choose the HMC proposal length (`mdSteps × timeStep`) so
  `τ_int` is O(10) rounds. NOTE mdSteps=20/ts=0.001 gives `τ_int≈362` for the C4
  torsion (`N_eff≈33` of 12000 frames) — unusable; mdSteps=100/ts=0.002 gives
  `τ_int≈10` (`N_eff≈1150`, 35× better). SHOULD prefer longer MD trajectories (larger
  `mdSteps`) over writing every round.
- **P3 (thinning).** SHALL thin the trace by `⌈τ_int⌉` (one retained frame per IACT)
  before computing any χ² statistic, OR use an `N_eff`-corrected statistic. Raw-frame
  χ² inflates by ~`τ_int` and rejects even the true hypothesis.
- **P4 (power / sizing).** SHALL size the run so the post-thinning independent count
  `N_indep` separates the flat and biased hypotheses at α: from a well-mixed
  calibration run, measure the effect size (e.g. `Hellinger(TORSIONAL, flat)` on a
  signal torsion) and require `N_indep` large enough that χ² rejects flat with power
  ≥ 0.99 while FIXMAN (flat) is not rejected. SHOULD target `N_indep ≳ 1000` (the
  calibration regime that gave `H(FIX,flat)=0.016`).
- **P5 (χ² gate).** With independent samples, apply χ² goodness-of-fit at α = 1e-4,
  (nbins − 1) dof (precedent `TestFixmanBoltzmann`). NOTE: where the run cannot afford
  `N_indep` sized for the χ² gate, a noise-floor-relative Hellinger bound
  (`H < k·√((nbins−1)/(8 N_indep))`) is an acceptable weaker surrogate for
  "flat within sampling noise"; it SHALL NOT be presented as a full-power χ² result.

### Tier 1 — statistical distribution recovery, idealized chains (`ROBOSAMPLE_SLOW_TESTS`) — Python

Python full-stack on the shipped `examples/fixman/*` prmtops (§3.3). No external
formula is used in Tier 1; the eq:24 oracle stays in Tier 0 (§3.3 NOTE). All
assertions run on independent samples per the protocol above.

- **T1.1 — C4 null / no-op control (SHOULD; recommended over drop).** The full-stack
  C4 metric is flat (§3.1 NOTE), so Fixman is a genuine no-op. If C4 is kept: SHALL,
  on independent samples of the central dihedral, assert TORSIONAL, FIXMAN, and
  FLEXIBLE marginals each **fail to reject flat** at α, and are pairwise close
  (`Hellinger` below a calibrated floor). This encodes the verified fact "Fixman does
  not introduce bias when the metric is already flat" and provides the
  FLEXIBLE-is-flat baseline and a cheap end-to-end pipeline smoke.
  - **NOTE (limited power — do not mistake for the primary check).** C4 does **not**
    test bias REMOVAL and **cannot** catch a Fixman SIGN error: with `det M` constant,
    `U_F` is constant for either sign, so a flipped-sign `calcFixman` is invisible
    here. Bias-removal and sign correctness are tested by Tier-0 (eq:24) and Tier-1
    C5+. C4 MAY be dropped if the reviewer prefers to exercise the pipeline via C5+.
  - **Recommendation:** repurpose C4 as this null control rather than drop it — it is
    the cheapest full-pipeline test, it converts a surprising empirical finding into a
    regression guard, and it anchors the FLEXIBLE baseline; its limitation is stated
    so it is not read as the correctness keystone.

- **Signal-torsion selection (PRECONDITION / guard for T1.2, T1.3).** SHALL identify
  signal-carrying torsions **empirically**: a torsion carries signal iff its TORSIONAL
  marginal (independent samples) is non-flat beyond the noise floor. Terminal / no-op
  torsions (TORSIONAL already flat) are EXCLUDED from the Fixman-flattening comparison
  and MAY be recorded as additional null controls. NOTE (vacuity handling —
  reconciled with the "do not re-run" constraint of the current implementation round,
  Open §7.2): if a chain yields NO signal torsion, it is a pure null control; the test
  SHALL still assert every marginal (TORSIONAL and FIXMAN) is flat and SHALL emit a
  loud recorded warning that the chain exercised no metric bias (so the vacuity is
  visible and traceable to §7.2), rather than silently passing as if bias removal had
  been demonstrated. Whether any full-stack chain in fact carries signal is the open
  question §7.2 and is not yet confirmed.
  - **NOTE (physical guidance, SHOULD).** `det M` is conformation-dependent for a
    torsion when rotating it reshapes the outboard mass distribution seen by another
    flexible hinge; this needs ≥ 2 coupled flexible torsions with non-collinear axes
    moving multi-atom groups on both effective sides of the axis. Interior torsions of
    the longer chains (C11/C15) carry the strongest signal; a torsion that moves only
    a single terminal atom, and the innermost torsion (no inboard hinge to reshape),
    are near/exact no-ops. The exact per-torsion signal map depends on how the full
    stack roots the tree and whether the free 6-DOF base couples into the internal
    marginals (Open §7) — hence the binding selection is empirical, not index-based.

- **T1.2 — C5 (Patriciu 2004).** SHALL on each signal torsion (if any exist); SHOULD
  across all. On each signal torsion, assert (i) FIXMAN is flat within noise at α,
  (ii) TORSIONAL is non-flat, (iii) `Hellinger(FIXMAN, flat) < Hellinger(TORSIONAL,
  flat)` strictly. NOTE (iii) is the sign-discriminating assertion: with a flipped
  Fixman sign, FIXMAN would be *more* biased than TORSIONAL and (i)/(iii) would fail —
  the check C4 cannot make. When no signal torsion exists, apply the vacuity handling
  above.
- **T1.3 — C11 / C15 capstone.** SHOULD / cost-gated. Same three assertions on the
  interior signal torsions (majority-based across the 8–12 torsions so a few noisy
  simultaneous marginals cannot make it flaky). C15 is the heaviest; C11 MAY be
  dropped if C4+C5+C15 already exercise the scaling.

### Tier 2 — proline dipeptide realism (`ROBOSAMPLE_SLOW_TESTS`, vacuum) — Python

Python full-stack on the shipped proline example (§3.2), vacuum ff99SB.

- **T2.0 — cyclic-path structural guard (always-on).** SHALL. Assert the molecule
  yields the expected loop-closure count (≥ 1 `DistanceConstraint`; exactly 2 for the
  shipped di-proline) and that `current_constraint_log_det()` is non-zero and
  conformation-dependent — otherwise the cyclic Fixman path is untested/vacuous.
- **T2.1 — reference recovery (slow).** Build FLEXIBLE and FIXMAN (torsional +
  tree-Fixman + loop-Fixman) in the **same** vacuum force field; sample backbone
  (φ, ψ) and proline χ torsions. SHALL apply the §5 statistical protocol (P1–P4:
  estimate IACT, thin to independent samples, size `N_indep`). Assert:
  - SHALL: FIXMAN closer to FLEXIBLE than TORSIONAL is
    (`Hellinger(FIX, FLEX) < Hellinger(TOR, FLEX)`), the qualitative Kandel 2016
    ordering, judged on the majority of INFORMATIVE torsions (those whose FLEXIBLE
    reference actually explored) and in the mean. NOTE the previous measurement was at
    `N_eff≈30` (autocorrelation-dominated and non-discriminating); the ordering SHALL
    be recomputed on independent samples.
  - SHOULD: `Hellinger(FIX, FLEX)` below an absolute threshold. The threshold SHALL be
    **recalibrated from a well-mixed run** (`τ_int` O(10)) and recorded in the test;
    the pre-revision `0.6` bound was set from an under-mixed run and is not
    authoritative (it has been removed pending recalibration).
  - NOTE (why not exact). Vacuum still carries 1-4/nonbonded coupling between
    constrained and flexible DOF — an extrinsic distortion Fixman cannot remove
    (Kandel 2016 §V). Bounded closeness with the correct ordering is the target, not
    exact FIXMAN=FLEXIBLE recovery.

## 6. Tolerances

- Tier 0: relative residuals `≤ 1e-6` engine-vs-engine (T0.3 block-vs-dense) and
  `≤ 1e-4` where a central FD is involved; T0.1 fit residual `≤ 1e-6` (exact geometry).
- `minEigD` bounded away from the singular-DOF lock for every chain configuration.
- Tier 1 χ²: critical value at α = 1e-4, `nbins − 1` dof, computed on **independent
  (thinned) samples only** (§5 P3); `N_indep` sized per §5 P4. Raw-per-round-frame χ²
  is invalid and SHALL NOT be used. The noise-floor Hellinger surrogate (§5 P5 NOTE)
  uses `k = 3` and a `0.12` absolute floor.
- Tier 2: the SHALL ordering carries no free tolerance; the SHOULD absolute Hellinger
  threshold is calibrated on independent samples from a well-mixed run and recorded
  in-test.

## 7. Open decisions for the reviewer / owner

1. **C4 full-stack role (RESOLVED here).** Repurpose as a null / no-op control
   (recommended, §5 T1.1) or drop. The eq:24 oracle is Tier-0-only; the pre-revision
   T1.1 "C4 TORSIONAL matches eq:24" is removed as physically impossible.
2. **Exact per-torsion signal map + free-vs-welded base — is any full-stack torsion
   actually biased? (OPEN — the one unresolved correctness question).** The
   rooted-tree factorization `ln det M = Σ_b ln D_b` predicts the innermost torsion is
   a metric no-op (no inboard hinge whose articulated inertia its rotation changes) and
   that non-innermost torsions carry signal; but whether the free 6-DOF base couples
   into the sampled internal marginals is unresolved — empirically it appears not to
   (full-stack C4, a single free-based torsion, is flat, whereas the C++ free-based
   2+2 build is not). **If it turns out that the full stack floats the base and every
   chain's torsions are flat (like C4), then the full-stack Fixman validation has no
   bias to remove and Tier-1 is vacuous — a real finding that would need engine
   investigation, not a test fix.** This governs which torsions the Tier-1 SHALL
   assertions cover but is deliberately NOT wired as a hard failure in the current
   round (the slow tier was not re-run, so signal presence is unconfirmed); §5's
   empirical selection + loud vacuity warning surfaces it either way. To resolve:
   instrument per-torsion TORSIONAL non-flatness in a well-mixed calibration run, or
   expose the per-body `D_b` so the signal map can be read deterministically. What
   would settle it: whether the full-stack torsional world welds the molecule base or
   floats it, and whether `calcFixman`'s `calcLogDetM` includes the free-root 6×6
   block in the sampled marginal.
3. **Proline case / force field.** Ace-Pro-Pro-Nme (shipped di-proline, 2 loops) vs
   Ace-Pro-Nme (1 loop); ff99SB vacuum vs Kandel's solvated GBSA. Binding invariant:
   FLEXIBLE and FIXMAN share one ff, and ≥ 1 conformation-dependent loop constraint
   fires.
4. **Tier-1 depth.** Whether C11 and C15 both run, or C15 alone as the capstone.

## 8. Non-goals

- Flexible cyclic-bead systems with a closed-form loop metric — no such reference
  exists in Jain/Patriciu.
- Branched acyclic peptide R_tor/R_Fix tables with `calcConstraintLogDet ≡ 0`.
- Solvent / implicit-solvent Fixman — extrinsic distortion outside a metric-correction
  validation.
