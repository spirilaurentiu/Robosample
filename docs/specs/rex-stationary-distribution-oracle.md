# Spec: Replica-exchange stationary-distribution oracle (REMC/Default driver-level)

Status: **draft** — researcher artifact, ready for review. Author: derived from a full read of the
label-swap REX driver (`src/Context.cpp`, `include/Context.hpp`, `include/ReplicaExchange.hpp`), the
parent REX spec (`docs/specs/replica-exchange-nonequilibrium-work.md`), the existing acceptance-algebra
tests (`tests/TestRexAcceptanceAlgebra.cpp`), and the OpenMM-free ensemble-validation machinery
(`tests/HmcDriver.hpp`, `tests/AnalyticForceBridge.hpp`, `tests/RobotBuilders.hpp`, `tests/StatTest.hpp`,
`tests/TestFixmanBoltzmann.cpp`, `tests/TestEnsembleValidation.cpp`).

Requirement keywords follow RFC 2119 (`styles/spec.md`): SHALL / SHOULD / MAY / NOTE.

Scope: this spec adds ONE new C++ gtest translation unit, `tests/TestRexStationaryDistribution.cpp`, plus
one small test-infrastructure mutator. It changes **no** engine source. It operationalizes invariant V2
("REMC marginal") of `docs/specs/replica-exchange-nonequilibrium-work.md` into a concrete, cheap,
closed-form oracle at the level of the exchange **driver's control structure**, not just the per-swap
acceptance algebra.

---

## Motivation

### Problem restatement (user's phrasing)

> Robosample's replica-exchange (REX) test coverage is currently *algebraic only* — it verifies the swap
> acceptance formula in isolation but never verifies that the REX *driver*, end to end, samples the
> correct distribution. We need a theory-derived behavioral oracle: a test proving the REX chain leaves
> the correct product-Boltzmann distribution invariant (detailed balance / stationary distribution) at
> the driver level.

### The gap in the codebase's vocabulary

`tests/TestRexAcceptanceAlgebra.cpp` (the "algebra" coverage) evaluates one swap move on fixed inputs: for
a chosen configuration pair and inverse temperatures it asserts the acceptance exponent `log α`
(`attemptREXSwap`'s `ETerm_equal` / `ETerm_nonequil` / `WTerm` branches, `src/Context.cpp:667-706`)
satisfies the detailed-balance identity `π(z)·P(z→Tz) = π(Tz)·P(Tz→z)` against a closed-form `π`, and that
a deliberately corrupted β-assignment breaks it. That suite never runs a Markov chain. It exercises none
of the driver machinery that composes swaps into a sampler:

- the inner GC-HMC sweep that propagates each replica between swap attempts (`src/Context.cpp:1112-1143`),
- the two mutually-inverse label maps `replica2ThermoIxs_` / `thermo2ReplicaIxs_` and their update on
  accept (`swapThermodynamicStates`, `src/Context.cpp:626-631`),
- the alternating-parity neighbour pairing (`prepareExchangePairs`, `src/Context.cpp:756-763`) driven by
  the dedicated `exchangeRound_` counter (`mixReplicas`, `src/Context.cpp:781-801`),
- the by-thermodynamic-state (not by-replica-identity) recording convention
  (`src/Context.cpp:1158-1163`, `writeOutputsCore`, `src/Context.cpp:1188-1197`).

A per-move acceptance identity is necessary but not sufficient for driver correctness. A driver that (a)
half-updates the inverse maps, (b) freezes the pairing parity (the B7 defect the parent spec fixes), (c)
records output indexed by replica identity instead of thermodynamic state, or (d) never actually attempts
a swap (the F1 "no swaps happen" regression), can pass every algebra test yet sample the wrong
distribution or fail to mix. No current test can distinguish those.

### What this oracle asserts

Let the outer chain run on a system whose per-temperature configurational Boltzmann marginal is known in
**closed form**. The oracle asserts, at the driver level, the two properties the algebra tests cannot
reach:

1. **Stationary distribution (marginal correctness).** Every thermodynamic state's sampled configurational
   marginal equals the analytic Boltzmann marginal at that state's temperature — measured via the existing
   `StatTest` chi-square / G machinery against a closed-form target, not self-consistency.
2. **Ergodic mixing over the ladder (swap-rate / round-trip).** Swaps actually occur, at a nonzero rate
   consistent with the closed-form energy-overlap, and replicas traverse the full temperature ladder
   (round trips), which the marginal check alone cannot detect.

### Binding constraints and the architectural limitation (stated up front)

- **The production `Context::RunREX` driver cannot be instantiated in this oracle.** It requires a fully
  OpenMM-backed `World`/`Context` (energy is a process-wide OpenMM singleton, `OpenMMContext::get()`;
  `openmmPotential`, `src/Context.cpp:1170`), which the analytic-bridge test suite deliberately does not
  build (`HmcDriver.hpp` header: "The C++ test suite never instantiates World"). Moreover the label-swap
  `RunREX` and its helpers "[have not] been compiled or run" (`include/Context.hpp:104-109`;
  `include/ReplicaExchange.hpp:11-15`). Therefore this oracle SHALL reimplement the REMC/Default **control
  structure** on the analytic bridge, mirroring the four production functions line-for-line (see Interface
  I2), exactly as `TestRexAcceptanceAlgebra.cpp` reimplements the acceptance switch and
  `TestFixmanBoltzmann.cpp` reimplements the World HMC kernel via `HmcDriver`. This is the suite's
  established convention for closed-form analytic oracles. The oracle validates the **algorithm** (label-
  swap REMC composed with inner GC-HMC preserves product-Boltzmann) and is the executable reference the
  production `Context::RunREX` SHALL match once it compiles against OpenMM. It is NOT a direct
  instantiation test of `Context::RunREX`; that follow-up is recorded in Open questions.
- Correctness precedence (`CLAUDE.md`): published theory > project specs > tests > implementation. The
  closed-form target is fixed by theory (Sugita–Okamoto, Chodera–Shirts, Kofke); the mirrored control
  structure is fixed by `src/Context.cpp`.
- Cost: the oracle SHALL be cheap enough for the slow tier — an analytic 1–3 DOF-per-body harmonic system,
  no OpenMM, no MD force field. It SHALL gate its high-N statistical assertions behind
  `ROBOSAMPLE_SLOW_TESTS` (`StatTest.hpp:31-37`), with an always-on smoke, matching every other ensemble
  test.

---

## Behavior

### The analytic system (closed-form configurational marginal)

**Reference system SR.** A `RobotModel` built with `RobotBuilders::buildForest` consisting of `n_c`
purely **translational** generalized coordinates rooted on Ground, together with the
`AnalyticForceBridge` isotropic harmonic anchor potential. Reference choice: `n_c = 6` from **two
`JointType::Cartesian` bodies** (each 3 translational DOFs, `RobotModel::jointNU(Cartesian) == 3`, no
quaternion), each carrying at least one real-mass atom via `RobotBuilders::attachAtoms`. A forest of six
`JointType::Slider` bodies (`jointNU(Slider) == 1`) is an acceptable equivalent; both give
configuration-independent mass metric and exactly-quadratic potential.

Claims about SR (proved in the derivation sketch):

- **C1.** With the anchor captured at the initial configuration (`q = 0`, `AnalyticForceBridge`
  constructor, `tests/AnalyticForceBridge.hpp:43-52`), the potential is **exactly quadratic and positive
  definite**: `U(q) = ½ qᵀ A q`, `A ≻ 0`, with no linear term (minimum at `q = 0`).
- **C2.** The mass metric `M` is **configuration-independent** (translations do not rotate the slide
  axes), so `det M` is constant; the GC-HMC configurational marginal is `∝ exp(−βU(q))` with **no Fixman
  term required** (`HmcDriver::useFixman = false`). NOTE: SR SHALL use only translational DOFs precisely
  so that `U` stays exactly quadratic; any rotational DOF (Torsion/BendStretch/Ball/Free) makes
  `U(q)` anharmonic through the Cartesian anchor and destroys the closed form.
- **C3.** At temperature `T_k` the Boltzmann distribution of the total potential energy is
  `U ~ Gamma(shape = n_c/2, scale = RT_k)` with `RT_k = k_B T_k`, `k_B = 0.0083144626` kJ/mol/K
  (`kBoltzmann_kJ`, `src/Context.cpp:21`; `kB`, `tests/TestEnsembleValidation.cpp:44`). Equivalently
  `2βU ~ χ²(n_c)`. This is the exact closed-form per-state marginal target and holds for **any** `A ≻ 0`
  (independent of masses/stiffnesses).

### The distribution the REMC chain SHALL leave invariant

The oracle runs `R = T` replicas, each a persistent copy of SR (its own `RobotState`, config `x_r`),
assigned to thermodynamic states by a permutation `σ` realized through the two inverse maps. The joint
target on `(x_1,…,x_R, σ)` is the product of per-temperature Boltzmann measures:

```
π(x, σ)  ∝  Π_r exp( −β_{σ(r)} U(x_r) ),     β_k = 1/(k_B T_k)
```

(Sugita–Okamoto eq:7 generalized-ensemble product weight; Chodera–Shirts eq:10 permutation-conditional
`π(S|X)`). Its per-state marginal — the distribution of "the configuration currently occupying state
`k`", i.e. `x_{thermo2ReplicaIxs_[k]}` — is exactly `Boltzmann(β_k)`, hence by C3 the observable

```
U at state k   ~   Gamma(n_c/2, RT_k).
```

This is the oracle's closed-form target for every state, with **no free parameters** and no reliance on a
long reference run.

### The mirrored REMC/Default control structure (one round)

The oracle reimplements, on SR, exactly the REMC/Default round of `Context::RunREX`
(`src/Context.cpp:1078-1164`):

1. **Inner Gibbs sweep.** For each thermodynamic state `k` in `[0, T)`, let `r = thermo2ReplicaIxs_[k]`;
   set replica `r`'s inner-move temperature to `RT_k` and run `nInner` `HmcDriver::move()` steps on its
   persistent `RobotState`. Mirrors `w.setTemperature(st.temperature)` + `generateSample()`
   (`src/Context.cpp:1118-1123`). The move is **persistent-configuration** (few HMC steps), never i.i.d.
   resampling — see the derivation-sketch NOTE on why exact resampling would blind the oracle to swap
   defects.
2. **Refresh reference potentials.** For each replica `r`, set `referencePotential = U(x_r)` (mirrors
   `src/Context.cpp:1149-1153`; here `U` is `AnalyticForceBridge::calcPotentialEnergy`, not
   `openmmPotential`).
3. **Mix (REMC only).** `mixReplicas(round)`: if `round % swapEvery == 0` and `T > 1`,
   `prepareExchangePairs(exchangeRound, oddity = 0)` builds neighbour pairs `startIdx = exchangeRound % 2`,
   `(startIdx, startIdx+1), (startIdx+2, startIdx+3), …`; for each state pair `(a, b)` call the swap; then
   `++exchangeRound`. **Parity SHALL come from the dedicated `exchangeRound` counter, not `round % 2`
   after the `swapEvery` gate** (the B7 revision-2 fix, `src/Context.cpp:789-797`). For `Default`, step 3
   is skipped entirely (`mixReplicas` returns for `RUN_TYPE::Default`, `src/Context.cpp:786`).
4. **Swap (REMC).** For state pair `(a=C, b=H)`: `X = thermo2ReplicaIxs_[a]`, `Y = thermo2ReplicaIxs_[b]`,
   `β_C = 1/(k_B T_a)`, `β_H = 1/(k_B T_b)`, `logα = −((β_H − β_C)(U(x_X) − U(x_Y)))` (mirrors the
   `RUN_TYPE::REMC` branch, `src/Context.cpp:648-649,667-669`); accept iff `logα ≥ 0 ∨ u < exp(logα)`
   (`src/Context.cpp:735-736`); on accept, **label-swap only**: `swap(replica2ThermoIxs_[X],
   replica2ThermoIxs_[Y])`, `swap(thermo2ReplicaIxs_[a], thermo2ReplicaIxs_[b])` (mirrors
   `swapThermodynamicStates`, `src/Context.cpp:626-631`); configs never move. REMC/Default never call
   `commitWorkAsFinal` (that is the F4-gated driven path, out of scope).
5. **Record by thermodynamic state.** After a burn-in of `equilRounds`, every `writeFreq` rounds and for
   each state `k`, record `U(x_{thermo2ReplicaIxs_[k]})` into state `k`'s accumulators (mirrors the
   by-state indexing of `src/Context.cpp:1158-1163`). Recording SHALL be indexed by **thermodynamic
   state**, never by replica-object identity — this is the convention that makes the per-state marginal
   equal `Gamma(n_c/2, RT_k)`.

### Derivation sketch

**Product-measure invariance.** The round is a composition of two `π`-invariant Markov kernels:

- *Inner Gibbs update.* For fixed `σ`, `HmcDriver::move()` is a GC-HMC kernel with stationary measure
  `exp(−β_{σ(r)} U(x_r))` for each replica (the standard HMC detailed-balance property already validated
  for this exact kernel by `tests/TestFixmanBoltzmann.cpp` and `tests/TestEnsembleValidation.cpp`).
  Updating each `x_r | σ` leaves `π` invariant (Gibbs update of the configuration block; Chodera–Shirts
  eq:10 conditions on the current permutation).
- *Label-swap update.* The swap proposes a transposition of the labels of the two replicas occupying
  states `a, b` — a deterministic, self-inverse (symmetric) proposal on `σ`. The joint-density ratio is
  `π(σ')/π(σ) = exp(−(β_b − β_a)(U(x_X) − U(x_Y)))`, so the Metropolis acceptance
  `min(1, π(σ')/π(σ))` has `logα = −((β_H − β_C)(U(x_X) − U(x_Y)))` — exactly the `RUN_TYPE::REMC`
  branch. This is the parallel-tempering swap of Swendsen–Wang (via Sugita–Okamoto eq:15 swap cost
  `Δ = (β_n − β_m)(E_i − E_j)`, eq:17 `min{1, exp(−Δ)}`; Chodera–Shirts eq:24 pair-swap acceptance). It
  is a Metropolis update on `σ` and leaves `π` invariant.

A composition of `π`-invariant kernels is `π`-invariant, so the product measure `π(x, σ)` is stationary
for the round, and (with the pairing giving an irreducible chain over permutations) it is the unique
stationary distribution. The per-state marginal is `Boltzmann(β_k)`; by C3 the recorded `U` at state `k`
is `Gamma(n_c/2, RT_k)`. QED sketch.

**Why the inner move SHALL be persistent, not i.i.d. resampling.** If the inner update were exact i.i.d.
resampling from `Boltzmann(β_{σ(r)})`, then at recording time the config at state `k` was just drawn from
`Boltzmann(β_k)` regardless of whether the preceding swap acceptance was correct — the swap becomes
invisible and a broken acceptance passes. A persistent, incompletely-decorrelating inner move (a few HMC
steps) carries configuration history across temperatures, so a wrong acceptance drives configs between
temperatures at wrong rates and biases at least one state's stationary marginal, which is detectable. This
is the same reason `TestRexAcceptanceAlgebra`'s single-replica NCMC case is needed alongside the
symmetric-swap case: the discriminating quantity must be exercised where it is not identically cancelled.

**Closed-form energy-overlap for the swap rate (secondary target).** For two harmonic replicas the mean
swap acceptance is a fixed function of `(n_c, T_a, T_b)` — Kofke's constant-heat-capacity model
(`kofke_2002`: eq:9 gamma-distributed energy for constant heat capacity, eq:1 swap acceptance
`min[1, exp(−(β_0−β_1)(U_1−U_0))]`, eq:10 exact acceptance as a 1-D overlap integral; eq:3/eq:11 show the
acceptance depends on temperatures only through their ratio, justifying a geometric ladder). The oracle
estimates this target to arbitrary precision by direct Monte-Carlo over two **independent** draws
`U_a ~ Gamma(n_c/2, RT_a)`, `U_b ~ Gamma(n_c/2, RT_b)` averaging `min(1, exp(−(β_b−β_a)(U_a−U_b)))` — a
closed-form-marginal oracle, not a self-consistency check.

**No-drive relation to the algebra tests.** At `n_c` translational DOFs with no BAT scaling, there is no
Jacobian and no work term: the driver reduces to plain parallel tempering. This oracle therefore has zero
overlap with the RENE/RENEMC/`WTerm`/Jacobian coverage of `TestRexAcceptanceAlgebra` and with V4/V5/V8/V10
of the parent spec (all of which are driven-path oracles). See Interface I4 for the explicit
non-overlap table.

---

## Invariants

Properties that SHALL hold after the oracle lands (each maps to a tagged check in Validation strategy).

- **INV-M1 (per-state stationary marginal).** Under `RUN_TYPE::REMC` and under `RUN_TYPE::Default`, the
  recorded `U` at each thermodynamic state `k` SHALL be statistically consistent with
  `Gamma(n_c/2, RT_k)` and SHALL be statistically inconsistent with the neighbour state's
  `Gamma(n_c/2, RT_{k±1})`. (Operationalizes V2 of the parent spec.)
- **INV-M2 (swap preserves the marginal).** The per-state marginals under `REMC` (swaps active) SHALL
  match the per-state marginals under `Default` (no swaps); the label-swap composition SHALL NOT bias any
  state's stationary distribution.
- **INV-M3 (map bijectivity).** After every accepted or rejected swap,
  `thermo2ReplicaIxs_[replica2ThermoIxs_[r]] == r` for all `r`, and `replica2ThermoIxs_` /
  `thermo2ReplicaIxs_` remain a permutation and its inverse. A half-updated map is a defect.
- **INV-M4 (swaps occur).** For `REMC` with `T > 1`, the total attempted-swap count SHALL be nonzero and
  the accepted/attempted ratio SHALL lie strictly in `(0, 1)`. Guards the F1 "no swaps happen" regression.
- **INV-M5 (ladder ergodicity).** For `REMC` with the geometric ladder, at least one replica SHALL
  complete a full cold↔hot round trip over the run, and every state SHALL be occupied by at least two
  distinct replicas. A parity-frozen pairing (the B7 defect) SHALL produce zero round trips. `Default`
  SHALL produce zero round trips (control).
- **INV-M6 (by-state recording).** The observable histogrammed for state `k` SHALL be the `U` of the
  replica currently occupying state `k` (`thermo2ReplicaIxs_[k]`), never the `U` of the fixed replica
  object `k`. Recording by replica identity biases the marginal into a temperature mixture.

---

## Interface

### Touch list

- **New file `tests/TestRexStationaryDistribution.cpp`** (C++ gtest). This is the entire deliverable's
  engine-observable surface. No production source changes.
- **Reused test infrastructure (unchanged):** `HmcDriver.hpp` (inner GC-HMC move), `AnalyticForceBridge.hpp`
  (harmonic `U`, closed-form), `RobotBuilders.hpp` (`buildForest`/`attachAtoms`, `BodySpec` with
  `JointType::Cartesian`/`Slider`), `StatTest.hpp` (`Histogram`, `chiSquareStatistic`, `chiSquareCritical`,
  `expectedFromWeights`, `MeanAccumulator`, `slowEnabled`), `TestHelpers.hpp` (`Rng`,
  `warnSlowTierSkipped`).
- **One small test-infrastructure addition:** `HmcDriver` SHALL gain a `void setRT(Real rt)` mutator (its
  `RT_` is private with no setter, `tests/HmcDriver.hpp:220`). The oracle keeps one `HmcDriver` per
  replica bound to that replica's persistent `RobotState` and sets its `RT_` to the current state's
  `RT_k` before each round's inner sweep — mirroring the production `w.setTemperature(st.temperature)` per
  round (`src/Context.cpp:1118`). This is test infra, not production code.
- **Shared closed-form helper (SHOULD):** the `Gamma(shape, scaleKT)` sub-sampled per-bin weight helper
  currently local to `tests/TestEnsembleValidation.cpp:59-75` (`gammaSubsampledWeights`) SHOULD be
  promoted to `StatTest.hpp` so this oracle and the ensemble suite share one implementation; otherwise the
  oracle MAY duplicate the tiny function in-file. Bin-centre evaluation is biased for a curved density at
  large N (documented at `tests/TestEnsembleValidation.cpp:112-114`); the sub-sampled form SHALL be used.

### Conventions at risk (the four production pieces the oracle mirrors verbatim)

The oracle SHALL mirror these line-for-line, each cross-referenced in-file to its `src/Context.cpp`
location, so a future refactor that diverges the driver from the mirror is caught:

- `swapThermodynamicStates` label-swap (`src/Context.cpp:626-631`) — INV-M3, INV-M6.
- `prepareExchangePairs` alternating parity (`src/Context.cpp:756-763`) driven by `exchangeRound`
  (`src/Context.cpp:789-797`) — INV-M5 (the B7 parity fix).
- `attemptREXSwap` `RUN_TYPE::REMC` branch and β definitions (`src/Context.cpp:648-649,667-669,735-736`)
  — INV-M1, INV-M2.
- by-thermodynamic-state recording (`src/Context.cpp:1158-1163`) — INV-M6.

SHOULD (maintainability, not blocking): the pure-integer pairing + label-swap logic SHOULD eventually be
extracted into a header included by both `Context` and this test so the mirror cannot silently drift; this
is a follow-up, not a precondition for landing the oracle.

### I4. Zero-overlap with existing coverage

| Concern | `TestRexAcceptanceAlgebra.cpp` (existing) | This oracle (new) |
| --- | --- | --- |
| Per-swap acceptance formula `logα` | Yes — asserts `π(z)P=π(Tz)P` on fixed inputs, catches β-swap, `correctionTerm` | No — **uses** the formula, does not re-derive it |
| Jacobian / `WTerm` / RENE / RENEMC | Yes — single-replica NCMC, `lnJac` sign | **Out of scope** (no BAT scaling; driven paths uncompiled) |
| Inner GC-HMC sweep + swap composition over many rounds | No | Yes — full REMC/Default chain |
| Label maps stay mutual inverses on accept | No | Yes — INV-M3 |
| Alternating-parity pairing / ladder mixing | No | Yes — INV-M5 round trips |
| By-state vs by-replica recording | No | Yes — INV-M6 |
| Sampled per-state marginal vs closed-form Boltzmann | No | Yes — INV-M1 chi-square vs `Gamma(n_c/2, RT_k)` |
| Swap rate vs closed-form energy-overlap | No | Yes — INV-M4 / secondary LEMMA |

### I5. Scope bound (recorded, non-negotiable)

- **In scope (runnable surface):** `RUN_TYPE::REMC` (parallel tempering) and `RUN_TYPE::Default`
  (independent replicas, no exchange). These are the fully-wired run-type branches
  (`include/ReplicaExchange.hpp:29`, `src/Context.cpp:786,667-669`).
- **Out of scope:** `RUN_TYPE::RENE` / `RUN_TYPE::REBASONTOP` (driven BAT-scaling round-loop is
  reviewed-on-paper, uncompiled — `include/Context.hpp:104-109`, `include/ReplicaExchange.hpp:11-15`) and
  `RUN_TYPE::RENEMC` (its driven round-loop is a Stage-2c TODO; `Context::RunREX` throws for it,
  `src/Context.cpp:1059-1065`). A stationary-distribution oracle for the driven run types SHALL be
  deferred until those round-loops compile and run; their per-swap acceptance algebra is already covered
  by `TestRexAcceptanceAlgebra.cpp`.

---

## Validation strategy

Oracles tagged PRECONDITION / INVARIANT / LEMMA per `CLAUDE.md`. Reference parameters (the coder MAY tune
within the stated bounds, but the discriminating structure is fixed): `n_c = 6`; `T = 4` states on a
geometric ladder, reference `{250, 325, 422.5, 549.25}` K (ratio 1.3) — geometric so the swap acceptance
is roughly uniform across the ladder (`kofke_2002` eq:3/eq:11); `swapEvery = 1`; `nInner ≈ 4`–`8` HMC
moves per replica per round with `HmcDriver` timestep tuned to `≈ 50`–`80 %` inner acceptance;
`equilRounds` a burn-in of `≥ 5×10³` rounds; recorded samples per state `N ≥ 10⁵` (full tier); histogram
of `U` per state with `≈ 24`–`30` bins over `[0, c·RT_k]`, `c ≈ 12`; goodness-of-fit at `α = 1e-4`
(`chiSquareCritical`). Smoke tier: `≈ 2×10³` rounds, always on.

- **O1 (PRECONDITION) — map bijectivity (INV-M3).** After every swap attempt, assert
  `thermo2ReplicaIxs_[replica2ThermoIxs_[r]] == r` for all `r` and that both arrays are permutations. A
  runtime guard in the mirrored swap, not a statistical test. Catches a half-updated `swapThermodynamicStates`.
- **O2 (PRECONDITION) — swaps occur (INV-M4).** For `REMC`, `T = 4`: total attempted-swap count `> 0` and
  `0 < accepted/attempted < 1`. Guards F1. Always on (cheap).
- **O3 (INVARIANT) — per-state marginal matches closed-form Boltzmann (INV-M1).** For `REMC`, for every
  state `k`: `chiSquareStatistic(U_hist_k, expectedFromWeights(gammaSubsampledWeights(·, n_c/2, RT_k), N))
  < chiSquareCritical(nbins−1, 1e-4)`. **Discrimination control (SHALL, not optional):** the same
  histogram SHALL **reject** the neighbour temperature, `chiSquareStatistic(U_hist_k, … RT_{k±1} …) >
  chiSquareCritical(nbins−1, 1e-4)` — mirroring the wrong-DOF control of
  `tests/TestEnsembleValidation.cpp:186-196`. Without the reject half the pass half has no power. Slow
  tier.
- **O4 (INVARIANT) — swap does not bias the marginal (INV-M2).** Run `Default` (no swaps) and `REMC`
  (swaps active) on identical SR and ladder; both SHALL satisfy O3. The `Default` run is the control that
  isolates inner-sampler correctness from the exchange logic: if `Default` fails O3 the inner GC-HMC (not
  the driver) is at fault. Slow tier.
- **O5 (LEMMA) — equipartition mean, exact value.** For each state `k`, `⟨U⟩_k = (n_c/2)·RT_k` within
  `4·stderr` (`MeanAccumulator`, `StatTest.hpp:190-215`; the `k = 4` two-sided false-fail convention of
  `StatTest.hpp:11`). Also `⟨2βU⟩_k = n_c` (chi-square mean). Cheap; MAY run in the smoke tier at reduced
  `N` with a looser band. Directly ties the marginal to equipartition.
- **O6 (LEMMA) — swap rate vs closed-form overlap.** For each neighbour pair `(a, b)`, the measured
  accepted/attempted ratio SHALL match the independent-Gamma Monte-Carlo estimate of
  `E[min(1, exp(−(β_b−β_a)(U_a−U_b)))]`, `U_a ~ Gamma(n_c/2, RT_a)`, `U_b ~ Gamma(n_c/2, RT_b)`, within the
  binomial standard error of the run's attempt count (`kofke_2002` eq:10 overlap integral estimated by
  sampling). Slow tier. Distinguishes a driver that swaps at the wrong rate (e.g. a subtly wrong β
  pairing) from a correct one, using a closed-form-marginal reference rather than self-consistency.
- **O7 (LEMMA) — ladder round trip, and the parity-frozen discriminator (INV-M5).** Track each replica's
  current state index each round; count a round trip when a replica reaches state `T−1` after having been
  at state `0` and returns to `0`. For `REMC`: total round trips `≥ 1` (SHOULD be several) and every state
  occupied by `≥ 2` distinct replicas. **Discriminator (SHALL):** a mutation that freezes the pairing
  parity (hold `exchangeRound` fixed so only `(0,1),(2,3)` are ever attempted, never `(1,2)`) SHALL
  partition `{0,1}` from `{2,3}` and produce **zero** round trips — proving O7 detects the B7 defect that
  O3 cannot (a parity-frozen chain keeps each intra-block state's marginal locally correct). `Default`
  SHALL produce zero round trips. Slow tier for the count magnitude; the parity mutation MAY run at
  reduced `N`.
- **O8 (LEMMA) — the marginal oracle has power (mutation).** At the same `N` used by O3, two deliberately
  corrupted drivers SHALL cause O3 to **fail** on at least one state (proving O3 is not vacuous), mirroring
  `TestRexAcceptanceAlgebra`'s paired correct/broken structure:
  1. **Swapped-β acceptance** — use `logα = −((β_H − β_C)(U(x_X) − U(x_Y)))` with `β_C`/`β_H` assigned to
     the wrong states (the `swapBetas` corruption of `tests/TestRexAcceptanceAlgebra.cpp:89-97`,
     lifted to the chain). This breaks detailed balance of the swap kernel, changes the stationary
     distribution, and biases at least one state's `U`-histogram.
  2. **By-replica recording (INV-M6 violation)** — histogram state `k`'s samples from the fixed replica
     object `k` (ignoring `thermo2ReplicaIxs_`). Because a replica wanders across temperatures under
     label swaps, its `U` is a temperature mixture and fails the single-`Gamma(n_c/2, RT_k)` fit.
  Both mutations SHALL be exercised as expected-fail assertions (the test asserts O3 would reject),
  keeping the corrupted paths out of the production-mirror code path.

NOTE on false-fail budget: O3/O4/O6/O7/O8 run at `α = 1e-4`; across the ~`4·T` distribution assertions the
expected suite false-fail rate stays well under 1 %, consistent with `StatTest.hpp:11-17`. Fixed seeds
(`TestHelpers.hpp` `Rng`, Rule 9) make each run reproducible.

---

## Consequences and trade-offs

- **The oracle validates the algorithm, not the production driver object.** Because `Context::RunREX`
  needs OpenMM and is uncompiled, the oracle mirrors its control structure rather than instantiating it.
  This is the suite's convention (algebra test, Fixman test) and is honest about what it proves: the
  label-swap REMC algorithm composed with GC-HMC preserves product-Boltzmann, and the mirror matches
  `src/Context.cpp` line-for-line. A future OpenMM-backed integration test (Open questions) closes the
  residual gap that the mirror could drift from the production driver. The line-referenced mirror plus the
  SHOULD to share the pairing/label-swap header bound that drift risk.
- **Pure-translational system forecloses richer geometry.** Keeping `U` exactly quadratic (translational
  DOFs, constant metric) is what buys the closed form. A metric-varying or torsional system would need
  Fixman-on and would lose the analytic `Gamma` target, reducing the oracle to self-consistency — which is
  exactly what this spec avoids. The driven-path oracles (deferred) are where torsional geometry belongs.
- **Slow tier cost.** `T = 4` replicas × `nInner` tiny HMC moves × `~10⁶` rounds is seconds on the
  analytic bridge (no OpenMM), well within the existing slow-tier budget of `TestFixmanBoltzmann` /
  `TestEnsembleValidation`.

---

## Open questions

- **OQ-1 (non-blocking) — production-driver integration test.** This oracle proves the REMC algorithm is
  correct and that the mirror matches `src/Context.cpp` today, but it does not execute `Context::RunREX`
  itself. Closing that gap requires an OpenMM-backed integration test that runs `Context::RunREX(REMC, …)`
  on a small real system and checks (a) the coordinate-swap `runREX` INVARIANT-EQUIV equivalence (parent
  spec, `include/Context.hpp:84-90`) and (b) the same per-state marginal against a trusted long single-
  temperature run. What is needed to proceed: `Context::RunREX` compiled and runnable against the OpenMM
  singleton (currently "not compiled or run", `include/Context.hpp:104-109`). Until then the mirror is the
  best available driver-level oracle; the parent spec's V2 remains the placeholder for the OpenMM version.
  This does not block landing the analytic oracle.
- **OQ-2 (non-blocking) — reference-parameter tuning.** The reference ladder (ratio 1.3, `n_c = 6`) is
  chosen for healthy swap acceptance and clear neighbour-`Gamma` separation, but the exact `nInner` /
  timestep that yields `50`–`80 %` inner acceptance and adequate O7 round-trip counts at the chosen `N` is
  an empirical calibration the coder SHALL fix by making the `Default` control (O4) and the round-trip
  count (O7) pass first; if round trips are too rare at ratio 1.3, narrow the ladder (more states or
  smaller ratio) rather than weaken the assertion. No theory question is open here — only calibration.

### References

- `sugita_okamoto_1999` — Sugita & Okamoto, *Chem. Phys. Lett.* (1999), DOI 10.1016/S0009-2614(99)01123-9.
  eq:7 product-of-Boltzmann generalized-ensemble weight; eq:15 swap cost `Δ=(β_n−β_m)(E_i−E_j)`; eq:17
  `min{1,exp(−Δ)}`.
- `chodera_2011_gibbs_replica_exchange` — Chodera & Shirts, *J. Chem. Phys.* (2011), DOI 10.1063/1.3660669.
  eq:10 permutation-conditional `π(S|X)`; eq:24 pair-swap Metropolis acceptance.
- `kofke_2002` — Kofke, *J. Chem. Phys.* (2002), DOI 10.1063/1.1507776. eq:1 swap acceptance; eq:9
  gamma-distributed energy for constant heat capacity; eq:10 exact acceptance overlap integral; eq:3/eq:11
  temperature-ratio scaling (geometric-ladder justification).
- Swendsen & Wang (1986) — origin of the replica-exchange / parallel-tempering swap move (indexed in
  `references/index.yaml` as a dependency of `sugita_okamoto_1999` / `kofke_2002`).
- `docs/specs/replica-exchange-nonequilibrium-work.md` — INV-1/INV-2, B6 (swap detail), B7 (pairing/parity
  fix), the label-swap object model, and V2 (the invariant this oracle operationalizes).
