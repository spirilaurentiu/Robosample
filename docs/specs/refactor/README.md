# Robosample Refactor - Spec Ticket Set

Status: draft, 2026-07-12. This directory holds the executable specs a coding
agent runs to refactor the `src/`+`include/` engine. It is downstream of the
recovery docs in [`../../architecture/`](../../architecture/): read
[`ARCHITECTURE.md`](../../architecture/ARCHITECTURE.md),
[`MODULES.md`](../../architecture/MODULES.md), and
[`TESTS.md`](../../architecture/TESTS.md) first - they define the layer map, the
split proposals (ticket IDs W1..., C1..., O1..., R1...), the invariants, and the test
triage that every ticket below references.

**Charter (from `.claude/agents/architect.md`).** Splits are **pure code motion**:
no logic change, no renames beyond file placement, no "improvements." Moving and
editing never share a commit. Every ticket is self-sufficient - a coding agent
receives one ticket and nothing else. Every ticket ends in machine-checkable exit
criteria. Public API is preserved unless a separate, human-approved API spec says
otherwise (target: zero exported-symbol delta on public symbols).

## Two-phase freeze

- **Phase A - tests frozen, source moves.** During all `SPLIT-###` execution,
  test files change in exactly one way: include-path fixes listed explicitly in
  the source ticket that caused them. After every ticket, the pass set, test
  count, and assertion counts match the Stage-0 baseline (`VERIFY.md`).
- **Phase B - source frozen, tests move.** After the source tree is stable and
  verified, `TEST-###` tickets run. The oracle inverts: frozen source + baselines
  verify the test refactor.

Source splits and test motion never share a phase. `DOC-###` execution happens
after the splits so Doxygen lands on final files (the DOC tickets are *written*
now, *executed* in the doc phase).

## Ticket inventory and execution order

Execute leaf-first (lower layers before higher) to avoid merge conflicts. IDs
match `MODULES.md`.

### Phase 0 - baselines (blocking prerequisite)
- `VERIFY.md` section 1 - capture the Stage-0 baselines (pass set x3, assertion counts,
  exported public symbols, representative-run outputs) **before any ticket runs.**

### Phase A - source SPLIT tickets, in order
1. **Leaf math/util:** `SPLIT-R1` (hinge_linalg), `SPLIT-R2` (robo_debug),
   `SPLIT-O1` (MTSIntegrator), `SPLIT-O2` (PeriodicBox fold).
2. **Bridge:** `SPLIT-DEDUP-FORCEREDUCER`, `SPLIT-O3` (ForceFactory), `SPLIT-O4`
   (AlchemyForceFactory), `SPLIT-O6` (GpuKinematics), `SPLIT-O5` (SystemBuilder).
3. **Dynamics:** `SPLIT-I1` (verletStep split), `SPLIT-R3` (RobotEngine TU split),
   `SPLIT-R4` (joint taxonomy -> JointKernels; gated on OQ-2),
   `SPLIT-DEDUP-CONSTRAINTROW`, `SPLIT-DEDUP-BATBENDSTRETCH`.
4. **World services (high-independence first):** `SPLIT-W1` (ModelBuilder),
   `SPLIT-W3` (ReactionReporter), `SPLIT-W5` (FixmanCorrection), `SPLIT-W2`
   (DockingMove), `SPLIT-W6` (GeometryFitter), `SPLIT-W7` (CartesianSolvent),
   `SPLIT-W8` (VelocityDistortion), `SPLIT-W4` (NcmcMove), `SPLIT-W9` (HmcMove - last).
5. **Workflow:** `SPLIT-C1` (StartupValidator), `SPLIT-C2` (OutputWriter),
   `SPLIT-C3` (GibbsSweep), `SPLIT-C4` (ReplicaExchangeDriver + SwapAcceptance),
   `SPLIT-C5` (DrivenRexDriver), then human-gated `SPLIT-C6` (LegacyCoordSwapRex -
   isolates the `runREX` oracle into its own TU; OQ-3 DECIDED keep, so C6 applies
   unless a human retires `runREX`). `runREX` retirement is the alternative to C6,
   not bundled.

The list above is a valid **total order**. It is not the only one: the merge
hazard is same-file editing, so tickets on disjoint source files may run
concurrently. The lane structure below is the authoritative DAG; `MODULES.md section 4`
carries the same graph. A single agent runs one lane serially; different lanes
run in parallel.

### Lane structure (parallel execution)

Four serial lanes, keyed by the god-file each edits, plus independent leaves and
two barriers. Within a lane, tickets serialize (each shifts the file's line
numbers, so downstream tickets re-anchor by symbol - the ranges are Stage-0
advisory, the symbol name is authoritative).

- **Lane E - RobotEngine:** `R1` -> `R2` -> `R3` -> `I1` -> `R4`. (`I1` edits
  `RobotIntegrator.hpp`; `R4` depends on `I1`, see CX-9. `R4` also gated on OQ-2.)
- **Lane B - Bridge/OpenMM:** `DEDUP-FORCEREDUCER` -> `O1` -> `O6` -> `O3` -> `O4` ->
  `O5`. (`O6` after `DEDUP-FORCEREDUCER`; `O5` last.)
- **Lane W - World:** `W1` -> `W3` -> `W5` -> `W2` -> `W6` -> `W7` -> `W8` -> `W4` ->
  `W9`. Longest chain (critical path). `W5` before `W2`/`W4` (CX-6 shared header).
- **Lane C - Context:** `C1` -> `C2` -> `C3` -> `C4` -> `C5` -> `C6` (human-gated;
  C6 depends on C4 having relocated `runREX`).
- **Leaves (own agent, any time):** `DEDUP-CONSTRAINTROW` (`Constraints`),
  `DEDUP-BATBENDSTRETCH` (`BatScaling`).

**Barriers - the only tickets that edit two god-files at once; they block the
lanes they straddle.** Run each either before both its lanes begin or after
both finish, never concurrently with them:

- `SPLIT-O2` (PeriodicBox fold) edits `OpenMMContext.cpp` **and** `Context.cpp`
  - barrier across Lane B and Lane C.
- `SPLIT-DEDUP-CTXFORWARDERS` (section 6.3/section 6.4 forwarders, CX-3) edits `Context.hpp`,
  `OpenMMContext.hpp`, **and** `Context.cpp` - barrier across Lane B and Lane C;
  it adds two public `Context` symbols, so it is human-approved (VERIFY section 3) and
  is NOT pure motion. Recommended slot: after Lane C completes.

**Cross-lane coordination the orchestrator (not any single agent) SHALL hold:**
1. **Versioned baseline.** `DEDUP-FORCEREDUCER` adds `TestForceReducer` to B1
   mid-stream; lanes asserting "pass set = B1" rebase onto B1' once it lands.
2. **Shared test-tree edits.** Disjoint source files can still share a test file
   (one that includes both `World.hpp` and `Context.hpp`); Phase-A include-fixes
   from two lanes may collide there. Serialize test-file edits or apply one
   post-merge include pass.
3. **Per-lane worktree/branch**, merged at ticket or lane completion.

### Doc phase - after all SPLITs verified
- `DOC-###` tickets, one per target module, calibration files first (one `.cpp`,
  one `.cu`, one test file -> human review -> exemplars attach to the rest), then
  dependency order (leaf helpers before public API).

### Phase B - test tickets, after source frozen
The `TEST-###` tickets cover the `TESTS.md section 5` seams. Full inventory:
- `TEST-001-slow-tier.md` - slow-tier optimization, gated by the mutation-power test.
- `TEST-002-shared-fixtures.md` - `RobotFixture` base + `TestPhysConstants.hpp`.
- `TEST-003-harmonic-bridge.md` - `HarmonicBridge<Policy>` consolidating the ad-hoc bridges.
- `TEST-004-sampling-harness.md` - `runHmcChain` / `runNcmcLoop` helpers.
- `TEST-005-split-oracle-monsters.md` - split the two oracle files by case family.
- `TEST-006-split-ncmcwork.md` - split `TestNCMCWork` into FAST algebra + SLOW work-chain.
- `TEST-007-mirror-module-structure.md` - mirror the module map; contract comments.
- `TEST-008-openmm-suite.md` - OpenMM's own suite (infra, outside the freeze).

### Behavioral specs (separate track, human-approved, NOT pure code motion)
- [`../rex-default-driver.md`](../rex-default-driver.md) - make label-swap `RunREX`
  REMC the default driver; `R = 1` as the degenerate case. Behavioral change; runs
  on its own, not inside a SPLIT ticket.
- [`../rex-stationary-distribution-oracle.md`](../rex-stationary-distribution-oracle.md)
  - the REX driver-level behavioral oracle (researcher artifact, pending review).

## SPLIT ticket template

Every SPLIT ticket SHALL use this shape (from the charter):

```markdown
# SPLIT-XX: <verb> <concept> from <file>

Status: draft (Phase A, source split - pure code motion). Executor: `coder`.
Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
<file> (<LOC>) mixes <N> responsibilities: <list>. This ticket extracts
<concept>. Invariants below come from ARCHITECTURE.md and are restated in full.

## Moves (exactly what)
- <symbol> (lines ~A-B) -> <target>.{hpp,cpp} [public | anonymous ns | detail::]
- ... (one line per symbol, with source line range and destination + visibility)

## Public API after this ticket
<target>.h exports: <exact list>. Nothing else. Header self-contained (IWYU),
forward-declare <X>, #pragma once.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- <origin file> gains #include "<target>" and loses the moved code; nothing else
  in it changes.
- Invariants that SHALL remain true: <the applicable INV-# from ARCHITECTURE.md>.
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- <test file>: <why it breaks under pure motion> -> <the one allowed include fix>.
  (If a test reaches a moved private helper, name the detail:: header it retargets to.)

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0.
- nm diff on public symbols: empty.
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- <target>.h <= 300 LOC; <target>.cpp <= 600 LOC (or a claimed, justified exception).
```

DOC and TEST tickets follow the shapes in `.claude/agents/architect.md` (Context
= module purpose/layer/ownership/invariants as *hypotheses to verify*; scope;
evidence pointers; calibration exemplar; exit criteria).

Executor and style per ticket class (current agent roster, `.claude/agents/`):

| Ticket | Executor | Style | Exit criteria |
| --- | --- | --- | --- |
| `SPLIT-###` | `coder` | `styles/spec.md` | VERIFY section 2 |
| `DOC-###` | `documenter` | `styles/reference.md` | VERIFY section 4 |
| `TEST-###` | `coder` | `styles/spec.md` | VERIFY section 5 |

The `coder` compiles and self-reviews each SPLIT/TEST ticket, then delegates the
Doxygen contract for any new or changed symbol to the `documenter`. The `reviewer`
performs the independent adversarial review before merge (Review exit: no confirmed
Blocking findings). Each ticket names its executor and style in its Status line.

---

## Cross-cutting findings - human decisions before execution

The six ticket authors surfaced issues that span tickets and need a human call
*before* any coder runs. Recorded here so they are not buried in individual
tickets.

**ACCEPTED (2026-07-13).** The human accepted every cross-cutting decision below.
Resolutions: CX-1 TU-split-first; CX-2 `GLOB_RECURSE` + subtree layout; CX-3 both
public-symbol changes (O2 deletion, CTXFORWARDERS additions) approved; CX-4 stretch-only
BAT dedup; CX-5 oracle stays independent; CX-6 the `engine_helpers.hpp` home with its
one predicted Phase-A test edit. CX-7/CX-8 are corrections, no decision. The
per-CX resolution is appended to each item.

### CX-1 - "Pure code motion" means TU-split, not class extraction (highest impact)
For `World` (W1-W9) and `Context` (C1, C2, C4, C5), the extracted methods depend
on the origin class's private state. A true standalone-class extraction would move
that state across a boundary - violating zero public-symbol nm-delta (these are
PyBind11-bound) and the token-identical-moved-block rule (`VERIFY section 2.7`). The
tickets therefore do a **multi-TU class split**: declarations stay on `World` /
`Context`; method *bodies* move into a new subtree TU as `World::` / `Context::`
definitions. `MODULES.md` frames C1-C5 / W1-W9 as *new classes* - that is the
*eventual* target, reached by a **separate, human-approved restructuring** after
the TU split, never inside a pure-motion ticket. **ACCEPTED (2026-07-13): TU-split-first.**
The standalone-class restructuring stays a separate, later, human-approved track;
`MODULES.md`'s new-class framing is the eventual target, not the Phase-A action.
Clean exceptions (genuine new modules/free functions, not TU splits): `C3`
(`robo::gibbs::stepWorldInGround`), `R1` (hinge_linalg), `R2` (robo_debug),
`O1`/`O3`/`O4`/`O6`, and the leaf math.

### CX-2 - CMake build change is a prerequisite (human-approved)
The `src`/`include` globs (`CMakeLists.txt`: `file(GLOB ROBOSAMPLE_CXX_SOURCE_FILES
.../src/*.cpp)` and `.../include/*.h*`, currently near :150/:152) match
**non-recursively**, so files placed in `src/world/**`, `src/bridge/**`, etc. will
not compile. The module subtree requires switching to `GLOB_RECURSE` (introduced by
the first split that creates a subdir - `W1`) or keeping files flat under `src/`.
The test glob (`file(GLOB _test_sources ... tests/Test*.cpp)`, near :602) already
picks up new `Test*.cpp`, but the oracle-fixture wiring (the `TestRoboticsOracle*`
`foreach`, near :621) hardcodes `TestRoboticsOracle`/`TestRoboticsOracleMolecule`
target names and SHALL gain the `TEST-005` split-binary names. Cite these globs by
symbol, not line: the line numbers drift. **ACCEPTED (2026-07-13): `GLOB_RECURSE`
+ subtree layout.** `SPLIT-W1` (the first ticket that creates a subdir) carries the
`GLOB_RECURSE` switch as its build-system prerequisite; the oracle-fixture `foreach`
gains the `TEST-005` split-binary names in Phase B.

### CX-3 - Public-symbol changes needing API approval (VERIFY section 3)
**ACCEPTED (2026-07-13): both changes approved.** These are the two sanctioned
non-empty public-symbol nm-deltas; every other ticket holds to zero.
- `SPLIT-O2` deletes `computePeriodicBoxVectors_Context` (zero in-scope callers).
  Deletion approved; no forwarder retained.
- `SPLIT-DEDUP-CTXFORWARDERS` adds two public `Context` symbols
  (`setSeparateForceGroups`/`setEnforcePeriodicBox`) so PyBind stops reaching the
  OpenMM singleton (ARCHITECTURE section 6.3). Addition approved.

### CX-4 - Numerical hazard: BAT dedup is not bitwise-safe
`SPLIT-DEDUP-BATBENDSTRETCH`: routing `theta0` through `readBatCoord` uses a
normalized-vector dot vs the current inline `dot(v1,v2)/(r*r2)` - **not
bitwise-identical**, so it would break the `VERIFY B3` bitwise-HMC-log baseline.
The ticket scopes the shared read to the stretch-`r` term only (bitwise-safe) and
leaves the bend inline; full sharing is a **separate human-approved numerical-change
ticket**, never folded into code motion.

### CX-5 - Do NOT dedup the linear-algebra oracle
`SPLIT-R1` extracts the engine's dense solver to `hinge_linalg`. `tests/RobotLinearAlgebra.hpp`
(`namespace robo_linalg`) is a **deliberately independent** re-implementation, and
`TestLinearAlgebraOracle` diffs the engine against it. Pointing the test at
`hinge_linalg` would make that comparison **tautological**. The oracle stays
independent. (Corrects the earlier `TESTS.md` note.)

### CX-6 - Shared file-local helpers need single-source homes (creates sequencing)
Anonymous-namespace helpers used by more than one extraction can't be shared across
TUs, so the first ticket to need them hoists them to a header:
- Quaternion converters `rotationToQuaternion`/`quatToRotation`/`safeLogSineSqr`
  (W2/W4/W5) -> a new `engine_helpers.hpp`. This creates a **W5 -> (W2, W4)** ordering
  dependency and the **only** predicted Phase-A test edit in the whole set: an
  `EngineHelpers.hpp` -> `engine_helpers.hpp` include retarget to avoid an ODR clash.
  Alternative `world/detail/` home would touch zero tests. **ACCEPTED (2026-07-13):
  `engine_helpers.hpp` home**, with its one predicted Phase-A test edit approved (the
  sole test edit in the whole SPLIT set; `VERIFY` criterion 9 expects exactly it).
- `spatialDot` (R3) -> engine-private `src/RobotEngine_internal.hpp`.
- `nmaDebugEnabled` (W8/W9), `nma_debug.hpp`.

### CX-7 - Corrected coupling assumption
`TestRexAcceptanceAlgebra` does **not** call `Context::attemptREXSwap` - it
reimplements the acceptance formulas in-file and never instantiates a `Context`.
So `C4` has no live test breakage; the retarget-to-`SwapAcceptance` guidance is a
forward note for a hypothetical future reproducer, not a predicted break. (Corrects
`TESTS.md section 3`.)

### CX-8 - Dead code found (move, do not delete)
`dockDebugEnabled` (`World.cpp:45`) has no call site. A code-motion ticket SHALL NOT
change symbol count, so `W2` moves it unchanged; removal is a reviewer decision.

### CX-9 - Ordering interactions to honor
`R4` and `I1` both touch the quaternion-drift block (`RobotIntegrator.hpp`
241-258). **Decided: I1 -> R4.** I1 extracts the drift loop (with the inline
`wHalf`/`advanceQuatExp` block) into `driftPositions`; R4 then lifts that block
out of `driftPositions` into `JointKernels::jointDriftQuat`. No line overlap -
R4 re-anchors to the post-I1 symbol, not the `verletStep` lines. `O5`
(SystemBuilder) lands after `O1/O3/O4/O6`. `W9` (HmcMove) last among World. `C5`
may surface a pre-existing latent compile defect in the uncompiled driven-REX TU
(OQ-5) - report, do not fix-forward.

### CX-10 - `C2`/`C4` both claim `Context.cpp:21` (`kBoltzmann_kJ`)
The top-of-file anonymous namespace holds two symbols: `boxFromReducedVectors`
(output-only) and `kBoltzmann_kJ` (REX-only). **Decided: C2 owns
`boxFromReducedVectors`, C4 owns `kBoltzmann_kJ`.** C2 (runs first) lifts only
the function and leaves the constant + namespace braces in place; C4 relocates
the constant to `RexInternal.hpp` and removes the now-empty namespace. C2 SHALL
NOT sweep the whole `namespace {...}` - burying `kBoltzmann_kJ` in
`OutputWriter.cpp` would break the REX TUs' build.
