# DOC-GibbsSweep: the single Gibbs-sweep primitive

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the one "run a sweep over the schedule" primitive that deduplicates the
  sweep loop written three times in `Context.cpp` (plain REX, driven RENE round, REMC
  subround) (ARCHITECTURE section 3; MODULES.md C3).
- **Layer:** Workflow (ARCHITECTURE section 2).
- **Ownership:** borrows the World schedule and the coordinate currency; owns nothing.
- **Invariants (HYPOTHESES):**
  - The sweep is exactly: `for each World in schedule: push coords -> generateSample ->
    pull coords` (ARCHITECTURE section 3). Per-atom Ground coordinates (nm) are the sole
    inter-world currency (INV-3, DOC-GeometryFitter). Document the loop's contract -
    ordering, the push/pull currency, and that it leaves no Worlds holding cross-sweep
    state.
  - It SHALL be behaviorally identical to the three inlined copies it replaces
    (pure-motion dedup); the three call sites are the evidence that the extracted
    primitive's contract matches each.

## 2. Scope

- **Files:** `workflow/GibbsSweep.{hpp,cpp}` (from the three duplicated loops in
  `Context.cpp` 457-490 / 977-1021 / 1112-1143).
- **Public symbols:** the sweep entry point.
- **Known gaps to close:** confirm the three former copies were identical up to the
  schedule they iterate; any behavioral difference between them is a finding (the
  dedup would not be pure motion).

## 3. Evidence pointers (tests exercising the module)

- TestGibbsWorlds - FAST contract (TESTS.md section 2). The Level-1 and REMC runs (VERIFY B3)
  exercise the sweep end to end.

## 4. Exit criteria

- Doxygen warning-free; the sweep primitive documented (order, currency, no cross-sweep
  state); comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-GibbsSweep-findings.md` present; any inter-copy difference recorded;
  no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
