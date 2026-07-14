# DOC-OutputWriter: run output emission (CSV / DCD / imaging)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `writeOutputs*`, `writeReactionRows`, the DCD scratch buffers, and
  whole-molecule periodic imaging - emit the moves CSV, per-replica energy CSV,
  per-replica DCD, and reaction-force CSV (ARCHITECTURE section 3, section 7; MODULES.md C2).
- **Layer:** Workflow (ARCHITECTURE section 2).
- **Ownership:** owns the DCD scratch buffers and the `DCDWriter` instances (DOC-
  DCDWriter); borrows coordinates/energies from the replicas per frame.
- **Invariants (HYPOTHESES):**
  - Periodic imaging is applied **here** (whole-molecule imaging using
    `PeriodicBox::reducedBoxVectors`, DOC-PeriodicBox), not in `DCDWriter` - document
    that boundary. The imaged coordinates are output-only and never fed back into
    sampling (a purely presentational transform).
  - The reaction rows come from `ReactionReporter` snapshots (DOC-ReactionReporter,
    INV-1 convention). CSV/DCD frame cadence is the run's write frequency; state it.

## 2. Scope

- **Files:** `workflow/OutputWriter.{hpp,cpp}` (from `Context.cpp`).
- **Public symbols:** `writeOutputs*`, `writeReactionRows`, the file-truncation/init
  surface.
- **Known gaps to close:** state that imaging is output-only (never re-enters
  sampling) and the per-file cadence; these are the contracts a reader needs, not the
  buffer bookkeeping.

## 3. Evidence pointers (tests exercising the module)

- No dedicated C++ test; the Level-1 run and the 4-replica REMC proof (VERIFY B3) diff
  the emitted CSV/DCD - those baselines are the behavioral oracle. Record the
  unit-test gap in findings.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; imaging-is-output-only and
  cadence stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-OutputWriter-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
