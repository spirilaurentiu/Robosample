# DOC-OpenMMSystemBuilder: OpenMM System/Context/Integrator construction

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the body of `initialize` - particle/box/force/integrator/platform
  construction that brings up the one OpenMM `System`/`Context`/`Integrator`
  (ARCHITECTURE section 3, section 7; MODULES.md O5). The ~235-line construction sequence.
- **Layer:** Bridge (ARCHITECTURE section 2).
- **Ownership:** constructs the OpenMM objects the `OpenMMContext` singleton then
  owns for the whole run (ARCHITECTURE section 4). The builder hands ownership to the
  singleton; recover the hand-off point.
- **Invariants (HYPOTHESES):** construction order matters - particles before forces,
  forces before integrator/context; box vectors come from `PeriodicBox`
  (`reducedBoxVectors`, DOC-PeriodicBox, single source after O2). Document the
  ordering contract and each input's role; verify the sequence against the call
  order in the source.

## 2. Scope

- **Files:** `bridge/OpenMMSystemBuilder.{hpp,cpp}` (extracted from
  `OpenMMContext.cpp::initialize`).
- **Public symbols:** the build entry point(s) and its inputs (topology, box, force
  config, platform selection).
- **Known gaps to close:** the mixed responsibilities were the reason for the split;
  document each construction phase's contract separately, not as one narrative.

## 3. Evidence pointers (tests exercising the module)

- TestAlchemy - the only real OpenMM `Context` (TESTS.md section 1). The Level-1 run and
  the 4-replica REMC proof run (VERIFY B3) exercise bring-up end to end.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; construction-order contract
  stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-OpenMMSystemBuilder-findings.md` present; no unmatched
  `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
