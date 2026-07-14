# DOC-ReactionReporter: per-body reaction-force snapshot and CSV

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `captureReactionSnapshot` + the five `reaction*_` members - an opt-in,
  orthogonal per-body spatial-force snapshot emitted as a reaction CSV at DCD cadence
  (ARCHITECTURE section 7; MODULES.md W3; MEMORY: reaction-force monitoring campaign - per
  frame, dynamical). High independence.
- **Layer:** Domain / world (ARCHITECTURE section 2).
- **Ownership:** owns its snapshot buffers; reads reaction forces from the solver
  (RobotEngine reaction methods); the CSV rows are consumed by OutputWriter.
- **Invariants (HYPOTHESES):**
  - The reported per-body spatial force follows **INV-1** convention (torque about
    body origin, net force, Ground frame) - the same convention as `bodyForceG`
    (ARCHITECTURE section 5). Document the reported quantity's definition and frame.
  - The snapshot is dynamical (captured per frame during the move), not a
    post-hoc recomputation. State the capture point (when in the move it is taken).

## 2. Scope

- **Files:** `world/ReactionReporter.{hpp,cpp}` (from `World.cpp`).
- **Public symbols:** `captureReactionSnapshot`, the accessors OutputWriter reads.
- **Known gaps to close:** state the reported force's convention (INV-1) and capture
  cadence/timing explicitly; it is opt-in, so the enable condition is part of the
  contract.

## 3. Evidence pointers (tests exercising the module)

- TestReactionForces - FAST contract (INV-1; TESTS.md section 2). Same INV-1 as the
  ForceReducer calibration exemplar - reuse that phrasing.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; INV-1 convention and cadence
  stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-ReactionReporter-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
