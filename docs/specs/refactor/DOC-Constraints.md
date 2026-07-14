# DOC-Constraints: loop-closure SHAKE/RATTLE + Fixman log-det

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** loop-closure distance constraints - SHAKE position projection, RATTLE
  velocity projection, and the loop-closure Fixman log-det - sharing one
  `assembleConstraintRow` (`G`-assembly) after the dedup (ARCHITECTURE section 7; MODULES.md
  dedup, section 2).
- **Layer:** Algorithms / dynamics (ARCHITECTURE section 2).
- **Ownership:** `ConstraintSet` owned by `World` by value (ARCHITECTURE section 4); operates
  on the borrowed `RobotState`. Note the managed forward-decl cycle with RobotEngine
  (ARCHITECTURE section 6.2) - document the interface, not the include.
- **Invariants (HYPOTHESES):**
  - **INV-6 constraint consistency** (ARCHITECTURE section 5): the `G` constraint Jacobian is
    assembled identically across SHAKE, RATTLE, and the Fixman log-det, so the
    correction matches the projection. `calcConstraintLogDet` returns **0 for acyclic
    molecules**. State both as contract; the dedup (`assembleConstraintRow`) is what
    makes the three uses provably identical - document the shared row's contract once.
- **Test-access coupling:** `ConstraintTestAccess.hpp` friends into private
  `solveSmallSpd`/`solveCoupling` (TESTS.md section 3); those retarget when Constraints
  splits. Treat them as internal - document per `documenter.md` section 4.1 (private helpers
  get local behavior, not external contract).

## 2. Scope

- **Files:** `dynamics/Constraints.{hpp,cpp}`.
- **Public symbols:** the SHAKE/RATTLE entry points, `calcConstraintLogDet`, the
  `ConstraintSet` surface. Private: `assembleConstraintRow`, `solveSmallSpd`,
  `solveCoupling`.
- **Known gaps to close:** state the acyclic->0 log-det contract and the shared-`G`
  guarantee explicitly (INV-6); they are the correctness link the Fixman correction
  relies on.

## 3. Evidence pointers (tests exercising the module)

- TestConstraints, TestConstraintSolver - FAST contract (INV-6; TESTS.md section 2/section 3).
  `ConstraintTestAccess.hpp` is the private-solver harness (TESTS.md section 3).

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; INV-6 (shared `G`, acyclic->0)
  stated as contract; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-Constraints-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
