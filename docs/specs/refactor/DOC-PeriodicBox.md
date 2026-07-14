# DOC-PeriodicBox: periodic-box value type

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the periodic-box value type and its reduced-box-vector math - the
  single source for box vectors after `SPLIT-O2` folds the duplicate
  `computePeriodicBoxVectors_Context` out of `OpenMMContext.cpp` (ARCHITECTURE section 7,
  MODULES.md O2).
- **Layer:** Domain data / model (ARCHITECTURE section 2).
- **Ownership:** value type, copied by value; no heap.
- **Invariants (HYPOTHESES):** `reducedBoxVectors` is the single implementation of
  the reduced-form computation (the whole point of O2 - one source, not two). Whole-
  molecule periodic imaging in `OutputWriter` reads these vectors; document the
  reduced-form contract (triclinic reduction convention) callers depend on. Verify
  the two former copies agreed before the fold - a disagreement is a finding, not
  something to paper over.

## 2. Scope

- **Files:** `model/PeriodicBox.hpp` (box value type + `reducedBoxVectors`).
- **Public symbols:** the box type and its vector/reduction accessors.
- **Known gaps to close:** state the reduced-vector convention precisely; it is
  consumed by OpenMM box setup and by DCD box records via OutputWriter - confirm both
  read the same convention.

## 3. Evidence pointers (tests exercising the module)

- TestPeriodicBoundary - FAST contract (TESTS.md section 2). Pins the imaging/reduction
  behavior.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; contracts behavioral;
  comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-PeriodicBox-findings.md` present; the O2 pre-fold-agreement check
  recorded; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
