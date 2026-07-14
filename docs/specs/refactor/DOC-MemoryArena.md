# DOC-MemoryArena: single-slab bump allocator

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the single-slab allocator backing `RobotState`'s per-step cache - one
  arena hands out the ~35 raw cache pointers (ARCHITECTURE section 4, MODULES.md section 5:
  single-concept, no split).
- **Layer:** Infrastructure / util (ARCHITECTURE section 2). No engine dependency.
- **Ownership:** the arena **owns** its slab; pointers it hands out are **borrowed**
  and valid only while the arena lives and is not reset. This lifetime rule is the
  central contract - recover it from `RobotState`'s allocation site and every
  accessor that dereferences an arena pointer.
- **Invariants (HYPOTHESES):** handed-out pointers stay valid until the arena is
  reset/destroyed; a reset invalidates all outstanding pointers. Document alignment
  guarantees only if a caller relies on them (verify).

## 2. Scope

- **Files:** `util/MemoryArena.hpp`.
- **Public symbols:** the arena type, its allocate/reset/reserve surface.
- **Known gaps to close:** no dedicated unit test. Recover the lifetime and
  invalidation contract from `RobotState` construction and use; record the coverage
  gap in findings.

## 3. Evidence pointers (tests exercising the module)

- None dedicated. Exercised indirectly through every `RobotState`-building test
  (TestKineticEnergy, TestMassMatrix, the oracles). Note this in findings.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; ownership/lifetime stated as
  contract; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-MemoryArena-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
