# DOC-ForceFactory: OpenMM force builders

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the `create*Force` builders and `addStandardExclusions` that assemble
  the OpenMM force objects during system construction (ARCHITECTURE section 7, MODULES.md
  O3).
- **Layer:** Bridge (ARCHITECTURE section 2).
- **Ownership:** the builders **construct** OpenMM force objects and transfer them to
  the OpenMM `System` (OpenMM takes ownership on `addForce`). Recover the transfer
  point; document who owns each force after construction.
- **Invariants (HYPOTHESES):** force-group assignment and the standard exclusion set
  are what the energy/force validation (the 1e-6 force-group band, VERIFY section 3)
  depends on. Document each builder's force-group placement and the exclusion policy
  as contract; verify against the OpenMMSystemBuilder call site.

## 2. Scope

- **Files:** `bridge/ForceFactory.{hpp,cpp}` (extracted from `OpenMMContext.cpp`).
- **Public symbols:** the `create*Force` builders and `addStandardExclusions`.
- **Known gaps to close:** state the force-group convention explicitly (it drives the
  separate-force-groups toggle and the validation band). Do not narrate the OpenMM
  API calls; state what forces exist and in which group after each builder runs.

## 3. Evidence pointers (tests exercising the module)

- TestAlchemy - the only test constructing a real OpenMM `Context` (TESTS.md section 1).
  Otherwise exercised by the Level-1 run and its force-group validation (VERIFY B3).

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; force-group/exclusion
  contracts stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-ForceFactory-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
