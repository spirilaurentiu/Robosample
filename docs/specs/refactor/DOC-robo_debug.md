# DOC-robo_debug: guarded NaN/Inf scanner

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the `ROBO_DEBUG`-guarded NaN/Inf scanner used to trap non-finite
  state during solver development (ARCHITECTURE section 7, MODULES.md R2, ~150 LOC).
- **Layer:** Infrastructure (ARCHITECTURE section 2). Diagnostic only; no effect on sampled
  distributions when the guard is off.
- **Ownership:** none - inspects borrowed buffers/state, produces diagnostics
  (abort/log). Confirm whether it aborts, throws, or logs, and on which stream.
- **Invariants (HYPOTHESES):** the scanner is compiled out entirely unless
  `ROBO_DEBUG` is defined; production sampling behavior is identical with it absent.
  Document this as the release-vs-debug behavioral difference (`documenter.md` section 6 for
  conditionally compiled code); verify against the build guards.

## 2. Scope

- **Files:** `math/robo_debug.hpp` (extracted from `RobotEngine.cpp` by `SPLIT-R2`).
- **Public symbols:** the scan entry points and any `ROBO_DEBUG` macros.
- **Known gaps to close:** no dedicated test exercises this module (guarded, dev
  only). Document from the guard sites and its call sites in the solver; where no
  call-site evidence exists, prefer omission over invention, or `@note Assumed:` +
  findings.

## 3. Evidence pointers (tests exercising the module)

- None dedicated. Call sites are inside the dynamics `.cpp` under `#ifdef
  ROBO_DEBUG`. Record the absence of test coverage in findings (a coverage gap, not
  a doc failure).

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; behavioral contracts only.
- Comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-robo_debug-findings.md` present (records the no-test gap); no
  `@note Assumed:` without a matching entry.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
