# DOC-StartupValidator: startup geometry sanity scan

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `checkStartupGeometry` - the O(N^2) clash / potential-energy sanity scan
  run during `Context::initialize` before sampling begins (ARCHITECTURE section 3, section 7;
  MODULES.md C1). Pure, high independence.
- **Layer:** Workflow (ARCHITECTURE section 2).
- **Ownership:** reads initial coordinates/topology; produces a pass/fail (or warning)
  result. Owns nothing persistent.
- **Invariants (HYPOTHESES):** it is a **precondition check**, not part of the sampled
  chain - it never mutates sampled state; a failure aborts/warns before any move.
  Document the pass/fail criterion (clash threshold, PE bound) and the failure action
  (abort vs warn) as contract; verify against the `initialize` call site.

## 2. Scope

- **Files:** `workflow/StartupValidator.{hpp,cpp}` (from `Context.cpp` ~237-395).
- **Public symbols:** `checkStartupGeometry` and any threshold accessors.
- **Known gaps to close:** state the exact criterion and the failure action; the
  O(N^2) cost is mechanism, not contract - do not document it as a guarantee.

## 3. Evidence pointers (tests exercising the module)

- No dedicated test found; TestStability is adjacent evidence for the clash criterion.
  Record the coverage gap in findings; recover the criterion from the source and the
  `initialize` caller.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; pass/fail criterion + failure
  action stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-StartupValidator-findings.md` present; the coverage gap recorded; no
  unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
