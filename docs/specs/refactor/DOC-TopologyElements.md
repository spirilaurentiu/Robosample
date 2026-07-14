# DOC-TopologyElements: SystemTopology SoA payload

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `SystemTopology`, the structure-of-arrays payload (~110 fields) that
  Python fills in place and hands to the engine - the Python<->engine boundary
  (ARCHITECTURE section 4, section 7).
- **Layer:** Domain data / model (ARCHITECTURE section 2).
- **Ownership:** owned by `Context` as a public member for the whole run
  (ARCHITECTURE section 4); filled in place from Python, then read by `World::buildModel`.
  The producer (Python `context.py`) and consumer (`ModelBuilder`) together define
  each field's contract.
- **Invariants (HYPOTHESES):** field arrays are parallel and indexed by a common
  atom/bond index; units follow INV-3 (nm). Document each field group's meaning,
  index space, unit, and which side (Python vs engine) writes it. Where a field's
  fill site and read site disagree on meaning, that is a finding.

## 2. Scope

- **Files:** `model/TopologyElements.hpp` (`SystemTopology`). The sub-struct split is
  **deferred** behind a decision record (MODULES.md section 5); document the god-struct as it
  stands, do not anticipate the split.
- **Public symbols:** `SystemTopology` and its field groups; any accessors/builders.
- **Known gaps to close:** this is the largest single-struct surface; the risk is
  transcribing field names instead of recovering meaning. Every field group gets a
  producer/consumer-verified contract or a findings entry - not a name restatement
  (`documenter.md` section 8 forbids restating the name).

## 3. Evidence pointers (tests exercising the module)

- TestBuilders, TestTransfer, TestAtomTransfer - FAST contract (TESTS.md section 2). Plus
  `ModelBuilder` call sites (build-time consumer) and Python `context.py` (producer).

## 4. Exit criteria

- Doxygen warning-free; every public field/symbol in scope documented with meaning +
  index space + unit + writer; contracts behavioral; comment-stripped diff empty
  (VERIFY section 4).
- `findings/DOC-TopologyElements-findings.md` present; producer/consumer
  disagreements recorded; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
