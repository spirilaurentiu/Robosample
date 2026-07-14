# DOC-LegacyCoordSwapRex: runREX coordinate-swap oracle (OQ-3)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).
**Conditional:** exists only if the retire-or-keep decision (OQ-3) keeps `runREX`
and its producer **SPLIT-C6** (`SPLIT-C6.md`) ran - C6 isolates `runREX` into
`workflow/rex/LegacyCoordSwapRex.{hpp,cpp}`. OQ-3 is **DECIDED: keep as oracle**
(ARCHITECTURE OQ-3), so this ticket is expected to apply. If a human instead
retires `runREX`, C6 and this DOC ticket are both dropped. If C6 is skipped but
`runREX` is kept, `runREX` lives in `ReplicaExchangeDriver` (per C4) and this
ticket's content folds into `DOC-ReplicaExchangeDriver` instead.

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `runREX` - the coordinate-swap replica-exchange driver, deliberately
  **frozen** as the INVARIANT-EQUIV differential oracle for the label-swap `RunREX`
  (ARCHITECTURE section 3, OQ-3; `Context.hpp:84-90` "do not extend"). It is not deprecated
  junk; it is the correctness baseline.
- **Layer:** Workflow (ARCHITECTURE section 2).
- **Ownership:** as the REX driver, but swaps **coordinates** between replicas rather
  than labels.
- **Invariants (HYPOTHESES):**
  - **INV-8** (ARCHITECTURE section 5): `runREX` and `RunREX` produce equivalent sampling;
    `runREX` **is** the differential oracle for `RunREX`. Document exactly this role -
    the contract is "equivalent-sampling reference," and the "do not extend" status is
    part of the contract.
  - It swaps coordinates (the inter-world currency, INV-3) between replicas; document
    the swap mechanism's observable effect as the reference behavior `RunREX` SHALL
    match.

## 2. Scope

- **Files:** `workflow/rex/LegacyCoordSwapRex.{hpp,cpp}` (from `Context.cpp`).
- **Public symbols:** `runREX` and its swap helpers.
- **Known gaps to close:** document the "frozen oracle / do not extend" status
  prominently as `@note`; a future extension would break its role as the differential
  baseline. State that the driver scripts historically call this at `NOF_REPLICAS = 1`
  (ARCHITECTURE section 1) as context, not as a recommendation.

## 3. Evidence pointers (tests exercising the module)

- `test_rex_label_swap_equivalence.py` - the label-swap-vs-coordinate-swap
  INVARIANT-EQUIV equivalence (Python; out of C++ scope but the oracle-relationship
  evidence, TESTS.md section 7).

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the frozen-oracle role (INV-8,
  "do not extend") stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-LegacyCoordSwapRex-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
