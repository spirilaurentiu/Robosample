# DOC-SwapAcceptance: REX swap acceptance algebra

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `attemptREXSwap` and its guards - the swap-acceptance algebra covering
  the four run-types (REMC / RENEMC / RENE / REBASONTOP), extracted with the driver
  (ARCHITECTURE section 3; MODULES.md C4). The theory-derived, tested acceptance formula.
- **Layer:** Workflow (ARCHITECTURE section 2).
- **Ownership:** stateless algebra over two `ThermodynamicState`s / `Replica`s
  (DOC-Replica); owns nothing.
- **Invariants (HYPOTHESES):**
  - **INV-8 REX detailed balance** (ARCHITECTURE section 5): the acceptance ratio satisfies
    detailed balance; the guards **reject on NaN/inf**. Document each run-type's
    acceptance criterion and the guard behavior (domain-error sentinel).
  - **INV-7-adjacent**: for the BAT-drive run-types, acceptance uses the shared anchor
    snapshot so the paired map is an exact involution (DOC-BatScaling); RENEMC uses an
    `ETerm_nonequil` term. Document each term's role per run-type; a sign or missing
    term biases acceptance (Critical).
  - The formula is Jacobian-sign-flip *blind* by construction where the test asserts it
    (TestRexAcceptanceAlgebra) - document the property the test defends.

## 2. Scope

- **Files:** `workflow/rex/SwapAcceptance.{hpp,cpp}` (`attemptREXSwap`,
  `checkInv7AndInv10Guards`, made public for the test - TESTS.md section 3).
- **Public symbols:** `attemptREXSwap`, the guard entry points, per-run-type
  acceptance helpers.
- **Known gaps to close:** document the four run-types' acceptance criteria separately
  and precisely (sign conventions, the nonequil term); this is the algebra the driver
  wires and the highest-consequence contract in workflow.

## 3. Evidence pointers (tests exercising the module)

- TestRexAcceptanceAlgebra (7 tests) - FAST contract, INV-8: detailed balance of the
  ratio, Jacobian sign-flip blindness, RENEMC `ETerm_nonequil`, domain-error sentinel
  (TESTS.md section 7). `test_rex_swap_acceptance_algebra.py` mirrors it (spec V8). These are
  the direct, theory-derived contract evidence.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; per-run-type acceptance +
  INV-8 guards stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-SwapAcceptance-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
