# DOC-Replica: Replica + ThermodynamicState value types

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).
Documented **before** the REX drivers that use them (leaf-first).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the `Replica` and `ThermodynamicState` value types underpinning
  label-swap replica exchange - a replica holds a configuration/label; a thermodynamic
  state holds a temperature/ensemble parameter; the driver swaps labels via two
  inverse index maps (ARCHITECTURE section 3; MODULES.md C4 / `rex/Replica.hpp`).
- **Layer:** Workflow (ARCHITECTURE section 2).
- **Ownership:** value/state types owned by the REX driver; document what each holds
  and its lifetime within a REX run.
- **Invariants (HYPOTHESES):**
  - The label-swap design keeps configurations in place and swaps **labels**, so the
    two index maps (`replica->state`, `state->replica`) are mutually inverse at all
    times (ARCHITECTURE section 3). Document that inverse-map invariant on the types that
    carry the mapping.
  - A `ThermodynamicState` defines the target Boltzmann distribution for a replica at
    its temperature (REX-spec INV-1/INV-2, TESTS.md section 7); document what parameters it fixes.

## 2. Scope

- **Files:** `workflow/rex/Replica.hpp`.
- **Public symbols:** `Replica`, `ThermodynamicState`, and the index-map accessors if
  co-located here.
- **Known gaps to close:** state the inverse-map invariant and the meaning of each
  field; these types are the vocabulary the driver/acceptance tickets build on.

## 3. Evidence pointers (tests exercising the module)

- TestRexAcceptanceAlgebra fixture (constructs states for the acceptance formula;
  TESTS.md section 7). The behavioral REX oracle (TESTS.md section 7, pending) will exercise these
  end to end.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the inverse-map and
  thermodynamic-state contracts stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-Replica-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
