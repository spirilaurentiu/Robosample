# DOC-NCMCProtocol: lambda protocol schedule

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `protocolLambda`, the single free function giving the NCMC switching
  schedule lambda(step) (ARCHITECTURE section 4, section 7; MODULES.md section 5: single-concept, no split).
- **Layer:** Algorithms / dynamics (ARCHITECTURE section 2).
- **Ownership:** stateless pure function.
- **Invariants (HYPOTHESES):** lambda runs the coupled<->decoupled endpoints the alchemy
  forces expect (DOC-AlchemyForceFactory); the schedule's endpoint values and
  monotonicity are the contract the NCMC work accounting assumes. A schedule/endpoint
  mismatch corrupts NCMC acceptance (MEMORY: NCMC campaign - work is path-dependent).
  Document lambda(0), lambda(N), and monotonicity as contract; verify against the `NcmcMove`
  caller.

## 2. Scope

- **Files:** `dynamics/NCMCProtocol.hpp` (header-only; moved, not split).
- **Public symbols:** `protocolLambda`.
- **Known gaps to close:** state the endpoint convention (which end is coupled) and
  the step->lambda mapping explicitly; this is the seam where the protocol and the alchemy
  force SHALL agree.

## 3. Evidence pointers (tests exercising the module)

- TestNCMCWork tests 1-8 - FAST algebra covering `protocolLambda` in isolation
  (TESTS.md section 2; the FAST/SLOW split of that file). These are the direct contract
  evidence.

## 4. Exit criteria

- Doxygen warning-free; `protocolLambda` documented; endpoint/monotonicity contract
  stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-NCMCProtocol-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
