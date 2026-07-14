# DOC-DrivenRexDriver: RENE / REBASONTOP driven-REX (uncompiled; OQ-5)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `runDrivenRound`, `driveReplica`, `runInterleavedRemcSubround` - the
  driven replica-exchange path (RENE / REBASONTOP) that interleaves BAT-drive/NMA
  nonequilibrium moves with exchange (ARCHITECTURE section 3; MODULES.md C5).
- **Layer:** Workflow (ARCHITECTURE section 2).
- **Ownership:** as ReplicaExchangeDriver, over the driven round structure.
- **Invariants (HYPOTHESES) - treat cautiously:**
  - **OQ-5 (flagged): this code is documented as reviewed-on-paper but never
    compiled/run** (ARCHITECTURE OQ-5), and RENEMC's round-loop throws (TESTS.md section 7).
    It may sit behind a feature flag until it has a test. **Do not write confident
    behavioral contracts.** Document the intended structure with explicit
    `@note Assumed:` (uncompiled, untested) and route the whole module to findings as
    an OPEN-QUESTION.
  - Where it invokes SwapAcceptance's RENE/REBASONTOP terms (INV-7-adjacent, INV-8) and the
    BAT drive (DOC-BatScaling), note the intended coupling but mark the end-to-end
    behavior as unverified.

## 2. Scope

- **Files:** `workflow/rex/DrivenRexDriver.{hpp,cpp}` (from `Context.cpp`).
- **Public symbols:** `runDrivenRound`, `driveReplica`, `runInterleavedRemcSubround`.
- **Known gaps to close:** the entire module is the OQ-5 gap - every contract here is
  `@note Assumed:` + findings until the code compiles and a test exists. Do not
  transcribe the on-paper design as established behavior.

## 3. Evidence pointers (tests exercising the module)

- None (uncompiled; RENEMC round-loop throws - TESTS.md section 7). Record the total absence
  of runnable evidence in findings; this is the defining characteristic of the ticket.

## 4. Exit criteria

- Doxygen warning-free; every public symbol carries at least an `@note Assumed:`
  (uncompiled/untested) contract; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-DrivenRexDriver-findings.md` present; OQ-5 recorded as an
  OPEN-QUESTION; **every** `@note Assumed:` has a matching findings entry.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
