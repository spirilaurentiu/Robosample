# DOC-MTSIntegrator: r-RESPA multiple-time-step integrator

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the r-RESPA multiple-time-step substep integrator (ARCHITECTURE section 7,
  MODULES.md O1), extracted from `OpenMMContext.cpp`.
- **Layer:** Bridge (ARCHITECTURE section 2) - it is an OpenMM integrator subclass/wrapper.
- **Ownership:** owned/held by the OpenMM `Context` per the OpenMM integrator
  contract. Recover which forces run at the inner vs outer substep from the
  construction/force-group call sites.
- **Invariants (HYPOTHESES):** the force-group split across substeps SHALL be
  reversible and preserve the sampled distribution (MTS was ruled out for the
  contact-world campaign, MEMORY: two-robot-contact - but it exists as an option).
  Document the substep decomposition (which force groups at which frequency) and the
  reversibility contract; verify against the force-group assignment in ForceFactory.

## 2. Scope

- **Files:** `bridge/MTSIntegrator.{hpp,cpp}` (extracted from `OpenMMContext.cpp`).
- **Public symbols:** the integrator type and its step/config surface.
- **Known gaps to close:** confirm whether this integrator is on any runnable
  production path or is dormant configuration; record which. Do not assert it is used
  if no call site selects it.

## 3. Evidence pointers (tests exercising the module)

- No dedicated test. Record the coverage gap in findings; recover the contract from
  the construction site and the OpenMM integrator base contract.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; substep decomposition and
  reversibility stated as contract; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-MTSIntegrator-findings.md` present; the no-test gap and
  used/dormant status recorded; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
