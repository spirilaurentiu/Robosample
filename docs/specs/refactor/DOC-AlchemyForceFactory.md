# DOC-AlchemyForceFactory: alchemical decoupling forces

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the alchemy/decoupling force construction and its `enableAlchemy`/`setAlchemicalLambda`
  control surface - the OpenMM side of NCMC alchemical switching (ARCHITECTURE section 7,
  MODULES.md O4). NCMC-only.
- **Layer:** Bridge (ARCHITECTURE section 2).
- **Ownership:** constructs decoupling forces and transfers them to the OpenMM
  `System`; holds handles/indices used later by `setAlchemicalLambda`. Recover the handle
  lifetime and who invalidates it.
- **Invariants (HYPOTHESES):** `setAlchemicalLambda` drives the same lambda the NCMC protocol
  (`protocolLambda`, DOC-NCMCProtocol) schedules; the lambda->force-scaling map here SHALL
  match the work-accounting the NCMC move assumes (the NCMC work is
  path-dependent - a mismatch corrupts the acceptance, MEMORY: NCMC campaign).
  Document the lambda semantics (endpoints lambda=0/lambda=1 meaning, monotonicity) as contract;
  verify against `AlchemyForceFactory`'s caller in `NcmcMove`.

## 2. Scope

- **Files:** `bridge/AlchemyForceFactory.{hpp,cpp}` (extracted from
  `OpenMMContext.cpp`).
- **Public symbols:** the decoupling-force builders and the `enableAlchemy`/`setAlchemicalLambda`
  surface.
- **Known gaps to close:** state precisely what lambda scales (which interactions
  decouple) and the coupled/decoupled endpoint convention; this is the contract NCMC
  correctness rests on, and it is not a mechanism detail.

## 3. Evidence pointers (tests exercising the module)

- TestAlchemy - real OpenMM `Context` (TESTS.md section 1). TestNCMCWork - NCMC work chains
  (SLOW; TESTS.md section 2) exercise the lambda path. Note TestNcmcExplicitSolvent is largely stubbed (3 permanent stubs, 6/8 skip at default tier)
  (TESTS.md section 3) - do not treat its skipped cases as evidence.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; lambda semantics stated as
  contract; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-AlchemyForceFactory-findings.md` present; no unmatched
  `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
