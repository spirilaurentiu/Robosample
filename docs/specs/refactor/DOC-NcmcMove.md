# DOC-NcmcMove: nonequilibrium candidate Monte Carlo move

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `ncmcMove`, `ncmcInnerGhmcStep`, and the trough teleport - the NCMC
  alchemical-switching move: decouple -> drive -> recouple with GHMC inner steps,
  accepting on accumulated nonequilibrium work (ARCHITECTURE section 7; MODULES.md W4).
- **Layer:** Domain / world (ARCHITECTURE section 2).
- **Ownership:** drives the World state through the lambda schedule; uses AlchemyForce
  Factory `setAlchemicalLambda` (DOC-AlchemyForceFactory) and `protocolLambda` (DOC-NCMCProtocol).
- **Invariants (HYPOTHESES):**
  - NCMC acceptance uses **accumulated protocol work**, which is path-dependent; the
    work accounting SHALL match the lambda schedule and the alchemy force's lambda->energy map.
    Document what work is accumulated, where, and the acceptance criterion. A
    mismatch between the switching path and the work integral biases acceptance
    (Critical; MEMORY: NCMC explicit-solvent campaign traced near-zero acceptance to
    integrator shadow work `dt^2*nDOF`, not reorganization).
  - The inner GHMC step is reversible and preserves the intermediate distribution at
    each lambda; state that requirement.
  - The trough teleport is a discrete proposal - document its proposal and acceptance
    correction separately.

## 2. Scope

- **Files:** `world/sampler/NcmcMove.{hpp,cpp}` (from `World.cpp`).
- **Public symbols:** `ncmcMove`, `ncmcInnerGhmcStep`, the trough-teleport entry.
- **Known gaps to close:** state the work-accumulation and acceptance contract
  precisely, and the shadow-work sensitivity noted above (context for the reviewer,
  not a claim the doc resolves).

## 3. Evidence pointers (tests exercising the module)

- TestNcmcTeleport (SLOW), TestNCMCWork tests 9-16 (SLOW work chains), TestNcmc
  ExplicitSolvent (largely stubbed (3 permanent stubs, 6/8 skip at default tier) - not evidence; TESTS.md section 3, D-T2). Contract for the work
  accounting comes from TestNCMCWork; record the explicit-solvent gap in findings.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the work-accounting and
  acceptance contract stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-NcmcMove-findings.md` present; the explicit-solvent coverage gap
  recorded; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
