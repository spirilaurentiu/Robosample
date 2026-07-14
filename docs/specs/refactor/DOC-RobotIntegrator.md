# DOC-RobotIntegrator: generalized-coordinate leapfrog HMC integrator

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the HMC integrator templated on the `Bridge` type - `stepTo` /
  `verletStep` running leapfrog in generalized coordinates (realizePosition ->
  realizeVelocity -> realizeArticulatedBodyInertias -> calcUDot -> position drift with
  quaternion exp-map -> SHAKE projection -> implicit-trapezoid velocity correction),
  plus `checkReversibility` (ARCHITECTURE section 2, section 3; MODULES.md I1). After `SPLIT-I1`,
  `driftPositions`/`velocityCorrector`/`cartesianSolventVerlet` are private helpers.
- **Layer:** Algorithms / dynamics (ARCHITECTURE section 2).
- **Ownership:** owns no persistent state; drives the borrowed `RobotState` through
  the realization stages via `RobotEngine`; obtains forces from the templated
  `Bridge`.
- **Invariants (HYPOTHESES):**
  - **INV-4 realization order** - the integrator is the canonical stage sequencer;
    document the exact order it enforces (this is where OQ-1's convention lives).
  - **INV-6 constraint consistency** - SHAKE position projection and the velocity
    correction use the same `G` as the Fixman log-det (Constraints).
  - **INV-9 quaternion double cover** - the exp-map drift and reversibility checks
    normalize quaternions and account for the double cover (ARCHITECTURE section 5).
  - It is templated on `Bridge` so OpenMM stays out of the header (ARCHITECTURE section 6.1);
    document the `Bridge` template requirements from all instantiations
    (`documenter.md` section 6), not from OpenMM alone.

## 2. Scope

- **Files:** `dynamics/RobotIntegrator.{hpp,cpp}` (the ~350/133-line integrator
  templates; I1 extracts the three privates).
- **Public symbols:** `stepTo`, `verletStep`, `checkReversibility`; the `Bridge`
  template parameter.
- **Known gaps to close:** document the leapfrog sequence as the **contract**
  (stages, projection, correction, reversibility guarantee), never as a
  step-by-step narration; a behavior-preserving reorder of internal helper calls SHALL
  not falsify the text.

## 3. Evidence pointers (tests exercising the module)

- TestIntegrator (drift/reversibility chains, **characterization** - label as such,
  TESTS.md section 3), TestReverseMobilizer, TestFreeJointKEPump (characterization),
  TestStability, TestMassScaleInvariance - contract/char (TESTS.md section 2/section 3).

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; INV-4/INV-6/INV-9 stated;
  characterization-test evidence labeled as such (not promoted to contract);
  comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-RobotIntegrator-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
