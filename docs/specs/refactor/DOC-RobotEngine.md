# DOC-RobotEngine: articulated-body solver (ABA recursions + mass ops)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the SimTK-free articulated-body solver - kinematics, articulated-body
  inertia recursion, forward dynamics (`calcUDot`), the mass-matrix operators
  (`multiplyBySqrtM`/`multiplyByMInv`/`calcKineticEnergy`/`calcLogDetM`), and reaction
  forces (ARCHITECTURE section 3, section 7). A **namespace of ~20 static methods** over
  `(const RobotModel&, RobotState&, ...)` - a structure-of-arrays driver, not an object
  (ARCHITECTURE section 4).
- **Layer:** Algorithms / dynamics (ARCHITECTURE section 2).
- **Ownership:** owns no state; reads the immutable model, writes state cache fields.
- **Invariants (HYPOTHESES):**
  - **INV-4 realization order** (ARCHITECTURE section 5): the recursions write caches in
    stage order position -> velocity -> ABA inertias -> udot. Document each method's
    stage - which caches it requires valid on entry and which it makes valid on exit.
    This is the write-side counterpart of the RobotState accessor contracts
    (DOC-RobotState); the two SHALL agree.
  - **INV-5 mass metric** (ARCHITECTURE section 5): momentum seeding and kinetic energy use
    the **same** mass-matrix operators; `calcLogDetM` is the Fixman kinetic term. State
    each operator's exact metric contract.
- After `SPLIT-R1`/`R2`/`R3`/`R4`: the dense solver (hinge_linalg), the NaN scanner
  (robo_debug), and the joint switches (JointKernels) have left; the class is one
  logical unit across `RobotEngine.cpp` + `RobotEngine_massops.cpp` +
  `RobotEngine_reaction.cpp` (cohesive TU split, same class - document at the header
  declaration, once per symbol).

## 2. Scope

- **Files:** `dynamics/RobotEngine.{hpp,cpp}`, `RobotEngine_massops.cpp`,
  `RobotEngine_reaction.cpp`.
- **Public symbols:** the kinematics/inertia/forward-dynamics recursions, the mass
  operators, the reaction-force methods.
- **Known gaps to close:** every method's stage precondition/postcondition (INV-4) -
  recovered from the realization sequence, the source of truth for DOC-RobotState.
  Keep the two tickets consistent; a disagreement is a finding in both.

## 3. Evidence pointers (tests exercising the module)

- TestRoboticsOracle, TestRoboticsOracleMolecule (differential vs Simbody, FAST\*),
  TestMassMatrix, TestKineticEnergy, TestInertia, TestMobilizerKinematics,
  TestSpatialAlgebra - contract (INV-4/INV-5; TESTS.md section 2/section 3). The oracles are the
  crown-jewel differential contract.

## 4. Exit criteria

- Doxygen warning-free; every public method documented with its stage in/out (INV-4)
  and, for mass ops, the metric (INV-5); comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-RobotEngine-findings.md` present; any stage disagreement with
  DOC-RobotState recorded; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
