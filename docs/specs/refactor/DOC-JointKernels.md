# DOC-JointKernels: per-joint-type kinematic kernels

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the centralized per-joint-type kernels - `X_FM`, `H_FM`, `HDot_FM`,
  `QDot`, `QDotDot`, and the quaternion drift branch - after `SPLIT-R4` gathers the
  four scattered joint switches into one place (ARCHITECTURE section 6.6, section 7; MODULES.md
  R4). Stateless free functions over `(const RobotModel&, RobotState&, ...)`
  (ARCHITECTURE section 4).
- **Layer:** Algorithms / dynamics (ARCHITECTURE section 2).
- **Ownership:** none - read the model, write cache fields on the state.
- **Invariants (HYPOTHESES):**
  - Joint-type dispatch is consistent across all six kernels (the R4 whole point:
    adding a joint touches one file, ARCHITECTURE section 6.6). Document the joint-type set
    and each kernel's contract per type.
  - **INV-9 quaternion double cover** (ARCHITECTURE section 5): free/quaternion joints
    normalize each step and reversibility accounts for the double cover. Document the
    normalization as a `@post` on the relevant kernels.
  - **OQ-2 (flagged) - UNVALIDATED guarantee.** `jointHDot_FM` for **BendStretch**,
    **SphericalCoords**, and **FreeLine** is "faithful on paper but NOT golden-tested"
    (ARCHITECTURE section 8, OQ-2). Do **not** write a confident contract for these three
    `HDot_FM` cases. Document them with an explicit `@note Assumed:` that the
    guarantee is untested, and **route the gap to findings** as a first-class
    OPEN-QUESTION. The other joint types' `HDot_FM` may be documented normally if call
    sites/tests agree.

## 2. Scope

- **Files:** `dynamics/JointKernels.{hpp,cpp}` (+ the switches R4 moved out of
  `calcQDot`/`calcQDotDot`/`verletStep` drift branch).
- **Public symbols:** `X_FM`, `H_FM`, `HDot_FM`, `QDot`, `QDotDot`, drift, per joint
  type.
- **Known gaps to close:** the OQ-2 three `HDot_FM` cases above; and confirm the
  post-R4 dispatch is the single source (a residual switch elsewhere is a finding).

## 3. Evidence pointers (tests exercising the module)

- TestJointKernels, TestMobilizer, TestReverseMobilizer, TestMobilizerKinematics,
  TestQuaternion - FAST contract (TESTS.md section 2). These are the golden coverage; the
  **absence** of a golden test for the three OQ-2 `HDot_FM` cases is exactly the gap.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; INV-9 stated; the OQ-2 three
  cases marked `@note Assumed:` (untested), not asserted; comment-stripped diff empty
  (VERIFY section 4).
- `findings/DOC-JointKernels-findings.md` present; OQ-2 recorded as an OPEN-QUESTION
  with the three joint types named; every `@note Assumed:` has a matching entry.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
