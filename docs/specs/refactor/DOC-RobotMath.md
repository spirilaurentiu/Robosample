# DOC-RobotMath: spatial-algebra value types

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the value PODs of spatial algebra - `Vec3`, `Rotation`, `Transform`,
  `Quat`, `SpatialVec`, `SpatialInertia`, `ArticulatedInertia`, `PhiMatrix`,
  `Inertia`, `MassProperties`. Breadth, not depth (ARCHITECTURE section 7).
- **Layer:** Infrastructure (ARCHITECTURE section 2). No downward engine dependency; depends
  only on `<cmath>`/std. Everything above uses it.
- **Ownership:** none - value types, copied by value, no heap, no aliasing concerns.
- **Invariants (HYPOTHESES):** Rotation/Transform compose as rigid-body frame maps;
  `Quat` operations respect the double cover (INV-9, ARCHITECTURE section 5) - the
  normalize/exp-map helpers here are what the joint code relies on. Verify against
  call sites; contradiction -> findings.

## 2. Scope

- **Files:** `include/robot_math.hpp` (867 LOC). If the optional split landed
  (MODULES.md section 1, sequenced last / may be deferred), document both
  `math/robot_vecmat.hpp` and `math/robot_spatial.hpp`; otherwise the single header.
  Confirm which files exist before starting.
- **Public symbols:** every exported type and its member/free operators.
- **Known gaps to close:** `Inertia` vs `SpatialInertia` vs `MassProperties` overlap
  (OQ-4). Do **not** assert they are interchangeable; document each type's actual
  contract at its call sites and route the overlap question to findings - a
  whole-repo usage scan is a human decision, not a doc claim.

## 3. Evidence pointers (tests exercising the module)

- TestVectorMath, TestTransform, TestOrientation, TestRotationConstruction,
  TestSpatialAlgebra, TestQuaternion, TestInertia - all FAST contract/algebra
  (TESTS.md section 2). These pin the composition and quaternion contracts.

## 4. Exit criteria

- Doxygen warning-free on touched files (warnings-as-errors).
- Every public symbol in scope documented; contracts state behavior, not the
  arithmetic that computes it (a rewrite SHALL NOT falsify them).
- Comment-stripped before/after diff empty (zero code change; VERIFY section 4).
- `findings/DOC-RobotMath-findings.md` present; the OQ-4 overlap recorded there; no
  `@note Assumed:` without a matching entry.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form.
Where an exemplar and `documenter.md` conflict, the spec wins and the conflict is a
finding.
