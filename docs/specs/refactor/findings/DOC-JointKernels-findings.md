# DOC-JointKernels findings

## OQ-2 (OPEN-QUESTION, first-class): untested HDot_FM cases
- `jointHDot_FM` for joint types **BendStretch**, **SphericalCoords**, and
  **FreeLine** is faithful to the Simbody derivation on paper but has NO golden
  test against a single-body Simbody reference in this port. These three feed the
  mobilizer-bias (Coriolis) acceleration, the historical home of silent
  energy-pumping bugs. Documented with an explicit `@note Assumed:` on
  `jointHDot_FM` (untested), not as a confident contract. Recommended oracle:
  energy conservation on a free body of each of the three types.
- The other seven joint types have `H_FM` constant in F, so `jointHDot_FM` is
  identically zero for them and the engine does not call it; documented as
  contract.

## Verified against call sites (dispatch is the single source)
- Callers of the kernels: `src/RobotEngine_kinematics.cpp`,
  `src/world/FixmanCorrection.cpp`, `include/RobotIntegrator.hpp`, and tests
  (`TestIntegrator`, `TestStability`, `TestEnsembleOrientation`,
  `tests/EngineHelpers.hpp`). No residual per-joint switch was found outside
  `JointKernels.{hpp,cpp}` for these six kernels; the post-R4 dispatch is the
  single source (verified by grep for the kernel names).
- Joint set documented from the `JointType` enum
  (`include/RobotModel.hpp:27`): Rigid, Torsion, Slider, Cylinder, BendStretch,
  Cartesian, Ball, SphericalCoords, FreeLine, Free.

## INV-9 (documented as contract)
- `advanceQuatExp` produces a unit quaternion by construction and is reversible
  under angular-velocity negation; `jointDriftQuat` delegates to it at the
  midpoint velocity. The parent-frame angular-velocity convention matches
  `quaternionDotFromAngVel` (robot_math). Documented as `@post` on both.

## Notes
- The `@note Assumed:` on `jointHDot_FM` is matched by this OQ-2 entry.
