# SPLIT-R4: centralize the joint taxonomy into JointKernels

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)

`ARCHITECTURE.md section 6.6` records a dependency defect: the per-joint-type `switch`
statements appear at FOUR sites instead of being centralized in `JointKernels`.
Adding a joint touches four files. The four sites are:

1. `src/JointKernels.cpp` - `jointX_FM` / `jointH_FM` / `jointHDot_FM` (already here).
2. `RobotEngine::calcQDot` - the `switch (m.bodyJoint[b])` mapping u -> qdot
   (`RobotEngine.cpp`, lines 738-786, cases Ball / FreeLine / Free / default).
3. `RobotEngine::calcQDotDot` - the `if/else` on `JointType` mapping udot -> qddot
   (`RobotEngine.cpp`, lines 1409-1451, branches Ball / FreeLine / Free / else).
4. The quaternion-drift branch - the per-quaternion-body `wHalf` construction plus
   `advanceQuatExp`, and `advanceQuatExp` itself (`RobotIntegrator.hpp`, lines 76-100).
   ORDERING: **R4 runs after I1** (and after R3). SPLIT-I1 relocates the drift loop
   (originally `verletStep` lines 241-258) into `RobotIntegrator::driftPositions`, so
   by the time R4 runs the `wHalf`/`advanceQuatExp` block lives inside
   `driftPositions`, not `verletStep`. R4 lifts it from there. Re-anchor by symbol
   (`driftPositions`), not by the original `verletStep` line numbers.

This ticket moves sites 2, 3, and 4 into `JointKernels`, leaving one file that owns
the taxonomy. It is pure code motion of already-correct logic - the `case` bodies
move byte-for-byte and are wrapped in a per-body entry point that `calcQDot` /
`calcQDotDot` / `verletStep` call once per body. Nothing about the math changes.

Because R4 relocates joint-type logic, it is GATED ON OQ-2 (`ARCHITECTURE.md section 9`):
three joint types - BendStretch, SphericalCoords, FreeLine - have `jointHDot_FM`
implementations that are "faithful on paper but NOT golden-tested" (`JointKernels.cpp`
lines 162-167). OQ-2 asks whether a golden test is a prerequisite for consolidation.
Per the decision policy this ticket carries the gap forward unchanged: R4 moves
existing logic and does not touch `jointHDot_FM`, so it neither closes nor widens the
gap. The untested-guarantee note SHALL be preserved verbatim on `jointHDot_FM`, and
this ticket SHALL NOT be marked complete until a human confirms OQ-2's "ship as pure
motion, gap carried forward" resolution (or supplies the golden test first).

Invariants from `ARCHITECTURE.md section 5`, restated in full:

- **INV-9 Quaternion double cover.** Free/quaternion joints normalize quaternions
  each step; reversibility checks account for the double cover. The drift branch
  moved here (site 4) contains the exact exponential-map advance and the
  `normalizeQuaternions` call ordering; both are preserved.
- **INV-4 Realization order.** `RobotState` caches are valid only in stage order
  position -> velocity -> articulated-body inertias -> udot. `calcQDot` reads u and
  X_FM (velocity stage); `calcQDotDot` reads u, udot, q, X_FM (post-udot). The
  centralized entry points read exactly the same inputs; call order is unchanged.

## Moves (exactly what)

New symbols in `include/JointKernels.hpp` + `src/JointKernels.cpp`
(`namespace robo`):

- `jointQDot` - per-body u -> qdot. Body: the `switch (jt)` from `calcQDot`
  (`RobotEngine.cpp` lines 738-786), cases Ball / FreeLine / Free / default, moved
  BYTE-FOR-BYTE. Signature carries exactly what the case bodies read:
  `void jointQDot(JointType jt, const Real* q, int qOff, const Real* u, int uOff,
  int dof, const Transform& X_FM, Real* qdotOut)`.
- `jointQDotDot` - per-body udot -> qddot. Body: the `if/else` from `calcQDotDot`
  (`RobotEngine.cpp` lines 1409-1451), branches Ball / FreeLine / Free / else, moved
  BYTE-FOR-BYTE. Signature: `void jointQDotDot(JointType jt, const Real* q, int qOff,
  const Real* u, const Real* udot, int uOff, int dof, const Transform& X_FM,
  Real* qddOut)`.
- `advanceQuatExp` - the exponential-map unit-quaternion advance, moved verbatim from
  `RobotIntegrator.hpp` lines 76-100 (it was itself moved there from RobotEngine.cpp;
  its sole callers are the drift branch and - after this ticket - `jointDriftQuat`).
- `jointDriftQuat` - per-quaternion-body position drift. Body: the `wHalf`
  construction + `advanceQuatExp` call, which SPLIT-I1 has relocated from
  `verletStep` into `RobotIntegrator::driftPositions` (originally `verletStep`
  lines 241-258); R4 lifts it from `driftPositions`, moved BYTE-FOR-BYTE (the
  FreeLine `R_FM * (u0,u1,0)` special case and the Ball/Free `(u0,u1,u2)` case). Signature carries the four inputs the block reads:
  `void jointDriftQuat(JointType jt, const Real* q0, int qOff, const Real* u0,
  const Real* udot0, int uOff, const Transform& X_FM, Real h, Real* qOut)`.

Rewrites at the origin sites (call, do not inline):

- `calcQDot` (kinematics TU after R3): the per-body loop body becomes
  `jointQDot(m.bodyJoint[b], q, qOff, u, uOff, dof, X_FM[b], qdotOut)`. The loop and
  its index setup are unchanged.
- `calcQDotDot` (kinematics TU): the per-body body becomes
  `jointQDotDot(jt, q, qOff, u, udot, uOff, dof, X_FM[b], qdd)`.
- `driftPositions` (`RobotIntegrator.hpp`, the I1 helper): the `for (bodyIx ...)`
  quaternion-drift loop keeps its `isQuaternionBody` guard and its four-slot
  writeback, but the body becomes `jointDriftQuat(m.bodyJoint[bodyIx], &q0[qOff],
  qOff ..., s.X_FM()[bodyIx], h, qNew.data())`. The surrounding scalar Taylor drift,
  the `normalizeQuaternions` call, and the Cartesian-solvent drift are untouched.
  `advanceQuatExp`'s definition leaves `RobotIntegrator.hpp` and its include of
  `JointKernels.hpp` supplies it.

The moved `case`/branch bodies SHALL be diffed against the four current sites and
shown byte-identical (whitespace/comment normalized) - this is the one ticket that
centralizes a switch, so that equivalence is the primary review artifact.

`jointHDot_FM` and its untested-guarantee note (`JointKernels.cpp` 155-199,
`JointKernels.hpp` 19-32) are NOT modified.

## Public API after this ticket

`JointKernels.hpp` exports in `namespace robo`, added to the existing three:
`jointQDot`, `jointQDotDot`, `jointDriftQuat`, `advanceQuatExp`. Header self-contained
(`#pragma once`, `#include "RobotModel.hpp"`, `#include "robot_math.hpp"`). These are
NEW public engine-internal symbols but NOT part of the `robo_bindings...so` Python
surface; the nm delta is on internal engine symbols only. Because section 3 defect-resolution
tickets MAY add a thin symbol, and these add no Python-visible export, criterion 4 is
met without a separate API spec - but the added symbols SHALL be listed in the ticket's
nm-delta rationale.

## Constraints

- Behavior-preserving code motion. Zero math change. The moved `case` bodies are
  byte-identical; only their host function (a per-body free function vs an inline
  switch) changes.
- `RobotEngine.hpp` is unchanged (`calcQDot`/`calcQDotDot` keep their signatures).
- `RobotIntegrator.hpp` loses `advanceQuatExp` and the inline `wHalf` block, gains
  `#include "JointKernels.hpp"` (it already transitively sees `robo` math types).
- Invariants that SHALL remain true: **INV-9**, **INV-4** (restated above). The
  quaternion normalization ordering and the velocity/udot read stages are preserved.
- Gated on OQ-2: do not merge without the human-confirmed "gap carried forward"
  resolution; `jointHDot_FM` and its note stay verbatim.
- One commit for the move; formatting fixes separate.

## Predicted test breakage (Phase-A include fixes only)

- `tests/TestJointKernels.cpp` already exercises `JointKernels` directly; it may add
  direct coverage of the new `jointQDot`/`jointQDotDot`, but under the Phase-A freeze
  it is NOT edited unless it fails to compile. It calls the pre-existing
  `jointX_FM`/`jointH_FM`/`jointHDot_FM`, whose signatures are unchanged -> no fix.
- The ~14 tests naming `calcQDot`/`calcQDotDot` (`TestIntegrator.cpp`,
  `TestMobilizer.cpp`, `TestMobilizerKinematics.cpp`, `HmcDriver.hpp`,
  `TestFreeJointKEPump.cpp`, `TestNcmc*.cpp`, `TestBiasForces.cpp`,
  `TestTwoRobotContact.cpp`, `TestConstraints.cpp`, `TestMassScaleInvariance.cpp`,
  `TestRoboticsOracle.cpp`, `TestAtomTransfer.cpp`, `TestNCMCWork.cpp`) call the
  unchanged `RobotEngine::calcQDot`/`calcQDotDot` public methods -> no fix.
- `advanceQuatExp` is named by `tests/TestEnsembleOrientation.cpp`,
  `tests/TestStability.cpp`, `tests/EngineHelpers.hpp`. It moves namespace-visible
  from `RobotIntegrator.hpp` to `JointKernels.hpp`. If any of these three includes
  `RobotIntegrator.hpp` solely for `advanceQuatExp`, the ONE allowed Phase-A fix is
  to add `#include "JointKernels.hpp"` (or rely on the transitive include, since
  `RobotIntegrator.hpp` will include `JointKernels.hpp`). Confirm at execution which
  of the three needs the include line; that is the only permitted test edit.

## Exit criteria (machine-checked)

- Full build passes; test suite / baseline outputs identical to Stage-0 (B1, B2, B3),
  both modes. Because the drift branch determines the deterministic
  HMC log (B3 bitwise), any non-identical output means a case body was not moved
  byte-for-byte - revert.
- Public-symbol nm diff vs B0: the four new `robo::joint*`/`advanceQuatExp` symbols
  are the only additions, all engine-internal (no Python-surface delta); rationale
  recorded in the ticket per section 3.
- The byte-identical case-body equivalence check (moved vs origin, four sites) passes.
- include-cycle script clean; `JointKernels.hpp` gains no upward edge (already depends
  only on `RobotModel`/`robot_math`).
- clang-tidy / clang-format / IWYU clean on `JointKernels.{hpp,cpp}`,
  `RobotIntegrator.hpp`, and the kinematics TU.
- `JointKernels.hpp` <= 300 LOC; `JointKernels.cpp` <= 600 LOC.
- Only the predicted include fix (at most one line, in at most one of the three named
  test files) is applied; no other test changed.
