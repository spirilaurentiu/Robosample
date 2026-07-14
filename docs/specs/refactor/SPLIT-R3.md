# SPLIT-R3: cohesive translation-unit split of RobotEngine.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)

After SPLIT-R1 (linear algebra removed) and SPLIT-R2 (NaN scanner removed),
`src/RobotEngine.cpp` still holds the whole `RobotEngine` class body - kinematics,
articulated-body dynamics, mass-matrix operators, and mobilizer reaction forces -
in one ~1000-LOC translation unit. This ticket splits that ONE `.cpp` into four
cohesive `.cpp` files that compile into the SAME `RobotEngine` class against the
SAME `include/RobotEngine.hpp`. No method signature, no member, and no header
changes: this is a file partition only.

Target files (per `MODULES.md`, layer `dynamics/`):
`RobotEngine_kinematics.cpp`, `RobotEngine_dynamics.cpp`, `RobotEngine_massops.cpp`,
`RobotEngine_reaction.cpp`. The original `RobotEngine.cpp` is retired (its
translation units are these four); `RobotEngine.hpp` is unchanged.

Invariants below come from `ARCHITECTURE.md section 5`, restated in full. A TU boundary
SHALL NOT cross any of them.

- **INV-1 Force->wrench convention.** Per-body spatial force in
  `RobotState::bodyForceG` is `(torque about the body origin, net force)` in the
  Ground frame. Both the host reduction (`ForceBridge.hpp`) and the CUDA
  `reduceForces` kernel SHALL produce this identically. (Relevant to the reaction TU,
  which consumes `bodyForceG`.)
- **INV-4 Realization order.** `RobotState` caches are valid only in stage order
  position -> velocity -> articulated-body inertias -> udot. The contract is enforced
  by convention and doc comments, not by types. The split preserves every method's
  precondition text; no method is reordered relative to its callers.
- **INV-5 Mass metric.** HMC momentum seeding and kinetic energy use the same mass
  matrix operators (`multiplyBySqrtM` / `multiplyByMInv` / `calcKineticEnergy`);
  `calcLogDetM` is the Fixman kinetic term. These land together in the massops TU so
  the shared metric stays one unit.

## Moves (exactly what)

Each method moves verbatim (definition body byte-for-byte) to its TU. Line ranges
are Stage-0 positions in the current `RobotEngine.cpp`; R1 and R2 run before R3 in
Lane E and remove ~600 lines above them, so re-anchor by symbol at execution.

`RobotEngine_kinematics.cpp`:
- `realizePosition` (lines 526-613)
- `realizeVelocity` (lines 615-728)
- `calcQDot` (lines 730-788)
- `calcQDotDot` (lines 1403-1453)
- `fillAtomPositionsFromBodies` (lines 1346-1375)
- `normalizeQuaternions` (lines 1380-1398)

`RobotEngine_dynamics.cpp`:
- `factorizeArticulatedInertias` (lines 798-988; includes the STEP-3 throwing gate
  867-901 and the `#if ROBO_DEBUG` diagnostic block 903-946)
- `seedArticulatedCentrifugal` (lines 995-1009)
- `realizeArticulatedBodyInertias` (lines 1015-1018)
- `calcUDot` (lines 1023-1082)

`RobotEngine_massops.cpp`:
- `multiplyByMInv` (lines 1087-1135)
- `multiplyBySqrtMInv` (lines 1140-1167)
- `multiplyBySqrtM` (lines 1183-1217)
- `calcLogDetM` (lines 1219-1243)
- `calcKineticEnergy` (lines 1249-1258)

`RobotEngine_reaction.cpp`:
- `calcMobilizerReactionForces` (lines 1282-1330)
- `findMobilizerReactionOnBodyAtMInGround` (lines 1332-1340)

Shared file-local helper - `spatialDot` (lines 320-322). It is used by the dynamics,
massops, and reaction TUs. Duplicating it into each TU is not pure motion. Move it
once to a NON-PUBLIC translation-unit-shared header `src/RobotEngine_internal.hpp`
(`namespace { }` is wrong across TUs - use an `inline` free function in
`namespace robo::detail`, or `static inline` in the internal header). Each of the
four TUs `#include "RobotEngine_internal.hpp"`. This header is engine-private (lives
under `src/`, not `include/`), holds only `spatialDot`, and depends only on
`robot_math.hpp`.

Per-TU includes (IWYU): every TU includes `RobotEngine.hpp`, `robot_math.hpp`,
`RobotEngine_internal.hpp`. Additionally:
- kinematics: `JointKernels.hpp` (jointX_FM/jointH_FM/jointHDot_FM), `RobotModel.hpp`.
- dynamics: `math/hinge_linalg.hpp` (invertDense/eigDecompAndTol) + the
  using-declarations from R1, `math/robo_debug.hpp` (the ROBO_DEBUG block),
  `<stdexcept>`/`<string>` (the throwing gate).
- massops: `math/hinge_linalg.hpp` (symSqrt/symSqrtInv/pseudoLogDet) + R1
  using-declarations, `<vector>` (`multiplyBySqrtM` scratch).
- reaction: `<vector>` (scratch `reacBo`/`atM`).

The R1 using-declarations (`using robo::detail::invertDense;` etc.) are replicated
in exactly the TUs that call each kernel (dynamics for invertDense/eigDecompAndTol;
massops for symSqrt/symSqrtInv/pseudoLogDet). This is not a rename - the call
expressions are unchanged.

Not in scope: the integrator templates already live in `RobotIntegrator.hpp`
(SPLIT-I1); `calcQDot`/`calcQDotDot` switch consolidation is SPLIT-R4. R3 moves
these methods as-is.

## Public API after this ticket

`include/RobotEngine.hpp` is UNCHANGED - same class, same 20 method declarations.
The four `.cpp` files define disjoint subsets of those methods. No symbol is added,
removed, or renamed. `RobotEngine_internal.hpp` is engine-private (not installed,
not on the public include path) and exports only `robo::detail::spatialDot`.

## Constraints

- Pure code motion. Zero logic changes. Zero symbol renames. Method bodies move
  byte-for-byte.
- `include/RobotEngine.hpp` does not change.
- Build system: replace the single `RobotEngine.cpp` entry in the engine object
  library target with the four new `.cpp` files. `spatialDot`'s single definition
  (internal header) prevents a duplicate-symbol link error.
- Invariants that SHALL remain true: **INV-1**, **INV-4**, **INV-5** (restated above).
  No method changes stage order or force convention; the split is orthogonal to all
  three.
- One commit for the move; formatting fixes separate.

## Predicted test breakage (Phase-A include fixes only)

**None.** Tests call `RobotEngine::` methods through `RobotEngine.hpp`, whose
declarations are unchanged; the linker binds each call to whichever of the four TUs
defines it. No test includes `RobotEngine.cpp` or `RobotEngine_internal.hpp`. Tests
that name `calcQDot`/`calcQDotDot` (`TestIntegrator.cpp`, `TestMobilizer.cpp`,
`TestJointKernels.cpp`, `HmcDriver.hpp`, and others) resolve the same declarations.
No include fix, no assertion-count change.

## Exit criteria (machine-checked)

- Full build passes; test suite / baseline outputs identical to Stage-0 (B1, B2, B3),
  both modes.
- Public-symbol nm diff vs B0: empty (same class, same defined methods - only their
  object-file location moved).
- include-cycle script clean; no upward layer edge introduced; `RobotEngine_internal.hpp`
  depends only on `robot_math.hpp`.
- clang-tidy / clang-format / IWYU clean on all four TUs and the internal header.
- Comment-stripped diff: the union of the four new `.cpp` files plus the internal
  header, compared against the pre-split `RobotEngine.cpp`, is empty apart from
  per-TU `#include` lines and the relocated `spatialDot`.
- Each new `.cpp` <= 600 LOC (dynamics is largest, ~330; well under cap); the
  internal header <= 300 LOC.
- No test file changed.
