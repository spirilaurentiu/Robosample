# SPLIT-W5: extract the FixmanCorrection (metric + orientation Jacobian) from World.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/World.cpp` (2619 LOC) is a god-object mixing ~12 responsibilities
(`ARCHITECTURE.md section 7`). This ticket extracts the **Fixman / coordinate-Jacobian
corrections** for torsional worlds: the compensating potential
`U_F = 1/2 RT (ln|M_tree| - ln det(G M^-^1 G^T) - ln|M_3N|)` (Spiridon & Minh 2017)
and the external-rotation Jacobian `-1/2 RT Sigma_b ln sin^2(gamma2_b)` (default OFF; wrong
for the engine's unit-quaternion roots). These are pure functions of
`(RobotModel, RobotState, ConstraintSet)` plus the constant `lnDetMCartesian_`.

High-independence extraction (pure functions, no per-step mutation beyond the
realize calls they already make). Sequenced early among the World splits
(`README.md section 4`). Invariants below come from `ARCHITECTURE.md` and are restated
in full.

## Moves (exactly what)
- `double World::calcFixman()` (`src/World.cpp` lines 1125-1155) ->
  `src/world/FixmanCorrection.cpp` [private World member definition relocated;
  declaration stays `World.hpp:652`].
- `double World::calcLogSineSqrGamma2() const` (lines 1218-1259) ->
  `src/world/FixmanCorrection.cpp` [private World member definition relocated;
  declaration stays `World.hpp:653`].
- File-scope helper `static bool freeRootAbsRotation(const RobotModel&, const
  RobotState&, int, Rotation&)` (lines 1177-1216, with its 1157-1176 doc banner)
  -> anonymous namespace of `src/world/FixmanCorrection.cpp` [internal linkage;
  used only by `calcLogSineSqrGamma2`, call site line 1250].
- File-local anonymous-namespace helper `Real safeLogSineSqr(Real)` (lines
  187-195) -> see the shared-converter note below [used only by
  `calcLogSineSqrGamma2`, line 1256].
- File-local anonymous-namespace helper `void rotationToQuaternion(const
  Rotation&, Real&, Real&, Real&, Real&)` (lines 156-184) -> see the shared-
  converter note below [used by `calcLogSineSqrGamma2` line 1254 **and** by W4's
  `ncmcApplyTroughTeleport`].

The `lnDetMCartesian_` member (`World.hpp:715`) stays declared on `World` (it is
written by `buildModel`/W1 at line 950 and read here); the inline
`currentConstraintLogDet()` accessor (`World.hpp:364-368`) is unchanged.

### Shared-converter note (cross-cutting; flag for human review)
`rotationToQuaternion` (W5, W4), `quatToRotation` (W2, W4) and `safeLogSineSqr`
(W5) are file-local statics shared across three World extractions. Duplicating
them per-`.cpp` is not pure motion (it changes the symbol count). The project
already documents a single-source target: `tests/EngineHelpers.hpp`'s banner
names `include/engine_helpers.hpp` as the intended home for exactly
`rotationToQuaternion` / `quatToRotation` / `safeLogSineSqr` (hoisted "verbatim so
the engine TUs and the unit-test TU bind the same symbol"). That hoist is **not
yet done** (no `include/engine_helpers.hpp` exists; `World.cpp` still owns the
copies). Because W5 is the first World ticket to need `rotationToQuaternion` +
`safeLogSineSqr`, it performs that hoist: move the three converters verbatim from
`World.cpp`'s anonymous namespace into `include/engine_helpers.hpp` (namespace
`EngineHelpers`, matching the existing `advanceQuatExp` home), and have
`FixmanCorrection.cpp` `#include` it. `World.cpp` (still owning `quatToRotation`
until W2/W4) also switches its remaining call sites to the hoisted symbol. W2 and
W4 then depend on this header existing. **This coupling is a sequencing
dependency across W2/W4/W5 and SHOULD be confirmed with a human before
execution** (alternative: a `world/detail/` header instead of the documented
`engine_helpers.hpp`).

## Public API after this ticket
No new `world/` public header for the Fixman methods (they stay private `World`
members). `include/engine_helpers.hpp` (created/extended here) exports
`EngineHelpers::rotationToQuaternion`, `EngineHelpers::quatToRotation`,
`EngineHelpers::safeLogSineSqr` (verbatim moves). `src/world/FixmanCorrection.cpp`
is a `World`-class TU (`#include "World.hpp"`) defining the two relocated
members; IWYU: `<cmath>`, `<algorithm>` (`std::clamp`, line 1255),
`engine_helpers.hpp`, `Constraints.hpp` (already via `World.hpp`), `World.hpp`.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames (the converters keep
  their exact names in `EngineHelpers`).
- `src/World.cpp` loses lines 1125-1155, 1177-1259, and (156-195 / 198-217 as the
  shared hoist proceeds); it gains `#include "engine_helpers.hpp"` for any
  converter call sites it still owns. Nothing else changes.
- Files land at `MODULES.md section 1` paths (`src/world/FixmanCorrection.cpp`;
  `include/engine_helpers.hpp` is a flat `include/` header, already covered by the
  `include/*.h*` glob). Recursive `src/` glob per SPLIT-W1.
- Invariants that SHALL remain true:
  - **INV-5 Mass metric.** HMC momentum seeding and kinetic energy use the same
    mass-matrix operators (`multiplyBySqrtM` / `multiplyByMInv` /
    `calcKineticEnergy`); `calcLogDetM` is the Fixman kinetic term. `calcFixman`
    consumes `RobotEngine::calcLogDetM` (line 1138) as `ln|M_tree|` after
    `realizeArticulatedBodyInertias` (line 1132) - this ordering and operator
    identity SHALL be preserved verbatim.
  - **INV-6 Constraint consistency.** The `G` (constraint Jacobian) assembly is
    identical across SHAKE, RATTLE, and the loop-closure Fixman log-det, so the
    correction matches the projection. `calcConstraintLogDet` returns 0 for
    acyclic molecules. `calcFixman` subtracts
    `constraints_.calcConstraintLogDet(model_, state_)` (line 1150) as `ln det(G
    M^-^1 G^T)`; the relocated body SHALL keep this exact term so acyclic molecules
    still yield a guaranteed no-op.
  - **INV-9 Quaternion double cover.** Free/quaternion joints normalize
    quaternions each step; reversibility checks account for the double cover.
    `calcLogSineSqrGamma2` reads the body's TRUE orientation via
    `freeRootAbsRotation` (NOT `X_GB[b].R()`) precisely because the per-block
    frame reset zeroes q; the relocated helper + `rotationToQuaternion` +
    `safeLogSineSqr` chain (lines 1250-1256) SHALL move verbatim so the
    pole-flooring artifact stays fixed and the term stays gated OFF for quaternion
    roots (`sampler_.useOrientationJacobian`).
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None from the `World` side: no test `#include`s `World.cpp`; no `friend`s. The
  Fixman methods are private and reached only through `generateSample` /
  `reinitialize` (W9), unchanged here. `TestTwoRobotContact.cpp` and
  `TestEquipartition.cpp` **re-derive** `calcFixman`'s formula (comments at
  `TestTwoRobotContact.cpp:283`, `TestEquipartition.cpp:217`) rather than calling
  the private method, so they are untouched.
- `tests/EngineHelpers.hpp` already carries its own copies of `rotationToQuaternion`
  / `quatToRotation` / `safeLogSineSqr` under `namespace EngineHelpers`. If the
  shared hoist lands `include/engine_helpers.hpp` with the **same** namespace +
  signatures the test header already declares, an ODR clash is possible where a
  test TU includes both. The one allowed Phase-A test edit is therefore: in the
  test files that pull `tests/EngineHelpers.hpp` for these three converters
  (`TestQuaternion.cpp`, `TestEnsembleOrientation.cpp`, `TestTransfer.cpp`,
  `TestGeometry.cpp`, `TeleportMove.hpp`, `HmcDriver.hpp`, `TestRoboticsOracleMolecule.cpp`),
  retarget the `#include "EngineHelpers.hpp"` to `#include "engine_helpers.hpp"`
  as the tests' banner already anticipates ("now deleted, both sides #include
  this"). This is the only permitted test change and it is an include-path
  retarget, not a logic edit. **If the human review chooses a `world/detail/`
  home instead of `engine_helpers.hpp`, this test edit is not needed and W5 makes
  zero test changes.**

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0.
- `nm` diff on public symbols: empty (`World::calcFixman` /
  `calcLogSineSqrGamma2` are private; the `EngineHelpers::*` converters are
  `inline` in a header - no new exported non-inline symbol).
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- `FixmanCorrection.cpp` <= 600 LOC (~120 expected).
- Comment-stripped before/after diff of every touched file is empty apart from the
  moved blocks and the added `#include`s.
