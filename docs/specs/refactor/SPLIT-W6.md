# SPLIT-W6: extract the GeometryFitter (Gibbs-block coordinate handoff) from World.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/World.cpp` (2619 LOC) is a god-object mixing ~12 responsibilities
(`ARCHITECTURE.md section 7`). This ticket extracts the **geometry re-fitting** concern:
`setAtomsLocationsInGround` (the Gibbs-block continuation entry point that copies
in per-atom Ground coordinates, rebuilds rigid-body frames, resets q, and retargets
loop-closure distances) and its worker `recomputeGeometry` (which rebuilds every
body's mass properties and the `X_PF`/`X_BM` joint transforms from the incoming
Cartesian geometry).

This is the mechanism behind INV-3: Cartesian coordinates are the sole inter-world
currency, and this pair is how one block hands the configuration to the next.
Sequenced after W1/W3/W5/W2 (`README.md section 4`). Invariants below come from
`ARCHITECTURE.md` and are restated in full.

## Moves (exactly what)
- `void World::setAtomsLocationsInGround(const std::vector<robo::Vec3>&)`
  (`src/World.cpp` lines 1035-1064, with its 1009-1034 doc banner) ->
  `src/world/GeometryFitter.cpp` [public World member definition relocated;
  declaration + the inline `getAtomsLocationsInGround` accessor stay
  `World.hpp:299-302`].
- `void World::recomputeGeometry(const robo::Vec3*)` (lines 1069-1120) ->
  `src/world/GeometryFitter.cpp` [private World member; decl `World.hpp:644`].
- File-local anonymous-namespace helper `void computeAllAtomFrames(const
  RobotModel::FrameGraph&, const Vec3*, Transform*, Transform*)` (lines 127-153)
  -> anonymous namespace of `src/world/GeometryFitter.cpp` [internal linkage; used
  only by `recomputeGeometry`, line 1070].
- File-local anonymous-namespace helper `void atomFrameFromGeometry(const Vec3*,
  int, int, int, int, Transform&, Transform*)` (lines 97-125) -> anonymous
  namespace of `src/world/GeometryFitter.cpp` [internal linkage; used only by
  `computeAllAtomFrames`, lines 135].

The `frameFlat_` / `xpcFlat_` scratch members (`World.hpp:740-741`) stay declared
on `World`; they are sized by `buildModel` (W1, lines 923-924) and consumed here.

## Public API after this ticket
No new public header. `World::setAtomsLocationsInGround` /
`getAtomsLocationsInGround` keep their existing public declarations in
`World.hpp`. `src/world/GeometryFitter.cpp` is a `World`-class TU (`#include
"World.hpp"`) defining the two relocated members and the two file-local helpers.
IWYU: `<algorithm>` (`std::min`, `std::fill`, `std::copy`), `<cmath>`
(`.norm()` uses; `robo::calcDihedralAngle`), `<vector>`, `robot_math.hpp` (via
`World.hpp`), `World.hpp`.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/World.cpp` loses lines 97-153 (the two helpers) and 1009-1120 (the two
  member bodies). Nothing else changes.
- Files land at `MODULES.md section 1` path (`src/world/GeometryFitter.cpp`); recursive
  `src/` glob per SPLIT-W1.
- Invariants that SHALL remain true:
  - **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
    between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
    currency. `setAtomsLocationsInGround` IS that handoff: it copies incoming
    Ground coordinates, rebuilds frames from the actual (not idealized) geometry,
    resets q=0 with identity quaternions, and retargets each loop-closure
    `restLength` to the carried-over distance (lines 1061-1063) - never to the
    force-field r0. The relocated body SHALL preserve this so Gibbs continuation
    and bit-for-bit reject-restore hold.
  - **INV-4 Realization order.** `RobotState` caches are valid only in stage order
    position -> velocity -> articulated-body inertias -> udot. `setAtomsLocationsInGround`
    ends by realizing position only (`RobotEngine::realizePosition`, line 1051)
    after the frame rebuild; the relocated body SHALL keep that single realize and
    not add velocity/inertia stages the callers do not expect.
  - **INV-1 Force->wrench convention.** `recomputeGeometry` refits `atomStation_B`
    and calls `bridge_.markStationsDirty()` (line 1119) so the fused CUDA path
    re-uploads stations before its next force reduction; preserving this call is
    required for the CUDA `reduceForces` kernel to keep producing the
    `(torque about Bo, net force)` wrench identically to the host path.
  - **INV-9 Quaternion double cover.** The q reset seeds identity quaternions at
    every `quaternionQStart` slot (lines 1046-1049); preserve verbatim so
    quaternion bodies re-initialize on the unit sphere.
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None. No test `#include`s `World.cpp`; no `friend`s. Tests that reference the
  geometry pipeline do so via **their own** hoisted copies of the frame helpers:
  `tests/EngineHelpers.hpp` and `RobotBuilders.hpp`/`TestGeometry.cpp` carry
  independent `atomFrameFromGeometry`/`computeAllAtomFrames`-equivalent code
  (`TestRoboticsOracleMolecule.cpp:91` calls out `atomFrameFromGeometry` as a
  `World.cpp`-private symbol it re-derives, not links). No test names
  `setAtomsLocationsInGround` as a private-helper reach-in; the four World.hpp
  includers use it only through its unchanged public declaration. No include-path
  fix is required by any test.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0.
- `nm` diff on public symbols: empty.
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- `GeometryFitter.cpp` <= 600 LOC (~200 expected).
- Comment-stripped before/after diff of `World.cpp` + `GeometryFitter.cpp` is
  empty apart from the moved blocks at their new location.
