# SPLIT-W1: extract the ModelBuilder (kinematic-tree construction) from World.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/World.cpp` (2619 LOC) is a god-object mixing ~12 responsibilities
(`ARCHITECTURE.md section 7`). This ticket extracts the **build-time** one: turning the
shared `SystemTopology` + this world's per-bond mobilities + this world's root
mobilities into an immutable `RobotModel` (body forest by union-find + BFS,
topological relabel, q/u index tables, quaternion-slot registry, z-matrix/BAT
rows, body-atom ranges, frame graph, mass-property init, loop-closure
constraints), plus the two in-place rebuild entry points that re-run it on a
root-mobility change.

This concept is the highest-independence World extraction: it runs once at
construction (and on `setRootMobility(ies)`), before any sampling, and touches no
per-step cache. It is sequenced **first** among the World splits
(`README.md section 4`). Invariants below come from `ARCHITECTURE.md` and are restated
in full.

## Moves (exactly what)
- `void World::buildModel(const SystemTopology&, const Selection&, const std::vector<JointType>&)`
  (`src/World.cpp` lines 601-968) -> `src/world/ModelBuilder.cpp` as a relocated
  `World::buildModel` definition [public API unchanged - declaration stays in
  `World.hpp:232-233`].
- `void World::setRootMobility(int, JointType)` (lines 983-994) ->
  `src/world/ModelBuilder.cpp` [public World member definition relocated;
  declaration stays `World.hpp:248`].
- `void World::setRootMobilities(const std::vector<JointType>&)` (lines 996-1007)
  -> `src/world/ModelBuilder.cpp` [public World member definition relocated;
  declaration stays `World.hpp:249`].
- File-local anonymous-namespace helper `struct DSU` (lines 79-91) ->
  anonymous namespace of `src/world/ModelBuilder.cpp` [internal linkage
  preserved]. Only `buildModel` uses it (call site line 635).
- File-local anonymous-namespace helper `bool isFlexible(JointType)` (lines
  93-95) -> anonymous namespace of `src/world/ModelBuilder.cpp` [internal linkage
  preserved]. Used only by `buildModel` (lines 641, 674, 853).

NOTE: `buildModel` calls `World::recomputeGeometry` (line 940) and
`RobotEngine::*`/`bridge_.markStationsDirty`. `recomputeGeometry` is extracted
later by **W6**; at W1 execution time it is still a `World` member declared in
`World.hpp`, so the relocated `buildModel` calls it unchanged. No W1<->W6 ordering
constraint results.

## Public API after this ticket
No new public header. `World`'s public surface is unchanged: `buildModel`,
`setRootMobility`, `setRootMobilities` keep their existing declarations and
signatures in `World.hpp`. `src/world/ModelBuilder.cpp` is a `World`-class
translation unit: it `#include "World.hpp"` and defines the three relocated
members plus the two file-local helpers in its own anonymous namespace. This is
the cohesive-TU-split pattern `MODULES.md` endorses for `RobotEngine` (R3).
Exception to the template's "target.h exports" clause is claimed: the concept
surface remains `World`'s existing public methods; a header exporting `DSU`/
`isFlexible` is undesirable (they are build-internal and SHALL keep internal
linkage). IWYU for the new `.cpp`: `<algorithm>` (`std::stable_sort`, line 757),
`<numeric>` (`std::iota` in `DSU`), `<queue>` (`std::queue`, line 714), `<set>`
is not needed here, `<stdexcept>` (`std::invalid_argument`, line 699), `<string>`
(`std::to_string`), plus `World.hpp`. `#pragma once` N/A (`.cpp`).

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/World.cpp` loses lines 79-95 (the two helpers) and 601-968 / 983-1007 (the
  three member bodies), and its includes for `<queue>` and `<numeric>` move with
  the code if no remaining `World.cpp` symbol needs them (verify: `<numeric>` is
  also used by `configureNcmc` `std::iota` at line 2167 - keep it in `World.cpp`;
  `<queue>` is used only by `buildModel` - it moves). Nothing else in `World.cpp`
  changes.
- Files land at their `MODULES.md section 1` target paths (`src/world/ModelBuilder.cpp`;
  no header). Because the `CMakeLists.txt` source glob
  `file(GLOB ROBOSAMPLE_CXX_SOURCE_FILES .../src/*.cpp)` is flat (non-recursive),
  this ticket switches it to `file(GLOB_RECURSE ...)` (and the matching
  `include/*.h*` glob) so `src/world/*.cpp` compiles. Cite the glob by symbol, not
  line - the line number drifts. This is build-system only; observable behavior is
  preserved. If an earlier (R/O-series) ticket already made the glob recursive, this
  is a no-op. (CX-2, ACCEPTED.)
- Invariants that SHALL remain true:
  - **Model-build order (`ARCHITECTURE.md section 3` Construction).** `buildModel` turns
    bonds + per-joint mobilities into a `RobotModel`: body forest via union-find
    + BFS, topological (parent<child) relabel, q/u index tables, static joint
    frames, mass properties, loop-closure constraints - in exactly this order.
    The relocated body preserves the sequence line-for-line (DSU join of welded
    bonds -> component->body map -> joint-edge BFS forest -> stable-sort relabel ->
    children ranges -> q/u/quaternion index tables -> z-matrix rows -> body-atom
    ranges -> frame graph -> `allocateFull` -> reference-geometry `recomputeGeometry`
    -> q/u reset -> `lnDetMCartesian_` -> loop-closure `DistanceConstraint`s).
  - **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
    between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
    currency. `buildModel` establishes topology only; it introduces no
    cross-sweep state.
  - **INV-6 Constraint consistency.** `calcConstraintLogDet` returns 0 for acyclic
    molecules. The relocated loop-closure loop (lines 957-967) SHALL keep pushing a
    `DistanceConstraint` for exactly the ring-closing bonds whose endpoints landed
    in different bodies - the same `G`-assembly source the Fixman log-det reads.
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None. No test `#include`s `World.cpp`; there are no `friend` declarations on
  `World`. `TestBatScalingJacobian.cpp`, `TestBatAnchorInvolution.cpp`,
  `TestTwoRobotContact.cpp`, and `TestRoboticsOracleMolecule.cpp` call
  `world.buildModel(...)` through the public declaration in `World.hpp`, which is
  unchanged; they keep `#include "World.hpp"` verbatim. `DSU`/`isFlexible` keep
  internal linkage and were never test-visible (`TestRoboticsOracleMolecule.cpp`
  carries its own local union-find, per its comment at lines 367-369). No
  include-path fix is required by any test.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0.
- `nm` diff on public symbols: empty (`World::buildModel`,
  `World::setRootMobility`, `World::setRootMobilities` keep the same mangled
  names; `DSU`/`isFlexible` are internal and excluded).
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- `src/world/ModelBuilder.cpp` <= 600 LOC (~400 expected). No header added.
- Comment-stripped before/after diff of `World.cpp` + `ModelBuilder.cpp` is empty
  apart from the moved blocks appearing at their new location.
