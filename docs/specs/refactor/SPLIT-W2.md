# SPLIT-W2: extract the DockingMove (ligand reposition + sphere sampling) from World.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/World.cpp` (2619 LOC) is a god-object mixing ~12 responsibilities
(`ARCHITECTURE.md section 7`). This ticket extracts the **docking** concern: the
auto-sized binding-sphere geometry, the rigid-kick proposal (uniform position in
the sphere + uniform SO(3) reorientation of each ligand), and the clash-free
initial-placement search. The move itself (accept/reject) lives in HMC (W9); this
extraction is the proposal machinery plus the docking configuration entry point.

High-independence: every symbol here is gated by `docking_` and touches only the
docking bookkeeping (`ligandGroups_`, `siteAtoms_`, `sampler_.sphereFactor`,
`dockingStuckCount_`) plus the RNG. Sequenced after W1/W3/W5 (`README.md section 4`).
Invariants below come from `ARCHITECTURE.md` and are restated in full.

## Moves (exactly what)
- `void World::configureDocking(std::vector<std::vector<int>>, std::vector<int>)`
  (`src/World.cpp` lines 294-299) -> `src/world/sampler/DockingMove.cpp` [public
  World member definition relocated; declaration stays `World.hpp:285`].
- `robo::Vec3 World::atomSetCentroid(const std::vector<int>&) const` (1282-1292)
  -> `src/world/sampler/DockingMove.cpp` [private World member; decl `World.hpp:665`].
- `robo::Vec3 World::atomSetMassCenter(const std::vector<int>&) const`
  (1294-1304) -> `DockingMove.cpp` [private; decl `World.hpp:666`].
- `double World::atomSetRadius(const std::vector<int>&, const robo::Vec3&) const`
  (1306-1316) -> `DockingMove.cpp` [private; decl `World.hpp:667`].
- `double World::groupSphereRadius(int) const` (1318-1322) -> `DockingMove.cpp`
  [private; decl `World.hpp:662`].
- `robo::Vec3 World::sampleUniformInSphere(double)` (1324-1329) ->
  `DockingMove.cpp` [private; decl `World.hpp:663`].
- `robo::Rotation World::sampleUniformRotation()` (1331-1340) -> `DockingMove.cpp`
  [private; decl `World.hpp:664`].
- `bool World::repositionLigands(bool)` (1342-1458) -> `DockingMove.cpp` [private;
  decl `World.hpp:661`].
- `int World::findGoodStartingPose()` (1477-1554) -> `DockingMove.cpp` [public
  World member definition relocated; declaration stays `World.hpp:312`].
- File-local anonymous-namespace helper `bool dockDebugEnabled()` (45-51) ->
  anonymous namespace of `DockingMove.cpp` [internal linkage]. NOTE: it appears
  **unused** in `World.cpp` (defined line 45, no call site found by grep) - move
  it with the docking concept regardless; flag the dead-code observation to the
  reviewer, do not delete it in a code-motion ticket.

`sampleUniformRotation` (line 1339) and (via W4 later) `ncmcApplyTroughTeleport`
use `quatToRotation`. Per the shared-converter note in **SPLIT-W5**, that
converter is hoisted to `include/engine_helpers.hpp`; `DockingMove.cpp`
`#include`s it. **W2 therefore depends on W5 having performed that hoist** (or on
the human-approved `world/detail/` alternative).

## Public API after this ticket
No new public header. `World`'s public docking surface (`configureDocking`,
`findGoodStartingPose`, `isDocking`) is unchanged. `src/world/sampler/
DockingMove.cpp` is a `World`-class TU (`#include "World.hpp"`) defining the ten
relocated members. IWYU: `<cmath>` (`std::sqrt`/`std::acos`/`std::cbrt`/`std::cos`/
`std::sin`, `std::cbrt` line 1327), `<cstdio>` (`std::fprintf`/`std::snprintf`),
`<cstdlib>` (`std::getenv` in `dockDebugEnabled`), `<stdexcept>`
(`std::runtime_error`, line 1553), `<vector>`, `engine_helpers.hpp`
(`quatToRotation`), `World.hpp`.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/World.cpp` loses lines 45-51, 294-299, 1282-1340, 1342-1458, 1477-1554.
  The `<cstdlib>` include stays in `World.cpp` only if another remaining symbol
  needs `std::getenv` (`nmaDebugEnabled`/`ncmcDebugEnabled` at lines 46-47/72-73
  still use it - keep `<cstdlib>` in `World.cpp`). Nothing else changes.
- Files land at `MODULES.md section 1` path (`src/world/sampler/DockingMove.cpp`);
  recursive `src/` glob per SPLIT-W1.
- Invariants that SHALL remain true:
  - **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
    between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
    currency. `repositionLigands` and `findGoodStartingPose` operate on
    `state_.atomPosG()` and re-fit through `setAtomsLocationsInGround` (lines 1455,
    1539) - the same Cartesian-coordinate currency; no cross-sweep state is
    introduced. `dockingStuckCount_` is a within-run counter, not per-replica
    carried state, and stays declared on `World`.
  - **INV-9 Quaternion double cover.** `sampleUniformRotation` uses Shoemake's
    uniform-quaternion -> SO(3) map (uniform on SO(3), lines 1331-1339); the
    reorientation SHALL remain Haar-uniform so the rigid kick is a symmetric
    Cartesian proposal (the docking move's reversibility rationale,
    `World.hpp:41-45` / `World.cpp:1273-1280`). Preserve the converter verbatim.
  - **Model-build order (`ARCHITECTURE.md section 3`).** `repositionLigands` calls
    `setAtomsLocationsInGround` (W6) to rebuild q/frames from the proposed pose;
    at W2 execution time that method is still a `World` member (W6 runs later), so
    the relocated body calls it unchanged - no W2<->W6 ordering constraint.
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None. No test `#include`s `World.cpp`; no `friend`s. The docking symbols are not
  named by any C++ test (grep of `tests/`: no `configureDocking` /
  `findGoodStartingPose` / `repositionLigands` call sites). `configureDocking` and
  `findGoodStartingPose` stay declared in the unchanged `World.hpp`. No test edit
  is required. (The `engine_helpers.hpp` include retarget, if any, is booked
  against SPLIT-W5, not here.)

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0.
- `nm` diff on public symbols: empty (`configureDocking`, `findGoodStartingPose`
  keep their mangled names; the eight private helpers and `dockDebugEnabled` are
  non-public).
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- `DockingMove.cpp` <= 600 LOC (~300 expected).
- Comment-stripped before/after diff of `World.cpp` + `DockingMove.cpp` is empty
  apart from the moved blocks at their new location.
