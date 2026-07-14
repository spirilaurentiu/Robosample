# SPLIT-W3: extract the ReactionReporter (per-body spatial-force snapshot) from World.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/World.cpp` (2619 LOC) is a god-object mixing ~12 responsibilities
(`ARCHITECTURE.md section 7`). This ticket extracts the **reaction-force reporter**
(reaction-force-monitoring campaign; its spec is not yet in the repo, so the
behavior below is recovered from source): an opt-in, orthogonal telemetry
concern that snapshots, per selected body, the sum of the OpenMM net applied
force (`bodyForceG`) and/or the static (u=0) mobilizer reaction, into a per-frame
CSV row buffer. A world that never enables it allocates nothing and does zero
extra work.

This is a high-independence extraction: reporter state is a self-contained set of
five members plus one value struct, written only through `setReactionReporter` /
`enableReactionReporter` and read only through `captureReactionSnapshot` /
`reactionSamples`. It is sequenced early among the World splits (`README.md section 4`).
Invariants below come from `ARCHITECTURE.md` and are restated in full.

## Moves (exactly what)
- `struct ReactionSample` (`include/World.hpp` lines 97-102, with its 82-96 doc
  banner) -> `include/robo/world/ReactionReporter.hpp` [public, top-level struct,
  same name/namespace]. `World.hpp` `#include`s the new header so the top-level
  name `ReactionSample` stays visible to every current includer.
- `void World::setReactionReporter(std::vector<int>)` (`src/World.cpp` lines
  403-413) -> `src/world/ReactionReporter.cpp` [public World member definition
  relocated; declaration stays `World.hpp:481`].
- `void World::enableReactionReporter(bool, bool, bool)` (lines 415-479) ->
  `src/world/ReactionReporter.cpp` [public World member definition relocated;
  declaration stays `World.hpp:520-522`].
- `void World::captureReactionSnapshot()` (lines 481-573) ->
  `src/world/ReactionReporter.cpp` [public World member definition relocated;
  declaration stays `World.hpp:547`].

The reporter member fields stay declared on `World` (private, `World.hpp`
695-699: `reactionReporter_`, `reactionIncludeOpenmm_`, `reactionIncludeReaction_`,
`reactionInterestingBodies_`, `reactionSamples_`) and the three inline accessors
(`isReactionReporter`, `reactionInterestingBodies`, `reactionSamples`,
`World.hpp:523-554`) are unchanged. Only the three out-of-line method bodies and
the value struct move.

## Public API after this ticket
`include/robo/world/ReactionReporter.hpp` exports: `struct ReactionSample { int
bodyIdx; int atomIdx; robo::Vec3 force; robo::Vec3 torque; };` and nothing else.
Header self-contained (IWYU): `#include "robot_math.hpp"` (for `robo::Vec3`);
`#pragma once`. `World`'s public surface (`setReactionReporter`,
`enableReactionReporter`, `captureReactionSnapshot`, `isReactionReporter`,
`reactionInterestingBodies`, `reactionSamples`) is unchanged. `src/world/
ReactionReporter.cpp` is a `World`-class TU (`#include "World.hpp"`) defining the
three relocated members; its IWYU set: `<set>` (`std::set`, line 434),
`<vector>`, `<stdexcept>` (`std::runtime_error`), `<cstddef>`, and `World.hpp`.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/World.cpp` loses lines 403-573; `include/World.hpp` loses the
  `ReactionSample` struct + banner (lines 82-102) and gains `#include
  "robo/world/ReactionReporter.hpp"`. Nothing else changes.
- New files land at their `MODULES.md section 1` target paths (`include/robo/world/
  ReactionReporter.hpp`, `src/world/ReactionReporter.cpp`); the recursive glob
  from SPLIT-W1 (`CMakeLists.txt`) covers them. If W1 has not yet run,
  switch lines 142/144 to `GLOB_RECURSE` here.
- Invariants that SHALL remain true:
  - **INV-1 Force->wrench convention.** Per-body spatial force in
    `RobotState::bodyForceG` is `(torque about the body origin, net force)` in the
    Ground frame. Both the host reduction (`ForceBridge.hpp`) and the CUDA
    `reduceForces` kernel SHALL produce this identically. `captureReactionSnapshot`
    reads `bodyForceG[b].linear`/`.angular` (lines 564-565) as force/torque about
    Bo in Ground; the relocated body SHALL preserve that read verbatim.
  - **INV-2 Virtual sites.** Virtual-site forces are already projected onto parent
    atoms by OpenMM; the reduction skips massless particles to avoid double
    counting. The reporter consumes the already-reduced `bodyForceG` and the
    representative-atom pick uses `bodyAtomsBeg` (line 555); no massless handling
    changes.
  - **INV-4 Realization order.** `RobotState` caches are valid only in stage order
    position -> velocity -> articulated-body inertias -> udot, enforced by
    convention. `captureReactionSnapshot`'s optional u=0 reaction path
    (`realizeArticulatedBodyInertias` -> save u -> zero -> `realizeVelocity` ->
    `calcUDot` -> `calcMobilizerReactionForces` -> restore u -> `realizeVelocity` ->
    `calcUDot`, lines 527-544) SHALL keep this exact sequence, including the
    save/zero/restore that makes the snapshot non-perturbing to the next
    `generateSample` (the no-perturbation property, spec Sec.3/6).
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None. No test `#include`s `World.cpp`; no `friend`s exist. `TestReactionForces.cpp`
  does not `#include "World.hpp"` (it exercises the reporter through the pybind /
  higher-level surface, not the C++ `World` type). `ReactionSample` is not named
  by any C++ test; even so, `World.hpp` re-includes `ReactionReporter.hpp`, so the
  top-level name stays visible to every current includer. No include-path fix is
  required by any test.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0.
- `nm` diff on public symbols: empty. `ReactionSample`'s definition relocates but
  its ODR identity and layout are unchanged.
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- `ReactionReporter.hpp` <= 300 LOC (~30 expected); `ReactionReporter.cpp`
  <= 600 LOC (~180 expected).
- Comment-stripped before/after diff of `World.cpp`, `World.hpp`, and the two new
  files is empty apart from the moved blocks at their new location and the added
  `#include`.
