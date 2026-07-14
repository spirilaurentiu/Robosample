# SPLIT-W4: extract the NcmcMove (alchemical decouple-move-recouple) from World.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/World.cpp` (2619 LOC) is a god-object mixing ~12 responsibilities
(`ARCHITECTURE.md section 7`). This ticket extracts the **per-molecule NCMC** move
(`docs/specs/ncmc-explicit-solvent/`): a lambda:1->0->1 alchemical
decouple-move-recouple proposal, its palindromic protocol schedule accessor, the
lambda=0 trough teleport, the Construction-II inner GHMC kernel, and the NCMC
configuration entry points. This is `MoveType::NcmcSwitch`, dispatched from
`generateSample` (W9) but otherwise a self-contained move regime.

Medium independence: it depends on the docking helpers (W2) for the teleport and
on the momentum/energy assembly (W9) at call time, but its own bodies form a
cohesive block at the tail of the file. Sequenced after W1/W3/W5/W2/W6/W7/W8,
immediately before W9 (`README.md section 4`). Invariants below come from
`ARCHITECTURE.md` and are restated in full.

## Moves (exactly what)
- `void World::configureNcmc(std::vector<int>, int, double)` (`src/World.cpp`
  lines 2150-2161) -> `src/world/sampler/NcmcMove.cpp` [**public** World member
  definition relocated (PyBind-bound, `PyBind11.cpp:81`); declaration stays
  `World.hpp:294`].
- `void World::configureNcmc(int, int, int, double)` overload (lines 2163-2170) ->
  `NcmcMove.cpp` [**public** World member definition relocated (PyBind-bound,
  `PyBind11.cpp:71`); decl `World.hpp:297`].
- `double World::protocolLambda(int) const` (lines 2172-2178) -> `NcmcMove.cpp`
  [private World member; decl `World.hpp:626`].
- `int World::ncmcTeleportRoot() const` (lines 2180-2201) -> `NcmcMove.cpp`
  [private World member; decl `World.hpp:627`].
- `void World::ncmcApplyTroughTeleport(int)` (lines 2203-2251) -> `NcmcMove.cpp`
  [private World member; decl `World.hpp:628`].
- `bool World::ncmcInnerGhmcStep(robo::Real, int, bool*)` (lines 2253-2393) ->
  `NcmcMove.cpp` [private World member; decl `World.hpp:643`].
- `bool World::ncmcMove()` (lines 2395-2619, end of file) -> `NcmcMove.cpp`
  [private World member; decl `World.hpp:625`].
- File-local anonymous-namespace helper `bool ncmcDebugEnabled()` (lines 64-77,
  with its `64-70` doc banner) -> anonymous namespace of
  `src/world/sampler/NcmcMove.cpp` [internal linkage; used **only** by
  `ncmcInnerGhmcStep` (call sites 2326, 2348) - no caller outside the moved set, so
  unlike W8's `nmaDebugEnabled` it needs no shared header and moves verbatim into
  `NcmcMove.cpp`'s own anonymous namespace].

The `generateSample` dispatch `if (sampler_.moveType == MoveType::NcmcSwitch)
return ncmcMove();` (`src/World.cpp:1560-1561`) **stays** in `World.cpp` (it is part
of `generateSample`, W9's concern) and calls the relocated `ncmcMove` through the
unchanged `World.hpp:625` declaration. The NCMC-related state stays declared on
`World`: `savedSolventPosG_` (`World.hpp:736`), `savedQ_` (line 737), `siteAtoms_`,
and the `sampler_.ncmc*` fields (the `RobotState` NCMC-solvent fields are the
concern-bleed of `OQ-6`; not relocated here). The `setNcmcTeleport` /
`setNcmcConstructionII` inline setters (`World.hpp:457-467`) are unchanged.

`ncmcApplyTroughTeleport` uses `quatToRotation` (line 2215) and
`rotationToQuaternion` (line 2224) - hoisted to `include/engine_helpers.hpp` by
**SPLIT-W5**; `NcmcMove.cpp` `#include`s it. It also calls the docking helpers
`atomSetCentroid`/`atomSetRadius`/`sampleUniformInSphere`/`sampleUniformRotation`
(lines 2216, 2219-2221), which **SPLIT-W2** relocates but keeps as `World`
members - so the relocated teleport calls them unchanged. **W4 depends on W2 and
W5 having run** (both precede it in the `README.md section 4` order).

## Public API after this ticket
No new public header. The two `configureNcmc` overloads keep their existing public
declarations in `World.hpp`. `src/world/sampler/NcmcMove.cpp` is a `World`-class
TU (`#include "World.hpp"`) defining the seven relocated members plus
`ncmcDebugEnabled`. IWYU: `<algorithm>` (`std::sort`/`std::unique`/`std::max`/
`std::copy`/`std::fill`), `<cmath>` (`std::isfinite`/`std::exp`), `<cstdio>`
(`std::fprintf`), `<cstdlib>` (`std::getenv`), `<limits>`
(`std::numeric_limits`, line 2338), `<numeric>` (`std::iota`, line 2167),
`<vector>`, `NCMCProtocol.hpp` (`robo::ncmc::protocolLambda`, line 2177; already
via `World.cpp`'s include), `engine_helpers.hpp`, `World.hpp`,
`RobotIntegrator.hpp`. The last is REQUIRED, not optional: `ncmcInnerGhmcStep`
(line 2320) and `ncmcMove` (line 2513) call `RobotEngine::stepTo`, a function
**template** whose definition lives in `RobotIntegrator.hpp` (`World.hpp` does not
include it). Without this include `NcmcMove.cpp` sees only the template
declaration and cannot instantiate `stepTo<Bridge>`, so it links only by
accidental COMDAT rescue from whichever sibling TU (`World.cpp` pre-W9, then
`HmcMove.cpp`) instantiates the same specialization - the exact fragile cross-TU
coupling the self-contained-TU rule forbids. IWYU/link SHALL be exercised on
`NcmcMove.cpp` as an ISOLATED TU, not relying on a sibling's instantiation.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/World.cpp` loses lines 64-77 (the helper + banner) and 2150-2620 (the file
  tail). The
  `#include "NCMCProtocol.hpp"` (`World.cpp:28`) moves to `NcmcMove.cpp` if no
  remaining `World.cpp` symbol needs it (verify: only `protocolLambda` uses
  `robo::ncmc::protocolLambda` - it moves). `<numeric>` stays in `World.cpp` only
  if another symbol needs it after W1 removed `buildModel` (`std::iota` here is the
  last user besides W1's already-moved `DSU`; move `<numeric>` with this code).
  Nothing else changes.
- Files land at `MODULES.md section 1` path (`src/world/sampler/NcmcMove.cpp`); recursive
  `src/` glob per SPLIT-W1.
- Invariants that SHALL remain true:
  - **INV-1 Force->wrench convention.** Per-body spatial force in
    `RobotState::bodyForceG` is `(torque about the body origin, net force)` in the
    Ground frame; the host reduction and the CUDA `reduceForces` kernel produce it
    identically. The NCMC reseed chains call `bridge_.evaluate` (the same
    reduction) but do not touch it - pure motion leaves the convention unchanged.
    Restated for parity with W9; no NCMC code alters the reduction.
  - **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
    between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
    currency. `ncmcMove` snapshots `savedQ_`/`savedSolventPosG_` at move start and
    restores them on reject (lines 2405-2416, 2606-2617); these are within-move
    rollback buffers, not cross-sweep carried state. On accept it emits the new pose
    through `RobotEngine::fillAtomPositionsFromBodies` (line 2599) - the same
    Cartesian-coordinate currency. Preserve the save/restore structure verbatim so
    no NCMC state leaks across sweeps.
  - **INV-4 Realization order.** `RobotState` caches are valid only in stage order
    position -> velocity -> articulated-body inertias -> udot. `ncmcApplyTroughTeleport`
    (lines 2243-2250), `ncmcInnerGhmcStep`'s per-substep reseed (lines 2313-2316)
    and its reject-flip reseed (lines 2380-2387) each re-run the exact
    realizePosition -> fillAtomPositions -> evaluate -> realizeVelocity ->
    realizeArticulatedBodyInertias -> calcUDot -> calcQDot -> calcQDotDot chain; this
    ordering is correctness-critical (the documented freeze root cause was a stale
    `a0`) and SHALL move verbatim.
  - **INV-5 Mass metric.** The move resamples momenta each block via `reinitialize`
    (W9) using `sqrt(M)`; NCMC does no explicit momentum flip precisely because
    the Maxwell-Boltzmann resample overwrites the sign (lines 2427-2433). The
    acceptance assembly `H_lambda` reuses `currentTotalEnergy`'s exact terms verbatim
    (F1/F2) - preserve so the inner GHMC stays pi_lambda-invariant.
  - **INV-6 Constraint consistency.** The teleport is applied only for an ACYCLIC
    Free-root region with `constraints_.empty()` (line 2465) - no
    constraint-manifold branch ambiguity to guard. Preserve that gate so the
    palindromic map stays its own reverse.
  - **INV-9 Quaternion double cover.** `ncmcApplyTroughTeleport` sets a Haar-random
    root orientation and CO-ROTATES the root angular+linear speeds by DeltaR = R_new
    R_old^T so KE is exactly preserved (M_ang(q) is orientation-dependent), reading
    and writing the quaternion via `quatToRotation`/`rotationToQuaternion` (lines
    2215-2240). Preserve verbatim, including the direct q/u write (no
    `setAtomsLocationsInGround`) so live momenta stay valid.
  - **NCMC acceptance separation (`docs/specs/ncmc-explicit-solvent/10-...`).**
    Construction I (endpoint-DeltaH) and Construction II (protocol-work) SHALL NOT
    mix; `ncmcMove`'s terminal branch (lines 2595-2597) selects
    `metropolis(0.0, work)` vs `metropolis(Hstart, Hend)` on
    `useMetropolizedInner`. Move verbatim.
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None. No test `#include`s `World.cpp`; no `friend`s. `TestNcmcExplicitSolvent.cpp`
  and `TestNcmcTeleport.cpp` do not `#include "World.hpp"`; they reproduce
  `ncmcInnerGhmcStep`'s acceptance against an OpenMM-free bridge
  (`TestNcmcExplicitSolvent.cpp:289` cites the private method as the target it
  mirrors, not links). `TestNCMCWork.cpp` and `TestFreeJointKEPump.cpp` also name
  `ncmcMove`/`protocolLambda`, but those references are the free function
  `robo::ncmc::protocolLambda` (from `NCMCProtocol.hpp`), not `World::` members,
  and neither `#include`s `World.hpp` - so they are unaffected. `configureNcmc`
  stays declared in the unchanged `World.hpp`. No include-path fix is required by
  any test.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0.
- `nm` diff on public symbols: empty (both `configureNcmc` overloads keep their
  mangled names; the five NCMC internals and `ncmcDebugEnabled` are non-public).
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- `NcmcMove.cpp` <= 600 LOC (~490 expected). No header added.
- Comment-stripped before/after diff of `World.cpp` + `NcmcMove.cpp` is empty
  apart from the moved blocks at their new location.
