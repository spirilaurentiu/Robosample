# SPLIT-C3: extract the GibbsSweep world-step primitive from Context.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/Context.cpp` (1335 LOC) writes the same per-World Gibbs-block step inline
three times, once in each sweep loop:

1. `runREX` - the coordinate-swap oracle sweep, iterating by replica index `r`
   over `worlds_` (`src/Context.cpp:457-490`; the step at **460-463**).
2. `runDrivenRound` - the driven-round equilibrium segment, iterating by
   thermodynamic-state index `k` over `st.worldIndexes`, skipping driven worlds
   (`src/Context.cpp:977-1021`; the step at **990-993**).
3. `RunREX` - the label-swap REMC/Default sweep, iterating by state index `k`
   (`src/Context.cpp:1112-1143`; the step at **1122-1125**).

The genuinely identical, duplicated unit inside all three loops is the INV-3
coordinate-currency exchange: push a coordinate buffer into a World, run one
sampling move, read the resulting coordinates back into the same buffer. Modulo
`w->` vs `w.` (site 1 dereferences a `unique_ptr`; sites 2-3 bind `World& w =
*worlds_[worldIx]`) and which buffer is passed (`replicaCoords_[r]` vs
`replicas_[r].atomsLocations`), the four lines are byte-identical:

```
w.setAtomsLocationsInGround(coords);
const bool accepted = w.generateSample();
const robo::Vec3* p = w.getAtomsLocationsInGround();
std::copy(p, p + systemTopology.numAtoms, coords.begin());
```

This ticket extracts exactly those four lines as the one primitive. It is pure
code motion: the primitive body is token-identical to each inline copy, and each
call site keeps its own surrounding scaffolding (schedule source, temperature and
per-world timestep/mdSteps/acceptRejectMode setters, the driven-world skip, the
BAT-anchor accumulation, the verbose print, and the `.moves.csv` logging).
Collapsing the *loops* is not in scope - their scaffolding differs and merging it
would require injecting site-specific behavior, which exceeds pure motion.

The primitive's entire contract is INV-3, restated in full:

- **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
  between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
  currency. The primitive is the operational statement of this invariant: the
  coordinate buffer is the only thing that crosses the World boundary in and out.
  `setAtomsLocationsInGround` re-fits every rigid-body frame and the `X_PF`/`X_BM`
  transforms from the pushed coordinates; `generateSample` runs one HMC move;
  `getAtomsLocationsInGround` reads the post-move coordinates. Nothing else
  persists across the call.

NOTE - occurrences that do NOT collapse (excluded deliberately, not missed):
the two pre-run initial-kick loops (`src/Context.cpp:421-429` in `runREX`,
`584-590` in `setupReplicaExchange`) also push/pull coordinates, but their core
op is `findGoodStartingPose()`, not `generateSample()` - a different operation,
not this primitive.

## Moves (exactly what)
- New free function -> `include/GibbsSweep.hpp` + `src/GibbsSweep.cpp` (target home
  `workflow/GibbsSweep`) [public free function in a `robo::gibbs` namespace]:

  ```
  namespace robo::gibbs {
  // Push coords into the World, run one sampling move, read the post-move
  // coordinates back into coords. Returns generateSample()'s accept flag.
  // INV-3: coords (nm, Ground frame, engine/OpenMM atom order) is the sole
  // inter-world currency; the World retains no per-replica state after this call.
  bool stepWorldInGround(World& w, std::vector<robo::Vec3>& coords, int numAtoms);
  }
  ```

  The definition is the four-line body above with `coords`/`numAtoms` in place of
  the site-specific buffer and `systemTopology.numAtoms`.

- Replace each inline copy with a call, leaving all surrounding lines untouched:
  - `runREX` (`src/Context.cpp:460-463`) ->
    `const bool accepted = robo::gibbs::stepWorldInGround(*w, replicaCoords_[r], systemTopology.numAtoms);`
    (`w` is the `unique_ptr` loop variable; `*w` yields the `World&`). Lines
    459 (`w->setTemperature`), 465-488 (verbose + `.moves.csv`) are unchanged.
  - `runDrivenRound` (`src/Context.cpp:990-993`) ->
    `const bool accepted = robo::gibbs::stepWorldInGround(w, replicas_[r].atomsLocations, systemTopology.numAtoms);`
    Lines 986-989 (schedule setters), 999 (`accumulateBatAnchorStats(w)`),
    1001-1015 (verbose) unchanged.
  - `RunREX` (`src/Context.cpp:1122-1125`) ->
    `const bool accepted = robo::gibbs::stepWorldInGround(w, replicas_[r].atomsLocations, systemTopology.numAtoms);`
    Lines 1118-1121 (schedule setters), 1127-1141 (verbose) unchanged.

Sequencing note: C3 runs before C4/C5. It edits the bodies of `runREX`,
`runDrivenRound`, and `RunREX` while they still live in `Context.cpp`; C4/C5 later
relocate those (now primitive-calling) definitions to their new TUs.

Oracle-coverage note: the `runDrivenRound` call site (line 990-993) sits in
driven-REX code that OQ-5 records as never compiled or run and TESTS.md section 7
records as throwing before its round loop. Its four-line rewrite is therefore NOT
exercised by the B1/B3 baselines - the coder verifies that site by inspection
(byte-identical body, `w`/buffer substitution only), not by a passing test. Flag any
divergence there as unverifiable-by-oracle, not as a passing change.

## Public API after this ticket
`include/GibbsSweep.hpp` exports exactly: `robo::gibbs::stepWorldInGround(World&,
std::vector<robo::Vec3>&, int)`. Nothing else. Header self-contained (IWYU):
forward-declare `class World;` is insufficient (the body calls `World` methods, so
the **`.cpp`** includes `"World.hpp"`); the **header** needs `#include
"robot_math.hpp"` (`robo::Vec3`) and `<vector>`, and forward-declares `class
World;`. `#pragma once`. `src/Context.cpp` gains `#include "GibbsSweep.hpp"`.

This is a new internal symbol, not part of the B0 public API contract; the nm
delta on **public** symbols is empty (the free function is engine-internal, not
bound to Python). `Context`'s class interface is unchanged.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames. The primitive body
  is token-identical to each inline copy; the three call sites lose exactly their
  four inline lines and gain one call each.
- Invariants that SHALL remain true: INV-3 (restated in full above). Because the
  primitive is the exact push/sample/pull triple and each caller keeps its own
  temperature/schedule/logging/BAT scaffolding, per-round trajectories are
  bit-identical - in particular the INVARIANT-EQUIV property (INV-8) between
  `runREX` and `RunREX` is preserved, since neither loop's World-visit order or
  RNG-consumption order changes.
- Build: compile `src/GibbsSweep.cpp` into the object library; the non-recursive
  `src/*.cpp` glob (`CMakeLists.txt`) reaches a flat `src/GibbsSweep.cpp`
  directly (or use `GLOB_RECURSE` for a subdir). One commit; formatting separate.

## Predicted test breakage (Phase-A include fixes only)
- None. `stepWorldInGround` is new and named by no test. No test includes
  `GibbsSweep.hpp`. The three drivers keep identical observable behavior, so the
  INVARIANT-EQUIV Python test `tests/test_rex_label_swap_equivalence.py`
  (asserts label-swap and coordinate-swap REMC produce identical per-round PE)
  still passes bit-for-bit, and the Level-1 example (B3) is unchanged. The only
  C++ test that `#include "Context.hpp"` (`tests/TestRoboticsOracleMolecule.cpp`)
  is unaffected - `Context.hpp` does not change.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0
  (B1/B2/B3), including `test_rex_label_swap_equivalence.py`'s per-round PE match.
- nm diff on public symbols vs B0: empty.
- include-cycle script clean vs B4 (`GibbsSweep.hpp` forward-declares `World`, so
  no new header cycle); clang-tidy / clang-format / IWYU clean on touched files.
- `GibbsSweep.hpp` <= 300 LOC (actual ~25); `GibbsSweep.cpp` <= 600 LOC (actual
  ~15). `Context.cpp` net change ~= 0 (four lines x three sites removed, three
  call lines added).
- Comment-stripped before/after diff: the new primitive at its location, the
  three one-line call replacements, and the added `#include`.
