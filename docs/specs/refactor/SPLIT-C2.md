# SPLIT-C2: extract the OutputWriter from Context.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/Context.cpp` (1335 LOC) mixes World lifecycle, OpenMM bring-up, startup
validation, Gibbs scheduling, two REX drivers, and output I/O. This ticket
extracts the output I/O: the per-round moves/energy CSV rows, the per-replica DCD
frame with whole-molecule periodic imaging, the reaction-force CSV rows, and the
DCD box-vector conversion helper.

Both REX drivers write outputs through the same core. `runREX` calls
`writeOutputs(replica, round, verbose)` (`src/Context.cpp:508`, wrapper at
1179-1181), and `RunREX` calls `writeOutputsCore(k, ...)` directly with a
`Replica`'s committed coordinates indexed by thermodynamic-state
(`src/Context.cpp:1161`). The shared core `writeOutputsCore`
(`src/Context.cpp:1188-1281`) does the periodic-imaging DCD scatter; the anonymous
-namespace helper `boxFromReducedVectors` (`src/Context.cpp:25-48`) converts
reduced lattice vectors to a CHARMM/DCD box and is used only there.

The behavioral contract to preserve verbatim: the `.moves.csv` schema is written
by the drivers (C4), not here; here it is the `.<idx>.csv` energy row
(`round,idx,T,pe`), the whole-molecule periodic imaging (COM folded in c->b->a
order onto the **output copy only**, `replicaCoords_`/`Replica` state left
unwrapped), the `atomsPrmtopIndex` permutation on scatter, and the
`.<idx>.reactions.csv` 10-column schema. The applicable invariant is INV-3.

- **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
  between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
  currency. (`writeOutputsCore` reads a coordinate buffer and writes files; the
  reaction-reporter re-sync in `runREX`, `src/Context.cpp:527-529`, is read-only
  and stays in the driver. Imaging mutates only `dcdScratch_`, never the
  authoritative coordinates - this is the property that keeps output side-effect
  free.)

## Moves (exactly what)
Cohesive translation-unit split (R3 pattern). The three writer methods stay
`Context` members; `boxFromReducedVectors` is a genuine free helper that moves as
a real symbol. `include/Context.hpp` is unchanged (the three method declarations
stay).

- Anonymous-namespace `boxFromReducedVectors(const std::vector<double>&)` (the
  function body only, `src/Context.cpp:25-48`) -> `src/workflow/OutputWriter.cpp`
  anonymous namespace [internal linkage, verbatim]. It is used only by
  `writeOutputsCore`, so it travels with it. **Do NOT move the sibling
  `kBoltzmann_kJ` (`src/Context.cpp:21`) or the enclosing `namespace {`/`}`
  braces (lines 20, 50):** that constant is read by REX code (`runREX`,
  `attemptREXSwap`) and is relocated by SPLIT-C4 to `RexInternal.hpp`. C2 runs
  before C4, so C2 lifts only `boxFromReducedVectors` out of the anonymous
  namespace and leaves the namespace in `Context.cpp` holding `kBoltzmann_kJ`;
  C4 empties and removes it. Burying `kBoltzmann_kJ` in `OutputWriter.cpp`'s
  anonymous namespace would make it unreachable by the REX TUs and break the
  build.
- `Context::writeOutputs(int, int, bool)` definition
  (`src/Context.cpp:1179-1181`) -> `OutputWriter.cpp` [Context member definition].
- `Context::writeOutputsCore(int, int, bool, const std::vector<robo::Vec3>&,
  double)` definition, including its leading comment (`src/Context.cpp:1183-1281`)
  -> `OutputWriter.cpp` [Context member definition, verbatim].
- `Context::writeReactionRows(int, int, const std::vector<ReactionSample>&)`
  definition, including its leading comment (`src/Context.cpp:1283-1309`) ->
  `OutputWriter.cpp` [Context member definition, verbatim].

The DCD/CSV **truncation** and DCD-writer construction in `Context::initialize`
(`src/Context.cpp:206-229`) stay in `initialize` (they are lifecycle setup, not
per-frame writing) and are out of scope for this ticket. Declarations stay at
`include/Context.hpp:273` (`writeOutputs`), `:283` (`writeOutputsCore`), `:288`
(`writeReactionRows`), all private.

The shared Context members these bodies read - `dcdWriters_`, `dcdScratch_`,
`systemTopology`, `writeCounter_` (touched only by the drivers), and the private
`openmmPotential(...)` helper - remain on `Context` and stay declared in
`Context.hpp`; the moved definitions call them across TU as ordinary member
access. `openmmPotential` (`src/Context.cpp:1170-1177`) stays in the residual
`Context.cpp` (also used by C4/C5).

## Public API after this ticket
No new exported symbol. `include/Context.hpp` byte-unchanged. A thin
`src/workflow/OutputWriter.hpp` MAY be omitted. The `.cpp` includes exactly what
the moved bodies reference; the IWYU exit check (VERIFY criterion 6) is
authoritative, so the list below is the expected result, not a guess to relax it:
`#include "Context.hpp"`, `#include "DCDWriter.hpp"`. `OpenMM.h` is needed only if a
moved body names an `OpenMM::` type (`openmmPotential` stays in `Context.cpp`, so it
likely is not). `PeriodicBox.hpp` is not required - the imaging math is inline, not
via `robo::pbc`. Standard headers the moved code uses: `<cmath>` (`std::floor`,
`std::sqrt`, `std::acos`, `std::max`/`min` in `boxFromReducedVectors`), `<cstdio>`
(`std::printf`), `<fstream>`, `<string>`, `<vector>`. Add `<algorithm>` only if a
moved body uses it. `#pragma once` if a header is added.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/Context.cpp` loses the `boxFromReducedVectors` definition (lines 25-48),
  1179-1181, 1183-1281, 1283-1309; nothing else changes. The anonymous-namespace
  braces (lines 20, 50) and `kBoltzmann_kJ` (line 21) STAY - SPLIT-C4 relocates
  the constant and removes the now-empty namespace. It keeps `#include "DCDWriter.hpp"`
  only if still used elsewhere (the DCD-writer construction in `initialize` uses
  `dcd::Writer`, so `DCDWriter.hpp` stays included by `Context.cpp` via
  `Context.hpp`). Move the `<cmath>` box-math dependency to the new TU.
- `Context`'s class interface is unchanged; the three methods remain private
  members, so the public-symbol nm set is untouched.
- Invariants that SHALL remain true: INV-3 (restated above). The periodic imaging
  writes only `dcdScratch_` and the output files; the authoritative coordinate
  buffers are never mutated - bit-identical to today.
- Build: compile `src/workflow/OutputWriter.cpp` into the object library; the
  non-recursive `src/*.cpp` glob (`CMakeLists.txt`) SHALL reach it (flat
  placement or `GLOB_RECURSE`). One commit for the move; formatting separate.

## Predicted test breakage (Phase-A include fixes only)
- None. `writeOutputs`/`writeOutputsCore`/`writeReactionRows`/
  `boxFromReducedVectors` are private/internal and named by no test (grep of
  `tests/`: no hit). Output files are validated only through the drivers by the
  Python suite (`test_rex_label_swap_equivalence.py` reads `.0.csv`/`.dcd`,
  `tests/test_openmm_potential_energy.py` etc.); those go through unchanged
  bindings and read byte-identical files, so B3 holds with no test edit.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0
  (B1/B2/B3), including every `.<idx>.csv`, `.dcd`, and `.reactions.csv` byte-for
  -byte (DCD frames tolerance-identical per B3).
- nm diff on public symbols vs B0: empty (all moved symbols were internal).
- include-cycle script clean vs B4; clang-tidy / clang-format / IWYU clean on the
  touched files.
- `OutputWriter.cpp` <= 600 LOC (actual ~150). `Context.cpp` shrinks by ~160 LOC.
- Comment-stripped before/after diff: only the moved blocks at their new location
  and the added `#include`s.
