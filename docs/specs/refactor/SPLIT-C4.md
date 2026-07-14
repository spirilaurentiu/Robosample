# SPLIT-C4: extract the ReplicaExchangeDriver and SwapAcceptance from Context.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/Context.cpp` (1335 LOC) carries ~700 lines of replica exchange entangled
with the rest of the workflow. Per the ARCHITECTURE.md section 1 design directive
(2026-07-12), **replica exchange is the intended center**: a run is `R` replicas
on a temperature ladder and single-replica is the degenerate `R = 1` case. The
extracted driver is therefore a **core** module, not an optional branch.

This ticket extracts two cohesive REX concerns as translation units:

- **ReplicaExchangeDriver** - the label-swap `RunREX` driver (`RUN_TYPE::REMC`/
  `Default`), the legacy coordinate-swap `runREX` oracle, replica/state setup, the
  label-swap and neighbour-pairing machinery, mixing, and the swap-matrix
  diagnostic.
- **SwapAcceptance** - the acceptance algebra `attemptREXSwap` (four run-types:
  REMC / RENEMC / RENE / REBASONTOP) and the driven-run precondition guard
  `checkInv7AndInv10Guards`.

The `Replica`/`ThermodynamicState`/`RUN_TYPE`/`ReplicaMixingScheme` value types
already live in `include/ReplicaExchange.hpp` (target `workflow/rex/Replica.hpp`)
and are **not** moved by this ticket.

The driven RENE/REBASONTOP round-loop (`runDrivenRound`/`driveReplica`/
`runInterleavedRemcSubround`) is extracted separately in **SPLIT-C5**. C4 moves
the RENE/REBASONTOP/RENEMC branches **of `attemptREXSwap`** (they are the
acceptance formulas, exercised directly by the algebra oracle and by C5's
round-loop), keeping them unchanged.

`runREX` stays. Per **OQ-3 (DECIDED: keep as oracle)**, the coordinate-swap
`runREX` is the deliberately frozen INVARIANT-EQUIV differential oracle for the
label-swap `RunREX` (`Context.hpp:84-90`, "Do NOT extend this method"). C4
relocates it unchanged alongside `RunREX`; it is not retired here (retirement is a
separate human-gated ticket, README section Phase A note).

Applicable invariants, restated in full:

- **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
  between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
  currency. (`RunREX` swaps **labels** - `swapThermodynamicStates` swaps the two
  inverse index maps, `src/Context.cpp:626-631` - never coordinates; `runREX`
  swaps `replicaCoords_` and is INV-3-non-compliant *by construction*, the exact
  property that makes it the differential oracle.)
- **INV-6 Constraint consistency.** The `G` (constraint Jacobian) assembly is
  identical across SHAKE, RATTLE, and the loop-closure Fixman log-det, so the
  correction matches the projection. `calcConstraintLogDet` returns 0 for acyclic
  molecules. (Manifested here as the "stored potential equals energy of stored
  coordinates" refresh before each swap, `src/Context.cpp:1149-1153`, 596-599.)
- **INV-7 BAT map/Jacobian agreement.** `applyBatScaling` and its Cartesian
  log-Jacobian read identical `(r,theta)` geometry; one anchor snapshot per round is
  shared by both swap partners so the paired map is an exact involution.
  (`checkInv7AndInv10Guards` enforces the Fixman-in-sampler precondition for
  driven run types, `src/Context.cpp:828-846`.)
- **INV-8 REX detailed balance.** `runREX` and `RunREX` produce equivalent
  sampling; `runREX` is the differential oracle for `RunREX`. Swap acceptance
  guards reject on NaN/inf. (The `forcedReject` guard, `src/Context.cpp:723-733`,
  SHALL move verbatim.)

## Moves (exactly what)
Cohesive translation-unit split (R3 pattern). Every method stays a `Context`
member - `include/Context.hpp` is unchanged, so the PyBind11-bound public symbols
(`runREX`, `RunREX`, `attemptREXSwap`, `checkInv7AndInv10Guards`,
`attemptedSwapsMatrix`/`acceptedSwapsMatrix`, the mixing setters) keep identical
signatures and the nm public delta is empty. Only **definitions** relocate.

Shared file-local constant:
- `kBoltzmann_kJ` (`src/Context.cpp:21`, anonymous namespace) is read by `runREX`
  (`src/Context.cpp:497-498`) and `attemptREXSwap` (`src/Context.cpp:648-649`),
  which land in two different TUs. C4 owns this constant. SPLIT-C2 runs earlier
  and deliberately leaves it in place (it moves only `boxFromReducedVectors` out
  of the anonymous namespace); C4 relocates `kBoltzmann_kJ` and then removes the
  now-empty `namespace {`/`}`. Move it to a shared internal header
  `src/workflow/rex/RexInternal.hpp` (`namespace { constexpr double kBoltzmann_kJ
  = 0.0083144626; }` or an `inline constexpr` in a `detail` namespace) included by
  both new TUs. Value unchanged - pure motion.

To `src/workflow/rex/ReplicaExchangeDriver.cpp` [Context member definitions,
verbatim]:
- `Context::runREX(int,int,int,bool)` incl. leading comment
  (`src/Context.cpp:397-535`) - the retained OQ-3 oracle. NOTE: the human-gated
  `SPLIT-C6` later isolates `runREX` from here into its own
  `workflow/rex/LegacyCoordSwapRex.{hpp,cpp}` TU; C4 co-locates it for now.
- `Context::setupReplicaExchange(RUN_TYPE)` (`src/Context.cpp:541-624`).
- `Context::swapThermodynamicStates(int,int)` (`src/Context.cpp:626-631`).
- `Context::prepareExchangePairs(int,int)` (`src/Context.cpp:756-763`).
- `Context::mixAllReplicas(int)` (`src/Context.cpp:765-779`).
- `Context::mixReplicas(int)` (`src/Context.cpp:781-801`).
- `Context::printSwapMatrix() const` (`src/Context.cpp:803-816`).
- `Context::RunREX(RUN_TYPE,int,int,int,bool)` incl. leading comment
  (`src/Context.cpp:1059-1168`).

To `src/workflow/rex/SwapAcceptance.cpp` [Context member definitions, verbatim]:
- `Context::attemptREXSwap(int,int)` incl. all four `switch` branches and the
  `forcedReject` NaN/inf guard (`src/Context.cpp:633-754`).
- `Context::checkInv7AndInv10Guards(RUN_TYPE) const`
  (`src/Context.cpp:823-880`).

Stays in the residual `Context.cpp`: `openmmPotential` (`src/Context.cpp:1170
-1177`), used by both new TUs and by C2/C5, remains a `Context` member. The REX
state members (`replicas_`, `thermodynamicStates_`, the two index maps, the swap
matrices, `rexRng_`/`rexUniform_`, mixing config, `batAnchorStats_`) stay declared
in `Context.hpp`; the moved definitions access them as ordinary members across TU.

## Public API after this ticket
No new exported symbol; `include/Context.hpp` byte-unchanged. PyBind11 bindings
(`src/PyBind11.cpp:588-656`: `run_rex`->`&Context::runREX`, `run_rex_label_swap`->
`&Context::RunREX`, `attempt_rex_swap`->`&Context::attemptREXSwap`,
`check_inv7_and_inv10_guards`, `attempted_swaps_matrix`/`accepted_swaps_matrix`)
resolve against the unchanged declarations - no PyBind edit. `RexInternal.hpp`
exports only the file-local `kBoltzmann_kJ` (internal linkage per TU). New TUs are
self-contained (IWYU): they carry the includes the moved bodies use -
`<algorithm>` (`std::swap`, `std::copy`), `<cmath>` (`std::exp`, `std::sqrt`,
`std::isnan`, `std::isinf`), `<cstdio>`, `<fstream>` (moves CSV in `runREX`),
`<limits>`, `<stdexcept>`, `<string>`, `<utility>`, `<vector>` - plus
`Context.hpp`, `ReplicaExchange.hpp`, `OpenMM.h` (only if referenced), and
`RexInternal.hpp`. `#pragma once` on the header.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames. `runREX` is moved
  **unchanged** (OQ-3 oracle); its INV-3-non-compliant coordinate swap is
  preserved exactly.
- `src/Context.cpp` loses lines 21 (the constant -> `RexInternal.hpp`), 397-535,
  541-624, 626-631, 633-754, 756-763, 765-779, 781-801, 803-816, 823-880,
  1059-1168. Nothing else changes. Includes now only used by the moved code
  (`<fstream>` if unused elsewhere) migrate to the new TUs.
- `Context`'s class interface is unchanged; public REX symbols keep their exact
  signatures, so the public-symbol nm set is untouched.
- Invariants that SHALL remain true: INV-3, INV-6, INV-7, INV-8 (restated above).
- Build: compile both new TUs into the object library; the non-recursive
  `src/*.cpp` glob (`CMakeLists.txt`) SHALL reach them (flat placement or
  `GLOB_RECURSE`). One commit for the move; formatting separate.

## Predicted test breakage (Phase-A include fixes only)
- None. `tests/TestRexAcceptanceAlgebra.cpp` does **not** instantiate a live
  `Context` or call `attemptREXSwap`/`checkInv7AndInv10Guards`; it **deliberately
  reimplements** the log-alpha formulas in-file (`TestRexAcceptanceAlgebra.cpp:15`
  "DELIBERATELY reimplements Context::attemptREXSwap's log-alpha formulas") and
  includes only `TestHelpers.hpp`/`robot_math.hpp` - it does not include
  `Context.hpp` and is unaffected by this move. `tests/test_rex_swap_acceptance_
  algebra.py` likewise reimplements the algebra in NumPy with no engine call.
  NOTE for the human: `attemptREXSwap`/`checkInv7AndInv10Guards` are `public`
  **only** so a reproducer *could* drive them directly (`Context.hpp:150-156`,
  162-169); no current test does. When a live-`Context` reproducer is later
  written (or the C++ algebra test is migrated to call the engine), it SHALL
  target the extracted **SwapAcceptance** TU's `Context::attemptREXSwap`/
  `checkInv7AndInv10Guards` - same symbols, same signatures, so the retarget is a
  no-op at the API level.
- `tests/test_rex_label_swap_equivalence.py` drives `run_rex`/`run_rex_label_swap`
  /`attempted_swaps_matrix` through the unchanged Python bindings; behavior is
  bit-identical (INVARIANT-EQUIV, INV-8), so it passes with no edit.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0
  (B1/B2/B3), including the 4-replica REMC swap-matrix proof run (B3) and
  `test_rex_label_swap_equivalence.py`'s per-round PE match.
- nm diff on public symbols vs B0: empty (all REX public symbols keep identical
  mangled names; internal churn allowed).
- include-cycle script clean vs B4; clang-tidy / clang-format / IWYU clean on
  touched files.
- `ReplicaExchangeDriver.cpp` <= 600 LOC (actual ~430 - if it exceeds the cap,
  claim the dispatch/driver exception, or split `runREX` into its own TU);
  `SwapAcceptance.cpp` <= 600 LOC (actual ~180); `RexInternal.hpp` <= 300 LOC.
  `Context.cpp` shrinks by ~610 LOC.
- Comment-stripped before/after diff: only the moved blocks at their new
  locations, the constant in `RexInternal.hpp`, and the added `#include`s.
