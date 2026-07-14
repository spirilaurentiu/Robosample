# SPLIT-O3: extract the standard force builders (ForceFactory) from OpenMMContext.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/OpenMMContext.cpp` (1165 LOC) inlines all OpenMM `Force` construction. This
ticket extracts the **standard** (non-alchemy) force builders and the shared
exclusion helper into a `ForceFactory`. These builders are pure functions of the
`SystemTopology` SoA plus the static predicate `OpenMMContext::isPeriodic`: they
read no `OpenMMContext` mutable member state, allocate an `OpenMM::*Force*`, and
return it. That property makes the extraction clean code motion into free
functions (the same hoisting pattern used for `robo::pbc` and the analytic force
bridge).

The alchemy/decoupling builders are NOT in scope here - they read alchemy member
state and move under SPLIT-O4. `isPeriodic` stays a public static on
`OpenMMContext`; `ForceFactory` calls it (or a shared copy) - see Constraints.

No ARCHITECTURE.md force->wrench invariant is exercised by these builders: they
define the potential (energy/force field), not the per-body reduction. INV-1/INV-2
live in the reduction path (ForceBridge / CUDA kernel), untouched here. The
correctness contracts these builders SHALL preserve are the ones already in their
comments: NBFIX table well-formedness checks
(`createCustomNonbondedForce`, lines 859-885), the CPU shared-neighbor-list
exclusion rule (`addStandardExclusions`, lines 954-964), the `*2.0` stiffness
convention on bonds/angles/Urey-Bradley, and the periodic->CutoffPeriodic collapse
for custom forces.

## Moves (exactly what)
Each moves to `bridge/ForceFactory.{hpp,cpp}` as a free function (namespace
`robo::forcefactory` or equivalent), taking `const SystemTopology&` and returning
the same pointer type. Declarations leave `include/OpenMMContext.hpp`.

- `createNonbondedForce` (`src/OpenMMContext.cpp` lines 755-819; decl
  `include/OpenMMContext.hpp:219`) -> `ForceFactory` [public free fn]
- `createGBSAOBCForce` (lines 821-857; decl `:220`) -> `ForceFactory` [public]
- `createCustomNonbondedForce` (lines 859-916; decl `:221-222`) -> `ForceFactory`
  [public] - calls `addStandardExclusions` and `isPeriodic`
- `addStandardExclusions` (lines 954-964; decl `:236-237`, already `static`) ->
  `ForceFactory` [public free fn]
- `createHarmonicBondForce` (lines 990-1001; decl `:238-239`) -> `ForceFactory`
- `createHarmonicAngleForce` (lines 1003-1015; decl `:240-241`) -> `ForceFactory`
- `createPeriodicTorsionForce` (lines 1017-1030; decl `:242-243`) -> `ForceFactory`
- `createImproperHarmonicTorsionForce` (lines 1032-1049; decl `:244-245`) ->
  `ForceFactory`
- `createCMAPTorsionForce` (lines 1051-1076; decl `:246-247`) -> `ForceFactory`
- `createUreyBradleyForce` (lines 1078-1088; decl `:248-249`) -> `ForceFactory`

Call sites inside `initialize()` (`src/OpenMMContext.cpp` lines 556-593:
`createNonbondedForce`, `createGBSAOBCForce`, `createCustomNonbondedForce`, and
the six bonded builders) change from member calls to
`robo::forcefactory::create*Force(systemTopology)`. `initialize()` itself is
extracted later (SPLIT-O5); after O3 it simply calls the new free functions.

`isPeriodic` (`include/OpenMMContext.hpp:146-149`, `static`) stays on
`OpenMMContext`. `createNonbondedForce` and `createCustomNonbondedForce`
reference it; they call `OpenMMContext::isPeriodic(...)`.

## Public API after this ticket
`bridge/ForceFactory.hpp` exports the ten free functions listed above (same
signatures, `const SystemTopology&` in, `OpenMM::*Force*` out). Header
self-contained: `#include "OpenMM.h"`, `#include "TopologyElements.hpp"`
(`SystemTopology`, `NonbondedMethod`, `ONE_4PI_EPS0`), and forward reference to
`OpenMMContext::isPeriodic` via `#include "OpenMMContext.hpp"` (or move
`isPeriodic` to a shared `NonbondedMethod` helper - see AMBIGUITY). `<set>`,
`<sstream>`, `<iomanip>`, `<string>`, `<vector>`, `<stdexcept>` as used.
`#pragma once`.

`OpenMMContext` loses ten **private** member functions from its declaration.
Public exported-symbol set is unchanged (these were private; nm public delta SHALL
be empty).

AMBIGUITY TO FLAG: `ForceFactory` including `OpenMMContext.hpp` only to reach the
static `isPeriodic` reintroduces a `bridge -> bridge` include the layer map would
rather avoid. Cleanest alternative: relocate `isPeriodic` to
`model/TopologyElements.hpp` beside `NonbondedMethod` (it is a pure predicate on
that enum), or to `ForceFactory` itself, leaving a forwarder on `OpenMMContext`.
Moving `isPeriodic` touches a public static symbol and therefore needs the same
human-approval flag as any API change; left as a decision point rather than
decided by fiat.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames (the free functions
  keep the exact names, minus the `OpenMMContext::` qualifier).
- `src/OpenMMContext.cpp` gains `#include "ForceFactory.hpp"` and loses the ten
  builder bodies; the `initialize()` call sites switch to the free functions.
  Nothing else changes.
- Invariants that SHALL remain true: no INV-# is force-reduction-adjacent here;
  the builders' internal contracts (NBFIX validation, exclusion mirroring,
  `*2.0` stiffness, periodic->CutoffPeriodic collapse) move verbatim and stay true.
- New `ForceFactory.cpp`/`.hpp` are auto-globbed (`CMakeLists.txt`) into
  `robosample_objects`; no CMake edit. One commit for the move; formatting
  separate.

## Predicted test breakage (Phase-A include fixes only)
- None. The moved functions are all **private** members of `OpenMMContext`; no
  test can name them. `tests/TestAlchemy.cpp` links the production
  `OpenMMContext.o` and calls only public `initialize()` /
  `computePotentialEnergyByGroup` / `ForceGroupEnergy` (`TestAlchemy.cpp:151,222-233`),
  which still resolve. No test includes `ForceFactory.hpp`. If a test target's
  link line lists object files explicitly rather than the whole object lib, the
  new `ForceFactory.o` joins that list - verify the test-link CMake (`CMakeLists.txt`
  around 594-616) links `robosample_objects` wholesale (it does via the object
  lib), so no per-test edit is needed.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0
  (B1/B2/B3). The Level-1 example and `TestAlchemy` energies stay within the
  existing 1e-6 force-group band.
- nm diff on public symbols vs B0: empty (unless the optional `isPeriodic` move is
  taken, which is a separate human-approved delta).
- include-cycle script clean vs B4; no new upward layer edge; clang-tidy /
  clang-format / IWYU clean on touched files.
- `ForceFactory.hpp` <= 300 LOC; `ForceFactory.cpp` <= 600 LOC (actual ~330).
  `OpenMMContext.cpp` shrinks by ~250 LOC.
- Comment-stripped before/after diff: only the moved builders at their new
  location, the `initialize()` call-site requalification, and the added `#include`.
