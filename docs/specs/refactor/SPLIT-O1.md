# SPLIT-O1: extract the MTSIntegrator class from OpenMMContext.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/OpenMMContext.cpp` (1165 LOC) mixes six OpenMM-adapter responsibilities:
the r-RESPA multiple-timestep integrator, the periodic-box math, the standard
force builders, the alchemy/decoupling force builders, the `initialize()` system
construction body, and the fused CUDA kinematics pipeline. This ticket extracts
the first: the file-local `MTSIntegrator` class - a self-contained
`OpenMM::CustomIntegrator` subclass with no dependency on `OpenMMContext` state.

`MTSIntegrator` is a reversible, symplectic r-RESPA Verlet: slow force groups get
fewer substeps than fast groups. It is constructed once, by `initialize()`, when
MTS is enabled. Extracting it isolates a leaf concept and shrinks the biggest
`.cpp` by ~46 LOC with no behavioral surface.

No ARCHITECTURE.md invariant governs this class directly. INV-1/INV-2 (force
reduction) and the box conventions live in other extractions; nothing here reads
per-body wrenches or box vectors. The only contract to preserve: the integrator
SHALL remain reversible + symplectic so it is safe inside HMC acceptance (the
existing class comment, `OpenMMContext.cpp:37-39`).

## Moves (exactly what)
- Header comment + `class MTSIntegrator : public OpenMM::CustomIntegrator`
  including its ctor and the private `createSubsteps` recursion
  (`src/OpenMMContext.cpp` lines 37-82) -> `bridge/MTSIntegrator.hpp` [public
  header, in-class definitions, or a `.hpp`/`.cpp` pair]. The class is currently a
  file-local (internal-linkage-adjacent) type; keep it a plain class in the new
  header. No `detail::`/anonymous-ns wrapping is required - it is referenced by
  name at the construction site.

The construction site stays in `initialize()`:
`std::make_unique<MTSIntegrator>(0.001, tiers)` (`src/OpenMMContext.cpp` line 599,
inside the block 596-604) is unchanged; it resolves against the new header.

## Public API after this ticket
`bridge/MTSIntegrator.hpp` exports: `class MTSIntegrator` (ctor
`MTSIntegrator(double stepSize, std::vector<std::pair<int,int>> groups)`; all
else private). Nothing else. Header self-contained (IWYU): `#include "OpenMM.h"`
(base class + `addPerDofVariable`/`addComputePerDof`/`addConstrain*`), `<algorithm>`
(`std::sort`, line 47), `<stdexcept>` (`std::invalid_argument`), `<string>`
(`std::to_string`), `<utility>` (`std::pair`), `<vector>`. `#pragma once`.

`OpenMMContext`'s public symbol set is unchanged: the MTS *control* surface
(`setMTS`, `getUseMTS`, `useMTS`, `mtsInnerSubsteps`, `kMtsSlowGroup`,
`kMtsFastGroup`) stays on `OpenMMContext` - only the integrator class moves.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/OpenMMContext.cpp` gains `#include "MTSIntegrator.hpp"` and loses lines
  37-82; nothing else in it changes. The `<algorithm>`/`<stdexcept>` includes it
  already carries stay (still used by `enableAlchemy`, `initialize`).
- Invariants that SHALL remain true: none of INV-1...INV-9 are touched. The class
  contract (reversible + symplectic r-RESPA, safe under HMC acceptance) is
  preserved bit-for-bit because the integrator string expressions move verbatim.
- New file is picked up automatically by the source glob (`CMakeLists.txt`,
  object lib `robosample_objects` line 321) if a `.cpp` is added; a header-only
  move needs no CMake change. One commit for the move; formatting fixes separate.

## Predicted test breakage (Phase-A include fixes only)
- None. `MTSIntegrator` is not named by any test (grep of `tests/`:
  `MTSIntegrator` appears only in `src/OpenMMContext.cpp` and the out-of-scope
  `Robosample/` tree). No test includes it, calls it, or links against it as a
  distinct symbol. `TestAlchemy.cpp` links the production `OpenMMContext.o` and
  uses only public `initialize()`/`ForceGroupEnergy`, so it is unaffected.

## Exit criteria (machine-checked)
- Full build passes (all reference configs); test suite / baseline outputs
  identical to Stage-0 (B1/B2/B3).
- nm diff on public symbols vs B0: empty (`MTSIntegrator` was never an exported
  public symbol of the API contract; internal churn is allowed).
- include-cycle script clean vs B4; clang-tidy / clang-format / IWYU clean on the
  two touched files.
- `MTSIntegrator.hpp` <= 300 LOC (actual ~50). `OpenMMContext.cpp` shrinks by ~46
  LOC.
- Comment-stripped before/after diff: the only non-empty hunks are the moved
  class at its new location and the added `#include`.
