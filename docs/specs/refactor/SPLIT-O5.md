# SPLIT-O5: extract the initialize() system-construction body (OpenMMSystemBuilder) from OpenMMContext.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`OpenMMContext::initialize()` (`src/OpenMMContext.cpp` lines 417-652, ~235 LOC) is
the OpenMM `System`/`Context`/`Integrator` construction procedure. It mixes six
concerns in sequence: particle creation, periodic-box setup, virtual-site
declaration + massless-particle safety net, GBSA usability gating, force-group
assignment and force wiring (delegating to the standard and alchemy builders), and
integrator + platform + context creation. This ticket extracts that procedure into
an `OpenMMSystemBuilder` so `OpenMMContext` retains only ownership of the built
objects and the runtime accessors.

**Sequencing - O5 runs LAST among the OpenMMContext splits.** `initialize()` calls
symbols the other five tickets relocate:
- MTSIntegrator (SPLIT-O1): construction at line 599.
- Standard force builders (SPLIT-O3): `createNonbondedForce`,
  `createGBSAOBCForce`, `createCustomNonbondedForce`, and the six bonded builders
  at lines 556-593.
- Alchemy factory (SPLIT-O4): the `alchemyEnabled` dispatch at lines 564-580
  (`createAlchemyDecouplingForces` / `createAlchemyCorrectionForce`).
- GpuKinematics (SPLIT-O6): the `ROBO_CUDA_KINEMATICS` env read at lines 642-648
  (and the `releaseGpuKinematics()`/`cudaKinematicsEnabled_` coupling).

O5 SHALL land only after O1, O3, O4, O6 are merged, so it extracts a body whose
callees already live in their target units. Extracting O5 first would force it to
carry those callees and conflict with the later tickets.

INV-1/INV-2 (force->wrench) are not exercised by construction - `initialize()`
builds the potential, not the per-body reduction. The contracts this body SHALL
preserve are its own, all present as comments that move verbatim: box vectors set
BEFORE context creation (PME grid), virtual sites declared before context, the
massless-real-particle guard (throws on violation, lines 483-501), GBSA<->periodic mutual
exclusion (lines 507-513), and the MTS/separate/single force-group assignment
policy (lines 529-552).

## Moves (exactly what)
Target `bridge/OpenMMSystemBuilder.{hpp,cpp}`. The builder produces the
`System` + `Integrator` + `Context` (and the `forceGroupLabels`, `hasVirtualSites`
side outputs) from a `SystemTopology` plus the caller's toggles
(`useMTS`, `mtsInnerSubsteps`, `separateForceGroups`, and the alchemy factory), and
`OpenMMContext::initialize()` becomes a thin driver that invokes it and stores the
results.

- Body of `OpenMMContext::initialize(const SystemTopology&)`
  (`src/OpenMMContext.cpp` lines 417-652) -> `OpenMMSystemBuilder` build routine.
  This includes:
  - particle loop (lines 421-424)
  - periodic-box block (lines 435-459) - reads already-reduced
    `systemTopology.boxVectors`; does NOT call the removed
    `computePeriodicBoxVectors_Context` (SPLIT-O2), so no dependency there
  - virtual-site declaration + massless guard (lines 461-501)
  - GBSA usability gating (lines 506-527)
  - `addForce` lambda + force-group policy (lines 529-552)
  - force wiring via the O3 free functions and the O4 alchemy factory
    (lines 554-594)
  - integrator selection (MTS/Verlet, lines 596-604)
  - platform selection `#if USE_CPU/REFERENCE/OPENCL/CUDA` (lines 606-624)
  - context creation + error handling (lines 626-635)
  - `initialized = true` and the CUDA env read (lines 637-648) - the env read
    sets `cudaKinematicsEnabled_`, an `OpenMMContext` member, so it either stays
    in the thin `initialize()` driver or the builder returns the flag.

The DECLARATION `OpenMMContext::initialize` (`include/OpenMMContext.hpp:24`) stays;
its body shrinks to: call `OpenMMSystemBuilder::build(...)`, move the returned
`system`/`integrator`/`context`/`forceGroupLabels`/`hasVirtualSites`/`numAtoms`
into the members, run the CUDA env read, return the success bool.

AMBIGUITY TO FLAG: `initialize()` writes many `OpenMMContext` members
(`system`, `context`, `integrator`, `numAtoms`, `forceGroupLabels`,
`hasVirtualSites`, `initialized`, `cudaKinematicsEnabled_`) and reads toggles
(`useMTS`, `mtsInnerSubsteps`, `separateForceGroups`) + the alchemy factory. Pure
motion therefore requires the builder to take those toggles + factory as inputs
and return the built objects as outputs (e.g. a `BuildResult` struct or out-params).
This is a mechanical parameter-threading, not a logic change, but it is the one
place where the extraction is not a literal cut-paste. The member SET and the
public `initialize` signature are unchanged. Confirm the input/output bundle shape
before implementing.

## Public API after this ticket
`OpenMMContext::initialize(const SystemTopology&) -> bool` is PRESERVED verbatim
(signature and observable effect). `bridge/OpenMMSystemBuilder.hpp` exports the
builder entry point (a free function or a `struct OpenMMSystemBuilder` with a
`build` method) taking `const SystemTopology&` + the toggles + the
`AlchemyForceFactory&` and returning the built `System`/`Integrator`/`Context` and
the `forceGroupLabels`/`hasVirtualSites` outputs. Header self-contained:
`#include "OpenMM.h"`, `#include "TopologyElements.hpp"`,
`#include "MTSIntegrator.hpp"`, `#include "ForceFactory.hpp"`,
`#include "AlchemyForceFactory.hpp"`, `<memory>`, `<vector>`, `<string>`.
`#pragma once`.

## Constraints
- Behavior-preserving relocation. Zero logic changes to construction order, the
  box/virtual-site/GBSA guards, the force-group policy, or platform selection.
  Zero renames. The only non-literal part is threading the members through the
  builder's inputs/outputs (a mechanical wiring, not a semantic change).
- `src/OpenMMContext.cpp` gains `#include "OpenMMSystemBuilder.hpp"` and reduces
  `initialize()` to the thin driver; the standard/alchemy/MTS includes it needed
  only for construction migrate to the builder.
- Invariants that SHALL remain true: box set before context (PME); virtual sites
  before context; massless-real-particle guard fires; GBSA excluded under
  periodic; force-group assignment policy identical; INV-1...INV-9 untouched.
- Ordering dependency: MERGE ONLY AFTER O1, O3, O4, O6. New builder files
  auto-globbed into `robosample_objects` (`CMakeLists.txt`); no CMake
  edit. One commit; formatting separate.

## Predicted test breakage (Phase-A include fixes only)
- None. `initialize()` stays a public `OpenMMContext` method with an unchanged
  signature. `tests/TestAlchemy.cpp` (`omm.initialize(s)`, line 233) and every
  other test that brings up OpenMM do so through this public call; they never name
  the extracted construction helpers. No test includes `OpenMMSystemBuilder.hpp`.

## Exit criteria (machine-checked)
- Full build passes in all reference configs; test suite / baseline outputs
  identical to Stage-0 (B1/B2/B3). The Level-1 example brings up an identical
  System (same force groups, same energies within the 1e-6 band).
- nm diff on public symbols vs B0: empty (public `initialize` preserved; the
  builder is new internal API).
- include-cycle script clean vs B4; no new upward layer edge; clang-tidy /
  clang-format / IWYU clean.
- `OpenMMSystemBuilder.hpp` <= 300 LOC; `OpenMMSystemBuilder.cpp` <= 600 LOC
  (actual ~240). `OpenMMContext.cpp` `initialize()` shrinks to a ~20-LOC driver;
  with O1-O4/O6 also landed, `OpenMMContext.cpp` drops from 1165 LOC to a residual
  well under the 600-LOC `.cpp` cap.
