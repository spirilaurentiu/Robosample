# SPLIT-O6: extract the fused CUDA robot-kinematics pipeline (GpuKinematics) from OpenMMContext.cpp - CUDA-only

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/OpenMMContext.cpp` contains the opt-in fused CUDA robot-kinematics pipeline
(GPU Cartesian-kinematics campaign; its spec is not yet in the repo, so the
behavior below is recovered from source): two nvrtc kernels and their host
driver. The whole unit is compiled only under `USE_CUDA` (guarded by
`#if USE_CUDA`). It is self-contained - a file-static device object plus the
`OpenMMContext` methods that drive it - and can move as one CUDA-only translation
unit, shrinking the biggest `.cpp` by ~230 LOC and isolating all nvrtc / CudaContext
internals behind one file.

The pipeline: kernel **K1** `pushPositions` writes `posG = X_GB*station + p`
straight into OpenMM's device `posq`; kernel **K2** `reduceForces` reduces the
device fixed-point force buffer into per-body spatial forces. Between and around
them, `syncInvOrder` rebuilds the atom->device-slot map whenever OpenMM's spatial
sort reorders atoms.

This ticket carries two ARCHITECTURE.md invariants because K2 IS the device half
of the force->wrench reduction:

- **INV-1 (force->wrench convention).** Per-body spatial force in
  `RobotState::bodyForceG` is `(torque about the body origin, net force)` in the
  Ground frame. Both the host reduction (`ForceBridge.hpp`) and the CUDA
  `reduceForces` kernel MUST produce this identically. In K2 the moment is
  `angular += r x f` with `r = X_GB.R * station` (about the body origin) and
  `linear += f` (`src/OpenMMContext.cpp:158-168`) - this convention moves verbatim.
- **INV-2 (virtual sites).** Virtual-site forces are already projected onto parent
  atoms by OpenMM; the reduction skips massless particles to avoid double
  counting. In K2 this is `if (isVirtual[a]) continue;`
  (`src/OpenMMContext.cpp:153`) - moves verbatim. K1 deliberately does NOT skip
  massless atoms (it writes every atom's `posq`; the comment at lines 108-113
  states why) - that asymmetry moves verbatim too.

Sequencing note: SPLIT-DEDUP-FORCEREDUCER runs BEFORE O6 (README section Phase-A step 2)
and validates the CUDA `reduceForces` kernel against the new host
`reduceAtomForcesToBodies` while the kernel still lives at
`src/OpenMMContext.cpp:132-170`. O6 then relocates that kernel; the dedup test
SHALL keep passing against the kernel at its new home.

## Moves (exactly what)
Target `bridge/GpuKinematics.{hpp,cpp}`, entirely under `#if USE_CUDA`. The
cleanest pure motion keeps the OpenMM-internal reach-through out of the header:
the anonymous-namespace device object and the `OpenMMContext::` CUDA method
DEFINITIONS both live in `GpuKinematics.cpp` (which `#include "OpenMMContext.hpp"`),
so the `OpenMMContext` public method DECLARATIONS stay in the header unchanged.

Device sub-object and helpers (currently `#if USE_CUDA` anonymous namespace):
- CUDA-only includes block (`src/OpenMMContext.cpp` lines 19-35: `CudaPlatform.h`,
  `ContextImpl.h`, `CudaContext.h`, `ComputeContext/Array/Kernel/Program.h`,
  `ContextSelector.h`, `<cstdlib>`, `<map>`, `<memory>`) -> `GpuKinematics.cpp`
- `kKinematicsKernelSource` (K1 `pushPositions` + K2 `reduceForces` nvrtc source
  string, lines 84-171) -> `GpuKinematics.cpp` [anonymous ns]
- `struct GpuKinematics` incl. `syncInvOrder` (lines 173-218) -> `GpuKinematics.cpp`
  [anonymous ns]
- `std::unique_ptr<GpuKinematics> gGpuKin` file-static (lines 220-222) ->
  `GpuKinematics.cpp` [anonymous ns, file-static]
- `getCudaComputeContext(OpenMM::Context*)` (lines 224-233) -> `GpuKinematics.cpp`
  [anonymous ns]

`OpenMMContext` member-method definitions (move the bodies; declarations stay in
`include/OpenMMContext.hpp`):
- `OpenMMContext::cudaKinematicsAvailable` (lines 238-244; decl `:168`) ->
  `GpuKinematics.cpp`
- `OpenMMContext::releaseGpuKinematics` (lines 246-258; decl `:294`) ->
  `GpuKinematics.cpp`
- `OpenMMContext::ensureKinematicsConstants` (lines 260-359; decl `:173-182`) ->
  `GpuKinematics.cpp`
- `OpenMMContext::pushBodyTransforms` (lines 361-387; decl `:185`) ->
  `GpuKinematics.cpp`
- `OpenMMContext::computeForcesAndEnergyOnDevice` (lines 389-399; decl `:189`) ->
  `GpuKinematics.cpp`
- `OpenMMContext::reduceForcesToBodies` (lines 401-415; decl `:193`) ->
  `GpuKinematics.cpp`

Left in place (unchanged): the `OpenMMContext` control surface
`setCudaKinematics`/`getCudaKinematics` (`:160-165`) and the member
`cudaKinematicsEnabled_` (`:293`); the `ROBO_CUDA_KINEMATICS` env read in
`initialize()` (`src/OpenMMContext.cpp:642-648`); and the `releaseGpuKinematics()`
call in `shutdown()` (`:27`, resolves to the relocated definition). `MTSIntegrator`
(lines 37-82) is SPLIT-O1, not part of this move.

## Public API after this ticket
`OpenMMContext`'s public CUDA-pipeline surface is unchanged (declarations stay in
`OpenMMContext.hpp`; only definitions move). `bridge/GpuKinematics.hpp` may be
empty of public API or expose nothing beyond an include guard - the entire device
type is anonymous-namespace file-static in the `.cpp`. If a header is warranted,
it forward-declares nothing OpenMM-typed and stays `#pragma once`.
`#else` no-op branches of each `OpenMMContext` CUDA method move with their
definitions so a non-CUDA build still links (the methods keep their `(void)`
argument casts).

## Constraints
- Pure code motion. Zero logic changes. Zero renames. The kernel source string,
  `syncInvOrder`, the `ContextSelector` scoping, and every `#if USE_CUDA`/`#else`
  guard move verbatim.
- Entire moved unit stays behind `#if USE_CUDA`. Under a non-CUDA build the new
  `.cpp` compiles to the no-op `#else` bodies exactly as today.
- `src/OpenMMContext.cpp` loses the moved bodies and the CUDA-only include block;
  it keeps the control setters and the env read. It does NOT gain a new include
  (the definitions live in `GpuKinematics.cpp`, which includes `OpenMMContext.hpp`).
- Invariants that SHALL remain true: **INV-1** (K2 wrench = `(rxf, f)` about body
  origin, Ground frame) and **INV-2** (K2 skips `isVirtual`; K1 writes all atoms),
  restated in full above. Also the reorder-safety contract of `syncInvOrder`
  (rebuild on live-order change, run before K1 and K2) moves verbatim.
- New `GpuKinematics.cpp` auto-globbed into `robosample_objects`
  (`CMakeLists.txt`); confirm CUDA-source handling matches the existing
  `OpenMMContext.cpp` compile flags (it is a `.cpp` compiled with the CUDA-enabled
  toolchain, not a `.cu`, so the existing glob suffices). One commit; formatting
  separate.

## Predicted test breakage (Phase-A include fixes only)
- None. No test names `GpuKinematics`, `gGpuKin`, `pushPositions`, or
  `reduceForces` as a symbol (grep of `tests/`: no hits). Tests drive the fused
  path only through the public `ForceBridge::evaluate` ->
  `OpenMMContext::pushBodyTransforms`/`reduceForcesToBodies` surface, whose
  declarations are unchanged. The SPLIT-DEDUP-FORCEREDUCER test that compares the
  host reducer against the CUDA kernel binds the kernel through
  `OpenMMContext::reduceForcesToBodies` (public), so it keeps resolving after the
  definition moves - no include fix.

## Exit criteria (machine-checked)
- Full build passes in `cuda-release` AND a non-CUDA config (`cpu-release`); test
  suite / baseline outputs identical to Stage-0 (B1/B2/B3). The DEDUP-FORCEREDUCER
  host-vs-CUDA equality test still passes against the relocated kernel.
- nm diff on public symbols vs B0: empty (only private/anonymous device symbols
  and `OpenMMContext` method definitions move; declarations unchanged).
- include-cycle script clean vs B4; clang-tidy / clang-format / IWYU clean.
- `GpuKinematics.hpp` <= 300 LOC; `GpuKinematics.cpp` <= 600 LOC (actual ~330).
  `OpenMMContext.cpp` shrinks by ~230 LOC and no longer includes the CUDA-internal
  OpenMM headers.
