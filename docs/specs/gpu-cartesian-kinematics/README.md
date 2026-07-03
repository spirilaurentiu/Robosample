# GPU Cartesian Kinematics + ABA Parallelization — Optimization Spec

Status: **01 + 02 implemented & validated** (fused position push and device force reduction — opt-in via
`ROBO_CUDA_KINEMATICS=1`, verified end-to-end on run_3SN6.lig.py; the integrator corrector hoist also
landed). **03 (ABA) is upgraded from sketch to grounded design** using the now-proven kernel
infrastructure and the constants taxonomy the stale-station bug validated; not yet implemented. Target:
`cuda` build (`USE_CUDA=ON`), with a mandatory CPU-only fallback (`USE_CPU`/`USE_REFERENCE`).

This spec plans a performance optimization of the two per-step hot loops that dominate Robosample's
inner sampling loop and of the articulated-body (ABA) recursions that surround them. It is written
against the current code, mapped in [Appendix A](#appendix-a--current-data-flow-as-mapped). It is an
*optimization* spec: it defines **what must stay correct** while the hot paths move to the GPU and
gain parallelism, and it declares — up front and explicitly — where numerical results are allowed to
change and by how much.

## 0. Scope and non-goals

**In scope**

1. Move the per-step Cartesian coordinate transform `posG[a] = X_GB[body(a)] * atomStation_B[a]` off
   the CPU and onto the GPU, writing the result **directly into OpenMM's device `posq` buffer**,
   eliminating the host `posCache_` copy and the host→device `setPositions` upload. See
   [`01-fused-position-push.md`](01-fused-position-push.md).
2. Move the per-step reduction of OpenMM per-atom forces into per-body spatial forces
   (`bodyForceG`) onto the GPU, so only `O(numBodies)` spatial forces cross PCIe per step instead of
   `O(numAtoms)` per-atom forces. See [`02-force-reduction.md`](02-force-reduction.md).
3. Parallelize the ABA recursions (`realizePosition`, `realizeVelocity`,
   `realizeArticulatedBodyInertias`, `calcUDot`, and the matrix-free operators) using the existing
   `bodyLevel` wavefront and the forest structure — portable OpenMP first, optional CUDA second. See
   [`03-aba-parallelization.md`](03-aba-parallelization.md).

**Non-goals**

- No change to the physics, the sampling algorithm, the acceptance test, the Fixman correction, or
  any joint kernel. The formulas in `RobotEngine.cpp` / `JointKernels.cpp` are preserved exactly
  (see §3, Numerical contract).
- No change to the OpenMM force-field evaluation itself (bonded/nonbonded kernels, PME, GBSA).
- No PGO/BOLT work; that is a separate opt-in facility (`nox -s build_optimized`).
- This spec does not port the ABA solve fully onto the GPU as a fused device-resident pipeline (that
  would let `bodyForceG` and `X_GB` never leave the device). That is called out as a **future phase**
  in [`03`](03-aba-parallelization.md) but is not specified here.

## 1. Why these targets

Per [Appendix A](#appendix-a--current-data-flow-as-mapped), a production run is
`rounds × replicas × worlds × mdSteps` Verlet steps — **hundreds of thousands to millions** — and each
Verlet step runs the force/derivative evaluation `evalDerivs` **1–11 times** (implicit-trapezoid
corrector, up to 10 sweeps). Every `evalDerivs` currently does, on the CPU and serially:

- a per-atom transform loop `posG[a] = X_GB[b].p() + X_GB[b].R()*station_B[a]` (`RobotEngine.cpp:606`),
- a per-atom `OpenMM::Vec3` repack into `posCache_` (`ForceBridge.hpp:37`),
- a host→device `Context::setPositions` upload (repacks again, splits to mixed precision, uploads),
- the OpenMM force eval (already on device),
- a device→host force download + a per-atom scatter into `bodyForceG` (`ForceBridge.hpp:56`).

For the target systems (**up to 1M atoms in up to 100k rigid bodies**) the per-atom traffic and the
per-atom serial CPU loops are the bottleneck around an already-GPU force eval. The transform is
embarrassingly parallel and its only per-step input, `X_GB`, is `O(numBodies)` — an order of magnitude
smaller than the `O(numAtoms)` positions it produces. This is the structural asymmetry the plan
exploits: **push `numBodies` transforms up, compute `numAtoms` positions on-device, pull `numBodies`
spatial forces down.**

## 2. The per-world invariant this plan is built on

Within the lifetime of one `World` (one Gibbs flexibility regimen, one `RobotModel`), the following are
**constant across all `mdSteps` steps** and can be uploaded to the device **once** at world setup /
first use (see `RobotModel.hpp`; confirmed constant-per-world in Appendix A):

- `atomStation_B[a]` — atom station in its body frame (`numAtoms` × `Vec3`).
- `atomBody[a]` — atom → body map (`numAtoms` × `int`).
- `atomMass[a]` (or a precomputed `isVirtual[a] = (atomMass[a]==0)` mask) — used to skip virtual sites.
- The body-sorted atom CSR (`bodyAtomsBeg/End`, `bodyAtoms`) — used by the force reduction.
- Topology for the ABA: `bodyParent`, `bodyLevel`, children CSR, `bodyQIndex/UIndex`, static frames.

The following **change every step** and are the only things uploaded per step:

- `X_GB[b]` — ground-to-body transforms (`numBodies` × `Transform`; 12 doubles each).

The following are the only things downloaded per step:

- `bodyForceG[b]` — per-body spatial force (`numBodies` × `SpatialVec`; 6 doubles each), plus the scalar
  potential energy already returned by OpenMM.

NOTE: This "upload constants once per world, stream `X_GB` per step" contract is the central design
requirement. Any implementation that re-uploads `atomStation_B`/`atomBody` per step has missed the point.

## 3. Numerical contract (read before touching any kernel)

The `optimizer` agent's default is bitwise-identical results. **This optimization cannot be bitwise**,
and the deviation is authorized and bounded here rather than discovered later.

**What is preserved exactly (NORMATIVE):**

- The mathematical formula at every site is unchanged: positions are
  `X_GB[b].p() + X_GB[b].R()*station_B[a]`; the body force is
  `BF[b].linear += f`, `BF[b].angular += (posG[a]-X_GB[b].p()) % f` summed over the body's atoms;
  virtual sites (`atomMass==0`) are excluded from the force reduction exactly as today
  (`ForceBridge.hpp:83`).
- The mixed-precision split written into `posq`/`posqCorrection` is **bit-for-bit the same operation
  OpenMM's own `setPositions` performs** (`CommonKernels.cpp:209-218`): `posq.{xyz} = (float)p`,
  `posqCorrection.{xyz} = (float)(p - (float)p)`. Given the same `double` position `p`, the two paths
  produce identical device buffers. See [`01`](01-fused-position-push.md) §4.
- Virtual-site positions are still produced by OpenMM (`computeVirtualSites`), not by the kernel.
- Detailed balance / the sampled distribution is unaffected: no formula, mass, or acceptance term
  changes.

**What is allowed to change (RECOMMENDED tolerance, reviewer-approvable):**

- The `double` value of `posG[a]` and of `bodyForceG[b]` may differ from the current CPU result by
  floating-point **association and FMA-contraction** differences between CPU (auto-vectorized,
  `-ffp-contract=fast`) and GPU (fused-multiply-add by default), plus, for the force reduction, the
  **summation order** of atoms within a body (a parallel/segmented reduction vs. the current
  sequential accumulation). Bound: relative difference at the `double` rounding level per operation
  (`~1e-15`), which for `posq` is almost always **absorbed by the `float` truncation** already inherent
  in mixed precision; per-body force differences are at the ULP level of the summands.
- Because MD is chaotic and the sampler is trajectory-level Metropolis, per-step trajectories will
  diverge over many steps even though the **step formula and the equilibrium distribution are
  unchanged**. This is the same class of non-reproducibility OpenMM already has between precision
  modes and GPU reorderings. The invariant to preserve is the *distribution*, not the *trajectory*.

**Verification obligation (NORMATIVE, deferred):** The user has directed that we run **no tests or
profiling** in this pass. This spec therefore ships a **verification harness design** in each concern
file (a one-shot A/B check: host-`setPositions` `posq` vs. kernel `posq`; CPU `bodyForceG` vs. GPU
`bodyForceG`, on ala-dipeptide) that MUST be run before this work is declared correct whenever the test
gate is re-enabled. Not running it now is a deliberate, recorded deviation — not evidence of
correctness.

## 4. CPU-only compatibility (NORMATIVE)

The engine must still build and run correctly with `USE_CUDA=OFF` (`USE_CPU`/`USE_REFERENCE`), where
OpenMM has no CUDA context and no device `posq`. Therefore:

- Every GPU path in [`01`](01-fused-position-push.md) and [`02`](02-force-reduction.md) is guarded by
  `#if USE_CUDA` **and** a runtime check that the active OpenMM platform is CUDA. The existing host
  path (`realizePosition`'s per-atom loop → `ForceBridge::setAtomPositionsInGround` →
  `Context::setPositions`; force download → `ForceBridge::getForcesFromOpenMM`) remains the fallback
  and remains the definition of correct behavior the GPU path is checked against.
- The ABA parallelization in [`03`](03-aba-parallelization.md) uses OpenMP, which is already linked on
  both CPU and CUDA builds (`CMakeLists.txt` `find_package(OpenMP REQUIRED)`), and degrades to the
  current serial loops at thread count 1. It requires no CUDA.
- No `.so`/Python API surface changes. The optimization is entirely below `ForceBridge` / `RobotEngine`.

## 5. Rollout order (implementation phasing — not architecture)

NOTE: This section is sequencing guidance only; the correctness requirements above do not depend on it.

1. **Phase 1 — Fused position push** ([`01`](01-fused-position-push.md)). Highest, cleanest win; touches
   only the write path. Includes the once-per-world constant upload and the per-step `X_GB` stream.
2. **Phase 2 — GPU force reduction** ([`02`](02-force-reduction.md)). Symmetric to Phase 1 on the read
   path; depends on the constant buffers Phase 1 already uploads.
3. **Phase 3 — ABA CPU parallelization** ([`03`](03-aba-parallelization.md) §A). OpenMP level-wavefront
   + forest parallelism; portable, no device dependency.
4. **Phase 4 — ABA CUDA (optional/future)** ([`03`](03-aba-parallelization.md) §B). One-thread-per-body
   dense algebra on the level wavefront; only justified for the large-system regime and only after
   Phases 1–3 are validated.

Phases 1 and 2 are independent of Phase 3/4 and of each other's correctness (each has its own fallback),
so they can land and be reverted independently.

Each phase is handed to `reviewer` before the next begins (per CLAUDE.md: optimizer hands every change
to reviewer; never self-certifies).

---

## Appendix A — current data flow (as mapped)

Grounding references (file:line) for the claims above; full maps are in the exploration notes.

- Per-step driver: `RobotEngine::verletStep` (`include/RobotIntegrator.hpp:102`); `evalDerivs` lambda
  (`:290-311`) runs `realizePosition → bridge.evaluate → realizeVelocity →
  realizeArticulatedBodyInertias → calcUDot → calcQDot/calcQDotDot`; corrector re-runs it up to 10×
  (`:326-372`). Move loop: `World::generateSample` (`src/World.cpp:1550`).
- Transform site: `RobotEngine::realizePosition` (`src/RobotEngine.cpp:532`), `X_GB` build at `:559-563`,
  per-atom `posG`/`stG` at `:602-611`. Duplicate in `fillAtomPositionsFromBodies` (`:1321-1337`).
- Position handoff: `ForceBridge::setAtomPositionsInGround` (`include/ForceBridge.hpp:37`), `posCache_`
  (`:150`), `ForceBridge::evaluate` (`:112`). `Context::setPositions` on the host at
  `OpenMMContext.hpp:115` / `OpenMMContext.cpp:297,309,350`.
- Force reduction: `ForceBridge::getForcesFromOpenMM` (`include/ForceBridge.hpp:56-109`); virtual-site
  skip at `:83`; scatter `BF[b].linear += f`, `BF[b].angular += r % f` with `r = posG[a]-X_GB[b].p()`.
- OpenMM is a **vendored fork built from source** (`.gitmodules` → `github.com/spirilaurentiu/openmm`,
  built via `cmake_modules/OpenMM.cmake`); CUDA platform selected at `OpenMMContext.cpp:272-276`
  (mixed precision). No direct device-buffer access exists today.
- SoA storage in one aligned `MemoryArena` per `RobotState` (`include/MemoryArena.hpp`,
  `include/RobotState.hpp`); tree invariant `bodyParent[b] < b`, `bodyLevel` reserved for wavefronts
  (`include/RobotModel.hpp:67`, `include/RobotEngine.hpp:20-21`); forest = molecules rooted at Ground.

## Appendix B — confirmed OpenMM device internals (for `01`/`02`)

From the vendored source (`openmm/platforms/cuda`, `openmm/platforms/common`):

- `CudaContext::getPosq()` → `CudaArray` of `float4` (single/mixed) or `double4` (double), length
  `paddedNumAtoms`. `getPosqCorrection()` → `float4`, exists only when `getUseMixedPrecision()`
  (`CudaContext.h:211,217`; allocated `posqCorrection.initialize<float4>(...)` at `CudaContext.cpp:237`).
- Internal atom order: `posq` is stored in reordered order; `order = cc.getAtomIndex()` maps device
  slot `i` → original particle index `order[i]` (`CommonKernels.cpp:183,188`). Device copy available via
  `getAtomIndexArray()` (`CudaContext.h:283`).
- Mixed-precision split reproduced by the kernel (`CommonKernels.cpp:209-218`).
- `setPositions` resets `posCellOffsets` and calls `reorderAtoms()` (`CommonKernels.cpp:223-225`).
- `reorderAtoms()` is a no-op unless a cutoff/neighbor-list is in use and ≥250 steps have passed
  (`ComputeContext.cpp:519-527`): **vacuum / implicit-solvent (no cutoff) never reorder → `order[i]==i`**;
  explicit-solvent (PME) reorders periodically. This distinction drives the two mapping regimes in
  [`01`](01-fused-position-push.md) §3.
- Virtual sites: `ContextImpl::computeVirtualSites()` (`ContextImpl.cpp:298`).
