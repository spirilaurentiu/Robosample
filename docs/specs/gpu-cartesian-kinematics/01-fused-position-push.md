# 01 — Fused GPU position push (`X_GB * station` → OpenMM `posq`)

Concern: replace the CPU transform + host repack + `setPositions` upload with a single CUDA kernel that
computes each atom's ground-frame position from the per-body transform and writes it **directly into
OpenMM's device `posq`** (and `posqCorrection` in mixed precision), so no per-atom position data crosses
PCIe. Read [`README.md`](README.md) §2 (per-world invariant) and §3 (numerical contract) first.

## 1. The operation being moved

Today (`RobotEngine.cpp:602-611`, then `ForceBridge.hpp:37-43`, then `OpenMMContext`/`CommonKernels.cpp:181`):

```
CPU:  for a: posG[a] = X_GB[b].p() + X_GB[b].R() * atomStation_B[a]   // b = atomBody[a]
CPU:  for a: posCache_[a] = OpenMM::Vec3(posG[a])                     // repack
H→D:  Context::setPositions(posCache_)                                // repack again + split + upload
```

Replace with:

```
H→D (per step):  upload X_GB[0..numBodies)                            // O(numBodies) transforms
GPU (per step):  for a in parallel: p = X_GB[atomBody[a]] * station_B[a]; write posq[slot(a)] (+corr)
GPU (per step):  computeVirtualSites()                               // OpenMM, unchanged
```

`atomStation_B`, `atomBody`, `isVirtual`, and the reorder mapping are device-resident constants uploaded
once per world (§2).

## 2. Device-resident state (NORMATIVE layout)

Owned by a new device-side helper (proposed: `CudaKinematicsBridge`, living beside `ForceBridge`, only
compiled under `USE_CUDA`). All allocated once per world, freed at world teardown.

**Constant for the world (uploaded once):**

| buffer | type | length | source |
|---|---|---|---|
| `d_station` | `double3` (or `Real3`) | `numAtoms` | `model.atomStation_B` |
| `d_atomBody` | `int` | `numAtoms` | `model.atomBody` |
| `d_isVirtual` | `uint8` | `numAtoms` | `model.atomMass[a]==0` |

**Streamed per step (uploaded each `evalDerivs`):**

| buffer | type | length | source |
|---|---|---|---|
| `d_XGB` | 12×`double` per body (R row-major + p) | `numBodies` | `state.X_GB()` |

NOTE: `X_GB` is a `robo::Transform` (a `Rotation`+`Vec3`); upload it in a packed POD form (`double[12]`)
via one `cudaMemcpyAsync` from the SoA arena. The arena is already contiguous and 64-byte aligned
(`MemoryArena.hpp`), so the source is a single flat span — no gather.

**Precision:** compute the transform in `double` on the device regardless of OpenMM's precision mode, to
match the CPU formula as closely as the hardware allows, then narrow at the write (§4). Uploading `X_GB`
as `double` is required for this; the per-step upload is still only `numBodies × 96 B`
(≈9.6 MB at 100k bodies, vs ≈24 MB of `double` positions at 1M atoms).

## 3. The atom→`posq`-slot mapping (the correctness crux)

`posq` is in OpenMM's internal reordered order; slot `i` holds particle `order[i] = getAtomIndex()[i]`.
The kernel must write atom `a`'s position into the slot `i` with `order[i] == a`, i.e. it needs the
**inverse** map `invOrder[a] = i`. Two regimes (per README Appendix B):

**Regime A — no atom reordering (vacuum / implicit solvent, no cutoff).** `order[i] == i` permanently;
`invOrder` is identity; atom `a` writes `posq[a]`. This covers ala-dipeptide and implicit-solvent FFAR1.

**Regime B — periodic reordering (explicit solvent, PME/cutoff).** `order` changes roughly every 250
steps inside OpenMM. Requirements (NORMATIVE):

- The kernel SHALL scatter using the **current** device `atomIndex` array (`getAtomIndexArray()`), read
  as the inverse map, refreshed whenever OpenMM reorders.
- Detecting a reorder: `ComputeContext::getAtomsWereReordered()` / the reorder happens inside
  `reorderAtoms()`. Because we now bypass `setPositions` (which is what normally *triggers* the
  post-set reorder), the driver of reordering becomes the force/neighbor-list path. The implementation
  SHALL, once per step, refresh a cached device `invOrder` from `atomIndex` when
  `getAtomsWereReordered()` is set (cheap: one `int` scatter kernel over `numAtoms`), else reuse the
  cached map.
- SHALL preserve the semantics OpenMM's own `setPositions` has around reordering: it zeroes
  `posCellOffsets` before reordering (`CommonKernels.cpp:223`). Since periodic-image cell offsets are
  only meaningful relative to a set of positions, writing a fresh full position set each step and
  letting the next force eval/reorder recompute offsets is consistent — but this interaction MUST be
  verified in the Regime B A/B check (§6) before explicit-solvent use.

RECOMMENDED: land Regime A first (it is the ala-dipeptide/implicit path and needs no `invOrder` logic),
gate Regime B behind the same runtime platform check plus a "hasCutoff" query, and fall back to the host
`setPositions` path for Regime B until its A/B check passes. This keeps explicit-solvent correctness on
the proven path while the vacuum/implicit path — the common case for the robotic worlds — gets the win.

## 4. The write (mixed-precision split — bit-exact with OpenMM)

For each real atom `a` (skip if `isVirtual[a]`; those are filled by `computeVirtualSites`), with
`b = atomBody[a]`, compute in `double`:

```
p.x = XGB[b].p.x + dot(XGB[b].R.row0, station[a])
p.y = XGB[b].p.y + dot(XGB[b].R.row1, station[a])
p.z = XGB[b].p.z + dot(XGB[b].R.row2, station[a])
i   = invOrder[a]                 // = a in Regime A
```

Then write exactly as OpenMM does (`CommonKernels.cpp:185-221`):

- **Double precision** (`getUseDoublePrecision()`): `posq[i].{x,y,z} = p.{x,y,z}` (leave `.w` charge
  untouched — see below).
- **Single precision** (neither double nor mixed): `posq[i].{x,y,z} = (float)p.{x,y,z}`.
- **Mixed precision** (default, `getUseMixedPrecision()`): `posq[i].{x,y,z} = (float)p.{x,y,z}` **and**
  `posqCorrection[i].{x,y,z} = (float)(p.{x,y,z} - (float)p.{x,y,z})`, `posqCorrection[i].w = 0`.

NORMATIVE: the kernel SHALL NOT touch `posq[i].w` (the per-atom charge, packed into `posq` by OpenMM at
context init, `hasAssignedPosqCharges`). Only `.xyz` are position. Write `.xyz` in place.

NOTE: Given the same `double p`, this write is byte-identical to what `setPositions` produces; the only
admissible difference vs. the current path is in `p` itself (CPU vs GPU `double` rounding of
`R*station + t`), per README §3. This is why the split is specified down to the cast.

## 5. Driver integration

- `ForceBridge::evaluate(s)` (`ForceBridge.hpp:112`) currently does
  `setAtomPositionsInGround(s); getForcesFromOpenMM(s);`. Under the CUDA fused path it becomes:
  `kin.pushPositionsToDevice(s.X_GB()); omm.computeVirtualSitesOnDevice(); getForces...` — no host
  `posCache_`, no `Context::setPositions`.
- The `X_GB` the kernel consumes is produced by `realizePosition`'s **body sweep** (`RobotEngine.cpp:559-563`),
  which stays on the CPU in Phase 1. Only the **per-atom loop** (`:602-611`) is what the kernel replaces.
  RECOMMENDED: keep computing `stG[a]` (the rotated station) — but see [`02`](02-force-reduction.md),
  which moves both the per-atom position loop's *consumer* (force reduction) to the device and lets
  `stG` be recomputed on-device from the same `X_GB` + `station`, so the CPU per-atom loop
  (`:602-611`) can be **removed entirely** on the CUDA path once Phases 1+2 are both in.
- `fillAtomPositionsFromBodies` (`RobotEngine.cpp:1321`, the accept path / DCD / saved positions) still
  needs host-side `posG` occasionally (rollback snapshots, trajectory output). This is **not** in the
  per-step hot loop (it runs at accept/reject and DCD cadence, not every corrector sweep). Leave it on
  the CPU. Only pull `posG` back to the host when actually needed (accept, DCD write), not every step.

NOTE: The `computeVirtualSites` call is required each step in any world that has virtual sites (e.g.
4-site water). In vacuum/implicit peptide worlds there are none and the call is a cheap no-op; keep it
unconditional to match OpenMM's `setPositions` contract rather than branching on a per-world flag.

## 6. Verification harness (deferred per README §3, but designed now)

A one-shot A/B equivalence check to run when the test gate is re-enabled:

- **Regime A (ala-dipeptide, vacuum):** for a fixed `X_GB`, run (a) the current host path
  (`setAtomPositionsInGround` + `Context::setPositions`) and (b) the fused kernel, each into a fresh
  context; download `posq`(+`posqCorrection`) from both and assert **bitwise equality** of the device
  buffers *when the input `posG` is forced identical* (feed the kernel the CPU-computed `posG` to
  isolate the write from the transform), then separately assert the *transform* agreement
  `max_a |p_gpu[a] - p_cpu[a]| < 1e-12` (double-rounding bound). Wire as a `ctest` case
  `TestFusedPositionPush` under the existing test harness.
- **Regime B (deca-alanine explicit solvent):** additionally step OpenMM enough to force at least one
  `reorderAtoms()`, then assert positions read back via `getState(Positions)` match the host path to the
  same `1e-12` bound, confirming the `invOrder`/cell-offset handling.

This harness is the gate that converts the declared tolerance (README §3) into evidence. It is **not**
run in this pass by the user's instruction; it is specified so that "not yet verified" is explicit.
