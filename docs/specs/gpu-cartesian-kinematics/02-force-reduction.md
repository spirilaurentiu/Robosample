# 02 — GPU per-atom force → per-body spatial force reduction

Concern: reduce OpenMM's per-atom Cartesian forces into per-body spatial forces (`bodyForceG`) **on the
device**, so the read path is symmetric to [`01`](01-fused-position-push.md): only `O(numBodies)` spatial
forces cross PCIe per step, not `O(numAtoms)` per-atom forces. Read [`README.md`](README.md) §3 first.

## 1. The operation being moved

Today (`ForceBridge::getForcesFromOpenMM`, `ForceBridge.hpp:56-109`): OpenMM forces are downloaded to the
host, then scattered per atom into per-body accumulators:

```
D→H:  forces = getState(Forces)                                  // O(numAtoms) download
CPU:  for b: BF[b] = 0;  for i in nu: mob[i] = 0
CPU:  for a:                                                     // scatter, sequential
        if atomMass[a]==0: continue                              // skip virtual sites
        b = atomBody[a];  f = forces[a]
        r = posG[a] - X_GB[b].p()                                // = stG[a] = X_GB[b].R()*station[a]
        BF[b].linear  += f
        BF[b].angular += r % f                                   // torque about body origin Bo, in Ground
```

Replace with a device segmented reduction over each body's atoms, producing `bodyForceG` on-device, then
download only `bodyForceG`:

```
GPU:  for b in parallel:                                          // one block/warp per body (CSR segment)
        acc_lin = 0; acc_ang = 0
        for a in [bodyAtomsBeg[b], bodyAtomsEnd[b]):
          if isVirtual[a]: continue
          f = deviceForce(a)                                      // read OpenMM device force buffer
          r = XGB[b].R() * station[a]                             // recompute stG on-device (= r)
          acc_lin += f;  acc_ang += r % f
        BF[b] = {angular: acc_ang, linear: acc_lin}
D→H:  download bodyForceG[0..numBodies)                           // O(numBodies) download
```

`mobilityForce[i]` is set to 0 here today (Fixman/bias torques would be added later); keep that zeroing
on the host — it is `O(nu)` and unrelated to the per-atom scatter.

## 2. Reading OpenMM's device force buffer (NORMATIVE)

OpenMM's CUDA force buffer is **64-bit fixed point**, not float: `CudaContext::getForce()` /
`getLongForceBuffer()` is a `long long` array of length `3*paddedNumAtoms`, laid out component-major
(`force[i]`, `force[i+padded]`, `force[i+2*padded]`), and converted to `double` newtons by
`scale = 1/(double)0x100000000` in `CommonKernels.cpp:289-295`:

```
forces[order[i]] = scale * Vec3(force[i], force[i+padded], force[i+2*padded])
```

Requirements:

- The reduction kernel SHALL read this fixed-point buffer and apply the **same** `scale` and the **same**
  `order`/`atomIndex` mapping used in `01` (device slot `i` ↔ particle `order[i]`), so that `deviceForce(a)`
  returns exactly what the host path's `forces[a]` would have been. Recompute nothing about the force
  itself — only its container changes.
- The kernel SHALL run **after** the force computation and any force reordering are complete for the step
  (i.e. after `ContextImpl::calcForcesAndEnergy` finishes), reading the finalized force buffer. It must
  not race the nonbonded/PME kernels.
- Virtual-site atoms are skipped exactly as today (`isVirtual[a]`, README §2). OpenMM redistributes
  virtual-site forces onto their parent real atoms as part of force computation, so by read time the real
  atoms already carry the full force — no change from current behavior.

## 3. The reduction (segmented, per body)

Use the body-sorted atom CSR (`bodyAtomsBeg/End`, `bodyAtoms`) already in `RobotModel` (README §2) so each
body owns a contiguous atom segment. RECOMMENDED mapping: **one warp (or block) per body**, warp-stride
over the segment, warp-reduce the six accumulators (`acc_lin.xyz`, `acc_ang.xyz`) in `double`, lane 0
writes `BF[b]`. This avoids global atomics entirely and gives deterministic-per-launch results.

NOTE on load balance: body atom counts are highly non-uniform (a solvent molecule = a few atoms; a large
rigid protein core = thousands). A fixed warp-per-body wastes lanes on tiny bodies and serializes huge
ones. RECOMMENDED for the large-system regime: a two-tier scheme — small bodies handled one-per-lane or
one-per-warp, large bodies (segment length above a threshold) handled by a whole block with a block-wide
reduction — selected from a per-world histogram of segment lengths computed once at upload. This is an
optimization detail, not a correctness requirement; a plain warp-per-body is correct and is the starting
point.

NORMATIVE: whatever the schedule, the six per-body sums SHALL be accumulated in `double`. The summation
*order* may differ from the host loop (that is the authorized tolerance, README §3), but the accumulator
*type* must not narrow to `float`.

## 4. `r` on-device, and retiring the CPU per-atom loop

`r = posG[a] - X_GB[b].p() = X_GB[b].R() * station[a]` — the rotated station, which `01`'s kernel already
computes as `stG`. Two options:

- **(a) Recompute** `r` inside this kernel from `d_XGB` + `d_station` (same inputs as `01`). Cheapest in
  memory; one extra mat-vec per atom.
- **(b) Cache** `stG[a]` from `01`'s position kernel into a device buffer and read it here.

RECOMMENDED: (a). It keeps the two kernels independent (either can run without the other's scratch) and
the mat-vec is trivial next to a force read. Either way, once **both** `01` and `02` are in on the CUDA
path, the entire host per-atom loop `RobotEngine.cpp:602-611` (which produced `posG`/`stG` for exactly
these two consumers) is dead on that path and SHALL be compiled out — leaving `X_GB` as the sole per-step
host→device payload and `bodyForceG` as the sole device→host payload.

## 5. Driver integration & fallback

- New device helper method (proposed `CudaKinematicsBridge::reduceForcesToBodies(...)` alongside `01`'s
  push), invoked from `ForceBridge::getForcesFromOpenMM` under `#if USE_CUDA` + runtime CUDA-platform
  check. Result `bodyForceG` is `cudaMemcpy`'d into the `RobotState` arena's `bodyForceG()` span.
- Fallback (non-CUDA or non-CUDA platform): the existing host `getForcesFromOpenMM` loop
  (`ForceBridge.hpp:56-109`) is unchanged and remains the correctness reference.
- `wantsAtomForces()` (`ForceBridge.hpp:82`): when some caller genuinely needs per-atom forces on the host
  (diagnostics, the reaction-force CSV campaign), keep a path that also downloads the per-atom force
  buffer. The default hot path does **not** download per-atom forces.

## 6. Verification harness (deferred per README §3, designed now)

`ctest` case `TestGpuForceReduction` (ala-dipeptide): after one OpenMM force eval, compute `bodyForceG`
both ways (host scatter vs. device segmented reduction) from the identical force buffer and assert
`max_b ||BF_gpu[b] - BF_cpu[b]|| < tol`, with `tol` at the ULP-scale of the per-body summand magnitudes
(a few×`1e-9` of the largest per-body force is a reasonable bound; the difference is pure summation-order
rounding, README §3). Include a body with many atoms (a whole rigid molecule) to exercise the reduction
tail, and a body with one atom to exercise the trivial segment. Not run in this pass; specified so the
tolerance is checkable later.
