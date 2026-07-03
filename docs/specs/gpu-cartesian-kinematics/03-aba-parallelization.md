# 03 — ABA recursion parallelization & residency

Concern: make the articulated-body (ABA) recursions faster **without changing their math** — first by
removing per-step redundancy the corrector exposes, then by exploiting the structured parallelism inside
each pass (forest + `bodyLevel` wavefront), and ultimately by keeping the whole corrector loop resident on
the device. Read [`README.md`](README.md) §3 (numerical contract) and §4 (CPU compatibility) first.

## 0. Status & what the implemented fused path (01/02) established

Phases **01 (fused position push)** and **02 (device force reduction)** are **implemented and validated**
(`src/OpenMMContext.cpp` kernels K1/K2, `include/ForceBridge.hpp` dispatch). Several facts learned there are
load-bearing for this spec and are stated once here:

- **The common-compute kernel infrastructure works and is the template for a device ABA.**
  `ComputeContext::compileProgram(src, defines)` → `createKernel` → `addArg`/`execute`, `ComputeArray`
  (`initialize<T>` / `upload` / `download`), and `ContextSelector` (run on OpenMM's stream) are all proven
  in K1/K2. A CUDA ABA (§B) reuses exactly this — no new build machinery, portable to CUDA+OpenCL.
- **The per-step host↔device budget is already `O(numBodies)`.** Per robotics step the fused path moves
  `X_GB` up (12·numBodies doubles) and `bodyForceG` down (6·numBodies doubles); atom stations (3·numAtoms)
  move up only once per round. A device ABA closes this to ~zero per step (§B).
- **The constants taxonomy is validated (the stale-station bug).** RobotModel data splits three ways, and
  conflating them corrupts results:
  - **Per-world constant** (built once in `buildModel`): topology — `bodyParent`, `bodyLevel`, children
    CSR, body-atom CSR, `bodyQIndex/UIndex`, `bodyJoint`, `atomBody`, `atomMass`.
  - **Per-transfer / per-round** (refit every coordinate handoff in `World::recomputeGeometry`): the body
    frames and mass properties — `X_PF`, `X_BM`, `bodyCom_B`, `bodyUnitInertia_B`, `bodyMass`,
    `atomStation_B`. **These are NOT world-lifetime constants** (this is what broke 01/02 until stations
    were re-uploaded per round — see README/campaign notes). A device ABA MUST re-upload them each round.
  - **Per-step varying**: `q`, `u`, `X_GB`, `X_FM`, `X_PB`, `Phi`, `Mk_G`, `V_*`, `A_GB`, `P`, `PPlus`,
    `D`/`DI`, `G`, `Z`, `zPlus`, `bodyForceG`, `mobilityForce`.
- **The ABA lives in body/SoA indexing, not OpenMM atom order.** K1/K2 needed an `invOrder` map because
  OpenMM spatially reorders its `posq`/force buffers; the ABA operates on the `RobotState` arena indexed by
  body `b` and generalized-speed offset `uOff`, which OpenMM never touches. So **the reorder/`invOrder`
  complexity does NOT apply to the ABA arrays** — only to the K1/K2 atom↔posq boundary. A device ABA is
  simpler here.

This file has three parts: **§0.5** a bitwise redundancy fix the corrector hoist exposed (do first — no
parallelism), **§A** portable OpenMP wavefront (Phase 3), **§B** a device-resident CUDA ABA (Phase 4).

## 0.5 The inertia-factorization hoist (bitwise, no parallelism, highest ROI) — IMPLEMENTED

Status: **landed & verified** (2026-07-03). `realizeArticulatedBodyInertias` split into
`factorizeArticulatedInertias` (P/PPlus/D/DI/G) + `seedArticulatedCentrifugal` (abcf), with a
backward-compatible wrapper for non-corrector callers; `verletStep` calls the factorization once in
`evalPos` and only the seed per `evalVel` sweep. Verified bitwise-identical (fused and host reproduce the
prior PE/KE/Fixman/H to every digit across rounds).

The velocity-Verlet driver (`RobotEngine::verletStep`, `include/RobotIntegrator.hpp`) now runs, per step:

- **`evalPos` once** — `realizePosition` + the OpenMM force eval (`bridge.evaluate`). Hoisted out of the
  corrector because both are pure functions of `q`, which the corrector never changes (landed already).
- **`evalVel` up to ~11×** — `realizeVelocity` + `realizeArticulatedBodyInertias` + `calcUDot` + `calcQDot`/
  `calcQDotDot`, iterated by the implicit-trapezoid corrector (converges `u`).

So **the corrector loop is now the hot ABA path**, and it re-runs `realizeArticulatedBodyInertias` every
sweep. But that pass computes two separable things:

- **The articulated-inertia factorization** — `P`, `PPlus`, `D`, `DI` (the per-body Jacobi eigensolve /
  `invertDense`), `G` (Kalman gain). These depend only on `X_GB`/`Phi`/`Mk_G` (i.e. on `q` and masses) —
  **NOT on `u`.** Constant across all corrector sweeps.
- **The articulated centrifugal seed** — the trailing loop `abcf[b] = P[b]·a_mob[b] + gyro[b]`. `a_mob`
  and `gyro` come from `realizeVelocity`, so this **is** velocity-dependent and must stay per-sweep.

**NORMATIVE opportunity:** split `realizeArticulatedBodyInertias` into `factorizeArticulatedInertias`
(P/PPlus/D/DI/G — position-only) and `seedArticulatedCentrifugal` (`abcf` — velocity-dependent). Call the
factorization once in `evalPos`; keep only the cheap `abcf` seed (and `calcUDot`) in `evalVel`. This cuts
the per-body Jacobi eigensolves from ~11×/step to 1×/step and is **bitwise-identical** (same `q` ⇒ same
P/D/G every sweep, so computing them once reproduces the per-sweep values exactly — the same argument that
justified the position hoist). NOTE: preserve the two load-bearing details when re-partitioning — the
child-sum uses the child's *outboard* quantity shifted by the child's own `Phi`
(`PPlus[c].shift(Phi[c].l())`), and `abcf` is seeded from the **mobilizer** Coriolis `a_mob`, not total
`a_tot`. This is the single highest-ROI ABA change and needs no parallelism; §A/§B build on top of it.

## 1. The dependency structure (what may and may not run concurrently)

From `RobotEngine.cpp`, with invariant `bodyParent[b] < b`:

- **Outward passes** — `realizePosition`, `realizeVelocity`, `calcUDot` pass 2, `multiplyBySqrtMInv`: body
  `b` reads its **parent**'s just-written result (`X_GB[b]=X_GB[p]·X_PB[b]`; `V_GB[b]=~Phi[b]·V_GB[p]+…`;
  `APlus=~Phi[b]·A_GB[p]`). A body may run only after its parent.
- **Inward passes** — the inertia factorization (§0.5), `calcUDot` pass 1, `multiplyByMInv`,
  `calcMobilizerReactionForces`: body `b` reads **all its children**'s results (`Pb += PPlus[c].shift(…)`;
  `z += Phi[c]·ZPlus[c]`). A body may run only after all its children.
- **Independent bodies:** (a) **different molecules are disjoint subtrees rooted at Ground** — fully
  independent for the whole pass; (b) **all bodies at the same `bodyLevel` are mutually independent** — the
  wavefront `bodyLevel` is reserved for (`RobotModel.hpp:67`, `RobotEngine.hpp:20-21`).
- **Fully independent (no inter-body dependency at all):** the per-atom transform loop (now SIMD, and on
  the fused path done by K1), the force reduction (K2 / host CSR loop), `calcKineticEnergy` (a reduction),
  `calcLogDetM` (independent per-body determinant + sum — `det M = ∏_b det D_b`, Jain 1997), the `abcf`
  seed, `calcQDot`/`calcQDotDot` (per-body `q`-blocks). These parallelize trivially.

NORMATIVE: any reordering SHALL preserve the two load-bearing details above (child outboard-shift; `a_mob`
seed).

## 2. §A — OpenMP level-wavefront + forest parallelism (Phase 3, portable)

Prioritize the **corrector-loop passes** (`realizeVelocity`, the `abcf` seed, `calcUDot`) — post-hoist those
are what run ~11×/step; `realizePosition` runs ~2×/step and the inertia factorization ~1×/step (§0.5).

### 2.1 Precompute the schedule once per world (NORMATIVE data, no math change)

At `World::buildModel`, alongside `bodyLevel`, build a **level CSR**: `levelBeg[L]`, `levelEnd[L]` into a
`bodyByLevel[]` array grouping bodies by level in increasing `b` order. Bodies are already relabeled into
level-major topological order (`World.cpp`), so `bodyByLevel` is likely the identity permutation or a
trivial slice — verify and, if so, iterate the existing arrays over `[levelBeg[L], levelEnd[L])` with no
extra indirection. Run-constant (per-world); compute once.

### 2.2 Parallelize each pass over the wavefront

Replace each outward `for (b=1..numBodies)` / inward `for (b=numBodies-1..1)` with a level loop that is
**serial across levels** and **parallel within a level**:

```
for (int L = 1; L <= maxLevel; ++L)                       // outward (reverse for inward)
    #pragma omp parallel for schedule(static)             // bodies at level L are independent
    for (int k = levelBeg[L]; k < levelEnd[L]; ++k) { int b = bodyByLevel[k]; /* body-b work */ }
```

Correctness rests on the wavefront property: within level `L`, every body's parent is at level `<L`
(finalized, outward) and every child at level `>L` (finalized, inward). Siblings at level `L` write disjoint
`b`, so no atomics/locks — the SoA is indexed by `b` and each task owns its `b`.

### 2.3 The fully-independent loops

`calcLogDetM`, the `abcf` seed, `calcQDot`, `calcQDotDot`, and the per-body dense algebra inside the inertia
factorization (`D=~H P H`, the Jacobi eigensolve, `symSqrt`/`symSqrtInv`) are per-body-local: a plain
`#pragma omp parallel for` over all bodies is correct. `calcKineticEnergy` / `calcLogDetM`'s scalar sum use
`reduction(+:…)` (or the serial-sum variant, §2.5 NOTE).

### 2.4 Cost model and gating (RECOMMENDED)

- Speedup ceiling is `numBodies / (critical-path length)`; the critical path is **tree depth** (`maxLevel`),
  not `numBodies`. Bushy systems (solvent/membrane/the 100k-body regime) parallelize well; a long unbranched
  backbone does not. Wavefront width = bodies per level is the available parallelism.
- OpenMP fork/join costs ~µs per region; for tiny worlds (`numBodies` ~ tens) the overhead dominates.
  RECOMMENDED: gate on a per-world threshold (`numBodies >= ROBO_ABA_OMP_MIN`, e.g. a few hundred, computed
  at `buildModel`); below it, the serial loop (via `if (parallel)` or `num_threads(1)`), one code path.
- The ABA and the OpenMM force eval run in strict sequence (never concurrently) so they share the full core
  count; respect `OMP_NUM_THREADS`, do not oversubscribe.
- NOTE: `-fopenmp-simd` is already enabled for Robosample's code (added for the transform loops); full
  OpenMP threading additionally needs `-fopenmp` linked to `robosample_objects` (a one-line CMake change) —
  `OpenMP::OpenMP_CXX` is already found/linked for the final targets.

### 2.5 Numerical effect

`schedule(static)` at a fixed thread count is deterministic, and the per-body child-sum order is **fixed by
the CSR** (`bodyChildren` order) — each body sums its own children sequentially inside its task — so the
inward passes are **bitwise-identical** to serial, and the outward passes (no cross-body reduction) also.
The only order-sensitive pieces are the global scalars (`calcKineticEnergy`, `calcLogDetM` sum) under
`reduction(+:)`; declare that under README §3 (bounded, ~`1e-15` rel), or use the serial-sum variant (NOTE)
to keep them bitwise too. **This is a stronger guarantee than the fused K1/K2 path**, whose float `posq` and
parallel force reduction admit ULP differences.

## 3. §B — device-resident CUDA ABA (Phase 4)

Now a grounded design (not a sketch), because 01/02 proved the infrastructure and the constants taxonomy.

### 3.1 Motivation — close the loop (the real prize)

Post-hoist, per step the fused path already keeps `posq` and forces on the device; the remaining host↔device
traffic is `X_GB` up and `bodyForceG` down, and the corrector (`evalVel`) runs ~11× on the **CPU** using the
once-downloaded `bodyForceG`. A device ABA lets:

- `bodyForceG` (already produced on-device by K2) feed a device `calcUDot` **without the download**,
- the inertia factorization + `realizeVelocity` + `calcUDot` corrector iterate **entirely on device**,
- `X_GB` be produced on-device by a device `realizePosition` (forward kinematics from `q`),

so that **nothing per-step crosses PCIe except at accept/DCD cadence.** The corrector's ~11 iterations then
cost no host round trips at all.

### 3.2 Shape

- **Kernels via the proven path** (`ComputeContext::compileProgram`, common-compute language, one
  `ComputeArray` per per-step slab). One thread (or small tile) per body; the per-body dense algebra
  (dof ≤ 6: `D=~H P H`, the `DI` inverse, `symSqrt`) is a natural per-thread workload and needs no `invOrder`
  (§0 — ABA arrays are in body indexing).
- **Schedule = the level CSR (§2.1).** One kernel launch per level (serial across levels, parallel within),
  outward and inward. Launch latency (~µs) × `2·maxLevel` per pass is the cost risk; for deep trees a single
  persistent kernel with grid-sync (cooperative groups) that walks levels in-kernel avoids per-level
  relaunch. Evaluate both.
- **Device state.** A device mirror of the per-step `RobotState` slabs (`q, u, X_GB, X_FM, X_PB, Phi, Mk_G,
  V_*, A_GB, P, PPlus, H, G, DI, Z, zPlus, bodyForceG, mobilityForce`) plus the topology (uploaded once per
  world) and the refit frames/masses (`X_PF, X_BM, bodyCom_B, bodyUnitInertia_B, bodyMass, atomStation_B` —
  **re-uploaded per round**, per the §0 taxonomy; getting this wrong is the stale-station failure mode). The
  64-byte-aligned `MemoryArena` maps cleanly to coalesced device arrays.
- **Boundary with K1/K2.** `bodyForceG` is the hand-off: K2 writes it on-device; a device `calcUDot` reads
  it in-place. A device `realizePosition` writes `X_GB` on-device, which K1 then consumes to write `posq`.
  With both, the per-step `X_GB`-up / `bodyForceG`-down transfers disappear.

### 3.3 Alternative for depth-bound trees (NOTE, not specified)

The level-wavefront span is the **tree depth**. For a single long macromolecular backbone (few branches,
deep chain) this is poor. The Divide-and-Conquer Algorithm (DCA; Featherstone 1999 / Anderson) reformulates
forward dynamics as an `O(log N)`-span binary assembly of subtrees — the canonical fix. It is a genuine
re-derivation of the dynamics (higher constant, needs a `researcher` spec + validation against the reference
recursion via the existing disasm-vs-refactor differential oracle, `[[robotics-oracle-campaign]]`). Only
pursue it if profiling shows the ABA is depth-bound on a real target system; the wavefront (§A/§3.2) already
saturates the GPU for the bushy solvent/membrane regime.

### 3.4 Numerical

Full ABA in device `double` differs from the CPU by FMA/association throughout — distribution-preserving,
trajectory-diverging (README §3), and the largest departure of the campaign. It needs its own A/B campaign
(device `udot`/`A_GB` vs host for a fixed `RobotState`) before use. Reusing the SKO/innovations factorization
(the current matrix-free `multiplyByMInv` / `calcLogDetM` / `multiplyBySqrtM` sweeps, `[[jain_1997]]`) means
`M⁻¹`, `√M` (momentum sampling), and `log det M` (Fixman) all come from the same device factorization — do
NOT form `M` and call a dense solver; the ABA *is* the sparse Cholesky, done matrix-free.

## 4. Decision gates & rollout (phasing — not architecture)

NOTE: sequencing only; the correctness requirements above do not depend on it.

1. **§0.5 inertia-factorization hoist** — do first. Bitwise, no parallelism, largest single ABA speedup
   (Jacobi eigensolves 11×→1×). Validate with the §5 equivalence test.
2. **§A OpenMP wavefront** — portable, bitwise for per-body work. Gate behind profiling showing the ABA is a
   material fraction of step time at scale, and behind `ROBO_ABA_OMP_MIN`.
3. **§B device-resident ABA** — only after §0.5/§A are reviewer-approved and profiling shows the
   `X_GB`/`bodyForceG` transfers or the CPU corrector dominate at the large-system scale. Not speculative
   (CLAUDE.md: correctness over performance).
4. **§3.3 DCA** — only if §B is depth-bound; researcher-owned, oracle-validated.

Each step is handed to `reviewer` before the next (CLAUDE.md: optimizer/coder never self-certifies).

## 5. Verification (designed now; run when the test gate is re-enabled)

- **§0.5 (bitwise):** `TestInertiaHoistEquivalence` — for a fixed `RobotState`, assert the corrector
  produces **bitwise-identical** `udot`/`A_GB`/`qdotdot` whether the inertia factorization runs per-sweep
  (old) or once-hoisted (new). Same reasoning as the position hoist.
- **§A (bitwise):** `TestAbaOmpEquivalence` — serial vs OpenMP per-body outputs (`X_GB`, `V_GB`, `P`,
  `PPlus`, `A_GB`, `udot`) **bitwise-identical** (§2.5), reduction scalars bounded (`1e-14` rel) or bitwise
  with the serial-sum variant. Thread counts 1, 2, 8.
- **§B (when attempted):** device-vs-host `udot`/`A_GB` A/B on ala-dipeptide (small) and a bushy system,
  tolerance per README §3; validate the per-round refit of frames/masses (the §0 taxonomy) by running ≥2
  rounds — the analogue of the stale-station regression test for 01/02.

Specified so the guarantees are checkable rather than asserted; none run in this pass per the user's
instruction.
