# 03 — ABA recursion parallelization

Concern: parallelize the articulated-body recursions without changing their math. The passes are a hard
**sequential chain** across the step (position → force → velocity → articulated inertias → `calcUDot`),
but **within** each pass there is real, structured parallelism: independent subtrees (the forest) and,
inside a tree, independent bodies at the same depth (`bodyLevel` wavefront). Read [`README.md`](README.md)
§3 (numerical contract) and §4 (CPU compatibility) first.

This file has two parts: **§A** portable OpenMP (Phase 3, specified for implementation) and **§B** a CUDA
device-resident ABA (Phase 4, future — sketched, not specified to build).

## 1. The dependency structure (what may and may not run concurrently)

From the map (`RobotEngine.cpp`), with invariant `bodyParent[b] < b`:

- **Outward passes** — `realizePosition` (`:532`), `realizeVelocity` (`:614`), `calcUDot` pass 2 (`:1039`),
  `multiplyBySqrtMInv` (`:1115`): body `b` reads its **parent**'s just-written result
  (`X_GB[b]=X_GB[p]*X_PB[b]`; `V_GB[b]=~Phi[b]*V_GB[p]+…`; `APlus=~Phi[b]*A_GB[p]`). A body may run only
  after its parent.
- **Inward passes** — `realizeArticulatedBodyInertias` (`:790`), `calcUDot` pass 1 (`:1017`),
  `multiplyByMInv` (`:1062`), `calcMobilizerReactionForces` (`:1257`): body `b` reads **all its
  children**'s results (`Pb += PPlus[c].shift(...)`; `z += Phi[c]*ZPlus[c]`). A body may run only after
  all its children.
- **Independent bodies:** any two bodies with no ancestor/descendant relationship are independent within a
  pass. In particular (a) **different molecules are disjoint subtrees rooted at Ground** — fully
  independent for the whole pass; (b) **all bodies at the same `bodyLevel` are mutually independent** —
  this is the wavefront the code comments reserve `bodyLevel` for (`RobotModel.hpp:67`,
  `RobotEngine.hpp:20-21`).
- **Fully independent (no inter-body dependency at all):** the per-atom position loop (`:606-611`, being
  moved to GPU in `01`), the force scatter (`ForceBridge.hpp:86`, moved in `02`), `calcKineticEnergy`
  (`:1228`, a reduction), `calcLogDetM` (`:1202`, independent per-body determinant + sum), the `abcf`
  seed loop (`:990`), `calcQDot`/`calcQDotDot` (per-body `q`-blocks). These parallelize trivially.

NORMATIVE: the parallelization SHALL preserve the two load-bearing details flagged in the map — the child
sum uses the child's *outboard* quantity shifted by the child's own `Phi`
(`PPlus[c].shift(Phi[c].l())`, `Phi[c]*ZPlus[c]`), and `abcf` is seeded from the **mobilizer** Coriolis
`a_mob`, not total `a_tot` (`RobotEngine.cpp:979-992`). Reordering iteration must not disturb these.

## 2. §A — OpenMP level-wavefront + forest parallelism (Phase 3, portable)

### 2.1 Precompute the schedule once per world (NORMATIVE data, no math change)

At `World::buildModel`, alongside the existing `bodyLevel`, build a **level CSR**: `levelBeg[L]`,
`levelEnd[L]` into a `bodyByLevel[]` array listing bodies grouped by (and contiguous within) level, in
increasing `b` order inside a level. Bodies are already relabeled into level-major topological order
(`World.cpp:599-634`), so `bodyByLevel` is likely the identity permutation or a trivial slice — verify and,
if so, iterate `bodyByLevel` as a plain `[levelBeg[L], levelEnd[L])` range over the existing arrays with no
extra indirection. This schedule is a **run-constant** (README §2); compute it once, never per step.

### 2.2 Parallelize each pass over the wavefront

Replace each `for (b = 1; b < numBodies; ++b)` (outward) / `for (b = numBodies-1; b >= 1; --b)` (inward)
with a level loop that is **serial across levels** and **parallel within a level**:

```
// outward
for (int L = 1; L <= maxLevel; ++L)
    #pragma omp parallel for schedule(static)          // bodies at level L are independent
    for (int k = levelBeg[L]; k < levelEnd[L]; ++k) { int b = bodyByLevel[k]; /* body-b work */ }

// inward
for (int L = maxLevel; L >= 1; --L)
    #pragma omp parallel for schedule(static)
    for (int k = levelBeg[L]; k < levelEnd[L]; ++k) { int b = bodyByLevel[k]; /* body-b work */ }
```

Correctness rests on the wavefront property: within level `L`, every body's parent is at level `<L`
(already finalized, outward) and every body's children are at level `>L` (already finalized, inward). No
two bodies in the same level write each other's inputs.

NOTE: the inward child-sum `Pb += PPlus[c].shift(...)` reads children at level `L+1` and writes `P[b]`/
`PPlus[b]` at level `L`; siblings at level `L` touch disjoint `b`, so the writes don't collide. No atomics
or locks are needed — the SoA is indexed by `b` and each task owns its `b`.

### 2.3 The fully-independent loops

`calcLogDetM`, `abcf` seed, `calcQDot`, `calcQDotDot`, and the per-body dense algebra inside
`realizeArticulatedBodyInertias` (the `D=~H P H` build, `invertDense` Jacobi eigensolve, `symSqrt`/
`symSqrtInv`) are per-body-local: a plain `#pragma omp parallel for` over all bodies (no level structure)
is correct. `calcKineticEnergy` and `calcLogDetM`'s scalar sum use `reduction(+:...)`.

### 2.4 Cost model and gating (RECOMMENDED)

- Speedup ceiling is `numBodies / (critical-path length)`. Critical path = **tree depth** = `maxLevel`,
  not `numBodies`. A long unbranched backbone (depth ≈ numBodies) parallelizes poorly; a bushy system
  (many molecules, many short sidechains — solvent, membrane, the target 100k-body regime) parallelizes
  well. Wavefront width = bodies per level is the available parallelism.
- Fixed OpenMP fork/join per level costs ~µs; for the small worlds (ala-dipeptide, `numBodies` ~ tens) the
  overhead dominates. RECOMMENDED: gate the parallel path on a per-world threshold
  (`numBodies >= ROBO_ABA_OMP_MIN`, e.g. a few hundred) computed once at `buildModel`; below it, take the
  current serial loop. A single `if (parallel)` branch around each level loop, or a `num_threads(1)`
  clause, keeps one code path.
- Thread count: default to OpenMP's, but the engine SHALL respect an env override (e.g. `OMP_NUM_THREADS`)
  and MUST NOT oversubscribe against OpenMM's own CPU threads when both run on CPU. Since the ABA and the
  OpenMM force eval run in strict sequence (never concurrently), they can share the full core count; no
  partitioning needed.

### 2.5 Numerical effect

`schedule(static)` with a fixed thread count gives a deterministic assignment, but per-body results still
combine in an order that can differ from the serial loop only where a body *sums over its children*
(inward passes). That child-sum order is **fixed by the CSR** (`bodyChildren` order), not by thread
scheduling — each body sums its own children sequentially inside its task — so the inward passes are in
fact **bitwise-identical** to serial. The outward passes have no cross-body reduction at all → also
bitwise-identical. The only genuinely order-sensitive reductions are the global scalars
(`calcKineticEnergy`, `calcLogDetM` sum), where `reduction(+:)` changes association: declare that under
README §3 (bounded, ~`1e-15` relative). **This is a stronger guarantee than Phases 1–2**: the ABA
per-body algebra is order-preserving under this schedule.

NOTE: if a `reduction(+:)` scalar difference is undesirable, compute the global sums by a serial pass over
the (already parallel-filled) per-body array instead of an OpenMP reduction — that restores bitwise
identity for those scalars too, at negligible cost (`O(numBodies)` adds, off the hot inner work).

## 3. §B — CUDA device-resident ABA (Phase 4, future — sketch only)

NOT specified for implementation here; recorded so Phase 1–3 choices don't foreclose it.

Motivation: for the 1M-atom / 100k-body regime, once `01`/`02` keep positions and forces on the device,
the remaining host↔device traffic is `X_GB` up and `bodyForceG` down each step. A device-resident ABA
would let `bodyForceG` feed `calcUDot` **without leaving the device**, and `X_GB` be produced on-device by
a device `realizePosition`, closing the loop so **nothing per-step crosses PCIe except at accept/DCD
cadence**.

Shape:

- One CUDA thread (or small tile) per body; dense per-body algebra (dof ≤ 6: the `D=~H P H` build, the
  `DI` inverse, `symSqrt`) is a natural per-thread workload.
- The level CSR from §2.1 is the launch schedule: **one kernel launch per level** (serial across levels,
  parallel within), outward and inward. Kernel-launch latency (~µs each) × `2·maxLevel` launches per pass
  is the main cost risk; for deep trees a single persistent-kernel + grid-sync (cooperative groups)
  variant that walks levels in-kernel avoids per-level relaunch. Evaluate both.
- Data already lives in the SoA arena; a device mirror of the per-step `RobotState` slabs
  (`X_GB, Phi, P, PPlus, H, G, DI, Z, zPlus, A_GB, …`) plus the constant `RobotModel` topology (uploaded
  once, README §2) is required. The 64-byte-aligned arena layout (`MemoryArena.hpp`) maps cleanly to
  coalesced device arrays.
- Numerical: this is the largest departure — full ABA in device `double`, order/FMA differences
  throughout. It is a distribution-preserving, trajectory-diverging change (README §3) and would need its
  own A/B campaign (compare device `udot` vs host `udot` for a fixed state) before use.

Decision gate for starting Phase 4: only after Phases 1–3 are reviewer-approved **and** profiling (when
re-enabled) shows the `X_GB`/`bodyForceG` transfers or the host ABA are a material fraction of step time
at the large-system scale. Do not start it speculatively (CLAUDE.md: correctness over performance, no
speculative optimization).

## 4. Verification (deferred per README §3, designed now)

- **§A bitwise cases:** `ctest` `TestAbaOmpEquivalence` — for a fixed `RobotState`, run each pass serial
  and OpenMP-parallel and assert **bitwise-identical** per-body outputs (`X_GB`, `V_GB`, `P`, `PPlus`,
  `A_GB`, `udot`) per §2.5, and bounded (`1e-14` rel) equality for the reduction scalars unless the
  serial-sum variant (§2.5 NOTE) is used, in which case bitwise. Run at thread counts 1, 2, 8.
- **§B (when attempted):** device-vs-host `udot`/`A_GB` A/B on ala-dipeptide and a bushy system, tolerance
  declared per README §3.

Neither is run in this pass, per the user's instruction; both are specified so the guarantees in §2.5 are
checkable rather than asserted.
