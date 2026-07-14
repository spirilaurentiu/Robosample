# SPLIT-W9: extract the HmcMove core (reinitialize + metropolis + generateSample) from World.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/World.cpp` (2619 LOC) is a god-object mixing ~12 responsibilities
(`ARCHITECTURE.md section 7`). This ticket extracts the **HMC move core** - the last and
most central World concern: the momentum draw + starting-Hamiltonian assembly
(`reinitialize`), the total-energy assembly at a proposed state
(`currentTotalEnergy`), the acceptance test (`metropolis`), and the move
dispatcher (`generateSample`, which routes NCMC/Cartesian/torsional/docking and
runs the constrained-Verlet trajectory).

`generateSample` touches every previously extracted service: it dispatches to
`ncmcMove` (W4), applies the docking kick via `repositionLigands` (W2), and
`reinitialize`/`currentTotalEnergy` fold in `calcFixman`/`calcLogSineSqrGamma2`
(W5), `setAtomsLocationsInGround` (W6), `drawSolventVelocities`/`calcSolventKE`
(W7) and `nmaKineticCorrection` (W8). It is therefore sequenced **last** among the
World splits (`README.md section 4`, `MODULES.md section 2`); all of W2/W4/W5/W6/W7/W8 are
prerequisites, because W9 runs last precisely so it depends on services already
extracted - each callee stays a `World` member with an unchanged declaration, so
the dispatch resolves without change, and sequencing W9 last keeps every sibling
ticket's moved-block diff off `HmcMove.cpp`. After this ticket, `World.cpp`'s
residual is identity/thermo/schedule setters + composition only.

**NOTE (scope).** `MODULES.md section 2` observes that `generateSample` (272 lines)
"collapses into a thin dispatch once W2/W4/W9 land" and the three move regimes
become a `Move` strategy. That strategy refactor is the **eventual** target, out
of scope here. W9 is pure code motion: it relocates the full `generateSample` body
verbatim, unchanged, dispatch branches and all.

Invariants below come from `ARCHITECTURE.md` and are restated in full.

## Moves (exactly what)
- `void World::reinitialize()` (`src/World.cpp` lines 1933-2064) -> `src/world/
  sampler/HmcMove.cpp` [private World member; decl `World.hpp:614`].
- `double World::currentTotalEnergy()` (lines 2102-2131) -> `HmcMove.cpp` [private
  World member; decl `World.hpp:616`].
- `bool World::metropolis(double, double)` (lines 2133-2148) -> `HmcMove.cpp`
  [private World member; decl `World.hpp:615`].
- `bool World::generateSample()` (lines 1559-1831) -> `HmcMove.cpp` [public World
  member definition relocated; declaration stays `World.hpp:305`].

No file-local anonymous helpers originate here; `reinitialize` /
`currentTotalEnergy` use `nmaDebugEnabled` from `include/robo/world/detail/
nma_debug.hpp` (hoisted by SPLIT-W8) - the include moves with the code out of
`World.cpp`. All state these methods read/write (`Hold_`, `savedPosG_`, `savedQ_`,
`generateSampleCalls_` (`World.hpp:684`), `equilPhase_`, `dockingStuckCount_`,
`sampler_`, `temperature_`/`RT_`/`beta_`, the RNG/distributions) stays declared on
`World`.

**PyBind / declaration stability.** Per CX-1 this is a multi-TU class split, not a
class extraction: the four declarations stay on `World`, only the bodies move as
`World::` definitions. `generateSample` is **public** but **not** directly
PyBind11-bound (grep of `src/PyBind11.cpp`: no binding); its live caller is the
Gibbs sweep in `Context.cpp` (`Context.cpp:461, 991, 1123`). Its declaration stays
at `World.hpp:305`, so the mangled symbol `World::generateSample()` is unchanged
and `Context.cpp` is not touched. `reinitialize` / `metropolis` /
`currentTotalEnergy` are private; their declarations stay at `World.hpp:614-616`.

## Public API after this ticket
No new public header. `World::generateSample` keeps its existing public
declaration in `World.hpp`. `src/world/sampler/HmcMove.cpp` is a `World`-class TU
(`#include "World.hpp"`) defining the four relocated members. IWYU: `<algorithm>`
(`std::copy`/`std::fill`/`std::min`), `<cmath>` (`std::isfinite`/`std::exp`/
`std::sqrt`), `<cstdio>` (`std::fprintf`), `<iostream>` (`std::cout` NMA trace,
lines 1995-2007, 2040-2042, 2109-2111), `<vector>`, `RobotIntegrator.hpp`
(`RobotEngine::stepTo`/`checkReversibility`, lines 1717, 1734; already a
`World.cpp` include), `detail/nma_debug.hpp`, `World.hpp`.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/World.cpp` loses lines 1559-1831, 1933-2064, 2102-2131, 2133-2148, and three
  includes complete their handoff at W9 (it is the last Lane-W ticket, so the last
  users of each leave `World.cpp`):
  - `#include "RobotIntegrator.hpp"` (`World.cpp:30`, templated `stepTo` /
    `checkReversibility`) -> `HmcMove.cpp`. After W4 relocated the NCMC `stepTo`
    calls (Stage-0 lines 2320, 2513), `generateSample`'s `checkReversibility` (line
    1717) and `stepTo` (line 1734) are the last non-comment users in `World.cpp`.
  - `#include <iostream>` (`World.cpp:20`, `std::cout` NMA trace) -> `HmcMove.cpp`.
    After W8 moved `nmaKineticCorrection`, the last `std::cout` in `World.cpp` is
    the NMA trace inside `reinitialize` (lines 1995-2007, 2040-2042) and
    `currentTotalEnergy` (lines 2109-2111). (W8 parked `<iostream>` in `World.cpp`
    "until W9"; W9 completes the handoff.)
  - `#include "robo/world/detail/nma_debug.hpp"` (added to `World.cpp` by W8) ->
    `HmcMove.cpp`. `reinitialize` (lines 1976, 2037) and `currentTotalEnergy` (line
    2108) are the last `nmaDebugEnabled` users; the include follows the code.
  `#include "NMA.hpp"` (`World.cpp:29`) is **not** a W9 concern: its only user is
  `setNMASoftModeFromHessian` (`computeRouteBNMA`, line 589), which W8 relocates.
  Nothing else in `World.cpp` changes.
- Files land at `MODULES.md section 1` path (`src/world/sampler/HmcMove.cpp`); recursive
  `src/` glob per SPLIT-W1.
- After this ticket `World.cpp` retains only: the constructor (lines 221-226),
  `add_sampler` (228-292), `setTemperature` (301-305), `setTimeStep`/`setMdSteps`/
  `setAcceptRejectMode` (307-317), `setBodyMassScale`/`setMassScaleByJoint`
  (369-387), `setReversibilityCheck` (389-401) - the residual World identity/
  thermo/schedule/composition surface (`MODULES.md section 1`).
- Invariants that SHALL remain true:
  - **INV-4 Realization order.** `RobotState` caches are valid only in stage order
    position -> velocity -> articulated-body inertias -> udot; the contract is
    enforced by convention and doc comments, not by types. This is the
    highest-risk invariant in the whole World split: `reinitialize`'s prelude
    (realizePosition -> realizeArticulatedBodyInertias, lines 1934-1935), its
    post-draw chain (evaluate -> realizeVelocity -> realizeArticulatedBodyInertias ->
    calcUDot -> calcQDot -> calcQDotDot, lines 2027-2032), and `currentTotalEnergy`'s
    (evaluate -> realizeVelocity -> calcKineticEnergy, lines 2103-2106) SHALL move
    verbatim in exactly this order. Any reordering silently corrupts the
    Hamiltonian.
  - **INV-5 Mass metric.** HMC momentum seeding and kinetic energy use the same
    mass-matrix operators (`multiplyBySqrtM` / `multiplyByMInv` /
    `calcKineticEnergy`); `calcLogDetM` is the Fixman kinetic term. `reinitialize`
    seeds `u = sqrt(RT)*M^{-1/2}*g` via `RobotEngine::multiplyBySqrtMInv` (line 2011)
    and `currentTotalEnergy` reads `calcKineticEnergy` (line 2106); the seeding
    operator and the KE operator SHALL stay the paired `sqrt(M^{-1})`/`1/2u^TMu`.
    Preserve the `ke_mix = ke - nmaCorr` and `+ keSolvent + fixman - 1/2RT logSineSqr`
    assembly (lines 2062, 2129) verbatim.
  - **INV-1 Force->wrench convention.** `reinitialize`/`currentTotalEnergy` call
    `bridge_.evaluate(state_)` (lines 2027, 2103), which reduces per-atom OpenMM
    forces to the per-body `(torque about Bo, net force)` Ground wrench; the
    host and CUDA reductions SHALL agree. No change here beyond relocation.
  - **INV-9 Quaternion double cover.** The constrained-Verlet trajectory
    (`RobotEngine::stepTo`, line 1734) advances quaternion joints with the exp-map
    and renormalizes each step; `generateSample`'s finite-q guard (lines 1738-1746)
    and reject-restore (lines 1810-1828) SHALL move verbatim so the double cover is
    handled and a rejected move restores the pre-move q bit-for-bit.
  - **INV-3 Stateless Worlds.** The move reads/writes only per-atom Ground
    coordinates as inter-world currency; the docking reject path restores from the
    saved Cartesian pose via `setAtomsLocationsInGround` (line 1815), never from a
    stale q. Preserve.
  - **HMC acceptance / detailed balance (`ARCHITECTURE.md section 3`, step 4).**
    `metropolis` accepts or rejects on the total-energy change `dH = Hnew - Hold`;
    on reject the saved `q` / positions are restored. Three properties SHALL move
    verbatim: (a) the burn-in / `AlwaysAccept` short-circuit (lines 2140-2142)
    returns true **before** drawing `uniform_`, so it SHALL NOT consume an RNG value
    - preserving the RNG stream; (b) the `dH <= 0` early-accept and the
    `uniform_(rng_) < exp(-beta*dH)` test (lines 2143-2147) are the exact Metropolis
    criterion over the `reinitialize`-assembled `Hold_` and the
    `currentTotalEnergy`-assembled `Hnew` - relocating `metropolis` apart from
    those two assemblers SHALL NOT perturb the paired energies; (c) the three
    reject-restore paths - Cartesian restores `savedPosG_` (lines 1580-1582),
    docking restores the pre-kick pose via `setAtomsLocationsInGround` (line 1815),
    torsional restores `savedQ_` then `realizePosition` + `fillAtomPositionsFromBodies`
    (lines 1825-1827) - leave the correct clean state for the next sweep. Restoring
    the wrong quantity breaks the reversibility the detailed-balance argument
    assumes. Move all three verbatim.
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None. No test `#include`s `World.cpp`; no `friend`s. The four World.hpp
  includers (`TestBatScalingJacobian`, `TestBatAnchorInvolution`,
  `TestTwoRobotContact`, `TestRoboticsOracleMolecule`) reach `generateSample`
  (where at all) through its unchanged public declaration; the C++ HMC-driver
  tests use `tests/HmcDriver.hpp` (its own analytic driver), not the private
  `reinitialize`/`metropolis`/`currentTotalEnergy`. No include-path fix is
  required by any test.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0 - in
  particular the deterministic HMC acceptance log (B3) is bit-identical, since
  this ticket relocates the acceptance math without touching it.
- `nm` diff on public symbols: empty (`generateSample` keeps its mangled name; the
  three helpers are private).
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- `HmcMove.cpp` <= 600 LOC (~450 expected: 273 + 132 + 30 + 16 moved-body lines plus
  includes); residual `World.cpp` <= 600 LOC.
- Comment-stripped before/after diff of `World.cpp` + `HmcMove.cpp` is empty apart
  from the moved blocks at their new location and the three relocated `#include`s.

---

## Findings for the orchestrator

- **No range overlap (clean interleave).** W9's Stage-0 ranges - `generateSample`
  1559-1831, `reinitialize` 1933-2064, `currentTotalEnergy` 2102-2131, `metropolis`
  2133-2148 - interleave with, but never overlap, the earlier Lane-W tickets:
  W7 (`setCartesianSolvent` / `drawSolventVelocities` / `calcSolventKE`, 1833-1931)
  sits between `generateSample` and `reinitialize`; W8 (`nmaKineticCorrection`,
  2066-2100) sits between `reinitialize` and `currentTotalEnergy`; W4
  (`configureNcmc...` / `ncmcMove`, 2150-2620) follows `metropolis`. No collision
  with W1 (601-968) or W8 (323-359, 575-596). Re-anchor by symbol.
- **`currentTotalEnergy` belongs to W9.** `MODULES.md section 2`'s W9 row names only
  `reinitialize` / `metropolis` / `generateSample`, but `currentTotalEnergy` is the
  end-of-leg energy assembler `generateSample` calls (line 1750) and shares the
  Fixman / solvent / NMA terms with `reinitialize`; W8 already treats it as a W9
  move (SPLIT-W8 shared-toggle note). Included here.
- **`NcmcMove.cpp` (W4) already lists `RobotIntegrator.hpp` as a required include**
  (SPLIT-W4 IWYU list) because `ncmcInnerGhmcStep` / `ncmcMove` instantiate
  `RobotEngine::stepTo<Bridge>`. W9 removes the include from `World.cpp` because
  `generateSample` is `World.cpp`'s last user of it.
