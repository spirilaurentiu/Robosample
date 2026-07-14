# SPLIT-W8: extract the VelocityDistortion (NMA soft-mode + BAT drive) from World.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/World.cpp` (2619 LOC) is a god-object mixing ~12 responsibilities
(`ARCHITECTURE.md section 7`). This ticket extracts the two **distortion** concerns that
sit around the momentum draw and the driven-exchange position map:

1. NMA Route-B velocity distortion - `setNMASoftModeFromHessian` (computes the
   softest internal-coordinate mode and stores the per-DOF scale factors) and
   `nmaKineticCorrection` (the `RT ln cosh(w*mu)` term that restores detailed
   balance for the biased-mixture momentum draw).
2. BAT-scaling drive (`docs/specs/replica-exchange-nonequilibrium-work.md`) -
   `previewBatScaling` (side-effect-free preview + D6 Jacobian) and
   `applyBatScalingDrive` (commit to this world's geometry, D7 hard SHALL guards).

Both write only the distortion members (`uScaleFactors_`, `nmaBias_`,
`lastDistortJacobianDetLog_`, `lastNScaled_`). Sequenced after
W1/W3/W5/W2/W6/W7, before W9 (`README.md section 4`). Invariants below come from
`ARCHITECTURE.md` and are restated in full.

## Moves (exactly what)
- `World::BatScalingResult World::previewBatScaling(const std::vector<robo::Vec3>&,
  double, const std::unordered_map<int,double>&, const std::unordered_map<int,double>&)
  const` (`src/World.cpp` lines 323-330) -> `src/world/sampler/
  VelocityDistortion.cpp` [public World member definition relocated; declaration
  stays `World.hpp:584-588`].
- `void World::applyBatScalingDrive(double, const std::unordered_map<int,double>&,
  const std::unordered_map<int,double>&)` (lines 332-359) ->
  `VelocityDistortion.cpp` [public World member definition relocated; declaration
  stays `World.hpp:598-600`].
- `double World::setNMASoftModeFromHessian(const std::vector<double>&, double,
  double)` (lines 575-596) -> `VelocityDistortion.cpp` [public World member
  definition relocated; declaration stays `World.hpp:561-563`].
- `double World::nmaKineticCorrection()` (lines 2066-2100) ->
  `VelocityDistortion.cpp` [private World member; decl `World.hpp:624`].
- File-local anonymous-namespace helper `bool nmaDebugEnabled()` (lines 56-62) ->
  see the shared-toggle note below [used by `nmaKineticCorrection` line 2093 **and**
  by W9's `reinitialize`/`currentTotalEnergy`].

`World::BatScalingResult` (`World.hpp:579-583`) **stays a public nested struct on
`World`** - it is not relocated. Tests reference the exact name
`World::BatScalingResult`; moving it would change that name and break them.
Likewise the distortion members (`uScaleFactors_`, `nmaBias_`,
`lastDistortJacobianDetLog_`, `lastNScaled_`, `World.hpp:725/731/707/708`) and the
`getDistortJacobianDetLog`/`getLastNScaled` inline accessors (`World.hpp:605-610`)
stay declared on `World`.

### Shared-toggle note (cross-cutting with W9)
`nmaDebugEnabled` is used by `nmaKineticCorrection` (W8) and by `reinitialize` /
`currentTotalEnergy` (W9). Because it is an anonymous-namespace symbol it cannot be
shared across two TUs. W8 hoists it verbatim to
`include/robo/world/detail/nma_debug.hpp` as an `inline bool nmaDebugEnabled()`;
`VelocityDistortion.cpp` and the (still-present) `World.cpp` both `#include` it.
When W9 moves `reinitialize`/`currentTotalEnergy` out of `World.cpp`, the include
follows the code. This is a small, verbatim, single-source hoist - no logic
change.

## Public API after this ticket
No new `world/` public header for the methods (they stay `World` members,
declarations unchanged). `include/robo/world/detail/nma_debug.hpp` exports
`inline bool nmaDebugEnabled()`. `src/world/sampler/VelocityDistortion.cpp` is a
`World`-class TU (`#include "World.hpp"`) defining the four relocated members.
IWYU: `<cmath>` (`std::sqrt`/`std::abs`/`std::exp`/`std::log1p`/`std::log`),
`<iostream>` (`std::cout` trace, lines 1995-2007, 2096), `<unordered_map>`,
`<vector>`, `NMA.hpp` (`robo::computeRouteBNMA`/`RouteBNMA`, line 589),
`BatScaling.hpp` (`robo::applyBatScaling`, line 328; already via `World.hpp`),
`detail/nma_debug.hpp`, `World.hpp`.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/World.cpp` loses lines 56-62 (the toggle, hoisted), 323-359, 575-596,
  2066-2100, and gains `#include "robo/world/detail/nma_debug.hpp"` (still needed
  by `reinitialize`/`currentTotalEnergy` until W9). The `<iostream>` include stays
  in `World.cpp` (still used by the `reinitialize` NMA trace at lines 1995-2007
  until W9). Nothing else changes.
- Files land at `MODULES.md section 1` paths; recursive `src/` glob per SPLIT-W1;
  `include/robo/world/detail/` is covered by the recursive `include/` glob.
- Invariants that SHALL remain true:
  - **INV-5 Mass metric.** HMC momentum seeding and kinetic energy use the same
    mass-matrix operators (`multiplyBySqrtM` / `multiplyByMInv` /
    `calcKineticEnergy`). `nmaKineticCorrection` forms `w = M^{1/2} u / sqrt(RT)` via
    `RobotEngine::multiplyBySqrtM` after `realizeArticulatedBodyInertias` (lines
    2082-2084) - the SAME `multiplyBySqrtM` operator the draw inverts. The
    correction `RT ln cosh(w*mu)` SHALL move verbatim so `ke_mix = ke - nmaCorr`
    keeps the Route-B mixture move detailed-balanced, and SHALL stay a strict no-op
    (returns 0) off `DistortOption::NMA` (lines 2072-2074) so ordinary HMC is
    untouched.
  - **INV-7 BAT map/Jacobian agreement.** `applyBatScaling` and its Cartesian
    log-Jacobian read identical `(r,theta)` geometry; one anchor snapshot per round is
    shared by both swap partners so the paired map is an exact involution.
    `previewBatScaling` is the single side-effect-free path that both the real
    drive and the finite-difference oracle call (line 351 reuses it), guaranteeing
    map and Jacobian read the same geometry; preserve that shared call so the
    involution property holds.
  - **INV-8 REX detailed balance.** The BAT drive feeds the nonequilibrium-work
    exchange; `applyBatScalingDrive`'s D7 hard SHALL - throw if `mdSteps != 0` or
    if `distortOption != ScaleBendStretch` (lines 337-346) - SHALL move verbatim:
    the exact two-endpoint work acceptance is valid only in the pure deterministic
    -scaling-map limit (no post-scale MD).
  - **INV-9 Quaternion double cover.** `setNMASoftModeFromHessian` positions the
    world at the minimized geometry via `setAtomsLocationsInGround` (line 587),
    which re-seeds identity quaternions; preserve the call so the Hessian is
    evaluated at q0=0 in the reset frames.
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None. No test `#include`s `World.cpp`; no `friend`s. `TestBatScalingJacobian.cpp`
  and `TestBatAnchorInvolution.cpp` (`#include "World.hpp"`) call
  `world.previewBatScaling(...)` and use the type `World::BatScalingResult`
  through the unchanged public declarations (`TestBatScalingJacobian.cpp:185-186,
  222`; `TestBatAnchorInvolution.cpp:99-100`); because `BatScalingResult` stays a
  nested `World` type and `previewBatScaling` stays public, both keep `#include
  "World.hpp"` verbatim. `setNMASoftModeFromHessian` / `nmaKineticCorrection` are
  not named by any test. No include-path fix is required by any test.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0.
- `nm` diff on public symbols: empty (`previewBatScaling`, `applyBatScalingDrive`,
  `setNMASoftModeFromHessian` keep their mangled names; `nmaKineticCorrection` is
  private; `nmaDebugEnabled` becomes `inline` - no exported non-inline symbol;
  `World::BatScalingResult` unchanged).
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- `VelocityDistortion.cpp` <= 600 LOC (~110 expected); `nma_debug.hpp` <= 300 LOC
  (~15 expected).
- Comment-stripped before/after diff of every touched file is empty apart from the
  moved blocks and the added `#include`.
