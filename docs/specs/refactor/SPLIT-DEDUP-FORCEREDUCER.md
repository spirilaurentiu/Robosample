# SPLIT-DEDUP-FORCEREDUCER: hoist the per-body force->wrench reduction into one shared function

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
The atom-force -> per-body spatial-wrench reduction exists **twice** with no shared
code (ARCHITECTURE.md section 6.5, defect 5):

- **Host loop:** `ForceBridge::getForcesFromOpenMM` in
  `include/ForceBridge.hpp` - the per-atom accumulation at lines 96-118 (inside the
  method spanning lines 66-119): for each real atom, `BF[b].linear += f` and
  `BF[b].angular += (posG[a] - X_GB[b].p()) x f`.
- **Device kernel:** `reduceForces` in `src/OpenMMContext.cpp` lines 132-170 (the
  K2 nvrtc kernel inside `kKinematicsKernelSource`): the same accumulation with
  `r = X_GB.R * station`, `linear += f`, `angular += r x f`, skipping virtual
  sites.

Two independent copies of the reduction are a divergence risk against INV-1 (a fix
to one silently skips the other) and neither is unit-testable in isolation. This
ticket - the highest-value bridge ticket - hoists the reduction into ONE host
function `reduceAtomForcesToBodies(...)` in `bridge/ForceReducer.{hpp,cpp}`,
routes the `ForceBridge` host path through it, and pins the CUDA kernel to it with
a differential test.

This is a **dependency-defect resolution ticket (VERIFY section 3)**: it may narrow
includes and adds a test that enters the B1 baseline BEFORE the dedup lands.

### Invariants this ticket is built to protect (restated in full)

- **INV-1 Force->wrench convention.** Per-body spatial force in
  `RobotState::bodyForceG` is `(torque about the body origin, net force)` in the
  Ground frame. Both the host reduction (`ForceBridge.hpp`) and the CUDA
  `reduceForces` kernel MUST produce this identically. Concretely: for body `b`
  with origin `X_GB[b].p()` in Ground, `bodyForceG[b].linear = Sigma_a f_a` and
  `bodyForceG[b].angular = Sigma_a (r_a - origin_b) x f_a`, sum over the body's real
  atoms, `f_a` the Ground-frame per-atom force. The host computes the lever arm as
  `posG[a] - X_GB[b].p()`; the kernel computes it as `X_GB[b].R * station_a`
  (algebraically equal - the station rotated into Ground is the atom offset from
  the body origin). The single hoisted function encodes exactly this convention so
  the two paths cannot drift.
- **INV-2 Virtual sites.** Virtual-site forces are already projected onto parent
  atoms by OpenMM; the reduction skips massless particles to avoid double counting.
  Host: `if (model_.atomMass[a] == 0) continue;` (`ForceBridge.hpp:107-109`).
  Kernel: `if (isVirtual[a]) continue;` (`OpenMMContext.cpp:153`). The hoisted
  function takes the skip predicate (mass==0 / isVirtual) and applies it once,
  identically. Re-reducing a virtual site's leftover force slot would double-count
  the M-site force (~1e3 kJ/mol/nm per water), a non-conservative kick - the exact
  bug the skip prevents (`ForceBridge.hpp:99-106`).

## Moves (exactly what)
- New free function `reduceAtomForcesToBodies(...)` -> `bridge/ForceReducer.{hpp,cpp}`
  [public free fn]. Signature captures the reduction inputs as POD/`robo::` spans:
  per-atom Ground forces, per-atom Ground positions (or body-frame stations +
  `X_GB`), `atomBody`, per-atom skip mask (mass==0), body origins (`X_GB[b].p()`),
  `numAtoms`, `numBodies`; output `robo::SpatialVec* bodyForceG`. It performs the
  INV-1/INV-2 accumulation and NOTHING else (no OpenMM types, no mobilityForce
  clearing, no atom-force caching - those stay in `ForceBridge`).
- `ForceBridge::getForcesFromOpenMM` host reduction loop
  (`include/ForceBridge.hpp` lines 96-118) -> replaced by a call to
  `reduceAtomForcesToBodies(...)`. The surrounding responsibilities STAY in
  `getForcesFromOpenMM`: the `evaluateForcesFromPositionsCache` call (67-68), the
  `BF[b]` zero-init (71-73), the `mobilityForce` clear (83-86, contract at 74-82),
  and the optional per-atom `atomForce` caching (92-93, 112-114). Whether the
  atom-force cache write stays a separate pass or is folded as an out-parameter of
  the reducer is an implementation choice - keep it observably identical.

The CUDA `reduceForces` kernel (`src/OpenMMContext.cpp:132-170`) is NOT rewritten
here (rewriting nvrtc to call a host function is impossible). It is instead PINNED
to `reduceAtomForcesToBodies` by the differential test below and stays the device
implementation of the same convention. (SPLIT-O6 later relocates that kernel to
`bridge/GpuKinematics.cpp`; the test follows it.)

## Public API after this ticket
`bridge/ForceReducer.hpp` exports `reduceAtomForcesToBodies(...)` (the one shared
host reduction). Header self-contained: `#include "RobotModel.hpp"` /
`"RobotState.hpp"` as needed for `robo::SpatialVec`/`robo::Vec3`/`robo::Transform`,
`#include "robot_math.hpp"`, `<vector>`. No OpenMM include (the reducer is
OpenMM-free; `ForceBridge` converts `OpenMM::Vec3` forces to `robo::Vec3` before
calling, or the reducer takes a raw `double`/`robo::Vec3` span). `#pragma once`.

`ForceBridge` and `OpenMMContext` public symbol sets are unchanged (the reduction
was inline, never an exported symbol). No public nm delta expected; if the reducer
signature is chosen to also serve as a public utility, that is internal API, not a
change to the B0 public contract.

## Constraints
- Behavior-preserving hoist. The host path result MUST be byte-identical to the
  current inline loop (same accumulation order per body, same skip, same lever-arm
  algebra). Zero logic change to the numbers `getForcesFromOpenMM` produces.
- Include edge (VERIFY section 3, checked vs B4): `ForceBridge.hpp` gains
  `#include "ForceReducer.hpp"`. `ForceReducer.hpp` introduces NO edge into OpenMM
  (it is OpenMM-free), so it does not deepen the section 6.1 solver->OpenMM inversion; it
  narrows toward resolving section 6.5. State the before/after edge explicitly in the PR.
- Invariants that SHALL remain true: **INV-1** and **INV-2**, restated in full
  above. The hoisted function is their single definition.
- New `ForceReducer.{hpp,cpp}` auto-globbed into `robosample_objects`
  (`CMakeLists.txt`); no CMake edit. One commit for the hoist; the test
  lands in its own commit BEFORE the hoist (see below); formatting separate.

## Required test (VERIFY section 3) - enters the B1 baseline BEFORE the dedup lands
Add a gtest (e.g. `tests/TestForceReducer.cpp`, auto-globbed by
`CMakeLists.txt`) that:

1. Builds a FIXED input: a small multi-body model with known per-atom Ground
   positions/stations, `X_GB` transforms (nontrivial rotation + translation), a
   fixed per-atom force vector, and at least one massless virtual site with a
   nonzero force in its own slot.
2. Computes per-body wrenches via the host `reduceAtomForcesToBodies`.
3. Computes per-body wrenches via the CUDA `reduceForces` kernel (driven through
   `OpenMMContext::reduceForcesToBodies`, i.e. the real device path) on the same
   input, under `USE_CUDA`.
4. Asserts the two per-body `(angular, linear)` wrenches are identical to a tight
   tolerance (mixed-precision fixed-point on device -> use the existing tight band,
   not exact-bit), AND that the virtual-site force is skipped by both (INV-2): the
   body wrench equals the sum over REAL atoms only.

This test SHALL be committed and pass (host==CUDA on the current, still-duplicated
code) BEFORE the dedup commit, so it enters B1 as a baseline guard. Under a
non-CUDA build the CUDA half is `GTEST_SKIP()`-gated; the host half still runs.
Per VERIFY section 3 this closes the INV-1 divergence risk between the two reductions.

## Predicted test breakage (Phase-A include fixes only)
- None to existing tests. The host reduction result is unchanged, so every test
  that exercises `ForceBridge::evaluate`/`getForcesFromOpenMM` (the OpenMM-backed
  tests: `TestAlchemy`, `TestNcmcExplicitSolvent`, `TestNCMCWork`,
  `TestReactionForces` via its analytic bridge, and the Level-1 example) produces
  identical numbers. `tests/ForceBridge.hpp` is a TEST-ONLY OpenMM-free stub
  (`class ForceBridge { void evaluate(RobotState&){} };`) and does NOT include the
  production reduction, so it is untouched. The one NEW file is
  `tests/TestForceReducer.cpp` (added by this ticket, not an edit to an existing
  test).

## Exit criteria (machine-checked)
- `TestForceReducer` present and passing in B1 (host==CUDA within tolerance,
  INV-2 skip verified) BEFORE the dedup commit; still passing after.
- Full build passes in `cuda-release` and `cpu-release`; pass set / assertion
  counts identical to B1/B2 (modulo the one added test, recorded in the baseline);
  representative outputs identical to B3.
- nm diff on public symbols vs B0: empty.
- include-cycle script clean vs B4; the stated `ForceBridge.hpp -> ForceReducer.hpp`
  edge is the only include change and introduces no cycle or OpenMM dependency;
  clang-tidy / clang-format / IWYU clean.
- `ForceReducer.hpp` <= 300 LOC; `ForceReducer.cpp` <= 600 LOC (actual small).
- Comment-stripped before/after diff of `ForceBridge.hpp`: the only non-empty
  hunk is the inline loop replaced by the `reduceAtomForcesToBodies` call plus the
  added `#include`.
