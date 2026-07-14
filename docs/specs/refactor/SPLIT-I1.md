# SPLIT-I1: split verletStep's body into drift / corrector / solvent helpers

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)

`include/RobotIntegrator.hpp` (639 LOC) carries the fixed-step velocity-Verlet driver
as function templates on the force-bridge type. `RobotEngine::verletStep`
(lines 102-452) is a single ~350-line template body doing three separable jobs
interleaved with shared local state: (1) the position drift (scalar Taylor +
quaternion exp-map + Cartesian-solvent drift), (2) the implicit-trapezoid velocity
corrector with its convergence loop, (3) the Cartesian-solvent velocity half-step.
This ticket extracts those three into named helper templates so `verletStep` reads as
a sequence of steps. It runs in the dynamics phase, before SPLIT-R3/R4 in the same
phase ordering does not matter for I1 (I1 touches only `RobotIntegrator.hpp`).

This is a behavior-preserving refactor, NOT byte-for-byte pure motion: the extracted
blocks share many `verletStep` locals (`h`, `q0`/`u0`/`udot0`/`qdd0`, `s`, `m`, the
solvent arrays `xs0`/`vs0`/`fs0`, the `evalPos`/`evalVel` lambdas), so extraction
threads that state through explicit parameters. The arithmetic in each block is
unchanged.

Invariants from `ARCHITECTURE.md section 5`, restated in full:

- **INV-9 Quaternion double cover.** Free/quaternion joints normalize quaternions
  each step; reversibility checks account for the double cover. The drift helper
  keeps the exact exp-map advance and the `normalizeQuaternions` call in the same
  position relative to the scalar Taylor drift.
- **INV-4 Realization order.** `RobotState` caches are valid only in stage order
  position -> velocity -> articulated-body inertias -> udot. `verletStep` establishes
  that order via `evalPos` (realizePosition + bridge + factorize) then `evalVel`
  (realizeVelocity + seed + calcUDot + qdots). The split preserves the single
  `evalPos` hoist and the per-sweep `evalVel` - the corrector helper receives
  `evalVel` by reference and calls it exactly where the current loop does.
- **INV-6 Constraint consistency.** The `G` assembly is identical across SHAKE,
  RATTLE, and the loop-closure Fixman log-det. `verletStep` calls
  `cset.enforcePositionConstraints` (after drift) and `cset.enforceVelocityConstraints`
  (after the corrector); both calls stay at their current boundaries between the new
  helpers, unchanged.

## Moves (exactly what)

New private static member function templates of `RobotEngine`, defined in
`RobotIntegrator.hpp` (declared in the header's integrator section; they need no
`RobotEngine.hpp` change if introduced as file-local `static` template helpers in
`RobotIntegrator.hpp` rather than class members - see Constraints). Each takes the
shared state it reads/writes by explicit parameter/reference.

- `driftPositions` - the scalar Taylor drift + quaternion exp-map drift.
  Body: `RobotIntegrator.hpp` lines 228-262 (the `for (i<nq)` Taylor loop, the
  `for (bodyIx ...)` quaternion loop with `wHalf`/`advanceQuatExp`, and the
  `normalizeQuaternions(m, s)` call). Reads `q0`, `u0`, `udot0`, `qdd0`, `h`, `m`,
  `s`; writes `q`. ORDERING (decided): **I1 runs before R4** (README section Phase A
  step 3, MODULES section 4). I1 extracts the whole drift loop - including the inline
  `wHalf`/`advanceQuatExp` quaternion block at lines 241-258 - into
  `driftPositions` verbatim. SPLIT-R4 then lifts that quaternion block out of
  `driftPositions` into `JointKernels::jointDriftQuat` and replaces it with a
  call. There is no line-range overlap: R4 re-anchors to the block's post-I1
  location inside `driftPositions`, not to `verletStep`.
- `cartesianSolventVerlet` - the Cartesian-solvent Verlet, which is TWO fragments:
  the position drift (lines 264-276, `x1 = x0 + h v0 + (h^2/2) a0`) and the velocity
  half (lines 442-448, `v1 = v0 + (h/2)(a0 + a1)`). Extract both into one helper with
  a `phase` selector, or two helpers `cartesianSolventDrift` / `cartesianSolventKick`;
  prefer the two-helper form so each is a single straight-line block moved verbatim.
  Reads `solvAtoms`, `solvInvM`, `xs0`/`vs0`/`fs0`, `posG`/`velG`/`frcG`, `h`.
- `velocityCorrector` - the implicit-trapezoid corrector.
  Body: `RobotIntegrator.hpp` lines 285-433 (the `u1_est` seed, the 10-iteration
  functional-iteration loop with the `velocitiesSane`/`evalVel` calls, the
  convergence bookkeeping, the dt-too-large `fprintf` guard). Reads `u0`, `udot0`,
  `h`, `nu`, the `velocitiesSane` and `evalVel` lambdas (by reference), `restorePreStep`
  (by reference); writes `u`, `*correctorConverged`. Returns the same
  `bool converged` / early-`false` outcomes verletStep currently produces.

`verletStep` after extraction is the driver: snapshot -> lambdas
(`forcesFinite`/`velocitiesSane`/`restorePreStep`/`evalPos`/`evalVel`) ->
`driftPositions` -> `cartesianSolventDrift` -> `refreshPos` +
`enforcePositionConstraints` -> `evalPos`/`evalVel` seed -> `velocityCorrector` ->
`enforceVelocityConstraints` + `realizeVelocity` -> `cartesianSolventKick` ->
`s.time += h`. Every call boundary matches a current line boundary.

`stepTo` (lines 454-469) and `checkReversibility` (lines 501-634) are unchanged.

## Public API after this ticket

Preferred: the three/four helpers are file-local `static` template helper functions
in `RobotIntegrator.hpp` (not class members), so `include/RobotEngine.hpp` needs NO
edit and no new public class member appears. They are templated on `Bridge` only
where they touch it (only `velocityCorrector` does, via `evalVel`); the drift and
solvent helpers are non-template and can be plain `inline` functions.

`RobotEngine.hpp` unchanged. No `robo_bindings...so` symbol delta (integrator templates
are header-only and instantiated per includer; helper templates add no exported
symbol).

## Constraints

- Behavior-preserving. Zero arithmetic change. The only new code is parameter
  threading; every moved expression is identical.
- Do NOT add private members to `RobotEngine` in `RobotEngine.hpp` - that would be a
  public-header edit. Use file-local helpers in `RobotIntegrator.hpp`.
- Invariants that SHALL remain true: **INV-9**, **INV-4**, **INV-6** (restated above).
  The two `cset.enforce*Constraints` calls and the `evalPos`/`evalVel` staging keep
  their exact positions.
- The `ROBO_CHECK` handling (the `#ifndef ROBO_CHECK` fallback at lines 45-48 and the
  end-of-file undef at 636-639) is unchanged; helpers that contain a `ROBO_CHECK` (via
  `evalPos`/`evalVel`) receive those lambdas by reference rather than re-expanding the
  macro.
- One commit for the extraction; formatting fixes separate.

## Predicted test breakage (Phase-A include fixes only)

**None.** `verletStep`/`stepTo`/`checkReversibility` keep their signatures; the new
helpers are file-local and unnamed by any test. The integrator tests
(`TestIntegrator.cpp`, `TestStability.cpp`, `TestEnsembleOrientation.cpp`,
`TestFreeJointKEPump.cpp`, `HmcDriver.hpp`, and the NCMC/contact suites) drive
`verletStep` through the same entry point with the same `AnalyticForceBridge`. No
include fix, no assertion-count change.

## Exit criteria (machine-checked)

- Full build passes; test suite / baseline outputs identical to Stage-0 (B1, B2, B3),
  both modes. B3 is bitwise for the deterministic HMC log - any divergence means an
  extracted block changed evaluation order; revert.
- `checkReversibility` residuals on the integrator tests match the baseline (the
  round-trip is the sharpest detector of a broken drift/corrector split).
- Public-symbol nm diff vs B0: empty.
- include-cycle script clean; no new include (the helpers use types already visible).
- clang-tidy / clang-format / IWYU clean on `RobotIntegrator.hpp`.
- `RobotIntegrator.hpp` shrinks; `verletStep` body drops from ~350 to a short driver.
  Header stays a template header (LOC cap is advisory for a header of function
  templates - claim the template-header exception if it exceeds 300 after the helper
  bodies are added, though net LOC should fall).
- No test file changed.
