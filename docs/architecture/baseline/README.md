# VERIFY Phase-0 baselines - capture record

Status: captured 2026-07-13, POST-split (the source split had already executed when
these were taken; the pre-split reference is git HEAD). Companion to
[`../../specs/refactor/VERIFY.md`](../../specs/refactor/VERIFY.md) section 1.

NOTE on timing: VERIFY.md section 1 specifies capturing these BEFORE any SPLIT ticket
runs. That did not happen; the splits ran first. So these artifacts describe the
post-split tree, and behavior-preservation is judged against git HEAD (the pre-split
commit) rather than a stored pre-split baseline.

## Captured

- **B0 - exported/defined symbols** (`symbols.txt`): `nm -C --defined-only` on
  `build/cuda-release/robo_bindings.cpython-312-x86_64-linux-gnu.so`, 11440 symbols.
  API-delta check against the accepted set (README CX-3): `Context::setSeparateForceGroups`
  and `Context::setEnforcePeriodicBox` present (+2); `computePeriodicBoxVectors_Context`
  absent (-1). No other intended public-symbol change.

- **B3 - representative run** (`run_l1_ala-dipeptide.txt`): the CLAUDE.md Level-1
  example (ala-dipeptide, seed 6000, 10 prod steps, `--validate true`). RESULT: PASS.
  Startup geometry clean (32 atoms, 0 hard clashes); 10 HMC moves (9 accept / 1
  reject) with energy conservation (dH ~ +-0.5 kJ/mol); force-group validation
  `result: PASS`, C++ PME PE -128.7423 vs OpenMM ref -128.7422, all six force groups
  matching within 1e-4..1e-6. Teardown emits benign OpenMM "Error deleting array:
  CUDA error" lines (destructors after context teardown; exit 0, validation passed).

## B1/B2 - gtest pass set (captured 2026-07-13)

`ctest -L robosample` (FAST tier, `ROBOSAMPLE_SLOW_TESTS` unset) run x3, byte-identical
each run (`passset.txt`, `passset_summary.txt`): 321 cases, 304 passed, 16 slow-tier
skips, 1 failed. `flaky.txt`: none (3 runs identical).

## Full exhaustive tier (2026-07-13) - GREEN

The full slow+fast suite (`ROBOSAMPLE_SLOW_TESTS=1 ctest -L robosample -j16`, openmm
excluded) runs in ~7:45 parallelized (the GPU RESOURCE_LOCK was scoped to `TestAlchemy`
alone by TEST-001, so the analytic-bridge tests run concurrently). Result: **321 tests,
0 failed** - every statistical oracle (Fixman slow tier, 3M-draw ensemble validation,
two-robot contact, robotics oracles) passes at the exhaustive tier. The 3 NCMC
explicit-solvent oracles are GATED (`GTEST_SKIP`) pending their in-flight campaign
(`known_failures.txt`), reviewer-confirmed byte-identical to HEAD (not a refactor
regression).

## Not captured here

- **B4/B5 - include graph / coverage map.** Not captured this pass.

## Build status

`cuda-release` compiles clean and links/install `robo_bindings.so`. One cross-lane
defect was fixed post-split: the I1-extracted free function `driftPositions`
(`RobotIntegrator.hpp`) called `normalizeQuaternions` unqualified; corrected to
`RobotEngine::normalizeQuaternions`.
