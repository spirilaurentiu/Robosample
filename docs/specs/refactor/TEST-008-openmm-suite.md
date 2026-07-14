# TEST-008: OpenMM's own test suite - off by default, all precisions, parallel

Status: final, 2026-07-12. Infrastructure spec, separate from the two-phase
freeze (it touches neither our source nor our tests). Executor: `coder`. Style:
`styles/spec.md`. It governs the build and CI wiring for the vendored OpenMM test
suite. Companion to `../../architecture/TESTS.md` and the root `CMakeLists.txt`
section "OpenMM's own test suite".

## Motivation

We plan to modify vendored OpenMM (a heavy-CUDA fork). Once we change OpenMM, its
own test suite is our regression oracle for those changes - it is the only
coverage that exercises the modified kernels directly. We do not own those tests
and SHALL NOT edit them; we own only how they are built and run.

The pre-existing wiring (`CMakeLists.txt`) had three gaps:

1. **Compiled unconditionally.** Under the `Tests` build type, ~100+ OpenMM test
   executables (53 CUDA + 40 serialization + 11 platform-independent) were always
   compiled, even for a dev who only touched our code. There was no off switch.
2. **Only GPU tests ran, and only in `mixed`.** `robosample_add_openmm_test`
   called `add_test` **only** in the `IS_GPU` branch and hard-coded
   `set(PRECISION mixed)`. The 11 platform-independent and 40 serialization tests
   (each with its own `main()`) were built but never registered - dead build
   cost. Numerical bugs that appear only in `single` or `double` were invisible.
3. **Serialized on the GPU.** Every GPU OpenMM test carried `RESOURCE_LOCK "gpu"`,
   so `ctest -j` ran them one at a time - the same lock our ASan sampling tests
   need, applied to tests that do not need it.

Original user framing: *"we will have to compile their tests ... i don't believe we
can run their tests in parallel since they want gpu resources, but i may be wrong
since they are only small tests ... they should be accessible from nox under some
special flag turned off by default (just like slow tests) ... ideally, openmm tests
should run in parallel for all precisions (mixed, float, double)."*

**Finding (parallelism, checked).** The user's doubt is settled in favor of
parallel. OpenMM's tests are small: the largest CUDA test builds a 10 000-particle
system, most use hundreds; each thin `.cpp` `#include`s a shared body
(`CudaTests.h` supplies `main()` + the `argv[1]` precision arg) and consumes tens
of MB of VRAM. On a 24 GB GPU dozens run concurrently. Separate processes get
separate CUDA contexts and separate memory, so sharing the GPU changes no result
- the only failure mode is VRAM exhaustion at very high `-j`, bounded by a `-j`
cap rather than a mutex. All 53 CUDA tests parse `argv[1]` as precision (the 3
that do not use `CudaTests.h` still read it as `CudaPrecision`), so all three
precisions are valid for every CUDA test.

## Behavior

- **`BUILD_OPENMM_TESTS` (CMake option, default `OFF`).** Gates compilation and
  registration of the entire OpenMM suite. Effective only under `BUILD_TESTING`
  (the `Tests` build type). OFF => the ~100+ executables are neither built nor
  registered; a normal `nox -s tests` pays nothing for them. This is the "special
  flag, off by default" analogue of the `ROBOSAMPLE_SLOW_TESTS` gate - but a
  build gate, not a runtime one, because the cost is compilation of third-party
  code, not sampling length.
- **All three precisions, in parallel.** Each GPU test registers three
  independent ctest cases - `<name>_single`, `<name>_mixed`, `<name>_double` -
  each labeled `openmm` and carrying **no** `RESOURCE_LOCK`. `ctest -j` runs them
  concurrently.
- **Non-GPU tests registered (build-gap fix).** CPU / reference / serialization /
  platform-independent tests register one case each (own `main()`, no precision
  arg), labeled `openmm`. They previously never ran.
- **Off by default in nox; opt-in session.** `nox -s tests` runs `ctest ... -LE
  openmm` (never OpenMM's suite). `nox -s openmm_tests` configures the Tests
  preset with `-DBUILD_OPENMM_TESTS=ON`, builds, and runs `ctest -L openmm` with a
  bounded default `-j 8` (override via posargs, e.g. `-- -j 16 -R Cuda_Nonbonded`).

## Invariants

- **Our tests are untouched.** `robosample_add_test` still stamps `LABELS
  "robosample"` and `RESOURCE_LOCK "gpu"` on our ASan sampling tests. They remain
  serialized on the GPU (they need it) and are unaffected by any OpenMM wiring.
- **No OpenMM source or test file is modified.** We do not own them. This ticket
  changes only `CMakeLists.txt` and `noxfile.py`.
- **The default developer and CI path is unchanged when the flag is off.** With
  `BUILD_OPENMM_TESTS=OFF`, target set, ctest case list, and `nox -s tests`
  behavior match the pre-change build exactly, minus the OpenMM targets that were
  never registered anyway.
- **The `.so` build is unaffected** - the OpenMM suite lives entirely inside
  `if(BUILD_TESTING)`, which excludes the Python-module config.
- **Label partition is total and disjoint:** every test carries exactly one of
  `robosample` or `openmm`; `-L openmm` and `-LE openmm` are exact complements.

## Interface

- CMake: `option(BUILD_OPENMM_TESTS ... OFF)`; the OpenMM test section wrapped in
  `if(BUILD_OPENMM_TESTS)`; `robosample_add_openmm_test` registers three precisions
  (GPU) or one case (non-GPU), all `LABELS "openmm"`, none locked.
- nox: `nox -s tests` gains `-LE openmm`; new `nox -s openmm_tests` session
  (configure with the option, build, `ctest -L openmm -j 8`, posargs override).
- Observable commands: `cmake --preset cuda-tests -DBUILD_OPENMM_TESTS=ON`;
  `nox -s openmm_tests`; `nox -s openmm_tests -- -j 16 --output-on-failure`.

## Validation strategy

- PRECONDITION: `cmake --preset cuda-tests` (flag off) configures clean; `ctest -N`
  lists zero `openmm`-labeled cases; the `robosample` case count matches the
  pre-change baseline. (Confirms off-by-default costs nothing.)
- INVARIANT: `cmake --preset cuda-tests -DBUILD_OPENMM_TESTS=ON` configures clean;
  `ctest -L openmm -N` lists `3 x (#GPU tests) + (#non-GPU tests)` cases;
  `ctest -LE openmm -N` is unchanged from the precondition. (Confirms the label
  partition and the three-precision expansion.)
- INVARIANT: on unmodified OpenMM, `nox -s openmm_tests` passes (0 or the tolerated
  ctest-exit-8 with a named failure investigated). Establishes the green baseline
  before any OpenMM fork change; a later OpenMM change is judged against it.
- INVARIANT: our own tests still serialize - `ctest -L robosample -j 8` shows one
  GPU test running at a time (the lock holds); `ctest -L openmm -j 8` shows many
  concurrent. Peak VRAM under the OpenMM run stays within the card (report it).
- LEMMA (parallelism is safe): the systems are <= 10 000 particles (evidence: the
  `numParticles` scan of `openmm/platforms/cuda/tests`), so N concurrent contexts
  fit in 24 GB for practical `-j`; correctness is context-local and independent of
  co-residency.

## Consequences and trade-offs

- Building the suite (flag on) compiles 100+ executables; ccache amortizes
  reconfigures. Toggling the flag on a shared build dir triggers one
  reconfigure/build of the OpenMM targets - expected, not a regression.
- `-j` is a VRAM ceiling knob, not a correctness knob. The default 8 is
  conservative for a 24 GB card; raise it when the card is larger or the tests
  smaller.

## Notes

- NOTE: the platform-independent (`openmm/tests`) and serialization tests run on
  whatever platform their `main()` selects (typically Reference); they are CPU
  work and parallelize freely. If a specific one needs an arg, that is an OpenMM
  concern, surfaced on first run - we do not patch their tests.
- NOTE: this spec is infrastructure and does not enter the Phase-A/Phase-B freeze.
  It MAY land at any time; it is sequenced before the OpenMM fork work it exists
  to guard.
