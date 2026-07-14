# TEST-007: Mirror the module map in the test tree

Status: draft, 2026-07-12. Phase-B test ticket (source frozen). Executor: `coder`. Style: `styles/spec.md`. Companion to
[`../../architecture/TESTS.md`](../../architecture/TESTS.md) section 5.9 (target
structure) and section 5.1-5.8 (the preceding TEST tickets), [`MODULES.md`](../../architecture/MODULES.md)
section 1 (the target directory layout this mirrors), the architect
test-documentation/target-structure sections, and [`VERIFY.md`](VERIFY.md) section 5.
Test-only motion; navigability + documentation discipline, zero behavior change.

## Context (why)

The 47 `Test*.cpp` files sit in a flat `tests/` directory with no correspondence
to the module map. `TESTS.md section 5.9` and the architect charter require the test
tree to **mirror the module map one-to-one**: `src/<module>/X.cpp <->
tests/<module>/XTest.cpp`, one fixture per class, behavior-sentence test names,
shared helpers under `tests/support/` held to production header discipline. This
ticket lands after TEST-002...006 have introduced `tests/support/` and split the
monoliths; it places every remaining test file into the mirrored directory and
brings each to the test-documentation standard.

The source module map is `MODULES.md section 1` - `workflow/`, `world/` (with
`world/sampler/`), `dynamics/`, `bridge/`, `model/`, `math/`, `io/`, `util/`. The
test tree SHALL reproduce it: e.g. `src/world/sampler/HmcMove.cpp <->
tests/world/sampler/HmcMoveTest.cpp`, `src/dynamics/RobotIntegrator.hpp <->
tests/dynamics/RobotIntegratorTest.cpp`.

## Scope

- **In:** relocate each `Test*.cpp` into the mirrored `tests/<module>/` directory
  and rename to `<ClassUnderTest>Test.cpp`; move the shared header helpers into
  `tests/support/`; bring every file to the file/fixture/test documentation
  standard; update CMake discovery to recurse.
- **Out:** any engine source; the split/consolidation tickets' *content* (this
  ticket relocates and documents the results of TEST-002...006, it does not
  re-split); any assertion or test-body change. **Renaming a suite or a `TEST`
  case is a behavior-observable change** to the pass set - see Constraints.

## Moves (exactly what)

### Directory mirror (one example per module; full mapping is 1:1 with `MODULES.md section 1`)
- `tests/dynamics/` - `RobotEngineTest.cpp`, `RobotIntegratorTest.cpp`,
  `ConstraintsTest.cpp`, `JointKernelsTest.cpp`, `BatScalingTest.cpp`,
  `NcmcProtocolTest.cpp` (the FAST half from TEST-006), the mass-matrix / kinetic
  / linear-algebra oracle tests.
- `tests/world/sampler/` - `HmcMoveTest.cpp`, `NcmcMoveTest.cpp`,
  `VelocityDistortionTest.cpp`, `CartesianSolventTest.cpp`, plus the ensemble /
  Boltzmann / equipartition / NCMC work-chain suites that exercise the sampler.
- `tests/world/` - `FixmanCorrectionTest.cpp`, `ReactionReporterTest.cpp`,
  `ModelBuilderTest.cpp`, `GeometryFitterTest.cpp`.
- `tests/bridge/` - `ForceReducerTest.cpp`, `AlchemyForceFactoryTest.cpp`,
  `MtsIntegratorTest.cpp`, and the oracle binaries from TEST-005 where they
  exercise bridge-level contracts.
- `tests/math/`, `tests/model/`, `tests/workflow/rex/` - the transform / spatial
  / linear-algebra tests, the periodic-box test, the REX acceptance-algebra test.
- `tests/support/` - the helpers now consolidated by the earlier tickets:
  `RobotFixture.hpp`, `TestPhysConstants.hpp` (TEST-002), `HarmonicBridge.hpp`
  (TEST-003), `SamplingHarness.hpp` (TEST-004), `RoboticsOracleRunners.hpp`
  (TEST-005), plus the pre-existing header-only helpers that already meet the
  discipline (`TestHelpers.hpp`, `RobotBuilders.hpp`, `HmcDriver.hpp`,
  `StatTest.hpp`, `AnalyticForceBridge.hpp`, `RoboticsOracleLoader.hpp`).

### Documentation standard (architect test-documentation section, `TESTS.md section 5.9`)
Applied to every relocated file, comment-only:
- **File level:** a `@file` block naming which module and which contracts the
  suite covers, and what it deliberately does not cover.
- **Fixture level:** what world the fixture builds and why that world is right for
  these contracts (the `RobotFixture`-derived fixtures from TEST-002).
- **Test level:** one comment per test naming the invariant it defends - a
  reference to the `ARCHITECTURE.md` invariant where one applies (e.g.
  `// verifies: INV-5 mass-metric`), Given/When/Then otherwise. Never restate the
  assertions.
- **Characterization tests labeled as such** (`TESTS.md section 3`: the integrator
  drift/reversibility chains, `TestStability`, `TestCyclicBoltzmann` pin current
  behavior, not intent).
- **Support utilities:** full Doxygen, production standard.

### CMake (`CMakeLists.txt`)
- The discovery glob `file(GLOB tests/Test*.cpp)` is flat. Change it to recurse
  the mirrored tree (`file(GLOB_RECURSE tests/*Test.cpp)` or an explicit
  per-directory add), preserving `CONFIGURE_DEPENDS`. The oracle fixture-wiring
  `foreach` (extended by TEST-005) matches by target name and is
  unaffected by the directory move. `tests/support/` holds only headers, so it is
  still not discovered as executables.

## Target structure (after this ticket)

```
tests/
  support/            # all shared fixtures/helpers/bridges/harness (header discipline)
  math/  model/  dynamics/  bridge/  world/  world/sampler/  workflow/  workflow/rex/  io/
    <ClassUnderTest>Test.cpp   # one fixture per class, mirroring src/<module>/<Class>
```

## Constraints

- **Behavior-observable identifiers are frozen.** The gtest **suite name** and
  **test name** of every case SHALL be unchanged - the pass set (`B1`) and
  assertion map (`B2`) key on `Suite.Case`, and gtest identity is independent of
  the file path. This ticket renames **files and directories**, never `TEST(...)`
  suite/case identifiers. A suite rename to match a class name is a separate,
  human-gated triage ticket, never bundled here.
- Test bodies are byte-identical to their post-TEST-006 state; only comments (the
  documentation standard) and file locations change. The comment-stripped diff of
  every relocated file SHALL be empty (mechanically proves zero code change).
- A test whose contract cannot be stated in one sentence is triage evidence
  (probably tests several things) and goes to the ticket's findings file, not
  forced into a one-liner (`TESTS.md section 5.9`, architect charter).
- Shared helpers under `tests/support/` are held to production header discipline:
  `#pragma once`, IWYU-clean, `@file` block, Doxygen on every public symbol.
- No test includes a `.cpp`. Test files respect the source size caps.
- One commit for the relocation + CMake glob update; the documentation pass is a
  separate commit (comment-only, so its comment-stripped diff is empty *and* its
  code diff is empty - pure comment addition).

## Predicted breakage

- Include paths inside relocated files that used flat `#include "TestHelpers.hpp"`
  resolve once the include directories cover `tests/support/` (they already
  resolve `tests/` via `robosample_objects` INCLUDE_DIRECTORIES,
  `CMakeLists.txt`); add `tests/support/` to that set if a relative include
  no longer resolves. This is the one predicted mechanical edit class.
- The flat `file(GLOB tests/Test*.cpp)` stops matching the moved files until it is
  made recursive; that is the predicted CMake edit. No test is dropped: the
  recursive glob discovers the same set at new paths.
- No behavioral breakage: source frozen, suite/case identifiers frozen, bodies
  byte-identical.

## Exit criteria (machine-checked; see VERIFY.md section 5)

- Build clean; **pass set identical to B1**; **assertion counts identical to B2**
  in both modes; three consecutive full-suite runs identical.
- **Coverage map (B5) not reduced.**
- The `Suite.Case` set from `--gtest_list_tests` equals the pre-ticket set
  exactly (no rename, add, or drop).
- Every relocated file lives under a `tests/<module>/` path mirroring its
  `src/<module>/` counterpart per `MODULES.md section 1`; one fixture per class under
  test.
- Comment-stripped diff of every relocated file is empty (relocation commit);
  code diff of the documentation commit is empty (comment-only).
- Every remaining test carries its contract comment; characterization tests are
  labeled; a `findings/TEST-007-findings.md` file is present (even if empty) for
  any test whose contract could not be stated in one sentence.
- **No test includes a `.cpp`.**
