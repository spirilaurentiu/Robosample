# TEST-002: Extract the shared build prologue into `tests/support/`

Status: draft, 2026-07-12. Phase-B test ticket (source frozen). Executor: `coder`. Style: `styles/spec.md`. Companion to
[`../../architecture/TESTS.md`](../../architecture/TESTS.md) section 4 (duplication) and
section 5.3/section 5.6 (the seams this executes), and [`VERIFY.md`](VERIFY.md) section 5 (exit
criteria). Test-only motion; touches no engine source.

## Context (why)

Every one of the 47 `Test*.cpp` files rebuilds `RobotModel` + `RobotState` +
force bridge in each test body - `TESTS.md section 4` records **zero gtest fixtures
across all 47 files**, the primary copy-paste source. Two concrete duplications
this ticket removes:

- The `RobotModel` -> `RobotState::allocateFull` -> `AnalyticForceBridge`
  construction prologue repeated across ~20 files. The canonical shape lives in
  `tests/RobotBuilders.hpp` (`buildForest`/`buildChain`/`buildSingle`/
  `buildBentTorsionChain`, all returning a `RobotModel`) followed by a
  hand-written `s.allocateFull(m); RobotEngine::realizePosition(m, s); ...`
  block and an `AnalyticForceBridge bridge(m, s, k)` construction
  (`tests/AnalyticForceBridge.hpp:44`). `TestEquipartition.cpp:48-50` and
  `TestFixmanBoltzmann.cpp:146` are representative call sites.
- The `kT300` / `kB` literal `0.0083144626 * 300.0`, duplicated across **9
  files** (`TESTS.md section 4` counts 10): `TestFixmanBoltzmann.cpp:61`,
  `TestMassScaleInvariance.cpp:46`, `TestEnsembleValidation.cpp:44,46`,
  `TestEquipartition.cpp:43,44`, `TestCyclicBoltzmann.cpp:53`,
  `TestEnsembleOrientation.cpp:80`, `TestNcmcTeleport.cpp:44`,
  `TestNcmcExplicitSolvent.cpp:72`, `TestTwoRobotContact.cpp` (uses `kT300`).

This ticket introduces the `tests/support/` module (`TESTS.md section 5.9`: shared
helpers held to production header discipline) with two headers: a `RobotFixture`
gtest base parameterized by a builder lambda, and a `TestPhysConstants.hpp`
holding the temperature literals.

## Scope

- **In:** two new headers under `tests/support/`; conversion of the test bodies
  that currently open-code the prologue or the literal to use them.
- **Out:** any engine source; `StatTest.hpp`'s slow gate (that is TEST-001); the
  six ad-hoc bridge subclasses (TEST-003); the sampling loops (TEST-004). This
  ticket moves the *construction* prologue and the *constant*, nothing else.

## Moves (exactly what)

### New: `tests/support/TestPhysConstants.hpp`
- One header-only `namespace rtest::phys` defining `kB = 0.0083144626`
  (kJ/mol/K), `kT300 = kB * 300.0`, and named reference temperatures the
  statistical tests already use (`TestEnsembleValidation.cpp:210` uses `T1`/`T2`
  for the two-temperature ratio; `TestEquipartition.cpp:70` sweeps `T`). Values
  SHALL be byte-identical to the literals they replace.
- Each of the 9 files above drops its local `constexpr double kB`/`kT300`
  (`TestFixmanBoltzmann.cpp:61` etc.) and includes this header.

### New: `tests/support/RobotFixture.hpp`
- A gtest fixture base, `rtest::RobotFixture`, that owns a `RobotModel m`, a
  `RobotState s`, a `robo::ConstraintSet cs`, and an
  `std::optional<AnalyticForceBridge> bridge`. Construction is driven by a
  builder callable the derived fixture supplies (a `std::function<RobotModel(Rng&)>`
  or a protected virtual), mirroring the existing lambda builders in
  `RobotBuilders.hpp`. `SetUp()` runs the fixed prologue exactly as the test
  bodies do today: build the model, `s.allocateFull(m)`, seed a valid config
  (`rtest::randomizeState`, `RobotBuilders.hpp:195`), `RobotEngine::realizePosition`
  + `fillAtomPositionsFromBodies`, then construct the `AnalyticForceBridge` with
  a caller-chosen stiffness `k` (default matching current call sites).
- The fixture holds a `Rng` seeded from a per-fixture constant so the derived
  suites keep their existing fixed seeds (`TestEquipartition.cpp:50` uses
  `0x9001`; each converted suite passes its own seed to preserve bit-identical
  streams).

## Target structure (after this ticket)

```
tests/support/
  TestPhysConstants.hpp     # kB, kT300, reference temperatures
  RobotFixture.hpp          # RobotModel+RobotState+ConstraintSet+AnalyticForceBridge base
```

`tests/support/` headers are header-only (no `.cpp`), so the CMake auto-glob
`file(GLOB tests/Test*.cpp)` (`CMakeLists.txt`) does not pick them up as test
executables and **no CMake change is required** - converted `Test*.cpp` files
`#include "support/RobotFixture.hpp"` and the existing include dirs resolve it.

## Constraints

- Pure test motion. The migrated prologue SHALL execute the **same engine calls
  in the same order** with the **same seeds and same stiffness `k`**, so every
  converted test draws a bit-identical sample stream. A fixture that changes the
  RNG seed, the realize order, or `k` is not a pure move - revert.
- No test includes a `.cpp` (`VERIFY.md section 5`).
- Headers held to production discipline: `#pragma once`, IWYU-clean, full
  Doxygen `@file` + fixture-rationale block (`TESTS.md section 5.9`, architect
  test-documentation section).
- Fixture and constant only. Bridges, loops, and gates are out of scope and
  SHALL NOT change here.
- One commit for the header addition + conversions; formatting fixes separate.

## Predicted breakage

None beyond the intended conversions. Source is frozen; the fixture reproduces
the existing prologue call-for-call, so the pass set and per-test assertion
counts are unchanged. Any converted suite whose assertion count moves indicates
the prologue was not reproduced faithfully and SHALL be reverted, not patched.

## Exit criteria (machine-checked; see VERIFY.md section 5)

- Build clean; **pass set identical to B1**; **test + assertion counts identical
  to B2** in both gated modes.
- **Coverage map (B5) not reduced.**
- **No test includes a `.cpp`.** Three consecutive full-suite runs identical.
- The 9 local `kB`/`kT300` definitions are gone; one `TestPhysConstants.hpp`
  remains and every former call site resolves to it.
- Every converted suite carries its fixture-rationale comment and per-test
  contract comment (`TESTS.md section 5.9`).
- `tests/support/RobotFixture.hpp` <= 300 LOC; `TestPhysConstants.hpp` <= 100 LOC.
