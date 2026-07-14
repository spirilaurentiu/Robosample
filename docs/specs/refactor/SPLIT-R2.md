# SPLIT-R2: extract the ROBO_DEBUG NaN scanner from RobotEngine.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)

`src/RobotEngine.cpp` (1456 LOC) carries ~150 lines of a compile-gated `ROBO_DEBUG`
NaN scanner (the `robodbg` namespace) inline in the dynamics translation unit. It is
a development aid, not part of the recursions, and it is the only reason
`<iostream>`/`<iomanip>` are pulled into this file. This ticket extracts it to a
guarded header so the dynamics `.cpp` files produced by SPLIT-R3 do not each
re-embed it. It runs after SPLIT-R1 (which removed the linear-algebra kernels from
the same anonymous namespace) and before SPLIT-R3.

The scanner reads only `RobotModel`/`RobotState` cache pointers and prints; it
changes no state. The single behavioral contract that SHALL survive: when
`ROBO_DEBUG == 0` every probe compiles to nothing (zero overhead). No `ARCHITECTURE.md`
INV governs a debug printer; the applicable constraint is the zero-overhead-when-off
property the file's own banner states (lines 150-158).

## Moves (exactly what)

From `src/RobotEngine.cpp`:

- the master switches `#define ROBO_DEBUG 1` and `#define ROBO_VERBOSE 2`
  (lines 157-158, with the DEBUG-LOGGING banner 148-156) -> `include/math/robo_debug.hpp`
- the `robodbg` namespace: `evalCount`, `stepCount`, `firstNanDumped` (lines
  163-165), `fin`/`finV`/`finS` (167-175), `dumpBody` (178-201), `scan` (207-311),
  all inside `#if ROBO_DEBUG ... #endif` (lines 160-313) -> `robo_debug.hpp`
- the `ROBO_CHECK(where)` macro, both branches (lines 314-317) -> `robo_debug.hpp`

The mutable counters (`evalCount`, `stepCount`, `firstNanDumped`) are file-static
globals today (internal linkage). In a header shared by several TUs they SHALL NOT
become multiply-defined. Define them once: declare `inline` (C++17
`inline` variables have a single definition across TUs) so `RobotEngine_dynamics.cpp`
and any other includer share one instance. This preserves the current single-counter
behavior (the "first NaN" latch is process-global today via one TU).

Stays in `RobotEngine.cpp` (references, not definitions): the inline
`#if ROBO_DEBUG { ... }` diagnostic block inside `factorizeArticulatedInertias`
(lines 903-946), which calls `robodbg::fin`, `robodbg::firstNanDumped`,
`robodbg::evalCount`, `robodbg::stepCount`. After R3 this block lives in
`RobotEngine_dynamics.cpp`; that TU includes `math/robo_debug.hpp`.

`RobotEngine.cpp` gains `#include "math/robo_debug.hpp"`.

Interaction with `RobotIntegrator.hpp`: that header already defines a fallback
`ROBO_CHECK` no-op guarded by `#ifndef ROBO_CHECK` (lines 45-48) and undefines it at
end of file (636-639). Moving the real `ROBO_CHECK` into `robo_debug.hpp` does not
change that: a TU that includes `robo_debug.hpp` gets the real macro; a TU that does
not (every integrator caller) keeps the `RobotIntegrator.hpp` no-op. The ODR note in
`RobotIntegrator.hpp` (lines 26-32) - "keeps all instantiations of a given Bridge
identical across TUs" - is preserved: the integrator bodies are unchanged.

## Public API after this ticket

`math/robo_debug.hpp` exports (only when `ROBO_DEBUG != 0`) in `namespace robodbg`:
`evalCount`, `stepCount`, `firstNanDumped` (`inline` variables), `fin`, `finV`,
`finS`, `dumpBody`, `scan`; and unconditionally the `ROBO_CHECK(where)` function-like
macro. Header self-contained (IWYU): `#pragma once`, `#include <cstdio>`,
`#include <iomanip>`, `#include <iostream>`, `#include <string>`,
`#include "RobotModel.hpp"`, `#include "RobotState.hpp"`, `#include "robot_math.hpp"`.

No public exported symbol delta (the scanner was internal-linkage file-static).

## Constraints

- Pure code motion. Zero logic changes. Zero symbol renames. The only edit beyond
  relocation is `static` file-globals -> `inline` variables, required for a header
  shared by multiple TUs after R3; this preserves the single-instance semantics the
  code has today and is not a behavioral change.
- `src/RobotEngine.cpp` gains `#include "math/robo_debug.hpp"`, loses the `robodbg`
  namespace, the two `#define`s, and the `ROBO_CHECK` macro definition; the
  `factorizeArticulatedInertias` diagnostic block (903-946) stays and now resolves
  `robodbg::` through the header. It may drop `<iostream>`/`<iomanip>` if no other
  use remains - verify with IWYU before removing.
- Build system: add `src/math/robo_debug.cpp` ONLY if any symbol needs a TU;
  with `inline` variables and `inline` functions the header is complete and no `.cpp`
  is required. Header-only is preferred (matches the guarded-aid intent).
- Zero-overhead-when-off: with `ROBO_DEBUG == 0` the `#if ROBO_DEBUG` body is
  excluded and `ROBO_CHECK` expands to `((void)0)`; the compiled engine is
  byte-identical to a build with the scanner absent. This SHALL be verified by a
  release build (which does not define `ROBO_DEBUG` beyond the header default).
- One commit for the move; formatting fixes separate.

## Predicted test breakage (Phase-A include fixes only)

**None.** No test names `robodbg` or `ROBO_CHECK`; the scanner is unreachable from
tests (internal-linkage today, and the header is included only by engine TUs). The
integrator tests rely on `RobotIntegrator.hpp`'s own `#ifndef ROBO_CHECK` no-op,
which is unaffected. No include fix, no assertion-count change.

## Exit criteria (machine-checked)

- Full build passes in `ROBO_DEBUG` on and off; test suite / baseline outputs
  identical to Stage-0 (B1, B2, B3).
- Public-symbol nm diff vs B0: empty.
- include-cycle script clean; `robo_debug.hpp` adds no upward layer edge (depends on
  `RobotModel`/`RobotState`/`robot_math`, all at or below the dynamics layer).
- clang-tidy / clang-format / IWYU clean on touched files.
- Comment-stripped diff of `RobotEngine.cpp` empty except the removed `robodbg`
  block/macros and the added `#include`.
- `robo_debug.hpp` <= 300 LOC (expected ~170).
- No test file changed.
