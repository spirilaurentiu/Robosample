# SPLIT-C1: extract the StartupValidator from Context.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/Context.cpp` (1335 LOC) mixes World lifecycle, OpenMM bring-up, startup
geometry validation, Gibbs scheduling (inlined three times), two REX drivers, and
output I/O. This ticket extracts the first standalone concern: the pre-sampling
startup health check `Context::checkStartupGeometry` - the O(N^2) clash / NaN /
potential-energy scan that refuses to start (or, under `ROBO_ALLOW_BAD_START`,
warns) on an unminimized or clashing input structure.

The check is a self-contained diagnostic run once inside `Context::initialize`
(call site `src/Context.cpp:197`). It reads only `systemTopology` fields, calls
the Context member `calcOpenMMPotentialEnergy()` once, and uses the single-source
minimum-image helper `robo::pbc::minimumImage` (`PeriodicBox.hpp`) plus
`OpenMMContext::isPeriodic`. It writes one `stderr` line and either returns or
throws `std::runtime_error`.

No ARCHITECTURE.md invariant INV-1...INV-9 governs the clash scan directly. The
behavioral contract to preserve verbatim is: the exact clash count, the
minimum-image distance used under a periodic box, the `stderr` diagnostic line
(same format string, `src/Context.cpp:341-348`), the identical throw/no-throw
decision (`anyNaN`, `nClash`, `peBad` thresholds `kClashNm = 0.08`,
`pe > 1.0e4`), and the `ROBO_ALLOW_BAD_START` downgrade. The one ambient
workflow invariant that SHALL remain true is INV-3.

- **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
  between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
  currency. (The startup check reads `systemTopology` coordinates only; it starts
  no sweep and holds no per-World state.)

## Moves (exactly what)
This is a cohesive translation-unit split (the R3 pattern already blessed for
`RobotEngine_*.cpp`). The method stays a `Context` member; only its **definition**
relocates. `Context`'s class declaration in `include/Context.hpp` is unchanged, so
no public or private symbol is added or removed.

- `Context::checkStartupGeometry()` definition - the leading comment block
  (`src/Context.cpp:232-236`) and the whole method body
  (`src/Context.cpp:237-395`) -> `src/workflow/StartupValidator.cpp` [Context
  member definition; internal linkage of the TU otherwise]. Moves **verbatim**:
  the `pairKey` lambda, the excluded-pair set assembly, the `isVirtual` scan, the
  `minImage` lambda that forwards to `robo::pbc::minimumImage`, the O(N^2) scan,
  the PE evaluation via `calcOpenMMPotentialEnergy()`, the message builder, and
  the `ROBO_ALLOW_BAD_START` branch.

The declaration stays at `include/Context.hpp:293` (`void checkStartupGeometry();`,
private). The call site `checkStartupGeometry();` at `src/Context.cpp:197` inside
`Context::initialize` is unchanged.

NOTE (deferred, not this ticket): promoting the body to a standalone free
function `robo::startup::check(const SystemTopology&, double initialPE)` would
require injecting the potential energy across the new boundary (the body calls
`calcOpenMMPotentialEnergy()` at `src/Context.cpp:338`). That parameter change is
a restructuring, not pure motion, and belongs to a separate human-approved API
ticket. This ticket keeps the body token-identical.

## Public API after this ticket
No new exported symbol. `include/Context.hpp` is byte-unchanged. If a thin header
`src/workflow/StartupValidator.hpp` is created, it exports nothing (the TU needs
no forward declarations beyond what `Context.hpp` already provides) and MAY be
omitted entirely - the definition compiles against `#include "Context.hpp"`,
`#include "OpenMMContext.hpp"`, and `#include "PeriodicBox.hpp"`. The new `.cpp`
is self-contained (IWYU): it carries the includes the moved body uses -
`<algorithm>` (`std::min`/`std::max`/`std::swap`), `<cmath>` (`std::isfinite`,
`std::sqrt`), `<cstdio>` (`std::fprintf`/`std::snprintf`), `<cstdlib>`
(`std::getenv`), `<limits>` if referenced, `<stdexcept>` (`std::runtime_error`),
`<string>`, `<unordered_set>`, `<utility>`, `<vector>` - plus `Context.hpp`,
`OpenMMContext.hpp`, `PeriodicBox.hpp`. `#pragma once` if a header is added.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/Context.cpp` loses lines 232-395; it keeps every include it still needs for
  the remaining code (it still uses `<unordered_set>` nowhere else - that include
  MAY migrate to the new TU). Nothing else in `Context.cpp` changes.
- `Context`'s class interface (`include/Context.hpp`) is unchanged; the method
  remains a private member, so the public-symbol nm set is untouched.
- Invariants that SHALL remain true: INV-3 (restated above). No INV-1...INV-9 logic
  is touched; the throw/no-throw decision and the `stderr` line are bit-identical.
- Build: the new `src/workflow/StartupValidator.cpp` SHALL be compiled into the
  object library. The source glob at `CMakeLists.txt` is non-recursive
  (`src/*.cpp`); either place the file directly under `src/` or change the glob to
  `GLOB_RECURSE`, then reconfigure. One commit for the move; formatting separate.

## Predicted test breakage (Phase-A include fixes only)
- None. `checkStartupGeometry` is private and named by no test (grep of `tests/`:
  no hit). The only C++ test that `#include "Context.hpp"` is
  `tests/TestRoboticsOracleMolecule.cpp`; it uses public construction/model APIs
  and is unaffected because `Context.hpp` does not change. The startup check is
  exercised indirectly through `initialize()` by the Level-1 example and the
  Python suite via the Python binding; behavior is bit-identical, so B1/B3 hold
  with no test edit.

## Exit criteria (machine-checked)
- Full build passes (all reference configs); test suite / baseline outputs
  identical to Stage-0 (B1/B2/B3), including the `[init] startup geometry check`
  `stderr` line byte-for-byte.
- nm diff on public symbols vs B0: empty (`checkStartupGeometry` was never a
  public symbol; internal churn allowed).
- include-cycle script clean vs B4; clang-tidy / clang-format / IWYU clean on the
  two touched files.
- `StartupValidator.cpp` <= 600 LOC (actual ~165). `Context.cpp` shrinks by ~165
  LOC.
- Comment-stripped before/after diff: the only non-empty hunks are the moved
  method at its new location and the added `#include` in the new TU.
