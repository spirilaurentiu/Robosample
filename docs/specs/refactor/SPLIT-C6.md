# SPLIT-C6: isolate the runREX coordinate-swap oracle into its own TU

Status: **human-gated (OQ-3).** Do not run without the human decision recorded in
`ARCHITECTURE.md section 9 OQ-3`. OQ-3 is DECIDED "keep as oracle," so the default action
is to isolate - not retire - `runREX`. This ticket is the producer that
`DOC-LegacyCoordSwapRex` documents; if a human instead decides to retire `runREX`,
this ticket and that DOC ticket are both dropped. Executor: `coder`. Style:
`styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)

`runREX` is the legacy coordinate-swap replica-exchange driver. Per OQ-3 it is
deliberately frozen as the INVARIANT-EQUIV differential oracle for the label-swap
`RunREX` (`Context.hpp:84-90`, "Do NOT extend this method"). SPLIT-C4 relocates it
**unchanged** into `ReplicaExchangeDriver.cpp` alongside `RunREX`, co-located with
the shared REX setup it calls. This ticket completes the `MODULES.md section 1` target
layout by moving `runREX` one step further, out of `ReplicaExchangeDriver` and
into its own `workflow/rex/LegacyCoordSwapRex.{hpp,cpp}` - a single-purpose TU
whose file name states its status: a frozen oracle, not production surface.

The move is motivated, not cosmetic: keeping the "do not extend" oracle in its own
translation unit makes the freeze visible at the file level and stops future REX
edits in `ReplicaExchangeDriver` from drifting into it.

This ticket runs **after SPLIT-C4** (which SHALL already have moved `runREX` and the
shared setup into `ReplicaExchangeDriver`). It is pure code motion.

Applicable invariant, restated in full:

- **INV-8 REX detailed balance.** `runREX` and `RunREX` produce equivalent
  sampling; `runREX` is the differential oracle for `RunREX`. Swap acceptance
  guards reject on NaN/inf. `runREX` swaps `replicaCoords_` directly and is
  INV-3-non-compliant *by construction* - the exact property that makes it the
  oracle. This move preserves the method byte-for-byte, so the INVARIANT-EQUIV
  relation with `RunREX` is unchanged.

## Moves (exactly what)

Cohesive translation-unit split (R3 pattern). `runREX` stays a `Context` member;
`include/Context.hpp` is unchanged, so the PyBind11-bound `run_rex` symbol keeps
its exact signature and the public nm delta is empty. Only the **definition**
relocates.

To `src/workflow/rex/LegacyCoordSwapRex.cpp` [Context member definition, verbatim]:

- `Context::runREX(int,int,int,bool)` including its leading comment - relocated by
  C4 into `ReplicaExchangeDriver.cpp` (originally `Context.cpp:397-535`; re-anchor
  by symbol, the C4-shifted location is authoritative). Moved BYTE-FOR-BYTE.

Stays in `ReplicaExchangeDriver.cpp` (shared by both drivers, NOT moved):
`setupReplicaExchange`, `swapThermodynamicStates`, `prepareExchangePairs`,
`mixAllReplicas`, `mixReplicas`, `printSwapMatrix`, and `RunREX`. `runREX` calls
`setupReplicaExchange` and the shared REX state members across the TU boundary as
ordinary member access - no signature change.

`RexInternal.hpp` (created by C4, holding `kBoltzmann_kJ`) is included by the new
TU; `runREX` reads that constant.

## Public API after this ticket

No new exported symbol; `include/Context.hpp` byte-unchanged. The PyBind11 binding
`run_rex -> &Context::runREX` (`PyBind11.cpp`) resolves against the unchanged
declaration - no PyBind edit. The new TU is self-contained (IWYU): it carries the
includes `runREX` uses - `<algorithm>`, `<fstream>` (the `.moves.csv` writes),
`<string>`, `<vector>` - plus `Context.hpp`, `ReplicaExchange.hpp`, and
`RexInternal.hpp`. `#pragma once` on the header if one is added; a thin header MAY
be omitted (the `.cpp` compiles against `Context.hpp`).

## Constraints

- Pure code motion. Zero logic changes. Zero symbol renames. `runREX` moves
  **unchanged** - the OQ-3 freeze is preserved exactly.
- `ReplicaExchangeDriver.cpp` loses the `runREX` definition and nothing else; the
  shared setup and `RunREX` stay. `src/workflow/rex/LegacyCoordSwapRex.cpp` gains
  it.
- Invariant that SHALL remain true: **INV-8** (restated above), and INVARIANT-EQUIV
  between `runREX` and `RunREX` - neither driver's visit order or RNG-consumption
  order changes, so the differential oracle still holds bit-for-bit.
- Build: compile the new TU into the object library; the non-recursive `src/*.cpp`
  glob (`CMakeLists.txt`) SHALL reach it (flat placement or `GLOB_RECURSE`, per
  CX-2). One commit for the move; formatting separate.

## Predicted test breakage (Phase-A include fixes only)

- None. `run_rex` runs through the unchanged Python binding;
  `tests/test_rex_label_swap_equivalence.py` drives `run_rex`/`run_rex_label_swap`
  and asserts INVARIANT-EQUIV per-round PE - bit-identical after this move, so it
  passes with no edit. No C++ test names `runREX` directly.

## Exit criteria (machine-checked)

- Full build passes; test suite / baseline outputs identical to Stage-0
  (B1/B2/B3), including `test_rex_label_swap_equivalence.py`'s per-round PE match.
- nm diff on public symbols vs B0: empty (`run_rex` keeps its exact mangled name).
- include-cycle script clean vs B4; clang-tidy / clang-format / IWYU clean on
  `LegacyCoordSwapRex.{hpp,cpp}` and `ReplicaExchangeDriver.cpp`.
- `LegacyCoordSwapRex.cpp` <= 600 LOC (actual ~140); header, if present, <= 300 LOC.
- Comment-stripped before/after diff: only the `runREX` block leaving
  `ReplicaExchangeDriver.cpp` and appearing in `LegacyCoordSwapRex.cpp`, plus the
  added `#include`s.
