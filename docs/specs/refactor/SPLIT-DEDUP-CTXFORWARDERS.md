# SPLIT-DEDUP-CTXFORWARDERS: add Context force-group / periodic-box forwarders

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

> **API-CHANGE TICKET - requires human approval before it runs.** This ticket adds
> **two public `Context` symbols**, so it is NOT pure code motion and criterion 4
> (empty public-symbol nm delta) is relaxed only under the human-approved flag
> (VERIFY section 3). It carries its own nm-delta rationale below. It changes include
> structure per ARCHITECTURE.md section 6.3.

## Context (why)
ARCHITECTURE.md section 6 defect 3: two `Context` bindings in `src/PyBind11.cpp` reach
the OpenMM singleton directly and **discard their `Context&` argument**:

- `set_separate_force_groups` (`src/PyBind11.cpp:681-687`): lambda
  `[](Context& /*ctx*/, bool enabled){ OpenMMContext::get().setSeparateForceGroups(enabled); }`.
- `set_enforce_periodic_box` (`src/PyBind11.cpp:688-698`): lambda
  `[](Context& /*ctx*/, bool enabled){ OpenMMContext::get().setEnforcePeriodicBox(enabled); }`.

Every other `Context` config toggle forwards through a `Context` method - `setMTS`
already does exactly this (`include/Context.hpp:252`, `src/Context.cpp:168-170`:
`OpenMMContext::get().setMTS(...)`). These two are the outliers: the bindings
bypass `Context` and couple the Python surface straight to the process singleton.
This ticket restores the pattern by adding two thin `Context` forwarders and
routing the bindings through them, so bindings call the object, not `get()`.

No ARCHITECTURE.md INV-1...INV-9 governs these toggles. The behavioral contract to
preserve: the two Python binding **names**, their `enabled`/`v` argument, their
docstrings, and their observable effect on the OpenMM context stay identical -
`set_enforce_periodic_box` still MUST default False under explicit solvent so the
robot engine receives whole molecules (the existing docstring,
`src/PyBind11.cpp:692-697`).

## Moves (exactly what)
Not a move - an addition plus a two-line binding rewrite.

- Add to `include/Context.hpp` (next to `setMTS`, `include/Context.hpp:252`), a
  thin forwarder pair mirroring `setMTS` verbatim in style:
  - `void setSeparateForceGroups(bool enabled);`
  - `void setEnforcePeriodicBox(bool enabled);`
- Add to `src/Context.cpp` (next to `Context::setMTS`, `src/Context.cpp:168-170`):
  - `void Context::setSeparateForceGroups(bool enabled) { OpenMMContext::get().setSeparateForceGroups(enabled); }`
  - `void Context::setEnforcePeriodicBox(bool enabled) { OpenMMContext::get().setEnforcePeriodicBox(enabled); }`

  These forward to the unchanged `OpenMMContext` setters
  (`include/OpenMMContext.hpp:94` `setSeparateForceGroups`, `:138`
  `setEnforcePeriodicBox`).
- Rewrite the two bindings in `src/PyBind11.cpp:681-698` to call the object,
  keeping the exact `.def` names, `py::arg`, and docstrings:
  - `.def("set_separate_force_groups", &Context::setSeparateForceGroups, py::arg("enabled"), "<unchanged docstring>")`
  - `.def("set_enforce_periodic_box", &Context::setEnforcePeriodicBox, py::arg("enabled"), "<unchanged docstring>")`

  The lambdas (and their `OpenMMContext::get()` reach-through) are deleted.

## Public API after this ticket
`Context` gains two public symbols: `Context::setSeparateForceGroups(bool)` and
`Context::setEnforcePeriodicBox(bool)`. The **Python** API is unchanged - the
binding names `set_separate_force_groups` / `set_enforce_periodic_box`, their
single `bool` argument, and their docstrings are identical; only the C++ callee
changes from a discarding lambda to the new `Context` method.

nm-delta rationale (criterion 4 relaxed, human-approved): the C++ public-symbol
set grows by exactly the two new `Context` methods. This is intended - it removes
a `Context`->singleton reach-through (ARCHITECTURE.md section 6.3) and aligns these
toggles with the established `setMTS` forwarder. No existing public symbol is
removed or renamed.

Include-edge change (checked against B4): `src/PyBind11.cpp` no longer needs
`OpenMMContext::get()` **for these two toggles**. It still includes
`OpenMMContext.hpp` for other uses (`OpenMMContext::get()` elsewhere in the file,
`OpenMMContext::ForceGroupEnergy`, `set_mts`), so the include is not removed - the
narrowing is one fewer reach-through, not one fewer include. `src/Context.cpp`
already includes `OpenMMContext.hpp` (via `Context.hpp`), so no new include.

## Constraints
- This is the ARCHITECTURE.md section 6.3 defect resolution: add a thin forwarding
  symbol; bindings call the object, not `get()`. It is human-approved and
  additive, never bundled into a code-motion split.
- Behavior identical: the forwarders are pure pass-throughs to the same
  `OpenMMContext` setters the lambdas called, so the OpenMM context sees the same
  calls in the same order.
- No header cycle: `Context.hpp` already includes `OpenMMContext.hpp`.
- One commit for the additive change; formatting fixes separate.

## Predicted test breakage (Phase-A include fixes only)
- None. Tests call the **Python** bindings, whose names and behavior are
  unchanged. Callers exercised: `tests/test_openmm_potential_energy.py`
  (`set_enforce_periodic_box`, `set_separate_force_groups`),
  `tests/test_ensemble_pe_ladder.py`, `tests/test_rex_label_swap_equivalence.py:81`,
  `tests/test_torsion_conformational.py`, and the drivers/validation under
  `python/robosample/` (`openmm_validation.py:379`, `run*.py`). All keep working
  unchanged; no test file is edited. `robo_bindings.pyi` (the generated stub,
  `python/robosample/robo_bindings.pyi:136,160`) already documents these binding
  names with a `bool` arg - no signature change, so the stub needs no regeneration
  for correctness (regenerate only if the doc phase requires it).

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0
  (B1/B2/B3). The Python toggles produce the same OpenMM-context behavior.
- nm diff on public symbols vs B0: exactly the two added `Context` methods, and
  nothing else - matching this ticket's declared human-approved delta. Any other
  delta fails the ticket.
- include-cycle script clean vs B4; the reach-through count from `PyBind11.cpp`
  into `OpenMMContext::get()` drops by two (these toggles); clang-tidy /
  clang-format / IWYU clean on the three touched files
  (`include/Context.hpp`, `src/Context.cpp`, `src/PyBind11.cpp`).
- Comment-stripped diff: the two added forwarder declarations, the two added
  definitions, and the two rewritten `.def` lines - nothing else.
