# DOC-OpenMMContext: the OpenMM adapter singleton (residual)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the process-singleton OpenMM adapter that, after O1-O6 extract
  builders/integrator/kinematics, **holds** the one `System`/`Context`/`Integrator`
  and exposes the runtime energy/force evaluation surface plus config toggles
  (ARCHITECTURE section 3, section 4, section 7).
- **Layer:** Bridge (ARCHITECTURE section 2).
- **Ownership:** process singleton via `get()`, lives the whole run; owns the OpenMM
  objects the builder handed it and (until O6) `gGpuKin` file-static CUDA state
  (ARCHITECTURE section 4).
- **Invariants (HYPOTHESES):**
  - The energy/force entry points return per-atom forces / PE in the engine's unit
    convention (INV-3, nm/kJ*mol^-^1) consumed by `ForceBridge`/`ForceReducer`
    (INV-1). Document what each entry point computes, its inputs (current `posq`),
    and completion semantics (host force path vs fused CUDA path).
  - Singleton reach-through defects (ARCHITECTURE section 6.3, section 6.4): `set_separate_force_
    groups`/`set_enforce_periodic_box` are called on `get()` from bindings, and
    `World::nDof()` reaches `getNumDegreesOfFreedom()`. Document the toggles' actual
    effect; record the reach-through coupling in findings (the thin-forwarder
    resolution is a separate ticket, not a doc claim).

## 2. Scope

- **Files:** `bridge/OpenMMContext.{hpp,cpp}` (residual after O1-O6).
- **Public symbols:** `get()`, the energy/force evaluation entry points, the config
  toggles, `getNumDegreesOfFreedom`, and remaining accessors.
- **Known gaps to close (flagged):** the **energy/force evaluation entry points are
  the least-documented runtime surface** while the config toggles are the most
  (ARCHITECTURE section 8). Invert that here: give the evaluation entry points full
  contracts (inputs, unit/frame of outputs, host-vs-CUDA path selection, when results
  are valid). Do not over-document the toggles.

## 3. Evidence pointers (tests exercising the module)

- TestAlchemy - the only real OpenMM `Context` (TESTS.md section 1); primary evidence for
  the evaluation entry points. Otherwise the Level-1 run (VERIFY B3).

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented, with the energy/force entry
  points fully contracted; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-OpenMMContext-findings.md` present; section 6.3/section 6.4 reach-through coupling
  recorded; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
