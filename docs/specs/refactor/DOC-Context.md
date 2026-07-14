# DOC-Context: thin run orchestrator (residual)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).
Sequenced **last** in workflow (it composes the extracted services).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the residual `Context` after C1-C5 extract validation / output / sweep /
  REX - a thin orchestrator that owns the Worlds, brings up the OpenMM singleton,
  accumulates BAT anchors, and delegates the run to the REX driver
  (ARCHITECTURE section 1, section 3, section 7; MODULES.md C-residual). Organized around replica exchange:
  `R = 1` is the degenerate case of the default REMC driver.
- **Layer:** Workflow (ARCHITECTURE section 2).
- **Ownership:** owns `SystemTopology` (public member; the Python<->engine payload),
  the `vector<unique_ptr<World>>` (heap, so `World&` stays valid across sweeps), and
  `rexRng_` (ARCHITECTURE section 4). Owns nothing OpenMM-typed directly (it uses the
  singleton via bring-up).
- **Invariants (HYPOTHESES):**
  - Worlds are heap-owned so references survive across sweeps (ARCHITECTURE section 4);
    document the lifetime guarantee callers rely on.
  - `addCartesianWorld` / `addRoboticWorld` / `addDockingWorld` each build a World via
    `World::buildModel` (DOC-ModelBuilder) - document the construction contract and the
    order requirement (`initialize` after all Worlds added).
  - Defect section 6.3: two bindings (`set_separate_force_groups`,
    `set_enforce_periodic_box`) reach `OpenMMContext::get()` and discard their
    `Context&`; the resolution adds thin forwarders on `Context`. Document the
    forwarders' observable effect; record the reach-through in findings.

## 2. Scope

- **Files:** `workflow/Context.{hpp,cpp}` (residual after C1-C5).
- **Public symbols:** the constructor, `addWorld*`, `initialize`, the run entry
  delegating to the REX driver, `SystemTopology` accessor, the new
  force-group/periodic-box forwarders.
- **Known gaps to close:** document the World-add -> initialize -> run lifecycle order
  as contract; note the section 6.3 forwarders in findings.

## 3. Evidence pointers (tests exercising the module)

- TestGibbsWorlds - FAST contract (TESTS.md section 2). The Level-1 run and 4-replica REMC
  proof (VERIFY B3) exercise the full lifecycle.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the lifecycle order and World
  ownership/lifetime stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-Context-findings.md` present; section 6.3 reach-through recorded; no unmatched
  `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
