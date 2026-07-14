# DOC-ForceBridge: robo:: <-> OpenMM adapter

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the sole adapter between `robo::` types and OpenMM - the boundary the
  templated integrator uses so OpenMM enters only the `World` translation unit and
  tests can substitute an analytic bridge (ARCHITECTURE section 2, section 7).
- **Layer:** Bridge (ARCHITECTURE section 2).
- **Ownership:** owned by `World` by value; an **observer** of `RobotModel` and an
  adapter to the `OpenMMContext` singleton (ARCHITECTURE section 4). Owns no OpenMM-typed
  member in its header after the section 6.1 resolution.
- **Invariants (HYPOTHESES):**
  - The `Bridge` concept is what `RobotIntegrator` is templated on (ARCHITECTURE section 2,
    section 6.1); document the interface contract each concrete bridge (OpenMM host, fused
    CUDA, analytic test bridge) SHALL satisfy - the force-evaluation call, its inputs
    (body transforms) and outputs (per-body wrenches, PE), and completion semantics.
  - The wrench output honors INV-1 (via `ForceReducer`); this bridge delegates the
    reduction, it does not re-implement it after the dedup.
  - Defect note (ARCHITECTURE section 6.1): the header SHALL no longer transitively include
    `OpenMMContext.hpp` after the forward-declaration resolution. Document the
    interface, not the include structure; but if the header still pulls the bridge
    concretely, record it as a finding.

## 2. Scope

- **Files:** `bridge/ForceBridge.{hpp,cpp}`.
- **Public symbols:** the bridge interface/concept and the OpenMM-host concrete
  implementation's public surface.
- **Known gaps to close:** the `Bridge` template concept's requirements are enforced
  by instantiation, not by a named interface. Document the semantic requirements the
  integrator relies on (`documenter.md` section 6 templates), evidenced by the three
  instantiations (OpenMM, CUDA, `AnalyticForceBridge`), not by one of them.

## 3. Evidence pointers (tests exercising the module)

- TestAnalyticForce, TestReactionForces, TestBiasForces - FAST contract (TESTS.md
  section 2/section 3). `AnalyticForceBridge` (`tests/`) is the substitute instantiation and is
  primary evidence for the concept's required surface.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the bridge-concept contract
  stated from all instantiations; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-ForceBridge-findings.md` present; the section 6.1 include state recorded; no
  unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio (ForceReducer is one of them - reuse its
INV-1 phrasing) once human-reviewed; match their form; spec wins any conflict;
conflict is a finding.
