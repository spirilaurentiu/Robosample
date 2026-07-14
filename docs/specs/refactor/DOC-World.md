# DOC-World: residual World (identity, thermo, composition)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).
Sequenced **after** all World-service tickets (it composes them).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the residual `World` after W1-W9 extract its ~12 responsibilities -
  identity, thermodynamic parameters, schedule setters, RNG ownership, and
  composition of the extracted services into one robot factorization + its sampler
  (ARCHITECTURE section 1, section 7; MODULES.md W-residual).
- **Layer:** Domain (ARCHITECTURE section 2). Owns `RobotModel`, `RobotState`, `ForceBridge`,
  `ConstraintSet`, `rng_` (ARCHITECTURE section 4).
- **Ownership:** owns the per-World objects by value (ARCHITECTURE section 4); the extracted
  services operate on them. A World is **stateless between sweeps** (INV-3).
- **Invariants (HYPOTHESES):**
  - **INV-3 stateless Worlds** (ARCHITECTURE section 5): no persistent per-replica state
    between sweeps; per-atom Ground coordinates are the sole inter-world currency.
    Document that the World's public surface carries no cross-sweep state.
  - Defect section 6.4: `World::nDof()` reaches `OpenMMContext::get().getNumDegreesOf
    Freedom()` - a hidden global not visible in the signature. Document the observable
    behavior; record the hidden coupling in findings (the thin-forwarder resolution is
    a separate ticket).

## 2. Scope

- **Files:** `world/World.{hpp,cpp}` (residual).
- **Public symbols:** identity/thermo/schedule setters, `nDof`, RNG accessor, the
  `generateSample`/`getAtomsLocationsInGround` facade forwarding to the samplers/
  fitter.
- **Known gaps to close:** the facade methods forward to extracted services - document
  each as the observable contract of the composed World, not by re-describing the
  service. Note section 6.4 in findings.

## 3. Evidence pointers (tests exercising the module)

- TestGibbsWorlds, TestBuilders - FAST contract (TESTS.md section 2). Plus every World-driving
  test as indirect evidence for the composed surface.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; INV-3 statelessness stated;
  comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-World-findings.md` present; section 6.4 hidden-global coupling recorded; no
  unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
