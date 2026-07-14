# DOC-NMA: Route-B normal-mode analysis

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** header-only Route-B normal-mode analysis over `RobotModel`/
  `RobotState` - `NMA::RouteBNMA` and its soft-mode products, used by velocity
  distortion (ARCHITECTURE section 4, section 7; MODULES.md section 5: single-concept, no split).
- **Layer:** Algorithms / dynamics (ARCHITECTURE section 2).
- **Ownership:** stateless free functions (`RouteBNMA` owns no state, ARCHITECTURE
  section 4); read the model/state, return modes.
- **Invariants (HYPOTHESES):** the soft modes feed HMC momentum distortion
  (VelocityDistortion); the mode basis and eigen-ordering convention are what the
  distortion relies on. Document the returned modes' definition (mass-weighting,
  ordering, normalization) as contract; verify against the VelocityDistortion caller.

## 2. Scope

- **Files:** `dynamics/NMA.hpp` (header-only; moved, not split).
- **Public symbols:** `RouteBNMA` and the mode-analysis surface.
- **Known gaps to close:** state the mass-weighting/normalization convention of the
  returned modes precisely - a distortion that assumes a different convention would
  bias momenta; this is contract, not detail.

## 3. Evidence pointers (tests exercising the module)

- TestNMALinearAlgebra - FAST contract (TESTS.md section 2). Uses the shared dense solver
  (see DOC-hinge_linalg).

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; mode convention stated;
  comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-NMA-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
