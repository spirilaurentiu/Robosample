# DOC-DockingMove: ligand repositioning move

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `repositionLigands`, sphere sampling, and `findGoodStartingPose` - the
  docking move that reseeds ligand poses (ARCHITECTURE section 7; MODULES.md W2). Gated by
  `docking_`; high independence.
- **Layer:** Domain / world (ARCHITECTURE section 2).
- **Ownership:** reads/writes ligand body transforms on the World state; uses the
  World RNG (`rng_`).
- **Invariants (HYPOTHESES):**
  - The move draws poses from a defined proposal (sphere sampling) and its acceptance
    SHALL respect the proposal's Jacobian so the move is a valid MCMC transition.
    Document the proposal distribution and how acceptance corrects for it; an
    unaccounted proposal bias is Critical.
  - `findGoodStartingPose` is initialization (pre-sampling), not part of the sampled
    chain - document it as setup, distinct from the reversible move. Verify which is
    which at the call site.

## 2. Scope

- **Files:** `world/sampler/DockingMove.{hpp,cpp}` (from `World.cpp`).
- **Public symbols:** `repositionLigands`, the sphere-sampling entry, and
  `findGoodStartingPose`.
- **Known gaps to close:** separate the reversible-move contract from the
  initialization contract; state the proposal distribution and its acceptance
  correction.

## 3. Evidence pointers (tests exercising the module)

- No dedicated test found. Record the coverage gap in findings; recover the contract
  from `generateSample`'s docking branch and the World RNG usage. Distinguish
  characterization-quality evidence from contract (`documenter.md` section 3.5).

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; proposal/acceptance vs
  initialization contracts separated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-DockingMove-findings.md` present; the no-test gap recorded; no
  unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
