# DOC-GeometryFitter: coordinate ingestion and frame re-fitting

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `setAtomsLocationsInGround` / `recomputeGeometry` - push per-atom
  Ground coordinates into a World and re-fit every rigid-body frame and the
  `X_PF`/`X_BM` transforms; and read coordinates back (ARCHITECTURE section 3, section 7; MODULES.md
  W6). This is the inter-world coordinate hand-off machinery.
- **Layer:** Domain / world (ARCHITECTURE section 2).
- **Ownership:** reads/writes the World's `RobotState`/`RobotModel` frames; borrows
  the caller's coordinate buffer.
- **Invariants (HYPOTHESES):**
  - **INV-3 stateless Worlds** (ARCHITECTURE section 5): per-atom Ground coordinates (nm) are
    the **sole** inter-world currency; a World carries no persistent per-replica state
    between sweeps. `setAtomsLocationsInGround` is the entry that re-establishes all
    per-body frames from those coordinates alone - document that it fully reconstructs
    the geometric state (no hidden carry-over).
  - After a fit, the model frames are consistent with the input coordinates; the
    round trip set->get is identity up to the fit. Verify against the Context sweep
    call sites (push -> sample -> pull, ARCHITECTURE section 3).

## 2. Scope

- **Files:** `world/GeometryFitter.{hpp,cpp}` (from `World.cpp`).
- **Public symbols:** `setAtomsLocationsInGround`, `getAtomsLocationsInGround`,
  `recomputeGeometry`.
- **Known gaps to close:** state the "reconstructs all frames from coordinates alone"
  contract (INV-3) explicitly - it is the guarantee that makes stateless Worlds
  correct, and it is convention-only in source.

## 3. Evidence pointers (tests exercising the module)

- TestAtomTransfer, TestTransfer, TestMobilizerKinematics, TestGeometry - FAST
  contract (TESTS.md section 2). The transfer tests are the direct INV-3 round-trip evidence.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; INV-3 reconstruction contract
  stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-GeometryFitter-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
