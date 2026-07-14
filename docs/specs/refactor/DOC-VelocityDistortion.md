# DOC-VelocityDistortion: NMA soft-mode / BAT-drive momentum coupling

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the momentum-distortion step of `reinitialize` - NMA soft-mode
  distortion and BAT-drive coupling that bias the seeded momenta toward slow
  collective directions (ARCHITECTURE section 3, section 7; MODULES.md W8).
- **Layer:** Domain / world (ARCHITECTURE section 2).
- **Ownership:** reads NMA modes (DOC-NMA) / BAT geometry (DOC-BatScaling); writes the
  momentum vector on the state.
- **Invariants (HYPOTHESES):**
  - The distortion is applied **after** `sqrt(M)` momentum seeding (INV-5) and SHALL
    preserve the momentum distribution required for HMC detailed balance - i.e. the
    distortion enters the Hamiltonian accounting so acceptance stays correct.
    Document what the distortion multiplies/rotates and how it is accounted for; a
    distortion unaccounted in the acceptance biases the distribution (Critical).
  - It uses the NMA mode convention (mass-weighting/normalization, DOC-NMA) and the
    BAT drive convention (INV-7, DOC-BatScaling) - verify agreement at the call
    boundary; a convention mismatch is a finding.

## 2. Scope

- **Files:** `world/sampler/VelocityDistortion.{hpp,cpp}` (from `World.cpp`).
- **Public symbols:** the soft-mode and BAT-drive momentum-coupling entry points.
- **Known gaps to close:** state precisely how the distortion is reflected in the
  acceptance/Hamiltonian; this is the seam where an unaccounted transform would break
  detailed balance.

## 3. Evidence pointers (tests exercising the module)

- TestNMALinearAlgebra (mode basis, FAST), TestBatScalingJacobian (BAT drive, FAST),
  TestBiasForces - contract/algebra (TESTS.md section 2). These cover the inputs; the
  momentum-accounting itself may lack a direct oracle - record the gap in findings.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the acceptance-accounting
  contract stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-VelocityDistortion-findings.md` present; any convention-mismatch or
  coverage gap recorded; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
