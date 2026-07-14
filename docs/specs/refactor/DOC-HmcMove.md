# DOC-HmcMove: the HMC move core (reinitialize + metropolis)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).
Sequenced **last** among World tickets (depends on the extracted services).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the HMC move core - `reinitialize` (seed momenta from `sqrt(M)`,
  optionally distort, draw solvent velocities, enforce velocity constraints, evaluate
  forces, assemble the starting Hamiltonian), `metropolis` (accept/reject on total
  energy, restore on reject), and the `generateSample` dispatch that collapses to a
  thin `Move`-strategy selector once W2/W4/W8 land (ARCHITECTURE section 3, section 7; MODULES.md
  W9).
- **Layer:** Domain / world (ARCHITECTURE section 2).
- **Ownership:** orchestrates the extracted services (VelocityDistortion,
  CartesianSolvent, FixmanCorrection, DockingMove, NcmcMove) and the integrator over
  the World state; uses the World RNG.
- **Invariants (HYPOTHESES):**
  - **INV-5** - momentum seeding uses `multiplyBySqrtM`; the KE that enters the
    Hamiltonian uses the same mass metric.
  - Detailed balance: `metropolis` accepts on the **total** energy change
    (PE + KE + Fixman + Jacobian), and on reject restores the saved `q`/positions
    exactly. Document the exact set of terms in the starting/ending Hamiltonian and
    the restore contract - an omitted term or an incomplete restore biases the
    distribution (Critical).
  - The starting Hamiltonian assembly ties together INV-5 (mass), INV-6 (constraint
    log-det), INV-7/INV-9 as applicable - document which terms are summed, citing the
    contributing modules, without re-deriving them.

## 2. Scope

- **Files:** `world/sampler/HmcMove.{hpp,cpp}` (from `World.cpp`:
  `reinitialize`/`metropolis`/`generateSample` branches).
- **Public symbols:** `reinitialize`, `metropolis`, `generateSample`, the `Move`
  dispatch.
- **Known gaps to close:** enumerate the exact Hamiltonian terms and the reject-restore
  set as contract; this is the acceptance seam and the highest-consequence contract in
  the World layer.

## 3. Evidence pointers (tests exercising the module)

- TestIntegrator (characterization), TestEquipartition, TestEnsembleValidation (SLOW),
  TestMassScaleInvariance (SLOW), TestKineticEnergy, TestMassMatrix - contract/char
  (INV-5, detailed balance; TESTS.md section 2/section 3). The ensemble/Boltzmann suite is the
  detailed-balance oracle.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the Hamiltonian-term set and
  reject-restore contract stated (INV-5, detailed balance); comment-stripped diff
  empty (VERIFY section 4).
- `findings/DOC-HmcMove-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
