# DOC-CartesianSolvent: Cartesian-solvent velocity draw / KE

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `setCartesianSolvent`, `drawSolventVelocities`, `calcSolventKE`, and
  the solvent save/restore - the Cartesian sub-integration of solvent atoms that ride
  alongside the generalized-coordinate solute during an HMC move (ARCHITECTURE section 7;
  MODULES.md W7).
- **Layer:** Domain / world (ARCHITECTURE section 2).
- **Ownership:** operates on the World's state solvent fields; borrows the solvent
  atom selection. Note OQ-6: NCMC-solvent fields live on `RobotState` (concern bleed,
  DOC-RobotState) - document usage, route the relocation question to findings.
- **Invariants (HYPOTHESES):**
  - **INV-5 mass metric** - the solvent velocity draw and `calcSolventKE` use the
    solvent masses consistently with the solute mass metric, so KE composes correctly
    into the Hamiltonian. Document the draw distribution (Maxwell-Boltzmann at the
    World temperature) and the KE definition as contract.
  - The solvent Verlet sub-step is reversible so it does not break HMC detailed
    balance; state that requirement.

## 2. Scope

- **Files:** `world/sampler/CartesianSolvent.{hpp,cpp}` (from `World.cpp`).
- **Public symbols:** `setCartesianSolvent`, `drawSolventVelocities`, `calcSolventKE`,
  save/restore.
- **Known gaps to close:** state the velocity-draw distribution and temperature, and
  the reversibility requirement of the sub-step; both are correctness-critical.

## 3. Evidence pointers (tests exercising the module)

- No strong dedicated test; TestNcmcExplicitSolvent is largely stubbed (3 permanent stubs, 6/8 skip at default tier) (TESTS.md section 3, D-T2) -
  do not treat skipped cases as evidence. TestEquipartition (KE sampling) is partial
  evidence for the draw. Record the coverage gap in findings.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; draw distribution + KE metric
  + reversibility stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-CartesianSolvent-findings.md` present; the OQ-6 field-bleed and the
  coverage gap recorded; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
