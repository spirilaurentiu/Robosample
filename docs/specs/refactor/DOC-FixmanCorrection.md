# DOC-FixmanCorrection: Fixman potential and Jacobian corrections

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `calcFixman`, `calcLogSineSqrGamma2`, and the Cartesian log-Jacobian
  (`lnDetMCartesian_`) - the Fixman/Jacobian corrections that make GC-HMC sample the
  correct Cartesian Boltzmann distribution (ARCHITECTURE section 7; MODULES.md W5). Pure
  functions, high independence.
- **Layer:** Domain / world (ARCHITECTURE section 2).
- **Ownership:** stateless; read the model/state, return a scalar correction.
- **Invariants (HYPOTHESES):**
  - **INV-5 mass metric** (ARCHITECTURE section 5): the Fixman kinetic term is `calcLogDetM`
    (RobotEngine); the correction here composes the mass-matrix log-det with the
    constraint log-det (`calcConstraintLogDet`, INV-6) and the Cartesian Jacobian.
    Document each term's role in the total correction and its sign convention.
  - **INV-6**: `calcConstraintLogDet` returns 0 for acyclic molecules, so the
    correction reduces correctly for trees. State the composition contract.
  - Fixman Tier-0 is validated (MEMORY: Fixman campaign - Tier 0 green); the slow-tier
    trouble was test methodology, not the engine. Document the correction's contract
    from the Tier-0 evidence and the caller (HmcMove reinitialize).

## 2. Scope

- **Files:** `world/FixmanCorrection.{hpp,cpp}` (from `World.cpp`).
- **Public symbols:** `calcFixman`, `calcLogSineSqrGamma2`, the Cartesian log-Jacobian
  accessor.
- **Known gaps to close:** state the sign convention and the exact set of terms
  summed into the reported correction; a sign or missing-term error is Critical
  (biases the distribution), so the contract SHALL be unambiguous.

## 3. Evidence pointers (tests exercising the module)

- TestFixmanBoltzmann (SLOW), TestFixmanIdealizedChains (FAST), TestCyclicBoltzmann
  (SLOW), TestConstraints (log-det) - contract (INV-5/INV-6; TESTS.md section 2). Tier-0
  algebraic checks are the trustworthy evidence per the Fixman campaign.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the term composition and sign
  convention stated (INV-5/INV-6); comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-FixmanCorrection-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
