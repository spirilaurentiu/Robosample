# DOC-hinge_linalg: dense n<=6 linear-algebra solver

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the pure numeric kernels the articulated-body solver needs on small
  dense blocks - `invertDense`, `jacobiSymEig`, `symSqrt`, `pseudoLogDet` - sized for
  a single joint's `n <= 6` generalized coordinates (ARCHITECTURE section 7, MODULES.md R1).
- **Layer:** Infrastructure (ARCHITECTURE section 2). Pure functions over dense buffers; no
  engine types, no state, no allocation beyond locals.
- **Ownership:** none - borrows caller buffers, writes results into caller-provided
  output. Confirm in/out direction and whether any operate in place.
- **Invariants (HYPOTHESES):** `symSqrt` feeds the mass-metric operators (INV-5) and
  `pseudoLogDet` feeds Fixman/constraint log-dets (INV-5/INV-6); the symmetric
  eigensolver's ordering/sign conventions are relied on downstream. Document the
  numerical contract callers depend on (definiteness assumptions, singular-value
  floor for the pseudo-log-det), not the iteration mechanism.

## 2. Scope

- **Files:** `math/hinge_linalg.{hpp,cpp}` (~450 LOC, extracted from
  `RobotEngine.cpp` by `SPLIT-R1`).
- **Public symbols:** `invertDense`, `jacobiSymEig`, `symSqrt`, `pseudoLogDet`, and
  any exported dispatch/size helpers.
- **Known gaps to close:** the preconditions on matrix definiteness/symmetry are
  enforced by convention at call sites, not by the code. State each as a `@pre` only
  if every caller guarantees it; if a caller can violate it, document the actual
  behavior on violation or route to findings.

## 3. Evidence pointers (tests exercising the module)

- TestLinearAlgebra, TestLinearAlgebraOracle, TestNMALinearAlgebra - FAST contract
  (TESTS.md section 2). NOTE: `tests/RobotLinearAlgebra.hpp` re-implements this solver and
  is retargeted to it by `SPLIT-R1` (TESTS.md section 3); treat that helper as a
  characterization mirror, not an independent contract.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; contracts are behavioral.
- Comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-hinge_linalg-findings.md` present; no `@note Assumed:` without a
  matching entry.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; the
spec wins any conflict and the conflict is a finding.
