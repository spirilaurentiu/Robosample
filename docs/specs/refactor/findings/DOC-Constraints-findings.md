# DOC-Constraints findings

## INV-6 (documented as contract)
- The shared `G` guarantee: `enforceVelocityConstraints` (RATTLE),
  `enforcePositionConstraints` (SHAKE), and `calcConstraintLogDet` all assemble
  the constraint Jacobian through the single private `assembleConstraintRow`
  helper, so the Fixman correction matches the projected subspace exactly.
  Documented as a `@note` on `calcConstraintLogDet` and on
  `assembleConstraintRow`.
- Acyclic -> 0: `calcConstraintLogDet` returns exactly 0 when there are no
  loop-closure constraints. Documented as a `@return` clause. Attested by
  `tests/TestCyclicBoltzmann.cpp`, `tests/TestConstraintSolver.cpp`, and the
  Fixman caller `src/world/FixmanCorrection.cpp`.

## Private test-access (per documenter.md 4.1)
- `solveSmallSpd`/`solveCoupling` are private, exposed to `ConstraintTestAccess`
  only for unit pinning (`tests/TestConstraintSolver.cpp`). Documented as
  internal helpers with local behavior, not as an external contract surface.

## Notes
- No `@note Assumed:` used; every documented pre/post is backed by the code or a
  contract test.
