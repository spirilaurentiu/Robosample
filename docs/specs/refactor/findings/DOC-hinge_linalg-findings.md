# DOC-hinge_linalg findings

## Verified against call sites
- `invertDense` return value (locked-direction count) is consumed by
  `src/RobotEngine_dynamics.cpp:84` and used to fire the CC1/CC4 fail-loud gate
  (`:125`). Documented as contract.
- `pseudoLogDet` feeds `calcLogDetM` (`src/RobotEngine_massops.cpp:193`); the
  null-direction->0 contribution is the consistency link with `invertDense`.
  Documented (INV-5/CC1).
- `symSqrt`/`symSqrtInv` feed the mass-metric sqrt operators
  (`src/RobotEngine_massops.cpp:107,149`). Documented (INV-5).

## Contradicted hypothesis
- Ticket section 1 hypothesis: "the symmetric eigensolver's ordering/sign
  conventions are relied on downstream." Partly CONTRADICTED. `invertDense`,
  `pseudoLogDet`, `symSqrt`, `symSqrtInv` sum over all eigenpairs and are
  order-agnostic; none relies on `jacobiSymEig`'s native ordering. The only
  consumer that needs an order (NMA mode ordering) sorts the output itself
  (`tests/TestNMALinearAlgebra.cpp:167` comment: "robo_linalg::jacobiSymEig
  -> unsorted; sort to compare"). Documented `jacobiSymEig` as returning
  eigenvalues in no particular order, with a @warning that callers needing an
  order must sort. No sign convention is relied on either (magnitudes taken).

## Notes
- `tests/RobotLinearAlgebra.hpp` is a characterization mirror of this solver
  (per ticket / TESTS.md 3), not independent contract evidence; no contract
  was inferred from it.
- No `@note Assumed:` used.
