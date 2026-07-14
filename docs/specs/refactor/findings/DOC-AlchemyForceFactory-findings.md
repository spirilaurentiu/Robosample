# DOC-AlchemyForceFactory findings

## Control-surface names (verified per task note)
The setter is `setAlchemicalLambda` and the enable is `enableAlchemy`. Confirmed
at every call site:
- production: `NcmcMove.cpp` calls `bridge_.enableAlchemy(...)` and
  `bridge_.setAlchemicalLambda(...)` (lines 60, 317, 383, 433, 507);
- tests: `TestAlchemy.cpp` calls `omm.enableAlchemy(...)` /
  `omm.setAlchemicalLambda(...)`.
`OpenMMContext` forwards `enableAlchemy` to the factory; `setAlchemicalLambda`
stays on `OpenMMContext` (drives the live Context via `lambda_inter`) and only
queries `AlchemyForceFactory::enabled()`.

## Lambda semantics (verified against contract tests + production caller)
- `lambda_inter == 1` reproduces the unmodified field. Verified: `TestAlchemy`
  "lambda=1 fidelity" (`E == Eplain`); `NcmcMove` brackets the protocol with
  `setAlchemicalLambda(1.0)` at start and end and evaluates acceptance at full H
  at lambda=1.
- `lambda_inter == 0` fully removes A<->rest coupling. Verified: `TestAlchemy`
  "lambda=0 decoupling" (`E == Eplain - S`).
- Vacuum/implicit path is exactly linear in lambda: `TestAlchemy` asserts
  `E == Eplain + (lambda-1)*S` across {1,0.75,0.5,0.25,0}. Documented.
- Explicit-solvent path additionally annihilates intra-A electrostatics for
  lambda<1 (charge scaled against everything through the PME charge offset) - a
  documented departure from pure decoupling, exact and unbiased at lambda=1.
  Stated in the contract as the ticket requires.

## Evidence excluded (per ticket)
`TestNcmcExplicitSolvent` is largely stubbed (permanent stubs; most cases skip at
the default tier). Its skipped cases were not treated as evidence.

No `@note Assumed:` entries. No ticket hypothesis contradicted.
