# DOC-SwapAcceptance findings (2026-07-13)

## Verified (per-run-type acceptance algebra, from SwapAcceptance.cpp)
- REMC: ETerm_equal = -(beta_H - beta_C)(refU_X - refU_Y). Documented.
- RENEMC: ETerm_nonequil, same PT form on driven-endpoint reference potentials, no
  Jacobian (INV-10 volume-preserving). Documented; acceptance path wired but its
  round-loop unimplemented (RunREX throws for RENEMC).
- RENE/REBASONTOP: WTerm = -(Work_X + Work_Y),
  Work_p = beta_target*U(x_p^tau) - beta_source*U(x_p^0) - lnJac_p. Documented.
- Guards: empty-states and RUN_TYPE::Default both throw std::logic_error;
  non-finite acceptance exponent (NaN/+-inf) forces an explicit automatic reject
  logged to stderr (INV-8 fail-loud). Documented.
- On accept, RENE/RENEMC/REBASONTOP commit both replicas' trial atomically
  (commitWorkAsFinal) before the label swap; REMC/Default do label swap only.
  Documented.

## Evidence-tier note
- tests/TestRexAcceptanceAlgebra.cpp is a contract test of the ALGEBRA: it
  reimplements the log-acceptance formulas inline (mirroring each switch branch)
  on a one-body analytic system; it does NOT instantiate Context or call
  attemptREXSwap. It defends: REMC detailed balance, the symmetric paired swap's
  blindness to a Jacobian sign flip (a pinned FINDING/limitation, not a guarantee),
  RENEMC no-Jacobian form, no-drive limit -> ETerm_equal, and the domain-error
  sentinel forcing reject. Its own header states it is "NOT compiled or run"
  (coordinator directive); the Python mirror test_rex_swap_acceptance_algebra.py
  checks the same closed forms. So the formulas are theory-derived and mirrored in
  tests, but no live-Context call site of attemptREXSwap is exercised.
