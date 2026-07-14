# DOC-DrivenRexDriver findings (2026-07-13)

## OPEN-QUESTION (OQ-5): entire module uncompiled and untested
- runDrivenRound, driveReplica, runInterleavedRemcSubround (DrivenRexDriver.cpp)
  are compiled-but-not-runtime-exercised: no run reaches them (RunREX dispatches
  RENE/REBASONTOP here, but the driven path has no end-to-end test and, per the
  coordinator directive 2026-07-12, was not compiled/run). RENEMC's round-loop is
  unimplemented (RunREX throws). There is NO runnable evidence for this module;
  that absence is its defining characteristic. Route to OPEN-QUESTIONS.

## Every @note Assumed: in the diff (all in include/Context.hpp) - matching entries
- runDrivenRound: intended round structure (equil sweep, pairing-before-drive, one
  frozen anchor snapshot per round, drive, WTerm swap, REBASONTOP interleave)
  documented from source; reviewed-on-paper, not build-confirmed. ASSUMED.
- driveReplica: intended B4 Q-scale-factor drive, WORK/WORK_Jacobian accumulation,
  and the domain_error -> WORK_Jacobian=-inf forced-reject mapping. Consuming side
  of the sentinel is unit-checked in TestRexAcceptanceAlgebra, but the try/catch
  itself is not exercised end-to-end. ASSUMED.
- runInterleavedRemcSubround: intended REMC-branch-borrow for D4 sub-rounds.
  ASSUMED.
- Related @note (uncompiled) entries carried on the trial state they touch:
  Replica WORK_* fields and commitWorkAsFinal (ReplicaExchange.hpp), Context BAT-
  anchor accumulate/snapshot/reset, RUN_TYPE RENE/REBASONTOP/RENEMC, RunREX driven
  dispatch, setInterleaveRemcEvery/setRebasontopSubrounds. All marked uncompiled.

## Un-stateable contracts
- No behavioral contract (accept-rate, detailed balance in situ, WORK correctness)
  can be stated for the driven path from evidence; only intended structure is
  documented, marked Assumed.
