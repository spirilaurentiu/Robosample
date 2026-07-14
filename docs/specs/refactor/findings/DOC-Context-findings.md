# DOC-Context findings (2026-07-13)


## Style: exemplar vs spec (documenter.md section 8) - recorded
- The approved calibration exemplar include/bridge/ForceReducer.hpp documents its
  contract in plain `//` block-comment prose, not Javadoc `/** */` with ordered
  tags. documenter.md section 8 mandates Javadoc. Per "spec wins any conflict; conflict
  is a finding," every symbol in this workflow batch (Context/GibbsSweep/
  StartupValidator/OutputWriter/ReplicaExchange + rex drivers) was written in
  Javadoc. Recommend a single house-style decision so the exemplars and the
  workflow docs agree (either convert the exemplars to Javadoc, or amend section 8).

## Verified
- Lifecycle order (construct -> add*World -> [pre-init config] -> initialize ->
  run entry) recovered from Context.cpp and the PyBind11 bindings; documented as
  contract on the class and on `initialize`.
- Worlds are heap-owned (`vector<unique_ptr<World>>`); `add*World` returns
  `World&` that stays valid across sweeps. Documented as the lifetime guarantee.

## Reach-through (defect section 6.3) - recorded, not fixed
- `setSeparateForceGroups` and `setEnforcePeriodicBox` (and `setMTS`) are thin
  forwarders onto the `OpenMMContext::get()` singleton (Context.cpp:130-140).
  These two forwarder declarations are ABSENT from git HEAD's Context.hpp; they
  were added by the uncommitted SPLIT/defect-resolution work, not by the
  documenter. The reach-through is real: Context holds no OpenMM object, it
  reaches the process-wide singleton. Documented the observable effect; the
  singleton coupling is the architectural note for OPEN-QUESTIONS.

## Contradicted ticket hypothesis
- Ticket section 3 claims "TestGibbsWorlds ... exercise the full lifecycle." FALSE at call-site
  level: tests/TestGibbsWorlds.cpp never constructs a `Context`, never calls
  `add*World`/`initialize`/`RunREX`. It pins engine-level (RobotEngine/RobotState)
  internal<->Cartesian transfer idempotence only. The actual end-to-end lifecycle
  oracle is the Level-1 run and the 4-replica REMC run (VERIFY B3), not a C++ unit
  test. Documentation of the lifecycle rests on source + Python-run evidence.

## Coverage
- No C++ unit test constructs a Context. Recorded as a gap.
