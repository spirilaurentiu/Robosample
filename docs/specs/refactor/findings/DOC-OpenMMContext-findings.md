# DOC-OpenMMContext findings

## Reach-through coupling (ARCHITECTURE 6.3 / 6.4, required record)
Confirmed at call sites; recorded here, not written as a doc claim:
- `Context::setSeparateForceGroups` / `Context::setEnforcePeriodicBox`
  (`src/Context.cpp:134,138`, bound in `src/PyBind11.cpp:682,686`) call the
  toggles on `OpenMMContext::get()` directly - the binding reaches through the
  singleton.
- `World::nDof()` (`include/World.hpp:320`) reaches
  `OpenMMContext::get().getNumDegreesOfFreedom()`.
The thin-forwarder resolution is a separate (non-doc) ticket.

## Documentation emphasis (per ticket)
Inverted the pre-split imbalance: the energy/force evaluation entry points
(`computePotentialEnergy`, `computePotentialEnergyByGroup`,
`evaluateForcesFromPositionsCache`, `integrateTrajectory`, plus the fused-CUDA
`computeForcesAndEnergyOnDevice` / `reduceForcesToBodies`) now carry full
contracts (inputs, nm/kJ units per INV-3, completion semantics, `initialize()`
precondition and its throw). Config toggles kept to one-line briefs.

## Candidate dead field (minor)
Private member `OpenMMContext::numAtoms` is written in `initialize()`
(`src/OpenMMContext.cpp:22`) and never read anywhere. Likely dead. Not a doc
claim; flagged for the split/cleanup track.

No `@note Assumed:` entries. No ticket hypothesis contradicted.
