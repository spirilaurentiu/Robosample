# DOC-GibbsSweep findings (2026-07-13)

## Scope of the dedup is narrower than the ticket framed it
- Ticket section 1 describes GibbsSweep as "the one 'run a sweep over the schedule'
  primitive that deduplicates the sweep loop written three times." What was
  actually extracted into `robo::gibbs::stepWorldInGround` (GibbsSweep.cpp) is the
  INNER per-World step only: `setAtomsLocationsInGround -> generateSample ->
  getAtomsLocationsInGround`. The three surrounding sweep LOOPS are NOT identical
  and were NOT collapsed into one function:
  - `runREX` (LegacyCoordSwapRex.cpp): iterates by replica index `r`, sets
    `temperatures_[r]`, writes a `.moves.csv` telemetry row.
  - `RunREX` REMC branch (ReplicaExchangeDriver.cpp): iterates by thermodynamic
    STATE (`thermo2ReplicaIxs_[k]`), applies per-state
    setTimeStep/setMdSteps/setAcceptRejectMode each position, no moves.csv.
  - `runDrivenRound` equilibrium segment (DrivenRexDriver.cpp): iterates by state,
    SKIPS driven worlds, and calls `accumulateBatAnchorStats`.
  The shared, byte-identical contract is exactly `stepWorldInGround`; the loops
  differ in iteration order, per-World runtime setters, telemetry, and driven-world
  skipping. Documented `stepWorldInGround`'s contract precisely (order, INV-3
  currency, no cross-sweep World state); did NOT claim the three loops are
  behaviorally identical, because they are not.

## Contradicted call-site hypothesis
- Ticket section 3: "TestGibbsWorlds ... exercise the sweep end to end." TestGibbsWorlds
  does not call `stepWorldInGround` (it operates on RobotEngine/RobotState, no
  World, no Context). The sweep primitive's only runtime evidence is the Level-1
  and REMC runs.
