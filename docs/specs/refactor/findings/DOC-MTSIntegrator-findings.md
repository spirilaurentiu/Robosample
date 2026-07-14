# DOC-MTSIntegrator findings

## Coverage gap (required record)
No dedicated test for `MTSIntegrator`. The contract (nested r-RESPA substep
program; reversible + symplectic; multiple-of-parent substep constraint;
force-group 0-31 range) is recovered from the constructor source and the OpenMM
`CustomIntegrator` base, and from the force-group split at its one construction
site (`OpenMMSystemBuilder::build`, slow group 0 = Nonbonded/GBSA/NBFIX, fast
group 1 = bonded).

## Used vs dormant (required record)
Dormant-but-reachable configuration.
- Reachable: `OpenMMContext::setMTS(enabled, innerSubsteps)` -> built in
  `OpenMMSystemBuilder::build` when `useMTS`. Exposed to Python as
  `Context.set_mts` (`src/Context.cpp:130`, bound `src/PyBind11.cpp:670`).
- Not selected by default: `innerSubsteps < 2` disables it; no production driver
  (`python/robosample/run.py`) calls `set_mts`. The `mm.MTSIntegrator`
  references in `python/robosample/simulate_openmm.py` / `roborun.py` are
  OpenMM's own Python integrator in helper scripts, not this C++ class.
- Scope: affects only the Cartesian world's on-device MD (`integrateTrajectory`);
  torsional worlds always get the full force sum. Documented.

No `@note Assumed:` entries.
