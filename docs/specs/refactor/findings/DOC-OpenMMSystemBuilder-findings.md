# DOC-OpenMMSystemBuilder findings

- Construction order verified against the source call order in
  `OpenMMSystemBuilder::build`: particles -> periodic box (before Context, for
  PME) -> virtual sites (before Context) + massless-real-particle safety net ->
  GBSA usability gating -> force-group assignment -> Forces (standard
  ForceFactory + AlchemyForceFactory) -> integrator (MTS or Verlet) -> platform
  -> Context. Documented as the ordering contract; matches the hypothesis.
- Ownership hand-off: `build` returns an owning `OpenMMSystemBuildResult`;
  `OpenMMContext::initialize` (`src/OpenMMContext.cpp:15-31`) moves each member
  into the singleton. Documented on the result struct.
- Partial-write-on-failure: `system`/`integrator`/`hasVirtualSites`/`numAtoms`/
  `forceGroupLabels` are populated even when `success == false` (built before the
  Context try/catch); `context`/`openMMVersion` valid only on success.
  Documented as a struct `@note`.
- Throw conditions documented from the source: non-9 box vectors under a
  periodic method, cutoff > half min box width, massless non-virtual-site
  particle, NCMC alchemy + NBFIX, >32 separate force groups. Context
  construction failure is reported via `.success`, not thrown.

No `@note Assumed:` entries.
