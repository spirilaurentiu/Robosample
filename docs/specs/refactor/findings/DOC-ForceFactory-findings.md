# DOC-ForceFactory findings

- Force-group placement: the `create*Force` builders do **not** set a force
  group. Verified against the call site: `OpenMMSystemBuilder::build`'s `addForce`
  lambda assigns the group (MTS slow=0/fast=1; else one-per-force under
  `separateForceGroups`; else group 0) after construction. The ticket hypothesis
  ("document each builder's force-group placement") is corrected: the placement
  is the builder-caller's contract, documented on the namespace and on
  `OpenMMSystemBuilder`. Documented the ownership-transfer point instead
  (`System::addForce` takes ownership).
- Exclusion policy: `addStandardExclusions` mirrors the main NonbondedForce's
  1-2/1-3 exclusions and 1-4 scaled pairs onto a custom force to satisfy the CPU
  platform's shared-neighbor-list rule; energy-neutral on every platform.
  Documented as contract.
- Naming quirk (not a bug): `createUreyBradleyForce` returns an
  `OpenMM::HarmonicBondForce` and `OpenMMSystemBuilder` labels it
  `"HarmonicBondForce"`, so two distinct force objects share that label in the
  force-group breakdown. Noted; behavior-preserving, not changed.

No `@note Assumed:` entries.
