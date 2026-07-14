# DOC-Replica findings (2026-07-13)

## File-location discrepancy (ticket scope vs reality)
- Ticket section 2 lists the scope file as workflow/rex/Replica.hpp. That file does not
  exist. Replica and ThermodynamicState (and the RUN_TYPE / ReplicaMixingScheme
  enums) are declared in include/ReplicaExchange.hpp. Documented there.

## Index maps are NOT co-located with the value types
- Ticket section 2 says "the index-map accessors if co-located here." They are not: the
  two inverse maps (replica2ThermoIxs_ / thermo2ReplicaIxs_) are private members of
  Context (include/Context.hpp), not fields of Replica/ThermodynamicState. The
  inverse-map invariant was documented on those Context members instead.

## Verified
- Replica holds committed coordinates + energies (INV-6: potential == energy of
  atomsLocations) and, for driven runs, a separate trial block. ThermodynamicState
  fixes a target temperature and the per-World schedule. Documented.
- potential == referencePotential for every run type (Fixman excluded
  unconditionally, D3); documented as such, not as a run-type-dependent split.

## OQ-5
- The WORK_* / referenceWORK_* trial fields and commitWorkAsFinal are only touched
  by the uncompiled driven path (OQ-5); marked uncompiled in the docs. No runtime
  evidence for their behavior beyond the inline-formula unit test.
