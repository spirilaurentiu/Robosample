# DOC-ReplicaExchangeDriver findings (2026-07-13)

## Verified
- Round structure recovered from ReplicaExchangeDriver.cpp: setup (identity maps,
  seeded replicas, docking pose-repair kick) -> per-round Gibbs sweep iterated by
  thermodynamic-state index -> refresh committed potentials (INV-6) -> mixReplicas
  (neighbour, parity from a dedicated exchangeRound_ counter) -> attemptREXSwap
  (label swap only) -> state-indexed output. Documented.
- Inverse-map invariant (replica2ThermoIxs_ / thermo2ReplicaIxs_ mutually inverse
  at all times) documented on the member declarations in Context.hpp.

## OQ-5 boundary (recorded)
- Runtime-verified surface is REMC + Default only. RENE/REBASONTOP dispatch to the
  uncompiled runDrivenRound (OQ-5). RENEMC throws std::logic_error in RunREX (its
  driven round-loop is unimplemented). Documented the boundary explicitly; did not
  write confident behavioral contracts for the driven path.

## Pending behavioral oracle (recorded gap)
- INV-8 evidence today is algebraic (attemptREXSwap acceptance formulas, tested)
  plus the label-swap-vs-coordinate-swap equivalence (Python
  test_rex_label_swap_equivalence.py). The end-to-end stationary-distribution /
  detailed-balance oracle for the full chain is a pending separate track
  (TESTS.md section 7). Documentation cites INV-8 as equivalence-to-runREX, NOT as a proven
  per-replica Boltzmann guarantee.
