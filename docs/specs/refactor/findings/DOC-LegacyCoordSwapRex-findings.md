# DOC-LegacyCoordSwapRex findings (2026-07-13)

## SPLIT-C6 confirmed
- runREX lives in workflow/rex/LegacyCoordSwapRex.cpp (SPLIT-C6 ran). Declaration
  in include/Context.hpp; documented there (single doc location).

## Verified
- runREX is the COORDINATE-swap driver: on an accepted exchange it swaps
  replicaCoords_[r] <-> replicaCoords_[r+1] directly (INV-3 non-compliant by
  construction). Exchange accept test is the REMC criterion
  (beta_a - beta_b)(E_a - E_b) >= 0 or u < exp(delta). Documented.
- Frozen-oracle "do not extend" status (INV-8 differential baseline for RunREX)
  documented prominently as @warning, per ticket.

## Observations (not bugs)
- runREX writes a per-move telemetry CSV (baseName.moves.csv); RunREX does not.
  A behavioral asymmetry between the two drivers' side outputs, harmless to the
  INVARIANT-EQUIV per-state energy/DCD comparison (moves.csv is not part of it).
- Only Python test_rex_label_swap_equivalence.py exercises the oracle relationship;
  no C++ test drives runREX. Recorded coverage note.
