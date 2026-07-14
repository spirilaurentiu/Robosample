# DOC-OutputWriter findings (2026-07-13)

## Verified
- Whole-molecule periodic imaging is applied in `writeOutputsCore` (OutputWriter.cpp)
  ONLY to the DCD scratch copy `dcdScratch_`; the sampled `coords`/`replicaCoords_`
  are never modified, so imaging is output-only and never re-enters sampling.
  Confirmed by tracing the buffers. Documented as contract.
- Cadence: CSV energy rows and DCD frames are written on production rounds at the
  run's `writeFreq` cadence (RunREX/runREX gate the call). The reaction CSV shares
  the same cadence (paired DCD frame index). Documented.
- `writeReactionRows` emits the fixed 10-column schema
  `frame,replica,body_idx,atom_idx,fx,fy,fz,tx,ty,tz`, header only when the file is
  empty, and no-ops on an empty row set (no empty header-only file). Documented.
- Output slot `idx` semantics: for runREX it is the replica index (fixed
  temperature slot); for RunREX it is the thermodynamic-state index. Both keep "a
  fixed slot is a fixed temperature." Documented.

## Coverage gap (recorded)
- No dedicated C++ unit test. The Level-1 run and the 4-replica REMC run (VERIFY B3)
  diff the emitted CSV/DCD as the behavioral oracle; those baselines are the only
  regression evidence.
