# DOC-ModelBuilder findings

Style: house `//` contract blocks per the approved calibration trio (see
DOC-World-findings.md). No code changed.

## State
world/ModelBuilder.cpp (buildModel, setRootMobility, setRootMobilities) and the
World.hpp declarations carried full contract documentation from SPLIT-W1. No
residual bare public symbol found; no edit required.

## Verified hypotheses
- buildModel is the sole writer of RobotModel; RobotModel is immutable after it
  returns (RobotModel mutated only inside buildModel/recomputeGeometry).
- Built-model post-condition (index tables bodyQIndex/bodyNQ/bodyUIndex/bodyNU,
  frames X_PF/X_BM, loop-closure DistanceConstraints) is stated as the produced
  structure, not the DSU/BFS mechanism (ModelBuilder.cpp banner + inline).
- Idempotent rebuild on root-mobility change: verified (buildModel clears its
  three append-only tables; setRootMobility(ies) re-run it). Documented at the
  World.hpp declarations (227-228) and ModelBuilder.cpp:420-432.
- Illegal-root behavior is a thrown std::invalid_argument, not silent Weld
  degradation (ModelBuilder.cpp:148-153) -- documented as contract.

## Assumed notes
None.
