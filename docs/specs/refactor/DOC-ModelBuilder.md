# DOC-ModelBuilder: kinematic-tree construction

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `buildModel` and its helpers - turn bonds + per-joint mobilities into
  a `RobotModel`: body forest via union-find + BFS, q/u index tables, static joint
  frames, mass properties, loop-closure constraints; plus root-mobility mutation
  (ARCHITECTURE section 3, section 7; MODULES.md W1). Build-time only, high independence.
- **Layer:** Domain / world (ARCHITECTURE section 2).
- **Ownership:** the sole **writer** of `RobotModel` (immutable thereafter,
  ARCHITECTURE section 4). Reads `SystemTopology` (borrowed from Context). Rebuilds on
  root-mobility change.
- **Invariants (HYPOTHESES):** the q/u index layout, body ordering, and static frames
  it produces are what the solver's mass metric (INV-5) and constraints (INV-6)
  consume; the DSU/BFS forest defines the kinematic tree. Document the produced
  model's structure as the post-condition (what `buildModel` guarantees about the
  returned model), not the union-find mechanism. `RobotModel`'s immutability boundary
  (DOC-RobotModel) begins where `buildModel` returns.

## 2. Scope

- **Files:** `world/ModelBuilder.{hpp,cpp}` (from `World.cpp` ~601-968 + root
  mobility + DSU + frame graph).
- **Public symbols:** `buildModel`, root-mobility mutation entry, and the build
  helpers made public.
- **Known gaps to close:** state the built model's guaranteed structure (index tables
  populated, frames set, constraints recorded) as `@post`; do not narrate BFS/DSU.

## 3. Evidence pointers (tests exercising the module)

- TestBuilders, TestGibbsWorlds, TestRoboticsOracleMolecule - FAST/FAST\* contract
  (TESTS.md section 2). The molecule oracle builds a model per fixture and is strong
  structural evidence.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the built-model post-condition
  stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-ModelBuilder-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
