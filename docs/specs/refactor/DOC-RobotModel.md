# DOC-RobotModel: immutable kinematic-tree topology

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the immutable articulated-robot topology produced by `buildModel` -
  the body forest, q/u index tables, static joint frames, mass properties, and
  loop-closure constraint records; plus joint-fact predicates the solver queries
  (ARCHITECTURE section 4, section 7; MODULES.md section 5: one concept, no split).
- **Layer:** Domain data / model (ARCHITECTURE section 2).
- **Ownership:** owned by `World` by value; **immutable after `buildModel`** and
  rebuilt only on root-mobility change (ARCHITECTURE section 4). Every solver function takes
  it as `const RobotModel&`. This const-after-build contract is central - recover it
  from `ModelBuilder` (the sole writer) and the solver (read-only consumers).
- **Invariants (HYPOTHESES):** the q/u index tables and body ordering define the
  generalized-coordinate layout the mass metric (INV-5) and constraints (INV-6)
  depend on; the joint-fact predicates classify joints consistently with
  `JointKernels` (INV related to section 6.6 taxonomy). Verify predicate meanings against
  the solver call sites.

## 2. Scope

- **Files:** `model/RobotModel.hpp`.
- **Public symbols:** the model type, its index-table/frame/mass accessors, and the
  joint-fact predicates.
- **Known gaps to close:** state the immutability boundary explicitly (what
  `buildModel` may set vs what is read-only thereafter). Do not document construction
  mechanism - that is `ModelBuilder`'s ticket.

## 3. Evidence pointers (tests exercising the module)

- TestBuilders, TestGibbsWorlds, TestRoboticsOracleMolecule - FAST/FAST\* contract
  (TESTS.md section 2). The molecule oracle builds a full model per fixture.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; the immutable-after-build
  contract stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-RobotModel-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
