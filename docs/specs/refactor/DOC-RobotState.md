# DOC-RobotState: per-step mutable solver cache (stage-validity contract)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).
**Highest-risk API in the doc set** (ARCHITECTURE section 8).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the per-step mutable cache the articulated-body solver reads and
  writes - ~35 raw cache pointers over one `MemoryArena` slab, exposed through
  single-letter accessors (`P()`, `Z()`, `G()`, ...) (ARCHITECTURE section 4, section 7, section 8;
  MODULES.md section 5: documented, not split).
- **Layer:** Domain data / model (ARCHITECTURE section 2).
- **Ownership:** owned by `World` by value; hands out **borrowed** raw cache pointers
  backed by its arena (see DOC-MemoryArena). Pointers are valid only while the state
  lives; document the lifetime.
- **Invariants (HYPOTHESES) - the core of this ticket:**
  - **INV-4 realization order** (ARCHITECTURE section 5, section 8; *source-only*, highest risk).
    Caches are valid **only** in stage order: position -> velocity ->
    articulated-body inertias -> udot. Reading a cache before its stage is computed is
    undefined. This contract is enforced by convention, not by types.

## 2. Scope

- **Files:** `model/RobotState.hpp`.
- **Public symbols:** every single-letter cache accessor and the state's
  construct/reset surface.
- **Known gaps to close (flagged):**
  1. **Stage-validity contract on the single-letter accessors (INV-4).** This is the
     genuine gap ARCHITECTURE section 8 names. For **each** accessor (`P()`, `Z()`, `G()`,
     ...) recover, from the solver's realization sequence
     (`realizePosition`/`realizeVelocity`/`realizeArticulatedBodyInertias`/`calcUDot`
     in RobotEngine/RobotIntegrator), the **earliest stage** at which its cache is
     valid, and state it as a `@pre`/`@warning`. The stage assignment SHALL come from
     the write sites, not from the accessor name.
  2. **NCMC-solvent field bleed (OQ-6).** `RobotState` carries NCMC-solvent fields
     (concern bleed, ARCHITECTURE section 7, OQ-6). Document them as they are used; the
     relocate-vs-accept question is a human decision -> route to findings, do not
     assert a target design.
- **OQ-1 boundary.** Whether a type-level stage guard should replace the convention
  is out of doc scope (behavioral change; ARCHITECTURE OQ-1). Document the convention
  as it exists; note OQ-1 in findings, do not implement or presuppose it.

## 3. Evidence pointers (tests exercising the module)

- TestKineticEnergy, TestMassMatrix - FAST contract (INV-5; TESTS.md section 2/section 3). These
  read caches at defined stages and are the strongest evidence for the stage each
  accessor requires. TestRoboticsOracle exercises per-body cache outputs across many
  cases (FAST\*).

## 4. Exit criteria

- Doxygen warning-free; **every** single-letter accessor carries an evidence-backed
  stage-validity contract (INV-4); comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-RobotState-findings.md` present; records OQ-1 (convention vs type
  guard), OQ-6 (NCMC field bleed), and any accessor whose valid stage the call sites
  do not agree on. No `@note Assumed:` without a matching entry.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
