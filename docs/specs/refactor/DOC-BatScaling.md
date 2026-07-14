# DOC-BatScaling: deterministic BAT position-scaling drive

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the deterministic bond/angle/torsion position-scaling drive used by
  RENE - `applyBend` + `applyStretch` + the Cartesian log-Jacobian, routed through a
  shared `readBatCoord` after the dedup (ARCHITECTURE section 7; MODULES.md, section 2).
- **Layer:** Algorithms / dynamics (ARCHITECTURE section 2).
- **Ownership:** operates on borrowed geometry; writes scaled coordinates. Stateless
  aside from its inputs.
- **Invariants (HYPOTHESES):**
  - **INV-7 BAT map/Jacobian agreement** (ARCHITECTURE section 5): `applyBatScaling` and its
    Cartesian log-Jacobian read **identical** `(r,theta)` geometry; one anchor snapshot
    per round is shared by both swap partners so the paired map is an exact involution
    (INV-9 in the BAT spec). State the map/Jacobian shared-geometry contract and the
    involution property as post-conditions; the shared `readBatCoord` is what makes
    them provably identical - document that once.
- Verify the involution against the anchor-snapshot call site in the REX round, not
  from within the scaling function alone.

## 2. Scope

- **Files:** `dynamics/BatScaling.{hpp,cpp}`.
- **Public symbols:** `applyBatScaling` and its Jacobian; private `applyBend`,
  `applyStretch`, `readBatCoord`.
- **Known gaps to close:** the god-function `applyBatScaling` (112 LOC) mixes bend +
  stretch + Jacobian - document the composite contract (what geometry it reads, what
  it writes, the log-Jacobian it returns), not the three internal steps as mechanism.

## 3. Evidence pointers (tests exercising the module)

- TestBatScalingJacobian, TestBatAnchorInvolution - FAST contract (INV-7; TESTS.md
  section 2/section 3). The anchor-involution test is the direct INV-7 involution oracle.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; INV-7 (shared geometry,
  involution, log-Jacobian) stated as contract; comment-stripped diff empty
  (VERIFY section 4).
- `findings/DOC-BatScaling-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
