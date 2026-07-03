# Two-robot contact world — NCMC complementarity (question 5, SECONDARY)

Brief. Terms from `00-...` §2 and `docs/specs/ncmc-explicit-solvent/`.

## 1. Position
The ballistic mixed move (`10-20-...`) is the primary line. A LOCAL soft-core NCMC
move is worth pursuing as a COMPLEMENTARY move for clashes the ballistic trajectory
cannot cross — hard overlaps where even the collapsed-`dt` Verlet diverges
(`forcesFinite`/`velocitiesSane` reject, `RobotIntegrator.hpp:157-188`). Soft-core
removes the `r⁻¹²` singularity during the switch so `dt` need not collapse
*inside the NCMC protocol*; this is orthogonal to, and composable with, the
contact world.

## 2. Recommendation (kept short)
- **RECOMMENDED:** a soft-core protocol whose Region A is restricted to the
  clashing INTERFACE atoms (not the whole robot), decoupling INTERmolecular
  nonbonded interactions and NOT annihilating INTRAmolecular electrostatics
  (decouple-not-annihilate). This targets the transition NCMC exists to pay for
  — a clash resolution — without distorting the internal landscape
  (`docs/specs/ncmc-explicit-solvent/30-region-and-protocol-policy.md`).
- **BLOCKING PRECONDITION:** Construction II currently rejects ~all inner GHMC
  substeps (`World.cpp:1928-2067`; `.claude` memory
  ncmc-explicit-solvent-acceptance). That defect SHALL be resolved (per the NCMC
  spec, not here) before local soft-core NCMC is evaluated — a soft-core protocol
  wrapped around an inner kernel that rejects everything cannot cross anything.
- **NOTE:** local soft-core inherits the Fixman/pitch-in-inner-accept correctness
  requirement (`docs/specs/ncmc-explicit-solvent/10-...` F1/F2). Under the mixed
  manifold the inner `H_λ` SHALL use the SAME `U_F` as `10-...` M2 (robot-only,
  plus the canceling constant); no new Fixman term is introduced by soft-core.

## 3. Not in scope here
Full soft-core protocol design, `λ`-schedule, and the Construction-II fix are the
subject of `docs/specs/ncmc-explicit-solvent/`; this file only records that a
LOCAL, interface-restricted, decouple-not-annihilate variant is the right
complementary move and is gated on that fix.
