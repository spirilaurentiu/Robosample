# DOC-RobotEngine findings

## INV-4 realization-order stage contracts (documented as contract)
- Each engine method now states which cache stages it requires valid on entry
  and which it makes valid on exit, recovered from the realization sequence
  position -> velocity -> articulated-body inertias -> udot:
  - realizePosition: pre q set; post position caches valid.
  - realizeVelocity: pre realizePosition + u; post velocity caches valid.
  - factorizeArticulatedInertias: pre realizePosition; post P/PPlus/D/DI/G (q-only).
  - seedArticulatedCentrifugal: pre realizeVelocity + factorize; post abcf.
  - realizeArticulatedBodyInertias: pre position+velocity; post inertia stage.
  - calcUDot: pre inertia stage + bodyForceG/mobilityForce; post udot/A_GB.
  - mass ops (multiplyByMInv/SqrtMInv/SqrtM, calcLogDetM, calcKineticEnergy):
    pre factorize/inertia stage; same mass metric (INV-5).
  - reaction forces: pre all four stages current.

## Cross-ticket consistency (INV-4 read-side vs write-side)
- These write-side stage contracts are the source of truth for the RobotState
  accessor (read-side) contracts (DOC-RobotState, a separate agent's ticket).
  DOC-RobotState is out of this agent's scope, so I could not diff the two here.
  No internal disagreement was found within RobotEngine. Flagged for the
  Architect to cross-check DOC-RobotState against these stages; a disagreement
  is a finding in both (per ticket).

## INV-5 mass metric (documented as contract)
- multiplyBySqrtMInv (momentum seeding), calcKineticEnergy, and calcLogDetM use
  the same mass metric; multiplyBySqrtM is the exact algebraic inverse of
  multiplyBySqrtMInv on the unlocked subspace. calcLogDetM uses the same
  null-space lock as the forward-dynamics pseudo-inverse (via pseudoLogDet), so
  the Fixman kinetic term stays consistent with the dynamics.

## Reaction-force convention (documented as contract, not a bug)
- calcMobilizerReactionForces does NOT subtract applied mobility (generalized
  joint) forces, matching Simbody; they are included in the reported reaction.
  Documented as a `@note`, not flagged as a defect (it is the intended ported
  behavior).

## Symbol-location note (coordination with DOC-RobotIntegrator)
- verletStep/stepTo/checkReversibility are DECLARED in RobotEngine.hpp but
  DEFINED as templates in RobotIntegrator.hpp. Per "document at the declaration,
  never in both", their contract Doxygen lives at the RobotEngine.hpp
  declaration (satisfying DOC-RobotIntegrator); the definitions in
  RobotIntegrator.hpp keep only implementation comments.

## Notes
- No `@note Assumed:` used.
