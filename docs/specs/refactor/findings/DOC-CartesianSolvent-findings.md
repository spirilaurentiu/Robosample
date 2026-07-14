# DOC-CartesianSolvent findings

Style: house `//` contract blocks (see DOC-World-findings.md). No code changed.

## State
world/sampler/CartesianSolvent.cpp (setCartesianSolvent, drawSolventVelocities,
calcSolventKE) and the World.hpp declarations (238, 597-598) carried full
contract documentation from SPLIT-W7. No residual edit required.

## Verified hypotheses
- Draw distribution: Maxwell-Boltzmann v_s ~ N(0, RT/m_s) per Cartesian
  component at the World temperature (CartesianSolvent.cpp:96-107). Documented.
- KE metric (INV-5): calcSolventKE = 1/2 sum_s m_s |v_s|^2, flat diagonal metric,
  composes additively into H as keSolvent (HmcMove.cpp:421-431, 454-462).
- Reversibility: the solvent draw mirrors the generalized-momentum draw (a Gibbs
  update of the velocity marginal), so the move needs no explicit flip
  (:96-99); the Verlet sub-step is the same reversible integrator. Stated.
- PRECONDITION P1: flagged atoms must sit in a 0-DOF Weld/Rigid body covered
  EXACTLY by the mask, else KE double-count / momentum double-draw; enforced by a
  throw (:37-85). Documented as contract.

## Flagged for the Architect
- OQ-6 (concern bleed): the NCMC/solvent (x,v) fields are stored on RobotState
  (state_.setCartSolvent / cartSolventAtoms / cartSolventInvMass), not on the
  CartesianSolvent concern. Usage documented; the relocation question is routed
  here, not resolved in the doc phase.

## Coverage gap
- No strong dedicated test. TestNcmcExplicitSolvent is largely stubbed
  (TESTS.md D-T2: 3 permanent stubs, 6/8 skip at default tier) -- not treated as
  evidence. TestEquipartition is partial evidence for the draw only. The draw
  distribution / KE / reversibility contracts are recovered from code + the
  reinitialize call site, not from a passing oracle. Recorded as a gap.

## Assumed notes
None.
