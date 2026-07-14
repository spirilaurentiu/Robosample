# DOC-HmcMove findings

Style: house `//` contract blocks (see DOC-World-findings.md). No code changed.

## State
world/sampler/HmcMove.cpp (reinitialize, currentTotalEnergy, metropolis,
generateSample) and the World.hpp private declarations (593-595) carried full
contract documentation from SPLIT-W9. No residual edit required.

## Verified hypotheses (acceptance seam -- highest-consequence, Critical if wrong)
- Hamiltonian term set is explicit and identical at start and end:
  H = pe + ke + keSolvent + fixman - 0.5*RT*logSineSqr - nmaCorr
  (reinitialize HmcMove.cpp:431; currentTotalEnergy :462). INV-5 (ke via the mass
  metric), INV-6 (fixman folds the constraint log-det), the NMA correction, the
  Cartesian-solvent KE, and the default-OFF orientation Jacobian are all present.
  Documented in-place; each off-path term is a strict 0 when its feature is off,
  so H matches plain HMC bit-for-bit off those paths.
- Momentum seeding (INV-5): u = sqrt(RT) * multiplyBySqrtMInv(g), g ~ N(0,I);
  never forms M/M^-1 (:306-385). The KE entering H uses the same metric.
- metropolis accepts on the TOTAL dH = Hnew - Hold_ (:466-481); equilPhase_ and
  AlwaysAccept short-circuit to accept. Verified.
- Reject-restore contract, mode-dependent (documented as such):
  * torsional: restore savedQ_, realizePosition, refill atom positions
    (:293-297);
  * docking: restore from the saved pre-kick CARTESIAN pose via
    setAtomsLocationsInGround (restoring q alone would reconstruct the KICKED
    pose, since the kick redefined the body frames), plus restore energy fields
    (:279-292);
  * Cartesian: restore savedPosG_ and peOld (:49-53).
  The restore is complete in each mode -- verified against the saved-state sets.

## Dispatch
- generateSample routes NcmcSwitch -> ncmcMove; else Cartesian (device MD) vs
  torsional/docking (constrained Verlet) branch (:28-54). Documented as the thin
  Move selector the ticket anticipated.

## Coverage
- TestEnsembleValidation / TestMassScaleInvariance (SLOW) and TestEquipartition /
  TestKineticEnergy / TestMassMatrix are the detailed-balance / INV-5 oracles;
  TestIntegrator is characterization. No gap beyond the SLOW-tier cost.

## Assumed notes
None.
