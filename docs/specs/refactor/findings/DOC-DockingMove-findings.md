# DOC-DockingMove findings

Style: house `//` contract blocks (see DOC-World-findings.md). No code changed.

## State
world/sampler/DockingMove.cpp (repositionLigands, sphere sampling,
findGoodStartingPose, configureDocking) and the World.hpp declarations
(263-291, 634-646) carried full contract documentation from SPLIT-W2. No
residual edit required.

## Verified hypotheses
- repositionLigands is a PURE PROPOSAL (relocate + reorient, no accept/reject);
  acceptance is judged in generateSample against Hold referenced to the pre-kick
  potential (World.hpp:634-640, HmcMove.cpp:56-132). Separation documented.
- Proposal distribution: uniform position in the per-ligand auto-sized sphere +
  uniform SO(3) reorientation -- a SYMMETRIC Cartesian proposal on the ligand's
  external q, hence Metropolised on dU alone (World.hpp:42-46 MoveType::RigidKick;
  DockingMove.cpp banner). The symmetry is why no explicit proposal-Jacobian term
  appears; documented as the correctness basis.
- findGoodStartingPose is INITIALIZATION (pre-sampling clash-free placement
  search), not part of the reversible chain; called by Context::runREX before
  round 0 (World.hpp:286-291). Documented as setup, distinct from the move.

## Coverage gap (ticket-flagged)
- No dedicated docking test found. The proposal/acceptance contract is recovered
  from generateSample's docking branch and the RNG usage, not from a passing
  oracle -- characterization-quality evidence for the symmetry claim, not a
  contract test. Recorded as a gap for the Architect.

## Assumed notes
None.
