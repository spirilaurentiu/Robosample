# DOC-ReactionReporter findings

Style: house `//` contract blocks (see DOC-World-findings.md). No code changed.

## State
world/ReactionReporter.cpp, include/robo/world/ReactionReporter.hpp
(ReactionSample), and the World.hpp declarations (448-533) carried full contract
documentation from SPLIT-W3. No residual edit required.

## Verified hypotheses
- INV-1 convention: the reported per-body row is a spatial force about the body
  origin Bo, in Ground -- the same (force, torque-about-Bo, Ground) convention as
  bodyForceG and the ForceReducer calibration exemplar
  (ReactionReporter.hpp:6-11). Verified; phrasing reused from the exemplar.
- Opt-in: a world that never calls enable/setReactionReporter allocates nothing
  and does zero extra work (banner + World.hpp:502-533). Enable condition is part
  of the contract.
- Capture timing: `captureReactionSnapshot()` is intended AFTER generateSample()
  resolves accept/reject, on the accepted q, at DCD cadence (World.hpp:509-530).
  Dynamical per-frame, not post-hoc -- verified against the reaction-force
  monitoring campaign (MEMORY).
- No-perturbation property: the optional u=0 reaction term saves/zeros/restores u
  and every u-derived cache, so a snapshot cannot perturb the next round
  (World.hpp:519-525). Cartesian worlds throw (integrator guard). Verified.

## Assumed notes
None.
