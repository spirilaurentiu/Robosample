# DOC-FixmanCorrection findings

Style: house `//` contract blocks (see DOC-World-findings.md). No code changed.

## State
world/FixmanCorrection.cpp (calcFixman, calcLogSineSqrGamma2) and the World.hpp
private declarations (625-632) carried full contract documentation from
SPLIT-W5. No residual edit required.

## Verified hypotheses (sign convention and term set -- Critical if wrong)
- calcFixman returns exactly `0.5 * RT_ * (lnDetM - lnDetZ - lnDetMCartesian_)`
  = 1/2 RT (ln|M_tree| - ln det(G M^-1 G^T) - ln|M_3N|), Spiridon & Minh 2017
  Eq. 3 (FixmanCorrection.cpp:119-121). Sign and the three summed terms are
  stated unambiguously in-place.
- INV-6: `calcConstraintLogDet` (lnDetZ) returns 0 for acyclic molecules, so the
  correction reduces to the tree term alone on trees (:115-116). Verified.
- INV-5: the kinetic-metric log-det is `RobotEngine::calcLogDetM` (the O(n) ABA
  D-determinant, :105); calcFixman composes it with the constraint log-det.
- calcLogSineSqrGamma2 (external-rotation Jacobian) is default OFF and documented
  as wrong for the engine's unit-quaternion roots (flat S^3 == Haar SO(3), no
  Jacobian); retained behind useOrientationJacobian for future Euler joints
  (:124-149). Consistent with the Tier-0-green Fixman campaign (MEMORY).

## Assumed notes
None.
