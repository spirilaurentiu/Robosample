# DOC-RobotIntegrator findings

## Symbol location
- The three public entry points (verletStep, stepTo, checkReversibility) are
  declared in `include/RobotEngine.hpp` and defined as templates in
  `include/RobotIntegrator.hpp`. Their contract Doxygen was written at the
  DECLARATION in RobotEngine.hpp (documented once, never in both). The private
  I1 helpers (driftPositions, cartesianSolventDrift, cartesianSolventKick,
  velocityCorrector) are declared+defined in RobotIntegrator.hpp and documented
  there as internal helpers (local behavior, per documenter.md 4.1). The ticket
  named the private set "driftPositions/velocityCorrector/cartesianSolventVerlet";
  the Cartesian-solvent helper is actually split into a drift and a kick
  (cartesianSolventDrift + cartesianSolventKick). Documented both.

## Invariants (documented as contract)
- INV-4: verletStep runs the realization stages in the canonical order; stated as
  a `@post` and cross-referenced to the RobotEngine stage methods. The exact
  helper-call sequence is documented as contract-level behavior (reversible
  leapfrog: drift + SHAKE projection + trapezoid correction), not as narration; a
  behavior-preserving reorder of the helper calls does not falsify the text.
- INV-6: SHAKE position projection and the velocity correction use the same `G`
  as the Fixman log-det (Constraints); stated on verletStep.
- INV-9: the quaternion position drift uses the exact exponential-map advance
  (unit by construction, reversible); stated on driftPositions and verletStep.

## Deliberate divergence from Simbody (documented as `@note`, not "matched away")
- The quaternion drift uses the exponential map rather than Simbody's
  linear-Taylor-plus-renormalize update. This is intentional: it bypasses the
  quaternion second derivative (qddot) and so is robust to a latent N/qddot
  inconsistency in the port's Free-joint kinematics that otherwise pumps kinetic
  energy. It reduces to the linear update as h -> 0, changing only the proposal,
  not the target distribution. Documented as a `@note` on driftPositions.
- OPEN gap (root cause of the above): the Free-joint quaternion kinematics
  (N, Ndot, qddot) are not golden-tested against a single free-body reference.
  Until that is closed the exp-map is the correct propagator. Recorded as a
  coverage gap (overlaps JointKernels OQ-2 for the q,u-dependent joints).

## Test-evidence labeling
- TestIntegrator (drift/reversibility chains) and TestFreeJointKEPump are
  CHARACTERIZATION per TESTS.md 3; they pin current behavior and were not
  promoted to contract. checkReversibility is documented as a configuration-local
  smoke test (its own `@note`), consistent with the characterization label.

## Notes
- No `@note Assumed:` used.
