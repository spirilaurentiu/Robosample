# DOC-RobotState findings

Executor: documenter. Source: include/RobotState.hpp. No code changed
(comment-stripped before/after diff empty).

## Stage-validity contract (INV-4) — verification result

The four-stage ladder (position -> velocity -> articulated-body inertias ->
udot) was recovered from the write sites, not the accessor names:

- Stage-P (realizePosition, src/RobotEngine_kinematics.cpp:70): writes X_GB,
  X_FM, X_PB, Phi, Mk_G, comG, H_FM, H, atomPosG, atomStationG.
- Stage-V (realizeVelocity, src/RobotEngine_kinematics.cpp:159): writes qdot
  (via calcQDot), V_FM, V_PB_G, V_GB, gyro, coriolisA, mobCoriolisA.
- Stage-A: factorizeArticulatedInertias (src/RobotEngine_dynamics.cpp:46)
  writes P, PPlus, G, DI (position-only); seedArticulatedCentrifugal
  (src/RobotEngine_dynamics.cpp:243) writes abCentrifugal and additionally
  requires Stage-V (reads a_mob, gyro). realizeArticulatedBodyInertias runs
  both.
- Stage-U (calcUDot, src/RobotEngine_dynamics.cpp:271): writes Z, zPlus, eps,
  A_GB, udot; reads the force inputs bodyForceG() and mobilityForce().
- qdotdot is written by calcQDotDot, which runs after Stage-U (reads udot).

No accessor's valid stage was contradicted across call sites. The integrator
(include/RobotIntegrator.hpp) and RobotEngine_reaction.cpp call the stages in
this order, and the FAST contract tests (TestKineticEnergy, TestMassMatrix,
TestRoboticsOracle*) read caches only after the corresponding stage. Contract
documented per accessor as @pre; agrees with the ticket hypothesis.

Input caches (not produced by a realization stage) documented as caller-owned
inputs instead of with an @pre stage: q, u (state), bodyForceG, mobilityForce
(force eval, written by ForceBridge), atomVelG, atomForceG (NCMC solvent).

## OQ-1 (convention vs type-level stage guard) — routed, not implemented

The stage-validity contract is enforced by call-order convention; the accessors
return raw mutable pointers with no compile-time stage tag. A type-level guard
(e.g. stage-tagged handles) would be a behavioral/API change and is out of doc
scope. Documented the convention as it exists; flagged here for the human OQ-1
decision. No presupposition of a target design.

## OQ-6 (NCMC-solvent field bleed) — documented as such, routed

RobotState carries solvent-relaxing-NCMC fields that are concern-bleed into the
per-step core cache. Each is marked "@note NCMC-solvent field (OQ-6)":
EnergySnapshot::keSolvent; atomVelG(); atomForceG(); setCartSolvent();
cartSolventAtoms(); cartSolventInvMass(); cartSolventMask(); wantsAtomForces();
and the private members cartSolventAtoms_/InvMass_/Mask_. All are guarded on
emptiness so the welded engine stays bit-identical. The relocate-vs-accept
decision is a human call; not asserted here.

## Dead accessor — finding (Medium, documentation/dead code)

BAT() (hBAT_, robo::Vec3[numZRows]) has no reader or writer anywhere in the
source tree (grep across src/include/tests/python: only the accessor definition
and the two allocation sites in allocateFull/allocateCompact). The slab is
reserved in both layouts but never accessed through the accessor. Documented
with an @warning; flagged here as a candidate for removal (out of doc scope;
Architect decision).

## Const-but-mutable accessors — documented, not a defect

Every pointer accessor is `const` yet returns a non-const `T*`; the solver
writes caches in place through these pointers. Documented at class level so the
apparent const-correctness quirk is not mistaken for a read-only guarantee.

## Test-coverage note

No dedicated RobotState unit test; exercised transitively by TestKineticEnergy,
TestMassMatrix, TestRoboticsOracle*, TestEnsembleOrientation (allocateFull).

## Assumed: notes

None. Every documented contract is backed by a write site or call site.
