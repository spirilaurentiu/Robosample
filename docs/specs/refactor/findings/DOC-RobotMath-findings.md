# DOC-RobotMath findings

## Scope decision (spec-vs-exemplar note)
- `documenter.md` 8 requires a `@param` with direction on every parameter. For
  the ~100 elementary value-type operators here (`Vec3::operator+`, `[]`, `+=`,
  scalar `*`, etc.) per-operator param docs would restate componentwise
  arithmetic already implied by the type. The ticket itself scopes this header
  "Breadth, not depth (ARCHITECTURE 7)". Resolution taken: documented each TYPE's
  contract (meaning, storage/index convention, ownership = trivially-copyable
  value) and gave full `@param`/`@return` to the semantically-loaded operations
  (`quaternionDotFromAngVel`, `crossMat`, `calcDihedralAngle`, the reexpress and
  spatial-inertia products). Elementary arithmetic operators inherit the type
  contract. Recorded here as the deliberate reading of the breadth-not-depth
  instruction.

## OQ-4: Inertia vs SpatialInertia vs MassProperties overlap (routed, not asserted)
- The three types were NOT documented as interchangeable. Their distinct
  call-site roles were documented as-is:
  - `Inertia`  - full second mass moment about a point; built from a point mass
    (`m(|p|^2 I - p p^T)`) or isotropic diagonal; additive. Used where a raw
    inertia tensor is accumulated.
  - `SpatialInertia` - the 6x6 spatial mass operator about the body origin
    (mass, com, UnitInertia); `operator*` yields spatial momentum. Delegates its
    block math to `ArticulatedInertia` (single source of truth).
  - `MassProperties` - storage of mass + com + per-unit-mass UnitInertia, with
    `reexpress` and `toSpatialInertia`.
  A whole-repo usage scan to decide whether any pair should be merged is a human
  decision (OQ-4), left OPEN. No doc claims they are equivalent.

## Verified conventions (documented as contract)
- `Vec3::operator%` is cross product, not modulo (verified by `crossMat` and the
  spatial-algebra tests). `SpatialVec` index `[0]`=angular, `[1]`=linear.
- `SymMat33` storage order `(xx,xy,yy,xz,yz,zz)`; `reexpress` = `~R I R`.
- INV-9: `quaternionDotFromAngVel` is the parent-frame map paired with the
  standard `R_FM` of `Rotation::fromQuaternion`; both `Quat::angVelToQdot` and
  `Rotation::convertAngVelToQuaternionDot` delegate to it. Documented the
  parent-frame precondition on the angular velocity as a `@warning` (INV-9).

## Notes
- No `@note Assumed:` used.
