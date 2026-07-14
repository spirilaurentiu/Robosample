# DOC-RobotModel findings

Source: include/RobotModel.hpp. No code changed (comment-stripped diff empty).

## Immutability boundary — verified

Documented the const-after-build contract: ModelBuilder (World::buildModel) is
the sole writer; every solver takes `const RobotModel&`; rebuilt only on a
root-mobility change (World::setRootMobility/setRootMobilities, confirmed at
src/world/ModelBuilder.cpp:437-456). Two documented exceptions to "set once at
build", both recovered from the field-group comments and consistent with the
mass-op call sites: bodyMass/bodyCom_B/bodyUnitInertia_B/atomStation_B are
recomputed per coordinate transfer; bodyMassScale is a sampling-time
run-constant. These are the immutability boundary the ticket asked to state.

## Joint-fact predicates — meanings verified against call sites (not names)

- jointHasConstantHFM: verified against RobotEngine::realizeVelocity
  (src/RobotEngine_kinematics.cpp:255), which branches on it to pick the cheap
  centripetal vs the general HDot_FM*u bias path. Documented accordingly.
- jointIsLegalRoot: verified against the World::buildModel root validation
  (per-molecule root joint from rootMobilities[mol]).
- jointNU/jointNQ/jointUsesQuaternion: verified against the q/u layout and the
  quaternion renormalizer; nq>nu only for Ball/FreeLine/Free.
Agrees with the ticket hypothesis (predicates consistent with JointKernels
taxonomy).

## SoA data fields

The ~40 public table members retain their inline index-space/unit annotations
(the original authors' evidence) and are covered by a struct-level contract
stating BFS/topological order, the q/u offset tables, and the per-transfer vs
fixed split. Documented at group granularity per the ticket; no field was
reduced to a name restatement.

## Minor (not a finding, noted)

jointIsLegalRoot carries redundant self-referential `// == Cartesian` /
`// == Rigid` inline comments; harmless, left untouched (documentation-only pass).

## Assumed: notes
None.
