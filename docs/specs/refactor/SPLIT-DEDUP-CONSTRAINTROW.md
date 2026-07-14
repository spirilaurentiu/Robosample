# SPLIT-DEDUP-CONSTRAINTROW: one assembleConstraintRow for the three G-assembly loops

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)

`ARCHITECTURE.md section 6` and `MODULES.md section 2` record that the constraint-Jacobian (`G`)
row assembly is written three times in the `ConstraintSet` code with no shared
function. The three sites build, per loop-closure constraint, the same thing: the
generalized-force row `jacobianT[con]` from the Cartesian bond vector, then its
mass-weighted image `minvJacobianT[con] = M^-1 jacobianT[con]`. The sites are:

1. `ConstraintSet::enforceVelocityConstraints` (RATTLE) - `src/Constraints.cpp`
   lines 57-76 (inner `for (con)` loop; also computes `rhs[con]`).
2. `ConstraintSet::calcConstraintLogDet` (loop-closure Fixman log-det) -
   `src/Constraints.cpp` lines 153-163.
3. `ConstraintSet::enforcePositionConstraints` (SHAKE) - `include/Constraints.hpp`
   lines 130-143 (inner `for (con)` loop, inside the header template; also computes
   `violation[con]`).

The three assemblies MUST stay identical - that identity is INV-6. Today it is
maintained by copy discipline; this ticket makes it structural by routing all three
through one private helper `ConstraintSet::assembleConstraintRow`.

Invariant from `ARCHITECTURE.md section 5`, restated in full:

- **INV-6 Constraint consistency.** The `G` (constraint Jacobian) assembly is
  identical across SHAKE, RATTLE, and the loop-closure Fixman log-det, so the
  correction matches the projection. `calcConstraintLogDet` returns 0 for acyclic
  molecules. This ticket's single `assembleConstraintRow` is the mechanism that
  makes "identical" hold by construction rather than by discipline.

## Moves (exactly what)

New private static helper on `ConstraintSet` (declared in `include/Constraints.hpp`
near `solveCoupling`/`solveSmallSpd`; defined in `src/Constraints.cpp`):

```
static void assembleConstraintRow(const RobotModel& model, RobotState& state,
                                  const DistanceConstraint& bond,
                                  std::vector<Vec3>& atomForceScratch,
                                  Real* jacobianTRow, Real* minvJacobianTRow);
```

The common body it captures (identical at all three sites):
- `const Vec3 vecAB = posGround[bond.atomA] - posGround[bond.atomB];`
- zero the `atomForce` scratch, then `atomForce[bond.atomA] = vecAB;
  atomForce[bond.atomB] = Vec3(0) - vecAB;`
- `mapAtomForcesToGeneralizedForces(model, state, atomForce.data(), jacobianTRow);`
- `RobotEngine::multiplyByMInv(model, state, jacobianTRow, minvJacobianTRow);`

`vecAB` is also needed by the site-specific pieces (site 1's `rhs`, site 3's
`violation`), so either return `vecAB` from the helper or have the caller recompute
it. To keep the site-specific arithmetic byte-identical, the helper SHALL return
`vecAB` (the callers then reuse the returned value for `rhs`/`violation`, preserving
the exact floating-point expression). This makes the change behavior-preserving:
`vecAB`, the atomForce fill, the `mapAtomForcesToGeneralizedForces` call, and the
`multiplyByMInv` call are literally the same operations in the same order.

Rewrites at the three sites - each inner-loop assembly block is replaced by a call to
`assembleConstraintRow`, keeping the surrounding loop, `rhs`/`violation` computation,
and solve calls unchanged:
- `enforceVelocityConstraints`: the block at lines 69-75 becomes the helper call;
  the `rhs[con] = dot(vecAB, velA - velB)` at line 67 uses the returned `vecAB`
  (the `velA`/`velB`/station computation at 61-66 stays above the call).
- `calcConstraintLogDet`: the block at lines 156-162 becomes the helper call.
- `enforcePositionConstraints` (header template): the block at lines 136-142 becomes
  the helper call; `violation[con] = dot(vecAB, vecAB) - restLength^2` at line 133
  uses the returned `vecAB`.

## Public API after this ticket

`include/Constraints.hpp` gains one PRIVATE static member `assembleConstraintRow`.
The public interface (`mapAtomForcesToGeneralizedForces`, `enforceVelocityConstraints`,
`calcConstraintLogDet`, `enforcePositionConstraints`, `applyVelSpaceIncrementToQ`,
`numConstraints`, `empty`) is unchanged. `solveCoupling`/`solveSmallSpd` stay private.
No `robo_bindings...so` symbol delta (the helper is engine-internal; a private member
adds no Python export).

Because a section 3 defect-resolution ticket MAY add a thin private symbol and this adds NO
public/Python symbol, VERIFY section 2 criterion 4 is met without a separate API spec.

## Constraints

- Behavior-preserving. Zero arithmetic change: the moved four statements are
  identical, and returning `vecAB` preserves the exact expression the site-specific
  `rhs`/`violation` lines consume.
- No include-structure change: `assembleConstraintRow` lives in the existing
  `Constraints.{hpp,cpp}`; the header template `enforcePositionConstraints` already
  sees `RobotEngine::multiplyByMInv` and `mapAtomForcesToGeneralizedForces`.
- Invariant that SHALL remain true: **INV-6** (restated above). After the dedup the
  three assemblies are the SAME function - the identity is no longer a review burden.
- One commit for the dedup; formatting fixes separate.

## Predicted test breakage (Phase-A include fixes only)

**None.**
- `tests/TestConstraints.cpp` defines its OWN test-local free function
  `assembleConstraintRow` (lines 123-207, called at 146/204/401) as an independent
  oracle - analogous to `RobotLinearAlgebra.hpp`. It is a different symbol (file-local
  free function vs `ConstraintSet` private static) in a different scope; the new
  private member does not collide with it and is not visible to the test (private +
  no `ConstraintTestAccess` forwarder for it). The test is NOT edited.
- `tests/TestConstraintSolver.cpp` reaches `solveSmallSpd`/`solveCoupling` through
  `tests/ConstraintTestAccess.hpp` (the `friend struct ConstraintTestAccess`
  forwarders). Those private statics are unchanged, so the friend window is unchanged.
  No include fix.
- The RATTLE/SHAKE/Fixman behavioral tests (`TestConstraints.cpp`,
  `TestCyclicBoltzmann.cpp`, `TestFixman*.cpp`) exercise the public methods, whose
  outputs are bit-identical. No assertion-count change.

## Exit criteria (machine-checked)

- Full build passes; test suite / baseline outputs identical to Stage-0 (B1, B2, B3),
  both modes. `TestConstraintSolver.cpp`'s direct `solveSmallSpd` ln|det| checks and
  the cyclic-Boltzmann Fixman tests are the sharpest detectors that the three
  assemblies stayed identical.
- Public-symbol nm diff vs B0: empty on Python surface; one added engine-internal
  private static, recorded in the section 3 nm-delta rationale.
- include-cycle script clean (no new edge).
- clang-tidy / clang-format / IWYU clean on `Constraints.{hpp,cpp}`.
- `Constraints.hpp` <= 300 LOC; `Constraints.cpp` <= 600 LOC (both shrink).
- No test file changed.
