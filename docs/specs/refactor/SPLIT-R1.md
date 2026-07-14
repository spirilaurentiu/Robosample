# SPLIT-R1: extract the dense n<=6 hinge linear algebra from RobotEngine.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)

`src/RobotEngine.cpp` (1456 LOC) mixes four responsibilities: the articulated-body
recursions (kinematics, inertia factorization, forward dynamics, mass operators,
reaction forces), a self-contained dense `n<=6` symmetric linear-algebra kernel set,
a `ROBO_DEBUG` NaN scanner, and file-local spatial helpers. This ticket extracts the
linear-algebra kernel set - the leaf with no dependency on `RobotModel`/`RobotState`,
only on `robo::Real`. It is the first ticket and lands before every other RobotEngine
split.

The extracted symbols implement the SINGLE source of truth for the null-space lock
(singular-DOF Fixman campaign, oracles CC1/CC4/N1; its spec is not yet in the repo,
so the behavior below is recovered from source): `invertDense` (dynamics
pseudo-inverse), `pseudoLogDet` (Fixman pseudo-determinant), and `symSqrtInv`
(sqrt(M) for NMA route B) all lock a null hinge direction on IDENTICAL `(d, V, tol)`
via `eigDecompAndTol`. That agreement is a correctness property; keeping these
routines in one unit preserves it.

Invariants below come from `ARCHITECTURE.md section 5` and are restated in full.

- **INV-5 Mass metric.** HMC momentum seeding and kinetic energy use the same mass
  matrix operators (`multiplyBySqrtM` / `multiplyByMInv` / `calcKineticEnergy`);
  `calcLogDetM` is the Fixman kinetic term. `symSqrt`/`symSqrtInv`/`pseudoLogDet`
  are the per-body kernels those operators call; their null-space convention SHALL
  stay identical to `invertDense`'s so the seed, the KE, and the Fixman log-det
  agree on which directions are locked.

## Moves (exactly what)

All symbols currently live in the anonymous namespace of `src/RobotEngine.cpp`
(the block spanning lines 68-521). They move to `include/math/hinge_linalg.hpp`
(declarations) + `src/math/hinge_linalg.cpp` (definitions), all inside
`namespace robo::detail`. The banner comments attached to each symbol move with it.

- `kNullLockAbs` (line 92) -> `hinge_linalg.hpp` [`robo::detail`], `constexpr`
- `nullLockTol` (lines 94-96) -> `hinge_linalg.{hpp,cpp}` [`robo::detail`]
- `eigDecompAndTol` (lines 104-117, incl. the shared-null-space-lock banner 75-91)
  -> `hinge_linalg.{hpp,cpp}` [`robo::detail`]
- `pseudoLogDet` (lines 129-145) -> `hinge_linalg.{hpp,cpp}` [`robo::detail`]
- `invertDense` (lines 348-382) -> `hinge_linalg.{hpp,cpp}` [`robo::detail`]
- `jacobiSymEig` (forward declaration line 73; definition lines 389-448)
  -> `hinge_linalg.{hpp,cpp}` [`robo::detail`]; the forward declaration is deleted
  (the header supplies it)
- `symSqrt` (lines 453-470) -> `hinge_linalg.{hpp,cpp}` [`robo::detail`]
- `symSqrtInv` (lines 493-518) -> `hinge_linalg.{hpp,cpp}` [`robo::detail`]

NOT moved (stay in the RobotEngine.cpp anonymous namespace):
- `spatialDot` (lines 320-322) - a `SpatialVec` helper, not a dense-block kernel;
  explicitly outside this ticket's symbol list. It stays until SPLIT-R3.
- the `robodbg` namespace (lines 148-317) - SPLIT-R2.

Internal-linkage change: these are `static`/anonymous-namespace today (internal
linkage). In `robo::detail` they gain external linkage. This is an internal-symbol
change only (VERIFY section 2 criterion 4 permits internal churn); no public symbol is
added or removed.

Call-site preservation: `RobotEngine.cpp` calls these unqualified (e.g.
`invertDense(D, dof, DI)` at line 836, `eigDecompAndTol` at line 872, `symSqrt` at
line 1154, `symSqrtInv` at line 1196, `pseudoLogDet` at line 1240). To keep every
call site byte-identical (zero renames), `RobotEngine.cpp` adds, after the include,
a file-scope `using robo::detail::invertDense; using robo::detail::eigDecompAndTol;
using robo::detail::pseudoLogDet; using robo::detail::symSqrt; using
robo::detail::symSqrtInv; using robo::detail::jacobiSymEig; using
robo::detail::nullLockTol;` (and `kNullLockAbs`). No call expression changes.

## Public API after this ticket

`math/hinge_linalg.hpp` exports in `robo::detail`: `kNullLockAbs`, `nullLockTol`,
`eigDecompAndTol`, `pseudoLogDet`, `invertDense`, `jacobiSymEig`, `symSqrt`,
`symSqrtInv`. Nothing else. Header self-contained (IWYU): `#pragma once`,
`#include <algorithm>`, `#include <cmath>`, `#include "robot_math.hpp"` (for
`robo::Real`). No `RobotModel`/`RobotState`/`ForceBridge` include. No forward
declarations needed.

The built `robo_bindings...so` public-symbol set is unchanged (these were internal).

## Constraints

- Pure code motion. Zero logic changes. Zero symbol renames. The `robo::detail`
  qualification is a namespace placement, not a rename; call sites are preserved by
  the using-declarations above.
- `src/RobotEngine.cpp` gains `#include "math/hinge_linalg.hpp"` and the
  using-declarations, and loses the eight moved symbols (and the forward
  declaration at line 73); nothing else in it changes. `spatialDot` and `robodbg`
  remain.
- Build system: add `src/math/hinge_linalg.cpp` to the engine object library
  target (the same target that compiles `RobotEngine.cpp`). `include/` is already
  on the compiler include path (headers are included bare), so `math/hinge_linalg.hpp`
  resolves without a new include directory.
- Invariant that SHALL remain true: **INV-5** (restated above). The null-space lock
  single-source-of-truth property - `invertDense`, `pseudoLogDet`, `symSqrtInv` all
  lock on the identical `(d, V, tol)` from `eigDecompAndTol` - SHALL survive the move
  unchanged (the three routines and `eigDecompAndTol` move together).
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)

**None.** No test includes the engine's anonymous-namespace kernels - they are
internal-linkage and unreachable. Every test that names `invertDense`,
`jacobiSymEig`, `symSqrt`, `symSqrtInv`, or `pseudoLogDet`
(`TestLinearAlgebra.cpp`, `TestLinearAlgebraOracle.cpp`, `TestMassMatrix.cpp`,
`TestRoboticsOracle.cpp`, `TestRoboticsOracleMolecule.cpp`, `TestNMALinearAlgebra.cpp`,
`TestFixmanIdealizedChains.cpp`) reaches them through the test-local independent
oracle `tests/RobotLinearAlgebra.hpp` (`namespace robo_linalg`), never through the
engine. That oracle is not touched by this ticket, so no include fix is needed and
no assertion count changes.

NOTE (conflict with the task directive; do NOT act on it inside R1). The
architecture note "RobotLinearAlgebra.hpp re-implements the engine's dense solver ->
after R1 point it at hinge_linalg" conflicts with that header's own charter
(`tests/RobotLinearAlgebra.hpp` lines 6-14): it is a deliberately independent,
hand-maintained oracle, and `TestLinearAlgebraOracle.cpp` diffs the engine against
it. Repointing it at `hinge_linalg` would make those oracle comparisons tautological
and delete the very independence the tests rely on. Repointing is therefore a
separate, human-approved Phase-B `TEST-###` decision, not part of this pure-motion
ticket. R1 leaves `RobotLinearAlgebra.hpp` unchanged.

## Exit criteria (machine-checked)

- Full build passes; test suite / baseline outputs identical to Stage-0 (VERIFY B1,
  B2, B3), both gated-off and gated-on modes.
- `nm -C --defined-only` diff on public symbols vs B0: empty. Internal-symbol churn
  (`invertDense` etc. now `robo::detail`, external linkage) is expected and allowed.
- include-cycle script clean (diff vs B4); `hinge_linalg.hpp` introduces no upward
  layer edge (it depends only on `robot_math.hpp`, the Infrastructure layer).
- clang-tidy / clang-format / IWYU clean on `RobotEngine.cpp`, `hinge_linalg.hpp`,
  `hinge_linalg.cpp`.
- Comment-stripped, whitespace-normalized before/after diff of `RobotEngine.cpp` is
  empty except the removed symbol blocks, the added `#include`, and the
  using-declarations.
- `hinge_linalg.hpp` <= 300 LOC; `hinge_linalg.cpp` <= 600 LOC (expected ~450).
- No test file changed.
