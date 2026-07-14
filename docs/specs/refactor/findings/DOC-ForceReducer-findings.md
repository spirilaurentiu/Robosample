# DOC-ForceReducer findings

Calibration exemplar. Documented at the header declaration
(`include/bridge/ForceReducer.hpp`), not at the `.cpp` definition, per the
"document at the declaration, never both" rule; `src/bridge/ForceReducer.cpp`
keeps its M-site derivation as an implementation comment.

- INV-1 stated as `@post` (linear = net force; angular = moment about the body
  origin; Ground frame). INV-2 (massless-slot skip to avoid double-counting the
  already-projected virtual-site force) stated as `@note`, as behavior not a
  loop guard.
- Ownership: borrows the per-atom arrays read-only; accumulates into
  caller-owned, caller-zeroed `bodyForceG`. Stated as `@pre` (zero-init) and
  `@pre` (in-range `atomBody`).
- Parity oracle cited: `tests/TestForceReducer.cpp` drives identical OpenMM
  per-atom forces through this host reducer and the CUDA `reduceForces` kernel
  and asserts per-body wrenches agree (INV-1), and independently sums real atoms
  only (INV-2). Confirmed at call sites: `ForceBridge::getForcesFromOpenMM` and
  the parity test are the only callers.

No `@note Assumed:` entries. No ticket hypothesis contradicted.
