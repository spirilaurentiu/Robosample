# DOC-ForceReducer: the single host force->wrench reduction

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).
**Calibration exemplar** (DOC-CALIBRATION section 1): fully documented during calibration;
this ticket covers only any residual symbols and records that.

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the single host implementation of the per-atom-force -> per-body
  spatial-wrench reduction, deduplicated from `ForceBridge.hpp` and the CUDA kernel
  by `SPLIT-DEDUP-FORCEREDUCER` (ARCHITECTURE section 6.5, MODULES.md section 2).
- **Layer:** Bridge (ARCHITECTURE section 2).
- **Ownership:** **borrows** OpenMM's per-atom force buffer (read-only), **writes**
  `RobotState::bodyForceG` (one wrench per body). Owns nothing.
- **Invariants (HYPOTHESES):**
  - **INV-1 force->wrench convention** (ARCHITECTURE section 5). `bodyForceG` is `(torque
    about the body origin, net force)` in the Ground frame. State this as the `@post`.
  - **INV-2 virtual sites** (ARCHITECTURE section 5). Massless particles are skipped because
    OpenMM already projected virtual-site forces onto parents; document the skip as
    behavior (avoids double counting), not as a loop guard.
  - Host<->device parity: this host reducer and the CUDA `reduceForces` kernel produce
    identical per-body wrenches on identical input - the dedup's whole purpose. The
    parity test (VERIFY section 3) enters the B1 baseline before the dedup lands; cite it.

## 2. Scope

- **Files:** `bridge/ForceReducer.{hpp,cpp}`.
- **Public symbols:** the reduction entry point(s).
- **Known gaps to close:** state the torque reference point (body origin) and the
  frame (Ground) explicitly - INV-1 is convention-only in source, so it SHALL read as
  a stated post-condition, not be left implicit.

## 3. Evidence pointers (tests exercising the module)

- TestReactionForces - FAST contract (INV-1; TESTS.md section 3). The host-vs-CUDA parity
  test added by the dedup ticket (VERIFY section 3) is the direct INV-1 oracle.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; INV-1/INV-2 stated as
  contract; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-ForceReducer-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

This module is a calibration exemplar. After human approval it is itself an attached
exemplar for the rest; keep its documented form authoritative. Spec wins any
conflict with an exemplar; conflict is a finding.
