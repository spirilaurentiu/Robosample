# DOC-GpuKinematics: fused nvrtc kinematics/force pipeline (CUDA)

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).
**Calibration exemplar** (DOC-CALIBRATION section 2): the CUDA calibration file; documented
during calibration. This ticket covers residual symbols and records that.

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the `gGpuKin` fused GPU pipeline - two nvrtc-JIT kernels that push
  per-body transforms into OpenMM's `posq` and reduce per-atom forces to per-body
  wrenches on-device, bypassing the host round trip (ARCHITECTURE section 3, section 7; MODULES.md
  O6). CUDA-only unit.
- **Layer:** Bridge (ARCHITECTURE section 2). `gGpuKin` is file-static CUDA state
  (ARCHITECTURE section 4).
- **Ownership:** owns device buffers, nvrtc modules, and the CUDA stream for the run;
  borrows OpenMM's `posq` (writes it) and force buffer (reads it). Recover buffer
  lifetimes from `gGpuKin` setup/teardown.
- **Invariants (HYPOTHESES):**
  - **INV-1** - the on-device `reduceForces` kernel produces per-body wrenches
    **identical** to the host `ForceReducer` (ARCHITECTURE section 5, section 6.5). This is the
    parity contract; cite the host-vs-CUDA parity test (VERIFY section 3).
  - CUDA launch contract per `documenter.md` section 5: recover, from the kernel index math
    **and** the host launch sites together - decomposition (thread<->body/atom
    mapping), launch config (grid/block, constraints), dynamic/static shared memory
    formula, `__syncthreads`/warp participation, the stream used and when results are
    valid (event/stream sync), and the write pattern/determinism of the reduction
    (atomics -> run-to-run bit variation, if any). Disagreement between index math and
    launch site is a finding, never a guess.
  - Runtime-selected element types, if the kernels are parameterized by a dtype tag
    (`documenter.md` section 5.4): document the tag<->buffer-match `@pre` at the dispatch
    boundary, the supported type set from evidence, and unsupported-tag behavior.

## 2. Scope

- **Files:** `bridge/GpuKinematics.{hpp,cpp}` (CUDA-only; nvrtc source strings +
  host launch wrappers).
- **Public symbols:** the pipeline setup/step/teardown surface and the host-side
  kernel launch wrappers.
- **Known gaps to close:** the kernels are nvrtc strings, not `.cu` - recover the
  launch contract from the host launch sites and the embedded index math. The
  Doxyfile CUDA settings (`documenter.md` section 9) may be needed; report in findings.

## 3. Evidence pointers (tests exercising the module)

- No dedicated CUDA unit test; exercised indirectly via `TestAlchemy` (the only real
  OpenMM `Context`, TESTS.md section 1) and the fused path in the Level-1 run. The
  host-vs-CUDA parity test (VERIFY section 3) is the INV-1 oracle. Record the direct-unit
  gap in findings.

## 4. Exit criteria

- Doxygen warning-free (CUDA settings applied); every public symbol and each
  `__global__` kernel documented with the full section 5 launch contract; INV-1 parity
  stated; comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-GpuKinematics-findings.md` present; any index-math/launch-site
  disagreement and the Doxyfile edit recorded; no unmatched `@note Assumed:`.

## 5. Calibration

This module is the CUDA calibration exemplar. After approval its documented form is
authoritative for later CUDA docs. Spec wins any conflict; conflict is a finding.
