# DOC-GpuKinematics findings

CUDA calibration exemplar. The two device kernels are runtime nvrtc source
strings (`kKinematicsKernelSource`), not compiled `.cu` code; the launch
contract is documented on that string definition and the host entry-point
contracts on their `OpenMMContext.hpp` declarations.

## Launch contract (recovered; index math and launch sites agree)
- Decomposition: K1 `pushPositions` one work item per atom
  (`for i = GLOBAL_ID; i < numAtoms; i += GLOBAL_SIZE`), writing
  `posq[invOrder[a]].xyz`; K2 `reduceForces` one work item per body, serially
  summing the body's CSR atom range.
- Launch config: `ComputeKernel::execute(numAtoms)` (K1) and
  `execute(numBodies)` (K2); block size chosen by OpenMM. Grid-stride form means
  the work-item count is the only launch parameter, and it matches the index
  math. No disagreement between index math and launch site.
- Shared memory: none. No `__syncthreads`, no warp primitives; work items
  independent.
- Stream/completion: enqueued on the `ComputeContext` stream (made current via
  `ContextSelector`); K2 output valid after the blocking `bodyForce.download`.
- Determinism: K1 writes each `posq.xyz` once (`.w` preserved); K2 accumulates
  in registers with no atomics and fixed traversal order -> bit-reproducible.
- INV-1 parity with the host `reduceAtomForcesToBodies`, INV-2 virtual-site skip
  (`isVirtual[a]`), both pinned by `tests/TestForceReducer.cpp`.
- No runtime dtype tag: `real`/`real4` fixed by the mixed-precision build
  (`HAS_POSQ_CORRECTION`); force buffer is `mm_long` fixed-point (scale 1/2^32).

## Coverage gap
No dedicated CUDA unit test. Exercised indirectly via `TestAlchemy` (the only
real OpenMM Context) and the Level-1 fused run; the host-vs-CUDA reduction is
pinned by `TestForceReducer`. The push/kinematics path (K1, invOrder sync,
reorder handling) has no direct unit test. Recorded per ticket.

## Doxyfile
No Doxyfile edit made or required. The kernels are nvrtc string literals, not
CUDA syntax the C++ compiler or Doxygen parses; the translation unit uses no
`__global__`/`__device__` attributes at C++ level, so Doxygen has nothing to
choke on. The section-9 `PREDEFINED`/`EXTENSION_MAPPING` additions were
unnecessary here.

No `@note Assumed:` entries.
