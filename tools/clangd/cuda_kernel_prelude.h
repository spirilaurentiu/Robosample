/*
 * clangd force-include prelude for OpenMM CUDA kernel sources.
 *
 * PROVENANCE (where OpenMM's real definitions come from)
 * -----------------------------------------------------
 * A kernel's definitions are never a header. OpenMM assembles them at runtime,
 * in C++, inside CudaContext::createModule: it string-builds the real/mixed
 * typedefs from the chosen precision, emits the compilationDefines map (SHFL,
 * BALLOT, SYNC_WARPS, make_real*) that the CudaContext ctor set from the CUDA
 * driver version, prepends CudaKernelSources::common (the common.cu text), then
 * the calling *Kernels.cpp add per-force #defines (NUM_ATOMS from the System;
 * code snippets like LOAD_ATOM1_PARAMETERS synthesized from a Force's params).
 * NVRTC then compiles the concatenation. None of it is persisted to a file.
 *
 * SCOPE (what this prelude does and deliberately does NOT do)
 * ----------------------------------------------------------
 * clangd parses the kernels in real CUDA mode (see .clangd), so CUDA builtins
 * (__device__, threadIdx, float4, atomicAdd, warp intrinsics) resolve natively.
 * This prelude adds only OpenMM's universal layer, identical for every kernel
 * and near-frozen (vendored OpenMM 8.5):
 *   - cuda_types.h: vector_types.h + the real/mixed aliases + compilationDefines,
 *   - #includes of the two files that ARE files: common.cu + vectorOps.cu.
 *
 * It does NOT reproduce the per-force code-snippet macros (LOAD_ATOM1_PARAMETERS,
 * FIRST_EXCLUSION_TILE, BIN_SHIFT, ...). Those are OpenMM's runtime codegen;
 * faking them would mean reimplementing CudaNonbondedUtilities in a header. The
 * residual squiggles they leave in nonbonded.cu / findInteractingBlocks.cu mark
 * exactly where that generated code is injected. Simpler kernels come out clean.
 *
 * common.cu / vectorOps.cu edited directly get cuda_types.h instead (see
 * .clangd), so they are not double-parsed through the includes below.
 *
 * Analysis aid only. NOT what NVRTC compiles; never include from the build.
 */
#ifndef ROBOSAMPLE_CUDA_KERNEL_PRELUDE_H
#define ROBOSAMPLE_CUDA_KERNEL_PRELUDE_H

/* The self-contained scalar layer (vector_types.h, real/mixed aliases, macros). */
#include "cuda_types.h"  // IWYU pragma: keep

/* The two real files. Paths are relative to this file (tools/clangd/).
 *   common.cu    -> KERNEL/RESTRICT/GLOBAL/mm_long/mm_ulong/realToFixedPoint/...
 *   vectorOps.cu -> dot/cross/normalize + the vector operator overloads.
 * Included from the tree so this prelude never drifts from the real preamble. */
#include "../../openmm/platforms/cuda/src/kernels/common.cu"     // IWYU pragma: keep
#include "../../openmm/platforms/cuda/src/kernels/vectorOps.cu"  // IWYU pragma: keep

#endif /* ROBOSAMPLE_CUDA_KERNEL_PRELUDE_H */
