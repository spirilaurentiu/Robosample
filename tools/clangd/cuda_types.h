#ifndef ROBOSAMPLE_CUDA_TYPES_H
#define ROBOSAMPLE_CUDA_TYPES_H

#include <vector_types.h>

using real = float;
using real2 = float2;
using real3 = float3;
using real4 = float4;
using mixed = double;
using mixed2 = double2;
using mixed3 = double3;
using mixed4 = double4;
using tileflags = unsigned int;

#define SYNC_WARPS __syncwarp();
#define SHFL(var, srcLane) __shfl_sync(0xffffffff, var, srcLane);
#define BALLOT(var) __ballot_sync(0xffffffff, var);
#define USE_MIXED_PRECISION 1
#define make_real2 make_float2
#define make_real3 make_float3
#define make_real4 make_float4
#define make_mixed2 make_double2
#define make_mixed3 make_double3
#define make_mixed4 make_double4

#endif /* ROBOSAMPLE_CUDA_TYPES_H */
