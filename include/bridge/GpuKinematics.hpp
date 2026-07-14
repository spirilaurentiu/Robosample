#pragma once

/**
 * @file GpuKinematics.hpp
 * @brief Translation-unit marker for the fused CUDA robot-kinematics pipeline
 *        (spec docs/specs/gpu-cartesian-kinematics); CUDA-only (`#if USE_CUDA`).
 *
 * This header declares no symbols. The pipeline's device sub-object - the K1
 * `pushPositions` and K2 `reduceForces` nvrtc kernels plus their host driver -
 * is anonymous-namespace, file-static state private to
 * `bridge/GpuKinematics.cpp` and never exposed. The host-side entry points that
 * drive it (`OpenMMContext::ensureKinematicsConstants`, `pushBodyTransforms`,
 * `computeForcesAndEnergyOnDevice`, `reduceForcesToBodies`,
 * `cudaKinematicsAvailable`, `releaseGpuKinematics`) are `OpenMMContext` member
 * definitions that live in `GpuKinematics.cpp`; their declarations - and their
 * contracts - stay in `OpenMMContext.hpp`, the only public surface of this unit.
 *
 * @see OpenMMContext.hpp for the pipeline's declared entry points and contracts.
 * @see bridge/GpuKinematics.cpp for the two device kernels and their launch
 *      contract (INV-1 on-device parity with the host `ForceReducer`).
 */
