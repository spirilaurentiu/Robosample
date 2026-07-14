#include "OpenMMContext.hpp"

#include "bridge/GpuKinematics.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

// Direct device-buffer access for the fused robot-kinematics pipeline. OpenMM is a
// vendored source build, so reaching CudaContext (via the public
// Context::getImpl().getPlatformData()) and driving kernels through the common
// ComputeContext abstraction is fully in-tree. All new code stays behind USE_CUDA.
#if USE_CUDA
#    include "../../openmm/platforms/cuda/include/CudaPlatform.h"
#    include "../../openmm/openmmapi/include/openmm/internal/ContextImpl.h"
#    include "../../openmm/platforms/cuda/include/CudaContext.h"
#    include "../../openmm/platforms/common/include/openmm/common/ComputeContext.h"
#    include "../../openmm/platforms/common/include/openmm/common/ComputeArray.h"
#    include "../../openmm/platforms/common/include/openmm/common/ComputeKernel.h"
#    include "../../openmm/platforms/common/include/openmm/common/ComputeProgram.h"
#    include "../../openmm/platforms/common/include/openmm/common/ContextSelector.h"
#    include <map>
#    include <memory>
#endif

// ===========================================================================
//  Fused CUDA robot-kinematics pipeline (spec docs/specs/gpu-cartesian-kinematics)
// ===========================================================================
#if USE_CUDA
namespace {

/**
 * @brief nvrtc source for the two fused-pipeline kernels, `pushPositions` (K1)
 *        and `reduceForces` (K2), compiled at runtime by
 *        `OpenMMContext::ensureKinematicsConstants`.
 *
 * @par Decomposition
 *      Both kernels are 1D grid-stride loops (`for (i = GLOBAL_ID; i < N; i +=
 *      GLOBAL_SIZE)`), so correctness is independent of grid/block shape.
 *      K1: one logical work item per atom (`N = numAtoms`); it writes that
 *      atom's Ground position into `posq[invOrder[a]].xyz`. K2: one logical work
 *      item per body (`N = numBodies`); it serially sums that body's wrench over
 *      the body's real atoms from the CSR range `[bodyAtomsBeg[b],
 *      bodyAtomsEnd[b])`.
 *
 * @par Launch configuration
 *      Launched through `OpenMM::ComputeKernel::execute(numAtoms)` (K1) and
 *      `execute(numBodies)` (K2); the block size is chosen by OpenMM's
 *      `ComputeContext`. The grid-stride form means the launch site's work-item
 *      count is the only launch parameter that matters and it matches the index
 *      math (no grid/block constraint, no divergence constraint).
 *
 * @par Shared memory
 *      None. Neither kernel uses static or dynamic shared memory, `__syncthreads`,
 *      or warp primitives; work items are fully independent.
 *
 * @par Stream and completion
 *      Enqueued on the `ComputeContext`'s current stream, made current by an
 *      `OpenMM::ContextSelector` in every host wrapper. K1's results are valid
 *      for the subsequent force evaluation on the same stream; K2's per-body
 *      output is valid after `reduceForcesToBodies` downloads `bodyForce`, a
 *      blocking device->host copy that synchronizes the stream.
 *
 * @par Write pattern and determinism
 *      K1 writes each atom's `posq.xyz` exactly once (`.w`, the packed charge, is
 *      left untouched); it deliberately writes every atom including massless
 *      ones (see the in-kernel note - a massless robot atom is still a real
 *      OpenMM particle). K2 accumulates into per-thread registers with no
 *      atomics and a fixed CSR traversal order, so each body's wrench is written
 *      once and is bit-for-bit reproducible run to run. The transforms are
 *      evaluated in `double`; rotations are row-major (matching
 *      `robo::Mat33::elems`).
 *
 * @note INV-1 parity: K2 produces per-body wrenches `(angular = moment about the
 *       body origin, linear = net force)` in Ground identical to the host
 *       `reduceAtomForcesToBodies`, and applies the same INV-2 virtual-site skip
 *       (`isVirtual[a]` -> parents already carry the site force). This on-device
 *       reduction is the whole reason the host reducer was deduplicated; the
 *       parity is pinned by tests/TestForceReducer.cpp.
 *
 * @note No runtime dtype tag: the buffer element types are fixed at compile
 *       time. `real`/`real4` follow OpenMM's mixed-precision build
 *       (`HAS_POSQ_CORRECTION` defined iff `getUseMixedPrecision()`); the force
 *       buffer is `mm_long` fixed-point, decoded with `scale = 1/2^32`.
 */
const char* const kKinematicsKernelSource = R"KERNSRC(
KERNEL void pushPositions(
        GLOBAL const double* RESTRICT station,   // 3*numAtoms, body-frame (nm)
        GLOBAL const int* RESTRICT atomBody,     // numAtoms
        GLOBAL const int* RESTRICT isVirtual,    // numAtoms (nonzero => skip)
        GLOBAL const double* RESTRICT xgb,       // 12*numBodies: 9 row-major R + 3 p
        GLOBAL const int* RESTRICT invOrder,     // numAtoms: atom -> device slot
        GLOBAL real4* RESTRICT posq,
#ifdef HAS_POSQ_CORRECTION
        GLOBAL real4* RESTRICT posqCorrection,
#endif
        int numAtoms) {
    // Write EVERY atom's position, exactly like Context::setPositions (the host path
    // pushes all atomPosG). Do NOT skip massless atoms here: some atoms are massless in
    // the robot model (skipped in the force reduction K2, matching getForcesFromOpenMM)
    // yet are real, force-bearing particles in OpenMM -- leaving their posq unwritten
    // keeps a stale value and corrupts the force field. True OpenMM virtual sites are
    // fixed up afterward by computeVirtualSites (called in pushBodyTransforms).
    for (int a = GLOBAL_ID; a < numAtoms; a += GLOBAL_SIZE) {
        int base = 12*atomBody[a];
        double sx = station[3*a], sy = station[3*a+1], sz = station[3*a+2];
        double x = xgb[base+9]  + xgb[base+0]*sx + xgb[base+1]*sy + xgb[base+2]*sz;
        double y = xgb[base+10] + xgb[base+3]*sx + xgb[base+4]*sy + xgb[base+5]*sz;
        double z = xgb[base+11] + xgb[base+6]*sx + xgb[base+7]*sy + xgb[base+8]*sz;
        int i = invOrder[a];
        posq[i].x = (real) x;                     // .w (charge) untouched
        posq[i].y = (real) y;
        posq[i].z = (real) z;
#ifdef HAS_POSQ_CORRECTION
        posqCorrection[i].x = (real) (x - (real) x);
        posqCorrection[i].y = (real) (y - (real) y);
        posqCorrection[i].z = (real) (z - (real) z);
#endif
    }
}

KERNEL void reduceForces(
        GLOBAL const mm_long* RESTRICT force,    // 3*paddedNumAtoms fixed point, component-major
        GLOBAL const double* RESTRICT station,
        GLOBAL const int* RESTRICT isVirtual,
        GLOBAL const double* RESTRICT xgb,
        GLOBAL const int* RESTRICT invOrder,
        GLOBAL const int* RESTRICT bodyAtomsBeg,
        GLOBAL const int* RESTRICT bodyAtomsEnd,
        GLOBAL const int* RESTRICT bodyAtoms,
        GLOBAL double* RESTRICT bodyForce,       // 6*numBodies: 3 angular + 3 linear
        int numBodies,
        int paddedNumAtoms) {
    const double scale = 1.0/(double) 0x100000000;
    for (int b = GLOBAL_ID; b < numBodies; b += GLOBAL_SIZE) {
        int base = 12*b;
        double r0 = xgb[base+0], r1 = xgb[base+1], r2 = xgb[base+2];
        double r3 = xgb[base+3], r4 = xgb[base+4], r5 = xgb[base+5];
        double r6 = xgb[base+6], r7 = xgb[base+7], r8 = xgb[base+8];
        double lx = 0, ly = 0, lz = 0, ax = 0, ay = 0, az = 0;
        for (int k = bodyAtomsBeg[b]; k < bodyAtomsEnd[b]; ++k) {
            int a = bodyAtoms[k];
            if (isVirtual[a]) continue;           // parents already carry the site force
            int i = invOrder[a];
            double fx = scale*(double) force[i];
            double fy = scale*(double) force[i+paddedNumAtoms];
            double fz = scale*(double) force[i+2*paddedNumAtoms];
            double sx = station[3*a], sy = station[3*a+1], sz = station[3*a+2];
            double rx = r0*sx + r1*sy + r2*sz;    // r = X_GB.R * station (about body origin)
            double ry = r3*sx + r4*sy + r5*sz;
            double rz = r6*sx + r7*sy + r8*sz;
            lx += fx; ly += fy; lz += fz;
            ax += ry*fz - rz*fy;                  // angular += r x f
            ay += rz*fx - rx*fz;
            az += rx*fy - ry*fx;
        }
        bodyForce[6*b+0] = ax; bodyForce[6*b+1] = ay; bodyForce[6*b+2] = az;
        bodyForce[6*b+3] = lx; bodyForce[6*b+4] = ly; bodyForce[6*b+5] = lz;
    }
}
)KERNSRC";

/**
 * @brief Per-world device state for the fused pipeline: the compiled kernels,
 *        the constant device arrays, and the per-step staging buffers.
 *
 * One instance exists at a time, held by the file-static `gGpuKin` and keyed to
 * one `(ComputeContext, worldToken, numAtoms, numBodies)` tuple. It borrows the
 * OpenMM `ComputeContext` (@ref cc, non-owning) and owns every `ComputeArray`
 * and `ComputeKernel` member; those are freed - with the CUDA context current -
 * when `gGpuKin` is reset by `OpenMMContext::releaseGpuKinematics` or replaced
 * for a different world. Constant arrays (@ref atomBody, @ref isVirtual, the
 * `bodyAtoms*` CSR) are uploaded once; @ref station is re-uploaded per round,
 * and @ref xgb / @ref bodyForce every step.
 */
struct GpuKinematics {
    OpenMM::ComputeContext* cc = nullptr;
    const void* worldToken = nullptr;
    int numAtoms = 0;
    int numBodies = 0;
    bool mixed = false;
    OpenMM::ComputeArray station, atomBody, isVirtual, invOrder;
    OpenMM::ComputeArray bodyAtomsBeg, bodyAtomsEnd, bodyAtoms;
    OpenMM::ComputeArray xgb, bodyForce; // per-step: uploaded / downloaded each call
    OpenMM::ComputeKernel pushKernel, reduceKernel;
    std::vector<int> invOrderHost;
    std::vector<int> lastOrder; // cached OpenMM slot->particle order to detect reorders
    std::vector<double> bodyForceHost;

    /**
     * @brief Refreshes @ref invOrder (atom -> device slot) whenever OpenMM's
     *        internal atom order has changed since the last call, re-uploading it
     *        only on an actual change.
     *
     * OpenMM spatially re-sorts atoms for the neighbor list (roughly every 250
     * steps under a cutoff), and a re-sort can occur inside another world's
     * on-device MD between this world's steps, after `getAtomsWereReordered()`
     * has been cleared. That flag is therefore unreliable; instead the live
     * `cc->getAtomIndex()` is compared against the cached @ref lastOrder and the
     * inverse map is rebuilt only when they differ.
     *
     * @pre  Call before K1 (`pushBodyTransforms`) and again before K2
     *       (`reduceForcesToBodies`): the force evaluation between them can itself
     *       trigger a re-sort.
     * @post @ref invOrder holds, for every real atom `a`, the device slot it
     *       occupies; padding slots and out-of-range particles are ignored (all
     *       slots are scanned because after a re-sort a real atom may sit at a
     *       slot index `>= numAtoms`).
     */
    void syncInvOrder() {
        const std::vector<int>& order = cc->getAtomIndex();
        if (order.size() == lastOrder.size()
            && std::equal(order.begin(), order.end(), lastOrder.begin())) {
            return; // order unchanged; invOrder still valid
        }
        // order maps device slot i -> original particle order[i]. Its length is
        // paddedNumAtoms (> numAtoms), and after a spatial reorder a REAL atom can sit in
        // a slot index >= numAtoms while padding entries (order[i] >= numAtoms) occupy low
        // slots. So we must scan ALL slots and invert only the real particles -- inverting
        // just the first numAtoms slots misses real atoms in high slots (leaving their posq
        // unwritten by K1 -> stale positions) and risks OOB on padding indices.
        invOrderHost.assign(static_cast<std::size_t>(numAtoms), -1);
        const int nSlots = static_cast<int>(order.size());
        for (int i = 0; i < nSlots; ++i) {
            const int a = order[static_cast<std::size_t>(i)];
            if (a >= 0 && a < numAtoms) {
                invOrderHost[static_cast<std::size_t>(a)] = i; // atom a lives in device slot i
            }
        }
        invOrder.upload(invOrderHost.data());
        lastOrder = order;
    }
};

// The device sub-object lives here (OpenMMContext is itself a singleton). Freed by
// OpenMMContext::releaseGpuKinematics() with the CUDA context current.
std::unique_ptr<GpuKinematics> gGpuKin;

OpenMM::ComputeContext* getCudaComputeContext(OpenMM::Context* ctx) {
    if (ctx == nullptr) {
        return nullptr;
    }
    auto* pd = reinterpret_cast<OpenMM::CudaPlatform::PlatformData*>(ctx->getImpl().getPlatformData());
    if (pd == nullptr || pd->contexts.empty()) {
        return nullptr;
    }
    return pd->contexts[0];
}

} // namespace
#endif // USE_CUDA

auto OpenMMContext::cudaKinematicsAvailable() const -> bool {
#if USE_CUDA
    return initialized && cudaKinematicsEnabled_;
#else
    return false;
#endif
}

void OpenMMContext::releaseGpuKinematics() {
#if USE_CUDA
    if (!gGpuKin) {
        return;
    }
    if (gGpuKin->cc != nullptr) {
        OpenMM::ContextSelector selector(*gGpuKin->cc);
        gGpuKin.reset();
    } else {
        gGpuKin.reset();
    }
#endif
}

void OpenMMContext::ensureKinematicsConstants(const void* worldToken,
                                              int numAtomsIn,
                                              int numBodiesIn,
                                              const double* station,
                                              const int* atomBody,
                                              const int* isVirtual,
                                              const int* bodyAtomsBeg,
                                              const int* bodyAtomsEnd,
                                              const int* bodyAtoms,
                                              bool stationsChanged) {
#if USE_CUDA
    ensureInitialized();
    OpenMM::ComputeContext* cc = getCudaComputeContext(context.get());
    if (cc == nullptr) {
        return; // not a CUDA platform: caller falls back to the host path
    }
    // atomStation_B is refit on every coordinate transfer (once per generateSample), so it
    // is constant only over one round's mdSteps -- NOT over the world's lifetime. atomBody,
    // the body-atom CSR and masses ARE constant. Build the constant arrays + compile the
    // kernels once per world (token-gated); re-upload the stations and refresh invOrder on
    // every call, since a stale station makes K1/K2 use the previous round's body frames
    // (rigid-fit error -> wrong forces -> blow-up).
    if (gGpuKin && gGpuKin->cc == cc && gGpuKin->worldToken == worldToken
        && gGpuKin->numAtoms == numAtomsIn && gGpuKin->numBodies == numBodiesIn) {
        if (stationsChanged) { // else the only per-step device traffic is X_GB up / bodyForce down
            OpenMM::ContextSelector selector(*cc);
            gGpuKin->station.upload(station);
        }
        return; // invOrder is refreshed in pushBodyTransforms/reduceForcesToBodies
    }
    OpenMM::ContextSelector selector(*cc);
    gGpuKin.reset(); // free any prior world's device arrays (context current)
    gGpuKin = std::make_unique<GpuKinematics>();
    GpuKinematics& g = *gGpuKin;
    g.cc = cc;
    g.worldToken = worldToken;
    g.numAtoms = numAtomsIn;
    g.numBodies = numBodiesIn;
    g.mixed = cc->getUseMixedPrecision();

    g.station.initialize<double>(*cc, static_cast<std::size_t>(3 * numAtomsIn), "robo_station");
    g.atomBody.initialize<int>(*cc, static_cast<std::size_t>(numAtomsIn), "robo_atomBody");
    g.isVirtual.initialize<int>(*cc, static_cast<std::size_t>(numAtomsIn), "robo_isVirtual");
    g.invOrder.initialize<int>(*cc, static_cast<std::size_t>(numAtomsIn), "robo_invOrder");
    g.bodyAtomsBeg.initialize<int>(*cc, static_cast<std::size_t>(numBodiesIn), "robo_bodyAtomsBeg");
    g.bodyAtomsEnd.initialize<int>(*cc, static_cast<std::size_t>(numBodiesIn), "robo_bodyAtomsEnd");
    g.bodyAtoms.initialize<int>(*cc, static_cast<std::size_t>(numAtomsIn), "robo_bodyAtoms");
    g.xgb.initialize<double>(*cc, static_cast<std::size_t>(12 * numBodiesIn), "robo_xgb");
    g.bodyForce.initialize<double>(*cc, static_cast<std::size_t>(6 * numBodiesIn), "robo_bodyForce");
    g.bodyForceHost.assign(static_cast<std::size_t>(6 * numBodiesIn), 0.0);

    g.station.upload(station);
    g.atomBody.upload(atomBody);
    g.isVirtual.upload(isVirtual);
    g.bodyAtomsBeg.upload(bodyAtomsBeg);
    g.bodyAtomsEnd.upload(bodyAtomsEnd);
    g.bodyAtoms.upload(bodyAtoms);
    g.syncInvOrder();

    std::map<std::string, std::string> defines;
    if (g.mixed) {
        defines["HAS_POSQ_CORRECTION"] = "1";
    }
    OpenMM::ComputeProgram program = cc->compileProgram(kKinematicsKernelSource, defines);
    g.pushKernel = program->createKernel("pushPositions");
    g.pushKernel->addArg(g.station);
    g.pushKernel->addArg(g.atomBody);
    g.pushKernel->addArg(g.isVirtual);
    g.pushKernel->addArg(g.xgb);
    g.pushKernel->addArg(g.invOrder);
    g.pushKernel->addArg(cc->getPosq());
    if (g.mixed) {
        g.pushKernel->addArg(cc->getPosqCorrection());
    }
    g.pushKernel->addArg(numAtomsIn);

    g.reduceKernel = program->createKernel("reduceForces");
    g.reduceKernel->addArg(cc->getLongForceBuffer());
    g.reduceKernel->addArg(g.station);
    g.reduceKernel->addArg(g.isVirtual);
    g.reduceKernel->addArg(g.xgb);
    g.reduceKernel->addArg(g.invOrder);
    g.reduceKernel->addArg(g.bodyAtomsBeg);
    g.reduceKernel->addArg(g.bodyAtomsEnd);
    g.reduceKernel->addArg(g.bodyAtoms);
    g.reduceKernel->addArg(g.bodyForce);
    g.reduceKernel->addArg(numBodiesIn);
    g.reduceKernel->addArg(cc->getPaddedNumAtoms());
#else
    (void) worldToken;
    (void) numAtomsIn;
    (void) numBodiesIn;
    (void) station;
    (void) atomBody;
    (void) isVirtual;
    (void) bodyAtomsBeg;
    (void) bodyAtomsEnd;
    (void) bodyAtoms;
#endif
}

void OpenMMContext::pushBodyTransforms(const double* xgbFlat) {
#if USE_CUDA
    if (!gGpuKin) {
        return;
    }
    GpuKinematics& g = *gGpuKin;
    OpenMM::ContextSelector selector(*g.cc);
    g.syncInvOrder(); // a reorder may have happened in another world's MD since our last step
    g.xgb.upload(xgbFlat);
    g.pushKernel->execute(g.numAtoms); // K1: posq <- X_GB*station

    // Writing posq directly bypasses the bookkeeping Context::setPositions does AFTER the
    // position write: zero the periodic-image cell offsets and re-sort atoms (reorderAtoms
    // self-limits to ~every 250 steps; between, OpenMM's displacement tracking keeps the
    // neighbor list valid, exactly as for the host setPositions path). This refreshes the
    // spatial decomposition after another world's on-device MD has moved the atoms.
    for (auto& off : g.cc->getPosCellOffsets()) {
        off = OpenMM::mm_int4(0, 0, 0, 0);
    }
    g.cc->reorderAtoms();
    if (hasVirtualSites) {
        context->computeVirtualSites();
    }
#else
    (void) xgbFlat;
#endif
}

auto OpenMMContext::computeForcesAndEnergyOnDevice() -> double {
#if USE_CUDA
    ensureInitialized();
    // Compute forces + energy on device from the posq we just wrote; leaves the
    // fixed-point force buffer populated for reduceForcesToBodies (no host download).
    potentialEnergy = context->getImpl().calcForcesAndEnergy(true, true);
    return potentialEnergy;
#else
    return 0.0;
#endif
}

void OpenMMContext::reduceForcesToBodies(double* bodyForceGFlat) {
#if USE_CUDA
    if (!gGpuKin) {
        return;
    }
    GpuKinematics& g = *gGpuKin;
    OpenMM::ContextSelector selector(*g.cc);
    g.syncInvOrder(); // the force calc between push and reduce can itself reorder atoms
    g.reduceKernel->execute(g.numBodies);
    g.bodyForce.download(g.bodyForceHost);
    std::copy(g.bodyForceHost.begin(), g.bodyForceHost.end(), bodyForceGFlat);
#else
    (void) bodyForceGFlat;
#endif
}
