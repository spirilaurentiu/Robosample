#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "OpenMM.h"
#include "TopologyElements.hpp"
#include "bridge/AlchemyForceFactory.hpp"

/**
 * @brief Process-singleton adapter that holds the one OpenMM
 *        System/Context/Integrator for the whole run and exposes the runtime
 *        energy/force evaluation surface plus configuration toggles.
 *
 * One OpenMM system serves the entire process. After construction the object
 * owns the OpenMM objects handed to it by `OpenMMSystemBuilder` and (under
 * `USE_CUDA`) the file-static `gGpuKin` device state. All evaluation entry
 * points require `initialize()` to have succeeded and throw
 * `std::runtime_error` otherwise.
 *
 * @note Units (INV-3): positions in nm, energies in kJ/mol, forces in
 *       kJ/mol/nm - the convention `ForceBridge`/`ForceReducer` consume.
 */
class OpenMMContext {
    public:
    /**
     * @brief Returns the process-wide singleton, constructing it on first call.
     * @return Reference to the single instance; lives for the whole process.
     */
    static auto get() -> OpenMMContext& {
        static OpenMMContext omm;
        return omm;
    }

    /**
     * @brief Builds the OpenMM System/Integrator/Context from @p systemTopology
     *        (via `OpenMMSystemBuilder`) and marks the singleton initialized.
     *
     * Must be called before any evaluation entry point. Reads the pre-set
     * toggles (`setMTS`, `setSeparateForceGroups`, alchemy enable) and, under
     * `USE_CUDA`, the `ROBO_CUDA_KINEMATICS` env var as the default for the fused
     * path.
     *
     * @param[in] systemTopology  SoA topology (particles, box, forces, method).
     *                            Borrowed for the duration of the call.
     * @return `true` on success; `false` iff OpenMM `Context` construction threw
     *         (the System/Integrator members are still populated, mirroring the
     *         builder's partial-write-on-failure).
     */
    auto initialize(const SystemTopology& systemTopology) -> bool;

    /// @brief Releases device kinematics state then the Context/System/Integrator
    ///        and resets to the uninitialized state. Safe to call when never
    ///        initialized.
    void shutdown() {
        releaseGpuKinematics(); // free device buffers BEFORE the CUDA context dies
        integrator.reset();
        context.reset();
        system.reset();
        forceGroupLabels.clear();
        hasVirtualSites = false;
        initialized = false;
    }

    /**
     * @brief Draws Maxwell-Boltzmann velocities on the live Context.
     * @param[in] temperature  Target temperature in Kelvin.
     * @param[in] seed         RNG seed for reproducibility.
     * @pre `initialize()` succeeded.
     */
    void setVelocitiesToTemperature(double temperature, int seed) {
        ensureInitialized();
        context->setVelocitiesToTemperature(temperature, seed);
    }

    // ---- NCMC alchemy (per-molecule intermolecular decoupling) -------------

    /**
     * @brief Records Region A (the decoupled atom set) on the owned
     *        `AlchemyForceFactory`; the decoupling forces are built later in
     *        `initialize()`, which holds the topology.
     *
     * Forwards to `AlchemyForceFactory::enableAlchemy`. Must be called before
     * `initialize()` to take effect. See @ref AlchemyForceFactory for the lambda
     * semantics and the endpoint contract.
     *
     * @param[in] atomIndices  Region A atom indices; stored ascending and
     *                         deduplicated (an arbitrary, possibly non-contiguous
     *                         set is allowed).
     */
    void enableAlchemy(const std::vector<int>& atomIndices);
    /// @brief Convenience overload building Region A as the contiguous range
    ///        `[atomBegin, atomEnd)`. @see enableAlchemy(const std::vector<int>&)
    void enableAlchemy(int atomBegin, int atomEnd);
    /**
     * @brief Sets the global `lambda_inter` on the live Context, driving the
     *        alchemical coupling.
     *
     * No-op when alchemy was never enabled. `lambda_inter == 1` reproduces the
     * unmodified force field (the endpoint HMC/NCMC acceptance relies on);
     * `lambda_inter == 0` fully removes the A<->rest coupling. See @ref
     * AlchemyForceFactory for what each path scales.
     *
     * @param[in] lambdaInter  Coupling parameter, conventionally in `[0, 1]`.
     * @pre `initialize()` succeeded (unless alchemy disabled, then no-op).
     */
    void setAlchemicalLambda(double lambdaInter);

    /// @brief Potential energy [kJ/mol] cached by the last evaluation/integration.
    /// @pre `initialize()` succeeded.
    [[nodiscard]] auto getPotentialEnergy() const -> double {
        ensureInitialized();
        return potentialEnergy;
    }
    /// @brief Kinetic energy [kJ/mol] cached by the last `integrateTrajectory`.
    /// @pre `initialize()` succeeded.
    [[nodiscard]] auto getKineticEnergy() const -> double {
        ensureInitialized();
        return kineticEnergy;
    }

    /**
     * @brief OpenMM's whole-system degree-of-freedom count, matching
     *        `StateDataReporter`'s convention.
     *
     * `3 * (particles with nonzero mass)` (massless virtual sites carry no
     * independent DOF), minus SHAKE/SETTLE distance constraints, minus 3 if a
     * `CMMotionRemover` is present. Consumed as the Cartesian world's `n_dof`,
     * because `RobotModel::nu` is not the physical DOF there (a Cartesian world
     * collapses every atom into a single `nu == 1` body).
     *
     * @return Physical degrees of freedom.
     * @pre `initialize()` succeeded.
     */
    [[nodiscard]] auto getNumDegreesOfFreedom() const -> int {
        ensureInitialized();
        int dof = 0;
        const int numParticles = system->getNumParticles();
        for (int i = 0; i < numParticles; ++i) {
            if (system->getParticleMass(i) > 0.0) {
                dof += 3;
            }
        }
        dof -= system->getNumConstraints();
        const int numForces = system->getNumForces();
        for (int i = 0; i < numForces; ++i) {
            if (dynamic_cast<const OpenMM::CMMotionRemover*>(&system->getForce(i)) != nullptr) {
                dof -= 3;
                break;
            }
        }
        return dof;
    }

    /// @brief One force group's labelled potential-energy contribution [kJ/mol],
    ///        returned by `computePotentialEnergyByGroup`.
    struct ForceGroupEnergy {
        int group = 0;
        std::string name;
        double energy = 0.0; // kJ/mol
    };

    /// @brief Toggle: when set, `initialize()` assigns each force its own group so
    ///        `computePotentialEnergyByGroup` can break energy down per force.
    ///        Set before `initialize()`.
    void setSeparateForceGroups(bool enabled) {
        separateForceGroups = enabled;
    }
    /// @brief Current separate-force-groups toggle.
    [[nodiscard]] auto getSeparateForceGroups() const -> bool {
        return separateForceGroups;
    }

    /**
     * @brief Multiple-timestep (r-RESPA) control. When enabled, `initialize()`
     *        builds an `MTSIntegrator` instead of a plain `VerletIntegrator`.
     *
     * Slow forces (Nonbonded, GBSA, NBFIX) go to force group 0 and evaluate once
     * per outer step; fast bonded forces go to group 1 and evaluate
     * @p innerSubsteps times. Affects only the Cartesian world's on-device MD
     * (`integrateTrajectory`); the torsional worlds always get the full force
     * sum. Must be called before `initialize()`.
     *
     * @param[in] enabled       Request MTS.
     * @param[in] innerSubsteps Fast-group substeps per outer step; a value `< 2`
     *                          disables MTS regardless of @p enabled.
     */
    void setMTS(bool enabled, int innerSubsteps) {
        useMTS = enabled && (innerSubsteps >= 2);
        mtsInnerSubsteps = innerSubsteps;
    }
    /// @brief Whether MTS is active (enabled and `innerSubsteps >= 2`).
    [[nodiscard]] auto getUseMTS() const -> bool {
        return useMTS;
    }

    /**
     * @brief Pushes positions to the live Context and places virtual sites.
     * @param[in] positions  Per-atom positions in OpenMM particle order (nm).
     *                       Borrowed.
     * @pre `initialize()` succeeded.
     */
    void setPositions(const std::vector<OpenMM::Vec3>& positions) {
        ensureInitialized();
        context->setPositions(positions);
        if (hasVirtualSites) {
            context->computeVirtualSites();
        }
    }
    /**
     * @brief Reads current positions from the live Context.
     * @param[out] out  Overwritten with per-atom positions in OpenMM particle
     *                  order (nm); wrapped into the primary box iff
     *                  `getEnforcePeriodicBox()`.
     * @pre `initialize()` succeeded.
     */
    void getPositions(std::vector<OpenMM::Vec3>& out) const {
        ensureInitialized();
        out = context->getState(OpenMM::State::Positions, enforcePeriodicBox).getPositions();
    }

    /**
     * @brief Toggle: whether OpenMM wraps coordinates into the primary box when
     *        state is pulled (positions/forces/energy `getState`). Default false.
     *
     * Must stay false under explicit solvent: the robot engine rebuilds each
     * molecule's internal frames from contiguous geometry, so a wrapped molecule
     * straddling a box face would yield a ~box-length "bond" and corrupt the
     * frame build. Energies and forces are unaffected (OpenMM always applies the
     * minimum image internally); this flag only changes the representation of
     * returned positions.
     */
    void setEnforcePeriodicBox(bool v) {
        enforcePeriodicBox = v;
    }
    /// @brief Current enforce-periodic-box toggle.
    [[nodiscard]] auto getEnforcePeriodicBox() const -> bool {
        return enforcePeriodicBox;
    }

    /// @brief True for nonbonded methods that require a periodic box
    ///        (CutoffPeriodic/Ewald/PME) and thus exclude GBSA.
    [[nodiscard]] static auto isPeriodic(NonbondedMethod m) -> bool {
        return m == NonbondedMethod::CutoffPeriodic || m == NonbondedMethod::Ewald
               || m == NonbondedMethod::PME;
    }

    // ---- CUDA robot kinematics pipeline (spec docs/specs/gpu-cartesian-kinematics) --

    /**
     * @brief Toggle: request the opt-in, CUDA-only fused path that computes
     *        `posG = X_GB*station + p` directly into OpenMM's device `posq` and
     *        reduces device forces to per-body wrenches on-device, so each robot
     *        step moves only `O(numBodies)` data up and down (no per-atom host
     *        round trip). Unavailable unless built with `USE_CUDA`.
     */
    void setCudaKinematics(bool enabled) {
        cudaKinematicsEnabled_ = enabled;
    }
    /// @brief Current fused-CUDA-kinematics toggle (not whether it can run).
    [[nodiscard]] auto getCudaKinematics() const -> bool {
        return cudaKinematicsEnabled_;
    }
    /// @brief True iff the fused CUDA path can run now: `USE_CUDA` build,
    ///        initialized, and toggle on. Always false on non-CUDA builds.
    [[nodiscard]] auto cudaKinematicsAvailable() const -> bool;
    /**
     * @brief Uploads the per-world constants (once per @p worldToken), rebuilds
     *        the atom->device-slot map, and compiles the two kernels on first use.
     *
     * All pointers are POD (no `robo::` types cross this boundary) and borrowed
     * for the call only. No-op / host-fallback on a non-CUDA platform.
     *
     * @param[in] worldToken     Opaque per-world identity; a change forces a full
     *                           rebuild. Borrowed as an identity value only.
     * @param[in] numAtoms       Atom count (sizes @p station, @p atomBody,
     *                           @p isVirtual, @p bodyAtoms).
     * @param[in] numBodies      Body count (sizes @p bodyAtomsBeg/End).
     * @param[in] station        Body-frame atom positions, `3*numAtoms` doubles
     *                           (xyz per atom, nm).
     * @param[in] atomBody       Body index per atom, `numAtoms`.
     * @param[in] isVirtual      Per-atom flag, `numAtoms`; nonzero => skipped in
     *                           the force reduction (INV-2).
     * @param[in] bodyAtomsBeg   CSR start offsets into @p bodyAtoms, `numBodies`.
     * @param[in] bodyAtomsEnd   CSR end offsets into @p bodyAtoms, `numBodies`.
     * @param[in] bodyAtoms      Body-sorted atom indices, `numAtoms`.
     * @param[in] stationsChanged  When the world is unchanged, re-upload @p station
     *                           only if true (stations are refit once per round).
     */
    void ensureKinematicsConstants(const void* worldToken,
                                   int numAtoms,
                                   int numBodies,
                                   const double* station,
                                   const int* atomBody,
                                   const int* isVirtual,
                                   const int* bodyAtomsBeg,
                                   const int* bodyAtomsEnd,
                                   const int* bodyAtoms,
                                   bool stationsChanged);
    /**
     * @brief K1: uploads @p xgbFlat, writes device `posq` (and `posqCorrection`
     *        under mixed precision), then places virtual sites on device.
     * @param[in] xgbFlat  Body transforms, `12*numBodies` doubles per body:
     *                     9 row-major rotation followed by 3 translation.
     *                     Borrowed. No-op unless the pipeline is set up.
     */
    void pushBodyTransforms(const double* xgbFlat);
    /**
     * @brief Computes forces and energy on device for the positions already in
     *        `posq`; leaves forces in the device fixed-point buffer for
     *        `reduceForcesToBodies`. No host download.
     * @return Potential energy [kJ/mol].
     * @pre `pushBodyTransforms` populated `posq` this step.
     */
    auto computeForcesAndEnergyOnDevice() -> double;
    /**
     * @brief K2: reduces the device force buffer to per-body wrenches (INV-1,
     *        identical to the host `reduceAtomForcesToBodies`) and downloads them.
     * @param[out] bodyForceGFlat  Per-body wrench, `6*numBodies` doubles per body:
     *                             3 angular (moment about the body origin) then
     *                             3 linear (net force), in Ground. Caller-sized.
     * @pre Runs after `computeForcesAndEnergyOnDevice`.
     */
    void reduceForcesToBodies(double* bodyForceGFlat);

    /**
     * @brief Sets @p positions, places virtual sites, and returns the total
     *        potential energy for that configuration.
     * @param[in] positions  Per-atom positions, OpenMM particle order (nm).
     *                       Borrowed.
     * @return Potential energy [kJ/mol]; also cached for `getPotentialEnergy()`.
     * @pre `initialize()` succeeded, else throws `std::runtime_error`.
     */
    auto computePotentialEnergy(const std::vector<OpenMM::Vec3>& positions) -> double;
    /**
     * @brief As `computePotentialEnergy`, additionally returning a per-force-group
     *        energy breakdown.
     * @param[in] positions  Per-atom positions, OpenMM particle order (nm).
     * @return `{total PE [kJ/mol], breakdown}`. The breakdown is per labelled
     *         group when `separateForceGroups` or MTS is active and labels exist;
     *         otherwise a single `{0, "All", total}` entry.
     * @pre `initialize()` succeeded, else throws `std::runtime_error`.
     */
    auto computePotentialEnergyByGroup(const std::vector<OpenMM::Vec3>& positions)
        -> std::pair<double, std::vector<ForceGroupEnergy>>;
    /**
     * @brief Sets @p positions, places virtual sites, and returns per-atom forces.
     * @param[in]  positions  Per-atom positions, OpenMM particle order (nm).
     * @param[out] outForces  Overwritten with per-atom forces (kJ/mol/nm), same
     *                        order. These are the input to `ForceReducer`; forces
     *                        left in massless virtual-site slots are handled by
     *                        the INV-2 skip downstream.
     * @pre `initialize()` succeeded, else throws `std::runtime_error`.
     */
    void evaluateForcesFromPositionsCache(const std::vector<OpenMM::Vec3>& positions,
                                          std::vector<OpenMM::Vec3>& outForces) const;
    /**
     * @brief Advances the live Context @p steps with the built integrator (Verlet
     *        or MTS), then caches the resulting potential and kinetic energy.
     * @param[in] steps                  Integration steps.
     * @param[in] timeStepInPicoseconds  Step size [ps]; set on the integrator.
     * @return `true` on success; `false` iff `integrator->step` threw (energies
     *         are still refreshed from the post-step state).
     * @pre `initialize()` succeeded, else throws `std::runtime_error`.
     */
    auto integrateTrajectory(int steps, double timeStepInPicoseconds) -> bool;

    private:
    OpenMMContext() = default;

    void ensureInitialized() const {
        if (!initialized) {
            throw std::runtime_error("OpenMM subsystem not initialized. Call OpenMMContext::initialize() "
                                     "before using any other functions.");
        }
    }

    std::unique_ptr<OpenMM::Context> context;
    std::unique_ptr<OpenMM::System> system;
    std::unique_ptr<OpenMM::Integrator> integrator;

    std::size_t numAtoms = 0;
    double potentialEnergy = 0;
    double kineticEnergy = 0;

    bool separateForceGroups = false;
    std::vector<std::pair<int, std::string>> forceGroupLabels;

    // NCMC per-molecule intermolecular decoupling: builders + control +
    // Region-A state, see bridge/AlchemyForceFactory.hpp (SPLIT-O4).
    AlchemyForceFactory alchemyFactory_;

    // r-RESPA multiple-timestep state.
    bool useMTS = false;
    int mtsInnerSubsteps = 4;

    bool enforcePeriodicBox = false;
    bool hasVirtualSites = false;
    bool initialized = false;

    // CUDA robot kinematics pipeline. The device sub-object (kernels + arrays) is a
    // file-static in OpenMMContext.cpp keyed to the singleton; released here.
    bool cudaKinematicsEnabled_ = false;
    void releaseGpuKinematics(); // no-op unless USE_CUDA
};