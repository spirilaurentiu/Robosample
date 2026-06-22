#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "OpenMM.h"
#include "TopologyElements.hpp"

// Thin Robosample-owned wrapper around a single OpenMM System/Context (a
// process-wide singleton). One OpenMM system for the whole run.
class OpenMMContext {
    public:
    static auto get() -> OpenMMContext& {
        static OpenMMContext omm;
        return omm;
    }

    auto initialize(const SystemTopology& systemTopology) -> bool;

    void shutdown() {
        integrator.reset();
        context.reset();
        system.reset();
        forceGroupLabels.clear();
        hasVirtualSites = false;
        initialized = false;
    }

    void setVelocitiesToTemperature(double temperature, int seed) {
        ensureInitialized();
        context->setVelocitiesToTemperature(temperature, seed);
    }

    [[nodiscard]] auto getPotentialEnergy() const -> double {
        ensureInitialized();
        return potentialEnergy;
    }
    [[nodiscard]] auto getKineticEnergy() const -> double {
        ensureInitialized();
        return kineticEnergy;
    }

    struct ForceGroupEnergy {
        int group = 0;
        std::string name;
        double energy = 0.0; // kJ/mol
    };

    void setSeparateForceGroups(bool enabled) {
        separateForceGroups = enabled;
    }
    [[nodiscard]] auto getSeparateForceGroups() const -> bool {
        return separateForceGroups;
    }

    // Multiple-timestep (r-RESPA) control. When enabled, initialize() builds an
    // MTSIntegrator instead of the single-rate VerletIntegrator: slow forces
    // (Nonbonded, GBSA, NBFIX) go to force group 0 and are evaluated once per
    // outer step; fast bonded forces go to group 1 and are evaluated
    // `innerSubsteps` times. Only affects the Cartesian world's on-device MD
    // (integrateTrajectory); the torsional worlds always get the full force sum.
    // MUST be called before initialize(). innerSubsteps < 2 disables MTS.
    void setMTS(bool enabled, int innerSubsteps) {
        useMTS = enabled && (innerSubsteps >= 2);
        mtsInnerSubsteps = innerSubsteps;
    }
    [[nodiscard]] auto getUseMTS() const -> bool {
        return useMTS;
    }

    void setPositions(const std::vector<OpenMM::Vec3>& positions) {
        ensureInitialized();
        context->setPositions(positions);
        if (hasVirtualSites) {
            context->computeVirtualSites();
        }
    }
    void getPositions(std::vector<OpenMM::Vec3>& out) const {
        ensureInitialized();
        out = context->getState(OpenMM::State::Positions, enforcePeriodicBox).getPositions();
    }

    // Whether OpenMM wraps coordinates into the primary box when we PULL state
    // (positions/forces/energy getState calls). For Robosample this MUST stay
    // false under explicit solvent: the robot engine consumes per-atom positions
    // and rebuilds each molecule's internal frames from contiguous geometry, so a
    // wrapped molecule that straddles a box face would yield a ~box-length "bond"
    // and corrupt the frame build. Energies/forces are unaffected by this flag
    // (OpenMM always applies the minimum image internally); wrapping is purely a
    // representation of the returned positions. Any periodic-image bookkeeping for
    // visualization is done as a WHOLE-MOLECULE rigid translation elsewhere, never
    // here. Default false.
    void setEnforcePeriodicBox(bool v) {
        enforcePeriodicBox = v;
    }
    [[nodiscard]] auto getEnforcePeriodicBox() const -> bool {
        return enforcePeriodicBox;
    }

    // True for the methods that require a periodic box (and thus exclude GBSA).
    [[nodiscard]] static auto isPeriodic(NonbondedMethod m) -> bool {
        return m == NonbondedMethod::CutoffPeriodic || m == NonbondedMethod::Ewald
               || m == NonbondedMethod::PME;
    }

    auto computePotentialEnergy(const std::vector<OpenMM::Vec3>& positions) -> double;
    auto computePotentialEnergyByGroup(const std::vector<OpenMM::Vec3>& positions)
        -> std::pair<double, std::vector<ForceGroupEnergy>>;
    void evaluateForcesFromPositionsCache(const std::vector<OpenMM::Vec3>& positions,
                                          std::vector<OpenMM::Vec3>& outForces) const;
    auto integrateTrajectory(int steps, double timeStepInPicoseconds) -> bool;
    auto computePeriodicBoxVectors_Context(double a_length,
                                           double b_length,
                                           double c_length,
                                           double alpha,
                                           double beta,
                                           double gamma)
        -> std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3>;

    private:
    OpenMMContext() = default;

    void ensureInitialized() const {
        if (!initialized) {
            throw std::runtime_error("OpenMM subsystem not initialized. Call OpenMMContext::initialize() "
                                     "before using any other functions.");
        }
    }

    [[nodiscard]] auto createNonbondedForce(const SystemTopology& systemTopology) -> OpenMM::NonbondedForce*;
    [[nodiscard]] auto createGBSAOBCForce(const SystemTopology& systemTopology) -> OpenMM::GBSAOBCForce*;
    [[nodiscard]] auto createCustomNonbondedForce(const SystemTopology& systemTopology)
        -> OpenMM::CustomNonbondedForce*;
    [[nodiscard]] auto createHarmonicBondForce(const SystemTopology& systemTopology)
        -> OpenMM::HarmonicBondForce*;
    [[nodiscard]] auto createHarmonicAngleForce(const SystemTopology& systemTopology)
        -> OpenMM::HarmonicAngleForce*;
    [[nodiscard]] auto createPeriodicTorsionForce(const SystemTopology& systemTopology)
        -> OpenMM::PeriodicTorsionForce*;
    [[nodiscard]] auto createImproperHarmonicTorsionForce(const SystemTopology& systemTopology)
        -> OpenMM::CustomTorsionForce*;
    [[nodiscard]] auto createCMAPTorsionForce(const SystemTopology& systemTopology)
        -> OpenMM::CMAPTorsionForce*;
    [[nodiscard]] auto createUreyBradleyForce(const SystemTopology& systemTopology)
        -> OpenMM::HarmonicBondForce*;

    std::unique_ptr<OpenMM::Context> context;
    std::unique_ptr<OpenMM::System> system;
    std::unique_ptr<OpenMM::Integrator> integrator;

    std::size_t numAtoms = 0;
    double potentialEnergy = 0;
    double kineticEnergy = 0;

    bool separateForceGroups = false;
    std::vector<std::pair<int, std::string>> forceGroupLabels;

    // r-RESPA multiple-timestep state.
    bool useMTS = false;
    int mtsInnerSubsteps = 4;
    static constexpr int kMtsSlowGroup = 0; // Nonbonded, GBSA, NBFIX
    static constexpr int kMtsFastGroup = 1; // bonds, angles, torsions, CMAP, UB

    bool enforcePeriodicBox = false;
    bool hasVirtualSites = false;
    bool initialized = false;
};