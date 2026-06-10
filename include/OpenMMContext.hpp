#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

#include "OpenMM.h"
#include "TopologyElements.hpp"
#include "Vec3.h"

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
            throw std::runtime_error("OPENMM subsystem not initialized. Call OPENMM::initialize() before "
                                     "using any other functions.");
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

    std::vector<std::vector<int>> atomMbxByWorld;
    SystemTopology systemTopology;

    std::unique_ptr<OpenMM::Context> context;
    std::unique_ptr<OpenMM::System> system;
    std::unique_ptr<OpenMM::Integrator> integrator;

    std::size_t numAtoms = 0;
    double potentialEnergy = 0, kineticEnergy = 0;

    bool enforcePeriodicBox = false;
    bool initialized = false;
};
