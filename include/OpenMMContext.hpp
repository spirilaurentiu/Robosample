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

    void evaluateForcesFromPositionsCache(const std::vector<OpenMMContext::Vec3>& positions,
                                          std::vector<OpenMMContext::Vec3>& outForces) const;

    auto integrateTrajectory(int steps, double timeStepInPicoseconds) -> bool;

    auto computePeriodicBoxVectors_Context(double a_length,
                                           double b_length,
                                           double c_length,
                                           double alpha,
                                           double beta,
                                           double gamma)
        -> std::tuple<OpenMMContext::Vec3, OpenMMContext::Vec3, OpenMMContext::Vec3>;

    private:
    OpenMMContext() = default;

    void ensureInitialized() const {
        if (!initialized) {
            throw std::runtime_error(
                "OPENMM subsystem not initialized. Call OpenMMContext::initialize() before "
                "using any other functions.");
        }
    }

    [[nodiscard]] auto createNonbondedForce(const SystemTopology& systemTopology)
        -> OpenMMContext::NonbondedForce*;
    [[nodiscard]] auto createGBSAOBCForce(const SystemTopology& systemTopology)
        -> OpenMMContext::GBSAOBCForce*;
    [[nodiscard]] auto createCustomNonbondedForce(const SystemTopology& systemTopology)
        -> OpenMMContext::CustomNonbondedForce*;
    [[nodiscard]] auto createHarmonicBondForce(const SystemTopology& systemTopology)
        -> OpenMMContext::HarmonicBondForce*;
    [[nodiscard]] auto createHarmonicAngleForce(const SystemTopology& systemTopology)
        -> OpenMMContext::HarmonicAngleForce*;
    [[nodiscard]] auto createPeriodicTorsionForce(const SystemTopology& systemTopology)
        -> OpenMMContext::PeriodicTorsionForce*;
    [[nodiscard]] auto createImproperHarmonicTorsionForce(const SystemTopology& systemTopology)
        -> OpenMMContext::CustomTorsionForce*;
    [[nodiscard]] auto createCMAPTorsionForce(const SystemTopology& systemTopology)
        -> OpenMMContext::CMAPTorsionForce*;
    [[nodiscard]] auto createUreyBradleyForce(const SystemTopology& systemTopology)
        -> OpenMMContext::HarmonicBondForce*;

    std::vector<std::vector<int>> atomMbxByWorld;
    SystemTopology systemTopology;

    std::unique_ptr<OpenMMContext::Context> context;
    std::unique_ptr<OpenMMContext::System> system;
    std::unique_ptr<OpenMMContext::Integrator> integrator;

    std::size_t numAtoms = 0;
    double potentialEnergy = 0, kineticEnergy = 0;

    bool enforcePeriodicBox = false;
    bool initialized = false;
};
