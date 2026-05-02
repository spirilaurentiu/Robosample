#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "Force.h"
#include "OpenMM.h"
#include "Scalar.h"
#include "TopologyElements.hpp"

class ForceGroup {
    public:
    void initialize(int forceGroupIndex,
                    const std::vector<int>& rigidBodies,
                    const SystemTopology& systemTopology,
                    const ForceFieldParams& ffParams,
                    const SimulationSettings& simSettings);

    [[nodiscard]] auto getForces() const -> const std::vector<OpenMM::Force*>& {
        return forces;
    }

    private:
    void registerForce(OpenMM::Force* force) {
        forces.push_back(force);
        force->setForceGroup(fg);
    }

    int fg = std::numeric_limits<int>::max();
    std::vector<OpenMM::Force*> forces;
};

class OPENMM {
    public:
    static auto initialize(const std::vector<std::vector<int>>& worlds,
                           const SystemTopology& systemTopology,
                           const ForceFieldParams& ffParams,
                           const SimulationSettings& simSettings) -> bool;

    static auto get() -> OPENMM& {
        static OPENMM omm;
        return omm;
    }

    static void shutdown() {
        OPENMM& omm = get();

        omm.integrator.reset();
        omm.context.reset();
        omm.system.reset();
        omm.initialized = false;
    }

    void setActiveForceGroup(int forceGroupIndex) {
        ensureInitialized();
        activeForceGroupIndex = forceGroupIndex;
    }

    void setVelocitiesToTemperature(SimTK::Real temperature, int seed) {
        ensureInitialized();
        context->setVelocitiesToTemperature(temperature, seed);
        context->setParameter(OpenMM::AndersenThermostat::Temperature(), temperature);
    }

    [[nodiscard]] auto getPotentialEnergy() const -> SimTK::Real {
        ensureInitialized();
        return potentialEnergy;
    }

    [[nodiscard]] auto getKineticEnergy() const -> SimTK::Real {
        ensureInitialized();
        return kineticEnergy;
    }

    void setPositions(const std::vector<SimTK::Vec3>& positions) {
        ensureInitialized();

        // Convert SimTK::Vec3 to OpenMM::Vec3
        for (std::size_t i = 0; i < positions.size(); ++i) {
            const SimTK::Vec3& coords = positions[i];
            ommAtomsPositionsCache[i] = OpenMM::Vec3(coords[0], coords[1], coords[2]);
        }

        // Set positions in OpenMM context
        context->setPositions(ommAtomsPositionsCache);
    }

    [[nodiscard]] auto getPositions() const -> const std::vector<SimTK::Vec3>& {
        ensureInitialized();
        return simbodyAtomsPositionsCache;
    }

    void negateVelocities() {
        ensureInitialized();
        const auto state = context->getState(OpenMM::State::Velocities);
        auto velocities = state.getVelocities();
        for (auto& v : velocities) {
            v = -v;
        }
        context->setVelocities(velocities);
    }

    [[nodiscard]] auto getVelocities() const -> const std::vector<SimTK::Vec3>& {
        ensureInitialized();
        return simbodyAtomsVelocitiesCache;
    }

    void updatePositionsCache(const std::vector<NonBondedMapping>& nonBondedMappings,
                              const SimTK::Vector_<SimTK::Vec3>& inclAtomPos_G);

    void evaluateEnergiesFromPositionCache(SimTK::Real& newPotentialEnergy, SimTK::Real& newKineticEnergy);

    void evaluateForcesFromPositionsCache(const std::vector<NonBondedMapping>& nonBondedMappings,
                                          const SimTK::Vector_<SimTK::Vec3>& inclAtomStation_G,
                                          SimTK::Vector_<SimTK::SpatialVec>& inclBodyForces_G) const;

    auto evaluatePotentialEnergyFromPositionsCache() const -> SimTK::Real;

    auto integrateTrajectory(int steps, SimTK::Real timeStepInPicoseconds) -> bool;

    auto computePeriodicBoxVectors_Context(SimTK::Real a_length,
                                           SimTK::Real b_length,
                                           SimTK::Real c_length,
                                           SimTK::Real alpha,
                                           SimTK::Real beta,
                                           SimTK::Real gamma)
        -> std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3>;

    private:
    OPENMM() = default;

    void ensureInitialized() const {
        if (!initialized) {
            throw std::runtime_error("OPENMM subsystem not initialized. Call OPENMM::initialize() before "
                                     "using any other functions.");
        }
    }

    std::unique_ptr<OpenMM::Context> context;
    std::unique_ptr<OpenMM::System> system;
    std::unique_ptr<OpenMM::Integrator> integrator;
    std::vector<ForceGroup> forceGroups;

    std::size_t numAtoms = 0;
    std::vector<OpenMM::Vec3> ommAtomsPositionsCache, ommAtomsPositionsCacheOld;
    std::vector<SimTK::Vec3> simbodyAtomsPositionsCache;
    std::vector<SimTK::Vec3> simbodyAtomsVelocitiesCache;
    SimTK::Real potentialEnergy = 0, kineticEnergy = 0;

    bool enforcePeriodicBox = false;
    bool initialized = false;
    int activeForceGroupIndex = -1;
};
