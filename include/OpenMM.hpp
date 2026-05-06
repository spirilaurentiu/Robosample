#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "Force.h"
#include "OpenMM.h"
#include "TopologyElements.hpp"

enum class OpenMMForceType : std::uint8_t {
    Nonbonded = 0,
    GBSAOBC,
    CustomNonbonded,
    HarmonicBond,
    HarmonicAngle,
    PeriodicTorsion,
    ImproperHarmonicTorsion,
    CMAPTorsion,
    UreyBradley
};

struct NonbondedParameters {
    int index = std::numeric_limits<int>::max();
    SimTK::Real charge = SimTK::NaN;
    SimTK::Real sigma = SimTK::NaN;
    SimTK::Real epsilon = SimTK::NaN;
};

struct HarmonicBondParameters {
    int index = std::numeric_limits<int>::max();
    int particle1 = std::numeric_limits<int>::max();
    int particle2 = std::numeric_limits<int>::max();
    SimTK::Real length = SimTK::NaN;
    SimTK::Real k = SimTK::NaN;
};

struct HarmonicAngleParameters {
    int index = std::numeric_limits<int>::max();
    int particle1 = std::numeric_limits<int>::max();
    int particle2 = std::numeric_limits<int>::max();
    int particle3 = std::numeric_limits<int>::max();
    SimTK::Real angle = SimTK::NaN;
    SimTK::Real k = SimTK::NaN;
};

struct PeriodicTorsionParameters {
    int index = std::numeric_limits<int>::max();
    int particle1 = std::numeric_limits<int>::max();
    int particle2 = std::numeric_limits<int>::max();
    int particle3 = std::numeric_limits<int>::max();
    int particle4 = std::numeric_limits<int>::max();
    int periodicity = -1;
    SimTK::Real phase = SimTK::NaN;
    SimTK::Real k = SimTK::NaN;
};

struct HarmonicImproperTorsionParameters {
    int index = std::numeric_limits<int>::max();
    int atom1GlobalIndex = std::numeric_limits<int>::max();
    int atom2GlobalIndex = std::numeric_limits<int>::max();
    int atom3GlobalIndex = std::numeric_limits<int>::max();
    int atom4GlobalIndex = std::numeric_limits<int>::max();
    SimTK::Real stiffnessInKJPerRadSq = SimTK::NaN;
    SimTK::Real nominalAngleInRad = SimTK::NaN;
};

struct CMAPTorsionParameters {
    int index = std::numeric_limits<int>::max();
    int mapIndex = std::numeric_limits<int>::max();

    // Torsion A atoms
    int torsionAAtom1GlobalIndex = std::numeric_limits<int>::max();
    int torsionAAtom2GlobalIndex = std::numeric_limits<int>::max();
    int torsionAAtom3GlobalIndex = std::numeric_limits<int>::max();
    int torsionAAtom4GlobalIndex = std::numeric_limits<int>::max();

    // Torsion B atoms
    int torsionBAtom1GlobalIndex = std::numeric_limits<int>::max();
    int torsionBAtom2GlobalIndex = std::numeric_limits<int>::max();
    int torsionBAtom3GlobalIndex = std::numeric_limits<int>::max();
    int torsionBAtom4GlobalIndex = std::numeric_limits<int>::max();
};

struct UreyBradleyParameters {
    int index = std::numeric_limits<int>::max();
    int atom1GlobalIndex = std::numeric_limits<int>::max();
    int atom3GlobalIndex = std::numeric_limits<int>::max();
    SimTK::Real stiffnessInKJPerNmSq = SimTK::NaN;
    SimTK::Real nominalLengthInNm = SimTK::NaN;
};

inline auto makePairKey(int particle1, int particle2) -> uint64_t {
    auto low = static_cast<uint64_t>(std::min(particle1, particle2));
    auto high = static_cast<uint64_t>(std::max(particle1, particle2));
    return (high << 32) | low;
}

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

    void setActiveForceGroup(int forceGroupIndex);

    void setVelocitiesToTemperature(SimTK::Real temperature, int seed) {
        ensureInitialized();
        context->setVelocitiesToTemperature(temperature, seed);

        // Update velocities cache
        const auto state = context->getState(OpenMM::State::Velocities);
        ommAtomsVelocitiesCache = state.getVelocities();
    }

    [[nodiscard]] auto getPotentialEnergy() const -> SimTK::Real {
        ensureInitialized();
        return potentialEnergy;
    }

    [[nodiscard]] auto getKineticEnergy() const -> SimTK::Real {
        ensureInitialized();
        return kineticEnergy;
    }

    [[nodiscard]] auto getPositions() const -> const std::vector<SimTK::Vec3>& {
        ensureInitialized();
        return simbodyAtomsPositionsCache;
    }

    void negateVelocities() {
        ensureInitialized();
        const auto state = context->getState(OpenMM::State::Velocities);
        auto velocities = state.getVelocities();
        for (auto& vel : velocities) {
            vel = -vel;
        }
        context->setVelocities(velocities);
    }

    [[nodiscard]] auto getVelocities() const -> const std::vector<SimTK::Vec3>& {
        ensureInitialized();
        return simbodyAtomsVelocitiesCache;
    }

    [[nodiscard]] auto getPositionsCache() const -> const std::vector<OpenMM::Vec3>& {
        ensureInitialized();
        return ommAtomsPositionsCache;
    }

    [[nodiscard]] auto getVelocitiesCache() const -> const std::vector<OpenMM::Vec3>& {
        ensureInitialized();
        return ommAtomsVelocitiesCache;
    }

    void updatePositionsCache(const std::vector<NonBondedMapping>& nonBondedMappings,
                              const SimTK::Vector_<SimTK::Vec3>& inclAtomPos_G);

    void evaluateEnergiesFromPositionCache(SimTK::Real& newPotentialEnergy, SimTK::Real& newKineticEnergy);

    void evaluateForcesFromPositionsCache(const std::vector<NonBondedMapping>& nonBondedMappings,
                                          const SimTK::Vector_<SimTK::Vec3>& inclAtomStation_G,
                                          SimTK::Vector_<SimTK::SpatialVec>& inclBodyForces_G) const;

    [[nodiscard]] auto evaluatePotentialEnergyFromPositionsCache() const -> SimTK::Real;

    auto integrateTrajectory(int steps, SimTK::Real timeStepInPicoseconds) -> bool;
    auto integrateTrajectory(std::vector<OpenMM::Vec3>& positions,
                             int direction,
                             std::vector<OpenMM::Vec3>& velocities,
                             int steps,
                             SimTK::Real timeStepInPicoseconds,
                             SimTK::Real& potentialEnergy,
                             SimTK::Real& kineticEnergy) -> bool;

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

    auto isSameRigidBody(int wIx, int src, int dst) -> bool {
        const bool sameMbx = atomMbxByWorld[wIx][src] == atomMbxByWorld[wIx][dst];
        const bool validRigidBodySize = numAtomsInRigidBodiesByWorld[wIx][atomMbxByWorld[wIx][src]] > 1;
        return sameMbx && validRigidBodySize;
    }

    [[nodiscard]] auto computeIntraRigidPairs(int wIx, int minBondedDistance)
        -> std::vector<std::pair<int, int>>;

    [[nodiscard]] auto createNonbondedForce(const std::vector<RoboAtom>& atoms,
                                            const std::vector<Scaling14>& scaling14s,
                                            const std::vector<Exclusion>& exclusions,
                                            const ForceFieldParams& ffParams,
                                            bool hasNBfix) -> OpenMM::NonbondedForce*;
    [[nodiscard]] auto createGBSAOBCForce(const std::vector<RoboAtom>& atoms,
                                          const ForceFieldParams& ffParams) -> OpenMM::GBSAOBCForce*;
    [[nodiscard]] auto createCustomNonbondedForce(const std::vector<RoboAtom>& atoms,
                                                  const std::vector<Scaling14>& scaling14s,
                                                  const std::vector<Exclusion>& exclusions,
                                                  const ForceFieldParams& ffParams,
                                                  bool useSwitchingFunction,
                                                  SimTK::Real switchingDistance)
        -> OpenMM::CustomNonbondedForce*;
    [[nodiscard]] auto createHarmonicBondForce() -> OpenMM::HarmonicBondForce*;
    [[nodiscard]] auto createHarmonicAngleForce() -> OpenMM::HarmonicAngleForce*;
    [[nodiscard]] auto createPeriodicTorsionForce() -> OpenMM::PeriodicTorsionForce*;
    [[nodiscard]] auto createImproperHarmonicTorsionForce() -> OpenMM::CustomTorsionForce*;
    [[nodiscard]] auto createCMAPTorsionForce() -> OpenMM::CMAPTorsionForce*;
    [[nodiscard]] auto createUreyBradleyForce() -> OpenMM::HarmonicBondForce*;

    std::vector<std::vector<int>> atomMbxByWorld;
    SystemTopology systemTopology;
    ForceFieldParams ffParams;
    SimulationSettings simSettings;

    std::vector<std::unordered_map<int, int>> numAtomsInRigidBodiesByWorld;
    int numWorlds = 0;

    std::unique_ptr<OpenMM::Context> context;
    std::unique_ptr<OpenMM::System> system;
    std::unique_ptr<OpenMM::Integrator> integrator;

    std::vector<std::vector<NonbondedParameters>> nonbondedParamsByWorld;
    std::vector<std::vector<HarmonicBondParameters>> harmonicBondParamsByWorld;
    std::vector<std::vector<HarmonicAngleParameters>> harmonicAngleParamsByWorld;
    std::vector<std::vector<PeriodicTorsionParameters>> periodicTorsionParamsByWorld;
    std::vector<std::vector<HarmonicImproperTorsionParameters>> harmonicImproperTorsionParamsByWorld;
    std::vector<std::vector<CMAPTorsionParameters>> cmapTorsionParamsByWorld;
    std::vector<std::vector<UreyBradleyParameters>> ureyBradleyParamsByWorld;

    int nonbondedForceIndex = std::numeric_limits<int>::max();
    int harmonicBondForceIndex = std::numeric_limits<int>::max();
    int harmonicAngleForceIndex = std::numeric_limits<int>::max();
    int periodicTorsionForceIndex = std::numeric_limits<int>::max();
    int improperHarmonicTorsionForceIndex = std::numeric_limits<int>::max();
    int cmapTorsionForceIndex = std::numeric_limits<int>::max();
    int ureyBradleyForceIndex = std::numeric_limits<int>::max();

    std::size_t numAtoms = 0;
    std::vector<OpenMM::Vec3> ommAtomsPositionsCache, ommAtomsPositionsCacheOld;
    std::vector<OpenMM::Vec3> ommAtomsVelocitiesCache;
    std::vector<SimTK::Vec3> simbodyAtomsPositionsCache;
    std::vector<SimTK::Vec3> simbodyAtomsVelocitiesCache;
    SimTK::Real potentialEnergy = 0, kineticEnergy = 0;

    bool enforcePeriodicBox = false;
    bool initialized = false;
};
