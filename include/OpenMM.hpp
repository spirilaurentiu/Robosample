#pragma once

#include <array>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>
#include <numeric>
#include <set>
#include <algorithm>

#include "Force.h"
#include "OpenMM.h"

#if USE_CPU
    #include "../Molmodel/src/gbsa/cpuObcInterface.h"
    #include "../openmm/platforms/cpu/include/CpuPlatform.h"
#elif USE_REFERENCE
    #include "../openmm/platforms/reference/include/ReferencePlatform.h"
#elif USE_OPENCL
    #include "../openmm/platforms/opencl/include/OpenCLPlatform.h"
#elif USE_CUDA
    #include "../openmm/platforms/cuda/include/CudaPlatform.h"
#endif

#include "TopologyElements.hpp"

enum NonbondedMethod : int {
    NoCutoff,
    CutoffNonPeriodic,
};

struct CMAPGrid {
	std::vector<SimTK::Real> energy;
	int size;
};

struct CMAPTorsion {
	int mapIndex;
	int a1, a2, a3, a4;
	int b1, b2, b3, b4;
};

struct UreyBradley {
    int a1, a3;
    SimTK::Real stiffnessInKJPerNmSq;
    SimTK::Real nominalLengthInNm;
};

struct Scaling14 {
    Scaling14() = default;
    Scaling14(
        int a1_, int a4_,
        SimTK::Real chargeProduct_, SimTK::Real epsilon_, SimTK::Real sigma_
    ) : a1(a1_), a4(a4_),
        chargeProduct(chargeProduct_), epsilon(epsilon_), sigma(sigma_)
    {}

    int a1, a4;
    SimTK::Real chargeProduct, epsilon, sigma;
};

struct Exclusion {
    Exclusion() = default;
    Exclusion(int a1_, int a2_) : a1(a1_), a2(a2_) {}

    int a1;
    int a2;
};

using CanonicalBond = std::pair<std::size_t, std::size_t>;
using CanonicalAngle = std::array<std::size_t, 3>;
using CanonicalTorsion = std::array<std::size_t, 4>;

[[nodiscard]] static inline std::pair<std::size_t, std::size_t> canonicalizeBond(std::size_t i, std::size_t j) noexcept {
    return { std::min(i, j), std::max(i, j) };
}

[[nodiscard]] static inline std::array<std::size_t, 3> canonicalizeAngle(std::size_t i, std::size_t j, std::size_t k) noexcept {
    return { std::min(i, k), j, std::max(i, k) };
}

[[nodiscard]] static inline std::array<std::size_t, 4> canonicalizeTorsion(std::size_t i, std::size_t j, std::size_t k, std::size_t l) noexcept {
    const std::array<std::size_t, 4> forward{i, j, k, l};
    const std::array<std::size_t, 4> reverse{l, k, j, i};
    return (forward < reverse) ? forward : reverse;
}

class ForceGroup {
public:
    void initialize(
        int forceGroupIndex,
        const std::vector<int>& rigidBodies,
        uint32_t seed,
		const std::vector<RoboAtom>& atoms,
		const std::vector<RoboBond>& bonds,
		const std::vector<RoboAngle>& angles,
		const std::vector<RoboPeriodicTorsion>& properPeriodicTorsions,
		const std::vector<RoboHarmonicImproperTorsion>& harmonicImproperTorsions,
		const std::vector<CMAPGrid>& cmapGrids,
		const std::vector<CMAPTorsion>& cmapTorsions,
        const std::vector<UreyBradley>& ureyBradleys,
		bool hasNBfix,
		int numTypes,
		const std::vector<SimTK::Real>& acoef,
		const std::vector<SimTK::Real>& bcoef,
		const std::vector<Exclusion>& exclusions,
        const std::vector<Scaling14>& scaling14s,
        bool useGBSAOBC2,
        SimTK::Real gbsaSolventDielectric,
        SimTK::Real gbsaSoluteDielectric,
        NonbondedMethod nonbondedMethod,
        SimTK::Real nonbondedCutoffInNm,
        SimTK::Real thermostatTemperature,
        SimTK::Real collisionFrequency
    );

    const std::vector<OpenMM::Force*>& getForces() const {
        return forces;
    }

private:
    void registerForce(OpenMM::Force* force) {
        forces.push_back(force);
        // force->setForceGroup(fg);
    }

    int fg = -1;
    std::vector<OpenMM::Force*> forces;
};

class OPENMM {
public:
    static bool initialize(
        uint32_t seed,
        const std::vector<std::vector<int>>& worlds,
		const std::vector<RoboAtom>& atoms,
		const std::vector<RoboBond>& bonds,
		const std::vector<RoboAngle>& angles,
		const std::vector<RoboPeriodicTorsion>& properPeriodicTorsions,
		const std::vector<RoboHarmonicImproperTorsion>& harmonicImproperTorsions,
		const std::vector<CMAPGrid>& cmapGrids,
		const std::vector<CMAPTorsion>& cmapTorsions,
        const std::vector<UreyBradley>& ureyBradleys,
        bool hasNBfix,
		int numTypes,
		const std::vector<SimTK::Real>& acoef,
		const std::vector<SimTK::Real>& bcoef,
        const std::vector<Exclusion>& exclusions,
        const std::vector<Scaling14>& scaling14s,
        bool useGBSAOBC2,
        SimTK::Real gbsaSolventDielectric,
        SimTK::Real gbsaSoluteDielectric,
        NonbondedMethod nonbondedMethod,
        SimTK::Real nonbondedCutoffInNm,
        SimTK::Real thermostatTemperature,
        SimTK::Real collisionFrequency
    );

    void setActiveForceGroup(int forceGroupIndex) {
        ensureInitialized();
        activeForceGroupIndex = forceGroupIndex;
    }

    void setVelocitiesToTemperature(SimTK::Real temperature, uint32_t seed);

    static OPENMM& get() {
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

    SimTK::Real getPotentialEnergy() const;
    SimTK::Real getKineticEnergy() const;

    void setPositions(const std::vector<SimTK::Vec3>& positions);
    const std::vector<SimTK::Vec3>& getPositions() const;

    bool integrateTrajectory(const SimTK::Vector_<SimTK::Vec3>& includedAtomPositionsInG, int steps, SimTK::Real timeStepInPicoseconds);

    void integrateTrajectory(const std::vector<SimTK::Vec3>& inPositions, std::vector<SimTK::Vec3>& outPositions, bool resetPositions, int steps, SimTK::Real timeStepInPicoseconds) {
        ensureInitialized();
        
        // Convert SimTK::Vec3 to OpenMM::Vec3
        if (resetPositions) {
            for (std::size_t i = 0; i < inPositions.size(); ++i) {
                const SimTK::Vec3& coords = inPositions[i];
                ommAtomsPositionsCache[i] = OpenMM::Vec3(coords[0], coords[1], coords[2]);
            }

            // Set positions in OpenMM context
            context->setPositions(ommAtomsPositionsCache);
        }

        // Integrate
        integrator->setIntegrationForceGroups(1 << activeForceGroupIndex);
		integrator->step(steps);

        // Return new positions
        const auto state = context->getState(OpenMM::State::Energy | OpenMM::State::Positions, enforcePeriodicBox, 1 << activeForceGroupIndex);
        const auto& pos = state.getPositions();
        outPositions.resize(pos.size());
        for (size_t i = 0, n = pos.size(); i < n; ++i) {
            const auto& c = pos[i];
            outPositions[i] = {c[0], c[1], c[2]};
        }

        potentialEnergy = state.getPotentialEnergy();
        kineticEnergy = state.getKineticEnergy();
    }

    void updatePositionsCache(const std::vector<NonBondedMapping>& nonBondedMappings, const SimTK::Vector_<SimTK::Vec3>& inclAtomPos_G);

    void evaluateForces(
        const std::vector<NonBondedMapping>& nonBondedMappings,
        const SimTK::Vector_<SimTK::Vec3>& inclAtomStation_G,
        SimTK::Vector_<SimTK::SpatialVec>& inclBodyForces_G) const;

    SimTK::Real evaluatePotentialEnergyFromPositionsCache() const;

    std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3> computePeriodicBoxVectors_Context(
        double a_length, double b_length, double c_length,
        double alpha, double beta, double gamma);

private:
    OPENMM() = default;

    void ensureInitialized() const {
        if (!initialized) {
            throw std::runtime_error("OPENMM subsystem not initialized. Call OPENMM::initialize() before using any other functions.");
        }
    }

	std::unique_ptr<OpenMM::Context> context;
	std::unique_ptr<OpenMM::System> system;
	std::unique_ptr<OpenMM::Integrator> integrator;
    std::vector<ForceGroup> forceGroups;
    
    std::size_t numAtoms = 0;
    std::vector<OpenMM::Vec3> ommAtomsPositionsCache, ommAtomsPositionsCacheOld;
    std::vector<SimTK::Vec3> simbodyAtomsPositionsCache;
    SimTK::Real potentialEnergy = 0, kineticEnergy = 0;
    
	bool enforcePeriodicBox = false;
    bool initialized = false;
    int activeForceGroupIndex = -1;
};
