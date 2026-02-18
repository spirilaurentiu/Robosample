#pragma once

#include "OpenMM.h"

#if USE_CPU
    #include "../Molmodel/src/gbsa/cpuObcInterface.h"
    #include "../openmm/platforms/cpu/include/CpuPlatform.h"
#elif USE_OPENCL
    #include "../openmm/platforms/opencl/include/OpenCLPlatform.h"
#elif USE_CUDA
    #include "../openmm/platforms/cuda/include/CudaPlatform.h"
#endif

#include "TopologyElements.hpp"
#include <fstream>

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

enum class ForceGroup : int {
    HarmonicBondForce = 0,
    HarmonicAngleForce,
    PeriodicTorsionForce,
    ImproperTorsionForce,
    NonbondedForce,
    CustomNonbondedForce,
    Thermostat,
    CMAPTorsion,
    GBSAOBC,
    UreyBradley,
    Count
};

static inline std::string to_string(ForceGroup fg) {
    switch (fg) {
        case ForceGroup::HarmonicBondForce: return "HarmonicBondForce";
        case ForceGroup::HarmonicAngleForce: return "HarmonicAngleForce";
        case ForceGroup::PeriodicTorsionForce: return "PeriodicTorsionForce";
        case ForceGroup::ImproperTorsionForce: return "ImproperTorsionForce";
        case ForceGroup::NonbondedForce: return "NonbondedForce";
        case ForceGroup::CustomNonbondedForce: return "CustomNonbondedForce";
        case ForceGroup::Thermostat: return "AndersenThermostat";
        case ForceGroup::CMAPTorsion: return "CMAPTorsionForce";
        case ForceGroup::GBSAOBC: return "GBSAOBCForce";
        case ForceGroup::UreyBradley: return "UreyBradleyForce";
        default: return "UnknownForceGroup";
    }
}

struct ForceRegistration {
    OpenMM::Force* force = nullptr;
    ForceGroup group = ForceGroup::Count;
};

using OpenMMEnergyComponents = std::unordered_map<std::string, SimTK::Real>;

class OPENMM {
public:
    static bool initialize(
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
        SimTK::Real collisionFrequency,
        bool testing);

    OpenMMEnergyComponents getEnergyComponents();

    SimTK::Real getPotentialEnergy(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets);
    
    void setVelocitiesToTemperature(SimTK::Real temperature, uint32_t seed);

    static OPENMM& get() {
        static OPENMM omm;
        return omm;
    }

    static void shutdown() {
        OPENMM& omm = get();
        omm.destroy();
        omm.initialized = false;
    }

    SimTK::Real getPotentialEnergy() const;
    SimTK::Real getKineticEnergy() const;

    void setPositions(const std::vector<SimTK::Vec3> &positions);
    const std::vector<SimTK::Vec3>& getPositions() const;

    bool integrateTrajectory(const SimTK::Vector_<SimTK::Vec3>& positions, int steps);
    void restorePositions();

    void getEnergyAndForces(
        bool positionsAlreadySet,
        const std::vector<NonBondedMapping>& nonBondedMappings,
        const SimTK::Vector_<SimTK::Vec3>& includedAtomStation_G,
        const SimTK::Vector_<SimTK::Vec3>& includedAtomPos_G,
        SimTK::Vector_<SimTK::SpatialVec>& includedBodyForces_G,
        SimTK::Real &energy);

    std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3> computePeriodicBoxVectors_Context(
        double a_length, double b_length, double c_length,
        double alpha, double beta, double gamma);

private:
    OPENMM() = default;
    void destroy();

    void ensureInitialized() const {
        SimTK_ASSERT_ALWAYS(initialized, "OPENMM::ensureInitialized(): OpenMM has not initialized.");
    }

    void registerForce(const ForceRegistration& fr);
    
	std::unique_ptr<OpenMM::Context> context;
	std::unique_ptr<OpenMM::System> system;
	std::unique_ptr<OpenMM::Integrator> integrator;

    std::vector<ForceRegistration> forceRegistry;
    
    std::size_t numAtoms = 0;
    std::vector<OpenMM::Vec3> ommAtomsPositionsCache, ommAtomsPositionsCacheOld;
    std::vector<SimTK::Vec3> simbodyAtomsPositionsCache;
    SimTK::Real potentialEnergy = 0, kineticEnergy = 0;
    
    bool testing = false;    
	bool enforcePeriodicBox = false;
    bool initialized = false;
};
