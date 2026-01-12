#pragma once

#include "OpenMM.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/internal/ThreadPool.h"

#if OPENMM_PLATFORM_CPU
    #include "../Molmodel/src/gbsa/cpuObcInterface.h"
    #include "../openmm/platforms/cpu/include/CpuPlatform.h"
#elif OPENMM_PLATFORM_OPENCL
    #include "../openmm/platforms/opencl/include/OpenCLPlatform.h"
#elif OPENMM_PLATFORM_CUDA
    #include "../openmm/platforms/cuda/include/CudaPlatform.h"
#endif

#include "TopologyElements.hpp"
#include <fstream>

struct OpenMMEnergyComponents {
    SimTK::Real totalEnergy = 0.0;
    SimTK::Real harmonicBondForce = 0.0;
    SimTK::Real harmonicAngleForce = 0.0;
    SimTK::Real periodicTorsionForce = 0.0;
    SimTK::Real nonbondedForce = 0.0;
    SimTK::Real andersenThermostat = 0.0;
    SimTK::Real gbsaObcForce = 0.0;
};

class OPENMM {
public:
    static bool initialize(uint32_t seed, SimTK::Real sqrtCoulombScale, SimTK::Real vdwGlobalScaleFactor, const std::vector<RoboAtom>& atoms, const std::vector<RoboBondStretch>& bonds, const std::vector<RoboBondBend>& angles, const std::vector<RoboBondTorsion>& torsions, bool testing);

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

    void integrateTrajectory(int steps);

    void getEnergyAndForces(
        const std::vector<std::size_t>& nax2daix,
        const std::vector<std::size_t>& nax2iax,
        const std::vector<std::size_t>& nax2ibx,
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
        SimTK_ASSERT_ALWAYS(initialized, "OpenMM subsystem not initialized.");
    }
    
	std::unique_ptr<OpenMM::Context> context;
	std::unique_ptr<OpenMM::System> system;
	std::unique_ptr<OpenMM::Integrator> integrator;
    
	OpenMM::State state;
    
    std::size_t numAtoms = 0;
    std::vector<OpenMM::Vec3> ommAtomsPositionsCache;
    std::vector<OpenMM::Vec3> ommForcesCache;
    std::vector<SimTK::Vec3> simbodyAtomsPositionsCache;
    SimTK::Real pe = 0, ke = 0;
    
    SimTK::Real coulomb14Scale = 1/1.2; // From Amber force fields
	SimTK::Real lj14Scale = 0.5; // From Amber force fields
    
    
    int groupHarmonicBondStretch = 0;
	int groupHarmonicAngleForce = 1;
	int groupPeriodicTorsionForce = 2;
	int groupNonbondedForce = 3;
	int groupThermostat = 4;
	int groupGBSAOBCForce = 5;
    
    bool testing = false;    
	bool enforcePeriodicBox = false;
    bool initialized = false;
};
