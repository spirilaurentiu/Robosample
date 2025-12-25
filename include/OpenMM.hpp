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

class OPENMM {
public:
    static bool initialize(uint32_t seed, SimTK::Real sqrtCoulombScale, SimTK::Real vdwGlobalScaleFactor, const std::vector<Atom>& atoms, const std::vector<BondStretch>& bonds, const std::vector<BondBend>& angles, const std::vector<BondTorsion>& torsions);
    void addForceGroup();

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

private:
    OPENMM() = default;
    void destroy();

    void ensureInitialized() const {
        if (!initialized) {
            throw std::runtime_error("Config singleton used before init()");
        }
    }
    
    //OpenMMPluginInterface refOpenMMPlugin; // __refOMM__
	std::unique_ptr<OpenMM::Platform> platform; // __refOMM__
	std::unique_ptr<OpenMM::Context> openMMContext; // __refOMM__
	std::unique_ptr<OpenMM::System> openMMSystem; // __refOMM__
    
	std::unique_ptr<OpenMM::NonbondedForce> ommNonbondedForce; // __refOMM__
	std::unique_ptr<OpenMM::GBSAOBCForce> ommGBSAOBCForce; // __refOMM__
	std::unique_ptr<OpenMM::HarmonicBondForce> ommHarmonicBondStretch; // __refOMM__
	std::unique_ptr<OpenMM::HarmonicAngleForce> ommHarmonicAngleForce; // __refOMM__
	std::unique_ptr<OpenMM::PeriodicTorsionForce> ommPeriodicTorsionForce; // __refOMM__
    
	std::unique_ptr<OpenMM::AndersenThermostat> openMMThermostat;
	std::unique_ptr<OpenMM::Integrator> openMMIntegrator;
    
	OpenMM::State state;
    
    std::size_t numAtoms = 0;
    std::vector<OpenMM::Vec3> ommAtomsPositionsCache;
    std::vector<SimTK::Vec3> simbodyAtomsPositionsCache;
    SimTK::Real pe = 0, ke = 0;
    
	bool enforcePeriodicBox = false;
    
    SimTK::Real coulomb14Scale = 1/1.2; // From Amber force fields
	SimTK::Real lj14Scale = 0.5; // From Amber force fields

    bool initialized = false;
};
