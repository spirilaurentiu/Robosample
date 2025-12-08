#include "OpenMM.hpp"

// #define __PBC__

bool OPENMM::initialize(uint32_t seed, SimTK::Real sqrtCoulombScale, SimTK::Real vdwGlobalScaleFactor, const std::vector<Atom>& atoms, const std::vector<BondStretch>& bonds, const std::vector<BondBend>& angles, const std::vector<BondTorsion>& torsions) {

    // Instantiate
    OPENMM& omm = get();

	// Set up variables
	omm.numAtoms = atoms.size();
	omm.ommAtomsPositions = std::vector<OpenMM::Vec3>(omm.numAtoms);

    // Allocate OpenMM forces
	omm.ommHarmonicBondStretch = std::make_unique<OpenMM::HarmonicBondForce>();
	omm.ommHarmonicAngleForce = std::make_unique<OpenMM::HarmonicAngleForce>();
	omm.ommPeriodicTorsionForce = std::make_unique<OpenMM::PeriodicTorsionForce>();
	// omm.ommGBSAOBCForce = std::make_unique<OpenMM::GBSAOBCForce>();
	omm.ommNonbondedForce = std::make_unique<OpenMM::NonbondedForce>();
	
	// Instantiate the thermostat with adjusted temperature
	//Real temperature = 300.0;
	//if(dumm->wantOpenMMIntegration){temperature = dumm->temperature;}
	omm.openMMThermostat = std::make_unique<OpenMM::AndersenThermostat>(300.0, 1);
	omm.openMMThermostat->setRandomNumberSeed(seed);
	
	// Allocate OpenMM system and add particles to it
	omm.openMMSystem = std::make_unique<OpenMM::System>();

#ifdef __PBC__ // _pbc_

	double angle_alpha = 1.5708;
	double angle_beta = 1.5708;
	double angle_gamma = 1.5708;

	double boxLength_X = 10; // Example box length in angstroms
	double boxLength_Y = 10; // Example box length in angstroms
	double boxLength_Z = 10; // Example box length in angstroms

	auto periodicBoxVectors = computePeriodicBoxVectors_Context(
		boxLength_X, boxLength_Y, boxLength_Z,
		angle_alpha, angle_beta, angle_gamma);

	OpenMM::Vec3 pbcVector_X = std::get<0>(periodicBoxVectors);
	OpenMM::Vec3 pbcVector_Y = std::get<1>(periodicBoxVectors);
	OpenMM::Vec3 pbcVector_Z = std::get<2>(periodicBoxVectors);

	openMMSystem->setDefaultPeriodicBoxVectors(pbcVector_X, pbcVector_Y, pbcVector_Z);
# endif
	
	// Nonbonded forces
	// ommNonbondedForce->setNonbondedMethod( OpenMM::NonbondedForce::NonbondedMethod( nonbondedMethod ) );
	// ommNonbondedForce->setCutoffDistance( nonbondedCutoff );
	omm.ommNonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::NonbondedMethod::NoCutoff);
	omm.ommNonbondedForce->setCutoffDistance(1.0); // in nm
	omm.ommNonbondedForce->setUseDispersionCorrection(false);
	// nonbondedForce->setUseSwitchingFunction( 0 );

	// Add atoms
	for (const auto& atom : atoms) {
		const SimTK::Real charge = atom.getChargeInE() * sqrtCoulombScale;
		const SimTK::Real sigma = atom.getSigmaInNm();
		const SimTK::Real epsilon = atom.getVdwWellDepthInKJ() * vdwGlobalScaleFactor;

		omm.openMMSystem->addParticle(atom.getMassInDaltons());
		omm.ommNonbondedForce->addParticle(charge, sigma, epsilon);
		// ommGBSAOBCForce->addParticle(charge, )

		// std::cout << "Atom " << atom.getUniqueAtomName() << ": charge = " << charge << " e, sigma = " << sigma << " nm, epsilon = " << epsilon << " kJ/mol" << std::endl;
	}

	// Add bonds
	std::vector<std::pair<int, int>> ommBonds;
	for (const auto& bond : bonds) {
		ommBonds.emplace_back(std::make_pair(bond.getParentAtomGlobalIndex(), bond.getChildAtomGlobalIndex()));
	}

	// Register all the 1-2 bonds between nonbond atoms for scaling.
	omm.ommNonbondedForce->createExceptionsFromBonds(ommBonds, omm.coulomb14Scale * sqrtCoulombScale, omm.lj14Scale * vdwGlobalScaleFactor);
	
	// GBSA
	// When it is called for the i'th time, it specifies the parameters for the i'th particle.
	//ommGBSAOBCForce->setSolventDielectric(80.0);
	//ommGBSAOBCForce->setSoluteDielectric(1.0);
	// Watch the units here. OpenMM works exclusively in MD (nm, kJ/mol). 
	// CPU GBSA uses Angstrom, kCal/mol.
	// for (auto atom : atoms) {
	// 	SimTK::Real charge = atom.getChargeInE();
	// 	getGbsaRadii(int numberOfAtoms, const int* atomicNumber, 
	// 				const int* numberOfCovalentPartners, 
	// 				const int* atomicNumberOfHCovalentPartner, 
	// 				RealOpenMM* gbsaRadii);
	// 	ommGBSAOBCForce.addParticle((worlds[0].forceField)->gbsaAtomicPartialCharges[nax],
	// 								(worlds[0].forceField)->gbsaRadii[nax]*OpenMM::NmPerAngstrom,
	// 								(worlds[0].forceField)->gbsaObcScaleFactors[nax]); 
	// }
	//// System takes over heap ownership of the force.
	//openMMSystem.addForce(ommGBSAOBCForce.get());
	//ommGBSAOBCForce.release();
				
	for (const auto& bond : bonds) {
		const SimTK::Real nominalLengthInNm = bond.getNominalLengthInNm();
		const SimTK::Real stiffnessInKJPerNmSq = bond.getStiffnessInKJPerNmSq();

		// force constants are expressed for the full quadratic form but OpenMM interprets them as the prefactor of 1/2 k(x-x0)^2, hence the factor of 2 here
		omm.ommHarmonicBondStretch->addBond(bond.getParentAtomGlobalIndex(), bond.getChildAtomGlobalIndex(), nominalLengthInNm, stiffnessInKJPerNmSq * 2);
	}

	// FORCES: ADD ANGLES (1-2-3)
	for (const auto& angle : angles) {
		const int a1num = angle.getGlobalIndex1();
		const int a2num = angle.getGlobalIndex2();
		const int a3num = angle.getGlobalIndex3();
		const SimTK::Real theta0 = angle.getNominalAngleInDeg() * SimTK::DuMM::Deg2Rad;
		const SimTK::Real forceKt = angle.getStiffnessInKJPerRadSq();

		// force constants are expressed for the full quadratic form but OpenMM interprets them as the prefactor of 1/2 k(x-x0)^2, hence the factor of 2 here
		omm.ommHarmonicAngleForce->addAngle(a1num, a2num, a3num, theta0, forceKt * 2);
	}

	// Add dihedrals. OpenMM does not distinguish between proper and improper dihedrals.
	for (const auto& t : torsions) {
		const int a1 = t.getGlobalIndex1();
		const int a2 = t.getGlobalIndex2();
		const int a3 = t.getGlobalIndex3();
		const int a4 = t.getGlobalIndex4();
		
		if (t.getPhaseInDegrees_1() != -1) {
			omm.ommPeriodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_1(), t.getPhaseInDegrees_1() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_1());
		}
		if (t.getPhaseInDegrees_2() != -1) {
			omm.ommPeriodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_2(), t.getPhaseInDegrees_2() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_2());
		}
		if (t.getPhaseInDegrees_3() != -1) {
			omm.ommPeriodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_3(), t.getPhaseInDegrees_3() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_3());
		}
		if (t.getPhaseInDegrees_4() != -1) {
			omm.ommPeriodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_4(), t.getPhaseInDegrees_4() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_4());
		}
		if (t.getPhaseInDegrees_5() != -1) {
			omm.ommPeriodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_5(), t.getPhaseInDegrees_5() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_5());
		}
	}

	// const int group0 = 0;
	// const int group1 = 1;
	// const int group2 = 2;
	// const int group3 = 3;
	// const int group4 = 4;

	// omm.ommHarmonicBondStretch->setForceGroup(group0);
	// omm.ommHarmonicAngleForce->setForceGroup(group1);
	// omm.ommPeriodicTorsionForce->setForceGroup(group2);
	// omm.ommNonbondedForce->setForceGroup(group3);
	// omm.openMMThermostat->setForceGroup(group4);

	omm.openMMSystem->addForce(omm.ommHarmonicBondStretch.get()); omm.ommHarmonicBondStretch.release();
	omm.openMMSystem->addForce(omm.ommHarmonicAngleForce.get()); omm.ommHarmonicAngleForce.release();
	omm.openMMSystem->addForce(omm.ommPeriodicTorsionForce.get()); omm.ommPeriodicTorsionForce.release();
	omm.openMMSystem->addForce(omm.ommNonbondedForce.get()); omm.ommNonbondedForce.release();

	// Get the thermostat
	omm.openMMSystem->addForce(omm.openMMThermostat.get()); omm.openMMThermostat.release();
		
	// Get the integrator
	omm.openMMIntegrator = std::make_unique<OpenMM::VerletIntegrator>(0.0007); // TODO should release?
		
    // Get the platform
    // By default, OpenMM builds a .so for each platform (CPU, OpenCL and CUDA)
    // When loading that .so, two functions get called
    // 1. registerPlatform() which does what you see below
    // 2. registerKernelFactories() which is used for Drude, Pme, Rpmd and other plugins (which we do not need as of right now)
#if OPENMM_PLATFORM_CPU
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        omm.platform = std::make_unique<OpenMM::CpuPlatform>();
        OpenMM::Platform::registerPlatform(omm.platform.get());
        platform.release();
    }
    constexpr auto PLATFORM_NAME = "CPU";

#elif OPENMM_PLATFORM_CUDA
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        omm.platform = std::make_unique<OpenMM::CudaPlatform>();
        OpenMM::Platform::registerPlatform(omm.platform.get());
        omm.platform.release();
    }
    constexpr auto PLATFORM_NAME = "CUDA";

#elif OPENMM_PLATFORM_OPENCL
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        omm.platform = std::make_unique<OpenMM::OpenCLPlatform>();
        OpenMM::Platform::registerPlatform(omm.platform.get());
        omm.platform.release();
    }
    constexpr auto PLATFORM_NAME = "OpenCL";
#endif

	bool allowReferencePlatform = true;
    // CREATE OPENMM CONTEXT based on PLATFORM
    try {
        auto& platform = OpenMM::Platform::getPlatformByName(PLATFORM_NAME);
        omm.openMMContext = std::make_unique<OpenMM::Context>(*omm.openMMSystem, *omm.openMMIntegrator, platform);
        const double speed = omm.openMMContext->getPlatform().getSpeed();

        if (speed <= 1 && !allowReferencePlatform) {
            std::cout << "ERROR: OpenMM not used: best available platform was " << PLATFORM_NAME << " with relative speed " << speed << std::endl;
            std::cout << "ERROR: Call setAllowOpenMMReference() if you want to use this anyway." << std::endl;
            return "";
        }

        std::cout << "NOTE: Created OpenMM context with " << PLATFORM_NAME << " platform with relative speed " << speed << std::endl;


    } catch (const std::exception& e) {
        // Could not create this platform so log and try the next one
        std::cout << "ERROR: OpenMM error during initialization: " << e.what() << std::endl;
        return "";
    }

// 	std::cout << "Robosample Context reference OpenMM loaded " << omm.openMMContext->getPlatform().getName() << std::endl;


// ////////////////////////////////////////////////////////////////////////////////
// 	bool checkPotentialEnergyManually = true;
// 	if (checkPotentialEnergyManually){
// 		for (const auto& atom : atoms) {
// 			omm.ommAtomsPositions[atom.getGlobalIndex()] = OpenMM::Vec3(atom.getXInNm(), atom.getYInNm(), atom.getZInNm());
// 		}

// 		omm.openMMContext->setPositions(omm.ommAtomsPositions);
		
// 		// calculate total potential energy manually
// 		omm.state = omm.openMMContext->getState(OpenMM::State::Energy);
// 		std::cout << "Total energy is " << omm.state.getPotentialEnergy() << std::endl;

// 		omm.state = omm.openMMContext->getState(OpenMM::State::Energy, omm.enforcePeriodicBox, 1<<group0);
// 		std::cout << "HarmonicBondForce is " << omm.state.getPotentialEnergy() << std::endl;

// 		omm.state = omm.openMMContext->getState(OpenMM::State::Energy, omm.enforcePeriodicBox, 1<<group1);
// 		std::cout << "HarmonicAngleForce is " << omm.state.getPotentialEnergy() << std::endl;

// 		omm.state = omm.openMMContext->getState(OpenMM::State::Energy, omm.enforcePeriodicBox, 1<<group2);
// 		std::cout << "PeriodicTorsionForce is " << omm.state.getPotentialEnergy() << std::endl;

// 		omm.state = omm.openMMContext->getState(OpenMM::State::Energy, omm.enforcePeriodicBox, 1<<group3);
// 		std::cout << "NonbondedForce is " << omm.state.getPotentialEnergy() << std::endl;

// 		omm.state = omm.openMMContext->getState(OpenMM::State::Energy, omm.enforcePeriodicBox, 1<<group4);
// 		std::cout << "AndersenThermostat is " << omm.state.getPotentialEnergy() << std::endl;
// 	}

// 	// for (int i = 0; i < 1000; i++) {
// 	// 	openMMIntegrator->step(1000);
// 	// 	const auto& state = openMMContext->getState(OpenMM::State::Energy);
// 	// 	std::cout << "Step " << i * 1000 << " energy is " << state.getPotentialEnergy() << std::endl;
// 	// }

// 	// for (int step = 0; step < 1000; step++) {
// 	// 	openMMIntegrator->step(1000);

// 	// 	int group = step % 4;
// 	// 	OpenMM::State state = openMMContext->getState(OpenMM::State::Energy, group);
// 	// 	std::cout << "Step " << step * 1000 << " energy (group " << group << ") is " << state.getPotentialEnergy() << std::endl;
// 	// }

	omm.initialized = true;
	return true;
}

SimTK::Real OPENMM::getPotentialEnergy(const SimTK::Compound::AtomTargetLocations& atomTargets) {

	ensureInitialized();
	
	// Convert SimTK::Vec3 to OpenMM::Vec3
	for (const auto& atomTarget : atomTargets) {
		const SimTK::Vec3& coords = atomTarget.second;
		ommAtomsPositions[atomTarget.first] = OpenMM::Vec3(coords[0], coords[1], coords[2]);
	}

	// Set positions in OpenMM context
	openMMContext->setPositions(ommAtomsPositions);
	
	// Get state with energy
    return openMMContext->getState(OpenMM::State::Energy).getPotentialEnergy();
}

void OPENMM::setVelocitiesToTemperature(SimTK::Real temperature, uint32_t seed) {
	ensureInitialized();

	// openMMThermostat->setDefaultTemperature(temperature);
	// openMMContext->setVelocitiesToTemperature(temperature);
}

SimTK::Real OPENMM::getPotentialEnergy() const {
	ensureInitialized();
	return openMMContext->getState(OpenMM::State::Energy).getPotentialEnergy();
}

SimTK::Real OPENMM::getKineticEnergy() const {
	ensureInitialized();
	return openMMContext->getState(OpenMM::State::Energy).getKineticEnergy();
}

void OPENMM::setPositions(const std::vector<SimTK::Vec3> &positions) {
	ensureInitialized();

	// Convert SimTK::Vec3 to OpenMM::Vec3
	for (std::size_t i = 0; i < positions.size(); ++i) {
		const SimTK::Vec3& coords = positions[i];
		ommAtomsPositions[i] = OpenMM::Vec3(coords[0], coords[1], coords[2]);

		// std::cout << "OPENMM::setPositions(): Setting position of atom " << i << " to (" << coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;
	}
	// Set positions in OpenMM context
	openMMContext->setPositions(ommAtomsPositions);
}

const std::vector<OpenMM::Vec3>& OPENMM::getPositions() const {
	ensureInitialized();

	// Get state with positions
	OpenMM::State state = openMMContext->getState(OpenMM::State::Positions, enforcePeriodicBox);
	return state.getPositions();
}

void OPENMM::integrateTrajectory(int steps) {
	ensureInitialized();

	// Prin kinetic and potential energy before step
	std::cout << "Before step: Potential Energy = " << getPotentialEnergy() << " kJ/mol, Kinetic Energy = " << getKineticEnergy() << " kJ/mol" << std::endl;

	openMMIntegrator->step(steps);

	// Print kinetic and potential energy after step
	std::cout << "After step: Potential Energy = " << getPotentialEnergy() << " kJ/mol, Kinetic Energy = " << getKineticEnergy() << " kJ/mol" << std::endl;
}

// dumm->setOpenMMvelocities
// OMM_calcKineticEnergy
// dumm->OMM_integrateTrajectory
// OMM_calcPotentialEnergy
// dumm->OMM_getPositions()
// dumm->OMM_setOpenMMPositions()

// OMM_*