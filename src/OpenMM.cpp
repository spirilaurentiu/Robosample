#include "OpenMM.hpp"

// #define __PBC__

void OPENMM::destroy() {
    ensureInitialized();

	// integrator.reset();
	// context.reset();
	// system.reset();
}

bool OPENMM::initialize(uint32_t seed,
						SimTK::Real sqrtCoulombScale,
						SimTK::Real vdwGlobalScaleFactor,
						const std::vector<RoboAtom>& atoms,
						const std::vector<RoboBondStretch>& bonds,
						const std::vector<RoboBondBend>& angles,
						const std::vector<RoboBondTorsion>& torsions,
						bool testing) {

    // Instantiate
    OPENMM& omm = get();

	// Set up variables
	omm.testing = testing;
	omm.numAtoms = atoms.size();
	omm.ommAtomsPositionsCache = std::vector<OpenMM::Vec3>(omm.numAtoms);
	omm.simbodyAtomsPositionsCache = std::vector<SimTK::Vec3>(omm.numAtoms);

    // Allocate OpenMM forces
	auto* harmonicBondStretch = new OpenMM::HarmonicBondForce();
	auto* harmonicAngleForce = new OpenMM::HarmonicAngleForce();
	auto* periodicTorsionForce = new OpenMM::PeriodicTorsionForce();
	auto* GBSAOBCForce = new OpenMM::GBSAOBCForce();
	auto* nonbondedForce = new OpenMM::NonbondedForce();
	
	// Instantiate the thermostat with adjusted temperature
	//Real temperature = 300.0;
	//if(dumm->wantOpenMMIntegration){temperature = dumm->temperature;}
	auto* thermostat = new OpenMM::AndersenThermostat(300.0, 1);
	thermostat->setDefaultTemperature(300.0);
	thermostat->setDefaultCollisionFrequency(1.0);
	thermostat->setRandomNumberSeed(seed);

	// Allocate OpenMM system and add particles to it
	omm.system = std::make_unique<OpenMM::System>();

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

	system->setDefaultPeriodicBoxVectors(pbcVector_X, pbcVector_Y, pbcVector_Z);
# endif
	
	// Nonbonded forces
	// nonbondedForce->setNonbondedMethod( OpenMM::NonbondedForce::NonbondedMethod( nonbondedMethod ) );
	// nonbondedForce->setCutoffDistance( nonbondedCutoff );
	nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::NonbondedMethod::NoCutoff);
	nonbondedForce->setCutoffDistance(1.0); // in nm
	nonbondedForce->setUseDispersionCorrection(false);
	// nonbondedForce->setUseSwitchingFunction( 0 );

	// Add atoms
	for (const auto& atom : atoms) {
		const SimTK::Real charge = atom.getChargeInE() * sqrtCoulombScale;
		const SimTK::Real sigma = atom.getSigmaInNm();
		const SimTK::Real epsilon = atom.getVdwWellDepthInKJ() * vdwGlobalScaleFactor;
		const SimTK::Real solventRadiusInNm = atom.getSolventRadiusInNm();
		const SimTK::Real screen = atom.getScreen();

		if (atom.isRoot()) {
			omm.system->addParticle(0.0); // massless root
		} else {
			omm.system->addParticle(atom.getMassInDaltons());
		}

		nonbondedForce->addParticle(charge, sigma, epsilon);
		GBSAOBCForce->addParticle(charge, solventRadiusInNm, screen);
	}

	// Add bonds
	std::vector<std::pair<int, int>> ommBonds;
	for (const auto& bond : bonds) {
		ommBonds.emplace_back(std::make_pair(bond.getParentAtomGlobalIndex(), bond.getChildAtomGlobalIndex()));
	}

	// Register all the 1-2 bonds between nonbond atoms for scaling.
	nonbondedForce->createExceptionsFromBonds(ommBonds, omm.coulomb14Scale * sqrtCoulombScale, omm.lj14Scale * vdwGlobalScaleFactor);
	
	for (const auto& bond : bonds) {
		const SimTK::Real nominalLengthInNm = bond.getNominalLengthInNm();
		const SimTK::Real stiffnessInKJPerNmSq = bond.getStiffnessInKJPerNmSq();

		// force constants are expressed for the full quadratic form but OpenMM interprets them as the prefactor of 1/2 k(x-x0)^2, hence the factor of 2 here
		harmonicBondStretch->addBond(bond.getParentAtomGlobalIndex(), bond.getChildAtomGlobalIndex(), nominalLengthInNm, stiffnessInKJPerNmSq * 2);
	}

	// FORCES: ADD ANGLES (1-2-3)
	for (const auto& angle : angles) {
		const int a1num = angle.getGlobalIndex1();
		const int a2num = angle.getGlobalIndex2();
		const int a3num = angle.getGlobalIndex3();
		const SimTK::Real theta0 = angle.getNominalAngleInDeg() * SimTK::DuMM::Deg2Rad;
		const SimTK::Real forceKt = angle.getStiffnessInKJPerRadSq();

		// Force constants are expressed for the full quadratic form but OpenMM interprets them as the prefactor of 1/2 k(x-x0)^2, hence the factor of 2 here
		harmonicAngleForce->addAngle(a1num, a2num, a3num, theta0, forceKt * 2);
	}

	// Add dihedrals. OpenMM does not distinguish between proper and improper dihedrals.
	for (const auto& t : torsions) {
		const int a1 = t.getGlobalIndex1();
		const int a2 = t.getGlobalIndex2();
		const int a3 = t.getGlobalIndex3();
		const int a4 = t.getGlobalIndex4();
		
		if (t.getPhaseInDegrees_1() != -1) {
			periodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_1(), t.getPhaseInDegrees_1() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_1());
		}
		if (t.getPhaseInDegrees_2() != -1) {
			periodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_2(), t.getPhaseInDegrees_2() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_2());
		}
		if (t.getPhaseInDegrees_3() != -1) {
			periodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_3(), t.getPhaseInDegrees_3() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_3());
		}
		if (t.getPhaseInDegrees_4() != -1) {
			periodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_4(), t.getPhaseInDegrees_4() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_4());
		}
		if (t.getPhaseInDegrees_5() != -1) {
			periodicTorsionForce->addTorsion(a1, a2, a3, a4, t.getPeriodicity_5(), t.getPhaseInDegrees_5() * SimTK::DuMM::Deg2Rad, t.getAmpInKJ_5());
		}
	}

	// Set force groups so that we can recover energy components later 
	if (omm.testing) {
		harmonicBondStretch->setForceGroup(omm.groupHarmonicBondStretch);
		harmonicAngleForce->setForceGroup(omm.groupHarmonicAngleForce);
		periodicTorsionForce->setForceGroup(omm.groupPeriodicTorsionForce);
		nonbondedForce->setForceGroup(omm.groupNonbondedForce);
		thermostat->setForceGroup(omm.groupThermostat);
		GBSAOBCForce->setForceGroup(omm.groupGBSAOBCForce);
	}

	// Add forces to system
	omm.system->addForce(harmonicBondStretch);
	omm.system->addForce(harmonicAngleForce);
	omm.system->addForce(periodicTorsionForce);
	omm.system->addForce(nonbondedForce);
	omm.system->addForce(thermostat);
	omm.system->addForce(GBSAOBCForce);

	// Get the integrator
	omm.integrator = std::make_unique<OpenMM::VerletIntegrator>(0.0007); // TODO should release?
		
    // Get the platform
	OpenMM::Platform* platform = nullptr;

#if OPENMM_PLATFORM_CPU
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        platform = new OpenMM::CpuPlatform();
        OpenMM::Platform::registerPlatform(platform);
    }
    constexpr auto PLATFORM_NAME = "CPU";

#elif OPENMM_PLATFORM_CUDA
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        platform = new OpenMM::CudaPlatform();
        OpenMM::Platform::registerPlatform(platform);
    }
    constexpr auto PLATFORM_NAME = "CUDA";

#elif OPENMM_PLATFORM_OPENCL
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        platform = new OpenMM::OpenCLPlatform();
        OpenMM::Platform::registerPlatform(platform);
    }
    constexpr auto PLATFORM_NAME = "OpenCL";
#endif

    try {
        omm.context = std::make_unique<OpenMM::Context>(*omm.system, *omm.integrator, *platform);

        const double speed = omm.context->getPlatform().getSpeed();
        std::cout << "Created OpenMM context with " << PLATFORM_NAME << " platform with relative speed " << speed << std::endl;

    } catch (const std::exception& e) {
        std::cout << "ERROR: OpenMM error during initialization: " << e.what() << std::endl;
        return false;
    }

	// All OpenMM components initialized successfully
	omm.initialized = true;

	// Set initial positions
	for (const auto& atom : atoms) {
		omm.ommAtomsPositionsCache[atom.getGlobalIndex()] = OpenMM::Vec3(atom.getXInNm(), atom.getYInNm(), atom.getZInNm());
	}
	omm.context->setPositions(omm.ommAtomsPositionsCache);

	return true;
}

OpenMMEnergyComponents OPENMM::getEnergyComponents() {
	ensureInitialized();

	SimTK_ASSERT_ALWAYS(testing, "OPENMM::getEnergyComponents() can only be called in testing mode.");

	OpenMMEnergyComponents components;

	// Total energy
	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox);
	components.totalEnergy = state.getPotentialEnergy();

	// Harmonic Bond Force
	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupHarmonicBondStretch);
	components.harmonicBondForce = state.getPotentialEnergy();

	// Harmonic Angle Force
	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupHarmonicAngleForce);
	components.harmonicAngleForce = state.getPotentialEnergy();

	// Periodic Torsion Force
	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupPeriodicTorsionForce);
	components.periodicTorsionForce = state.getPotentialEnergy();

	// Nonbonded Force
	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupNonbondedForce);
	components.nonbondedForce = state.getPotentialEnergy();

	// Thermostat
	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupThermostat);
	components.andersenThermostat = state.getPotentialEnergy();

	// GBSA OBC Force
	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupGBSAOBCForce);
	components.gbsaObcForce = state.getPotentialEnergy();

	return components;
}

SimTK::Real OPENMM::getPotentialEnergy(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets) {

	ensureInitialized();
	
	// Convert SimTK::Vec3 to OpenMM::Vec3
	int i = 0;
	for (const auto& topology : atomTargets) {
		for (const auto& atomTarget : topology) {
			const SimTK::Vec3& coords = atomTarget.second;
			ommAtomsPositionsCache[i] = OpenMM::Vec3(coords[0], coords[1], coords[2]);
			i++;
		}
	}

	// Set positions in OpenMM context
	context->setPositions(ommAtomsPositionsCache);
	
	// Get state with energy
    return context->getState(OpenMM::State::Energy).getPotentialEnergy();
}

void OPENMM::setVelocitiesToTemperature(SimTK::Real temperature, uint32_t seed) {
	ensureInitialized();
	context->setVelocitiesToTemperature(temperature, seed);
	context->setParameter(OpenMM::AndersenThermostat::Temperature(), temperature);
}

SimTK::Real OPENMM::getPotentialEnergy() const {
	ensureInitialized();
	return pe;
}

SimTK::Real OPENMM::getKineticEnergy() const {
	ensureInitialized();
	return ke;
}

void OPENMM::setPositions(const std::vector<SimTK::Vec3> &positions) {
	ensureInitialized();

	// Convert SimTK::Vec3 to OpenMM::Vec3
	for (std::size_t i = 0; i < positions.size(); ++i) {
		const SimTK::Vec3& coords = positions[i];
		ommAtomsPositionsCache[i] = OpenMM::Vec3(coords[0], coords[1], coords[2]);
	}

	// Set positions in OpenMM context
	context->setPositions(ommAtomsPositionsCache);
}

const std::vector<SimTK::Vec3>& OPENMM::getPositions() const {
	ensureInitialized();
	return simbodyAtomsPositionsCache;
}

void OPENMM::integrateTrajectory(int steps) {
	ensureInitialized();

	integrator->step(steps);

	// Cache the results
	state = context->getState(OpenMM::State::Positions | OpenMM::State::Energy | OpenMM::State::Forces, enforcePeriodicBox);

	// Save the energies
	pe = state.getPotentialEnergy();
	ke = state.getKineticEnergy();

	// Save the positions
	for (std::size_t i = 0; i < simbodyAtomsPositionsCache.size(); ++i) {
		simbodyAtomsPositionsCache[i] = SimTK::Vec3(
			state.getPositions()[i][0],
			state.getPositions()[i][1],
			state.getPositions()[i][2]
		);
	}

	// Save the forces
	ommForcesCache = state.getForces();
}

void OPENMM::getEnergyAndForces(
        const std::vector<std::size_t>& nax2daix,
        const std::vector<std::size_t>& nax2iax,
        const std::vector<std::size_t>& nax2ibx,
        const SimTK::Vector_<SimTK::Vec3>& includedAtomStation_G,
        const SimTK::Vector_<SimTK::Vec3>& includedAtomPos_G,
        SimTK::Vector_<SimTK::SpatialVec>& includedBodyForces_G,
        SimTK::Real &energy)
{
	ensureInitialized();
	SimTK_ASSERT_ALWAYS(nax2iax.size() == nax2ibx.size(),
		"OPENMM::getEnergyAndForces(): nax2iax and nax2ibx must have the same size.");
	SimTK_ASSERT_ALWAYS(includedAtomPos_G.size() == numAtoms,
		"OPENMM::getEnergyAndForces(): includedAtomPos_G size must match number of atoms.");

	// std::cout << "\tOPENMM::getEnergyAndForces(): numAtoms = " << numAtoms << ", nax2iax.size() = " << nax2iax.size() << std::endl;

	// Set positions in OpenMM context
	for (std::size_t nax = 0; nax < nax2iax.size(); ++nax) {
		const std::size_t dAIx = nax2daix[nax];
		const std::size_t iax = nax2iax[nax];
		const SimTK::Vec3& pos_G = includedAtomPos_G[iax];

		// std::cout << "\tnax " << nax << ": dAIx = " << dAIx << ", iax = " << iax
		// 	<< ", pos_G = (" << pos_G[0] << ", " << pos_G[1] << ", " << pos_G[2] << ")" << std::endl;

		ommAtomsPositionsCache[dAIx] = OpenMM::Vec3(pos_G[0], pos_G[1], pos_G[2]);
	}
	context->setPositions(ommAtomsPositionsCache);

	// bool checkPotentialEnergyManually = true;
	// if (checkPotentialEnergyManually){
	// 	const int groupHarmonicBondStretch = 0;
	// 	const int groupHarmonicAngleForce = 1;
	// 	const int groupPeriodicTorsionForce = 2;
	// 	const int groupNonbondedForce = 3;
	// 	const int groupThermostat = 4;
	// 	const int groupGBSAOBCForce = 5;
		
	// 	// calculate total potential energy manually
	// 	state = context->getState(OpenMM::State::Energy);
	// 	std::cout << "Total energy is " << state.getPotentialEnergy() << std::endl;

	// 	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupHarmonicBondStretch);
	// 	std::cout << "HarmonicBondForce is " << state.getPotentialEnergy() << std::endl;

	// 	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupHarmonicAngleForce);
	// 	std::cout << "HarmonicAngleForce is " << state.getPotentialEnergy() << std::endl;

	// 	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupPeriodicTorsionForce);
	// 	std::cout << "PeriodicTorsionForce is " << state.getPotentialEnergy() << std::endl;

	// 	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupNonbondedForce);
	// 	std::cout << "NonbondedForce is " << state.getPotentialEnergy() << std::endl;

	// 	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupThermostat);
	// 	std::cout << "AndersenThermostat is " << state.getPotentialEnergy() << std::endl;

	// 	state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<groupGBSAOBCForce);
	// 	std::cout << "GBSAOBCForce is " << state.getPotentialEnergy() << std::endl;
	// }

	// SimTK_ASSERT_ALWAYS(false, "checkpoint set positions");


	// Get state with energy and forces
	state = context->getState(OpenMM::State::Positions | OpenMM::State::Energy | OpenMM::State::Forces, enforcePeriodicBox);
	pe = state.getPotentialEnergy();
	ommForcesCache = state.getForces();

	// Return energy
	energy += pe;

	// Map forces from atoms to bodies
    for (std::size_t nax = 0; nax < nax2ibx.size(); ++nax)
    {
		const std::size_t dAIx = nax2daix[nax];
		const std::size_t iax = nax2iax[nax];
		const std::size_t ibx = nax2ibx[nax];

    	const SimTK::Vec3 simForce(ommForcesCache[dAIx][0], ommForcesCache[dAIx][1], ommForcesCache[dAIx][2]);
    	includedBodyForces_G[ibx] += SimTK::SpatialVec(includedAtomStation_G[iax] % simForce, simForce);
    }
}

std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3> OPENMM::computePeriodicBoxVectors_Context(
	double a_length, double b_length, double c_length,
    double alpha, double beta, double gamma)
{
	ensureInitialized();
	
	const double TOL = 1e-6;

	// // Convert angles from degrees to radians
	// alpha = SimTK::Deg2Rad * alpha;
	// beta  = SimTK::Deg2Rad * beta;
	// gamma = SimTK::Deg2Rad * gamma;
    
		// Compute the box vectors
    OpenMM::Vec3 a(a_length, 0.0, 0.0);

    OpenMM::Vec3 b(b_length * std::cos(gamma),
           b_length * std::sin(gamma),
           0.0);

    double cx = c_length * std::cos(beta);
    double cy = c_length * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
    double cz = std::sqrt(c_length * c_length - cx * cx - cy * cy);

    OpenMM::Vec3 c(cx, cy, cz);

    // Zero out small components
    for (int i = 0; i < 3; i++) {
        if (std::abs(a[i]) < TOL) a[i] = 0.0;
        if (std::abs(b[i]) < TOL) b[i] = 0.0;
        if (std::abs(c[i]) < TOL) c[i] = 0.0;
    }

    // Reduced form (OpenMM requirement)
    if (b[1] != 0.0)
        c -= b * std::round(c[1] / b[1]);
    if (a[0] != 0.0)
        c -= a * std::round(c[0] / a[0]);
    if (a[0] != 0.0)
        b -= a * std::round(b[0] / a[0]);

    return std::make_tuple(a, b, c);
}