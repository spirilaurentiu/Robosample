#include "OpenMM.hpp"

// #define __PBC__

void OPENMM::destroy() {
    ensureInitialized();

	// integrator.reset();
	// context.reset();
	// system.reset();
}

void addExclusionsToSet(const std::vector<std::set<int>>& bonded12, std::set<int>& exclusions, int baseParticle, int fromParticle, int currentLevel) {
    for (int i : bonded12[fromParticle]) {
        if (i != baseParticle)
            exclusions.insert(i);
        if (currentLevel > 0)
            addExclusionsToSet(bonded12, exclusions, baseParticle, i, currentLevel-1);
    }
}

bool OPENMM::initialize(
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
        bool testing) {

	OpenMM::NonbondedForce::NonbondedMethod nonbondedForceMethod;
	OpenMM::CustomNonbondedForce::NonbondedMethod customNonbondedForceMethod;
	OpenMM::GBSAOBCForce::NonbondedMethod gbsaForceMethod;

	switch (nonbondedMethod)
	{
	case NonbondedMethod::NoCutoff:
		nonbondedForceMethod = OpenMM::NonbondedForce::NoCutoff;
		customNonbondedForceMethod = OpenMM::CustomNonbondedForce::NoCutoff;
		gbsaForceMethod = OpenMM::GBSAOBCForce::NoCutoff;
		break;
	case NonbondedMethod::CutoffNonPeriodic:
		nonbondedForceMethod = OpenMM::NonbondedForce::CutoffNonPeriodic;
		customNonbondedForceMethod = OpenMM::CustomNonbondedForce::CutoffNonPeriodic;
		gbsaForceMethod = OpenMM::GBSAOBCForce::CutoffNonPeriodic;
		break;
	default:
		SimTK_ASSERT_ALWAYS(false, "OPENMM::initialize: Unsupported nonbonded method.");
	}

    // Instantiate
    OPENMM& omm = get();

	// Set up variables
	omm.testing = testing;
	omm.numAtoms = atoms.size();
	omm.ommAtomsPositionsCache = std::vector<OpenMM::Vec3>(omm.numAtoms);
	omm.ommAtomsPositionsCacheOld = std::vector<OpenMM::Vec3>(omm.numAtoms);
	omm.simbodyAtomsPositionsCache = std::vector<SimTK::Vec3>(omm.numAtoms);

	// Allocate OpenMM system and add particles to it
	omm.system = std::make_unique<OpenMM::System>();

	// Instantiate the thermostat with adjusted temperature
	auto* thermostat = new OpenMM::AndersenThermostat(thermostatTemperature, collisionFrequency);
	thermostat->setRandomNumberSeed(seed);
	omm.registerForce({thermostat, ForceGroup::Thermostat});
	
	// Nonbonded forces
	auto* nonbondedForce = new OpenMM::NonbondedForce();
	nonbondedForce->setNonbondedMethod(nonbondedForceMethod);
	nonbondedForce->setCutoffDistance(nonbondedCutoffInNm);

	if (nonbondedForceMethod == OpenMM::NonbondedForce::CutoffNonPeriodic) {
		nonbondedForce->setUseDispersionCorrection(true);
	} else {
		nonbondedForce->setUseDispersionCorrection(false);
	}

	// 1-4 VdW Correction (Bond-based so we can target specific pairs)
	OpenMM::CustomNonbondedForce* customNonbonded = nullptr;
	if (hasNBfix) {
		customNonbonded = new OpenMM::CustomNonbondedForce("(a/r6)^2-b/r6; r6=r^6; a=acoef(type1, type2); b=bcoef(type1, type2);");
		customNonbonded->addTabulatedFunction("acoef", new OpenMM::Discrete2DFunction(numTypes, numTypes, acoef));
		customNonbonded->addTabulatedFunction("bcoef", new OpenMM::Discrete2DFunction(numTypes, numTypes, bcoef));
		customNonbonded->addPerParticleParameter("type");

		customNonbonded->setNonbondedMethod(customNonbondedForceMethod);
		customNonbonded->setCutoffDistance(nonbondedCutoffInNm);

		customNonbonded->setUseSwitchingFunction(nonbondedForce->getUseSwitchingFunction());
		customNonbonded->setSwitchingDistance(nonbondedForce->getSwitchingDistance());
	}

	// GBSA OBC Force
	// implicitSolventKappa
	OpenMM::GBSAOBCForce* GBSAOBCForce = nullptr;
	if (useGBSAOBC2) {
		GBSAOBCForce = new OpenMM::GBSAOBCForce();
		GBSAOBCForce->setSolventDielectric(gbsaSolventDielectric);
		GBSAOBCForce->setSoluteDielectric(gbsaSoluteDielectric);
		GBSAOBCForce->setNonbondedMethod(gbsaForceMethod);
		GBSAOBCForce->setCutoffDistance(nonbondedCutoffInNm);

		/*
		 * The following comment is copy-pasted from the OpenMM documentation for GBSAOBCForce (note that we only support CutoffNonPeriodic):
		 * When using GBSAOBCForce, the System should also include a NonbondedForce, and both objects must specify
		 * identical charges for all particles. Otherwise, the results will not be correct. Furthermore, if the
		 * nonbonded method is set to CutoffNonPeriodic or CutoffPeriodic, you should call setReactionFieldDielectric(1.0)
		 * on the NonbondedForce to turn off the reaction field approximation, which does not produce correct results
		 * when combined with GBSA.
		 */
		if (nonbondedForceMethod == OpenMM::NonbondedForce::CutoffNonPeriodic) {
			nonbondedForce->setReactionFieldDielectric(1.0);
			nonbondedForce->setUseDispersionCorrection(false);
		}
	}

	// Add atoms
	for (const auto& atom : atoms) {
		const SimTK::Real charge = atom.physics.chargeInE;
		const SimTK::Real sigma = atom.physics.sigmaInNm;
		const SimTK::Real epsilon = atom.physics.vdwWellDepthInKJ;

		bool finiteCharge = std::isfinite(charge);
		if (!finiteCharge) {
			std::string errorMsg = "Atom " + atom.identity.uniqueAtomName + " has non-finite charge: " + std::to_string(charge);
			SimTK_ASSERT_ALWAYS(finiteCharge, errorMsg.c_str());
		}

		bool finiteSigma = std::isfinite(sigma);
		if (!finiteSigma) {
			std::string errorMsg = "Atom " + atom.identity.uniqueAtomName + " has non-finite sigma: " + std::to_string(sigma);
			SimTK_ASSERT_ALWAYS(finiteSigma, errorMsg.c_str());
		}

		bool finiteEpsilon = std::isfinite(epsilon);
		if (!finiteEpsilon) {
			std::string errorMsg = "Atom " + atom.identity.uniqueAtomName + " has non-finite epsilon: " + std::to_string(epsilon);
			SimTK_ASSERT_ALWAYS(finiteEpsilon, errorMsg.c_str());
		}

		if (atom.connectivity.root) {
			omm.system->addParticle(0.0); // massless root
		} else {
			omm.system->addParticle(atom.physics.massInDaltons);
		}

		if (hasNBfix) {
			nonbondedForce->addParticle(charge, 1.0, 0.0);
			customNonbonded->addParticle({static_cast<SimTK::Real>(atom.identity.nonbondedIndex)});
		} else {
			nonbondedForce->addParticle(charge, sigma, epsilon);
		}

		if (useGBSAOBC2) {
			const SimTK::Real solventRadiusInNm = atom.physics.solventRadiusInNm;
			bool validSolventRadius = std::isfinite(solventRadiusInNm) && solventRadiusInNm > 0.0;
			if (!validSolventRadius) {
				std::string errorMsg = "Atom " + atom.identity.uniqueAtomName + " has invalid solvent radius: " + std::to_string(solventRadiusInNm);
				SimTK_ASSERT_ALWAYS(validSolventRadius, errorMsg.c_str());
			}

			const SimTK::Real screen = atom.physics.screen;
			bool validScreen = std::isfinite(screen) && screen >= 0.0 && screen <= 1.0;
			if (!validScreen) {
				std::string errorMsg = "Atom " + atom.identity.uniqueAtomName + " has invalid screen value: " + std::to_string(screen);
				SimTK_ASSERT_ALWAYS(validScreen, errorMsg.c_str());
			}
			
			GBSAOBCForce->addParticle(charge, solventRadiusInNm, screen);
		}
	}

	if (GBSAOBCForce) {
		omm.registerForce({GBSAOBCForce, ForceGroup::GBSAOBC});
	}

	for (const auto& scaling14 : scaling14s) {
		bool validChargeProduct = std::isfinite(scaling14.chargeProduct);
		if (!validChargeProduct) {
			std::string errorMsg = "Scaling 1-4 between atoms " + std::to_string(scaling14.a1) + " and " + std::to_string(scaling14.a4) + " has non-finite charge product: " + std::to_string(scaling14.chargeProduct);
			SimTK_ASSERT_ALWAYS(validChargeProduct, errorMsg.c_str());
		}

		bool validSigma = std::isfinite(scaling14.sigma);
		if (!validSigma) {
			std::string errorMsg = "Scaling 1-4 between atoms " + std::to_string(scaling14.a1) + " and " + std::to_string(scaling14.a4) + " has non-finite sigma: " + std::to_string(scaling14.sigma);
			SimTK_ASSERT_ALWAYS(validSigma, errorMsg.c_str());
		}

		bool validEpsilon = std::isfinite(scaling14.epsilon);
		if (!validEpsilon) {
			std::string errorMsg = "Scaling 1-4 between atoms " + std::to_string(scaling14.a1) + " and " + std::to_string(scaling14.a4) + " has non-finite epsilon: " + std::to_string(scaling14.epsilon);
			SimTK_ASSERT_ALWAYS(validEpsilon, errorMsg.c_str());
		}

		nonbondedForce->addException(scaling14.a1, scaling14.a4, scaling14.chargeProduct, scaling14.sigma, scaling14.epsilon);
		if (hasNBfix) {
			customNonbonded->addExclusion(scaling14.a1, scaling14.a4);
		}
	}

	for (const auto& exclusion : exclusions) {
		nonbondedForce->addException(exclusion.a1, exclusion.a2, 0.0, 0.1, 0.0);
		if (hasNBfix) {
			customNonbonded->addExclusion(exclusion.a1, exclusion.a2);
		}
	}
	
	// Register non-bonded forces
	omm.registerForce({nonbondedForce, ForceGroup::NonbondedForce});
	if (hasNBfix) {
		omm.registerForce({customNonbonded, ForceGroup::CustomNonbondedForce});
	}
	
	// Add bonds
	// TODO Force constants are expressed for the full quadratic form but OpenMM interprets them as the prefactor of 1/2 k(x-x0)^2, hence the factor of 2 here
	auto* harmonicBondStretch = new OpenMM::HarmonicBondForce();
	for (const auto& bond : bonds) {
		const auto& g = bond.globalIndices;
		harmonicBondStretch->addBond(g[0], g[1], bond.nominalLengthInNm, bond.stiffnessInKJPerNmSq * 2);
	}
	omm.registerForce({harmonicBondStretch, ForceGroup::HarmonicBondForce});

	// Add angles
	// TODO Force constants are expressed for the full quadratic form but OpenMM interprets them as the prefactor of 1/2 k(x-x0)^2, hence the factor of 2 here
	auto* harmonicAngleForce = new OpenMM::HarmonicAngleForce();
	for (const auto& angle : angles) {
		const auto& g = angle.globalIndices;
		harmonicAngleForce->addAngle(g[0], g[1], g[2], angle.nominalAngleInDeg * SimTK::DuMM::Deg2Rad, angle.stiffnessInKJPerRadSq * 2);
	}
	omm.registerForce({harmonicAngleForce, ForceGroup::HarmonicAngleForce});

	// Define torsions
	// Note that OpenMM only supports only periodic torsions by default (AMBER style), so we need to handle improper harmonic torsions differently
	// For AMBER style torsions, OpenMM does not care if they are proper or improper
	auto* periodicTorsionForce = new OpenMM::PeriodicTorsionForce();
	for (const auto& t : properPeriodicTorsions) {
		const auto& g = t.globalIndices;
		for (int i = 0; i < t.numTerms; ++i) {
			const auto& term = t.terms[i];
			periodicTorsionForce->addTorsion(g[0], g[1], g[2], g[3], term.periodicity, term.phaseDeg * SimTK::DuMM::Deg2Rad, term.amplitudeKJ);
		}
	}
	omm.registerForce({periodicTorsionForce, ForceGroup::PeriodicTorsionForce});

	// Periodic torsions are handles with a custom...
	if (!harmonicImproperTorsions.empty()) {
		// "0.5 * k * (theta - theta0)^2"
		std::stringstream ss;
		ss << std::setprecision(17) << "k*min(dtheta, 2*" << SimTK::Pi << "-dtheta)^2; dtheta = abs(theta-theta0)";

		auto* improperTorsionForce = new OpenMM::CustomTorsionForce(ss.str());
		improperTorsionForce->addPerTorsionParameter("k");
		improperTorsionForce->addPerTorsionParameter("theta0");

		for (const auto& t : harmonicImproperTorsions) {
			const int a1 = t.globalIndices[0];
			const int a2 = t.globalIndices[1];
			const int a3 = t.globalIndices[2];
			const int a4 = t.globalIndices[3];

			std::vector<double> params = { t.stiffnessInKJPerRadSq, t.nominalAngleInRad }; 
			improperTorsionForce->addTorsion(a1, a2, a3, a4, params);
		}
		omm.registerForce({improperTorsionForce, ForceGroup::ImproperTorsionForce});
	}

	// Add correction map torsions (CMAPs) force
	if (!cmapGrids.empty()) {
		auto* cmapTorsionForce = new OpenMM::CMAPTorsionForce();
		cmapTorsionForce->setUsesPeriodicBoundaryConditions(false);
		
		for (const auto& grid : cmapGrids) {
			cmapTorsionForce->addMap(grid.size, grid.energy);
		}
		for (const auto& cmapTorsion : cmapTorsions) {
			cmapTorsionForce->addTorsion(
				cmapTorsion.mapIndex,
				cmapTorsion.a1, cmapTorsion.a2, cmapTorsion.a3, cmapTorsion.a4,
				cmapTorsion.b1, cmapTorsion.b2, cmapTorsion.b3, cmapTorsion.b4
			);
		}
		omm.registerForce({cmapTorsionForce, ForceGroup::CMAPTorsion});
	}

	// Add Urey-Bradley Potential
	if (!ureyBradleys.empty()) {
		auto* ubForce = new OpenMM::HarmonicBondForce();
		for (const auto& term : ureyBradleys) {
			ubForce->addBond(term.a1, term.a3, term.nominalLengthInNm, term.stiffnessInKJPerNmSq * 2);
		}
		omm.registerForce({ubForce, ForceGroup::UreyBradley});
	}

	// Get the integrator
	omm.integrator = std::make_unique<OpenMM::VerletIntegrator>(0.0007); // TODO should release?
		
#if OPENMM_PLATFORM_CPU
    OpenMM::Platform* platform = new OpenMM::CpuPlatform();
    OpenMM::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "CPU";

#elif OPENMM_PLATFORM_CUDA
    OpenMM::Platform* platform = new OpenMM::CudaPlatform();
    OpenMM::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "CUDA";

#elif OPENMM_PLATFORM_OPENCL
    OpenMM::Platform* platform = new OpenMM::OpenCLPlatform();
    OpenMM::Platform::registerPlatform(platform);
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
		omm.ommAtomsPositionsCache[atom.identity.globalIndex] = OpenMM::Vec3(atom.position[0], atom.position[1], atom.position[2]);
	}
	omm.context->setPositions(omm.ommAtomsPositionsCache);

	return true;
}

OpenMMEnergyComponents OPENMM::getEnergyComponents() {
	ensureInitialized();
	SimTK_ASSERT_ALWAYS(testing, "OPENMM::getEnergyComponents() can only be called in testing mode.");

	OpenMMEnergyComponents components;

	// Total energy
	// Intentional value capture: relies on C++17 guaranteed copy elision.
	const auto state = context->getState(OpenMM::State::Energy, enforcePeriodicBox);
	components["TotalEnergy"] = state.getPotentialEnergy();

	// Individual forces
	for (const auto& f : forceRegistry) {
		const std::string key = to_string(f.group);
		const auto state = context->getState(OpenMM::State::Energy, enforcePeriodicBox, 1<<static_cast<int>(f.group));
		components[key] = state.getPotentialEnergy();
	}
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
	return potentialEnergy;
}

SimTK::Real OPENMM::getKineticEnergy() const {
	ensureInitialized();
	return kineticEnergy;
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

bool OPENMM::integrateTrajectory(const SimTK::Vector_<SimTK::Vec3>& includedAtomPositionsInG, int steps) {
	ensureInitialized();
	
	// TODO This loop transforms from SimTK::Vec3 to OpenMM::Vec3
    // When you have the time and consider that this is worth the effort, you may force OpenMM to accept SimTK::Vec3 directly
    // 1. Context::setPositions which calls
    // 2. ContextImpl::setPositions which calls (for CUDA):
    // 3. CudaUpdateStateDataKernel::setPositions
	ommAtomsPositionsCache.resize(includedAtomPositionsInG.size());
	for (size_t i = 0, n = includedAtomPositionsInG.size(); i < n; ++i) {
        const auto& c = includedAtomPositionsInG[i];
        ommAtomsPositionsCache[i] = {c[0], c[1], c[2]};
    }

	// Cache old positions in case we need to restore them later
	ommAtomsPositionsCacheOld = ommAtomsPositionsCache;

	// Set the new positions to GPU
	context->setPositions(ommAtomsPositionsCache);

	// Try to integrate
	bool success = true;
	try {
		integrator->step(steps);
	} catch (const std::exception& e) {
		// // Restore old positions in case of integration failure
		// restorePositions();
		success = false;
	}

	// Intentional value capture: relies on C++17 guaranteed copy elision.
	const auto state = context->getState(OpenMM::State::Positions, enforcePeriodicBox);
	const auto& pos = state.getPositions();
	
	// Get new positions
	simbodyAtomsPositionsCache.resize(pos.size());
	for (size_t i = 0, n = pos.size(); i < n; ++i) {
		const auto& c = pos[i];
		simbodyAtomsPositionsCache[i] = {c[0], c[1], c[2]};
	}

	return success;
}

void OPENMM::restorePositions() {
	// ensureInitialized();
	// ommAtomsPositionsCache = ommAtomsPositionsCacheOld;
	// context->setPositions(ommAtomsPositionsCache);
}

void OPENMM::getEnergyAndForces(
	bool positionsAlreadySet,
    const std::vector<NonBondedMapping>& nonBondedMappings,
    const SimTK::Vector_<SimTK::Vec3>& includedAtomStation_G,
    const SimTK::Vector_<SimTK::Vec3>& includedAtomPos_G,
    SimTK::Vector_<SimTK::SpatialVec>& includedBodyForces_G,
    SimTK::Real &energy)
{
	ensureInitialized();

	// Set positions in OpenMM context only when requested
	// This is to prevent setting the positions again after we integrated with OpenMM which would be wasteful
	if (!positionsAlreadySet) {
		for (const auto& mapping : nonBondedMappings) {
			const std::size_t dAIx = mapping.dummAtomIndex;
			const std::size_t iax = mapping.includedAtomIndex;
			const SimTK::Vec3& pos_G = includedAtomPos_G[iax];

			ommAtomsPositionsCache[dAIx] = OpenMM::Vec3(pos_G[0], pos_G[1], pos_G[2]);
		}
		context->setPositions(ommAtomsPositionsCache);
	}

	// Get state with energy and forces
	// Intentional value capture: relies on C++17 guaranteed copy elision.
	const auto state = context->getState(OpenMM::State::Energy | OpenMM::State::Forces, enforcePeriodicBox);

	potentialEnergy = state.getPotentialEnergy();
	kineticEnergy = state.getKineticEnergy();
	const auto& forces = state.getForces();

	// Return only the potential energy for now
	energy += potentialEnergy;

	// Map forces from atoms to bodies
    for (const auto& mapping : nonBondedMappings)
    {
		const std::size_t dAIx = mapping.dummAtomIndex;
		const std::size_t iax = mapping.includedAtomIndex;
		const std::size_t ibx = mapping.bodyIndex;

    	const SimTK::Vec3 simForce(forces[dAIx][0], forces[dAIx][1], forces[dAIx][2]);
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

void OPENMM::registerForce(const ForceRegistration& fr) {
	forceRegistry.push_back(fr);
	system->addForce(fr.force);
	fr.force->setForceGroup(static_cast<int>(fr.group));
}