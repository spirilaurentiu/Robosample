#include "OpenMM.hpp"

#include <chrono>
#include <vector>

void ForceGroup::initialize(
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
) {
	std::cout << "Initializing OpenMM force group " << forceGroupIndex << ":" << std::endl;

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

	// Save force group index
	fg = forceGroupIndex;

	// Instantiate the thermostat with adjusted temperature
	auto* thermostat = new OpenMM::AndersenThermostat(thermostatTemperature, collisionFrequency);
	thermostat->setRandomNumberSeed(seed);
	registerForce(thermostat);
	
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

	// Compute rigid body sizes - only multi-atom bodies cause terms to be skipped
	std::unordered_map<int, int> rigidBodySizes;
	for (int rb : rigidBodies) {
		rigidBodySizes[rb]++;
	}

	auto inSameMultiAtomRigidBody = [&](int a, int b) {
		return rigidBodies[a] == rigidBodies[b] && rigidBodySizes.at(rigidBodies[a]) > 1;
	};

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
		registerForce(GBSAOBCForce);
	}

	int numAddedExceptions = 0;

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

		// // If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely omitted from force and energy calculations
		// if (rigidBodies[scaling14.a1] == rigidBodies[scaling14.a4]) {
		// 	nonbondedForce->addException(scaling14.a1, scaling14.a4, 0.0, 0.1, 0.0);
		// 	++numAddedExceptions;
		// } else {
		// 	nonbondedForce->addException(scaling14.a1, scaling14.a4, scaling14.chargeProduct, scaling14.sigma, scaling14.epsilon);
		// }

		nonbondedForce->addException(scaling14.a1, scaling14.a4, scaling14.chargeProduct, scaling14.sigma, scaling14.epsilon);

		if (hasNBfix) {
			customNonbonded->addExclusion(scaling14.a1, scaling14.a4);
		}
	}

	// Exclude 1-2 and 1-3 interactions
	for (const auto& exclusion : exclusions) {
		// If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely omitted from force and energy calculations
		nonbondedForce->addException(exclusion.a1, exclusion.a2, 0.0, 0.1, 0.0);
		++numAddedExceptions;

		if (hasNBfix) {
			customNonbonded->addExclusion(exclusion.a1, exclusion.a2);
		}
	}









	// std::vector<std::pair<int,int>> intraRigidPairs;

	// // Build adjacency list
	// std::vector<std::vector<int>> adj(atoms.size());
	// for (const auto& b : bonds) {
	// 	int a = b.globalIndices[0];
	// 	int c = b.globalIndices[1];
	// 	adj[a].push_back(c);
	// 	adj[c].push_back(a);
	// }

	// // Group atoms by rigid body
	// std::unordered_map<int, std::vector<int>> rbAtoms;
	// for (int i = 0; i < rigidBodies.size(); ++i)
	// 	rbAtoms[rigidBodies[i]].push_back(i);

	// for (const auto& [rb, atomList] : rbAtoms) {

	// 	for (int a : atomList) {

	// 		std::queue<std::pair<int,int>> q;
	// 		std::unordered_set<int> visited;

	// 		q.push({a,0});
	// 		visited.insert(a);

	// 		while (!q.empty()) {
	// 			auto [v,depth] = q.front();
	// 			q.pop();

	// 			if (depth == 4) continue;

	// 			for (int nb : adj[v]) {
	// 				if (rigidBodies[nb] != rb) continue;
	// 				if (visited.insert(nb).second)
	// 					q.push({nb, depth+1});
	// 			}
	// 		}

	// 		for (int b : atomList) {
	// 			if (b <= a) continue;
	// 			if (!visited.count(b))
	// 				intraRigidPairs.emplace_back(a,b);
	// 		}
	// 	}
	// }

	// for (const auto& pair : intraRigidPairs) {
	// 	// If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely omitted from force and energy calculations
	// 	nonbondedForce->addException(pair.first, pair.second, 0.0, 0.1, 0.0);
	// 	++numAddedExceptions;

	// 	if (hasNBfix) {
	// 		customNonbonded->addExclusion(pair.first, pair.second);
	// 	}
	// }








	std::cout << "\tAdded " << numAddedExceptions << " nonbonded exceptions." << std::endl;
	
	// Register non-bonded forces
	registerForce(nonbondedForce);
	if (hasNBfix) {
		registerForce(customNonbonded);
	}
	
	// Add bonds
	int numAddedBonds = 0;
	auto* harmonicBondStretch = new OpenMM::HarmonicBondForce();

	for (const auto& bond : bonds) {
		const auto& g = bond.globalIndices;

		// Don't add bonds between atoms in the same rigid body
		if (inSameMultiAtomRigidBody(g[0], g[1])) continue;

		harmonicBondStretch->addBond(g[0], g[1], bond.nominalLengthInNm, bond.stiffnessInKJPerNmSq * 2);
		++numAddedBonds;
	}

	registerForce(harmonicBondStretch);
	std::cout << "\tAdded " << numAddedBonds << " bonds." << std::endl;

	// Add angles
	int numAddedAngles = 0;
	auto* harmonicAngleForce = new OpenMM::HarmonicAngleForce();

	for (const auto& angle : angles) {
		const auto& g = angle.globalIndices;

		// Skip angles where all three atoms are in the same rigid body
		if (inSameMultiAtomRigidBody(g[0], g[1]) && inSameMultiAtomRigidBody(g[1], g[2])) continue;

		harmonicAngleForce->addAngle(g[0], g[1], g[2], angle.nominalAngleInDeg * SimTK::DuMM::Deg2Rad, angle.stiffnessInKJPerRadSq * 2);
		++numAddedAngles;
	}
	
	registerForce(harmonicAngleForce);
	std::cout << "\tAdded " << numAddedAngles << " angles." << std::endl;

	// Define torsions
	// Note that OpenMM only supports only periodic torsions by default (AMBER style), so we need to handle improper harmonic torsions differently
	// For AMBER style torsions, OpenMM does not care if they are proper or improper
	int numAddedPeriodicTorsions = 0;
	auto* periodicTorsionForce = new OpenMM::PeriodicTorsionForce();

	for (const auto& t : properPeriodicTorsions) {
		const auto& g = t.globalIndices;
		for (int i = 0; i < t.numTerms; ++i) {
			const auto& term = t.terms[i];

			// Skip propers where the middle two atoms are in the same rigid body
			if (!t.improper && inSameMultiAtomRigidBody(g[1], g[2])) continue;

			// Skip improper torsions where all four atoms are in the same rigid body
			if (t.improper && inSameMultiAtomRigidBody(g[0], g[1]) && inSameMultiAtomRigidBody(g[1], g[2]) && inSameMultiAtomRigidBody(g[2], g[3])) continue;

			periodicTorsionForce->addTorsion(g[0], g[1], g[2], g[3], term.periodicity, term.phaseDeg * SimTK::DuMM::Deg2Rad, term.amplitudeKJ);
			++numAddedPeriodicTorsions;
		}
	}
	std::cout << "\tAdded " << numAddedPeriodicTorsions << " periodic torsions." << std::endl;
	registerForce(periodicTorsionForce);

	// Harmonic improper torsions (CHARMM) are handled using a CustomTorsionForce since OpenMM does not support them natively
	if (!harmonicImproperTorsions.empty()) {
		std::stringstream ss;
		ss << std::setprecision(17) << "k*min(dtheta, 2*" << SimTK::Pi << "-dtheta)^2; dtheta = abs(theta-theta0)";

		auto* improperTorsionForce = new OpenMM::CustomTorsionForce(ss.str());
		improperTorsionForce->addPerTorsionParameter("k");
		improperTorsionForce->addPerTorsionParameter("theta0");

		int numAddedImproperTorsions = 0;

		for (const auto& t : harmonicImproperTorsions) {
			const int a1 = t.globalIndices[0];
			const int a2 = t.globalIndices[1];
			const int a3 = t.globalIndices[2];
			const int a4 = t.globalIndices[3];

			// Improper torsions where all four atoms are in the same rigid body won't change, so we can skip them
			if (inSameMultiAtomRigidBody(a1, a2) && inSameMultiAtomRigidBody(a2, a3) && inSameMultiAtomRigidBody(a3, a4)) continue;


			std::vector<double> params = { t.stiffnessInKJPerRadSq, t.nominalAngleInRad }; 
			improperTorsionForce->addTorsion(a1, a2, a3, a4, params);
			++numAddedImproperTorsions;
		}
		registerForce(improperTorsionForce);
		std::cout << "\tAdded " << numAddedImproperTorsions << " harmonic improper torsions." << std::endl;
	}

	// Add correction map torsions (CMAPs) force
	if (!cmapGrids.empty()) {
		auto* cmapTorsionForce = new OpenMM::CMAPTorsionForce();
		cmapTorsionForce->setUsesPeriodicBoundaryConditions(false);
		
		for (const auto& grid : cmapGrids) {
			cmapTorsionForce->addMap(grid.size, grid.energy);
		}

		int numAddedCmapTorsions = 0;

		for (const auto& cmapTorsion : cmapTorsions) {
			const bool phi_fixed =
				inSameMultiAtomRigidBody(cmapTorsion.a1, cmapTorsion.a2) &&
				inSameMultiAtomRigidBody(cmapTorsion.a2, cmapTorsion.a3) &&
				inSameMultiAtomRigidBody(cmapTorsion.a3, cmapTorsion.a4);
			const bool psi_fixed =
				inSameMultiAtomRigidBody(cmapTorsion.b1, cmapTorsion.b2) &&
				inSameMultiAtomRigidBody(cmapTorsion.b2, cmapTorsion.b3) &&
				inSameMultiAtomRigidBody(cmapTorsion.b3, cmapTorsion.b4);
			if (phi_fixed && psi_fixed) continue;

			cmapTorsionForce->addTorsion(
				cmapTorsion.mapIndex,
				cmapTorsion.a1, cmapTorsion.a2, cmapTorsion.a3, cmapTorsion.a4,
				cmapTorsion.b1, cmapTorsion.b2, cmapTorsion.b3, cmapTorsion.b4
			);
			++numAddedCmapTorsions;
		}

		std::cout << "\tAdded " << numAddedCmapTorsions << " CMAP torsions." << std::endl;
		registerForce(cmapTorsionForce);
	}

	// Add Urey-Bradley Potential
	if (!ureyBradleys.empty()) {
		auto* ubForce = new OpenMM::HarmonicBondForce();
		int numAddedUreyBradleys = 0;

		for (const auto& term : ureyBradleys) {
			if (inSameMultiAtomRigidBody(term.a1, term.a3)) continue;

			ubForce->addBond(term.a1, term.a3, term.nominalLengthInNm, term.stiffnessInKJPerNmSq * 2);
			++numAddedUreyBradleys;
		}

		std::cout << "\tAdded " << numAddedUreyBradleys << " Urey-Bradley terms." << std::endl;
		registerForce(ubForce);
	}
}

bool OPENMM::initialize(
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
        SimTK::Real collisionFrequency)
{
    // Instantiate
    OPENMM& omm = get();

	// Set up variables
	omm.numAtoms = atoms.size();
	omm.ommAtomsPositionsCache = std::vector<OpenMM::Vec3>(omm.numAtoms);
	omm.ommAtomsPositionsCacheOld = std::vector<OpenMM::Vec3>(omm.numAtoms);
	omm.simbodyAtomsPositionsCache = std::vector<SimTK::Vec3>(omm.numAtoms);

	// Allocate OpenMM system and add particles to it
	omm.system = std::make_unique<OpenMM::System>();

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
	}

	// Add a force group for each world
	for (int forceGroupIndex = 0; forceGroupIndex < worlds.size(); ++forceGroupIndex) {
		omm.forceGroups.emplace_back();
		omm.forceGroups.back().initialize(
			forceGroupIndex,
			worlds[forceGroupIndex],
			seed,
			atoms,
			bonds,
			angles,
			properPeriodicTorsions,
			harmonicImproperTorsions,
			cmapGrids,
			cmapTorsions,
			ureyBradleys,
			hasNBfix,
			numTypes,
			acoef,
			bcoef,
			exclusions,
			scaling14s,
			useGBSAOBC2,
			gbsaSolventDielectric,
			gbsaSoluteDielectric,
			nonbondedMethod,
			nonbondedCutoffInNm,
			thermostatTemperature,
			collisionFrequency
		);

		for (const auto& force : omm.forceGroups.back().getForces()) {
			omm.system->addForce(force);
		}

		break;
	}

	// Create the integrator
	omm.integrator = std::make_unique<OpenMM::VerletIntegrator>(0.0007);
		
#if USE_CPU
    OpenMM::Platform* platform = new OpenMM::CpuPlatform();
    OpenMM::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "CPU";

#elif USE_REFERENCE
	OpenMM::Platform* platform = new OpenMM::ReferencePlatform();
	OpenMM::Platform::registerPlatform(platform);
	constexpr auto PLATFORM_NAME = "Reference";

#elif USE_OPENCL
    OpenMM::Platform* platform = new OpenMM::OpenCLPlatform();
    OpenMM::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "OpenCL";

#elif USE_CUDA
    OpenMM::Platform* platform = new OpenMM::CudaPlatform();
    OpenMM::Platform::registerPlatform(platform);
    constexpr auto PLATFORM_NAME = "CUDA";
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

bool OPENMM::integrateTrajectory(const SimTK::Vector_<SimTK::Vec3>& includedAtomPositionsInG, int steps, SimTK::Real timeStepInPicoseconds) {
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

	ommAtomsPositionsCacheOld = ommAtomsPositionsCache;
	context->setPositions(ommAtomsPositionsCache);
	integrator->setStepSize(timeStepInPicoseconds);

	// Try to integrate
	bool success = true;
	try {
		// integrator->setIntegrationForceGroups(1 << activeForceGroupIndex);
		integrator->step(steps);
	} catch (const std::exception& e) {
		// // Restore old positions in case of integration failure
		ommAtomsPositionsCache = ommAtomsPositionsCacheOld;
		context->setPositions(ommAtomsPositionsCache);
		success = false;
	}

	// Intentional value capture: relies on C++17 guaranteed copy elision.
	// , 1 << activeForceGroupIndex
	const auto state = context->getState(OpenMM::State::Positions | OpenMM::State::Energy | OpenMM::State::Velocities, enforcePeriodicBox);
	const auto& pos = state.getPositions();
	const auto& velocities = state.getVelocities();

	// DO NOT UNCOMMENT THIS - THIS IS A BUG AND BREAK SOMETHING SOMEHWERE Store energies
	potentialEnergy = state.getPotentialEnergy();
	kineticEnergy = state.getKineticEnergy();

	// std::cout << "\tIntegrated " << steps << " steps " << " at " << timeStepInPicoseconds << " ps with force group " << activeForceGroupIndex << ". Potential energy: " << state.getPotentialEnergy() << " kJ/mol, Kinetic energy: " << state.getKineticEnergy() << " kJ/mol" << std::endl;
	
	// Get new positions
	simbodyAtomsPositionsCache.resize(pos.size());
	for (size_t i = 0, n = pos.size(); i < n; ++i) {
		const auto& c = pos[i];
		simbodyAtomsPositionsCache[i] = {c[0], c[1], c[2]};
	}

	return success;
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
	// Intentional value capture: relies on C++17 guaranteed copy elision
	// OpenMM already evaluates the potential energy internally in most kernels (the energy reduction is cheap once forces are computed)
	// In that case retrieving the energy is essentially free.
	// , 1 << activeForceGroupIndex
	const auto state = context->getState(OpenMM::State::Energy | OpenMM::State::Forces, enforcePeriodicBox);

	// auto start = std::chrono::high_resolution_clock::now();

	// const auto state = context->getState(
	// 	OpenMM::State::Energy | OpenMM::State::Forces,
	// 	enforcePeriodicBox,
	// 	1 << activeForceGroupIndex
	// );

	// auto end = std::chrono::high_resolution_clock::now();
	// auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	// std::cout << "\tgetState(" << activeForceGroupIndex << ") took " << duration << " us" << std::endl;

	potentialEnergy = state.getPotentialEnergy();
	kineticEnergy = state.getKineticEnergy();
	const auto& forces = state.getForces();

	// std::cout << "\tForce group " << activeForceGroupIndex << " energy: " << potentialEnergy << " kJ/mol" << std::endl;

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
