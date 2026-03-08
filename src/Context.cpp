#include "Context.hpp"
#include "Compound.h"
#include "OpenMM.hpp"
#include "Replica.hpp"
#include "World.hpp"
#include "Sampler.hpp"
#include "common.h"
#include "readAmberInput.hpp"

#include <cstddef>
#include <iostream>
#include <sys/stat.h>
#include <sys/sysinfo.h>

#include <stack>

/*!
 * <!-- Constructor: sets temperatures, random engine and checks for CUDA_ROOT -->
*/
Context::Context(const std::string& baseName_arg, uint32_t seed, uint32_t threads, uint32_t nofRoundsTillReblock, RUN_TYPE runType, uint32_t swapFreq, uint32_t swapFixmanFreq, bool testing)
{
	// Set the base name of the simulation
	std::cout << "Context with base name: " << baseName + "_" + std::to_string(seed) << std::endl << std::flush;
	this->baseName = baseName_arg + "_" + std::to_string(seed);

	// Alert user of CUDA environment variables
#if USE_CUDA
	if (SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT").empty()){
		std::cerr << cwar_prefix << "CUDA_ROOT not set." << std::endl;
	} else {
		std::cout << cinf_prefix << "CUDA_ROOT set to " << SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT") << std::endl;
	}
#endif

	// Use a random seed if none is provided
	if (seed == 0) {
		std::random_device rd;
		this->seed = rd();
	}else{
		this->seed = seed;
	}

	// Set the random seed
	randomEngine = buildRandom32(seed);

	// Set the number of threads
	if (threads < 0) {
		std::cerr << "Invalid number of threads (negative value). Default number (0) of threads will be used." << std::endl;
		numThreads = 0;
	} else {
		numThreads = threads;
	}

	this->roundsTillReblock = nofRoundsTillReblock;
	this->runType = runType;
	this->swapEvery = swapFreq;
	this->swapFixman = swapFixmanFreq;

	// foutU = std::string(baseName + "_U.bin");
	// foutUDot = std::string(baseName + "_U_dot.bin");
	// foutTorque = std::string(baseName + "_torque.bin");

	// Run in testing mode
	this->testing = testing;
}

void Context::setVerbose(bool verbose){
	this->verbose = verbose;
}

/*!
 * <!--  -->
*/
bool Context::setOutput(const std::string& outDir) {

	// Set the log filename
	std::string logFilename = outDir + "/log." + std::to_string(seed);

	// Open the log file
	logFile = std::ofstream(logFilename);
	if ( !logFile.is_open() ) {
		std::cerr << cerr_prefix << "Failed to open log file " << logFilename << "." << std::endl;
		return false;
	}

	// Set the directory where the logs and the trajectories are stored
	if( !SimTK::Pathname::fileExists(outDir + "/pdbs") ){
		const int err = mkdir((outDir + "/pdbs").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (err == -1){
			std::cerr << cerr_prefix << "Failed to create " << outDir + "/pdbs" << "." << std::endl;
			return false;
		}
	}

	setOutputDir(outDir);

	return true;
}

void Context::loadAmberSystem(
	const std::vector<int>& roots_,
	const std::vector<RoboAtom>& atoms_,
	const std::vector<RoboBond>& bonds_,
	const std::vector<RoboAngle>& angles_,
	const std::vector<RoboPeriodicTorsion>& properPeriodicTorsions_,
	const std::vector<RoboHarmonicImproperTorsion>& harmonicImproperTorsions_,
	const std::vector<TopologyRange>& topologyRanges)
{
	// Copy the data
	roots = roots_;
	atoms = atoms_;
	bonds = bonds_;
	angles = angles_;
	properPeriodicTorsions = properPeriodicTorsions_;
	harmonicImproperTorsions = harmonicImproperTorsions_;\

	numMolecules = roots.size();

	// Construct a Compound for every atom
	// Since we iterate the list of atoms, we can also validate the global indices
	// We expect them to be contiguous from 0 to N-1
	int expectedGlobalIndex = 0;
	for(auto& atom : atoms) {

		// Create the SimTK::Compound
		atom.createSingleAtom();

		// Validate the continuity of the rest of the sequence
		if (atom.identity.globalIndex  != expectedGlobalIndex) {
			std::string errorMsg = "Sequence Error: Global indices are not contiguous. "
				+ std::to_string(expectedGlobalIndex) + " expected, but found "
				+ std::to_string(atom.identity.globalIndex) + " for atom `"
				+ atom.identity.uniqueAtomName + "`.";

			SimTK_ASSERT_ALWAYS(false, errorMsg.c_str());
		}

		expectedGlobalIndex++;
	}

	// Add new topologies
	topologies.reserve(numMolecules);
	for(std::size_t molIx = 0; molIx < roots.size(); molIx++) {

		// New empty topology
		SimTK::Compound::Name name = "MOL_" + std::to_string(molIx);
		Topology topology(name, SimTK::CompoundSystem::CompoundIndex(molIx), roots[molIx]);

		// Set spans
		const auto atomRangeBegin = topologyRanges[molIx].getRange(TopologyRangeType::ATOM).first;
		const auto atomRangeEnd = topologyRanges[molIx].getRange(TopologyRangeType::ATOM).second;
		topology.setAtoms(safe_subspan(atoms, atomRangeBegin, atomRangeEnd));

		const auto bondRangeBegin = topologyRanges[molIx].getRange(TopologyRangeType::BOND).first;
		const auto bondRangeEnd = topologyRanges[molIx].getRange(TopologyRangeType::BOND).second;
		topology.setBonds(safe_subspan(bonds, bondRangeBegin, bondRangeEnd));

		const auto angleRangeBegin = topologyRanges[molIx].getRange(TopologyRangeType::ANGLE).first;
		const auto angleRangeEnd = topologyRanges[molIx].getRange(TopologyRangeType::ANGLE).second;
		topology.setAngles(safe_subspan(angles, angleRangeBegin, angleRangeEnd));

		const auto properPeriodicTorsionRangeBegin = topologyRanges[molIx].getRange(TopologyRangeType::PERIODIC_TORSION).first;
		const auto properPeriodicTorsionRangeEnd = topologyRanges[molIx].getRange(TopologyRangeType::PERIODIC_TORSION).second;
		topology.setPeriodicTorsions(safe_subspan(properPeriodicTorsions, properPeriodicTorsionRangeBegin, properPeriodicTorsionRangeEnd));

		const auto improperHarmonicTorsionRangeBegin = topologyRanges[molIx].getRange(TopologyRangeType::IMPROPER_HARMONIC_TORSION).first;
		const auto improperHarmonicTorsionRangeEnd = topologyRanges[molIx].getRange(TopologyRangeType::IMPROPER_HARMONIC_TORSION).second;
		topology.setImproperHarmonicTorsions(safe_subspan(harmonicImproperTorsions, improperHarmonicTorsionRangeBegin, improperHarmonicTorsionRangeEnd));

		// Set root atom
		// Its compound atom index is set inside the next loop
		RoboAtom& rootAtom = atoms[roots[molIx]];
		topology.setBaseAtom(*rootAtom.compoundSingleAtom, SimTK::Transform());
		topology.convertInboardBondCenterToOutboard();

		// std::cout << cinf_prefix << "Set root atom " << rootAtom.identity.uniqueAtomName
		// 		  << " for molecule " << molIx
		// 		  << std::endl << std::flush;

		// Add non-ring closing bonds first
		for(auto& bond : topology.updBonds()) {

			if (bond.ringClosing) continue;

			RoboAtom& parent = atoms[bond.globalIndices[0]];
			RoboAtom& child = atoms[bond.globalIndices[1]];

			// Get next available bond center ID for the parent
			const int parentNofBonds = parent.connectivity.numBondsInvolved;
			const int parentNofFreebonds = parent.connectivity.numAvailableBonds;
			const int parentNextAvailBondCenter = parentNofBonds - parentNofFreebonds + 1;

			// Cook the parentBondCenterPathName = RESNAME + RESID + _ATOMNAME + bond int(next)
			const SimTK::Compound::BondCenterPathName parentBondCenterPathName = parent.identity.uniqueAtomName + "/bond" + std::to_string(parentNextAvailBondCenter);

			// std::cout << cinf_prefix << "Bonding child atom " << child.identity.uniqueAtomName << " cAIx " << child.identity.compoundAtomIndex
			// 		  << " to parent bond center " << parentBondCenterPathName << " cAIx " << parent.identity.compoundAtomIndex
			// 		  << " for molecule " << molIx
			// 		  << std::endl << std::flush;

			// Actual bonding with default mobility (torsion)
			topology.bondAtom(*child.compoundSingleAtom, parentBondCenterPathName, bond.nominalLengthInNm, 109.47 * SimTK::Deg2Rad, SimTK::BondMobility::Default);

			// Set the local compound atom index for this child
			// 1 is child, 0 is for parent
			// SimTK::Compound::AtomIndex childCAIx = topology.getBondAtomIndex(SimTK::Compound::BondIndex(topology.getNumBonds() - 1), 1);
			// child.setCompoundAtomIndex(childCAIx);
			topology.setAtomMass(child.identity.compoundAtomIndex, child.physics.massInDaltons);

			// Set the local compound atom index for the parent if it is the root
			if(bond.globalIndices[0] == roots[molIx]) {
				// parent.setCompoundAtomIndex(parentCAIx);
				topology.setAtomMass(parent.identity.compoundAtomIndex, parent.physics.massInDaltons);
			}

			// // Handle ions
			// if (internCoords.getRoot(molIx).second == -1) {
			// 	atoms[internCoords.getRoot(molIx).first].setMoleculeIndex(molIx);
			// 	atoms[internCoords.getRoot(molIx).first].setCompoundAtomIndex(SimTK::Compound::AtomIndex(0));
			// }	

			parent.connectivity.numAvailableBonds--;
			child.connectivity.numAvailableBonds--;
		}

		// Add ring closing bonds
		for(const auto& bond : topology.updBonds()) {

			if (!bond.ringClosing) continue;

			RoboAtom& parent = atoms[bond.globalIndices[0]];
			RoboAtom& child = atoms[bond.globalIndices[1]];

			// Molmodel expects bond center names to not have been taken yet
			const int childNofBonds = child.connectivity.numBondsInvolved;
			const int childNofFreebonds = child.connectivity.numAvailableBonds;
			const int childNextAvailBondCenter = childNofBonds - childNofFreebonds + 1;
			const SimTK::Compound::BondCenterPathName bondCenterName1 = child.identity.uniqueAtomName + "/bond" + std::to_string(childNextAvailBondCenter);

			// Molmodel expects bond center names to not have been taken yet
			const int parentNofBonds = parent.connectivity.numBondsInvolved;
			const int parentNofFreebonds = parent.connectivity.numAvailableBonds;
			const int parentNextAvailBondCenter = parentNofBonds - parentNofFreebonds + 1;
			const SimTK::Compound::BondCenterPathName bondCenterName2 = parent.identity.uniqueAtomName + "/bond" + std::to_string(parentNextAvailBondCenter);

			// std::cout << cinf_prefix << "Adding ring closing bond between bond centers "
			// 		  << bondCenterName1 << " cAIx " << child.identity.compoundAtomIndex
			// 		  << " and " << bondCenterName2 << " cAIx " << parent.identity.compoundAtomIndex
			// 		  << " for molecule " << molIx
			// 		  << std::endl << std::flush;

			topology.addRingClosingBond(bondCenterName1, bondCenterName2, bond.nominalLengthInNm, 109.47 * SimTK::Deg2Rad, SimTK::BondMobility::Rigid);

			parent.connectivity.numAvailableBonds--;
			child.connectivity.numAvailableBonds--;
		}

		// Make sure all bonds have been used
		// Usually, it's from missing ring closing bonds
		std::vector<std::size_t> unsatisfiedAtomsIndices;
		for (auto& atom : topology.getAtoms()) {
			if (atom.connectivity.numAvailableBonds != 0) {
				unsatisfiedAtomsIndices.push_back(atom.identity.globalIndex);
			}
		}

		if (!unsatisfiedAtomsIndices.empty()) {
			std::string error_msg = "Not all bonds have been satisfied when building topology for molecule " + std::to_string(molIx) + ". ";
			for (const auto& atomIx : unsatisfiedAtomsIndices) {
				const RoboAtom& atom = atoms[atomIx];
				error_msg += "\tAtom " + atom.identity.uniqueAtomName + " has " + std::to_string(atom.connectivity.numAvailableBonds) + " unsatisfied bonds.";
			}

			SimTK_ASSERT_ALWAYS(unsatisfiedAtomsIndices.empty(), error_msg.c_str());
		}

		// Define the biotype of the atom
		for (auto& atom: topology.getAtoms()) {
			// It calls SimTK::Biotype::defineBiotype and checks if it already exists
			// This expects three parameters: the unique name of the atom (combination of residue name, residue id, force field atom type name and global id), the residue name and the force field atom type name
			// However, we skip the residue name e.g N from ALA_1 and N from ALA_2 are different atoms because of the valence (N terminus vs backbone valence)
			// Thus, we rely solely on the unique atom name and a custom atom type name (see exactly in the Python wrapper how this is set)
			topology.setAtomBiotype(atom.identity.uniqueAtomName.c_str(), "", atom.identity.chargedAtomTypeName.c_str());
			atom.biotypeIndex = topology.getAtomBiotypeIndex(atom.identity.compoundAtomIndex);
		}

		// Get coordinates
		SimTK::Compound::AtomTargetLocations atomTargets;
		for (const auto& a : topology.getAtoms()) {
			atomTargets.insert(std::make_pair(a.identity.compoundAtomIndex, a.position));
		}
		atomTargetLocationsCache.emplace_back(atomTargets);

		topology.setTopLevelTransform(SimTK::Transform(SimTK::Rotation(), rootAtom.position));

		// Match the topology to the input coordinates
		topology.loadIndicesMaps();
		topology.matchAtomTargetLocations(atomTargets);

		topologies.push_back(topology);
	}
}

OpenMMEnergyComponents Context::initializeOpenMM(
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
    const std::vector<Scaling14>& scaling14s)
{
	OPENMM::initialize(
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
		solventDielectric,
		soluteDielectric,
		nonbondedMethod,
		nonbondedCutoffInNm,
		300.0, // thermostatTemperature
		1.0, // collisionFrequency
		testing
	);
	
	return OPENMM::get().getEnergyComponents();
}

/*! <!--  --> */
void Context::Initialize() {

	// // Initialize the Z matrix
	// int firstWIx = 0;
	// SimTK::State& lastAdvancedState = worlds[firstWIx].updIntegrator().updAdvancedState();

	// // Get coordinates from source
	// const auto& firstWorldsAtomsLocations = worlds[firstWIx].getAtomsLocationsInGround(lastAdvancedState);

	// // Get Z-matrix indexes table	
	// calcZMatrixTable();
	// PrintZMatrixTable();
	// reallocZMatrixBAT();
	// calcZMatrixBAT(firstWIx, firstWorldsAtomsLocations);
	// PrintZMatrixBAT();
	// PrintZMatrixMobods(firstWIx, lastAdvancedState);
	
	// for(int k = 0; k < nofReplicas; k++){
	// 	replicas[k].reallocZMatrixBAT();
	// 	replicas[k].calcZMatrixBAT(firstWorldsAtomsLocations);
	// 	//replicas[k].PrintZMatrixBAT();
	// }
}

/*!
 * <!-- Add flexibilities and CompoundSystem model -->
*/
void Context::addWorld(
	bool fixmanTorque,
	int samplesPerRound,
	ROOT_MOBILITY rootMobility,
	const std::vector<std::vector<BondFlexibility>>& rollFlexibilities,
	bool useOpenMM,
	bool visual,
	SimTK::Real visualizerFrequency)
{
	// Create new world and add its index
	worldIndexes.push_back(worldIndexes.size());
	Span<Topology> t {topologies};
	worlds.emplace_back(worldIndexes.back(), t, testing, visual, visualizerFrequency);

	// Set force field scale factor.
	if (useAmberForceFieldScaleFactors) {
		worlds.back().setAmberForceFieldScaleFactors();
	} else {
		worlds.back().setGlobalForceFieldScaleFactor(globalForceFieldScaleFactor);
	}

	// Set the nonbonded method and cutoff
	worlds.back().updForceField().setNonbondedMethod(nonbondedMethod);
	worlds.back().updForceField().setNonbondedCutoff(nonbondedCutoffInNm);

	// If requested, add Fixman torque as an additional force subsystem
	if (fixmanTorque) {
		worlds.back().addFixmanTorque();
		worlds.back().updFixmanTorque()->setScaleFactor(1);
	}

	// Set the number of threads for DuMM
	// TODO
	if (numThreads == 1) {
		worlds.back().updForceField().setUseMultithreadedComputation(false);
	} else {
		worlds.back().updForceField().setNumThreadsRequested(numThreads);
	}

	// Set GBSA scaling and VdW mixing rule
	worlds.back().setGbsaGlobalScaleFactor(gbsaGlobalScaleFactor);
	worlds.back().updForceField().setVdwMixingRule(SimTK::DuMMForceFieldSubsystem::LorentzBerthelot); // DuMMForceFieldSubsystem::WaldmanHagler

	// Set how many times to run sample_iteration()
	worlds.back().setSamplesPerRound(samplesPerRound);

	// Set temperatures for sampler and Fixman torque is applied to this world
	worlds.back().setTemperature(tempIni);

	// Set seed for random number generators
	worlds.back().setSeed(randomEngine());

	// @TODO what does it do? how does it work if we have multiple molecules?
	// Propagate root mobility
	worlds.back().setRootMobility(rootMobility);

	// Store the number of worlds
	nofWorlds = worlds.size();

	// Prepare the world for OpenMM
	if (useOpenMM) {
		worlds.back().updForceField().setUseOpenMMAcceleration(true);
	}

	// Generate DuMM parameters: DuMM atom types, charged atom types, bond types, angle types and torsion types
	worlds.back().generateDummParams(atoms, bonds, angles, properPeriodicTorsions, harmonicImproperTorsions);

	// Allocate root mobilities - TODO
	rootMobilitiesStr.push_back({});
	for(unsigned int molIx = 0; molIx < topologies.size(); molIx++){
		rootMobilitiesStr.back().push_back("Rigid");
	}

	// Helper to create a normalized bond key (always min-first) e.g., (3,5) and (5,3) both become (3,5)
	auto make_bond_key = [](int i, int j) {
		return std::make_pair(std::min(i, j), std::max(i, j));
	};

	// This is perfectly collision-safe because std::map uses lexicographical comparison, not hashing like std::unordered_map 
	// Key: {atom_i, atom_j}, Value: {mobility, was_satisfied}
    struct FlexStatus { SimTK::BondMobility::Mobility mobility; bool satisfied{false}; };

    std::map<std::pair<int, int>, FlexStatus> flex_map;
	for (const auto& flex_list : rollFlexibilities) {
        for (const auto& flex : flex_list) {
            flex_map[make_bond_key(flex.globalIndex1, flex.globalIndex2)] = {flex.mobility, false};
        }
    }

	// Apply user bond flexibilities
	for (auto& topology : topologies) {
        for (auto& bond : topology.getBonds()) {
            const int p = bond.globalIndices[0];
            const int c = bond.globalIndices[1];
            const auto key = make_bond_key(p, c);

            SimTK::BondMobility::Mobility mobility = SimTK::BondMobility::Mobility::Rigid;

			// Check if user specified flexibility for this bond
            auto it = flex_map.find(key);
            if (it != flex_map.end()) {
                it->second.satisfied = true;

				if (bond.ringClosing && it->second.mobility != SimTK::BondMobility::Rigid) {
					std::cout << "\tWARNING: Custom bond mobility (" 
						<< SimTK::BondMobility::getBondMobilityName(it->second.mobility)
						<< ") cannot be applied to the ring-closing bond between atoms " 
						<< atoms[p].identity.uniqueAtomName << " (index " << p << ") and "
						<< atoms[c].identity.uniqueAtomName << " (index " << c << "). "
						<< "Ring-closing bonds must be 'Rigid'; defaulting to Rigid mobility." 
						<< std::endl;
				} else {
					mobility = it->second.mobility;
				}
            }

			// This add a new bond mobibility to the bonds lists
			// Each bond has precisely one bond mobility per world
            bond.addBondMobility(mobility);

			// However, this is a trick
			// A robot (molecule) is not stored in a Compound, but in CompoundSystem which builds it from a topology (which inherits from Compound)
			// We call setBondMobility (which is inherited from Compound) to temporarily set the bond mobility in the topology
			// Later, when we call CompoundSystem::modelOneCompound, the bond mobilities are read from the topology and used to build the robot accordingly
			const SimTK::Compound::AtomName parentAtomName = atoms[p].identity.uniqueAtomName;
			const SimTK::Compound::AtomName childAtomName = atoms[c].identity.uniqueAtomName;
			topology.setBondMobility(mobility, parentAtomName, childAtomName);
            
			// Print status of bond flexibility setting
			bool verbose_local = false;
			if (verbose_local) {
				const RoboAtom& parentAtom = atoms[p];
				const RoboAtom& childAtom = atoms[c];
				std::cout << "Setting bond flexibility for atoms "
						  << parentAtom.identity.uniqueAtomName << " (BAT index " << p << ") and "
						  << childAtom.identity.uniqueAtomName << " (BAT index " << c << ") to "
						  << SimTK::BondMobility::getBondMobilityName(mobility)
						  << std::endl;
			}
        }
    }

    // Input sanitization: check that all user-specified flexibilities were satisfied
    for (auto const& entry : flex_map) {
        const FlexStatus& status = entry.second;
        
        if (!status.satisfied) {
            const std::pair<int, int>& indices = entry.first;
            std::string error_msg = "Error: User-specified bond flexibility for atoms " + 
                                	 atoms[indices.first].identity.uniqueAtomName + " (BAT index " + std::to_string(indices.first) + ") and " +
									 atoms[indices.second].identity.uniqueAtomName + " (BAT index " + std::to_string(indices.second) + ") " +
									 "was not found in the system.";
            
            // Condition set to 'false' so the assert actually triggers when it hits this block
            SimTK_ASSERT_ALWAYS(false, error_msg.c_str());
        }
    }

	// Let DuMM model this robot
	worlds.back().modelTopologies();
	worlds.back().setAtomTargetLocationsToState(atomTargetLocationsCache);

	// Check that this world is not overconstrained
	const auto& state = worlds.back().getIntegrator().getAdvancedState();
	worlds.back().getCompoundSystem().realize(state, SimTK::Stage::Acceleration);

	if (worlds.back().isOverconstrained(state, std::cout)) {
		const std::string msg = "World " + std::to_string(worldIndexes.back()) + " is overconstrained.";
		throw std::runtime_error(msg);
	}

	// Find bodies for roll mobilities
	std::vector<std::vector<SimTK::MobilizedBodyIndex>> mobodLocks;

	for (const auto& flex_list : rollFlexibilities) {
		mobodLocks.emplace_back();

		for (const auto& flex : flex_list) {
			const auto aIx1 = atoms[flex.globalIndex1].identity.compoundAtomIndex;
			const auto aIx2 = atoms[flex.globalIndex2].identity.compoundAtomIndex;

			// Sanity check
			const auto topoIx1 = atoms[flex.globalIndex1].identity.moleculeIndex;
			const auto topoIx2 = atoms[flex.globalIndex2].identity.moleculeIndex;
			if (topoIx1 != topoIx2) {
				std::string error_msg = "Error: Atoms " + atoms[flex.globalIndex1].identity.uniqueAtomName + " (BAT index " + std::to_string(flex.globalIndex1) + ") and " +
									   atoms[flex.globalIndex2].identity.uniqueAtomName + " (BAT index " + std::to_string(flex.globalIndex2) + ") " +
									   "are in different molecules, but a bond flexibility was specified between them.";
				throw std::runtime_error(error_msg);
			}

			// Get mobilized bodies
			const SimTK::MobilizedBodyIndex mbx1 = topologies[topoIx1].getAtomMobilizedBodyIndex(aIx1);
			const SimTK::MobilizedBodyIndex mbx2 = topologies[topoIx2].getAtomMobilizedBodyIndex(aIx2);
			if (mbx1 == mbx2) {
				std::string error_msg = "Error: Atoms " + atoms[flex.globalIndex1].identity.uniqueAtomName + " (BAT index " + std::to_string(flex.globalIndex1) + ") and " +
									   atoms[flex.globalIndex2].identity.uniqueAtomName + " (BAT index " + std::to_string(flex.globalIndex2) + ") " +
									   "are in the same mobilized body, but a bond flexibility was specified between them.";
				throw std::runtime_error(error_msg);
			}

			// Assign mobilized bodies
			const SimTK::MobilizedBodyIndex mbx2Parent = worlds.back().getMatterSubsystem().getMobilizedBody(mbx2).getParentMobilizedBody();
			if (mbx2Parent == mbx1) {
				mobodLocks.back().push_back(mbx2);
			} else {
				mobodLocks.back().push_back(mbx1);
			}
		}
	}

	worlds.back().setMobodLocks(mobodLocks);
}

/** Add task spaces */
void Context::addTaskSpacesLS(void)
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		worlds[worldIx].addTaskSpaceLS();
	}
}

/** Add rod constraints */
void Context::addConstraints(void)
{
		// for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		// 	worlds[worldIx].addRodConstraint(worlds[worldIx].updIntegrator().updAdvancedState());
		// }
}

// Print status
void Context::printStatus(void){
	// for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
	// 	SimTK::VerletIntegrator& checkIntegrator = *worlds[worldIx].integ;
	// 	const SimTK::State& checkState = checkIntegrator.getState();
	// 	const SimTK::Stage& checkStage = checkState.getSystemStage();
	// 	std::cout << "Context world " << worldIx << " integ state stage "
	// 		<< checkStage << std::endl << std::flush;
	// 	std::cout << "Context world " << worldIx << " integ state nof Subsystems "
	// 		<< checkState.getNumSubsystems() << ":" << std::endl << std::flush;
	// 	for(int i = 0; i < checkState.getNumSubsystems(); i++){
	// 		std::cout
	// 			<< " Subsystem Name: "
	// 			<< checkState.getSubsystemName(SimTK::SubsystemIndex(i))
	// 			<< " Stage: "
	// 			<< checkState.getSubsystemStage(SimTK::SubsystemIndex(i))
	// 			<< " Version: "
	// 			<< checkState.getSubsystemVersion(SimTK::SubsystemIndex(i))
	// 			<< std::endl << std::flush;
	// 	}
	// 	//SimTK::State& checkAdvState = checkIntegrator.updAdvancedState();
	// 	//const SimTK::Stage& checkAdvStage = checkAdvState.getSystemStage();
	// 	//std::cout << "Context world " << worldIx << " integ advState stage "
	// 	//	<< checkAdvStage << std::endl << std::flush;


	// 	// CompoundSystem <- MolecularMechanicsSystem <- MultibodySystem <- System
	// 	SimTK::CompoundSystem& compoundSystem = *(worlds[worldIx].getCompoundSystem());
	// 	std::cout << "Context world " << worldIx << " compoundSystem nof compounds "
	// 		<< compoundSystem.getNumCompounds() << std::endl;
	// 	std::cout << "Context world " << worldIx << " System Topology realized "
	// 		<< compoundSystem.getNumRealizationsOfThisStage(SimTK::Stage::Topology)
	// 		<< " times.\n" << std::flush;

	// 	// Matter
	// 	////const SimTK::System& checkSystem = (worlds[worldIx].matter)->getSystem();
	// 	SimTK::SimbodyMatterSubsystem& matter = *(worlds[worldIx].matter);
	// 	std::cout << "Context world " << worldIx
	// 		<< " matter nofBodies " << matter.getNumBodies()
	// 		<< " nofConstraints " << matter.getNumConstraints()
	// 		<< "\n" << std::flush;

	// 	// GeneralForceSubsystem
	// 	SimTK::GeneralForceSubsystem& gfs = *(worlds[worldIx].forces);
	// 	std::cout << "Context world " << worldIx
	// 		<< " gfs nofForces " << gfs.getNumForces()
	// 		<< "\n" << std::flush;

	// 	SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[worldIx].forceField);
	// 	std::cout << "Context world " << worldIx
	// 		<< " dumm nofThreads " << dumm.getNumThreadsRequested()
	// 		<< " useOpenMM " << dumm.getUseOpenMMAcceleration()
	// 		<< " " << dumm.isUsingOpenMM()
	// 		<< "\n" << std::flush;


	// }
}

// Print thermodynamics
void Context::printThermodynamics()
{
	// for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
	// 	std::cout << "World " << worldIx << " temperature = "
	// 		<< worlds[worldIx].getTemperature()
	// 		<< std::endl;
	// 	if(worlds[worldIx].isUsingFixmanTorque()){
	// 		std::cout << "World " << worldIx
	// 		<< " FixmanTorque temperature = "
	// 		<< worlds[worldIx].updFixmanTorque()->getTemperature()
	// 		<< std::endl;
	// 	}
	// 	for (int samplerIx = 0; samplerIx < worlds[worldIx].getNofSamplers(); samplerIx++){
	// 		std::cout << "World " << worldIx << " Sampler " << samplerIx
	// 			<< " temperature = " << worlds[worldIx].updSampler(samplerIx)->getTemperature()
	// 			<< " initial const state PE: " << std::setprecision(20)
	// 			//<< worlds[worldIx].forces->getMultibodySystem().calcPotentialEnergy(worlds[worldIx].updIntegrator().updAdvancedState())
	// 			//<< worlds[worldIx].forces->getMultibodySystem().calcPotentialEnergy(updAdvancedState(worldIx, samplerIx))
	// 			<< " useFixmanPotential = "
	// 			<< pHMC(worlds[worldIx].updSampler(samplerIx))->isUsingFixmanPotential()
	// 			<< std::endl;
	// 	}

	// }
}

// Print Simbody related information
void Context::PrintSimbodyMobods(){
	// for(std::size_t worldIx = 0; worldIx < nofWorlds; worldIx++){
	// 	std::cout << "Context::PrintSimbodyMobods world " << worldIx << "\n";
	// 	for(std::size_t molIx = 0; molIx < numMolecules; molIx++){
	// 		std::cout << "Context::PrintSimbodyMobods molecule " << molIx << "\n";
	// 		const Topology& topology = worlds[worldIx].getTopology(molIx);

	// 		for(std::size_t i = 0; i < topology.getNumAtoms(); i++){
	// 			SimTK::Compound::AtomIndex aIx
	// 				= (topology.subAtomList[i]).identity.compoundAtomIndex;
	// 			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(aIx);
	// 			std::cout << "i aIx mbx " << i << " " << aIx << " "
	// 				<< mbx << std::endl << std::flush;
	// 		}
	// 	}
	// }
}

// Print DuMM atoms stations in mobilized body frame
void Context::checkAtomStationsThroughDumm()
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for (int samplerIx = 0;
			samplerIx < worlds[worldIx].getNofSamplers();
			samplerIx++){
			(worlds[worldIx].updSampler(samplerIx))->checkAtomStationsThroughDumm();
		}
	}
}

// Another way to do it is setting the number of rounds
int Context::getRequiredNofRounds()
{
	return requiredNofRounds;
}

void Context::setRequiredNofRounds(int argNofRounds)
{
	requiredNofRounds = argNofRounds;
}

int Context::getNofRoundsTillReblock()
{
	return roundsTillReblock;
}

void Context::setNofRoundsTillReblock(int nofRoundsTillReblock)
{
	this->roundsTillReblock = nofRoundsTillReblock;
}

void Context::updNofRoundsTillReblock(int nofRoundsTillReblock)
{
	this->roundsTillReblock = nofRoundsTillReblock;
}

// Adaptive Gibbs blocking: TODO: consider moving in World
void Context::allocateReblockQsCache()
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		QsCache.push_back(std::vector<std::vector<SimTK::Real>>(roundsTillReblock));
		//std::cout << "Context::AddWorld QsCache size " << QsCache.size() << std::endl;
	}
}

// TODO This seems wrong !!!
void Context::allocateReblockQsCacheQVectors(){
	SimTK_ASSERT_ALWAYS(false, "Context::allocateReblockQsCacheQVectors: Not implemented yet.");
	
	// // Adaptive Gibbs blocking: // TODO generalized coord may not always be Real
	// if(QsCache[0][0].size() == 0){
	// 	std::size_t worldIx = 0;
	// 	for(auto& world : worlds) {
	// 		int nQs = world.getCompoundSystem().getMatterSubsystem().getSystem().getDefaultState().getNQ();
	// 		//std::cout << "World " << worldIx  << " has " << nQs << " Qs" << std::endl;
	// 		//std::cout << "Context::realizeTopology QsCache[" << worldIx << "] size " << QsCache[worldIx].size() << std::endl;

	// 		for(int t = 0; t < roundsTillReblock; t++) { // TODO use insert (why use insert?)
	// 			for(int qi = 0; qi < nQs; qi++){
	// 				QsCache[worldIx][t].push_back(0);
	// 			}

	// 		//std::cout << "Context::realizeTopology QsCache[" << worldIx << "]["<< t << "] size " << QsCache[worldIx][t].size() << std::endl;
	// 		}

	// 		worldIx++;
	// 	}
	// }
}

// Return the world index in position 'which'. To be used when rotationg
std::size_t Context::getWorldIndex(std::size_t which) const
{
	return worldIndexes[which];
}

// --- Arrange different mixing parameters ---
void Context::initializeMixingParamters(){assert(!"Not implemented"); throw std::exception();}
//------------

// --- Mix ---
void Context::RotateWorlds(){assert(!"Not implemented"); throw std::exception();}
//------------

// 2D roundsTillReblock; 3D nofQs
SimTK::Real Context::Pearson(std::vector<std::vector<SimTK::Real>> inputVector, int QIx1, int QIx2)
{
	if(inputVector.size() < 1){
		std::cout << "Context::Pearson: Too few entries in the input vector" << std::endl;
		return std::numeric_limits<SimTK::Real>::min();
	}

	SimTK::Real miu0 = 0, miu1 = 0;
	SimTK::Real sqMiu0 = 0, sqMiu1 = 0;
	SimTK::Real crossMiu = 0;
	SimTK::Real var0 = 0, var1 = 0;
	SimTK::Real stdev0 = 0, stdev1 = 0;
	SimTK::Real result;

	// Get averages
	// std::cout << "Context::Pearson: inputVector " << std::endl;
	for(const auto& in : inputVector){
		if(in.size() < 2){
			std::cout << std::setprecision(1) << std::fixed;
			std::cout << "Context::Pearson: Too few Qs" << std::endl;

			return std::numeric_limits<SimTK::Real>::min();
		}

		// for(unsigned int j = 0; j < in.size(); j++){
		//     std::cout << in[j] << " ";
		// }
		// std::cout << std::endl;

		miu0 += in[QIx1];
		miu1 += in[QIx2];

		sqMiu0 += (in[QIx1] * in[QIx1]);
		sqMiu1 += (in[QIx2] * in[QIx2]);

		crossMiu += (in[QIx1] * in[QIx2]);
	}

	miu0 /= static_cast<SimTK::Real>(inputVector.size());
	miu1 /= static_cast<SimTK::Real>(inputVector.size());

	sqMiu0 /= static_cast<SimTK::Real>(inputVector.size());
	sqMiu1 /= static_cast<SimTK::Real>(inputVector.size());

	crossMiu /= static_cast<SimTK::Real>(inputVector.size());

	var0 = sqMiu0 - (miu0 * miu0);
	var1 = sqMiu1 - (miu1 * miu1);

	stdev0 = std::sqrt(var0);
	stdev1 = std::sqrt(var1);

	result = (crossMiu - (miu0 * miu1)) / (stdev0 * stdev1);

	return result;
}

/**
 *  Pass compounds to the new world
 */
void Context::passTopologiesToNewWorld(int newWorldIx)
{
	// Go through all the molecules
	for (auto& topology : topologies) {

		// Aquire the CompoundSystem
		topology.setMultibodySystem(worlds[newWorldIx].updCompoundSystem());

		// Reset mobilized body indices in Compound
		for (const auto& atom : topology.getAtoms()) {
			SimTK::Compound::AtomIndex aIx = atom.identity.compoundAtomIndex;
			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, worlds[newWorldIx].getForceField());
			topology.setAtomMobilizedBodyIndex(aIx, mbx);
		}
	}
}


////////////////////////
// REX
////////////////////////

// Set the number of replicas. This could be dangerous
void Context::setNofReplicas(const size_t& argNofReplicas)
{
	this->nofReplicas = argNofReplicas;
}

/*!
 * <!-- Adds a replica to the vector of Replica objects and sets the coordinates
 * of the replica's atomsLocations -->
*/
void Context::addReplica()
{
	if (atoms.empty()) {
		throw std::runtime_error("Context::addReplica: No atoms defined in the system.");
	}

	// Add replica and the vector of worlds
	replicas.emplace_back(Replica(nofReplicas, atoms, roots, topologies, zMatrixTable));
	nofReplicas++;

	// Set replicas coordinates
	// This also updates coordinate buffers needed to write DCD files
	replicas.back().setAtomsLocationsInGround(atomTargetLocationsCache);
	replicas.back().set_WORK_AtomsLocationsInGround(atomTargetLocationsCache);

	// Retrieve atom locations from replica and check they are consistent with the input atom positions
	// This serves to check that global to prmtop index mapping is consistent
	// Errors here mean that the DCD file will be wrong
	const auto replicaX = replicas.back().getX();
	const auto replicaY = replicas.back().getY();
	const auto replicaZ = replicas.back().getZ();

	for (const auto& atom : atoms) {
		const std::size_t prmtopIndex = atom.identity.prmtopIndex;
		const SimTK::Vec3 atomLocation(
			replicaX[atom.identity.globalIndex],
			replicaY[atom.identity.globalIndex],
			replicaZ[atom.identity.globalIndex]
		);
		const auto diff = (atomLocation - atom.position).norm();

		if (std::fabs(diff) > 10e-9) {
			const std::string msg = "Global to prmtop index mapping is inconsistent for atom " + atom.identity.uniqueAtomName +
							  " (global index " + std::to_string(atom.identity.globalIndex) +
							  ", prmtop index " + std::to_string(prmtopIndex) + "). ";
			throw std::runtime_error(msg);
		}
	}
}

void Context::addThermodynamicState(
	SimTK::Real T,
	const std::vector<AcceptRejectMode>& acceptRejectModes,
	const std::vector<int>& rexDistortOptions,
	const std::vector<std::string>& rexDistortArgs,
	const std::vector<int>& rexFlowOptions,
	const std::vector<int>& rexWorkOptions,
	const std::vector<IntegratorType>& rexIntegrators,
	const std::vector<int>& argWorldIndexes,
	const std::vector<SimTK::Real>& timestepsInThisReplica,
	const std::vector<int>& mdstepsInThisReplica)
{
	// Allocate and construct
	thermodynamicStates.emplace_back(
		ThermodynamicState(nofThermodynamicStates,
			T,
			argWorldIndexes,
			timestepsInThisReplica,
			mdstepsInThisReplica,
			atoms,			
			zMatrixTable
		)
	);
	
	// Set temperature
	thermodynamicStates.back().setTemperature(T); // seems redundant

	// Set the sampling methods
	thermodynamicStates.back().setAcceptRejectModes(acceptRejectModes);

	// Set non-equilibrium params
	thermodynamicStates.back().setDistortOptions(rexDistortOptions);
	thermodynamicStates.back().setDistortArgs(rexDistortArgs);
	thermodynamicStates.back().setFlowOptions(rexFlowOptions);
	thermodynamicStates.back().setWorkOptions(rexWorkOptions);
	thermodynamicStates.back().appendLog(baseName + ".repl" + std::to_string(nofThermodynamicStates) + ".csv");
	thermodynamicStates.back().appendDCDReporter(baseName + ".repl" + std::to_string(nofThermodynamicStates) + ".dcd", atoms.size(), topologies.size());

	// Set integrating method
	thermodynamicStates.back().setIntegrators(rexIntegrators);

	// Done
	nofThermodynamicStates++;

}

// Get the number of replicas
const size_t& Context::getNofReplicas() const
{
	return this->nofReplicas;
}

// Set the number of thermodynamic states
// Also allocates the matrix of attempted and accepted swaps
void Context::allocateSwapMatrices()
{
	// Allocate the number of attempted swaps
	nofAttemptedSwapsMatrix.resize(nofThermodynamicStates);
	for(size_t i = 0; i < nofThermodynamicStates; i++){
		nofAttemptedSwapsMatrix[i].resize(nofThermodynamicStates);
	}

	// Fill with zeros
	std::fill(nofAttemptedSwapsMatrix.begin(), nofAttemptedSwapsMatrix.end(), std::vector<int>(nofThermodynamicStates, 0));

	// Allocate the number of accepted swaps
	nofAcceptedSwapsMatrix.resize(nofThermodynamicStates);
	for(size_t i = 0; i < nofThermodynamicStates; i++){
		nofAcceptedSwapsMatrix[i].resize(nofThermodynamicStates);
	}

	// Fill with zeros
	std::fill(nofAcceptedSwapsMatrix.begin(), nofAcceptedSwapsMatrix.end(), std::vector<int>(nofThermodynamicStates, 0));

}

// Get the number of replicas
const size_t& Context::getNofThermodynamicStates() const
{
	return nofThermodynamicStates;
}

// Set the intial mapping between replicas and thermoStates
void Context::loadReplica2ThermoIxs()
{
	// Set index of replicas the same as those of the thermodynamic states
	for(size_t thermoState_k = 0; thermoState_k < nofThermodynamicStates; thermoState_k++){

		replica2ThermoIxs.insert(
			std::pair<int, int>
			(thermoState_k, thermoState_k));

	}

	for(size_t thermoState_k = 0; thermoState_k < nofThermodynamicStates; thermoState_k++){

		thermo2ReplicaIxs.insert(
			std::pair<int, int>
			(thermoState_k, thermoState_k));

	}

	// Make thermodynamic state to point to replica's BAT
	for(size_t thermoState_k = 0; thermoState_k < nofThermodynamicStates; thermoState_k++){
		thermodynamicStates[thermoState_k].setZMatrixBATPointer(
			(replicas[thermoState_k].getZMatrixBATPointer())
		);		
	}

}

void Context::setThermostatesNonequilibrium(){

	// Set index of replicas the same as those of the thermodynamic states
	for(size_t thermoState_k = 0;
	thermoState_k < nofThermodynamicStates;
	thermoState_k++){
		
		std::vector<int> distortOptions =
			thermodynamicStates[thermoState_k].getDistortOptions();
		
		for(auto distOpt : distortOptions){
			
			if(distOpt != 0){
				thermodynamicStates[thermoState_k].setNonequilibrium(1);
				std::cout << "THERMO " << thermoState_k << " nonequil" << std::endl;
			}
		}
		
	}
}

void Context::PrintReplicaMaps(){

	std::cout << "Replica -> Thermo:\n";
	for(const auto& elem : replica2ThermoIxs){
		std::cout << elem.first << " " << elem.second << "\n";
	}

	std::cout << "Thermo -> Replica:\n";
	for(const auto& elem : thermo2ReplicaIxs){
		std::cout << elem.first << " " << elem.second << "\n";
	}

}

/*!
 * <!-- Get Fixman potential already calculated from replica -->
*/
SimTK::Real Context::getFixman(int replica_i)
{
    return replicas[replica_i].getFixman();
}

// Calculate Fixman potential of replica I in replica J's back world. Uj(X_i)
SimTK::Real Context::calcFixman_IinJ(int replica_i, int replica_j)
{
	const SimTK::Real U_j = replicas[replica_j].getFixman();

	if (replica_j == replica_i) {
		// Same replica
		return U_j;

	} else if (U_j <= SimTK::Eps) {
		// Fully flexible world
		return U_j;

	}else{

		// Get replica i thermodynamic state
       	const int thermoState_j = replica2ThermoIxs[replica_j];

		// Get replica i back world
		const int world_j_front = thermodynamicStates[thermoState_j].getWorldIndexes().front();
       	const int world_j_back = thermodynamicStates[thermoState_j].getWorldIndexes().back();

		// Pass compounds to the new world
		passTopologiesToNewWorld(world_j_back);
		   
		// Transfer coordinates from replica i to back world of replica j
		const auto& X_i = replicas[replica_i].getAtomsLocationsInGround();
        worlds[world_j_back].setAtomTargetLocationsToState(X_i);

		// Transfer buffer coordinates of replica i back to back world
		// Get thermoState corresponding to this replica
		const int thermoIx = replica2ThermoIxs[replica_j];

		// Get worlds indexes of this thermodynamic state
		const std::vector<int>& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();

		// Transfer coordinates
		const auto& coords = replicas[replica_j].getAtomsLocationsInGround();
		worlds[worldIndexes.back()].setAtomTargetLocationsToState(coords);

		passTopologiesToNewWorld(world_j_front);

		// Calculate Fixman in replica i back world
		return worlds[world_j_back].calcFixman();
	}
}

void Context::swapThermodynamicStates(int replica_i, int replica_j){

	// Get replicas' thermodynamic states indexes
	int thermoState_i = replica2ThermoIxs[replica_i];
	int thermoState_j = replica2ThermoIxs[replica_j];

	// Record this swap
	nofAcceptedSwapsMatrix[thermoState_i][thermoState_j] += 1;
	nofAcceptedSwapsMatrix[thermoState_j][thermoState_i] += 1;

	// Swap thermodynamic states
	int temp = replica2ThermoIxs[replica_i];
	replica2ThermoIxs[replica_i] = replica2ThermoIxs[replica_j];
	replica2ThermoIxs[replica_j] = temp;

	// Mirror this operation in the reverse map
	temp = thermo2ReplicaIxs[thermoState_i];
	thermo2ReplicaIxs[thermoState_i] = thermo2ReplicaIxs[thermoState_j];
	thermo2ReplicaIxs[thermoState_j] = temp;

	// Swap the BAT pointers too
	thermodynamicStates[thermoState_i].setZMatrixBATPointer(
		(replicas[replica_j].getZMatrixBATPointer())
	);
}

void Context::swapPotentialEnergies(int replica_i, int replica_j)
{
	// Exchange potential energies (not necessary)
	SimTK::Real tempE = replicas[replica_i].getPotentialEnergy();
	replicas[replica_i].setPotentialEnergy(replicas[replica_j].getPotentialEnergy());
	replicas[replica_j].setPotentialEnergy(tempE);
}

void Context::swapReferencePotentialEnergies(int replica_i, int replica_j)
{
	// Exchange reference potential energies (not necessary)
	SimTK::Real tempE = replicas[replica_i].getReferencePotentialEnergy();
	replicas[replica_i].setReferencePotentialEnergy(replicas[replica_j].getReferencePotentialEnergy());
	replicas[replica_j].setReferencePotentialEnergy(tempE);
}

/*! <!-- restoreReplica --> */\
void Context::rewindReplica(void)
{
	throw std::runtime_error("Context::rewindReplica: Not implemented yet.");

	// Return to equilibrium worlds coordinates
	// - no need because it is restored in RunREX

	// Return to equilibrium worlds energies
	// - no need because it is restored in RunREX

}

/*! <!-- Attempt swap between replicas r_i and r_j
 * Code inspired from OpenmmTools
 * Chodera JD and Shirts MR. Replica exchange and expanded ensemble simulations
 * as Gibbs multistate: Simple improvements for enhanced mixing. J. Chem. Phys.
 * , 135:194110, 2011. DOI:10.1063/1.3660669
 *  replica_i and replica_j are variable
 * --> */
bool Context::attemptREXSwap(int thermoState_C, int thermoState_H)
{
	// Extract information only from Replica and ThermodynamicState objects
	// Do not use World objects

	// Get replicas' thermodynamic states indexes
	//int thermoState_C = replica2ThermoIxs[replica_X];
	//int thermoState_H = replica2ThermoIxs[replica_Y];

	const int replica_X = thermo2ReplicaIxs[thermoState_C];
	const int replica_Y = thermo2ReplicaIxs[thermoState_H];
	
	// Record this attempt
	nofAttemptedSwapsMatrix[thermoState_C][thermoState_H] += 1;
	nofAttemptedSwapsMatrix[thermoState_H][thermoState_C] += 1;

	// Convenient vars (Ballard-Jarzinski nomenclature)
	const SimTK::Real beta_C = thermodynamicStates[thermoState_C].getBeta();
	const SimTK::Real beta_H = thermodynamicStates[thermoState_H].getBeta();
	
	const SimTK::Real U_Xset = replicas[replica_X].getPotentialEnergy(); // last equil potential
	const SimTK::Real U_Yset = replicas[replica_Y].getPotentialEnergy(); // last equil potential

	const SimTK::Real refU_Xset = replicas[replica_X].getReferencePotentialEnergy(); // last equil reference potential
	const SimTK::Real refU_Yset = replicas[replica_Y].getReferencePotentialEnergy(); // last equil reference potential

	const SimTK::Real W_X = replicas[replica_X].getWORK(); // work without Jacobian
	const SimTK::Real W_Y = replicas[replica_Y].getWORK(); // work without Jacobian

	const SimTK::Real U_Xtau = replicas[replica_X].get_WORK_PotentialEnergy_New(); // last non-equil potential
	const SimTK::Real U_Ytau = replicas[replica_Y].get_WORK_PotentialEnergy_New(); // last non-equil potential

	const SimTK::Real refU_Xtau = replicas[replica_X].get_WORK_ReferencePotentialEnergy_New(); // last non-equil potential
	const SimTK::Real refU_Ytau = replicas[replica_Y].get_WORK_ReferencePotentialEnergy_New(); // last non-equil potential

	const SimTK::Real lnJac_X = replicas[replica_X].get_WORK_Jacobian(); // non-equil Jacobian
	const SimTK::Real lnJac_Y = replicas[replica_Y].get_WORK_Jacobian(); // non-equil Jacobian

	// Reduced potentials X0
	const SimTK::Real uC_Xset = beta_C * U_Xset; // Replica i reduced potential in state i
	const SimTK::Real uH_Yset = beta_H * U_Yset; // Replica j reduced potential in state j
	const SimTK::Real uH_Xset = beta_H * U_Xset; // Replica i reduced potential in state j
	const SimTK::Real uC_Yset = beta_C * U_Yset; // Replica j reduced potential in state i

	const SimTK::Real ref_uC_Xset = beta_C * refU_Xset; // Replica i reduced reference potential in state i
	const SimTK::Real ref_uH_Yset = beta_H * refU_Yset; // Replica j reduced reference potential in state j
	const SimTK::Real ref_uH_Xset = beta_H * refU_Xset; // Replica i reduced reference potential in state j
	const SimTK::Real ref_uC_Yset = beta_C * refU_Yset; // Replica j reduced reference potential in state i

	// Reduced potential Xtau
	const SimTK::Real uC_Xtau = beta_C * U_Xtau; // Replica i reduced potential in state i
	const SimTK::Real uH_Ytau = beta_H * U_Ytau; // Replica j reduced potential in state j
	const SimTK::Real uH_Xtau = beta_H * U_Xtau; // Replica i reduced potential in state j
	const SimTK::Real uC_Ytau = beta_C * U_Ytau; // Replica j reduced potential in state i

	const SimTK::Real ref_uC_Xtau = beta_C * refU_Xtau; // Replica i reduced reference potential in state i
	const SimTK::Real ref_uH_Ytau = beta_H * refU_Ytau; // Replica j reduced reference potential in state j
	const SimTK::Real ref_uH_Xtau = beta_H * refU_Xtau; // Replica i reduced reference potential in state j
	const SimTK::Real ref_uC_Ytau = beta_C * refU_Ytau; // Replica j reduced reference potential in state i

	// Get Fixman potential for the work coordinates
	const SimTK::Real Fix_Xtau = replicas[replica_X].get_WORK_Fixman();
	const SimTK::Real Fix_Ytau = replicas[replica_Y].get_WORK_Fixman();

	// Get Fixman potential for the equilibrium coordinates
	const SimTK::Real Fix_Xset = replicas[replica_X].getFixman();
	const SimTK::Real Fix_Yset = replicas[replica_Y].getFixman();

	// Get reduced Fixman potentials
	const SimTK::Real fixH_Xtau = beta_H * Fix_Xtau;
	const SimTK::Real fixC_Ytau = beta_C * Fix_Ytau;

	const SimTK::Real fixC_Xset = beta_C * Fix_Xset;
	const SimTK::Real fixH_Yset = beta_H * Fix_Yset;

	// Include the Fixman term if indicated
	SimTK::Real Fix_ii = 0, Fix_jj = 0, Fix_ij = 0, Fix_ji = 0;
	if (swapFixman) {
        if (thermoState_C == 0){
			std::cout << "Swap between " << thermoState_C << " and " << thermoState_H << " ";

            // Replica i reduced Fixman potential in state i
            Fix_ii = beta_C * calcFixman_IinJ(replica_X, replica_X);

            // Replica j reduced Fixman potential in state j
            Fix_jj = beta_H * calcFixman_IinJ(replica_Y, replica_Y);

            // Replica i reduced Fixman potential in state j
            Fix_ij = beta_H * calcFixman_IinJ(replica_X, replica_Y);

            // Replica j reduced Fixman potential in state i
            Fix_ji = beta_C * calcFixman_IinJ(replica_Y, replica_X); 
        }else{
            Fix_ii = Fix_jj = Fix_ij = Fix_ji = 0;
        }

        std::cout << "Uii Ujj Uij Uji " << Fix_ii << " " << Fix_jj << " " << Fix_ij << " " << Fix_ji << std::endl;
	}

	// LogP energy equilibrium
	const SimTK::Real ETerm_equil = -1.0 * ((ref_uH_Xset - ref_uC_Xset) + (ref_uC_Yset - ref_uH_Yset));

	// LogP energy non-equilibrium
	const SimTK::Real ETerm_nonequil = -1.0 * ((ref_uH_Xtau - ref_uC_Xtau) + (ref_uC_Ytau - ref_uH_Ytau));

	// LogP work
	// Get work from X replica
	//SimTK::Real Work_X = (ref_uH_Xtau - ref_uC_Xset) + (fixH_Xtau - fixC_Xset) - lnJac_X; // variant 1
	const SimTK::Real Work_X = (ref_uH_Xtau - ref_uC_Xset) - lnJac_X;                           // variant 2

	// Get work from Y replica
	//SimTK::Real Work_Y = (ref_uC_Ytau - ref_uH_Yset) + (fixC_Ytau - fixH_Yset) - lnJac_Y; // variant 1
	const SimTK::Real Work_Y = (ref_uC_Ytau - ref_uH_Yset) - lnJac_Y;                           // variant 2

	// Get total work
	const SimTK::Real WTerm = -1.0 * (Work_X + Work_Y);

	// // CORRECTION TERM FOR REBAS : probability of choosing
	// const SimTK::Real miu_C = qScaleFactorsMiu.at(thermoState_C);
	// const SimTK::Real miu_H = qScaleFactorsMiu.at(thermoState_H);
	// const SimTK::Real std_C = qScaleFactorsStd.at(thermoState_C);
	// const SimTK::Real std_H = qScaleFactorsStd.at(thermoState_H);
	
	const SimTK::Real s_X = qScaleFactors.at(thermoState_C);
	const SimTK::Real s_Y = qScaleFactors.at(thermoState_H);
	const SimTK::Real s_X_1 = 1.0 / s_X;
	const SimTK::Real s_Y_1 = 1.0 / s_Y;

	// Correction term is 1 for now
	const SimTK::Real qC_s_X = 1.0, qH_s_Y = 1.0, qH_s_X_1 = 1.0, qC_s_Y_1 = 1.0;
	const SimTK::Real correctionTerm = (qH_s_X_1 * qC_s_Y_1) / (qC_s_X * qH_s_Y);

	// Just for DEBUG
    // cout << "th_C, th_H, re_X, re_Y, beta_C, beta_H, U_Xset, U_Yset, U_Xtau, U_Ytau "
    //      << thermoState_C << " " << thermoState_H << " "
    //      << replica_X << " " << replica_Y << " "
    //      << beta_C << " " << beta_H << " "
    //      << refU_Xset << " "
    //      << refU_Yset << " "
    //      << refU_Xtau << " "
    //      << refU_Ytau << " "
	// 	 << Work_X << " "
	// 	 << Work_Y << " "
    //      << endl;

	const bool printTerms = false, printWithoutText = false;
	if (printTerms){
		std::cout << "thermoIxs " << thermoState_C << " " << thermoState_H << std::endl;
		std::cout << "replicaIxs " << replica_X << " " << replica_Y << std::endl;
		std::cout << "bibjwiwj " << beta_C << " " << beta_H << " " << std::endl;
		std::cout << "LiiLjj " << uC_Xtau << " " << uH_Ytau << " "
							   << uH_Xtau << " " << uC_Ytau << std::endl;
		std::cout << "EiiEjj " << uC_Xset << " " << uH_Yset << " "
							   << uH_Xset << " " << uC_Yset << std::endl;
		std::cout << "Transferred E i j " << W_X << " " << W_Y << std::endl;
		std::cout << "ETerm " << ETerm_equil << std::endl;
		std::cout << "ETerm_noneq " << ETerm_nonequil << std::endl;
		std::cout << "WTerm " << WTerm << std::endl;
		std::cout << "correctionTerm s_i s_f " << correctionTerm 
			<< " " << s_X << " " << s_Y << " " << s_X_1 << " " << s_Y_1
			<< " " << qC_s_X << " " << qH_s_Y << " " << qH_s_X_1 << " " << qC_s_Y_1
			<< std::endl;
	}
	if(printWithoutText){

		std::stringstream rexDetStream;
		rexDetStream.str("");

		rexDetStream 
			<< "REXdetails" << ", " << thermoState_C << ", " << thermoState_H << ", "
			<< replica_X << ", " << replica_Y << ", "
			<< beta_C << ", " << beta_H << ", "

			<< uC_Xset << ", " << uH_Yset << ", " << uH_Xset << ", " << uC_Yset << ", "
			<< uC_Xtau << ", " << uH_Ytau << ", " << uH_Xtau << ", " << uC_Ytau << ", "

			<< ref_uC_Xset << ", " << ref_uH_Yset << ", " << ref_uH_Xset << ", " << ref_uC_Yset << ", "
			<< ref_uC_Xtau << ", " << ref_uH_Ytau << ", " << ref_uH_Xtau << ", " << ref_uC_Ytau << ", "
			
			<< lnJac_X << ", " << lnJac_Y << ", "
			<< Work_X << ", " << Work_Y << ", "
			<< s_X << ", " << s_Y << ", " << s_X_1 << ", " << s_Y_1 << ", "
			<< qC_s_X << ", " << qH_s_Y << ", " << qH_s_X_1 << ", " << qC_s_Y_1 << ", "
			<< ETerm_equil << ", " << WTerm << ", " << correctionTerm << ", "   
		;

		std::cout << rexDetStream.str();		
	}

	// ----------------------------------------------------------------
	// EVALUATE
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real log_p_accept = 1.0;

	// Calculate log_p_accept
	if( runType == RUN_TYPE::REMC) {
		log_p_accept = ETerm_equil ;
	} else if (runType == RUN_TYPE::RENEMC) {
		log_p_accept = ETerm_nonequil + std::log(correctionTerm) ;
	} else if (runType == RUN_TYPE::RENE || runType == RUN_TYPE::REBASONTOP){
		log_p_accept = WTerm + std::log(correctionTerm) ;
	}

	// Draw from uniform distribution
	SimTK::Real unifSample = uniformRealDistribution(randomEngine);

	bool testingMode = false; 
	if (testingMode){
		std::cerr << "WARNING: REX EXCHANGE IN TESTING MODE " << std::endl;

		enum TestingWay : int {
			ALWAYS_ACCEPT = 0,
			ALWAYS_REJECT
		};
		
		const TestingWay testingWay = TestingWay::ALWAYS_ACCEPT;

		if(testingWay == TestingWay::ALWAYS_ACCEPT) {
			log_p_accept = 1.0;
		} else if (testingWay == TestingWay::ALWAYS_REJECT) {
			log_p_accept = -1.0; // std::exp(log_p_accept) = 0.3678794411714424
			unifSample = 1.0;
		}
	}

	// Accept
	if (log_p_accept >= 0.0 || unifSample < std::exp(log_p_accept)) {
		if (runType == RUN_TYPE::RENE || runType == RUN_TYPE::REBASONTOP) {
			replicas[replica_X].incrementWorldsNofSamples();
			replicas[replica_Y].incrementWorldsNofSamples();
			thermodynamicStates[thermoState_C].incrementWorldsNofSamples();
			thermodynamicStates[thermoState_H].incrementWorldsNofSamples();

			bool onlyNonequilWorlds = true;
			for(int wIx = 0; wIx < nofWorlds; wIx++) {
				if(worlds[wIx].getSampler(0)->getDistortOption() == 0) {
					onlyNonequilWorlds = false;
					break;
				}
			}

			if (onlyNonequilWorlds) {
				replicas[replica_X].incrementNofSamples();
				replicas[replica_Y].incrementNofSamples();
				thermodynamicStates[thermoState_C].incrementNofSamples();
				thermodynamicStates[thermoState_H].incrementNofSamples();
			}
			
			// Calculate replica BAT
			//replicas[replica_X].calcZMatrixBAT_WORK();
			//replicas[replica_Y].calcZMatrixBAT_WORK();
			// Calculate thermodynamic states BAT stats
			//thermodynamicStates[thermoState_C].calcZMatrixBATStats();
			//thermodynamicStates[thermoState_H].calcZMatrixBATStats();
		}

		// if (runType == RUN_TYPE::RENE || runType == RUN_TYPE::RENEMC || runType == RUN_TYPE::REBASONTOP)
		{
			
			// Update replicas coordinates from work generated coordinates
			set_WORK_CoordinatesAsFinal(replica_X);
			set_WORK_CoordinatesAsFinal(replica_Y);

			// Update replica's energy from work last potential energy
			set_WORK_PotentialAsFinal(replica_X);
			set_WORK_PotentialAsFinal(replica_Y);
		}

		// Swap thermodynamic states
		swapThermodynamicStates(replica_X, replica_Y);
		//swapPotentialEnergies(replica_X, replica_Y);
		//swapReferencePotentialEnergies(replica_X, replica_Y);

		// std::cout << "1" <<", " << unifSample << std::endl << std::endl;

		// Swap accepted
		return true;
	} else {
		// rewindReplica();

		// Return to equilibrium worlds coordinates
		// - no need because it is restored in RunREX
		// Return to equilibrium worlds energies
		// - no need because it is restored in RunREX
		// Don't swap thermodynamics states nor energies

		// std::cout << "0" <<", " << unifSample << std::endl << std::endl;

		// Swap rejected
		return false;
	}
}

/*!
 * <!--	Get printing REX swap details -->
*/
void Context::getMsg_RexDetHeader(std::stringstream& rexDetHeader) {
	rexDetHeader << "REXdetails"
		<< ", "<< "thermoState_C" << ", "<< "thermoState_H" << ", "<< "replica_X" << ", "<< "replica_Y" << ", "<< "beta_C" << ", "<< "beta_H"                                                                                                                                                                                                                                                                                           << ", "<< "uC_X0" << ", "<< "uH_Y0" << ", "<< "uH_X0" << ", "<< "uC_Y0"
		<< ", "<< "uC_Xtau" << ", "<< "uH_Ytau" << ", "<< "uH_Xtau" << ", "<< "uC_Ytau"                                                                                                                                                                                                                                                                                                                                                 << ", "<< "ref_uC_X0" << ", "<< "ref_uH_Y0" << ", "<< "ref_uH_X0" << ", "<< "ref_uC_Y0"
		<< ", "<< "ref_uC_Xtau" << ", "<< "ref_uH_Ytau" << ", "<< "ref_uH_Xtau" << ", "<< "ref_uC_Ytau"
		<< ", "<< "lnJac_X" << ", "<< "lnJac_Y"
		<< ", "<< "W_X" << ", "<< "W_Y"
		<< ", "<< "s_X" << ", "<< "s_Y" << ", "<< "s_X_1" << ", "<< "s_Y_1"
		<< ", "<< "qC_s_X" << ", "<< "qH_s_Y" << ", "<< "qH_s_X_1" << ", "<< "qC_s_Y_1"
		<< ", "<< "ETerm_equil" << ", "<< "WTerm" << ", "<< "correctionTerm"
		<< ", "<< "acc" << ", "<< "unif";
}

// Exhange all replicas
void Context::mixAllReplicas(int nSwapAttempts)
{
	std::uniform_int_distribution<std::size_t> randReplicaDistrib(0, nofReplicas-1);

	// Try nSwapAttempts to swap between random replicas
	for(size_t swap_k = 0; swap_k < nSwapAttempts; swap_k++){

		// Get two random replicas
		//auto replica_i = randReplicaDistrib(randomEngine);
		//auto replica_j = randReplicaDistrib(randomEngine);
		const std::size_t thermoState_i = randReplicaDistrib(randomEngine);
		const std::size_t thermoState_j = randReplicaDistrib(randomEngine);

		std::cout << "Attempt to swap thermoStates " << thermoState_i << " and " << thermoState_j << std::endl;

		// Attempt to swap
		attemptREXSwap(thermoState_i, thermoState_j);
	}
}

void Context::prepareExchangePairs(int rexRound, int oddity)
{
    // No need to fill exchangePairs with -1 if you use the list for the actual loop
    exchangePairList.clear();
    
    const int K = static_cast<int>(nofThermodynamicStates);
    if (K < 2) return;

    // // 4-step cyclic scheme determines if we start at index 0 (even) or 1 (odd)
    // int phase = (rexRound * 2 + oddity) % 4;
    // int startIdx = (phase == 1 || phase == 2) ? 1 : 0;

	// 2-step cyclic scheme (simpler)
	const int startIdx = (rexRound + oddity) % 2;
    for (int thIx = startIdx; thIx + 1 < K; thIx += 2)
    {
        exchangePairList.emplace_back(thIx, thIx + 1);
        
        // Only update the lookup table if other parts of the code actually use it
        exchangePairs[thIx] = thIx + 1;
        exchangePairs[thIx + 1] = thIx;
    }
}


void Context::mixReplicas(int mixi, int oddity)
{
    if ((mixi % swapEvery) != 0) return;

	// 1. Unified Guard Clause
    // If it's DEFAULT, we don't mix. If it's only 1 replica, we can't mix.
    if (runType == RUN_TYPE::DEFAULT || nofReplicas <= 1) return;

    // 2. Perform the swaps
    for (const auto& [thermoState_i, thermoState_j] : exchangePairList) {
        attemptREXSwap(thermoState_i, thermoState_j);
    }
}

// Load replica's atomLocations into it's front world
int Context::restoreReplicaCoordinatesToFrontWorld(int whichReplica)
{
	const int thermoIx = replica2ThermoIxs[whichReplica];
	const auto& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();
	const int currWorldIx = worldIndexes.front();

	const auto& coords = replicas[whichReplica].getAtomsLocationsInGround();
	worlds[currWorldIx].setAtomTargetLocationsToState(coords);

	return currWorldIx;
}

/*!
 * <!-- Load replica's atomLocations into it's back world -->
*/
void Context::restoreReplicaCoordinatesToBackWorld(int whichReplica)
{
	const int thermoIx = replica2ThermoIxs[whichReplica];
	const auto& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();
	const int backWorldIx = worldIndexes.back();

	const auto& coords = replicas[whichReplica].getAtomsLocationsInGround();
	worlds[backWorldIx].setAtomTargetLocationsToState(coords);
}

// Stores replica's front world's coordinates into it's atomsLocations
// This should always be a fully flexible world
void Context::storeReplicaCoordinatesFromFrontWorld(int whichReplica)
{
	const int thermoIx = replica2ThermoIxs[whichReplica];
	const auto& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();
	const int frontWorldIx = worldIndexes.front();

	const auto& coords = worlds[frontWorldIx].getAtomTargetLocationsCache();
	replicas[whichReplica].setAtomsLocationsInGround(coords);
}



// Store first world coordinates into replica's work coords buffer
void Context::store_WORK_CoordinatesFromFrontWorld(int whichReplica)
{
	const int thermoIx = replica2ThermoIxs[whichReplica];
	const auto& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();
	const int frontWorldIx = worldIndexes.front();

	const auto& coords = worlds[frontWorldIx].getAtomTargetLocationsCache();
	replicas[whichReplica].set_WORK_AtomsLocationsInGround(coords);
}

// Store front world potential energy into work last energy buffer of the
// replica 
void Context::store_WORK_ReplicaEnergyFromFrontWorldFull(int replicaIx)
{
	const int thermoIx = replica2ThermoIxs[replicaIx];
	const auto& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();
	const int frontWorldIx = worldIndexes.front();

	const SimTK::Real potentialEnergy = worlds[frontWorldIx].getSampler(0)->getCurrentEnergy().potential;
	replicas[replicaIx].set_WORK_PotentialEnergy_New(potentialEnergy);
}

// Store any WORK Jacobians contribution from back world
void Context::store_WORK_JacobianFromBackWorld(int replicaIx)
{
	const int thermoIx = replica2ThermoIxs[replicaIx];
	const auto& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();
	const int backWorldIx = worldIndexes.back();

	const SimTK::Real jac = worlds[backWorldIx].getSampler(0)->getDistortJacobianDetLog();
	replicas[replicaIx].set_WORK_Jacobian(jac);
}

// Get energy of the back world and store it in replica thisReplica
void Context::storeReplicaEnergyFromBackWorld(int replicaIx)
{
	const int thermoIx = replica2ThermoIxs[replicaIx];
	const auto& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();
	const int backWorldIx = worldIndexes.back();

	const SimTK::Real energy = worlds[backWorldIx].getSampler(0)->getCurrentEnergy().potential + worlds[backWorldIx].getSampler(0)->getCurrentEnergy().fixman;
	replicas[replicaIx].setPotentialEnergy(energy);
}

// Get ennergy of the front world and store it in replica thisReplica
void Context::storeReplicaEnergyFromFrontWorldFull(int replicaIx)
{
	const int thermoIx = replica2ThermoIxs[replicaIx];
	const auto& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();
	const int frontWorldIx = worldIndexes.front();

	const SimTK::Real energy = worlds[frontWorldIx].getSampler(0)->getCurrentEnergy().potential;
	replicas[replicaIx].setPotentialEnergy(energy);
}

// Get Fixman of the back world and store it in replica thisReplica
void Context::storeReplicaFixmanFromBackWorld(int replicaIx)
{
	const int thermoIx = replica2ThermoIxs[replicaIx];
	const auto& worldIndexes = thermodynamicStates[thermoIx].getWorldIndexes();
	const int backWorldIx = worldIndexes.back();

	const SimTK::Real U = worlds[backWorldIx].getSampler(0)->getCurrentEnergy().fixman;
	replicas[replicaIx].setFixman(U);
}

// Update replicas coordinates from work generated coordinates
void Context::set_WORK_CoordinatesAsFinal(int replicaIx)
{
	replicas[replicaIx].updAtomsLocationsInGround_FromWORK();
}

// Update replica's energy from work last potential energy
void Context::set_WORK_PotentialAsFinal(int replicaIx)
{
	replicas[replicaIx].setPotentialEnergy_FromWORK();
}

/*!
 * <!-- Set all of a replica's worlds' paramters -->
*/
void Context::initializeReplica(int thisReplica)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thisThermoStateIx = replica2ThermoIxs[thisReplica];

	// Get this world indexes from the corresponding thermoState
	const std::vector<int>& replicaWorldIxs = thermodynamicStates[thisThermoStateIx].getWorldIndexes();
	size_t replicaNofWorlds = replicaWorldIxs.size();

	// Set temperature for all of this replica's worlds
	SimTK::Real T = thermodynamicStates[thisThermoStateIx].getTemperature();
	
	// Get thermodynamic state from map
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].setTemperature( T );
		worlds[replicaWorldIxs[i]].setBoostTemperature( T );
	}

	// Set samplers parameters for this replica
	std::vector<SimTK::Real> replicaTimesteps = thermodynamicStates[thisThermoStateIx].getTimesteps();
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setTimestep(replicaTimesteps[i], false);
	}
	
	std::vector<int> replicaMdsteps = thermodynamicStates[thisThermoStateIx].getMdsteps();
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setMDStepsPerSample(replicaMdsteps[i]);
	}

	std::cout << "initialTss set to ";
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		std::cout << worlds[replicaWorldIxs[i]].getSampler(0)->getTimestep() << " " ;
	}
	std::cout << std::endl;

	std::cout << "initialMDSs set to ";
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		std::cout << worlds[replicaWorldIxs[i]].getSampler(0)->getMDStepsPerSample() << " " ;
	}
	std::cout << std::endl;

	// Let the thermodynamic state compute and store its equil-nonequil partitioning
	thermodynamicStates[thisThermoStateIx].computeNonequilPartitioning();	
	thermodynamicStates[thisThermoStateIx].printPartitioning(std::cout);
}

/*!
 * <!--	 -->
*/
void Context::setReplicaExchangePairs(unsigned int startingFrom)
{
	assert((startingFrom <= 1) &&
	"Replica exchange scheme has to start from 0 or 1.");

	int thermoState_i = 0;
	int thermoState_j = 1;

	// Odd scheme implies 0-N exchange
	if(startingFrom == 1){
		exchangePairs[0] = exchangePairs.size() - 1;
	}

	// Go through neighboring thermodynamic states
	for(size_t thermoState_k = startingFrom;
	thermoState_k < (nofThermodynamicStates - 1);
	thermoState_k += 2)
	{
		
		// Get thermodynamic states
		thermoState_i = thermoState_k;
		thermoState_j = thermoState_k + 1;

		// Get replicas corresponding to the thermodynamic states
		int replica_i = thermo2ReplicaIxs[thermoState_i];
		int replica_j = thermo2ReplicaIxs[thermoState_j];

		// Set the vector of exchange pairs
        exchangePairs[replica_i] = replica_j;

	}
}

/*! <!--	 -->
*/
const int Context::getThermoPair(int replicaIx)
{
	assert((exchangePairs.size() > 0) &&
	"Replica exchange pairs not set.");
	
	return exchangePairs[replicaIx];
}

// Prepare Q, U, and tau altering function parameters
void Context::PrepareNonEquilibriumParams_Q(){

	if(nofThermodynamicStates == 0){return;}
	
	// Initialize a vector of scalingFactors for scaling Qs (non-equil)
	qScaleFactorsEven.resize(nofThermodynamicStates, 1.0);
	qScaleFactorsOdd.resize(nofThermodynamicStates, 1.0);
	qScaleFactorsMiu.resize(nofThermodynamicStates, 1.0);
	qScaleFactorsStd.resize(nofThermodynamicStates, 0.0);
	qScaleFactors.resize(nofThermodynamicStates, 1.0);

	// Set the even scale factors equal to the sqrt(Ti/Tj)
	// and distribute it according the some distribution
	for(size_t thermoIx = 0; thermoIx < nofThermodynamicStates - 1; thermoIx += 4){
		// s_i = T_j
		qScaleFactorsEven.at(thermoIx)     = thermodynamicStates[thermoIx + 1].getTemperature();
		qScaleFactorsEven.at(thermoIx + 1) = thermodynamicStates[thermoIx].getTemperature();

		// s_i /= T_i
		qScaleFactorsEven.at(thermoIx)     /= thermodynamicStates[thermoIx].getTemperature();
		qScaleFactorsEven.at(thermoIx + 1) /= thermodynamicStates[thermoIx + 1].getTemperature();

		// s_i = sqrt(s_i)
		qScaleFactorsEven.at(thermoIx) = std::sqrt(qScaleFactorsEven.at(thermoIx));
		qScaleFactorsEven.at(thermoIx + 1) = std::sqrt(qScaleFactorsEven.at(thermoIx + 1));
	}

	// Set the odd scale factors equal to the sqrt(Ti/Tj)
	// and distribute it according the some distribution
	for(size_t thermoIx = 1; thermoIx < nofThermodynamicStates - 1; thermoIx += 4){

		// s_i = T_j
		qScaleFactorsOdd.at(thermoIx)     = thermodynamicStates[thermoIx + 1].getTemperature();
		qScaleFactorsOdd.at(thermoIx + 1) = thermodynamicStates[thermoIx].getTemperature();

		// s_i /= T_i
		qScaleFactorsOdd.at(thermoIx)     /= thermodynamicStates[thermoIx].getTemperature();
		qScaleFactorsOdd.at(thermoIx + 1) /= thermodynamicStates[thermoIx + 1].getTemperature();

		// s_i = sqrt(s_i)
		qScaleFactorsOdd.at(thermoIx) = std::sqrt(qScaleFactorsOdd.at(thermoIx));
		qScaleFactorsOdd.at(thermoIx + 1) = std::sqrt(qScaleFactorsOdd.at(thermoIx + 1));
	}

	for(size_t thermoIx = 0; thermoIx < nofThermodynamicStates; thermoIx++){
		std::cout << "ScaleFactor even for thermoState " << thermoIx << " "
			<< qScaleFactorsEven.at(thermoIx) << std::endl;
	}
	for(size_t thermoIx = 0; thermoIx < nofThermodynamicStates; thermoIx++){
		std::cout << "ScaleFactor odd for thermoState " << thermoIx << " "
			<< qScaleFactorsOdd.at(thermoIx) << std::endl;
	}

}

/*!
 * <!--	Set world distort parameters -->
*/
void Context::setWorldDistortParameters(int whichWorld, SimTK::Real scaleFactor)
{
		// Set the scaling factor
		HMCSampler *worldFirstSampler = (worlds[whichWorld].updSampler(0));
		worldFirstSampler->setBendStretchStdevScaleFactor(
			scaleFactor);
}

// Set thermodynamic and simulation parameters for one replica
void Context::setReplicasWorldsParameters(int thisReplica, bool alwaysAccept, bool adaptTimestep, int mixi)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thisThermoStateIx = replica2ThermoIxs[thisReplica];

	// Get this world indexes from the corresponding thermoState
	std::vector<int> replicaWorldIxs = 
		thermodynamicStates[thisThermoStateIx].getWorldIndexes();
	size_t replicaNofWorlds = replicaWorldIxs.size();

	// -------------
	// Set temperature for all of this replica's worlds
	// Get thermodynamic state from map
	// =============
	SimTK::Real T = thermodynamicStates[thisThermoStateIx].getTemperature();

	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].setTemperature( T );
		worlds[replicaWorldIxs[i]].setBoostTemperature( T );
	}

	//std::cout << "Temperature set to " << T << std::endl << std::flush;

	// -------------
	// Set sampling parameters
	// =============
	// Set sampler names
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setAcceptRejectMode(
			thermodynamicStates[thisThermoStateIx].getAcceptRejectModes()[i]
		);
	}

	// -------------
	// Set simulation parameters
	// =============

	// Set integrator
	for(std::size_t i = 0; i < replicaNofWorlds; i++){

		worlds[replicaWorldIxs[i]].updSampler(0)->setIntegratorType(
			thermodynamicStates[thisThermoStateIx].getIntegrators()[i]
		);
	}

	// Set timestep and nof MD steps
	const std::vector<SimTK::Real>& replicaTimesteps =
		thermodynamicStates[thisThermoStateIx].getTimesteps();
	const std::vector<int>& replicaMdsteps =
		thermodynamicStates[thisThermoStateIx].getMdsteps();

	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		AcceptRejectMode acceptRejectMode = thermodynamicStates[thisThermoStateIx].getAcceptRejectModes()[i];
		int MDStepsPerSample = replicaMdsteps[i];
		SimTK::Real timestep = replicaTimesteps[i];
		bool adaptiveTimestep = adaptTimestep;

		if (alwaysAccept) {
			acceptRejectMode = AcceptRejectMode::AlwaysAccept;
			MDStepsPerSample /= 10;
			timestep /= 10;
			adaptiveTimestep = false;
		}

		worlds[replicaWorldIxs[i]].updSampler(0)->setAcceptRejectMode(acceptRejectMode);
		worlds[replicaWorldIxs[i]].updSampler(0)->setMDStepsPerSample(MDStepsPerSample);
		worlds[replicaWorldIxs[i]].updSampler(0)->setTemperature(T);

		worlds[replicaWorldIxs[i]].updSampler(0)->setTimestep(timestep, adaptiveTimestep);
		if (worlds[replicaWorldIxs[i]].updSampler(0)->integratorType == IntegratorType::OMMVV) {
			worlds[replicaWorldIxs[i]].updSampler(0)->dumm->setOpenMMTimestep(timestep / 10);
		}
	}


	// SET NON_EQUIL PARAMS ------------------------- 
	// Non-equilibrium params change with every replica / thermoState
	for(std::size_t worldCnt = 0; worldCnt < replicaNofWorlds; worldCnt++){

			std::string how;
			bool randSignOpt = false;
			if(thermodynamicStates[thisThermoStateIx].getDistortOptions()[worldCnt] != 0){
				if(  thermodynamicStates[thisThermoStateIx].getDistortArgs().size()  ){
					how = thermodynamicStates[thisThermoStateIx].getDistortArgs()[worldCnt];
				}else{
					how = std::string("deterministic");
				}
			}

		// Send DISTORT_OPTION from the input to the sampler
		// std::cout <<"STUDY_Context::setReplicasWorldsParameters"
		// 	<<" world "<< replicaWorldIxs[worldCnt]
		// 	<<" DISTORT_OPT "<< thermodynamicStates[thisThermoStateIx].getDistortOptions()[worldCnt]
		// 	<<" DISTORT_ARGS "<< thermodynamicStates[thisThermoStateIx].getDistortArgs()[worldCnt]
		// 	<<" QScaleFactor "<< qScaleFactorsMiu.at(thisThermoStateIx)
		// 	<< std::endl << std::flush;

		worlds[replicaWorldIxs[worldCnt]].updSampler(0)->setDistortOption(
			thermodynamicStates[thisThermoStateIx].getDistortOptions()[worldCnt]);

		// Perturb and set scale Q scale factor
		worlds[replicaWorldIxs[worldCnt]].updSampler(0)->setBendStretchStdevScaleFactor(
			perturbScalingFactor( how, qScaleFactorsMiu.at(thisThermoStateIx), randSignOpt) );

	} // _end_ Non-equil parameters

	# pragma region REBAS_TEST
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setReplica( thisReplica );
		worlds[replicaWorldIxs[i]].updSampler(0)->setThermodynamicState( thisThermoStateIx );
	}	
	# pragma endregion REBAS_TEST


	// Print info
	// std::cout << "Timesteps set to ";
	// for(std::size_t i = 0; i < replicaNofWorlds; i++){
	// 	std::cout
	// 		<< worlds[replicaWorldIxs[i]].getSampler(0)->getTimestep()
	// 		<< " " ;
	// }
	// std::cout << std::endl;
	// std::cout << "Mdsteps set to ";
	// for(std::size_t i = 0; i < replicaNofWorlds; i++){
	// 	std::cout 
	// 		<< worlds[replicaWorldIxs[i]].getSampler(0)->getMDStepsPerSample()
	// 		<< " " ;
	// }
	// std::cout << std::endl;
	// =============

}

// TODO turn strings into enum
SimTK::Real Context::perturbScalingFactor(
	std::string how, SimTK::Real scalefactor, bool randSignOpt)
{

	// Deterministic
	if (how == "deterministic"){
		// Do nothing
	}
	
	// Truncated normal
	if (how == "Gauss"){
		SimTK::Real scaleFactorStd = 0.3;
		SimTK::Real  leftLimit = -5;
		SimTK::Real rightLimit = +5;

		std::cout << "SFdistrib Gauss "
			<< scaleFactorStd << " " << leftLimit << " " << rightLimit
			<< std::endl;

		worlds[0].updSampler(0)->convoluteVariable(
			scalefactor, "truncNormal", scaleFactorStd, leftLimit, rightLimit);
	}

	// Uniform distribution
	if (how == "uniform"){
		scalefactor = 
			worlds[0].updSampler(0)->uniformRealDistributionRandTrunc(
				0.8, 1.25); //0.625, 1.600);
	}

	// Assign a random direction: stretch or compress
	if (how == "Bernoulli"){

		SimTK::Real randDir = worlds[0].updSampler(0)->uniformRealDistribution_m1_1(randomEngine);
		scalefactor = (randDir > 0) ? scalefactor : (1.0/scalefactor) ;

		// std::cout <<"STUDY_Context::perturbScalingFactor"
		// 	<<" randDir " << randDir
		// 	<<" scalefactor " << scalefactor
		// 	<< std::endl << std::flush;
	}

	// Assign a random sign (optional)
	if(randSignOpt){
		SimTK::Real randSign;
		SimTK::Real randUni_m1_1 =
			worlds[0].updSampler(0)->uniformRealDistribution_m1_1(randomEngine);
		randSign = (randUni_m1_1 > 0) ? 1 : -1 ;
		scalefactor *= randSign;
	}
	
	return scalefactor;
}

// Set nonequilibrium parameters for one replica
void Context::updWorldsDistortOptions(int thisReplica)
{

	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thisThermoStateIx = replica2ThermoIxs[thisReplica];

	// Get this world indexes from the corresponding thermoState
	std::vector<int> replicaWorldIxs = 
		thermodynamicStates[thisThermoStateIx].getWorldIndexes();
	size_t replicaNofWorlds = replicaWorldIxs.size();

	// SET NON_EQUIL PARAMS ------------------------- 
	// Non-equilibrium params change with every replica / thermoState

	for(std::size_t i = 0; i < replicaNofWorlds; i++){

		// Send DISTORT_OPTION from the input to the sampler
		worlds[replicaWorldIxs[i]].updSampler(0)->setDistortOption(
			thermodynamicStates[thisThermoStateIx].getDistortOptions()[i]
		);		

		// Set scale Q scale factor
		setWorldDistortParameters(replicaWorldIxs[i],
			qScaleFactors.at(thisThermoStateIx));

	}

}

/*!
 * <!-- Rewind back world -->
*/
void Context::RewindBackWorld(int thisReplica)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thisThermoStateIx = replica2ThermoIxs[thisReplica];

	// Get this world indexes from the corresponding thermoState
	std::vector<int> replicaWorldIxs = 
		thermodynamicStates[thisThermoStateIx].getWorldIndexes();

	// == TRANSFER == coordinates from last world to current
	// TODO: eliminate in the last iteration
	int frontIx = replicaWorldIxs.front();
	int backIx = replicaWorldIxs.back();
	if(replicaWorldIxs.size() > 1) {
		transferCoordsFromWorldToWorld(frontIx, backIx);
	}

	// == ROTATE == worlds indices (translate from right to left)
	std::rotate(replicaWorldIxs.begin(),
		replicaWorldIxs.begin() + 1,
		replicaWorldIxs.end());

}

// Run front world, rotate and transfer
int Context::RunFrontWorldAndRotate(std::vector<int> & worldIxs)
{
	bool validated = false;

	int frontWorldIx = -1;
	int backWorldIx = -1;

	// == SAMPLE == from the front world
	frontWorldIx = worldIxs.front();
	validated = RunWorld(frontWorldIx, "");

	// Write pdbs every world
	//writePdbs(nofRounds, frontWorldIx);
	
	// == ROTATE == worlds indices (translate from right to left)
	std::rotate(worldIxs.begin(),
		worldIxs.begin() + 1,
		worldIxs.end());

	// == TRANSFER == coordinates from back world to front
	frontWorldIx = worldIxs.front();
	backWorldIx = worldIxs.back();

	if(worldIxs.size() > 1) {

		// spacedcout("[YDIRBUG]");
		// std::cout << "Transfer from world " << backWorldIx << " to " << frontWorldIx ;
		// spacedcout("[YDIRBUG]"); ceol;

		transferCoordsFromWorldToWorld(backWorldIx, frontWorldIx);


		#ifdef PRINTALOT 
			if(validated){
				std::cout << std::endl;
			}else{
				std::cout << " invalid sample." << std::endl;
			}
		#endif

	}

	return worldIxs.front();

}

/**
 * Update the scale factors
 */ 
void Context::updThermostatesQScaleFactors(int mixi)
{

	// Prepare non-equilibrium scale factors
	if(mixi % 2){ // odd batch
		qScaleFactorsMiu = qScaleFactorsOdd;
	}else{ // even batch
		qScaleFactorsMiu = qScaleFactorsEven;
	}

	// Get scaling factor
	qScaleFactors = qScaleFactorsMiu;

	// Random sign for the scaling factors
	// bool randSignOpt = false;
	// The names of the probability distributions operators
	//std::vector<std::string> how;
	// Go through all thermodynamic states which should correspond to
	// scale factors
	// for(size_t thermoIx = 0; thermoIx < qScaleFactors.size(); thermoIx++){
	// 	std::vector<int>& worldIxs =  thermodynamicStates[thermoIx].updWorldIndexes();
	// 	size_t thermoNofWorlds = worldIxs.size();
	// 	for(std::size_t worldCnt = 0; worldCnt < thermoNofWorlds; worldCnt++){
	// 		if(thermodynamicStates[thermoIx].getDistortOptions()[worldCnt] != 0){
	// 			if(  thermodynamicStates[thermoIx].getDistortArgs().size()  ){
	// 				how = split(thermodynamicStates[thermoIx].getDistortArgs()[worldCnt], "_");
	// 			}else{
	// 				how = {"deterministic"};
	// 			}
	// 			break;
	// 		}
	// 		// Distribute scale factor
	// 		if(qScaleFactors.at(thermoIx) != 1){ // This is questionable
	// 			qScaleFactors.at(thermoIx) = perturbScalingFactor( how, qScaleFactorsMiu.at(thermoIx), randSignOpt);
	// 			std::cout <<"STUDY_Context::updQScaleFactors"
	// 				<<" thermoIx "<< thermoIx
	// 				<<" qScaleFactors.at(thermoIx) "<< qScaleFactors.at(thermoIx)
	// 				<< std::endl;
	// 		}
	// 	} // _end_ for worldCnt
	// } // _end_ for thermoIx

}

// Print to log and write pdbs
void Context::REXLog(int mixi, int replicaIx)
{
	// Write energy and geometric features to logfile
	if(printFreq || pdbRestartFreq){
		if( !(mixi % printFreq) ){

			for(auto wIx: worldIndexes){
				PrintToLog(replicaIx, wIx, 0);
			}

		}
		// Write pdb
		if( pdbRestartFreq != 0){
			if((mixi % pdbRestartFreq) == 0){
				writePdbs(mixi, replica2ThermoIxs[replicaIx]);
			}
		}
	} // wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
}

// rexnewfunc
void Context::incrementNofSamples(void){

	for(size_t rk = 0; rk < nofReplicas; rk++){
		replicas[rk].incrementNofSamples();
	}

	for(size_t tk = 0; tk < nofThermodynamicStates; tk++){
		thermodynamicStates[tk].incrementNofSamples();
	}

}

/*!
 * <!--  -->
*/
void Context::transferQStatistics(int thermoIx, int srcStatsWIx, int destStatsWIx)
{
	auto* sampler = worlds[destStatsWIx].updSampler(0);

	sampler->set_BMps_means(thermodynamicStates[thermoIx].getBMps_means(srcStatsWIx));
	sampler->set_PFrs_means(thermodynamicStates[thermoIx].getPFrs_means(srcStatsWIx));

	sampler->set_dBMps(thermodynamicStates[thermoIx].get_dBMps(srcStatsWIx));
	sampler->set_dPFrs(thermodynamicStates[thermoIx].get_dPFrs(srcStatsWIx));
	sampler->setPreviousQs(thermodynamicStates[thermoIx].getCurrentQs(srcStatsWIx));
	sampler->setQmeans(thermodynamicStates[thermoIx].getQmeans(srcStatsWIx));
	sampler->setQdiffs(thermodynamicStates[thermoIx].getQdiffs(srcStatsWIx));
	sampler->setQvars(thermodynamicStates[thermoIx].getQvars(srcStatsWIx));
}

/*!
 * <!-- Run a particular world -->
*/
bool Context::RunWorld(int whichWorld, const std::string& header)
{
	// Prepare output
	std::stringstream worldOutStream;
	worldOutStream.str(""); // empty

	// == SAMPLE == from the current world
	bool validated = false;
	const int numSamples = worlds[whichWorld].getSamplesPerRound();
	const int distortOption = worlds[whichWorld].getSampler(0)->getDistortOption();

	// Equilibrium world
	if(distortOption == 0) {

		// Generate samples
		// std::cout << "[EQ] World " << whichWorld
		// 	<< " generating " << numSamples << " samples." << std::endl;
		validated = worlds[whichWorld].generateSamples(numSamples, worldOutStream, header, verbose);
		// std::cout << "[EQ] World " << whichWorld
		// 	<< " generated " << numSamples << " samples." << std::endl;

		// size_t wIx = 1; // We want the U and UDot of the torsional dynamics world
		// if (whichWorld == wIx) {
		// 	SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[wIx].updForceField());
		// 	SimTK::SimbodyMatterSubsystem& matter = *(worlds[wIx].matter);
			
		// 	// auto& someState = worlds[whichWorld].updIntegrator().updAdvancedState();
		// 	// SimTK::Vector udot = someState.getUDot();
        //     // SimTK::Vector torques;
        //     // torques.resize(someState.getNU());
        //     // matter.multiplyByM(someState, udot, torques);

		// 	// std::cout << "Torsional dynamics world UDot: " << std::endl;
		// 	// std::cout << torques.sum() << std::endl;
		// 	// std::cout << "Torsional dynamics world torques: " << std::endl;

		// 	// UDotCache = std::vector<SimTK::Real>(torques.size());
		// 	// for (size_t i = 0; i < torques.size(); i++) {
		// 	// 	UDotCache[i] = torques[i];
		// 	// 	UCache.push_back(0); // U is not used in this case, so we just push 0
		// 	// }

		// 	// const auto& backupU = someState.getU();
		// 	// someState.updU() = 0; // SET VELOCITIES TO ZERO
		// 	// worlds[wIx].updSampler(0)->system->realize(someState, SimTK::Stage::Acceleration);

		// 	// // print all bonds
		// 	// for (const auto& bond : orderedBonds[0]) {
		// 	// 	std::cout << "BondLink: " << bond.first << " - " << bond.second << std::endl;
		// 	// }

		// 	// Iterate molecules
		// 	for(size_t topoIx = 0; topoIx < numMolecules; topoIx++){

		// 		// Get molecule and it's bonds
		// 		Topology& topology = topologies[topoIx];
		// 		const std::vector<BOND>& BONDS = orderedBonds[topoIx];

		// 		// Iterate molecule's bonds
		// 		for(size_t BOIx = 0; BOIx < BONDS.size(); BOIx++){

		// 			// Get current bond
		// 			const BOND& currBOND = BONDS[BOIx];
		// 			size_t boIx = BONDS_to_bonds[topoIx][BOIx];
		// 			BondLink& bond = bonds[boIx];

		// 			// Get bond's atoms
		// 			Atom& childAtom  = atoms[currBOND.first];
		// 			Atom& parentAtom = atoms[currBOND.second];

		// 			SimTK::Compound::AtomIndex child_cAIx = childAtom.identity.compoundAtomIndex;
		// 			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.identity.compoundAtomIndex;

		// 			SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
		// 			SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);

		// 			SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);
		// 			SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);

		// 			const SimTK::MobilizedBody& childMobod = matter.getMobilizedBody(childMbx);
		// 			const SimTK::MobilizedBody& parentMobod = matter.getMobilizedBody(parentMbx);



		// 			// int min_aix = std::min(currBOND.first, currBOND.second);
		// 			// int max_aix = std::max(currBOND.first, currBOND.second);
		// 			// std:: cout << "bond " << min_aix << " " << max_aix << " has childMbx: " << childMbx << " and parentMbx " << parentMbx << std::endl;

		// 			// AtomIndex1.push_back(max_aix);
		// 			// if (min_aix == 234 && max_aix == 244) {
		// 			// 	std::cout << "bond 234 244 has mbx: " << childMbx << " and parentMbx " << parentMbx << std::endl;
		// 			// }
		// 			// if (min_aix == 238 && max_aix == 241) {
		// 			// 	std::cout << "bond 234 244 has mbx: " << childMbx << " and parentMbx " << parentMbx << std::endl;
		// 			// }
		// 			// if (min_aix == 266 && max_aix == 269) {
		// 			// 	std::cout << "bond 234 244 has mbx: " << childMbx << " and parentMbx " << parentMbx << std::endl;
		// 			// }
		// 			// if (min_aix == 472 && max_aix == 474) {
		// 			// 	std::cout << "bond 234 244 has mbx: " << childMbx << " and parentMbx " << parentMbx << std::endl;
		// 			// }
		// 			// if (min_aix == 551 && max_aix == 554) {
		// 			// 	std::cout << "bond 234 244 has mbx: " << childMbx << " and parentMbx " << parentMbx << std::endl;
		// 			// }
		// 			// if (min_aix == 640 && max_aix == 642) {
		// 			// 	std::cout << "bond 234 244 has mbx: " << childMbx << " and parentMbx " << parentMbx << std::endl;
		// 			// }


		// 			if ((childMbx != parentMbx)) { //  || true
		// 				if (!binaryFileIsInitialized) {
		// 					int min_aix = std::min(currBOND.first, currBOND.second);
		// 					AtomIndex0.push_back(min_aix);

		// 					int max_aix = std::max(currBOND.first, currBOND.second);
		// 					AtomIndex1.push_back(max_aix);

		// 					// if (min_aix == 234 && max_aix == 244) {
		// 					// 	std::cout << "bond 234 244 has mbx: " << childMbx << std::endl;
		// 					// }

		// 					// if (min_aix == 238 && max_aix == 241) {
		// 					// 	std::cout << "bond 234 244 has mbx: " << childMbx << std::endl;
		// 					// }

		// 					// if (min_aix == 266 && max_aix == 269) {
		// 					// 	std::cout << "bond 234 244 has mbx: " << childMbx << std::endl;
		// 					// }

		// 					// if (min_aix == 472 && max_aix == 474) {
		// 					// 	std::cout << "bond 234 244 has mbx: " << childMbx << std::endl;
		// 					// }

		// 					// if (min_aix == 551 && max_aix == 554) {
		// 					// 	std::cout << "bond 234 244 has mbx: " << childMbx << std::endl;
		// 					// }

		// 					// if (min_aix == 640 && max_aix == 642) {
		// 					// 	std::cout << "bond 234 244 has mbx: " << childMbx << std::endl;
		// 					// }



		// 					// if (min_aix == 4 && max_aix == 6) {
		// 					// 	std::cout << "childMbx - 2: " << childMbx - 2 << std::endl;
		// 					// }
		// 				}

		// 				// const SimTK::Transform X_GP = parentMobod.getBodyTransform(someState); // Transform from G to P
		// 				// const SimTK::Transform X_PG = ~X_GP; // Transform from P to G

		// 				// const SimTK::Transform X_GB = childMobod.getBodyTransform(someState); // Transform from G to B
		// 				// const SimTK::Transform X_BG = ~X_GB; // Transform from B to G

		// 				// const SimTK::Inertia I_PB_P = childMobod.calcBodyInertiaAboutAnotherBodyStation(someState, parentMobod, Vec3(0, 0, 0)); // Inertia expressed in P
		// 				// const SimTK::Vec3 b_PB_P = childMobod.findBodyAngularAccelerationInAnotherBody(someState, parentMobod); // In P

		// 				// const SimTK::Vec3 Torque_P = I_PB_P * b_PB_P; // Torque in P

		// 				// const SimTK::Transform& X_PM = parentMobod.getInboardFrame(someState); // Mobilizer frame M, expressed in P
		// 				// const SimTK::UnitVec3 pinAxis_G = X_PM.R().z(); // z-axis of frame M, expressed in P

		// 				// const SimTK::Real u_dot = dot(Torque_P, pinAxis_G); // Angular acceleration projected onto pin axis

		// 				// UCache.push_back(0);
		// 				// UDotCache.push_back(u_dot);
		// 			}
		// 		} // every bond
		// 	} // every molecule

		// 	// Print two rows for atom indices
		// 	// The last column in each row is whether it was accepted
		// 	// Since this is not simulation, it will write 0
		// 	if (!binaryFileIsInitialized) {
		// 		initializeBinaryFile(foutU, AtomIndex0.size() + 1);
		// 		writeRowToBinaryFile(foutU, AtomIndex0, false, false);
		// 		writeRowToBinaryFile(foutU, AtomIndex1, false, false);

		// 		initializeBinaryFile(foutUDot, AtomIndex0.size() + 1);
		// 		writeRowToBinaryFile(foutUDot, AtomIndex0, false, false);
		// 		writeRowToBinaryFile(foutUDot, AtomIndex1, false, false);

		// 		initializeBinaryFile(foutTorque, AtomIndex0.size() + 1);
		// 		writeRowToBinaryFile(foutTorque, AtomIndex0, false, false);
		// 		writeRowToBinaryFile(foutTorque, AtomIndex1, false, false);

		// 		binaryFileIsInitialized = true;
		// 	}

		// 	writeRowToBinaryFile(foutU, worlds[wIx].updSampler(0)->UCache, true, validated);
		// 	writeRowToBinaryFile(foutUDot, worlds[wIx].updSampler(0)->UDotCache, true, validated);
		// 	writeRowToBinaryFile(foutTorque, worlds[wIx].updSampler(0)->TorqueCache, true, validated);

		// 	// someState.updU() = backupU; // Restore velocities
		// }
	// Non-equilibrium world
	} else if (distortOption != 0) {

		// // Generate samples
		
		// // drl
		// #ifdef __DRILLING__ // SCALEQ

        //     // Get drl data
        //     const std::vector<std::vector<double>>& drl_bon_Energies = worlds[whichWorld].getEnergies_drl_bon();
        //     const std::vector<std::vector<double>>& drl_ang_Energies = worlds[whichWorld].getEnergies_drl_ang();
        //     const std::vector<std::vector<double>>& drl_tor_Energies = worlds[whichWorld].getEnergies_drl_tor();
        //     const std::vector<std::vector<double>>& drl_n14_Energies = worlds[whichWorld].getEnergies_drl_n14();
        //     const std::vector<std::vector<double>>& drl_vdw_Energies = worlds[whichWorld].getEnergies_drl_vdw();
        //     const std::vector<std::vector<double>>& drl_cou_Energies = worlds[whichWorld].getEnergies_drl_cou();

        //     // validated = worlds[whichWorld].generateSamples(numSamples, worldOutStream){

        //         //warn("under drilling conditions");

        //         // Update Robosample bAtomList
        //         SimTK::State& currentAdvancedState = worlds[whichWorld].updIntegrator().updAdvancedState();
        //         worlds[whichWorld].updateAtomListsFromSimbody(currentAdvancedState); // Update Robosample bAtomList
        //         // ''''''''''''''''''''
        //         // coutspaced("SCALING_BAT init:"); ceolf;
        //         // replicas[0].calcZMatrixBAT( worlds[whichWorld].getAtomsLocationsInGround( worlds[whichWorld].updIntegrator().updAdvancedState() ));
        //         // thermodynamicStates[0].PrintZMatrixBAT();
        //         // ''''''''''''''''''''

        //         // // Reinitialize the sampler
        //         // worlds[whichWorld].updSampler(0)->reinitialize(worldOutStream, verbose);

        //         SimTK::Real pe_beforeScale = worlds[whichWorld].getForces().getMultibodySystem().calcPotentialEnergy(worlds[whichWorld].getIntegrator().getAdvancedState());

        //         if(false && ((whichWorld == 3)
        //                 //&& (std::abs(worlds[whichWorld].updSampler(0)->QScaleFactor - 1.0) > 0.00001)
        //         )){ 
        //             scout("[SCALING_PES]: before") <<" " << pe_beforeScale << eolf;
        //             scout("drl_bon_E"); ceol; PrintCppVector(drl_bon_Energies, 6, "bonE", "bonE");
        //             scout("drl_ang_E"); ceol; PrintCppVector(drl_ang_Energies, 6, "angE", "angE");
        //             scout("drl_tor_E"); ceol; PrintCppVector(drl_tor_Energies, 6, "torE", "torE");
        //             scout("drl_n14_E"); ceol; PrintCppVector(drl_n14_Energies, 6, "n14E", "n14E");
        //             scout("drl_vdw_E"); ceol; PrintCppVector(drl_vdw_Energies, 6, "vdwE", "vdwE");
        //             scout("drl_cou_E"); ceol; PrintCppVector(drl_cou_Energies, 6, "couE", "couE");
        //             std::cout<<std::flush;
        //         } // __end__ choose a world to print drilling

        //         auto runSamplingLoop = [&](SimTK::State& state) {
        //             for (int sampleIx = 0; sampleIx < numSamples; ++sampleIx) {
        //                 if (verbose) {
        //                     worldOutStream << header << " ";
        //                     worlds[whichWorld].updSampler(0)->getMsg_InitialParams(worldOutStream);
        //                 }
        //                 validated = worlds[whichWorld].updSampler(0)->sample_iteration(worldOutStream, verbose) && validated;
        //                 if (verbose) {worldOutStream << std::endl;}
        //             }
        //         };              

        //         // GENERATE the requested number of samples
        //         if (worlds[whichWorld].getIsRollFlexibilities()) {
        //             for (int mobIntIx = 1; mobIntIx < worlds[whichWorld].getMatterSubsystem().getNumBodies(); ++mobIntIx) {
        //                 worlds[whichWorld].lockAllMobilizers();
        //                 const SimTK::MobilizedBody& mobod = worlds[whichWorld].getMatterSubsystem().getMobilizedBody(SimTK::MobilizedBodyIndex(mobIntIx));
        //                 mobod.unlock(currentAdvancedState);
        //                 runSamplingLoop(currentAdvancedState);
        //             }
        //         } else {
        //             runSamplingLoop(currentAdvancedState);
        //         }

        //         SimTK::Real pe_afterScale = worlds[whichWorld].getForces().getMultibodySystem().calcPotentialEnergy(worlds[whichWorld].getIntegrator().getAdvancedState());

        //         // ''''''''''''''''''''
        //         // coutspaced("SCALING_BAT after:"); ceolf;
        //         // replicas[0].calcZMatrixBAT( worlds[whichWorld].getAtomsLocationsInGround( worlds[whichWorld].updIntegrator().updAdvancedState() ));
        //         // thermodynamicStates[0].PrintZMatrixBAT();
        //         // ''''''''''''''''''''
        //         if(false && ((whichWorld == 3)
        //                 //&& (std::abs(worlds[whichWorld].updSampler(0)->QScaleFactor - 1.0) > 0.00001)
        //         )){
        //             scout("[SCALING_PES]: after") <<" " << pe_afterScale << eolf;
        //             scout("drl_bon_E"); ceol; PrintCppVector(drl_bon_Energies, 6, "bonE", "bonE");
        //             scout("drl_ang_E"); ceol; PrintCppVector(drl_ang_Energies, 6, "angE", "angE");
        //             scout("drl_tor_E"); ceol; PrintCppVector(drl_tor_Energies, 6, "torE", "torE");
        //             scout("drl_n14_E"); ceol; PrintCppVector(drl_n14_Energies, 6, "n14E", "n14E");
        //             scout("drl_vdw_E"); ceol; PrintCppVector(drl_vdw_Energies, 6, "vdwE", "vdwE");
        //             scout("drl_cou_E"); ceol; PrintCppVector(drl_cou_Energies, 6, "couE", "couE");
        //             std::cout<<std::flush;
        //         } // __end__ choose a world to print drilling
                
        //     // }

		// #else

		// 	validated = worlds[whichWorld].generateSamples(numSamples, worldOutStream, header, verbose); // =

		// #endif

	}

	// Print geometry to output stream too
	# pragma region REBAS_TEST
		worldOutStream << " ";

						for (const auto& distanceIx : distanceIxs) {
							if( distanceIx[0] == whichWorld ) {
								worldOutStream << std::fixed << std::setprecision(3) << Distance(distanceIx[0], distanceIx[1], 0, distanceIx[2], distanceIx[3]) << " ";
							}
						}

						for (const auto& angleIx : angleIxs){
							if( angleIx[0] == whichWorld ) {
								worldOutStream << std::fixed << std::setprecision(3) << Roboangle(angleIx[0], angleIx[1], 0, angleIx[2], angleIx[3], angleIx[4]) << " ";
							}
						}

						for (const auto& dihedralIx : dihedralIxs){
							if( dihedralIx[0] == whichWorld ) {
								worldOutStream << std::fixed << std::setprecision(3) << Dihedral(dihedralIx[0], dihedralIx[1], 0, dihedralIx[2], dihedralIx[3], dihedralIx[4], dihedralIx[5]) << " ";

								// std::cout <<"STUDY_Context::RunWorld" 
								// <<" | "<< dihedralIx[0] <<" "<< dihedralIx[1] <<" "<< dihedralIx[2] <<" "<< dihedralIx[3] <<" "<< dihedralIx[4] <<" "<< dihedralIx[5]
								// <<" | "<< atoms[dihedralIx[0]].getInName() <<" "<< atoms[dihedralIx[1]].getInName() <<" "<< atoms[dihedralIx[2]].getInName() <<" "<< atoms[dihedralIx[3]].getInName()
								// <<" | "<< Dihedral(dihedralIx[0], dihedralIx[1], 0, dihedralIx[2], dihedralIx[3], dihedralIx[4], dihedralIx[5])
								// << std::endl;

							}
						}
	# pragma endregion REBAS_TEST

	// Print the world output stream
	if (verbose) {
		std::cout << worldOutStream.str() << std::flush;
	}

	return validated;
}


/*! <!--  -->*/
void Context::RunReplicaWorldRange(int replicaIx, int startWorldCnt, int nofWorldsCounted, bool isNonEquilibrium)
{
	// Get thermodynamic state and its' worlds
	Replica& replica = replicas[replicaIx];
	const int thermoIx = replica2ThermoIxs[replicaIx];
	ThermodynamicState& thermoState = thermodynamicStates[thermoIx];
	const std::vector<int>& thermoWorldIxs = thermoState.updWorldIndexes();
	const std::vector<int>& distortOpts = thermoState.getDistortOptions();
	
	replica.updWORK() = 0.0;
	replica.upd_WORK_Jacobian() = 0.0;
	
	// Loop through all worlds within the thermodynamic schedule
	// This is an index into the thermodynamic state’s world list
	for(const int worldScheduleIndex : thermoState.updWorldIndexes()) {

		// Get the physical world object in the simulation
		const int worldIndex = thermoWorldIxs[worldScheduleIndex];
		World& currWorld = worlds[worldIndex];
		const int distortIx = distortOpts[worldScheduleIndex];

		// Don't transfer if it's the first world in the range
		const bool firstWorldInRange = (worldScheduleIndex == startWorldCnt);

		// Also don't transfer if it's the first world in the range and it's also the first world overall
		const bool firstWorldOverall = (worldScheduleIndex == 0);

		if (!firstWorldInRange && !firstWorldOverall) {
			transferCoordsFromWorldToWorld(thermoWorldIxs[worldScheduleIndex - 1], worldIndex);
		}

		// Header
		std::string headerToRunWorld = "REX";
		headerToRunWorld += ", " + std::to_string(replicaIx);
		headerToRunWorld += ", " + std::to_string(thermoIx);
		headerToRunWorld += ", " + std::to_string(worldIndex);
				
		// Run
		std::cout << "Running world " << worldIndex <<
			" withdistortIx=" << distortIx <<
			" for replica " << replicaIx <<
			" at thermodynamic state " << thermoIx <<
			" with temperature=" << thermodynamicStates[thermoIx].getTemperature() << std::endl;
			
		bool validated = RunWorld(worldIndex, headerToRunWorld);
		// std::cout << "World " << worldIndex << " validated: " << validated << std::endl;

		// Transfer coordinates
		const bool isEquilibrium = (distortIx == 0);
		const bool intoWORK = !isEquilibrium;
		transferCoordsFromWorldToReplica(worldIndex, replicaIx, intoWORK);

		if (isEquilibrium) {
			// Calculate Q statistics
			if (worlds[worldIndex].getSampler(0)->getAcc()) {
				thermoState.calcQStats(worldIndex, currWorld.getBMps(), currWorld.getPFrs(), currWorld.getAdvancedQs(), currWorld.getNofSamples());
			} else {
				thermoState.calcQStats(worldIndex, currWorld.getBMps(), currWorld.getPFrs(), SimTK::Vector(currWorld.getNQs(), SimTK::Real(0)), currWorld.getNofSamples());
			}
		}

		// Increment the nof samples for replica and thermostate
		replica.incrementWorldsNofSamples(1);
		thermoState.incrementWorldsNofSamples(1);
	}
}

void Context::writeLog(int mixi, int replicaIx) {
	// Check if we want to write to log

	int thermoIx = replica2ThermoIxs[replicaIx];
	const auto& coords = replicas[replicaIx].getAtomsLocationsInGround();
	const auto& replica = replicas[replicaIx];

	int whichLog = replica2ThermoIxs[replicaIx];
	if (!thermodynamicStates[whichLog].logFile.is_open()) {
		return;
	}

	// Print stats for each world
	for (const auto wIx : worldIndexes) {
		SimTK::State& currentAdvancedState = worlds[wIx].updIntegrator().updAdvancedState();
		const auto& sampler = pHMC((worlds[wIx].getSampler(0)));

		const auto temperature = sampler->getTemperature();
		const auto NU = currentAdvancedState.getNU();
		const auto acceptedSteps = sampler->acceptedSteps;
		const auto pe_o = sampler->getPreviousEnergy().potential;
		const auto pe_n = sampler->getCurrentEnergy().potential;
		const auto ke_o = sampler->getPreviousEnergy().kinetic;
		const auto ke_n = sampler->getCurrentEnergy().kinetic;
		const auto fix_o = sampler->getPreviousEnergy().fixman;
		const auto fix_n = sampler->getCurrentEnergy().fixman;
		const auto timestep = sampler->getTimestep();
		const auto mdstep = sampler->getMDStepsPerSample();
		const SimTK::Real acc = sampler->numAccepted_period / static_cast<SimTK::Real>(sampler->numSamples_period);

		sampler->numAccepted_period = 0;
		sampler->numSamples_period = 0;

		// Write to log
		// round_ix replica_ix temperature world_ix NU accepted_steps pe_o pe_set ke_o ke_n fix_o fix_n fix_set timestep mdstep acc
		thermodynamicStates[whichLog].logFile
				<< std::fixed << std::setprecision(0) << mixi << ","
                << replicaIx << ","
                << std::fixed << std::setprecision(3) << temperature << ","
                << std::fixed << std::setprecision(0) << wIx << ","
                << NU << ","
                << acceptedSteps << ","
                << std::fixed << std::setprecision(2)
				<< pe_o << ","
				<< pe_n << ","
                << ke_o << ","
                << ke_n << ","
                << fix_o << ","
                << fix_n << ","
				<< timestep << ","
				<< mdstep << ","
				<< acc << std::endl;
	}
}

std::stringstream worldOutStream_TEST;

/*!
 * <!-- Run replica exchange protocol -->
*/
void Context::RunREX(int equilRounds, int prodRounds)
{
	// They all start with replica 0 coordinates
	// TODO does not work in debug
	for (int worldIx = 0; worldIx < worlds.size(); worldIx++) {
		World& world = worlds[worldIx];

		// Add this worlds BAT coordinates to it's samplers
		addSubZMatrixBATsToWorld(worldIx, 0);
		//scout("Context::initializeFromFile PrintSubZMatrixBAT: ") << eol;
		//world.updSampler(0)->PrintSubZMatrixBAT();
	}

	// Set a vector of replica pairs for exchanges
	exchangePairs.resize(nofReplicas);

	// Consider renaming
	loadReplica2ThermoIxs();
	PrintReplicas();

	// Initialize non-equilibrium parameters
	PrepareNonEquilibriumParams_Q();
	setThermostatesNonequilibrium();

	// Allocate space for swap matrices
	allocateSwapMatrices();

	// Initialize replicas
	for (size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++){
		initializeReplica(replicaIx);

	}
	
	// Allocate and calculate initial Q statistics for all thermodynamic states
	for (auto& thermoState : thermodynamicStates) {
		for (const auto worldIx : thermoState.getWorldIndexes()) {
			World& currWorld = worlds[worldIx];
			thermoState.calcQStats(worldIx, currWorld.getBMps(), currWorld.getPFrs(), SimTK::Vector(currWorld.getNQs(), SimTK::Real(0)), currWorld.getNofSamples());
		}
	}

	// Print a header =========================================================
	std::stringstream rexOutput;
	rexOutput.str("");
	rexOutput << "REX, " << "replicaIx" << ", " << "thermoIx" << ", " << "wIx" ;
	worlds[0].getSampler(0)->getMsg_Header(rexOutput);
	rexOutput << std::endl;
	getMsg_RexDetHeader(rexOutput);
	std::cout << rexOutput.str() << std::endl;

	std::cout << rexOutput.str() << std::endl;
	// ------------------------------------------------------------------------

	// REPLICA EXCHANGE MAIN LOOP -------------------------------------------->
	std::size_t mixIndex = 0;
	for (std::size_t cycleIndex = 0; cycleIndex < equilRounds + prodRounds; cycleIndex++) {

		// SIMULATE EACH REPLICA --------------------------------------------->
		for (std::size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++) { 

			// Get thermodynamic state and its' worlds
			Replica& replica = replicas[replicaIx];
			const int thermoIx = replica2ThermoIxs[replicaIx];
			ThermodynamicState& thermoState = thermodynamicStates[thermoIx];
			const std::vector<int>& thermoWorldIxs = thermoState.updWorldIndexes();
			const std::vector<int>& distortOpts = thermoState.getDistortOptions();
			const Partitioning& wPart = thermoState.getNonequilPartitioning();

			// Update BAT map for all the replica's world
			updSubZMatrixBATsToAllWorlds(replicaIx);

			// Update simulation parameters
			setReplicasWorldsParameters(replicaIx, false, true, mixIndex);

			transferCoordsFromReplicaToWorld(replicaIx, 0);
			transferQStatistics(thermoIx, thermoWorldIxs[wPart.nofEquilibriumWorlds - 1], thermoWorldIxs[wPart.nofEquilibriumWorlds - 1]);

			// Simulate this replica
			RunReplicaWorldRange(replicaIx, 0, wPart.nofEquilibriumWorlds, false);

			replica.incrementNofSamples(1);
			thermoState.incrementNofSamples(1);

			// SimTK::State& state = worlds[thermoWorldIxs[wPart.N2_wCnt]].integ->updAdvancedState();
			// replica.calcZMatrixBAT( worlds[thermoWorldIxs[wPart.N2_wCnt]].getAtomsLocationsInGround( state ));
		}

		if(runType == RUN_TYPE::REBASONTOP) {
			runType = RUN_TYPE::REMC;

			for(int remcIx = 0; remcIx < 6; remcIx++){
				prepareExchangePairs(mixIndex, 0);
        		mixReplicas(mixIndex, 0);
				mixIndex++;
				//PrintNofAcceptedSwapsMatrix();
			}
			PrintNofAcceptedSwapsMatrix();
			runType = RUN_TYPE::REBASONTOP;
    	}

		// @@@@@@@@@@ LOOP THROUGH REPLICAS (NON-EQUILIBRIUM) ------------------- RUN A ----------->

		// Update work scale factors
		prepareExchangePairs(mixIndex, 0);
		updThermostatesQScaleFactors(mixIndex); // depends on exchange pairs, so needs to be updated after prepareExchangePairs !!!

		for (int replicaIx = 0; replicaIx < nofReplicas; replicaIx++){ // BY_REPLICA
			const int thermoIx = replica2ThermoIxs[replicaIx];
			Replica& replica = replicas[replicaIx];
			ThermodynamicState& thermoState = thermodynamicStates[thermoIx];
			const auto& thermoWorldIxs = thermoState.updWorldIndexes();
			const auto& distortOpts = thermoState.getDistortOptions();
			const size_t thermoNofWorlds = thermoWorldIxs.size();
			const Partitioning& wPart = thermoState.getNonequilPartitioning();

			// // Update BAT map for all the replica's world
			// updSubZMatrixBATsToAllWorlds(replicaIx);

			if(wPart.nofNonequilibriumWorlds){
				setReplicasWorldsParameters(replicaIx, false, true, mixIndex);
				transferCoordsFromReplicaToWorld(replicaIx, thermoWorldIxs[wPart.N1_wCnt]);
				transferQStatistics(thermoIx, thermoWorldIxs[wPart.N2_wCnt], thermoWorldIxs[wPart.N1_wCnt]);
				RunReplicaWorldRange(replicaIx, wPart.N1_wCnt, thermoNofWorlds, true);

				// Write log and DCD
				//writeReplicaLogAndDCD(mixi, replicaIx, printFreq);

				replica.incrementNofSamples(1);
				thermoState.incrementNofSamples(1);

				// SimTK::State& state = worlds.back().integ->updAdvancedState();
				// replica.calcZMatrixBAT( worlds.back().getAtomsLocationsInGround( state ));
			}

		} // _end_ loop through replicas (NON-EQUILIBRIUM) RUN A

		mixReplicas(mixIndex, 0);
		mixIndex++;
    	//PrintNofAcceptedSwapsMatrix();


		// @@@@@@@@@@ LOOP THROUGH REPLICAS (NON-EQUILIBRIUM) ------------------- RUN B ----------->

		// Update work scale factors
		prepareExchangePairs(mixIndex, 0);
		updThermostatesQScaleFactors(mixIndex); // depends on exchange pairs, so needs to be updated after prepareExchangePairs !!!

		for (int replicaIx = 0; replicaIx < nofReplicas; replicaIx++){ // BY_REPLICA
			const int thermoIx = replica2ThermoIxs[replicaIx];
			Replica& replica = replicas[replicaIx];
			ThermodynamicState& thermoState = thermodynamicStates[thermoIx];
			const auto& thermoWorldIxs = thermoState.updWorldIndexes();
			const auto& distortOpts = thermoState.getDistortOptions();
			const size_t thermoNofWorlds = thermoWorldIxs.size();
			assert((thermoWorldIxs.size() == distortOpts.size()));
			const Partitioning& wPart = thermoState.getNonequilPartitioning();

			// Update BAT map for all the replica's world
			updSubZMatrixBATsToAllWorlds(replicaIx);

			if(wPart.nofNonequilibriumWorlds){
				setReplicasWorldsParameters(replicaIx, false, true, mixIndex);
				transferCoordsFromReplicaToWorld(replicaIx, thermoWorldIxs[wPart.N1_wCnt]);
				transferQStatistics(thermoIx, thermoWorldIxs[wPart.N2_wCnt], thermoWorldIxs[wPart.N1_wCnt]);
				RunReplicaWorldRange(replicaIx, wPart.N1_wCnt, thermoNofWorlds, true);

				// // Write log and DCD
				// if ((mixIndex + 1) % printFreq == 0) {
				// 	writeLog(mixIndex, replicaIx);
				// 	REXLog(mixIndex, replicaIx);
				// 	std::cout << std::flush;

				// 	// Write DCDs
				// 	int whichDCD = replica2ThermoIxs[replicaIx];

				// 	// TODO cache these
				// 	std::vector<SimTK::Real> x (replicas[replicaIx].getX().size(), 0.0);
				// 	std::vector<SimTK::Real> y (replicas[replicaIx].getY().size(), 0.0);
				// 	std::vector<SimTK::Real> z (replicas[replicaIx].getZ().size(), 0.0);

				// 	for (const auto& atom : atoms) {
				// 		const std::size_t prmtopIndex = atom.identity.prmtopIndex;
				// 		x[prmtopIndex] = replicas[replicaIx].getX()[atom.identity.globalIndex] * 10;
				// 		y[prmtopIndex] = replicas[replicaIx].getY()[atom.identity.globalIndex] * 10;
				// 		z[prmtopIndex] = replicas[replicaIx].getZ()[atom.identity.globalIndex] * 10;
				// 	}
						
				// 	thermodynamicStates[whichDCD].writeDCD(x, y, z);
				// }

				replica.incrementNofSamples(1);
				thermoState.incrementNofSamples(1);

				// SimTK::State& state = worlds.back().integ->updAdvancedState();
				// replica.calcZMatrixBAT( worlds.back().getAtomsLocationsInGround( state ));
			}

		} // _end_ loop through replicas (NON-EQUILIBRIUM) RUN B

		mixReplicas(mixIndex, 0);
		mixIndex++;
    	// PrintNofAcceptedSwapsMatrix();

		if ((mixIndex + 1) % printFreq == 0) {
			for (int replicaIx = 0; replicaIx < nofReplicas; replicaIx++) {
				writeLog(mixIndex, replicaIx);
				REXLog(mixIndex, replicaIx);
			}
			std::cout << std::flush;

			// Write DCDs
			for (int replicaIx = 0; replicaIx < nofReplicas; replicaIx++) {
				
				// TODO cache these
				std::vector<SimTK::Real> x (replicas[replicaIx].getX().size(), 0.0);
				std::vector<SimTK::Real> y (replicas[replicaIx].getY().size(), 0.0);
				std::vector<SimTK::Real> z (replicas[replicaIx].getZ().size(), 0.0);
				
				for (const auto& atom : atoms) {
					const std::size_t prmtopIndex = atom.identity.prmtopIndex;
					x[prmtopIndex] = replicas[replicaIx].getX()[atom.identity.globalIndex] * 10;
					y[prmtopIndex] = replicas[replicaIx].getY()[atom.identity.globalIndex] * 10;
					z[prmtopIndex] = replicas[replicaIx].getZ()[atom.identity.globalIndex] * 10;
				}
				
				const int whichDCD = replica2ThermoIxs[replicaIx];
				thermodynamicStates[whichDCD].writeDCD(x, y, z);
			}
		}

		this->nofRounds++; 
	}

	//PrintNofAttemptedSwapsMatrix();
	PrintNofAcceptedSwapsMatrix();
	//PrintReplicaMaps();

	// foutU.close();
	// foutUDot.close();

	OPENMM::get().shutdown();
}

void Context::transferCoordsFromWorldToWorld(int sourceWorldIndex, int destinationWorldIndex) {
	const auto& atomTargetLocations = worlds[sourceWorldIndex].getAtomTargetLocationsCache();
	worlds[destinationWorldIndex].setAtomTargetLocationsToState(atomTargetLocations);
}

void Context::transferCoordsFromReplicaToWorld(int sourceReplicaIndex, int destinationWorldIndex) {
	const auto& atomTargetLocations = replicas[sourceReplicaIndex].getAtomsLocationsInGround();
	worlds[destinationWorldIndex].setAtomTargetLocationsToState(atomTargetLocations);
}

void Context::transferCoordsFromWorldToReplica(int sourceWorldIndex, int destinationReplicaIndex, bool intoWORK) {

	const World& world = worlds[sourceWorldIndex];
	Replica& replica = replicas[destinationReplicaIndex];

	const auto& atomTargetLocations = worlds[sourceWorldIndex].getAtomTargetLocationsCache();
	const auto potential = world.getSampler(0)->currentEnergy.potential;
	const auto fixman = world.getSampler(0)->currentEnergy.fixman;

	if(intoWORK){
		replica.updWORK() += world.getWork();  // TODO merge with Jacobians
		replica.upd_WORK_Jacobian() += world.getSampler(0)->getDistortJacobianDetLog();		
		replica.upd_WORK_AtomsLocationsInGround(atomTargetLocations);
		replica.set_WORK_PotentialEnergy_New(potential);
		replica.set_WORK_Fixman(world.getSampler(0)->currentEnergy.fixman);
		replica.set_WORK_ReferencePotentialEnergy_New(potential);
	}else{
		replica.setAtomsLocationsInGround(atomTargetLocations);
		replica.setPotentialEnergy(potential);
		replica.setFixman(fixman);
		replica.setReferencePotentialEnergy(potential);		
	}
}

void Context::initializeBinaryFile(const std::string &filename, uint32_t num_columns) {
    std::ofstream ofs(filename, std::ios::binary | std::ios::trunc); // overwrite file
    if (!ofs) {
        throw std::runtime_error("Cannot create file: " + filename);
    }

	std::cout << "Created binary file " << filename 
			  << " with " << num_columns << " columns." << std::endl;

    uint32_t num_rows = 0; // initially zero rows
    ofs.write(reinterpret_cast<const char *>(&num_columns), sizeof(num_columns));
    ofs.write(reinterpret_cast<const char *>(&num_rows), sizeof(num_rows));
    ofs.close();
}

// Append a row and update the row count in header
void Context::writeRowToBinaryFile(const std::string &filename, const std::vector<SimTK::Real> &row, bool has_acceptance, bool accepted) {
    std::fstream file(filename, std::ios::binary | std::ios::in | std::ios::out);
    if (!file) {
        throw std::runtime_error("File does not exist: " + filename);
    }

    uint32_t num_columns = 0, num_rows = 0;
    file.read(reinterpret_cast<char *>(&num_columns), sizeof(num_columns));
    file.read(reinterpret_cast<char *>(&num_rows), sizeof(num_rows));

    if (row.size() != (num_columns - 1)) { // -1 for acceptance column
		std::cout << "Row size: " << row.size() << " does not match expected number of columns: " << (num_columns - 1) << std::endl;
        throw std::runtime_error("Row size does not match expected number of columns.");
    }

    // Increment row count
    ++num_rows;
    file.seekp(sizeof(num_columns), std::ios::beg);
    file.write(reinterpret_cast<const char *>(&num_rows), sizeof(num_rows));

    // Seek to end to append
    file.seekp(0, std::ios::end);

    // Write row data
    file.write(reinterpret_cast<const char *>(row.data()), row.size() * sizeof(SimTK::Real));

    // Append acceptance column
    // SimTK::Real value_to_write = has_acceptance ? (accepted ? 1.0 : 0.0)
    //                                        : std::numeric_limits<SimTK::Real>::quiet_NaN();
	SimTK::Real value_to_write = (has_acceptance & accepted) ? 1.0 : 0.0;
    file.write(reinterpret_cast<const char *>(&value_to_write), sizeof(SimTK::Real));

    file.close();
}



/*!
 * <!--	zmatrixbat_ 
 * very inefficient so far -->
*/
void Context::setSubZmatrixBATStatsToSamplers(int thermoIx, int whichWorld)
{

	//scout("Context::setSubZmatrixBATStatsToSamplers") << eol;
	//PrintZMatrixTableAndBAT();

	// BAT stats containers to send to samplers
	std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATmeans;
	std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATdiffs;
	std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATvars;
	std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATvars_Alien;

	size_t zMatCnt = 0;

	for (const auto& zRow : zMatrixTable) {

		// Get Compound AtomIndex
		int childPositionInZMat = 0;
		RoboAtom& atom0 = atoms[zRow[childPositionInZMat]];
		SimTK::Compound::AtomIndex cAIx = atom0.identity.compoundAtomIndex;

		// Get statistics from thermodynamic state
		std::vector<SimTK::Real>& BATMeans = thermodynamicStates[thermoIx].getBATMeansRow(zMatCnt);
		std::vector<SimTK::Real>& BATDiffs = thermodynamicStates[thermoIx].getBATDiffsRow(zMatCnt);
		std::vector<SimTK::Real>& BATVars = thermodynamicStates[thermoIx].getBATVarsRow(zMatCnt);

		// Get statistics from thermostate exchange pair
		int alienThermoIx = getThermoPair(thermoIx);
		std::vector<SimTK::Real>& BATVars_Alien = thermodynamicStates[alienThermoIx].getBATVarsRow(zMatCnt);

		// Insert entry into sampler stats containers
		inBATmeans.insert({cAIx, BATMeans});
		inBATdiffs.insert({cAIx, BATDiffs});
		inBATvars.insert({cAIx, BATVars});
		inBATvars_Alien.insert({cAIx, BATVars_Alien});

		// Increment Z Matrix row
		zMatCnt++;

	} // ZMatrix row

	assert((inBATmeans.size() != 0) && 
		"Context BATmeans size is 0.");
	assert((inBATdiffs.size() != 0) && 
		"Context BATdiffs size is 0.");
	assert((inBATvars.size() != 0) && 
		"Context BATvars size is 0.");
	assert((inBATvars_Alien.size() != 0) && 
		"Context BATvars_Alien size is 0.");

	// scout("Context::setSubZmatrixBATStatsToSamplers") << eol;
	// for (const auto& [key, value] : inBATmeans) {
	// 	std::cout << "cAIx: " << key << " ";
	// 	std::cout << "BAT: ";
	// 	for (const auto& val : value) {
	// 		std::cout << val << " ";
	// 	}
	// 	std::cout << std::endl;
	// }

	// Set samplers BAT stats
	pHMC(worlds[whichWorld].updSampler(0))->setSubZMatrixBATStats(
		inBATmeans, inBATdiffs, inBATvars, inBATvars_Alien);

}

/*!
 * <!--  -->
*/
SimTK::Real Context::calcReplicaTransferedEnergy(int replicaIx)
{
	// Get thermoState corresponding to this replica
	int thisThermoStateIx = replica2ThermoIxs[replicaIx];

	// Get this world indexes from the corresponding thermoState
	std::vector<int> replicaWorldIxs = 
		thermodynamicStates[thisThermoStateIx].getWorldIndexes();

	// Get nof worlds in this replica
	size_t replicaNofWorlds = replicaWorldIxs.size();

	// Accumulate energy transfer here
	SimTK::Real deltaEnergy = 0;

	// Accumulate heat from equilibrium worlds and
	// work from perturbation kernels of nonequil worlds
	for(std::size_t worldIx = 0; worldIx < replicaNofWorlds; worldIx++){
			deltaEnergy += ( worlds[worldIx].getWorkOrHeat() );
	}

	return deltaEnergy;

}

/*!
 * <!-- Gather work contributions from all the worlds -->
*/
SimTK::Real Context::calcReplicaWork(int replicaIx)
{
	// Get thermoState corresponding to this replica
	int thisThermoStateIx = replica2ThermoIxs[replicaIx];

	// Get this world indexes from the corresponding thermoState
	std::vector<int> replicaWorldIxs = 
		thermodynamicStates[thisThermoStateIx].getWorldIndexes();

	// Get nof worlds in this replica
	size_t replicaNofWorlds = replicaWorldIxs.size();

	// Accumulate energy transfer here
	SimTK::Real Work = 0;

	// Accumulate heat from equilibrium worlds and
	// work from perturbation kernels of nonequil worlds
	for(std::size_t worldIx = 0; worldIx < replicaNofWorlds; worldIx++){
			Work += ( worlds[worldIx].getWork() );
	}

	return Work;
}


/*!
 * <!-- Print info about all the replicas and thermo states -->
*/
void Context::PrintReplicas()
{
	for (size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++){
		replicas[replicaIx].Print();
	}

	for(size_t thermoStateIx = 0;
	thermoStateIx < nofThermodynamicStates;
	thermoStateIx++){
		thermodynamicStates[thermoStateIx].Print();
	}
}

void Context::PrintNofAcceptedSwapsMatrix() {
	const std::size_t M = nofAcceptedSwapsMatrix.size();

	// Compute maximum width
	std::size_t maxWidth = 0;
	for (std::size_t i = 0; i < M; ++i) {
		for (std::size_t j = 0; j < M; ++j) {
			std::ostringstream oss;
			oss << nofAcceptedSwapsMatrix[i][j];
			maxWidth = std::max(maxWidth, oss.str().size());
		}
	}

	// Print aligned matrix
	for (std::size_t i = 0; i < M; ++i) {
		std::cout << "RSM";
		for (std::size_t j = 0; j < M; ++j) {
			std::cout << " " << std::setw(maxWidth) << nofAcceptedSwapsMatrix[i][j];
		}
		std::cout << '\n';
	}
}

void Context::PrintNofAttemptedSwapsMatrix(){

	size_t M = nofAttemptedSwapsMatrix.size();

	std::cout << "Number of attempted swaps matrix:\n";
	for(size_t i = 0; i < M; i++){
		for(size_t j = 0; j < M; j++){
			std::cout << nofAttemptedSwapsMatrix[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void Context::writePdbs(int someIndex, int thermodynamicStateIx)
{
	//for(int world_i = 0; world_i < this->nofWorlds; world_i++){

		// Update bAtomList in Topology
		const SimTK::State& pdbState =
				worlds[worldIndexes.front()].updIntegrator().updAdvancedState();
			worlds[worldIndexes.front()].updateAtomListsFromSimbody(pdbState);

		// const SimTK::State& pdbState =
		// 	worlds[world_i].updIntegrator().updAdvancedState();
		// worlds[world_i].updateAtomListsFromCompound(pdbState);

		// Write
		for(int mol_i = 0; mol_i < numMolecules; mol_i++){
			topologies[mol_i].writeAtomListPdb(
				outputDir,
				"/pdbs/sb."
					+ pdbPrefix + "." + std::to_string(mol_i) + "."
					+ "s" + std::to_string(thermodynamicStateIx) + ".",
					//+ "w" + std::to_string(world_i) + ".",
					".pdb",
					10,
					someIndex);
		}

	//}

}

void Context::randomizeWorldIndexes()
{
	// Random int for random world order
	std::uniform_int_distribution<std::size_t>
		randWorldDistrib(1, nofWorlds-1); // TODO between 1 and nOfWorlds-1?

	if(getNofWorlds() >= 3){

		// Swap world indeces between vector position 2 and random
		auto randVecPos = randWorldDistrib(randomEngine);
		//std::cout << "Swapping position 1 with "
		//	<< randVecPos << std::endl;

		auto secondWorldIx = worldIndexes[1];
		auto randWorldIx = worldIndexes[randVecPos];

		worldIndexes[1] = randWorldIx;
		worldIndexes[randVecPos] = secondWorldIx;

	}
}

/*!
 * <!--  -->
*/
void Context::passThroughBonds_template(int whichWorld)
{
	// Iterate molecules
	for(const auto& topology : topologies) {
		
		// Iterate molecule's bonds
		for(const auto& bond: topology.getBonds()) {

			RoboAtom& childAtom  = atoms[bond.globalIndices[1]];
			RoboAtom& parentAtom = atoms[bond.globalIndices[0]];

			// Get Compound atom indexes
			SimTK::Compound::AtomIndex child_cAIx = childAtom.identity.compoundAtomIndex;
			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.identity.compoundAtomIndex;
			
			// Get DuMM atom indezes
			SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
			SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);
		}
	}
}

// Print to log and write pdbs
void Context::RunLog(int round)
{
	// Write energy and geometric features to logfile
	if(printFreq || pdbRestartFreq){
		if( !(round % getPrintFreq()) ){
			for(auto wIx: worldIndexes){
				PrintToLog(0, wIx, 0);
			}
		}

		// Write pdb
		if( pdbRestartFreq != 0){
			if((round % pdbRestartFreq) == 0){
				writePdbs(round);
			}
		}
	}
}


/** Analysis related functions **/
void Context::addDistance(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t aIx1, std::size_t aIx2)
{
	distanceIxs.push_back({ whichWorld, whichCompound, aIx1, aIx2 });
}

// Get distances
void Context::addDistances(const std::vector<std::size_t>& distanceIx)
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for(unsigned int ai = 0; ai < distanceIx.size() / 2; ai++){
			addDistance(worldIx, 0,
				distanceIx[2*ai + 0], distanceIx[2*ai + 1]);
		}
	}
}

void Context::addAngle(std::size_t whichWorld, std::size_t whichCompound,
	std::size_t aIx1, std::size_t aIx2, std::size_t aIx3)
{
	angleIxs.push_back({ whichWorld, whichCompound, aIx1, aIx2, aIx3 });
}

// Get dihedrals. TODO : only adds to the first Topology
void Context::addAngles(const std::vector<std::size_t>& angleIx)
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for(unsigned int ai = 0; ai < angleIx.size() / 3; ai++){
			addAngle(worldIx, 0,
				angleIx[3*ai + 0], angleIx[3*ai + 1],
				angleIx[3*ai + 2]);
		}
	}
}

void Context::addDihedral(std::size_t whichWorld, std::size_t whichCompound,
	std::size_t aIx1, std::size_t aIx2, std::size_t aIx3, std::size_t aIx4)
{
	dihedralIxs.push_back({ whichWorld, whichCompound, aIx1, aIx2, aIx3, aIx4 });
}

// Get dihedrals. TODO : only adds to the first Topology
void Context::addDihedrals(const std::vector<std::size_t>& dihedralIx)
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for(unsigned int ai = 0; ai < dihedralIx.size() / 4; ai++){
			addDihedral(worldIx, 0,
				dihedralIx[4*ai + 0], dihedralIx[4*ai + 1],
				dihedralIx[4*ai + 2], dihedralIx[4*ai + 3]);

			// std::cout <<"Context::addDihedrals"
			// 	<<" | "<< dihedralIx[4*ai + 0] <<" "<< dihedralIx[4*ai + 1] <<" "<< dihedralIx[4*ai + 2] <<" "<< dihedralIx[4*ai + 3]
			// 	<<" | "<< atoms[dihedralIx[4*ai + 0]].getNumber() <<" "<< atoms[dihedralIx[4*ai + 1]].getNumber() <<" "<< atoms[dihedralIx[4*ai + 2]].getNumber() <<" "<< atoms[dihedralIx[4*ai + 3]].getNumber()
			// 	<<" | "<< atoms[dihedralIx[4*ai + 0]].getInName() <<" "<< atoms[dihedralIx[4*ai + 1]].getInName() <<" "<< atoms[dihedralIx[4*ai + 2]].getInName() <<" "<< atoms[dihedralIx[4*ai + 3]].getInName()
			// 	<< std::endl;
		
		}
	}
}

// --- Printing functions --

// Print energy information
void Context::PrintSamplerDataToLog(std::size_t whichWorld, std::size_t whichSampler)
{
	SimTK::State& currentAdvancedState = worlds[whichWorld].updIntegrator().updAdvancedState();

	const auto NU = currentAdvancedState.getNU();
	const auto acceptedSteps = pHMC((worlds[whichWorld].samplers[0]))->acceptedSteps;
	const auto pe_o = pHMC((worlds[whichWorld].samplers[0]))->getPreviousEnergy().potential;
	const auto ke_o = pHMC((worlds[whichWorld].samplers[0]))->getPreviousEnergy().kinetic;
	const auto ke_n = pHMC((worlds[whichWorld].samplers[0]))->getCurrentEnergy().kinetic;
	const auto fix_o = pHMC((worlds[whichWorld].samplers[0]))->getPreviousEnergy().fixman;
	const auto fix_n = pHMC((worlds[whichWorld].samplers[0]))->getCurrentEnergy().fixman;

	// Write to a file instead of stdout
	logFile
		<< std::fixed << std::setprecision(3)
		<< worlds[whichWorld].updSampler(whichSampler)->getTemperature() << " "
		<< std::fixed << std::setprecision(0)
		<< whichWorld << " "
		<< NU << " "
		<< acceptedSteps << " "
		<< std::fixed << std::setprecision(2) << pe_o << " "
		<< ke_o << " "
		<< ke_n << " "
		<< fix_o << " "
		<< fix_n << " ";
}

// Print geometric parameters during simulation
void Context::PrintGeometry(SetupReader& setupReader, std::size_t whichWorld)
{
	if(setupReader.get("GEOMETRY")[0] == "TRUE"){
		// Get distances indeces
		std::vector<int> distanceIx(setupReader.get("DISTANCE").size());
		for(unsigned int i = 0; i < setupReader.get("DISTANCE").size(); i++){
			distanceIx.emplace_back(atoi(setupReader.get("DISTANCE")[i].c_str()));
		}

		// Get distances
		for(size_t ai = 0; ai < (setupReader.get("DISTANCE").size() / 2); ai++){
			/*
			std::cout << std::setprecision(4)
			<< Distance(whichWorld, 0, 0,
				distanceIx[2*ai + 0], distanceIx[2*ai + 1]) << " ";
			*/
			printf("%.2f ", Distance(whichWorld, 0, 0,
				 distanceIx[2*ai + 0], distanceIx[2*ai + 1]));

		}

		// Get dihedrals indeces
		std::vector<int> dihedralIx(setupReader.get("DIHEDRAL").size());
		for(unsigned int i = 0; i < setupReader.get("DIHEDRAL").size(); i++){
			dihedralIx.emplace_back(atoi(setupReader.get("DIHEDRAL")[i].c_str()));
		}
		// Get dihedrals
		for(size_t ai = 0; ai < (setupReader.get("DIHEDRAL").size() / 4); ai++){
			/*
			std::cout << std::setprecision(4)
			<< Dihedral(whichWorld, 0, 0,
				dihedralIx[4*ai + 0], dihedralIx[4*ai + 1],
				dihedralIx[4*ai + 2], dihedralIx[4*ai + 3]) << " ";
			*/

			printf("%.2f ", Dihedral(whichWorld, 0, 0,
				dihedralIx[4*ai + 0], dihedralIx[4*ai + 1],
				dihedralIx[4*ai + 2], dihedralIx[4*ai + 3]));

		}
		//std::cout << std::endl;
		printf("\n");
	}else{
		//std::cout << std::endl;
		printf("\n");
	}
}

void Context::PrintGeometryToLog(std::size_t whichWorld, std::size_t whichSampler)
{
	PrintDistancesToLog(whichWorld, whichSampler);
	PrintAnglesToLog(whichWorld, whichSampler);
	PrintDihedralsQsToLog(whichWorld, whichSampler);
}

void Context::PrintDistancesToLog(std::size_t whichWorld, std::size_t whichSampler)
{
	for (const auto& distanceIx : distanceIxs) {
		if( distanceIx[0] == whichWorld ) {
			logFile << std::fixed << std::setprecision(3) << Distance(distanceIx[0], distanceIx[1], 0, distanceIx[2], distanceIx[3]) << " ";
		}
	}
}

void Context::PrintAnglesToLog(std::size_t whichWorld, std::size_t whichSampler)
{
	for (const auto& angleIx : angleIxs){
		if( angleIx[0] == whichWorld ) {
			logFile << std::fixed << std::setprecision(3) << Roboangle(angleIx[0], angleIx[1], 0, angleIx[2], angleIx[3], angleIx[4]) << " ";
		}
	}
}

void Context::PrintDihedralsToLog(std::size_t whichWorld, std::size_t whichSampler)
{
	for (const auto& dihedralIx : dihedralIxs){
		if( dihedralIx[0] == whichWorld ) {
			logFile << std::fixed << std::setprecision(3) << Dihedral(dihedralIx[0], dihedralIx[1], 0, dihedralIx[2], dihedralIx[3], dihedralIx[4], dihedralIx[5]) << " ";
		}
	}
}

void Context::PrintDihedralsQsToLog(std::size_t whichWorld, std::size_t whichSampler)
{
	for (const auto& dihedralIx : dihedralIxs){
		if( dihedralIx[0] == whichWorld ) {

			// std::cout << "Context::PrintDihedralsQs w c s a1 a2 a3 a4 ="
			//     << " " << dihedralIx[0] << " " << dihedralIx[1] << " " << 0
			//     << " " << (worlds[whichWorld].getTopology(0)).getAtomName( SimTK::Compound::AtomIndex(dihedralIx[2]))
			//     << " " << (worlds[whichWorld].getTopology(0)).getAtomName( SimTK::Compound::AtomIndex(dihedralIx[3]))
			//     << " " << (worlds[whichWorld].getTopology(0)).getAtomName( SimTK::Compound::AtomIndex(dihedralIx[4]))
			//     << " " << (worlds[whichWorld].getTopology(0)).getAtomName( SimTK::Compound::AtomIndex(dihedralIx[5]))
			//     << std::endl;

			logFile << std::fixed << std::setprecision(3) << Dihedral(
				dihedralIx[0], dihedralIx[1], 0,
				dihedralIx[2], dihedralIx[3],
				dihedralIx[4], dihedralIx[5]) << " ";

			// const Topology& topology = worlds[whichWorld].getTopology(dihedralIx[1]);
			// SimTK::State& currentAdvancedState = worlds[whichWorld].updIntegrator().updAdvancedState();
			// SimTK::MobilizedBodyIndex mbx3 = topology.getAtomMobilizedBodyIndex(
			//     SimTK::Compound::AtomIndex(dihedralIx[4]) );
			// SimTK::MobilizedBody::Pin& mobod3 = (SimTK::MobilizedBody::Pin&) (worlds[whichWorld].uptMatterSubsystem().updMobilizedBody(mbx3));

			// //std::cout << mbx3 << std::endl ;
			// //std::cout << currentAdvancedState.getQ() << std::endl;
			// //fprintf(logFile, "%.3f ", currentAdvancedState.getQ()[mbx3] );
			// fprintf(logFile, "%.3f ", mobod3.getQ(currentAdvancedState) );

		}
	}
}

void Context::PrintFreeE2EDist(std::size_t whichWorld, int whichCompound)
{
	// const Topology& topology = worlds[whichWorld].getTopology(whichCompound);
	// SimTK::State& currentAdvancedState = worlds[whichWorld].updIntegrator().updAdvancedState();

	// for (const auto& distanceIx : distanceIxs) {
	// 	if( distanceIx[0] == whichWorld ) {

	// 		logFile << std::fixed << std::setprecision(3) << 
	// 			Distance(distanceIx[0], distanceIx[1], 0, distanceIx[2], distanceIx[3]) << std::endl;

	// 		SimTK::MobilizedBodyIndex mbx1 = topology.getAtomMobilizedBodyIndex(
	// 			SimTK::Compound::AtomIndex(distanceIx[2]) );
	// 		SimTK::MobilizedBodyIndex mbx2 = topology.getAtomMobilizedBodyIndex(
	// 			SimTK::Compound::AtomIndex(distanceIx[3]) );
	// 		SimTK::MobilizedBody& mobod1 = worlds[whichWorld].uptMatterSubsystem().updMobilizedBody(mbx1);
	// 		SimTK::MobilizedBody& mobod2 = worlds[whichWorld].uptMatterSubsystem().updMobilizedBody(mbx2);
	// 		SimTK::Transform X_PF1 = mobod1.getInboardFrame(currentAdvancedState);
	// 		SimTK::Transform X_PF2 = mobod2.getInboardFrame(currentAdvancedState);
	// 		//SimTK::Transform X_BM1 = mobod1.getOutboardFrame(currentAdvancedState);
	// 		//SimTK::Transform X_BM2 = mobod2.getOutboardFrame(currentAdvancedState);
	// 		SimTK::Transform X_FM1 = mobod1.getMobilizerTransform(currentAdvancedState);
	// 		SimTK::Transform X_FM2 = mobod2.getMobilizerTransform(currentAdvancedState);

	// 		SimTK::Transform deltaX_PF = X_PF2.p() - X_PF1.p();
	// 		logFile << std::fixed << std::setprecision(3) << ((-1 * X_FM1.p()) + deltaX_PF.p() + X_FM2.p()).norm() << std::endl;

	// 		//std::cout << "X_PF1:" << std::endl << X_PF1 << std::endl;
	// 		//std::cout << "X_FM1:" << std::endl << X_FM1 << std::endl;
	// 		//std::cout << "X_BM1:" << std::endl << X_BM1 << std::endl;
	// 		//std::cout << "X_PF2:" << std::endl << X_PF2 << std::endl;
	// 		//std::cout << "X_FM2:" << std::endl << X_FM2 << std::endl;
	// 		//std::cout << "X_BM2:" << std::endl << X_BM2 << std::endl;

	// 		//SimTK::Vec3 a1pos = X_PF1.R() * X_FM1.p();
	// 		//SimTK::Vec3 a2pos = X_PF2.R() * X_FM2.p();
	// 		//fprintf(logFile, "%.3f ",
	// 		//    (a1pos - a2pos).norm() );

	// 	}
	// }


}

void Context::PrintToLog(std::size_t whichReplica,
	std::size_t whichWorld, std::size_t whichSampler)
{
	logFile << whichReplica << " ";

	PrintSamplerDataToLog(whichWorld, whichSampler);

	PrintGeometryToLog(whichWorld, whichSampler);

	logFile << "\n";
}

// Write intial pdb for reference
// TODO: what's the deal with mc_step
void Context::writeInitialPdb()
{

	// - we need this to get compound atoms
	int currentWorldIx = worldIndexes.front();
	SimTK::State& advancedState = worlds[currentWorldIx].updIntegrator().updAdvancedState();

	constexpr int mc_step = -1;

	// Pass compounds to the new world
	passTopologiesToNewWorld(currentWorldIx);

	// 
	worlds[currentWorldIx].updateAtomListsFromSimbody(advancedState);
	std::cout << "Writing pdb initial" << mc_step << ".pdb" << std::endl;

	// 
	for(unsigned int mol_i = 0; mol_i < topologies.size(); mol_i++){
		topologies[mol_i].writeAtomListPdb(getOutputDir(),
		"/pdbs/sb." + getPdbPrefix() + ".", ".pdb", 10, mc_step);
	}

}

// Write final pdb for reference
void Context::writeFinalPdb()
{

	// Update bAtomList in Topology
	const SimTK::State& pdbState =
			worlds[worldIndexes.front()].updIntegrator().updAdvancedState();
		worlds[worldIndexes.front()].updateAtomListsFromSimbody(pdbState);

	// Write
	for(unsigned int mol_i = 0; mol_i < numMolecules; mol_i++){
		topologies[mol_i].writeAtomListPdb(
			getOutputDir(),
			"/pdbs/final."
			+ getPdbPrefix() + std::to_string(mol_i) + ".",
			".pdb",
			10,
			getRequiredNofRounds());
	}

}

// Get / set pdb files writing frequency
int Context::getPdbRestartFreq()
{
	return this->pdbRestartFreq;
}

//
void Context::setPdbRestartFreq(int argFreq)
{
	this->pdbRestartFreq = argFreq;
}

const std::string& Context::getRestartDir() const
{
	return this->restartDir;
}

void Context::setRestartDir(const std::string& argRestartDir)
{
	this->restartDir = argRestartDir;
}

// Get / set printing frequency
int Context::getPrintFreq()
{
	return this->printFreq;
}

// 
void Context::setPrintFreq(int argFreq)
{
	this->printFreq = argFreq;
}

std::string Context::getOutputDir()
{
	return this->outputDir;
}

void Context::setOutputDir(std::string arg)
{
	this->outputDir = arg;
}

void Context::setPdbPrefix(const std::string& argPdbPrefix)
{
	this->pdbPrefix = argPdbPrefix;
}

std::string Context::getPdbPrefix()
{
	return this->pdbPrefix;
}

SimTK::Real Context::Roboangle(std::size_t whichWorld,
	std::size_t whichCompound, std::size_t whichSampler,
	int a1, int a2, int a3)
{

	SimTK::State& state = worlds[whichWorld].updIntegrator().updAdvancedState();

	Topology& topology = worlds[whichWorld].updTopology(whichCompound);

	const SimTK::DuMMForceFieldSubsystem& dumm = worlds[whichWorld].getForceField();
	const SimTK::SimbodyMatterSubsystem& matter = worlds[whichWorld].getMatterSubsystem();

	int cAIx_1 = atoms[a1].identity.compoundAtomIndex;
	int cAIx_2 = atoms[a2].identity.compoundAtomIndex;
	int cAIx_3 = atoms[a3].identity.compoundAtomIndex;

	SimTK::Vec3 a1pos, a2pos, a3pos, a4pos;
	a1pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_1)),
		dumm, matter, state);
	a2pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_2)),
		dumm, matter, state);
	a3pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_3)),
		dumm, matter, state);

	return bAngle(a1pos, a2pos, a3pos);

}


SimTK::Real Context::Dihedral(std::size_t whichWorld,
	std::size_t whichCompound, std::size_t whichSampler,
	int a1, int a2, int a3, int a4)
{

	SimTK::State& state = worlds[whichWorld].updIntegrator().updAdvancedState();

	Topology& topology = worlds[whichWorld].updTopology(whichCompound);

	SimTK::DuMMForceFieldSubsystem& dumm = worlds[whichWorld].updForceField();

	SimTK::SimbodyMatterSubsystem& matter = worlds[whichWorld].updMatterSubsystem();

	SimTK::Vec3 a1pos, a2pos, a3pos, a4pos;

	// std::cout <<"Context::Dihedral" 
	// 	<<" | "<< whichWorld <<" "<< whichCompound <<" "<< whichSampler <<" "<< a1 <<" "<< a2 <<" "<< a3 <<" "<< a4
	// 	<<" | "<< atoms[a1].getInName() <<" "<< atoms[a2].getInName() <<" "<< atoms[a3].getInName() <<" "<< atoms[a4].getInName()
	// << std::endl;

	int cAIx_1 = atoms[a1].identity.compoundAtomIndex;
	int cAIx_2 = atoms[a2].identity.compoundAtomIndex;
	int cAIx_3 = atoms[a3].identity.compoundAtomIndex;
	int cAIx_4 = atoms[a4].identity.compoundAtomIndex;

	a1pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_1)),
		dumm, matter, state);
	a2pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_2)),
		dumm, matter, state);
	a3pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_3)),
		dumm, matter, state);
	a4pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_4)),
		dumm, matter, state);

	return bDihedral(a1pos, a2pos, a3pos, a4pos);

}

SimTK::Real Context::Distance(
	std::size_t whichWorld, std::size_t whichCompound, std::size_t whichSampler,
	int a1, int a2)
{

	SimTK::State& state = worlds[whichWorld].updIntegrator().updAdvancedState();

	Topology& topology = worlds[whichWorld].updTopology(whichCompound);

	SimTK::DuMMForceFieldSubsystem& dumm = worlds[whichWorld].updForceField();

	SimTK::SimbodyMatterSubsystem& matter = worlds[whichWorld].updMatterSubsystem();

	SimTK::Vec3 a1pos, a2pos;

	int cAIx_1 = atoms[a1].identity.compoundAtomIndex;
	int cAIx_2 = atoms[a2].identity.compoundAtomIndex;

	a1pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_1)),
		dumm, matter, state);
	a2pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_2)),
		dumm, matter, state);

	return (a1pos - a2pos).norm();

}

// Realize Topology Stage for all the Worlds
void Context::realizeTopology() {
	for(auto& world : worlds) {
		world.updCompoundSystem().realizeTopology();
	}

}

/**
 *  Pass compounds to the new world
 */
void Context::areAllDuMMsTheSame(void)
{

	// std::cout << "=== areAllDuMMsTheSame ===\n";

	// // All molecules
	// for(std::size_t molIx = 0; molIx < numMolecules; molIx++){

	// 	// All Compound atoms
	// 	for(std::size_t k = 0; k < topologies[molIx].getNumAtoms(); k++){

	// 		// Get atom index
	// 		SimTK::Compound::AtomIndex aIx =
	// 			(topologies[molIx].subAtomList[k]).identity.compoundAtomIndex;

	// 		std::cout << aIx << " ";

	// 		// All worlds
	// 		for(int worldIx = 0; worldIx < getNofWorlds(); worldIx++){

	// 			// Get mobod
	// 			SimTK::MobilizedBodyIndex mbx =
	// 				topologies[molIx].getAtomMobilizedBodyIndexThroughDumm(aIx,
	// 				*(worlds[worldIx].forceField) );

	// 			SimTK::DuMM::AtomIndex dAIx = topologies[molIx].getDuMMAtomIndex(aIx);

	// 			std::cout << dAIx << " ";	
	// 		}
	// 		std::cout << std::endl;

	// 	}

	// 	// TODO Restante DANGER
    //     //c.setTopLevelTransform(compoundTransform * c.getTopLevelTransform());
	// }
}




// Teodor's membrane

/** Implicit membrane mimicked by half-space contacts */
void Context::addContactImplicitMembrane(const float memZWidth, const SetupReader& setupReader){
	// // Before adding the membrane, we add the contacts and join them
	// // to the appropiate Contact Cliques

	// // Each of the flags are formatted as such:
	// // CONTACTS_X  int1 int2 , int3 int4 int5
	// // X is one of {0,1,2,3}, the ints are atom indices (0-Based)
	// // and the comma separates the topologies. In the example given
	// // the contacts are set on atoms int1, int2 for topology 0,
	// // and on atoms int3, int4 and int5 for topology 1.
	// // If the user wishes to skip a topology, then they'd 
	// // input "-1" as the only atom index.

	// std::vector<std::vector<std::vector<int>>> cliqueAtomIxs;

	
	// for (int contactCliqueIx = 0; contactCliqueIx < 4; contactCliqueIx++){

	// 	// Empty vector of prmtop atom indexes
	// 	cliqueAtomIxs.push_back({});
	// 	cliqueAtomIxs[contactCliqueIx].push_back({});

	// 	// Get values for this contactCliqueIx
	// 	std::string contactClique_key = "CONTACTS_";
	// 	contactClique_key.append( std::to_string(contactCliqueIx) );
	// 	const std::vector<std::string>& contactClique_vals = setupReader.get(contactClique_key);
		
	// 	int cur_topology = 0;
	// 	if (contactClique_vals.size() > 2) {

	// 		// Get atom indexes for this clique
	// 		for (const auto& value : contactClique_vals){
				
	// 			if (value == ",") { //TODO: This does not account for 'int1,'. Fix this.
	// 				cliqueAtomIxs[contactCliqueIx].push_back({});
	// 				cur_topology++;
	// 			}
	// 			else {
	// 				cliqueAtomIxs[contactCliqueIx][cur_topology].push_back(std::stoi(value));
	// 			}
	// 		}

	// 		// Check
	// 		if(cur_topology != topologies.size()){
	// 			std::cout << "[WARNING] " 
	// 				<< "Number of topologies in CONTACT_ keys don't match the actual number of topologies\n";
	// 		}

	// 		// Add contact atom indexes for all worlds
	// 		for(unsigned int worldIx = 0; worldIx < getNofWorlds(); worldIx++){
	// 			for (int topologyIx = 0;topologyIx < cliqueAtomIxs[contactCliqueIx].size(); topologyIx++){
					
	// 				worlds[worldIx].addContacts(
	// 						cliqueAtomIxs[contactCliqueIx][topologyIx],
	// 						topologyIx,
	// 						SimTK::ContactCliqueId(contactCliqueIx));
	// 			}
	// 		}
	// 	}
	// }

	// // Add membrane to all worlds.
	// for(unsigned int worldIx = 0; worldIx < getNofWorlds(); worldIx++){
	// 	worlds[worldIx].addMembrane(memZWidth);
	// }

	// // Print
	// std::cout << "\n########## MEMBRANE STATS ##########\n";
	// std::cout << "Atom cliques are: \n";
	// for (int contactClique=0; contactClique<4; contactClique++){
	// 	for (const auto& topologyIx : cliqueAtomIxs[contactClique]) {
	// 	for (int atomIx : topologyIx) {
	// 		std::cout << atomIx << " ";
	// 	}
	// 	std::cout << " / ";
	// }
	// std::cout << std::endl;
	// }
	// std::cout << "########## MEMBRANE STATS ##########\n\n";

	// // TODO: Do we need this here (looks like World's buissiness)
	// realizeTopology();
}

void Context::setNumThreads(int threads) {
	if (threads < 0) {
		std::cerr << "Invalid number of threads (negative value). Default number (0) of threads will be used." << std::endl;
		numThreads = 0;
	} else {
		numThreads = threads;
	}
}

void Context::setNonbonded(NonbondedMethod method, SimTK::Real cutoffInNm) {

	SimTK_ASSERT_ALWAYS(cutoffInNm >= 0, "Context::setNonbonded: Cutoff distance cannot be negative.");

	nonbondedMethod = method;
	nonbondedCutoffInNm = cutoffInNm;
}

void Context::setGBSAOptions(bool useGBSAOBC2, SimTK::Real solventDielectric, SimTK::Real soluteDielectric) {
	this->useGBSAOBC2 = useGBSAOBC2;
	this->solventDielectric = solventDielectric;
	this->soluteDielectric = soluteDielectric;

	if (useGBSAOBC2) {
		gbsaGlobalScaleFactor = 1.0;
		// gbsaGlobalScaleFactor = 1.0 / 1.2; // GBSA OBC2 scaling factor
	} else {
		gbsaGlobalScaleFactor = 0.0; // No scaling
	}
}

void Context::setForceFieldScaleFactors(SimTK::Real globalScaleFactor) {
	useAmberForceFieldScaleFactors = false;

	if (globalScaleFactor < 0 && globalScaleFactor > 1) {
		std::cerr << "Invalid force field scale factor (valid range is 0 to 1). Default value of 1.0 will be used instead." << std::endl;
		globalForceFieldScaleFactor = 1.0;
	} else {
		globalForceFieldScaleFactor = globalScaleFactor;
	}
}



// ===========================================================================
// ===========================================================================
// ZMatrix BAT
// ===========================================================================
// ===========================================================================

/*!
 * <!--	zmatrixbat_ -->
*/
void Context::addZMatrixTableRow(const std::vector<int>& newRow) {
	// Add bounds checking if needed
	zMatrixTable.push_back(newRow);
}

/*!
 * <!--	zmatrixbat_ -->
*/
int Context::getZMatrixTableEntry(int rowIndex, int colIndex) const {
	// Add bounds checking if needed
	return zMatrixTable[rowIndex][colIndex];
}

/*!
 * <!--	zmatrixbat_ -->
*/
void Context::setZMatrixTableEntry(int rowIndex, int colIndex, int value) {
	// Add bounds checking if needed
	zMatrixTable[rowIndex][colIndex] = value;
}

/*!
 * <!--	zmatrixbat_ -->
*/
void Context::PrintZMatrixTable() const {
	for (const auto& row : zMatrixTable) {
		scout("ZMatrinxTableEntry: ");
		for (int value : row) {
			std::cout << std::setw(6) << value <<" "; 
		}
		std::cout << std::endl; 
	}
}

/*!
 * <!--	zmatrixbat_ -->
*/
void Context::setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value) {

	// Set the value at the specified position
	zMatrixBAT[rowIndex][colIndex] = value;
}

/*!
 * <!--	zmatrixbat_ -->
*/
const std::vector<SimTK::Real>& Context::getZMatrixBATRow(size_t rowIndex) const {

	assert(rowIndex < zMatrixBAT.size());

	// Check if the indices are within bounds
	return zMatrixBAT[rowIndex];

}

/*!
 * <!--	zmatrixbat_ -->
*/
std::vector<SimTK::Real>& Context::updZMatrixBATRow(size_t rowIndex) {

	assert(rowIndex < zMatrixBAT.size());

	// Check if the indices are within bounds
	return zMatrixBAT[rowIndex];

}


/*!
 * <!-- zmatrixbat_ Get Z-matrix indexes table -->
*/
void
Context::calcZMatrixTable(void)
{
	// // Iterate molecules
	// for(size_t topoIx = 0; topoIx < numMolecules; topoIx++){

	// 	// Get molecule and it's bonds 
	// 	const std::vector<BOND>& BONDS = orderedBonds[topoIx];

	// 	if (BONDS.size() == 0) {
	// 		continue;
	// 	}

	// 	std::cout << "BONDS " << BONDS[0].first << " " << BONDS[0].second << std::endl;
	// 	addZMatrixTableRow(std::vector<int> {
	// 		BONDS[0].first,
	// 		BONDS[0].second,
	// 		-1, -2
	// 	});

	// 	addZMatrixTableRow(std::vector<int> {
	// 		BONDS[1].first,
	// 		BONDS[1].second,
	// 		BONDS[internCoords.findBondByFirst(topoIx, BONDS[1].second)].second,
	// 		-1
	// 	});

	// 	// Iterate molecule's bonds
	// 	for(size_t BOIx = 2; BOIx < BONDS.size(); BOIx++){

	// 		// Get current bond
	// 		const BOND& currBOND = BONDS[BOIx];

	// 		// Get bond's atoms
	// 		int childNo = currBOND.first;
	// 		int parentNo = currBOND.second;
	// 		int gparentNo = BONDS[internCoords.findBondByFirst(topoIx, parentNo)].second;
			
	// 		// Not all gparents have ggparents
	// 		int ggparentNo = -3;
	// 		int ggparentBONDIx = internCoords.findBondByFirst(topoIx, gparentNo);
	// 		if(ggparentBONDIx >= 0){
	// 			ggparentNo = BONDS[internCoords.findBondByFirst(topoIx, gparentNo)].second;
	// 		}

	// 		addZMatrixTableRow(std::vector<int> {childNo, parentNo, gparentNo, ggparentNo});
	// 	} // every bond
	// } // every molecule	
}


/*!
 * <!-- zmatrixbat_ Allocate Z Matrix BAT -->
*/
void Context::reallocZMatrixBAT(void){

	zMatrixBAT.resize(zMatrixTable.size());
	for (auto& row : zMatrixBAT) {
		row.resize(3, SimTK::NaN);
	}

}

/*!
 * <!-- zmatrixbat_ Calculate Z-matrix -->
*/
void
Context::calcZMatrixBAT(
	int wIx,
	const std::vector< std::vector<
		std::pair <RoboAtom *, SimTK::Vec3 > > >&
		otherWorldsAtomsLocations)
{

	// // Iterate molecules
	// int allCnt = 0;

	// int topoIx = 0;

	// // Get locations of this molecule
	// std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
	// for (int i = 0; i < topologies.size(); i++) {
	// 	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> temp;
	// 	std::cout << "Context::calcZMatrixBAT" << std::endl;
	// 	worlds[wIx].extractAtomTargets(i, otherWorldsAtomsLocations, temp);
	// 	atomTargets.insert(temp.begin(), temp.end());
	// }

	// int rowCnt = 0;
	// SimTK::Real bondLength, bondBend, bondTorsion;
	// for (const auto& row : zMatrixTable) {
		
	// 	bondLength = SimTK::NaN;
	// 	bondBend = SimTK::NaN;
	// 	bondTorsion = SimTK::NaN;

	// 	SimTK::Compound::AtomIndex a0_cAIx, a1_cAIx;
		
	// 	// Calculate bond length
	// 	a0_cAIx = atoms[row[0]].identity.compoundAtomIndex;
	// 	a1_cAIx = atoms[row[1]].identity.compoundAtomIndex;
	// 	SimTK::Vec3 a0loc = findAtomTarget(atomTargets, a0_cAIx);
	// 	SimTK::Vec3 a1loc = findAtomTarget(atomTargets, a1_cAIx);

	// 	SimTK::Vec3 v_a0a1 = a0loc - a1loc;
	// 	SimTK::Real bondLength = std::sqrt(SimTK::dot(v_a0a1, v_a0a1));

	// 	if(row[2] >= 0){

	// 		SimTK::Compound::AtomIndex a2_cAIx;
	// 		a2_cAIx = atoms[row[2]].identity.compoundAtomIndex;
	// 		SimTK::Vec3 a2loc = findAtomTarget(atomTargets, a2_cAIx);

	// 		// Calculate angle
	// 		UnitVec3 v1(v_a0a1);
	// 		UnitVec3 v2(a2loc - a1loc);

	// 		Real dotProduct = SimTK::dot(v1, v2);
	// 		assert(dotProduct < 1.1);
	// 		assert(dotProduct > -1.1);
	// 		if (dotProduct > 1.0) dotProduct = 1.0;
	// 		if (dotProduct < -1.0) dotProduct = -1.0;
	// 		bondBend = std::acos(dotProduct);

	// 		if(row[3] >= 0){
	// 			SimTK::Compound::AtomIndex
	// 				a3_cAIx = atoms[row[3]].identity.compoundAtomIndex;
	// 			SimTK::Vec3 a3loc = findAtomTarget(atomTargets, a3_cAIx);

	// 			bondTorsion = bDihedral(a0loc, a1loc, a2loc, a3loc);
	// 		}

	// 	} // angle
		
	// 	setZMatrixBATValue(rowCnt, 0, bondLength);
	// 	setZMatrixBATValue(rowCnt, 1, bondBend);
	// 	setZMatrixBATValue(rowCnt, 2, bondTorsion);

	// 	if(row[3] == -2){
	// 		topoIx++;
	// 	}

	// 	rowCnt++;

	// } // every zMatrix row		

}


/*!
 * <!--	zmatrixbat_ -->
*/
SimTK::Real Context::getZMatrixBATValue(size_t rowIndex, size_t colIndex) const {
	// Check if the indices are within bounds
	if (rowIndex < zMatrixBAT.size() && colIndex < zMatrixBAT[0].size()) {
		// Return the value at the specified position
		return zMatrixBAT[rowIndex][colIndex];
	} else {
		// Indices are out of bounds, handle this case accordingly
		return SimTK::NaN;
	}
}

/*!
 * <!--	zmatrixbat_ -->
*/
void Context::PrintZMatrixBAT() const {
	for (const auto& row : zMatrixBAT) {
		for (SimTK::Real value : row) {
			std::cout << std::setw(6) << value << " ";
		}
		std::cout << std::endl;
	}
}

/*!
 * <!--	zmatrixbat_ -->
*/
void Context::PrintZMatrixTableAndBAT() const 
{

	size_t zMatCnt = 0;

	for (const auto& row : zMatrixTable) {

		scout("ZMatrixBATEntry: ");

		// Print indexes
		for (int value : row) {
			std::cout << std::setw(9) << value <<" "; 
		}

		// Print BAT values
		const std::vector<SimTK::Real>& BATrow = getZMatrixBATRow(zMatCnt);
		
		for (SimTK::Real BATvalue : BATrow) {
			std::cout << std::setw(6) << BATvalue << " ";
		}

		ceol;

		zMatCnt++;
	}

	for(int rk = 0; rk < nofReplicas; rk++){
		scout("Replica ") << rk <<" " <<"BAT " << eol;
		replicas[rk].PrintZMatrixBAT();
	}

	for(size_t tk = 0;
	tk < nofThermodynamicStates;
	tk++){
		scout("ThermoState ") << tk <<" " <<"BAT " << eol;
		thermodynamicStates[tk].PrintZMatrixBAT();
	}

}

/*!
 * <!--	zmatrixbat_ -->
*/
void Context::addZMatrixBATRow(const std::vector<SimTK::Real>& newRow) {
	zMatrixBAT.push_back(newRow);
}

/*!
 * <!-- zmatrixbat_ BAT JAcobian -->
*/
SimTK::Real
Context::calcInternalBATJacobianLog(void)
	{

		// Get log of the Cartesian->BAT Jacobian
		SimTK::Real logJacBAT = 0.0;

		for(size_t zCnt = 0; zCnt < zMatrixBAT.size(); zCnt++){

				// Get bond term
				SimTK::Real currBond = zMatrixBAT[zCnt][0];
				
				if(currBond != SimTK::NaN){
				
					logJacBAT += 4.0 * std::log(currBond);
				}

				// Get the angle term
				SimTK::Real currAngle = zMatrixBAT[zCnt][1];

				if(currAngle != SimTK::NaN){

					logJacBAT += 2.0 * std::log(std::sin(currAngle));
					
				}

		}

		return logJacBAT;

	}


/*!
 * <!-- zmatrixbat_ Get BAT coordinates modifyable by a selected world -->
*/
void
Context::addSubZMatrixBATsToWorld(
	int wIx
	, int replicaIx)
{
	
	// // Get world
	// World& world = worlds[wIx];

	// // Get generalized coordinates
	// const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// // Iterate ZMatrix and bonds
	// int prevMolIx = -1;
	// size_t zMatCnt = 0;
	// for (const auto& row : zMatrixTable) {

	// 	std::vector<SimTK::Real>& BATrow = updZMatrixBATRow(zMatCnt);
	// 	std::vector<SimTK::Real>& replicaBATrow = replicas[replicaIx].updZMatrixBATRow(zMatCnt);

	// 	assert( (std::abs(BATrow[0] - replicaBATrow[0]) < 0.00001 || std::isnan(BATrow[0])) &&
	// 			(std::abs(BATrow[1] - replicaBATrow[1]) < 0.00001 || std::isnan(BATrow[1])) &&
	// 			(std::abs(BATrow[2] - replicaBATrow[2]) < 0.00001 || std::isnan(BATrow[2])) &&
	// 			"BATrow and replicaBATrow not equal.");

	// 	// scout("BATrow_ref") <<" ";
	// 	// for (SimTK::Real BATvalue : BATrow) {
	// 	// 	std::cout << std::setw(6) << BATvalue << " ";
	// 	// }
	// 	// for (SimTK::Real BATvalue : BATrow_ref) {
	// 	// 	std::cout << std::setw(6) << BATvalue << " ";
	// 	// }
	// 	// ceol;

	// 	// Get bond's atoms
	// 	Atom& childAtom  = atoms[row[0]];
	// 	Atom& parentAtom = atoms[row[1]];

	// 	// Get molecule
	// 	int childMolIx = childAtom.getMoleculeIndex();
	// 	int parentMolIx = parentAtom.getMoleculeIndex();
	// 	assert((childMolIx == parentMolIx) &&
	// 		"Atoms from different molecules");
	// 	Topology& topology = topologies[childMolIx];

	// 	// Iterate BONDS and get bond // ======================================
	// 	int BOIx;
	// 	if(prevMolIx != childMolIx){
	// 		BOIx = 0;
	// 	}

	// 	// Get bond
	// 	const std::vector<BOND>& BONDS = allBONDS[childMolIx];
	// 	const BOND& currBOND = BONDS[BOIx];
	// 	size_t boIx = BONDS_to_bonds[childMolIx][BOIx];
	// 	BondLink& bond = bonds[boIx];

	// 	// Insert reference to BAT entry into samplers if the bond is flexible
	// 	//scout(" ") << MobilityStr [ bond.getBondMobility(wIx) ] <<" ";
	// 	if(bond.getBondMobility(wIx) != SimTK::BondMobility::Rigid){

	// 		// Get Molmodel indexes
	// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.identity.compoundAtomIndex;
	// 		SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);

	// 		// Get mbx
	// 		SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
	// 		SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);

	// 		// Insert key and value into the map
	// 		for(size_t sami = 0; sami < worlds[wIx].samplers.size(); sami++){

	// 			//(pHMC((worlds[wIx].samplers[sami]))->subZMatrixBATs_ref).insert({childMbx, BATrow});
	// 			(pHMC((worlds[wIx].samplers[sami]))->subZMatrixBATs_ref).insert({childMbx, replicaBATrow});
	// 		}

	// 	} // is bond rigid

	// 	BOIx++;

	// 	if(prevMolIx != childMolIx){
	// 		prevMolIx = childMolIx;
	// 	} // every BOND -------------------------------------------------------

	// 	zMatCnt++;
	// }

	// //world.samplers[0].variableBATs = worldBATs;
		
}

/*!
 * <!-- zmatrixbat_ Get BAT coordinates modifyable by a selected world -->
*/
void
Context::updSubZMatrixBATsToWorld(
	  int wIx
	, int replicaIx)
{
	
	// // Get world
	// World& world = worlds[wIx];

	// // Get generalized coordinates
	// const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// // Iterate ZMatrix and bonds
	// int prevMolIx = -1;
	// size_t zMatCnt = 0;
	// for (const auto& row : zMatrixTable) {

	// 	std::vector<SimTK::Real>& replicaBATrow = replicas[replicaIx].updZMatrixBATRow(zMatCnt);

	// 	// Get bond's atoms
	// 	Atom& childAtom  = atoms[row[0]];
	// 	Atom& parentAtom = atoms[row[1]];

	// 	// Get molecule
	// 	int childMolIx = childAtom.getMoleculeIndex();
	// 	int parentMolIx = parentAtom.getMoleculeIndex();
	// 	assert((childMolIx == parentMolIx) &&
	// 		"Atoms from different molecules");
	// 	Topology& topology = topologies[childMolIx];

	// 	// Iterate BONDS and get bond // ======================================
	// 	int BOIx;
	// 	if(prevMolIx != childMolIx){
	// 		BOIx = 0;
	// 	}

	// 	// Get bond
	// 	const std::vector<BOND>& BONDS = allBONDS[childMolIx];
	// 	const BOND& currBOND = BONDS[BOIx];
	// 	size_t boIx = BONDS_to_bonds[childMolIx][BOIx];
	// 	BondLink& bond = bonds[boIx];

	// 	//scout(" ") << MobilityStr [ bond.getBondMobility(wIx) ] <<" ";
	// 	if(bond.getBondMobility(wIx) != SimTK::BondMobility::Rigid){

	// 		// Get Molmodel indexes
	// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.identity.compoundAtomIndex;
	// 		SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);

	// 		// Get mbx
	// 		SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
	// 		SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);

	// 		// Insert key and value into the map
	// 		for(size_t sami = 0; sami < worlds[wIx].samplers.size(); sami++){

	// 			(pHMC((worlds[wIx].samplers[sami]))->subZMatrixBATs_ref).at(childMbx) = replicaBATrow;

	// 		}

	// 	}

	// 	BOIx++;

	// 	if(prevMolIx != childMolIx){
	// 		prevMolIx = childMolIx;
	// 	} // every BOND -------------------------------------------------------

	// 	zMatCnt++;
	// }
		
}


/*!
 * <!-- Go through the vector of worlds and if their equilibrium worlds run and 
 * rotate. -->
*/
void Context::updSubZMatrixBATsToAllWorlds(int replicaIx)
{
	// // Get thermoState corresponding to this replica
	// int thisThermoStateIx = replica2ThermoIxs[replicaIx];

	// // Get this world indexes from the corresponding thermoState
	// std::vector<int>& replicaWorldIxs = 
	// 	thermodynamicStates[thisThermoStateIx].updWorldIndexes();

	// // Get nof worlds in this replica
	// size_t replicaNofWorlds = replicaWorldIxs.size();

	// // Set BAT map for all replica's worlds
	// int currFrontWIx = -1;
	// for(std::size_t worldCnt = 0; worldCnt < replicaNofWorlds; worldCnt++){

	// 	updSubZMatrixBATsToWorld( replicaWorldIxs[worldCnt], replicaIx);

	// } // every world

}


/*!
 * <!-- Print BAT coordinates from a selected world -->
 */
void
Context::PrintWorldSubZMatrixBATs(
	int wIx)
{
	
	// // Get world
	// World& world = worlds[wIx];

	// // Get generalized coordinates
	// const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// // Iterate ZMatrix and bonds
	// int prevMolIx = -1;
	// size_t zMatCnt = 0;
	// for (const auto& row : zMatrixTable) {

	// 	std::vector<SimTK::Real>& BATrow = updZMatrixBATRow(zMatCnt);

	// 	// Get bond's atoms
	// 	Atom& childAtom  = atoms[row[0]];
	// 	Atom& parentAtom = atoms[row[1]];

	// 	// Get molecule
	// 	int childMolIx = childAtom.getMoleculeIndex();
	// 	int parentMolIx = parentAtom.getMoleculeIndex();
	// 	assert((childMolIx == parentMolIx) &&
	// 		"Atoms from different molecules");
	// 	Topology& topology = topologies[childMolIx];

	// 	// Iterate BONDS and get bond // ======================================
	// 	int BOIx;
	// 	if(prevMolIx != childMolIx){
	// 		BOIx = 0;
	// 	}

	// 	// Get bond
	// 	const std::vector<BOND>& BONDS = allBONDS[childMolIx];
	// 	const BOND& currBOND = BONDS[BOIx];
	// 	size_t boIx = BONDS_to_bonds[childMolIx][BOIx];
	// 	BondLink& bond = bonds[boIx];

	// 	//scout(" ") << MobilityStr [ bond.getBondMobility(wIx) ] <<" ";
	// 	if(bond.getBondMobility(wIx) != SimTK::BondMobility::Rigid){

	// 		// Get Molmodel indexes
	// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.identity.compoundAtomIndex;
	// 		SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);

	// 		// Get mbx
	// 		SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
	// 		SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);

	// 		// Insert key and value into the map
	// 		for(size_t sami = 0; sami < worlds[wIx].samplers.size(); sami++){

	// 			std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>&>&
	// 				variableBATs = pHMC((worlds[wIx].samplers[sami]))->updSubZMatrixBATsRef();
					
	// 			scout("WorldBAT ") << wIx <<" "; 
	// 			for(auto varBAT : variableBATs.at(childMbx)){
	// 				cout << varBAT <<" ";
	// 			}
	// 			ceol;
				

	// 		}

	// 	}

	// 	BOIx++;

	// 	if(prevMolIx != childMolIx){
	// 		prevMolIx = childMolIx;
	// 	} // every BOND -------------------------------------------------------

	// 	zMatCnt++;
	// }
		
}


/*!
 * <!-- zmatrixbat_ Relationship BAT - mobod transforms -->
*/
void Context::PrintZMatrixMobods(int wIx, SimTK::State& someState)
{

	// // Get world
	// World& world = worlds[wIx];

	// // Get generalized coordinates
	// SimTK::Vector qVector = someState.getQ();
	// std::cout << "Q= " << qVector << eol;

	// const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// // Iterate ZMatrix and bonds
	// int prevMolIx = -1;
	// size_t zMatCnt = 0;
	// for (const auto& row : zMatrixTable) {

	// 	scout("ZMatrixBATEntry: ");

	// 	// Print indexes
	// 	for (int value : row) {
	// 		std::cout << std::setw(6) << value <<" "; 
	// 	}

	// 	// Print BAT values
	// 	const std::vector<SimTK::Real>& BATrow = getZMatrixBATRow(zMatCnt);
	// 	for (SimTK::Real BATvalue : BATrow) {
	// 		std::cout << std::setw(6) << BATvalue << " ";
	// 	}

	// 	// Get bond's atoms
	// 	Atom& childAtom  = atoms[row[0]];
	// 	Atom& parentAtom = atoms[row[1]];

	// 	// Get molecule
	// 	int childMolIx = childAtom.getMoleculeIndex();
	// 	int parentMolIx = parentAtom.getMoleculeIndex();
	// 	assert((childMolIx == parentMolIx) &&
	// 		"Atoms from different molecules");
	// 	Topology& topology = topologies[childMolIx];


	// 	// Get current bond
	// 	int BOIx;
	// 	if(prevMolIx != childMolIx){
	// 		BOIx = 0;
	// 	}

	// 	const std::vector<BOND>& BONDS = allBONDS[childMolIx];
	// 	const BOND& currBOND = BONDS[BOIx];
	// 	size_t boIx = BONDS_to_bonds[childMolIx][BOIx];
	// 	BondLink& bond = bonds[boIx];

	// 	//scout(" ") << MobilityStr [ bond.getBondMobility(wIx) ] <<" ";
	// 	if(bond.getBondMobility(wIx) != SimTK::BondMobility::Rigid){

	// 		// Get Molmodel indexes
	// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.identity.compoundAtomIndex;
	// 		SimTK::Compound::AtomIndex parent_cAIx = parentAtom.identity.compoundAtomIndex;

	// 		SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
	// 		SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);

	// 		// Get child-parent mobods
	// 		SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
			
	// 		SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);
	// 		const SimTK::MobilizedBody &childMobod = world.getMatterSubsystem().getMobilizedBody(childMbx);
	// 		SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);
	// 		const SimTK::MobilizedBody &parentMobod = world.getMatterSubsystem().getMobilizedBody(parentMbx);

	// 		childMobod.getFirstQIndex(someState);

	// 		scout(" ") << childMbx <<" " << parentMbx <<" ";

	// 		scout("| ")
	// 			<< childMobod.getQAsVector(someState) <<" |" ;

	// 		//scout("| ")
	// 		//	<< parentMobod.getQAsVector(someState) <<" |";
	// 		// if(mbx != parentMbx){
	// 		// 	// Get default transforms
	// 		// 	const SimTK::Transform& X_PF = mobod.getInboardFrame(someState);
	// 		// 	const SimTK::Transform& X_BM = mobod.getOutboardFrame(someState);
	// 		// 	const SimTK::Transform& X_FM = mobod.getMobilizerTransform(someState);
	// 		// 	PrintTransform(X_PF, 6, "X_PF");
	// 		// 	PrintTransform(X_BM, 6, "X_BM");
	// 		// 	PrintTransform(X_FM, 6, "X_FM");
	// 		// }else{
	// 		// 	//ceol;
	// 		// }
	// 	}

	// 	ceol;

	// 	BOIx++;

	// 	if(prevMolIx != childMolIx){
	// 		prevMolIx = childMolIx;
	// 	}

	// 	zMatCnt++;
	// }

}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ZMatrix BAT
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------



// ===========================================================================
// TRANSFORMERS LAB
// ===========================================================================

/*!
 * <!--  -->
*/
void Context::Print_TRANSFORMERS_Work(void)
{

		// scout("Transformers table by atom: no dAIx topoIx cAIx mobods") << eol;
		// size_t cnt = 0;
		// for(auto &atom : atoms){

		// 	size_t topoIx = atom.getMoleculeIndex();
		// 	Topology& topology = topologies[topoIx];

		// 	const SimTK::Compound::AtomIndex cAIx = atom.identity.compoundAtomIndex;

		// 	// Get dumm atom index (set in modelOneCompound Step 1)
		// 	SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(cAIx);

		// 	// Get vector of worlds mbxs which contain this atom
		// 	std::vector<SimTK::MobilizedBodyIndex> worldsMbxs;

		// 	size_t wIx = 0;
		// 	for(auto &world : worlds){
		// 		SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
		// 		SimTK::MobilizedBodyIndex mbx = dumm.getAtomBody(dAIx);
		// 		worldsMbxs.push_back(mbx);
		// 		wIx++;
		// 	}

		// 	// What do we have so far
		// 	scout("atomEntry: ") 
		// 		<< cnt <<" "
		// 		<< dAIx <<" "
		// 		<< topoIx <<" "
		// 		<< cAIx <<" ";

		// 	scout(" ");	
		// 	for(auto mbx: worldsMbxs){
		// 		std::cout << mbx <<" ";
		// 	}
		// 	ceol;	

		// 	// Increase atom number
		// 	cnt++;
		// } // every atom

		// // --------------------------------------------------------------------

		// scout("Transformers table by bond: allCnt topoIx BOIx boIx BOchild BOparent child_cAIx parent_cAIx child_dAIx parent_dAIx childMbx parentMbx flex flexStr") << eol;

		// const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

		// assert(allBONDS.size() == numMolecules 
		// 	&& "internal coordinates nof molecules wrong");

		// // Counter for all bonds
		// size_t allCnt = 0;

		// // Iterate molecules
		// for(size_t topoIx = 0; topoIx < numMolecules; topoIx++){


		// 	// Get molecule and it's bonds
		// 	Topology& topology = topologies[topoIx];
		// 	const std::vector<BOND>& BONDS = allBONDS[topoIx];

		// 	// Iterate molecule's bonds
		// 	for(size_t BOIx = 0; BOIx < BONDS.size(); BOIx++){

		// 		// Get current bond
		// 		const BOND& currBOND = BONDS[BOIx];
		// 		size_t boIx = BONDS_to_bonds[topoIx][BOIx];
		// 		BondLink& bond = bonds[boIx];

		// 		// Get bond's atoms
		// 		Atom& childAtom  = atoms[currBOND.first];
		// 		Atom& parentAtom = atoms[currBOND.second];

		// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.identity.compoundAtomIndex;
		// 		SimTK::Compound::AtomIndex parent_cAIx = parentAtom.identity.compoundAtomIndex;

		// 		SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
		// 		SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);



		// 		// Is it a base atom
		// 		// const SimTK::Compound::SingleAtom &parentCompoundAtom = parentAtom.getSingleAtom();
		// 		// SimTK::Compound::AtomIndex parCAIx = parentAtom.identity.compoundAtomIndex;
		// 		// const SimTK::Compound::AtomPathName parAtomPathName = parentCompoundAtom.getAtomName(parCAIx);
		// 		if(parentAtom.isRoot() == true){

		// 			scout("bondEntry: -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ");

		// 			// Iterate worlds and get child-parent mobods
		// 			size_t wIx = 0;
		// 			for(auto &world : worlds){
		// 				SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
		// 				SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);
		// 				const SimTK::MobilizedBody &parentMobod = world.getMatterSubsystem().getMobilizedBody(parentMbx);
		// 				const SimTK::MobilizedBody &grandMobod = parentMobod.getParentMobilizedBody();
		// 				SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();

		// 				std::cout << parentMbx <<" " << grandMbx <<" "
		// 					<< getMobility(rootMobilitiesStr[wIx][topoIx]) <<" "
		// 					<< rootMobilitiesStr[wIx][topoIx] <<" ";
						
		// 				wIx++;
		// 			} ceol;
			
		// 		}


		// 		// What we have so far
		// 		scout("bondEntry: ") << allCnt <<" " << topoIx <<" " << BOIx <<" " << boIx <<" ";
		// 		BONDS[BOIx].Print();

		// 		scout(" ") << child_cAIx <<" " << parent_cAIx <<" " << child_dAIx <<" " << parent_dAIx <<" ";

		// 		// Get vector of worlds mbxs which contain this atom
		// 		std::vector<std::pair<
		// 			SimTK::MobilizedBodyIndex, SimTK::MobilizedBodyIndex>> worldsMbxBonds;

		// 		// Iterate worlds and get child-parent mobods
		// 		size_t wIx = 0;
		// 		for(auto &world : worlds){
		// 			SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
		// 			SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);
		// 			SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);

		// 			worldsMbxBonds.push_back(std::pair<SimTK::MobilizedBodyIndex, SimTK::MobilizedBodyIndex>{childMbx, parentMbx});

		// 			wIx++;
		// 		}

		// 		// Iterate worlds and print
		// 		wIx = 0;	
		// 		for(auto mbxPair: worldsMbxBonds){
		// 			std::cout << mbxPair.first <<" " << mbxPair.second <<" " << bond.getBondMobility(wIx) <<" " << MobilityStr[ bond.getBondMobility(wIx) ] <<" ";
		// 			wIx++;
		// 		} ceol;

		// 		allCnt++;

		// 	} // every bond

		// } // every molecule
}

// void Context::PrintUDot(void)
// {
// 	size_t wIx = 1; // We want the U and UDot of the torsional dynamics world
// 	const auto& currentAdvancedState = worlds[wIx].updIntegrator().updAdvancedState();
// 	SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[wIx].updForceField());

// 	const auto& U = worlds[wIx].updSampler(0)->UCache;
// 	const auto& UDot = worlds[wIx].updSampler(0)->UDotCache;

// 	// std::cout << "U.size() = " << U.size() << std::endl;
// 	// std::cout << "UDot.size() = " << UDot.size() << std::endl;

// 	const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();
// 	assert(allBONDS.size() == numMolecules && "internal coordinates nof molecules wrong");

// 	// Iterate molecules
// 	for(size_t topoIx = 0; topoIx < numMolecules; topoIx++){

// 		// Get molecule and it's bonds
// 		Topology& topology = topologies[topoIx];
// 		const std::vector<BOND>& BONDS = allBONDS[topoIx];

// 		// Iterate molecule's bonds
// 		for(size_t BOIx = 0; BOIx < BONDS.size(); BOIx++){

// 			// Get current bond
// 			const BOND& currBOND = BONDS[BOIx];
// 			size_t boIx = BONDS_to_bonds[topoIx][BOIx];
// 			BondLink& bond = bonds[boIx];

// 			// Get bond's atoms
// 			Atom& childAtom  = atoms[currBOND.first];
// 			Atom& parentAtom = atoms[currBOND.second];

// 			SimTK::Compound::AtomIndex child_cAIx = childAtom.identity.compoundAtomIndex;
// 			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.identity.compoundAtomIndex;

// 			SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
// 			SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);

// 			SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);
// 			SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);

// 			if (childMbx != parentMbx) {
// 				int min_aix = std::min(currBOND.first, currBOND.second);
// 				int max_aix = std::max(currBOND.first, currBOND.second);
// 				std::string key = std::to_string(min_aix) + "-" + std::to_string(max_aix);
// 				// std::cout << "key " << key << " childMbx " << childMbx - 2 << std::endl;

// 				// MobilizedBodyIndex starts from 1, so we subtract 1 to match our indexing
// 				// I think the first one is the ground, so we subtract 2 to get the correct index
// 				const auto& u = U[childMbx - 2];
// 				uCache[key].push_back(u);

// 				const auto& uDot = UDot[childMbx - 2];
// 				uDotCache[key].push_back(uDot);
// 			}
// 		} // every bond
// 	} // every molecule

// }
// TRANSFORMERS LAB
// ===========================================================================

//////////////////////////////////////////////////////
//-------------         Q Stats        ---------------
//////////////////////////////////////////////////////

/*!
 * <!--  -->
*/
void Context::reserveThermostatsQs(void)
{
	// // Iterate thermodynamic states
	// for(size_t thIx = 0; thIx < nofThermodynamicStates; thIx++){
	// 	const std::vector<int> & thermoWorldIxs = thermodynamicStates[thIx].getWorldIndexes();
	// 	thermodynamicStates[thIx].allocateQStats(thermoWorldIxs.size());
	// }
}

/*!
 * <!--  -->
*/
void Context::setThermostatesQs(void)
{

	// // Iterate thermodynamic states
	// for(size_t thIx = 0; thIx < nofThermodynamicStates; thIx++){
	// 	const std::vector<int> & thermoWorldIxs = thermodynamicStates[thIx].getWorldIndexes();
	// 	// Iterate worlds
	// 	for(const auto worldIx : thermoWorldIxs){
	// 		SimTK::State& worldCurrentState = worlds[worldIx].updIntegrator().updAdvancedState();
	// 		//int NQ = (worlds[worldIx].getSimbodyMatterSubsystem())->getNQ(worldCurrentState);
	// 		thermodynamicStates[thIx].setWorldQs(worldIx, (getWorld(0).getSimbodyMatterSubsystem())->getQ(worldCurrentState));
	// 	}
	// }

}


/*!
 * <!--  -->
*/
void Context::printQStats(int thIx)
{
	thermodynamicStates[thIx].printQStats();
}
