#include "Context.hpp"
#include "World.hpp"
#include "Sampler.hpp"
#include "readAmberInput.hpp"

#include <sys/stat.h>
#include <sys/sysinfo.h>

#include <stack>

/*!
 * <!-- Constructor: sets temperatures, random engine and checks for CUDA_ROOT
 * -->
*/
Context::Context(const std::string& baseName_arg, uint32_t seed, uint32_t threads, uint32_t nofRoundsTillReblock, RUN_TYPE runType, uint32_t swapFreq, uint32_t swapFixmanFreq)
{
	// Set the base name of the simulation
	std::cout << "Context with base name: " << baseName + "_" + std::to_string(seed) << std::endl << std::flush;
	this->baseName = baseName_arg + "_" + std::to_string(seed);

	// Alert user of CUDA environment variables
#if OPENMM_PLATFORM_CUDA
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

	foutU = std::string(baseName + "_U.bin");
	foutUDot = std::string(baseName + "_U_dot.bin");
	foutTorque = std::string(baseName + "_torque.bin");
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

std::vector<TopologyRange> Context::findMoleculeRnages() const {
	std::vector<TopologyRange> ranges;

	// Create molecule index-based spans for atoms, bonds, angles and torsions
	for (std::size_t molIx = 0; molIx < roots.size(); molIx++) {

		TopologyRange r;

		std::size_t begin = std::numeric_limits<std::size_t>::max();
		std::size_t end   = std::numeric_limits<std::size_t>::min();
		bool found = false;

		// Search atom
		for (std::size_t i = 0; i < atoms.size(); ++i) {
			if (atoms[i].getMoleculeIndex() != molIx) continue;
			begin = std::min(begin, i);
			end   = std::max(end, i);
			found = true;
		}

		SimTK_ASSERT_ALWAYS(found && end >= begin, std::string("No atoms found for molecule index " + std::to_string(molIx)).c_str());
		SimTK_ASSERT_ALWAYS(end < atoms.size(), std::string("Atom index out of range for molecule " + std::to_string(molIx)).c_str());

		r.atomRange.begin = begin;
		r.atomRange.end = end + 1; // end is exclusive

		// Search bonds
		begin = std::numeric_limits<std::size_t>::max();
		end   = std::numeric_limits<std::size_t>::min();
		found = false;

		for (std::size_t i = 0; i < bonds.size(); ++i) {
			if (bonds[i].getMoleculeIndex() != molIx) continue;
			begin = std::min(begin, i);
			end   = std::max(end, i);
			found = true;
		}

		SimTK_ASSERT_ALWAYS(found && end >= begin, std::string("No bonds found for molecule index " + std::to_string(molIx)).c_str());
		SimTK_ASSERT_ALWAYS(end < bonds.size(), std::string("Bond index out of range for molecule " + std::to_string(molIx)).c_str());

		r.bondRange.begin = begin;
		r.bondRange.end = end + 1; // end is exclusive

		// // Search angles
		// begin = angles.size();
		// end   = 0;
		// found = false;

		// for (std::size_t i = 0; i < angles.size(); ++i) {
		// 	begin = std::min(begin, i);
		// 	end   = std::max(end, i);
		// 	found = true;
		// }

		// SimTK_ASSERT_ALWAYS(found && end >= begin, std::string("No angles found for molecule index " + std::to_string(molIx)).c_str());
		// SimTK_ASSERT_ALWAYS(end < angles.size(), std::string("Angle index out of range for molecule " + std::to_string(molIx)).c_str());

		// r.angleRange.begin = begin;
		// r.angleRange.end = end + 1; // end is exclusive

		// // Search torsions
		// begin = torsions.size();
		// end   = 0;
		// found = false;

		// for (std::size_t i = 0; i < torsions.size(); ++i) {
		// 	begin = std::min(begin, i);
		// 	end   = std::max(end, i);
		// 	found = true;
		// }

		// SimTK_ASSERT_ALWAYS(found && end >= begin, std::string("No torsions found for molecule index " + std::to_string(molIx)).c_str());
		// SimTK_ASSERT_ALWAYS(end < torsions.size(), std::string("Torsion index out of range for molecule " + std::to_string(molIx)).c_str());

		// r.torsionRange.begin = begin;
		// r.torsionRange.end = end + 1; // end is exclusive

		ranges.push_back(r);
	}

	return ranges;
}

void Context::loadAmberSystem(const std::vector<int>& inRoots, const std::vector<Atom>& inAtoms, const std::vector<BondStretch>& inBonds, const std::vector<BondBend>& inAngles, const std::vector<BondTorsion>& inTorsions) {
	roots = inRoots;
	atoms = inAtoms;
	bonds = inBonds;
	angles = inAngles;
	torsions = inTorsions;

	numMolecules = roots.size();

	// Construct a Compound for every atom
	for(auto& atom : atoms) {
		atom.createSingleAtom();
	}

	// Create spans based on molecule index
	const std::vector<TopologyRange> ranges = findMoleculeRnages();

	// Add new topologies
	topologies.reserve(numMolecules);
	for(std::size_t molIx = 0; molIx < roots.size(); molIx++) {

		// New empty topology
		SimTK::Compound::Name name = "MOL_" + std::to_string(molIx);
		Topology topology(name, SimTK::CompoundSystem::CompoundIndex(molIx), roots[molIx]);

		// Set spans
		topology.setAtoms(Span<Atom>(atoms.begin() + ranges[molIx].atomRange.begin, atoms.begin() + ranges[molIx].atomRange.end));
		topology.setBonds(Span<BondStretch>(bonds.begin() + ranges[molIx].bondRange.begin, bonds.begin() + ranges[molIx].bondRange.end));
		// topology.setAngles(Span<BondAngle>(angles.begin() + ranges[molIx].angleRange.begin, angles.begin() + ranges[molIx].angleRange.end));
		// topology.setTorsions(Span<BondTorsion>(torsions.begin() + ranges[molIx].torsionRange.begin, torsions.begin() + ranges[molIx].torsionRange.end));

		// Set root atom. Its compound atom index is set inside the next loop
		int rootAtomGlobalIndex = roots[molIx];
		Atom& rootAtom = atoms[rootAtomGlobalIndex];
		topology.setBaseAtom(rootAtom.getSingleAtom(), Transform());
		topology.convertInboardBondCenterToOutboard();

		// Add non-ring closing bonds first
		for(auto& bond : topology.updBonds()) {

			if (bond.isRingClosing()) continue;

			Atom& parent = atoms[bond.getParentAtomGlobalIndex()];
			Atom& child = atoms[bond.getChildAtomGlobalIndex()];

			// Get next available BondCenter id
			int parentNofBonds = parent.getNumBondsInvolved();
			int parentNofFreebonds = parent.getNumAvailableBonds();
			int parentNextAvailBondCenter = parentNofBonds - parentNofFreebonds + 1;

			// Cook the parentBondCenterPathName = RESNAME + RESID + _ATOMNAME + bond int(next)
			SimTK::Compound::BondCenterPathName parentBondCenterPathName = parent.getUniqueAtomName() + "/bond" + std::to_string(parentNextAvailBondCenter);

			// Actual bonding with default mobility (torsion)
			topology.bondAtom(child.getSingleAtom(), parentBondCenterPathName, 0.149, 0); // SimTK::BondMobility::Mobility = SimTK::BondMobility::Default

			// Set the local compound atom index for this child
			// 1 is child, 0 is for parent
			SimTK::Compound::AtomIndex childCAIx = topology.getBondAtomIndex(Compound::BondIndex(topology.getNumBonds() - 1), 1);
			child.setCompoundAtomIndex(childCAIx);
			topology.setAtomMass(childCAIx, child.getMassInDaltons());

			std::cout << "STEP_1: Set compound atom index for atom " << child.getGlobalIndex() << " to " << childCAIx << std::endl;

			// Set the local compound atom index for the parent if it is the root
			if(bond.getParentAtomGlobalIndex() == rootAtomGlobalIndex) {
				SimTK::Compound::AtomIndex parentCAIx = topology.getBondAtomIndex(Compound::BondIndex(topology.getNumBonds() - 1), 0);
				parent.setCompoundAtomIndex(parentCAIx);
				topology.setAtomMass(parentCAIx, parent.getMassInDaltons());
			}

			// // Handle ions
			// if (internCoords.getRoot(molIx).second == -1) {
			// 	atoms[internCoords.getRoot(molIx).first].setMoleculeIndex(molIx);
			// 	atoms[internCoords.getRoot(molIx).first].setCompoundAtomIndex(SimTK::Compound::AtomIndex(0));
			// }						

			parent.decrementAvailableBonds();
			child.decrementAvailableBonds();
		}

		// Add ring closing bonds
		for(const auto& bond : topology.updBonds()) {

			if (!bond.isRingClosing()) continue;

			Atom& parent = atoms[bond.getParentAtomGlobalIndex()];
			Atom& child = atoms[bond.getChildAtomGlobalIndex()];

			SimTK::Compound::BondCenterPathName bondCenterName1;
			if (child.getGlobalIndex() == rootAtomGlobalIndex) {
				bondCenterName1 = child.getUniqueAtomName() + "/bond" + std::to_string(child.getNumAvailableBonds());
			} else {
				bondCenterName1 = child.getUniqueAtomName() + "/bond" + std::to_string(child.getNumBondsInvolved() - child.getNumAvailableBonds() + 1);
			}

			SimTK::Compound::BondCenterPathName bondCenterName2;
			if (parent.getGlobalIndex() == rootAtomGlobalIndex) {
				bondCenterName2 = parent.getUniqueAtomName() + "/bond" + std::to_string(parent.getNumAvailableBonds());
			} else {
				bondCenterName2 = parent.getUniqueAtomName() + "/bond" + std::to_string(parent.getNumBondsInvolved() - parent.getNumAvailableBonds() + 1);
			}

			topology.addRingClosingBond(bondCenterName1, bondCenterName2, 0.14, 109*Deg2Rad, BondMobility::Rigid);					

			parent.decrementAvailableBonds();
			child.decrementAvailableBonds();
		}

		// Define the biotype of the atom
		for (auto& atom: topology.getAtoms()) {
			// It calls SimTK::Biotype::defineBiotype and checks if it already exists
			// topology.setAtomBiotype(atom.getUniqueAtomName().c_str(), atom.getResidueName().c_str(), atom.getAtomName().c_str());
			topology.setAtomBiotype(atom.getUniqueAtomName().c_str(), "", atom.getChargedAtomName().c_str());
			atom.setBiotypeIndex(topology.getAtomBiotypeIndex(atom.getCompoundAtomIndex()));
		}

		// Get coordinates
		Compound::AtomTargetLocations atomTargets;
		for (const auto& a : topology.getAtoms()) {
			SimTK::Vec3 atomCoords(a.getX(), a.getY(), a.getZ());
			atomTargets.insert(std::make_pair(a.getCompoundAtomIndex(), atomCoords));
		}

		// Top level shift corresponds to the root atom's position
		SimTK::Vec3 topLevelShift = SimTK::Vec3(rootAtom.getX(), rootAtom.getY(), rootAtom.getZ());
		topology.setTopLevelTransform(Transform(Rotation(), topLevelShift));

		bool flipAllChirality = true;
		topology.matchDefaultBondLengths(atomTargets);
		topology.matchDefaultAtomChirality(atomTargets, 0.01, flipAllChirality);
		topology.matchDefaultBondAngles(atomTargets);
		topology.matchDefaultDirections(atomTargets);
		topology.matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
		topology.matchDefaultTopLevelTransform(atomTargets);

		topology.loadIndicesMaps();

		topologies.push_back(topology);
	}

	// Add a check that all bonds have been satisfied

	// Context topologies to all the worlds
	// Iterate worlds
	for(size_t wCnt = 0; wCnt < worlds.size(); wCnt++){

		scout("World ") << wCnt << eol;

		// Get world and its force field
		World& world = worlds[wCnt];

		// Pass current topology to the current world
		world.topologies = &topologies;
	}
}

/*! <!--  --> */
void Context::Initialize() {

	worlds[0].setDuMMAtomIndexes(); // REVISE

	// @TODO take a look
	for (int worldIx = 0; worldIx < worlds.size(); worldIx++) {
		worlds[worldIx].setMyContext(this);
	}

	// Initialize the Z matrix
	int firstWIx = 0;
	SimTK::State& lastAdvancedState = worlds[firstWIx].integ->updAdvancedState();

	// Get coordinates from source
	const auto& firstWorldsAtomsLocations = worlds[firstWIx].getAtomsLocationsInGround(lastAdvancedState);

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

	// Initialize OpenMM
	OMMRef_initialize();
	//OMMRef_calcPotential(true, true);
}


#ifndef tracerefOMM
#define tracerefOMM(msg) std::cout<<__FILE__<<":"<<__LINE__<<" __refOMM__ "<<msg<<std::endl<<std::flush;
#endif

const double TOL = 1e-6;

std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3> computePeriodicBoxVectors_Context(
	double a_length, double b_length, double c_length,
    double alpha, double beta, double gamma)
{
    
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



/*! __refOMM__ <!-- Initialize OpenMM --> */
std::string Context::OMMRef_initialize(void)
{
	// Allocate OpenMM forces
	ommNonbondedForce = std::make_unique<OpenMM::NonbondedForce>();
	//auto ommGBSAOBCForce = std::make_unique<OpenMM::GBSAOBCForce>(); // TODO
	ommHarmonicBondStretch = std::make_unique<OpenMM::HarmonicBondForce>();
	ommHarmonicAngleForce = std::make_unique<OpenMM::HarmonicAngleForce>();
	ommPeriodicTorsionForce = std::make_unique<OpenMM::PeriodicTorsionForce>();

	// Instantiate the thermostat with adjusted temperature
	//Real temperature = 300.0;
	//if(dumm->wantOpenMMIntegration){temperature = dumm->temperature;}
	openMMThermostat = std::make_unique<OpenMM::AndersenThermostat>(300.0, 1);
	openMMThermostat->setRandomNumberSeed(seed);
	
	// Allocate OpenMM system and add particles to it
	openMMSystem = std::make_unique<OpenMM::System>();

    // ----------------------------------------------
    // PBC - Periodic Boundary Conditions __begin__
    // ----------------------------------------------
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
	ommNonbondedForce->setNonbondedMethod( OpenMM::NonbondedForce::NonbondedMethod( nonbondedMethod ) );
	ommNonbondedForce->setCutoffDistance( nonbondedCutoff );
	// nonbondedForce->setUseSwitchingFunction( 0 );

	// Scale charges by sqrt of scale factor so that products of charges scale linearly.
	const Real sqrtCoulombScale = std::sqrt(1.0);

	// Add atoms
	for (auto atom : atoms) {   
		SimTK::Real charge = atom.getChargeInE() * sqrtCoulombScale;
		const SimTK::Real sigma = atom.getSigmaInNm();
		const SimTK::Real epsilon = atom.getVdwWellDepthInKJ() * DuMM::KJ2Kcal * vdwGlobalScaleFactor;

		openMMSystem->addParticle(atom.getMassInDaltons());
		ommNonbondedForce->addParticle(charge, sigma, epsilon);
	}

	// Add bonds
	std::vector<std::pair<int, int>> ommBonds;
	for (auto bond : bonds) {
		ommBonds.emplace_back(std::make_pair(bond.getParentAtomGlobalIndex(), bond.getChildAtomGlobalIndex()));
	}

	// Register all the 1-2 bonds between nonbond atoms for scaling.
	// World::setAmberForceFieldScaleFactors(0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.8333333333, 1.0);
	ommNonbondedForce->createExceptionsFromBonds(ommBonds, 0.8333333333, 0.5);
	
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
				
	for (auto bond : bonds) {
		SimTK::Real nominalLengthInNm = bond.getNominalLengthInNm();
		SimTK::Real stiffnessInKJPerNmSq = bond.getStiffnessInKJPerNmSq();

		ommHarmonicBondStretch->addBond(bond.getParentAtomGlobalIndex(), bond.getChildAtomGlobalIndex(), nominalLengthInNm, stiffnessInKJPerNmSq);
	}

	// FORCES: ADD ANGLES (1-2-3)
	for (const auto& dummAngle : angles) {
		const int a1num = dummAngle.getGlobalIndex1();
		const int a2num = dummAngle.getGlobalIndex2();
		const int a3num = dummAngle.getGlobalIndex3();
		SimTK::Real theta0 = dummAngle.getNominalAngleInDeg() * DuMM::Deg2Rad;
		SimTK::Real forceKt = dummAngle.getStiffnessInKJPerRadSq() * DuMM::Kcal2KJ * 2.0;

		ommHarmonicAngleForce->addAngle(a1num, a2num, a3num, theta0, forceKt);
	}

	// Add dihedrals. OpenMM does not distinguish between proper and improper dihedrals.
	for (const auto& t : torsions) {
		int a1 = t.getGlobalIndex1();
		int a2 = t.getGlobalIndex2();
		int a3 = t.getGlobalIndex3();
		int a4 = t.getGlobalIndex4();
		int per = t.getPeriodicity();
		SimTK::Real phaseInRad = t.getPhaseInDegrees() * DuMM::Deg2Rad;
		SimTK::Real k_kJ = t.getAmpInKJ();

		ommPeriodicTorsionForce->addTorsion(a1, a2, a3, a4, per, phaseInRad, k_kJ);
	}

	openMMSystem->addForce(ommHarmonicBondStretch.get()); ommHarmonicBondStretch.release();
	openMMSystem->addForce(ommHarmonicAngleForce.get()); ommHarmonicAngleForce.release();
	openMMSystem->addForce(ommPeriodicTorsionForce.get()); ommPeriodicTorsionForce.release();
	openMMSystem->addForce(ommNonbondedForce.get()); ommNonbondedForce.release();

	// Get the thermostat
	openMMSystem->addForce(openMMThermostat.get()); openMMThermostat.release();
		
	// Get the integrator
	openMMIntegrator = std::make_unique<OpenMM::VerletIntegrator>(0.0007); // TODO should release?
		
    // Get the platform
    // By default, OpenMM builds a .so for each platform (CPU, OpenCL and CUDA)
    // When loading that .so, two functions get called
    // 1. registerPlatform() which does what you see below
    // 2. registerKernelFactories() which is used for Drude, Pme, Rpmd and other plugins (which we do not need as of right now)
#if OPENMM_PLATFORM_CPU
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        platform = std::make_unique<OpenMM::CpuPlatform>();
        OpenMM::Platform::registerPlatform(platform.get());
        platform.release();
    }
    constexpr auto PLATFORM_NAME = "CPU";

#elif OPENMM_PLATFORM_CUDA
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        platform = std::make_unique<OpenMM::CudaPlatform>();
        OpenMM::Platform::registerPlatform(platform.get());
        platform.release();
    }
    constexpr auto PLATFORM_NAME = "CUDA";

#elif OPENMM_PLATFORM_OPENCL
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        platform = std::make_unique<OpenMM::OpenCLPlatform>();
        OpenMM::Platform::registerPlatform(platform.get());
        platform.release();
    }
    constexpr auto PLATFORM_NAME = "OpenCL";
#endif

	bool allowReferencePlatform = true;
    // CREATE OPENMM CONTEXT based on PLATFORM
    try {
        auto& platform = OpenMM::Platform::getPlatformByName(PLATFORM_NAME);
        openMMContext = std::make_unique<OpenMM::Context>(*openMMSystem, *openMMIntegrator, platform);
        const double speed = openMMContext->getPlatform().getSpeed();

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

    std::cout << "Robosample Context reference OpenMM loaded " << openMMContext->getPlatform().getName() << std::endl;


////////////////////////////////////////////////////////////////////////////////
if("checkPotential"){
	std::vector<OpenMM::Vec3> atomsPositions = std::vector<OpenMM::Vec3>(atoms.size());
	int aCnt = -1;
	for (auto atom : atoms){
		aCnt++;
		atomsPositions[aCnt] = OpenMM::Vec3(atom.getX(), atom.getY(), atom.getZ());
	}
	openMMContext->setPositions(atomsPositions);
	OpenMM::State openMMState = openMMContext->getState(
		(true?OpenMM::State::Forces:0) | (true?OpenMM::State::Energy:0)
	);
	std::cout << "\nROBO_OpenMM POTENTIAL " << openMMState.getPotentialEnergy() << std::endl << std::flush;
}
/////////////////////////////////////////////////////////////////////////////////////////////



    return openMMContext->getPlatform().getName();
	
}


/*! __refOMM__
 * <!-- Calculate OpenMM energy -->
*/
SimTK::Real Context::OMMRef_calcPotential(const SimTK::Compound::AtomTargetLocations& atomTargets, bool wantEnergy, bool wantForces)
{
	if (!openMMContext) {
        std::cerr << "ERROR: OpenMM Context is not initialized!" << std::endl;
        return 0.0;
    }

    if (atoms.size() != openMMContext->getSystem().getNumParticles()) {
        std::cerr << "ERROR: Mismatch between atoms vector size (" << atoms.size()
                  << ") and OpenMM system particles (" 
                  << openMMContext->getSystem().getNumParticles() << ")!" << std::endl;
        return 0.0;
    }

	SimTK::Real refPotential = 0.0;
	std::vector<OpenMM::Vec3> ommAtomsPositions = std::vector<OpenMM::Vec3>(atoms.size());
	
	// Convert SimTK::Vec3 to OpenMM::Vec3
	for (const auto& atom : atoms) {
		SimTK::Compound::AtomIndex atomIndex = atom.getCompoundAtomIndex();
		const SimTK::Vec3& coords = atomTargets.at(atomIndex);
		ommAtomsPositions[atom.getGlobalIndex()] = OpenMM::Vec3(coords[0], coords[1], coords[2]);
	}

    // Pass the converted positions to OpenMM
    openMMContext->setPositions(ommAtomsPositions);

    // Ask for energy, forces, or both.
    OpenMM::State openMMState = openMMContext->getState(
        (wantForces?OpenMM::State::Forces:0) | (wantEnergy?OpenMM::State::Energy:0)
    );

    // std::cout << "Energy: " << openMMState.getPotentialEnergy() << std::endl;

    if (wantForces) {
        const std::vector<OpenMM::Vec3>& openMMForces = openMMState.getForces();
		int aCnt = -1;
		for (auto atom : atoms){
			aCnt++;
            const OpenMM::Vec3& ommForce = openMMForces[aCnt];
        }

    }

    if (wantEnergy){
		refPotential += openMMState.getPotentialEnergy();
	}


	//openMMState.getEnergies_drl_bon();
	if (verbose) {
		//std::cout << "Robosample reference OpenMM energy " << refPotential << std::endl;
	}

	return refPotential;
}

/*!
 * <!-- Add flexibilities and CompoundSystem model -->
*/
void Context::addWorld(bool fixmanTorque, int samplesPerRound, ROOT_MOBILITY rootMobility, const std::vector<BOND_FLEXIBILITY>& flexibilities, bool useOpenMM, bool visual, SimTK::Real visualizerFrequency) {

	// Create new world and add its index
	worldIndexes.push_back(worldIndexes.size());
	worlds.emplace_back(worldIndexes.back(), numMolecules, visual, visualizerFrequency);

	// @TODO should only pass natoms to worlds
	worlds.back().setMyContext(this);

	// Set force field scale factor.
	if (useAmberForceFieldScaleFactors) {
		worlds.back().setAmberForceFieldScaleFactors();
	} else {
		worlds.back().setGlobalForceFieldScaleFactor(globalForceFieldScaleFactor);
	}

	// Set the nonbonded method and cutoff
	worlds.back().updForceField()->setNonbondedMethod(nonbondedMethod);
	worlds.back().updForceField()->setNonbondedCutoff(nonbondedCutoff);

	// If requested, add Fixman torque as an additional force subsystem
	if (fixmanTorque) {
		worlds.back().addFixmanTorque();
		worlds.back().updFixmanTorque()->setScaleFactor(1);
	}

	// Set the number of threads for DuMM
	if (numThreads == 1) {
		worlds.back().updForceField()->setUseMultithreadedComputation(false);
	} else {
		worlds.back().updForceField()->setNumThreadsRequested(numThreads);
	}

	// Set GBSA scaling and VdW mixing rule
	worlds.back().setGbsaGlobalScaleFactor(gbsaGlobalScaleFactor);
	worlds.back().updForceField()->setVdwMixingRule(DuMMForceFieldSubsystem::LorentzBerthelot); // DuMMForceFieldSubsystem::WaldmanHagler

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
		worlds.back().forceField->setUseOpenMMAcceleration(true);
	}

	// Generate DuMM parameters: DuMM atom types, charged atom types, bond types, angle types and torsion types
	worlds.back().generateDummParams(atoms, bonds, angles, torsions);

	// Allocate root mobilities
	rootMobilitiesStr.push_back({});
	for(unsigned int molIx = 0; molIx < topologies.size(); molIx++){
		rootMobilitiesStr.back().push_back("Rigid");
	}

	// Set root mobilities
	for (const auto& flex : flexibilities) {

		if(flex.i == -1){

			assert("Set root mobilities: atom index fault." &&
				((flex.j >= 0) && (flex.j < atoms.size())));

			Atom& userRootAtom = atoms[flex.j];

			int molIx = userRootAtom.getMoleculeIndex();
			if(molIx >= numMolecules){
				std::cerr << "Set root mobilities: molecule index fault." << std::endl;
				break;
			}

			std::cout << "Set root mobilities -1=" << flex.i << " molecule " << molIx <<" at atom " << flex.j <<" to " << flex.mobility << std::endl;

		} // found a root mobility
	} // every flexibility

	// Set flexibilities
	for (auto& topology: topologies) {
		for (auto& bond : topology.getBonds()) {
			// Default is Rigid
			BondMobility::Mobility mobility = BondMobility::Mobility::Rigid;

			for (const auto& flex : flexibilities) {
				if ((bond.getParentAtomGlobalIndex() == flex.i && bond.getChildAtomGlobalIndex() == flex.j) || (bond.getParentAtomGlobalIndex() == flex.j && bond.getChildAtomGlobalIndex() == flex.i)) {
					mobility = flex.mobility;
					break;
				}
			}

			bond.addBondMobility(mobility);
		} // every bond


	}
	// for (auto& bond : orderedBonds) {
	// 	Topology& topology = topologies[bond.getMoleculeIndex()];
	// 	SimTK::Compound::BondIndex compoundBondIx = bond.getBondIndex();

	// 	// Check if the user set the bond mobility
	// 	BondMobility::Mobility mobility = BondMobility::Mobility::Rigid;
	// 	for (const auto& flex : flexibilities) {
	// 		if (bond.isThisMe(flex.i, flex.j)) {
	// 			mobility = flex.mobility;
	// 			break;
	// 		}
	// 	}

	// 	bond.addBondMobility(mobility);
	// 	topology.setBondMobility(mobility, compoundBondIx);
	// }

	// 
	worlds.back().AllocateCoordBuffers(atoms.size());

	// desk_mass_related
	// for (auto& atom : atoms) {
	// 	SimTK::DuMM::AtomIndex dAIx = atom.updDuMMAtomIndex();
	// 	SimTK::mdunits::Mass atomMass = atom.getMass();
	// 	worlds.back().forceField->setDuMMAtomMass(dAIx, atomMass);
	// } // every atom

	worlds.back().topologies = &topologies;
	for(std::size_t topologyIx = 0; topologyIx < topologies.size(); topologyIx++) {

		// Add topologies to CompoundSystem and add it to the visualizer's vector of molecules
		worlds.back().adoptTopology(topologyIx);

		// Was "Cartesian"
		worlds.back().compoundSystem->modelOneCompound(SimTK::CompoundSystem::CompoundIndex(topologyIx), topologies[topologyIx].updAtomFrameCache(), "Rigid");

		// Initialize the cache
		for (const auto& atom : topologies[topologyIx].getAtoms()) {
			SimTK::Compound::AtomIndex aIx = atom.getCompoundAtomIndex();
			worlds.back().atomTargetLocaltionsCache.insert(std::make_pair(atom.getCompoundAtomIndex(), atom.getCoords()));
		}
	}

	// This is one to many map
	worlds.back().loadMbx2AIxMap();

	// Allocate whatever needed Simbody dependent vectors from World here
	worlds.back().allocateStatsContainers();

	// // print the inertia tensors of all bodies
	// // PrintSimbodyMobods();
	// for(std::size_t worldIx = 0; worldIx < nofWorlds; worldIx++){
	// 	//TRACE("Context::PrintSimbodyMobods world " << worldIx);

	// 	for(std::size_t molIx = 0; molIx < numMolecules; molIx++){

	// 		//TRACE("Context::PrintSimbodyMobods molecule " << molIx);
	// 		const Topology& topology = worlds[worldIx].getTopology(molIx);

	// 		for(std::size_t i = 0; i < topology.getNumAtoms(); i++){
	// 			SimTK::Compound::AtomIndex aIx = (topology.subAtomList[i]).getCompoundAtomIndex();
	// 			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(aIx);

	// 			//TRACE("i = " << i << "; aIx = " << aIx << "; mbx = " << mbx << ";");
	// 		}
	// 	}
	// }
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
		for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
			worlds[worldIx].addRodConstraint(
				worlds[worldIx].integ->updAdvancedState()
			);
		}
}

// Print status
void Context::printStatus(void){
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		if (worlds[worldIx].integ  == nullptr ){
			std::cout << "Context: integrator is null" << std::endl;
			break;
		}
		SimTK::VerletIntegrator& checkIntegrator = *worlds[worldIx].integ;
		const SimTK::State& checkState = checkIntegrator.getState();
		const SimTK::Stage& checkStage = checkState.getSystemStage();
		std::cout << "Context world " << worldIx << " integ state stage "
			<< checkStage << std::endl << std::flush;
		std::cout << "Context world " << worldIx << " integ state nof Subsystems "
			<< checkState.getNumSubsystems() << ":" << std::endl << std::flush;
		for(int i = 0; i < checkState.getNumSubsystems(); i++){
			std::cout
				<< " Subsystem Name: "
				<< checkState.getSubsystemName(SimTK::SubsystemIndex(i))
				<< " Stage: "
				<< checkState.getSubsystemStage(SimTK::SubsystemIndex(i))
				<< " Version: "
				<< checkState.getSubsystemVersion(SimTK::SubsystemIndex(i))
				<< std::endl << std::flush;
		}
		//SimTK::State& checkAdvState = checkIntegrator.updAdvancedState();
		//const SimTK::Stage& checkAdvStage = checkAdvState.getSystemStage();
		//std::cout << "Context world " << worldIx << " integ advState stage "
		//	<< checkAdvStage << std::endl << std::flush;


		// CompoundSystem <- MolecularMechanicsSystem <- MultibodySystem <- System
		SimTK::CompoundSystem& compoundSystem = *(worlds[worldIx].getCompoundSystem());
		std::cout << "Context world " << worldIx << " compoundSystem nof compounds "
			<< compoundSystem.getNumCompounds() << std::endl;
		std::cout << "Context world " << worldIx << " System Topology realized "
			<< compoundSystem.getNumRealizationsOfThisStage(SimTK::Stage::Topology)
			<< " times.\n" << std::flush;

		// Matter
		////const SimTK::System& checkSystem = (worlds[worldIx].matter)->getSystem();
		SimTK::SimbodyMatterSubsystem& matter = *(worlds[worldIx].matter);
		std::cout << "Context world " << worldIx
			<< " matter nofBodies " << matter.getNumBodies()
			<< " nofConstraints " << matter.getNumConstraints()
			<< "\n" << std::flush;

		// GeneralForceSubsystem
		SimTK::GeneralForceSubsystem& gfs = *(worlds[worldIx].forces);
		std::cout << "Context world " << worldIx
			<< " gfs nofForces " << gfs.getNumForces()
			<< "\n" << std::flush;

		SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[worldIx].forceField);
		std::cout << "Context world " << worldIx
			<< " dumm nofThreads " << dumm.getNumThreadsRequested()
			<< " useOpenMM " << dumm.getUseOpenMMAcceleration()
			<< " " << dumm.isUsingOpenMM()
			<< "\n" << std::flush;


	}
}

// Print thermodynamics
void Context::printThermodynamics()
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		std::cout << "World " << worldIx << " temperature = "
			<< worlds[worldIx].getTemperature()
			<< std::endl;
		if(worlds[worldIx].isUsingFixmanTorque()){
			std::cout << "World " << worldIx
			<< " FixmanTorque temperature = "
			<< worlds[worldIx].updFixmanTorque()->getTemperature()
			<< std::endl;
		}
		for (int samplerIx = 0; samplerIx < worlds[worldIx].getNofSamplers(); samplerIx++){
			std::cout << "World " << worldIx << " Sampler " << samplerIx
				<< " temperature = " << worlds[worldIx].updSampler(samplerIx)->getTemperature()
				<< " initial const state PE: " << std::setprecision(20)
				//<< worlds[worldIx].forces->getMultibodySystem().calcPotentialEnergy(worlds[worldIx].integ->updAdvancedState())
				//<< worlds[worldIx].forces->getMultibodySystem().calcPotentialEnergy(updAdvancedState(worldIx, samplerIx))
				<< " useFixmanPotential = "
				<< pHMC(worlds[worldIx].updSampler(samplerIx))->isUsingFixmanPotential()
				<< std::endl;
		}

	}
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
	// 				= (topology.subAtomList[i]).getCompoundAtomIndex();
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

World& Context::getWorld(std::size_t whichWorld)
{
	return worlds[whichWorld];
}

const World& Context::getWorld(std::size_t whichWorld) const
{
	return worlds[whichWorld];
}

std::size_t Context::getNofWorlds() const
{
	return nofWorlds;
}

std::vector<World>& Context::getWorlds() {
	return worlds;
}

const std::vector<World>& Context::getWorlds() const {
	return worlds;
}

/////////////////////////
// --- Mixing parameters ---
/////////////////////////

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
	// Adaptive Gibbs blocking: // TODO generalized coord may not always be Real
	if(QsCache[0][0].size() == 0){
		std::size_t worldIx = 0;
		for(auto& world : worlds) {
			int nQs = world.getCompoundSystem()->getMatterSubsystem().getSystem().getDefaultState().getNQ();
			//std::cout << "World " << worldIx  << " has " << nQs << " Qs" << std::endl;
			//std::cout << "Context::realizeTopology QsCache[" << worldIx << "] size " << QsCache[worldIx].size() << std::endl;

			for(int t = 0; t < roundsTillReblock; t++) { // TODO use insert (why use insert?)
				for(int qi = 0; qi < nQs; qi++){
					QsCache[worldIx][t].push_back(0);
				}

			//std::cout << "Context::realizeTopology QsCache[" << worldIx << "]["<< t << "] size " << QsCache[worldIx][t].size() << std::endl;
			}

			worldIx++;
		}
	}
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
		topology.setMultibodySystem(*worlds[newWorldIx].compoundSystem);

		// Reset mobilized body indices in Compound
		for (const auto& atom : topology.getAtoms()) {
			SimTK::Compound::AtomIndex aIx = atom.getCompoundAtomIndex();
			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, *(worlds[newWorldIx].forceField));
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
void Context::addReplica(int index)
{
	// Add replica and the vector of worlds
	replicas.emplace_back(Replica(index
		, atoms
		, roots
		, topologies
		, zMatrixTable));

	// Set replicas coordinates

    // std::vector<std::vector<std::pair <Atom *, SimTK::Vec3>>> referenceAtomsLocationsFromFile;
	// replicas.back().setAtomsLocationsInGround(referenceAtomsLocationsFromFile);
	// replicas.back().set_WORK_AtomsLocationsInGround(referenceAtomsLocationsFromFile);

	const auto& referenceAtomsLocations = worlds[0].getCurrentAtomsLocationsInGround();
	replicas.back().setAtomsLocationsInGround(referenceAtomsLocations);
	replicas.back().set_WORK_AtomsLocationsInGround(referenceAtomsLocations);

	// Increment nof replicas
	nofReplicas++;
}

void Context::addThermodynamicState(
	int index,
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
		ThermodynamicState(index,
			T,
			argWorldIndexes,
			timestepsInThisReplica,
			mdstepsInThisReplica,
			atoms,			
			zMatrixTable
			//, zMatrixBAT
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
	thermodynamicStates.back().appendLog(baseName + ".repl" + std::to_string(index) + ".csv");
	thermodynamicStates.back().appendDCDReporter(baseName + ".repl" + std::to_string(index) + ".dcd", atoms.size(), topologies.size());

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

/*!
 * <!-- Calculate Fixman potential of replica J in replica I's back world. Ui(X_j) -->
*/
SimTK::Real Context::calcFixman(int replica_i, int replica_j)
{
	SimTK_ASSERT_ALWAYS(true, "Context::calcFixman: Not implemented yet.");

	// SimTK::Real Ui = replicas[replica_i].getFixman();

	// if (replica_i == replica_j){ // same replica
	// 	return Ui;
	// }else if(Ui <= 0.0000001){ // fully flexible world
	// 	return Ui;
	// }else{
	//     //std::cout << "CALC FIXMAN" << replica_i << " replica " << replica_j << "\n" << std::flush;
	// 	// Get replica i thermodynamic state
    //     int thermoState_i = replica2ThermoIxs[replica_i];

	// 	// Get replica i back world
    //     int world_i_front = thermodynamicStates[thermoState_i].getWorldIndexes().front();
    //     int world_i_back = thermodynamicStates[thermoState_i].getWorldIndexes().back();

	// 	// Get coordinates from replica j
	// 	const std::vector<std::vector<
    //         std::pair <Atom *, SimTK::Vec3>>>&
	// 	X_j = replicas[replica_j].getAtomsLocationsInGround();

    //     // Pass compounds to the new world
    //     passTopologiesToNewWorld(world_i_back);

	// 	// Transfer coordinates from j to i
    //     SimTK::State& state = (worlds[world_i_back].integ)->updAdvancedState();
    //     worlds[world_i_back].setAtomsLocationsInGround_REFAC(state, X_j);

	// 	// Calculate Fixman in replica i back world
	// 	SimTK::Real Fixman = worlds[world_i_back].calcFixman();

	// 	// Transfer buffer coordinates of replica i
	// 	// back to back world
	// 	restoreReplicaCoordinatesToBackWorld(replica_i);

	// 	passTopologiesToNewWorld(world_i_front);

	// 	// Return
	// 	return Fixman;
	// }

}


// Calculate Fixman potential of replica J in replica I's back world. Ui(X_j)
SimTK::Real Context::calcFixman_JinI(int replica_i, int replica_j)
{
	SimTK_ASSERT_ALWAYS(true, "Context::calcFixman_JinI: Not implemented yet.");

	// SimTK::Real U_i = replicas[replica_i].getFixman();

	// if (replica_i == replica_j){ // same replica
	// 	return U_i;

	// }else if(U_i <= 0.0000001){ // fully flexible world
	// 	return U_i;

	// }else{

	// 	// Get replica i thermodynamic state
    //    	int thermoState_i = replica2ThermoIxs[replica_i];

	// 	// Get replica i back world
	// 	int world_i_front = thermodynamicStates[thermoState_i].getWorldIndexes().front();
    //    	int world_i_back = thermodynamicStates[thermoState_i].getWorldIndexes().back();

	// 	// Get coordinates from replica j
	// 	const std::vector<std::vector<
    //         std::pair <Atom *, SimTK::Vec3>>>&
	// 	X_j = replicas[replica_j].getAtomsLocationsInGround();

    //     // Pass compounds to the new world
    //     passTopologiesToNewWorld(world_i_back);

	// 	// Transfer coordinates from j to i
    //     SimTK::State& state = (worlds[world_i_back].integ)->updAdvancedState();
    //     worlds[world_i_back].setAtomsLocationsInGround_REFAC(state, X_j);

	// 	// Calculate Fixman in replica i back world
	// 	SimTK::Real Fixman = worlds[world_i_back].calcFixman();

	// 	// Transfer buffer coordinates of replica i back to back world
	// 	restoreReplicaCoordinatesToBackWorld(replica_i);

	// 	passTopologiesToNewWorld(world_i_front);

	// 	// Return
	// 	return Fixman;
	// }

}

// Calculate Fixman potential of replica I in replica J's back world. Uj(X_i)
SimTK::Real Context::calcFixman_IinJ(int replica_i, int replica_j)
{
	SimTK::Real U_j = replicas[replica_j].getFixman();

	if (replica_j == replica_i){ // same replica
		return U_j;

	}else if(U_j <= 0.0000001){ // fully flexible world
		return U_j;

	}else{

		// Get replica i thermodynamic state
       	int thermoState_j = replica2ThermoIxs[replica_j];

		// Get replica i back world
		int world_j_front = thermodynamicStates[thermoState_j].getWorldIndexes().front();
       	int world_j_back = thermodynamicStates[thermoState_j].getWorldIndexes().back();

		// Get coordinates from replica j
		const SimTK::Compound::AtomTargetLocations X_i = replicas[replica_i].getAtomsLocationsInGround();

        // Pass compounds to the new world
        passTopologiesToNewWorld(world_j_back);

		// Transfer coordinates from j to i
        SimTK::State& state = (worlds[world_j_back].integ)->updAdvancedState();
        worlds[world_j_back].setAtomsLocationsInGround_REFAC(state, X_i);

		// Calculate Fixman in replica i back world
		SimTK::Real Fixman = worlds[world_j_back].calcFixman();

		// Transfer buffer coordinates of replica i back to back world
		restoreReplicaCoordinatesToBackWorld(replica_j);

		passTopologiesToNewWorld(world_j_front);

		// Return
		return Fixman;
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
	//assert(!"Not implemented");

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
bool Context::attemptREXSwap(int replica_X, int replica_Y)
{
	bool returnValue = false;

	#pragma region convienent_vars

	// Get replicas' thermodynamic states indexes
	int thermoState_C = replica2ThermoIxs[replica_X];
	int thermoState_H = replica2ThermoIxs[replica_Y];

	// Record this attempt
	nofAttemptedSwapsMatrix[thermoState_C][thermoState_H] += 1;
	nofAttemptedSwapsMatrix[thermoState_H][thermoState_C] += 1;

	// For useful functions
	auto genericSampler = worlds[0].updSampler(0);

	// ----------------------------------------------------------------
	// Convenient vars (Ballard-Jarzinski nomenclature)
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real beta_C = thermodynamicStates[thermoState_C].getBeta();
	SimTK::Real beta_H = thermodynamicStates[thermoState_H].getBeta();
	
	SimTK::Real U_X0 = replicas[replica_X].getPotentialEnergy(); // last equil potential
	SimTK::Real U_Y0 = replicas[replica_Y].getPotentialEnergy(); // last equil potential

	SimTK::Real refU_X0 = replicas[replica_X].getReferencePotentialEnergy(); // last equil reference potential
	SimTK::Real refU_Y0 = replicas[replica_Y].getReferencePotentialEnergy(); // last equil reference potential

	SimTK::Real W_X = replicas[replica_X].getWORK(); // work without Jacobian
	SimTK::Real W_Y = replicas[replica_Y].getWORK(); // work without Jacobian

	SimTK::Real U_Xtau = replicas[replica_X].get_WORK_PotentialEnergy_New(); // last non-equil potential
	SimTK::Real U_Ytau = replicas[replica_Y].get_WORK_PotentialEnergy_New(); // last non-equil potential

	SimTK::Real refU_Xtau = replicas[replica_X].get_WORK_ReferencePotentialEnergy_New(); // last non-equil potential
	SimTK::Real refU_Ytau = replicas[replica_Y].get_WORK_ReferencePotentialEnergy_New(); // last non-equil potential

	SimTK::Real lnJac_X = replicas[replica_X].get_WORK_Jacobian(); // non-equil Jacobian
	SimTK::Real lnJac_Y = replicas[replica_Y].get_WORK_Jacobian(); // non-equil Jacobian

	#pragma endregion convienent_vars

	// ----------------------------------------------------------------
	// Reduced potentials X0
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real uC_X0 = beta_C * U_X0; // Replica i reduced potential in state i
	SimTK::Real uH_Y0 = beta_H * U_Y0; // Replica j reduced potential in state j
	SimTK::Real uH_X0 = beta_H * U_X0; // Replica i reduced potential in state j
	SimTK::Real uC_Y0 = beta_C * U_Y0; // Replica j reduced potential in state i

	SimTK::Real ref_uC_X0 = beta_C * refU_X0; // Replica i reduced reference potential in state i
	SimTK::Real ref_uH_Y0 = beta_H * refU_Y0; // Replica j reduced reference potential in state j
	SimTK::Real ref_uH_X0 = beta_H * refU_X0; // Replica i reduced reference potential in state j
	SimTK::Real ref_uC_Y0 = beta_C * refU_Y0; // Replica j reduced reference potential in state i

	// ----------------------------------------------------------------
	// Reduced potential Xtau
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real uC_Xtau = beta_C * U_Xtau; // Replica i reduced potential in state i
	SimTK::Real uH_Ytau = beta_H * U_Ytau; // Replica j reduced potential in state j
	SimTK::Real uH_Xtau = beta_H * U_Xtau; // Replica i reduced potential in state j
	SimTK::Real uC_Ytau = beta_C * U_Ytau; // Replica j reduced potential in state i

	SimTK::Real ref_uC_Xtau = beta_C * refU_Xtau; // Replica i reduced reference potential in state i
	SimTK::Real ref_uH_Ytau = beta_H * refU_Ytau; // Replica j reduced reference potential in state j
	SimTK::Real ref_uH_Xtau = beta_H * refU_Xtau; // Replica i reduced reference potential in state j
	SimTK::Real ref_uC_Ytau = beta_C * refU_Ytau; // Replica j reduced reference potential in state i

	// Get Fixman potential for the work coordinates
	SimTK::Real Fix_Xtau = replicas[replica_X].get_WORK_Fixman();
	SimTK::Real Fix_Ytau = replicas[replica_Y].get_WORK_Fixman();

	// Get Fixman potential for the equilibrium coordinates
	SimTK::Real Fix_X0 = replicas[replica_X].getFixman();
	SimTK::Real Fix_Y0 = replicas[replica_Y].getFixman();

	// Get reduced Fixman potentials
	SimTK::Real fixH_Xtau = beta_H * Fix_Xtau;
	SimTK::Real fixC_Ytau = beta_C * Fix_Ytau;

	SimTK::Real fixC_X0 = beta_C * Fix_X0;
	SimTK::Real fixH_Y0 = beta_H * Fix_Y0;

	// Include the Fixman term if indicated
	SimTK::Real Fix_ii = 0, Fix_jj = 0, Fix_ij = 0, Fix_ji = 0;
	if (swapFixman){

        if (thermoState_C == 0){
			std::cout << "Swap between " << thermoState_C << " and "
				<< thermoState_H << " ";

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

        std::cout << "Uii Ujj Uij Uji " << Fix_ii << " " << Fix_jj
            << " " << Fix_ij << " " << Fix_ji << std::endl;
	}

	// ----------------------------------------------------------------
	// LOGP ENERGY EQUILIBRIUM
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real ETerm_equil    = ref_uH_X0 - ref_uC_X0;
				ETerm_equil   += ref_uC_Y0 - ref_uH_Y0;
				ETerm_equil = -1.0 * ETerm_equil;

	// ----------------------------------------------------------------
	// LOGP ENERGY NON-EQUILIBRIUM
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real ETerm_nonequil  = ref_uH_Xtau - ref_uC_Xtau;
			    ETerm_nonequil += ref_uC_Ytau - ref_uH_Ytau;
				ETerm_nonequil = -1.0 * ETerm_nonequil;

	// ----------------------------------------------------------------
	// LOGP WORK
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Get work from X replica
	//SimTK::Real Work_X = (ref_uH_Xtau - ref_uC_X0) + (fixH_Xtau - fixC_X0) - lnJac_X; // variant 1
	SimTK::Real Work_X = (ref_uH_Xtau - ref_uC_X0) - lnJac_X;                           // variant 2

	// Get work from Y replica
	//SimTK::Real Work_Y = (ref_uC_Ytau - ref_uH_Y0) + (fixC_Ytau - fixH_Y0) - lnJac_Y; // variant 1
	SimTK::Real Work_Y = (ref_uC_Ytau - ref_uH_Y0) - lnJac_Y;                           // variant 2

	// Get total work
	SimTK::Real WTerm = -1.0 * (Work_X + Work_Y);

	// ----------------------------------------------------------------
	// CORRECTION TERM FOR REBAS : probability of choosing
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	#pragma region correctionTerm

	SimTK::Real correctionTerm = 1.0;
	SimTK::Real miu_C = qScaleFactorsMiu.at(thermoState_C);
	SimTK::Real miu_H = qScaleFactorsMiu.at(thermoState_H);
	SimTK::Real std_C = qScaleFactorsStd.at(thermoState_C);
	SimTK::Real std_H = qScaleFactorsStd.at(thermoState_H);
	
	SimTK::Real s_X = qScaleFactors.at(thermoState_C);
	SimTK::Real s_Y = qScaleFactors.at(thermoState_H);
	SimTK::Real s_X_1 = 1.0 / s_X;
	SimTK::Real s_Y_1 = 1.0 / s_Y;

	// Correction term is 1 for now
	SimTK::Real qC_s_X = 1.0, qH_s_Y = 1.0, qH_s_X_1 = 1.0, qC_s_Y_1 = 1.0;
	correctionTerm = (qH_s_X_1 * qC_s_Y_1) / (qC_s_X * qH_s_Y);

	#pragma endregion correctionTerm

	// ----------------------------------------------------------------
	// PRINT
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	bool printTerms = false, printWithoutText = true;
	if (printTerms){
		std::cout << "thermoIxs " << thermoState_C << " " << thermoState_H << std::endl;
		std::cout << "replicaIxs " << replica_X << " " << replica_Y << std::endl;
		std::cout << "bibjwiwj " << beta_C << " " << beta_H << " " << std::endl;
		std::cout << "LiiLjj " << uC_Xtau << " " << uH_Ytau << " "
							   << uH_Xtau << " " << uC_Ytau << std::endl;
		std::cout << "EiiEjj " << uC_X0 << " " << uH_Y0 << " "
							   << uH_X0 << " " << uC_Y0 << std::endl;
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
			<< "REXdetails " << ", " << thermoState_C << ", " << thermoState_H << ", "
			<< replica_X << ", " << replica_Y << ", "
			<< beta_C << ", " << beta_H << ", "

			<< uC_X0 << ", " << uH_Y0 << ", " << uH_X0 << ", " << uC_Y0 << ", "
			<< uC_Xtau << ", " << uH_Ytau << ", " << uH_Xtau << ", " << uC_Ytau << ", "

			<< ref_uC_X0 << ", " << ref_uH_Y0 << ", " << ref_uH_X0 << ", " << ref_uC_Y0 << ", "
			<< ref_uC_Xtau << ", " << ref_uH_Ytau << ", " << ref_uH_Xtau << ", " << ref_uC_Ytau << ", "
			
			<< lnJac_X << ", " << lnJac_Y << ", "
			<< W_X << ", " << W_Y << ", "
			<< ", " << s_X << ", " << s_Y << ", " << s_X_1 << ", " << s_Y_1 << ", "
			<< ", " << qC_s_X << ", " << qH_s_Y << ", " << qH_s_X_1 << ", " << qC_s_Y_1 
			<< ", " << ETerm_equil << ", " << WTerm << ", " << correctionTerm << ", "   
		;

		std::cout << rexDetStream.str();		

	}

	// ----------------------------------------------------------------
	// EVALUATE
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real log_p_accept = 1.0;

	// Calculate log_p_accept
	if(runType == RUN_TYPE::REMC){

		log_p_accept = ETerm_equil ;

	}else if(runType == RUN_TYPE::RENEMC){

		log_p_accept = ETerm_nonequil + std::log(correctionTerm) ;

	}else if( runType == RUN_TYPE::RENE){

		log_p_accept = WTerm + std::log(correctionTerm) ;

	}

	// Draw from uniform distribution
	SimTK::Real unifSample = uniformRealDistribution(randomEngine);

	bool testingMode = false; 

	if(testingMode){
		# pragma region REBAS_TEST
		std::cerr << "WARNING: REX EXCHANGE IN TESTING MODE " << std::endl;

		enum TestingWay {
			ALWAYS_ACCEPT,
			ALWAYS_REJECT};
		
		TestingWay testingWay = TestingWay::ALWAYS_ACCEPT;  				// ALWAYS_REJECT

		if(testingWay == TestingWay::ALWAYS_ACCEPT){
			log_p_accept = 1.0;

		}else if(testingWay == TestingWay::ALWAYS_REJECT){
			log_p_accept = -1.0; // std::exp(log_p_accept) = 0.3678794411714424
			unifSample = 1.0;
		}
			
		# pragma endregion REBAS_TEST
	}

	// Accept
	if((log_p_accept >= 0.0) || (unifSample < std::exp(log_p_accept))){

		if(runType == RUN_TYPE::RENE){

			replicas[replica_X].incrementWorldsNofSamples();
			replicas[replica_Y].incrementWorldsNofSamples();
			thermodynamicStates[thermoState_C].incrementWorldsNofSamples();
			thermodynamicStates[thermoState_H].incrementWorldsNofSamples();

			bool onlyNonequilWorlds = true;
			for(int wIx = 0; wIx < nofWorlds; wIx++){
				if(worlds[wIx].getSampler(0)->getDistortOption() == 0){
					onlyNonequilWorlds = false;
					break;
				}
			}

			if(onlyNonequilWorlds){
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

		if((runType == RUN_TYPE::RENE) || (runType == RUN_TYPE::RENEMC)){
			
			// Update replicas coordinates from work generated coordinates
			set_WORK_CoordinatesAsFinal(replica_X);
			set_WORK_CoordinatesAsFinal(replica_Y);

			// Update replica's energy from work last potential energy
			set_WORK_PotentialAsFinal(replica_X);
			set_WORK_PotentialAsFinal(replica_Y);
		}

		// Swap thermodynamic states
		swapThermodynamicStates(replica_X, replica_Y);
		swapPotentialEnergies(replica_X, replica_Y);
		swapReferencePotentialEnergies(replica_X, replica_Y);

		std::cout << "1" 
		<<", " << unifSample
		<< std::endl << std::endl;

		returnValue = true;

	// Reject
	}else{

		rewindReplica();

		// Return to equilibrium worlds coordinates
		// - no need because it is restored in RunREX
		// Return to equilibrium worlds energies
		// - no need because it is restored in RunREX
		// Don't swap thermodynamics states nor energies

		std::cout << "0"
		<<", " << unifSample 
		<< std::endl << std::endl;

		returnValue = false;
	}

	return returnValue;

}

/*!
 * <!--	Get printing REX swap details -->
*/
void Context::getMsg_RexDetHeader(
	std::stringstream& rexDetHeader)
{

	rexDetHeader << "REXdetails"
		<< ", "<< "thermoState_C" << ", "<< "thermoState_H" << ", "<< "replica_X" << ", "<< "replica_Y" << ", "<< "beta_C" << ", "<< "beta_H"                                                                                                                                                                                                                                                                                           << ", "<< "uC_X0" << ", "<< "uH_Y0" << ", "<< "uH_X0" << ", "<< "uC_Y0"
		<< ", "<< "uC_Xtau" << ", "<< "uH_Ytau" << ", "<< "uH_Xtau" << ", "<< "uC_Ytau"                                                                                                                                                                                                                                                                                                                                                 << ", "<< "ref_uC_X0" << ", "<< "ref_uH_Y0" << ", "<< "ref_uH_X0" << ", "<< "ref_uC_Y0"
		<< ", "<< "ref_uC_Xtau" << ", "<< "ref_uH_Ytau" << ", "<< "ref_uH_Xtau" << ", "<< "ref_uC_Ytau"
		<< ", "<< "lnJac_X" << ", "<< "lnJac_Y"
		<< ", "<< "W_X" << ", "<< "W_Y"
		<< ", "<< "s_X" << ", "<< "s_Y" << ", "<< "s_X_1" << ", "<< "s_Y_1"
		<< ", "<< "qC_s_X" << ", "<< "qH_s_Y" << ", "<< "qH_s_X_1" << ", "<< "qC_s_Y_1"
		<< ", "<< "ETerm_equil" << ", "<< "WTerm" << ", "<< "correctionTerm"
		<< ", "<< "acc" << ", "<< "unif"
	;

}

void Context::mixReplicas(int mixi)
{
	if((mixi % swapEvery) == 0){

		int startFrom = mixi % 2;

		for(int thermoState_i = startFrom; thermoState_i <= (nofThermodynamicStates - 2); thermoState_i += 2){

			int thermoState_j = thermoState_i + 1;

			bool swapped = attemptREXSwap(thermo2ReplicaIxs[thermoState_i], thermo2ReplicaIxs[thermoState_j]);
		}

	}
}

// Load replica's atomLocations into it's front world
int Context::restoreReplicaCoordinatesToFrontWorld(int whichReplica)
{

	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[whichReplica];

	//std::cout <<  "Context::restoreReplicaCoordinatesToFrontWorld thermoIx " << thermoIx << std::endl << std::flush;

	// Get worlds indexes of this thermodynamic state
	std::vector<int> worldIndexes =
		thermodynamicStates[thermoIx].getWorldIndexes();

	//std::cout <<  " worldIndexes[1] " << worldIndexes[1] << std::flush;

	// Set thoermoState front world from replica coordinate buffer
	// Will use worlds own integrator's advanced state
	int currWorldIx;
	currWorldIx = worldIndexes.front();

	SimTK::State& state = worlds[currWorldIx].integ->updAdvancedState();

	//worlds[currWorldIx].setAtomsLocationsInGround(state,
	//	replicas[whichReplica].getAtomsLocationsInGround());
	// std::cout << "Context::restoreReplicaCoordinatesToFrontWorld" << std::endl;
	state = setAtoms_SP_NEW(currWorldIx, state, replicas[whichReplica].getAtomsLocationsInGround());		

	return currWorldIx;

}

/*!
 * <!-- Load replica's atomLocations into it's back world -->
*/
void Context::restoreReplicaCoordinatesToBackWorld(int whichReplica)
{

	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[whichReplica];

	// Get worlds indexes of this thermodynamic state
	std::vector<int> worldIndexes =
		thermodynamicStates[thermoIx].getWorldIndexes();

	// Set thoermoState front world from replica coordinate buffer
	// Will use worlds own integrator's advanced state
	SimTK::State& state =
		(worlds[worldIndexes.back()].integ)->updAdvancedState();

	worlds[worldIndexes.back()].setAtomsLocationsInGround_REFAC(state,
		replicas[whichReplica].getAtomsLocationsInGround());


}

// Stores replica's front world's coordinates into it's atomsLocations
// This should always be a fully flexible world
void Context::storeReplicaCoordinatesFromFrontWorld(int whichReplica)
{

	//std::cout <<  "storeReplicaCoordinatesFromFrontWorld " << whichReplica << ": ";

	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[whichReplica];

	//std::cout <<  " thermoIx " << thermoIx;

	// Get worlds indexes of this thermodynamic state
	std::vector<int> worldIndexes =
		thermodynamicStates[thermoIx].getWorldIndexes();

	//std::cout <<  " worldIndexes[1] " << worldIndexes[1];

	// Update replica atomsLocations from back
	replicas[whichReplica].updAtomsLocationsInGround(
		worlds[worldIndexes.front()].getCurrentAtomsLocationsInGround()
	);

	//std::cout << " worldIndexes.front() " << worldIndexes.front();

	//std::cout << std::endl;
}

// Store first world coordinates into replica's work coords buffer
void Context::store_WORK_CoordinatesFromFrontWorld(int whichReplica)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[whichReplica];

	//std::cout <<  " thermoIx " << thermoIx;

	// Get worlds indexes of this thermodynamic state
	std::vector<int> worldIndexes =
		thermodynamicStates[thermoIx].getWorldIndexes();

	// Update replica atomsLocations from back
	replicas[whichReplica].upd_WORK_AtomsLocationsInGround(
		worlds[worldIndexes.front()].getCurrentAtomsLocationsInGround()
	);

}

// Store front world potential energy into work last energy buffer of the
// replica 
void Context::store_WORK_ReplicaEnergyFromFrontWorldFull(int replicaIx)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[replicaIx];

	// Get the index of the front world
	int frontWorldIx =
		thermodynamicStates[thermoIx].getWorldIndexes().front();

	// Get the front world energy
	SimTK::Real energy =
		//worlds[frontWorldIx].CalcFullPotentialEnergyIncludingRigidBodies();
		worlds[frontWorldIx].calcPotentialEnergy();


	// Set this replica's energy
	replicas[replicaIx].set_WORK_PotentialEnergy_New(energy);

}

// Store any WORK Jacobians contribution from back world
void Context::store_WORK_JacobianFromBackWorld(int replicaIx)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[replicaIx];

    // Get the index of the back world
	int backWorldIx =
    thermodynamicStates[thermoIx].getWorldIndexes().back();

    // Set this replica's WORK Jacobians potential
	
    SimTK::Real jac = (worlds[backWorldIx].updSampler(0))->getDistortJacobianDetLog();
	replicas[replicaIx].set_WORK_Jacobian(jac);
}

// Get energy of the back world and store it in replica thisReplica
void Context::storeReplicaEnergyFromBackWorld(int replicaIx)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[replicaIx];

	// Get the index of the back world
	int backWorldIx =
		thermodynamicStates[thermoIx].getWorldIndexes().back();

	// Get the back world energy
	SimTK::Real energy =
		pHMC((worlds[backWorldIx].samplers[0]))->pe_set +
		pHMC((worlds[backWorldIx].samplers[0]))->fix_set;

	// Set this replica's energy
	replicas[replicaIx].setPotentialEnergy(energy);
}

// Get ennergy of the front world and store it in replica thisReplica
void Context::storeReplicaEnergyFromFrontWorldFull(int replicaIx)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[replicaIx];

	// Get the index of the front world
	int frontWorldIx =
		thermodynamicStates[thermoIx].getWorldIndexes().front();

	// Get the front world energy
	SimTK::Real energy =
		//worlds[frontWorldIx].CalcFullPotentialEnergyIncludingRigidBodies(); // DOESN'T WORK with OPENMM
		worlds[frontWorldIx].calcPotentialEnergy();

	// Add the Fixman potential to the energy (DANGEROUS)
	//energy += pHMC((worlds[backWorldIx].samplers[0]))->fix_set;

	// Set this replica's energy
	replicas[replicaIx].setPotentialEnergy(energy);

}

// Get Fixman of the back world and store it in replica thisReplica
void Context::storeReplicaFixmanFromBackWorld(int replicaIx)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[replicaIx];

    // Get the index of the back world
	int backWorldIx =
    thermodynamicStates[thermoIx].getWorldIndexes().back();

    // Set this replica's Fixman potential
    SimTK::Real U = pHMC((worlds[backWorldIx].samplers[0]))->fix_set;
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
	//std::cout << "iniTemperature set to " << T << std::endl << std::flush;

	// -------------
	// Set samplers parameters for this replica
	// =============
	std::vector<SimTK::Real> replicaTimesteps =
		thermodynamicStates[thisThermoStateIx].getTimesteps();
	std::vector<int> replicaMdsteps =
		thermodynamicStates[thisThermoStateIx].getMdsteps();

	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setTimestep(replicaTimesteps[i], false);
	}

	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setMDStepsPerSample(
			replicaMdsteps[i]);
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
		transferCoordinates_WorldToWorld(frontIx, backIx);
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

		transferCoordinates_WorldToWorld(backWorldIx, frontWorldIx);

		//SimTK::Real cumulDiff_Cart = checkTransferCoordinates_Cart(backWorldIx, frontWorldIx);
		//SimTK::Real cumulDiff_BAT = checkTransferCoordinates_BAT(backWorldIx, frontWorldIx);
		// if(cumulDiff_BAT > 0.001){
		// 	std::cout << "\nBad reconstruction " << cumulDiff_BAT << std::endl;
		// }

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
		std::cout << "[EQ] World " << whichWorld
			<< " generating " << numSamples << " samples." << std::endl;
		validated = worlds[whichWorld].generateSamples(numSamples, worldOutStream, header, verbose);
		std::cout << "[EQ] World " << whichWorld
			<< " generated " << numSamples << " samples." << std::endl;

		// size_t wIx = 1; // We want the U and UDot of the torsional dynamics world
		// if (whichWorld == wIx) {
		// 	SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[wIx].updForceField());
		// 	SimTK::SimbodyMatterSubsystem& matter = *(worlds[wIx].matter);
			
		// 	// auto& someState = (worlds[whichWorld]).integ->updAdvancedState();
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

		// 			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
		// 			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();

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

		// Generate samples
		
		// drl
		#ifdef __DRILLING__ // SCALEQ

            // Get drl data
            const std::vector<std::vector<double>>& drl_bon_Energies = worlds[whichWorld].getEnergies_drl_bon();
            const std::vector<std::vector<double>>& drl_ang_Energies = worlds[whichWorld].getEnergies_drl_ang();
            const std::vector<std::vector<double>>& drl_tor_Energies = worlds[whichWorld].getEnergies_drl_tor();
            const std::vector<std::vector<double>>& drl_n14_Energies = worlds[whichWorld].getEnergies_drl_n14();
            const std::vector<std::vector<double>>& drl_vdw_Energies = worlds[whichWorld].getEnergies_drl_vdw();
            const std::vector<std::vector<double>>& drl_cou_Energies = worlds[whichWorld].getEnergies_drl_cou();

            // validated = worlds[whichWorld].generateSamples(numSamples, worldOutStream){

                //warn("under drilling conditions");

                // Update Robosample bAtomList
                SimTK::State& currentAdvancedState = (worlds[whichWorld]).integ->updAdvancedState();
                (worlds[whichWorld]).updateAtomListsFromSimbody(currentAdvancedState); // Update Robosample bAtomList
                // ''''''''''''''''''''
                // coutspaced("SCALING_BAT init:"); ceolf;
                // replicas[0].calcZMatrixBAT( (worlds[whichWorld]).getAtomsLocationsInGround( (worlds[whichWorld]).integ->updAdvancedState() ));
                // thermodynamicStates[0].PrintZMatrixBAT();
                // ''''''''''''''''''''

                // Reinitialize the sampler
                validated = (worlds[whichWorld]).updSampler(0)->reinitialize(currentAdvancedState, worldOutStream, verbose);

                SimTK::Real pe_beforeScale = (worlds[whichWorld]).forces->getMultibodySystem().calcPotentialEnergy((worlds[whichWorld]).integ->updAdvancedState());

                if(false && ((whichWorld == 3)
                        //&& (std::abs((worlds[whichWorld]).updSampler(0)->QScaleFactor - 1.0) > 0.00001)
                )){ 
                    scout("[SCALING_PES]: before") <<" " << pe_beforeScale << eolf;
                    scout("drl_bon_E"); ceol; PrintCppVector(drl_bon_Energies, 6, "bonE", "bonE");
                    scout("drl_ang_E"); ceol; PrintCppVector(drl_ang_Energies, 6, "angE", "angE");
                    scout("drl_tor_E"); ceol; PrintCppVector(drl_tor_Energies, 6, "torE", "torE");
                    scout("drl_n14_E"); ceol; PrintCppVector(drl_n14_Energies, 6, "n14E", "n14E");
                    scout("drl_vdw_E"); ceol; PrintCppVector(drl_vdw_Energies, 6, "vdwE", "vdwE");
                    scout("drl_cou_E"); ceol; PrintCppVector(drl_cou_Energies, 6, "couE", "couE");
                    std::cout<<std::flush;
                } // __end__ choose a world to print drilling

                auto runSamplingLoop = [&](SimTK::State& state) {
                    for (int sampleIx = 0; sampleIx < numSamples; ++sampleIx) {
                        if (verbose) {
                            worldOutStream << header << " ";
                            (worlds[whichWorld]).updSampler(0)->getMsg_InitialParams(worldOutStream);
                        }
                        validated = (worlds[whichWorld]).updSampler(0)->sample_iteration(state, worldOutStream, verbose) && validated;
                        if (verbose) {worldOutStream << std::endl;}
                    }
                };              

                // GENERATE the requested number of samples
                if ((worlds[whichWorld]).getIsRollFlexibilities()) {
                    for (int mobIntIx = 1; mobIntIx < (worlds[whichWorld]).matter->getNumBodies(); ++mobIntIx) {
                        (worlds[whichWorld]).lockAllMobilizers();
                        const SimTK::MobilizedBody& mobod = (worlds[whichWorld]).matter->getMobilizedBody(SimTK::MobilizedBodyIndex(mobIntIx));
                        mobod.unlock(currentAdvancedState);
                        runSamplingLoop(currentAdvancedState);
                    }
                } else {
                    runSamplingLoop(currentAdvancedState);
                }

                SimTK::Real pe_afterScale = (worlds[whichWorld]).forces->getMultibodySystem().calcPotentialEnergy((worlds[whichWorld]).integ->updAdvancedState());

                // ''''''''''''''''''''
                // coutspaced("SCALING_BAT after:"); ceolf;
                // replicas[0].calcZMatrixBAT( (worlds[whichWorld]).getAtomsLocationsInGround( (worlds[whichWorld]).integ->updAdvancedState() ));
                // thermodynamicStates[0].PrintZMatrixBAT();
                // ''''''''''''''''''''
                if(false && ((whichWorld == 3)
                        //&& (std::abs((worlds[whichWorld]).updSampler(0)->QScaleFactor - 1.0) > 0.00001)
                )){
                    scout("[SCALING_PES]: after") <<" " << pe_afterScale << eolf;
                    scout("drl_bon_E"); ceol; PrintCppVector(drl_bon_Energies, 6, "bonE", "bonE");
                    scout("drl_ang_E"); ceol; PrintCppVector(drl_ang_Energies, 6, "angE", "angE");
                    scout("drl_tor_E"); ceol; PrintCppVector(drl_tor_Energies, 6, "torE", "torE");
                    scout("drl_n14_E"); ceol; PrintCppVector(drl_n14_Energies, 6, "n14E", "n14E");
                    scout("drl_vdw_E"); ceol; PrintCppVector(drl_vdw_Energies, 6, "vdwE", "vdwE");
                    scout("drl_cou_E"); ceol; PrintCppVector(drl_cou_Energies, 6, "couE", "couE");
                    std::cout<<std::flush;
                } // __end__ choose a world to print drilling
                
            // }

		#else

			validated = worlds[whichWorld].generateSamples(numSamples, worldOutStream, header, verbose); // =

		#endif

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
void Context::RunReplicaRefactor_SIMPLE(int mixi, int replicaIx)
{
	# pragma region CONVENIENT_VARS
		// Get thermodynamic state and its' worlds
		Replica& replica = replicas[replicaIx];
		int thermoIx = replica2ThermoIxs[replicaIx];
		ThermodynamicState& thermoState = thermodynamicStates[thermoIx];
		std::vector<int>& thermoWorldIxs = thermoState.updWorldIndexes();
		std::vector<int> & distortOpts = thermoState.getDistortOptions();
		size_t thermoNofWorlds = thermoWorldIxs.size();
		assert((thermoWorldIxs.size() == distortOpts.size()));
	# pragma endregion CONVENIENT_VARS

	replica.updWORK() = 0.0;
	replica.upd_WORK_Jacobian() = 0.0;

	// Loop through all worlds
	for(std::size_t thWCnt = 0; thWCnt < thermoNofWorlds; thWCnt++){

		int wIx  = thermoWorldIxs[thWCnt];
		World& currWorld = worlds[wIx];
		HMCSampler* sampler_p = worlds[wIx].samplers[0].get();
		int distortIx = distortOpts[thWCnt];

		// Transfer coordinates to the next world
		if(thWCnt == 0){
			//std::cout << "Transfer coordinates from replica " << replicaIx << " thermoState " << thermoIx << " to world " << thermoWorldIxs.front() << std::endl;
			transferCoordinates_ReplicaToWorld(replicaIx, thermoWorldIxs.front());
			transferQStatistics(thermoIx, thermoWorldIxs.back(), thermoWorldIxs.front());
		}else{
			std::cout << "Transfer coordinates from world " << thermoWorldIxs[thWCnt - 1] << " to world " << wIx << std::endl;
			transferCoordinates_WorldToWorld(thermoWorldIxs[thWCnt - 1], wIx);
			std::cout << "Transfer Q stats from world " << thermoWorldIxs[thWCnt - 1] << " to world " << wIx << std::endl;
			transferQStatistics(thermoIx, thermoWorldIxs[thWCnt - 1], wIx);
			std::cout << "done" << std::endl;
		}

		// Header
		std::string headerToRunWorld = "REX";
					headerToRunWorld += ", " + std::to_string(replicaIx);
					headerToRunWorld += ", " + std::to_string(thermoIx);
					headerToRunWorld += ", " + std::to_string(wIx);
				
		// Run
		bool validated = true;

		if(false || (wIx == 2) //&& (std::abs(sampler_p->QScaleFactor - 1.0) > 0.00001)
		){
		    //std::cout<<"BMps_means "; PrintCppVector(thermoState.getBMps_means(wIx));
			//worlds[wIx].PrintBATFromSimbody(); // BENDSTRETCH
		}

		std::cout << "Running world " << wIx << " (distortIx=" << distortIx << ") for replica " << replicaIx << " at thermodynamic state " << thermoIx << std::endl;

		validated = RunWorld(wIx, headerToRunWorld ) && validated;

		std::cout << "World " << wIx << " validated: " << validated << std::endl;

		if(MEMDEBUG){stdcout_memdebug("Context::RunReplicaRefactor_SIMPLE 6.5");}

		// Calculate Q statistics
		if(sampler_p->getAcc() == true){
			thermoState.calcQStats(wIx, currWorld.getBMps(), currWorld.getPFrs(), currWorld.getAdvancedQs(), currWorld.getNofSamples());
		}else{
			thermoState.calcQStats(wIx, currWorld.getBMps(), currWorld.getPFrs(), SimTK::Vector(currWorld.getNQs(), SimTK::Real(0)), currWorld.getNofSamples());
		}

		// ======================== EQUILIBRIUM ======================
		if(distortIx == 0){

			replica.updAtomsLocationsInGround(currWorld.getCurrentAtomsLocationsInGround());

			replica.setPotentialEnergy(currWorld.calcPotentialEnergy());

			replica.setFixman(sampler_p->fix_set);

			replica.setReferencePotentialEnergy(OMMRef_calcPotential(replica.getAtomsLocationsInGround(), true, true));

		// ======================== NON-EQUILIBRIUM ======================
		}else{

			replica.updWORK() += currWorld.getWork();  // TODO merge with Jacobians
			replica.upd_WORK_Jacobian() += sampler_p->getDistortJacobianDetLog();

			replica.upd_WORK_AtomsLocationsInGround(currWorld.getCurrentAtomsLocationsInGround());

			replica.set_WORK_PotentialEnergy_New(currWorld.calcPotentialEnergy());
		
			replica.set_WORK_Fixman(sampler_p->fix_set);

			replica.set_WORK_ReferencePotentialEnergy_New(OMMRef_calcPotential(replica.get_WORK_AtomsLocationsInGround(), true, true));

		} // __end__ Non/Equilibrium =======================================

		if(MEMDEBUG){stdcout_memdebug("Context::RunReplicaRefactor_SIMPLE 6.6");}

		# pragma region REBAS_TEST
		const SimTK::State& pdbState = currWorld.integ->updAdvancedState();
		currWorld.updateAtomListsFromSimbody(pdbState);

		// std::string pdbMiddle = pdbPrefix + "." + std::to_string(0) + "." + "s" + std::to_string(thermoIx) + "." + "w" + std::to_string(wIx) + ".";
		// std::cout << "Writing " << pdbMiddle << std::endl;

		// Write pdb
		if(pdbRestartFreq){
			if((mixi % pdbRestartFreq) == 0){
				if(wIx == 0){
					for(int mol_i = 0; mol_i < numMolecules; mol_i++){
						topologies[mol_i].writeAtomListPdb(
							outputDir,
							"/pdbs/sb." + pdbPrefix + "." + std::to_string(mol_i) + "." + "s" + std::to_string(thermoIx) + "." + "w" + std::to_string(wIx) + ".",
								".pdb",
								10,
								mixi);
					}
				}
			}
		}
		# pragma endregion REBAS_TEST

		// Increment the nof samples for replica and thermostate
		replica.incrementWorldsNofSamples(1);
		thermoState.incrementWorldsNofSamples(1);

	} // __end__ Loop through all worlds

	if ((mixi + 1) % printFreq == 0) {
		writeLog(mixi + 1, replicaIx);
		REXLog(mixi + 1, replicaIx);
		std::cout << std::flush;

		int whichDCD = replica2ThermoIxs[replicaIx];
		auto [x, y, z] = replicas[replicaIx].getCoordinates();

		// Convert from nm to Angstrom
		for (auto& coord : x) coord *= 10;
		for (auto& coord : y) coord *= 10;
		for (auto& coord : z) coord *= 10;

		thermodynamicStates[whichDCD].writeDCD(x, y, z);
	}

	replica.incrementNofSamples(1);
	thermoState.incrementNofSamples(1);

	SimTK::State& state = worlds.back().integ->updAdvancedState();
	// replica.calcZMatrixBAT( worlds.back().getAtomsLocationsInGround( state ));

	transferCoordinates_WorldToWorld(thermoWorldIxs.back(), thermoWorldIxs.front());

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
		SimTK::State& currentAdvancedState = worlds[wIx].integ->updAdvancedState();
		const auto& sampler = pHMC((worlds[wIx].samplers[0]));

		const auto temperature = sampler->getTemperature();
		const auto NU = currentAdvancedState.getNU();
		const auto acceptedSteps = sampler->acceptedSteps;
		const auto pe_o = sampler->pe_o;
		const auto pe_n = sampler->pe_n;
		const auto pe_set = sampler->pe_set;
		const auto ke_o = sampler->ke_o;
		const auto ke_n = sampler->ke_n;
		const auto ke_set = sampler->ke_set;
		const auto fix_o = sampler->fix_o;
		const auto fix_n = sampler->fix_n;
		const auto fix_set = sampler->fix_set;
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
                << pe_set << ","
                << ke_o << ","
                << ke_n << ","
				<< ke_set << ","
                << fix_o << ","
                << fix_n << ","
                << fix_set << ","
				<< timestep << ","
				<< mdstep << ","
				<< acc << std::endl;
	}
}


/*!
 * <!-- Run replica exchange protocol -->
*/
void Context::RunREX(int equilRounds, int prodRounds)
{
    if(MEMDEBUG){stdcout_memdebug("Context::RunREX 1");}


	// desk_mass_related
	for (int worldIx = 0; worldIx < worlds.size(); worldIx++) {
		for (auto& atom : atoms) {
			
			World& currWorld = worlds[worldIx];
			//SimTK::DuMMForceFieldSubsystem dumm = *(currWorld.forceField);

			SimTK::DuMM::AtomIndex dAIx = atom.getDuMMAtomIndex();
			SimTK::mdunits::Mass atomMass = atom.getMassInDaltons();

			currWorld.forceField->setDuMMAtomMass(dAIx, atomMass);
		}	
	}

	// Is this necesary =======================================================
	realizeTopology();
	
    if(MEMDEBUG){stdcout_memdebug("Context::RunREX 2");}

	// Allocate space for swap matrices
	allocateSwapMatrices();

	// Initialize replicas
	for (size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++){

		// Set intial parameters
		initializeReplica(replicaIx);

	} // ======================================================================

    if(MEMDEBUG){stdcout_memdebug("Context::RunREX 3");}

	// Print a header =========================================================
	std::stringstream rexOutput;
	rexOutput.str("");

	rexOutput << "REX, " << "replicaIx" << ", " << "thermoIx" << ", " << "wIx" ;

	worlds[0].getSampler(0)->getMsg_Header(rexOutput);
	rexOutput << std::endl;

	getMsg_RexDetHeader(rexOutput);

	std::cout << rexOutput.str() << std::endl;
	// ------------------------------------------------------------------------


	// Internal debug studies
	bool givenTsMode = false;

	// Useful vars
	int nofMixes = requiredNofRounds;
	int currFrontWIx = -1;

	// REPLICA EXCHANGE MAIN LOOP -------------------------------------------->
	for(size_t mixi = 0; mixi < equilRounds + prodRounds; mixi++) {

		//std::cout << " REX batch " << mixi << std::endl;

		// Reset replica exchange pairs vector
		if(runType != RUN_TYPE::DEFAULT){                                                 // (9) + (10)
			if(replicaMixingScheme == ReplicaMixingScheme::neighboring){
				setReplicaExchangePairs(mixi % 2);
			}
		}

		// Update work scale factors
		updThermostatesQScaleFactors(mixi);

		//Print_TRANSFORMERS_Work(); // BENDSTRETCH_5
		// if (mixi >= equilRounds) {
		// 	PrintUDot();
		// }

    	if(MEMDEBUG){stdcout_memdebug("Context::RunREX 4");}

		// SIMULATE EACH REPLICA --------------------------------------------->
		for (size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++){ 

			// Update BAT map for all the replica's world
			updSubZMatrixBATsToAllWorlds(replicaIx);

	    	if(MEMDEBUG){stdcout_memdebug("Context::RunREX 5");}

			// Load the front world
			currFrontWIx = restoreReplicaCoordinatesToFrontWorld(replicaIx);           // (1)

			// Set thermo and simulation parameters for the worlds in this replica
			if (mixi < equilRounds) {
				setReplicasWorldsParameters(replicaIx, true, false, mixi);
			} else {
				setReplicasWorldsParameters(replicaIx, false, true, mixi);
			}

	    	if(MEMDEBUG){stdcout_memdebug("Context::RunREX 6");}

			// ======================== SIMULATE ======================
			//RunReplicaRefactor(mixi, replicaIx);
			RunReplicaRefactor_SIMPLE(mixi, replicaIx);
				
			// Copy the new timestep and mdstep if we should be adapting
			if (mixi >= equilRounds) {
				std::vector<SimTK::Real> newTimesteps(worlds.size());
				std::vector<int> newMDSteps(worlds.size());

				for (std::size_t i = 0; i < worlds.size(); i++) {
					newTimesteps[i] = worlds[i].getSampler(0)->getTimestep();
					newMDSteps[i] = worlds[i].getSampler(0)->getMDStepsPerSample();
				}

				int thisThermoStateIx = replica2ThermoIxs[replicaIx];
				thermodynamicStates[thisThermoStateIx].setTimesteps(newTimesteps);
				thermodynamicStates[thisThermoStateIx].setMdsteps(newMDSteps);
			}

	    	if(MEMDEBUG){stdcout_memdebug("Context::RunREX 7");}

			// for(const auto wIx : worldIndexes){ // @@@@@@@@@@@@@
			// 	std::cout << "BMps thIx " << replica2ThermoIxs[replicaIx];
			// 	std::cout << " wIx " << wIx << " nq " << worlds[wIx].getNQs() << " wN " << worlds[wIx].getNofSamples() <<" : ";
			// 	worlds[wIx].PrintXBMps();
			// 	std::cout << std::endl;
			// }
			// printQStats(replica2ThermoIxs[replicaIx]); // @@@@@@@@@@@@@

		} // end replicas simulations

		// Mix replicas
		if((runType != RUN_TYPE::DEFAULT) && (nofReplicas != 1)){
			
			mixReplicas(mixi); // check this
			                                               
			PrintNofAcceptedSwapsMatrix();
		}else{

			// float unifSampleDummy = 1.0;
			// std::cout << "0"
			// <<", " << unifSampleDummy
			// << endl << endl;

			PrintNofAcceptedSwapsMatrix();

		}

		this->nofRounds++; 

	} // end rounds

	//PrintNofAttemptedSwapsMatrix();
	PrintNofAcceptedSwapsMatrix();
	//PrintReplicaMaps();

	// foutU.close();
	// foutUDot.close();

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
		Atom& atom0 = atoms[zRow[childPositionInZMat]];
		SimTK::Compound::AtomIndex cAIx = atom0.getCompoundAtomIndex();

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

void Context::PrintNofAcceptedSwapsMatrix(){

	size_t M = nofAcceptedSwapsMatrix.size();

	//std::cout << "Number of accepted swaps matrix:\n";
	for(size_t i = 0; i < M; i++){
		std::cout << "RSM" ;
		for(size_t j = 0; j < M; j++){
			std::cout << ", " << nofAcceptedSwapsMatrix[i][j] ;
		}
		std::cout << "\n";
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

const int Context::getSwapEvery(){
	return swapEvery;
}

void Context::setSwapEvery(const int& n){
	swapEvery = n;
}

void Context::writePdbs(int someIndex, int thermodynamicStateIx)
{

	//for(int world_i = 0; world_i < this->nofWorlds; world_i++){

		// Update bAtomList in Topology
		const SimTK::State& pdbState =
				worlds[worldIndexes.front()].integ->updAdvancedState();
			worlds[worldIndexes.front()].updateAtomListsFromSimbody(pdbState);

		// const SimTK::State& pdbState =
		// 	worlds[world_i].integ->updAdvancedState();
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
 * <!-- Coordinate transfer -->
*/
void Context::transferCoordinates_WorldToWorld(int srcWIx, int destWIx)
{
	// Get advanced states of the integrators
	SimTK::State& lastAdvancedState = worlds[srcWIx].integ->updAdvancedState();
	SimTK::State& currentAdvancedState = worlds[destWIx].integ->updAdvancedState();

	// // Get BAT coordinates
	// calcZMatrixBAT(srcWIx, otherWorldsAtomsLocations);

	//PrintZMatrixTableAndBAT();
	//PrintZMatrixMobods(srcWIx, lastAdvancedState);
	
	// Pass compounds to the new world
	passTopologiesToNewWorld(destWIx);

	// Focus on destination world
	World& destWorld = worlds[destWIx];
	SimTK::State& someState = currentAdvancedState;

	// New setAtomsLocations
	currentAdvancedState = setAtoms_SP_NEW(destWIx, someState, worlds[srcWIx].getAtomsLocationsInGround(lastAdvancedState));

	// SimTK::Real cumulDiff_Cart = checkTransferCoordinates_Cart(srcWIx, destWIx);
	// SimTK::Real cumulDiff_BAT = checkTransferCoordinates_BAT(srcWIx, destWIx, false);
	// scout("checkTransfer:") <<" " << cumulDiff_Cart <<" " << cumulDiff_BAT << eol;
	// if(cumulDiff_BAT > 0.001){
	// 	std::cout << "\nBad reconstruction " << cumulDiff_BAT << std::endl;
	// }
}

/*!
 * <!-- Check coordinate transfer -->
*/
SimTK::Real Context::checkTransferCoordinates_Cart(int srcWIx, int destWIx)
{
	// // Get advanced states of the integrators
	// SimTK::State& srcAdvancedState = worlds[srcWIx].integ->updAdvancedState();
	// SimTK::State& destAdvancedState = worlds[destWIx].integ->updAdvancedState();

	// // Get coordinates from source World
	// const SimTK::Compound::AtomTargetLocations& srcWorldsAtomsLocations = worlds[srcWIx].getAtomsLocationsInGround(srcAdvancedState);

	// // Get coordinates from destintaion World
	// const SimTK::Compound::AtomTargetLocations& destWorldsAtomsLocations = worlds[destWIx].getAtomsLocationsInGround(destAdvancedState);

	// SimTK_ASSERT_ALWAYS(
	// 	srcWorldsAtomsLocations.size() == destWorldsAtomsLocations.size(),
	// 	"Context::checkTransferCoordinates_Cart: Source and destination worlds have different number of atoms.");

	// // Get the max
	// SimTK::Real vMaxComp = 0;
	// SimTK::Real currDiff;d
	// SimTK::Real cumulDiff = 0;

	// for (const auto& atom : atoms) {
	// 	const SimTK::Compound::AtomIndex aIx = atom.getCompoundAtomIndex();
	// 	const SimTK::Vec3& srcPos = srcWorldsAtomsLocations.at(aIx);
	// 	const SimTK::Vec3& destPos = destWorldsAtomsLocations.at(aIx);

	// 	for (int dim = 0; dim < 3; dim++) {
	// 		currDiff = std::abs(srcPos[dim] - destPos[dim]);
	// 		if (currDiff > vMaxComp) {
	// 			vMaxComp = currDiff;
	// 		}
	// 		cumulDiff += currDiff;
	// 	}
	// }

	// scout("Check coords transfer Cart vMaxComp cumulDiff ") 
	// 	<< std::setprecision(10) << std::fixed
	// 	<< vMaxComp << " " << cumulDiff << std::endl;

	// return cumulDiff;
}

/*!
 * <!-- Check coordinate transfer -->
*/
SimTK::Real Context::checkTransferCoordinates_BAT(int srcWIx, int destWIx, bool wantJacobian)
{

	// // Get advanced states of the integrators
	// SimTK::State& srcAdvancedState = worlds[srcWIx].integ->updAdvancedState();
	// SimTK::State& destAdvancedState = worlds[destWIx].integ->updAdvancedState();

	// // Get coordinates from source World
	// const std::vector<std::vector<std::pair<
	// 	Atom *, SimTK::Vec3> > >&
	// 	srcWorldsAtomsLocations =
	// worlds[srcWIx].getAtomsLocationsInGround(srcAdvancedState);

	// // Get coordinates from destintaion World
	// const std::vector<std::vector<std::pair<
	// 	Atom *, SimTK::Vec3> > >&
	// 	destWorldsAtomsLocations =
	// worlds[destWIx].getAtomsLocationsInGround(destAdvancedState);

	// // Iterate molecules
	// int allCnt = 0;
	// int topoIx = 0;

	// // Get locations of this molecule
	// std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> srcAtomTargets;
	// std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> destAtomTargets;
	// for (int i = 0; i < topologies.size(); i++) {
	// 	std::cout << "Context::checkTransferCoordinates_BAT" << std::endl;
	// 	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> temp;
	// 	worlds[srcWIx].extractAtomTargets(i, srcWorldsAtomsLocations, temp);
	// 	srcAtomTargets.insert(temp.begin(), temp.end());
	// }
	// for (int i = 0; i < topologies.size(); i++) {
	// 	std::cout << "Context::checkTransferCoordinates_BAT" << std::endl;
	// 	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> temp;
	// 	worlds[destWIx].extractAtomTargets(i, destWorldsAtomsLocations, temp);
	// 	destAtomTargets.insert(temp.begin(), temp.end());
	// }

	// int rowCnt = 0;
	// SimTK::Real bondLength_src, bondBend_src, bondTorsion_src;
	// SimTK::Real bondLength_des, bondBend_des, bondTorsion_des;
	// SimTK::Real cosTheta_src, cosTheta_des;
	// SimTK::Real vMaxComp = 0;
	// SimTK::Real currDiff;
	// SimTK::Real cumulDiff = 0;

	// for (const auto& row : zMatrixTable) {
		
	// 	bondLength_src = SimTK::NaN; bondLength_des = SimTK::NaN;
	// 	cosTheta_src = SimTK::NaN; cosTheta_des = SimTK::NaN;
	// 	bondBend_src = SimTK::NaN; bondBend_des = SimTK::NaN;
	// 	bondTorsion_src = SimTK::NaN; bondTorsion_des = SimTK::NaN;

	// 	SimTK::Compound::AtomIndex a0_cAIx, a1_cAIx;
		
	// 	// Calculate source bond length
	// 	a0_cAIx = atoms[row[0]].getCompoundAtomIndex();
	// 	a1_cAIx = atoms[row[1]].getCompoundAtomIndex();
		
	// 	SimTK::Vec3 a0loc = findAtomTarget(srcAtomTargets, a0_cAIx);
	// 	SimTK::Vec3 a1loc = findAtomTarget(srcAtomTargets, a1_cAIx);

	// 	SimTK::Vec3 v_a0a1 = a0loc - a1loc;
	// 	bondLength_src = std::sqrt(SimTK::dot(v_a0a1, v_a0a1));

	// 	if(row[2] >= 0){

	// 		SimTK::Compound::AtomIndex a2_cAIx;
	// 		a2_cAIx = atoms[row[2]].getCompoundAtomIndex();
	// 		SimTK::Vec3 a2loc = findAtomTarget(srcAtomTargets, a2_cAIx);

	// 		// Calculate source angle
	// 		UnitVec3 v1(v_a0a1);
	// 		UnitVec3 v2(a2loc - a1loc);

	// 		cosTheta_src = SimTK::dot(v1, v2);
	// 		assert(cosTheta_src < 1.1);
	// 		assert(cosTheta_src > -1.1);
	// 		if (cosTheta_src > 1.0) cosTheta_src = 1.0;
	// 		if (cosTheta_src < -1.0) cosTheta_src = -1.0;
	// 		bondBend_src = std::acos(cosTheta_src);

	// 		if(row[3] >= 0){
	// 			SimTK::Compound::AtomIndex
	// 				a3_cAIx = atoms[row[3]].getCompoundAtomIndex();
	// 			SimTK::Vec3 a3loc = findAtomTarget(srcAtomTargets, a3_cAIx);

	// 			bondTorsion_src = bDihedral(a0loc, a1loc, a2loc, a3loc);
	// 		}

	// 	} // angle

	// 	// Calculate dest bond length
	// 	a0loc = findAtomTarget(destAtomTargets, a0_cAIx);
	// 	a1loc = findAtomTarget(destAtomTargets, a1_cAIx);

	// 	v_a0a1 = a0loc - a1loc;
	// 	bondLength_des = std::sqrt(SimTK::dot(v_a0a1, v_a0a1));

	// 	if(row[2] >= 0){

	// 		SimTK::Compound::AtomIndex a2_cAIx;
	// 		a2_cAIx = atoms[row[2]].getCompoundAtomIndex();
	// 		SimTK::Vec3 a2loc = findAtomTarget(destAtomTargets, a2_cAIx);

	// 		// Calculate source angle
	// 		UnitVec3 v1(v_a0a1);
	// 		UnitVec3 v2(a2loc - a1loc);

	// 		cosTheta_des = SimTK::dot(v1, v2);
	// 		assert(cosTheta_des < 1.1);
	// 		assert(cosTheta_des > -1.1);
	// 		if (cosTheta_des > 1.0) cosTheta_des = 1.0;
	// 		if (cosTheta_des < -1.0) cosTheta_des = -1.0;
	// 		bondBend_des = std::acos(cosTheta_des);

	// 		if(row[3] >= 0){
	// 			SimTK::Compound::AtomIndex
	// 				a3_cAIx = atoms[row[3]].getCompoundAtomIndex();
	// 			SimTK::Vec3 a3loc = findAtomTarget(destAtomTargets, a3_cAIx);

	// 			bondTorsion_des = bDihedral(a0loc, a1loc, a2loc, a3loc);
	// 		}

	// 	} // angle

	// 	if(wantJacobian){
	// 		// Jacobians
	// 		SimTK::Real sinTheta_src = std::sqrt(1.0 - (cosTheta_src*cosTheta_src));
	// 		SimTK::Real sinTheta_des = std::sqrt(1.0 - (cosTheta_des*cosTheta_des));
	// 		SimTK::Real cosPhi_src = std::cos(bondTorsion_src);
	// 		SimTK::Real sinPhi_src = std::sqrt(1.0 - (cosPhi_src*cosPhi_src));
	// 		SimTK::Real cosPhi_des = std::cos(bondTorsion_des);
	// 		SimTK::Real sinPhi_des = std::sqrt(1.0 - (cosPhi_des*cosPhi_des));

	// 		std::vector<std::vector<SimTK::Real>> Jacobian_src(3, std::vector<SimTK::Real>(3));
	// 		std::vector<std::vector<SimTK::Real>> Jacobian_des(3, std::vector<SimTK::Real>(3));

	// 		Jacobian_src[0][0] = sinPhi_src * cosTheta_src;
	// 		Jacobian_src[0][1] = bondLength_src * cosPhi_src * cosTheta_src;
	// 		Jacobian_src[0][2] = -bondLength_src * sinPhi_src * sinTheta_src;

	// 		Jacobian_src[1][0] = sinPhi_src * sinTheta_src;
	// 		Jacobian_src[1][1] = bondLength_src * cosPhi_src * sinTheta_src;
	// 		Jacobian_src[1][2] = -bondLength_src * sinPhi_src * cosTheta_src;

	// 		Jacobian_src[2][0] = cosPhi_src;
	// 		Jacobian_src[2][1] = -bondLength_src * cosPhi_src;
	// 		Jacobian_src[2][2] = 0;		

	// 		Jacobian_des[0][0] = sinPhi_des * cosTheta_des;
	// 		Jacobian_des[0][1] = bondLength_des * cosPhi_des * cosTheta_des;
	// 		Jacobian_des[0][2] = -bondLength_des * sinPhi_des * sinTheta_des;

	// 		Jacobian_des[1][0] = sinPhi_des * sinTheta_des;
	// 		Jacobian_des[1][1] = bondLength_des * cosPhi_des * sinTheta_des;
	// 		Jacobian_des[1][2] = -bondLength_des * sinPhi_des * cosTheta_des;

	// 		Jacobian_des[2][0] = cosPhi_des;
	// 		Jacobian_des[2][1] = -bondLength_des * cosPhi_des;
	// 		Jacobian_des[2][2] = 0;	

	// 		std::vector<SimTK::Real> cartV_src(3, 0.0);
	// 		std::vector<SimTK::Real> BATV_src{bondLength_src, bondBend_src, bondTorsion_src};
	// 		std::vector<SimTK::Real> cartV_des(3, 0.0);
	// 		std::vector<SimTK::Real> BATV_des{bondLength_des, bondBend_des, bondTorsion_des};

	// 		for (size_t m_I = 0; m_I < Jacobian_src.size(); ++m_I) {
	// 			for (size_t m_J = 0; m_J < Jacobian_src[m_I].size(); ++m_J) {
	// 				cartV_src[m_I] += Jacobian_src[m_I][m_J] * BATV_src[m_J];
	// 			}
	// 		}

	// 		for (size_t m_I = 0; m_I < Jacobian_des.size(); ++m_I) {
	// 			for (size_t m_J = 0; m_J < Jacobian_des[m_I].size(); ++m_J) {
	// 				cartV_des[m_I] += Jacobian_des[m_I][m_J] * BATV_des[m_J];
	// 			}
	// 		}

	// 		// scout("checkTransfer BAT") 
	// 		// << " " << getZMatrixTableEntry(rowCnt, 0) << " " << getZMatrixTableEntry(rowCnt, 1)
	// 		// << " " << getZMatrixTableEntry(rowCnt, 2) << " " << getZMatrixTableEntry(rowCnt, 3)
	// 		// 	<< " " << bondLength_src << " " << bondBend_src << " " << bondTorsion_src
	// 		// 	<< " " << bondLength_des << " " << bondBend_des << " " << bondTorsion_des
	// 		// 	<< " " << bondLength_des - bondLength_src << " " << bondBend_des - bondBend_src << " " << bondTorsion_des - bondTorsion_src
	// 		// 	<< " " << cartV_src[0] << " " << cartV_src[1] << " " << cartV_src[2]
	// 		// 	<< " " << cartV_des[0] << " " << cartV_des[1] << " " << cartV_des[2];
	// 		// ceol;

	// 		currDiff = std::abs(cartV_src[0] - cartV_des[0]);
	// 		if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
	// 		cumulDiff += currDiff;
	// 		if(vMaxComp < currDiff){vMaxComp = currDiff;}

	// 		currDiff = std::abs(cartV_src[1] - cartV_des[1]);
	// 		if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
	// 		cumulDiff += currDiff;
	// 		if(vMaxComp < currDiff){vMaxComp = currDiff;}

	// 		currDiff = std::abs(cartV_src[2] - cartV_des[2]);
	// 		if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
	// 		cumulDiff += currDiff;
	// 		if(vMaxComp < currDiff){vMaxComp = currDiff;}

	// 	}else{
	// 		currDiff = std::abs(bondLength_src - bondLength_des);
	// 		if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
	// 		cumulDiff += currDiff;
	// 		if(vMaxComp < currDiff){vMaxComp = currDiff;}

	// 		currDiff = std::abs(bondBend_src - bondBend_des);
	// 		if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
	// 		cumulDiff += currDiff;
	// 		if(vMaxComp < currDiff){vMaxComp = currDiff;}

	// 		currDiff = std::abs(bondTorsion_src - bondTorsion_des);
	// 		if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
	// 		cumulDiff += currDiff;
	// 		if(vMaxComp < currDiff){vMaxComp = currDiff;}

	// 		// scout("checkTransfer BAT") 
	// 		// << " " << getZMatrixTableEntry(rowCnt, 0) << " " << getZMatrixTableEntry(rowCnt, 1)
	// 		// << " " << getZMatrixTableEntry(rowCnt, 2) << " " << getZMatrixTableEntry(rowCnt, 3)
	// 		// 	<< " " << bondLength_src << " " << bondBend_src << " " << bondTorsion_src
	// 		// 	<< " " << bondLength_des << " " << bondBend_des << " " << bondTorsion_des
	// 		// 	<< " " << bondLength_des - bondLength_src << " " << bondBend_des - bondBend_src << " " << bondTorsion_des - bondTorsion_src
	// 		// 	;
	// 		// ceol;			
	// 	}

	// 	if(row[3] == -2){
	// 		topoIx++;
	// 	}

	// 	rowCnt++;

	// } // every zMatrix row

	// scout("Check coords transfer BAT vMaxComp cumulDiff ")
	// 	<< std::setprecision(10) << std::fixed
	// 	<< vMaxComp << " " << cumulDiff << std::endl;

	// return cumulDiff;

}


void Context::transferCoordinates_ReplicaToWorld(int replicaIx, int destWIx)
{
	std::cout << "Context::transferCoordinates_ReplicaToWorld" << std::endl;
	SimTK::State& state = worlds[destWIx].integ->updAdvancedState();	
	state = setAtoms_SP_NEW(destWIx, state, replicas[replicaIx].getAtomsLocationsInGround());	
}


// SP_NEW_TRANSFER ============================================================

/*! <!-- -->
*/
SimTK::State& Context::setAtoms_CompoundsAndDuMM(int destWIx, SimTK::State& someState, const SimTK::Compound::AtomTargetLocations& atomTargets)
{
	// Get destination world
	World& destWorld = worlds[destWIx];

	// Arrays of Transforms
	SimTK::Transform G_X_T;

	// Loop through molecules/topologies
	for(std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++)
	{

		// // 0. VISUALIZER
		// ///////////////////////////////////////////////////////////
		// // Set the decorator
		// if (destWorld.visual == true) {
		// 	destWorld.paraMolecularDecorator->setAtomTargets(
		// 		otherWorldsAtomsLocations[topoIx]);
		// }

		// 1. COMPOUND MATCHDEFAULT
		///////////////////////////////////////////////////////////

		// Convenient vars
		Topology& currTopology = topologies[topoIx];
		int currNAtoms = currTopology.getNumAtoms();

		// Use Molmodel's Compound match functions to set the new conf
		G_X_T = destWorld.setAtoms_Compound_Match(topoIx, atomTargets);

		// 2.1 MORE COMPOUND FOR DUMM
		///////////////////////////////////////////////////////////

		// Get locations for DuMM
		std::vector<SimTK::Vec3> locaationsInMobods(currNAtoms, SimTK::Vec3(0));

		// Set CompoundAtom frameInMobilizedFrame and get loc in mobod
		destWorld.setAtoms_Compound_FramesAndLocsInMobods(topoIx, atomTargets, locaationsInMobods);

		// 2.2 DUMM
		///////////////////////////////////////////////////////////

		// Set atoms' stations on body
		destWorld.setAtoms_SetDuMMStations(topoIx, locaationsInMobods);
	} // every Topology

	return someState;

}

/*! <!--  -->
*/
void Context::setAtoms_XPF_XBM(int wIx) {
	// Get world's Simbody
	SimTK::State& someState = worlds[wIx].integ->updAdvancedState();
	SimTK::SimbodyMatterSubsystem& matter = *worlds[wIx].matter;
	SimTK::DuMMForceFieldSubsystem &dumm = *worlds[wIx].forceField;

	// Iterate molecules
	for (auto& topology : topologies) {
		for (auto& bond : topology.getBonds()) {

			// Get Molmodel compound atom indices
			SimTK::Compound::AtomIndex childAIx = atoms[bond.getChildAtomGlobalIndex()].getCompoundAtomIndex();
			SimTK::Compound::AtomIndex parentAIx = atoms[bond.getParentAtomGlobalIndex()].getCompoundAtomIndex();

			// Get atoms' mobods
			SimTK::MobilizedBodyIndex childAtomMbx = topology.getAtomMobilizedBodyIndexThroughDumm(childAIx, dumm);
			SimTK::MobilizedBody& childAtomMobod = matter.updMobilizedBody(childAtomMbx);
			SimTK::MobilizedBodyIndex parentAtomMbx = topology.getAtomMobilizedBodyIndexThroughDumm(parentAIx, dumm);
			SimTK::MobilizedBody& parentAtomMobod = matter.updMobilizedBody(parentAtomMbx);

			// Get parent body of the child atom's body
			const SimTK::MobilizedBody& parentMobod =  childAtomMobod.getParentMobilizedBody();
			SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

			// Set mobods transforms on flexible joints
			if(bond.getBondMobility(wIx) != SimTK::BondMobility::Mobility::Rigid) {
				// Create transforms for child default inboard frame (XPF) and default outboard frame (XBM)
				SimTK::Transform XPF, XBM;

				// Bound to Ground
				if(parentMobod.isGround()){
					SimTK::Transform G_X_T = topology.getTopLevelTransform();
					SimTK::Transform T_X_base = topology.getTopTransform(SimTK::Compound::AtomIndex(0));

					XPF = G_X_T * T_X_base;
					XBM = SimTK::Transform();
				} else {
					// Get parent-child BondCenters relationship
					SimTK::Transform X_parentBC_childBC = topology.getDefaultBondCenterFrameInOtherBondCenterFrame(childAIx, parentAIx);
					SimTK::Transform X_childBC_parentBC = ~X_parentBC_childBC;

					// Get parent-child BC transform
					SimTK::Transform X_parentAtom_BCpar = topology.calcDefaultBondCenterFrameInParentAtomFrame(parentAIx, childAIx);
					SimTK::Transform X_childAtom_BCchi = topology.calcDefaultBondCenterFrameInChildAtomFrame(parentAIx, childAIx);
					SimTK::Transform X_BCchi_childAtom = ~X_childAtom_BCchi;

					// Get Top frame
					SimTK::Transform T_X_root = topology.getTopTransform(childAIx);

					// Get Top to parent frame
					const std::pair<int, SimTK::Compound::AtomIndex>& topoAtomPair = worlds[wIx].getMobodRootAtomIndex(parentMbx);
					SimTK::Compound::AtomIndex parentMobodAIx = topoAtomPair.second;

					//SimTK::Compound::AtomIndex parentRootAIx = worlds[wIx].getMbx2aIx()[parentMbx];
					SimTK::Compound::AtomIndex parentRootAIx = parentMobodAIx;
					
					// Origin of the parent mobod
					SimTK::Transform T_X_Proot = topology.getTopTransform(parentRootAIx);
					SimTK::Transform Proot_X_T = ~T_X_Proot;
					SimTK::Transform Proot_X_root = Proot_X_T * T_X_root;

					// // Print parent-child BC transforms
					// std::string bondMbxs = std::to_string(int(parentAtomMbx)) + ":" + std::to_string(int(childAtomMbx));
					// SimTK::Test::PrintTransform(X_parentAtom_BCpar, 6, "parAt_BC:" + bondMbxs, "X_parAt_BC:" + bondMbxs);
					// SimTK::Test::PrintTransform(X_parentBC_childBC, 6, "parBC_chiBC:" + bondMbxs, "X_parBC_chiBC:" + bondMbxs);
					// SimTK::Test::PrintTransform(X_BCchi_childAtom, 6, "BC_chiAt:" + bondMbxs, "BC_chiAt:" + bondMbxs);
					// SimTK::Test::PrintTransform(Proot_X_root, 6, "Proot_X_root:" + bondMbxs, "Proot_X_root:" + bondMbxs);

					// Get inboard dihedral angle
					SimTK::Angle inboardBondDihedralAngle = topology.bgetDefaultInboardDihedralAngle(childAIx);
					SimTK::Transform InboardDihedral_XAxis = SimTK::Rotation(inboardBondDihedralAngle, SimTK::XAxis);
					SimTK::Transform InboardDihedral_ZAxis = SimTK::Rotation(inboardBondDihedralAngle, SimTK::ZAxis);

					// Get inboard bond length
					SimTK::Real inboardBondlength = topology.bgetDefaultInboardBondLength(childAIx);
					SimTK::Transform InboardLength_mZAxis = SimTK::Transform(Rotation(), Vec3(0, 0, -inboardBondlength));

					// Samuel Flores' terminology
					SimTK::Transform M_X_pin = SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis);

					// Get the old PxFxMxB transform
					SimTK::Transform oldX_PB = Proot_X_root;

					// B_X_Ms
					SimTK::Transform B_X_M = X_to_Z; // aka M_X_pin
					SimTK::Transform B_X_M_anglePin = X_parentBC_childBC;
					SimTK::Transform B_X_M_pin 		= X_parentBC_childBC * X_to_Z;
					SimTK::Transform B_X_M_univ 	= X_parentBC_childBC * Y_to_Z;

					// P_X_Fs = old P_X_B * B_X_M
					SimTK::Transform P_X_F 			= Proot_X_root * B_X_M;
					SimTK::Transform P_X_F_anglePin = Proot_X_root * B_X_M_anglePin;
					SimTK::Transform P_X_F_pin 		= Proot_X_root * B_X_M_pin;
					SimTK::Transform P_X_F_univ 	= Proot_X_root * B_X_M;

					//Spherical ===============================================================
					SimTK::Real bondBend = getZMatrixBATValue(6, 1);
					SimTK::Transform XXX;
					//SimTK::Transform XXX(SimTK::Rotation(-1.0 * (bondBend - (SimTK::Pi / 2.0)), SimTK::YAxis));
					SimTK::Transform XXXorthospherical;
					SimTK::Transform XXXinv = ~XXX;

					// Proot -> root -> parentBC -> chilBC=X -> Z
					SimTK::Transform P_X_F_spheric = SimTK::Transform() * Proot_X_root  * X_parentBC_childBC * X_to_Y * Y_to_Z * XXX;
					
					// Z -> X=childBC -> parentBC
					SimTK::Transform M_X_B_spheric = SimTK::Transform() * Z_to_Y * Y_to_X * X_childBC_parentBC;

					SimTK::Transform B_X_M_spheric = ~M_X_B_spheric;

					// OrthoSpherical ==========================================================
					SimTK::Transform P_X_F_orthospheric = X_parentAtom_BCpar; // BAT from Compound
					SimTK::Transform M_X_B_orthospheric = X_parentBC_childBC * X_BCchi_childAtom; // BAT from Compound
					SimTK::Transform B_X_M_orthospheric = ~M_X_B_orthospheric; // X_childAtom_BC * X_childBC_parentBC;

					switch (bond.getBondMobility(wIx)) {
						case SimTK::BondMobility::Mobility::AnglePin:
						case SimTK::BondMobility::Mobility::Slider:
						case SimTK::BondMobility::Mobility::BendStretch:
							XPF = P_X_F_anglePin;
							XBM = B_X_M_anglePin;
							break;

						case SimTK::BondMobility::Mobility::Torsion:
						case SimTK::BondMobility::Mobility::Cylinder:
							XPF = P_X_F_pin;
							XBM = B_X_M_pin;
							break;

						case SimTK::BondMobility::Mobility::BallM:
						case SimTK::BondMobility::Mobility::Rigid:
						case SimTK::BondMobility::Mobility::Translation:
							XPF = P_X_F;
							XBM = B_X_M;
							break;

						case SimTK::BondMobility::Mobility::Spherical:
							XPF = P_X_F_spheric;
							XBM = B_X_M_spheric;
							break;

						case SimTK::BondMobility::Mobility::OrthoSpherical:
							XPF = P_X_F_orthospheric;
							XBM = B_X_M_orthospheric;
							break;

						default:
							warn("Warning: unknown mobility");
							XPF = P_X_F_anglePin;
							XBM = B_X_M_anglePin;
							break;
					}

				}

				childAtomMobod.setDefaultInboardFrame(XPF);
				childAtomMobod.setDefaultOutboardFrame(XBM);
			}

			// Set mobods X_PFs and X_BMs for atoms
			SimTK::Transform G_X_T = topology.getTopLevelTransform();
			
			if(atoms[bond.getChildAtomGlobalIndex()].isRoot()){
				SimTK::Transform T_X_base = topology.getTopTransform(parentAIx);
				SimTK::Transform G_X_base = G_X_T * T_X_base;
				parentAtomMobod.setDefaultInboardFrame(G_X_base);
				parentAtomMobod.setDefaultOutboardFrame(SimTK::Transform());
			}else if(atoms[bond.getParentAtomGlobalIndex()].isRoot()){
				SimTK::Transform T_X_base = topology.getTopTransform(parentAIx);
				SimTK::Transform G_X_base = G_X_T * T_X_base;
				childAtomMobod.setDefaultInboardFrame(G_X_base);
				childAtomMobod.setDefaultOutboardFrame(SimTK::Transform());				
			}
		}
	}
}

/*!
 * <!--  -->
*/
SimTK::State&
Context::setAtoms_MassProperties(
	int wIx
)
{
	SimTK::State& someState = worlds[wIx].integ->updAdvancedState();
	
	SimTK::SimbodyMatterSubsystem& currMatter = *worlds[wIx].matter;
	SimTK::DuMMForceFieldSubsystem &dumm = *worlds[wIx].forceField;
	SimTK::CompoundSystem& compoundSystem = *worlds[wIx].compoundSystem;

	// Set mass properties for mobilized bodies
	// Loop through mobilized bodies
	for (SimTK::MobilizedBodyIndex mbx(1); mbx < currMatter.getNumBodies(); ++mbx)
	{
		SimTK::MobilizedBody& mobod = currMatter.updMobilizedBody(mbx);
		DuMM::ClusterIndex clusterIx = dumm.bgetMobodClusterIndex(mbx);
		SimTK::MassProperties massProperties = dumm.calcClusterMassProperties(clusterIx);
		mobod.setDefaultMassProperties(massProperties);
	}

	// Recover the modified state (may not be necessary)
	compoundSystem.realizeTopology();

	return compoundSystem.updDefaultState();

}


/*!
 * <!--  -->
*/
SimTK::Transform Context::calc_XFM(int wIx, const Topology& topology, SimTK::Compound::AtomIndex& childAIx, SimTK::Compound::AtomIndex& parentAIx, SimTK::BondMobility::Mobility mobility, const SimTK::State& someState) const
{
	// Get world's forcefield and matter
	SimTK::DuMMForceFieldSubsystem &dumm = *worlds[wIx].forceField;
	SimTK::SimbodyMatterSubsystem& matter = *worlds[wIx].matter;

	// Get body and parentBody
	SimTK::MobilizedBodyIndex childMbx = topology.getAtomMobilizedBodyIndexThroughDumm(childAIx, dumm);
	const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(childMbx);
	const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
	SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

	
	//
	SimTK::Real bondBend = getZMatrixBATValue(6, 1);
	//SimTK::Transform XXX(SimTK::Rotation(-1.0 * (bondBend - (SimTK::Pi / 2.0)), SimTK::YAxis));
	SimTK::Transform XXX;
	SimTK::Transform XXXinv = ~XXX;
	SimTK::Transform X_FMspherical = SimTK::Transform()
		* XXXinv
	;

	SimTK::Transform XXX_orthospherical;
	SimTK::Transform XXXinv_orthospherical = ~XXX_orthospherical;
	SimTK::Transform X_FMorthospherical = SimTK::Transform()
		* XXXinv_orthospherical
	;

	// Return
	if(mobility == SimTK::BondMobility::Mobility::Spherical){
		return X_FMspherical;
	}else if(mobility == SimTK::BondMobility::Mobility::OrthoSpherical){
		return X_FMorthospherical;
	}else{
		return Transform();
	}

}


/*!
 * <!--  -->
*/
SimTK::State&
Context::setAtoms_XFM(
	int wIx,
	SimTK::State& someState)
{

	// Get world's forcefield and matter
	SimTK::DuMMForceFieldSubsystem &dumm = *worlds[wIx].forceField;
	SimTK::SimbodyMatterSubsystem& matter = *worlds[wIx].matter;

	// Iterate molecules
	for(const auto& topology : topologies) {

		// Iterate molecule's bonds
		for(const auto& bond : topology.getBonds()) {
			Atom& childAtom  = atoms[bond.getChildAtomGlobalIndex()];
			Atom& parentAtom = atoms[bond.getParentAtomGlobalIndex()];
			//Atom& gparentAtom = atoms[gparentNo];
			//Atom& ggparentAtom = atoms[ggparentNo];			

			int childTopoIx = childAtom.getMoleculeIndex();
			int parentTopoIx = parentAtom.getMoleculeIndex();

			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();
			//SimTK::Compound::AtomIndex gparent_cAIx = gparentAtom.getCompoundAtomIndex();
			//SimTK::Compound::AtomIndex ggparent_cAIx = ggparentAtom.getCompoundAtomIndex();			

			SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
			SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);

			// Get atoms' mobods
			SimTK::MobilizedBodyIndex childMbx = topology.getAtomMobilizedBodyIndexThroughDumm(child_cAIx, dumm);
			SimTK::MobilizedBody& mobod = matter.updMobilizedBody(childMbx);
			const SimTK::MobilizedBody& constParentMobod =  mobod.getParentMobilizedBody();
			SimTK::MobilizedBodyIndex parentMbx = constParentMobod.getMobilizedBodyIndex();
			SimTK::MobilizedBody& parentMobod = matter.updMobilizedBody(parentMbx);

			if(parentMobod.getMobilizedBodyIndex() != 0){ // parent not Ground

				// Calculate
				SimTK::Transform X_FM = calc_XFM(wIx, topology, child_cAIx, parent_cAIx, bond.getBondMobility(wIx), someState);

				//PrintTransform(X_FM, 6, "X_FMreceived");

				//worlds[wIx].compoundSystem->realize(someState, SimTK::Stage::Position);				
				//PrintTransform(mobod.getMobilizerTransform(someState), 6, "X_FMcurrent");

				if(bond.getBondMobility(wIx) == SimTK::BondMobility::Mobility::Spherical)
				{

					// //mobod.getMatterSubsystem().getSystem().getSystemGuts().getVersion();
					// //mobod.getMatterSubsystem().getSystem().getVersion();
					// mobod.getMatterSubsystem().getSystem().realize(someState, SimTK::Stage::Position);
					// //worlds[wIx].compoundSystem->realize(someState, SimTK::Stage::Position);
					// mobod.setQToFitTransform(someState, X_FM);
					// //someState.updQ()[1] = 0.1;
					// worlds[wIx].compoundSystem->realize(
					// 	someState, SimTK::Stage::Position);
					// //PrintTransform(mobod.getMobilizerTransform(someState),
					// //	6, "X_FMafter");
					// scout("mobodQ= ") << mobod.getQAsVector(someState) << eol;
					// scout("stateQ= ") << someState.updQ() << eol;
				
				}
				
				if(bond.getBondMobility(wIx) == SimTK::BondMobility::Mobility::OrthoSpherical){
				}
			}
		} // every bond
	} // every molecule

	return someState;

}


/*! <!-- Transfer geometry to a world
 --> */
SimTK::State& Context::setAtoms_SP_NEW(int destWIx, SimTK::State& someState, const SimTK::Compound::AtomTargetLocations& atomTargets) {

	// Get destination world
	World& destWorld = worlds[destWIx];

	// Get dumm force field and matter
	SimTK::DuMMForceFieldSubsystem &dumm = *destWorld.forceField;
	SimTK::SimbodyMatterSubsystem& matter = *destWorld.matter;

	// Match Compound and DuMM coordinates
	someState = setAtoms_CompoundsAndDuMM(destWIx, someState, atomTargets);

	// Set default child mobod inboard (X_PF) and outboard (X_BM) frames
	// This method only uses internal coordinates per molecule bonds info
	setAtoms_XPF_XBM(destWIx);

	// Recover the modified state (may not be necessary)
	destWorld.compoundSystem->realizeTopology();
	someState = destWorld.compoundSystem->updDefaultState();

	// Set every mobod's mass properties
	someState = setAtoms_MassProperties(destWIx);

	// Set X_FMs
	someState = setAtoms_XFM(destWIx, someState);

	// Realize position
	destWorld.compoundSystem->realize(someState, SimTK::Stage::Position);

	return someState;
}

// SP_NEW_TRANSFER ------------------------------------------------------------

/*!
 * <!--  -->
*/
void Context::passThroughBonds_template(int whichWorld)
{
	// Iterate molecules
	for(const auto& topology : topologies) {
		
		// Iterate molecule's bonds
		for(const auto& bond: topology.getBonds()) {

			Atom& childAtom  = atoms[bond.getChildAtomGlobalIndex()];
			Atom& parentAtom = atoms[bond.getParentAtomGlobalIndex()];

			// Get Compound atom indexes
			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();
			
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


//------------
//------------


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
	SimTK::State& currentAdvancedState = worlds[whichWorld].integ->updAdvancedState();

	const auto NU = currentAdvancedState.getNU();
	const auto acceptedSteps = pHMC((worlds[whichWorld].samplers[0]))->acceptedSteps;
	const auto pe_o = pHMC((worlds[whichWorld].samplers[0]))->pe_o;
	const auto pe_set = pHMC((worlds[whichWorld].samplers[0]))->pe_set;
	const auto ke_o = pHMC((worlds[whichWorld].samplers[0]))->ke_o;
	const auto ke_n = pHMC((worlds[whichWorld].samplers[0]))->ke_n;
	const auto fix_o = pHMC((worlds[whichWorld].samplers[0]))->fix_o;
	const auto fix_n = pHMC((worlds[whichWorld].samplers[0]))->fix_n;
	const auto fix_set = pHMC((worlds[whichWorld].samplers[0]))->fix_set;

	// Write to a file instead of stdout
	logFile
		<< std::fixed << std::setprecision(3)
		<< worlds[whichWorld].updSampler(whichSampler)->getTemperature() << " "
		<< std::fixed << std::setprecision(0)
		<< whichWorld << " "
		<< NU << " "
		<< acceptedSteps << " "
		<< std::fixed << std::setprecision(2) << pe_o << " "
		<< pe_set << " "
		<< ke_o << " "
		<< ke_n << " "
		<< fix_o << " "
		<< fix_n << " "
		<< fix_set << " ";
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
			// SimTK::State& currentAdvancedState = worlds[whichWorld].integ->updAdvancedState();
			// SimTK::MobilizedBodyIndex mbx3 = topology.getAtomMobilizedBodyIndex(
			//     SimTK::Compound::AtomIndex(dihedralIx[4]) );
			// SimTK::MobilizedBody::Pin& mobod3 = (SimTK::MobilizedBody::Pin&) (worlds[whichWorld].matter->updMobilizedBody(mbx3));

			// //std::cout << mbx3 << std::endl ;
			// //std::cout << currentAdvancedState.getQ() << std::endl;
			// //fprintf(logFile, "%.3f ", currentAdvancedState.getQ()[mbx3] );
			// fprintf(logFile, "%.3f ", mobod3.getQ(currentAdvancedState) );

		}
	}
}

void Context::PrintFreeE2EDist(std::size_t whichWorld, int whichCompound)
{
	const Topology& topology = worlds[whichWorld].getTopology(whichCompound);
	SimTK::State& currentAdvancedState = worlds[whichWorld].integ->updAdvancedState();

	for (const auto& distanceIx : distanceIxs) {
		if( distanceIx[0] == whichWorld ) {

			logFile << std::fixed << std::setprecision(3) << 
				Distance(distanceIx[0], distanceIx[1], 0, distanceIx[2], distanceIx[3]) << std::endl;

			SimTK::MobilizedBodyIndex mbx1 = topology.getAtomMobilizedBodyIndex(
				SimTK::Compound::AtomIndex(distanceIx[2]) );
			SimTK::MobilizedBodyIndex mbx2 = topology.getAtomMobilizedBodyIndex(
				SimTK::Compound::AtomIndex(distanceIx[3]) );
			SimTK::MobilizedBody& mobod1 = worlds[whichWorld].matter->updMobilizedBody(mbx1);
			SimTK::MobilizedBody& mobod2 = worlds[whichWorld].matter->updMobilizedBody(mbx2);
			SimTK::Transform X_PF1 = mobod1.getInboardFrame(currentAdvancedState);
			SimTK::Transform X_PF2 = mobod2.getInboardFrame(currentAdvancedState);
			//SimTK::Transform X_BM1 = mobod1.getOutboardFrame(currentAdvancedState);
			//SimTK::Transform X_BM2 = mobod2.getOutboardFrame(currentAdvancedState);
			SimTK::Transform X_FM1 = mobod1.getMobilizerTransform(currentAdvancedState);
			SimTK::Transform X_FM2 = mobod2.getMobilizerTransform(currentAdvancedState);

			SimTK::Transform deltaX_PF = X_PF2.p() - X_PF1.p();
			logFile << std::fixed << std::setprecision(3) << ((-1 * X_FM1.p()) + deltaX_PF.p() + X_FM2.p()).norm() << std::endl;

			//std::cout << "X_PF1:" << std::endl << X_PF1 << std::endl;
			//std::cout << "X_FM1:" << std::endl << X_FM1 << std::endl;
			//std::cout << "X_BM1:" << std::endl << X_BM1 << std::endl;
			//std::cout << "X_PF2:" << std::endl << X_PF2 << std::endl;
			//std::cout << "X_FM2:" << std::endl << X_FM2 << std::endl;
			//std::cout << "X_BM2:" << std::endl << X_BM2 << std::endl;

			//SimTK::Vec3 a1pos = X_PF1.R() * X_FM1.p();
			//SimTK::Vec3 a2pos = X_PF2.R() * X_FM2.p();
			//fprintf(logFile, "%.3f ",
			//    (a1pos - a2pos).norm() );

		}
	}


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
	SimTK::State& advancedState = worlds[currentWorldIx].integ->updAdvancedState();

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
			worlds[worldIndexes.front()].integ->updAdvancedState();
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

	SimTK::State& state = worlds[whichWorld].integ->updAdvancedState();

	Topology& topology = worlds[whichWorld].updTopology(whichCompound);

	SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[whichWorld].forceField);

	SimTK::SimbodyMatterSubsystem& matter = *(worlds[whichWorld].matter);

	int cAIx_1 = atoms[a1].getCompoundAtomIndex();
	int cAIx_2 = atoms[a2].getCompoundAtomIndex();
	int cAIx_3 = atoms[a3].getCompoundAtomIndex();

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

	SimTK::State& state = worlds[whichWorld].integ->updAdvancedState();

	Topology& topology = worlds[whichWorld].updTopology(whichCompound);

	SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[whichWorld].forceField);

	SimTK::SimbodyMatterSubsystem& matter = *(worlds[whichWorld].matter);

	SimTK::Vec3 a1pos, a2pos, a3pos, a4pos;

	// std::cout <<"Context::Dihedral" 
	// 	<<" | "<< whichWorld <<" "<< whichCompound <<" "<< whichSampler <<" "<< a1 <<" "<< a2 <<" "<< a3 <<" "<< a4
	// 	<<" | "<< atoms[a1].getInName() <<" "<< atoms[a2].getInName() <<" "<< atoms[a3].getInName() <<" "<< atoms[a4].getInName()
	// << std::endl;

	int cAIx_1 = atoms[a1].getCompoundAtomIndex();
	int cAIx_2 = atoms[a2].getCompoundAtomIndex();
	int cAIx_3 = atoms[a3].getCompoundAtomIndex();
	int cAIx_4 = atoms[a4].getCompoundAtomIndex();

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

	SimTK::State& state = worlds[whichWorld].integ->updAdvancedState();

	Topology& topology = worlds[whichWorld].updTopology(whichCompound);

	SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[whichWorld].forceField);

	SimTK::SimbodyMatterSubsystem& matter = *(worlds[whichWorld].matter);

	SimTK::Vec3 a1pos, a2pos;

	int cAIx_1 = atoms[a1].getCompoundAtomIndex();
	int cAIx_2 = atoms[a2].getCompoundAtomIndex();

	a1pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_1)),
		dumm, matter, state);
	a2pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(cAIx_2)),
		dumm, matter, state);

	return (a1pos - a2pos).norm();

}

// Writeble reference to a samplers advanced state
SimTK::State& Context::updAdvancedState(std::size_t whichWorld, std::size_t whichSampler)
{
	return (pHMC(worlds[whichWorld].updSampler(whichSampler))->updTimeStepper()->updIntegrator()).updAdvancedState();
}

// Realize Topology Stage for all the Worlds
void Context::realizeTopology() {
	for(auto& world : worlds) {
		world.getCompoundSystem()->realizeTopology();
	}

}

// Realize Topology Stage for all the Worlds
void Context::realizePosition() {
	for(auto& world : worlds) {
		SimTK::State& someState = world.integ->updAdvancedState();
		world.getCompoundSystem()->realize(someState, SimTK::Stage::Position);
	}

}

SimTK::Real Context::getPotentialEnergy(std::size_t world, std::size_t sampler) const {
	return pHMC((worlds[world].samplers[sampler]))->pe_o;
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
	// 			(topologies[molIx].subAtomList[k]).getCompoundAtomIndex();

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


	// Before adding the membrane, we add the contacts and join them
	// to the appropiate Contact Cliques

	// Each of the flags are formatted as such:
	// CONTACTS_X  int1 int2 , int3 int4 int5
	// X is one of {0,1,2,3}, the ints are atom indices (0-Based)
	// and the comma separates the topologies. In the example given
	// the contacts are set on atoms int1, int2 for topology 0,
	// and on atoms int3, int4 and int5 for topology 1.
	// If the user wishes to skip a topology, then they'd 
	// input "-1" as the only atom index.

	std::vector<std::vector<std::vector<int>>> cliqueAtomIxs;

	
	for (int contactCliqueIx = 0; contactCliqueIx < 4; contactCliqueIx++){

		// Empty vector of prmtop atom indexes
		cliqueAtomIxs.push_back({});
		cliqueAtomIxs[contactCliqueIx].push_back({});

		// Get values for this contactCliqueIx
		std::string contactClique_key = "CONTACTS_";
		contactClique_key.append( std::to_string(contactCliqueIx) );
		const std::vector<std::string>& contactClique_vals = setupReader.get(contactClique_key);
		
		int cur_topology = 0;
		if (contactClique_vals.size() > 2) {

			// Get atom indexes for this clique
			for (const auto& value : contactClique_vals){
				
				if (value == ",") { //TODO: This does not account for 'int1,'. Fix this.
					cliqueAtomIxs[contactCliqueIx].push_back({});
					cur_topology++;
				}
				else {
					cliqueAtomIxs[contactCliqueIx][cur_topology].push_back(std::stoi(value));
				}
			}

			// Check
			if(cur_topology != topologies.size()){
				std::cout << "[WARNING] " 
					<< "Number of topologies in CONTACT_ keys don't match the actual number of topologies\n";
			}

			// Add contact atom indexes for all worlds
			for(unsigned int worldIx = 0; worldIx < getNofWorlds(); worldIx++){
				for (int topologyIx = 0;topologyIx < cliqueAtomIxs[contactCliqueIx].size(); topologyIx++){
					
					worlds[worldIx].addContacts(
							cliqueAtomIxs[contactCliqueIx][topologyIx],
							topologyIx,
							SimTK::ContactCliqueId(contactCliqueIx));
				}
			}
		}
	}

	// Add membrane to all worlds.
	for(unsigned int worldIx = 0; worldIx < getNofWorlds(); worldIx++){
		worlds[worldIx].addMembrane(memZWidth);
	}

	// Print
	std::cout << "\n########## MEMBRANE STATS ##########\n";
	std::cout << "Atom cliques are: \n";
	for (int contactClique=0; contactClique<4; contactClique++){
		for (const auto& topologyIx : cliqueAtomIxs[contactClique]) {
		for (int atomIx : topologyIx) {
			std::cout << atomIx << " ";
		}
		std::cout << " / ";
	}
	std::cout << std::endl;
	}
	std::cout << "########## MEMBRANE STATS ##########\n\n";


	// TODO: Do we need this here (looks like World's buissiness)
	realizeTopology();

}

void Context::setNumThreads(int threads) {
	if (threads < 0) {
		std::cerr << "Invalid number of threads (negative value). Default number (0) of threads will be used." << std::endl;
		numThreads = 0;
	} else {
		numThreads = threads;
	}
}

void Context::setNonbonded(int method, SimTK::Real cutoff) {

	nonbondedMethod = method;

	if (cutoff < 0) {
		std::cerr << "Invalid cutoff requested (negative value). Default cutoff of 1.2 nm will be used instead." << std::endl;
		nonbondedCutoff = 1.2;
    }else{
		nonbondedCutoff = cutoff;
    }
}

void Context::setGBSA(SimTK::Real globalScaleFactor) {
	if (globalScaleFactor < 0 && globalScaleFactor > 1) {
		std::cerr << "Invalid GBSA scale factor (valid range is 0 to 1). Default value of 0.0 will be used instead." << std::endl;
		gbsaGlobalScaleFactor = 0.0;
	}else{
		gbsaGlobalScaleFactor = globalScaleFactor;
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
		std::pair <Atom *, SimTK::Vec3 > > >&
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
	// 	a0_cAIx = atoms[row[0]].getCompoundAtomIndex();
	// 	a1_cAIx = atoms[row[1]].getCompoundAtomIndex();
	// 	SimTK::Vec3 a0loc = findAtomTarget(atomTargets, a0_cAIx);
	// 	SimTK::Vec3 a1loc = findAtomTarget(atomTargets, a1_cAIx);

	// 	SimTK::Vec3 v_a0a1 = a0loc - a1loc;
	// 	SimTK::Real bondLength = std::sqrt(SimTK::dot(v_a0a1, v_a0a1));

	// 	if(row[2] >= 0){

	// 		SimTK::Compound::AtomIndex a2_cAIx;
	// 		a2_cAIx = atoms[row[2]].getCompoundAtomIndex();
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
	// 				a3_cAIx = atoms[row[3]].getCompoundAtomIndex();
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
	// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
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
	// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
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
	// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
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
	// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
	// 		SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();

	// 		SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
	// 		SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);

	// 		// Get child-parent mobods
	// 		SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
			
	// 		SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);
	// 		const SimTK::MobilizedBody &childMobod = world.matter->getMobilizedBody(childMbx);
	// 		SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);
	// 		const SimTK::MobilizedBody &parentMobod = world.matter->getMobilizedBody(parentMbx);

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

		// 	const SimTK::Compound::AtomIndex cAIx = atom.getCompoundAtomIndex();

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

		// 		SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
		// 		SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();

		// 		SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
		// 		SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);



		// 		// Is it a base atom
		// 		// const SimTK::Compound::SingleAtom &parentCompoundAtom = parentAtom.getSingleAtom();
		// 		// SimTK::Compound::AtomIndex parCAIx = parentAtom.getCompoundAtomIndex();
		// 		// const SimTK::Compound::AtomPathName parAtomPathName = parentCompoundAtom.getAtomName(parCAIx);
		// 		if(parentAtom.isRoot() == true){

		// 			scout("bondEntry: -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ");

		// 			// Iterate worlds and get child-parent mobods
		// 			size_t wIx = 0;
		// 			for(auto &world : worlds){
		// 				SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
		// 				SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);
		// 				const SimTK::MobilizedBody &parentMobod = world.matter->getMobilizedBody(parentMbx);
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
// 	const auto& currentAdvancedState = worlds[wIx].integ->updAdvancedState();
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

// 			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
// 			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();

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
	// 		SimTK::State& worldCurrentState = worlds[worldIx].integ->updAdvancedState();
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
