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
Context::Context(SimTK::Real Ti, SimTK::Real Tf, uint32_t argSeed)
{
	// Use a random seed if none is provided
	if (argSeed == 0) {
		std::random_device rd;
		this->seed = rd();
	}else{
		this->seed = argSeed;
	}

	// Set the random seed
	randomEngine = buildRandom32(seed);

	// Set the temperature range
	tempIni = Ti;
	tempFin = Tf;

	// Check the temperature range
	SimTK_ASSERT_ALWAYS(tempIni <= tempFin,
		"Context::Context: Initial temperature must be less than or equal to the final temperature.");
	SimTK_ASSERT_ALWAYS(tempIni > 0,
		"Context::Context: Initial temperature must be greater than 0.");
	SimTK_ASSERT_ALWAYS(tempFin > 0,
		"Context::Context: Final temperature must be greater than 0.");

	// Alert user of CUDA environment variables
	if constexpr (OPENMM_PLATFORM_CUDA) {
		if (SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT").empty()){
			std::cerr << cwar_prefix << "CUDA_ROOT not set." << std::endl;
		} else {
			std::cerr << cinf_prefix << "CUDA_ROOT set to " << SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT") << std::endl;
		}
	}
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

/*!
 * <!-- Print atom's Compound and DuMM indexes -->
*/
void Context::PrintAtomsDebugInfo(void){

	for(size_t prmIx = 0; prmIx < natoms; prmIx++){

		bSpecificAtom& bAtom  = atoms[prmIx];

		int topoIx = bAtom.getMoleculeIndex();

		Topology& topology = topologies[topoIx];

		SimTK::Compound::AtomIndex cAIx = bAtom.getCompoundAtomIndex();

		SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(cAIx);
		
		scout("prmIx cAIx dAIx  ") << prmIx << " " << cAIx << " " << dAIx << std::endl;
	}
}

/*! <!-- Read flexibilities from a file and load it into
 * our vector of BOND_FLEXIBILITY  -->
*/
 std::vector<BOND_FLEXIBILITY>& Context::readFlexibility(
	std::string flexFileFN,
 	std::vector<BOND_FLEXIBILITY>& flexibilities) {

	// Get flexible bonds from file. Numbering starts at 0 in prmtop
	std::ifstream flexF(flexFileFN);
	while(flexF.good()){

		// Get a line
		std::string line;
		std::getline(flexF, line);
		if(!line.empty()){

			// Comment
			if(line.at(0) == '#'){continue;}

			// Get words
			std::istringstream iss(line);
			std::string word;
			std::vector<std::string> lineWords;

			while(iss >> word){
				lineWords.push_back(std::move(word));
			}

			// Check the line
			if(lineWords.size() >= 3 ){
				int index_1 = std::stoi(lineWords[0]);
				int index_2 = std::stoi(lineWords[1]);
				std::string mobility =  lineWords[2];

				flexibilities.push_back({index_1, index_2, mobilityMap[mobility]});
			}
			else{
				scout("Bad flex file format");
			}
		}
	}

	return flexibilities;
 }

/*!
 * <!-- Initialize Context 
 * 1.  Setup general input-output parameters
 * 2.  Construct topologies based on what's read from an AmberReader
 * 3.  Add Worlds 
 * 4.  Add contacts: (Add membrane)
 * 5.  Add samplers
 * 6.  Replica exchange setup
 * 7.  Non-equilibrium setup
 * 8.  BAT and Z-matrix
 * 9.  Binding site
 * 10. Geometry calculations
 * 11  Task spaces
 * 12. Constraints -->
*/
bool Context::initializeFromFile(const std::string &inpFN)
{
	/*
	// __new__
	void __newFunction__(void){
		std::cout << "__newFunction__\n";
	}
	// __end__ __new__
	*/
	//SimTK::Test::__newFunction__();

	// Input
	setupReader.ReadSetup(inpFN);
	setupReader.dump(true);
	if ( !CheckInputParameters(setupReader) ) {return false;}

	// Output
	if (!setOutput(setupReader.get("OUTPUT_DIR")[0])) {return false;}
	this->molDir = GetMoleculeDirectoryShort(setupReader.get("MOLECULES")[0]);
	this->pdbPrefix = (molDir + setupReader.get("SEED")[0]);
	this->restartDir = setupReader.get("RESTART_DIR")[0];
	if(restartDir.length() == 0){restartDir = ".";}
	std::cout << "Molecule directory: " << molDir << "\nRestart directory: " << restartDir << "\n";

	// Rounds
	this->requiredNofRounds = std::stoi(setupReader.get("ROUNDS")[0]);
	this->roundsTillReblock = std::stoi((setupReader.get("ROUNDS_TILL_REBLOCK"))[0]);
	this->pdbRestartFreq = std::stoi(setupReader.get("WRITEPDBS")[0]);
	this->printFreq = std::stoi(setupReader.get("PRINT_FREQ")[0]);

	// // Run type (Needed in Run RunREX, RunReplicaRefactor and attemptREXSwap)
	// RUN_TYPE whatRunType = setRunType(setupReader.get("RUN_TYPE")[0]);
	// std::cout << "Run type " << RUN_TYPE_Str[int(whatRunType)] << std::endl;

	// Construct topologies based on what's read from an AmberReader
	// Build Gmolmodel and Molmodel Graphs
	std::string prmtop = setupReader.get("MOLECULES")[0] + "/" + setupReader.get("PRMTOP")[0];
	std::string inpcrd = restartDir + "/" + setupReader.get("INPCRD")[0] + ".s0" + ".rst7";

	loadAmberSystem(prmtop, inpcrd);
	//scout("Context PrintAtoms.\n"); PrintAtoms();

	// Get Z-matrix indexes table based on InternalCoordinates BONDS	
	calcZMatrixTable(); // PrintZMatrixTable();
	reallocZMatrixBAT();

	// Set nonbonded (needed for Add World)
	setNonbonded(
		std::stoi(setupReader.get("NONBONDED_METHOD")[0]),
		std::stod(setupReader.get("NONBONDED_CUTOFF")[0]));

	// Set GBSA
	setGBSA(std::stod(setupReader.get("GBSA")[0]));

    if(MEMDEBUG){stdcout_memdebug("Context::initializeFromFile 1");}

	// desk_mass_related
	for (auto& topology : (topologies)){

		// Iterate through atoms
		for (auto& atom : topology.subAtomList) {

			// Set masses
			SimTK::Compound::AtomIndex cAIx = atom.getCompoundAtomIndex();
			SimTK::mdunits::Mass atomMass = atom.getMass();
			topology.setAtomMass(cAIx, atomMass);

		}
	}

	// Add Worlds to the  Every World instantiates a:
	// CompoundSystem, SimbodyMatterSubsystem, GeneralForceSubsystem,
	// DuMMForceSubsystem, Integrator, TimeStepper and optionally:
	// DecorationSubsystem, Visualizer, VisuzlizerReporter, ParaMolecularDecorator
	for (int worldIx = 0; worldIx < setupReader.get("WORLDS").size(); worldIx++){

		std::string flexFileFN = setupReader.get("MOLECULES")[0] + "/" + setupReader.get("FLEXFILE")[worldIx];
		std::vector<BOND_FLEXIBILITY> flexibilities = {};
		readFlexibility(flexFileFN, flexibilities);

		addWorld(
			setupReader.get("FIXMAN_TORQUE")[worldIx] == "TRUE",
			std::stoi(setupReader.get("SAMPLES_PER_ROUND")[worldIx]),
			ROOT_MOBILITY::WELD,
			flexibilities,
			setupReader.get("OPENMM")[worldIx] == "TRUE",
			setupReader.get("VISUAL")[worldIx] == "TRUE",
			setupReader.get("VISUAL")[worldIx] == "TRUE" ? SimTK::Real(std::stod(setupReader.get("TIMESTEPS")[worldIx])) : 0);

		flexibilities.clear();
	}

	// desk_mass_related
	// Set DuMM atom indexes into atoms of bAtomList
	worlds[0].setDuMMAtomIndexes(); // REVISE

	#ifdef __DRILLING__
		scout("Correspondence cAIx dAIx\n");
		int ix = -1;
		for(auto atom: atoms){
			ix++;
			spacedcout("aIx cAIx dAIx", ix, atom.getCompoundAtomIndex(), atom.getDuMMAtomIndex(), "\n");
		}
	#endif

	// Victor - take a look
	for (int worldIx = 0; worldIx < worlds.size(); worldIx++) {
		worlds[worldIx].setMyContext(this);
	}

	// // Get how much available memory we have
	// struct sysinfo info;
	// if (sysinfo(&info) == 0) {
    //     std::cout << "Free memory: " << info.freeram / (1024 * 1024) << " MB" << std::endl;
    // }
	// Add contacts: (Add membrane).
	// if (setupReader.get("MEMBRANE")[0] != "ERROR_KEY_NOT_FOUND"){
	// 	addContactImplicitMembrane(std::stof(setupReader.get("MEMBRANE")[0]), setupReader);
	// }

	// Add samplers
	for (int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		worlds[worldIx].addSampler(
			SamplerName::HMC,
			sampleGenerator[setupReader.get("SAMPLERS")[worldIx]],
			integratorName[setupReader.get("INTEGRATORS")[worldIx]],
			thermostateName[setupReader.get("THERMOSTAT")[worldIx]],
			std::stod(setupReader.get("TIMESTEPS")[worldIx]),
			std::stoi(setupReader.get("MDSTEPS")[worldIx]),
			setupReader.find("MDSTEPS_STD") ? std::stoi(setupReader.get("MDSTEPS_STD")[worldIx]) : 0,
			std::stof(setupReader.get("BOOST_TEMPERATURE")[worldIx]),
			std::stoi(setupReader.get("BOOST_MDSTEPS")[worldIx]),
			setupReader.find("DISTORT_OPTION") ? std::stoi(setupReader.get("DISTORT_OPTION")[worldIx]) : 0,
			setupReader.find("FLOW_OPTION") ? std::stoi(setupReader.get("FLOW_OPTION")[worldIx]) : 0,
			setupReader.find("WORK_OPTION") ? std::stoi(setupReader.get("WORK_OPTION")[worldIx]) : 0,
			setupReader.get("FIXMAN_POTENTIAL")[worldIx] == "TRUE"
		);
	} // every world

    if(MEMDEBUG){stdcout_memdebug("Context::initializeFromFile 2");}

	// -- REX: REPLICA EXCHANGE --

	// Run type (Needed in Run RunREX, RunReplicaRefactor and attemptREXSwap)
	RUN_TYPE whatRunType = setRunType(setupReader.get("RUN_TYPE")[0]);
	std::string whatRunTypeStr = RUN_TYPE_Str[int(whatRunType)];
	std::cout << "Run type " << whatRunTypeStr << std::endl;

	if(	(whatRunTypeStr == "REMC")   || 
		(whatRunTypeStr == "RENEMC") ||
		(whatRunTypeStr == "RENE") 
		|| true)
	{

		//int NRepl = this->nofReplicas;

		SetupReader rexReader;

		// Storage for thermodynamic state temperatures
		std::vector<SimTK::Real> temperatures;

		// Storage for sampling details
		std::vector<std::vector<std::string>> rexSamplers;
		std::vector<std::vector<int>> rexDistortOptions;
		std::vector<std::vector<std::string>> rexDistortArgs;
		std::vector<std::vector<int>> rexFlowOptions;
		std::vector<std::vector<int>> rexWorkOptions;

		// Storage for each replica simulation parameters
		std::vector<std::vector<std::string>> rexIntegrators;
		std::vector<std::vector<SimTK::Real>> rexTimesteps;
		std::vector<std::vector<int>> rexWorldIndexes;
		std::vector<std::vector<int>> rexMdsteps;
		std::vector<std::vector<int>> rexSamplesPerRound;

		// Read REX config file
		size_t nofReplicas = rexReader.readREXConfigFile(
			setupReader.get("REX_FILE")[0],
			temperatures,

			rexSamplers,
			rexDistortOptions,
			rexDistortArgs,
			rexFlowOptions,
			rexWorkOptions,
			rexIntegrators,

			rexTimesteps,
			rexWorldIndexes,
			rexMdsteps,
			rexSamplesPerRound);

		// Add replicas
		for(int k = 0; k < nofReplicas; k++){
			inpcrd = restartDir + "/" + setupReader.get("INPCRD")[0] + ".s" + std::to_string(k) + ".rst7";
			inpcrdFNs.push_back(inpcrd);
		}

		for(int k = 0; k < nofReplicas; k++){
			loadAtomsCoordinates(prmtop, inpcrdFNs[k]);
			std::cout << "Loaded coordinates of replica " << k << " from " << inpcrdFNs[k] << std::endl;
			addReplica(k);
		}

		// Add thermodynamic states
		for(int k = 0; k < temperatures.size(); k++){
			addThermodynamicState(k, temperatures[k],

				rexSamplers[k],
				rexDistortOptions[k],
				rexDistortArgs[k],
				rexFlowOptions[k],
				rexWorkOptions[k],
				rexIntegrators[k],

				rexWorldIndexes[k], rexTimesteps[k], rexMdsteps[k]);
		}

		// Set a vector of replica pairs for exchanges
	    exchangePairs.resize(nofReplicas);

		// Consider renaming
		loadReplica2ThermoIxs();

		PrintReplicas();

		// How many Gibbs rounds until replica swpas occurs
		setSwapEvery(std::stoi(setupReader.get("REX_SWAP_EVERY")[0]));

		setSwapFixman(std::stoi(setupReader.get("REX_SWAP_FIXMAN")[0]));

	}

	// END REPLICA EXCHANGE -------------------

    if(MEMDEBUG){stdcout_memdebug("Context::initializeFromFile 3");}

	// Non-equilibrium
	PrepareNonEquilibriumParams_Q();
	setThermostatesNonequilibrium();

	// BAT and Z-matrix
	int firstWIx = 0;
	SimTK::State& lastAdvancedState = worlds[firstWIx].integ->updAdvancedState();

	// Get coordinates from source
	const std::vector<std::vector<std::pair<
	bSpecificAtom *, SimTK::Vec3> > >& firstWorldsAtomsLocations =
		worlds[firstWIx].getAtomsLocationsInGround(lastAdvancedState);

	calcZMatrixBAT(firstWIx, firstWorldsAtomsLocations);
	//PrintZMatrixMobods(firstWIx, lastAdvancedState);

    if(MEMDEBUG){stdcout_memdebug("Context::initializeFromFile 4");}

	for(int k = 0; k < nofReplicas; k++){
		replicas[k].reallocZMatrixBAT();
		replicas[k].calcZMatrixBAT(firstWorldsAtomsLocations);
		//replicas[k].PrintZMatrixBAT();
	}
	
	// They all start with replica 0 coordinates
	// TODO does not work in debug
	for (int worldIx = 0; worldIx < worlds.size(); worldIx++) {
		World& world = worlds[worldIx];

		// Add this worlds BAT coordinates to it's samplers
		addSubZMatrixBATsToWorld(worldIx, 0);
		//scout("Context::initializeFromFile PrintSubZMatrixBAT: ") << eol;
		//world.updSampler(0)->PrintSubZMatrixBAT();
	}

	// Binding site
	// TODO If there is any problem here, it might be because this block was above the previous one
	if (setupReader.get("BINDINGSITE_ATOMS")[0] != "ERROR_KEY_NOT_FOUND" &&
		setupReader.get("BINDINGSITE_MOLECULES")[0] != "ERROR_KEY_NOT_FOUND" &&
		setupReader.get("SPHERE_RADIUS")[0] != "ERROR_KEY_NOT_FOUND") {
		// Set binding site TopologyIx, AtomIx and Sphere Radius
		// Generate Amber-style Atom Lists
		std::vector<std::vector<int>> amberAtomIXs;
		amberAtomIXs.push_back({});
		int cur_topology = 0;
		for (const auto& value : setupReader.get("BINDINGSITE_ATOMS")) {
			if (value == ",") {
				amberAtomIXs.push_back({});
				cur_topology++;
			}
			else {
				amberAtomIXs[cur_topology].push_back(std::stoi(value));
			}
		}
		std::vector<int> topologyIXs;
		for (const auto& value : setupReader.get("BINDINGSITE_MOLECULES")) {
			topologyIXs.push_back(std::stoi(value));
		}
		float sphere_radius = std::stod(setupReader.get("SPHERE_RADIUS")[0]);
		std::cout << "Robosample Sphere Radius: " << sphere_radius << std::endl;
		for (auto& world : worlds) {
			world.setTopologyIXs(topologyIXs);
			world.setAmberAtomIXs(amberAtomIXs);
			HMCSampler* sampler_p = pHMC(world.updSampler(0));
			sampler_p->setSphereRadius(sphere_radius);
		}
	}

	// Geometry calculations
	if(setupReader.get("GEOMETRY")[0] == "TRUE"){
		std::vector<std::size_t> distanceIx;
		std::vector<std::size_t> angleIx;
		std::vector<std::size_t> dihedralIx;
		for(unsigned int i = 0; i < setupReader.get("DISTANCE").size(); i++){
			distanceIx.push_back(atoi(setupReader.get("DISTANCE")[i].c_str()));
		}
		for(unsigned int i = 0; i < setupReader.get("ANGLE").size(); i++){
			angleIx.push_back(atoi(setupReader.get("ANGLE")[i].c_str()));
		}
		for(unsigned int i = 0; i < setupReader.get("DIHEDRAL").size(); i++){
			dihedralIx.push_back(atoi(setupReader.get("DIHEDRAL")[i].c_str()));
		}

		addDistances(distanceIx);
		addAngles(angleIx);
		addDihedrals(dihedralIx);
	}

	// Task spaces
	bool usingTaskSpace = false;
	if(usingTaskSpace){
		addTaskSpacesLS();
	}

	// Constraints
	//context.addConstraints();

	//scout("PrintAtomsDebugInfo") << eol;
	//PrintAtomsDebugInfo();

    if(MEMDEBUG){stdcout_memdebug("Context::initializeFromFile 4");}

	// scout("BiotypeIndex check\n");
	// int ix = -1;
	// for(auto atom: atoms){
	// 	ix++;
	// 	spacedcout("atom vs Compound BiotypeIndex" ,
	// 		atom.getBiotypeIndex(),
	// 		(topologies[atom.getMoleculeIndex()]).getAtomBiotypeIndex(atom.getCompoundAtomIndex()),
	// 	"\n");
	// }

    return true;
}

/*!
 * <!-- Find molecules -->
*/
std::vector<int> Context::findMolecules(const AmberReader& reader) {
	// build adjacency list
	std::vector<std::vector<int>> adjacency(natoms);
	const int nbonds = reader.getNumberBonds();
	for (int i = 0; i < nbonds; i++) {
		const int b0 = reader.getBondsAtomsIndex1(i);
		const int b1 = reader.getBondsAtomsIndex2(i);
		adjacency[b0].push_back(b1);
	}

	// map atoms to molecules
	// moleculeAtomIndices[i][j] points to the jth atom (in prmtop order, [0, natoms)) in the ith molecule
	std::vector<std::vector<int>> moleculeAtomIndices;
	std::vector<bool> visited(natoms, false);
	int numMols = -1;

	// this reduces to counting the number of connected components from an adjacency list
	for (int i = 0; i < natoms; i++) {

		// build graph of the current molecule using bfs
		std::stack<int> nodesToVisit;
		nodesToVisit.push(i);

		// new molecule found
		if (!visited[i]) {
			numMols++;
			moleculeAtomIndices.push_back(std::vector<int>());
		}

		// non recursive dfs
		while (!nodesToVisit.empty()) {
			int atom = nodesToVisit.top();
			nodesToVisit.pop();

			if (visited[atom]) {
				continue;
			}

			visited[atom] = true;
			moleculeAtomIndices[numMols].push_back(atom);
			for (auto to_visit : adjacency[atom]) {
				nodesToVisit.push(to_visit);
			}
		}
	}

	// sort the atoms in each molecule by their index in the prmtop file
	std::vector<int> moleculesBegin;
	for (auto& molecule : moleculeAtomIndices) {
		moleculesBegin.push_back(*std::min(molecule.begin(), molecule.end()));
	}

	return moleculesBegin;
}

/*!
 * <!-- Set coordinates for bSpecific atom list from a Amber rst7 file -->
*/
void Context::loadAtomsCoordinates(const std::string& prmtop, const std::string& inpcrdFN) {

	// Load Amber files
	AmberReader reader;
	reader.readAmberFiles(inpcrdFN, prmtop);

	// Check
	//natoms = reader.getNumberAtoms();
	//atoms.size();

	// Iterate through atoms and set as much as possible from amberReader
	for(int aCnt = 0; aCnt < natoms; aCnt++) {

		// Set coordinates in nm (AMBER uses Angstroms)
		atoms[aCnt].setX(reader.getAtomsXcoord(aCnt) / 10.0);
		atoms[aCnt].setY(reader.getAtomsYcoord(aCnt) / 10.0);
		atoms[aCnt].setZ(reader.getAtomsZcoord(aCnt) / 10.0);
		atoms[aCnt].setCartesians(
			reader.getAtomsXcoord(aCnt) / 10.0,
			reader.getAtomsYcoord(aCnt) / 10.0,
			reader.getAtomsZcoord(aCnt) / 10.0 );

	}
}

/*! <!-- Set Gmolmodel atoms properties from a reader: number, name, element,
 * initial name, force field type, charge, coordinates, mass, LJ parameters.
 * 1-to-1 correspondence between prmtop and Gmolmodel.
 * This does not set anything in Compund or DuMM. -->
*/
void Context::loadAtoms(const AmberReader& reader) {
	// Alloc memory for atoms and bonds list
	natoms = reader.getNumberAtoms();
	atoms.resize(natoms);

	// Iterate through atoms and set as much as possible from amberReader
	for(int aCnt = 0; aCnt < natoms; aCnt++) {
		// Assign an index like in prmtop
		atoms[aCnt].setNumber(aCnt);

		// TODO this should disappear
		atoms[aCnt].setDummAtomClassIndex(SimTK::DuMM::AtomClassIndex(aCnt));
		atoms[aCnt].setChargedAtomTypeIndex(SimTK::DuMM::ChargedAtomTypeIndex(aCnt));

		// This is the name of the atom in the .prmtop file
		// Examples: "O1", "C1", "C2", "H1", "H10"
		const std::string initialName = reader.getAtomsName(aCnt);
		
		// Set element
		const int atomicNumber = reader.getAtomicNumber(aCnt);
 		atoms[aCnt].setAtomicNumber(atomicNumber);
		atoms[aCnt].setElem(elementCache.getSymbolByAtomicNumber(atomicNumber));

		// // // Assign a "unique" name. The generator is however limited.
		// // // Examples: AAAA, AAAB, AAAC, AAAD etc
		//atoms[aCnt].generateName(aCnt);

		// Store the initial name from prmtop
		// Examples: "O1", "C1", "C2", "H1", "H10"
		atoms[aCnt].setInName(initialName);

		// Set force field atom type
		// Examples: "O1", "C1", "C2", "H1", "H10"
		//atoms[i].setFfType(initialName);
		atoms[aCnt].setFfType(reader.getAtomsType(aCnt));

		// Set charge as it is used in Amber
		constexpr SimTK::Real chargeMultiplier = 18.2223;
		atoms[aCnt].setCharge(reader.getAtomsCharge(aCnt) / chargeMultiplier);

		// Set coordinates in nm (AMBER uses Angstroms)
		atoms[aCnt].setX(reader.getAtomsXcoord(aCnt) / 10.0);
		atoms[aCnt].setY(reader.getAtomsYcoord(aCnt) / 10.0);
		atoms[aCnt].setZ(reader.getAtomsZcoord(aCnt) / 10.0);
		atoms[aCnt].setCartesians(
			reader.getAtomsXcoord(aCnt) / 10.0,
			reader.getAtomsYcoord(aCnt) / 10.0,
			reader.getAtomsZcoord(aCnt) / 10.0 );

		// Set mass
		SimTK::Real mass = reader.getAtomsMass(aCnt);
		atoms[aCnt].setMass(mass);

		// Set Lennard-Jones parameters
		atoms[aCnt].setVdwRadius(reader.getAtomsRVdW(aCnt));
		atoms[aCnt].setLJWellDepth(reader.getAtomsEpsilon(aCnt));

		// Set residue name and index
		atoms[aCnt].setResidueName(reader.getResidueLabel(aCnt));
		atoms[aCnt].residueIndex = reader.getResidueIndex(aCnt);

		// Assign an unique atom name
		atoms[aCnt].setName(
			atoms[aCnt].getResidueName() + std::to_string(atoms[aCnt].residueIndex) + "_" +
			//atoms[aCnt].getFftype() + "_" +
			atoms[aCnt].getInName() + " " +
			std::to_string(atoms[aCnt].getNumber())
		);

		// Save
		std::cout << std::fixed << std::setprecision(6);
		std::cout <<"STUDY_Context::loadAtoms" <<" atomicNumber "<< atomicNumber <<" mass "<< mass << std::endl;
		elementCache.addElement(atomicNumber, mass);
	}
}

/*! <!-- Build bBond list and BSpecificAtoms neighbors
 *  1-to-1 correspondence between prmtop and Gmolmodel. -->
 **/
void Context::loadBonds(const AmberReader& reader) {
	
	assert( (!atoms.empty()) &&
		"Context::loadBonds: atom list empty.");

	// Allocate memory for bonds list
	nbonds = reader.getNumberBonds();
	bonds.reserve(nbonds);

	// Iterate bonds and get atom indeces
	// This establishes a 1-to-1 correspondence between prmtop and Gmolmodel
	for(int bCnt = 0; bCnt < nbonds; bCnt++) {
		bBond bond;
		bond.setIndex(bCnt);
		bond.i = reader.getBondsAtomsIndex1(bCnt);
		bond.j = reader.getBondsAtomsIndex2(bCnt);

		bond.setForceEquil( reader.getBondsEqval(bCnt) );
		bond.setForceK( reader.getBondsForceK(bCnt) );

		bonds.push_back(bond);

		atoms[bond.i].addNeighborIndex(bonds[bCnt].j);
		atoms[bond.i].addBondIndex(bCnt);

		atoms[bond.j].addNeighborIndex(bonds[bCnt].i);
		atoms[bond.j].addBondIndex(bCnt);
	}

	// Assign the number of bonds an atom has and set the number of freebonds
	// equal to the number of bonds for now
	for(auto& atom : atoms) {
		atom.setNbonds(atom.getNumBonds());
		atom.setFreebonds(atom.getNumBonds());
	}
}

/*! <!-- Set dumm angle list to be used in World generateDummParams. --> **/
void Context::loadAngles(const AmberReader& reader) {
	dummAngles.reserve(reader.getNumberAngles());

	for (int i = 0; i < reader.getNumberAngles(); i++) {
		DUMM_ANGLE angle;
		angle.first = reader.getAnglesAtomsIndex1(i);
		angle.second = reader.getAnglesAtomsIndex2(i);
		angle.third = reader.getAnglesAtomsIndex3(i);

		angle.k = reader.getAnglesForceK(i);
		angle.equil = reader.getAnglesEqval(i);
		angle.equil = static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * angle.equil));

		dummAngles.push_back(angle);
	}
}

/*! <!-- Set dumm torsion list to be used in World generateDummParams. --> **/
void Context::loadTorsions(const AmberReader& reader) {
	// There are multiple torsions defined with the same four indices
	// This vector shows us where each torsion begins (first) and how long it is (second)
	std::vector<std::pair<int, int>> pairStartAndLens = reader.getPairStartAndLen();

	auto checkBond = [&](int a1, int a2) {
		for(int i = 0; i < nbonds; i++){
			if( (bonds[i]).isThisMe(a1, a2) ) {
				return true;
			}
		}
		return false;
	};

	for(const auto& tor : pairStartAndLens){
		int t = tor.first;
		int periodicity = tor.second;

		// Fill this torsion
		DUMM_TORSION dummTorsion;
		dummTorsion.num = periodicity;

		// Check if a quad of atom indices is a normal dihedral or an improper dihedral, by checking if consecutive atoms are bonded
		// TODO: amberReader getDihedralsAtomsIndex4() is negative if this is an improper, but we never look at it and even abs() it
		dummTorsion.first = reader.getDihedralsAtomsIndex1(tor.first);
		dummTorsion.second = reader.getDihedralsAtomsIndex2(tor.first);
		dummTorsion.third = reader.getDihedralsAtomsIndex3(tor.first);
		dummTorsion.fourth = reader.getDihedralsAtomsIndex4(tor.first);

		if (checkBond(dummTorsion.first, dummTorsion.second) &&
			checkBond(dummTorsion.second, dummTorsion.third) &&
			checkBond(dummTorsion.third, dummTorsion.fourth)) {
			dummTorsion.improper = false;
		} else {
			dummTorsion.improper = true;
		}

		// Define parameters
		if (periodicity >= 1) {
			dummTorsion.period[0] = static_cast<int>(reader.getDihedralsPeriod(t));
			dummTorsion.k[0] = reader.getDihedralsForceK(t);
			dummTorsion.phase[0] = static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * reader.getDihedralsPhase(t)));
		}
		if (periodicity >= 2) {
			dummTorsion.period[1] = static_cast<int>(reader.getDihedralsPeriod(t + 1));
			dummTorsion.k[1] = reader.getDihedralsForceK(t + 1);
			dummTorsion.phase[1] = static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * reader.getDihedralsPhase(t + 1)));
		}
		if (periodicity >= 3) {
			dummTorsion.period[2] = static_cast<int>(reader.getDihedralsPeriod(t + 2));
			dummTorsion.k[2] = reader.getDihedralsForceK(t + 2);
			dummTorsion.phase[2] = static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * reader.getDihedralsPhase(t + 2)));
		}
		if (periodicity >= 4) {
			dummTorsion.period[3] = static_cast<int>(reader.getDihedralsPeriod(t + 3));
			dummTorsion.k[3] = reader.getDihedralsForceK(t + 3);
			dummTorsion.phase[3] = static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * reader.getDihedralsPhase(t + 3)));
		}

		dummTorsions.push_back(dummTorsion);
	}
}


void Context::PrintBond(bBond& bond)
{
	bond.Print();

	scout(" ") << atoms[bond.i].getName() <<" ";
	scout(" ") << atoms[bond.j].getName() <<" ";

	ceol;	
}

void Context::PrintBonds(void){

	for(size_t cnt = 0; cnt < nbonds; cnt++) {
		bonds[cnt].Print();

		scout(" ") << atoms[bonds[cnt].i].getName() <<" ";
		scout(" ") << atoms[bonds[cnt].j].getName() <<" ";

		ceol;
	}

}

int Context::checkBonds(void)
{
	for(size_t cnt = 0; cnt < nbonds; cnt++) {

		bBond& bond = bonds[cnt];
		
		if(bond.getMoleculeIndex() < 0){
			std::cerr << "Bond did not set its molecule index." << std::endl;
			PrintBond(bond);
			return 1;
		}

	}

	return 0;
}

/*! <!-- Create a Comopund for each atom -->
*/
void Context::setAtomsCompounds() {
	for(auto& atom : atoms) {
		const std::string& currAtomName = atom.getName();
		const int atomicNumber = atom.getAtomicNumber();
		const int mass = atom.getMass();

		// Create Compound
		atom.setAtomCompound(
			elementCache.getElement(atomicNumber, mass));

	}
}

/*! <!-- Define biotypes for each atom -->
*/
void Context::addBiotypes() {

	// Set a residue name
	std::string resName = "MOL0"; // TODO: delete

	// Define atoms' Biotypes and BiotypeIndexes with their indeces and names
	int aCnt = -1;
	for(auto& atom : atoms) {
		aCnt++;

		// Define a new biotype
		SimTK::BiotypeIndex biotypeIndex = SimTK::Biotype::defineBiotype(
			elementCache.getElement(atom.getAtomicNumber(), atom.getMass()),
			atom.getNBonds(),
			atom.getResidueName().c_str(),
			atom.getName().c_str(),
			SimTK::Ordinality::Any
		);

		atom.setBiotypeIndex(biotypeIndex);

		// Assign atom's biotype as a composed name: name + force field type
		std::string biotype = atom.getResidueName() + atom.getName() + atom.getFftype();
		atom.setBiotype(biotype);

	}
}

/*!
 * <!-- Load Amber System. Will load:
 * 1.  bSpecificAtom list
 * 2.  bBond list
 * 3.  DUMM_ANGLE list for generateParams 
 * 4.  DUMM_TORSIONlist for generateParams
 * 5.  Compounds for every atoms first
 * 6.  Biotypes to use in force field
 * 7.  BAT and Z-matrix
 * 8.  Molmodel Compound Graph
 * 9.  Match Compounds configurations
 * 10. Geometry calculations
 * 11  Task spaces
 * 12. Constraints -->
*/
void Context::loadAmberSystem(const std::string& prmtop, const std::string& inpcrd) {
	
	// Load Amber files
	AmberReader reader;
	reader.readAmberFiles(inpcrd, prmtop);

	// Read Amber and load Gmolmodel lists
	loadAtoms(reader);     // bSpecificAtom list
	loadBonds(reader);     // bBond list
	loadAngles(reader);    // DUMM_ANGLE list for generateParams
	loadTorsions(reader);  // DUMM_TORSIONlist for generateParams	

	// Calculate InternalCoordinates BONDS / BAT graphs
	calc_Gmolmodel_Graph(); 

	// Construct a Compound for every atom
	setAtomsCompounds();

	// Biotype will be used to look up molecular
	// force field specific parameters for an atom type
	addBiotypes();

	// // Calculate InternalCoordinates BONDS / BAT graphs
	// calc_Gmolmodel_Graph();

	// bBonds to BAT / InternalCoordinates bonds map
	load_BONDS_to_bonds( internCoords.getBonds() );

	// InternalCoordinates -> Molmodel Compound graphs with bondAtom
	build_Molmodel_AcyclicGraphs();

	// Close rings
	addRingClosingBonds_All(); 

	//PrintAtoms();
	//PrintBonds();
	if(checkBonds() != 0){
		std::cerr << "[ERROR] " << "bonds checks failed. Exiting" << std::endl;
		exit(1);
	}

	// Generate Topologies sub_array views (also sort bonds - not BONDS)
	generateTopologiesSubarrays();

	// Match Compounds configurations to atoms Cartesian coords
	matchDefaultConfigurations();

	// Allocate coordinate buffers
	Xs.resize(natoms);
	Ys.resize(natoms);
	Zs.resize(natoms);
}

/*!
 * <!-- This must be called after the worlds are added and flexibility vector
 * is loaded -->
*/
/* void Context::setRootMobilitiesFromFlexFiles(void)
{

	// Allocate space for root mobilities
	for(std::size_t worldIx = 0; worldIx < worlds.size(); worldIx++){
		rootMobilities.push_back({});
		
		for(unsigned int molIx = 0; molIx < getNofMolecules(); molIx++) {
			rootMobilities.back().push_back("Rigid");
		
		}
	}

	// Set root mobilities
	for(std::size_t worldIx = 0; worldIx < worlds.size(); worldIx++) {
		
			for (const auto& flex : worlds[worldIx].getFlexibilities()) {

				if(flex.i == -1){

					assert("Set root mobilities: atom index fault." &&
						((flex.j >= 0) && (flex.j < atoms.size())));

					bSpecificAtom& userRootAtom = atoms[flex.j];

					int molIx = userRootAtom.getMoleculeIndex();
					if(molIx >= getNofMolecules()){
						std::cerr << "Set root mobilities: molecule index fault." << std::endl;
						break;
					}

					rootMobilities[worldIx][molIx] = flex.mobility;

				} // found a root mobility
		} // every flexibility
	} // every world

	// Allocate space for root mobilities
	for(std::size_t worldIx = 0; worldIx < worlds.size(); worldIx++){
		
		for(unsigned int molIx = 0; molIx < getNofMolecules(); molIx++) {
			scout("Set root mobility for world ") << worldIx <<" "
				<< " molecule " << molIx <<" to "
				<< rootMobilities[worldIx][molIx] << eol;
		
		}
	}
} */

void Context::appendDCDReporter(const std::string& filename) {
	traj.createTrajectory(filename, "dcd", natoms, topologies.size());
	wantDCD = true;
}

/*!
 * <!-- Main run function -->
*/
void Context::Run() {
	if(    (getRunType() == RUN_TYPE::DEFAULT) 
		|| (getRunType() == RUN_TYPE::REMC)
		|| (getRunType() == RUN_TYPE::RENEMC)
		|| (getRunType() == RUN_TYPE::RENE)){
		RunREX();

	}else{
		std::cout << "[WARNING] " << "Unknown run type. Running default.\n" ;
		RunREX();
	}

}

std::string Context::GetMoleculeDirectoryShort(const std::string& path) const
{
	std::size_t lastSlashPos = path.find_last_of("/");
	return path.substr(lastSlashPos + 1, path.length());
}

// Check input
bool Context::CheckInputParameters(const SetupReader& setupReader) {

	// Context specific parameters
	if ( !SimTK::Pathname::fileExists(setupReader.get("OUTPUT_DIR")[0]) ) {
		std::cerr << cerr_prefix << setupReader.get("OUTPUT_DIR")[0] << " does not exist" << std::endl;
		return false;
	}

	if ( std::stoi((setupReader.get("ROUNDS"))[0]) < 0 ) {
		std::cerr << cerr_prefix << std::stoi((setupReader.get("ROUNDS"))[0]) << " must be positive" << std::endl;
		return false;
	}

	if ( std::stoi((setupReader.get("ROUNDS_TILL_REBLOCK"))[0]) < 0 ) {
		std::cerr << cerr_prefix << std::stoi((setupReader.get("ROUNDS_TILL_REBLOCK"))[0]) << " must be positive" << std::endl;
		return false;
	}
	
	if ( std::stoll((setupReader.get("SEED"))[0]) < 0 ) {
		std::cerr << cerr_prefix << std::stoll((setupReader.get("SEED"))[0]) << " must be positive" << std::endl;
		return false;
	}

	if(setupReader.find("RANDOM_WORLD_ORDER")){
		if(setupReader.get("RANDOM_WORLD_ORDER").size() == 0){
			std::cerr << cerr_prefix << "Please specify RANDOM_WORLD_ORDER with TRUE/FALSE" << std::endl;
			return false;
		} else {
			isWorldsOrderRandom = ((setupReader.get("RANDOM_WORLD_ORDER")[0] == "TRUE") ? true : false);
		}
	}else{
		std::cerr << cerr_prefix << "Please specify if the world order should be random through RANDOM_WORLD_ORDER keyword" << std::endl;
		return false;
	}

	// Set these numbers first
	std::size_t inpNofWorlds = setupReader.get("WORLDS").size();
	std::size_t inpNofMols = setupReader.get("MOLECULES").size();
	std::size_t inpNofEmbeddedTopologies = inpNofWorlds * inpNofMols;

	if(setupReader.get("ROOTS").size() != inpNofEmbeddedTopologies){
		std::cerr << cerr_prefix << "Must have the same no. of root atoms as the no. of Topologies = nofWorlds x nofMolecules.\n" << std::endl;
		return false;
	}

	if(setupReader.get("ROOT_MOBILITY").size() != setupReader.get("ROOTS").size()){
		std::cerr << cerr_prefix << "Must have the same no. of root mobilities as the no. of root atoms.\n" << std::endl;
		return false;
	}

	// World Samplers specific parameters
	if(setupReader.get("SAMPLERS").size() != inpNofWorlds){
		std::cerr << cerr_prefix << "Must have the same no. of samplers as the no. of worlds.\n" << std::endl;
		return false;
	}

	if(inpNofWorlds > setupReader.get("TIMESTEPS").size()){
		std::cerr << cerr_prefix << "Must have the at least same no. of timesteps as the no. of worlds.\n" << std::endl;
		return false;
	}else{
		for(std::size_t worldIx = 0; worldIx < inpNofWorlds; worldIx++){
			if( std::abs(std::stod(setupReader.get("TIMESTEPS")[worldIx])) <= 0.0000001){
				std::cout << "Warning: timestep for world " << worldIx << " is too small.\n";
			}
		}
	}

	// Molecule specific parameters
	// TODO please stop throwing in constructors
	for(std::size_t molIx = 0; molIx < inpNofMols; molIx++) {
		const auto prmtop = setupReader.get("MOLECULES")[molIx] + "/" + setupReader.get("PRMTOP")[molIx];
		//const auto inpcrd = setupReader.get("MOLECULES")[molIx] + "/" + setupReader.get("INPCRD")[molIx] + ".rst7";

		if(!SimTK::Pathname::fileExists(prmtop)) {
			std::cerr << cerr_prefix << "Molecule " + std::to_string(molIx) + " prmtop not found at " + prmtop << std::endl;
			return false;
		}

		// if(!SimTK::Pathname::fileExists(inpcrd)) {
		// 	std::cerr << cerr_prefix << "Molecule " + std::to_string(molIx) + " inpcrd not found at " + inpcrd << std::endl;
		// 	return false;
		// }
	}

	// Topology specific paramters
	for(std::size_t worldIx = 0; worldIx < inpNofWorlds; worldIx++) {
		for(std::size_t molIx = 0; molIx < inpNofMols; molIx++) {

			if(!SimTK::Pathname::fileExists(setupReader.get("MOLECULES")[molIx] + std::string("/") + setupReader.get("FLEXFILE")[molIx]) ){
				std::cerr << cerr_prefix << "world " + std::to_string(worldIx) + " molecule " + std::to_string(molIx) + " flex not found" << std::endl;
				return false;
			}

		}
	}

    if(setupReader.find("REX_SWAP_FIXMAN")){
		if(setupReader.get("REX_SWAP_FIXMAN").size() == 0){
			std::cerr << cerr_prefix << "The DISTORT_OPTION key is present. Please specify a value" << std::endl;
			return false;
		}else{
            swapFixman = std::stoi(setupReader.get("REX_SWAP_FIXMAN")[0]);
		}
	}

	// // Set the type of simulation and temperatures
	// if ( setupReader.get("RUN_TYPE")[0] == "REMC" ) {
	// 	setRunType(RUN_TYPE::REMC);
	// } else if ( setupReader.get("RUN_TYPE")[0] == "RENEMC" ) {
	// 	setRunType(RUN_TYPE::RENEMC);
	// }else if ( setupReader.get("RUN_TYPE")[0] == "RENE" ) {
	// 	setRunType(RUN_TYPE::RENE);
	// } else {
	// 	tempIni = std::stof(setupReader.get("TEMPERATURE_INI")[0]);
	// 	tempFin = std::stof(setupReader.get("TEMPERATURE_FIN")[0]);
	// 	setRunType(RUN_TYPE::DEFAULT);
	// }
	// std::cout << "[CheckInput] " << "Im doing " << static_cast<int>(runType) << std::endl;

	// If we got here we can set global variables
	// Reserve memory
	nofWorlds = inpNofWorlds;
	nofMols = inpNofMols;
	nofEmbeddedTopologies = inpNofEmbeddedTopologies;

	worlds.reserve(nofWorlds);
	worldIndexes.reserve(nofWorlds);

	return true;
}

// Input molecular files TODO : merge with loadCoordinatesFile
bool Context::loadTopologyFile(std::string topologyFilename)
{
	// function args were std::size_t whichWorld, int whichMolecule, std::string topologyFilename
	std::ifstream file(topologyFilename);
	if(!file){
		std::cout << topologyFilename << " not found." << std::endl;
		return false;
	}
	//topFNs[whichWorld].push_back(topologyFilename);
	topFNs.push_back(topologyFilename);

	nofMols = topFNs.size();

	return true;
}

// Load inpcrd / rst7 file. Input only provides a prefix
bool Context::loadCoordinatesFile(std::string coordinatesFilename)
{
	// function args were std::size_t whichWorld, int whichMolecule, std::string coordinatesFilename
	std::ifstream file(coordinatesFilename) ;
	if(!file){
		std::cout << coordinatesFilename << " not found." << std::endl;
		return false;
	}
	//crdFNs[whichWorld].push_back(coordinatesFilename);
	crdFNs.push_back(coordinatesFilename);
	return true;
}

void Context::PrintCoordinates(
    const std::vector<std::vector
    <std::pair <bSpecificAtom *, SimTK::Vec3>>>& atomsLocations)
{
	for(auto& topology : atomsLocations){
		for(auto& atomCoordinates : topology){
			std::cout
			//<< (atomCoordinates.first)->inName
			<< (atomCoordinates.first)->getCompoundAtomIndex() << " "
			<< atomCoordinates.second[0] << " "
			<< atomCoordinates.second[1] << " "
			<< atomCoordinates.second[2] << std::endl;
		}
	}
}

// TODO merge with loadFlexibleBondsSpecs
bool Context::loadRigidBodiesSpecs(std::size_t whichWorld, int, std::string RBSpecsFN)
{
	// function args were :std::size_t whichWorld, int whichMolecule, std::string RBSpecsFN
	std::ifstream file(RBSpecsFN);
	if(!file){
		std::cout << RBSpecsFN << " not found." << std::endl;
		return false;
	}
	rbSpecsFNs[whichWorld].push_back(RBSpecsFN);

	// Update nofEmbeddedTopologies
	int s = 0;
	for(int i = 0; i < rbSpecsFNs.size(); i++){
		s += rbSpecsFNs[i].size();
	}
	nofEmbeddedTopologies = s;

	return true;
}

// Add flexibility filename to whichWorld row of the flexibility filenames
// matrix
bool Context::loadFlexibleBondsSpecs(
	std::size_t whichWorld,
	std::string flexSpecsFN)
{
	// function args were std::size_t whichWorld, int whichMolecule, std::string flexSpecsFN
	std::ifstream file(flexSpecsFN);
	if(!file){
		std::cout << flexSpecsFN << " not found." << std::endl;
		return false;
	}
	flexSpecsFNs[whichWorld].push_back(flexSpecsFN);
	return true;
}

/*!
 * <!--  -->
*/
void Context::setRegimen (std::size_t whichWorld, int, std::string regimen)
{
	// function args were std::size_t whichWorld, int whichMolecule, std::string regimen
	regimens[whichWorld].push_back(regimen);
}

/*!
 * <!-- Add flexibilities and CompoundSystem model -->
*/
void Context::addWorld(
	bool fixmanTorque,
	int samplesPerRound,
	ROOT_MOBILITY rootMobility,
	const std::vector<BOND_FLEXIBILITY>& flexibilities,
	bool useOpenMM,
	bool visual, SimTK::Real visualizerFrequency) {

	// Create new world and add its index
	worldIndexes.push_back(worldIndexes.size());
	worlds.emplace_back(worldIndexes.back(), nofMols, visual, visualizerFrequency);

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

	rbSpecsFNs.push_back(std::vector<std::string>());
	flexSpecsFNs.push_back(std::vector<std::string>());
	regimens.push_back(std::vector<std::string>());

	// Set seed for random number generators
	worlds.back().setSeed(randomEngine());

	// Propagate root mobility
	worlds.back().setRootMobility(rootMobility);

	// Store the number of worlds
	nofWorlds = worlds.size();

	// Prepare the world for OpenMM
	if (useOpenMM) {
		worlds.back().forceField->setUseOpenMMAcceleration(true);
	}

	// Generate DuMM parameters: DuMM atom types, charged atom types, bond types, angle types and torsion types
	worlds.back().generateDummParams(atoms, bonds, dummAngles, dummTorsions, elementCache);

	// Define the inverse map
	std::map<BondMobility::Mobility, std::string> inverseMobilityMap;

	// Populate the inverse map by swapping keys and values
	for (const auto& pair : mobilityMap) {
		inverseMobilityMap[pair.second] = pair.first;
	}

	// Allocate root mobilities
	rootMobilities.push_back({});
	for(unsigned int molIx = 0; molIx < topologies.size(); molIx++){
		rootMobilities.back().push_back("Rigid");
	}

	for (const auto& flex : flexibilities) {

		if(flex.i == -1){

			assert("Set root mobilities: atom index fault." &&
				((flex.j >= 0) && (flex.j < atoms.size())));

			bSpecificAtom& userRootAtom = atoms[flex.j];

			int molIx = userRootAtom.getMoleculeIndex();
			if(molIx >= getNofMolecules()){
				std::cerr << "Set root mobilities: molecule index fault." << std::endl;
				break;
			}

			std::cout << "Set root mobilities -1=" << flex.i << " molecule " << molIx <<" at atom " << flex.j <<" to " << flex.mobility << std::endl;

			(rootMobilities.back())[molIx] = inverseMobilityMap[flex.mobility];

		} // found a root mobility
	} // every flexibility

	// Set flexibilities
	for (auto& bond : bonds) {
		Topology& topology = topologies[bond.getMoleculeIndex()];
		SimTK::Compound::BondIndex compoundBondIx = bond.getBondIndex();

		// Check if the user set the bond mobility
		BondMobility::Mobility mobility = BondMobility::Mobility::Rigid;
		for (const auto& flex : flexibilities) {
			if (bond.isThisMe(flex.i, flex.j)) {
				mobility = flex.mobility;
				break;
			}
		}

		bond.addBondMobility(mobility);
		topology.setBondMobility(mobility, compoundBondIx);
	}

	// 
	worlds.back().AllocateCoordBuffers(natoms);

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
		modelOneEmbeddedTopology(topologyIx, worldIndexes.back());

		// This is many to one map
		topologies[topologyIx].loadAIx2MbxMap();
	}

	// This is one to many map
	worlds.back().loadMbx2AIxMap();

	// Allocate whatever needed Simbody dependent vectors from World here
	worlds.back().allocateStatsContainers();

}

/** Load molecules based on loaded filenames **/
void Context::AddMolecules(
	int requestedNofMols,
	SetupReader& setupReader
){
	assert(!"Deprecated function. Use AddMoleculesFromSetup instead.");

	// // const int molIx = 0;
	// // auto prmtop = setupReader.get("MOLECULES")[molIx] + "/" + setupReader.get("PRMTOP")[molIx];
	// // auto inpcrd = setupReader.get("MOLECULES")[molIx] + "/" + setupReader.get("INPCRD")[molIx] + ".rst7";

	// // readAmberInput reader;
	// // reader.readAmberFiles(inpcrd, prmtop);
	// // loadAtoms(reader);
	// // loadBonds(reader);
	// // loadAngles(reader);
	// // loadTorsions(reader);
	// // setAtomCompoundTypes();
	// // addBiotypes();



	// topologies.reserve(requestedNofMols);
	// moleculeCount = -1;

	// std::vector<std::string> argRoots = setupReader.get("ROOTS");
	// //std::vector<std::string> argRootMobilities = setupReader.get("ROOT_MOBILITY");

	// //std::vector<readAmberInput> amberReader(requestedNofMols);
	// amberReader.resize(requestedNofMols);

	// // Iterate through topology filenames vector
	// //for(unsigned int molIx = 0; molIx < nofMols; molIx++){
	// for(unsigned int molIx = 0; molIx < requestedNofMols; molIx++){

	// 	// Add filenames to Context filenames vectors
	// 	// This has to be called before Worlds constructors so that
	// 	// reserve will be called for molecules and topologies
	// 	//int nofMols = static_cast<int>(setupReader.get("MOLECULES").size());

	// 	std::string topFN =
	// 		setupReader.get("MOLECULES")[molIx] + std::string("/")
	// 		+ setupReader.get("PRMTOP")[molIx];

	// 	std::string crdFN =
	// 		setupReader.get("MOLECULES")[molIx] + std::string("/")
	// 		+ setupReader.get("INPCRD")[molIx] + ".rst7";

	// 	// make amber reader create multiple molecule
	// 	// each molecule then generate one topology
	// 	// amber: vector<molecule info>

	// 	loadTopologyFile( topFN );

	// 	loadCoordinatesFile( crdFN );

	// 	// Initialize an input reader
	// 	//readAmberInput amberReader;
	// 	amberReader[molIx].readAmberFiles(crdFNs[molIx], topFNs[molIx]);

	// 	// Keep track of the number of molecules
	// 	moleculeCount++; // Used for unique names of molecules // from world
	// 	std::string moleculeName = "MOL" + std::to_string(moleculeCount); 
	// 	roots.emplace_back(std::stoi(argRoots[molIx])); // from world
	// 	//rootMobilities.emplace_back("Pin"); // TODO: move to setflexibilities
	// 	topologies.emplace_back(Topology{moleculeName}); // TODO is this ok?

	// 	// // topologies[molIx].setAtomList(atoms);

	// 	// Set Gmolmodel atoms properties from a reader: number, name, element
	// 	// initial name, force field type, charge, coordinates, mass,
	// 	// LJ parameters
	// 	// Also resizes atomFrameCache
	// 	topologies[molIx].SetGmolAtomPropertiesFromReader(&amberReader[molIx]);

	// 	// Set bonds properties from reader: bond indeces, atom neighbours
	// 	topologies[molIx].SetGmolBondingPropertiesFromReader(&amberReader[molIx]);

	// 	// Set atoms Molmodel types (Compound::SingleAtom derived) based on
	// 	// their valence // from world
	// 	topologies[molIx].SetGmolAtomsCompoundTypes();

	// 	// topologies[molIx].load(amberReader[molIx]);

	// 	// Add Biotype indeces and Biotype names representing Biotypes
	// 	topologies[molIx].bAddBiotypes(//&amberReader[molIx]);
	// 	);
	// 	// topologies[molIx].BAT();

	// 	// Build Robosample graph and Compound graph.
	// 	// It also asigns atom indexes in Compound
	// 	// This is done only once and it needs
	// 	topologies[molIx].buildGraphAndMatchCoords(
	// 		//std::stoi(
	// 			roots.back()
	// 		//)
	// 		);

	// 	// Helper function for calc MBAT determinant
	// 	topologies[molIx].loadTriples();

	// 	// Map of Compound atom indexes to Robosample atom indexes
	// 	topologies[molIx].loadCompoundAtomIx2GmolAtomIx();
	// 	//std::cout << "Topology " << molIx << " info\n";
	// 	//topologies[molIx].printMaps();

	// }


}


// ============================================================================
// ============================================================================
// ==========================   SINGLE PRMTOP    ==============================
// ============================================================================
// ============================================================================
/*!
 * <!-- Sets the root atom for a topology -->
*/
void Context::setRootAtom(Topology& topology, int rootAmberIx)
{
	// Set Topology base atom
	bSpecificAtom& rootAtom = atoms[rootAmberIx];

	rootAtom.setIsRoot(true);

	const SimTK::Compound::SingleAtom& compoundRootAtom
		= rootAtom.getSingleAtom();

	topology.baseAtomNumber = rootAmberIx;

	Transform location = Transform();
	topology.setBaseAtom(compoundRootAtom, location);

	topology.setAtomBiotype(rootAtom.getName(), rootAtom.getResidueName(),
							rootAtom.getName());
	topology.convertInboardBondCenterToOutboard();
	topology.baseSetFlag = 1;

	scout("rootAmberIx root ") 
	// << rootAtom.getName() <<" " << rootAtom.getInName() 
	<< " " << rootAmberIx 
	<< eol;	
}

/*!
 * <!--  -->
*/
void Context::load_BONDS_to_bonds(
	const std::vector<std::vector<BOND>>& BATbonds)
{

	BONDS_to_bonds.resize(BATbonds.size());

	for(unsigned int molIx = 0; molIx < BATbonds.size(); molIx++){
		BONDS_to_bonds[molIx].resize(BATbonds[molIx].size(), -111111);
	}

	bonds_to_BONDS.resize(bonds.size());

	// bBonds and BAT bonds equivalence
	scout(" bonds BATbonds euqivalence ") << eol;
	bool found = false;

	for(unsigned int molIx = 0; molIx < BATbonds.size(); molIx++){

		for(size_t cnt = 0; cnt < BATbonds[molIx].size(); cnt++){

			found = false;
			for(size_t bbCnt = 0; bbCnt < bonds.size(); bbCnt++){

				if(	(BATbonds[molIx][cnt].first  == bonds[bbCnt].j) &&
					(BATbonds[molIx][cnt].second == bonds[bbCnt].i)){
					found = true;

				}else if((BATbonds[molIx][cnt].first  == bonds[bbCnt].i) &&
						(BATbonds[molIx][cnt].second == bonds[bbCnt].j)){
					found = true;

				}

				if(found){
					BONDS_to_bonds[molIx][cnt] = bbCnt;
					bonds_to_BONDS[bbCnt] = std::make_pair(molIx, cnt);
					break;
				}

			}
		}
	}

	// scout("BONDS_to_bonds");
	// for(unsigned int molIx = 0; molIx < BATbonds.size(); molIx++){
	// 	for(size_t cnt = 0; cnt < BATbonds[molIx].size(); cnt++){
	
	// 		cout << molIx << " " << cnt << " " << BONDS_to_bonds[molIx][cnt] ;
	// 			ceol;
	// 		scout("bBond "); bonds[BONDS_to_bonds[molIx][cnt]].Print();
	// 			ceol;
	// 		scout("BATbonds "); BATbonds[molIx][cnt].Print();
	// 			ceol;
	// 		ceol;
	
	// 	}
	// }

}

/*!
 * <!--  -->
*/
void Context::reset_BONDS_to_bonds(
	const std::vector<std::vector<BOND>>& BATbonds)
{

	// BONDS_to_bonds.resize(BATbonds.size());
	// for(unsigned int molIx = 0; molIx < BATbonds.size(); molIx++){
	// 	BONDS_to_bonds[molIx].resize(BATbonds[molIx].size(), -111111);
	// }

	// bBonds and BAT bonds equivalence
	scout(" reset bonds BATbonds euqivalence ") << eol;
	bool found = false;

	for(unsigned int molIx = 0; molIx < BATbonds.size(); molIx++){

		for(size_t cnt = 0; cnt < BATbonds[molIx].size(); cnt++){

			found = false;
			for(size_t bbCnt = 0; bbCnt < bonds.size(); bbCnt++){

				if(	(BATbonds[molIx][cnt].first  == bonds[bbCnt].j) &&
					(BATbonds[molIx][cnt].second == bonds[bbCnt].i)){
					found = true;

				}else if((BATbonds[molIx][cnt].first  == bonds[bbCnt].i) &&
						(BATbonds[molIx][cnt].second == bonds[bbCnt].j)){
					found = true;

				}

				if(found){
					BONDS_to_bonds[molIx][cnt] = bbCnt;
					bonds_to_BONDS[bbCnt] = std::make_pair(molIx, cnt);
					break;
				}

			}
		}
	}

	// scout("BONDS_to_bonds");
	// for(unsigned int molIx = 0; molIx < BATbonds.size(); molIx++){
	// 	for(size_t cnt = 0; cnt < BATbonds[molIx].size(); cnt++){
	
	// 		cout << molIx << " " << cnt << " " << BONDS_to_bonds[molIx][cnt] ;
	// 			ceol;
	// 		scout("bBond "); bonds[BONDS_to_bonds[molIx][cnt]].Print();
	// 			ceol;
	// 		scout("BATbonds "); BATbonds[molIx][cnt].Print();
	// 			ceol;
	// 		ceol;
	
	// 	}
	// }

}


/*! <!-- 
 * A) Get child and parent bSpecificAtoms.
 * B) Cook a parentBondCenterPathName = RESNAME + RESID + _ATOMNAME + bond int(next)
 * C) Call Compound::bondAtom(childAtomCompund, BondCenterPathName
 --> */
void Context::buildAcyclicGraph(
	Topology& topology,
	int rootAmberIx,
	int molIx)
{

	const std::vector<BOND>& theseBonds = internCoords.getMoleculeBonds(molIx);

	// Iterate bonds from internal coordinates
	std::size_t bCnt = 0;
	std::vector<BOND>::const_iterator bIt;
	for(bIt = theseBonds.begin(); bIt != theseBonds.end(); bIt++, bCnt++){

		// scout("internCoord bond first second ");
		// bIt->Print();
		// ceol;

		// ===================== GET ATOMS IN THIS BOND ===================
		// Child
		int childAmberIx = bIt->first;
		bSpecificAtom& child = (atoms[childAmberIx]);
		const SimTK::Compound::SingleAtom& childCompoundAtom
			= child.getSingleAtom();

		// Parent
		int parentAmberIx = bIt->second;
		bSpecificAtom& parent = (atoms[parentAmberIx]);
		const SimTK::Compound::SingleAtom& parentCompoundAtom
			= parent.getSingleAtom();

		// Set a molecule identifier
		child.setMoleculeIndex(molIx);
		parent.setMoleculeIndex(molIx);

		// Print
		// scout("bSpecificAtoms parent-child ") << parentAmberIx << " " <<  childAmberIx
		// 	<< eol;

		// ====================== PARENT BOND CENTER ======================
		
		// Convenient vars
		std::stringstream parentBondCenterPathName;
		std::string parentBondCenterPathNameStr = "";
		int parentNextAvailBondCenter = -111111;
		std::string parentNextAvailBondCenterStr = "-111111";
		int parentNofBonds = parent.getNBonds();
		int parentNofFreebonds = parent.getFreebonds();

		// Get next available BondCenter id
		parentNextAvailBondCenter = parentNofBonds - parentNofFreebonds + 1;
		//if((parentNofBonds != 1)){ // this decides if it's bond or bond1
			parentNextAvailBondCenterStr =
				std::to_string(parentNextAvailBondCenter);
		//}

		// Cook the parentBondCenterPathName = RESNAME + RESID + _ATOMNAME + bond int(next)
		parentBondCenterPathName << parent.getName()
			<< "/bond" 
			<< parentNextAvailBondCenterStr;
		parentBondCenterPathNameStr = parentBondCenterPathName.str();

		// ========================================= BOND
		// @@@@@@@@@@@@@@@   BOND    @@@@@@@@@@@@@@@
		// -----------------------------------------
		scout("STUDY Bonding ")
			<< "child " << child.getName() <<" " << child.getInName()
			<<" " << child.getNumber() <<" "
			<< "to parent " << parent.getName() <<" " << parent.getInName() <<" "
			<< parent.getNumber() <<" "
			<< "with bond center name " << parentBondCenterPathNameStr <<" "
			<< eol;

		topology.bondAtom(child.getSingleAtom(),
				(parentBondCenterPathNameStr).c_str(), 0.149, 0);

		// Set the final Biotype
		topology.setAtomBiotype(child.getName(),
								child.getResidueName().c_str(),
								child.getName());

		// Set bSpecificAtom atomIndex to the last atom added to bond
		child.setCompoundAtomIndex(topology.getBondAtomIndex(
			Compound::BondIndex(topology.getNumBonds() - 1), 1));

		// TODO Not so sure about this TODO check if necessary !!!!
		if(parentAmberIx == rootAmberIx){
			parent.setCompoundAtomIndex(topology.getBondAtomIndex(
				Compound::BondIndex(topology.getNumBonds() - 1), 0));
		}

		#ifdef __DRILLING__
			// Print compound atom indices for child and parent
			spacedcout("chiNo", "chi_cAIx", "parNo", "par_cAIx",
				child.getNumber(), child.getCompoundAtomIndex(),
				parent.getNumber(), parent.getCompoundAtomIndex());
			ceol;
		#endif

		// Set bBond Molmodel Compound::BondIndex
		bBond& bond = bonds[ BONDS_to_bonds[molIx][bCnt] ];
		int currentCompoundBondIndex = topology.getNumBonds() - 1;

		bond.setBondIndex(Compound::BondIndex(currentCompoundBondIndex));

		// Not sure where is useful
		bond.setVisited(1);

		// Also set it's molecule index
		bond.setMoleculeIndex(molIx);

		// ====================== RECORD MODIFICATION =====================

		// Decrease freebonds
		parent.decrFreebonds();
		child.decrFreebonds();

	}

	// Handle ions
	if (internCoords.getRoot(molIx).second == -1) {
		atoms[internCoords.getRoot(molIx).first].setMoleculeIndex(molIx);
		atoms[internCoords.getRoot(molIx).first].setCompoundAtomIndex(SimTK::Compound::AtomIndex(0));
	}

	// Add to the list of topologies
	// topologies.push_back(topology);

}

void Context::closeARingWithThisBond(Topology& topology, bBond& bond, int molIx)
{

	// ===================== GET ATOMS IN THIS BOND ===================
	// Child
	//int childAmberIx = bIt->first;
	int childAmberIx = bond.i;
	bSpecificAtom& child = (atoms[childAmberIx]);
	const SimTK::Compound::SingleAtom& childCompoundAtom = child.getSingleAtom();

	// Parent
	//int parentAmberIx = bIt->second;
	int parentAmberIx = bond.j;
	bSpecificAtom& parent = (atoms[parentAmberIx]);
	const SimTK::Compound::SingleAtom& parentCompoundAtom = parent.getSingleAtom();

	// Set a molecule identifier
	child.setMoleculeIndex(molIx);
	parent.setMoleculeIndex(molIx);

	// Bond
	bSpecificAtom &leftNode  = atoms[bond.i];
	bSpecificAtom &rightNode = atoms[bond.j];

	std::stringstream sbuff;
	if (leftNode.getNumber() == topology.baseAtomNumber) {
		sbuff << leftNode.getName() << "/bond" << leftNode.getFreebonds();
	} else {
		sbuff << leftNode.getName() << "/bond" << leftNode.getNBonds() - leftNode.getFreebonds() + 1;
	}

	std::stringstream otsbuff;
	if (rightNode.getNumber() == topology.baseAtomNumber) {
		otsbuff << rightNode.getName() << "/bond" << rightNode.getFreebonds();
	} else {
		otsbuff << rightNode.getName() << "/bond" << rightNode.getNBonds() - rightNode.getFreebonds() + 1;
	}

	// ======================== ACTUAL BONDING ========================
	// scout("Bonding ring ")
	// 	<< child.getName() <<" " << child.getInName()
	// 	<<" " << child.getNumber() <<" " << sbuff.str() <<" "
	// 	<< "to " << parent.getName() <<" " << parent.getInName() <<" "
	// 	<< parent.getNumber() <<" "
	// 	<< "with bond center name " << otsbuff.str() <<" "
	// 	<< eolf;

	topology.addRingClosingBond(
			(sbuff.str()).c_str(),
			(otsbuff.str()).c_str(),
			0.14,
			109*Deg2Rad,
			BondMobility::Rigid);

	// Set bBond Molmodel Compound::BondIndex
	int currentCompoundBondIndex = topology.getNumBonds() - 1;
	bond.setBondIndex(Compound::BondIndex(currentCompoundBondIndex));
	bond.setAsRingClosing();

	// Set the final Biotype
	topology.setAtomBiotype(child.getName(),
							child.getResidueName().c_str(),
							child.getName());
	topology.setAtomBiotype(parent.getName(),
							parent.getResidueName().c_str(),
							parent.getName());									

	// Not sure where is useful
	bond.setVisited(1);

	// Also set it's molecule index
	child.setMoleculeIndex(molIx);
	parent.setMoleculeIndex(molIx);
	bond.setMoleculeIndex(molIx);

	// ====================== RECORD MODIFICATION =====================

	// Decrease freebonds
	parent.decrFreebonds();
	child.decrFreebonds();

}

/*!
 * <!--  -->
*/
void Context::addRingClosingBonds(
	Topology& topology,
	int rootAmberIx,
	int molIx)
{

	for (const auto& b : internCoords.getRingClosingBonds()) {
		for (int bCnt = 0; bCnt < bonds.size(); bCnt++) {
			if (bonds[bCnt].isThisMe(b.first, b.second)) {

				// check if the bond is in this topology with internCoords.getMoleculeBonds( molIx )
				// scout("Ring closing bond: ") << bCnt << " " << b.second << " " << b.first << "; ";
				bool found = false;
				for (const auto& molBond : internCoords.getMoleculeBonds(molIx)) {
					if ((molBond.first == b.first && molBond.second == b.second) ||
						(molBond.second == b.first && molBond.second == b.first)) {
						// std::cout << "Ring closing bond: " << bCnt << " " << b.second << " " << b.first << "; ";
						found = true;
						break;
					}
				}

				if (!found) {
					continue;
				}

				// std::cout << "Ring closing bond: " << bCnt << " " << b.second << " " << b.first << "; ";

				bBond& bond = bonds[bCnt];

				closeARingWithThisBond(topology, bond, molIx);

				break;
			}
		}
	}


}


/*!
 * <!-- Pass newly created topologies to the worlds -->
*/
void Context::passTopologiesToWorlds(void){

	// Iterate worlds
	for(size_t wCnt = 0; wCnt < worlds.size(); wCnt++){

		scout("World ") << wCnt << eol;

		// Get world and its force field
		World& world = worlds[wCnt];

		// Pass current topology to the current world
		world.topologies = &topologies;
	}
}

/*!
 * <!--
 * ========================================================================
 * ======== (0) Read atoms and bonds from all the molecules ===============
 * ========================================================================
  -->
*/

/*! <!-- __fill__ --> */
void Context::build_GmolGraph_MolmodelAcyclicGraph(
){

	// ========================================================================
	// ======== Construct a Compound for every atom ===========================
	// ========================================================================		
	setAtomsCompounds();

	// Biotype will be used to look up molecular
	// force field specific parameters for an atom type
	addBiotypes();

	// ========================================================================
	// ======== (1) Get BAT graphs ============================================
	// ========================================================================

	// Find a root in the unvisited atoms and build BAT graphs
	nofMols = 0;
	while( internCoords.computeRoot( getAtoms() )){ // find a root

		nofMols++;
		//internCoords.PrintRoot();
		
		// Compute the new molecule's BAT coordinates
		internCoords.computeBAT( getAtoms() );
		internCoords.updateVisited(atoms);
		//internCoords.PrintBAT();

	}

	internCoords.computeLevelsAndOffsets( getAtoms() );

	// ========================================================================
	// ======== (2) BAT bonds to bonds ========================================
	// ========================================================================
	//const std::vector<std::vector<BOND>>& BATbonds = internCoords.getBonds();

	load_BONDS_to_bonds( internCoords.getBonds() );

	// ========================================================================
	// ======== (2) Build graphs with bondAtom ================================
	// ========================================================================
	topologies.reserve(nofMols);
	moleculeCount = -1;

	for(unsigned int molIx = 0; molIx < nofMols; molIx++){

		// Add an empty topology
		std::string moleculeName = "MOL" + std::to_string(++moleculeCount);
		Topology topology(moleculeName);

		// --------------------------------------------------------------------
		//  (1) findARoot 
		// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		const int rootAmberIx = internCoords.getRoot( molIx ).first;
		topology.bSpecificAtomRootIndex = rootAmberIx;
		setRootAtom( topology, rootAmberIx );

		// --------------------------------------------------------------------
		// (2) buildAcyclicGraph
		// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		buildAcyclicGraph(topology, rootAmberIx, molIx);

		//addRingClosingBonds_SP_NEW(topology, rootAmberIx, molIx);

		// --------------------------------------------------------------------
		// (4) Add new topology 
		// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		// Add to the list of topologies
		topologies.push_back(topology);

	}
}




/*! <!-- Get BAT graphs --> */
void Context::calc_Gmolmodel_Graph()
{

	// Find a root in the unvisited atoms and build BAT graphs
	nofMols = 0;
	while( internCoords.computeRoot( getAtoms() )){ // find a root

		nofMols++;
		internCoords.PrintRoot();
		
		// Compute the new molecule's BAT coordinates
		internCoords.computeBAT( getAtoms() );
		internCoords.updateVisited(atoms);
		//internCoords.PrintBAT();

	}

	internCoords.computeLevelsAndOffsets( getAtoms() );

}

/*! <!-- Build Molmodel graphs with bondAtom --> */
void Context::build_Molmodel_AcyclicGraphs(void)
{

	topologies.reserve(nofMols);
	moleculeCount = -1;

	std::cout << "STUDY Context::build_Molmodel_AcyclicGraphs nofMols " <<" "<< nofMols <<std::endl<<std::flush;

	for(unsigned int molIx = 0; molIx < nofMols; molIx++){

		// Add an empty topology
		std::string moleculeName = "MOL" + std::to_string(++moleculeCount);
		Topology topology(moleculeName);

		// --------------------------------------------------------------------
		//  (1) findARoot 
		// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		const int rootAmberIx = internCoords.getRoot( molIx ).first;
		topology.bSpecificAtomRootIndex = rootAmberIx;

		std::cout << "STUDY Context::build_Molmodel_AcyclicGraphs rootAmberIx " <<" "<< rootAmberIx <<std::endl<<std::flush;

		setRootAtom( topology, rootAmberIx );

		// --------------------------------------------------------------------
		// (2) buildAcyclicGraph
		// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		buildAcyclicGraph(topology, rootAmberIx, molIx);

		// --------------------------------------------------------------------
		// (4) Add new topology 
		// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		// Add to the list of topologies
		topologies.push_back(topology);

	}
}


/*! <!-- __fill__ --> */
void Context::new_build_GmolGraph_MolmodelAcyclicGraph()
{

	calc_Gmolmodel_Graph();

	load_BONDS_to_bonds( internCoords.getBonds() );

	build_Molmodel_AcyclicGraphs();

}


























/*!
 * <!-- Subarray view of atoms -->
*/
void Context::addRingClosingBonds_All(void)
{
	// Variant 1
	/* //for(int molIx = nofMols - 1; molIx > -1; molIx--){
	for(int molIx = 0; molIx < nofMols; molIx++){

		Topology& currTopology = topologies[molIx];

		addRingClosingBonds_SP_NEW(
			currTopology,
			internCoords.getRoot( molIx ).first,
			molIx);

	} */

	// Variant 2 [BIGPROT]
	/* for (int bbCnt = 0; bbCnt < bonds.size(); bbCnt++) {

		bBond& bond = bonds[bbCnt];

		if(bond.isVisited() == 0){

			int molIx = -1, BOCnt = -1;

			molIx = bonds_to_BONDS[bbCnt].first;
			BOCnt = bonds_to_BONDS[bbCnt].second;

			Topology& topology = topologies[molIx];
			scout("Ring closing bond ") ; PrintBond(bond);
			closeARingWithThisBond(topology, bond, molIx);

		} // if is visited
	} // every bBond */


	// Variant 3
	for (int bbCnt = 0; bbCnt < bonds.size(); bbCnt++) {

		bBond& bond = bonds[bbCnt];

		if(bond.isVisited() == 0){

			bSpecificAtom& left = atoms[bond.i];
			bSpecificAtom& right = atoms[bond.j];

			int molIx = left.getMoleculeIndex();

			bond.setMoleculeIndex(molIx);

			Topology& topology = topologies[molIx];

			//scout("Ring closing bond ") ; PrintBond(bond);

			closeARingWithThisBond(topology, bond, molIx);

		} // if is visited
	} // every bBond

}

/*!
 * <!-- Subarray view of atoms -->
*/
void Context::generateSubAtomLists(void){
	
	// Resize subAtomList
	subAtomLists.resize( getNofMolecules() );

	// Detect begin and start indeces
	std::vector<std::pair<size_t, size_t>> molRanges;
	molRanges.reserve(getNofMolecules());

	std::vector<size_t> molStarts;
	int nextMolStart = -1;

	// Fill a vector of molecule start indeces
	std::set<size_t> uniqueMolIx; // unique molecule indeces

	for(size_t atomCnt = 0; atomCnt < atoms.size(); atomCnt++){

		bool inserted =
			uniqueMolIx.insert( atoms[atomCnt].getMoleculeIndex() ).second;

		if ( inserted ) {
			molStarts.push_back(atomCnt);
		}
	}

	// Also add a last element
	molStarts.push_back(atoms.size());

	// Build ranges of molecules
	for (size_t i = 0; i < molStarts.size() - 1; i++) {
		molRanges.push_back(std::make_pair(molStarts[i], molStarts[i + 1]));
	}

	// Sort by molecule index
	std::sort(molRanges.begin(), molRanges.end(),[this]
		(const std::pair<std::size_t, std::size_t>& lhs,
		 const std::pair<std::size_t, std::size_t>& rhs) {
			return atoms[lhs.first].getMoleculeIndex() < atoms[rhs.first].getMoleculeIndex();
	});

	// Print
	// scout("subAtoms ranges\n");
	// for(size_t cnt = 0; cnt < molRanges.size(); cnt++){
	// 	cout << molRanges[cnt].first << " " << molRanges[cnt].second << eol;
	// }

	std::cout << molRanges.size() << " " << getNofMolecules() << std::endl;

	assert("atoms molecule indexes different than nofMols" &&
		(molRanges.size() == getNofMolecules()));

	// Set atom sublists for every Compound
	int sIx = -1;
	if( !(molRanges.empty()) ){
		for( size_t cnt = 0; cnt < molRanges.size(); cnt ++ ){

			sIx++;
			subAtomLists[sIx].set_view(
				atoms.begin() + molRanges[cnt].first,
				atoms.begin() + molRanges[cnt].second);

		}
	}

	// Check
	// for(const auto& view:subAtomLists){
	// 	scout("New molecule") << eol;
	// 	for(auto it = view.begin(); it != view.end(); it++){
	// 		it->Print(0);
	// 	}
	// }


}

/*!
 * <!-- Subarray view of bonds -->
*/
void Context::generateSubBondLists(void){
	
	// Resize subAtomList
	subBondLists.resize( getNofMolecules() );

	// Detect begin and start indeces
	std::vector<std::pair<size_t, size_t>> molRanges;
	molRanges.reserve(getNofMolecules());

	std::vector<size_t> molStarts;
	int nextMolStart = -1;

	// Fill a vector of molecule start indeces
	std::set<size_t> uniqueMolIx; // unique molecule indeces

	for(size_t bondCnt = 0; bondCnt < bonds.size(); bondCnt++){

		bool inserted =
			uniqueMolIx.insert( bonds[bondCnt].getMoleculeIndex() ).second;

		if ( inserted ) {
			molStarts.push_back(bondCnt);
		}
	}

	// Also add a last element
	molStarts.push_back(bonds.size());

	// Build ranges of molecules
	for (size_t i = 0; i < molStarts.size() - 1; i++) {
		molRanges.push_back(std::make_pair(molStarts[i], molStarts[i + 1]));
	}

	// Sort by molecule index
	std::sort(molRanges.begin(), molRanges.end(),[this]
		(const std::pair<std::size_t, std::size_t>& lhs,
		 const std::pair<std::size_t, std::size_t>& rhs) {
			return bonds[lhs.first].getMoleculeIndex() < bonds[rhs.first].getMoleculeIndex();
	});

	// scout("subBonds ranges\n");
	// for(size_t cnt = 0; cnt < molRanges.size(); cnt++){
	// 	cout << molRanges[cnt].first << " " << molRanges[cnt].second << eol;
	// }

	// Set atom sublists for every Compound
	int sIx = -1;
	if( !(molRanges.empty()) ){
		for( size_t cnt = 0; cnt < molRanges.size(); cnt ++ ){

			sIx++;
			subBondLists[sIx].set_view(
				bonds.begin() + molRanges[cnt].first,
				bonds.begin() + molRanges[cnt].second);
		}
	}

	// // Check
	// for(const auto& view:subBondLists){
	// 	scout("New molecule") << eol;
	// 	for(auto it = view.begin(); it != view.end(); it++){
	// 		it->Print();
	// 		ceol;
	// 	}
	// }

}

/*!
 * <!--
 * ========================================================================
 * ======== (3) Generate subarray views for atoms and bonds ===============
 * ========================================================================
 * -->
*/
void Context::generateTopologiesSubarrays(void){

	// Atoms - should never be sorted though
	// std::sort(atoms.begin(), atoms.end(), [](
	// 	const bSpecificAtom& lhs, const bSpecificAtom& rhs){
	// 		return lhs.getMoleculeIndex() < rhs.getMoleculeIndex();
	// 	}
	// );
	
	// Generate array_views for atoms in every topology
	generateSubAtomLists();

	// scout("Loaded bonds") << eol;
	// PrintBonds();
	
	// SORT the bonds after molecule:
	// ATENTION: this changes all maps containg bonds indeces
	std::sort(bonds.begin(), bonds.end(), [](
		const bBond& lhs, const bBond& rhs){
			return lhs.getMoleculeIndex() < rhs.getMoleculeIndex();
		}
	);

	// Keep a map of BONDS to bonds after sorting
	reset_BONDS_to_bonds( internCoords.getBonds() );

	// Allocate space for bond mappings
	bondMapping.reserve(bonds.size());

	// Create bond mappings
	for (const auto& bond : bonds) {
		
		// Get Robosample bond index
		const int rIx = bond.getIndex();

		// Get Molmodel bond index
		const SimTK::Compound::BondIndex molmodelBIx = bond.getBondIndex();

		// Map Robosample bond index to Molmodel bond index
		bondMapping[rIx] = molmodelBIx;
	}

	// Pass bond mappings to all topologies (all molecules)
	for (auto& topology : topologies) {
		topology.setBondMappings(bondMapping);
	}

	// Generate array_views for bonds and bonds in every topology
	generateSubBondLists();

	for(unsigned int molIx = 0; molIx < nofMols; molIx++){

		Topology& topology = topologies[molIx];

		topology.setSubAtomList(
			subAtomLists[molIx].begin(),
			subAtomLists[molIx].end(),
			elementCache);

		topology.generateAIx2TopXMaps();

		topology.setSubBondList(
			subBondLists[molIx].begin(),
			subBondLists[molIx].end());

	}

	// Set offset for all subBondLists
	for(unsigned int molIx = 0; molIx < nofMols; molIx++){

		Topology& topology = topologies[molIx];

		auto minNumberAtomIt = std::min_element(topology.subAtomList.begin(), topology.subAtomList.end(),
			[](const bSpecificAtom& lhs, const bSpecificAtom& rhs) {
				return lhs.getNumber() < rhs.getNumber();
			});
		
		topology.subBondList.set_offset( minNumberAtomIt->getNumber() );
	}

}


/*!
 * <!-- Assign Compound coordinates by matching bAtomList coordinates -->
*/
void Context::matchDefaultConfigurationFromAtomsCoords(Topology& topology, int molIx)
{
	// Convenient vars
	std::map<Compound::AtomIndex, SimTK::Vec3> atomTargets;
	array_view<std::vector<bSpecificAtom>::iterator>& topoSubAtomList =
		topology.subAtomList;

	// For every atom in this topology deposit coords in atomTargets
	for(int ix = 0; ix < topology.getNumAtoms(); ++ix){
		
		SimTK::Vec3 atomCoords(
			topoSubAtomList[ix].getX(),
			topoSubAtomList[ix].getY(),
			topoSubAtomList[ix].getZ());

		atomTargets.insert(std::pair<Compound::AtomIndex, SimTK::Vec3> (
			topoSubAtomList[ix].getCompoundAtomIndex(),
			atomCoords));

	}

	// Match coordinates

	// Optimize superpose in matchDefaultTopLevelTransform
	SimTK::Vec3 topLevelShift = SimTK::Vec3(
		atoms[internCoords.getRoot(molIx).first].getX(),
		atoms[internCoords.getRoot(molIx).first].getY(),
		atoms[internCoords.getRoot(molIx).first].getZ()
	);
	topology.setTopLevelTransform(Transform(Rotation(), topLevelShift));
	//PrintTransform(topology.getTopLevelTransform(), 3, topology.getName(), "Top");

	//std::cout << "Match start." << "\n" << std::flush;
	bool flipAllChirality = true;
	topology.matchDefaultBondLengths(atomTargets);
	//std::cout << "matchDefaultBondLengths done. " << "\n" << std::flush;
	topology.matchDefaultAtomChirality(atomTargets, 0.01, flipAllChirality);
	//std::cout << "matchDefaultAtomChirality done. " << "\n" << std::flush;
	topology.matchDefaultBondAngles(atomTargets);
	//std::cout << "matchDefaultBondAngles done. " << "\n" << std::flush;
	topology.matchDefaultDirections(atomTargets);
	//std::cout << "matchDefaultDirections done. " << "\n" << std::flush;
	topology.matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
	//std::cout << "matchDefaultDefaultDihedralAngles done. " << "\n" << std::flush;
	topology.matchDefaultTopLevelTransform(atomTargets);
	//std::cout << "matchDefaultDefaultTopLevelTransform done. " << "\n" << std::flush;
	//PrintTransform(topology.getTopLevelTransform(), 3, topology.getName(), "Top");

}

/*!
 * <!-- Match Compounds configurations to atoms Cartesian coords -->
*/
void Context::matchDefaultConfigurations(void){

	// --------------------------------------------------------------------
	// matchDefaultConfiguration
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	// Assign Compound coordinates by matching bAtomList coordinates
	for(unsigned int molIx = 0; molIx < nofMols; molIx++){

		Topology& topology = topologies[molIx];

		matchDefaultConfigurationFromAtomsCoords(topology, molIx);
		// PrintAtoms();
	}

	// Helper function for calc MBAT determinant // TO BE DELETED
	for(unsigned int molIx = 0; molIx < nofMols; molIx++){

		Topology& topology = topologies[molIx];

		topology.loadTriples_SP_NEW();

	}

	// Map of Compound atom indexes to Robosample atom indexes
	for(unsigned int molIx = 0; molIx < nofMols; molIx++){

		Topology& topology = topologies[molIx];

		topology.loadCompoundAtomIx2GmolAtomIx();

	}
	// ========================================================================
	// ======== (3) Rest from AddMolecules ====================================
	// ========================================================================

	// Context topologies to all the worlds
	passTopologiesToWorlds();

}



/*!
 * <!-- Long print of all atoms properties -->
*/
void Context::PrintAtoms(void){
	for(const auto& atom : atoms) {
		atom.Print(0);
	}
}

// TODO: move to world
/** It calls DuMMs defineAtomClass. These Molmodel functions contain
information regarding the force field parameters. **/
void Context::updDummAtomClasses(
	std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
	, int worldIx			
)
{
	// Accumulate DuMM parameters in these vectors
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allBondsACIxs;
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allAnglesACIxs;
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allDihedralsACIxs;
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allImpropersACIxs;


		SimTK::DuMM::AtomClassIndex aCIx;
		std::string atomClassName;
		// Iterate through AtomClasses map and put AtomClasses in Dumm
		std::map<AtomClassParams, AtomClassId>::const_iterator it;
		for (	it = aClassParams2aClassId.begin();
			it != aClassParams2aClassId.end(); ++it){

			const AtomClassParams& atomParams = it->first;
			const AtomClassId& atomClassId = it->second;

			aCIx = atomClassId.dummAtomClassIndex;
			atomClassName = atomClassId.name;

			/* std::cout << "Context::transferAtomClasses "
				<< aCIx << " " << atomClassName ;
			atomParams.dump(); */

			// Define an AtomClass
			worlds[worldIx].forceField->defineAtomClass(aCIx, atomClassName.c_str(),
				atomParams.atomicNumber,
				atomParams.valence,
				atomParams.vdwRadius,
				atomParams.LJWellDepth
			);
		}
}

// ============================================================================
// NEW WAY TO ADD PARAMS
// ============================================================================
/*!
 * <!-- It calls DuMMs defineAtomClass, defineChargedAtomTye and
 * setBiotypeChargedAtomType for every atom. These Molmodel functions contain
 * information regarding the force field parameters.-->
*/
void Context::generateDummAtomClasses(AmberReader& amberReader)
{
	// Iterate worlds
	for(size_t wCnt = 0; wCnt < worlds.size(); wCnt++){
		
		// Convenient vars
		World& world = worlds[wCnt];
		SimTK::DuMMForceFieldSubsystem& dumm = *(world.forceField);
		std::vector<bool> founditInDuMM(atoms.size(), false);

		// scout("World ") << wCnt << eol;

		// ========================================================================
		// ======== (1) DuMM atom classes =========================================
		// ========================================================================

		// Iterate atoms
		for(size_t aCnt = 0; aCnt < atoms.size(); aCnt++){

			SimTK::DuMM::AtomClassIndex dummAtomClassIndex;
			std::string atomClassName;

			// Define an AtomClassparam for this atom
			AtomClassParams atomClassParams(
				atoms[aCnt].getAtomicNumber(),
				atoms[aCnt].getNBonds(),
				atoms[aCnt].getVdwRadius() / 10.0, // nm
				atoms[aCnt].getLJWellDepth() * 4.184 // kcal to kJ
			);

			std::string str(atoms[aCnt].getFftype());
			const SimTK::String simtkstr(str);

			founditInDuMM[aCnt] = dumm.hasAtomClass(simtkstr);
			// std::cout << "Context::generateDummAtomClasses_SP_NEW "
			// 	<< "world " << wCnt << " atom " << aCnt << " "
			// 	<< atoms[aCnt].getFftype() << " "
			// 	<< founditInDuMM[aCnt] << std::endl << std::flush;

			// Define the AtomClass
			if (!founditInDuMM[aCnt]){ // We don't have this AtomClass
				
				// Get an AtomClass index from DuMM
				dummAtomClassIndex = dumm.getNextUnusedAtomClassIndex();

				// Define an AtomClass name
				atomClassName = atoms[aCnt].getFftype();

				// Define an AtomClass
				dumm.defineAtomClass(dummAtomClassIndex,
					atomClassName.c_str(),
					atomClassParams.atomicNumber,
					atomClassParams.valence,
					atomClassParams.vdwRadius,
					atomClassParams.LJWellDepth
				);

				// scout("Added atom aCnt atomClassIndex ") 
				// 	<< aCnt <<" " << dummAtomClassIndex <<" "
				// 	<< "|" << atomClassName <<"| "
				// 	<< 	atomClassParams.atomicNumber <<" "
				// 	<< atomClassParams.valence <<" "
				// 	<< atomClassParams.vdwRadius <<" "
				// 	<< atomClassParams.LJWellDepth <<" "
				// 	<< eol;				

			}else{ // We have this AtomClass

				// Get AtomClass index from DuMM
				dummAtomClassIndex = dumm.getAtomClassIndex(
					atoms[aCnt].getFftype());
			}			

			// Insert AtomClass index in Gmolmodel atom list too
			atoms[aCnt].setDummAtomClassIndex(dummAtomClassIndex);

		} // every atom


		// Define ChargedAtomTypeIndeces
		SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
		std::string chargedAtomTypeName;

		// ========================================================================
		// ======== (2) DuMM charged atom types ===================================
		// ========================================================================

		// Iterate atoms
		for(size_t aCnt = 0; aCnt < atoms.size(); aCnt++){

			// Get a ChargedAtomType index
			chargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex();
			atoms[aCnt].setChargedAtomTypeIndex(chargedAtomTypeIndex);

			// Define a chargedAtomType name
			chargedAtomTypeName =  atoms[aCnt].getResidueName(); // taken from Amber
			chargedAtomTypeName += atoms[aCnt].getBiotype();

			// Define a ChargedAtomType (AtomClass with a charge)
			dumm.defineChargedAtomType(
				chargedAtomTypeIndex,
				chargedAtomTypeName.c_str(),
				atoms[aCnt].getDummAtomClassIndex(),
				atoms[aCnt].charge
			);

			// scout("SP_NEW_LAB Defined chargedAtomType ") << chargedAtomTypeName  <<" "
			// 	<< chargedAtomTypeIndex <<" "
			// 	<< "|" << chargedAtomTypeName <<"| "
			// 	<< atoms[aCnt].getDummAtomClassIndex() <<" "
			// 	<< atoms[aCnt].charge
			// 	<< eol;

			// Associate a ChargedAtomTypeIndex with a Biotype index
			dumm.setBiotypeChargedAtomType(
				atoms[aCnt].getChargedAtomTypeIndex(),
				atoms[aCnt].getBiotypeIndex()
			);

			// scout("SP_NEW_LAB setBiotypeChargedAtomType biotypeIx chargedATIx ")
			// 	<< atoms[aCnt].getBiotypeIndex()  <<" "
			// 	<< atoms[aCnt].getChargedAtomTypeIndex() <<" "
			// 	<< eol;

		} // every atom

	} // every world
}

/*!
 * <!-- Calls DuMM defineBondStretch to define bonds parameters. -->
*/
void Context::bAddDummBondParams(AmberReader& amberReader)
{

	// Iterate worlds
	for(size_t wCnt = 0; wCnt < worlds.size(); wCnt++){

		// Get world and its force field
		World& world = worlds[wCnt];
		SimTK::DuMMForceFieldSubsystem& dumm = *(world.forceField);

		// Keep track of inserted AtomClass pairs		
		std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allBondsACIxs;
	
		// scout("Dumm bonds") << eol;

		// Iterate through bonds and define their parameters
		for(size_t bCnt = 0; bCnt < bonds.size(); bCnt++){

			// Get atoms
			int atomNumber1 = bonds[bCnt].i;
			int atomNumber2 = bonds[bCnt].j;
			bSpecificAtom& atom1 = atoms[atomNumber1];
			bSpecificAtom& atom2 = atoms[atomNumber2];

			// Generate a pair of atom classes for this bond
			std::vector<SimTK::DuMM::AtomClassIndex> thisBondACIxs;

			thisBondACIxs.push_back( SimTK::DuMM::AtomClassIndex(
				atom1.getDummAtomClassIndex()) );
			thisBondACIxs.push_back( SimTK::DuMM::AtomClassIndex(
				atom2.getDummAtomClassIndex()) );

			// Check if we already have this bond
			bool foundit = false;
			for(auto& row:allBondsACIxs){
				if ( IsTheSameBond (thisBondACIxs, row) ){
					foundit = true;
					break;
				}
			}

			// Add bond to Dumm
			if (  !foundit ){ // bond was not found

				dumm.defineBondStretch_KA(
					atom1.getDummAtomClassIndex(),
					atom2.getDummAtomClassIndex(),
					bonds[bCnt].getForceK(),  //k1
					bonds[bCnt].getForceEquil()   //equil1
				);				

				// Print
				Topology& topology1 = topologies[ atom1.getMoleculeIndex() ];
				Topology& topology2 = topologies[ atom2.getMoleculeIndex() ];

				const SimTK::Compound::AtomIndex cAIx1 = atom1.getCompoundAtomIndex(); 
				const SimTK::Compound::AtomIndex cAIx2 = atom2.getCompoundAtomIndex(); 
				const SimTK::DuMM::AtomIndex dAIx1 = topology1.getDuMMAtomIndex(cAIx1);
				const SimTK::DuMM::AtomIndex dAIx2 = topology2.getDuMMAtomIndex(cAIx2);

				// LAB BEGIN
				//SimTK::DuMM::IncludedAtomIndex inclAIx = a.getIncludedAtomIndex();
				// dumm.getIncludedAtomIndexOfDummAtom(dAIx);
				// SimTK::DuMM::AtomIndex	 dAIx = getAtomIndexOfNonbondAtom (SimTK::DuMM::NonbondAtomIndex nbDAIx);
				// SimTK::DuMM::AtomIndex dAIx = dumm.getAtomIndexOfIncludedAtom(inclDAIx);
				// SimTK::DuMM::IncludedAtomIndex inclDAIx = dumm.getIncludedAtomIndexOfNonbondAtom(nbDAIx);
				// LAB END

				// scout("bond ") << bonds[bCnt].getMoleculeIndex() <<" "
				// 	<< atomNumber1 <<" " << atomNumber2 <<" "
				// 	<< atom1.getInName() <<" " << atom2.getInName() <<" "
				// 	<< cAIx1 <<" " << cAIx2 <<" "
				// 	<< bonds[bCnt].getForceK() <<" "
				// 	<< bonds[bCnt].getForceEquil() <<" "
				// 	<< eol;

				// Put the entry in our map too
				allBondsACIxs.push_back(thisBondACIxs);

			}
		} // evvery bond
	} // every world

}

/*!
 * <!-- Calls DuMM defineBondBend to define angle parameters. -->
*/
void Context::bAddDummAngleParams(AmberReader& amberReader)
{

	// Iterate worlds
	for(size_t wCnt = 0; wCnt < worlds.size(); wCnt++){

		// scout("World ") << wCnt << eol;

		// Get world and its force field
		World& world = worlds[wCnt];
		SimTK::DuMMForceFieldSubsystem& dumm = *(world.forceField);

		// Keep track of inserted AtomClass pairs
		std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allAnglesACIxs;

		// Iterate angles and define their parameters
		for(int angCnt = 0; angCnt < amberReader.getNumberAngles(); angCnt++){

			// Gmolmodel Atom indeces
			int a1 = amberReader.getAnglesAtomsIndex1(angCnt);
			int a2 = amberReader.getAnglesAtomsIndex2(angCnt);
			int a3 = amberReader.getAnglesAtomsIndex3(angCnt);

			// Generate a triple of atom class indexes for this angle
			std::vector<SimTK::DuMM::AtomClassIndex> thisAngleACIxs;
			thisAngleACIxs.push_back( SimTK::DuMM::AtomClassIndex(
				atoms[a1].getDummAtomClassIndex()) );
			thisAngleACIxs.push_back( SimTK::DuMM::AtomClassIndex(
				atoms[a2].getDummAtomClassIndex()) );
			thisAngleACIxs.push_back( SimTK::DuMM::AtomClassIndex(
				atoms[a3].getDummAtomClassIndex()) );

			// Check if we already have this angle
			bool foundit = false;
			for(auto& row:allAnglesACIxs){
				if ( IsTheSameAngle (thisAngleACIxs, row) ){
					foundit = true;
					break;
				}
			}

			// Add angle to Dumm
			if (  !foundit ){ // angle was not found
			
				dumm.defineBondBend_KA(
					atoms[a1].getDummAtomClassIndex(),
					atoms[a2].getDummAtomClassIndex(),
					atoms[a3].getDummAtomClassIndex(),
					amberReader.getAnglesForceK(angCnt),
					static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE 
						* amberReader.getAnglesEqval(angCnt))) // TODO 32 vs 64 bit
				);

				// Put the entry in our map too
				allAnglesACIxs.push_back(thisAngleACIxs);

				// scout("angle ")
				// 	<< a1 <<" " << a2 <<" " << a3 <<" "
				// 	<< atoms[a1].getInName() <<" "
				// 	<< atoms[a2].getInName() <<" "
				// 	<< atoms[a3].getInName() <<" "
				// 	<< amberReader.getAnglesForceK(angCnt) <<" "
				// 	<< static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE 
				// 		* amberReader.getAnglesEqval(angCnt)))
				// 	<< eol;
			}
		}		

	} // every world
}

/*!
 * <!-- Check if a1 and a2 are bonded. -->
*/
bool Context::checkBond(int a1, int a2)
{
	for(int i = 0; i < nbonds; i++){
		if( (bonds[i]).isThisMe(a1, a2) )
		{
			return true;
		}
	}
	return false;
}

/*!
 * <!-- Calls DuMM defineBondTorsion. -->
*/
void Context::bAddDummTorsionParams(AmberReader& amberReader)
{

	// Iterate worlds
	for(size_t wCnt = 0; wCnt < worlds.size(); wCnt++){

		// scout("World ") << wCnt << eol;

		// Get world and its force field
		World& world = worlds[wCnt];
		SimTK::DuMMForceFieldSubsystem& dumm = *(world.forceField);

		// Keep track of inserted AtomClass pairs
		std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allDihedralsACIxs;
		std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allImpropersACIxs;

		// Amber reader dihedrals vector
		std::vector<std::pair<int, int>> pairStartAndLens =
			amberReader.getPairStartAndLen();

		for(unsigned int index = 0; index < pairStartAndLens.size(); index++){

			// Get start and len of this dihedral
			int torsCnt        = pairStartAndLens[index].first;
			int numberOf = pairStartAndLens[index].second;

			// Get Amber indeces
			int amber_aIx_1 = amberReader.getDihedralsAtomsIndex1(torsCnt);
			int amber_aIx_2 = amberReader.getDihedralsAtomsIndex2(torsCnt);
			int amber_aIx_3 = amberReader.getDihedralsAtomsIndex3(torsCnt);
			int amber_aIx_4 = amberReader.getDihedralsAtomsIndex4(torsCnt);


			// Check if a quad of atom indices is a normal dihedral
			// or an improper dihedral, by checking if consecutive
			// atoms are bonded 

			bool dihedral=false;
			bool improper=true;

			if (checkBond(amber_aIx_1, amber_aIx_2) &&
				checkBond(amber_aIx_2, amber_aIx_3) &&
				checkBond(amber_aIx_3, amber_aIx_4))
			{
				dihedral = true;
				improper = false;
			}

			// Get AtomClass indeces 
			SimTK::DuMM::AtomClassIndex aCIx1 =
				atoms[amber_aIx_1].getDummAtomClassIndex();
			SimTK::DuMM::AtomClassIndex aCIx2 =
				atoms[amber_aIx_2].getDummAtomClassIndex();
			SimTK::DuMM::AtomClassIndex aCIx3 =
				atoms[amber_aIx_3].getDummAtomClassIndex();
			SimTK::DuMM::AtomClassIndex aCIx4 =
				atoms[amber_aIx_4].getDummAtomClassIndex();

				
			// Generate a quad of atom class indexes for this dihedral, 
			// regardless of whether it's a torsion or an improper
			std::vector<SimTK::DuMM::AtomClassIndex> thisDihedralACIxs;
			thisDihedralACIxs.push_back(aCIx1);
			thisDihedralACIxs.push_back(aCIx2);
			thisDihedralACIxs.push_back(aCIx3);
			thisDihedralACIxs.push_back(aCIx4);

			if (dihedral){

				// If it is a normal dihedral, check if we have it in our
				// dihedral list
				bool foundit = false;

				for(auto& row:allDihedralsACIxs)
				{
					if ( IsTheSameTorsion (thisDihedralACIxs, row))
					{
						foundit = true;	break;
					}
				}

				if (!foundit){ // dihedral was not found

					// scout("dihedral ")
					// 		<< amber_aIx_1 <<" " << amber_aIx_2 <<" "
					// 		<< amber_aIx_3 <<" " << amber_aIx_4 <<" "
					// 		<< atoms[amber_aIx_1].getInName() <<" "
					// 		<< atoms[amber_aIx_2].getInName() <<" "
					// 		<< atoms[amber_aIx_3].getInName() <<" "
					// 		<< atoms[amber_aIx_4].getInName() <<" ";

					// Define the dihedrals
					if(numberOf == 1){
						dumm.defineBondTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt)))
						);

						// scout("") << amberReader.getDihedralsForceK(torsCnt) <<" "
						// 	<< static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE *
						// 		amberReader.getDihedralsPhase(torsCnt))) <<" ";

					}
					else if(numberOf == 2){
						dumm.defineBondTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt))),
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt + 1)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt + 1),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt+1)))
						);

						// scout("") << amberReader.getDihedralsForceK(torsCnt + 1) <<" "
						// 	<< static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE *
						// 		amberReader.getDihedralsPhase(torsCnt + 1))) <<" ";

					}
					else if(numberOf == 3){
						dumm.defineBondTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt))),
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt + 1)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt + 1),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt+1))),
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt + 2)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt + 2),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt+2)))
						);

						// scout("") << amberReader.getDihedralsForceK(torsCnt + 2) <<" "
						// 	<< static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE *
						// 		amberReader.getDihedralsPhase(torsCnt + 2))) <<" ";

					}else if (numberOf == 4){
						dumm.defineBondTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt))),
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt + 1)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt + 1),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt+1))),
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt + 2)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt + 2),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt+2))),
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt+ 3)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt+3),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt+3)))
						);

						// scout("") << amberReader.getDihedralsForceK(torsCnt + 3) <<" "
						// 	<< static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE *
						// 		amberReader.getDihedralsPhase(torsCnt + 3))) <<" ";
					}
					
					//ceol;

					// Add the dihedral to the list of impropers.
					allDihedralsACIxs.push_back(thisDihedralACIxs);

				} // END of not foundit
		
			}
			
			if (improper){

				// If it is an improper dihedral, we check if it exitsts, without
				// checking for the reverse (order matters for impropers)

				bool foundit = false;
				for(auto& row:allImpropersACIxs){
					if ((thisDihedralACIxs == row))
						{
						foundit = true;	
						break;
						}
					}
				
				if (!foundit){ // improper was not found

					// scout("improper ")
					// 		<< amber_aIx_1 <<" " << amber_aIx_2 <<" "
					// 		<< amber_aIx_3 <<" " << amber_aIx_4 <<" "
					// 		<< atoms[amber_aIx_1].getInName() <<" "
					// 		<< atoms[amber_aIx_2].getInName() <<" "
					// 		<< atoms[amber_aIx_3].getInName() <<" "
					// 		<< atoms[amber_aIx_4].getInName() <<" ";

					// Define the dihedrals
					if(numberOf == 1){
						dumm.defineAmberImproperTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt)))
						);

						// scout("") << amberReader.getDihedralsForceK(torsCnt) <<" "
						// 	<< static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE *
						// 		amberReader.getDihedralsPhase(torsCnt))) <<" ";

					}
					else if(numberOf == 2){
						dumm.defineAmberImproperTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt))),
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt + 1)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt + 1),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt+1)))
						);

						// scout("") << amberReader.getDihedralsForceK(torsCnt + 1) <<" "
						// 	<< static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE *
						// 		amberReader.getDihedralsPhase(torsCnt + 1))) <<" ";

					}
					else if(numberOf == 3){
						dumm.defineAmberImproperTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt))),
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt + 1)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt + 1),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt+1))),
							static_cast<int>(amberReader.getDihedralsPeriod(torsCnt + 2)), // TODO wants int, returns double
							amberReader.getDihedralsForceK(torsCnt + 2),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader.getDihedralsPhase(torsCnt+2)))
						);

						// scout("") << amberReader.getDihedralsForceK(torsCnt + 2) <<" "
						// 	<< static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE *
						// 		amberReader.getDihedralsPhase(torsCnt + 2))) <<" ";

					}

					//ceol;

					// Add the improper to the list of impropers.
					allImpropersACIxs.push_back(thisDihedralACIxs);
				}
			} // improper

		} // torsions


	} // every world

}

/*!
 * <!-- It calls DuMMs defineAtomClass, defineChargedAtomTye and
 * setBiotypeChargedAtomType for every atom. These Molmodel functions contain
 * information regarding the force field parameters.-->
*/
void Context::addDummParams(AmberReader& amberReader)
{

	// DANGER: BONDS WERE SORTED =>
	// amberReader and bonds have no longer the same indeces

	// Add atom types and charged atmo types
	generateDummAtomClasses(amberReader);

	// Add bonds parameters
	bAddDummBondParams(amberReader);

	// Add angle parameters
	bAddDummAngleParams(amberReader);

	// Add torsion parameters
	bAddDummTorsionParams(amberReader);

}


// ============================================================================
// MODEL
// ============================================================================
/*!
 * <!-- Set all flexibilities for all the worlds to Rigid. -->
*/
void Context::initializeFlexibility(void)
{

	// Iterate worlds
	for(size_t worldIx = 0; worldIx < worlds.size(); worldIx++){

		// Set all Compounds' bonds mobilities to Rigid as a default
		for(size_t topCnt = 0; topCnt < topologies.size(); topCnt++){

			Topology& topology = topologies[topCnt];

			for(size_t topBCnt = 0; topBCnt < topology.getNumBonds(); topBCnt++){
				topology.setBondMobility(BondMobility::Rigid,
					SimTK::Compound::BondIndex(topBCnt));
			}
		}

		// Set all Gmolmodel bonds Rigid also
		for(size_t bCnt = 0; bCnt < bonds.size(); bCnt++){
			bonds[bCnt].addBondMobility(SimTK::BondMobility::Rigid);
		}

	}

}

/*!
 * <!-- Set flexibilities. -->
*/
void Context::setFlexibility(
	std::string argRegimen,
	std::string flexFN,
	int whichWorld)
{

	// Get flexible bonds from file. Numbering starts at 0 in prmtop
	std::string line;
	std::ifstream flexF(flexFN);

	int molIx = -1;

	while(flexF.good()){

		// Get a line
		std::getline(flexF, line);
		if(!line.empty()){

			// Comment
			if(line.at(0) == '#'){continue;}

			// Get words
			std::istringstream iss(line);
			std::string word;
			std::vector<std::string> lineWords;

			while(iss >> word){
				lineWords.push_back(std::move(word));
			}

			// Check the line
			if(lineWords.size() >= 3 ){

				int index_1 = std::stoi(lineWords[0]);
				int index_2 = std::stoi(lineWords[1]);
				std::string mobility =  lineWords[2];

				// Add root mobilities
				// if(index_1 == -1){
				// 	molIx++;
				// 	if(molIx < rootMobilities[whichWorld].size()){
				// 		rootMobilities[whichWorld][molIx] = mobility;
						
				// 		// scout("Added rootMobility world ")
				// 		// 	<< whichWorld << " mol " << molIx <<" "
				// 		// 	<< rootMobilities[whichWorld][molIx] <<" "
				// 		// 	<< eol;

				// 	}else{
				// 		std::cerr << "[WARNING] Too many root mobilities\n";
				// 	}
				// }

				// Iterate Gmolmodel bonds
				for(size_t bCnt = 0; bCnt < bonds.size(); bCnt++){

					bBond& bond = bonds[bCnt];

					if(bond.isThisMe(index_1, index_2) ){

						// Get this bond's topology and Comopund index
						Topology& topology = topologies[bond.getMoleculeIndex()];
						SimTK::Compound::BondIndex compoundBondIx = bond.getBondIndex();

						// Get the bond mobility
						auto it = mobilityMap.find(mobility);
						BondMobility::Mobility bondMobility = (it != mobilityMap.end()) ? it->second : BondMobility::Rigid;
						std::cout << "Setting bond mobility " << mobility << " " << bondMobility << std::endl;

						// Set the bond mobility
						bond.setBondMobility(bondMobility, whichWorld);
						topology.setBondMobility(bondMobility, compoundBondIx);

					} // found the bond

				} // search every bond
			} // more than two words per line

			else{
				scout("Bad flex file format");
			} // less than two words per line

		} // line is not empty
	} // every line

}

/*!
 * <!-- 
 * Careful: topology index comes first -->
*/
void Context::modelOneEmbeddedTopology(
	int whichTopology,
	int whichWorld
	//,std::string rootMobilizer
	)
{

	// Add a root mobilizer
	//this->rootMobilities.push_back(rootMobilizer);

	// Call Molmodel to model Compound

	std::vector<SimTK::Transform>& atomFrameCache =
		topologies[whichTopology].atomFrameCache;

	worlds[whichWorld].compoundSystem->modelOneCompound(
		SimTK::CompoundSystem::CompoundIndex(whichTopology),
		atomFrameCache,
		SimTK::String(rootMobilities[whichWorld][whichTopology])
		);

	// Get the forcefield within this world
	SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[whichWorld].forceField);

	// For every atom
	size_t topoNatoms = topologies[whichTopology].getNumAtoms();
	for(std::size_t aCnt = 0; aCnt < topoNatoms; aCnt++){

		// Get atom's mobod
		SimTK::Compound::AtomIndex aIx =
			(topologies[whichTopology].subAtomList[aCnt]).getCompoundAtomIndex();

		SimTK::MobilizedBodyIndex mbx =
			topologies[whichTopology].getAtomMobilizedBodyIndex(aIx);

		// SimTK::MobilizedBodyIndex mbxCheck =
		// 	topologies[whichTopology].getAtomMobilizedBodyIndexThroughDumm(
		// 		aIx,
		// 		dumm);

		// scout("world topology atom aIx")
		// 	<< whichWorld <<" " << whichTopology <<" "
		// 	<< aCnt <<" " << aIx <<" "
		// 	<< mbx <<" "
		// 	// << mbxCheck <<" " 
			// << eol;

	}

}


// ============================================================================
// WORK
// ============================================================================

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

// Print Molmodel related information
void Context::PrintMolmodelAndDuMMTypes(){
	for(std::size_t worldIx = 0; worldIx < nofWorlds; worldIx++){
		std::cout << "Context::PrintMolmodelAndDuMMTypes world " << worldIx << "\n";
		for(std::size_t molIx = 0; molIx < nofMols; molIx++){
			std::cout << "Context::PrintMolmodelAndDuMMTypes molecule " << molIx << "\n";
			const Topology& topology = worlds[worldIx].getTopology(molIx);
			SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[worldIx].forceField);
			topology.PrintMolmodelAndDuMMTypes(dumm);
		}
	}
}

// Print Simbody related information
void Context::PrintSimbodyMobods(){
	for(std::size_t worldIx = 0; worldIx < nofWorlds; worldIx++){
		std::cout << "Context::PrintSimbodyMobods world " << worldIx << "\n";
		for(std::size_t molIx = 0; molIx < nofMols; molIx++){
			std::cout << "Context::PrintSimbodyMobods molecule " << molIx << "\n";
			const Topology& topology = worlds[worldIx].getTopology(molIx);

			for(std::size_t i = 0; i < topology.getNumAtoms(); i++){
				SimTK::Compound::AtomIndex aIx
					= (topology.subAtomList[i]).getCompoundAtomIndex();
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(aIx);
				std::cout << "i aIx mbx " << i << " " << aIx << " "
					<< mbx << std::endl << std::flush;
			}
		}
	}
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

int Context::getNofMolecules()
{
	return nofMols;
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
	for(std::size_t molIx = 0; molIx < nofMols; molIx++){

		// Aquire the CompoundSystem
		topologies[molIx].setMultibodySystem(*worlds[newWorldIx].compoundSystem);

		int presumNAtoms = topologies[molIx].getNumAtoms();

		// Reset mobilized body indeces in Compound
		for(std::size_t k = 0; k < presumNAtoms; k++){

			// Get atom index
			SimTK::Compound::AtomIndex aIx =
				(topologies[molIx].subAtomList[k]).getCompoundAtomIndex();

			// Get mobod
			SimTK::MobilizedBodyIndex mbx =
				topologies[molIx].getAtomMobilizedBodyIndexThroughDumm(aIx,
				*(worlds[newWorldIx].forceField) );

			// Set the corresponding mobod
			topologies[molIx].setAtomMobilizedBodyIndex(aIx, mbx);

		}

	}
}


////////////////////////
// REX
////////////////////////

// Set the number of replicas. This could be dangerous
void Context::setNofReplicas(const size_t& argNofReplicas)
{
	nofReplicas = argNofReplicas;
}

// Adds a replica to the vector of Replica objects and
// does not set the coordinates of the replica's atomsLocations !
//void Context::addReplica(int index)
//{
//	replicas.emplace_back(Replica(index));
//	nofReplicas++;
//	replicas.back().Print();
//}

//// Adds a replica to the vector of Replica objects and sets the coordinates
//// of the replica's atomsLocations
//void Context::addReplica(int index,
//		std::vector<int>& argWorldIndexes)
//{
//	// Add replica and the vector of worlds
//	replicas.emplace_back(Replica(index, argWorldIndexes));
//
//	// Set replicas coordinates
//	std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>
//		referenceAtomsLocations =
//		worlds[0].getCurrentAtomsLocationsInGround();
//
//	replicas.back().setAtomsLocationsInGround(referenceAtomsLocations);
//
//	nofReplicas++;
//}

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
		, internCoords
		, zMatrixTable));

	// Set replicas coordinates

    // std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>> referenceAtomsLocationsFromFile;
	// replicas.back().setAtomsLocationsInGround(referenceAtomsLocationsFromFile);
	// replicas.back().set_WORK_AtomsLocationsInGround(referenceAtomsLocationsFromFile);

	const auto& referenceAtomsLocations = worlds[0].getCurrentAtomsLocationsInGround();
	replicas.back().setAtomsLocationsInGround(referenceAtomsLocations);
	replicas.back().set_WORK_AtomsLocationsInGround(referenceAtomsLocations);

	// Increment nof replicas
	nofReplicas++;
}

void Context::addThermodynamicState(int index,
		SimTK::Real T,
		std::vector<std::string>& rexSamplers,
		std::vector<int>& rexDistortOptions,
		std::vector<std::string>& rexDistortArgs,
		std::vector<int>& rexFlowOptions,
		std::vector<int>& rexWorkOptions,
		std::vector<std::string>& rexIntegrators,
		std::vector<int>& argWorldIndexes,
		std::vector<SimTK::Real>& timestepsInThisReplica,
		std::vector<int>& mdstepsInThisReplica)
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
	thermodynamicStates.back().setSamplers(rexSamplers);

	// Set non-equilibrium params
	thermodynamicStates.back().setDistortOptions(rexDistortOptions);
	thermodynamicStates.back().setDistortArgs(rexDistortArgs);
	thermodynamicStates.back().setFlowOptions(rexFlowOptions);
	thermodynamicStates.back().setWorkOptions(rexWorkOptions);

	// Set integrating method
	thermodynamicStates.back().setIntegrators(rexIntegrators);

	// Done
	nofThermodynamicStates++;
}

// Get the number of replicas
const size_t& Context::getNofReplicas() const
{
	return nofReplicas;
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
	std::fill(nofAttemptedSwapsMatrix.begin(), nofAttemptedSwapsMatrix.end(),
		vector<int>(nofThermodynamicStates, 0));

	// Allocate the number of accepted swaps
	nofAcceptedSwapsMatrix.resize(nofThermodynamicStates);
	for(size_t i = 0; i < nofThermodynamicStates; i++){
		nofAcceptedSwapsMatrix[i].resize(nofThermodynamicStates);
	}

	// Fill with zeros
	std::fill(nofAcceptedSwapsMatrix.begin(), nofAcceptedSwapsMatrix.end(),
		vector<int>(nofThermodynamicStates, 0));

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
	SimTK::Real Ui = replicas[replica_i].getFixman();

	if (replica_i == replica_j){ // same replica
		return Ui;
	}else if(Ui <= 0.0000001){ // fully flexible world
		return Ui;
	}else{
	    //std::cout << "CALC FIXMAN" << replica_i << " replica " << replica_j << "\n" << std::flush;
		// Get replica i thermodynamic state
        int thermoState_i = replica2ThermoIxs[replica_i];

		// Get replica i back world
        int world_i_front = thermodynamicStates[thermoState_i].getWorldIndexes().front();
        int world_i_back = thermodynamicStates[thermoState_i].getWorldIndexes().back();

		// Get coordinates from replica j
		const std::vector<std::vector<
            std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		X_j = replicas[replica_j].getAtomsLocationsInGround();

        // Pass compounds to the new world
        passTopologiesToNewWorld(world_i_back);

		// Transfer coordinates from j to i
        SimTK::State& state = (worlds[world_i_back].integ)->updAdvancedState();
        worlds[world_i_back].setAtomsLocationsInGround_REFAC(state, X_j);

		// Calculate Fixman in replica i back world
		SimTK::Real Fixman = worlds[world_i_back].calcFixman();

		// Transfer buffer coordinates of replica i
		// back to back world
		restoreReplicaCoordinatesToBackWorld(replica_i);

		passTopologiesToNewWorld(world_i_front);

		// Return
		return Fixman;
	}

}


// Calculate Fixman potential of replica J in replica I's back world. Ui(X_j)
SimTK::Real Context::calcFixman_JinI(int replica_i, int replica_j)
{
	SimTK::Real U_i = replicas[replica_i].getFixman();

	if (replica_i == replica_j){ // same replica
		return U_i;

	}else if(U_i <= 0.0000001){ // fully flexible world
		return U_i;

	}else{

		// Get replica i thermodynamic state
       	int thermoState_i = replica2ThermoIxs[replica_i];

		// Get replica i back world
		int world_i_front = thermodynamicStates[thermoState_i].getWorldIndexes().front();
       	int world_i_back = thermodynamicStates[thermoState_i].getWorldIndexes().back();

		// Get coordinates from replica j
		const std::vector<std::vector<
            std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		X_j = replicas[replica_j].getAtomsLocationsInGround();

        // Pass compounds to the new world
        passTopologiesToNewWorld(world_i_back);

		// Transfer coordinates from j to i
        SimTK::State& state = (worlds[world_i_back].integ)->updAdvancedState();
        worlds[world_i_back].setAtomsLocationsInGround_REFAC(state, X_j);

		// Calculate Fixman in replica i back world
		SimTK::Real Fixman = worlds[world_i_back].calcFixman();

		// Transfer buffer coordinates of replica i back to back world
		restoreReplicaCoordinatesToBackWorld(replica_i);

		passTopologiesToNewWorld(world_i_front);

		// Return
		return Fixman;
	}

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
		const std::vector<std::vector<
            std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		X_i = replicas[replica_i].getAtomsLocationsInGround();

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

/*! <!-- restoreReplica --> */\
void Context::rewindReplica(void){
	assert(!"Not implemented");

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

	SimTK::Real W_X = replicas[replica_X].getWORK(); // work without Jacobian
	SimTK::Real W_Y = replicas[replica_Y].getWORK(); // work without Jacobian

	SimTK::Real U_Xtau = replicas[replica_X].get_WORK_PotentialEnergy_New(); // last non-equil potential
	SimTK::Real U_Ytau = replicas[replica_Y].get_WORK_PotentialEnergy_New(); // last non-equil potential

	SimTK::Real lnJac_X = replicas[replica_X].get_WORK_Jacobian(); // non-equil Jacobian
	SimTK::Real lnJac_Y = replicas[replica_Y].get_WORK_Jacobian(); // non-equil Jacobian

	std::cout <<" STUDY relationships"
		<<" "<< "U_X0 U_Y0 W_X W_Y U_Xtau U_Ytau lnJac_X lnJAC_Y"
		<<" "<< U_X0 <<" "<< U_Y0 <<" "<< W_X <<" "<< W_Y <<" "<< U_Xtau <<" "<< U_Ytau <<" "<< lnJac_X <<" "<< lnJac_Y
		//<<" "<< (U_X0 + W_X) - U_Xtau <<" "<< (U_Y0 + W_Y) - U_Ytau
		<< std::endl;

	assert( ((std::abs((U_X0 + W_X) - U_Xtau) < 0.00001) &&
			 (std::abs((U_Y0 + W_Y) - U_Ytau) < 0.00001)) &&
			"Work seems incorrectly calculated");

	if( ((std::abs((U_X0 + W_X) - U_Xtau) > 0.00001) ||
		 (std::abs((U_Y0 + W_Y) - U_Ytau) > 0.00001)) ){
		std::cerr << "Work seems incorrectly calculated\n";
	}

	#pragma endregion convienent_vars

	// ----------------------------------------------------------------
	// Reduced potentials X0
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real uC_X0 = beta_C * U_X0; // Replica i reduced potential in state i
	SimTK::Real uH_Y0 = beta_H * U_Y0; // Replica j reduced potential in state j
	SimTK::Real uH_X0 = beta_H * U_X0; // Replica i reduced potential in state j
	SimTK::Real uC_Y0 = beta_C * U_Y0; // Replica j reduced potential in state i

	// ----------------------------------------------------------------
	// Reduced potential Xtau
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real uC_Xtau = beta_C * U_Xtau; // Replica i reduced potential in state i
	SimTK::Real uH_Ytau = beta_H * U_Ytau; // Replica j reduced potential in state j
	SimTK::Real uH_Xtau = beta_H * U_Xtau; // Replica i reduced potential in state j
	SimTK::Real uC_Ytau = beta_C * U_Ytau; // Replica j reduced potential in state i

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
	SimTK::Real ETerm_equil    = uH_X0 - uC_X0;
				ETerm_equil   += uC_Y0 - uH_Y0;
				ETerm_equil = -1.0 * ETerm_equil;

	// ----------------------------------------------------------------
	// LOGP ENERGY NON-EQUILIBRIUM
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real ETerm_nonequil  = uH_Xtau - uC_Xtau;
			    ETerm_nonequil += uC_Ytau - uH_Ytau;
				ETerm_nonequil = -1.0 * ETerm_nonequil;

	// ----------------------------------------------------------------
	// LOGP WORK
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Get work from X replica
	SimTK::Real Work_X = (uH_Xtau - uC_X0) + (fixH_Xtau - fixC_X0) - lnJac_X;

	// Get work from Y replica
	SimTK::Real Work_Y = (uC_Ytau - uH_Y0) + (fixC_Ytau - fixH_Y0) - lnJac_Y;

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
	bool printTerms = true, printWithoutText = true;
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
			<< beta_C << ", "
			<< beta_H << ", "
			<< uC_X0 << ", " << uH_Y0 << ", " << uH_X0 << ", " << uC_Y0 << ", "
			<< uC_Xtau << ", " << uH_Ytau << ", " << uH_Xtau << ", " << uC_Ytau << ", "
			<< lnJac_X << ", "
			<< lnJac_Y << ", "
			<< W_X << ", "
			<< W_Y << ", "
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
	if(getRunType() == RUN_TYPE::REMC){

		log_p_accept = ETerm_equil ;

	}else if(getRunType() == RUN_TYPE::RENEMC){

		log_p_accept = ETerm_nonequil + std::log(correctionTerm) ;

	}else if( getRunType() == RUN_TYPE::RENE){

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
		
		TestingWay testingWay = TestingWay::ALWAYS_REJECT;  				// ALWAYS_REJECT

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

		if( getRunType() == RUN_TYPE::RENE){

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

		// Update replicas coordinates from work generated coordinates
		set_WORK_CoordinatesAsFinal(replica_X);
		set_WORK_CoordinatesAsFinal(replica_Y);

		// Update replica's energy from work last potential energy
		set_WORK_PotentialAsFinal(replica_X);
		set_WORK_PotentialAsFinal(replica_Y);

		// Swap thermodynamic states
		swapThermodynamicStates(replica_X, replica_Y);
		swapPotentialEnergies(replica_X, replica_Y);

		std::cout << "1" 
		<<", " << unifSample 
		<< endl << endl;

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
		<< endl << endl;

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
	rexDetHeader 
		<< "REXdetails" << ", " << "thermoState_C" << ", " << "thermoState_H" << ", "
		<< "replica_X" << ", " << "replica_Y" << ", "
		<< "betaC" << ", " << "beatH" << ", "
		<< "eC_X0" << ", " << "eH_Y0" << ", " << "eH_X0" << ", " << "eC_Y0" << ", "
		<< "lC_Xtau" << ", " << "lH_Ytau" << ", " << "lH_Xtau" << ", " << "lC_Ytau" << ", "
		<< "wJ_X" << ", " << "wJ_Y" << ", " << "wE_X" << ", " << "wE_Y" << ", "
		<< ", " << "s_X" << ", " << "s_Y" << ", " << "s_X_1" << ", " << "s_Y_1" << ", "
		<< ", " << "qC_s_X" << ", " << "qH_s_Y" << ", " << "qH_s_X_1" << ", " << "qC_s_Y_1" << ", "
		<< "ETerm_equil" << ", " << "WTerm" << ", " << "correctionTerm" << ", " << "swapped"
	;
}

// Mix neighboring replicas // TODO: delete
// Thermodyanmic states are fixed; replicas are variables
void Context::mixNeighboringReplicas(unsigned int startingFrom)
{
	assert((startingFrom <= 1) &&
	"Replica exchange scheme has to start from 0 or 1.");

	int thermoState_i = 0;
	int thermoState_j = 1;

	// Go through neighboring thermodynamic states
	for(size_t thermoState_k = startingFrom;
	thermoState_k < (nofThermodynamicStates - 1);
	thermoState_k += 2){

		// Get thermodynamic states
		thermoState_i = thermoState_k;
		thermoState_j = thermoState_k + 1;

		// Get replicas corresponding to the thermodynamic states
		int replica_i = thermo2ReplicaIxs[thermoState_i];
		int replica_j = thermo2ReplicaIxs[thermoState_j];

		// Attempt to swap
		bool swapped = false;
		swapped = attemptREXSwap(replica_i, replica_j);

	}


}


// Exhange all replicas
void Context::mixAllReplicas(int nSwapAttempts)
{
	std::uniform_int_distribution<std::size_t> randReplicaDistrib(0, nofReplicas-1);

	// Try nSwapAttempts to swap between random replicas
	for(size_t swap_k = 0; swap_k < nSwapAttempts; swap_k++){

		// Get two random replicas
		auto replica_i = randReplicaDistrib(randomEngine);
		auto replica_j = randReplicaDistrib(randomEngine);

		std::cout << "Attempt to swap replicas " << replica_i
			<< " and " << replica_j << std::endl;

		// Attempt to swap
		bool swapped = false;

		swapped = attemptREXSwap(replica_i, replica_j);

	}
}

/*!
 * <!--  -->
*/
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

	SimTK::State& state =
		(worlds[currWorldIx].integ)->updAdvancedState();

	//worlds[currWorldIx].setAtomsLocationsInGround(state,
	//	replicas[whichReplica].getAtomsLocationsInGround());
	state = setAtoms_SP_NEW(currWorldIx, state,
		replicas[whichReplica].getAtomsLocationsInGround());		

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

/*!
 * <!--	 -->
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
void Context::setReplicasWorldsParameters(int thisReplica)
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
		worlds[replicaWorldIxs[i]].updSampler(0)->setSampleGenerator(
			thermodynamicStates[thisThermoStateIx].getSamplers()[i]
		);
	}

	// -------------
	// Set simulation parameters
	// =============

	// Set integrator
	for(std::size_t i = 0; i < replicaNofWorlds; i++){

		worlds[replicaWorldIxs[i]].updSampler(0)->setIntegratorName(
			thermodynamicStates[thisThermoStateIx].getIntegrators()[i]
		);
	}

	// Set timestep and nof MD steps
	const std::vector<SimTK::Real>& replicaTimesteps =
		thermodynamicStates[thisThermoStateIx].getTimesteps();
	const std::vector<int>& replicaMdsteps =
		thermodynamicStates[thisThermoStateIx].getMdsteps();

	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setTimestep(replicaTimesteps[i], false);
	}

	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setMDStepsPerSample(
			replicaMdsteps[i]);
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

RUN_TYPE Context::getRunType(void) const
{
	return runType;
}

void Context::setRunType(RUN_TYPE runTypeArg)
{
	this->runType = runTypeArg;
}

RUN_TYPE Context::setRunType(const std::string& runTypeArgStr)
{
	if(runTypeArgStr == "Normal"){
		this->runType = RUN_TYPE::DEFAULT;
	}else if(runTypeArgStr == "SimulatedTempering") {
		this->runType = RUN_TYPE::DEFAULT;
	}else if(runTypeArgStr == "REMC"){
		this->runType = RUN_TYPE::REMC;
	}else if(runTypeArgStr == "RENEMC"){
		this->runType = RUN_TYPE::RENEMC;
	}else if(runTypeArgStr == "RENE"){
		this->runType = RUN_TYPE::RENE;
	}else{
		this->runType = RUN_TYPE::DEFAULT;
	}

	return this->runType;

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
		// const SimTK::Vector & BMps = worlds[srcStatsWIx].getBMps();
	worlds[destStatsWIx].updSampler(0)->set_dBMps(thermodynamicStates[thermoIx].get_dBMps(srcStatsWIx));

	worlds[destStatsWIx].updSampler(0)->setPreviousQs(thermodynamicStates[thermoIx].getCurrentQs(srcStatsWIx));
	worlds[destStatsWIx].updSampler(0)->setQmeans(thermodynamicStates[thermoIx].getQmeans(srcStatsWIx));
	worlds[destStatsWIx].updSampler(0)->setQdiffs(thermodynamicStates[thermoIx].getQdiffs(srcStatsWIx));
	worlds[destStatsWIx].updSampler(0)->setQvars(thermodynamicStates[thermoIx].getQvars(srcStatsWIx));
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
		validated = worlds[whichWorld].generateSamples(
			numSamples, worldOutStream, header);

	// Non-equilibrium world
	} else if (distortOption != 0) {

		// Generate samples
		
		// drl
		#ifdef __DRILLING__ // SCALEQ

			// Get drl data
			const std::vector<std::vector<double>>& drl_bon_Energies = worlds[whichWorld].getEnergies_drl_bon();
			const std::vector<std::vector<double>>& drl_ang_Energies = worlds[whichWorld].getEnergies_drl_ang();
			const std::vector<std::vector<double>>& drl_tor_Energies = worlds[whichWorld].getEnergies_drl_tor();

			// validated = worlds[whichWorld].generateSamples(numSamples, worldOutStream){

				//warn("under drilling conditions");

				// Update Robosample bAtomList
				SimTK::State& currentAdvancedState = (worlds[whichWorld]).integ->updAdvancedState();
				(worlds[whichWorld]).updateAtomListsFromCompound(currentAdvancedState);

				// ''''''''''''''''''''
				// coutspaced("SCALING_BAT init:"); ceolf;
				// replicas[0].calcZMatrixBAT( (worlds[whichWorld]).getAtomsLocationsInGround( (worlds[whichWorld]).integ->updAdvancedState() ));
				// thermodynamicStates[0].PrintZMatrixBAT();
				// ''''''''''''''''''''

				// Reinitialize the sampler
				validated = (worlds[whichWorld]).updSampler(0)->reinitialize(currentAdvancedState,
					worldOutStream);

				SimTK::Real pe_beforeScale = (worlds[whichWorld]).forces->getMultibodySystem().calcPotentialEnergy((worlds[whichWorld]).integ->updAdvancedState());
				// scout("[SCALING_PES]: before") <<" " << pe_beforeScale << eolf;
				// scout("drl_bon_E"); ceol; PrintCppVector(drl_bon_Energies);
				// scout("drl_ang_E"); ceol; PrintCppVector(drl_ang_Energies);
				// scout("drl_tor_E"); ceol; PrintCppVector(drl_tor_Energies);

				// GENERATE the requested number of samples
				for(int k = 0; k < numSamples; k++) {

					worldOutStream << header << " ";
					(worlds[whichWorld]).updSampler(0)->getMsg_InitialParams(worldOutStream);

					validated = (worlds[whichWorld]).updSampler(0)->sample_iteration(
						currentAdvancedState, worldOutStream) && validated;

					worldOutStream << "\n";

				}

				SimTK::Real pe_afterScale = (worlds[whichWorld]).forces->getMultibodySystem().calcPotentialEnergy((worlds[whichWorld]).integ->updAdvancedState());

				// ''''''''''''''''''''
				// coutspaced("SCALING_BAT after:"); ceolf;
				// replicas[0].calcZMatrixBAT( (worlds[whichWorld]).getAtomsLocationsInGround( (worlds[whichWorld]).integ->updAdvancedState() ));
				// thermodynamicStates[0].PrintZMatrixBAT();
				// ''''''''''''''''''''

				// scout("[SCALING_PES]: after") <<" " << pe_afterScale << eolf;
				// scout("drl_bon_E"); ceol; PrintCppVector(drl_bon_Energies);
				// scout("drl_ang_E"); ceol; PrintCppVector(drl_ang_Energies);
				// scout("drl_tor_E"); ceol; PrintCppVector(drl_tor_Energies);

			// }

		#else

			validated = worlds[whichWorld].generateSamples(
				numSamples, worldOutStream, header); // =

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
	worldOutStream << "\n";
	std::cout << worldOutStream.str();

	return validated;
}

/*!
 * <!-- Run a vector of worlds -->
*/ 
void Context::RunWorlds(std::vector<int>& specificWIxs, int replicaIx)
{
	int thermoIx = replica2ThermoIxs[replicaIx];
	bool validated = true;
	for(std::size_t spWCnt = 0; spWCnt < specificWIxs.size() - 1; spWCnt++){ // -1 so we can transfer

		// Run
		int srcStatsWIx  = specificWIxs[spWCnt];

		//std::cout << "REX, " << replicaIx << ", " << thermoIx << " , " << srcStatsWIx;
		std::string headerToRunWorld = "REX, " + std::to_string(replicaIx)
					+ ", " + std::to_string(thermoIx)
					+ " , " + std::to_string(srcStatsWIx);
		validated = RunWorld(srcStatsWIx, headerToRunWorld ) && validated;

		// Calculate Q statistics ^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&
		// worlds[srcStatsWIx].PrintXBMps(); // @@@@@@@@@@@@@
		// const SimTK::Vector & BMps = worlds[srcStatsWIx].getBMps();
		// for(int mbx = 1; mbx < worlds[srcStatsWIx].matter->getNumBodies(); mbx++){
		// 	std::cout <<" " << BMps[mbx] ;
		// }
		// ^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&

		if( pHMC((worlds[srcStatsWIx].samplers[0]))->getAcc() == true){
			thermodynamicStates[thermoIx].calcQStats(
				srcStatsWIx, worlds[srcStatsWIx].getBMps(), worlds[srcStatsWIx].getAdvancedQs(), worlds[srcStatsWIx].getNofSamples());
		}else{
			thermodynamicStates[thermoIx].calcQStats(
				srcStatsWIx, worlds[srcStatsWIx].getBMps(), SimTK::Vector(worlds[srcStatsWIx].getNQs(), SimTK::Real(0)), worlds[srcStatsWIx].getNofSamples());
		}
		
		// Transfer coordinates to the next world
		int destStatsWIx = specificWIxs[spWCnt + 1];
		transferCoordinates_WorldToWorld(srcStatsWIx, destStatsWIx);
		transferQStatistics(thermoIx, srcStatsWIx, destStatsWIx);

		// // Calculate replica BAT and BAT stats
		// World& currWorld = worlds[specificWIxs[spWCnt]];
		// SimTK::State& currState = currWorld.integ->updAdvancedState();
		// replicas[replicaIx].calcZMatrixBAT( currWorld.getAtomsLocationsInGround( currState ));

	}

	// Run the last world
	int srcStatsWIx = specificWIxs.back();

	//std::cout << "REX, " << replicaIx << ", " << thermoIx << ", " << specificWIxs.back();
	std::string headerToRunWorld = "REX, " + std::to_string(replicaIx) 
				+ ", " + std::to_string(thermoIx)
				+ ", " + std::to_string(srcStatsWIx);
	validated = RunWorld(srcStatsWIx, headerToRunWorld) && validated;

	// Calculate Q statistics ^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&
	// worlds[srcStatsWIx].PrintXBMps(); // @@@@@@@@@@@@@
	// const SimTK::Vector & BMps = worlds[srcStatsWIx].getBMps();
	// for(int mbx = 1; mbx < worlds[srcStatsWIx].matter->getNumBodies(); mbx++){
	// 	std::cout <<" " << BMps[mbx] ;
	// } // @@@@@@@@@@@@@

	if( pHMC((worlds[srcStatsWIx].samplers[0]))->getAcc() == true){
		thermodynamicStates[thermoIx].calcQStats(
			srcStatsWIx, worlds[srcStatsWIx].getBMps(), worlds[srcStatsWIx].getAdvancedQs(), worlds[srcStatsWIx].getNofSamples());
	}else{
		thermodynamicStates[thermoIx].calcQStats(
			srcStatsWIx, worlds[srcStatsWIx].getBMps(), SimTK::Vector(worlds[srcStatsWIx].getNQs(), SimTK::Real(0)), worlds[srcStatsWIx].getNofSamples());
	}

	if(true){
		World& currWorld = worlds[specificWIxs.back()];
		SimTK::State& currState = currWorld.integ->updAdvancedState();
		replicas[replicaIx].calcZMatrixBAT( currWorld.getAtomsLocationsInGround( currState ));
	}

	// #ifndef PRINTALOT
	// #define PRINTALOT
	// #endif

	#ifdef PRINTALOT 
		if(validated){
			std::cout << std::endl;
		}else{
			std::cout << " invalid sample." << std::endl;
		}
	#endif
}

/*!
 * <!--  -->
*/
void Context::RunReplicaRefactor(
	int mixi,
	int replicaIx)
{

	assert(!"Context::RunReplicaRefactor function will be removed.");

	/*
	// MOVE THIS INTO THERMODYNAMIC STATES
	// Get equilibrium and non-equilibrium worlds
	int thermoIx = replica2ThermoIxs[replicaIx];
	std::vector<int>& replicaWorldIxs = thermodynamicStates[thermoIx].updWorldIndexes();
	size_t replicaNofWorlds = replicaWorldIxs.size();
	std::vector<int> & distortOpts = thermodynamicStates[thermoIx].getDistortOptions();

	std::vector<int> equilWIxs;
	std::vector<int> nonequilWIxs;
	equilWIxs.reserve(distortOpts.size());
	nonequilWIxs.reserve(distortOpts.size());

	for(int dIx = 0; dIx < distortOpts.size(); dIx++){
		if(distortOpts[dIx] == 0){
			equilWIxs.push_back(dIx);
		}else{
			nonequilWIxs.push_back(dIx);
		}
	}

	// Set world samples incrementation
	int wSamIncr = 0;
	int samIncr = 0;
	(getRunType() == RUN_TYPE::RENE) ? wSamIncr = replicaNofWorlds - 1 : replicaNofWorlds;
	(getRunType() == RUN_TYPE::RENE) ? samIncr = 1 : samIncr = 1; // Always increment replica's nof samples

	//scout("\nEquil world indexes") ;if(equilWIxs.size() > 0){PrintCppVector(equilWIxs);}scout("\nNonequil world indexes") ;if(nonequilWIxs.size() > 0){PrintCppVector(nonequilWIxs);}ceolf;

	if((equilWIxs.size() == 0) && (nonequilWIxs.size() == 0)){
		warn("No worlds");
		return;
	}

	// Transfer Q statistics ^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&
	transferQStatistics(thermoIx, replicaWorldIxs.back(), replicaWorldIxs.front());
	// ^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&^&

	if(!equilWIxs.empty()){ // Equilibrium

		RunWorlds(equilWIxs, replicaIx);
		
		replicas[replicaIx].updAtomsLocationsInGround(worlds[equilWIxs.back()].getCurrentAtomsLocationsInGround());
	
		replicas[replicaIx].setPotentialEnergy(worlds[equilWIxs.back()].calcPotentialEnergy());

		REXLog(mixi, replicaIx);

		if(!nonequilWIxs.empty()){ // Non-Equilibrium

			transferCoordinates_WorldToWorld(equilWIxs.back(), nonequilWIxs.front());
			transferQStatistics(thermoIx, equilWIxs.back(), nonequilWIxs.front());

			RunWorlds(nonequilWIxs, replicaIx);
		
			replicas[replicaIx].setWORK( calcReplicaWork(replicaIx) ); // possibly not correct
		
			replicas[replicaIx].upd_WORK_AtomsLocationsInGround(worlds[nonequilWIxs.back()].getCurrentAtomsLocationsInGround());
		
			replicas[replicaIx].set_WORK_PotentialEnergy_New(worlds[nonequilWIxs.back()].calcPotentialEnergy());
		
			SimTK::Real jac = (worlds[nonequilWIxs.back()].updSampler(0))->getDistortJacobianDetLog();
			replicas[replicaIx].set_WORK_Jacobian(jac);

			//std::cout << "Set Jacobian for replica " << replicaIx <<" to " << jac << std::endl;
		
			SimTK::Real fix_set_back = pHMC((worlds[nonequilWIxs.back()].samplers[0]))->fix_set;
			replicas[replicaIx].setFixman(fix_set_back);

			transferCoordinates_WorldToWorld(nonequilWIxs.back(), equilWIxs.front());

		}else{

			replicas[replicaIx].upd_WORK_AtomsLocationsInGround(worlds[equilWIxs.back()].getCurrentAtomsLocationsInGround()); // Victor bugfix

			transferCoordinates_WorldToWorld(equilWIxs.back(), equilWIxs.front());

		}

		// Increment the nof samples for replica and thermostate
		replicas[replicaIx].incrementWorldsNofSamples(wSamIncr);
		thermodynamicStates[thermoIx].incrementWorldsNofSamples(wSamIncr);
		replicas[replicaIx].incrementNofSamples(samIncr);
		thermodynamicStates[thermoIx].incrementNofSamples(samIncr);		

		return;
	}
		
	if(!nonequilWIxs.empty()){ // Only Non-Equilibrium

		RunWorlds(nonequilWIxs, replicaIx);
	
		replicas[replicaIx].setWORK( calcReplicaWork(replicaIx) ); // possibly not correct
	
		replicas[replicaIx].upd_WORK_AtomsLocationsInGround(worlds[nonequilWIxs.back()].getCurrentAtomsLocationsInGround());
	
		replicas[replicaIx].set_WORK_PotentialEnergy_New(worlds[nonequilWIxs.back()].calcPotentialEnergy());
	
		SimTK::Real jac = (worlds[nonequilWIxs.back()].updSampler(0))->getDistortJacobianDetLog();
		replicas[replicaIx].set_WORK_Jacobian(jac);
	
		SimTK::Real fix_set_back = pHMC((worlds[nonequilWIxs.back()].samplers[0]))->fix_set;
		replicas[replicaIx].setFixman(fix_set_back);

		transferCoordinates_WorldToWorld(nonequilWIxs.back(), nonequilWIxs.front());

	}

	// Increment the nof samples for replica and thermostate
	replicas[replicaIx].incrementWorldsNofSamples(wSamIncr);
	thermodynamicStates[thermoIx].incrementWorldsNofSamples(wSamIncr);

	replicas[replicaIx].incrementNofSamples(0); // Only non-equilibrium
	thermodynamicStates[thermoIx].incrementNofSamples(0); // Only non-equilibrium
	*/

}


/*! <!--  -->
*/
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

	replica.updWORK() = 0.0;  // TODO merge with Jacobians
	replica.upd_WORK_Jacobian() = 0.0;

	// Main loop
	for(std::size_t thWCnt = 0; thWCnt < thermoNofWorlds; thWCnt++){

		int wIx  = thermoWorldIxs[thWCnt];
		int distortIx = distortOpts[thWCnt];
		World& currWorld = worlds[wIx];
		HMCSampler* sampler_p = worlds[wIx].samplers[0].get();

		// Transfer coordinates to the next world
		if(thWCnt == 0){
			transferCoordinates_ReplicaToWorld(replicaIx, thermoWorldIxs.front()); //state = setAtoms_SP_NEW(thermoWorldIxs.front(), state, replica.getAtomsLocationsInGround());
			transferQStatistics(thermoIx, thermoWorldIxs.back(), thermoWorldIxs.front());
		}else{
			transferCoordinates_WorldToWorld(thermoWorldIxs[thWCnt - 1], wIx);
			transferQStatistics(thermoIx, thermoWorldIxs[thWCnt - 1], wIx);
		}

		// Header
		std::string headerToRunWorld = "REX";
					headerToRunWorld += ", " + std::to_string(replicaIx);
					headerToRunWorld += ", " + std::to_string(thermoIx);
					headerToRunWorld += ", " + std::to_string(wIx);

		// Run
		bool validated = true;
		validated = RunWorld(wIx, headerToRunWorld ) && validated;

		// Calculate Q statistics
		if( sampler_p->getAcc() == true){
			thermoState.calcQStats(wIx, currWorld.getBMps(), currWorld.getAdvancedQs(), currWorld.getNofSamples());
		}else{
			thermoState.calcQStats(wIx, currWorld.getBMps(), SimTK::Vector(currWorld.getNQs(), SimTK::Real(0)), currWorld.getNofSamples());
		}

		// ======================== EQUILIBRIUM ======================
		if(distortIx == 0){

			replica.updAtomsLocationsInGround(currWorld.getCurrentAtomsLocationsInGround());

			replica.setPotentialEnergy(currWorld.calcPotentialEnergy());

			replica.setFixman(sampler_p->fix_set); // DID_I_REALLY_FORGET

		// ======================== NON-EQUILIBRIUM ======================
		}else{

			replica.updWORK() += currWorld.getWork();  // TODO merge with Jacobians
			replica.upd_WORK_Jacobian() += sampler_p->getDistortJacobianDetLog();

			replica.upd_WORK_AtomsLocationsInGround(currWorld.getCurrentAtomsLocationsInGround());

			replica.set_WORK_PotentialEnergy_New(currWorld.calcPotentialEnergy());
		
			replica.set_WORK_Fixman(sampler_p->fix_set); // DID_I_REALLY_FORGET

		} // __end__ Non/Equilibrium =======================================

		# pragma region REBAS_TEST
		const SimTK::State& pdbState = currWorld.integ->updAdvancedState();
		currWorld.updateAtomListsFromSimbody(pdbState);

		// std::string pdbMiddle = pdbPrefix + "." + std::to_string(0) + "." + "s" + std::to_string(thermoIx) + "." + "w" + std::to_string(wIx) + ".";
		// std::cout << "Writing " << pdbMiddle << std::endl;

		// Write pdb
		if(pdbRestartFreq){
			if((mixi % pdbRestartFreq) == 0){
				if(wIx == 0){
					for(int mol_i = 0; mol_i < getNofMolecules(); mol_i++){
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

	} // __end__ Main loop

	replica.incrementNofSamples(1);
	thermoState.incrementNofSamples(1);

	SimTK::State& state = worlds.back().integ->updAdvancedState();
	replica.calcZMatrixBAT( worlds.back().getAtomsLocationsInGround( state ));

	transferCoordinates_WorldToWorld(thermoWorldIxs.back(), thermoWorldIxs.front());
	
}


/*!
 * <!-- Run replica exchange protocol -->
*/
void Context::RunREX()
{

    if(MEMDEBUG){stdcout_memdebug("Context::RunREX 1");}


	// desk_mass_related
	for (int worldIx = 0; worldIx < worlds.size(); worldIx++) {
		for (auto& atom : atoms) {
			
			World& currWorld = worlds[worldIx];
			//SimTK::DuMMForceFieldSubsystem dumm = *(currWorld.forceField);

			SimTK::DuMM::AtomIndex dAIx = atom.updDuMMAtomIndex();
			SimTK::mdunits::Mass atomMass = atom.getMass();

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
	for(size_t mixi = 0; mixi < nofMixes; mixi++){

		//std::cout << " REX batch " << mixi << std::endl;

		// Reset replica exchange pairs vector
		if(getRunType() != RUN_TYPE::DEFAULT){                                                 // (9) + (10)
			if(replicaMixingScheme == ReplicaMixingScheme::neighboring){
				setReplicaExchangePairs(mixi % 2);
			}
		}

		// Update work scale factors
		updThermostatesQScaleFactors(mixi);

		//Print_TRANSFORMERS_Work();

    	if(MEMDEBUG){stdcout_memdebug("Context::RunREX 4");}

		// SIMULATE EACH REPLICA --------------------------------------------->
		for (size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++){

			// Update BAT map for all the replica's world
			updSubZMatrixBATsToAllWorlds(replicaIx);

	    	if(MEMDEBUG){stdcout_memdebug("Context::RunREX 5");}

			// Load the front world
			currFrontWIx = restoreReplicaCoordinatesToFrontWorld(replicaIx);           // (1)

			// Set thermo and simulation parameters for the worlds in this replica
			setReplicasWorldsParameters(replicaIx);

	    	if(MEMDEBUG){stdcout_memdebug("Context::RunREX 6");}

			// ======================== SIMULATE ======================
			//RunReplicaRefactor(mixi, replicaIx);
			RunReplicaRefactor_SIMPLE(mixi, replicaIx);

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
		if((getRunType() != RUN_TYPE::DEFAULT) && (nofReplicas != 1)){
			
			mixReplicas(mixi); // check this
			                                               
			PrintNofAcceptedSwapsMatrix();
		}

		this->nofRounds++; 

		if (wantDCD) {
			writeDCDs();
		}

	} // end rounds

	//PrintNofAttemptedSwapsMatrix();
	PrintNofAcceptedSwapsMatrix();
	//PrintReplicaMaps();

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
		bSpecificAtom& atom0 = atoms[zRow[childPositionInZMat]];
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
		for(int mol_i = 0; mol_i < getNofMolecules(); mol_i++){
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

void Context::writeDCDs()
{
	// Write the DCD from the first world
	auto& world = worlds[worldIndexes.front()];

	// Update bAtomList in Topology
	const SimTK::State& pdbState = world.integ->updAdvancedState();
	world.updateAtomListsFromSimbody(pdbState);

	for (int i = 0; i < natoms; i++) {
		Xs[i] = atoms[i].getX() * 10;
		Ys[i] = atoms[i].getY() * 10;
		Zs[i] = atoms[i].getZ() * 10;
	}

	// Actual write
	traj.appendTimestep("dcd", Xs, Ys, Zs);
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

	// Get coordinates from source
	const std::vector<std::vector<std::pair<
		bSpecificAtom *, SimTK::Vec3> > >&
		otherWorldsAtomsLocations =
	worlds[srcWIx].getAtomsLocationsInGround(lastAdvancedState);

	// Get BAT coordinates
	calcZMatrixBAT(srcWIx, otherWorldsAtomsLocations);

	//PrintZMatrixTableAndBAT();
	//PrintZMatrixMobods(srcWIx, lastAdvancedState);
	
	// Pass compounds to the new world
	passTopologiesToNewWorld(destWIx);

	// Focus on destination world
	World& destWorld = worlds[destWIx];
	SimTK::State& someState = currentAdvancedState;

	// New setAtomsLocations
	currentAdvancedState = setAtoms_SP_NEW(destWIx, someState, otherWorldsAtomsLocations);

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

	// Get advanced states of the integrators
	SimTK::State& srcAdvancedState = worlds[srcWIx].integ->updAdvancedState();
	SimTK::State& destAdvancedState = worlds[destWIx].integ->updAdvancedState();

	// Get coordinates from source World
	const std::vector<std::vector<std::pair<
		bSpecificAtom *, SimTK::Vec3> > >&
		srcWorldsAtomsLocations =
	worlds[srcWIx].getAtomsLocationsInGround(srcAdvancedState);

	// Get coordinates from destintaion World
	const std::vector<std::vector<std::pair<
		bSpecificAtom *, SimTK::Vec3> > >&
		destWorldsAtomsLocations =
	worlds[destWIx].getAtomsLocationsInGround(destAdvancedState);

	if(srcWorldsAtomsLocations.size() != destWorldsAtomsLocations.size()){
		std::cout << 
		"Context::checkTransferCoordinates worldsAtomsLocations not the same size"
		<< std::endl;
		return 99999999;
	}
	for(size_t toIx = 0; toIx < srcWorldsAtomsLocations.size(); toIx++){
		if(srcWorldsAtomsLocations[toIx].size() != destWorldsAtomsLocations[toIx].size()){
			std::cout << 
			"Context::checkTransferCoordinates worldsAtomsLocations " << toIx
			<< " not the same size" << std::endl;
			return 99999999;
		}
	}

	// Get the max
	SimTK::Real vMaxComp = 0;
	SimTK::Real currDiff;
	SimTK::Real cumulDiff = 0;
	for(size_t toIx = 0; toIx < srcWorldsAtomsLocations.size(); toIx++){
		for(size_t atIx = 0; atIx < srcWorldsAtomsLocations[toIx].size(); atIx++){

			currDiff = std::abs((srcWorldsAtomsLocations[toIx][atIx]).second[0] - 
					 			(destWorldsAtomsLocations[toIx][atIx]).second[0]);
			if(currDiff > vMaxComp){vMaxComp = currDiff;}
			cumulDiff += currDiff;

			currDiff = std::abs((srcWorldsAtomsLocations[toIx][atIx]).second[1] - 
					 			(destWorldsAtomsLocations[toIx][atIx]).second[1]);
			if(currDiff > vMaxComp){vMaxComp = currDiff;}
			cumulDiff += currDiff;

			currDiff = std::abs((srcWorldsAtomsLocations[toIx][atIx]).second[2] - 
					 			(destWorldsAtomsLocations[toIx][atIx]).second[2]);
			if(currDiff > vMaxComp){vMaxComp = currDiff;}
			cumulDiff += currDiff;

		}
	}

	// for(size_t toIx = 0; toIx < srcWorldsAtomsLocations.size(); toIx++){
	// 	for(size_t atIx = 0; atIx < srcWorldsAtomsLocations[toIx].size(); atIx++){		
	// 		std::cout << "transfer " << toIx << " " << atIx
	// 			<< " " << double((srcWorldsAtomsLocations[toIx][atIx]).second[0])
	// 			<< " " << double((srcWorldsAtomsLocations[toIx][atIx]).second[1])
	// 			<< " " << double((srcWorldsAtomsLocations[toIx][atIx]).second[2])
	// 			<< " " << double((destWorldsAtomsLocations[toIx][atIx]).second[0])
	// 			<< " " << double((destWorldsAtomsLocations[toIx][atIx]).second[1])
	// 			<< " " << double((destWorldsAtomsLocations[toIx][atIx]).second[2])
	// 			<< " " << double((destWorldsAtomsLocations[toIx][atIx]).second[0]) - double((srcWorldsAtomsLocations[toIx][atIx]).second[0])
	// 			<< " " << double((destWorldsAtomsLocations[toIx][atIx]).second[1]) - double((srcWorldsAtomsLocations[toIx][atIx]).second[1])
	// 			<< " " << double((destWorldsAtomsLocations[toIx][atIx]).second[2]) - double((srcWorldsAtomsLocations[toIx][atIx]).second[2])									
	// 		<< std::endl;
	// 	}
	// }

	scout("Check coords transfer Cart vMaxComp cumulDiff ") 
		<< std::setprecision(10) << std::fixed
		<< vMaxComp << " " << cumulDiff << std::endl;

	return cumulDiff;

}

/*!
 * <!-- Check coordinate transfer -->
*/
SimTK::Real Context::checkTransferCoordinates_BAT(int srcWIx, int destWIx, bool wantJacobian)
{

	// Get advanced states of the integrators
	SimTK::State& srcAdvancedState = worlds[srcWIx].integ->updAdvancedState();
	SimTK::State& destAdvancedState = worlds[destWIx].integ->updAdvancedState();

	// Get coordinates from source World
	const std::vector<std::vector<std::pair<
		bSpecificAtom *, SimTK::Vec3> > >&
		srcWorldsAtomsLocations =
	worlds[srcWIx].getAtomsLocationsInGround(srcAdvancedState);

	// Get coordinates from destintaion World
	const std::vector<std::vector<std::pair<
		bSpecificAtom *, SimTK::Vec3> > >&
		destWorldsAtomsLocations =
	worlds[destWIx].getAtomsLocationsInGround(destAdvancedState);

	// Iterate molecules
	int allCnt = 0;
	int topoIx = 0;

	// Get locations of this molecule
	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> srcAtomTargets;
	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> destAtomTargets;
	for (int i = 0; i < topologies.size(); i++) {
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> temp;
		worlds[srcWIx].extractAtomTargets(i, srcWorldsAtomsLocations, temp);
		srcAtomTargets.insert(temp.begin(), temp.end());
	}
	for (int i = 0; i < topologies.size(); i++) {
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> temp;
		worlds[destWIx].extractAtomTargets(i, destWorldsAtomsLocations, temp);
		destAtomTargets.insert(temp.begin(), temp.end());
	}

	int rowCnt = 0;
	SimTK::Real bondLength_src, bondBend_src, bondTorsion_src;
	SimTK::Real bondLength_des, bondBend_des, bondTorsion_des;
	SimTK::Real cosTheta_src, cosTheta_des;
	SimTK::Real vMaxComp = 0;
	SimTK::Real currDiff;
	SimTK::Real cumulDiff = 0;

	for (const auto& row : zMatrixTable) {
		
		bondLength_src = SimTK::NaN; bondLength_des = SimTK::NaN;
		cosTheta_src = SimTK::NaN; cosTheta_des = SimTK::NaN;
		bondBend_src = SimTK::NaN; bondBend_des = SimTK::NaN;
		bondTorsion_src = SimTK::NaN; bondTorsion_des = SimTK::NaN;

		SimTK::Compound::AtomIndex a0_cAIx, a1_cAIx;
		
		// Calculate source bond length
		a0_cAIx = atoms[row[0]].getCompoundAtomIndex();
		a1_cAIx = atoms[row[1]].getCompoundAtomIndex();
		
		SimTK::Vec3 a0loc = findAtomTarget(srcAtomTargets, a0_cAIx);
		SimTK::Vec3 a1loc = findAtomTarget(srcAtomTargets, a1_cAIx);

		SimTK::Vec3 v_a0a1 = a0loc - a1loc;
		bondLength_src = std::sqrt(SimTK::dot(v_a0a1, v_a0a1));

		if(row[2] >= 0){

			SimTK::Compound::AtomIndex a2_cAIx;
			a2_cAIx = atoms[row[2]].getCompoundAtomIndex();
			SimTK::Vec3 a2loc = findAtomTarget(srcAtomTargets, a2_cAIx);

			// Calculate source angle
			UnitVec3 v1(v_a0a1);
			UnitVec3 v2(a2loc - a1loc);

			cosTheta_src = SimTK::dot(v1, v2);
			assert(cosTheta_src < 1.1);
			assert(cosTheta_src > -1.1);
			if (cosTheta_src > 1.0) cosTheta_src = 1.0;
			if (cosTheta_src < -1.0) cosTheta_src = -1.0;
			bondBend_src = std::acos(cosTheta_src);

			if(row[3] >= 0){
				SimTK::Compound::AtomIndex
					a3_cAIx = atoms[row[3]].getCompoundAtomIndex();
				SimTK::Vec3 a3loc = findAtomTarget(srcAtomTargets, a3_cAIx);

				bondTorsion_src = bDihedral(a0loc, a1loc, a2loc, a3loc);
			}

		} // angle

		// Calculate dest bond length
		a0loc = findAtomTarget(destAtomTargets, a0_cAIx);
		a1loc = findAtomTarget(destAtomTargets, a1_cAIx);

		v_a0a1 = a0loc - a1loc;
		bondLength_des = std::sqrt(SimTK::dot(v_a0a1, v_a0a1));

		if(row[2] >= 0){

			SimTK::Compound::AtomIndex a2_cAIx;
			a2_cAIx = atoms[row[2]].getCompoundAtomIndex();
			SimTK::Vec3 a2loc = findAtomTarget(destAtomTargets, a2_cAIx);

			// Calculate source angle
			UnitVec3 v1(v_a0a1);
			UnitVec3 v2(a2loc - a1loc);

			cosTheta_des = SimTK::dot(v1, v2);
			assert(cosTheta_des < 1.1);
			assert(cosTheta_des > -1.1);
			if (cosTheta_des > 1.0) cosTheta_des = 1.0;
			if (cosTheta_des < -1.0) cosTheta_des = -1.0;
			bondBend_des = std::acos(cosTheta_des);

			if(row[3] >= 0){
				SimTK::Compound::AtomIndex
					a3_cAIx = atoms[row[3]].getCompoundAtomIndex();
				SimTK::Vec3 a3loc = findAtomTarget(destAtomTargets, a3_cAIx);

				bondTorsion_des = bDihedral(a0loc, a1loc, a2loc, a3loc);
			}

		} // angle

		if(wantJacobian){
			// Jacobians
			SimTK::Real sinTheta_src = std::sqrt(1.0 - (cosTheta_src*cosTheta_src));
			SimTK::Real sinTheta_des = std::sqrt(1.0 - (cosTheta_des*cosTheta_des));
			SimTK::Real cosPhi_src = std::cos(bondTorsion_src);
			SimTK::Real sinPhi_src = std::sqrt(1.0 - (cosPhi_src*cosPhi_src));
			SimTK::Real cosPhi_des = std::cos(bondTorsion_des);
			SimTK::Real sinPhi_des = std::sqrt(1.0 - (cosPhi_des*cosPhi_des));

			std::vector<std::vector<SimTK::Real>> Jacobian_src(3, std::vector<SimTK::Real>(3));
			std::vector<std::vector<SimTK::Real>> Jacobian_des(3, std::vector<SimTK::Real>(3));

			Jacobian_src[0][0] = sinPhi_src * cosTheta_src;
			Jacobian_src[0][1] = bondLength_src * cosPhi_src * cosTheta_src;
			Jacobian_src[0][2] = -bondLength_src * sinPhi_src * sinTheta_src;

			Jacobian_src[1][0] = sinPhi_src * sinTheta_src;
			Jacobian_src[1][1] = bondLength_src * cosPhi_src * sinTheta_src;
			Jacobian_src[1][2] = -bondLength_src * sinPhi_src * cosTheta_src;

			Jacobian_src[2][0] = cosPhi_src;
			Jacobian_src[2][1] = -bondLength_src * cosPhi_src;
			Jacobian_src[2][2] = 0;		

			Jacobian_des[0][0] = sinPhi_des * cosTheta_des;
			Jacobian_des[0][1] = bondLength_des * cosPhi_des * cosTheta_des;
			Jacobian_des[0][2] = -bondLength_des * sinPhi_des * sinTheta_des;

			Jacobian_des[1][0] = sinPhi_des * sinTheta_des;
			Jacobian_des[1][1] = bondLength_des * cosPhi_des * sinTheta_des;
			Jacobian_des[1][2] = -bondLength_des * sinPhi_des * cosTheta_des;

			Jacobian_des[2][0] = cosPhi_des;
			Jacobian_des[2][1] = -bondLength_des * cosPhi_des;
			Jacobian_des[2][2] = 0;	

			std::vector<SimTK::Real> cartV_src(3, 0.0);
			std::vector<SimTK::Real> BATV_src{bondLength_src, bondBend_src, bondTorsion_src};
			std::vector<SimTK::Real> cartV_des(3, 0.0);
			std::vector<SimTK::Real> BATV_des{bondLength_des, bondBend_des, bondTorsion_des};

			for (size_t m_I = 0; m_I < Jacobian_src.size(); ++m_I) {
				for (size_t m_J = 0; m_J < Jacobian_src[m_I].size(); ++m_J) {
					cartV_src[m_I] += Jacobian_src[m_I][m_J] * BATV_src[m_J];
				}
			}

			for (size_t m_I = 0; m_I < Jacobian_des.size(); ++m_I) {
				for (size_t m_J = 0; m_J < Jacobian_des[m_I].size(); ++m_J) {
					cartV_des[m_I] += Jacobian_des[m_I][m_J] * BATV_des[m_J];
				}
			}

			// scout("checkTransfer BAT") 
			// << " " << getZMatrixTableEntry(rowCnt, 0) << " " << getZMatrixTableEntry(rowCnt, 1)
			// << " " << getZMatrixTableEntry(rowCnt, 2) << " " << getZMatrixTableEntry(rowCnt, 3)
			// 	<< " " << bondLength_src << " " << bondBend_src << " " << bondTorsion_src
			// 	<< " " << bondLength_des << " " << bondBend_des << " " << bondTorsion_des
			// 	<< " " << bondLength_des - bondLength_src << " " << bondBend_des - bondBend_src << " " << bondTorsion_des - bondTorsion_src
			// 	<< " " << cartV_src[0] << " " << cartV_src[1] << " " << cartV_src[2]
			// 	<< " " << cartV_des[0] << " " << cartV_des[1] << " " << cartV_des[2];
			// ceol;

			currDiff = std::abs(cartV_src[0] - cartV_des[0]);
			if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
			cumulDiff += currDiff;
			if(vMaxComp < currDiff){vMaxComp = currDiff;}

			currDiff = std::abs(cartV_src[1] - cartV_des[1]);
			if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
			cumulDiff += currDiff;
			if(vMaxComp < currDiff){vMaxComp = currDiff;}

			currDiff = std::abs(cartV_src[2] - cartV_des[2]);
			if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
			cumulDiff += currDiff;
			if(vMaxComp < currDiff){vMaxComp = currDiff;}

		}else{
			currDiff = std::abs(bondLength_src - bondLength_des);
			if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
			cumulDiff += currDiff;
			if(vMaxComp < currDiff){vMaxComp = currDiff;}

			currDiff = std::abs(bondBend_src - bondBend_des);
			if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
			cumulDiff += currDiff;
			if(vMaxComp < currDiff){vMaxComp = currDiff;}

			currDiff = std::abs(bondTorsion_src - bondTorsion_des);
			if(SimTK::isNaN(currDiff)){currDiff = 0.0;}
			cumulDiff += currDiff;
			if(vMaxComp < currDiff){vMaxComp = currDiff;}

			// scout("checkTransfer BAT") 
			// << " " << getZMatrixTableEntry(rowCnt, 0) << " " << getZMatrixTableEntry(rowCnt, 1)
			// << " " << getZMatrixTableEntry(rowCnt, 2) << " " << getZMatrixTableEntry(rowCnt, 3)
			// 	<< " " << bondLength_src << " " << bondBend_src << " " << bondTorsion_src
			// 	<< " " << bondLength_des << " " << bondBend_des << " " << bondTorsion_des
			// 	<< " " << bondLength_des - bondLength_src << " " << bondBend_des - bondBend_src << " " << bondTorsion_des - bondTorsion_src
			// 	;
			// ceol;			
		}

		if(row[3] == -2){
			topoIx++;
		}

		rowCnt++;

	} // every zMatrix row

	scout("Check coords transfer BAT vMaxComp cumulDiff ")
		<< std::setprecision(10) << std::fixed
		<< vMaxComp << " " << cumulDiff << std::endl;

	return cumulDiff;

}


void Context::transferCoordinates_ReplicaToWorld(int replicaIx, int destWIx)
{
	SimTK::State& state = worlds[destWIx].integ->updAdvancedState();	
	state = setAtoms_SP_NEW(destWIx, state, replicas[replicaIx].getAtomsLocationsInGround());	
}


// SP_NEW_TRANSFER ============================================================

/*! <!-- -->
*/
SimTK::State& Context::setAtoms_CompoundsAndDuMM(
	int destWIx,
	SimTK::State& someState,
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>> &
		otherWorldsAtomsLocations)
{
	// Get destination world
	World& destWorld = worlds[destWIx];

	// Arrays of Transforms
	SimTK::Transform G_X_T;

	// Loop through molecules/topologies
	for(std::size_t topoIx = 0; topoIx < otherWorldsAtomsLocations.size();
	topoIx++)
	{

		// 0. VISUALIZER
		///////////////////////////////////////////////////////////
		// Set the decorator
		if (destWorld.visual == true) {
			destWorld.paraMolecularDecorator->setAtomTargets(
				otherWorldsAtomsLocations[topoIx]);
		}

		// 1. COMPOUND MATCHDEFAULT
		///////////////////////////////////////////////////////////

		// Convenient vars
		Topology& currTopology = (topologies)[topoIx];
		int currNAtoms = currTopology.getNumAtoms();

		// Convert input coordinates to Compound-friendly datatype
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
		destWorld.extractAtomTargets(
			topoIx, otherWorldsAtomsLocations, atomTargets);

		// Use Molmodel's Compound match functions to set the new conf
		G_X_T = destWorld.setAtoms_Compound_Match(
			topoIx, atomTargets);

		//matchDefaultConfiguration(topoIx, atomTargets);

		// 2.1 MORE COMPOUND FOR DUMM
		///////////////////////////////////////////////////////////

		// Get locations for DuMM
		SimTK::Vec3 locs[currNAtoms];

		// Set CompoundAtom frameInMobilizedFrame and get loc in mobod
		destWorld.setAtoms_Compound_FramesAndLocsInMobods(
			topoIx, atomTargets, locs);

		// 2.2 DUMM
		///////////////////////////////////////////////////////////

		// Set atoms' stations on body
		destWorld.setAtoms_SetDuMMStations(topoIx, locs);
	} // every Topology

	return someState;

}

/*! <!--  -->
*/
SimTK::State&
Context::setAtoms_XPF_XBM( int wIx){
	// Get world's Simbody
	SimTK::State& someState = worlds[wIx].integ->updAdvancedState();
	SimTK::SimbodyMatterSubsystem& matter = *worlds[wIx].matter;
	SimTK::DuMMForceFieldSubsystem &dumm = *worlds[wIx].forceField;

	// All the bonds from all Compounds
	const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// Iterate molecules
	size_t allBONDSCnt = 0; // counter for all bonds
	for(size_t topoIx = 0; topoIx < getNofMolecules(); topoIx++){

		// Get molecule and it's bonds
		Topology& topology = topologies[topoIx];
		const std::vector<BOND>& BONDS = allBONDS[topoIx];

		// Iterate molecule's bonds
		for(size_t BOIx = 0; BOIx < BONDS.size(); BOIx++){

			// Get current bond
			const BOND& currBOND = BONDS[BOIx];

			// Get bond's atoms
			int childNo = currBOND.first;
			int parentNo = currBOND.second;

			bSpecificAtom& childAtom  = atoms[currBOND.first];
			bSpecificAtom& parentAtom = atoms[currBOND.second];

			// int childTopoIx = childAtom.getMoleculeIndex();
			// int parentTopoIx = parentAtom.getMoleculeIndex();

			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();

			// SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
			// SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);

			// Get atoms' mobods
			SimTK::MobilizedBodyIndex childAtomMbx = topology.getAtomMobilizedBodyIndexThroughDumm(child_cAIx, dumm);
			SimTK::MobilizedBody& childAtomMobod = matter.updMobilizedBody(childAtomMbx);
			SimTK::MobilizedBodyIndex parentAtomMbx = topology.getAtomMobilizedBodyIndexThroughDumm(parent_cAIx, dumm);
			SimTK::MobilizedBody& parentAtomMobod = matter.updMobilizedBody(parentAtomMbx);

			const SimTK::MobilizedBody& parentMobod_const =  childAtomMobod.getParentMobilizedBody();
			SimTK::MobilizedBodyIndex parentMbx = parentMobod_const.getMobilizedBodyIndex();
			SimTK::MobilizedBody& parentMobod = matter.updMobilizedBody(parentMbx);
		
			// Set mobods transforms on flexible joints
			size_t boIx = BONDS_to_bonds[topoIx][BOIx];
			bBond& bond = bonds[boIx];
			if(bond.getBondMobility(wIx) != SimTK::BondMobility::Mobility::Rigid){
				
				std::vector<SimTK::Transform> mobodTs =
					calc_XPF_XBM(
						wIx, topology,
						child_cAIx, parent_cAIx,
						bond.getBondMobility(wIx),
						someState);
				
				childAtomMobod.setDefaultInboardFrame(mobodTs[0]);
				childAtomMobod.setDefaultOutboardFrame(mobodTs[1]);

			}

			// Set mobods X_PFs and X_BMs for atoms
			SimTK::Transform G_X_T = topology.getTopLevelTransform();
			
			if(atoms[parentNo].getIsRoot()){
				SimTK::Transform T_X_base = topology.getTopTransform_FromMap(parent_cAIx);
				SimTK::Transform G_X_base = G_X_T * T_X_base;
				parentAtomMobod.setDefaultInboardFrame(G_X_base);
				parentAtomMobod.setDefaultOutboardFrame(SimTK::Transform());
			}else if(atoms[childNo].getIsRoot()){
				SimTK::Transform T_X_base = topology.getTopTransform_FromMap(parent_cAIx);
				SimTK::Transform G_X_base = G_X_T * T_X_base;
				childAtomMobod.setDefaultInboardFrame(G_X_base);
				childAtomMobod.setDefaultOutboardFrame(SimTK::Transform());				
			}

			// Next bond
			allBONDSCnt++;

		} // every bond

	} // every molecule

	// ------------------------------------------------------------------------


	// Return
	return someState;

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
SimTK::Transform
Context::calc_XFM(
	int wIx,
	Topology& topology,	
	SimTK::Compound::AtomIndex& childAIx,
	SimTK::Compound::AtomIndex& parentAIx,
	SimTK::BondMobility::Mobility mobility,
	const SimTK::State& someState) const
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

	// New
	SimTK::Transform X_FMspherical = SimTK::Transform()
		* XXXinv
	;

	// Return
	if(mobility == SimTK::BondMobility::Mobility::Spherical){
		return X_FMspherical;
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

	// Get bonds
	const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// Iterate molecules
	int allCnt = 0;
	for(size_t topoIx = 0; topoIx < getNofMolecules(); topoIx++){

		// Get molecule and it's bonds
		Topology& topology = topologies[topoIx];
		const std::vector<BOND>& BONDS = allBONDS[topoIx];

		// Iterate molecule's bonds
		for(size_t BOIx = 0; BOIx < BONDS.size(); BOIx++){

			// Get current bond
			const BOND& currBOND = BONDS[BOIx];
			size_t boIx = BONDS_to_bonds[topoIx][BOIx];
			bBond& bond = bonds[boIx];

			// Get bond's atoms
			int childNo = currBOND.first;
			int parentNo = currBOND.second;
			//int gparentNo = BONDS[internCoords.findBondByFirst(topoIx, parentNo)].second;
			//int ggparentNo = BONDS[internCoords.findBondByFirst(topoIx, gparentNo)].second;			

			bSpecificAtom& childAtom  = atoms[currBOND.first];
			bSpecificAtom& parentAtom = atoms[currBOND.second];
			//bSpecificAtom& gparentAtom = atoms[gparentNo];
			//bSpecificAtom& ggparentAtom = atoms[ggparentNo];			

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
				SimTK::Transform X_FM =
				calc_XFM(wIx, topology,
					child_cAIx, parent_cAIx,
					bond.getBondMobility(wIx),
					someState);

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

			}

			// Next bond
			allCnt++;

		} // every bond

	} // every molecule

	return someState;

}

/*! <!-- __no_desk__. 
 * --> */
std::vector<SimTK::Transform>
Context::calc_XPF_XBM(
	int wIx, Topology& topology,
	SimTK::Compound::AtomIndex& childAIx,
	SimTK::Compound::AtomIndex& parentAIx,
	SimTK::BondMobility::Mobility mobility,
	const SimTK::State& someState)
{
	// Get world's forcefield and matter
	SimTK::DuMMForceFieldSubsystem &dumm = *worlds[wIx].forceField;
	SimTK::SimbodyMatterSubsystem& matter = *worlds[wIx].matter;

	// Get child atom's body
	SimTK::MobilizedBodyIndex childAtomMbx = topology.getAtomMobilizedBodyIndexThroughDumm(childAIx, dumm);
	const SimTK::MobilizedBody& childAtomMobod = matter.getMobilizedBody(childAtomMbx);

	// Get parent atom's body
	SimTK::MobilizedBodyIndex parentAtomMbx = topology.getAtomMobilizedBodyIndexThroughDumm(parentAIx, dumm);
	const SimTK::MobilizedBody& parentAtomMobod = matter.getMobilizedBody(parentAtomMbx);

	// Get parent body of the child atom's body
	const SimTK::MobilizedBody& parentMobod =  childAtomMobod.getParentMobilizedBody();
	SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

	// Bound to Ground
	if(parentMobod.isGround()){
		SimTK::Transform G_X_T = topology.getTopLevelTransform();
		SimTK::Transform T_X_base =
			topology.getTopTransform_FromMap(SimTK::Compound::AtomIndex(0));

		return std::vector<SimTK::Transform>
			{G_X_T * T_X_base,
			SimTK::Transform()};
	}

	// Get parent-child BondCenters relationship
	SimTK::Transform X_parentBC_childBC =
	  topology.getDefaultBondCenterFrameInOtherBondCenterFrame(
		childAIx, parentAIx);
	SimTK::Transform X_childBC_parentBC = ~X_parentBC_childBC;

	// Get Top frame
	SimTK::Transform T_X_root = topology.getTopTransform_FromMap(childAIx);

	// Get Top to parent frame

	const std::pair<int, SimTK::Compound::AtomIndex>& topoAtomPair = worlds[wIx].getMobodRootAtomIndex(parentMbx);
	SimTK::Compound::AtomIndex parentMobodAIx = topoAtomPair.second;

	//SimTK::Compound::AtomIndex parentRootAIx = worlds[wIx].getMbx2aIx()[parentMbx];
	SimTK::Compound::AtomIndex parentRootAIx = parentMobodAIx;
	
	// Origin of the parent mobod
	SimTK::Transform T_X_Proot = topology.getTopTransform_FromMap(parentRootAIx);
	SimTK::Transform Proot_X_T = ~T_X_Proot;
	SimTK::Transform Proot_X_root = Proot_X_T * T_X_root;
	
	// Get inboard dihedral angle
	SimTK::Angle inboardBondDihedralAngle =
		topology.bgetDefaultInboardDihedralAngle(childAIx);
	SimTK::Transform InboardDihedral_XAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::XAxis);
	SimTK::Transform InboardDihedral_ZAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::ZAxis);

	// Get inboard bond length
	SimTK::Real inboardBondlength = topology.bgetDefaultInboardBondLength(childAIx);
	SimTK::Transform InboardLength_mZAxis
		= SimTK::Transform(Rotation(), Vec3(0, 0, -inboardBondlength));

	// Samuel Flores' terminology
	SimTK::Transform M_X_pin =
		SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis);

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
	SimTK::Transform XXXinv = ~XXX;

	// Proot -> root -> parentBC -> chilBC=X -> Z
	SimTK::Transform P_X_F_spheric = SimTK::Transform()
		* Proot_X_root 
		* X_parentBC_childBC
		* X_to_Y
		* Y_to_Z
		* XXX
	;

	// Z -> X=childBC -> parentBC
	SimTK::Transform M_X_B_spheric = SimTK::Transform()
		* Z_to_Y
		* Y_to_X
		* X_childBC_parentBC
	;

	SimTK::Transform B_X_M_spheric = ~M_X_B_spheric;
	// SimTK::Transform B_X_M_spheric =
	// 	X_parentBC_childBC
	// 	* X_to_Z
	// 	* XXXinv
	// ;

	// ------------------------------------------------------------------------

	bool anglePin_OR = (
		   (mobility == SimTK::BondMobility::Mobility::AnglePin)
		|| (mobility == SimTK::BondMobility::Mobility::Slider)
		|| (mobility == SimTK::BondMobility::Mobility::BendStretch)
	);

	if( (anglePin_OR)
		//&& ((atom->neighbors).size() == 1)
	){
		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};

	}else if( (anglePin_OR)
		//&& ((atom->neighbors).size() != 1)
	){
		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};

	}else if(	(mobility == SimTK::BondMobility::Mobility::Torsion)
			||	(mobility == SimTK::BondMobility::Mobility::Cylinder)
	){
		return std::vector<SimTK::Transform> {P_X_F_pin, B_X_M_pin};

	}else if((mobility == SimTK::BondMobility::Mobility::BallM)
	|| (mobility == SimTK::BondMobility::Mobility::Rigid)
	|| (mobility == SimTK::BondMobility::Mobility::Translation) // Cartesian
	){
		return std::vector<SimTK::Transform> {P_X_F, B_X_M};

	// Spherical
	}else if((mobility == SimTK::BondMobility::Mobility::Spherical)
	){
		return std::vector<SimTK::Transform> {
			P_X_F_spheric,
			B_X_M_spheric
		};

	}else{
		warn("Warning: unknown mobility");
		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};
	}

}


/*! <!-- Transfer geometry to a world
 --> */
SimTK::State& Context::setAtoms_SP_NEW(
	int destWIx, SimTK::State& someState,
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>> & otherWorldsAtomsLocations)
{
	// Get destination world
	World& destWorld = worlds[destWIx];

	// Get dumm force field and matter
	SimTK::DuMMForceFieldSubsystem &dumm = *destWorld.forceField;
	SimTK::SimbodyMatterSubsystem& matter = *destWorld.matter;

	// Match Compound and DuMM coordinates
	someState = setAtoms_CompoundsAndDuMM(destWIx, someState, otherWorldsAtomsLocations);

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

	//std::cout << someState.getSystemStage() << std::endl;
	return someState;
	
}


/*!
 * <!--  -->
*/
SimTK::Compound::AtomIndex
Context::getChemicalParent_IfIAmRoot(
	int wIx,
	int atomNo,
	SimTK::DuMMForceFieldSubsystem &dumm
	)
{

	SimTK::Compound::AtomIndex aIx = atoms[atomNo].getCompoundAtomIndex();

	assert(!"Implementation not ready");
}

// SP_NEW_TRANSFER ------------------------------------------------------------

/*!
 * <!--  -->
*/
void Context::passThroughBonds_template(int whichWorld)
{

	//worlds[whichWorld].newFunction();

	const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// Iterate molecules
	int allCnt = 0;
	for(size_t topoIx = 0; topoIx < getNofMolecules(); topoIx++){

		// Get molecule and it's bonds
		Topology& topology = topologies[topoIx];
		const std::vector<BOND>& BONDS = allBONDS[topoIx];

		// Iterate molecule's bonds
		for(size_t BOIx = 0; BOIx < BONDS.size(); BOIx++){

			// Get current bond
			const BOND& currBOND = BONDS[BOIx];
			size_t boIx = BONDS_to_bonds[topoIx][BOIx];
			bBond& bond = bonds[boIx];

			// Get bond's atoms
			int childNo = currBOND.first;
			int parentNo = currBOND.second;
			//int gparentNo = BONDS[internCoords.findBondByFirst(topoIx, parentNo)].second;
			//int ggparentNo = BONDS[internCoords.findBondByFirst(topoIx, gparentNo)].second;			

			bSpecificAtom& childAtom  = atoms[currBOND.first];
			bSpecificAtom& parentAtom = atoms[currBOND.second];
			//bSpecificAtom& gparentAtom = atoms[gparentNo];
			//bSpecificAtom& ggparentAtom = atoms[ggparentNo];

			int childTopoIx = childAtom.getMoleculeIndex();
			int parentTopoIx = parentAtom.getMoleculeIndex();

			// Get Compound atom indexes
			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();
			//SimTK::Compound::AtomIndex gparent_cAIx = gparentAtom.getCompoundAtomIndex();
			//SimTK::Compound::AtomIndex ggparent_cAIx = ggparentAtom.getCompoundAtomIndex();			

			// Get DuMM atom indezes
			SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
			SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);

			scout("newFunc") <<" " << parentNo <<" " << childNo <<" "  << parent_dAIx <<" " << child_dAIx << eol; 	
		}
	}

}

// Go through all the worlds and generate samples
void Context::RunOneRound()
{

	for(std::size_t worldIx = 0; worldIx < getNofWorlds(); worldIx++){	

		// Rotate worlds indices (translate from right to left)
		if(isWorldsOrderRandom){
			randomizeWorldIndexes();
		}

		// Non-equilibrium parameters
		if( worlds[worldIx].getSampler(0)->getDistortOpt() < 0){

			// Determine the scale factor from the boost temperature
			SimTK::Real Tj = worlds[worldIx].getSampler(0)->getBoostTemperature();
			SimTK::Real Ti = worlds[worldIx].getSampler(0)->getTemperature();
			SimTK::Real qScaleFactor = std::sqrt(Tj / Ti);
			
			// Distribute the scale factor
			std::string how = "deterministic";
			bool randSignOpt = false;
			qScaleFactor = perturbScalingFactor(how, qScaleFactor, randSignOpt);

			// Send the scale factor to the world's samplers
			if(qScaleFactor != 1){
				setWorldDistortParameters(worldIx, qScaleFactor);
			}

		}


		// Run front world
		RunFrontWorldAndRotate(worldIndexes);

	} // END iteration through worlds

}

// Print to log and write pdbs
void Context::RunLog(int round)
{
	if (wantDCD) {
		writeDCDs();
	}

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

// Normal run
void Context::Run(int, SimTK::Real Ti, SimTK::Real Tf)
{

	// Initialize world indeces
	std::size_t currentWorldIx = worldIndexes.front();
	std::size_t lastWorldIx = 0;

	// Write an initial pdb
	for(int moli = 0; moli < this->nofMols; moli++){
		topologies[moli].writeAtomListPdb(outputDir,
			"/pdbs/ini." + std::to_string(moli) + "." + setupReader.get("SEED")[0] + ".", ".pdb",
			10, 0);
		std::cout << "XXXXXXXXXXXXXXXXXXXx" << std::endl;
	}

	// Write pdb at time 0
    writeInitialPdb();

	if( std::abs(Tf - Ti) < SimTK::TinyReal){ // Don't heat

		bool givenTsMode = false;

		// Main loop: iterate through rounds
		for(int round = 0; round < requiredNofRounds; round++){

			RunOneRound();

			RunLog(round);

			this->nofRounds++;

		}

	}else{// if Ti != Tf heating protocol
		SimTK::Real Tincr = (Tf - Ti) / static_cast<SimTK::Real>(requiredNofRounds);
		SimTK::Real currT = Ti;
		for(int round = 0; round < requiredNofRounds; round++){ // Iterate rounds

			// Set current temperature
			currT += Tincr;
			for (auto& world : worlds) {
				world.setTemperature(currT);
			}
			std::cout << "T= " << currT << std::endl;

			RunOneRound();

			RunLog(round);

			this->nofRounds++;
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
	for(unsigned int mol_i = 0; mol_i < nofMols; mol_i++){
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

// Write pdb
void Context::writePdb(std::size_t)
{
	// function args were std::size_t whichWorld
	assert(!"Not implemented"); throw std::exception();
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

/** Print the number of threads each World got **/
void Context::PrintNumThreads() {
	std::size_t worldIx = 0;
	for(auto& world : worlds) {
		std::cout << "World " << worldIx  << " requested "
			<< world.updForceField()->getNumThreadsRequested()
			<< " and got "
			<< world.updForceField()->getNumThreadsInUse()
			<< std::endl;

		worldIx++;
	}
}

SimTK::Real Context::getPotentialEnergy(std::size_t world, std::size_t sampler) const {
	return pHMC((worlds[world].samplers[sampler]))->pe_o;
}


// TESTING FUNCTIONS

/**
 *  Pass compounds to the new world
 */
void Context::areAllDuMMsTheSame(void)
{

	std::cout << "=== areAllDuMMsTheSame ===\n";

	// All molecules
	for(std::size_t molIx = 0; molIx < nofMols; molIx++){

		// All Compound atoms
		for(std::size_t k = 0; k < topologies[molIx].getNumAtoms(); k++){

			// Get atom index
			SimTK::Compound::AtomIndex aIx =
				(topologies[molIx].subAtomList[k]).getCompoundAtomIndex();

			std::cout << aIx << " ";

			// All worlds
			for(int worldIx = 0; worldIx < getNofWorlds(); worldIx++){

				// Get mobod
				SimTK::MobilizedBodyIndex mbx =
					topologies[molIx].getAtomMobilizedBodyIndexThroughDumm(aIx,
					*(worlds[worldIx].forceField) );

				SimTK::DuMM::AtomIndex dAIx = topologies[molIx].getDuMMAtomIndex(aIx);

				std::cout << dAIx << " ";	
			}
			std::cout << std::endl;

		}

		// TODO Restante DANGER
        //c.setTopLevelTransform(compoundTransform * c.getTopLevelTransform());
	}
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
	// Get bonds
	const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// Tag bonds that are checked 
	std::vector<int> taggedBonds; // TODO delete

	// Iterate molecules
	int allCnt = 0;
	for(size_t topoIx = 0; topoIx < getNofMolecules(); topoIx++){

		// Get molecule and it's bonds
		Topology& topology = topologies[topoIx]; // TODO delete
		const std::vector<BOND>& BONDS = allBONDS[topoIx];

		if (BONDS.size() == 0) {
			continue;
		}

		std::cout << "BONDS " << BONDS[0].first << " " << BONDS[0].second << std::endl;
		addZMatrixTableRow(std::vector<int> {
			BONDS[0].first,
			BONDS[0].second,
			-1, -2
		});

		addZMatrixTableRow(std::vector<int> {
			BONDS[1].first,
			BONDS[1].second,
			BONDS[internCoords.findBondByFirst(topoIx, BONDS[1].second)].second,
			-1
		});

		// Iterate molecule's bonds
		for(size_t BOIx = 2; BOIx < BONDS.size(); BOIx++){

			// Get current bond
			const BOND& currBOND = BONDS[BOIx];

			// Get bond's atoms
			int childNo = currBOND.first;
			int parentNo = currBOND.second;
			int gparentNo = BONDS[internCoords.findBondByFirst(topoIx, parentNo)].second;
			
			// Not all gparents have ggparents
			int ggparentNo = -3;
			int ggparentBONDIx = internCoords.findBondByFirst(topoIx, gparentNo);
			if(ggparentBONDIx >= 0){
				ggparentNo = BONDS[internCoords.findBondByFirst(topoIx, gparentNo)].second;
			}

			bSpecificAtom& childAtom = atoms[childNo]; // TODO delete
			bSpecificAtom& parentAtom = atoms[parentNo]; // TODO delete
			int childTopoIx = childAtom.getMoleculeIndex(); // TODO delete
			int parentTopoIx = parentAtom.getMoleculeIndex(); // TODO delete

			addZMatrixTableRow(std::vector<int> {childNo, parentNo, gparentNo, ggparentNo});

			// Next bond
			allCnt++;

		} // every bond

	} // every molecule	
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
		std::pair <bSpecificAtom *, SimTK::Vec3 > > >&
		otherWorldsAtomsLocations)
{

	// Iterate molecules
	int allCnt = 0;

	int topoIx = 0;

	// Get locations of this molecule
	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
	for (int i = 0; i < topologies.size(); i++) {
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> temp;
		worlds[wIx].extractAtomTargets(i, otherWorldsAtomsLocations, temp);
		atomTargets.insert(temp.begin(), temp.end());
	}

	int rowCnt = 0;
	SimTK::Real bondLength, bondBend, bondTorsion;
	for (const auto& row : zMatrixTable) {
		
		bondLength = SimTK::NaN;
		bondBend = SimTK::NaN;
		bondTorsion = SimTK::NaN;

		SimTK::Compound::AtomIndex a0_cAIx, a1_cAIx;
		
		// Calculate bond length
		a0_cAIx = atoms[row[0]].getCompoundAtomIndex();
		a1_cAIx = atoms[row[1]].getCompoundAtomIndex();
		SimTK::Vec3 a0loc = findAtomTarget(atomTargets, a0_cAIx);
		SimTK::Vec3 a1loc = findAtomTarget(atomTargets, a1_cAIx);

		SimTK::Vec3 v_a0a1 = a0loc - a1loc;
		SimTK::Real bondLength = std::sqrt(SimTK::dot(v_a0a1, v_a0a1));

		if(row[2] >= 0){

			SimTK::Compound::AtomIndex a2_cAIx;
			a2_cAIx = atoms[row[2]].getCompoundAtomIndex();
			SimTK::Vec3 a2loc = findAtomTarget(atomTargets, a2_cAIx);

			// Calculate angle
			UnitVec3 v1(v_a0a1);
			UnitVec3 v2(a2loc - a1loc);

			Real dotProduct = SimTK::dot(v1, v2);
			assert(dotProduct < 1.1);
			assert(dotProduct > -1.1);
			if (dotProduct > 1.0) dotProduct = 1.0;
			if (dotProduct < -1.0) dotProduct = -1.0;
			bondBend = std::acos(dotProduct);

			if(row[3] >= 0){
				SimTK::Compound::AtomIndex
					a3_cAIx = atoms[row[3]].getCompoundAtomIndex();
				SimTK::Vec3 a3loc = findAtomTarget(atomTargets, a3_cAIx);

				bondTorsion = bDihedral(a0loc, a1loc, a2loc, a3loc);
			}

		} // angle
		
		setZMatrixBATValue(rowCnt, 0, bondLength);
		setZMatrixBATValue(rowCnt, 1, bondBend);
		setZMatrixBATValue(rowCnt, 2, bondTorsion);

		if(row[3] == -2){
			topoIx++;
		}

		rowCnt++;

	} // every zMatrix row		

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

		for(size_t zCnt = 0; zCnt = zMatrixBAT.size(); zCnt++){

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
	
	// Get world
	World& world = worlds[wIx];

	// Get generalized coordinates
	const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// Iterate ZMatrix and bonds
	int prevMolIx = -1;
	size_t zMatCnt = 0;
	for (const auto& row : zMatrixTable) {

		std::vector<SimTK::Real>& BATrow = updZMatrixBATRow(zMatCnt);
		std::vector<SimTK::Real>& replicaBATrow = replicas[replicaIx].updZMatrixBATRow(zMatCnt);

		assert( (std::abs(BATrow[0] - replicaBATrow[0]) < 0.00001 || std::isnan(BATrow[0])) &&
				(std::abs(BATrow[1] - replicaBATrow[1]) < 0.00001 || std::isnan(BATrow[1])) &&
				(std::abs(BATrow[2] - replicaBATrow[2]) < 0.00001 || std::isnan(BATrow[2])) &&
				"BATrow and replicaBATrow not equal.");

		// scout("BATrow_ref") <<" ";
		// for (SimTK::Real BATvalue : BATrow) {
		// 	std::cout << std::setw(6) << BATvalue << " ";
		// }
		// for (SimTK::Real BATvalue : BATrow_ref) {
		// 	std::cout << std::setw(6) << BATvalue << " ";
		// }
		// ceol;

		// Get bond's atoms
		bSpecificAtom& childAtom  = atoms[row[0]];
		bSpecificAtom& parentAtom = atoms[row[1]];

		// Get molecule
		int childMolIx = childAtom.getMoleculeIndex();
		int parentMolIx = parentAtom.getMoleculeIndex();
		assert((childMolIx == parentMolIx) &&
			"Atoms from different molecules");
		Topology& topology = topologies[childMolIx];

		// Iterate BONDS and get bond // ======================================
		int BOIx;
		if(prevMolIx != childMolIx){
			BOIx = 0;
		}

		// Get bond
		const std::vector<BOND>& BONDS = allBONDS[childMolIx];
		const BOND& currBOND = BONDS[BOIx];
		size_t boIx = BONDS_to_bonds[childMolIx][BOIx];
		bBond& bond = bonds[boIx];

		// Insert reference to BAT entry into samplers if the bond is flexible
		//scout(" ") << MobilityStr [ bond.getBondMobility(wIx) ] <<" ";
		if(bond.getBondMobility(wIx) != SimTK::BondMobility::Rigid){

			// Get Molmodel indexes
			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
			SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);

			// Get mbx
			SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
			SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);

			// Insert key and value into the map
			for(size_t sami = 0; sami < worlds[wIx].samplers.size(); sami++){

				//(pHMC((worlds[wIx].samplers[sami]))->subZMatrixBATs_ref).insert({childMbx, BATrow});
				(pHMC((worlds[wIx].samplers[sami]))->subZMatrixBATs_ref).insert({childMbx, replicaBATrow});
			}

		} // is bond rigid

		BOIx++;

		if(prevMolIx != childMolIx){
			prevMolIx = childMolIx;
		} // every BOND -------------------------------------------------------

		zMatCnt++;
	}

	//world.samplers[0].variableBATs = worldBATs;
		
}

/*!
 * <!-- zmatrixbat_ Get BAT coordinates modifyable by a selected world -->
*/
void
Context::updSubZMatrixBATsToWorld(
	  int wIx
	, int replicaIx)
{
	
	// Get world
	World& world = worlds[wIx];

	// Get generalized coordinates
	const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// Iterate ZMatrix and bonds
	int prevMolIx = -1;
	size_t zMatCnt = 0;
	for (const auto& row : zMatrixTable) {

		std::vector<SimTK::Real>& replicaBATrow = replicas[replicaIx].updZMatrixBATRow(zMatCnt);

		// Get bond's atoms
		bSpecificAtom& childAtom  = atoms[row[0]];
		bSpecificAtom& parentAtom = atoms[row[1]];

		// Get molecule
		int childMolIx = childAtom.getMoleculeIndex();
		int parentMolIx = parentAtom.getMoleculeIndex();
		assert((childMolIx == parentMolIx) &&
			"Atoms from different molecules");
		Topology& topology = topologies[childMolIx];

		// Iterate BONDS and get bond // ======================================
		int BOIx;
		if(prevMolIx != childMolIx){
			BOIx = 0;
		}

		// Get bond
		const std::vector<BOND>& BONDS = allBONDS[childMolIx];
		const BOND& currBOND = BONDS[BOIx];
		size_t boIx = BONDS_to_bonds[childMolIx][BOIx];
		bBond& bond = bonds[boIx];

		//scout(" ") << MobilityStr [ bond.getBondMobility(wIx) ] <<" ";
		if(bond.getBondMobility(wIx) != SimTK::BondMobility::Rigid){

			// Get Molmodel indexes
			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
			SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);

			// Get mbx
			SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
			SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);

			// Insert key and value into the map
			for(size_t sami = 0; sami < worlds[wIx].samplers.size(); sami++){

				(pHMC((worlds[wIx].samplers[sami]))->subZMatrixBATs_ref).at(childMbx) = replicaBATrow;

			}

		}

		BOIx++;

		if(prevMolIx != childMolIx){
			prevMolIx = childMolIx;
		} // every BOND -------------------------------------------------------

		zMatCnt++;
	}
		
}


/*!
 * <!-- Go through the vector of worlds and if their equilibrium worlds run and 
 * rotate. -->
*/
void Context::updSubZMatrixBATsToAllWorlds(int replicaIx)
{
	// Get thermoState corresponding to this replica
	int thisThermoStateIx = replica2ThermoIxs[replicaIx];

	// Get this world indexes from the corresponding thermoState
	std::vector<int>& replicaWorldIxs = 
		thermodynamicStates[thisThermoStateIx].updWorldIndexes();

	// Get nof worlds in this replica
	size_t replicaNofWorlds = replicaWorldIxs.size();

	// Set BAT map for all replica's worlds
	int currFrontWIx = -1;
	for(std::size_t worldCnt = 0; worldCnt < replicaNofWorlds; worldCnt++){

		updSubZMatrixBATsToWorld( replicaWorldIxs[worldCnt], replicaIx);

	} // every world

}


/*!
 * <!-- Print BAT coordinates from a selected world -->
 */
void
Context::PrintWorldSubZMatrixBATs(
	int wIx)
{
	
	// Get world
	World& world = worlds[wIx];

	// Get generalized coordinates
	const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// Iterate ZMatrix and bonds
	int prevMolIx = -1;
	size_t zMatCnt = 0;
	for (const auto& row : zMatrixTable) {

		std::vector<SimTK::Real>& BATrow = updZMatrixBATRow(zMatCnt);

		// Get bond's atoms
		bSpecificAtom& childAtom  = atoms[row[0]];
		bSpecificAtom& parentAtom = atoms[row[1]];

		// Get molecule
		int childMolIx = childAtom.getMoleculeIndex();
		int parentMolIx = parentAtom.getMoleculeIndex();
		assert((childMolIx == parentMolIx) &&
			"Atoms from different molecules");
		Topology& topology = topologies[childMolIx];

		// Iterate BONDS and get bond // ======================================
		int BOIx;
		if(prevMolIx != childMolIx){
			BOIx = 0;
		}

		// Get bond
		const std::vector<BOND>& BONDS = allBONDS[childMolIx];
		const BOND& currBOND = BONDS[BOIx];
		size_t boIx = BONDS_to_bonds[childMolIx][BOIx];
		bBond& bond = bonds[boIx];

		//scout(" ") << MobilityStr [ bond.getBondMobility(wIx) ] <<" ";
		if(bond.getBondMobility(wIx) != SimTK::BondMobility::Rigid){

			// Get Molmodel indexes
			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
			SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);

			// Get mbx
			SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
			SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);

			// Insert key and value into the map
			for(size_t sami = 0; sami < worlds[wIx].samplers.size(); sami++){

				std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>&>&
					variableBATs = pHMC((worlds[wIx].samplers[sami]))->updSubZMatrixBATsRef();
					
				scout("WorldBAT ") << wIx <<" "; 
				for(auto varBAT : variableBATs.at(childMbx)){
					cout << varBAT <<" ";
				}
				ceol;
				

			}

		}

		BOIx++;

		if(prevMolIx != childMolIx){
			prevMolIx = childMolIx;
		} // every BOND -------------------------------------------------------

		zMatCnt++;
	}
		
}


/*!
 * <!-- zmatrixbat_ Relationship BAT - mobod transforms -->
*/
void Context::PrintZMatrixMobods(int wIx, SimTK::State& someState)
{

	// Get world
	World& world = worlds[wIx];

	// Get generalized coordinates
	SimTK::Vector qVector = someState.getQ();
	std::cout << "Q= " << qVector << eol;

	const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

	// Iterate ZMatrix and bonds
	int prevMolIx = -1;
	size_t zMatCnt = 0;
	for (const auto& row : zMatrixTable) {

		scout("ZMatrixBATEntry: ");

		// Print indexes
		for (int value : row) {
			std::cout << std::setw(6) << value <<" "; 
		}

		// Print BAT values
		const std::vector<SimTK::Real>& BATrow = getZMatrixBATRow(zMatCnt);
		for (SimTK::Real BATvalue : BATrow) {
			std::cout << std::setw(6) << BATvalue << " ";
		}

		// Get bond's atoms
		bSpecificAtom& childAtom  = atoms[row[0]];
		bSpecificAtom& parentAtom = atoms[row[1]];

		// Get molecule
		int childMolIx = childAtom.getMoleculeIndex();
		int parentMolIx = parentAtom.getMoleculeIndex();
		assert((childMolIx == parentMolIx) &&
			"Atoms from different molecules");
		Topology& topology = topologies[childMolIx];


		// Get current bond
		int BOIx;
		if(prevMolIx != childMolIx){
			BOIx = 0;
		}

		const std::vector<BOND>& BONDS = allBONDS[childMolIx];
		const BOND& currBOND = BONDS[BOIx];
		size_t boIx = BONDS_to_bonds[childMolIx][BOIx];
		bBond& bond = bonds[boIx];

		//scout(" ") << MobilityStr [ bond.getBondMobility(wIx) ] <<" ";
		if(bond.getBondMobility(wIx) != SimTK::BondMobility::Rigid){

			// Get Molmodel indexes
			SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
			SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();

			SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
			SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);

			// Get child-parent mobods
			SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
			
			SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);
			const SimTK::MobilizedBody &childMobod = world.matter->getMobilizedBody(childMbx);
			SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);
			const SimTK::MobilizedBody &parentMobod = world.matter->getMobilizedBody(parentMbx);

			childMobod.getFirstQIndex(someState);

			scout(" ") << childMbx <<" " << parentMbx <<" ";

			scout("| ")
				<< childMobod.getQAsVector(someState) <<" |" ;

			//scout("| ")
			//	<< parentMobod.getQAsVector(someState) <<" |";
			// if(mbx != parentMbx){
			// 	// Get default transforms
			// 	const SimTK::Transform& X_PF = mobod.getInboardFrame(someState);
			// 	const SimTK::Transform& X_BM = mobod.getOutboardFrame(someState);
			// 	const SimTK::Transform& X_FM = mobod.getMobilizerTransform(someState);
			// 	PrintTransform(X_PF, 6, "X_PF");
			// 	PrintTransform(X_BM, 6, "X_BM");
			// 	PrintTransform(X_FM, 6, "X_FM");
			// }else{
			// 	//ceol;
			// }
		}

		ceol;

		BOIx++;

		if(prevMolIx != childMolIx){
			prevMolIx = childMolIx;
		}

		zMatCnt++;
	}

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

		scout("Transformers table by atom: no dAIx topoIx cAIx mobods") << eol;
		size_t cnt = 0;
		for(auto &atom : atoms){

			size_t topoIx = atom.getMoleculeIndex();
			Topology& topology = topologies[topoIx];

			const SimTK::Compound::AtomIndex cAIx = atom.getCompoundAtomIndex();

			// Get dumm atom index (set in modelOneCompound Step 1)
			SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(cAIx);

			// Get vector of worlds mbxs which contain this atom
			std::vector<SimTK::MobilizedBodyIndex> worldsMbxs;

			size_t wIx = 0;
			for(auto &world : worlds){
				SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
				SimTK::MobilizedBodyIndex mbx = dumm.getAtomBody(dAIx);
				worldsMbxs.push_back(mbx);
				wIx++;
			}

			// What do we have so far
			scout("atomEntry: ") 
				<< cnt <<" "
				<< dAIx <<" "
				<< topoIx <<" "
				<< cAIx <<" ";

			scout(" ");	
			for(auto mbx: worldsMbxs){
				std::cout << mbx <<" ";
			}
			ceol;	

			// Increase atom number
			cnt++;
		} // every atom

		// --------------------------------------------------------------------

		scout("Transformers table by bond: allCnt topoIx BOIx boIx BOchild BOparent child_cAIx parent_cAIx child_dAIx parent_dAIx childMbx parentMbx flex flexStr") << eol;

		const std::vector<std::vector<BOND>> &allBONDS = internCoords.getBonds();

		assert(allBONDS.size() == getNofMolecules() 
			&& "internal coordinates nof molecules wrong");

		// Counter for all bonds
		size_t allCnt = 0;

		// Iterate molecules
		for(size_t topoIx = 0; topoIx < getNofMolecules(); topoIx++){


			// Get molecule and it's bonds
			Topology& topology = topologies[topoIx];
			const std::vector<BOND>& BONDS = allBONDS[topoIx];

			// Iterate molecule's bonds
			for(size_t BOIx = 0; BOIx < BONDS.size(); BOIx++){

				// Get current bond
				const BOND& currBOND = BONDS[BOIx];
				size_t boIx = BONDS_to_bonds[topoIx][BOIx];
				bBond& bond = bonds[boIx];

				// Get bond's atoms
				bSpecificAtom& childAtom  = atoms[currBOND.first];
				bSpecificAtom& parentAtom = atoms[currBOND.second];

				SimTK::Compound::AtomIndex child_cAIx = childAtom.getCompoundAtomIndex();
				SimTK::Compound::AtomIndex parent_cAIx = parentAtom.getCompoundAtomIndex();

				SimTK::DuMM::AtomIndex child_dAIx = topology.getDuMMAtomIndex(child_cAIx);
				SimTK::DuMM::AtomIndex parent_dAIx = topology.getDuMMAtomIndex(parent_cAIx);



				// Is it a base atom
				// const SimTK::Compound::SingleAtom &parentCompoundAtom = parentAtom.getSingleAtom();
				// SimTK::Compound::AtomIndex parCAIx = parentAtom.getCompoundAtomIndex();
				// const SimTK::Compound::AtomPathName parAtomPathName = parentCompoundAtom.getAtomName(parCAIx);
				if(parentAtom.getIsRoot() == true){

					scout("bondEntry: -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ");

					// Iterate worlds and get child-parent mobods
					size_t wIx = 0;
					for(auto &world : worlds){
						SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
						SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);
						const SimTK::MobilizedBody &parentMobod = world.matter->getMobilizedBody(parentMbx);
						const SimTK::MobilizedBody &grandMobod = parentMobod.getParentMobilizedBody();
						SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();

						std::cout << parentMbx <<" " << grandMbx <<" "
							<< getMobility(rootMobilities[wIx][topoIx]) <<" "
							<< rootMobilities[wIx][topoIx] <<" ";
						
						wIx++;
					} ceol;
			
				}


				// What we have so far
				scout("bondEntry: ") << allCnt <<" " << topoIx <<" " << BOIx <<" " << boIx <<" ";
				BONDS[BOIx].Print();

				scout(" ") << child_cAIx <<" " << parent_cAIx <<" " << child_dAIx <<" " << parent_dAIx <<" ";

				// Get vector of worlds mbxs which contain this atom
				std::vector<std::pair<
					SimTK::MobilizedBodyIndex, SimTK::MobilizedBodyIndex>> worldsMbxBonds;

				// Iterate worlds and get child-parent mobods
				size_t wIx = 0;
				for(auto &world : worlds){
					SimTK::DuMMForceFieldSubsystem& dumm = *(world.updForceField());
					SimTK::MobilizedBodyIndex childMbx = dumm.getAtomBody(child_dAIx);
					SimTK::MobilizedBodyIndex parentMbx = dumm.getAtomBody(parent_dAIx);

					worldsMbxBonds.push_back(std::pair<SimTK::MobilizedBodyIndex, SimTK::MobilizedBodyIndex>{childMbx, parentMbx});

					wIx++;
				}

				// Iterate worlds and print
				wIx = 0;	
				for(auto mbxPair: worldsMbxBonds){
					std::cout << mbxPair.first <<" " << mbxPair.second <<" " << bond.getBondMobility(wIx) <<" " << MobilityStr[ bond.getBondMobility(wIx) ] <<" ";
					wIx++;
				} ceol;

				allCnt++;

			} // every bond

		} // every molecule

}
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
 * <!-- Calculate Q statistics -->
*/
void Context::calcQStats(int thIx)
{

	// Get world indexes
	const std::vector<int> & thermoWorldIxs = thermodynamicStates[thIx].getWorldIndexes();

	// Iterate worlds
	for(const auto worldIx : thermoWorldIxs){

		// Get world's Qs
		SimTK::State& worldCurrentState = worlds[worldIx].integ->updAdvancedState();
		int NQ = (worlds[worldIx].getSimbodyMatterSubsystem())->getNQ(worldCurrentState);
		const SimTK::Vector & worldQs = (getWorld(worldIx).getSimbodyMatterSubsystem())->getQ(worldCurrentState);

		// Get Q statistics
		bool found = thermodynamicStates[thIx].calcQStats(worldIx, worlds[worldIx].getBMps(), worldQs, worlds[worldIx].getNofSamples());
		if(!found){
			warn("Context::calcQStats: World not " + std::to_string(worldIx) + " found. Q statistics not calculated...");
		}

	}


}

/*!
 * <!--  -->
*/
void Context::printQStats(int thIx)
{
	thermodynamicStates[thIx].printQStats();
}


//////////////////////////////////////////////////////
//-------------         Q Stats        ---------------
//////////////////////////////////////////////////////