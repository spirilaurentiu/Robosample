#include "Context.hpp"
#include "World.hpp"
#include "Sampler.hpp"
#include "readAmberInput.hpp"

#include <sys/stat.h>
#include <sys/sysinfo.h>

bool Context::initializeFromFile(const std::string &file)
{
	// Read input into a SetupReader object
	setupReader.ReadSetup(file);
	setupReader.dump(true);

	// Set the log filename
	std::string logFilename = CreateLogfilename(setupReader.get("OUTPUT_DIR")[0], std::stoi(setupReader.get("SEED")[0]));
	logFile = std::ofstream(logFilename);
	if ( !logFile.is_open() ) {
		std::cerr << cerr_prefix << "Failed to open log file " << logFilename << std::endl;
		return false;
	}

	// Set the directory where the logs and the trajectories are stored
	if ( !CreateOutputDirectory(setupReader.get("OUTPUT_DIR")[0]) ) {
		return false;
	}
	setOutputDir(setupReader.get("OUTPUT_DIR")[0]);

	if ( !CheckInputParameters(setupReader) ) {
		return false;
	}

	// Get molecules directory
	std::string molDir = GetMoleculeDirectoryShort(setupReader.get("MOLECULES")[0]);
	std::cout << "Molecule directory: " << molDir << std::endl << std::flush;
	setPdbPrefix(molDir + setupReader.get("SEED")[0]);

	// Alert user of CUDA environment variables
	if(SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT").empty()){
		std::cout << "CUDA_ROOT not set." << std::endl;
	}else{
		std::cout << "CUDA_ROOT set to " << SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT") << std::endl;
	}

	// Requested nof Worlds in input. We'll try to construct them all
	// but we're not sure if we'll going to succesed.
	const auto requestedNofMols = nofMols;
	std::cout << "Requested " << nofMols << " molecules in " << nofWorlds << " worlds" << std::endl;

	/////////// Add Worlds to context ////////////
	// Add Worlds to the  Every World instantiates a:
	// CompoundSystem, SimbodyMatterSubsystem, GeneralForceSubsystem,
	// DuMMForceSubsystem, Integrator, TimeStepper and optionally:
	// DecorationSubsystem, Visualizer, VisuzlizerReporter,
	//  ParaMolecularDecorator

	// Deal with visualizer before adding worlds.
	std::vector<double> visualizerFrequencies;
	int i = -1;
	for(auto ts : setupReader.get("TIMESTEPS")){
		i++;
		if (setupReader.get("VISUAL")[i] == "TRUE"){
			visualizerFrequencies.push_back(std::stod(ts));
		}else{
			visualizerFrequencies.push_back(0);
		}
	}

	// Add Worlds
	addEmptyWorlds(setupReader.get("WORLDS").size(), visualizerFrequencies);
	
	// Get how much available memory we have
	struct sysinfo info;
	if (sysinfo(&info) == 0) {
        std::cout << "Free memory: " << info.freeram / (1024 * 1024) << " MB" << std::endl;
    }

	// Request threads
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		setNumThreadsRequested(worldIx, std::stoi(setupReader.get("THREADS")[worldIx]));
	}
	PrintNumThreads();

	// Set force field scale factor.
	if(setupReader.get("FFSCALE")[0] == "AMBER"){
		setAmberForceFieldScaleFactors();
	}else{
		setGlobalForceFieldScaleFactor(std::stod(setupReader.get("FFSCALE")[0]));
	}

	// Set GBSA scale factor
	setGbsaGlobalScaleFactor(std::stod(setupReader.get("GBSA")[0]));

	// Use OpenMM if possible
	if(setupReader.get("OPENMM")[0] == "TRUE"){
		setUseOpenMMAcceleration(true);
	}

	if(setupReader.get("OPENMM_CalcOnlyNonbonded")[0] == "TRUE"){
		setUseOpenMMCalcOnlyNonBonded(true);
	} else {
		setUseOpenMMCalcOnlyNonBonded(false);
	}

	for(unsigned int worldIx = 0; worldIx < getNofWorlds(); worldIx++){
		// Only NoCutoff (0) and CutoffNonPeriodic(1) methods are supported. Additional 3 methods available in
		// OpenMM could be integrated as well.
		if(setupReader.get("NONBONDED_METHOD")[worldIx] == "1"){
			setNonbondedMethod(worldIx, 1);
			setNonbondedCutoff(worldIx, std::stod( setupReader.get("NONBONDED_CUTOFF")[worldIx] ));
		}

		if(setupReader.get("INTEGRATORS")[worldIx] == "OMMVV"){
			setUseOpenMMIntegration(worldIx,
				std::stod( setupReader.get("BOOST_TEMPERATURE")[worldIx]),
				std::stod(setupReader.get("TIMESTEPS")[worldIx]));
		}
		std::cout<<"SETTING INTEGRATOR in ROBOSAMPLE "<<std::endl << std::flush;
	}


    // Set Lennard-Jones mixing rule
	setVdwMixingRule(DuMMForceFieldSubsystem::LorentzBerthelot);

	// Add molecules based on the setup reader
	// amber -> robo
	AddMolecules(requestedNofMols, setupReader);

	//std::cout << "OS memory 2.\n" << exec("free") << std::endl;
	int finalNofMols = getNofMolecules();
	if(requestedNofMols != finalNofMols){
		std::cerr << cerr_prefix << "Something went wrong while adding the world" << std::endl;
		return false;
	}
	std::cout << "Added " << finalNofMols << " molecules" << std::endl;

	// Loads parameters into DuMM
	addDummParams(finalNofMols, setupReader);

	// Adopts compound by the CompoundSystem and loads maps of indexes
	model(finalNofMols, setupReader);

	// Allocate space for containers that keep statistics if we're doing any
	allocWorldsStatsContainers();

	// Adaptive Gibbs blocking
	setNofRoundsTillReblock(std::stoi((setupReader.get("ROUNDS_TILL_REBLOCK"))[0]));

	// Only now we can allocate memory for reblocking Q vectors
	//allocateReblockQsCacheQVectors();

	// Add Fixman torque (Additional ForceSubsystem) if required
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		if(setupReader.get("FIXMAN_TORQUE")[worldIx] == "TRUE"){
			addFixmanTorque(worldIx);
		}
	}

	// Is this necessary?
	realizeTopology();


	// Add membrane.
	bool haveMembrane = (setupReader.get("MEMBRANE")[0] != "ERROR_KEY_NOT_FOUND");
	if (haveMembrane){
		float memZWidth = std::stof(setupReader.get("MEMBRANE")[0]);
		addContactImplicitMembrane(memZWidth, setupReader);
	
	}

	// Add empty samplers to the worlds.
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		BaseSampler *p = addSampler(worldIx,
			to_upper(setupReader.get("SAMPLERS")[worldIx]),
			setupReader.get("INTEGRATORS")[worldIx]);
	}

	//std::cout << "OS memory 3\n" << exec("free") << std::endl;
	//////////////////////
	// Thermodynamics
	//////////////////////

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

		for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
			getWorld(worldIx)->setTopologyIXs(topologyIXs);
			getWorld(worldIx)->setAmberAtomIXs(amberAtomIXs);
			HMCSampler* sampler_p = pHMC(updWorld(worldIx)->updSampler(0));
			sampler_p->setSphereRadius(sphere_radius);
		}
	}

	// Set thermostats to the samplers
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {

		for (int samplerIx = 0;
		samplerIx < getWorld(worldIx)->getNofSamplers();
		samplerIx++) {
			HMCSampler* sampler_p = pHMC(updWorld(worldIx)->updSampler(samplerIx));
			sampler_p->setThermostat(
				setupReader.get("THERMOSTAT")[worldIx]);
		}

		setTemperature(worldIx, std::stof(setupReader.get("TEMPERATURE_INI")[worldIx]));
	}

	// Set the guidance Hamiltonian parameters
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		for (int samplerIx = 0;
		samplerIx < getWorld(worldIx)->getNofSamplers();
		samplerIx++) {
			HMCSampler* sampler_p = pHMC(updWorld(worldIx)->updSampler(samplerIx));
			sampler_p->setBoostTemperature(std::stof(setupReader.get("BOOST_TEMPERATURE")[worldIx]));
		}
	}

	//////////////////////
	// MD parameters
	//////////////////////
	// Set timesteps
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		for (int samplerIx = 0;
		samplerIx < getWorld(worldIx)->getNofSamplers();
		samplerIx++) {

			// Set timesteps
			setTimestep(worldIx, samplerIx, std::stof(setupReader.get("TIMESTEPS")[worldIx]));

			// Activate Fixman potential if needed
			if(setupReader.get("FIXMAN_POTENTIAL")[worldIx] == "TRUE"){
				useFixmanPotential(worldIx, samplerIx);
			}
		}
	}

	// Set the nunmber of MD steps
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for (int samplerIx = 0;
		samplerIx < getWorld(worldIx)->getNofSamplers();
		samplerIx++) {
			setNofMDStepsPerSample(worldIx,
				samplerIx,
				std::stoi(setupReader.get("MDSTEPS")[worldIx]));

			HMCSampler* sampler_p = pHMC(updWorld(worldIx)->updSampler(samplerIx));

			if(setupReader.find("MDSTEPS_STD")){
				sampler_p->setMDStepsPerSampleStd(std::stoi(setupReader.get("MDSTEPS_STD")[worldIx]));
			}else{
				sampler_p->setMDStepsPerSampleStd(0);
			}

		}
	}

	// Guidance Hamiltonian MD
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for (int samplerIx = 0; samplerIx < getWorld(worldIx)->getNofSamplers(); samplerIx++) {
			HMCSampler* sampler_p =
			pHMC(updWorld(worldIx)->updSampler(samplerIx));
			sampler_p->setBoostMDSteps(std::stoi(setupReader.get("BOOST_MDSTEPS")[worldIx]));
		}
	}

	//////////////////////
	// Non-equilibrium parameters
	//////////////////////

	// Q distortin parameters
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for (int samplerIx = 0; samplerIx < getWorld(worldIx)->getNofSamplers(); samplerIx++) {
			HMCSampler* sampler_p = pHMC(updWorld(worldIx)->updSampler(samplerIx));
			sampler_p->setDistortOption(std::stoi(setupReader.get("DISTORT_OPTION")[worldIx]));
		}
	}

	//////////////////////
	// Simulation parameters
	//////////////////////
	// Set the seeds for reproducibility. Samplers have to be here already.
	// Let the user set one seed only. TODO: better algorithm (Victor).
	if( setupReader.find("SEED") ){
		if( !(setupReader.get("SEED").empty()) ){
			for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
				setSeed(worldIx, 0, std::stoi(setupReader.get("SEED")[0]) + worldIx );
			}
		}
	}

	// Initialize samplers
	/** Set simulation temperature,
	velocities to desired temperature, variables that store the configuration
	and variables that store the energies, both needed for the
	acception-rejection step. Also realize velocities and initialize
	the timestepper. **/
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		
		for (int samplerIx = 0;
		samplerIx < getWorld(worldIx)->getNofSamplers();
		samplerIx++){

			initializeSampler(worldIx, samplerIx);
		
		}
	}

	// Set the number of samples per round
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		setNofSamplesPerRound(worldIx, std::stoi(setupReader.get("SAMPLES_PER_ROUND")[worldIx]));
	}

	// Simulation parameters
	int currentWorldIx = 0;
	int round_mcsteps = 0;

	setRequiredNofRounds(std::stoi(setupReader.get("ROUNDS")[0]));

	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		round_mcsteps += getNofSamplesPerRound(worldIx);
	}

	// Set pdb writing frequency
	setPdbRestartFreq( std::stoi(setupReader.get("WRITEPDBS")[0]) );

	// Get atom indeces for geometry calculations
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

	// // Write intial pdb for reference
	// if(setupReader.get("WRITEPDBS")[0] == "TRUE"){
	// 	writeInitialPdb();
	// }

	// Get output printing frequency
	setPrintFreq( std::stoi(setupReader.get("PRINT_FREQ")[0]) );

	// Realize topology for all the Worlds
	realizeTopology();

	// U Scale Factors uses maps stored in Topology
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		(updWorld(worldIx))->setUScaleFactorsToMobods();
	}

	// Realize topology for all the Worlds
	realizeTopology();

	//std::cout << "OS memory 4.\n" << exec("free") << std::endl;
	// Check atom stations for debug purposes
	checkAtomStationsThroughDumm();

	// Load/store Mobilized bodies joint types in samplers
	loadMbxsToMobilities();


	// Setup task spaces
	addTaskSpacesLS();

	// -- Setup REX --
	std::string runType = setupReader.get("RUN_TYPE")[0];
	if((runType[0] == 'R') || (runType[1] == 'E')){

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

		// Add replicas
		for(int k = 0; k < nofReplicas; k++){
			addReplica(k);
		}

		// Consider renaming
		loadReplica2ThermoIxs();

		PrintReplicas();

		// How many Gibbs rounds until replica swpas occurs
		setSwapEvery(std::stoi(setupReader.get("REX_SWAP_EVERY")[0]));

		setSwapFixman(std::stoi(setupReader.get("REX_SWAP_FIXMAN")[0]));

	}

	PrepareNonEquilibriumParams_Q();
		
	//std::cout << "OS memory 5.\n" << exec("free") << std::endl;
	// -- Run --
	
	setThermostatesNonequilibrium();


	/* // Add constraints
	//context.addConstraints();

	PrintInitialRecommendedTimesteps();

	if(setupReader.get("RUN_TYPE")[0] == "Normal"){
		setRunType(0);
		Run(getRequiredNofRounds(),
			std::stof(setupReader.get("TEMPERATURE_INI")[0]),
			std::stof(setupReader.get("TEMPERATURE_FIN")[0]));
			
	}else if(setupReader.get("RUN_TYPE")[0] == "SimulatedTempering") {
			setRunType(0);
			RunSimulatedTempering(getRequiredNofRounds(),
				std::stof(setupReader.get("TEMPERATURE_INI")[0]),
				std::stof(setupReader.get("TEMPERATURE_FIN")[0]));

	}else if(setupReader.get("RUN_TYPE")[0] == "REMC"){
		setRunType(1);
		RunREX();

	}else if(setupReader.get("RUN_TYPE")[0] == "RENEMC"){
		setRunType(2);
		RunREX();

	}else if(setupReader.get("RUN_TYPE")[0] == "RENE"){
		setRunType(3);
		RunREX();

	}else{
		setRunType(0);
		Run(getRequiredNofRounds(),
			std::stof(setupReader.get("TEMPERATURE_INI")[0]),
			std::stof(setupReader.get("TEMPERATURE_FIN")[0]));
	} */

    return true;
}

#include <stack>

std::vector<int> Context::findMolecules(const readAmberInput& reader) {
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

void Context::loadBonds(const readAmberInput& reader) {
	// Allocate memory for bonds list
	nbonds = reader.getNumberBonds();
	bonds.reserve(nbonds);

	// Iterate bonds and get atom indeces
	// This establishes a 1-to-1 correspondence between prmtop and Gmolmodel
	for(int i = 0; i < nbonds; i++) {
		bBond bond;
		bond.setIndex(i);
		bond.i = reader.getBondsAtomsIndex1(i);
		bond.j = reader.getBondsAtomsIndex2(i);
		bonds.push_back(bond);

		// BAD! this will break if we invalidate the bonds vector
		atoms[bond.i].addNeighbor(&atoms[bonds[i].j]);
		atoms[bond.i].addBond(&bonds[i]);

		atoms[bond.j].addNeighbor(&atoms[bonds[i].i]);
		atoms[bond.j].addBond(&bonds[i]);
	}

	// Assign the number of bonds an atom has and set the number of freebonds
	// equal to the number of bonds for now
	for(auto& atom : atoms) {
		atom.setNbonds(atom.bondsInvolved.size());
		atom.setFreebonds(atom.bondsInvolved.size());
	}
}

void Context::loadAtoms(const readAmberInput& reader) {
	// Alloc memory for atoms and bonds list
	natoms = reader.getNumberAtoms();
	atoms.reserve(natoms);

	// Iterate through atoms and set as much as possible from amberReader
	for(int i = 0; i < natoms; i++) {
		bSpecificAtom atom;

		// Assign an index like in prmtop
		atom.setNumber(i);

		// This is the name of the atom in the .prmtop file
		// Examples: "O1", "C1", "C2", "H1", "H10"
		const std::string initialName = reader.getAtomsName(i);
		
		// Set element
		const int atomicNumber = reader.getAtomicNumber(i);
 		atom.setAtomicNumber(atomicNumber);
		atom.setElem(elementCache.getSymbolByAtomicNumber(atomicNumber));

		// Assign a "unique" name. The generator is however limited.
		// Examples: AAAA, AAAB, AAAC, AAAD etc
		atom.generateName(i);

		// Store the initial name from prmtop
		// Examples: "O1", "C1", "C2", "H1", "H10"
		atom.setInName(initialName);

		// Set force field atom type
		// Examples: "O1", "C1", "C2", "H1", "H10"
		atom.setFfType(initialName);

		// Set charge as it is used in Amber
		constexpr SimTK::Real chargeMultiplier = 18.2223;
		atom.setCharge(reader.getAtomsCharge(i) / chargeMultiplier);

		// Set coordinates in nm (AMBER uses Angstroms)
		atom.setX(reader.getAtomsXcoord(i) / 10.0);
		atom.setY(reader.getAtomsYcoord(i) / 10.0);
		atom.setZ(reader.getAtomsZcoord(i) / 10.0);
		atom.setCartesians(
			reader.getAtomsXcoord(i) / 10.0,
			reader.getAtomsYcoord(i) / 10.0,
			reader.getAtomsZcoord(i) / 10.0 );

		// Set mass
		const int mass = reader.getAtomsMass(i);
		atom.setMass(mass);

		// Set Lennard-Jones parameters
		atom.setVdwRadius(reader.getAtomsRVdW(i));
		atom.setLJWellDepth(reader.getAtomsEpsilon(i));

		// Set residue name and index
		atom.setResidueName("UNK");
		atom.residueIndex = 1;

		// Add element to cache
		elementCache.addElement(atomicNumber, mass);

		atoms.push_back(atom);
	}
}

void Context::setAtomCompoundTypes() {
	// Angle TetrahedralAngle = 109.47 * Deg2Rad;

	// Set Gmolmodel name and element and inboard length
	for(auto& atom : atoms) {
		
		// bond centers must have been already set
		const std::string& currAtomName = atom.getName();
		const int atomicNumber = atom.getAtomicNumber();
		const int mass = atom.getMass();
		atom.setAtomCompoundType(elementCache.getElement(atomicNumber, mass));
		
		// // Add BondCenters
		// const int currAtomNBonds = atom.getNBonds();
		// if(currAtomNBonds > 0){
		// 	if (currAtomNBonds == 1){
		// 		atom.addFirstBondCenter("bond1", currAtomName);
		// 	} else {
		// 		atom.addFirstTwoBondCenters("bond1", "bond2", currAtomName, UnitVec3(1, 0, 0), UnitVec3(-0.5, 0.866025, 0.0));
		// 		if (currAtomNBonds > 2) {
		// 			atom.addLeftHandedBondCenter("bond3", currAtomName, TetrahedralAngle, TetrahedralAngle);
		// 		}
		// 		if (currAtomNBonds > 3){
		// 			atom.addRightHandedBondCenter("bond4", currAtomName, TetrahedralAngle, TetrahedralAngle);
		// 		}
		// 	}

		// 	// Set the inboard BondCenter
		// 	atom.setInboardBondCenter("bond1");
		// 	atom.setDefaultInboardBondLength(0.19);
		// }
	}
}

void Context::addBiotypes() {
	// Iterate atoms and define Biotypes with their indeces and names
	for(auto& atom : atoms) {
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

void Context::loadAmberSystem(const std::string& prmtop, const std::string& inpcrd) {
	readAmberInput reader;
	reader.readAmberFiles(inpcrd, prmtop);

	// look at prmtop and create a mapping between atoms and molecule index. for each molecule:
	// - create the final atomlist via internalcoordinates. this is where we read from files
	// - add dumm parameters and classes now because we need molecule names
	// - reorder all lists to bat (batomlist, bonds, flexibilitie etc)
	// - pass to topology ????

	loadAtoms(reader);
	loadBonds(reader);
	setAtomCompoundTypes();
	addBiotypes();

	// addDummParams()

	// int num_atoms = 0;
	// const auto moleculeAtomIndices = findMolecules(reader);
	// for (const auto& indices : moleculeAtomIndices) {

	// 	std::cout << indices << std::endl;
	// }

	return;

	// // molecules are stored as a vector of atoms
	// // bAtomList[moleculeIndex][atomIndex] where atomIndex is in range [0, atomsPerMolecule[moleculeIndex])
	// std::vector<std::vector<bSpecificAtom>> bAtomList(numMols + 1);
	// for (int i = 0; i < numMols; i++) {
	// 	bAtomList[i].resize(atomsPerMolecule[i]);
	// }

	// std::vector<bSpecificAtom> bAtomList;
	// bAtomList.resize(natoms);

	// std::vector<bBond> bonds;
	// bonds.resize(nbonds);

	// for (int i = 0; i < natoms; i++) {
	// 	bSpecificAtom atom;

	// 	// Assign an index like in prmtop
	// 	atom.setNumber(i);

	// 	// N, H1, OX etc
	// 	auto atomType = reader.getAtomsType(i);
	// 	ROBOSAMPLE::trim(atomType);
	// 	atom.setInName(atomType);
	// 	atom.setFfType(atomType);

 	// 	atom.setAtomicNumber( reader.getAtomicNumber(i) );

	// 	// Assign a "unique" name. The generator is however limited.
	// 	// AAAA, AAAB, AAAC, AAAD, AAAF etc - why do we need this?
	// 	atom.setName(GetUniqueName(i));

	// 	std::cout << i << " " << atom.getAtomicNumber() << " " << atom.getInName() << " " << atom.getName() << std::endl;

	// 	// Set charge as it is used in Amber 						????????????????????????
	// 	SimTK::Real chargeMultiplier = 18.2223;
	// 	atom.setCharge(reader.getAtomsCharge(i) / chargeMultiplier);

	// 	// Set coordinates in nm
	// 	atom.setX(reader.getAtomsXcoord(i) / 10.0);
	// 	atom.setY(reader.getAtomsYcoord(i) / 10.0);
	// 	atom.setZ(reader.getAtomsZcoord(i) / 10.0);

	// 	// Set mass
	// 	atom.setMass(reader.getAtomsMass(i));

	// 	// Set Lennard-Jones parameters
	// 	atom.setVdwRadius(reader.getAtomsRVdW(i));
	// 	atom.setLJWellDepth(reader.getAtomsEpsilon(i));

	// 	// Set residue name and index
	// 	atom.residueName = std::string("UNK");
	// 	atom.residueIndex = 1;

	// 	switch(atom.getAtomicNumber()) {
	// 	case 1:
	// 		atom.setElem("H");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name, Element(1, "Hydrogen", "H", atom.getMass()) );
	// 		break;
	// 	case 6:
	// 		atom.setElem("C");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(6, "Carbon", "C", atom.getMass()));
	// 		break;
	// 	case 7:
	// 		atom.setElem("N");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name, Element(7, "Nitrogen", "N", atom.getMass()));
	// 		break;
	// 	case 8:
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name, Element(8, "Oxygen", "O", atom.getMass()));
	// 		break;
	// 	case 9:
	// 		atom.setElem("F");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name, Element(9, "Fluorine", "F", atom.getMass()));
	// 		break;
	// 	case 11:
	// 		atom.setElem("Na");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(11, "Sodium", "Na", atom.getMass()));
	// 		break;
	// 	case 12:
	// 		atom.setElem("Mg");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(12, "Magnesium", "Mg", atom.getMass()));
	// 		break;
	// 	case 15:
	// 		atom.setElem("P");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(15, "Phosphorus", "P", atom.getMass()));
	// 		break;
	// 	case 16:
	// 		atom.setElem("S");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(16, "Sulfur", "S", atom.getMass()));
	// 		break;
	// 	case 17:
	// 		atom.setElem("Cl");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(17, "Chlorine", "Cl", atom.getMass()));
	// 		break;
	// 	case 19:
	// 		atom.setElem("K");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(19, "Potassium", "K", atom.getMass()));
	// 		break;
	// 	case 20:
	// 		atom.setElem("Ca");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(20, "Calcium", "Ca", atom.getMass()));
	// 		break;
	// 	case 23:
	// 		atom.setElem("V");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(23, "Vanadium", "V", atom.getMass()));
	// 		break;
	// 	case 25:
	// 		atom.setElem("Mn");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(25, "Manganese", "Mn", atom.getMass()));
	// 		break;
	// 		atom.setElem("O");
	// 	case 26:
	// 		atom.setElem("Fe");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(26, "Iron", "Fe", atom.getMass()));
	// 		break;
	// 	case 30:
	// 		atom.setElem("Zn");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(30, "Zinc", "Zn", atom.getMass()));
	// 		break;
	// 	case 35:
	// 		atom.setElem("Br");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name,  Element(35, "Bromine", "Br", atom.getMass()));
	// 		break;
	// 	case 53:
	// 		atom.setElem("I");
	// 		atom.compoundSingleAtom = new
	// 			Compound::SingleAtom(atom.name, Element(53, "Iodine", "I", atom.getMass()));
	// 		break;
	// 	default:
	// 		std::cerr << "Atom with atomic number " << atom.getAtomicNumber() << " not available." << std::endl;
	// 		// return false;
	// 	}

	// 	// Add bond centers
	// 	Angle TetrahedralAngle = 109.47 * Deg2Rad;

	// 	// Add BondCenters
	// 	for(int i = 0; i < (natoms); i++) {
	// 		(atom.compoundSingleAtom)->setCompoundName("SingleAtom");

	// 		if(atom.nbonds > 0){
	// 			if(atom.nbonds == 1){ // One bond only
	// 				(atom.compoundSingleAtom)->addFirstBondCenter("bond1",
	// 					atom.name);
	// 			}else if(atom.nbonds > 1){ // multiple bonds
	// 				(atom.compoundSingleAtom)->addFirstTwoBondCenters("bond1", "bond2",
	// 					atom.name, UnitVec3(1, 0, 0), UnitVec3(-0.5, 0.866025, 0.0));
	// 				if(atom.nbonds > 2){
	// 					(atom.compoundSingleAtom)->addLeftHandedBondCenter("bond3",
	// 						atom.name, TetrahedralAngle, TetrahedralAngle);
	// 				}
	// 				if(atom.nbonds > 3){
	// 					(atom.compoundSingleAtom)->addRightHandedBondCenter("bond4",
	// 						atom.name, TetrahedralAngle, TetrahedralAngle);
	// 				}
	// 			}

	// 			// Set the inboard BondCenter
	// 			(atom.compoundSingleAtom)->setInboardBondCenter("bond1");
	// 			(atom.compoundSingleAtom)->setDefaultInboardBondLength(0.19);

	// 		} // bonded atoms iterator
	// 	} // atoms iterator

	// 	// // Add atom to the list
	// 	// int whichMol = whichMolecule[i];
	// 	// bAtomList[whichMol].push_back(atom);
	// 	bAtomList.push_back(atom);
	// }

	// // Iterate bonds and get atom indeces
	// // This establishes a 1-to-1 correspondence between prmtop and Gmolmodel
	// for(int i=0; i<nbonds; i++){
	// 	bBond bond;
	// 	bond.setIndex(i);
	// 	bond.i = reader.getBondsAtomsIndex1(i);
	// 	bond.j = reader.getBondsAtomsIndex2(i);
	// 	bonds.push_back(bond);
	// }

	// // Assign neighbors and bonds involved for each atom
	// // which translates into pushing bSpecificAtom * and bBond *
	// // into their apropriate vectors
	// for(int i = 0; i < nbonds; i++){
	// 	(bAtomList[ bonds[i].i  ]).addNeighbor( &(bAtomList[ bonds[i].j  ]) );
	// 	(bAtomList[ bonds[i].i  ]).addBond( &(bonds[i]) );

	// 	(bAtomList[ bonds[i].j  ]).addNeighbor( &(bAtomList[ bonds[i].i  ]) );
	// 	(bAtomList[ bonds[i].j  ]).addBond( &(bonds[i]) );
	// }

	// // Assign the number of bonds an atom has and set the number of freebonds
	// // equal to the number of bonds for now
	// for (int i = 0; i < natoms; i++) {
	// 	bAtomList[i].nbonds = bAtomList[i].bondsInvolved.size();
	// 	bAtomList[i].freebonds = bAtomList[i].bondsInvolved.size();
	// }
}

/** 
 * Main run function
*/
void Context::Run() {
	if(getRunType() == RunType::Default) {
		Run(getRequiredNofRounds(), tempIni, tempFin);

	}else if(  (getRunType() == RunType::REMC)
			|| (getRunType() == RunType::RENEMC)
			|| (getRunType() == RunType::RENE)){
		RunREX();

	}else{
		std::cout << "[WARNING] " << "Unknown run type. Running default.\n" ;
		Run(getRequiredNofRounds(), tempIni, tempFin);

	}
}

/** 
 * Main run function
*/
void Context::Run(int rounds) {
	if(getRunType() == RunType::Default) {
		Run(rounds, tempIni, tempFin);

	}else if(  (getRunType() == RunType::REMC)
			|| (getRunType() == RunType::RENEMC)
			|| (getRunType() == RunType::RENE)){

		RunREX();
	}else{
		std::cout << "[WARNING] " << "Unknown run type. Running default.\n" ;
		Run(rounds, tempIni, tempFin);

	}
}

bool Context::CreateOutputDirectory(const std::string& outDir)
{
    if( !SimTK::Pathname::fileExists(outDir + "/pdbs") ){
		const int err = mkdir((outDir + "/pdbs").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (err == -1){
			std::cerr << cerr_prefix << "Failed to create " << outDir + "/pdbs" << std::endl;
			return false;
		}
	}
    
	return true;
}

std::string Context::CreateLogfilename(const std::string& outDir, long long int seed) const
{
	return outDir + std::string("/log.") + std::to_string(seed);
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
	
	if ( std::stoi((setupReader.get("SEED"))[0]) < 0 ) {
		std::cerr << cerr_prefix << std::stoi((setupReader.get("SEED"))[0]) << " must be positive" << std::endl;
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
		const auto inpcrd = setupReader.get("MOLECULES")[molIx] + "/" + setupReader.get("INPCRD")[molIx] + ".rst7";

		if(!SimTK::Pathname::fileExists(prmtop)) {
			std::cerr << cerr_prefix << "Molecule " + std::to_string(molIx) + " prmtop not found at " + prmtop << std::endl;
			return false;
		}

		if(!SimTK::Pathname::fileExists(inpcrd)) {
			std::cerr << cerr_prefix << "Molecule " + std::to_string(molIx) + " inpcrd not found at " + inpcrd << std::endl;
			return false;
		}
	}

	// Topology specific paramters
	for(std::size_t worldIx = 0; worldIx < inpNofWorlds; worldIx++) {
		for(std::size_t molIx = 0; molIx < inpNofMols; molIx++) {

			if(!SimTK::Pathname::fileExists(setupReader.get("MOLECULES")[molIx] + std::string("/") + setupReader.get("RBFILE")[molIx]) ){
				std::cerr << cerr_prefix << "world " + std::to_string(worldIx) + " molecule " + std::to_string(molIx) + " rb not found" << std::endl;
				return false;
			}

			if(!SimTK::Pathname::fileExists(setupReader.get("MOLECULES")[molIx] + std::string("/") + setupReader.get("FLEXFILE")[molIx]) ){
				std::cerr << cerr_prefix << "world " + std::to_string(worldIx) + " molecule " + std::to_string(molIx) + " flex not found" << std::endl;
				return false;
			}

		}
	}

	// Normal modes options
	NDistortOpt.resize(inpNofWorlds, 0);
	if(setupReader.find("DISTORT_OPTION")){
		for(std::size_t worldIx = 0; worldIx < inpNofWorlds; worldIx++){
			NDistortOpt[worldIx] = std::stoi(setupReader.get("DISTORT_OPTION")[worldIx]);
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

	// Set the type of simulation and temperatures
	if ( setupReader.get("RUN_TYPE")[0] == "REMC" ) {
		setRunType(RunType::REMC);
	} else if ( setupReader.get("RUN_TYPE")[0] == "RENEMC" ) {
		setRunType(RunType::RENEMC);
	}else if ( setupReader.get("RUN_TYPE")[0] == "RENE" ) {
		setRunType(RunType::RENE);
	} else {
		tempIni = std::stof(setupReader.get("TEMPERATURE_INI")[0]);
		tempFin = std::stof(setupReader.get("TEMPERATURE_FIN")[0]);

		setRunType(RunType::Default);
	}
	std::cout << "[CheckInput] " << "Im doing " << static_cast<int>(runType) << std::endl;

	// If we got here we can set global variables
	// Reserve memory
	nofWorlds = inpNofWorlds;
	nofMols = inpNofMols;
	nofEmbeddedTopologies = inpNofEmbeddedTopologies;

	worlds.reserve(nofWorlds);
	worldIndexes.reserve(nofWorlds);

	return true;
}

void Context::reserveWorldsAndTopologies( int inpNofWorlds, int inpNofMols,
	int inpNofEmbeddedTopologies)
{
	nofWorlds = inpNofWorlds;
	nofMols = inpNofMols;
	nofEmbeddedTopologies = inpNofEmbeddedTopologies;

	worlds.reserve(nofWorlds);
	worldIndexes.reserve(nofWorlds);
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
bool Context::loadFlexibleBondsSpecs(std::size_t whichWorld, int, std::string flexSpecsFN)
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

void Context::setRegimen (std::size_t whichWorld, int, std::string regimen)
{
	// function args were std::size_t whichWorld, int whichMolecule, std::string regimen
	regimens[whichWorld].push_back(regimen);
}

// Add a number of empty worlds
// Each world initializes the following objects:
//  - CompoundSystem
//  - SimbodyMatterSubsystem, GeneralForceSubsystem, DecorationSubsystem,
//		Visualizer, Visualizer::Reporter, DuMMForceFieldSubsystem,
//  - Integrator with a TimeStepper on top
void Context::addEmptyWorlds(std::size_t argNofWorlds,
	std::vector<double> visualizerFrequencies)
{
	for(unsigned int worldIx = 0;
		worldIx < argNofWorlds;
		worldIx++){
		if(visualizerFrequencies[worldIx] > 0){
			addWorld(true, visualizerFrequencies[worldIx]);
		}else{
			addWorld(false);
		}
	}

	if(argNofWorlds != nofWorlds){
		std::cerr << cerr_prefix << "Something went wrong while adding the world\n";
		throw std::exception();
		std::exit(1);
	}

	std::cout << "Added " << nofWorlds << " empty worlds" << std::endl;

	//
	//allocateReblockQsCache();

}

// Add an empty world to the context
World * Context::addWorld(bool visual, SimTK::Real visualizerFrequency){

	// Increment worldIndexes
	worldIndexes.push_back(worldIndexes.size());

	// Call World constructor
	worlds.emplace_back(worldIndexes.back(), nofMols, visual,
		visualizerFrequency);

	// Add another row in the matrix of flexibility filenames
	rbSpecsFNs.push_back(std::vector<std::string>());
	// Add another row in the matrix of rigidity filenames
	flexSpecsFNs.push_back(std::vector<std::string>());
	// Add another row in the matrix of regimen filenames
	regimens.push_back(std::vector<std::string>());

	// Variable World specific parameters
	nofSamplesPerRound.push_back(1);
	nofMDStepsPerSample.push_back(1);
	timesteps.push_back(0.002); // ps
	nofBoostStairs.push_back(0);

	// Store the number of worlds
	nofWorlds = worlds.size();

	//
	return &worlds.back();
}

// Use a SetupReader Object to read worlds information from a file
void Context::LoadWorldsFromSetup(SetupReader&)
{
	// function args were SetupReader& setupReader
	assert(!"Not implemented.");
	throw std::exception();
}

void Context::modelOneEmbeddedTopology(int whichTopology,
int whichWorld,
std::string rootMobilizer)
{
		this->rootMobilities.push_back(rootMobilizer);

		(updWorld(whichWorld))->compoundSystem->modelOneCompound(
			SimTK::CompoundSystem::CompoundIndex(whichTopology),
			rootMobilizer);

		SimTK::DuMMForceFieldSubsystem& dumm = *((updWorld(whichWorld))->forceField);

		for(std::size_t k = 0; k < topologies[whichTopology].getNumAtoms(); k++){
			SimTK::Compound::AtomIndex aIx =
				(topologies[whichTopology].bAtomList[k]).getCompoundAtomIndex();
			SimTK::MobilizedBodyIndex mbx =
				topologies[whichTopology].getAtomMobilizedBodyIndex(aIx);
			//std::cout << "k aIx mbx " << k << " " << aIx << " " << mbx;

			SimTK::MobilizedBodyIndex mbxCheck =
				topologies[whichTopology].getAtomMobilizedBodyIndexThroughDumm(aIx,
				dumm);

			//std::cout << " mbxCheck " << mbxCheck ;
			//std::cout << std::endl << std::flush;

		}
}

/** Load molecules based on loaded filenames **/
void Context::AddMolecules(
	int requestedNofMols,
	SetupReader& setupReader
){

	topologies.reserve(requestedNofMols);
	moleculeCount = -1;

	std::vector<std::string> argRoots = setupReader.get("ROOTS");
	std::vector<std::string> argRootMobilities = setupReader.get("ROOT_MOBILITY");

	//std::vector<readAmberInput> amberReader(requestedNofMols);
	amberReader.resize(requestedNofMols);

	// Iterate through topology filenames vector
	//for(unsigned int molIx = 0; molIx < nofMols; molIx++){
	for(unsigned int molIx = 0; molIx < requestedNofMols; molIx++){

		// Add filenames to Context filenames vectors
		// This has to be called before Worlds constructors so that
		// reserve will be called for molecules and topologies
		//int nofMols = static_cast<int>(setupReader.get("MOLECULES").size());

		std::string topFN =
			setupReader.get("MOLECULES")[molIx] + std::string("/")
			+ setupReader.get("PRMTOP")[molIx];

		std::string crdFN =
			setupReader.get("MOLECULES")[molIx] + std::string("/")
			+ setupReader.get("INPCRD")[molIx] + ".rst7";

		// make amber reader create multiple molecule
		// each molecule then generate one topology
		// amber: vector<molecule info>

		loadTopologyFile( topFN );

		loadCoordinatesFile( crdFN );

		// Initialize an input reader
		//readAmberInput amberReader;
		amberReader[molIx].readAmberFiles(crdFNs[molIx], topFNs[molIx]);

		// Keep track of the number of molecules
		moleculeCount++; // Used for unique names of molecules // from world
		std::string moleculeName = "MOL" + std::to_string(moleculeCount); 
		roots.emplace_back(argRoots[molIx]); // from world
		//rootMobilities.emplace_back("Pin"); // TODO: move to setflexibilities
		topologies.emplace_back(Topology{moleculeName}); // TODO is this ok? 

		// Set Gmolmodel atoms properties from a reader: number, name, element
		// initial name, force field type, charge, coordinates, mass,
		// LJ parameters
		topologies[molIx].SetGmolAtomPropertiesFromReader(&amberReader[molIx]); 

		// Set bonds properties from reader: bond indeces, atom neighbours
		topologies[molIx].SetGmolBondingPropertiesFromReader(&amberReader[molIx]);

		// Set atoms Molmodel types (Compound::SingleAtom derived) based on
		// their valence // from world
		topologies[molIx].SetGmolAtomsCompoundTypes();

		// topologies[molIx].load(amberReader[molIx]);

		// Add Biotype indeces and Biotype names representing Biotypes
		topologies[molIx].bAddBiotypes(&amberReader[molIx]);

		// topologies[molIx].BAT();

		// Build Robosample graph and Compound graph.
		// It also asigns atom indexes in Compound
		// This is done only once and it needs
		topologies[molIx].buildGraphAndMatchCoords(
			std::stoi(roots.back()));

		// Helper function for calc MBAT determinant
		topologies[molIx].loadTriples();

		// Map of Compound atom indexes to Robosample atom indexes
		topologies[molIx].loadCompoundAtomIx2GmolAtomIx();
		//std::cout << "Topology " << molIx << " info\n";
		//topologies[molIx].printMaps();

	}


}

/* // Loads parameters into DuMM, adopts compound by the CompoundSystem
// and loads maps of indexes
void Context::addDummParams(
	int requestedNofMols,
	SetupReader& setupReader
){
	//std::vector<std::string> argRoots = setupReader.get("ROOTS");

	// Accumulate DuMM parameters in these vectors
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allBondsACIxs;
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allAnglesACIxs;
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allDihedralsACIxs;

	// Iterate through molecules
	for(unsigned int molIx = 0; molIx < requestedNofMols; molIx++){

		// Initialize an input reader
		readAmberInput amberReader;
		amberReader.readAmberFiles(crdFNs[molIx], topFNs[molIx]);
		
		// Pass current topology to the current world
		updWorld(0)->topologies = &topologies;

		// Add parameters in DuMM
		updWorld(0)->generateDummParams(molIx, &amberReader,
			allBondsACIxs, allAnglesACIxs, allDihedralsACIxs);

		for(unsigned int worldIx = 1; worldIx < nofWorlds; worldIx++){

			// Pass current topology to the current world
			(updWorld(worldIx))->topologies = &topologies;

			// Add parameters in DuMM
			(updWorld(worldIx))->transferDummParams(molIx, &amberReader);
		}

	}

} */

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
			(updWorld(worldIx))->forceField->defineAtomClass(aCIx, atomClassName.c_str(),
				atomParams.atomicNumber,
				atomParams.valence,
				atomParams.vdwRadius,
				atomParams.LJWellDepth
			);
		}
}

// Loads parameters into DuMM, adopts compound by the CompoundSystem
// and loads maps of indexes
void Context::addDummParams(
	int requestedNofMols,
	SetupReader& setupReader
){
	//std::vector<std::string> argRoots = setupReader.get("ROOTS");

	// Get an input reader
	std::vector<readAmberInput> amberReader(requestedNofMols);

	// Load Amber info for all molecules
	for(unsigned int molIx = 0; molIx < requestedNofMols; molIx++){
		(amberReader[molIx]).readAmberFiles(crdFNs[molIx], topFNs[molIx]);
	}

	// Accumulate DuMM parameters in these vectors
	std::map<AtomClassParams, AtomClassId> aClassParams2aClassId;
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allBondsACIxs;
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allAnglesACIxs;
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allDihedralsACIxs;
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allImpropersACIxs;


	// Load DuMM parameters for the first world
	for(unsigned int molIx = 0; molIx < requestedNofMols; molIx++){
		std::cout << "Context::addDummParams WORLD " << 0 << " topology " << molIx << std::endl << std::flush;

		// Pass current topology to the current world
		updWorld(0)->topologies = &topologies;

		// Add parameters in DuMM
		updWorld(0)->generateDummParams(molIx, &amberReader[molIx],
			aClassParams2aClassId,
			allBondsACIxs, allAnglesACIxs, allDihedralsACIxs, allImpropersACIxs);
	}

	// Load DuMM params for the rest of the worlds
	for(unsigned int worldIx = 1; worldIx < nofWorlds; worldIx++){

		// Accumulate DuMM parameters in these vectors
		//aClassParams2aClassId = std::map<AtomClassParams, AtomClassId>();
		allBondsACIxs = (std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>());
		allAnglesACIxs = (std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>());
		allDihedralsACIxs = (std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>());
		allImpropersACIxs = (std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>());


		updDummAtomClasses(aClassParams2aClassId, worldIx);

		// Iterate through molecules
		for(unsigned int molIx = 0; molIx < requestedNofMols; molIx++){
			std::cout << "Context::addDummParams WORLD " << worldIx << " topology " << molIx << std::endl << std::flush;

			// Pass current topology to the current world
			(updWorld(worldIx))->topologies = &topologies;

			// Add parameters in DuMM
			(updWorld(worldIx))->transferDummParams(molIx, &amberReader[molIx],
			aClassParams2aClassId,
			allBondsACIxs, allAnglesACIxs, allDihedralsACIxs, allImpropersACIxs);
		}

	}
	
}


// Loads parameters into DuMM, adopts compound by the CompoundSystem
// and loads maps of indexes
void Context::model(
	int requestedNofMols,
	SetupReader& setupReader
){
	std::vector<std::string> argRootMobilities = setupReader.get("ROOT_MOBILITY");

	// Iterate through molecules
	for(unsigned int molIx = 0; molIx < requestedNofMols; molIx++){

		for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){

			// TODO check if count rbfile == count flexfile
			loadRigidBodiesSpecs( worldIx, molIx,
				setupReader.get("MOLECULES")[molIx] + std::string("/")
				+ setupReader.get("RBFILE")[(requestedNofMols * worldIx) + molIx]
			);

			loadFlexibleBondsSpecs( worldIx, molIx,
				setupReader.get("MOLECULES")[molIx] + std::string("/")
				+ setupReader.get("FLEXFILE")[(requestedNofMols * worldIx) + molIx]
			);

			setRegimen( worldIx, molIx,
				setupReader.get("WORLDS")[worldIx] ); // TODO: delete from Topology

			std::cout << " Context::AddMolecule for world "<< worldIx << " " << std::endl;
			std::cout << " Context::AddMolecule molIx "<< molIx << " " << std::endl;
			std::cout << " Context::AddMolecule topFNs[molIx] "<< topFNs[molIx]
				<< " " << crdFNs[molIx] << " " << rbSpecsFNs[worldIx][molIx]
				<< std::endl << std::flush;

			//(updWorld(worldIx))->AllocateCoordBuffers(molIx); // TODO: remove

			// Set BondFlexibilities in Compound
			std::cout << "Context setting flexibility for mol "
				<< molIx << " world " << worldIx
				<< " regimenFN " << regimens[worldIx][molIx]
				<< " flexSpecsFNs " << flexSpecsFNs[worldIx][molIx]
				<< std::endl;

			topologies[molIx].setFlexibility(regimens[worldIx][molIx],
				flexSpecsFNs[worldIx][molIx], worldIx);

			// Set UScale factors. TODO: move in World
			topologies[molIx].setUScaleFactorsToBonds(flexSpecsFNs[worldIx][molIx]);

			// Add topology by CompoundSystem and add it to the
			//Visualizer's vector of molecules
			(updWorld(worldIx))->adoptTopology(molIx);

			// Calls modelOneCompound from CompoundSystem
			// amber (robo) -> dumm
			modelOneEmbeddedTopology(molIx, worldIx,
				argRootMobilities[(requestedNofMols * worldIx) + molIx]);

			// Realize Topology Stage involvs all the SubSystems
			//(updWorld(worldIx))->getCompoundSystem()->realizeTopology();

			topologies[molIx].loadAIx2MbxMap();
			(updWorld(worldIx))->loadMbx2AIxMap();
		}

	}

	// Realize topology for all the worlds all subsystems
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		(updWorld(worldIx))->getCompoundSystem()->realizeTopology();
	}

}

// Allocate space for containers that keep statistics if we're doing any
void Context::allocWorldsStatsContainers()
{

	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		worlds[worldIx].allocateStatsContainers();
	}

}


// Load/store Mobilized bodies joint types in samplers
void Context::loadMbxsToMobilities()
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for (int samplerIx = 0; samplerIx < getWorld(worldIx)->getNofSamplers(); samplerIx++){
			std::cout << "Loading mbx2mobility" << std::endl;

			// Pass compounds to the new world
			passTopologiesToNewWorld(worldIx);

			(updWorld(worldIx)->updSampler(samplerIx))->loadMbx2mobility(worldIx);
		}
	}
}

void Context::modelTopologies(std::vector<std::string> GroundToCompoundMobilizerTypes)
{

	// Model molecules
	std::cout << "Context::modelTopologies nof embedded Topologies "
		<< nofEmbeddedTopologies << "\n" << std::flush;

	for ( std::size_t molIx = 0; molIx < topologies.size(); molIx++){

		for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){

			std::cout << "Model molecule " << molIx
				<< " embedded in world " << worldIx << std::endl;

			this->rootMobilities.push_back(
				GroundToCompoundMobilizerTypes[(nofMols * worldIx) + molIx]);

			(updWorld(worldIx))->compoundSystem->modelOneCompound(
				SimTK::CompoundSystem::CompoundIndex(molIx),
				rootMobilities[(nofMols * worldIx) + molIx]);

			for(std::size_t k = 0; k < topologies[molIx].getNumAtoms(); k++){
				SimTK::Compound::AtomIndex aIx =
					(topologies[molIx].bAtomList[k]).getCompoundAtomIndex();
				SimTK::MobilizedBodyIndex mbx =
					topologies[molIx].getAtomMobilizedBodyIndex(aIx);
				//std::cout << "k aIx " << k << " " << aIx
				//	<< " " << mbx << std::endl << std::flush;
			}
		}

	}

// ZONE


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
		if ((updWorld(worldIx))->integ  == nullptr ){
			std::cout << "Context: integrator is null" << std::endl;
			break;
		}
		SimTK::VerletIntegrator& checkIntegrator = *(updWorld(worldIx))->integ;
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
		SimTK::CompoundSystem& compoundSystem = *((updWorld(worldIx))->getCompoundSystem());
		std::cout << "Context world " << worldIx << " compoundSystem nof compounds "
			<< compoundSystem.getNumCompounds() << std::endl;
		std::cout << "Context world " << worldIx << " System Topology realized "
			<< compoundSystem.getNumRealizationsOfThisStage(SimTK::Stage::Topology)
			<< " times.\n" << std::flush;

		// Matter
		////const SimTK::System& checkSystem = ((updWorld(worldIx))->matter)->getSystem();
		SimTK::SimbodyMatterSubsystem& matter = *((updWorld(worldIx))->matter);
		std::cout << "Context world " << worldIx
			<< " matter nofBodies " << matter.getNumBodies()
			<< " nofConstraints " << matter.getNumConstraints()
			<< "\n" << std::flush;

		// GeneralForceSubsystem
		SimTK::GeneralForceSubsystem& gfs = *((updWorld(worldIx))->forces);
		std::cout << "Context world " << worldIx
			<< " gfs nofForces " << gfs.getNumForces()
			<< "\n" << std::flush;

		SimTK::DuMMForceFieldSubsystem& dumm = *((updWorld(worldIx))->forceField);
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
			<< getWorld(worldIx)->getTemperature()
			<< std::endl;
		if(isUsingFixmanTorque(worldIx)){
			std::cout << "World " << worldIx
			<< " FixmanTorque temperature = "
			<< updWorld(worldIx)->updFixmanTorque()->getTemperature()
			<< std::endl;
		}
		for (int samplerIx = 0; samplerIx < getWorld(worldIx)->getNofSamplers(); samplerIx++){
			std::cout << "World " << worldIx << " Sampler " << samplerIx
				<< " temperature = " << updWorld(worldIx)->updSampler(samplerIx)->getTemperature()
				<< " initial const state PE: " << std::setprecision(20)
				//<< (updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((updWorld(worldIx))->integ->updAdvancedState())
				//<< (updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy(updAdvancedState(worldIx, samplerIx))
				<< " useFixmanPotential = "
				<< pHMC(updWorld(worldIx)->updSampler(samplerIx))->isUsingFixmanPotential()
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
			SimTK::DuMMForceFieldSubsystem& dumm = *((updWorld(worldIx))->forceField);
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
					= (topology.bAtomList[i]).getCompoundAtomIndex();
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(aIx);
				std::cout << "i aIx " << i << " " << aIx << " "
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
			samplerIx < getWorld(worldIx)->getNofSamplers();
			samplerIx++){
			(updWorld(worldIx)->updSampler(samplerIx))->checkAtomStationsThroughDumm();
		}
	}
}

// Get world
World * Context::getWorld() {
	return &worlds.back();
}

// Get a specific world
World * Context::getWorld(std::size_t which) {
	return &worlds[which];
}

// Get the last mutable world
World * Context::updWorld(){
	return &worlds.back();
}

// Get a mutable specific world
World * Context::updWorld(std::size_t which) {
	return &worlds[which];
}

std::size_t Context::getNofWorlds() const
{
	return nofWorlds;
}

SimTK::DuMMForceFieldSubsystem * Context::updForceField(std::size_t whichWorld)
{
	return worlds[whichWorld].updForceField();
}

int Context::getNofMolecules()
{
	return nofMols;
}

// Set mixing rule for Lennard-Jones
void Context::setVdwMixingRule(SimTK::DuMMForceFieldSubsystem::VdwMixingRule mixingRule){
	for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
		(updWorld(worldIx))->updForceField()->setVdwMixingRule(mixingRule);
	}
}

/////////////////////////
// --- Thermodynamics ---
/////////////////////////

// Get/set the main temperature (acc/rej temperature for MC)
SimTK::Real Context::getTemperature(std::size_t whichWorld) const
{
	return worlds[whichWorld].temperature;
}

void  Context::setTemperature(std::size_t whichWorld,
	float someTemperature)
{
	std::cout << " Context::setTemperature for world " 
		<< whichWorld << " " << someTemperature << std::endl;
	worlds[whichWorld].setTemperature(someTemperature);
}

// Set a temperature for all the worlds
void  Context::setTemperature(float someTemperature){
	for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
		worlds[worldIx].setTemperature(someTemperature);
	}
}


// If HMC, get/set the guidance Hamiltonian temperature
SimTK::Real Context::getGuidanceTemperature(std::size_t, std::size_t)
{
	// function args were std::size_t whichWorld, std::size_t whichSampler
	assert(!"Not implemented"); throw std::exception();

	return SimTK::NaN;
}

void Context::setGuidanceTemperature(std::size_t, std::size_t, SimTK::Real)
{
	// function args were std::size_t whichWorld,
	// std::size_t whichSampler, float someTemperature
	assert(!"Not implemented"); throw std::exception();
}
//------------

/////////////////////////
// --- Simulation parameters ---
/////////////////////////
/* 
* Add a sampler
*/
BaseSampler * Context::addSampler(
	std::size_t whichWorld,
	SamplerName whichSampler,
	IntegratorName whichIntegrator)
{
	assert(!"Not implemented");
/* 	// We only use HMCSampler for now. Later we'll add LAHMC and Girolami
	if( !samplerName.empty() ){

		// Add HMCSampler
		BaseSampler *p = worlds[whichWorld].addSampler(whichSampler);

		// Set the chain generation method (ex. Markov Cahin Monte Carlo)
		pHMC(p)->setSampleGenerator(whichSampler);

		// Set the integration method
		pHMC(p)->setIntegratorName(whichIntegrator);	

		return p;
	}else{
		// Replace with a macro
		std::cerr << cerr_prefix << "Context No sampler specified.\n";throw std::exception();std::exit(1);
	}	 */
}

/* 
* Add a sampler
*/
BaseSampler * Context::addSampler(
	std::size_t whichWorld,
	std::string samplerName,
	std::string integratorName)
{

	// We only use HMCSampler for now. Later we'll add LAHMC and Girolami
	if( !samplerName.empty() ){

		// Add HMCSampler
		BaseSampler *p = worlds[whichWorld].addSampler(SamplerName::HMC);
		
		// Set the chain generation method (ex. Markov Cahin Monte Carlo)
		pHMC(p)->setSampleGenerator(samplerName);

		// Set the integration method
		pHMC(p)->setIntegratorName(integratorName);

		return p;

	}else{
		// Replace with a macro
		std::cerr << "Context No sampler specified.\n";throw std::exception();std::exit(1);
	}

}

void Context::initializeSampler(std::size_t whichWorld,
	std::size_t whichSampler)
{

	World& thisWorld = worlds[whichWorld];
	SimTK::State& worldAdvancedState = thisWorld.integ->updAdvancedState();
	HMCSampler *poHMC = worlds[whichWorld].updSampler(whichSampler);
	poHMC->initialize( worldAdvancedState );

	// auto compoundSystem = worlds[whichWorld].getCompoundSystem();
	// auto forces = worlds[whichWorld].getGeneralForceSubsystem();
	// auto matter = worlds[whichWorld].getSimbodyMatterSubsystem();
	// worlds[whichWorld].updSampler(whichSampler)->initializeTaskSpace(*compoundSystem, *forces, *matter);
	// worlds[whichWorld].updSampler(whichSampler)->initializeTaskSpace(
	// *worlds[whichWorld].getCompoundSystem(),
	// *worlds[whichWorld].getGeneralForceSubsystem(),
	// *worlds[whichWorld].getSimbodyMatterSubsystem());
}


// Amber like scale factors.
void Context::setAmberForceFieldScaleFactors(std::size_t whichWorld)
{
	worlds[whichWorld].setAmberForceFieldScaleFactors();
}

// Amber like scale factors.
void Context::setAmberForceFieldScaleFactors()
{
	for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
		worlds[worldIx].setAmberForceFieldScaleFactors();
	}
}

// Set a global scaling factor for the forcefield
void Context::setGlobalForceFieldScaleFactor(
	std::size_t whichWorld, SimTK::Real globalScaleFactor)
{
	worlds[whichWorld].setGlobalForceFieldScaleFactor(globalScaleFactor);
}

void Context::setGlobalForceFieldScaleFactor(SimTK::Real globalScaleFactor)
{
	for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
		worlds[worldIx].setGlobalForceFieldScaleFactor(globalScaleFactor);
	}
}

// Set GBSA implicit solvent scale factor
void Context::setGbsaGlobalScaleFactor(std::size_t whichWorld, SimTK::Real gbsaGlobalScaleFactor)
{
	worlds[whichWorld].setGbsaGlobalScaleFactor(gbsaGlobalScaleFactor);
}

void Context::setGbsaGlobalScaleFactor(SimTK::Real gbsaGlobalScaleFactor){
	for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
		worlds[worldIx].setGbsaGlobalScaleFactor(gbsaGlobalScaleFactor);
	}
}

// If HMC, get/set the number of MD steps
int Context::getNofMDStepsPerSample(std::size_t whichWorld, std::size_t whichSampler)
{
   return pHMC(worlds[whichWorld].updSampler(whichSampler))->getMDStepsPerSample();
}

void Context::setNofMDStepsPerSample(
	std::size_t whichWorld,
	std::size_t whichSampler,
	int MDStepsPerSample)
{
   nofMDStepsPerSample[whichWorld] = MDStepsPerSample;

   pHMC(worlds[whichWorld].updSampler(whichSampler))->setMDStepsPerSample(MDStepsPerSample);

}

// If HMC, get/set timestep forMD
SimTK::Real Context::getTimestep(std::size_t whichWorld, std::size_t whichSampler) const
{
	return pHMC(worlds[whichWorld].getSampler(whichSampler))->getTimeStepper()->getIntegrator().getPredictedNextStepSize();
}

void Context::setTimestep(std::size_t whichWorld, std::size_t whichSampler, SimTK::Real argTimestep)
{
	//worlds[whichWorld].updSampler(whichSampler)->updTimeStepper()->updIntegrator().setFixedStepSize(argTimestep);
	pHMC(worlds[whichWorld].updSampler(whichSampler))->setTimestep(argTimestep);
}

// Use Fixman torque as an additional force subsystem
void Context::addFixmanTorque(std::size_t whichWorld)
{
	worlds[whichWorld].addFixmanTorque();
}

bool Context::isUsingFixmanTorque(std::size_t whichWorld) const
{
	return worlds[whichWorld].isUsingFixmanTorque();
}

void Context::setFixmanTorqueScaleFactor(std::size_t whichWorld, SimTK::Real scaleFactor)
{
	std::cout << "Context::setFixmanTorqueScaleFactor: ( (FixmanTorque *) (worlds["
	<< whichWorld << "]->updFixmanTorque()) )->setScaleFactor(" << scaleFactor << ") "<< std::endl;
	( (FixmanTorque *) (worlds[whichWorld].updFixmanTorque()) )->setScaleFactor(scaleFactor);
}

void Context::setFixmanTorqueTemperature(std::size_t whichWorld, SimTK::Real argTemperature)
{
	std::cout << "Context::setFixmanTemperature: ( (FixmanTorque *) (worlds["
	<< whichWorld << "]->updFixmanTorque()) )->setTemperature(" << argTemperature << ") "<< std::endl;
	( (FixmanTorque *) (worlds[whichWorld].updFixmanTorque()) )->setTemperature(argTemperature);
}

// Use Fixman potential
void Context::useFixmanPotential(std::size_t whichWorld, std::size_t whichSampler)
{
	pHMC(worlds[whichWorld].updSampler(whichSampler))->useFixmanPotential();
}

bool Context::isUsingFixmanPotential(std::size_t whichWorld, std::size_t whichSampler)
{
	return pHMC(worlds[whichWorld].updSampler(whichSampler))->isUsingFixmanPotential();
}


//------------

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

// This seems wrong !!!
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

// Get the number of samples returned by the sampler in one round
SimTK::Real Context::getNofSamplesPerRound(std::size_t whichWorld)
{
	return nofSamplesPerRound[whichWorld];
}

// Set the number of samples returned by the sampler in one round
void Context::setNofSamplesPerRound(std::size_t whichWorld, SimTK::Real MCStepsPerRound)
{
	nofSamplesPerRound[whichWorld] = MCStepsPerRound;
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

// SImulate Tempering
void Context::setNofBoostStairs(std::size_t whichWorld, int howManyStairs)
{
	nofBoostStairs[whichWorld] = howManyStairs;
}

int Context::getNofBoostStairs(std::size_t whichWorld)
{
	return nofBoostStairs[whichWorld];
}

// Simulated Tempering
void Context::RunSimulatedTempering(int, SimTK::Real, SimTK::Real) {

	std::size_t currentWorldIx = worldIndexes.front();
	std::size_t lastWorldIx = 0;

	// Write the initial Default Configuration of the first Compound of the first World
	PdbStructure  pdb(worlds[0].getTopology(0));
	std::ostringstream sstream;
	sstream << "pdbs/sb_" << (updWorld(worldIndexes.back())->getTopology(0)).getName() <<"_ini"<<".pdb";
	std::string ofilename = sstream.str();
	std::filebuf fb;
	std::cout<<"Writing pdb file: "<<ofilename<<std::endl;
	fb.open(ofilename.c_str(), std::ios::out);
	std::ostream os(&fb);
	pdb.write(os); // automatically multiplies by ten (nm to A)
	fb.close();

	// Simulated Tempering specifics
	//SimTK::Real iniBoostBeta = updWorld(worldIndexes.front())->updSampler(0)->getBeta(); // Intial beta
	//SimTK::Real finBoostBeta = 1.0 / (updWorld(worldIndexes.front())->updSampler(0)->getBoostTemperature() * SimTK_BOLTZMANN_CONSTANT_MD);
	//SimTK::Real iniBoostT = updWorld(worldIndexes.front())->updSampler(0)->getTemperature();
	//SimTK::Real finBoostT = updWorld(worldIndexes.front())->updSampler(0)->getBoostTemperature();
	//SimTK::Real dBoostT = (finBoostT - iniBoostT) / this->nofBoostStairs[0];

	// Main
	for(int round = 0; round < requiredNofRounds; round++){ // Iterate rounds
		for(std::size_t worldIx = 0; worldIx < getNofWorlds(); worldIx++){ // Iterate worlds

			// Rotate worlds indeces (translate from right to left)
			std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());

			// Get indeces
			currentWorldIx = worldIndexes.front();
			lastWorldIx = worldIndexes.back();

			// Transfer coordinates from last world to current
			SimTK::State& lastAdvancedState = updWorld(lastWorldIx)->integ->updAdvancedState();
			SimTK::State& currentAdvancedState = updWorld(currentWorldIx)->integ->updAdvancedState();

			if(worldIndexes.size() > 1) {
				// DANGER ZONE
				const std::vector<std::vector<std::pair<
					bSpecificAtom *, SimTK::Vec3> > >&
					otherWorldsAtomsLocations = updWorld(worldIndexes.back())->getAtomsLocationsInGround(lastAdvancedState);

					// Pass compounds to the new world
					//passTopologiesToNewWorld(currentWorldIx);

					currentAdvancedState = updWorld(currentWorldIx)->setAtomsLocationsInGround(
							currentAdvancedState, otherWorldsAtomsLocations);
				// SAFE ZONE
				//	currentAdvancedState = updWorld(currentWorldIx)->setAtomsLocationsInGround(
				//			currentAdvancedState,
				//			updWorld(worldIndexes.back())->getAtomsLocationsInGround(lastAdvancedState));
				// ZONE

			}

			// Check if reconstructions is done correctly
			//double backSetE = pMC(updWorld(lastWorldIx)->updSampler(0))->getSetPE();
			//double backCalcE = updWorld(lastWorldIx)->forceField->CalcFullPotEnergyIncludingRigidBodies(lastAdvancedState);
			//double currOldE = pMC(updWorld(currentWorldIx)->updSampler(0))->getOldPE();
			//double currCalcE = updWorld(currentWorldIx)->forceField->CalcFullPotEnergyIncludingRigidBodies(currentAdvancedState);

			// Set old potential energy of the new world
			pHMC(updWorld(currentWorldIx)->updSampler(0))->setOldPE(
					pHMC(updWorld(lastWorldIx)->updSampler(0))->getSetPE() );

			// Reinitialize current sampler
			updWorld(currentWorldIx)->updSampler(0)->reinitialize(currentAdvancedState);

			// Update
			for(int k = 0; k < int(getNofSamplesPerRound(currentWorldIx)); k++){ 
				updWorld(currentWorldIx)->updSampler(0)->sample_iteration(currentAdvancedState);
			} // END for samples

		} // for i in worlds

		// Print energy and geometric feat	ures
		if( !(round % getPrintFreq()) ){
			/* PrintSamplerDataToLog(worldIndexes.back());
			PrintDistancesToLog(worldIndexes.back());
			PrintAnglesToLog(worldIndexes.back());
			PrintDihedralsQsToLog(worldIndexes.back());
			fprintf(logFile, "\n"); */
			PrintToLog(0, worldIndexes.back(), 0);
		}

		// Write pdb
		if( pdbRestartFreq != 0){
			if(((round) % pdbRestartFreq) == 0){

				const SimTK::State& pdbState = updWorld(worldIndexes.front())->integ->updAdvancedState();
				updWorld(worldIndexes.front())->updateAtomListsFromCompound(pdbState);

				for(int mol_i = 0; mol_i < getNofMolecules(); mol_i++){
					updWorld(worldIndexes.front())->getTopology(mol_i).writeAtomListPdb(getOutputDir(),
						"/pdbs/sb." +
						getPdbPrefix() + "." + std::to_string(mol_i) + ".",
						".pdb", 10, round);
				}
			}
		} // if write pdbs

		this->nofRounds++;
	} // for i in rounds
}

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
		topologies[molIx].setMultibodySystem(
			*((updWorld(newWorldIx))->compoundSystem) );

		// Reset mobilized body indeces in Compound
		for(std::size_t k = 0; k < topologies[molIx].getNumAtoms(); k++){

			// Get atom index
			SimTK::Compound::AtomIndex aIx =
				(topologies[molIx].bAtomList[k]).getCompoundAtomIndex();

			// Get mobod
			SimTK::MobilizedBodyIndex mbx =
				topologies[molIx].getAtomMobilizedBodyIndexThroughDumm(aIx,
				*(worlds[newWorldIx].forceField) );

			// SimTK::MobilizedBodyIndex mbx =
			// 	topologies[molIx].getAtomMobilizedBodyIndexFromMap(aIx,
			// 	newWorldIx);

			// Set the corresponding mobod
			topologies[molIx].setAtomMobilizedBodyIndex(aIx, mbx);

		}

		// TODO Restante DANGER
        //c.setTopLevelTransform(compoundTransform * c.getTopLevelTransform());
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

// Adds a replica to the vector of Replica objects and sets the coordinates
// of the replica's atomsLocations
void Context::addReplica(int index)
{
	// Add replica and the vector of worlds
	replicas.emplace_back(Replica(index
		//argWorldIndexes,
		//timestepsInThisReplica,
		//mdstepsInThisReplica
	));

	///* // EXPERIMENTAL
    std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>
		referenceAtomsLocationsFromFile;

    // Iterate through molecules
    for(unsigned int molIx = 0; molIx < nofMols; molIx++){

        // Coordinate file prefix
        //std::string crdPrefix = crdFNs[molIx].substr(0, crdFNs[molIx].find("."));
        std::string crdPrefix = crdFNs[molIx].substr(0, crdFNs[molIx].find_last_of('.'));
        //std::string crdPrefix = crdFNs[molIx];

		// Read files
		readAmberInput amberReader;
		std::string crdFN = crdPrefix + ".s" + std::to_string(index) + ".rst7";
		std::cout << "Context::addReplica: " << "loading " << crdFN << std::endl << std::flush;
		amberReader.readAmberFiles(crdFN,  topFNs[0]);

		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> currentTopologyInfo;
		currentTopologyInfo.reserve( topologies[molIx].getNAtoms() );

		// Add each atom's location from file
		int i = -1;
		for (auto& atom : topologies[molIx].bAtomList) {
            i++;

			SimTK::Vec3 location(
                amberReader.getAtomsXcoord(i) / 10.0,
                amberReader.getAtomsYcoord(i) / 10.0,
                amberReader.getAtomsZcoord(i) / 10.0
            );

			currentTopologyInfo.emplace_back(&atom, location);
		}

		// Add topology to reference
		referenceAtomsLocationsFromFile.emplace_back(currentTopologyInfo);

    }//*/

	// Set replicas coordinates
	std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>
		referenceAtomsLocations =
		worlds[0].getCurrentAtomsLocationsInGround();

	replicas.back().setAtomsLocationsInGround(referenceAtomsLocationsFromFile);

	replicas.back().set_WORK_AtomsLocationsInGround(referenceAtomsLocationsFromFile);

    /* // EXPERIMENTAL
	replicas.back().PrintCoordinates();
    std::cout << "Context::addReplica world[0] reference coordinates\n" << std::flush;
	for(auto molecule : referenceAtomsLocations){
        for (auto atomPair : molecule){
            std::cout << atomPair.first->atomIndex << " " << atomPair.second << std::endl;
        }
	}
    */

	// Done
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
			mdstepsInThisReplica
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
	for(size_t thermoState_k = 0;
	thermoState_k < nofThermodynamicStates;
	thermoState_k++){

		replica2ThermoIxs.insert(
			std::pair<int, int>
			(thermoState_k, thermoState_k));

	}

	for(size_t thermoState_k = 0;
	thermoState_k < nofThermodynamicStates;
	thermoState_k++){

		thermo2ReplicaIxs.insert(
			std::pair<int, int>
			(thermoState_k, thermoState_k));

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
			if(distOpt == -1){
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

// Get Fixman potential already calculated from replica
SimTK::Real Context::getFixman(int replica_i)
{
    return replicas[replica_i].getFixman();
}

// Calculate Fixman potential of replica J in replica I's back world. Ui(X_j)
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

		/* TEST if replica j buffer has the same coordinates as its back world
		int thermoState_j = replica2ThermoIxs[replica_j];
		int world_j_back = thermodynamicStates[thermoState_j].getWorldIndexes().back();
		int world_j_front = thermodynamicStates[thermoState_j].getWorldIndexes().front();
		const std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>
        X_j_back =  worlds[world_j_back].getCurrentAtomsLocationsInGround();
        PrintCoordinates(X_j_back);
        //PrintCoordinates(X_j);
        */

        // Pass compounds to the new world
        passTopologiesToNewWorld(world_i_back);

		// Transfer coordinates from j to i
        SimTK::State& state = (worlds[world_i_back].integ)->updAdvancedState();
        worlds[world_i_back].setAtomsLocationsInGround(state, X_j);

        /* TEST if Xj is the same as world i back
		const std::vector<std::vector<
            std::pair <bSpecificAtom *, SimTK::Vec3>>>&
        X_i_back =  worlds[world_i_back].getCurrentAtomsLocationsInGround();
        PrintCoordinates(X_j);
        PrintCoordinates(X_i_back);
        */

		// Calculate Fixman in replica i back world
		SimTK::Real Fixman = worlds[world_i_back].calcFixman();

		// Transfer buffer coordinates of replica i
		// back to back world
		//transferCoordinates(world_i_front, world_i_back);
		restoreReplicaCoordinatesToBackWorld(replica_i);

		/* TEST if replica i buffer is the same as replica i back world
		const std::vector<std::vector<
            std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		X_i = replicas[replica_i].getAtomsLocationsInGround();
        const std::vector<std::vector<
            std::pair <bSpecificAtom *, SimTK::Vec3>>>&
        X_i_back =  worlds[world_i_back].getCurrentAtomsLocationsInGround();
        PrintCoordinates(X_i);
        PrintCoordinates(X_i_back);
        */

		passTopologiesToNewWorld(world_i_front);
        //restoreReplicaCoordinatesToFrontWorld(replica_i);

		// Return
        //std::cout << "CALC_FIXMAN return " << Fixman << std::endl << std::flush;
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
        worlds[world_i_back].setAtomsLocationsInGround(state, X_j);

		// Calculate Fixman in replica i back world
		SimTK::Real Fixman = worlds[world_i_back].calcFixman();

		// Transfer buffer coordinates of replica i
		// back to back world
		//transferCoordinates(world_i_front, world_i_back);
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
	    //std::cout << "CALC FIXMAN" << replica_j << " replica " << replica_i << "\n" << std::flush;

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
        worlds[world_j_back].setAtomsLocationsInGround(state, X_i);

		// Calculate Fixman in replica i back world
		SimTK::Real Fixman = worlds[world_j_back].calcFixman();

		// Transfer buffer coordinates of replica i
		// back to back world
		//transferCoordinates(world_j_front, world_j_back);
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
}

void Context::swapPotentialEnergies(int replica_i, int replica_j)
{
	// Exchange potential energies (not necessary)
	SimTK::Real tempE = replicas[replica_i].getPotentialEnergy();
	replicas[replica_i].setPotentialEnergy(replicas[replica_j].getPotentialEnergy());
	replicas[replica_j].setPotentialEnergy(tempE);
}

// Attempt swap between replicas r_i and r_j
// Code inspired from OpenmmTools
// Chodera JD and Shirts MR. Replica exchange and expanded ensemble simulations
// as Gibbs multistate: Simple improvements for enhanced mixing. J. Chem. Phys.
//, 135:194110, 2011. DOI:10.1063/1.3660669
// replica_i and replica_j are variable
bool Context::attemptREXSwap(int replica_X, int replica_Y)
{
	bool returnValue = false;

	// Get replicas' thermodynamic states indexes
	int thermoState_C = replica2ThermoIxs[replica_X];
	int thermoState_H = replica2ThermoIxs[replica_Y];

	// Record this attempt
	nofAttemptedSwapsMatrix[thermoState_C][thermoState_H] += 1;
	nofAttemptedSwapsMatrix[thermoState_H][thermoState_C] += 1;

	// For useful functions
	auto genericSampler = updWorld(0)->updSampler(0);

	// ----------------------------------------------------------------
	// INITIAL PE x,y 0
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Replica i reduced potential in state i
	SimTK::Real eC_X0 = thermodynamicStates[thermoState_C].getBeta()
		* replicas[replica_X].getPotentialEnergy();

	// Replica j reduced potential in state j
	SimTK::Real eH_Y0 = thermodynamicStates[thermoState_H].getBeta()
		* replicas[replica_Y].getPotentialEnergy();

	// Replica i reduced potential in state j
	SimTK::Real eH_X0 = thermodynamicStates[thermoState_H].getBeta()
		* replicas[replica_X].getPotentialEnergy();

	// Replica j reduced potential in state i
	SimTK::Real eC_Y0 = thermodynamicStates[thermoState_C].getBeta()
		* replicas[replica_Y].getPotentialEnergy();

	// ----------------------------------------------------------------
	// LAST PE x,y tau
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Replica i reduced potential in state i
	SimTK::Real lC_Xtau = thermodynamicStates[thermoState_C].getBeta()
		* (replicas[replica_X].getPotentialEnergy() + replicas[replica_X].getTransferedEnergy());

	// Replica j reduced potential in state j
	SimTK::Real lH_Ytau = thermodynamicStates[thermoState_H].getBeta()
		* (replicas[replica_Y].getPotentialEnergy() + replicas[replica_Y].getTransferedEnergy());

	// Replica i reduced potential in state j
	SimTK::Real lH_Xtau = thermodynamicStates[thermoState_H].getBeta()
		* (replicas[replica_X].getPotentialEnergy() + replicas[replica_X].getTransferedEnergy());

	// Replica j reduced potential in state i
	SimTK::Real lC_Ytau = thermodynamicStates[thermoState_C].getBeta()
		* (replicas[replica_Y].getPotentialEnergy() + replicas[replica_Y].getTransferedEnergy());
	// ========================================================================

	// Include the Fixman term if indicated
	SimTK::Real Uii = 0, Ujj = 0, Uij = 0, Uji = 0;
	if (swapFixman){

        if (thermoState_C == 0){
			std::cout << "Swap between " << thermoState_C << " and "
				<< thermoState_H << " ";

            // Replica i reduced Fixman potential in state i
            Uii = thermodynamicStates[thermoState_C].getBeta()
                * calcFixman_IinJ(replica_X, replica_X);

            // Replica j reduced Fixman potential in state j
            Ujj = thermodynamicStates[thermoState_H].getBeta()
                * calcFixman_IinJ(replica_Y, replica_Y);

            // Replica i reduced Fixman potential in state j
            Uij = thermodynamicStates[thermoState_H].getBeta()
                * calcFixman_IinJ(replica_X, replica_Y);

            // Replica j reduced Fixman potential in state i
            Uji = thermodynamicStates[thermoState_C].getBeta()
                * calcFixman_IinJ(replica_Y, replica_X); 
        }else{
            Uii = Ujj = Uij = Uji = 0;
        }

        std::cout << "Uii Ujj Uij Uji " << Uii << " " << Ujj
            << " " << Uij << " " << Uji << std::endl;
	}

	// ----------------------------------------------------------------
	// LOGP ENERGY EQUILIBRIUM
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real ETerm_equil    = eH_X0 - eC_X0;
	ETerm_equil               += eC_Y0 - eH_Y0;
	ETerm_equil = -1.0 * ETerm_equil;

	// ----------------------------------------------------------------
	// LOGP ENERGY NON-EQUILIBRIUM
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real
	ETerm_nonequil  = lH_Xtau - lC_Xtau;
	ETerm_nonequil += lC_Ytau - lH_Ytau;
	ETerm_nonequil = -1.0 * ETerm_nonequil;

	// ----------------------------------------------------------------
	// LOGP WORK
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real Work_X = lH_Xtau - eC_X0 - replicas[replica_X].get_WORK_Jacobian();
	SimTK::Real Work_Y = lC_Ytau - eH_Y0 - replicas[replica_Y].get_WORK_Jacobian();
	SimTK::Real WTerm = -1.0 * (Work_X + Work_Y);

	// ----------------------------------------------------------------
	// CORRECTION TERM FOR NON-EQUILIBRIUM
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real correctionTerm = 1.0;
	SimTK::Real miu_C = qScaleFactorsMiu.at(thermoState_C);
	SimTK::Real miu_H = qScaleFactorsMiu.at(thermoState_H);
	SimTK::Real std_C = qScaleFactorsStd.at(thermoState_C);
	SimTK::Real std_H = qScaleFactorsStd.at(thermoState_H);
	
	SimTK::Real s_X = qScaleFactors.at(thermoState_C);
	SimTK::Real s_Y = qScaleFactors.at(thermoState_H);
	SimTK::Real s_X_1 = 1.0 / s_X;
	SimTK::Real s_Y_1 = 1.0 / s_Y;

	SimTK::Real qC_s_X = 1.0, qH_s_Y = 1.0, qH_s_X_1 = 1.0, qC_s_Y_1 = 1.0;

	correctionTerm = (qH_s_X_1 * qC_s_Y_1) / (qC_s_X * qH_s_Y);

	// ----------------------------------------------------------------
	// PRINT
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	bool printTerms = false, printWithoutText = true;
	if (printTerms){
		std::cout << "thermoIxs " << thermoState_C << " " << thermoState_H << std::endl;
		std::cout << "replicaIxs " << replica_X << " " << replica_Y << std::endl;
		std::cout << "bibjwiwj "
			<< thermodynamicStates[thermoState_C].getBeta() << " "
			<< thermodynamicStates[thermoState_H].getBeta() << " "
			<< std::endl;
		std::cout << "LiiLjj " << lC_Xtau << " " << lH_Ytau << " "
			<< lH_Xtau << " " << lC_Ytau << std::endl;
		std::cout << "EiiEjj " << eC_X0 << " " << eH_Y0 << " "
			<< eH_X0 << " " << eC_Y0 << std::endl;
		std::cout << "Transferred E i j " << replicas[replica_X].getTransferedEnergy()
			<< " " << replicas[replica_Y].getTransferedEnergy() << std::endl;
		std::cout << "ETerm " << ETerm_equil << std::endl;
		std::cout << "WTerm " << WTerm << std::endl;
		std::cout << "correctionTerm s_i s_f " << correctionTerm 
			<< " " << s_X << " " << s_Y << " " << s_X_1 << " " << s_Y_1
			<< " " << qC_s_X << " " << qH_s_Y << " " << qH_s_X_1 << " " << qC_s_Y_1
			<< std::endl;
	}
	if(printWithoutText){
		std::cout << "REXdetails " << thermoState_C << " " << thermoState_H << " " 
			<< replica_X << " " << replica_Y << " " 
			<< thermodynamicStates[thermoState_C].getBeta() << " "
			<< thermodynamicStates[thermoState_H].getBeta() << " "
			<< eC_X0 << " " << eH_Y0 << " " << eH_X0 << " " << eC_Y0 << " "
			<< lC_Xtau << " " << lH_Ytau << " " << lH_Xtau << " " << lC_Ytau << " "
			<< replicas[replica_X].get_WORK_Jacobian() << " " 
			<< replicas[replica_Y].get_WORK_Jacobian() << " "
			<< replicas[replica_X].getTransferedEnergy() << " "
			<< replicas[replica_Y].getTransferedEnergy() << " "
			<< " " << s_X << " " << s_Y << " " << s_X_1 << " " << s_Y_1 << " "
			<< " " << qC_s_X << " " << qH_s_Y << " " << qH_s_X_1 << " " << qC_s_Y_1 << " "
			<< ETerm_equil << " " << WTerm << " " << correctionTerm << " "
		<< std::endl;
	}

	// ----------------------------------------------------------------
	// EVALUATE
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	SimTK::Real log_p_accept = 1.0;
	/* if(getRunType() > 0){
		log_p_accept = ETerm_nonequil + std::log(correctionTerm);
	}else{
		log_p_accept = ETerm_equil ;
	} */

	if(getRunType() == RunType::REMC){

		log_p_accept = ETerm_equil ;

	}else if(getRunType() == RunType::RENEMC){

		log_p_accept = ETerm_nonequil + std::log(correctionTerm) ;

	}else if( getRunType() == RunType::RENE){

		log_p_accept = WTerm + std::log(correctionTerm) ;

	}

	// Draw from uniform distribution
	SimTK::Real unifSample = uniformRealDistribution(randomEngine);

	// Accept or reject
	if((log_p_accept >= 0.0) || (unifSample < std::exp(log_p_accept))){

		// Update replicas coordinates from work generated coordinates
		set_WORK_CoordinatesAsFinal(replica_X);
		set_WORK_CoordinatesAsFinal(replica_Y);

		// Update replica's energy from work last potential energy
		set_WORK_PotentialAsFinal(replica_X);
		set_WORK_PotentialAsFinal(replica_Y);

		// Swap thermodynamic states
		swapThermodynamicStates(replica_X, replica_Y);
		swapPotentialEnergies(replica_X, replica_Y);

		std::cout << "swapped\n" << endl;

		returnValue = true;

	}else{

		// Return to equilibrium worlds coordinates
		// - no need because it is restored in RunRENS

		// Return to equilibrium worlds energies
		// - no need because it is restored in RunRENS

		// Don't swap thermodynamics states nor energies

		std::cout << "left\n" << endl;

		returnValue = false;
	}

	return returnValue;

}

// Exhange all replicas
void Context::mixAllReplicas(int nSwapAttempts)
{
	// Get a random number generator
	//std::random_device rd; // obtain a random number
	//std::mt19937 randomEngine(rd()); // seed the generator
	std::uniform_int_distribution<std::size_t>
		randReplicaDistrib(0, nofReplicas-1);

	// Try nSwapAttempts to swap between random replicas
	for(size_t swap_k = 0; swap_k < nSwapAttempts; swap_k++){

		// Get two random replicas
		auto replica_i = randReplicaDistrib(randomEngine);
		auto replica_j = randReplicaDistrib(randomEngine);
		std::cout << "Attempt to swap replicas " << replica_i
			<< " and " << replica_j << std::endl;

		// Attempt to swap
		bool swapped = false;

/* 		if( getRunType() == 1 ){ // TODO make it enum
			swapped = attemptREXSwap(replica_i, replica_j);

		}else if ( getRunType() == 2 ){
			swapped = attemptRENSSwap(replica_i, replica_j);
		} */

		swapped = attemptREXSwap(replica_i, replica_j);

	}
}

// Mix neighboring replicas
// Thermodyanmic states are fixed; replicas are variables
void Context::mixNeighboringReplicas(unsigned int startingFrom)
{
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

		/* std::cout << "mixNeighboringReplicas thermoStates "
			<< thermoState_i << " " << thermoState_j << "\n";
		std::cout << "Attempt to swap replicas " << replica_i
			<< " and " << replica_j << std::endl << std::flush; */

		// Attempt to swap
		bool swapped = false;
		/* if( getRunType() == 1 ){ // TODO make it enum
			swapped = attemptREXSwap(replica_i, replica_j);

		}else if ( getRunType() == 2 ){
			swapped = attemptRENSSwap(replica_i, replica_j);
		} */

		swapped = attemptREXSwap(replica_i, replica_j);

	}
}

// Mix replicas - hold it for the moment
void Context::mixReplicas()
{
	assert(!"Not implemented"); throw std::exception();
//	if(replicaMixingScheme == ReplicaMixingScheme::neighboring){
//		mixNeighboringReplicas();
//	}else{
//		mixAllReplicas(nofReplicas*nofReplicas*nofReplicas);
//	}
}

// Load replica's atomLocations into it's front world
int Context::restoreReplicaCoordinatesToFrontWorld(int whichReplica)
{

	//std::cout <<  "Context::restoreReplicaCoordinatesToFrontWorld " << whichReplica << ": " << std::flush;

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

	worlds[currWorldIx].setAtomsLocationsInGround(state,
		replicas[whichReplica].getAtomsLocationsInGround());

	//std::cout << "Context::restoreReplicaCoordinatesToFrontWorld worldIndexes.front() " << worldIndexes.front() << std::endl << std::flush;

	return currWorldIx;

}

// Load replica's atomLocations into it's back world
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

	worlds[worldIndexes.back()].setAtomsLocationsInGround(state,
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
		worlds[frontWorldIx].CalcPotentialEnergy();


	// Set this replica's energy
	replicas[replicaIx].set_WORK_PotentialEnergy_New(energy);

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
		worlds[frontWorldIx].CalcPotentialEnergy();

	// Add the Fixman potential to the energy (DANGEROUS)
	//energy += pHMC((worlds[backWorldIx].samplers[0]))->fix_set;

	// Set this replica's energy
	replicas[replicaIx].setPotentialEnergy(energy);

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


void Context::initializeReplica(int thisReplica)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thisThermoStateIx = replica2ThermoIxs[thisReplica];

	// Get this world indexes from the corresponding thermoState
	std::vector<int> replicaWorldIxs = thermodynamicStates[thisThermoStateIx].getWorldIndexes();
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
		worlds[replicaWorldIxs[i]].updSampler(0)->setTimestep(
			replicaTimesteps[i]);
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
	// =============

	// // For some reason I run the worlds a few times
	// for(size_t ri = 0; ri < 1; ri++){
	// 	for(std::size_t worldIx = 0; worldIx < replicaNofWorlds; worldIx++){

	// 		int frontIx = -1;
	// 		int backIx = -1;

	// 		// -------------
	// 		// SAMPLE from the current world
	// 		frontIx = replicaWorldIxs.front();
	// 		//std::cout << "Sample world " << front << "\n";
	// 		int accepted = worlds[frontIx].generateSamples(0);
	// 		// =============

	// 		// -------------
	// 		// ROTATE
	// 		///*print*/std::cout << "Rotate from";/*print*/
	// 		///*print*/for(int k = 0; k < replicaNofWorlds; k++){std::cout << " " << replicaWorldIxs[k];}/*print*/

	// 		// Rotate worlds indices (translate from right to left)
	// 	   	std::rotate(replicaWorldIxs.begin(),
	// 			replicaWorldIxs.begin() + 1,
	// 			replicaWorldIxs.end());

	// 		///*print*/std::cout << " to";/*print*/
	// 		///*print*/for(int k = 0; k < replicaNofWorlds; k++){std::cout << " " << replicaWorldIxs[k];}/*print*/
	// 		///*print*/std::cout << "\n";/*print*/
	// 		// =============

	// 		// -------------
	// 		// TRANSFER coordinates from last world to current
	// 		// TODO: eliminate in the last iteration
	// 		frontIx = replicaWorldIxs.front();
	// 		backIx = replicaWorldIxs.back();

	// 		if(replicaNofWorlds > 1) {
	// 			//std::cout << "Transfer from world " << backIx
	// 			//	<< " to " << frontIx << std::endl;

	// 			transferCoordinates(backIx, frontIx);
	// 		}
	// 		// =============

	// 	} // END iteration through worlds
	// } // END iteration through rounds


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

// Set world distort parameters
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

	std::cout << "Temperature set to " << T << std::endl << std::flush;

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
		worlds[replicaWorldIxs[i]].updSampler(0)->setTimestep(
			replicaTimesteps[i]);
	}

	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setMDStepsPerSample(
			replicaMdsteps[i]);
	}

	// Print info
	std::cout << "Timesteps set to ";
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		std::cout
			<< worlds[replicaWorldIxs[i]].getSampler(0)->getTimestep()
			<< " " ;
	}
	std::cout << std::endl;

	std::cout << "Mdsteps set to ";
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		std::cout 
			<< worlds[replicaWorldIxs[i]].getSampler(0)->getMDStepsPerSample()
			<< " " ;
	}
	std::cout << std::endl;
	// =============

}

// TODO turn strings into enum
SimTK::Real Context::distributeScalingFactor(
	std::vector<std::string> how, SimTK::Real scalefactor, bool randSignOpt)
{

	// Deterministic
	if (std::find(how.begin(), how.end(), "deterministic") != how.end()){
		// Do nothing
	}
	
	// Truncated normal
	if (std::find(how.begin(), how.end(), "Gauss") != how.end()){
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
	if (std::find(how.begin(), how.end(), "uniform") != how.end()){
		scalefactor = 
			worlds[0].updSampler(0)->uniformRealDistributionRandTrunc(
				0.8, 1.25); //0.625, 1.600);
	}

	// Assign a random direction: stretch or compress
	if (std::find(how.begin(), how.end(), "Bernoulli") != how.end()){
		std::cout << "SFdistrib Bernoulli "
			<< std::endl;

		SimTK::Real randDir =
			worlds[0].updSampler(0)->uniformRealDistribution_m1_1(
				randomEngine);
		scalefactor = (randDir > 0) ? scalefactor : (1.0/scalefactor) ;
	}

	// Assign a random sign (optional)
	if(randSignOpt){
		SimTK::Real randSign;
		SimTK::Real randUni_m1_1 =
			worlds[0].updSampler(0)->uniformRealDistribution_m1_1(
				randomEngine);
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

// Run a particular world
bool Context::RunWorld(int whichWorld)
{
	// == SAMPLE == from the current world
	bool validated = false;

	// Equilibrium world
	if(NDistortOpt[whichWorld] == 0){

		validated = worlds[whichWorld].generateSamples(
			int(nofSamplesPerRound[whichWorld]));

	// Non-equilibrium world
	}else if(NDistortOpt[whichWorld] == -1){

		// Generate samples
		validated = worlds[whichWorld].generateSamples(
			int(nofSamplesPerRound[whichWorld]));

	}
	return validated;

}

// Rewind back world
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
		transferCoordinates(frontIx, backIx);
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
	validated = RunWorld(frontWorldIx);

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
		std::cout << "Transfer from world " << backWorldIx << " to " << frontWorldIx ;
		transferCoordinates(backWorldIx, frontWorldIx);

		if(validated){
			std::cout << std::endl;
		}else{
			std::cout << " invalid sample." << std::endl;
		}

	}

	return worldIxs.front();

}

// Go through all of this replica's worlds and generate samples
int Context::RunReplicaAllWorlds(int thisReplica, int howManyRounds)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thisThermoStateIx = replica2ThermoIxs[thisReplica];

	// Get this world indexes from the corresponding thermoState
	std::vector<int> &replicaWorldIxs = 
		thermodynamicStates[thisThermoStateIx].updWorldIndexes();

	// Get nof worlds in this replica
	size_t replicaNofWorlds = replicaWorldIxs.size();

	// Set thermo and simulation parameters for the worlds in this replica
	setReplicasWorldsParameters(thisReplica);

	// Go through the requested nof rounds
	int frontWIx = -1;
	for(size_t ri = 0; ri < howManyRounds; ri++){

		// Go through each world of this replica
		for(std::size_t worldCnt = 0; worldCnt < replicaNofWorlds; worldCnt++){

			// Run front world
			frontWIx = RunFrontWorldAndRotate(replicaWorldIxs);

		} // END iteration through worlds

	} // END iteration through rounds

	return frontWIx;

}

/* vector<string> split(const string& i_str, const string& i_delim)
{
    vector<string> result;

    size_t found = i_str.find(i_delim);
    size_t startIndex = 0;

    while(found != string::npos)
    {
        result.push_back(string(i_str.begin()+startIndex, i_str.begin()+found));
        startIndex = found + i_delim.size();
        found = i_str.find(i_delim, startIndex);
    }
    if(startIndex != i_str.size())
        result.push_back(string(i_str.begin()+startIndex, i_str.end()));
    return result;
} */

/**
 * Update the scale factors
 */ 
void Context::updQScaleFactors(int mixi)
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
	bool randSignOpt = false;

	// The names of the probability distributions operators
	std::vector<std::string> how;

	// Go through all thermodynamic states which should correspond to
	// scale factors
	for(size_t thermoIx = 0; thermoIx < qScaleFactors.size(); thermoIx++){

		// Get thermoState world indexes 
		std::vector<int>& worldIxs = 
			thermodynamicStates[thermoIx].updWorldIndexes();
		size_t thermoNofWorlds = worldIxs.size();

		// Go through all the worlds this thermostate should go through
		for(std::size_t worldCnt = 0; worldCnt < thermoNofWorlds; worldCnt++){

			// Assign distribution from the first nonequilibrium world FOR NOW
			if(thermodynamicStates[thermoIx].getDistortOptions()[worldCnt] != 0){

				// Set the way the scaling factors are distributed
				if(  thermodynamicStates[thermoIx].getDistortArgs().size()  ){

					how = split(thermodynamicStates[thermoIx].getDistortArgs()[
						worldCnt], "_");

				}else{
					
					how = {"deterministic"};
				}

				break;
			}

			// Distribute scale factor
			if(qScaleFactors.at(thermoIx) != 1){ // This is questionable
				qScaleFactors.at(thermoIx) = distributeScalingFactor(
					how, qScaleFactorsMiu.at(thermoIx), randSignOpt);
			}

		}
	}

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


// Run replica exchange protocol
void Context::RunREX()
{

	// Is this necesary =======================================================
	realizeTopology(); 

	// Allocate space for swap matrices
	allocateSwapMatrices();

	// Run each replica one time initially
	// std::cout << " REX batch " << 0 << std::endl;
	for (size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++){
	// 	std::cout << "REX replica " << replicaIx << std::endl;
	//
	// 		// Set intial parameters
	 		initializeReplica(replicaIx);
	//
	// 		// Copy coordinates from replica to front world
	// 		restoreReplicaCoordinatesToFrontWorld(replicaIx);
	//
	// 		// Iterate this replica's worlds
	// 		RunReplicaAllWorlds(replicaIx, swapEvery);
	//
	// 		// Copy coordinates from front world to replica
	// 		storeReplicaCoordinatesFromFrontWorld(replicaIx);
	//
	// 		// Store energy
	// 		storeReplicaEnergyFromFrontWorldFull(replicaIx);
	// 		storeReplicaFixmanFromBackWorld(replicaIx);
	//
	// 		/* PrintToLog(replicaIx, worldIndexes.front(), 0);
	// 		writePdbs(0,	replica2ThermoIxs[replicaIx]); */
	//
	// 		// Write energy and geometric features to logfile
	// 		REXLog(0, replicaIx);
	//
	} // ===================================================================

	// PrintNofAcceptedSwapsMatrix();

	bool givenTsMode = false;

	int nofMixes = int(requiredNofRounds / swapEvery);
	int currFrontWIx = -1;

	// REPLICA EXCHANGE MAIN LOOP -------------------------------------------->
	for(size_t mixi = 0; mixi < nofMixes; mixi++){

		std::cout << " REX batch " << mixi << std::endl;
		nofRounds = mixi;

		updQScaleFactors(mixi);

		// SIMULATE EACH REPLICA --------------------------------------------->
		for (size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++){
			std::cout << "REX replica " << replicaIx << std::endl;

// std::cout << "WorkCoordinate   after (1.0) replica " << replicaIx << " work " << std::endl;		// DEBUG !!!!!!!!!!!!!!!!!!!!!!
// replicas[replicaIx].Print_WORK_Coordinates();													// DEBUG !!!!!!!!!!!!!!!!!!!!!!
// std::cout << "================ after (1.0) replica " << replicaIx << " work \n" << std::flush;	// DEBUG !!!!!!!!!!!!!!!!!!!!!!
// std::cout << "Context::RunREX OMM (1.0) world " << 0 << std::endl;									// DEBUG !!!!!!!!!!!!!!!!!!!!!!
// (worlds[0].getSampler(0))->OMM_PrintLocations();													// DEBUG !!!!!!!!!!!!!!!!!!!!!!
// std::cout << std::flush;																			// DEBUG !!!!!!!!!!!!!!!!!!!!!!
// std::cout << "Context::RunREX Sim (1.0) world " << 0 << std::endl;									// DEBUG !!!!!!!!!!!!!!!!!!!!!!
// PrintCoordinates(worlds[0].getCurrentAtomsLocationsInGround());										// DEBUG !!!!!!!!!!!!!!!!!!!!!!
// std::cout << "Context::RunREX (1.0) world " << 0 << "================\n" << std::flush;				// DEBUG !!!!!!!!!!!!!!!!!!!!!!

			// ========================== LOAD ========================
			// Load the front world
			currFrontWIx = restoreReplicaCoordinatesToFrontWorld(replicaIx);                  // (1)

			// Set non-equilibrium parameters: get scale factors
			updWorldsDistortOptions(replicaIx);

			// Set thermo and simulation parameters for the worlds in this replica
			setReplicasWorldsParameters(replicaIx);

			// ----------------------------------------------------------------
			// EQUILIBRIUM
			// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

					// ======================== SIMULATE ======================
					// Includes worlds rotations
					currFrontWIx = RunReplicaEquilibriumWorlds(replicaIx, swapEvery);         // (2) + (3) + (4)

					// Write energy and geometric features to logfile
					REXLog(mixi, replicaIx);

					// ========================= UNLOAD =======================
					// Deposit front world coordinates into the replica
					storeReplicaCoordinatesFromFrontWorld(replicaIx);                         // (4.5)

					// Deposit energy terms
					storeReplicaEnergyFromFrontWorldFull(replicaIx);

			// ----------------------------------------------------------------
			// NON-EQUILIBRIUM
			// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

					// ======================== SIMULATE ======================
					// Includes worlds rotations                                              // (5)
					currFrontWIx = RunReplicaNonequilibriumWorlds(replicaIx, swapEvery);

					// ========================= UNLOAD =======================
					replicas[replicaIx].setTransferedEnergy( calcReplicaWork(replicaIx) );

					// Deposit work coordinates into the replica
					store_WORK_CoordinatesFromFrontWorld(replicaIx);                          // (8)

					// Deposit energy terms
					store_WORK_ReplicaEnergyFromFrontWorldFull(replicaIx);

					// Store any transformation Jacobians contribution
					store_WORK_JacobianFromBackWorld(replicaIx);			

					// Store Fixman if required
					storeReplicaFixmanFromBackWorld(replicaIx);

					if(currFrontWIx != 0){std::cout << "=== RUN FIRST WORLD NOT 0 === " << currFrontWIx << std::endl;}

		} // end replicas simulations

		// Mix replicas
		if(getRunType() != RunType::Default){                                                 // (9) + (10)
			if(replicaMixingScheme == ReplicaMixingScheme::neighboring){
				if ((mixi % 2) == 0){
					mixNeighboringReplicas(0);
				}else{
					mixNeighboringReplicas(1);
				}
			}else{
				mixAllReplicas(nofReplicas*nofReplicas*nofReplicas);
			}

			PrintNofAcceptedSwapsMatrix();
		} 

	} // end rounds

	//PrintNofAttemptedSwapsMatrix();
	PrintNofAcceptedSwapsMatrix();
	//PrintReplicaMaps();

}

/**
 * Go through the vector of worlds and if their equilibrium worlds run and 
 * rotate. 
*/
int Context::RunReplicaEquilibriumWorlds(int replicaIx, int swapEvery)
{
	// Get thermoState corresponding to this replica
	int thisThermoStateIx = replica2ThermoIxs[replicaIx];

	// Get this world indexes from the corresponding thermoState
	std::vector<int>& replicaWorldIxs = 
		thermodynamicStates[thisThermoStateIx].updWorldIndexes();

	// Get nof worlds in this replica
	size_t replicaNofWorlds = replicaWorldIxs.size();

	/* std::cout << "Context::RunReplicaEquilibriumWorlds replicaWorldIxs before ";
	PrintCppVector(replicaWorldIxs, " | ", "|\n"); */

	// Run the equilibrium worlds
	int currFrontWIx = -1;
	for(std::size_t worldCnt = 0; worldCnt < replicaNofWorlds; worldCnt++){

		if( thermodynamicStates[thisThermoStateIx].getDistortOptions()[worldCnt]
		== 0){

			// Run front world
			currFrontWIx = RunFrontWorldAndRotate(replicaWorldIxs);
		
		}
	} // END iteration through worlds

	/* std::cout << "Context::RunReplicaEquilibriumWorlds replicaWorldIxs after ";
	PrintCppVector(replicaWorldIxs, " | ", "|\n"); */

	return currFrontWIx;

}

int Context::RunReplicaNonequilibriumWorlds(int replicaIx, int swapEvery)
{
	/* std::cout << "Context::RunReplicaNonequilibriumWorlds\n"; */
	// Get thermoState corresponding to this replica
	int thisThermoStateIx = replica2ThermoIxs[replicaIx];

	// Get this world indexes from the corresponding thermoState
	std::vector<int>& replicaWorldIxs = 
		thermodynamicStates[thisThermoStateIx].updWorldIndexes();

	// Get nof worlds in this replica
	size_t replicaNofWorlds = replicaWorldIxs.size();

	/* std::cout << "Context::RunReplicaNonquilibriumWorlds replicaWorldIxs before ";
	PrintCppVector(replicaWorldIxs, " | ", "|\n"); */

	// Run the non-equilibrium worlds
	int frontWIx = -1;
	for(std::size_t worldCnt = 0; worldCnt < replicaNofWorlds; worldCnt++){

		if(thermodynamicStates[thisThermoStateIx].getDistortOptions()[worldCnt]
		!= 0){

			// Run front world
			frontWIx = RunFrontWorldAndRotate(replicaWorldIxs);
		
		}
	} // END iteration through worlds

	/* std::cout << "Context::RunReplicaNonquilibriumWorlds replicaWorldIxs before ";
	PrintCppVector(replicaWorldIxs, " | ", "|\n"); */

	return frontWIx;
}

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


// Print info about all the replicas and thermo states
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

	std::cout << "Number of accepted swaps matrix:\n";
	for(size_t i = 0; i < M; i++){
		for(size_t j = 0; j < M; j++){
			std::cout << nofAcceptedSwapsMatrix[i][j] << " ";
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
			worlds[worldIndexes.front()].updateAtomListsFromCompound(pdbState);

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

void Context::randomizeWorldIndexes()
{
	// Random int for random world order
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<std::size_t>
		randWorldDistrib(1, nofWorlds-1); // TODO between 1 and nOfWorlds-1?

	if(getNofWorlds() >= 3){

		// Swap world indeces between vector position 2 and random
		auto randVecPos = randWorldDistrib(gen);
		//std::cout << "Swapping position 1 with "
		//	<< randVecPos << std::endl;

		auto secondWorldIx = worldIndexes[1];
		auto randWorldIx = worldIndexes[randVecPos];

		worldIndexes[1] = randWorldIx;
		worldIndexes[randVecPos] = secondWorldIx;

	}
}

void Context::transferCoordinates(int src, int dest)
{
	// Get advanced states of the integrators
	SimTK::State& lastAdvancedState = updWorld(src)->integ->updAdvancedState();
	SimTK::State& currentAdvancedState = updWorld(dest)->integ->updAdvancedState();

	// Get coordinates from source
	const std::vector<std::vector<std::pair<
		bSpecificAtom *, SimTK::Vec3> > >&
		otherWorldsAtomsLocations =
	updWorld(src)->getAtomsLocationsInGround(lastAdvancedState);

	// Pass compounds to the new world
	passTopologiesToNewWorld(dest);

	// Set coordinates to destination (reconstruct)
	currentAdvancedState = updWorld(dest)->setAtomsLocationsInGround(
		currentAdvancedState, otherWorldsAtomsLocations);
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
			std::vector<std::string> how = { "deterministic", "Bernoulli"};
			bool randSignOpt = false;
			qScaleFactor = distributeScalingFactor(how, qScaleFactor, randSignOpt);

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

// Print recommended timesteps
void Context::PrintInitialRecommendedTimesteps(void)
{
	std::cout << "Initial recommended timesteps ";
	for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
		std::cout << "World " << worldIx << " " ;
;
		worlds[worldIx].PrintInitialRecommendedTimesteps();
	}
	std::cout << std::endl;
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
			setTemperature(currT);
			std::cout << "T= " << currT << std::endl;

			RunOneRound();

			RunLog(round);

			this->nofRounds++;
		}
	}


}

// Set number of threads
void Context::setNumThreadsRequested(std::size_t which, int howMany)
{
	std::cout << "Robosample requested " << howMany << " threads " << std::endl;
	if (howMany == 1){
		worlds[which].updForceField()->setUseMultithreadedComputation(false);
	}else{
		worlds[which].updForceField()->setNumThreadsRequested(howMany);
	}
}

void Context::setUseOpenMMAcceleration(bool arg)
{
	for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
		worlds[worldIx].updForceField()->setUseOpenMMAcceleration(arg);
	}
}


void Context::setUseOpenMMIntegration(std::size_t which, Real temperature, Real stepsize)
{
	worlds[which].updForceField()->setUseOpenMMIntegration(true);
	worlds[which].updForceField()->setDuMMTemperature(temperature);
	worlds[which].updForceField()->setOpenMMstepsize(stepsize);
}

void Context::setNonbondedMethod(std::size_t which, int methodInx)
{
    if (methodInx >= 0 && methodInx <= 1 ){
        worlds[which].updForceField()->setNonbondedMethod(methodInx);
    }else{
        std::cout<< "Invalid nonbonded method. (0 = nocutoff; 1 = cutoffNonPeriodic). Default NoCutoff method will be used." << endl;
    }
}

void Context::setUseOpenMMCalcOnlyNonBonded(bool arg)
{

    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
        worlds[worldIx].updForceField()->setUseOpenMMCalcOnlyNonBonded(arg);
    }
}

void Context::setNonbondedCutoff(std::size_t which, Real cutoffNm)
{
    if (cutoffNm >= 0 ){
        worlds[which].updForceField()->setNonbondedCutoff(cutoffNm);
    }else{
        std::cout<< "Negative cutoff requested. Default cutoff = 2.0 nm will be used instead" << endl;
    }
}


/** Get/Set seed for reproducibility. **/
void Context::setSeed(std::size_t whichWorld, std::size_t whichSampler, uint32_t argSeed)
{
	worlds[whichWorld].updSampler(whichSampler)->setSeed(argSeed);
}

uint32_t Context::getSeed(std::size_t whichWorld, std::size_t whichSampler) const
{
	return worlds[whichWorld].getSampler(whichSampler)->getSeed();
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
			logFile << std::fixed << std::setprecision(3) 
				<< Dihedral(dihedralIx[0], dihedralIx[1], 0,
				dihedralIx[2], dihedralIx[3], dihedralIx[4], dihedralIx[5]) << " ";
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
	SimTK::State& advancedState =
		(updWorld(currentWorldIx))->integ->updAdvancedState();

	constexpr int mc_step = -1;

	// Pass compounds to the new world
	passTopologiesToNewWorld(currentWorldIx);

	// 
	(updWorld(currentWorldIx))->updateAtomListsFromCompound(advancedState);
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
		worlds[worldIndexes.front()].updateAtomListsFromCompound(pdbState);

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
void Context::WritePdb(std::size_t)
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

std::string Context::getPdbPrefix()
{
	return this->pdbPrefix;
}

void Context::setPdbPrefix(std::string arg)
{
	this->pdbPrefix = arg;
}

SimTK::Real Context::Roboangle(std::size_t whichWorld,
	std::size_t whichCompound, std::size_t whichSampler,
	int a1, int a2, int a3)
{

	SimTK::State& state = worlds[whichWorld].integ->updAdvancedState();

	Topology& topology = worlds[whichWorld].updTopology(whichCompound);

	SimTK::DuMMForceFieldSubsystem& dumm = *(worlds[whichWorld].forceField);

	SimTK::SimbodyMatterSubsystem& matter = *(worlds[whichWorld].matter);

	SimTK::Vec3 a1pos, a2pos, a3pos, a4pos;
	a1pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)),
		dumm, matter, state);
	a2pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)),
		dumm, matter, state);
	a3pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)),
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
	a1pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)),
		dumm, matter, state);
	a2pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)),
		dumm, matter, state);
	a3pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)),
		dumm, matter, state);
	a4pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a4)),
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

	a1pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)),
		dumm, matter, state);
	a2pos = topology.calcAtomLocationInGroundFrameThroughSimbody(
		SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)),
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
				(topologies[molIx].bAtomList[k]).getCompoundAtomIndex();

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
						
						(updWorld(worldIx))->addContacts(
								cliqueAtomIxs[contactCliqueIx][topologyIx],
								topologyIx,
								SimTK::ContactCliqueId(contactCliqueIx));
					}
				}
			}
		}

		// Add membrane to all worlds.
		for(unsigned int worldIx = 0; worldIx < getNofWorlds(); worldIx++){
			(updWorld(worldIx))->addMembrane(memZWidth);
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



