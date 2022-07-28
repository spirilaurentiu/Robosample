#include "Context.hpp"

// Includes to get the structure of additional classes

#include "World.hpp"
#include "Sampler.hpp"

// Helper classes for REX

//////////////////////////
// CLASS THERMODYNAMICSTATE
//////////////////////////

class ThermodynamicState{
  public:
	ThermodynamicState();
	ThermodynamicState(int index);

	ThermodynamicState(int index, const SimTK::Real& T
		,std::vector<int>& argWorldIndexes
		,std::vector<SimTK::Real>& argTimesteps
		,std::vector<int>& argMdsteps
	);

	~ThermodynamicState(){}

	const int getIndex(void){return myIndex;}
	void setIndex(const int& someIndex){myIndex = someIndex;}

	void setTemperature(const SimTK::Real& T);
	const SimTK::Real& getTemperature();
	const SimTK::Real& getBeta();

	const std::vector<int> getWorldIndexes(void);
	const std::vector<SimTK::Real> getTimesteps(void);
	const std::vector<int> getMdsteps(void);

	void Print(void);

  private:

	// Index
	int myIndex;

	// Temperature
	SimTK::Real temperature;
	SimTK::Real beta;

	// Worlds related parameters //R2TNEW
	std::vector<int> worldIndexes; //R2TNEW
	std::vector<SimTK::Real> timesteps; //R2TNEW
	std::vector<int> mdsteps; //R2TNEW

};

ThermodynamicState::ThermodynamicState()
{
	myIndex = 0;
	temperature = 300;
	SimTK::Real RT = temperature *
		static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	beta = 1.0 / RT;
}

ThermodynamicState::ThermodynamicState(int index)
{
	myIndex = index;
	temperature = 300;
	SimTK::Real RT = temperature *
		static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	beta = 1.0 / RT;
}

ThermodynamicState::ThermodynamicState(int index, const SimTK::Real& T
	,std::vector<int>& argWorldIndexes
	,std::vector<SimTK::Real>& argTimesteps
	,std::vector<int>& argMdsteps
	)
{
	// Own index
	myIndex = index;

	// Temperature related
	temperature = T;
	SimTK::Real RT = temperature *
		static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	beta = 1.0 / RT;

	// Worlds related parameters
	worldIndexes = argWorldIndexes;
	timesteps = argTimesteps;
	mdsteps = argMdsteps;
}

void ThermodynamicState::setTemperature(const SimTK::Real& T)
{
	temperature = T;
	SimTK::Real RT = temperature *
		static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	beta = 1.0 / RT;
}

const SimTK::Real& ThermodynamicState::getTemperature()
{
	return temperature;
}

const SimTK::Real& ThermodynamicState::getBeta()
{
	return beta;
}

//
const std::vector<int> ThermodynamicState::getWorldIndexes(void)
{
	return worldIndexes;
}

const std::vector<SimTK::Real> ThermodynamicState::getTimesteps(void)
{
	return timesteps;
}

const std::vector<int> ThermodynamicState::getMdsteps(void)
{
	return mdsteps;
}

void ThermodynamicState::Print(void)
{
	std::cout << "ThermodynamicState::Print index T "
		<< myIndex << " " << temperature
		<< std::endl;

	std::cout << "ThermodynamicState::Print index worldIndexes " << myIndex;
	for(auto worldIndex : worldIndexes){
		std::cout << " " << worldIndex;
	}
	std::cout << std::endl;
}


//////////////////////////
// CLASS REPLICA
//////////////////////////

class Replica{
  public:
	Replica();
	Replica(int index);
	Replica(int index,
		std::vector<int>& argWorldIndexes);
	Replica(int index,
		std::vector<int>& argWorldIndexes,
		std::vector<SimTK::Real>& argTimesteps,
		std::vector<int>& argMdsteps);
	~Replica();

	const	std::vector<
		std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>&
		getAtomsLocationsInGround();

	// Reserve memory and set values
	void setAtomsLocationsInGround(std::vector<
		std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>&);

	// This assumes allocation has been done already
	void updAtomsLocationsInGround(std::vector<
		std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>);

	// Load atomLocations coordinates into the front world
	void restoreCoordinates(void){}

	// Stores coordinates from front world into atomsLocations
	void storeCoordinates(void){}

	const SimTK::Real getPotentialEnergy(void){return potential;}
	void setPotentialEnergy(const SimTK::Real& somePotential){potential = somePotential;}

	void Print(void){}
	void PrintCoordinates(void);

  private:

	int myIndex;

	//// Worlds related parameters //R2TOLD
	//std::vector<int> worldIndexes; //R2TOLD
	//std::vector<SimTK::Real> timesteps; //R2TOLD
	//std::vector<int> mdsteps; //R2TOLD

	// Replica configurations
	std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>
		atomsLocations;

	// Replica potential energy
	SimTK::Real potential; // TODO: turn into a vector for worlds
};

Replica::Replica()
{
	myIndex = 0;
}

Replica::Replica(int argIndex)
{
	myIndex = argIndex;
}

Replica::~Replica()
{
}

// Get coordinates from this replica
const	std::vector<
	std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>&
	Replica::getAtomsLocationsInGround()
{
	return atomsLocations;
}

// Set the coordinates of this replica
// Also allocate memory
void Replica::setAtomsLocationsInGround(
	std::vector< // topologies
	std::vector< // one topology
	std::pair <bSpecificAtom *, SimTK::Vec3>>>& // atom location
		otherAtomsLocations
	)
{

	for (auto& topology : otherAtomsLocations){
		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> currentTopologyInfo;
		currentTopologyInfo.reserve(topology.size());

		for(auto& otherAtom : topology){
			bSpecificAtom* batom = otherAtom.first;
			SimTK::Vec3 location = otherAtom.second;
			std::pair<bSpecificAtom*, SimTK::Vec3> atomLocationPair(batom, location);
			currentTopologyInfo.push_back(atomLocationPair);
		}

		atomsLocations.emplace_back(currentTopologyInfo);

	}

	//atomsLocations = otherAtomsLocations;
}

// Update the coordinates of this replica
void Replica::updAtomsLocationsInGround(
	std::vector< // topologies
	std::vector< // one topology
	std::pair <bSpecificAtom *, SimTK::Vec3>>> // atom location
		otherAtomsLocations
	)
{
	int i = -1;
	for (auto& topology : otherAtomsLocations){
		i += 1;

		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>&
			myTopology = atomsLocations[i];

		int j = -1;
		for(auto& otherAtom : topology){
			j += 1;

			std::pair<bSpecificAtom *, SimTK::Vec3>&
			myAtom = myTopology[j];

			myAtom.first  = otherAtom.first;
			myAtom.second = otherAtom.second;

		}

	}

}

void Replica::PrintCoordinates(void)
{
	for(auto& topology : atomsLocations){
		for(auto& atomCoordinates : topology){
			std::cout
			//<< (atomCoordinates.first)->inName
			<< (atomCoordinates.first)->getCompoundAtomIndex() << " "
			<< atomCoordinates.second
			<< std::endl;
		}
	}
}
//////////////////////////
// CLASS CONTEXT
//////////////////////////

//Primitive error handling
void Context::throwAndExit(std::string errMsg, int errCode){
		std::cerr << "Context: " << errMsg;
		throw std::exception();
		std::exit(errCode);
}


// Default constructor
Context::Context(const SetupReader& setupReader, std::string logFilename)
: thermodynamicStates(), replicas()
{
	nofWorlds = 0;
	nofMols = 0;
	nofEmbeddedTopologies = 0;

	isWorldsOrderRandom = false;

	CheckInputParameters(setupReader);

	BUFSIZE = 1024 * 1048576; // 1048576
	buffer = std::make_unique<char[]>(BUFSIZE);
	logFile = fopen(logFilename.c_str(), "w+");
	if ( setvbuf(logFile, &buffer[0], _IOFBF, BUFSIZE) != 0){
	   perror("setvbuf()");
	}

	nofReplicas = 0;
	nofThermodynamicStates = 0;
	replicaMixingScheme = ReplicaMixingScheme::neighboring;
	swapEvery = 1;

}

// Destructor
Context::~Context(){
	fclose(logFile);

	//for(int i = 0; i < replicas.size(); i++){
	//	replicas[i].thermodynamicStateIx.clear();
	//}
}

// Check input
void Context::CheckInputParameters(const SetupReader& setupReader) {

	// Context specific parameters
	assert(SimTK::Pathname::fileExists(setupReader.get("OUTPUT_DIR")[0]));

	assert( std::stoi((setupReader.get("ROUNDS"))[0]) >= 0);
	assert( std::stoi((setupReader.get("ROUNDS_TILL_REBLOCK"))[0]) > 0);

	// TODO: if seed not provided, generate one and warn
	assert( std::stoi((setupReader.get("SEED"))[0]) > 0);

	if(setupReader.find("RANDOM_WORLD_ORDER")){
		if(setupReader.get("RANDOM_WORLD_ORDER").size() == 0){
			std::cerr << "Please specify RANDOM_WORLD_ORDER with TRUE/FALSE\n";
			throw std::exception();
			std::exit(1);
		}else{
			isWorldsOrderRandom = ((setupReader.get("RANDOM_WORLD_ORDER")[0] == "TRUE") ? true : false);
		}
	}else{
		std::cerr << "Please specify if the world order should be random through RANDOM_WORLD_ORDER keyword\n";
		throw std::exception();
		std::exit(1);
	}

	// Set these numbers first
	std::size_t inpNofWorlds = setupReader.get("WORLDS").size();
	std::size_t inpNofMols = setupReader.get("MOLECULES").size();
	std::size_t inpNofEmbeddedTopologies = inpNofWorlds * inpNofMols;

	if(setupReader.get("ROOTS").size() != inpNofEmbeddedTopologies){
		std::cerr << "Must have the same no. of root atoms as the no. of Topologies = nofWorlds x nofMolecules.\n";
		throw std::exception();
		std::exit(1);
	}

	if(setupReader.get("ROOT_MOBILITY").size() != setupReader.get("ROOTS").size()){
		std::cerr << "Must have the same no. of root mobilities as the no. of root atoms.\n";
		throw std::exception();
		std::exit(1);
	}

	// World Samplers specific parameters
	if(setupReader.get("SAMPLER").size() != inpNofWorlds){
		std::cerr << "Must have the same no. of samplers as the no. of worlds.\n";
		throw std::exception();
		std::exit(1);
	}

	if(inpNofWorlds > setupReader.get("TIMESTEPS").size()){
		std::cerr << "Must have the at least same no. of timesteps as the no. of worlds.\n";
		throw std::exception();
		std::exit(1);
	}else{
		for(std::size_t worldIx = 0; worldIx < inpNofWorlds; worldIx++){
			if( std::abs(std::stod(setupReader.get("TIMESTEPS")[worldIx])) <= 0.0000001){
				std::cout << "Warning: timestep for world " << worldIx << " is too small.\n";
			}
		}
	}

	// Molecule specific parameters
	for(std::size_t molIx = 0; molIx < inpNofMols; molIx++){
		if(!SimTK::Pathname::fileExists(
			setupReader.get("MOLECULES")[molIx] + std::string("/")
			+ setupReader.get("PRMTOP")[molIx]) ){
			throwAndExit("Molecule " + std::to_string(molIx) + " prmtop not found\n", 1);

		}
		if(!SimTK::Pathname::fileExists(
			setupReader.get("MOLECULES")[molIx] + std::string("/")
			+ setupReader.get("INPCRD")[molIx] + ".rst7") ){
			throwAndExit("Molecule " + std::to_string(molIx) + " inpcrd not found\n", 1);
			}
	}

	// Topology specific paramters
	for(std::size_t worldIx = 0; worldIx < inpNofWorlds; worldIx++){
		for(std::size_t molIx = 0; molIx < inpNofMols; molIx++){
			if(!SimTK::Pathname::fileExists(
				setupReader.get("MOLECULES")[molIx] + std::string("/")
				+ setupReader.get("RBFILE")[molIx]) ){
				throwAndExit("world " + std::to_string(worldIx) +
					" molecule " + std::to_string(molIx) +
					" rb not found\n", 1);
			}
			if(!SimTK::Pathname::fileExists(
				setupReader.get("MOLECULES")[molIx] + std::string("/")
				+ setupReader.get("FLEXFILE")[molIx]) ){
				throwAndExit("world " + std::to_string(worldIx) +
					" molecule " + std::to_string(molIx) +
					" flex not found\n", 1);
			}

		}
	}

	// Normal modes options
	NMAOption.resize(inpNofWorlds, 0);
	if(setupReader.find("NMA_OPTION")){
		if(setupReader.get("RANDOM_WORLD_ORDER").size() == 0){
			std::cerr << "The NMA_OPTION key is present. Please specify a value.\n";
			throw std::exception();
			std::exit(1);
		}else{
			for(std::size_t worldIx = 0; worldIx < inpNofWorlds; worldIx++){
				NMAOption[worldIx] = std::stoi(setupReader.get("NMA_OPTION")[worldIx]);
			}
		}
	}

	// If we got here we can set global variables
	// Reserve memory
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
			AddWorld(true, visualizerFrequencies[worldIx]);
		}else{
			AddWorld(false);
		}
	}

	if(argNofWorlds != nofWorlds){
		std::cerr << "Something went wrong while adding the world\n";
		throw std::exception();
		std::exit(1);
	}

	std::cout << "Added " << nofWorlds << " empty worlds" << std::endl;

	//
	//allocateReblockQsCache();

}

// Add an empty world to the context
World * Context::AddWorld(bool visual, SimTK::Real visualizerFrequency){

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

void Context::modelOneEmbeddedTopology(int whichTopology, int whichWorld, std::string rootMobilizer)
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

		loadTopologyFile( topFN );

		loadCoordinatesFile( crdFN );

		// Initialize an input reader
		readAmberInput amberReader;
		amberReader.readAmberFiles(crdFNs[molIx], topFNs[molIx]);

		moleculeCount++; // Used for unique names of molecules // from world
		std::string moleculeName = "MOL" + std::to_string(moleculeCount); // from world
		roots.emplace_back(argRoots[molIx]); // from world
		//rootMobilities.emplace_back("Pin"); // TODO: move to setflexibilities
		topologies.emplace_back(Topology{moleculeName}); // TODO is this ok? // from world
		// Set atoms properties from a reader: number, name, element, initial
		// name, force field type, charge, coordinates, mass, LJ parameters
		topologies[molIx].SetGmolAtomPropertiesFromReader(&amberReader); // from world
		// Set bonds properties from reader: bond indeces, atom neighbours // from world
		topologies[molIx].SetGmolBondingPropertiesFromReader(&amberReader); // from world
		// Set atoms Molmodel types (Compound::SingleAtom derived) based on // from world
		// their valence // from world
		topologies[molIx].SetGmolAtomsMolmodelTypesTrial(); // from world

		// Add Biotypes // from world
		topologies[molIx].bAddBiotypes(&amberReader); // from world

		// Build Robosample graph and Compound graph.
		// It also asigns atom indexes in Compound
		// This is done only once and it needs
		topologies[molIx].buildGraphAndMatchCoords(
			std::stoi(argRoots[0]));

		topologies[molIx].loadTriples();

		// Map of Compound atom indexes to Robosample atom indexes
		topologies[molIx].loadCompoundAtomIx2GmolAtomIx();

	}


}

// Loads parameters into DuMM, adopts compound by the CompoundSystem
// and loads maps of indexes
void Context::addDummParams(
	int requestedNofMols,
	SetupReader& setupReader
){
	//std::vector<std::string> argRoots = setupReader.get("ROOTS");

	// Iterate through molecules
	for(unsigned int molIx = 0; molIx < requestedNofMols; molIx++){

		// Initialize an input reader
		readAmberInput amberReader;
		amberReader.readAmberFiles(crdFNs[molIx], topFNs[molIx]);

		for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){

			// Pass current topology to the current world
			(updWorld(worldIx))->topologies = &topologies;

			// Add parameters in DuMM
			(updWorld(worldIx))->AddDummParams(molIx, &amberReader);
			//topologies[molIx].bAddDummParams(&amberReader,
			//	*((updWorld(worldIx))->forceField));
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


			// Add entry to flexibility filenames matrix
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
			modelOneEmbeddedTopology(molIx, worldIx, argRootMobilities[worldIx]);

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

// Load/store Mobilized bodies joint types in samplers
void Context::loadMbxsToMobilities(void)
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
void Context::printThermodynamics(void)
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
				//<< (context.updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((updWorld(worldIx))->integ->updAdvancedState())
				//<< (context.updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy(context.updAdvancedState(worldIx, samplerIx))
				<< " useFixmanPotential = "
				<< pHMC(updWorld(worldIx)->updSampler(samplerIx))->isUsingFixmanPotential()
				<< std::endl;
		}

	}
}

// Print Molmodel related information
void Context::PrintMolmodelAndDuMMTypes(void){
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
void Context::PrintSimbodyMobods(void){
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
void Context::checkAtomStationsThroughDumm(void)
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
SimTK::Real Context::getTemperature(std::size_t whichWorld) const {
	return worlds[whichWorld].temperature;
}

void  Context::setTemperature(std::size_t whichWorld, float someTemperature){
	std::cout << " Context::setTemperature for world "<< whichWorld << " " << someTemperature << std::endl;
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
	// function args were std::size_t whichWorld, std::size_t whichSampler, float someTemperature
	assert(!"Not implemented"); throw std::exception();
}
//------------

/////////////////////////
// --- Simulation parameters ---
/////////////////////////

BaseSampler * Context::addSampler(
	std::size_t whichWorld,
	SamplerName whichSampler)
{
	return worlds[whichWorld].addSampler(whichSampler);
}

BaseSampler * Context::addSampler(
	std::size_t whichWorld,
	std::string samplerName)
{

	if(samplerName == "VV"){
		BaseSampler *p = worlds[whichWorld].addSampler(SamplerName::HMC);
		pHMC(p)->setAlwaysAccept(true);
		return p;
	}else if(samplerName == "HMC"){
		BaseSampler *p = worlds[whichWorld].addSampler(SamplerName::HMC);
		pHMC(p)->setAlwaysAccept(false);
		return p;
	}

}

void Context::initializeSampler(std::size_t whichWorld,
	std::size_t whichSampler)
{
	worlds[whichWorld].updSampler(whichSampler)->initialize(
		worlds[whichWorld].integ->updAdvancedState());
}

// Amber like scale factors.
void Context::setAmberForceFieldScaleFactors(std::size_t whichWorld){
	worlds[whichWorld].setAmberForceFieldScaleFactors();
}

// Amber like scale factors.
void Context::setAmberForceFieldScaleFactors(void){
	for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
		worlds[worldIx].setAmberForceFieldScaleFactors();
	}
}

// Set a global scaling factor for the forcefield
void Context::setGlobalForceFieldScaleFactor(
	std::size_t whichWorld, SimTK::Real globalScaleFactor){
	worlds[whichWorld].setGlobalForceFieldScaleFactor(globalScaleFactor);
}

void Context::setGlobalForceFieldScaleFactor(SimTK::Real globalScaleFactor){
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
int Context::getNofMDStepsPerSample(std::size_t whichWorld, std::size_t whichSampler){
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
int Context::getNofRounds()
{
	return nofRounds;
}

void Context::setNofRounds(int argNofRounds)
{
	nofRounds = argNofRounds;
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
void Context::allocateReblockQsCache(void)
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		QsCache.push_back(std::vector<std::vector<SimTK::Real>>(roundsTillReblock));
		//std::cout << "Context::AddWorld QsCache size " << QsCache.size() << std::endl;
	}
}

// This seems wrong !!!
void Context::allocateReblockQsCacheQVectors(void){
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
int Context::getNofSamplesPerRound(std::size_t whichWorld)
{
	return nofSamplesPerRound[whichWorld];
}

// Set the number of samples returned by the sampler in one round
void Context::setNofSamplesPerRound(std::size_t whichWorld, int MCStepsPerRound)
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

// -- Main ---
void Context::Run(void)
{
	assert(!"Not implemented"); throw std::exception();
	throw std::exception();
}

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
	for(int round = 0; round < nofRounds; round++){ // Iterate rounds
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
				std::vector<std::vector<std::pair<
					bSpecificAtom *, SimTK::Vec3> > >
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
			for(int k = 0; k < getNofSamplesPerRound(currentWorldIx); k++){ // Iterate through samples
				updWorld(currentWorldIx)->updSampler(0)->sample_iteration(currentAdvancedState);
			} // END for samples

		} // for i in worlds

		// Print energy and geometric features
		if( !(round % getPrintFreq()) ){
			PrintSamplerData(worldIndexes.back());
			PrintDistances(worldIndexes.back());
			PrintDihedralsQs(worldIndexes.back());
			fprintf(logFile, "\n");
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


// Pass compounds to the new world
void Context::passTopologiesToNewWorld(int newWorldIx)
{
	for(std::size_t molIx = 0; molIx < nofMols; molIx++){
		topologies[molIx].setMultibodySystem(
			*((updWorld(newWorldIx))->compoundSystem) );

		// Reset mobilized body indeces in Compound
		for(std::size_t k = 0; k < topologies[molIx].getNumAtoms(); k++){
			SimTK::Compound::AtomIndex aIx =
				(topologies[molIx].bAtomList[k]).getCompoundAtomIndex();
		//	//SimTK::MobilizedBodyIndex mbx =
		//	//	worlds[newWorldIx].getMobilizedBodyIndex(aIx);
			SimTK::MobilizedBodyIndex mbx =
				topologies[molIx].getAtomMobilizedBodyIndexFromMap(aIx, newWorldIx);
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
    std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>> referenceAtomsLocationsFromFile;

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

	// Done
	nofThermodynamicStates++;
}

// Get the number of replicas
const size_t& Context::getNofReplicas(void) const
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
const size_t& Context::getNofThermodynamicStates(void) const
{
	return nofThermodynamicStates;
}

// Set the intial mapping between replicas and thermoStates
void Context::loadReplica2ThermoIxs(void)
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

void Context::PrintReplicaMaps(void){

	std::cout << "Replica -> Thermo:\n";
	for(const auto& elem : replica2ThermoIxs){
		std::cout << elem.first << " " << elem.second << "\n";
	}

	std::cout << "Thermo -> Replica:\n";
	for(const auto& elem : thermo2ReplicaIxs){
		std::cout << elem.first << " " << elem.second << "\n";
	}

}

// Attempt swap between replicas r_i and r_j
// Code inspired from OpenmmTools
// Chodera JD and Shirts MR. Replica exchange and expanded ensemble simulations
// as Gibbs multistate: Simple improvements for enhanced mixing. J. Chem. Phys.
//, 135:194110, 2011. DOI:10.1063/1.3660669
bool Context::attemptSwap(int replica_i, int replica_j)
{

	bool returnValue = false;

	// Get replicas' thermodynamic states indexes
	int thermoState_i = replica2ThermoIxs[replica_i];
	int thermoState_j = replica2ThermoIxs[replica_j];

	// Record this attempt
	nofAttemptedSwapsMatrix[thermoState_i][thermoState_j] += 1;
	nofAttemptedSwapsMatrix[thermoState_j][thermoState_i] += 1;

	// Replica i reduced potential of state i
	SimTK::Real Eii = thermodynamicStates[thermoState_i].getBeta()
		* replicas[replica_i].getPotentialEnergy();

	// Replica j reduced potential of state j
	SimTK::Real Ejj = thermodynamicStates[thermoState_j].getBeta()
		* replicas[replica_j].getPotentialEnergy();

	// Replica i reduced potential of state j
	SimTK::Real Eij = thermodynamicStates[thermoState_j].getBeta()
		* replicas[replica_i].getPotentialEnergy();

	// Replica j reduced potential of state i
	SimTK::Real Eji = thermodynamicStates[thermoState_i].getBeta()
		* replicas[replica_j].getPotentialEnergy();

	SimTK::Real log_p_accept = -1.0 * (Eij + Eji) + Eii + Ejj;

	SimTK::Real unifSample = uniformRealDistribution(randomEngine);

	//std::cout << "Replicas energies = "
	//	<< replicas[replica_i].getPotentialEnergy() << " "
	//	<< replicas[replica_j].getPotentialEnergy() << " "
	//	<< std::endl;
	//std::cout << "log_p_accept = " << log_p_accept << std::endl;

	// Acceptance criterion
	if((log_p_accept >= 0.0) || (unifSample < std::exp(log_p_accept))){

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

		// Exchange potential energies (not necessary)
		SimTK::Real tempE = replicas[replica_i].getPotentialEnergy();
		replicas[replica_i].setPotentialEnergy(replicas[replica_j].getPotentialEnergy());
		replicas[replica_j].setPotentialEnergy(tempE);

		//std::cout << "Swap done." << std::endl;

		returnValue = true;

		//PrintReplicaMaps();

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
		attemptSwap(replica_i, replica_j);
	}
}

// Mix neighboring replicas
void Context::mixNeighboringReplicas(void)
{
	// Go through neighboring thermodynamic states
	for(size_t thermoState_k = 0;
	thermoState_k < nofThermodynamicStates - 1;
	thermoState_k++){
		// Get thermodynamic states
		int thermoState_i = thermoState_k;
		int thermoState_j = thermoState_k + 1;

		// Get replicas corresponding to the thermodynamic states
		int replica_i = thermo2ReplicaIxs[thermoState_i];
		int replica_j = thermo2ReplicaIxs[thermoState_j];

		//std::cout << "Attempt to swap replicas " << replica_i
		//	<< " and " << replica_j << std::endl;

		// Attempt to swap
		attemptSwap(replica_i, replica_j);
	}
}

// Mix replicas
void Context::mixReplicas(void)
{
	if(replicaMixingScheme == ReplicaMixingScheme::neighboring){
		mixNeighboringReplicas();
	}else{
		mixAllReplicas(nofReplicas*nofReplicas*nofReplicas);
	}
}


// Load replica's atomLocations into it's front world
void Context::restoreReplicaCoordinatesToFront(int whichReplica)
{

	//std::cout <<  "Context::restoreReplicaCoordinatesToFront " << whichReplica << ": " << std::flush;

	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[whichReplica];

	//std::cout <<  "Context::restoreReplicaCoordinatesToFront thermoIx " << thermoIx << std::endl << std::flush;

	// Get worlds indexes of this thermodynamic state
	std::vector<int> worldIndexes =
		thermodynamicStates[thermoIx].getWorldIndexes();

	//std::cout <<  " worldIndexes[1] " << worldIndexes[1] << std::flush;

	// Set thoermoState front world from replica coordinate buffer
	// Will use worlds own integrator's advanced state
	SimTK::State& state =
		(worlds[worldIndexes.front()].integ)->updAdvancedState();

	worlds[worldIndexes.front()].setAtomsLocationsInGround(state,
		replicas[whichReplica].getAtomsLocationsInGround());

	//std::cout << "Context::restoreReplicaCoordinatesToFront worldIndexes.front() " << worldIndexes.front() << std::endl << std::flush;

}

// Stores replica's front world's coordinates into it's atomsLocations
// This should always be a fully flexible world
void Context::storeReplicaCoordinatesFromFront(int whichReplica)
{

	//std::cout <<  "storeReplicaCoordinatesFromFront " << whichReplica << ": ";

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


// Get ennergy of the back world and store it in replica thisReplica
void Context::storeReplicaEnergyFromBack(int replicaIx)
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

// Get ennergy of the back world and store it in replica thisReplica
void Context::storeReplicaEnergyFromFrontFull(int replicaIx)
{
	// Get thermoState corresponding to this replica
	// KEYWORD = replica, VALUE = thermoState
	int thermoIx = replica2ThermoIxs[replicaIx];

	// Get the index of the back world
	int frontWorldIx =
		thermodynamicStates[thermoIx].getWorldIndexes().front();

	// Get the back world energy
	SimTK::Real energy =
		worlds[frontWorldIx].CalcFullPotentialEnergyIncludingRigidBodies();

	// Set this replica's energy
	replicas[replicaIx].setPotentialEnergy(energy);
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

	std::cout << "iniTemperature set to " << T << std::endl << std::flush;

	// -------------
	// Set samplers parameters for this replica
	// =============
	std::vector<SimTK::Real> replicaTimesteps = thermodynamicStates[thisThermoStateIx].getTimesteps();
	std::vector<int> replicaMdsteps = thermodynamicStates[thisThermoStateIx].getMdsteps();

	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setTimestep(
			replicaTimesteps[i]);
	}
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setMDStepsPerSample(
			replicaMdsteps[i]);
	}


	std::cout << "iniTimesteps set to ";
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		std::cout << worlds[replicaWorldIxs[i]].getSampler(0)->getTimestep() << " " ;
	}
	std::cout << std::endl;
	std::cout << "iniMdsteps set to ";
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		std::cout << worlds[replicaWorldIxs[i]].getSampler(0)->getMDStepsPerSample() << " " ;
	}
	std::cout << std::endl;
	// =============


	for(size_t ri = 0; ri < 1; ri++){
		for(std::size_t worldIx = 0; worldIx < replicaNofWorlds; worldIx++){


			int frontIx = -1;
			int backIx = -1;

			// -------------
			// SAMPLE from the current world
			frontIx = replicaWorldIxs.front();
			//std::cout << "Sample world " << front << "\n";
			int accepted = worlds[frontIx].generateSamples(
				0,
				NMAOption[frontIx]);

			// =============

			// -------------
			// ROTATE
			///std::cout << "Rotate from";
			///for(int k = 0; k < replicaNofWorlds; k++){std::cout << " " << replicaWorldIxs[k];}

			// Rotate worlds indices (translate from right to left)
		   	std::rotate(replicaWorldIxs.begin(),
				replicaWorldIxs.begin() + 1,
				replicaWorldIxs.end());

			std::cout << " to";
			for(int k = 0; k < replicaNofWorlds; k++){std::cout << " " << replicaWorldIxs[k];}
			//std::cout << "\n";
			// =============


			// -------------
			// TRANSFER coordinates from last world to current
			// TODO: eliminate in the last iteration
			frontIx = replicaWorldIxs.front();
			backIx = replicaWorldIxs.back();

			if(replicaNofWorlds > 1) {
				//std::cout << "Transfer from world " << backIx
				//	<< " to " << frontIx << std::endl;

				transferCoordinates(backIx, frontIx);
			}
			// =============

		} // END iteration through worlds
	} // END iteration through rounds

}

// Go through all of this replica's worlds and generate samples
void Context::RunReplica(int thisReplica, int howManyRounds)
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

	std::cout << "Temperature set to " << T << std::endl << std::flush;

	// -------------
	// Set samplers parameters for this replica
	// =============
	std::vector<SimTK::Real> replicaTimesteps = thermodynamicStates[thisThermoStateIx].getTimesteps();
	std::vector<int> replicaMdsteps = thermodynamicStates[thisThermoStateIx].getMdsteps();

	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setTimestep(
			replicaTimesteps[i]);
	}
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		worlds[replicaWorldIxs[i]].updSampler(0)->setMDStepsPerSample(
			replicaMdsteps[i]);
	}

	std::cout << "Timesteps set to ";
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		std::cout << worlds[replicaWorldIxs[i]].getSampler(0)->getTimestep() << " " ;
	}
	std::cout << std::endl;
	std::cout << "Mdsteps set to ";
	for(std::size_t i = 0; i < replicaNofWorlds; i++){
		std::cout << worlds[replicaWorldIxs[i]].getSampler(0)->getMDStepsPerSample() << " " ;
	}
	std::cout << std::endl;
	// =============


	for(size_t ri = 0; ri < howManyRounds; ri++){
		for(std::size_t worldIx = 0; worldIx < replicaNofWorlds; worldIx++){

			int frontIx = -1;
			int backIx = -1;

			// -------------
			// SAMPLE from the current world
			frontIx = replicaWorldIxs.front();
			//std::cout << "Sample world " << front << "\n";
			int accepted = worlds[frontIx].generateSamples(
				nofSamplesPerRound[frontIx],
				NMAOption[frontIx]);

			// =============

			// -------------
			// ROTATE
			///*print*/std::cout << "Rotate from";/*print*/
			///*print*/for(int k = 0; k < replicaNofWorlds; k++){std::cout << " " << replicaWorldIxs[k];}/*print*/

			// Rotate worlds indices (translate from right to left)
		   	std::rotate(replicaWorldIxs.begin(),
				replicaWorldIxs.begin() + 1,
				replicaWorldIxs.end());

			///*print*/std::cout << " to";/*print*/
			///*print*/for(int k = 0; k < replicaNofWorlds; k++){std::cout << " " << replicaWorldIxs[k];}/*print*/
			///*print*/std::cout << "\n";/*print*/
			// =============


			// -------------
			// TRANSFER coordinates from last world to current
			// TODO: eliminate in the last iteration
			frontIx = replicaWorldIxs.front();
			backIx = replicaWorldIxs.back();

			if(replicaNofWorlds > 1) {
				//std::cout << "Transfer from world " << backIx
				//	<< " to " << frontIx << std::endl;

				transferCoordinates(backIx, frontIx);
			}
			// =============

		} // END iteration through worlds
	} // END iteration through rounds

}

// Run replica exchange protocol
void Context::RunREX(void)
{

	// Is this necesary
	realizeTopology();

	// Allocate space for swap matrices
	allocateSwapMatrices();

	// Run each replica one time initially
	std::cout << " REX batch " << 0 << std::endl;
	for (size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++){
		std::cout << "REX replica " << replicaIx << std::endl;

            //std::cout << "Replica front world coordinates:\n";
			//replicas[replicaIx].PrintCoordinates();
			initializeReplica(replicaIx); // EXPERIMENTAL
			restoreReplicaCoordinatesToFront(replicaIx); // EXPERIMENTAL

			// Iterate this replica's worlds
			RunReplica(replicaIx, swapEvery);

			// This should always be a fully flexible world
			storeReplicaCoordinatesFromFront(replicaIx);

			// Store energy
			storeReplicaEnergyFromFrontFull(replicaIx);

			PrintToLog(worldIndexes.front());

			writePdbs(0,	thermo2ReplicaIxs[replicaIx]);
	}

	// Mix replicas
	mixReplicas();

	PrintNofAcceptedSwapsMatrix();

	// Main loop
	int nofMixes = int(nofRounds / swapEvery);
	for(size_t mixi = 1; mixi < nofMixes; mixi++){
		std::cout << " REX batch " << mixi << std::endl;

		// Run each replica serially
		for (size_t replicaIx = 0; replicaIx < nofReplicas; replicaIx++){
			std::cout << "REX replica " << replicaIx << std::endl;

			//std::cout << "Replica front world coordinates:\n";
			//replicas[replicaIx].PrintCoordinates();
			restoreReplicaCoordinatesToFront(replicaIx);

			// Iterate this replica's worlds
			RunReplica(replicaIx, swapEvery);

			// This should always be a fully flexible world
			storeReplicaCoordinatesFromFront(replicaIx);

			// STORE energy
			storeReplicaEnergyFromFrontFull(replicaIx);

			// Write energy and geometric features to logfile
			if( !(mixi % getPrintFreq()) ){
				PrintToLog(worldIndexes.front());
			}

			// Write pdb
			if( pdbRestartFreq != 0){
				int thermoStateIx = thermo2ReplicaIxs[replicaIx];
				if((mixi % pdbRestartFreq) == 0){
					writePdbs(mixi,
					thermoStateIx);
				}
			}

		} // end replicas simulations

		// Mix replicas
		mixReplicas();

		PrintNofAcceptedSwapsMatrix();

	} // end rounds

	//PrintNofAttemptedSwapsMatrix();
	PrintNofAcceptedSwapsMatrix();
	//PrintReplicaMaps();


}

// Print info about all the replicas and thermo states
void Context::PrintReplicas(void)
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

void Context::PrintNofAcceptedSwapsMatrix(void){

	size_t M = nofAcceptedSwapsMatrix.size();

	std::cout << "Number of accepted swaps matrix:\n";
	for(size_t i = 0; i < M; i++){
		for(size_t j = 0; j < M; j++){
			std::cout << nofAcceptedSwapsMatrix[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void Context::PrintNofAttemptedSwapsMatrix(void){

	size_t M = nofAttemptedSwapsMatrix.size();

	std::cout << "Number of attempted swaps matrix:\n";
	for(size_t i = 0; i < M; i++){
		for(size_t j = 0; j < M; j++){
			std::cout << nofAttemptedSwapsMatrix[i][j] << " ";
		}
		std::cout << "\n";
	}
}

const int Context::getSwapEvery(void){
	return swapEvery;
}

void Context::setSwapEvery(const int& n){
	swapEvery = n;
}

void Context::writePdbs(int someIndex, int thermodynamicStateIx)
{

	// Update bAtomList in Topology
	const SimTK::State& pdbState =
	worlds[worldIndexes.front()].integ->updAdvancedState();
	worlds[worldIndexes.front()].updateAtomListsFromCompound(pdbState);

	// Write
	for(int mol_i = 0; mol_i < getNofMolecules(); mol_i++){
		topologies[mol_i].writeAtomListPdb(
			outputDir,
			"/pdbs/sb."
			    + pdbPrefix + "." + std::to_string(mol_i) + "."
			    + "s" + std::to_string(thermodynamicStateIx) + ".",
			".pdb",
			10,
			someIndex);
	}
}

void Context::randomizeWorldIndexes(void)
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
	std::vector<std::vector<std::pair<
		bSpecificAtom *, SimTK::Vec3> > >
		otherWorldsAtomsLocations =
	updWorld(src)->getAtomsLocationsInGround(lastAdvancedState);

	// Pass compounds to the new world
	passTopologiesToNewWorld(dest);

	// Set coordinates to destination (reconstruct)
	currentAdvancedState = updWorld(dest)->setAtomsLocationsInGround(
		currentAdvancedState, otherWorldsAtomsLocations);
}

// Go through all the worlds and generate samples
void Context::RunOneRound(void)
{

	for(std::size_t worldIx = 0; worldIx < getNofWorlds(); worldIx++){

		// Rotate worlds indices (translate from right to left)
		if(isWorldsOrderRandom){
			randomizeWorldIndexes();
		}

		// Rotate the vector of world indexes
	   	std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());

		// Get indeces
		int currentWorldIx = worldIndexes.front();
		int lastWorldIx = worldIndexes.back();

		// Transfer coordinates from last world to current
		if(worldIndexes.size() > 1) {

			std::cout << "Transfer from world " << lastWorldIx
				<< " to " << currentWorldIx << std::endl;

			transferCoordinates(lastWorldIx, currentWorldIx);
		}

		// Generate samples from the current world
		int accepted = worlds[currentWorldIx].generateSamples(
			nofSamplesPerRound[currentWorldIx],
			NMAOption[currentWorldIx]);

	} // END iteration through worlds
}


// Normal run
void Context::Run(int, SimTK::Real Ti, SimTK::Real Tf)
{

	// Initialize world indeces
	std::size_t currentWorldIx = worldIndexes.front();
	std::size_t lastWorldIx = 0;

	// Write an initial pdb
    topologies[0].writeAtomListPdb(
        outputDir,
        "/pdbs/ini.",
        ".pdb",
        10,
        0);

    writeInitialPdb();

	if( std::abs(Tf - Ti) < SimTK::TinyReal){ // Don't heat

		// Main loop: iterate through rounds
		for(int round = 0; round < nofRounds; round++){

			RunOneRound();

			// Write energy and geometric features to logfile
			if( !(round % getPrintFreq()) ){
				PrintToLog(worldIndexes.back());
				PrintToLog(worldIndexes.front());
			}

			// Write pdb
			if( pdbRestartFreq != 0){
				if((round % pdbRestartFreq) == 0){
					writePdbs(round);
				}
			}

		}

	}else{// if Ti != Tf heating protocol
		SimTK::Real Tincr = (Tf - Ti) / static_cast<SimTK::Real>(nofRounds);
		SimTK::Real currT = Ti;
		for(int round = 0; round < nofRounds; round++){ // Iterate rounds

			// Set current temperature
			currT += Tincr;
			setTemperature(currT);
			std::cout << "T= " << currT << std::endl;

			RunOneRound();

			// Write energy and geometric features to logfile
			if( !(round % getPrintFreq()) ){
				PrintToLog(worldIndexes.back());
				PrintToLog(worldIndexes.front());
			}

			// Write pdb
			if( pdbRestartFreq != 0){
				if((round % pdbRestartFreq) == 0){
					writePdbs(round);
				}
			}

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
void Context::addDistance(std::size_t whichWorld, std::size_t whichCompound, int aIx1, int aIx2)
{
	// TODO some are 64 bit, some 32 bit. What to do?
	std::vector<int> tempV;
	tempV.push_back(static_cast<int>(whichWorld));
	tempV.push_back(static_cast<int>(whichCompound));
	tempV.push_back(aIx1);
	tempV.push_back(aIx2);

	distanceIxs.push_back(tempV);
}

// Get distances
void Context::addDistances(std::vector<int> distanceIx)
{
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for(unsigned int ai = 0; ai < distanceIx.size() / 2; ai++){
			addDistance(worldIx, 0, distanceIx[2*ai + 0], distanceIx[2*ai + 1]);
		}
	}
}

void Context::addDihedral(std::size_t whichWorld, std::size_t whichCompound, int aIx1, int aIx2, int aIx3, int aIx4)
{
	// TODO some are 64 bit, some 32 bit. What to do?
	std::vector<int> tempV;
	tempV.push_back(static_cast<int>(whichWorld));
	tempV.push_back(static_cast<int>(whichCompound));
	tempV.push_back(aIx1);
	tempV.push_back(aIx2);
	tempV.push_back(aIx3);
	tempV.push_back(aIx4);

	dihedralIxs.push_back(tempV);
}

// Get dihedrals. TODO : only adds to the first Topology
void Context::addDihedrals(std::vector<int> dihedralIx){
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
void Context::PrintSamplerData(std::size_t whichWorld)
{

	SimTK::State& currentAdvancedState = worlds[whichWorld].integ->updAdvancedState();
/*    std::cout << currentAdvancedState.getNU() << ' '
		<< worlds[whichWorld].updSampler(0)->getAcceptedSteps() << ' '
		<< std::setprecision(4) << std::fixed
		<< worlds[whichWorld].updSampler(0)->getOldPE() << ' '
		<< worlds[whichWorld].updSampler(0)->getSetPE() << ' '
		<< worlds[whichWorld].updSampler(0)->getLastAcceptedKE() << ' '
		<< worlds[whichWorld].updSampler(0)->getProposedKE() << ' '
		<< worlds[whichWorld].updSampler(0)->getOldFixman() << ' '
		<< worlds[whichWorld].updSampler(0)->getSetFixman() << ' '
		<< worlds[whichWorld].updSampler(0)->getProposedFixman() << ' '
		;
*/
/*    // Use printf for faster output
	printf("%d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f "
		, currentAdvancedState.getNU()
		, worlds[whichWorld].updSampler(0)->getAcceptedSteps()
		, worlds[whichWorld].updSampler(0)->getOldPE()
		, worlds[whichWorld].updSampler(0)->getSetPE()
		, worlds[whichWorld].updSampler(0)->getLastAcceptedKE()
		, worlds[whichWorld].updSampler(0)->getProposedKE()
		, worlds[whichWorld].updSampler(0)->getOldFixman()
		, worlds[whichWorld].updSampler(0)->getSetFixman()
		, worlds[whichWorld].updSampler(0)->getProposedFixman()
	);
*/

/*    // Avoid get function calls
	printf("%d %d %.2f %.2f %.2f %.2f %.2f %.2f"
		, currentAdvancedState.getNU()
		, (worlds[whichWorld].samplers[0])->acceptedSteps
		, (worlds[whichWorld].samplers[0])->pe_o
		, (worlds[whichWorld].samplers[0])->pe_set
		, (worlds[whichWorld].samplers[0])->ke_proposed
		, (worlds[whichWorld].samplers[0])->ke_n
		, (worlds[whichWorld].samplers[0])->fix_o
		, (worlds[whichWorld].samplers[0])->fix_set
	);
*/

	// Write to a file instead of stdout
	fprintf(logFile, "%d %d %.2f %.2f %.2f %.2f %.12f %.12f %.12f "
		, currentAdvancedState.getNU()
		, pHMC(worlds[whichWorld].samplers[0])->acceptedSteps
		, pHMC((worlds[whichWorld].samplers[0]))->pe_o
		, pHMC((worlds[whichWorld].samplers[0]))->pe_set
		, pHMC((worlds[whichWorld].samplers[0]))->ke_proposed
		, pHMC((worlds[whichWorld].samplers[0]))->ke_n
		, pHMC((worlds[whichWorld].samplers[0]))->fix_o
		, pHMC((worlds[whichWorld].samplers[0]))->fix_n
		, pHMC((worlds[whichWorld].samplers[0]))->fix_set
	);
	fflush(logFile);

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

void Context::PrintGeometry(std::size_t whichWorld)
{
	// TODO same 32 vs 64 bit thing. See the many other function below. Might use a vector<struct>, not a vector<vector>
	for (const auto& distanceIx : distanceIxs) {
		if( distanceIx[0] == static_cast<int>(whichWorld)){
			fprintf(logFile, "%.3f ",
				Distance(distanceIx[0], distanceIx[1], 0, distanceIx[2], distanceIx[3]));
		}
	}

	for (const auto& dihedralIx : dihedralIxs){
		if( dihedralIx[0] == static_cast<int>(whichWorld)){
			fprintf(logFile, "%.3f ", Dihedral(dihedralIx[0], dihedralIx[1], 0, dihedralIx[2], dihedralIx[3], dihedralIx[4], dihedralIx[5]));
		}
	}
}

void Context::PrintDistances(std::size_t whichWorld)
{
	for (const auto& distanceIx : distanceIxs) {
		if( distanceIx[0] == static_cast<int>(whichWorld)){
			fprintf(logFile, "%.3f ", Distance(distanceIx[0], distanceIx[1], 0, distanceIx[2], distanceIx[3]) );
		}
	}
}

void Context::PrintDihedrals(std::size_t whichWorld)
{
	for (const auto& dihedralIx : dihedralIxs){
		if( dihedralIx[0] == static_cast<int>(whichWorld)){
			fprintf(logFile, "%.3f ", Dihedral(dihedralIx[0], dihedralIx[1], 0, dihedralIx[2], dihedralIx[3], dihedralIx[4], dihedralIx[5]) );
		}
	}
}

void Context::PrintDihedralsQs(std::size_t whichWorld)
{
	for (const auto& dihedralIx : dihedralIxs){
		if( dihedralIx[0] == static_cast<int>(whichWorld)){

			// std::cout << "Context::PrintDihedralsQs w c s a1 a2 a3 a4 ="
			//     << " " << dihedralIx[0] << " " << dihedralIx[1] << " " << 0
			//     << " " << (worlds[whichWorld].getTopology(0)).getAtomName( SimTK::Compound::AtomIndex(dihedralIx[2]))
			//     << " " << (worlds[whichWorld].getTopology(0)).getAtomName( SimTK::Compound::AtomIndex(dihedralIx[3]))
			//     << " " << (worlds[whichWorld].getTopology(0)).getAtomName( SimTK::Compound::AtomIndex(dihedralIx[4]))
			//     << " " << (worlds[whichWorld].getTopology(0)).getAtomName( SimTK::Compound::AtomIndex(dihedralIx[5]))
			//     << std::endl;

			fprintf(logFile, "%.3f ", Dihedral(
				dihedralIx[0], dihedralIx[1], 0,
				dihedralIx[2], dihedralIx[3],
				dihedralIx[4], dihedralIx[5]) );

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
		if( distanceIx[0] == static_cast<int>(whichWorld)){

			fprintf(logFile, "%.3f ",
				Distance(distanceIx[0], distanceIx[1], 0,
				   distanceIx[2], distanceIx[3]) );

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
			fprintf(logFile, "%.3f ",
				((-1 * X_FM1.p()) + deltaX_PF.p() + X_FM2.p()).norm() );

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

void Context::PrintToLog(int worldIx)
{
	PrintSamplerData(worldIx);
	PrintDistances(worldIx);
	PrintDihedralsQs(worldIx);
	fprintf(logFile, "\n");
}

// Write intial pdb for reference
// TODO: what's the deal with mc_step
void Context::writeInitialPdb(void)
{

	// - we need this to get compound atoms
	int currentWorldIx = worldIndexes.front();
	SimTK::State& advancedState = (updWorld(currentWorldIx))->integ->updAdvancedState();

	constexpr int mc_step = -1;

	// Pass compounds to the new world
	passTopologiesToNewWorld(currentWorldIx);

	(updWorld(currentWorldIx))->updateAtomListsFromCompound(advancedState);
	std::cout << "Writing pdb initial" << mc_step << ".pdb" << std::endl;
	for(unsigned int mol_i = 0; mol_i < topologies.size(); mol_i++){
		topologies[mol_i].writeAtomListPdb(getOutputDir(),
		"/pdbs/sb." + getPdbPrefix() + ".", ".pdb", 10, mc_step);
	}
}

// Write final pdb for reference
void Context::writeFinalPdb(void)
{
	for(unsigned int mol_i = 0; mol_i < nofMols; mol_i++){
		topologies[mol_i].writeAtomListPdb(
			getOutputDir(), "/pdbs/final." + getPdbPrefix() + ".", ".pdb", 10,
			getNofRounds());
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

SimTK::Real Context::Distance(std::size_t whichWorld, std::size_t whichCompound, std::size_t, int a1, int a2)
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


///////////////////////////
// CTYPES
///////////////////////////

//
//extern "C"
//{
//	Context* makeContext
//}

//------------



