#pragma once

#include "Robo.hpp"
#include "Sampler.hpp"
#include "SetupReader.hpp"
#include "World.hpp"
#include "ThermodynamicState.hpp"
#include "Replica.hpp"
#include "SetupReader.hpp"

#include <fstream>

class Sampler;
class World;

enum class ReplicaMixingScheme : int {
	all = 0,
	neighboring = 1
};

enum class RunType : int {
	Default = 0,
	SimulatedTempering = 1,
	REX = 2,
	RENS = 3
};

class Context{
public:
	bool initializeFromFile(const std::string& file);
	void run();
	void run(int steps);

	// Input functions
	bool loadTopologyFile(std::string topologyFilename);
	bool loadCoordinatesFile(std::string coordinatesFilename);
	void PrintCoordinates(const std::vector<std::vector
        <std::pair <bSpecificAtom *,
		SimTK::Vec3>>>& atomsLocations);

	bool loadRigidBodiesSpecs(std::size_t whichWorld, int whichMolecule, std::string RBSpecsFN);
	bool loadFlexibleBondsSpecs(std::size_t whichWorld, int whichMolecule, std::string FlexSpecsFN);
	void setRegimen (std::size_t whichWorld, int whichMolecule, std::string regimen);

	/** Load molecules based on loaded filenames. One molecule
	creates a topology object and based on amberReader forcefield
	 parameters - defines Biotypes; - adds BAT parameters to DuMM **/
	void AddMolecules(
		int requestedNofMols,
		SetupReader& setupReader
		//std::vector<std::string> argRoots,
		//std::vector<std::string> argRootMobilities
	);

	/** It calls DuMMs defineAtomClass. These Molmodel functions contain
	information regarding the force field parameters. **/
	void updDummAtomClasses(
		std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
		, int worldIx		
	);

	// Loads parameters into DuMM
	void addDummParams(
		int requestedNofMols,
		SetupReader& setupReader
	);

	// Adopts compound by the CompoundSystem
	// and loads maps of indexes
	void model(
		int requestedNofMols,
		SetupReader& setupReader
	);

	void modelOneEmbeddedTopology(int whichTopology, int whichWorld, std::string rootMobilizer);
	void modelTopologies(std::vector<std::string> GroundToCompoundMobilizerTypes);

	// Add task spaces
	void addTaskSpacesLS(void);

	/** Add constraints */
	void addConstraints(void);

	void realizeTopology();
	void realizePosition();

	void LoadWorldsFromSetup(SetupReader&);

	void passTopologiesToNewWorld(int newWorldIx);

	int getNofMolecules();
	//------------

	// --- Thermodynamics ---
	// Set mixing rule for Lennard-Jones
	void setVdwMixingRule(SimTK::DuMMForceFieldSubsystem::VdwMixingRule mixingRule);

	// Get/set the main temperature (acc/rej temperature for MC)
	SimTK::Real getTemperature(std::size_t whichWorld) const;
	void  setTemperature(std::size_t whichWorld, float someTemperature);
	void setTemperature(float someTemperature);

	// If HMC, get/set the guidance Hamiltonian temperature
	SimTK::Real getGuidanceTemperature(std::size_t whichWorld, std::size_t whichSampler);
	void  setGuidanceTemperature(std::size_t whichWorld, std::size_t whichSampler, SimTK::Real someTemperature);
	//------------

	// --- Simulation parameters ---

	// Add sampler
	BaseSampler* addSampler(std::size_t whichWorld,
		std::string samplerName, std::string integratorName);

	BaseSampler* addSampler(std::size_t whichWorld,
		SamplerName whichSampler, IntegratorName whichIntegrator);

	// Samplers have to set parameters after Simbody subsystems generation
	void initializeSampler(std::size_t whichWorld, std::size_t whichSampler);

	// Amber like scale factors.
	void setAmberForceFieldScaleFactors(std::size_t whichWorld);
	void setAmberForceFieldScaleFactors(void);

	// Set a global scaling factor for the forcefield
	void setGlobalForceFieldScaleFactor(std::size_t whichWorld, SimTK::Real);
	void setGlobalForceFieldScaleFactor(SimTK::Real);

	// Set GBSA implicit solvent scale factor
	void setGbsaGlobalScaleFactor(std::size_t whichWorld, SimTK::Real);
	void setGbsaGlobalScaleFactor(SimTK::Real);

	// If HMC, get/set the number of MD steps
	int getNofMDStepsPerSample(std::size_t whichWorld, std::size_t whichSampler);
	void setNofMDStepsPerSample(std::size_t whichWorld, std::size_t whichSampler, int MDStepsPerSample);

	// If HMC, get/set timestep forMD
	SimTK::Real getTimestep(std::size_t whichWorld, std::size_t whichSampler) const;
	void setTimestep(std::size_t whichWorld, std::size_t whichSampler, SimTK::Real timeStep);

	// Use Fixman torque as an additional force subsystem
	void useFixmanPotential(std::size_t whichWorld, std::size_t whichSampler);
	bool isUsingFixmanPotential(std::size_t whichWorld, std::size_t whichSampler);

	void addFixmanTorque(std::size_t whichWorld);
	bool isUsingFixmanTorque(std::size_t whichWorld) const;

	void setFixmanTorqueScaleFactor(std::size_t whichWorld, SimTK::Real scaleFactor);
	void setFixmanTorqueTemperature(std::size_t whichWorld, SimTK::Real temperature);
	//------------

	// --- Mixing parameters ---
	// Another way to do it is setting the number of rounds
	int getRequiredNofRounds();
	void setRequiredNofRounds(int argNofRounds);

	int getNofRoundsTillReblock();
	void setNofRoundsTillReblock(int nofRoundsTillReblock);
	void updNofRoundsTillReblock(int nofRoundsTillReblock);

	SimTK::Real getNofSamplesPerRound(std::size_t whichWorld);
	void setNofSamplesPerRound(std::size_t whichWorld, SimTK::Real MCStepsPerRound);

	std::size_t getWorldIndex(std::size_t which) const;

	// Adaptive Gibbs blocking: TODO: consider moving in World
	void allocateReblockQsCache(void);
	void allocateReblockQsCacheQVectors(void);

	// --- Arrange different mixing parameters ---
	void initializeMixingParamters();
	//------------

	void addEmptyWorlds(std::size_t NofWorlds, std::vector<double> visualizerFrequencies);
	World * addWorld(bool visual, SimTK::Real visualizerFrequency = 0.0015);
	//World * AddWorld(World *, bool visual);

	// Load/store Mobilized bodies joint types in samplers
	void loadMbxsToMobilities(void);

	World * getWorld();
	World * getWorld(std::size_t which);

	World * updWorld();
	World * updWorld(std::size_t which);

	// Returns the size of the worlds vector
	std::size_t getNofWorlds() const;

	SimTK::DuMMForceFieldSubsystem * updForceField(std::size_t whichWorld);

	// Writeble reference to a samplers advanced state
	SimTK::State& updAdvancedState(std::size_t whichWorld, std::size_t whichSampler);

	void RotateWorlds();
	//------------

	// --- Main ---
	void randomizeWorldIndexes(void);
	void transferCoordinates(int src, int dest);

	// Go through all the worlds and generate samples
	void RunOneRound(void);
	void Run(int howManyRounds, SimTK::Real Ti, SimTK::Real Tf);
	void RunSimulatedTempering(int howManyRounds, SimTK::Real Ti, SimTK::Real Tf);
	void setNofBoostStairs(std::size_t whichWorld, int howManyStairs);
	int getNofBoostStairs(std::size_t whichWorld);
	void setNumThreadsRequested(std::size_t which, int howMany);
	void setUseOpenMMAcceleration(bool arg);
	void setUseOpenMMIntegration(std::size_t which, Real temperature, Real stepsize);
	void setUseOpenMMCalcOnlyNonBonded(bool arg);
	void setNonbondedMethod(std::size_t whichWorld, int methodInx);
	void setNonbondedCutoff(std::size_t whichWorld, Real cutoffNm);

	SimTK::Real Pearson(std::vector<std::vector<SimTK::Real>> someVector,
		int QIx1, int QIx2); // 2D roundsTillReblock; 3D nofQs

	/** Print the number of threads each World got **/
	void PrintNumThreads();

	/** Get/Set seed for reproducibility. **/
	void setSeed(std::size_t whichWorld, std::size_t whichSampler, uint32_t seed);
	uint32_t getSeed(std::size_t whichWorld, std::size_t whichSampler) const;

	//------------

	/** Analysis related functions **/
	void addDistance(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t aIx1, std::size_t aIx2);
	void addAngle(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t aIx1, std::size_t aIx2, std::size_t aIx3);
	void addDihedral(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t aIx1, std::size_t aIx2, std::size_t aIx3, std::size_t aIx4);

	void addDistances(const std::vector<std::size_t>& distanceIx);
	void addAngles(const std::vector<std::size_t>& angleIx);
	void addDihedrals(const std::vector<std::size_t>& dihedralIx);

	// Allocate space for containers that keep statistics if we're doing any
	void allocWorldsStatsContainers(void);

	// --- Output ---
	void printThermodynamics(void);
	void printStatus(void);

	// Print Molmodel related information
	void PrintMolmodelAndDuMMTypes(void);

	// Print DuMM atoms stations in mobilized body frame
	void checkAtomStationsThroughDumm(void);

	// Print Simbody related information
	void PrintSimbodyMobods(void);
	void PrintGeometry(SetupReader&, std::size_t whichWorld);
	void PrintFreeE2EDist(std::size_t whichWorld, int whichCompound);

	void PrintGeometryToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintDistancesToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintAnglesToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintDihedralsToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintDihedralsQsToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintSamplerDataToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintToLog(std::size_t whichReplica,
		std::size_t whichWorld, std::size_t whichSampler);

	// Write intial/final pdb for reference
	void writeInitialPdb(void);
	void writeFinalPdb(void);

	//
	void WritePdb(std::size_t whichWorld);

	//
	void writePdbs(int someIndex, int thermodynamicStateIx = 0);

	SimTK::Real Dihedral(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t whichSampler, int a1, int a2, int a3, int a4);
	SimTK::Real Roboangle(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t whichSampler, int a1, int a2, int a3);
	SimTK::Real Distance(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t whichSampler, int a1, int a2);

	int getPdbRestartFreq();
	void setPdbRestartFreq(int argFreq);
	int getPrintFreq();
	void setPrintFreq(int argFreq);

	std::string getOutputDir();
	void setOutputDir(std::string arg);
	std::string getPdbPrefix();
	void setPdbPrefix(std::string arg);
	//------------

	//////////////////////////////////
	//// REPLICA EXCHANGE FUNCTIONS //
	//////////////////////////////////

	void setNofReplicas(const size_t& argNofReplicas);
	const size_t& getNofReplicas(void) const;
	//void setNofThermodynamicStates(const size_t& argNofThermodynamicStates);
	const size_t& getNofThermodynamicStates(void) const;

	void allocateSwapMatrices(void);

	// Add one replica
	void addReplica(int index);

	// Add one thermodynamic state
	void addThermodynamicState(int index, SimTK::Real T,

		std::vector<std::string>& rexSamplers,
		std::vector<int>& rexDistortOptions,
		std::vector<std::string>& rexDistortArgs,
		std::vector<int>& rexFlowOptions,
		std::vector<int>& rexWorkOptions,
		std::vector<std::string>& rexIntegrators,

		std::vector<int>& argWorldIndexes,
		std::vector<SimTK::Real>& timestepsInThisReplica,
		std::vector<int>& mdstepsInThisReplica);

	// Prepare Q, U, and tau altering function parameters
	void PrepareNonEquilibriumParams_Q(void);

	// Set thermodynamic states nonequilibrium flag
	void setThermostatesNonequilibrium(void);

	// Set the intial mapping between replicas and thermoStates
	void loadReplica2ThermoIxs(void);

	// Get Fixman potential already calculated from replica
	SimTK::Real getFixman(int replica_i);

	// Calculate Fixman potential of replica J in replica I's back world. Ui(X_j)
	SimTK::Real calcFixman(int replica_i, int replica_j);
	// Calculate Fixman potential of replica J in replica I's back world. Ui(X_j)
	SimTK::Real calcFixman_JinI(int replica_i, int replica_j);
	// Calculate Fixman potential of replica I in replica J's back world. Uj(X_i)
	SimTK::Real calcFixman_IinJ(int replica_i, int replica_j);

	const int& getSwapFixman(void){return swapFixman;}
	void setSwapFixman(const int argSwapFixman){swapFixman = argSwapFixman;}

	SimTK::Real calcReplicaTransferedEnergy(int replicaIx);
	SimTK::Real calcReplicaWork(int replicaIx);


	// SWaps replicas thermodynamic states
	void swapThermodynamicStates(int replica_i, int replica_j);

	// Swap replicas' potential energies
	void swapPotentialEnergies(int replica_i, int replica_j);

	// Exchanges thermodynamic states between replicas
	bool attemptREXSwap(int replica_i, int replica_j);

	bool attemptRENSSwap(int replica_i, int replica_j);

	const int getSwapEvery(void);
	void setSwapEvery(const int& n);

	// We hold this for the moment 
	void mixAllReplicas(int nSwapAttempts);

	// startingFrom argument is for alternating odd and even neighbors
	void mixNeighboringReplicas(unsigned int startingFrom);

	// Mix replicas
	void mixReplicas(void);

	// ========================================================================
	// Configuration manipulation functions between worlds and replicas
	// This can be quite costly since they imply transfer between worlds

	// Load replica's atomLocations into it's front world. Returns world index
	int restoreReplicaCoordinatesToFrontWorld(int whichReplica);

	// Load replica's atomLocations into it's back world
    void restoreReplicaCoordinatesToBackWorld(int whichReplica);

	// Stores replica's front world's coordinates into it's atomsLocations

	// This should always be a fully flexible world
	void storeReplicaCoordinatesFromFrontWorld(int whichReplica);

	// Store work world coordinates into the replica
	void store_WORK_CoordinatesFromFrontWorld(int replicaIx);

	// Store work world energy into the replica 
	void store_WORK_ReplicaEnergyFromFrontWorldFull(int replicaIx);

	// ========================================================================
	// Energy manipulation functions between worlds and replicas
	// This can be quite costly - energy calculation (O^2)

	// Get ennergy of the back world and store it in replica thisReplica
	void storeReplicaEnergyFromBackWorld(int thisReplica);

    	// Get ennergy of the front world and store it in replica thisReplica
	void storeReplicaEnergyFromFrontWorldFull(int thisReplica);

	// Store any WORK Jacobians contribution from back world
	void store_WORK_JacobianFromBackWorld(int replicaIx);

	// Get Fixman of the back world and store it in replica thisReplica
    void storeReplicaFixmanFromBackWorld(int replicaIx);

	// Update replicas coordinates from work generated coordinates
	void set_WORK_CoordinatesAsFinal(int replicaIx);

	// Update replica's energy from work last potential energy
	void set_WORK_PotentialAsFinal(int replicaIx);

	// ------------------------------------------------------------------------

	void initializeReplica(int whichReplica);

	// Reset worlds parameters according to thermodynamic state
	void setReplicasWorldsParameters(int thisReplica);

	// Given a scaling factor
	SimTK::Real distributeScalingFactor(
		std::vector<std::string> how, SimTK::Real sf,
		bool randSignOpt = false);

	// Set world distort parameters
	void setWorldDistortParameters(int whichWorld, SimTK::Real scaleFactor);

	// Set nonequilibrium parameters for one replica
	void updWorldsDistortOptions(int thisReplica);
	void updQScaleFactors(int mixi);

	int getRunType(void){return runType;}
	void setRunType(int runTypeArg){this->runType = runTypeArg;}

	// Run a particular world
	int RunWorld(int whichWorld);

	// Rewind back world
	void RewindBackWorld(int thisReplica);

	// Run front world, rotate and transfer. Return worldIxs.front
	int RunFrontWorldAndRotate(std::vector<int> & worldIxs);

	// Go through all of this replica's worlds and generate samples
	int RunReplicaAllWorlds(int whichReplica, int howManyRounds);

	// Print to log and write pdbs
	void RunLog(int roundi);
	void REXLog(int mixi, int replicaIx);

	void RunREX();

	// Helper Functions for RENS

	int RunReplicaEquilibriumWorlds(int replicaIx, int swapEvery);
	int RunReplicaNonequilibriumWorlds(int replicaIx, int swapEvery);

	// RENS
	void RunRENS(void);

	void PrintReplicas(void);
	void PrintReplicaMaps(void);
	void PrintNofAttemptedSwapsMatrix(void);
	void PrintNofAcceptedSwapsMatrix(void);

	SimTK::Real getPotentialEnergy(std::size_t world, std::size_t sampler) const;

private:
	bool CreateOutputDirectory(const std::string& outDir);
	std::string CreateLogfilename(const std::string& outDir, long long int seed) const;
	std::string GetMoleculeDirectoryShort(const std::string& path) const;
	bool CheckInputParameters(const SetupReader& setupReader);
	void reserveWorldsAndTopologies(int inpNofWorlds, int inpNofMols,
		int inpNofEmbeddedTopologies);

	std::vector<int> TopologyIXs;
	std::vector<std::vector<int>> AmberAtomIXs;
	std::vector<World> worlds;
	std::vector<readAmberInput> amberReader;

	std::vector<int> worldIndexes;
	// Molecules files
	std::vector<std::string> topFNs;
	std::vector<std::string> crdFNs;
	std::vector<std::vector<std::string>> rbSpecsFNs;
	std::vector<std::vector<std::string>> flexSpecsFNs;
	std::vector<std::vector<std::string>> regimens;
	std::vector<std::string> rootMobilities; // WORLD CONFLICT

	// Nof molecules
	int moleculeCount = -1;

	// Molecules (topologies<-Compounds) objects
	//std::vector<bMoleculeReader *> moleculeReaders;
	std::vector<Topology> topologies;
	std::vector<std::string> roots;
	//std::vector<std::string> rootMobilities;

	/** Vectors of Cartesian coordinates **/
	std::vector<SimTK::Real> Xs;
	std::vector<SimTK::Real> Ys;
	std::vector<SimTK::Real> Zs;

	// WORLD END

	// Simulation parameters
	int requiredNofRounds = -1;
	int nofRounds = -1;
	//int total_mcsteps;

	std::size_t nofWorlds = 0;
	bool isWorldsOrderRandom = false;
	std::vector<SimTK::Real> nofSamplesPerRound;
	std::vector<int> nofMDStepsPerSample;
	std::vector<SimTK::Real> timesteps;

	std::vector<int> nofBoostStairs;

	std::size_t nofMols = 0;
	std::size_t nofEmbeddedTopologies = 0; // nofWorlds x nofMols

	int pdbRestartFreq = false;
	int printFreq = -1;

	std::string outputDir;
	std::string pdbPrefix;

	// Geometric features analysis
	// First two integers specifiy the world and the Compound. The rest
	// specifies atom indeces
	std::vector< std::vector<std::size_t> > distanceIxs;
	std::vector< std::vector<std::size_t> > angleIxs;
	std::vector< std::vector<std::size_t> > dihedralIxs;

	SimTK::Real geom1[PRINT_BUFFER_SIZE];
	SimTK::Real geom2[PRINT_BUFFER_SIZE];
	SimTK::Real geom3[PRINT_BUFFER_SIZE];

	// Output
	std::ofstream logFile;

	// Adaptive Gibbs blocking variables
	int roundsTillReblock;
	std::vector<std::vector<std::vector<SimTK::Real>>>
		QsCache; // 1D nofWorlds; 2D roundsTillReblock; 3D nofQs

	// Normal mode analysis
	std::vector<int> NDistortOpt;

	////////////////////////
	//// REPLICA EXCHANGE //
	////////////////////////
	int runType;

	std::vector<ThermodynamicState> thermodynamicStates;
	std::vector<Replica> replicas;

	// Mapping between replicas and thermodynamic states indexes
	// KEYWORD = replica, VALUE = thermoState
	std::map<int, int> replica2ThermoIxs;

	// Mapping between replicas and thermodynamic states indexes
	// KEYWORD = thermoState, VALUE = replica
	std::map<int, int> thermo2ReplicaIxs;

	// Counter matrix of accepted swaps
	std::vector<std::vector<int>> nofAttemptedSwapsMatrix;
	std::vector<std::vector<int>> nofAcceptedSwapsMatrix;

	size_t nofReplicas = 0;
	size_t nofThermodynamicStates = 0;
	ReplicaMixingScheme replicaMixingScheme = ReplicaMixingScheme::neighboring;

	int swapFixman = 1;

	std::random_device rd;
	std::mt19937 randomEngine;
	std::uniform_real_distribution<SimTK::Real> uniformRealDistribution =
		    std::uniform_real_distribution<SimTK::Real>(SimTK::Zero, SimTK::One);

	int swapEvery = 1;

	// Non-equilibrium parameters
	std::vector<SimTK::Real> qScaleFactorsEven;
	std::vector<SimTK::Real> qScaleFactorsOdd;
	std::vector<SimTK::Real> qScaleFactorsMiu;
	std::vector<SimTK::Real> qScaleFactorsStd;

	std::vector<SimTK::Real> qScaleFactors;

	std::string cerr_prefix = "[ERROR] ";

	RunType run_type = RunType::Default;
	SimTK::Real tempIni = 0,
		tempFin = 0;

	SetupReader setupReader;
};

