#pragma once

#include "Robo.hpp"
#include "Sampler.hpp"
#include "SetupReader.hpp"
#include "World.hpp"
#include "ThermodynamicState.hpp"
#include "Replica.hpp"
#include "SetupReader.hpp"
#include "TrajectoryObject.hpp"

#include <fstream>

class Sampler;
class World;

enum class ReplicaMixingScheme : int {
	all = 0,
	neighboring = 1
};

enum class RUN_TYPE : int {
	DEFAULT = 0,
	REMC,
	RENEMC,
	RENE
};

class Context {
public:
	int natoms = std::numeric_limits<int>::min();
	std::vector<bSpecificAtom> atoms;
	std::string baseName;

	/**
	 * @brief Initialize simulation variables.
	 * @param Ti Initial temperature. Cannot be 0.
	 * @param Tf Final temperature. Cannot be 0. Must be greater than Ti.
	 * @param seed Seed to use for random number generation. If 0, a random seed is used.
	 * @param threads Number of threads to use.
	 * @param nofRoundsTillReblock Number of rounds until reblocking.
	 * @param runType Type of simulation to run.
	*/
	Context(std::string baseName, uint32_t seed, uint32_t threads, uint32_t nofRoundsTillReblock, RUN_TYPE runType, uint32_t swapFreq, uint32_t swapFixmanFreq);

	/**
	 * @brief Print atom's Compound and DuMM indexes.
	*/
	void PrintAtomsDebugInfo(void);

 	std::vector<BOND_FLEXIBILITY>& readFlexibility(
		std::string flexFileFN,
		std::vector<BOND_FLEXIBILITY>& flexibilities);

	/**	
	* @brief Read all parameters from an input file
	* @param filename input file name
	* @return succes of the function
	*/	
	bool initializeFromFile(const std::string& filename);

	void setNonbonded(int method, SimTK::Real cutoff);
	void setGBSA(SimTK::Real globalScaleFactor);
	void setForceFieldScaleFactors(SimTK::Real globalScaleFactor);

	bool setOutput(const std::string& outDir);

	void loadAmberSystem(const std::string& prmtop, const std::string& inpcrd);
	void initialize();

	int BAT2Amber(int batIndex) {
		return internCoords.BAT2amber(batIndex);
	}

	void appendLog(const std::string& filename);
	void appendDCDReporter(const std::string& filename);

	// Input functions
	bool loadCoordinatesFile(std::string coordinatesFilename);
	void PrintCoordinates(const std::vector<std::vector
        <std::pair <bSpecificAtom *,
		SimTK::Vec3>>>& atomsLocations);

	bool loadRigidBodiesSpecs(std::size_t whichWorld, int whichMolecule, std::string RBSpecsFN);
	bool loadFlexibleBondsSpecs(std::size_t whichWorld, std::string FlexSpecsFN);
	void setRegimen (std::size_t whichWorld, int whichMolecule, std::string regimen);

	// ============================================================================
	// ============================================================================
	// ==========================   SINGLE PRMTOP    ==============================
	// ============================================================================
	// ============================================================================

	/**  */
	void setRootAtom(Topology& topology, int molIx);

	/**  */
	void load_BONDS_to_bonds(const std::vector<std::vector<BOND>>& BATbonds);

	/** If bonds are resorted */
	void reset_BONDS_to_bonds(const std::vector<std::vector<BOND>>& BATbonds);

	/**  */
	void buildAcyclicGraph(
		Topology& topology,
		int rootAmberIx,
		int molIx);

	void closeARingWithThisBond(Topology& topology, bBond& bond, int molIx);

	/**  */
	void addRingClosingBonds(
		Topology& topology,
		int rootAmberIx,
		int molIx
	);

	void addRingClosingBonds_All(void);

	/**  */
	void generateSubAtomLists(void);

	/**  */
	void generateSubBondLists(void);

	/** Pass Context topologies to all the worlds */
	void passTopologiesToWorlds(void);

	/** Read atoms and bonds from all the molecules */
	void readMolecules(void);
	
	/**  */
	void constructTopologies();

	/**  */
	void generateTopologiesSubarrays(void);

	/** Assign Compound coordinates by matching bAtomList coordinates */
	void matchDefaultConfiguration(Topology& topology, int molIx);

	/** Match Compounds configurations to atoms Cartesian coords */
	void matchDefaultConfigurations(void);

	// ------------- PARAMETERS -------------

	/** Long print of all atoms properties */
	void PrintAtoms(void);

	/** It calls DuMMs defineAtomClass. These Molmodel functions contain
	information regarding the force field parameters. **/
	void updDummAtomClasses(
		std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
		, int worldIx		
	);

	// Add Dumm params for single prmtop
	void generateDummAtomClasses(readAmberInput& amberReader);
	void bAddDummBondParams(readAmberInput& amberReader);
	void bAddDummAngleParams(readAmberInput& amberReader);
	bool checkBond(int a1, int a2);
	void bAddDummTorsionParams(readAmberInput& amberReader);

	void addDummParams(readAmberInput& amberReader);

	// -------------------------

	// ---------
	/** Set all flexibilities for all the worlds to Rigid. */
	void initializeFlexibility(void);

	/** Set flexibilities. */
	void setFlexibility(
		std::string argRegimen,
		std::string flexFN,
		int whichWorld);

	void modelOneEmbeddedTopology(int whichTopology,
		int whichWorld
		//,std::string rootMobilizer
		);
	// ---------

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

	// --- Simulation parameters ---

	//------------

	// --- Mixing parameters ---
	// Another way to do it is setting the number of rounds
	int getRequiredNofRounds();
	void setRequiredNofRounds(int argNofRounds);

	std::size_t getWorldIndex(std::size_t which) const;

	// Adaptive Gibbs blocking: TODO: consider moving in World
	void allocateReblockQsCache(void);
	void allocateReblockQsCacheQVectors(void);

	// --- Arrange different mixing parameters ---
	void initializeMixingParamters();
	//------------

	void addWorld(
		bool fixmanTorque,
		int samplesPerRound,
		ROOT_MOBILITY rootMobility,
		const std::vector<BOND_FLEXIBILITY>& flexibilities,
		bool useOpenMM = true,
		bool visual = false,
		SimTK::Real visualizerFrequency = 0);

	std::size_t getNofWorlds() const;
	World& getWorld(std::size_t which);
	const World& getWorld(std::size_t which) const;

	std::vector<World>& getWorlds();
	const std::vector<World>& getWorlds() const;

	// Writeble reference to a samplers advanced state
	SimTK::State& updAdvancedState(std::size_t whichWorld, std::size_t whichSampler);

	void RotateWorlds();
	//------------

	// --- Main ---
	void randomizeWorldIndexes(void);
	void transferCoordinates(int src, int dest);
	SimTK::Real checkTransferCoordinates_Cart(int srcWIx, int destWIx);
	SimTK::Real checkTransferCoordinates_BAT(int srcWIx, int destWIx, bool wantJacobian = false);
	
	// Relationship BAT - mobod transforms
	void PrintZMatrixMobods(int wIx, SimTK::State& someState);

	// TRANSFER ============================================================

	SimTK::State&
	setAtoms_CompoundsAndDuMM(
		int destWIx,
		SimTK::State& someState,
		const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>> &
			otherWorldsAtomsLocations);

	SimTK::State&
	setAtoms_XPF_XBM(
		int wIx
	);

	SimTK::State&
	setAtoms_MassProperties(
		int wIx
	);

	SimTK::Transform
	calc_XFM(
		int wIx,
		Topology& topology,	
		SimTK::Compound::AtomIndex& childAIx,
		SimTK::Compound::AtomIndex& parentAIx,
		SimTK::BondMobility::Mobility mobility,
		const SimTK::State& someState) const;

	SimTK::State&
	setAtoms_XFM(
		int wIx,
		SimTK::State& someState);

	std::vector<SimTK::Transform>
	calc_XPF_XBM(
		int wIx,
		Topology& topology,
		SimTK::Compound::AtomIndex& childNo,
		SimTK::Compound::AtomIndex& parentNo,
		SimTK::BondMobility::Mobility mobility,
		const SimTK::State& someState
	);

	SimTK::State&
	setAtoms_SP_NEW(
		int destWIx,
		SimTK::State& someState,
		const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>> &
			otherWorldsAtomsLocations);

	SimTK::Compound::AtomIndex
	getChemicalParent_IfIAmRoot(
		int wIx,
		int atomNo,
		SimTK::DuMMForceFieldSubsystem &dumm
	);

	// X axis to Z axis switch
	const SimTK::Transform X_to_Z 
		=  SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis);
	const SimTK::Transform Z_to_X = ~X_to_Z;

	// Y axis to Z axis switch
	const SimTK::Transform Y_to_Z =
		SimTK::Transform(SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::XAxis));
	const SimTK::Transform Z_to_Y = ~Y_to_Z;

	// X axis to X axis switch
	const SimTK::Transform Y_to_X =
		SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::ZAxis);
	const SimTK::Transform X_to_Y = ~Y_to_X;

	// TRANSFER ------------------------------------------------------------

	// Drilling drl
	void passThroughBonds_template(int whichWorld);

	// Go through all the worlds and generate samples
	void RunOneRound(void);
	void Run(int howManyRounds, SimTK::Real Ti, SimTK::Real Tf);

	SimTK::Real Pearson(std::vector<std::vector<SimTK::Real>> someVector,
		int QIx1, int QIx2); // 2D roundsTillReblock; 3D nofQs

	/** Print the number of threads each World got **/
	void PrintNumThreads();

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
	void addThermodynamicState(int index,
		SimTK::Real T,
		SimTK::Real boostT,
		std::vector<AcceptRejectMode>& rexSamplers,
		std::vector<int>& rexDistortOptions,
		std::vector<std::string>& rexDistortArgs,
		std::vector<int>& rexFlowOptions,
		std::vector<int>& rexWorkOptions,
		std::vector<int>& argWorldIndexes,
		std::vector<SimTK::Real>& timestepsInThisReplica,
		std::vector<int>& mdSteps,
		std::vector<int>& boostMDSteps);

	/**
	* @brief zmatrixbat_
	* @param
	*/	
	void setReplicaExchangePairs(unsigned int startingFrom);

	/**
	* @brief zmatrixbat_
	* @param
	*/
	const int getThermoPair(int replicaIx);


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

	SimTK::Real calcReplicaTransferedEnergy(int replicaIx);
	SimTK::Real calcReplicaWork(int replicaIx);


	// SWaps replicas thermodynamic states
	void swapThermodynamicStates(int replica_i, int replica_j);

	// Swap replicas' potential energies
	void swapPotentialEnergies(int replica_i, int replica_j);

	// Exchanges thermodynamic states between replicas
	void getMsg_RexDetHeader(std::stringstream& rexDetHeader);
	bool attemptREXSwap(int replica_i, int replica_j);

	// startingFrom argument is for alternating odd and even neighbors
	void mixNeighboringReplicas(unsigned int startingFrom);

	// We hold this for the moment 
	void mixAllReplicas(int nSwapAttempts);

	// Mix replicas
	void mixReplicas(int mixi);

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
	void setReplicasWorldsParameters(int thisReplica, bool alwaysAccept, bool adaptTimestep, int mixi);

	// Given a scaling factor
	SimTK::Real distributeScalingFactor(
		std::vector<std::string> how, SimTK::Real sf,
		bool randSignOpt = false);

	// Set world distort parameters
	void setWorldDistortParameters(int whichWorld, SimTK::Real scaleFactor);

	// Set nonequilibrium parameters for one replica
	void updWorldsDistortOptions(int thisReplica);
	void updQScaleFactors(int mixi);

	RUN_TYPE getRunType(void) const;
	void setRunType(RUN_TYPE runTypeArg);

	// Run a particular world
	bool RunWorld(int whichWorld);

	// Rewind back world
	void RewindBackWorld(int thisReplica);

	// Run front world, rotate and transfer. Return worldIxs.front
	int RunFrontWorldAndRotate(std::vector<int> & worldIxs);
	void RunWorlds(std::vector<int>& specificWIxs, int replicaIx);

	// Print to log and write pdbs
	void RunLog(int roundi);
	void REXLog(int mixi, int replicaIx);
	void writeLog(int mixi, int replicaIx);

	void incrementNofSamples(void);


	void RunREX();

	void RunREXNew(int equilRounds, int prodRounds);

	// Helper Functions for REX

	int RunReplicaEquilibriumWorlds(int replicaIx, int swapEvery);

	void setSubZmatrixBATStatsToSamplers(int thermoIx, int worldCnt);
	int RunReplicaNonequilibriumWorlds(int replicaIx, int swapEvery);

	// Go through all of this replica's worlds and generate samples
	int RunReplicaAllWorlds(int mixi, int replicaIx, int swapEvery);

	// Transfer Q statistics
	void transferQStatistics(int thermoIx, int srcStatsWIx, int destStatsWIx);
	void RunReplicaRefactor(int mixi, int replicaIx);

	void PrintReplicas(void);
	void PrintReplicaMaps(void);
	void PrintNofAttemptedSwapsMatrix(void);
	void PrintNofAcceptedSwapsMatrix(void);

	SimTK::Real getPotentialEnergy(std::size_t world, std::size_t sampler) const;


	//////////////////////////////////
	/////     TEST FUNCTIONS     /////
	//////////////////////////////////
	void areAllDuMMsTheSame(void);

	void PrintBond(bBond& bond);
	void PrintBonds(void);
	int checkBonds(void);

	// Transformers
	void Print_TRANSFORMERS_Work(void);

	// Function to find and return the value for a given AtomIndex
	SimTK::Vec3
	findAtomTarget(
		const std::map<SimTK::Compound::AtomIndex, SimTK::Vec3>& atomTargets,
		SimTK::Compound::AtomIndex searchIndex)
	{
		auto it = atomTargets.find(searchIndex);

		if (it != atomTargets.end()) {
			return it->second;
		} else {
			return SimTK::Vec3(SimTK::NaN);
		}
	}

private:
	std::string GetMoleculeDirectoryShort(const std::string& path) const;
	bool CheckInputParameters(const SetupReader& setupReader);

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
	std::vector<std::vector<std::string>> rootMobilities;

	// Nof molecules
	int moleculeCount = -1;

	// Molecules (topologies<-Compounds) objects
	//std::vector<bMoleculeReader *> moleculeReaders;
	std::vector<Topology> topologies;
	std::vector<int> roots;
	//std::vector<std::string> rootMobilities;
	InternalCoordinates internCoords;

	// WORLD END

	// Simulation parameters
	int requiredNofRounds = -1;
	int nofRounds = -1;
	//int total_mcsteps;

	std::size_t nofWorlds = 0;
	bool isWorldsOrderRandom = false;

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

	////////////////////////
	//// REPLICA EXCHANGE //
	////////////////////////

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
	std::string cwar_prefix = "[WARNING] ";
	std::string cinf_prefix = "[INFO] ";

	RUN_TYPE runType = RUN_TYPE::DEFAULT;
	SimTK::Real tempIni = 0,
		tempFin = 0;

	SetupReader setupReader;

	// Every molecule has an array_view for atoms and bonds
	std::vector<array_view<std::vector<bSpecificAtom>::iterator>>
		subAtomLists;
	std::vector<array_view<std::vector<bBond>::iterator>>
		subBondLists;

	
	int nbonds = std::numeric_limits<int>::min();
	std::vector<bBond> bonds;
	std::vector<std::vector<int>> BONDS_to_bonds; // correspondence
	std::vector<std::pair<int, int>> bonds_to_BONDS;

	std::vector<DUMM_ANGLE> dummAngles;
	std::vector<DUMM_TORSION> dummTorsions;
	
	ELEMENT_CACHE elementCache;

	void loadAtoms(const readAmberInput& reader);

	void loadAtomsCoordinates(const std::string& prmtop, const std::string& inpcrdFN);
	
	void loadBonds(const readAmberInput& reader);
	void loadAngles(const readAmberInput& reader);
	void loadTorsions(const readAmberInput& reader);
	void setAtomCompoundTypes();
	void addBiotypes();
	std::vector<bSpecificAtom>& getAtoms() {
        return atoms;
    }

	uint32_t seed = 0;
	int numThreads = 0;
	int nonbondedMethod = 0; // 0 = NoCutoff, 1 = CutoffNonPeriodic
	SimTK::Real nonbondedCutoff = 1.2; // 1.2 nm, not used by default (no cutoff)
	SimTK::Real gbsaGlobalScaleFactor = 0.0; // Default is 0 (vacuum). Use 1 for water

	bool useAmberForceFieldScaleFactors = true;
	SimTK::Real globalForceFieldScaleFactor = 1.0; // Used in place of Amber scaling (not default)

	// Random number generator
	Random32 randomEngine;

public:
	/** Implicit membrane mimicked by half-space contacts */
	void addContactImplicitMembrane(const float memZWidth, const SetupReader& setupReader);

	std::map<std::string, BondMobility::Mobility> mobilityMap = {
		{ "Free", BondMobility::Free },
		{ "Pin", BondMobility::Torsion },
		{ "Cartesian", BondMobility::Translation },
		{ "Rigid", BondMobility::Rigid },
		{ "Weld", BondMobility::Rigid },
		{ "Slider", BondMobility::Slider },
		{ "AnglePin", BondMobility::AnglePin },
		{ "BendStretch", BondMobility::BendStretch },
		{ "Spherical", BondMobility::Spherical }
	};

    std::vector<std::string> MobilityStr {
		"Zero",
        "Free",
        "Torsion",
        "Rigid", 
        "BallF", 
        "BallM", 
        "Cylinder", 
        "Translation", 
        "FreeLine", 
        "LineOrientationF", 
        "LineOrientationM", 
        "UniversalM", 
        "Spherical", 
        "AnglePin",
        "BendStretch",
        "Slider"   
    };

	BondMobility::Mobility getMobility(const std::string& mobilityStr) {
		// Assume MobilityStr is a vector defined elsewhere in your code
		// std::vector<std::string> MobilityStr = { ... };

		auto it = std::find(MobilityStr.begin(), MobilityStr.end(), mobilityStr);
		
		if (it != MobilityStr.end()) {
			// If the string is found, return the corresponding enum value
			return static_cast<BondMobility::Mobility>(std::distance(MobilityStr.begin(), it) + 1);
		} else {
			// If the string is not found, return the default value
			return BondMobility::Default;
		}
	}	

private:
	std::vector<SimTK::Real> Xs, Ys, Zs;

    /** @name Z Matrix and BAT functions
	*/

    /**@{**/
	
	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

	std::vector<std::vector<int>> zMatrixTable;
	std::vector<std::vector<SimTK::Real>> zMatrixBAT;

	/**	
	* @brief
	* @param
	* @return
	*/
    // zmatrixbat_ Function to add a new row to the zMatrixTable
    void addZMatrixTableRow(const std::vector<int>& newRow) ;

    // zmatrixbat_ Getter for a specific entry
    int getZMatrixTableEntry(int rowIndex, int colIndex) const ;

    // zmatrixbat_ Setter for a specific entry
    void setZMatrixTableEntry(int rowIndex, int colIndex, int value) ;

    // zmatrixbat_ Print function for the zMatrixTable
    void PrintZMatrixTable() const ;

    // zmatrixbat_ Setter for a specific entry
    void setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value) ;

    // zmatrixbat_ Function to get a given row
    const std::vector<SimTK::Real>& getZMatrixBATRow(size_t rowIndex) const;

    // zmatrixbat_ Function to get a given row
    std::vector<SimTK::Real>& updZMatrixBATRow(size_t rowIndex) ;


	// Get Z-matrix indexes table 
	void
	calcZMatrixTable(void);

	// Allocate Z Matrix BAT
	void reallocZMatrixBAT(void);

	//
	void
	calcZMatrixBAT(	int wIx,
	const std::vector< std::vector<
	std::pair <bSpecificAtom *, SimTK::Vec3 > > >&
		otherWorldsAtomsLocations);

    // zmatrixbat_ Function to get the value for a given row and column in zMatrixBAT
    SimTK::Real getZMatrixBATValue(size_t rowIndex, size_t colIndex) const ;

    // zmatrixbat_ Function to print the zMatrixBAT
    void PrintZMatrixBAT() const ;

	// zmatrixbat_ 
	void PrintZMatrixTableAndBAT() const;

    // zmatrixbat_  Function to add a new row to the zMatrixBAT
    void addZMatrixBATRow(const std::vector<SimTK::Real>& newRow);

	// WORK Q PERTURB BEND STRETCH ============================================

	/**
	* @brief zmatrixbat_ Get log of the Cartesian->BAT Jacobian
	* @param
	*/
	SimTK::Real
	calcInternalBATJacobianLog(void);

	/**
	* @brief zmatrixbat_ Add BAT coordinates 
	* @param
	*/
	void
	addSubZMatrixBATsToWorld(
		int wIx, int replicaIx);

	/**
	* @brief zmatrixbat_ Get BAT coordinates modifyable by a selected world
	* @param
	*/
	void
	updSubZMatrixBATsToWorld(
		int wIx, int replicaIx);

	/**
	* @brief zmatrixbat_ Set BAT coordinates modifyable to all worlds of a replica
	* @param
	*/
	void updSubZMatrixBATsToAllWorlds(int replicaIx);

	/**
	* @brief zmatrixbat_ Print BAT coordinates
	* @param
	*/
	void
	PrintWorldSubZMatrixBATs(
		int wIx);

	// WORK Q PERTURB BEND STRETCH --------------------------------------------

	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

	/**@}**/

	// Molmodel to Gmolmodel (and inverse) bond mappings
	// bondMapping[gmolmodelBondIndex] = molmodelBondIndex (inverse mapping)
	// bondMapping[molmodelBondIndex] = gmolmodelBondIndex (normal mapping)
	std::unordered_map<int, SimTK::Compound::BondIndex> bondMapping;

	// Pairs of replicas to be exchanged
	std::vector<int> exchangePairs;

	std::map<std::string, AcceptRejectMode> sampleGenerator = {
		{ "EMPTY", AcceptRejectMode::AlwaysAccept },
		{ "MC", AcceptRejectMode::MetropolisHastings },
	};

	std::map<std::string, IntegratorName> integratorName = {
		{ "EMPTY", IntegratorName::None },
		{ "VV", IntegratorName::Verlet },
		{ "VERLET", IntegratorName::Verlet },
		{ "EULER", IntegratorName::Euler },
		{ "EULER2", IntegratorName::Euler2 },
		{ "CPODES", IntegratorName::Cpodes },
		{ "RUNGEKUTTA", IntegratorName::RungeKutta },
		{ "RUNGEKUTTA2", IntegratorName::RungeKutta2 },
		{ "RUNGEKUTTA3", IntegratorName::RungeKutta3 },
		{ "RUNGEKUTTAFELDBERG", IntegratorName::RungeKuttaFeldberg },
		{ "BENDSTRETCH", IntegratorName::BendStretch },
		{ "OMMVV", IntegratorName::OMMVV },
		{ "BOUND_WALK", IntegratorName::BoundWalk },
		{ "BOUND_HMC", IntegratorName::BoundHMC },
		{ "STATIONS_TASK", IntegratorName::StationsTask },
		{ "NOF_INTEGRATORS", IntegratorName::NofIntegrators },
	};

	std::map<std::string, ThermostatName> thermostateName = {
		{ "NONE", ThermostatName::NONE },
		{ "ANDERSEN", ThermostatName::ANDERSEN },
		{ "BERENDSEN", ThermostatName::BERENDSEN },
		{ "LANGEVIN", ThermostatName::LANGEVIN },
		{ "NOSE_HOOVER", ThermostatName::NOSE_HOOVER },
	};


	//////////////////////////////////
	//---         Q Stats        -----
	//////////////////////////////////
	void reserveThermostatsQs(void);
	
	void setThermostatesQs(void);

	void calcQStats(int thIx);

	void printQStats(int thIx);

	//////////////////////////////////
	//---         Q Stats        -----
	//////////////////////////////////

};
