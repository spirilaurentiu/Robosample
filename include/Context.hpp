#pragma once

#include "Robo.hpp"
#include "Sampler.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"
#include "ThermodynamicState.hpp"
#include "Replica.hpp"
#include "TrajectoryObject.hpp"
#include "OpenMM.hpp"

class Sampler;
class World;

enum class ReplicaMixingScheme : int {
	all = 0,
	neighboring = 1
};

enum class RUN_TYPE : int {
	DEFAULT = 0,
	REMC, // Replica Exchange Monte Carlo
	RENEMC, // Replica Exchange Non-Equilibrium Monte Carlo
	RENE, // Replica Exchange Non-Equilibrium
	REBASONTOP // Replica
};

enum class TopologyRangeType : int {
	ATOM = 0,
	BOND,
	ANGLE,
	PERIODIC_TORSION,
	IMPROPER_HARMONIC_TORSION,
	COUNT
};

class TopologyRange {
	// Stores [begin, end) pairs for each type
    std::array<std::pair<int, int>, (int)TopologyRangeType::COUNT> ranges;

public:
	TopologyRange(std::vector<int> startCounts) {
		ranges[(int)TopologyRangeType::ATOM] = { startCounts[0], startCounts[0] };
		ranges[(int)TopologyRangeType::BOND] = { startCounts[1], startCounts[1] };
		ranges[(int)TopologyRangeType::ANGLE] = { startCounts[2], startCounts[2] };
		ranges[(int)TopologyRangeType::PERIODIC_TORSION] = { startCounts[3], startCounts[3] };
		ranges[(int)TopologyRangeType::IMPROPER_HARMONIC_TORSION] = { startCounts[4], startCounts[4] };
	}

	void close(std::vector<int> endCounts) {
		ranges[(int)TopologyRangeType::ATOM].second = endCounts[0];
		ranges[(int)TopologyRangeType::BOND].second = endCounts[1];
		ranges[(int)TopologyRangeType::ANGLE].second = endCounts[2];
		ranges[(int)TopologyRangeType::PERIODIC_TORSION].second = endCounts[3];
		ranges[(int)TopologyRangeType::IMPROPER_HARMONIC_TORSION].second = endCounts[4];
	}

	const std::pair<int, int>& getRange(TopologyRangeType type) const {
		return ranges[(int)type];
	}
};

class Context{

    /** @name Constructor **/
    /**@{**/
	/**@}**/

	std::string baseName;
	bool verbose = false;

public:

// vector<vector<ATOM>> for each molecule
	// void createSystem(std::vector<ATOM> atoms, std::vector<BOND> bonds);

	/**
	 * @brief Initialize simulation variables.
	 * @param baseName Base name for output files.
	 * @param Ti Initial temperature. Cannot be 0.
	 * @param Tf Final temperature. Cannot be 0. Must be greater than Ti.
	 * @param seed Seed to use for random number generation. If 0, a random seed is used.
	 * @param threads Number of threads to use.
	 * @param nofRoundsTillReblock Number of rounds until reblocking.
	 * @param runType Type of simulation to run.
	*/
	Context(const std::string& baseName, uint32_t seed, uint32_t threads, uint32_t nofRoundsTillReblock, RUN_TYPE runType, uint32_t swapFreq, uint32_t swapFixmanFreq, bool testing);

	void setVerbose(bool verbose);
	void setNumThreads(int threads);
	void setGBSAOptions(bool useGBSAOBC2, SimTK::Real solventDielectric, SimTK::Real soluteDielectric);
	void setForceFieldScaleFactors(SimTK::Real globalScaleFactor);
	bool setOutput(const std::string& outDir);
	
	void setNofRoundsTillReblock(int nofRoundsTillReblock);
	void setRequiredNofRounds(int argNofRounds);

	void setNonbonded(NonbondedMethod method, SimTK::Real cutoffInNm);

	void loadAmberSystem(
		const std::vector<int>& roots_,
		const std::vector<RoboAtom>& atoms_,
		const std::vector<RoboBond>& bonds_,
		const std::vector<RoboAngle>& angles_,
		const std::vector<RoboPeriodicTorsion>& properPeriodicTorsions_,
		const std::vector<RoboHarmonicImproperTorsion>& harmonicImproperTorsions_,
		const std::vector<TopologyRange>& topologyRanges,
		const ZMatrix& _zMatrix
	);

	void initializeOpenMM(
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
        const std::vector<Scaling14>& scaling14s
	);
	SimTK::Real calculatePotentialEnergy(int worldIndex);

	const std::string& getAtomNameByPrmtopIndex(int prmtopIndex) const {
		for (const auto& atom : atoms) {
			if (atom.identity.prmtopIndex == prmtopIndex) {
				return atom.identity.uniqueAtomName;
			}
		}
		throw std::runtime_error("Atom with specified prmtop index not found.");
	}

	void Initialize();

	void addWorld(
		bool fixmanTorque,
		int samplesPerRound,
		ROOT_MOBILITY rootMobility,
		const std::vector<std::vector<BondFlexibility>>& flexibilities,
		bool useOpenMM = true,
		bool visual = false,
		SimTK::Real visualizerFrequency = 0
	);

	// Add task spaces
	void addTaskSpacesLS(void);

	/** Add constraints */
	void addConstraints(void);

	void passTopologiesToNewWorld(int newWorldIx);

	// --- Simulation parameters ---

	//------------

	// --- Mixing parameters ---
	// Another way to do it is setting the number of rounds
	int getRequiredNofRounds();

	int getNofRoundsTillReblock();

	void updNofRoundsTillReblock(int nofRoundsTillReblock);

	std::size_t getWorldIndex(std::size_t which) const;

	// Adaptive Gibbs blocking: TODO: consider moving in World
	void allocateReblockQsCache(void);
	void allocateReblockQsCacheQVectors(void);

	// --- Arrange different mixing parameters ---
	void initializeMixingParamters();
	//------------

	std::size_t getNofWorlds() const {
		return worlds.size();
	}

	World& getWorld(std::size_t which) {
		return worlds[which];
	}

	const World& getWorld(std::size_t which) const {
		return worlds[which];
	}

	std::vector<World>& getWorlds() {
		return worlds;
	}

	const std::vector<World>& getWorlds() const {
		return worlds;
	}

	void RotateWorlds();
	//------------

	// --- Main ---
	void randomizeWorldIndexes(void);

	// Relationship BAT - mobod transforms
	void PrintZMatrixMobods(int wIx, SimTK::State& someState);

	// Drilling drl
	void passThroughBonds_template(int whichWorld);

	SimTK::Real Pearson(std::vector<std::vector<SimTK::Real>> someVector, int QIx1, int QIx2); // 2D roundsTillReblock; 3D nofQs

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

	// Print DuMM atoms stations in mobilized body frame
	void checkAtomStationsThroughDumm(void);

	// Print Simbody related information
	void PrintSimbodyMobods(void);
	void PrintFreeE2EDist(std::size_t whichWorld, int whichCompound);

	void PrintGeometryToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintDistancesToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintAnglesToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintDihedralsToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintDihedralsQsToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintSamplerDataToLog(std::size_t whichWorld, std::size_t whichSampler);
	void PrintToLog(std::size_t whichReplica,
		std::size_t whichWorld, std::size_t whichSampler);

    /** @name Write coordinates **/
    /**@{**/

	/**	
	* @brief Write pdb file.
	*/
	void writeInitialPdb(void);
	void writeFinalPdb(void);
	void writePdbs(int someIndex, int thermodynamicStateIx = 0);

	// Output helpers
	int getPdbRestartFreq();
	void setPdbRestartFreq(int argFreq);

	const std::string& getRestartDir() const;
	void setRestartDir(const std::string& argRestartDir);

	void setPdbPrefix(const std::string& argPdbPrefix);
	std::string getPdbPrefix();

	int getPrintFreq();
	void setPrintFreq(int argFreq);

	std::string getOutputDir();
	void setOutputDir(std::string arg);

	/**@}**/

	SimTK::Real Dihedral(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t whichSampler, int a1, int a2, int a3, int a4);
	SimTK::Real Roboangle(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t whichSampler, int a1, int a2, int a3);
	SimTK::Real Distance(std::size_t whichWorld, std::size_t whichCompound,
		std::size_t whichSampler, int a1, int a2);

	//////////////////////////////////
	//// REPLICA EXCHANGE FUNCTIONS //
	//////////////////////////////////

	void setNofReplicas(const size_t& argNofReplicas);
	const size_t& getNofReplicas(void) const;
	//void setNofThermodynamicStates(const size_t& argNofThermodynamicStates);
	const size_t& getNofThermodynamicStates(void) const;

	void allocateSwapMatrices(void);

	// Add one replica
	void addReplica();

	// Add one thermodynamic state
	void addThermodynamicState(
		SimTK::Real T,
		const std::vector<AcceptRejectMode>& acceptRejectModes,
		const std::vector<int>& rexDistortOptions,
		const std::vector<std::string>& rexDistortArgs,
		const std::vector<int>& rexFlowOptions,
		const std::vector<int>& rexWorkOptions,
		const std::vector<IntegratorType>& rexIntegrators,
		const std::vector<int>& argWorldIndexes,
		const std::vector<SimTK::Real>& timestepsInThisReplica,
		const std::vector<int>& mdstepsInThisReplica
	);

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
	void swapReferencePotentialEnergies(int replica_i, int replica_j);

	// Exchanges thermodynamic states between replicas
	void getMsg_RexDetHeader(std::stringstream& rexDetHeader);
	void rewindReplica();
	bool attemptREXSwap(int thermoState_C, int thermoState_H);

	int getSwapEvery() const { return swapEvery; }
	void setSwapEvery(int n) { swapEvery = n; }

	// Mix replicas
	void mixAllReplicas(int nSwapAttempts);
	
	void prepareExchangePairs(int rexRound, int oddity);
	void mixReplicas(int mixi, int oddity = 0);

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
	SimTK::Real perturbScalingFactor(
		std::string how, SimTK::Real sf,
		bool randSignOpt = false);

	// Set world distort parameters
	void setWorldDistortParameters(int whichWorld, SimTK::Real scaleFactor);

	// Set nonequilibrium parameters for one replica
	void updWorldsDistortOptions(int thisReplica);
	void updThermostatesQScaleFactors(int mixi);

	// Rewind back world
	void RewindBackWorld(int thisReplica);

	// Run front world, rotate and transfer. Return worldIxs.front
	int RunFrontWorldAndRotate(std::vector<int> & worldIxs);

	// Print to log and write pdbs
	void RunLog(int roundi);
	void REXLog(int mixi, int replicaIx);

	void writeLog(int mixi, int replicaIx);

	void incrementNofSamples(void);

    /** @name Replica exchange **/
    /**@{**/

	// Run a particular world
	bool RunWorld(int whichWorld, const std::string& header);
	void RunReplicaWorldRange(int replicaIx, int startWorldCnt, int nofWorldsCounted, bool isNonEquilibrium);	
	/**	
	* @brief Main function
	* @param
	* @return
	*/
	void RunREX(int equilRounds, int prodRounds);

	void transferCoordsFromWorldToWorld(int sourceWorldIndex, int destinationWorldIndex);
	void transferCoordsFromWorldToReplica(int sourceWorldIndex, int destinationReplicaIndex, bool intoWORK);
	void transferCoordsFromReplicaToWorld(int sourceReplicaIndex, int destinationWorldIndex);

	void setSubZmatrixBATStatsToSamplers(int thermoIx, int worldCnt);

	// Transfer Q statistics
	void transferQStatistics(int thermoIx, int srcStatsWIx, int destStatsWIx);

	void PrintReplicas(void);
	void PrintReplicaMaps(void);
	void PrintNofAttemptedSwapsMatrix(void);
	void PrintNofAcceptedSwapsMatrix(void);

	// Transformers
	void Print_TRANSFORMERS_Work(void);

	std::vector<SimTK::Real> AtomIndex0, AtomIndex1, UDotCache, UCache;
	bool binaryFileIsInitialized = false;

	std::string foutU, foutUDot, foutTorque;

	void initializeBinaryFile(const std::string &filename, uint32_t num_columns);
	void writeRowToBinaryFile(const std::string &filename, const std::vector<SimTK::Real> &row, bool has_acceptance, bool accepted);

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

	// Run in testing mode
	bool testing = false;

	std::vector<int> TopologyIXs;
	std::vector<std::vector<int>> AmberAtomIXs;
	std::vector<World> worlds;

	std::vector<int> worldIndexes;
	std::vector<std::vector<std::string>> rootMobilitiesStr;

	int moleculeCount = -1;

	// Simulation parameters
	int requiredNofRounds = -1;
	int nofRounds = -1;

	std::size_t nofWorlds = 0;
	bool isWorldsOrderRandom = false;

	int pdbRestartFreq = false;
	int printFreq = -1;

	std::string molDir;
	std::string outputDir;
	std::string restartDir;
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

	std::size_t nofReplicas = 0;
	std::size_t nofThermodynamicStates = 0;
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
	SimTK::Real tempIni = 0, tempFin = 0;

	std::vector<RoboAtom> atoms;
	std::vector<RoboBond> bonds;
	std::vector<RoboAngle> angles;
	std::vector<RoboPeriodicTorsion> properPeriodicTorsions;
	std::vector<RoboHarmonicImproperTorsion> harmonicImproperTorsions;

	int numMolecules = 0;

	std::vector<Topology> topologies;
	std::vector<int> roots;
	
	uint32_t seed = 0;
	int numThreads = 0;
	NonbondedMethod nonbondedMethod = NonbondedMethod::NoCutoff;
	SimTK::Real nonbondedCutoffInNm = 1.2; // 1.2 nm, not used by default (no cutoff)
	SimTK::Real vdwGlobalScaleFactor = 1.0; // Default is 1.0

	bool useGBSAOBC2 = false;
	SimTK::Real solventDielectric = 78.5;
	SimTK::Real soluteDielectric = 1.0;
	SimTK::Real gbsaGlobalScaleFactor = 0.0; // Default is 0 (vacuum). Use 1 for implicit solvent (water)

	bool useAmberForceFieldScaleFactors = true;
	SimTK::Real globalForceFieldScaleFactor = 1.0; // Used in place of Amber scaling (not default)

	// Random number generator
	Random32 randomEngine;

public:
	// /** Implicit membrane mimicked by half-space contacts */
	// void addContactImplicitMembrane(const float memZWidth, const SetupReader& setupReader);

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
        "Slider",
		"OrthoSpherical"
    };

	SimTK::BondMobility::Mobility getMobility(const std::string& mobilityStr) {
		// Assume MobilityStr is a vector defined elsewhere in your code
		// std::vector<std::string> MobilityStr = { ... };

		auto it = std::find(MobilityStr.begin(), MobilityStr.end(), mobilityStr);
		
		if (it != MobilityStr.end()) {
			// If the string is found, return the corresponding enum value
			return static_cast<SimTK::BondMobility::Mobility>(std::distance(MobilityStr.begin(), it) + 1);
		} else {
			// If the string is not found, return the default value
			return SimTK::BondMobility::Default;
		}
	}	

	
	/**	
	* @brief Get Z-matrix indexes table
	* @param
	*/
	void calcZMatrixTable(void);

	/**	
	* @brief Allocate Z Matrix BAT
	* @param
	*/
	void reallocZMatrixBAT(void);

	
	/**	
	* @brief Print function for the zMatrixTable
	* @param
	*/
	void PrintZMatrixTable() const ;

	
	/**	
	* @brief Function to print the zMatrixBAT
	* @param
	*/
	void PrintZMatrixBAT() const ;

	/**	
	* @brief
	* @param
	*/
	void PrintZMatrixTableAndBAT() const;


private:
	ZMatrix zMatrix;

	std::vector<SimTK::Compound::AtomTargetLocations> atomTargetLocationsCache;

    /** @name Z Matrix and BAT functions
	*/

    /**@{**/
	
	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

	std::vector<std::vector<int>> zMatrixTable;
	std::vector<std::vector<SimTK::Real>> zMatrixBAT;

	/**	
	* @brief Add a new row to the zMatrixTable
	* @param
	*/
    void addZMatrixTableRow(const std::vector<int>& newRow) ;

	/**	
	* @brief Getter for a specific entry
	* @param
	* @return
	*/
    int getZMatrixTableEntry(int rowIndex, int colIndex) const ;

	/**	
	* @brief Setter for a specific entry
	* @param
	*/
    void setZMatrixTableEntry(int rowIndex, int colIndex, int value) ;

	/**	
	* @brief Setter for a specific entry
	* @param
	*/
    void setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value) ;

	/**	
	* @brief Function to get a given row
	* @param
	*/
    const std::vector<SimTK::Real>& getZMatrixBATRow(size_t rowIndex) const;

	/**	
	* @brief Function to get a given row
	* @param
	* @return
	*/
    std::vector<SimTK::Real>& updZMatrixBATRow(size_t rowIndex) ;



	/**	
	* @brief 
	* @param
	*/
	void calcZMatrixBAT(int wIx,
		const std::vector< std::vector<
			std::pair <RoboAtom *, SimTK::Vec3 > > >&
				otherWorldsAtomsLocations);

	/**	
	* @brief Function to get the value for a given row and column in zMatrixBAT
	* @param
	* @return
	*/
    SimTK::Real getZMatrixBATValue(size_t rowIndex, size_t colIndex) const ;

	/**	
	* @brief Function to add a new row to the zMatrixBAT
	* @param
	*/
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

	// Explicit pairs (replica_i, replica_j)
	std::vector<std::pair<int, int>> exchangePairList;

	// Quick lookup: exchangePairs[i] = j or -1
    std::vector<int> exchangePairs;

	std::map<std::string, AcceptRejectMode> acceptRejectModes = {
		{ "EMPTY", AcceptRejectMode::AlwaysAccept },
		{ "MC", AcceptRejectMode::MetropolisHastings },
	};

	std::map<std::string, IntegratorType> integratorTypes = {
		{ "EMPTY", IntegratorType::EMPTY },
		{ "VV", IntegratorType::VERLET },
		{ "VERLET", IntegratorType::VERLET },
		{ "EULER", IntegratorType::EULER },
		{ "EULER2", IntegratorType::EULER2 },
		{ "CPODES", IntegratorType::CPODES },
		{ "RUNGEKUTTA", IntegratorType::RUNGEKUTTA },
		{ "RUNGEKUTTA2", IntegratorType::RUNGEKUTTA2 },
		{ "RUNGEKUTTA3", IntegratorType::RUNGEKUTTA3 },
		{ "RUNGEKUTTAFELDBERG", IntegratorType::RUNGEKUTTAFELDBERG },
		{ "BENDSTRETCH", IntegratorType::BENDSTRETCH },
		{ "OMMVV", IntegratorType::OMMVV },
		{ "BOUND_WALK", IntegratorType::BOUND_WALK },
		{ "BOUND_HMC", IntegratorType::BOUND_HMC },
		{ "STATIONS_TASK", IntegratorType::STATIONS_TASK },
		{ "NOF_INTEGRATORS", IntegratorType::NOF_INTEGRATORS },
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

	void printQStats(int thIx);

	//////////////////////////////////
	//---         Q Stats        -----
	//////////////////////////////////

};
