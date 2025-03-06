#pragma once

#include "Robo.hpp"
#include "Sampler.hpp"
#include "SetupReader.hpp"
#include "World.hpp"
#include "ThermodynamicState.hpp"
#include "Replica.hpp"
#include "SetupReader.hpp"
#include "TrajectoryObject.hpp"

#include "OpenMM.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/common/windowsExportCommon.h"

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

const std::vector<std::string> RUN_TYPE_Str = {"DEFAULT", "REMC", "RENEMC", "RENE"};

const std::unordered_map<std::string, RUN_TYPE>
RUN_TYPE_MAP{
	{"DEFAULT", RUN_TYPE::DEFAULT},
	{"REMC", RUN_TYPE::REMC},
	{"RENEMC", RUN_TYPE::RENEMC},
	{"RENE", RUN_TYPE::RENE}
};

//==============================================================================
//                           CLASS Context
//==============================================================================
/** 
 * This defines the Context class.
 
Topology 0:                                                      :       
          :                                                      :
Position 0:                                                      :
          :          ┌─────────────────────────┐                 :
                     │                         │                 :
                     │       INITIALIZE        │                 :
                     │                         │                 :
                     └────────────┬────────────┘                 : 
                                  │                              :
                                  │                              :
                            ┌─────▼─────┐                        :
                            │    RUN    │                        :     
                            └─────┬─────┘                        :
                                  │                              :
**/
class Context{

    /** @name Constructor **/
    /**@{**/
	/**@}**/

	std::string baseName;
	bool verbose = true;

public:
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
	Context(const std::string& baseName, 
			uint32_t seed,
			uint32_t threads,
			uint32_t nofRoundsTillReblock,
			RUN_TYPE runType,
			uint32_t swapFreq,
			uint32_t swapFixmanFreq);

	/**
	 * @brief Print atom's Compound and DuMM indexes.
	*/
	void PrintAtomsDebugInfo(void);

	void setVerbose(bool verbose);

    /** @name System and simulation setup.
	 * 1. Read input file.
	 * 2. Read Amber files and construct topologies.
	 * 3. Add worlds.
	 * 4. Add samplers to worlds.
	 * 5. Read REX input and build replicas/thermostates. **/
    /**@{**/

	/**	
	* @brief Read input file.
	* @param var
	* @return
	*/
	bool setOutput(const std::string& outDir);
	bool CheckInputParameters(const SetupReader& setupReader);
	std::string GetMoleculeDirectoryShort(const std::string& path) const;
	
	void setNofRoundsTillReblock(int nofRoundsTillReblock);
	void setRequiredNofRounds(int argNofRounds);

	void setNonbonded(int method, SimTK::Real cutoff);

	RUN_TYPE getRunType(void) const;
	void setRunType(RUN_TYPE runTypeArg);
	RUN_TYPE setRunType(const std::string& runTypeArgStr);



	/**	
	* @brief Read Amber files.
	* @param var
	*/
	void loadAmberSystem(const std::string& prmtop, const std::string& inpcrd);

	/**
	 * @brief Set Molmodel atom masses from out list of atoms
	*/
	void setAtomMasses();

	/**
	 * @brief Calc Z matrix and BAT load replicas and nonequil
	*/
	void Initialize();

	/**	
	* @brief Construct topologies.
	* @param var
	*/


	//void constructTopologies

	/**	
	* @brief Add worlds.
	* @param var
	* @return
	*/

 	std::vector<BOND_FLEXIBILITY>& readFlexibility(
		std::string flexFileFN,
		std::vector<BOND_FLEXIBILITY>& flexibilities);
		
	void addWorld(
		bool fixmanTorque,
		int samplesPerRound,
		ROOT_MOBILITY rootMobility,
		const std::vector<BOND_FLEXIBILITY>& flexibilities,
		bool useOpenMM = true,
		bool visual = false,
		SimTK::Real visualizerFrequency = 0);

	/**
	 * @brief Load coordinates into atoms, add replicas and pass coordinates.
	 * @param prmtop Amber prmtop file.
	 * @param restartDir Restart directory.
	 * @param nofReplicas Number of replicas.
	*/
	bool addReplicasAndLoadCoordinates(const std::string& name, const std::string& prmtop, const std::string& restartDir, int nofReplicas);

	/**	
	* @brief
		* Initialize Context 
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
		* 12. Constraints
	* @param filename input file name
	* @return succes of the function
	*/	
	bool initializeFromFile(const std::string& filename);


	/**@}**/


	void setNumThreads(int threads);
	void setGBSA(SimTK::Real globalScaleFactor);
	void setForceFieldScaleFactors(SimTK::Real globalScaleFactor);

	//void setRootMobilitiesFromFlexFiles(void);

	// Input functions
	bool loadTopologyFile(std::string topologyFilename);
	bool loadCoordinatesFile(std::string coordinatesFilename);
	void PrintCoordinates(const std::vector<std::vector
        <std::pair <bSpecificAtom *,
		SimTK::Vec3>>>& atomsLocations);

	bool loadRigidBodiesSpecs(std::size_t whichWorld, int whichMolecule, std::string RBSpecsFN);
	bool loadFlexibleBondsSpecs(std::size_t whichWorld, std::string FlexSpecsFN);
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

	/**  */
	void new_build_GmolGraph_MolmodelAcyclicGraph(); // TODO delete

	/**	
	* @brief Calculate BAT graphs
	*/
	void calc_Gmolmodel_Graph();

	/**	
	* @brief Build Molmodel Compound / graphs
	*/
	void build_Molmodel_AcyclicGraphs();


	/** @brief __fill__ */
	void generateTopologiesSubarrays(void);

	/** Assign Compound coordinates by matching bAtomList coordinates */
	void matchDefaultConfigurationFromAtomsCoords(Topology& topology, int molIx);

	void matchDefaultConfiguration(int molIx, std::map<Compound::AtomIndex, SimTK::Vec3> atomTargets){assert(!"Not implemented");}

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
	void generateDummAtomClasses(AmberReader& amberReader);
	void bAddDummBondParams(AmberReader& amberReader);
	void bAddDummAngleParams(AmberReader& amberReader);
	bool checkBond(int a1, int a2);
	void bAddDummTorsionParams(AmberReader& amberReader);

	void addDummParams(AmberReader& amberReader);

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

	void passTopologiesToNewWorld(int newWorldIx);

	int getNofMolecules();
	//------------

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
	void transferCoordinates_WorldToWorld(int src, int dest);
	SimTK::Real checkTransferCoordinates_Cart(int srcWIx, int destWIx);
	SimTK::Real checkTransferCoordinates_BAT(int srcWIx, int destWIx, bool wantJacobian = false);

	void transferCoordinates_ReplicaToWorld(int replicaIx, int destWIx);

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

    /** @name Write coordinates **/
    /**@{**/

	/**	
	* @brief Write pdb file.
	*/
	void writeInitialPdb(void);
	void writeFinalPdb(void);
	void writePdb(std::size_t whichWorld);
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

	// Write dcd
	void writeDCDs();

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
	void addReplica(int index);

	// Add one thermodynamic state
	void addThermodynamicState(
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
	void swapReferencePotentialEnergies(int replica_i, int replica_j);

	// Exchanges thermodynamic states between replicas
	void getMsg_RexDetHeader(std::stringstream& rexDetHeader);
	void rewindReplica(void);
	bool attemptREXSwap(int replica_i, int replica_j);

	const int getSwapEvery(void);
	void setSwapEvery(const int& n);

	// StartingFrom argument is for alternating odd and even neighbors
	void mixNeighboringReplicas(unsigned int startingFrom);

	// Mix replicas
	void mixAllReplicas(int nSwapAttempts);
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
	void RunWorlds(std::vector<int>& specificWIxs, int replicaIx);
	void RunReplicaRefactor_SIMPLE(int mixi, int replicaIx);	
	/**	
	* @brief Main function
	* @param
	* @return
	*/
	void RunREX(int equilRounds, int prodRounds);

	void Run();

	/**@}**/


	void setSubZmatrixBATStatsToSamplers(int thermoIx, int worldCnt);

	// Transfer Q statistics
	void transferQStatistics(int thermoIx, int srcStatsWIx, int destStatsWIx);

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

    /**	
	* @brief Initialize OpenMM with general data
	* @param allowReferencePlatform Allow OpenMM reference platform
	* @return OpenMM context's platform's name
	*/ 
	std::string OMMRef_initialize(void);

	SimTK::Real OMMRef_calcPotential(
		const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>> & atomsLocations,
		bool wantEnergy, bool wantForces);

private:

	//OpenMMPluginInterface refOpenMMPlugin; // __refOMM__
	std::unique_ptr<OpenMM::Platform> platform; // __refOMM__
	std::unique_ptr<OpenMM::Context> openMMContext; // __refOMM__
	std::unique_ptr<OpenMM::System> openMMSystem; // __refOMM__

	std::unique_ptr<OpenMM::NonbondedForce> ommNonbondedForce; // __refOMM__
	//std::unique_ptr<OpenMM::GBSAOBCForce> ommGBSAOBCForce; // __refOMM__
	std::unique_ptr<OpenMM::HarmonicBondForce> ommHarmonicBondStretch; // __refOMM__
	std::unique_ptr<OpenMM::HarmonicAngleForce> ommHarmonicAngleForce; // __refOMM__
	std::unique_ptr<OpenMM::PeriodicTorsionForce> ommPeriodicTorsionForce; // __refOMM__

	std::unique_ptr<OpenMM::AndersenThermostat> openMMThermostat;
	std::unique_ptr<OpenMM::Integrator> openMMIntegrator;

	mutable OpenMM::State openMMState;

	// end __refOMM__

	std::vector<int> TopologyIXs;
	std::vector<std::vector<int>> AmberAtomIXs;
	std::vector<World> worlds;
	std::vector<AmberReader> amberReader;

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

	int natoms = std::numeric_limits<int>::min();
	std::vector<bSpecificAtom> atoms;

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

	std::vector<int> findMolecules(const AmberReader& reader);

	void loadAtomsCoordinates(const std::string& prmtop, const std::string& inpcrdFN);
	
	void loadAtoms(const AmberReader& reader);
	void loadBonds(const AmberReader& reader);
	void loadAngles(const AmberReader& reader);
	void loadTorsions(const AmberReader& reader);

	void setAtomsCompounds();
	void addBiotypes();
	std::vector<bSpecificAtom>& getAtoms() {
        return atoms;
    }

	uint32_t seed = 0;
	int numThreads = 0;
	int nonbondedMethod = 0; // 0 = NoCutoff, 1 = CutoffNonPeriodic, 2 = CutoffPeriodic .. TODO: implement enum
	SimTK::Real nonbondedCutoff = 1.2; // 1.2 nm, not used by default (no cutoff)
	SimTK::Real vdwGlobalScaleFactor = 1.0; // Default is 1.0
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

	std::vector<std::string> inpcrdFNs;

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
			std::pair <bSpecificAtom *, SimTK::Vec3 > > >&
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

	// Pairs of replicas to be exchanged
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

	void calcQStats(int thIx);

	void printQStats(int thIx);

	//////////////////////////////////
	//---         Q Stats        -----
	//////////////////////////////////

};
