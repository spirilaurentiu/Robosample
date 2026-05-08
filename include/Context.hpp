#pragma once

#include "OpenMM.hpp"
#include "Replica.hpp"
#include "Sampler.hpp"
#include "ThermodynamicState.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"
#include "bgeneral.hpp"

class Context {
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
    Context(const std::string& baseName,
            uint32_t seed,
            uint32_t nofRoundsTillReblock,
            RUN_TYPE runType,
            uint32_t swapFreq,
            uint32_t swapFixmanFreq,
            bool testing);

    void setVerbose(bool verbose);
    void setGBSAOptions(bool useGBSAOBC2, SimTK::Real solventDielectric, SimTK::Real soluteDielectric);
    bool setOutput(const std::string& outDir);

    void setNofRoundsTillReblock(int nofRoundsTillReblock);
    void setRequiredNofRounds(int argNofRounds);

    void setNonbonded(NonbondedMethod method, SimTK::Real cutoffInNm);

    void loadAmberSystem(const SystemTopology& systemTopology,
                         const ForceFieldParams& ffParams,
                         const SimulationSettings& simSettings,
                         const ZMatrix& zMatrix);
    auto initializeOpenMM() -> bool;

    SimTK::Real calculatePotentialEnergy(int worldIndex);

    auto getAtomNameByPrmtopIndex(int prmtopIndex) const -> const std::string& {
        for (const auto& atom : systemTopology.atoms) {
            if (atom.identity.prmtopIndex == prmtopIndex) {
                return atom.identity.uniqueAtomName;
            }
        }
        throw std::runtime_error("Atom with specified prmtop index not found.");
    }

    auto validateContext() -> bool;

    void addWorld(bool fixmanTorque,
                  int samplesPerRound,
                  const std::vector<std::vector<BondFlexibility>>& flexibilities);

    // Add task spaces
    void addTaskSpacesLS();

    /** Add constraints */
    void addConstraints();

    void passTopologiesToNewWorld(int newWorldIx);

    // --- Simulation parameters ---

    //------------

    // --- Mixing parameters ---
    // Another way to do it is setting the number of rounds
    int getRequiredNofRounds();

    int getNofRoundsTillReblock();

    void updNofRoundsTillReblock(int nofRoundsTillReblock);

    // Adaptive Gibbs blocking: TODO: consider moving in World
    void allocateReblockQsCache();
    void allocateReblockQsCacheQVectors();

    // --- Arrange different mixing parameters ---
    void initializeMixingParameters();
    //------------

    auto getNofWorlds() const -> std::size_t {
        return worlds.size();
    }

    auto getWorld(std::size_t which) -> World& {
        return worlds[which];
    }

    auto getWorld(std::size_t which) const -> const World& {
        return worlds[which];
    }

    auto getWorlds() -> std::vector<World>& {
        return worlds;
    }

    auto getWorlds() const -> const std::vector<World>& {
        return worlds;
    }

    // Relationship BAT - mobod transforms
    void PrintZMatrixMobods(int wIx, SimTK::State& someState);

    auto Pearson(std::vector<std::vector<SimTK::Real>> someVector,
                 int QIx1,
                 int QIx2) -> SimTK::Real; // 2D roundsTillReblock; 3D nofQs

    //------------

    /** Analysis related functions **/
    void addDistance(std::size_t whichWorld, std::size_t whichCompound, std::size_t aIx1, std::size_t aIx2);
    void addAngle(std::size_t whichWorld,
                  std::size_t whichCompound,
                  std::size_t aIx1,
                  std::size_t aIx2,
                  std::size_t aIx3);
    void addDihedral(std::size_t whichWorld,
                     std::size_t whichCompound,
                     std::size_t aIx1,
                     std::size_t aIx2,
                     std::size_t aIx3,
                     std::size_t aIx4);
    void addDistances(const std::vector<std::size_t>& distanceIx);
    void addAngles(const std::vector<std::size_t>& angleIx);
    void addDihedrals(const std::vector<std::size_t>& dihedralIx);

    // --- Output ---
    void printThermodynamics();

    // Print DuMM atoms stations in mobilized body frame
    void checkAtomStationsThroughDumm();

    /**
     * @brief Write pdb file.
     */
    void writeInitialPdb();
    void writeFinalPdb();
    void writePdbs(int someIndex, int thermodynamicStateIx = 0);

    // Output helpers
    auto getPdbRestartFreq() -> int;
    void setPdbRestartFreq(int argFreq);

    auto getRestartDir() const -> const std::string&;
    void setRestartDir(const std::string& argRestartDir);

    void setPdbPrefix(const std::string& argPdbPrefix);
    auto getPdbPrefix() -> std::string;

    auto getOutputDir() -> std::string;
    void setOutputDir(std::string arg);

    //////////////////////////////////
    //// REPLICA EXCHANGE FUNCTIONS //
    //////////////////////////////////

    void allocateSwapMatrices();

    // Add one replica
    void addReplica();

    // Add one thermodynamic state
    void addThermodynamicState(SimTK::Real T,
                               const std::vector<AcceptRejectMode>& acceptRejectModes,
                               const std::vector<int>& rexDistortOptions,
                               const std::vector<std::string>& rexDistortArgs,
                               const std::vector<int>& rexFlowOptions,
                               const std::vector<int>& rexWorkOptions,
                               const std::vector<IntegratorType>& rexIntegrators,
                               const std::vector<int>& argWorldIndexes,
                               const std::vector<SimTK::Real>& timestepsInThisReplica,
                               const std::vector<int>& mdstepsInThisReplica);

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
    void PrepareNonEquilibriumParams_Q();

    // Set thermodynamic states nonequilibrium flag
    void setThermostatesNonequilibrium();

    // Set the initial mapping between replicas and thermoStates
    void loadReplica2ThermoIxs();

    // Get Fixman potential already calculated from replica
    SimTK::Real getFixman(int replica_i);

    // Calculate Fixman potential of replica I in replica J's back world. Uj(X_i)
    SimTK::Real calcFixman_IInJ(int replica_i, int replica_j);

    const int& getSwapFixman() {
        return swapFixman;
    }
    void setSwapFixman(const int argSwapFixman) {
        swapFixman = argSwapFixman;
    }

    SimTK::Real calcReplicaTransferredEnergy(int replicaIx);
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

    int getSwapEvery() const {
        return swapEvery;
    }
    void setSwapEvery(int n) {
        swapEvery = n;
    }

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
    SimTK::Real perturbScalingFactor(std::string how, SimTK::Real sf, bool randSignOpt = false);

    // Set world distort parameters
    void setWorldDistortParameters(int whichWorld, SimTK::Real scaleFactor);

    // Set nonequilibrium parameters for one replica
    void updWorldsDistortOptions(int thisReplica);
    void updThermostatesQScaleFactors(int mixi);

    void writeLog(int mixi, int replicaIx);
    void writeDCD(int replicaIx);

    void incrementNofSamples();

    /** @name Replica exchange **/
    /**@{**/

    // Run a particular world
    auto RunWorld(int whichWorld, const std::string& header, bool shouldPrint) -> bool;
    void RunReplicaWorldRange(int replicaIx,
                              int startWorldCnt,
                              int nofWorldsCounted,
                              bool isNonEquilibrium,
                              bool shouldPrint);
    /**
     * @brief Main function
     * @param
     * @return
     */
    void RunREX(int numEquilibrationRounds, int numProductionRounds, int writeFrequency, bool writeToStdio);

    void setSubZmatrixBATStatsToSamplers(int thermoIx, int worldCnt);

    // Transfer Q statistics
    void transferQStatistics(int thermoIx, int srcStatsWIx, int destStatsWIx);

    void PrintReplicas();
    void PrintReplicaMaps();
    void PrintNofAttemptedSwapsMatrix();
    void PrintNofAcceptedSwapsMatrix();

    // Transformers
    void Print_TRANSFORMERS_Work();

    std::vector<SimTK::Real> AtomIndex0, AtomIndex1, UDotCache, UCache;
    bool binaryFileIsInitialized = false;

    std::string foutU, foutUDot, foutTorque;

    void initializeBinaryFile(const std::string& filename, uint32_t num_columns);
    void writeRowToBinaryFile(const std::string& filename,
                              const std::vector<SimTK::Real>& row,
                              bool has_acceptance,
                              bool accepted);

    private:
    // Run in testing mode
    bool testing = false;

    std::vector<World> worlds;
    std::vector<std::size_t> worldIndices;

    int moleculeCount = -1;

    // Simulation parameters
    int requiredNofRounds = -1;

    std::size_t nofWorlds = 0;

    int pdbRestartFreq = 0;

    std::string molDir;
    std::string outputDir;
    std::string restartDir;
    std::string pdbPrefix;

    // Geometric features analysis
    // First two integers specify the world and the Compound. The rest
    // specifies atom indices
    std::vector<std::vector<std::size_t>> distanceIxs;
    std::vector<std::vector<std::size_t>> angleIxs;
    std::vector<std::vector<std::size_t>> dihedralIxs;

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

    int nofReplicas = 0;
    int nofThermodynamicStates = 0;
    ReplicaMixingScheme replicaMixingScheme = ReplicaMixingScheme::Neighboring;

    int swapFixman = 1;
    int swapEvery = 1;

    std::uniform_real_distribution<SimTK::Real> uniformRealDistribution =
        std::uniform_real_distribution<SimTK::Real>(SimTK::Zero, SimTK::One);

    // Non-equilibrium parameters
    std::vector<SimTK::Real> qScaleFactorsEven;
    std::vector<SimTK::Real> qScaleFactorsOdd;
    std::vector<SimTK::Real> qScaleFactorsMiu;
    std::vector<SimTK::Real> qScaleFactorsStd;

    std::vector<SimTK::Real> qScaleFactors;

    std::string cerr_prefix = "[ERROR] ";
    std::string cwar_prefix = "[WARNING] ";
    std::string cinf_prefix = "[INFO] ";

    RUN_TYPE runType = RUN_TYPE::Default;
    SimTK::Real tempIni = 0, tempFin = 0;

    SystemTopology systemTopology;
    ForceFieldParams ffParams;
    SimulationSettings simSettings;
    ZMatrix zMatrix;

    int numMolecules = 0;

    std::vector<Topology> topologies;

    uint32_t seed = 0;
    NonbondedMethod nonbondedMethod = NonbondedMethod::NoCutoff;
    SimTK::Real nonbondedCutoffInNm = 1.2;  // 1.2 nm, not used by default (no cutoff)
    SimTK::Real vdwGlobalScaleFactor = 1.0; // Default is 1.0

    bool useGBSAOBC2 = false;
    SimTK::Real solventDielectric = 78.5;
    SimTK::Real soluteDielectric = 1.0;
    SimTK::Real gbsaGlobalScaleFactor = 0.0; // Default is 0 (vacuum). Use 1 for implicit solvent (water)

    // Random number generator
    Random32 randomEngine;

    std::vector<SimTK::Real> dcdXBuffer;
    std::vector<SimTK::Real> dcdYBuffer;
    std::vector<SimTK::Real> dcdZBuffer;

    public:
    // /** Implicit membrane mimicked by half-space contacts */
    // void addContactImplicitMembrane(const float memZWidth, const SetupReader& setupReader);

    /**
     * @brief Get Z-matrix indexes table
     * @param
     */
    void calcZMatrixTable();

    /**
     * @brief Allocate Z Matrix BAT
     * @param
     */
    void reallocZMatrixBAT();


    /**
     * @brief Print function for the zMatrixTable
     * @param
     */
    void PrintZMatrixTable() const;


    /**
     * @brief Function to print the zMatrixBAT
     * @param
     */
    void PrintZMatrixBAT() const;

    /**
     * @brief
     * @param
     */
    void PrintZMatrixTableAndBAT() const;


    private:
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
    void addZMatrixTableRow(const std::vector<int>& newRow);

    /**
     * @brief Getter for a specific entry
     * @param
     * @return
     */
    int getZMatrixTableEntry(int rowIndex, int colIndex) const;

    /**
     * @brief Setter for a specific entry
     * @param
     */
    void setZMatrixTableEntry(int rowIndex, int colIndex, int value);

    /**
     * @brief Setter for a specific entry
     * @param
     */
    void setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value);

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
    std::vector<SimTK::Real>& updZMatrixBATRow(size_t rowIndex);


    /**
     * @brief
     * @param
     */
    void calcZMatrixBAT(
        int wIx,
        const std::vector<std::vector<std::pair<RoboAtom*, SimTK::Vec3>>>& otherWorldsAtomsLocations);

    /**
     * @brief Function to get the value for a given row and column in zMatrixBAT
     * @param
     * @return
     */
    SimTK::Real getZMatrixBATValue(size_t rowIndex, size_t colIndex) const;

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
    SimTK::Real calcInternalBATJacobianLog();

    /**
     * @brief zmatrixbat_ Add BAT coordinates
     * @param
     */
    void addSubZMatrixBATsToWorld(int wIx, int replicaIx);

    /**
     * @brief zmatrixbat_ Get BAT coordinates modifiable by a selected world
     * @param
     */
    void updSubZMatrixBATsToWorld(int wIx, int replicaIx);

    /**
     * @brief zmatrixbat_ Set BAT coordinates modifiable to all worlds of a replica
     * @param
     */
    void updSubZMatrixBATsToAllWorlds(int replicaIx);

    /**
     * @brief zmatrixbat_ Print BAT coordinates
     * @param
     */
    void PrintWorldSubZMatrixBATs(int wIx);

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

    //////////////////////////////////
    //---         Q Stats        -----
    //////////////////////////////////
    void reserveThermostatsQs();

    void setThermostatesQs();

    void printQStats(int thIx);

    //////////////////////////////////
    //---         Q Stats        -----
    //////////////////////////////////
};
