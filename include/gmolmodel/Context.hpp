#ifndef __CONTEXT_HPP__
#define __CONTEXT_HPP__

#include "Robo.hpp"
#include "Sampler.hpp"
#include "SetupReader.hpp"
#include "World.hpp"

class Sampler;
class World;

class Context{

public:
	// Context(const SetupReader& setupReader, World *, std::string logFilenameArg);
	Context(const SetupReader& setupReader, std::string logFilenameArg);
	~Context();

	void printStatus(void);
	void throwAndExit(std::string errMsg, int errCode);

	World * AddWorld(bool visual, SimTK::Real visualizerFrequency = 0.0015);
	//World * AddWorld(World *, bool visual);

	World * getWorld();
	World * getWorld(std::size_t which);

	World * updWorld();
	World * updWorld(std::size_t which);

	// Returns the size of the worlds vector
	std::size_t getNofWorlds() const;

	SimTK::DuMMForceFieldSubsystem * updForceField(std::size_t whichWorld);

	// Writeble reference to a samplers advanced state
	SimTK::State& updAdvancedState(std::size_t whichWorld, std::size_t whichSampler);

	// --- Use a SetupReader Object to read worlds information from a file ---
	bool loadTopologyFile(/*std::size_t whichWorld, int whichMolecule,*/ std::string topologyFilename);
	bool loadCoordinatesFile(/*std::size_t whichWorld, int whichMolecule,*/ std::string coordinatesFilename);
	bool loadRigidBodiesSpecs(std::size_t whichWorld, int whichMolecule, std::string RBSpecsFN);
	bool loadFlexibleBondsSpecs(std::size_t whichWorld, int whichMolecule, std::string FlexSpecsFN);
	void setRegimen (std::size_t whichWorld, int whichMolecule, std::string regimen);

	/** Load molecules based on loaded filenames **/
	void AddMolecules(std::vector<std::string> argRoots,
		std::vector<std::string> argRootMobilities);
	void modelTopologies(std::vector<std::string> GroundToCompoundMobilizerTypes);

	void realizeTopology();

	void LoadWorldsFromSetup(SetupReader&);

	int getNofMolecules();
	//------------

	// --- Thermodynamics ---a
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
	BaseSampler* addSampler(std::size_t whichWorld, SamplerName whichSampler);
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
	int getNofRounds();
	void setNofRounds(int nofRounds);

	int getNofSamplesPerRound(std::size_t whichWorld);
	void setNofSamplesPerRound(std::size_t whichWorld, int MCStepsPerRound);

	std::size_t getWorldIndex(std::size_t which) const;

	// --- Arrange different mixing parameters ---
	void initializeMixingParamters();
	//------------

	// --- Mix ---
	void RotateWorlds();
	//------------

	// --- Main ---
	void Run(SetupReader&);
	void Run(int howManyRounds, SimTK::Real Ti, SimTK::Real Tf);
	void RunSimulatedTempering(int howManyRounds, SimTK::Real Ti, SimTK::Real Tf);
	void setNofBoostStairs(std::size_t whichWorld, int howManyStairs);
	int getNofBoostStairs(std::size_t whichWorld);
	void setNumThreadsRequested(std::size_t which, int howMany);


	void setUseOpenMMAcceleration(bool arg);
    void setUseOpenMMCalcOnlyNonBonded(bool arg);
    void setNonbondedMethod(std::size_t whichWorld, int methodInx);
    void setNonbondedCutoff(std::size_t whichWorld, Real cutoffNm);

	SimTK::Real Pearson(std::vector<std::vector<SimTK::Real>> someVector, int QIx1, int QIx2); // 2D roundsTillReblock; 3D nofQs

	/** Print the number of threads each World got **/
	void PrintNumThreads();

	/** Get/Set seed for reproducibility. **/
	void setSeed(std::size_t whichWorld, std::size_t whichSampler, uint32_t seed);
	uint32_t getSeed(std::size_t whichWorld, std::size_t whichSampler) const;

	//------------

	/** Analysis related functions **/
	void addDistance(std::size_t whichWorld, std::size_t whichCompound, int aIx1, int aIx2);
	void addDihedral(std::size_t whichWorld, std::size_t whichCompound, int aIx1, int aIx2, int aIx3, int aIx4);

	// --- Printing functions ---

	// Print Molmodel related information
	void PrintMolmodelAndDuMMTypes(void);

	// Print Simbody related information
	void PrintSimbodyMobods(void);
	void PrintSamplerData(std::size_t whichWorld);
	void PrintGeometry(SetupReader&, std::size_t whichWorld);
	void PrintGeometry(std::size_t whichWorld);
	void PrintDistances(std::size_t whichWorld);
	void PrintDihedrals(std::size_t whichWorld);
	void PrintDihedralsQs(std::size_t whichWorld);
	void PrintFreeE2EDist(std::size_t whichWorld, int whichCompound);
	void WritePdb(std::size_t whichWorld);
	SimTK::Real Dihedral(std::size_t whichWorld, std::size_t whichCompound, std::size_t whichSampler, int a1, int a2, int a3, int a4);
	SimTK::Real Distance(std::size_t whichWorld, std::size_t whichCompound, std::size_t whichSampler, int a1, int a2);

	int getPdbRestartFreq();
	void setPdbRestartFreq(int argFreq);
	int getPrintFreq();
	void setPrintFreq(int argFreq);

	std::string getOutputDir();
	void setOutputDir(std::string arg);
	std::string getPdbPrefix();
	void setPdbPrefix(std::string arg);
	//------------

public:
	std::vector<std::size_t> worldIndexes;

private:
	void ValidateSetupReader(const SetupReader& setupReader);

	std::vector<World> worlds;

	// Molecules files
	std::vector<std::string> topFNs;
	std::vector<std::string> crdFNs;
	std::vector<std::vector<std::string>> rbSpecsFNs;
	std::vector<std::vector<std::string>> flexSpecsFNs;
	std::vector<std::vector<std::string>> regimens;
	std::vector<std::string> rootMobilities;

	// Simulation parameters
	int nofRounds;
	//int total_mcsteps;

	std::size_t nofWorlds;
	bool isWorldsOrderRandom;
	std::vector<int> nofSamplesPerRound;
	std::vector<int> nofMDStepsPerSample;
	std::vector<SimTK::Real> timesteps;

	std::vector<int> nofBoostStairs;

	std::size_t nofMols;
	std::size_t nofTopologies; // nofWorlds x nofMols
	//
	bool reproducible;
	int pdbRestartFreq;
	int printFreq;

	// 
	std::string outputDir;
	std::string pdbPrefix;

	// Geometric features analysis
	// First two integers specifiy the world and the Compound. The rest
	// specifies atom indeces
	std::vector< std::vector<int> > distanceIxs;
	std::vector< std::vector<int> > dihedralIxs;

	SimTK::Real geom1[PRINT_BUFFER_SIZE];
	SimTK::Real geom2[PRINT_BUFFER_SIZE];
	SimTK::Real geom3[PRINT_BUFFER_SIZE];

	// Output
	unsigned int BUFSIZE;
	std::unique_ptr<char[]> buffer;
	FILE *logFile;

	// Adaptive Gibbs blocking variables
	int roundsTillReblock;
	std::vector<std::vector<std::vector<SimTK::Real>>> QsCache; // 1D nofWorlds; 2D roundsTillReblock; 3D nofQs

	// Normal mode analysis
	std::vector<int> NMAOption;

};

#endif //__CONTEXT_HPP__





