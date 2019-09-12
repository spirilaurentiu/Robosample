#ifndef __CONTEXT_HPP__
#define __CONTEXT_HPP__

#include "Robo.hpp"
#include "Sampler.hpp"

class Sampler;
class World;

class Context{

public:
    Context(World *, std::string logFilenameArg);
    Context(std::string logFilenameArg);
    ~Context();

    World * AddWorld(bool visual, SimTK::Real visualizerFrequency = 0.0015);
    //World * AddWorld(World *, bool visual);

    World * getWorld(void) const;
    World * getWorld(int which) const;

    World * updWorld(void);
    World * updWorld(int which);

    unsigned int getNofWorlds(void);

    SimTK::DuMMForceFieldSubsystem * updForceField(int whichWorld);

    // Writeble reference to a samplers advanced state
    SimTK::State& updAdvancedState(int whichWorld, int whichSampler);

    // --- Use a SetupReader Object to read worlds information from a file ---
    bool loadTopologyFile(int whichWorld, int whichMolecule, std::string topologyFilename);
    bool loadCoordinatesFile(int whichWorld, int whichMolecule, std::string coordinatesFilename);
    bool loadRigidBodiesSpecs(int whichWorld, int whichMolecule, std::string RBSpecsFN);
    bool loadFlexibleBondsSpecs(int whichWorld, int whichMolecule, std::string FlexSpecsFN);
    void setRegimen (int whichWorld, int whichMolecule, std::string regimen);

    /** Load molecules based on loaded filenames **/
    void AddMolecules(std::vector<std::string> argRoots);
    void modelTopologies(void);

    void realizeTopology(void);

    void LoadWorldsFromSetup(SetupReader&);

    int getNofMolecules(void);
    //------------

    // --- Thermodynamics ---a
    // Get/set the main temperature (acc/rej temperature for MC)
    float getTemperature(int whichWorld);
    void  setTemperature(int whichWorld, float someTemperature);

    // If HMC, get/set the guidance Hamiltonian temperature
    float getGuidanceTemperature(int whichWorld, int whichSampler);
    void  setGuidanceTemperature(int whichWorld, int whichSampler, float someTemperature);
    //------------

    // --- Simulation parameters ---
    int addSampler(int whichWorld, SamplerName whichSampler);
    void initializeSampler(int whichWorld, int whichSampler);

    // Amber like scale factors.
    void setAmberForceFieldScaleFactors(int whichWorld);

    // Set a global scaling factor for the forcefield
    void setGlobalForceFieldScaleFactor(int whichWorld, SimTK::Real);

    // Set GBSA implicit solvent scale factor
    void setGbsaGlobalScaleFactor(int whichWorld, SimTK::Real);

    // If HMC, get/set the number of MD steps
    int getNofMDStepsPerSample(int whichWorld, int whichSampler);
    void setNofMDStepsPerSample(int whichWorld, int whichSampler, int MDStepsPerSample);

    // If HMC, get/set timestep forMD
    const float getTimestep(int whichWorld, int whichSampler);
    void setTimestep(int whichWorld, int whichSampler, float timeStep);

    // Use Fixman torque as an additional force subsystem
    void useFixmanPotential(int whichWorld, int whichSampler);
    bool isUsingFixmanPotential(int whichWorld, int whichSampler);

    void addFixmanTorque(int whichWorld);
    bool isUsingFixmanTorque(int whichWorld);

    void setFixmanTorqueScaleFactor(int whichWorld, double scaleFactor);
    void setFixmanTorqueTemperature(int whichWorld, double temperature);
    //------------

    // --- Mixing parameters ---
    // Another way to do it is setting the number of rounds
    int getNofRounds(void);
    void setNofRounds(int nofRounds);

    int getNofSamplesPerRound(int whichWorld);
    void setNofSamplesPerRound(int whichWorld, int MCStepsPerRound);

    int getWorldIndex(int which);

    // --- Arrange different mixing parameters ---
    void initializeMixingParamters(void);
    //------------

    // --- Mix ---
    void RotateWorlds(void);
    //------------

    // --- Main ---
    void Run(SetupReader&);
    void Run(int howManyRounds, float Ti, float Tf);
    void setNumThreadsRequested(int which, int howMany);
    void setUseOpenMMAcceleration(bool arg);

    /** Print the number of threads each World got **/
    void PrintNumThreads(void);

    /** Get/Set seed for reproducibility. **/
    void setSeed(int whichWorld, int whichSampler, unsigned long long int);
    unsigned long long int getSeed(int whichWorld, int whichSampler);

    //------------

    /** Analysis related functions **/
    void addDistance(int whichWorld, int whichCompound, int aIx1, int aIx2);
    void addDihedral(int whichWorld, int whichCompound, int aIx1, int aIx2, int aIx3, int aIx4);

    // --- Printing functions ---
    void PrintSamplerData(unsigned int whichWorld);
    void PrintGeometry(SetupReader&, int whichWorld);
    void PrintGeometry(int whichWorld);
    void PrintDistances(int whichWorld);
    void PrintDihedrals(int whichWorld);
    void PrintDihedralsQs(int whichWorld);
    void PrintFreeE2EDist(int whichWorld, int whichCompound);
    void WritePdb(int whichWorld);
    SimTK::Real Dihedral(int whichWorld, int whichCompound, int whichSampler, int a1, int a2, int a3, int a4);
    SimTK::Real Distance(int whichWorld, int whichCompound, int whichSampler, int a1, int a2);

    int getPdbRestartFreq(void);
    void setPdbRestartFreq(int argFreq);
    int getPrintFreq(void);
    void setPrintFreq(int argFreq);

    std::string getOutputDir(void);
    void setOutputDir(std::string arg);
    std::string getPdbPrefix(void);
    void setPdbPrefix(std::string arg);
    //------------

public:
    std::vector<int> worldIndexes;

private:
    std::vector<World *> worlds;

    // Molecules files
    std::vector<std::vector<std::string>> topFNs;
    std::vector<std::vector<std::string>> crdFNs;
    std::vector<std::vector<std::string>> rbSpecsFNs;
    std::vector<std::vector<std::string>> flexSpecsFNs;
    std::vector<std::vector<std::string>> regimens;

    // Simulation parameters
    int nofRounds;
    //int total_mcsteps;
    std::vector<int> nofSamplesPerRound;
    std::vector<int> nofMDStepsPerSample;
    std::vector<float> timesteps;

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
    char *buffer;
    FILE *logFile;

};

#endif //__CONTEXT_HPP__





