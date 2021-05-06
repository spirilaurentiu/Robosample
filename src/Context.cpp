#include "Context.hpp"

// Includes to get the structure of additional classes

#include "World.hpp"
#include "Sampler.hpp"

// Default constructor
Context::Context(const SetupReader& setupReader, std::string logFilename){
    ValidateSetupReader(setupReader);

    BUFSIZE = 1024 * 1048576; // 1048576
    buffer = std::make_unique<char[]>(BUFSIZE);
    logFile = fopen(logFilename.c_str(), "w+");
    if ( setvbuf(logFile, &buffer[0], _IOFBF, BUFSIZE) != 0){
       perror("setvbuf()");
    }

    // Adaptive Gibbs blocki
    roundsTillReblock = std::stoi((setupReader.get("ROUNDS_TILL_REBLOCK"))[0]);

}

// // Constructor
// Context::Context(const SetupReader& setupReader, World *inp_p_world, std::string logFilename){
//     ValidateSetupReader(setupReader);
    
//     //total_mcsteps = 0;
//     worlds.emplace_back(std::unique_ptr<World>(inp_p_world)); // TODO make_unique converts ptr to int, might declare ctor explicit
//     worldIndexes.push_back(0);

//     topFNs.push_back(std::vector<std::string>());
//     crdFNs.push_back(std::vector<std::string>());
//     rbSpecsFNs.push_back(std::vector<std::string>());
//     flexSpecsFNs.push_back(std::vector<std::string>());
//     regimens.push_back(std::vector<std::string>());

//     nofSamplesPerRound.push_back(1);
//     nofMDStepsPerSample.push_back(1);
//     timesteps.push_back(0.002); // ps
//     nofBoostStairs.push_back(0);

//     pdbRestartFreq = 0;
//     outputDir = "pdbs";
//     pdbPrefix = "x";

//     BUFSIZE = 1024 * 1048576; // 1048576
//     buffer = std::make_unique<char[]>(BUFSIZE);
//     logFile = fopen(logFilename.c_str(), "w+");
//     if ( setvbuf(logFile, &buffer[0], _IOFBF, BUFSIZE) != 0){
//        perror("setvbuf()");
//     }

// }

// Add an empty world to the context
World * Context::AddWorld(bool visual, SimTK::Real visualizerFrequency){

    worldIndexes.push_back(worldIndexes.size());
    worlds.emplace_back(worldIndexes.back(), visual, visualizerFrequency);

    topFNs.push_back(std::vector<std::string>());
    crdFNs.push_back(std::vector<std::string>());
    rbSpecsFNs.push_back(std::vector<std::string>());
    flexSpecsFNs.push_back(std::vector<std::string>());
    regimens.push_back(std::vector<std::string>());

    nofSamplesPerRound.push_back(1);
    nofMDStepsPerSample.push_back(1);
    timesteps.push_back(0.002); // ps
    nofBoostStairs.push_back(0);

    // Adaptive Gibbs blocking
    QsCache.push_back(std::vector<std::vector<SimTK::Real>>(roundsTillReblock));
    //std::cout << "Context::AddWorld QsCache size " << QsCache.size() << std::endl;

    return &worlds.back();
}

// Destructor
Context::~Context(){
    fclose(logFile);
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
    return worlds.size();
}

SimTK::DuMMForceFieldSubsystem * Context::updForceField(std::size_t whichWorld)
{
    return worlds[whichWorld].updForceField();
}

// Input molecular files
bool Context::loadTopologyFile(std::size_t whichWorld, int, std::string topologyFilename)
{
    // function args were std::size_t whichWorld, int whichMolecule, std::string topologyFilename
    std::ifstream file(topologyFilename);
    if(!file){
        std::cout << topologyFilename << " not found." << std::endl;
        return false;
    }
    topFNs[whichWorld].push_back(topologyFilename);
    return true;
}

bool Context::loadCoordinatesFile(std::size_t whichWorld, int, std::string coordinatesFilename)
{
    // function args were std::size_t whichWorld, int whichMolecule, std::string coordinatesFilename
    std::ifstream file(coordinatesFilename);
    if(!file){
        std::cout << coordinatesFilename << " not found." << std::endl;
        return false;
    }
    crdFNs[whichWorld].push_back(coordinatesFilename);
    return true;
}

bool Context::loadRigidBodiesSpecs(std::size_t whichWorld, int, std::string RBSpecsFN)
{
    // function args were std::size_t whichWorld, int whichMolecule, std::string RBSpecsFN
    std::ifstream file(RBSpecsFN);
    if(!file){
        std::cout << RBSpecsFN << " not found." << std::endl;
        return false;
    }
    rbSpecsFNs[whichWorld].push_back(RBSpecsFN);
    return true;
}

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

/** Load molecules based on loaded filenames **/
void Context::AddMolecules(std::vector<std::string> argRoots)
{
    // TODO assert that the filenames vectors are not empty
    // Iterate through Worlds
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
    	std::cout << " Context::AddMolecule for world "<< worldIx << " " << std::endl;
        // Iterate through topology filenames vector
        for(unsigned int molIx = 0; molIx < topFNs[worldIx].size(); molIx++){
    	    std::cout << " Context::AddMolecule molIx "<< molIx << " " << std::endl;
    	    std::cout << " Context::AddMolecule topFNs[worldIx][molIx] "<< topFNs[worldIx][molIx] 
		<< " " << crdFNs[worldIx][molIx] << " " << rbSpecsFNs[worldIx][molIx] 
		//<< " " << flexSpecsFNs[worldIx][molIx] 
		//<< " " << regimens[worldIx][molIx] 
		<< std::endl << std::flush;
            // Initialize an input reader
            readAmberInput amberReader;
            amberReader.readAmberFiles(crdFNs[worldIx][molIx], topFNs[worldIx][molIx]);

            // Add the molecule to the World
            (updWorld(worldIx))->AddMolecule(&amberReader,
                    rbSpecsFNs[worldIx][molIx], flexSpecsFNs[worldIx][molIx],
                    regimens[worldIx][molIx], argRoots[worldIx]);
        }
    }
}

void Context::modelTopologies(std::vector<std::string> GroundToCompoundMobilizerTypes)
{
    // Model molecules
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
        this->rootMobilities.push_back(GroundToCompoundMobilizerTypes[worldIx]);
        (updWorld(worldIx))->modelTopologies(GroundToCompoundMobilizerTypes[worldIx]);
    }

}

int Context::getNofMolecules()
{
    return worlds[0].getNofMolecules();
}

// Set mixing rule for Lennard-Jones
void Context::setVdwMixingRule(SimTK::DuMMForceFieldSubsystem::VdwMixingRule mixingRule){
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
        (updWorld(worldIx))->updForceField()->setVdwMixingRule(mixingRule);
    }
}


// Use a SetupReader Object to read worlds information from a file
void Context::LoadWorldsFromSetup(SetupReader&)
{
    // function args were SetupReader& setupReader
    assert(!"Not implemented.");
    throw std::exception();
}

// --- Thermodynamics ---a
// Get/set the main temperature (acc/rej temperature for MC)
SimTK::Real Context::getTemperature(std::size_t whichWorld) const {
    return worlds[whichWorld].temperature;
}

void  Context::setTemperature(std::size_t whichWorld, float someTemperature){
    std::cout << " Context::setTemperature for world "<< whichWorld << " " << someTemperature << std::endl;
    worlds[whichWorld].setTemperature(someTemperature);
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

// --- Simulation parameters ---

BaseSampler * Context::addSampler(std::size_t whichWorld, SamplerName whichSampler)
{
    return worlds[whichWorld].addSampler(whichSampler);
}

void Context::initializeSampler(std::size_t whichWorld, std::size_t whichSampler)
{
    worlds[whichWorld].updSampler(whichSampler)->initialize( worlds[whichWorld].integ->updAdvancedState());
}

// Amber like scale factors.
void Context::setAmberForceFieldScaleFactors(std::size_t whichWorld)
{
    worlds[whichWorld].setAmberForceFieldScaleFactors();
}

// Set a global scaling factor for the forcefield
void Context::setGlobalForceFieldScaleFactor(std::size_t whichWorld, SimTK::Real globalScaleFactor)
{
    worlds[whichWorld].setGlobalForceFieldScaleFactor(globalScaleFactor);
}

// Set GBSA implicit solvent scale factor
void Context::setGbsaGlobalScaleFactor(std::size_t whichWorld, SimTK::Real gbsaGlobalScaleFactor)
{
    worlds[whichWorld].setGbsaGlobalScaleFactor(gbsaGlobalScaleFactor);
}

// If HMC, get/set the number of MD steps
int Context::getNofMDStepsPerSample(std::size_t whichWorld, std::size_t whichSampler){
   return pHMC(worlds[whichWorld].updSampler(whichSampler))->getMDStepsPerSample();
}

void Context::setNofMDStepsPerSample(std::size_t whichWorld, std::size_t whichSampler, int MDStepsPerSample)
{
   nofMDStepsPerSample[whichWorld] = MDStepsPerSample; // RE
   pHMC(worlds[whichWorld].updSampler(whichSampler))->setMDStepsPerSample(MDStepsPerSample); // NEW

   
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
    pMC(worlds[whichWorld].updSampler(whichSampler))->useFixmanPotential();
}

bool Context::isUsingFixmanPotential(std::size_t whichWorld, std::size_t whichSampler)
{
    return pMC(worlds[whichWorld].updSampler(whichSampler))->isUsingFixmanPotential();
}


//------------

// --- Mixing parameters ---

// Another way to do it is setting the number of rounds
int Context::getNofRounds()
{
    return nofRounds;
}

void Context::setNofRounds(int argNofRounds)
{
    nofRounds = argNofRounds;
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
void Context::Run(SetupReader&)
{
    // function args were SetupReader& setupReader
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
    // function args were int howManyRounds, SimTK::Real Ti, SimTK::Real Tf, bool isWorldOrderRandom

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
                currentAdvancedState = updWorld(currentWorldIx)->setAtomsLocationsInGround(
                        currentAdvancedState,
                        updWorld(lastWorldIx)->getAtomsLocationsInGround(lastAdvancedState));
            }

            // Check if reconstructions is done correctly
            //double backSetE = pMC(updWorld(lastWorldIx)->updSampler(0))->getSetPE();
            //double backCalcE = updWorld(lastWorldIx)->forceField->CalcFullPotEnergyIncludingRigidBodies(lastAdvancedState);
            //double currOldE = pMC(updWorld(currentWorldIx)->updSampler(0))->getOldPE();
            //double currCalcE = updWorld(currentWorldIx)->forceField->CalcFullPotEnergyIncludingRigidBodies(currentAdvancedState);

            // Set old potential energy of the new world
            pMC(updWorld(currentWorldIx)->updSampler(0))->setOldPE(
                    pMC(updWorld(lastWorldIx)->updSampler(0))->getSetPE() );

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
        SimTK::State& pdbState = (updWorld(worldIndexes.front()))->integ->updAdvancedState();
        if( getPdbRestartFreq() != 0){
            if(((round) % getPdbRestartFreq()) == 0){
                (updWorld(worldIndexes.front()))->updateAtomListsFromCompound(pdbState);
                for(int mol_i = 0; mol_i < getNofMolecules(); mol_i++){
                    ((updWorld(worldIndexes.front()))->getTopology(mol_i)).writeAtomListPdb(getOutputDir(),
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


// Main
void Context::Run(int, SimTK::Real Ti, SimTK::Real Tf, bool isWorldOrderRandom)
{
    // function args were int howManyRounds, SimTK::Real Ti, SimTK::Real Tf, bool isWorldOrderRandom

    // Initialize world indeces
    std::size_t currentWorldIx = worldIndexes.front();
    std::size_t lastWorldIx = 0;

    // Random int for random world order
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<std::size_t> randWorldDistrib(1, getNofWorlds()-1); // TODO between 1 and nOfWorlds-1?

    // Main loop: iterate through rounds
    if( std::abs(Tf - Ti) < SimTK::TinyReal){ // Don't heat
        for(int round = 0; round < nofRounds; round++){

            // std::cout << "Entering round " << round + 1 << "/" << nofRounds << "\n";

            // Iterate through worlds
            for(std::size_t worldIx = 0; worldIx < getNofWorlds(); worldIx++){

                // Rotate worlds indices (translate from right to left)
                if(isWorldOrderRandom){
                    if(getNofWorlds() >= 3){
                        // Swap world indeces between vector position 2 and random
                        auto randVecPos = randWorldDistrib(gen);
                        //std::cout << "Swapping position 1 with " << randVecPos << std::endl;
                        auto secondWorldIx = worldIndexes[1];
                        auto randWorldIx = worldIndexes[randVecPos];

                        worldIndexes[1] = randWorldIx;
                        worldIndexes[randVecPos] = secondWorldIx;

                    }
                }

               	std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());

                // Get indeces
                currentWorldIx = worldIndexes.front();
                lastWorldIx = worldIndexes.back();

                // Transfer coordinates from last world to current
                SimTK::State& lastAdvancedState = updWorld(lastWorldIx)->integ->updAdvancedState();
                SimTK::State& currentAdvancedState = updWorld(currentWorldIx)->integ->updAdvancedState();

                if(worldIndexes.size() > 1) { // It also loads bAtomList
                	currentAdvancedState = updWorld(currentWorldIx)->setAtomsLocationsInGround(
                        	currentAdvancedState,
                        	updWorld(lastWorldIx)->getAtomsLocationsInGround(lastAdvancedState));
                }else{ // Just load bAtomList
                    updWorld(currentWorldIx)->updateAtomListsFromCompound(currentAdvancedState);
                }

                // Set old potential energy of the new world via OpenMM
                auto OldPE = updWorld(currentWorldIx)->updSampler(0)->forces->getMultibodySystem().calcPotentialEnergy(currentAdvancedState);
		        pMC(updWorld(currentWorldIx)->updSampler(0))->setOldPE(OldPE);
			    //updWorld(currentWorldIx)->forceField->CalcFullPotEnergyIncludingRigidBodies(currentAdvancedState));

                std::cout << "World " << currentWorldIx << ", NU " << currentAdvancedState.getNU() << ":\n";

                // Reinitialize current sampler
                updWorld(currentWorldIx)->updSampler(0)->reinitialize(currentAdvancedState);

                // Make the requested number of samples
                for(int k = 0; k < getNofSamplesPerRound(currentWorldIx); k++) {
                    auto accepted = updWorld(currentWorldIx)->updSampler(0)->sample_iteration(currentAdvancedState);
                    if (accepted) {

                        // CONTACT DEBUG
                        int numForces = updWorld(currentWorldIx)->contactForces->getNumContactForces(currentAdvancedState);
                        SimTK::Real dissEnergy = updWorld(currentWorldIx)->contactForces->getDissipatedEnergy(currentAdvancedState);
                        bool hasDefaultForceGenerator = updWorld(currentWorldIx)->contactForces->hasDefaultForceGenerator();

                        const MultibodySystem & mbs = updWorld(currentWorldIx)->contactForces->getMultibodySystem();
                        int nofMobods = mbs.getMatterSubsystem().getNumBodies();

                        const ContactTrackerSubsystem & cts = updWorld(currentWorldIx)->contactForces->getContactTrackerSubsystem();
                        int ctsNofSurfaces = cts.getNumSurfaces();
                        

                        std::cout << "CONTACT INFO:"
                            << " #forces= " << numForces
                            << " dissEnergy= " << dissEnergy
                            << " hasDefaultForceGenerator= " << hasDefaultForceGenerator
                            << " #mobods= " << nofMobods 
                            << " ctsNofSurfaces= " << ctsNofSurfaces
                        << std::endl;
                        // CONTACT DEBUG enD
                    }
                }
    
            } // END iteration through worlds

            // Print energy and geometric features
            if( !(round % getPrintFreq()) ){
                PrintSamplerData(worldIndexes.back());
                PrintDistances(worldIndexes.back());
                PrintDihedralsQs(worldIndexes.back());
                fprintf(logFile, "\n");
                PrintSamplerData(worldIndexes.front());
                PrintDistances(worldIndexes.front());
                PrintDihedralsQs(worldIndexes.front());
                fprintf(logFile, "\n");
            }

            // std::cout << "\n";
    
            // Write pdb
            SimTK::State& pdbState = (updWorld(worldIndexes.front()))->integ->updAdvancedState();
            if( getPdbRestartFreq() != 0){
                if(((round) % getPdbRestartFreq()) == 0){
                    (updWorld(worldIndexes.front()))->updateAtomListsFromCompound(pdbState);
                    for(int mol_i = 0; mol_i < getNofMolecules(); mol_i++){
                        ((updWorld(worldIndexes.front()))->getTopology(mol_i)).writeAtomListPdb(getOutputDir(),
                                                                                                "/pdbs/sb." +
                                                                                                getPdbPrefix() + "." + std::to_string(mol_i) + ".",
                                                                                                ".pdb", 10, round);
                    }
                }
            }
    
        } // END main loop

    }else{// if Ti != Tf heating protocol
        SimTK::Real Tincr = (Tf - Ti) / static_cast<SimTK::Real>(nofRounds);
        SimTK::Real currT = Ti;
        for(int round = 0; round < nofRounds; round++){ // Iterate rounds
    
            // Set current temperature
            currT += Tincr;
            std::cout << "T= " << currT << std::endl;

            for(std::size_t worldIx = 0; worldIx < getNofWorlds(); worldIx++){ // Iterate worlds
    
                // Rotate worlds indeces (translate from right to left)
                std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());
    
                // Get indeces
                currentWorldIx = worldIndexes.front();
                lastWorldIx = worldIndexes.back();

                // Set temperatures
                updWorld(lastWorldIx)->setTemperature(currT);
                updWorld(currentWorldIx)->setTemperature(currT);

                if(isUsingFixmanTorque(worldIx)){
                    setFixmanTorqueTemperature(worldIx, currT);
                }
    
                // Transfer coordinates from last world to current
                SimTK::State& lastAdvancedState = updWorld(worldIndexes.back())->integ->updAdvancedState();
                SimTK::State& currentAdvancedState = updWorld(currentWorldIx)->integ->updAdvancedState();

                if(worldIndexes.size() > 1) {
                    currentAdvancedState = updWorld(currentWorldIx)->setAtomsLocationsInGround(
                            currentAdvancedState,
                            updWorld(worldIndexes.back())->getAtomsLocationsInGround(lastAdvancedState));
                }else{ // Load bAtomList though
                    updWorld(currentWorldIx)->updateAtomListsFromCompound(currentAdvancedState);
                }
    
                // Set old potential energy of the new world
                const auto oldPE = pMC(updWorld(worldIndexes.back())->updSampler(0))->getSetPE();
                pMC(updWorld(currentWorldIx)->updSampler(0))->setOldPE(oldPE);
    
                // Reinitialize current sampler
                updWorld(currentWorldIx)->updSampler(0)->reinitialize(currentAdvancedState);
    
                // Update
                for(int k = 0; k < getNofSamplesPerRound(currentWorldIx); k++){ // Iterate through samples
                    //updWorld(currentWorldIx)->updSampler(0)->propose(currentAdvancedState);
                    updWorld(currentWorldIx)->updSampler(0)->sample_iteration(currentAdvancedState);

                } // END for samples
    
            } // for i in worlds
    
            // Print energy and geometric features
            if( !(round % getPrintFreq()) ){
                PrintSamplerData(worldIndexes.back());
                PrintDistances(worldIndexes.front());
                PrintDihedralsQs(worldIndexes.back());
                fprintf(logFile, "\n");
            }
    
            // Write pdb
            SimTK::State& pdbState = (updWorld(worldIndexes.front()))->integ->updAdvancedState();
            if( getPdbRestartFreq() != 0){
                if(((round) % getPdbRestartFreq()) == 0){
                    (updWorld(worldIndexes.front()))->updateAtomListsFromCompound(pdbState);
                    for(int mol_i = 0; mol_i < getNofMolecules(); mol_i++){
                        ((updWorld(worldIndexes.front()))->getTopology(mol_i)).writeAtomListPdb(getOutputDir(),
                                                                                                "/pdbs/sb." +
                                                                                                getPdbPrefix() + "." + std::to_string(mol_i) + ".",
                                                                                                ".pdb", 10, round);
                    }
                }
            } // if write pdbs
    
        } // for i in rounds

    } // if heating protocol


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

void Context::setUseOpenMMCalcOnlyNonBonded(bool arg)
{
	
     for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
	worlds[worldIx].updForceField()->setUseOpenMMCalcOnlyNonBonded(arg);
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
    // TODO some are 64 bit, some 32 bit. What do?
    std::vector<int> tempV;
    tempV.push_back(static_cast<int>(whichWorld));
    tempV.push_back(static_cast<int>(whichCompound));
    tempV.push_back(aIx1);
    tempV.push_back(aIx2);

    distanceIxs.push_back(tempV);
}

void Context::addDihedral(std::size_t whichWorld, std::size_t whichCompound, int aIx1, int aIx2, int aIx3, int aIx4)
{
    // TODO some are 64 bit, some 32 bit. What do?
    std::vector<int> tempV;
    tempV.push_back(static_cast<int>(whichWorld));
    tempV.push_back(static_cast<int>(whichCompound));
    tempV.push_back(aIx1);
    tempV.push_back(aIx2);
    tempV.push_back(aIx3);
    tempV.push_back(aIx4);

    dihedralIxs.push_back(tempV);
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
    fprintf(logFile, "%d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f "
        , currentAdvancedState.getNU()
        , pHMC(worlds[whichWorld].samplers[0])->acceptedSteps
        , pMC((worlds[whichWorld].samplers[0]))->pe_o
        , pMC((worlds[whichWorld].samplers[0]))->pe_set
        , pHMC((worlds[whichWorld].samplers[0]))->ke_proposed
        , pHMC((worlds[whichWorld].samplers[0]))->ke_n
        , pMC((worlds[whichWorld].samplers[0]))->fix_o
        , pMC((worlds[whichWorld].samplers[0]))->fix_n
        , pMC((worlds[whichWorld].samplers[0]))->fix_set
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

            fprintf(logFile, "%.3f ", Dihedral(dihedralIx[0], dihedralIx[1], 0, dihedralIx[2], dihedralIx[3], dihedralIx[4], dihedralIx[5]) );

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


SimTK::Real Context::Dihedral(std::size_t whichWorld, std::size_t whichCompound, std::size_t, int a1, int a2, int a3, int a4)
{
    // function args were std::size_t whichWorld, std::size_t whichCompound, std::size_t whichSampler, int a1, int a2, int a3, int a4

    //SimTK::State& currentAdvancedState = (updWorld(whichWorld))->integ->updAdvancedState();

    SimTK::State& currentAdvancedState = worlds[whichWorld].integ->updAdvancedState();

    //SimTK::State& currentAdvancedState = (worlds[whichWorld].updSampler(whichSampler)->updTimeStepper()->updIntegrator()).updAdvancedState();

    const Topology& topology = worlds[whichWorld].getTopology(whichCompound);
    SimTK::Vec3 a1pos, a2pos, a3pos, a4pos;
    a1pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
    a2pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));
    a3pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)));
    a4pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a4)));

    //std::cout << " poss: " << a1pos << ' ' << a2pos << ' ' << a3pos << ' ' << a4pos << ' ';
    //std::cout << " dih: "  << bDihedral(a1pos, a2pos, a3pos, a4pos) << '|' ;

    return bDihedral(a1pos, a2pos, a3pos, a4pos);

}

SimTK::Real Context::Distance(std::size_t whichWorld, std::size_t whichCompound, std::size_t, int a1, int a2)
{
    // function args were std::size_t whichWorld, std::size_t whichCompound, std::size_t whichSampler, int a1, int a2

    SimTK::State& currentAdvancedState = worlds[whichWorld].integ->updAdvancedState();

    const Topology& topology = worlds[whichWorld].getTopology(whichCompound);
    SimTK::Vec3 a1pos, a2pos;
    a1pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
    a2pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));

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

void Context::ValidateSetupReader(const SetupReader& setupReader) {
    const std::size_t nofWorlds = setupReader.get("WORLDS").size();

    for(std::size_t worldIx = 0; worldIx < nofWorlds; worldIx++){
        for(std::size_t molIx = 0; molIx < setupReader.get("MOLECULES").size(); molIx++){
            assert(SimTK::Pathname::fileExists(
                    setupReader.get("MOLECULES")[molIx] + std::string("/")
                + setupReader.get("PRMTOP")[molIx]) );
            assert(SimTK::Pathname::fileExists(
                    setupReader.get("MOLECULES")[molIx] + std::string("/")
                + setupReader.get("INPCRD")[molIx]) );
            assert(SimTK::Pathname::fileExists(
                    setupReader.get("MOLECULES")[molIx] + std::string("/")
                + setupReader.get("RBFILE")[molIx]) );
            assert(SimTK::Pathname::fileExists(
                    setupReader.get("MOLECULES")[molIx] + std::string("/")
                + setupReader.get("FLEXFILE")[molIx]) );

        }
    }

    assert(SimTK::Pathname::fileExists(setupReader.get("OUTPUT_DIR")[0]));
}






//------------



