#include "Context.hpp"

// Includes to get the structure of additional classes

#include "World.hpp"
#include "Sampler.hpp"

// Default constructor
Context::Context(std::string logFilename){
    BUFSIZE = 1048576;
    buffer = new char[BUFSIZE];
    logFile = fopen(logFilename.c_str(), "w+");
    if ( setvbuf(logFile, buffer, _IOFBF, BUFSIZE) != 0){
       perror("setvbuf()");
    }
}

// Constructor
Context::Context(World *inp_p_world, std::string logFilename){
    //total_mcsteps = 0;
    worlds.push_back(inp_p_world);
    worldIndexes.push_back(0);

    topFNs.push_back(std::vector<std::string>());
    crdFNs.push_back(std::vector<std::string>());
    rbSpecsFNs.push_back(std::vector<std::string>());
    flexSpecsFNs.push_back(std::vector<std::string>());
    regimens.push_back(std::vector<std::string>());

    nofSamplesPerRound.push_back(1);
    nofMDStepsPerSample.push_back(1);
    timesteps.push_back(0.002); // ps

    pdbRestartFreq = 0;
    outputDir = "pdbs";
    pdbPrefix = "x";

    BUFSIZE = 1048576;
    buffer = new char[BUFSIZE];
    logFile = fopen(logFilename.c_str(), "w+");
    if ( setvbuf(logFile, buffer, _IOFBF, BUFSIZE) != 0){
       perror("setvbuf()");
    }

}

// Add an empty world to the context
World * Context::AddWorld(bool visual, SimTK::Real visualizerFrequency){

    worldIndexes.push_back(worldIndexes.size());
    World * inp_p_world = new World(worldIndexes.back(), visual, visualizerFrequency);
    worlds.push_back(inp_p_world);

    topFNs.push_back(std::vector<std::string>());
    crdFNs.push_back(std::vector<std::string>());
    rbSpecsFNs.push_back(std::vector<std::string>());
    flexSpecsFNs.push_back(std::vector<std::string>());
    regimens.push_back(std::vector<std::string>());

    nofSamplesPerRound.push_back(1);
    nofMDStepsPerSample.push_back(1);
    timesteps.push_back(0.002); // ps

    return worlds.back();
}

// Destructor
Context::~Context(){
    // Delete each world
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
         delete worlds[worldIx];
    }
    worlds.clear();
    worldIndexes.clear();
   
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
        topFNs[worldIx].clear();
        crdFNs[worldIx].clear();
        rbSpecsFNs[worldIx].clear();
        flexSpecsFNs[worldIx].clear();
        regimens[worldIx].clear();
    }
    topFNs.clear();
    crdFNs.clear();
    rbSpecsFNs.clear();
    flexSpecsFNs.clear();
    regimens.clear();
 
    nofSamplesPerRound.clear();
    nofMDStepsPerSample.clear();
    timesteps.clear();

    fclose(logFile);
    delete buffer;
 
}

// Get world
World * Context::getWorld(void) const{
    return worlds.back();
}

// Get a specific world
World * Context::getWorld(int which) const{
    return worlds[which];
}

// Get the last mutable world
World * Context::updWorld(void){
    return worlds.back();
}

// Get a mutable specific world
World * Context::updWorld(int which){
    return worlds[which];
}

unsigned int Context::getNofWorlds(void)
{
    return worlds.size();
}



SimTK::DuMMForceFieldSubsystem * Context::updForceField(int whichWorld)
{
    return worlds[whichWorld]->updForceField();
}


// Input molecular files
bool Context::loadTopologyFile(int whichWorld, int whichMolecule, std::string topologyFilename)
{
    std::ifstream file(topologyFilename);
    if(!file){
        std::cout << topologyFilename << " not found." << std::endl;
        return false;
    }
    topFNs[whichWorld].push_back(topologyFilename);
    return true;
}

bool Context::loadCoordinatesFile(int whichWorld, int whichMolecule, std::string coordinatesFilename)
{
    std::ifstream file(coordinatesFilename);
    if(!file){
        std::cout << coordinatesFilename << " not found." << std::endl;
        return false;
    }
    crdFNs[whichWorld].push_back(coordinatesFilename);
    return true;
}

bool Context::loadRigidBodiesSpecs(int whichWorld, int whichMolecule, std::string RBSpecsFN)
{
    std::ifstream file(RBSpecsFN);
    if(!file){
        std::cout << RBSpecsFN << " not found." << std::endl;
        return false;
    }
    rbSpecsFNs[whichWorld].push_back(RBSpecsFN);
    return true;
}

bool Context::loadFlexibleBondsSpecs(int whichWorld, int whichMolecule, std::string flexSpecsFN)
{
    std::ifstream file(flexSpecsFN);
    if(!file){
        std::cout << flexSpecsFN << " not found." << std::endl;
        return false;
    }
    flexSpecsFNs[whichWorld].push_back(flexSpecsFN);
    return true;
}

void Context::setRegimen (int whichWorld, int whichMolecule, std::string regimen)
{
    regimens[whichWorld].push_back(regimen);
}

/** Load molecules based on loaded filenames **/
void Context::AddMolecules(std::vector<std::string> argRoots)
{
    // TODO assert that the filenames vectors are not empty
    // Iterate through Worlds
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
        // Iterate through topology filenames vector
        for(unsigned int molIx = 0; molIx < topFNs[worldIx].size(); molIx++){
            // Initialize an input reader
            readAmberInput *amberReader = new readAmberInput();
            amberReader->readAmberFiles(crdFNs[worldIx][molIx], topFNs[worldIx][molIx]);

            // Add the molecule to the World
            (updWorld(worldIx))->AddMolecule(amberReader,
                    rbSpecsFNs[worldIx][molIx], flexSpecsFNs[worldIx][molIx],
                    regimens[worldIx][molIx], argRoots[worldIx]);

            // Delete the reader
            delete amberReader; 
        }
    }
}

void Context::modelTopologies(void)
{
    // Model molecules
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
        (updWorld(worldIx))->ModelTopologies();
    }

}

int Context::getNofMolecules(void)
{
    return worlds[0]->getNofMolecules();
}


// Use a SetupReader Object to read worlds information from a file
void Context::LoadWorldsFromSetup(SetupReader& setupReader)
{
    assert(!"Not implemented.");
}

// --- Thermodynamics ---a
// Get/set the main temperature (acc/rej temperature for MC)
float Context::getTemperature(int whichWorld){
    return worlds[whichWorld]->temperature;
}

void  Context::setTemperature(int whichWorld, float someTemperature){
    std::cout << " Context::setTemperature for world "<< whichWorld << " " << someTemperature << std::endl;
    worlds[whichWorld]->setTemperature(someTemperature);
}

    // If HMC, get/set the guidance Hamiltonian temperature
float Context::getGuidanceTemperature(int whichWorld, int whichSampler)
{
   assert(!"Not implemented"); 
}

void Context::setGuidanceTemperature(int whichWorld, int whichSampler, float someTemperature)
{
   assert(!"Not implemented"); 
}
//------------

// --- Simulation parameters ---

int Context::addSampler(int whichWorld, SamplerName whichSampler)
{
    worlds[whichWorld]->addSampler(whichSampler);
}

void Context::initializeSampler(int whichWorld, int whichSampler)
{
    worlds[whichWorld]->updSampler(whichSampler)->initialize( worlds[whichWorld]->integ->updAdvancedState());
}

// Amber like scale factors.
void Context::setAmberForceFieldScaleFactors(int whichWorld)
{
    worlds[whichWorld]->setAmberForceFieldScaleFactors();
}

// Set a global scaling factor for the forcefield
void Context::setGlobalForceFieldScaleFactor(int whichWorld, SimTK::Real globalScaleFactor)
{
    worlds[whichWorld]->setGlobalForceFieldScaleFactor(globalScaleFactor);
}

// Set GBSA implicit solvent scale factor
void Context::setGbsaGlobalScaleFactor(int whichWorld, SimTK::Real gbsaGlobalScaleFactor)
{
    worlds[whichWorld]->setGbsaGlobalScaleFactor(gbsaGlobalScaleFactor);
}

// If HMC, get/set the number of MD steps
int Context::getNofMDStepsPerSample(int whichWorld, int whichSampler){
   return pHMC(worlds[whichWorld]->updSampler(whichSampler))->getMDStepsPerSample();
}

void Context::setNofMDStepsPerSample(int whichWorld, int whichSampler, int MDStepsPerSample)
{
   nofMDStepsPerSample[whichWorld] = MDStepsPerSample; // RE
   pHMC(worlds[whichWorld]->updSampler(whichSampler))->setMDStepsPerSample(MDStepsPerSample); // NEW
}

// If HMC, get/set timestep forMD
const float Context::getTimestep(int whichWorld, int whichSampler)
{
    return pHMC(worlds[whichWorld]->updSampler(whichSampler))->getTimeStepper()->getIntegrator().getPredictedNextStepSize();

}

void Context::setTimestep(int whichWorld, int whichSampler, float argTimestep)
{
    //worlds[whichWorld]->updSampler(whichSampler)->updTimeStepper()->updIntegrator().setFixedStepSize(argTimestep);
    pHMC(worlds[whichWorld]->updSampler(whichSampler))->setTimestep(argTimestep);
}

// Use Fixman torque as an additional force subsystem
void Context::addFixmanTorque(int whichWorld)
{
    worlds[whichWorld]->addFixmanTorque();
}

bool Context::isUsingFixmanTorque(int whichWorld)
{
    return worlds[whichWorld]->isUsingFixmanTorque();
}

void Context::setFixmanTorqueScaleFactor(int whichWorld, double scaleFactor)
{
    std::cout << "Context::setFixmanTorqueScaleFactor: ( (FixmanTorque *) (worlds[" 
    << whichWorld << "]->updFixmanTorque()) )->setScaleFactor(" << scaleFactor << ") "<< std::endl;
    ( (FixmanTorque *) (worlds[whichWorld]->updFixmanTorque()) )->setScaleFactor(scaleFactor);
}

void Context::setFixmanTorqueTemperature(int whichWorld, double argTemperature)
{
    std::cout << "Context::setFixmanTemperature: ( (FixmanTorque *) (worlds[" 
    << whichWorld << "]->updFixmanTorque()) )->setTemperature(" << argTemperature << ") "<< std::endl;
    ( (FixmanTorque *) (worlds[whichWorld]->updFixmanTorque()) )->setTemperature(argTemperature);
}

// Use Fixman potential
void Context::useFixmanPotential(int whichWorld, int whichSampler)
{
    pMC(worlds[whichWorld]->updSampler(whichSampler))->useFixmanPotential();
}

bool Context::isUsingFixmanPotential(int whichWorld, int whichSampler)
{
    return pMC(worlds[whichWorld]->updSampler(whichSampler))->isUsingFixmanPotential();
}


//------------

// --- Mixing parameters ---

// Another way to do it is setting the number of rounds
int Context::getNofRounds(void)
{
    return nofRounds;
}

void Context::setNofRounds(int argNofRounds)
{
    nofRounds = argNofRounds;
}

// Get the number of samples returned by the sampler in one round
int Context::getNofSamplesPerRound(int whichWorld)
{
    return nofSamplesPerRound[whichWorld];
}

// Set the number of samples returned by the sampler in one round
void Context::setNofSamplesPerRound(int whichWorld, int MCStepsPerRound)
{
    nofSamplesPerRound[whichWorld] = MCStepsPerRound;
}

// Return the world index in position 'which'. To be used when rotationg
int Context::getWorldIndex(int which)
{
    return worldIndexes[which];
}

// --- Arrange different mixing parameters ---
void Context::initializeMixingParamters(void){assert(!"Not implemented");}
//------------

// --- Mix ---
void Context::RotateWorlds(void){assert(!"Not implemented");}
//------------

// -- Main ---
void Context::Run(SetupReader& setupReader)
{
    assert(!"Not implemented");
} 

// Main
void Context::Run(int howManyRounds, float Ti, float Tf)
{

    int currentWorldIx = worldIndexes.front();
    int lastWorldIx = 0;

    // TODO move or delete
    // Write the initial Default Configuration of the first Compound
    // of the first World
    PdbStructure  pdb(worlds[0]->getTopology(0));
    std::ostringstream sstream;
    sstream << "pdbs/sb_" << ((updWorld(worldIndexes.back()))->getTopology(0)).getName() <<"_ini"<<".pdb";
    std::string ofilename = sstream.str();
    std::filebuf fb;
    std::cout<<"Writing pdb file: "<<ofilename<<std::endl;
    fb.open(ofilename.c_str(), std::ios::out);
    std::ostream os(&fb);
    pdb.write(os); // automatically multiplies by ten (nm to A)
    fb.close();
    //

    if( std::abs(Tf - Ti) < SimTK::TinyReal){
        for(int round = 0; round < nofRounds; round++){ // Iterate rounds

            //std::cout << "round " << round << std::endl;

// TIME START -----------------
//std::chrono::steady_clock::time_point start0 = std::chrono::steady_clock::now();
// TIME START

            for(unsigned int worldIx = 0; worldIx < getNofWorlds(); worldIx++){ // Iterate worlds

                // Rotate worlds indeces (translate from right to left)
                std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());

                //std::cout << "world " << worldIndexes.front() << std::endl;

                // Get indeces
                currentWorldIx = worldIndexes.front();
                lastWorldIx = worldIndexes.back();

                // Transfer coordinates from last world to current
                SimTK::State& lastAdvancedState = (updWorld(worldIndexes.back()))->integ->updAdvancedState();
                SimTK::State& currentAdvancedState = (updWorld(currentWorldIx))->integ->updAdvancedState();

// TIME STOP
//std::chrono::steady_clock::time_point end0 = std::chrono::steady_clock::now();
//std::cout << "Context::Run end0 - start0 "
//              << std::chrono::duration_cast<std::chrono::microseconds>(end0 - start0).count()
//              << " us.\n";
// TIME STOP ===================
                if(worldIndexes.size() > 1) {
                    currentAdvancedState = (updWorld(currentWorldIx))->setAtomsLocationsInGround(
                            currentAdvancedState,
                            (updWorld(worldIndexes.back()))->getAtomsLocationsInGround(lastAdvancedState));
                }

                double backSetE = pMC(updWorld(worldIndexes.back())->updSampler(0))->getSetPE();
                double backCalcE = updWorld(worldIndexes.back())->forceField->CalcFullPotEnergyIncludingRigidBodies(
                        lastAdvancedState);
                double currOldE = pMC(updWorld(currentWorldIx)->updSampler(0))->getOldPE();
                double currCalcE = updWorld(currentWorldIx)->forceField->CalcFullPotEnergyIncludingRigidBodies(
                        currentAdvancedState);

                // Check if reconstructions is done correctly
                //assert(std::abs(backCalcE - currCalcE) < 0.00001);
                //assert((backCalcE - currCalcE) > SimTK::SignificantReal);
                TRACE("Energy after reconstruction "); TRACE(std::to_string(std::abs(backCalcE - currCalcE)).c_str());

  // TIME STOP
//std::chrono::steady_clock::time_point end1_0 = std::chrono::steady_clock::now();
//std::cout << "Context::Run end1_0 - start0 "
//              << std::chrono::duration_cast<std::chrono::microseconds>(end1_0 - start0).count()
//              << " us.\n";
// TIME STOP ===================

                // Set old potential energy of the new world
                pMC((updWorld(currentWorldIx))->updSampler(0))->setOldPE(
                    pMC((updWorld(worldIndexes.back()))
                    ->updSampler(0))->getSetPE() );


/*                std::cout << "RunPe backSet backCalc currCalc currOld "
                          << backSetE << " " << backCalcE << " "
                          << currCalcE << " " << currOldE
                          << std::endl;*/
                if(std::abs(backCalcE - currCalcE) > 100.0) {
                    std::cout << "Writing Compound pdbs for round " << round << std::endl;
                    (updWorld(worldIndexes.back()))->updateAtomLists(lastAdvancedState);
                    ((updWorld(worldIndexes.back()))->getTopology(0)).writeAtomListPdb(
                            getOutputDir(), std::string("/pdbs/sb.")
                                            + std::to_string(worldIx) + std::string("."), ".0.pdb", 10, round);

                    (updWorld(currentWorldIx))->updateAtomLists(currentAdvancedState);
                    ((updWorld(currentWorldIx))->getTopology(0)).writeAtomListPdb(
                            getOutputDir(), std::string("/pdbs/sb.")
                                            + std::to_string(worldIx) + std::string("."), ".1.pdb", 10, round);

                    exit(1);
                }

// TIME STOP
//std::chrono::steady_clock::time_point end1_1 = std::chrono::steady_clock::now();
//std::cout << "Context::Run end1_1 - start0 "
//              << std::chrono::duration_cast<std::chrono::microseconds>(end1_1 - start0).count()
//              << " us.\n";
// TIME STOP ===================

    
                // Reinitialize current sampler
                updWorld(currentWorldIx)->updSampler(0)->reinitialize(currentAdvancedState);
// TIME STOP
//std::chrono::steady_clock::time_point end1_2 = std::chrono::steady_clock::now();
//std::cout << "Context::Run end1_2 - start0 "
//              << std::chrono::duration_cast<std::chrono::microseconds>(end1_2 - start0).count()
//              << " us.\n";
// TIME STOP ===================

    
                // Update
                for(int k = 0; k < getNofSamplesPerRound(currentWorldIx); k++){ // Iterate through samples
                    updWorld(currentWorldIx)->updSampler(0)->propose(currentAdvancedState);
                    updWorld(currentWorldIx)->updSampler(0)->update(currentAdvancedState);
                            // , getNofMDStepsPerSample(currentWorldIx, 0)); RE

// TIME STOP
//std::chrono::steady_clock::time_point end1_3 = std::chrono::steady_clock::now();
//std::cout << "Context::Run end1_3 - start0 "
//              << std::chrono::duration_cast<std::chrono::microseconds>(end1_3 - start0).count()
//              << " us.\n";
// TIME STOP ===================

                } // END for samples
    
            } // for i in worlds

            // Print energy and geometric features
            if( !(round % getPrintFreq()) ){
                PrintSamplerData(worldIndexes.back());
                PrintDistances(worldIndexes.front());
                PrintDihedralsQs(worldIndexes.back());
                fprintf(logFile, "\n");
            }
    
// TIME STOP
//std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
//std::cout << "Context::Run end3 - start0 "
//              << std::chrono::duration_cast<std::chrono::microseconds>(end3 - start0).count()
//              << " us.\n";
// TIME STOP ===================

            // Write pdb
            SimTK::State& pdbState = (updWorld(worldIndexes.front()))->integ->updAdvancedState();
            if( getPdbRestartFreq() != 0){
                if(((round) % getPdbRestartFreq()) == 0){
                    (updWorld(worldIndexes.front()))->updateAtomLists(pdbState);
                    for(int mol_i = 0; mol_i < getNofMolecules(); mol_i++){
                        ((updWorld(worldIndexes.front()))->getTopology(mol_i)).writeAtomListPdb(getOutputDir(),
                                                                                                "/pdbs/sb." +
                                                                                                getPdbPrefix() + ".",
                                                                                                ".pdb", 10, round);
                    }
                }
            } // if write pdbs
    
// TIME STOP
//std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();
//std::cout << "Context::Run end4 - start0 "
//              << std::chrono::duration_cast<std::chrono::microseconds>(end4 - start0).count()
//              << " us.\n";
// TIME STOP ===================

        } // for i in rounds

    }else{// if Ti != Tf heating protocol
        SimTK::Real Tincr = (Tf - Ti) / nofRounds;
        SimTK::Real currT = Ti;
        for(int round = 0; round < nofRounds; round++){ // Iterate rounds
    
            // Set current temperature
            currT += Tincr;
            std::cout << "T= " << currT << std::endl;


            for(unsigned int worldIx = 0; worldIx < getNofWorlds(); worldIx++){ // Iterate worlds
    
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
                SimTK::State& lastAdvancedState = (updWorld(worldIndexes.back()))->integ->updAdvancedState();
                SimTK::State& currentAdvancedState = (updWorld(currentWorldIx))->integ->updAdvancedState();

                if(worldIndexes.size() > 1) {
                    currentAdvancedState = (updWorld(currentWorldIx))->setAtomsLocationsInGround(
                            currentAdvancedState,
                            (updWorld(worldIndexes.back()))->getAtomsLocationsInGround(lastAdvancedState));
                }
    
                // Set old potential energy of the new world
                pMC((updWorld(currentWorldIx))->updSampler(0))->setOldPE(
                        pMC((updWorld(worldIndexes.back()))
                        ->updSampler(0))->getSetPE() );
    
                // Reinitialize current sampler
                updWorld(currentWorldIx)->updSampler(0)->reinitialize(currentAdvancedState);
    
                // Update
                for(int k = 0; k < getNofSamplesPerRound(currentWorldIx); k++){ // Iterate through samples
                    updWorld(currentWorldIx)->updSampler(0)->propose(currentAdvancedState);
                    updWorld(currentWorldIx)->updSampler(0)->update(currentAdvancedState);
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
                    (updWorld(worldIndexes.front()))->updateAtomLists(pdbState);
                    for(int mol_i = 0; mol_i < getNofMolecules(); mol_i++){
                        ((updWorld(worldIndexes.front()))->getTopology(mol_i)).writeAtomListPdb(getOutputDir(),
                                                                                                "/pdbs/sb." +
                                                                                                getPdbPrefix() + ".",
                                                                                                ".pdb", 10, round);
                    }
                }
            } // if write pdbs
    
        } // for i in rounds

    } // if heating protocol


}

// Set number of threads
void Context::setNumThreadsRequested(int which, int howMany)
{
    worlds[which]->updForceField()->setNumThreadsRequested(howMany);
}

void Context::setUseOpenMMAcceleration(bool arg)
{
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
        worlds[worldIx]->updForceField()->setUseOpenMMAcceleration(arg);
    }
}

/** Get/Set seed for reproducibility. **/
void Context::setSeed(int whichWorld, int whichSampler, unsigned long long int argSeed)
{
    worlds[whichWorld]->updSampler(whichSampler)->setSeed(argSeed);
}

unsigned long long int Context::getSeed(int whichWorld, int whichSampler)
{
    return worlds[whichWorld]->updSampler(whichSampler)->getSeed();
}



    //------------
//------------


/** Analysis related functions **/
void Context::addDistance(int whichWorld, int whichCompound, int aIx1, int aIx2)
{
    std::vector<int> tempV;
    tempV.push_back(whichWorld);
    tempV.push_back(whichCompound);
    tempV.push_back(aIx1);
    tempV.push_back(aIx2);

    distanceIxs.push_back(tempV);
}

void Context::addDihedral(int whichWorld, int whichCompound, int aIx1, int aIx2, int aIx3, int aIx4)
{
    std::vector<int> tempV;
    tempV.push_back(whichWorld);
    tempV.push_back(whichCompound);
    tempV.push_back(aIx1);
    tempV.push_back(aIx2);
    tempV.push_back(aIx3);
    tempV.push_back(aIx4);

    dihedralIxs.push_back(tempV);
}


// --- Printing functions --

// Print energy information
void Context::PrintSamplerData(unsigned int whichWorld)
{

    SimTK::State& currentAdvancedState = worlds[whichWorld]->integ->updAdvancedState();
/*    std::cout << currentAdvancedState.getNU() << ' '
        << worlds[whichWorld]->updSampler(0)->getAcceptedSteps() << ' '
        << std::setprecision(4) << std::fixed
        << worlds[whichWorld]->updSampler(0)->getOldPE() << ' '
        << worlds[whichWorld]->updSampler(0)->getSetPE() << ' '
        << worlds[whichWorld]->updSampler(0)->getLastAcceptedKE() << ' '
        << worlds[whichWorld]->updSampler(0)->getProposedKE() << ' '
        << worlds[whichWorld]->updSampler(0)->getOldFixman() << ' '
        << worlds[whichWorld]->updSampler(0)->getSetFixman() << ' '
        << worlds[whichWorld]->updSampler(0)->getProposedFixman() << ' '
        ;
*/
/*    // Use printf for faster output
    printf("%d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f "
        , currentAdvancedState.getNU() 
        , worlds[whichWorld]->updSampler(0)->getAcceptedSteps() 
        , worlds[whichWorld]->updSampler(0)->getOldPE() 
        , worlds[whichWorld]->updSampler(0)->getSetPE() 
        , worlds[whichWorld]->updSampler(0)->getLastAcceptedKE() 
        , worlds[whichWorld]->updSampler(0)->getProposedKE() 
        , worlds[whichWorld]->updSampler(0)->getOldFixman() 
        , worlds[whichWorld]->updSampler(0)->getSetFixman() 
        , worlds[whichWorld]->updSampler(0)->getProposedFixman() 
    );
*/

/*    // Avoid get function calls
    printf("%d %d %.2f %.2f %.2f %.2f %.2f %.2f"
        , currentAdvancedState.getNU()
        , (worlds[whichWorld]->samplers[0])->acceptedSteps
        , (worlds[whichWorld]->samplers[0])->pe_o
        , (worlds[whichWorld]->samplers[0])->pe_set
        , (worlds[whichWorld]->samplers[0])->ke_proposed
        , (worlds[whichWorld]->samplers[0])->ke_n
        , (worlds[whichWorld]->samplers[0])->fix_o
        , (worlds[whichWorld]->samplers[0])->fix_set
    );
*/

    // Write to a file instead of stdout
    fprintf(logFile, "%d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f "
        , currentAdvancedState.getNU()
        , pHMC(worlds[whichWorld]->samplers[0])->acceptedSteps
        , pMC((worlds[whichWorld]->samplers[0]))->pe_o
        , pMC((worlds[whichWorld]->samplers[0]))->pe_set
        , pHMC((worlds[whichWorld]->samplers[0]))->ke_proposed
        , pHMC((worlds[whichWorld]->samplers[0]))->ke_n
        , pMC((worlds[whichWorld]->samplers[0]))->fix_o
        , pMC((worlds[whichWorld]->samplers[0]))->fix_n
        , pMC((worlds[whichWorld]->samplers[0]))->fix_set
    );

}

// Print geometric parameters during simulation
void Context::PrintGeometry(SetupReader& setupReader, int whichWorld)
{
    if(setupReader.get("GEOMETRY")[0] == "TRUE"){
        // Get distances indeces
        int distanceIx[setupReader.get("DISTANCE").size()];
        for(unsigned int i = 0; i < setupReader.get("DISTANCE").size(); i++){
            distanceIx[i] = atoi(setupReader.get("DISTANCE")[i].c_str());
        }
        // Get distances
        for(int ai = 0; ai < (setupReader.get("DISTANCE").size() / 2); ai++){
            /*
            std::cout << std::setprecision(4) 
            << this->Distance(whichWorld, 0, 0, 
                distanceIx[2*ai + 0], distanceIx[2*ai + 1]) << " ";
            */
            
            printf("%.2f ", this->Distance(whichWorld, 0, 0,
                 distanceIx[2*ai + 0], distanceIx[2*ai + 1]));

        }

        // Get dihedrals indeces
        int dihedralIx[setupReader.get("DIHEDRAL").size()];
        for(unsigned int i = 0; i < setupReader.get("DIHEDRAL").size(); i++){
            dihedralIx[i] = atoi(setupReader.get("DIHEDRAL")[i].c_str());
        }
        // Get dihedrals
        for(int ai = 0; ai < (setupReader.get("DIHEDRAL").size() / 4); ai++){
            /*
            std::cout << std::setprecision(4) 
            << this->Dihedral(whichWorld, 0, 0, 
                dihedralIx[4*ai + 0], dihedralIx[4*ai + 1], 
                dihedralIx[4*ai + 2], dihedralIx[4*ai + 3]) << " ";
            */

            printf("%.2f ", this->Dihedral(whichWorld, 0, 0,
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

void Context::PrintGeometry(int whichWorld)
{

    for ( unsigned int i = 0; i < distanceIxs.size(); i++){
        if( distanceIxs[i][0] == whichWorld){
            fprintf(logFile, "%.3f ", 
                this->Distance(distanceIxs[i][0], distanceIxs[i][1], 0, 
                   distanceIxs[i][2], distanceIxs[i][3]) );
        }
    }

    for ( unsigned int i = 0; i < dihedralIxs.size(); i++){
        if( dihedralIxs[i][0] == whichWorld){
            fprintf(logFile, "%.3f ", this->Dihedral(dihedralIxs[i][0], dihedralIxs[i][1], 0,
                dihedralIxs[i][2], dihedralIxs[i][3], 
                dihedralIxs[i][4], dihedralIxs[i][5]) );
        }
    }
}

void Context::PrintDistances(int whichWorld)
{

    for ( unsigned int i = 0; i < distanceIxs.size(); i++){
        if( distanceIxs[i][0] == whichWorld){
            fprintf(logFile, "%.3f ", 
                this->Distance(distanceIxs[i][0], distanceIxs[i][1], 0, 
                   distanceIxs[i][2], distanceIxs[i][3]) );
        }
    }
}

void Context::PrintDihedrals(int whichWorld)
{
    for ( unsigned int i = 0; i < dihedralIxs.size(); i++){
        if( dihedralIxs[i][0] == whichWorld){
            fprintf(logFile, "%.3f ", this->Dihedral(dihedralIxs[i][0], dihedralIxs[i][1], 0,
                dihedralIxs[i][2], dihedralIxs[i][3], 
                dihedralIxs[i][4], dihedralIxs[i][5]) );
        }
    }
}

void Context::PrintDihedralsQs(int whichWorld)
{
    for ( unsigned int i = 0; i < dihedralIxs.size(); i++){
        if( dihedralIxs[i][0] == whichWorld){
            //fprintf(logFile, "%.3f ", this->Dihedral(dihedralIxs[i][0], dihedralIxs[i][1], 0,
            //    dihedralIxs[i][2], dihedralIxs[i][3], 
            //    dihedralIxs[i][4], dihedralIxs[i][5]) );

            const Topology& topology = worlds[whichWorld]->getTopology(dihedralIxs[i][1]);
            SimTK::State& currentAdvancedState = worlds[whichWorld]->integ->updAdvancedState();

            SimTK::MobilizedBodyIndex mbx3 = topology.getAtomMobilizedBodyIndex(
                SimTK::Compound::AtomIndex(dihedralIxs[i][4]) );
            SimTK::MobilizedBody::Pin& mobod3 = (SimTK::MobilizedBody::Pin&) (worlds[whichWorld]->matter->updMobilizedBody(mbx3));
           
            //std::cout << mbx3 << std::endl ;
            //std::cout << currentAdvancedState.getQ() << std::endl; 
            //fprintf(logFile, "%.3f ", currentAdvancedState.getQ()[mbx3] );
            fprintf(logFile, "%.3f ", mobod3.getQ(currentAdvancedState) );

        }
    }
}

void Context::PrintFreeE2EDist(int whichWorld, int whichCompound)
{
    const Topology& topology = worlds[whichWorld]->getTopology(whichCompound);
    SimTK::State& currentAdvancedState = worlds[whichWorld]->integ->updAdvancedState();

    for ( unsigned int i = 0; i < distanceIxs.size(); i++){
        if( distanceIxs[i][0] == whichWorld){

            fprintf(logFile, "%.3f ", 
                this->Distance(distanceIxs[i][0], distanceIxs[i][1], 0, 
                   distanceIxs[i][2], distanceIxs[i][3]) );

            SimTK::MobilizedBodyIndex mbx1 = topology.getAtomMobilizedBodyIndex(
                SimTK::Compound::AtomIndex(distanceIxs[i][2]) );
            SimTK::MobilizedBodyIndex mbx2 = topology.getAtomMobilizedBodyIndex(
                SimTK::Compound::AtomIndex(distanceIxs[i][3]) );
            SimTK::MobilizedBody& mobod1 = worlds[whichWorld]->matter->updMobilizedBody(mbx1);
            SimTK::MobilizedBody& mobod2 = worlds[whichWorld]->matter->updMobilizedBody(mbx2);
            SimTK::Transform X_PF1 = mobod1.getInboardFrame(currentAdvancedState);
            SimTK::Transform X_PF2 = mobod2.getInboardFrame(currentAdvancedState);
            SimTK::Transform X_BM1 = mobod1.getOutboardFrame(currentAdvancedState);
            SimTK::Transform X_BM2 = mobod2.getOutboardFrame(currentAdvancedState);
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
int Context::getPdbRestartFreq(void)
{
    return this->pdbRestartFreq;
}

// 
void Context::setPdbRestartFreq(int argFreq)
{
    this->pdbRestartFreq = argFreq;
}

// Write pdb
void Context::WritePdb(int whichWorld){assert(!"Not implemented");}


// Get / set printing frequency
int Context::getPrintFreq(void)
{
    return this->printFreq;
}

// 
void Context::setPrintFreq(int argFreq)
{
    this->printFreq = argFreq;
}

std::string Context::getOutputDir(void)
{
    return this->outputDir;
}

void Context::setOutputDir(std::string arg)
{
    this->outputDir = arg;
}

std::string Context::getPdbPrefix(void)
{
    return this->pdbPrefix;
}

void Context::setPdbPrefix(std::string arg)
{
    this->pdbPrefix = arg;
}


SimTK::Real Context::Dihedral(int whichWorld, int whichCompound, int whichSampler, int a1, int a2, int a3, int a4){

    //SimTK::State& currentAdvancedState = (updWorld(whichWorld))->integ->updAdvancedState();

    SimTK::State& currentAdvancedState = worlds[whichWorld]->integ->updAdvancedState();

    //SimTK::State& currentAdvancedState = (worlds[whichWorld]->updSampler(whichSampler)->updTimeStepper()->updIntegrator()).updAdvancedState();

    const Topology& topology = worlds[whichWorld]->getTopology(whichCompound);
    SimTK::Vec3 a1pos, a2pos, a3pos, a4pos;
    a1pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
    a2pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));
    a3pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)));
    a4pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a4)));

    //std::cout << " poss: " << a1pos << ' ' << a2pos << ' ' << a3pos << ' ' << a4pos << ' ';
    //std::cout << " dih: "  << bDihedral(a1pos, a2pos, a3pos, a4pos) << '|' ;

    return bDihedral(a1pos, a2pos, a3pos, a4pos);

}

SimTK::Real Context::Distance(int whichWorld, int whichCompound, int whichSampler, int a1, int a2){

    SimTK::State& currentAdvancedState = worlds[whichWorld]->integ->updAdvancedState();

    const Topology& topology = worlds[whichWorld]->getTopology(whichCompound);
    SimTK::Vec3 a1pos, a2pos;
    a1pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
    a2pos = topology.calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));

    return (a1pos - a2pos).norm();

}

// Writeble reference to a samplers advanced state
SimTK::State& Context::updAdvancedState(int whichWorld, int whichSampler)
{
    return (pHMC(worlds[whichWorld]->updSampler(whichSampler))->updTimeStepper()->updIntegrator()).updAdvancedState();
}

// Realize Topology Stage for all the Worlds
void Context::realizeTopology() {
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++) {
        (worlds[worldIx]->getCompoundSystem())->realizeTopology();
    }
}

/** Print the number of threads each World got **/
void Context::PrintNumThreads(void) {
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++) {
        std::cout << "World " << worldIx  << " requested "
            << worlds[worldIx]->updForceField()->getNumThreadsRequested()
            << " and got "
            << worlds[worldIx]->updForceField()->getNumThreadsInUse()
            << std::endl;
    }
}




//------------




