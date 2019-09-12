/**@file
Implementation of Sampler class. **/

#include "Sampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

Sampler::Sampler(SimTK::CompoundSystem *argCompoundSystem,
                 SimTK::SimbodyMatterSubsystem *argMatter,
                 //Topology *argResidue,
                 SimTK::Compound *argResidue,
                 SimTK::DuMMForceFieldSubsystem *argDumm,
                 SimTK::GeneralForceSubsystem *argForces,
                 SimTK::TimeStepper *argTimeStepper)
{
    this->compoundSystem = argCompoundSystem;
    this->matter = argMatter;
    this->residue = argResidue;
    this->dumm = argDumm;
    this->forces = argForces;
    this->timeStepper = argTimeStepper;
    this->system = &(matter->getSystem());

}

// Destructor
Sampler::~Sampler(){
    ;
}

// Get set the seed
unsigned long long int Sampler::getSeed(void)
{
    return this->seed;
}

/** Store the value of the seed internally and also feed it to the random
number generator **/
void Sampler::setSeed(unsigned long long int argSeed)
{
    if(argSeed == 0){
        std::chrono::system_clock::time_point tp
            = std::chrono::system_clock::now();
        std::chrono::system_clock::duration dtn = tp.time_since_epoch();
        long clockSeed = dtn.count();
        randomEngine.seed( clockSeed );
        seed = clockSeed;
    }else{
        randomEngine.seed( argSeed );
        seed = argSeed;
    }

}

// Compute mass matrix determinant
SimTK::Real Sampler::calcMassDeterminant(SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Real D0 = 1.0;
    matter->calcDetM(state, V, DetV, &D0);
    return D0;
}

// Compute mass matrix determinant

SimTK::Real Sampler::calcMassDeterminant(const SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Real D0 = 1.0;
    matter->calcDetM(state, V, DetV, &D0);
    return D0;
}

/** Returns the number of MC trials done by this integrator. **/
int Sampler::getNofSamples(void){
    return nofSamples;
}


// Update - to be implemented by every specific sampler

//void Sampler::update(SimTK::State& somState){}

// TODO move
void Sampler::PrintSimbodyStateCache(SimTK::State& someState){
    std::cout << " System Stage: " << someState.getSystemStage() << std::endl;
    for(int i = 0; i < someState.getNumSubsystems(); i++){
        std::cout
            << " Subsystem Name: "
            << someState.getSubsystemName(SimTK::SubsystemIndex(i))
            << " Stage: "
            << someState.getSubsystemStage(SimTK::SubsystemIndex(i))
            << " Version: "
            << someState.getSubsystemVersion(SimTK::SubsystemIndex(i))
            << std::endl;
    }
}

void Sampler::initialize(SimTK::State& someState) {
    // Sampling
    int nofSamples = 0;
}

void Sampler::reinitialize(SimTK::State& someState) {
    // Sampling
    int nofSamples = 0;
}

// Getter for macroscopic temperature
SimTK::Real Sampler::getTemperature() const {
    return temperature;
}

/** Setter for macroscopic temperature. Also sets the RT **/
void Sampler::setTemperature(SimTK::Real temperature) {
    Sampler::temperature = temperature;
    RT = this->temperature * SimTK_BOLTZMANN_CONSTANT_MD;
}

SimTK::Real Sampler::getRT() const {
    return RT;
}

SimTK::Real Sampler::generateRandomNumber(GmolRandDistributionType distributionType) {
    if(distributionType == UNIFORM){
        return uniformRealDistribution(randomEngine);
    }else if(distributionType == NORMAL){
        return gaurand(randomEngine);
    }
}


