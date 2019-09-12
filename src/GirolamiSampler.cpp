/**@file
Implementation of GirolamiSampler class. **/

#include "GirolamiSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

//** Constructor **/
GirolamiSampler::GirolamiSampler(SimTK::CompoundSystem *argCompoundSystem,
                                     SimTK::SimbodyMatterSubsystem *argMatter,
                                     SimTK::Compound *argResidue,
                                     SimTK::DuMMForceFieldSubsystem *argDumm,
                                     SimTK::GeneralForceSubsystem *argForces,
                                     SimTK::TimeStepper *argTimeStepper)
    : HamiltonianMonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
    , MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
    , Sampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
{
    this->useFixmanTorque = true;  
}

/** Destructor **/
GirolamiSampler::~GirolamiSampler()
{
}

/** Seed the random number generator. Set simulation temperature,
velocities to desired temperature, variables that store the configuration
and variables that store the energies, both needed for the
acception-rejection step. Also realize velocities and initialize
the timestepper. **/
void GirolamiSampler::initialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman)
{
    // Seed the random number generator
    randomEngine.seed( std::time(0) );

    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    timeStepper->initialize(compoundSystem->getDefaultState());

    // Set the simulation temperature
    setTemperature(argTemperature); // Needed for Fixman

    // Store the configuration
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

    // Store potential energies
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    // Store Fixman potential
    this->useFixman = argUseFixman;
    if(useFixman){
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }

    // Initialize velocities to temperature
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);

    // Store kinetic energies
    setProposedKE(matter->calcKineticEnergy(someState));
    setLastAcceptedKE(getProposedKE());

    // Store total energies
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman();
    this->etot_set = this->etot_proposed;
  
}

/** Same as initialize **/
void GirolamiSampler::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature)
{
    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

    // Set the simulation temperature
    setTemperature(argTemperature); // Needed for Fixman

    // Store the configuration
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

    // Store potential energies
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    // Store Fixman potential
    if(useFixman){
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }

    // Initialize velocities to temperature
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);

    // Store kinetic energies
    setProposedKE(matter->calcKineticEnergy(someState));
    setLastAcceptedKE(getProposedKE());

    // Store total energies
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman();
    this->etot_set = this->etot_proposed;

}

/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. TODO: break in two functions:
initializeVelocities and propagate/integrate **/
void GirolamiSampler::propose(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    // Seed the random number generator every move
    //randomEngine.seed(4294653137UL); // for reproductibility

    // Initialize configuration - not necessary unless we modify the
    // configuration in addition to velocities
    system->realize(someState, SimTK::Stage::Position);
    int t = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        TVector[t] = SetTVector[t];
        t++;
    }

    // TODO: change the names from Old to Proposed and Set to lastAccepted
    setOldPE(getSetPE());
    setOldFixman(getSetFixman());

    // Initialize velocities according to the Maxwell-Boltzmann distribution
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int i=0; i < nu; ++i){
        V[i] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);

    // TODO: Implement this in a different function
    //SimTK::Real temperatureBoost = 2.58; // sqrt(2000/300) : brings temperature from 300 to 1000
    //SqrtMInvV *= sqrtRT * temperatureBoost; // Set stddev according to temperature

    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;

    // TODEL
    // TODO: Implement this in a different function
////
////    // Unlock all mobilizers
////    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
////        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
////        mobod.unlock(someState);
////    }
////
////    // Randomly choose what mobods to lock and store the chosen ones into
////    // a binary array
////    int BinaryArray[matter->getNumBodies()];
////    for(int i = 0; i < matter->getNumBodies(); i++){
////        BinaryArray[i] = 1;
////    }
////    for(int i = 0; i < 3; i++){
////        SimTK::Real rand_no = uniformRealDistribution(randomEngine);
////        BinaryArray[int(std::floor(matter->getNumBodies() * rand_no))] = 0;
////    }
////    
////    // Lock the chosen mobilizers
////    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
////        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
////        if(BinaryArray[int(mbx) - 1] == 0){
////            SimTK::Motion::Level motionLevel = SimTK::Motion::Level::Position;
////            mobod.lock(someState, motionLevel);
////        }
////    }
    // END TODEL

    // Store the proposed kinetic energy
    system->realize(someState, SimTK::Stage::Velocity);
    setProposedKE(matter->calcKineticEnergy(someState));

    // Store the proposed total energy
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman();

    // TODEL
////    std::cout << "Qs and Us before stepTo:" << std::endl;
////    PrintBigMat(someState.getQ(), someState.getNQ(), 3, "Q");
////    PrintBigMat(someState.getU(), someState.getNU(), 3, "U");
    // END TODEL

    // Integrate (propagate trajectory)
    this->timeStepper->stepTo(someState.getTime() + (timestep*nosteps));

    // TODEL
////    std::cout << "Qs and Us after stepTo:" << std::endl;
////    PrintBigMat(someState.getQ(), someState.getNQ(), 3, "Q");
////    PrintBigMat(someState.getU(), someState.getNU(), 3, "U");
    // END TODEL

}

/** Main function that contains all the 3 steps of HMC.
Implements the acception-rejection step and sets the state of the
compound to the appropriate conformation wether it accepted or not. **/
void GirolamiSampler::update(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);

    // Do a trial move
    propose(someState, timestep, nosteps);

    // Get needed energies
    SimTK::Real pe_o  = getOldPE();
    if(useFixman){
        SimTK::Real fix_o = getOldFixman();
    }
    SimTK::Real ke_proposed  = getProposedKE();

    if(useFixman){
        fix_n = calcFixman(someState);
    }else{
        fix_n = 0.0;
    }
    SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM

    system->realize(someState, SimTK::Stage::Velocity);
    SimTK::Real ke_n = matter->calcKineticEnergy(someState);

    SimTK::Real etot_proposed, etot_n;
    if(useFixman){
        etot_n = pe_n + ke_n + fix_n;
        etot_proposed = pe_o + ke_proposed + fix_o;
    }else{
        etot_n = pe_n + ke_n;
        etot_proposed = pe_o + ke_proposed;
    }

    etot_proposed;
    etot_n;


    std::cout<<std::setprecision(5)<<std::fixed;
    std::cout << "pe_o " << pe_o + getREP() << " ke_o " << ke_proposed << " fix_o " << fix_o << " rep " << getREP()
        << " pe_n " << pe_n  + getREP() << " ke_n " << ke_n << " fix_n " << fix_n
        //<< " rand_no " << rand_no << " RT " << RT << " exp(-(etot_n - etot_proposed) " << exp(-(etot_n - etot_proposed) / RT)
        //<< " etot_n " << etot_n  + getREP() << " etot_proposed " << etot_proposed + getREP()
        ;

    // Apply Metropolis criterion
////    if(1){ // Always accept // TODO
    //assert(!std::isnan(pe_n));
    if ( (!std::isnan(pe_n)) || (etot_n < etot_proposed) || (rand_no < exp(-(etot_n - etot_proposed)/RT)) ){ // Accept
        std::cout << " acc 1 " ;
        setSetTVector(someState);
        //sendConfToEvaluator(); // OPENMM
        setSetPE(pe_n);
        setSetFixman(fix_n);
        setLastAcceptedKE(ke_n);
        //this->etot_set = getSetPE() + getSetFixman() + getLastAcceptedKE(); // TODO:seems wrong
        this->etot_set = getSetPE() + getSetFixman() + getProposedKE(); // TODO
    }else{ // Reject
        std::cout << " acc 0 " ;
        assignConfFromSetTVector(someState);
    }

    std::cout << " pe_os " << getSetPE() + getREP() << " ke_os " << getLastAcceptedKE() << " fix_os " << getSetFixman()
        //<< " pe_n " << pe_n << " ke_n " << ke_n << " fix_n " << fix_n
        << std:: endl;
    //std::cout << "Number of times the force field was evaluated: " << dumm->getForceEvaluationCount() << std::endl;

    // Keep track of how many MC trials have been done 
    ++nofSamples;


}








