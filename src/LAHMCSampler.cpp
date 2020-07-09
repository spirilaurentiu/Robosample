/**@file
Implementation of LAHMCSampler class. **/

#include "LAHMCSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

//** Constructor **/
LAHMCSampler::LAHMCSampler(SimTK::CompoundSystem *argCompoundSystem
                                     ,SimTK::SimbodyMatterSubsystem *argMatter
                                     ,SimTK::Compound *argResidue
                                     ,SimTK::DuMMForceFieldSubsystem *argDumm
                                     ,SimTK::GeneralForceSubsystem *argForces
                                     ,SimTK::TimeStepper *argTimeStepper
                                     )
    : Sampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper),
    MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper),
    HMCSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
{
    this->useFixman = false;  
    this->fix_n = this->fix_o = 0.0;
    this->logSineSqrGamma2_n = this->logSineSqrGamma2_o = 0.0;
    this->residualEmbeddedPotential = 0.0;
    nofSamples = 0;
    this->alwaysAccept = false;
    this->timestep = 0.002; // ps
    this->temperature = 0.0;
    this->boostT = this->temperature;
    MDStepsPerSample = 0;
    proposeExceptionCaught = false;
    adaptTimestep = false;
}

/** Destructor **/
LAHMCSampler::~LAHMCSampler()
{
}

/** Calculate sqrt(M) using Eigen. For debug purposes. **/
//void LAHMCSampler::calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper)
//{
//    assert("!Not implemented");
//}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular matrix and it is computed 
multipling a set of orthonormal vectors with the sqrt(MInv). **/
//void LAHMCSampler::calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
//{
//    int nu = someState.getNU();
//    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvU: passed matrix doesn't have nu x nu size.");
//
//    SimTK::Vector V(nu);
//    SimTK::Vector SqrtMInvV(nu);
//
//    // Zero the matrix and the temporary vector
//    for (int i=0; i < nu; ++i){
//        V[i] = 0;
//        for (int j=0; j < nu; ++j){
//            SqrtMInv(i, j) = 0;
//        }
//    }
//
//    // Calculate the values inside the matrix
//    for (int i=0; i < nu; ++i){
//        V[i] = 1;
//        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
//        for (int j=0; j < nu; ++j){
//            SqrtMInv(i, j) = SqrtMInvV[j];
//        }
//        V[i] = 0;
//    }
//
//
//}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
This is lower triangular matrix and it is computed by multipling a set of
 orthonormal vectors with the sqrt(MInv) and transpose it. **/
//void LAHMCSampler::calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
//{
//    int nu = someState.getNU();
//    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvL: passed matrix doesn't have nu x nu size.");
//
//    SimTK::Vector V(nu);
//    SimTK::Vector SqrtMInvV(nu);
//
//    // Zero the matrix and the temporary vector
//    for (int i=0; i < nu; ++i){
//        V[i] = 0;
//        for (int j=0; j < nu; ++j){
//            SqrtMInv(i, j) = 0;
//        }
//    }
//
//    // Calculate the values inside the matrix
//    for (int i=0; i < nu; ++i){
//        V[i] = 1;
//        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
//        for (int j=0; j < nu; ++j){
//            SqrtMInv(j, i) = SqrtMInvV[j];
//        }
//        V[i] = 0;
//    }
//
//}

/** Stores the accepted kinetic energy. This should be set right after a
move is accepted. It's a component of the total energy stored. **/
void LAHMCSampler::setLastAcceptedKE(SimTK::Real inpKE)
{
    this->ke_lastAccepted = inpKE;
}

/** Sets the proposed kinetic energy before the proposal. This should be
set right after the velocities are initialized. **/
void LAHMCSampler::setProposedKE(SimTK::Real inpKE)
{
    this->ke_proposed = inpKE;
}

/** Get/set the TimeStepper that manages the integrator **/
//const SimTK::TimeStepper * LAHMCSampler::getTimeStepper(void)
//{
//    return timeStepper;
//}

//SimTK::TimeStepper * LAHMCSampler::updTimeStepper(void)
//{
//    return timeStepper;
//}

//void LAHMCSampler::setTimeStepper(SimTK::TimeStepper * someTimeStepper)
//{
//    timeStepper = someTimeStepper;
//}

/** Seed the random number generator. Set simulation temperature,
velocities to desired temperature, variables that store the configuration
and variables that store the energies, both needed for the
acception-rejection step. Also realize velocities and initialize
the timestepper. **/
void LAHMCSampler::initialize(SimTK::State& someState )
{
    // After an event handler has made a discontinuous change to the
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    timeStepper->initialize(compoundSystem->getDefaultState());

    // Set the simulation temperature
//r    setTemperature(argTemperature); // Needed for Fixman

    int nu = someState.getNU();

    // Randomize configuration
    //if(randomizeConformation == true){
    //    system->realize(someState, SimTK::Stage::Position);
    //    int nq = someState.getNQ();
    //    SimTK::Vector QV(nq);
    //    for (int j=7; j < nq; ++j){
    //        QV[j] = uniformRealDistribution_mpi_pi(randomEngine);
    //    }
    //    someState.updQ() = QV;
    //}
    //

    // Store the configuration
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

    // Store potential energies
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    // Store Fixman potential
//r    this->useFixman = argUseFixman;
    if(useFixman){
        std::cout << "Look Ahead Hamiltonian Monte Carlo sampler: using Fixman potential." << std::endl;
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());

        setOldLogSineSqrGamma2( ((Topology *)residue)->calcLogSineSqrGamma2(someState));
        setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());

        setOldLogSineSqrGamma2(0.0);
        setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
    }

    // Initialize velocities to temperature
    // TODO Shouldn't be here
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
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();
    this->etot_set = this->etot_proposed;

  
}

/** Same as initialize **/
//r void LAHMCSampler::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature) 
void LAHMCSampler::reinitialize(SimTK::State& someState)
{
     // After an event handler has made a discontinuous change to the
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

    // Set the simulation temperature
//r    setTemperature(argTemperature); // Needed for Fixman

    // Store the configuration
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

//std::cout << "reinitialize "
//<< dumm->CalcFullPotEnergyIncludingRigidBodies(someState) << std::endl;

    // Store potential energies
    //setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    // Store Fixman potential
    if(useFixman){
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());

        setOldLogSineSqrGamma2( ((Topology *)residue)->calcLogSineSqrGamma2(someState));
        setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());

        setOldLogSineSqrGamma2(0.0);
        setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
    }

    // Initialize velocities to temperature
    // TODO Shouldn't be here
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
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();
    this->etot_set = this->etot_proposed;

}

/** Get/Set the timestep for integration **/
//float LAHMCSampler::getTimestep(void)
//{
//    return timestep;
//    //return timeStepper->updIntegrator().getPredictedNextStepSize();
//}

//void LAHMCSampler::setTimestep(float argTimestep)
//{
//    if(argTimestep <= 0){
//        adaptTimestep = true;
//        this->timestep = 0.0002;
//	this->timestepIncr = 0.0001;
//    }else{
//        adaptTimestep = false;
//    	timeStepper->updIntegrator().setFixedStepSize(argTimestep);
//    	this->timestep = argTimestep;
//    }
//}

/** Get/Set boost temperature **/
//SimTK::Real LAHMCSampler::getBoostTemperature(void)
//{
//    return this->boostT;
//}
//
//void LAHMCSampler::setBoostTemperature(SimTK::Real argT)
//{
//    this->boostT = argT;
//    this->boostFactor = std::sqrt(this->boostT / this->temperature);
//    this->unboostFactor = 1 / boostFactor;
//    std::cout << "HMC: boost temperature: " << this->boostT << std::endl;
//    std::cout << "HMC: boost velocity scale factor: " << this->boostFactor << std::endl;
//}
//
//void LAHMCSampler::setBoostMDSteps(int argMDSteps)
//{
//    this->boostMDSteps = argMDSteps;
//    std::cout << "HMC: boost MD steps: " << this->boostMDSteps << std::endl;
//
//}
/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. TODO: break in two functions:
initializeVelocities and propagate/integrate **/
void LAHMCSampler::propose(SimTK::State& someState)
{
    // Initialize configuration - not necessary unless we modify the
    // configuration in addition to velocities
    system->realize(someState, SimTK::Stage::Position);

    int t = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        TVector[t] = SetTVector[t];
        t++;
    }

    // TODO: change the names from Old to Proposed and Set to lastAccepted
    pe_o = pe_set;
    fix_o = fix_set;
    logSineSqrGamma2_o = logSineSqrGamma2_set;

    // Initialize velocities according to the Maxwell-Boltzmann distribution
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int i=0; i < nu; ++i){
        V[i] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);

    SqrtMInvV *= (sqrtRT); // Set stddev according to temperature

    // Raise the temperature
    someState.updU() = SqrtMInvV;

    // Realize velocity
    system->realize(someState, SimTK::Stage::Velocity);

    // Store the proposed energies
    // setProposedKE(matter->calcKineticEnergy(someState));
    this->ke_proposed = matter->calcKineticEnergy(someState);
    //std::cout << " ke before timestepping " << this->ke_proposed << std::endl;
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();



    try {
        // IF NOT PRINT EVERY STEP // Integrate (propagate trajectory)
	if(adaptTimestep){
		if( !(nofSamples % 30) ){ // Do it only so often
			//std::cout << "adapt ts " << acceptedSteps << ' ' << nofSamples << std::endl;
			// Compute average acceptance in the buffer
			int sum = 0, sqSum = 0;
			for (std::list<int>::iterator p = acceptedStepsBuffer.begin(); p != acceptedStepsBuffer.end(); ++p){
        			sum += (int)*p;
				sqSum += (((int)*p) * ((int)*p));
			}
    			float acceptance = float(sum) / float(acceptedStepsBufferSize);
			//float sqAcceptance = sqSum / acceptedStepsBufferSize;
			//float stdAcceptance = std::sqrt(sqAcceptance - (acceptance*acceptance));
			float stdAcceptance = 0.03;

			// Increase or decrease timestep
			if(acceptance > (0.651 + stdAcceptance)){
				timestep += timestepIncr;
			}else if(acceptance < (0.651 - stdAcceptance)){
				if(timestep > 0.0001){
					timestep -= timestepIncr;
				}
			}
			

		}
        	this->timeStepper->stepTo(someState.getTime() + (timestep*MDStepsPerSample));
	}else{
        	this->timeStepper->stepTo(someState.getTime() + (timestep*MDStepsPerSample));
	}
        system->realize(someState, SimTK::Stage::Position);

    }catch(const std::exception&){
        proposeExceptionCaught = true;
        //std::cout << "HMC_stepTo_or_realizePosition_threw_exception";

        int i = 0;
        for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
            const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
            mobod.setQToFitTransform(someState, SetTVector[i]);
            i++;
        }
        system->realize(someState, SimTK::Stage::Position);


    }

}

/** Main function that contains all the 3 steps of HMC.
Implements the acception-rejection step and sets the state of the
compound to the appropriate conformation wether it accepted or not. **/
bool LAHMCSampler::update(SimTK::State& someState, SimTK::Real newBeta)
{
    bool acc;
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);

    // Get new Fixman potential
    if(useFixman){
        fix_n = calcFixman(someState);
        logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
    }else{
        fix_n = 0.0;
        logSineSqrGamma2_n = 0.0;
    }

    // Get new kinetic energy
    system->realize(someState, SimTK::Stage::Velocity);
    ke_n = matter->calcKineticEnergy(someState);
    //std::cout << " ke after timestepping " << this->ke_n << std::endl;

    // Get new potential energy
    if ( getThermostat() == ANDERSEN ){
        pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
        //pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // ELIZA FULL
	    //detmbat_n = ((Topology *)residue)->calcLogDetMBAT(someState);
	    logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
    }
    else{
        pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
        //pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // ELIZA FULL
	    //detmbat_n = ((Topology *)residue)->calcLogDetMBAT(someState);
        logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
    }

    // Calculate total energy
    if(useFixman){
        etot_n = pe_n + ke_n + fix_n - (0.5 * RT * logSineSqrGamma2_n);
        etot_proposed = pe_o + ke_proposed + fix_o - (0.5 * RT * logSineSqrGamma2_o);
    }else{
        etot_n = pe_n + ke_n;
        etot_proposed = pe_o + ke_proposed;
    }

    PrintDetailedEnergyInfo(someState);


    // Decide and get a new sample
    if ( getThermostat() == ANDERSEN ){ // MD with Andersen thermostat
        acc = true;
        std::cout << " acc" << std::endl;
        setSetTVector(someState);
        pe_set = pe_n;
        fix_set = fix_n;
        logSineSqrGamma2_set = logSineSqrGamma2_n;
        ke_lastAccepted = ke_n;
        etot_set = pe_set + fix_set + ke_proposed + logSineSqrGamma2_set;

        //setBeta(newBeta);

        ++acceptedSteps;
        acceptedStepsBuffer.push_back(1);
        acceptedStepsBuffer.pop_front();
    }else{ // Apply Metropolis correction
        if ( (proposeExceptionCaught == false) &&
                (!std::isnan(pe_n)) && ((etot_n < etot_proposed) ||
             (rand_no < exp(-(etot_n - etot_proposed) * this->beta)))) { // Accept based on full energy
             //(rand_no < exp(-(pe_n - pe_o) * this->beta)))) { // Accept based on potential energy
             acc = true;
             std::cout << " acc" << std::endl;
             setSetTVector(someState);
             pe_set = pe_n;
             fix_set = fix_n;
             logSineSqrGamma2_set = logSineSqrGamma2_n;
             ke_lastAccepted = ke_n;
             etot_set = pe_set + fix_set + ke_proposed + logSineSqrGamma2_set;

             //setBeta(newBeta);

             ++acceptedSteps;
             acceptedStepsBuffer.push_back(1);
             acceptedStepsBuffer.pop_front();

        } else { // Reject
             acc = false;
             std::cout << " nacc" << std::endl;
             assignConfFromSetTVector(someState);
             proposeExceptionCaught = false;
             acceptedStepsBuffer.push_back(0);
             acceptedStepsBuffer.pop_front();
        }
    }

    ++nofSamples;
    return acc;

}

int LAHMCSampler::getMDStepsPerSample() const {
    return MDStepsPerSample;
}

void LAHMCSampler::setMDStepsPerSample(int mdStepsPerSample) {
    MDStepsPerSample = mdStepsPerSample;
}

/** Print detailed energy information **/
void LAHMCSampler::PrintDetailedEnergyInfo(SimTK::State& someState)
{

    std::cout << std::setprecision(5) << std::fixed;
    std::cout << "pe_o " << pe_o << " pe_n " << pe_n
        << " pe_nB " << getPEFromEvaluator(someState)
        << " ke_prop " << ke_proposed << " ke_n " << ke_n
        << " fix_o " << fix_o << " fix_n " << fix_n << " "
        << " logSineSqrGamma2_o " << logSineSqrGamma2_o << " logSineSqrGamma2_n " << logSineSqrGamma2_n << " "
        //<< " detmbat_n " << detmbat_n //<< " detmbat_o " << detmbat_o << " "
        << " ts " << timestep  << " exp(bdE) " << exp(-(etot_n - etot_proposed) / RT)
        << " etot_n " << etot_n  << " etot_proposed " << etot_proposed
        //<< std::endl
        ;
}


