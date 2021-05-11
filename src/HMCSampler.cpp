/**@file
Implementation of HMCSampler class. **/

#include "HMCSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"
#include "World.hpp"

#define CHECK_IF_NAN(n) \
	if(std::isnan(n)) { \
		std::cout << "\t[WARNING] invalid sample: " << #n << " is nan!\n"; \
		return false; \
	}

//** Constructor **/
HMCSampler::HMCSampler(World* argWorld, SimTK::CompoundSystem *argCompoundSystem,
	SimTK::SimbodyMatterSubsystem *argMatter,
	//SimTK::Compound *argResidue,
	std::vector<Topology> &argTopologies,
	SimTK::DuMMForceFieldSubsystem *argDumm,
	SimTK::GeneralForceSubsystem *argForces,
	SimTK::TimeStepper *argTimeStepper) :
		Sampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper)
		//, MonteCarloSampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper)
{
	// BEGIN SAMPLER
	assert(argCompoundSystem != nullptr);
	assert(argMatter != nullptr);
	assert(argDumm != nullptr);
	assert(argForces != nullptr);
	assert(argTimeStepper != nullptr);

	this->system = &argMatter->getSystem();

	//this->residue = argResidue;
	assert(topologies.size() > 0);
	this->residue = &topologies[0];

	// Set total number of atoms and dofs
	natoms = 0;
	for (const auto& topology: topologies){
		natoms += topology.getNumAtoms();
	}

	int ThreeFrom3D = 3;
	ndofs = natoms * ThreeFrom3D;
	// END SAMPLER

	TVector = std::vector<SimTK::Transform>(matter->getNumBodies());
	SetTVector = std::vector<SimTK::Transform>(matter->getNumBodies());
	proposeExceptionCaught = false;

	for(int i = 0; i < acceptedStepsBufferSize; i++){ 
		acceptedStepsBuffer.push_back(0);
	}

	this->useFixman = false;  
	this->fix_n = this->fix_o = 0.0;
	this->logSineSqrGamma2_n = this->logSineSqrGamma2_o = 0.0;
	this->residualEmbeddedPotential = 0.0;
	nofSamples = 0;
	this->alwaysAccept = false;

	this->prevPrevAcceptance = SimTK::NaN;
	this->prevAcceptance = SimTK::NaN;
	this->acceptance = SimTK::NaN;

	
	this->prevPrevTimestep = SimTK::NaN; // ps
	this->prevTimestep = SimTK::NaN; // ps
	this->timestep = 0;

	this->temperature = 0.0;
	this->boostT = this->temperature;
	MDStepsPerSample = 0;
	proposeExceptionCaught = false;
	shouldAdaptTimestep = false;

	dR.resize(ndofs, 0);
	dRdot.resize(ndofs, 0);
}

/** Destructor **/
HMCSampler::~HMCSampler()
{
}


// Use Fixman potential
void HMCSampler::useFixmanPotential(void)
{
    useFixman = true;
}

// Return true if use Fixman potential
bool HMCSampler::isUsingFixmanPotential(void) const
{
    return useFixman;
}

// Get Fixman potential
SimTK::Real HMCSampler::getOldFixman(void) const
{
    return this->fix_o;
}

// Set set Fixman potential
void HMCSampler::setSetFixman(SimTK::Real argFixman)
{
    this->fix_set = argFixman;
}

// Get set Fixman potential
SimTK::Real HMCSampler::getSetFixman(void) const
{
    return this->fix_set;
}

// Get the stored potential energy
SimTK::Real HMCSampler::getOldPE(void) const
{
    return this->pe_o;
}

// Set stored potential energy
void HMCSampler::setOldPE(SimTK::Real argPE)
{
    this->pe_o = argPE;
}

// Get the set potential energy
SimTK::Real HMCSampler::getSetPE(void) const
{
    return this->pe_set;
}

// Set set potential energy
void HMCSampler::setSetPE(SimTK::Real argPE)
{
    this->pe_set = argPE;
}

// Set/get Residual Embedded Potential
void HMCSampler::setREP(SimTK::Real inp)
{
    this->residualEmbeddedPotential = inp;
}

SimTK::Real HMCSampler::getREP(void) const
{
    return this->residualEmbeddedPotential;
}

// Set a thermostat
void HMCSampler::setThermostat(ThermostatName argThermostat){
    this->thermostat = argThermostat;
}
// Set a thermostat
void HMCSampler::setThermostat(std::string thermoName){
    thermoName.resize(thermoName.size());
    std::transform(thermoName.begin(), thermoName.end(), thermoName.begin(), ::tolower);

    if(thermoName == "none"){
        this->thermostat = ThermostatName::NONE;
    }else if(thermoName == "andersen"){
        this->thermostat = ThermostatName::ANDERSEN;
    }else if(thermoName == "berendsen"){
        this->thermostat = ThermostatName::BERENDSEN;
    }else if(thermoName == "langevin"){
        this->thermostat = ThermostatName::LANGEVIN;
    }else if(thermoName == "nose_hoover"){
        this->thermostat = ThermostatName::NOSE_HOOVER;
    }else{
        std::cerr << "Invalid argument: " << thermoName << '\n';
    }
}

// Set a thermostat
void HMCSampler::setThermostat(const char *argThermostat){
    setThermostat(std::string(argThermostat));
}

// Get the name of the thermostat
ThermostatName HMCSampler::getThermostat(void) const{
    return thermostat;
}

////////////////////////////////////
// FIMAN POTENTIAL RELATED
////////////////////////////////////


// Get the potential energy from an external source as far as the sampler
// is concerned - OPENMM has to be inserted here
SimTK::Real HMCSampler::getPEFromEvaluator(SimTK::State& someState){
        return forces->getMultibodySystem().calcPotentialEnergy(someState);
	// Eliza's potential energy's calculations including rigid bodies
	// internal energy
    //return dumm->CalcFullPotEnergyIncludingRigidBodies(someState);

}

SimTK::Real HMCSampler::calcFixman(SimTK::State& someState){
    int nu = someState.getNU();
    SimTK::Vector V(nu);

    system->realize(someState, SimTK::Stage::Position);
    matter->realizeArticulatedBodyInertias(someState); // Move in calcDetM ?

    // Get detM
    SimTK::Vector DetV(nu);
    SimTK::Real D0 = 1.0;

    // TODO: remove the request for Dynamics stage cache in SImbody files
    //std::cout << "MonteCarloSampler::calcFixman Stage: "<< matter->getStage(someState) << std::endl;
    matter->calcDetM(someState, V, DetV, &D0);

    assert(RT > SimTK::TinyReal);
    SimTK::Real result = 0.5 * RT * ( std::log(D0) - ((Topology *)residue)->calcLogDetMBATInternal(someState) );

    if(SimTK::isInf(result)){
        result = 0.0;
    }
    
    return result;
}

// Compute Fixman potential numerically
SimTK::Real HMCSampler::calcNumFixman(SimTK::State&){
    // args were SimTK::State& someState
    assert(!"Not implemented.");
    return 0;
    /*
    // Get M
    int nu = someState.getNU();
    SimTK::Matrix M(nu, nu);
    matter->calcM(someState, M);

    // Eigen M determinant
    Eigen::MatrixXd EiM(nu, nu);
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            EiM(i, j) = M(i, j);
        }
    }
    SimTK::Real EiDetM = EiM.determinant();
    assert(RT > SimTK::TinyReal);
    SimTK::Real result = 0.5 * RT * std::log(EiDetM);
    return result;
    */
}

// Set old Fixman potential
void HMCSampler::setOldFixman(SimTK::Real argFixman)
{
    this->fix_o = argFixman;
}

// Set/get External MBAT contribution potential
void HMCSampler::setOldLogSineSqrGamma2(SimTK::Real argX){
    this->logSineSqrGamma2_o = argX;
}

SimTK::Real HMCSampler::getOldLogSineSqrGamma2(void) const
{
    return this->logSineSqrGamma2_o;
}

// Set/get External MBAT contribution potential
SimTK::Real HMCSampler::getSetLogSineSqrGamma2(void) const
{
    return this->logSineSqrGamma2_set;
}

void HMCSampler::setSetLogSineSqrGamma2(SimTK::Real argX){
    this->logSineSqrGamma2_set = argX;
}

void HMCSampler::setProposedLogSineSqrGamma2(SimTK::Real argFixman)
{
    this->logSineSqrGamma2_n = argFixman;
}

////////////////////////////////////
// CONFIGURATION RELATED
////////////////////////////////////

// Stores the configuration into an internal vector of transforms TVector
// Get the stored configuration
SimTK::Transform * HMCSampler::getTVector(void)
{
    return &TVector[0];
}

// Stores the configuration into an internal vector of transforms TVector
void HMCSampler::setTVector(const SimTK::State& someState)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    TVector[i] = mobod.getMobilizerTransform(someState);
    i++;
  }
}

// Stores the configuration into an internal vector of transforms TVector
void HMCSampler::setTVector(SimTK::Transform *inpTVector)
{
    // TODO pointer parameter is bad
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    TVector[i] = inpTVector[i];
    i++;
  }
}

// Stores the set configuration into an internal vector of transforms TVector
void HMCSampler::setSetTVector(const SimTK::State& someState)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    SetTVector[i] = mobod.getMobilizerTransform(someState);
    i++;
  }
}

// Stores the configuration into an internal vector of transforms TVector
// Get the stored configuration
SimTK::Transform * HMCSampler::getSetTVector(void)
{
    return &SetTVector[0];
}

// Restores configuration from the internal set vector of transforms TVector
void HMCSampler::assignConfFromSetTVector(SimTK::State& someState)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    system->realize(someState, SimTK::Stage::Position);
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    mobod.setQToFitTransform(someState, SetTVector[i]);

    system->realize(someState, SimTK::Stage::Position);


    i++;
  }
  system->realize(someState, SimTK::Stage::Position);
}

/** Calculate sqrt(M) using Eigen. For debug purposes. **/
void HMCSampler::calcNumSqrtMUpper(SimTK::State&, SimTK::Matrix&) const
{
	// function signature was SimTK::State& someState, SimTK::Matrix& SqrtMUpper
    assert("!Not implemented");
}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular matrix and it is computed 
multipling a set of orthonormal vectors with the sqrt(MInv). **/
void HMCSampler::calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv) const
{
    const int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvU: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = SqrtMInvV[j];
        }
        V[i] = 0;
    }


}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
This is lower triangular matrix and it is computed by multipling a set of
 orthonormal vectors with the sqrt(MInv) and transpose it. **/
void HMCSampler::calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv) const
{
    const int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvL: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(j, i) = SqrtMInvV[j];
        }
        V[i] = 0;
    }

}

/** Stores the accepted kinetic energy. This should be set right after a
move is accepted. It's a component of the total energy stored. **/
void HMCSampler::setLastAcceptedKE(SimTK::Real inpKE)
{
    this->ke_lastAccepted = inpKE;
}

/** Sets the proposed kinetic energy before the proposal. This should be
set right after the velocities are initialized. **/
void HMCSampler::setProposedKE(SimTK::Real inpKE)
{
    this->ke_proposed = inpKE;
}

/** Get/set the TimeStepper that manages the integrator **/
const SimTK::TimeStepper * HMCSampler::getTimeStepper()
{
    return timeStepper;
}

SimTK::TimeStepper * HMCSampler::updTimeStepper()
{
    return timeStepper;
}

void HMCSampler::setTimeStepper(SimTK::TimeStepper * someTimeStepper)
{
    timeStepper = someTimeStepper;
}

/** Put generalized velocities scale factors into a fixed size array to avoid
 searching for them into a map every time the velocities are intialized **/
void HMCSampler::loadUScaleFactors(SimTK::State& someState)
{
    int nu = someState.getNU();
    UScaleFactors.resize(nu, 1);

    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        // const int mnu = mobod.getNumU(someState);
		const float scaleFactor = world->getMobodUScaleFactor(mbx);

        //std::cout << "RED ZONE mbx scaleFactor uIxes " << int(mbx) << ' ' << scaleFactor;
		for(SimTK::UIndex uIx = mobod.getFirstUIndex(someState); uIx < mobod.getFirstUIndex(someState) + mobod.getNumU(someState); uIx++ ){
        	//std::cout << ' ' << int(uIx) ;
			UScaleFactors[int(uIx)] = scaleFactor;
        }
        //std::cout << '\n';
    }
}

/** Seed the random number generator. Set simulation temperature,
velocities to desired temperature, variables that store the configuration
and variables that store the energies, both needed for the
acception-rejection step. Also realize velocities and initialize
the timestepper. **/
void HMCSampler::initialize(SimTK::State& someState)
{
    // After an event handler has made a discontinuous change to the
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    timeStepper->initialize(compoundSystem->getDefaultState());

    const int nu = someState.getNU();

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
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SetTVector[mbx - 1] = TVector[mbx - 1] = mobod.getMobilizerTransform(someState);
    }

    // Initialize QsBuffer with zeros
    int nq = matter->getNQ(someState);
    int totSize = QsBufferSize * nq;
    for(int i = 0; i < totSize; i++){ 
        //QsBuffer.push_back(SimTK::Vector(nq, SimTK::Real(0)));
        QsBuffer.push_back(SimTK::Real(0));
    }

    // Store potential energies
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    // Store Fixman potential
    if(useFixman){
        std::cout << "Hamiltonian Monte Carlo sampler: using Fixman potential.\n";
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

    loadUScaleFactors(someState);
}

/** Same as initialize **/
//r void HMCSampler::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature) 
void HMCSampler::reinitialize(SimTK::State& someState)
{
     // After an event handler has made a discontinuous change to the
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

    // Store the configuration
	// In our case, a MobilizedBody is an atom
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

    // Store potential energies
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

    loadUScaleFactors(someState);
}

/** Get/Set the timestep for integration **/
SimTK::Real HMCSampler::getTimestep() const
{
    return timestep;
}

void HMCSampler::setTimestep(SimTK::Real argTimestep)
{
    if(argTimestep <= 0){
        shouldAdaptTimestep = true;
    	timeStepper->updIntegrator().setFixedStepSize(-argTimestep);
        this->timestep = -argTimestep;
    }else{
        shouldAdaptTimestep = false;
    	timeStepper->updIntegrator().setFixedStepSize(argTimestep);
    	this->timestep = argTimestep;
    }
}

/** Get/Set boost temperature **/
SimTK::Real HMCSampler::getBoostTemperature()
{
    return this->boostT;
}

void HMCSampler::setBoostTemperature(SimTK::Real argT)
{
    this->boostT = argT;
    this->boostFactor = std::sqrt(this->boostT / this->temperature);
    this->unboostFactor = 1 / boostFactor;
    std::cout << "HMC: boost temperature: " << this->boostT << std::endl;
    std::cout << "HMC: boost velocity scale factor: " << this->boostFactor << std::endl;
}

void HMCSampler::setBoostMDSteps(int argMDSteps)
{
    this->boostMDSteps = argMDSteps;
    std::cout << "HMC: boost MD steps: " << this->boostMDSteps << std::endl;

}

/** Store configuration as a set of Transforms **/
void HMCSampler::storeOldConfigurationAndPotentialEnergies(SimTK::State& someState){ // func
    system->realize(someState, SimTK::Stage::Position);

    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        TVector[mbx - 1] = SetTVector[mbx - 1];
    }

    // TODO: change the names from Old to Proposed and Set to lastAccepted
    pe_o = pe_set;
    fix_o = fix_set;
    logSineSqrGamma2_o = logSineSqrGamma2_set;
} // func

/** Initialize velocities according to the Maxwell-Boltzmann
distribution.  Coresponds to R operator in LAHMC **/
void HMCSampler::initializeVelocities(SimTK::State& someState){
	// Check if we can use our cache
    const int nu = someState.getNU();
	if (nu != RandomCache.nu) {
		// Rebuild the cache
		// We also get here if the cache is not initialized
		RandomCache.V.resize(nu);
		RandomCache.SqrtMInvV.resize(nu);
		RandomCache.sqrtRT = std::sqrt(RT);
		RandomCache.nu = nu;

		// we don't get to use multithreading here
		RandomCache.FillWithGaussian();
	} else {
		// wait for random number generation to finish (should be done by this stage)
		RandomCache.task.wait();
	}

	// V[i] *= UScaleFactors[i] - note that V is already populated with random numbers 
	std::transform(RandomCache.V.begin(), RandomCache.V.end(), // apply an operation on this
		UScaleFactors.begin(), // and this
		RandomCache.V.begin(), // and store here
		std::multiplies<SimTK::Real>()); // this is the operation
	
	// Scale by square root of the inverse mass matrix
	matter->multiplyBySqrtMInv(someState, RandomCache.V, RandomCache.SqrtMInvV);

	// Set stddev according to temperature
	RandomCache.SqrtMInvV *= (RandomCache.sqrtRT);

	// Raise the temperature
	someState.updU() = RandomCache.SqrtMInvV;

	// Realize velocity
	system->realize(someState, SimTK::Stage::Velocity);

	// ask for a number of random numbers and check if we are done the next time we hit this function
	RandomCache.task = std::async(std::launch::async, RandomCache.FillWithGaussian);
}

/** Store the proposed energies **/
void HMCSampler::calcProposedKineticAndTotalEnergy(SimTK::State& someState){

    // Store proposed kinetic energy
    // setProposedKE(matter->calcKineticEnergy(someState));
    this->ke_proposed = matter->calcKineticEnergy(someState);

    // Store proposed total energy
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();
}

void HMCSampler::adaptTimestep(SimTK::State&)
{
	// function signature was SimTK::State& someState

	std::cout << std::endl;
	//std::cout << "Adapt: nofSamples= " << nofSamples << std::endl;
	if( (nofSamples % acceptedStepsBufferSize) == (acceptedStepsBufferSize-1) ){ // Do it only so often
		std::cout << "Adapt BEGIN: ts= " << timestep;

		// SimTK::Real stdAcceptance = 0.03;
		SimTK::Real idealAcceptance = 0.651;
		// SimTK::Real smallestAcceptance = 0.0001;
		// SimTK::Real timestepIncr = 0.001;

		// Compute acceptance in the buffer
		int sum = std::accumulate(acceptedStepsBuffer.begin(), acceptedStepsBuffer.end(), 0);

		SimTK::Real prevPrevPrevAcceptance = prevPrevAcceptance; // to reset
		prevPrevAcceptance = prevAcceptance;
		prevAcceptance = acceptance;
    	acceptance = sum / acceptedStepsBufferSize;

		std::cout << " ppAcc= " << prevPrevAcceptance << " pAcc= " << prevAcceptance << " acc= " << acceptance << std::endl;

		if( !SimTK::isNaN(prevPrevAcceptance) ){ // Passed first two initial evaluations
			// Calculate gradients
			SimTK::Real a_n, a_n_1, a_n_2, t_n, t_n_1, t_n_2;
			a_n = acceptance; a_n_1 = prevAcceptance; a_n_2 = prevPrevAcceptance;
			t_n = timestep; t_n_1 = prevTimestep; t_n_2 = prevPrevTimestep;

			SimTK::Real F_n =   std::abs(a_n   - idealAcceptance);
			SimTK::Real F_n_1 = std::abs(a_n_1 - idealAcceptance);
			SimTK::Real F_n_2 = std::abs(a_n_2 - idealAcceptance);

			SimTK::Real dF_n   = F_n     - F_n_1;
			SimTK::Real dF_n_1 = F_n_1   - F_n_2;
			SimTK::Real dt_n = (t_n   - t_n_1);
			SimTK::Real dt_n_1 = (t_n_1   - t_n_2);
			SimTK::Real gradF_n   = dF_n / dt_n;
			SimTK::Real gradF_n_1 = dF_n_1 / dt_n_1;

			// Calculate gamma
			SimTK::Real gamma;
			SimTK::Real num, denom;
			num = std::abs(dt_n * (gradF_n - gradF_n_1));
			denom = (gradF_n - gradF_n_1) * (gradF_n - gradF_n_1);
			gamma = num / denom;

			// Increase or decrease timestep
			SimTK::Real newTimestep = t_n - (gamma * gradF_n);
			std::cout << "Adapt data: "
				<< " " << F_n  << " " << F_n_1  << " " << F_n_2 << " " << dF_n << " " << dF_n_1 
				<< " " << dt_n << " " << dt_n_1 << " " << gradF_n << " " << gradF_n_1
				<< " " << num << " " << denom << " " << gamma
				<< std::endl;
			if( (std::abs(prevTimestep - timestep) > 0.0001) || (newTimestep < 0.00089) || (SimTK::isNaN(newTimestep)) || ((acceptance - prevAcceptance) < 0.001)){
				std::cout << "newTimestep rejected\n"; 
    				SimTK::Real r = uniformRealDistribution(randomEngine);
				if(acceptance < idealAcceptance){
					timestep = timestep - (0.3*r * timestep);
				}else{
					timestep = timestep + (0.3*r * timestep);
				}

				acceptance = prevAcceptance;
				prevAcceptance = prevPrevAcceptance;
				prevPrevAcceptance = prevPrevPrevAcceptance;
			} else { // Reset acceptances
				std::cout << "newTimestep accepted\n"; 
				prevPrevTimestep = prevTimestep;
				prevTimestep = timestep;
				timestep = newTimestep;
			}

    			timeStepper->updIntegrator().setFixedStepSize(timestep);
			std::cout << "Adapt END: "  << " ppTs= " << prevPrevTimestep << " pTs= " << prevTimestep << " ts= " << timestep << std::endl;
		} else { // Alter the intial timesteps to get a valid dF_n next time

		//SimTK::ArticulatedInertia abi = matter->getArticulatedBodyInertia(someState, SimTK::MobilizedBodyIndex(2));
	    	//const SimTK::MobilizedBody& mobod2 = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(2));
		//const SimTK::MassProperties mp2 = mobod2.getBodyMassProperties(someState);
		//const SimTK::Inertia i2 = mp2.calcInertia();
		//SimTK::Mat33 mi2 = i2.toMat33();
		//std::cout << std::endl;
		//std::cout << mi2(0,0) << " " << mi2(0,1) << " " << mi2(0,2) << std::endl;
		//std::cout << mi2(1,0) << " " << mi2(1,1) << " " << mi2(1,2) << std::endl;
		//std::cout << mi2(2,0) << " " << mi2(2,1) << " " << mi2(2,2) << std::endl;
		//std::cout << "det(2)^1/5 " << std::pow(SimTK::det(mi2), 0.2) << std::endl;

			prevPrevTimestep = prevTimestep;
			prevTimestep = timestep;
			if( SimTK::isNaN(prevPrevTimestep) ){
				timestep = 0.0009; // smaller than H vibration
    				timeStepper->updIntegrator().setFixedStepSize(timestep);
			}else{
				timestep = (prevTimestep + prevPrevTimestep) / 2;
    				timeStepper->updIntegrator().setFixedStepSize(timestep);
			}
			std::cout << "Adapt END: ppTs= " << prevPrevTimestep << " pTs= " << prevTimestep << " ts= " << timestep << std::endl;
		}

	} // is time to adapt 
}

/**  **/
void HMCSampler::adaptWorldBlocks(SimTK::State& someState){

    std::cout << std::setprecision(10) << std::fixed;
	int nq = someState.getNQ();
	// int totSize = QsBufferSize * nq;
	if( (nofSamples % QsBufferSize) == (QsBufferSize-1) ){
		std::cout << "Adapt blocks BEGIN: \n";
		//std::cout << "Print by row: \n";
		//for(int confi = 0; confi < QsBufferSize; confi++){
		//	for(int qi = 0; qi < nq; qi++){
		//		SimTK::Real x;
		//		int pos;
		//
		//		std::list<SimTK::Real>::iterator it = QsBuffer.begin();
		//		pos = (confi * nq) + qi;
		//		std::advance(it, pos);
		//		std::cout << *it << ' ';
		//	}
		//	std::cout << std::endl;
		//}

		//std::cout << "Print by column: \n";
		std::vector<std::vector<SimTK::Real>> QsBufferVec(nq, vector<SimTK::Real>(QsBufferSize));
		for(int qi = 0; qi < nq; qi++){

			std::vector<SimTK::Real> thisQ(QsBufferSize);

			//std::cout << "QsBuffer= \n";
			for(int confi = 0; confi < QsBufferSize; confi++){
				thisQ[confi] = QsBuffer[(confi * nq) + qi];
				//std::cout << QsBuffer[(confi * nq) + qi] << ' ';
			}
			//std::cout << std::endl;

			//std::cout << "thisQ= \n";
			//for(int confi = 0; confi < QsBufferSize; confi++){
			//	std::cout << thisQ[confi] << ' ';
			//}
			//std::cout << std::endl;

			QsBufferVec[qi] = thisQ;

		}

		//std::cout << "QsBufferVec= \n";
		//for(int qi = 0; qi < nq; qi++){
		//	for(int confi = 0; confi < QsBufferSize; confi++){
		//		std::cout << QsBufferVec[qi][confi] << ' ';
		//	}
		//	std::cout << std::endl;
		//}
		
		std::cout << "Compute correlation matrix: \n";
		SimTK::Real corr;
		for(int qi = 0; qi < nq; qi++){
			//for(int qj = 0; (qj < nq) && (qj < qi); qj++){
			for(int qj = 0; (qj < nq); qj++){
				if(qIndex2jointType[SimTK::QIndex(qi)] == JointType::ANGULAR360){
					if(qIndex2jointType[SimTK::QIndex(qj)] == JointType::ANGULAR360){
						corr = circCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}else if(qIndex2jointType[SimTK::QIndex(qj)] == JointType::LINEAR){
						corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}else{
						corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}
				}else if(qIndex2jointType[SimTK::QIndex(qi)] == JointType::LINEAR){
					if(qIndex2jointType[SimTK::QIndex(qj)] == JointType::LINEAR){
						corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}else if(qIndex2jointType[SimTK::QIndex(qj)] == JointType::ANGULAR360){
						corr = circCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}else{
						corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}
				}else{
					corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
					std::cout << corr << ' ';
				}
			}
			std::cout << std::endl;
		}

		std::cout << "Adapt blocks END: \n";
	}

}


/** Apply the L operator **/
void HMCSampler::integrateTrajectory(SimTK::State& someState){
    try {
        this->timeStepper->stepTo(someState.getTime() + (timestep*MDStepsPerSample));

        system->realize(someState, SimTK::Stage::Position);

    }catch(const std::exception&){
		std::cout << "\t[WARNING] propose exception caught!\n";
        proposeExceptionCaught = true;

        for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
            const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
            mobod.setQToFitTransform(someState, SetTVector[mbx - 1]);
        }

        system->realize(someState, SimTK::Stage::Position);
    }
}

/**  **/
void HMCSampler::integrateTrajectoryOneStepAtATime(SimTK::State& someState
	//, std::function<void()> initialFunc
	//, std::function<void()> loopFunc
	//, std::function<void()> finalFunc
	){

	//// TODO: Intial function BEGIN ////
	//// Push initial R and Rdots
	//SimTK::Vec3 atomR;
	//SimTK::Vec3 atomRdot;
	//SimTK::Compound::AtomIndex aIx; 
	//for(int i = 0; i < topologies.size(); i++){
	//	for(int j = 0; j < (topologies[i])->getNumAtoms(); j++){
	//		aIx = ((topologies[i])->bAtomList[j]).getCompoundAtomIndex();
	//		atomR = (topologies[i])->calcAtomLocationInGroundFrame(someState, aIx);
	//		R.push_back(atomR[0]);
	//		R.push_back(atomR[1]);
	//		R.push_back(atomR[2]);
	//	}
	//}
	//// Initial function END ////

	//// Main loop ////
	for(int k = 0; k < MDStepsPerSample; k++){
		this->timeStepper->stepTo(someState.getTime() + (timestep));
                system->realize(someState, SimTK::Stage::Position);

		//// Calculate generalized quantities
    		//int nu = someState.getNU();
    		//SimTK::Vector MU(nu);
    		//for (int i=0; i < nu; ++i){
        	//	MU[i] = 1;
		//}
		//matter->multiplyByM(someState, someState.getU(), MU);
		//std::cout << "MU= " << MU << std::endl;
		
    		//int nu = someState.getNU();
    		//SimTK::Vector U(nu);
		//U = someState.getU();
		//std::cout << "U= " << U << std::endl;

		/////////////////////////


		//// TODO:Loop function BEGIN ////
		//	// Push current R and Rdot
		//	for(int i = 0; i < topologies.size(); i++){
		//		for(int j = 0; j < (topologies[i])->getNumAtoms(); j++){
		//			aIx = ((topologies[i])->bAtomList[j]).getCompoundAtomIndex();
		//			atomR = (topologies[i])->calcAtomLocationInGroundFrame(someState, aIx);
		//			R.push_back(atomR[0]);
		//			R.push_back(atomR[1]);
		//			R.push_back(atomR[2]);
		//			atomRdot = (topologies[i])->calcAtomVelocityInGroundFrame(someState, aIx);
		//			Rdot.push_back(atomRdot[0]);
		//			Rdot.push_back(atomRdot[1]);
		//			Rdot.push_back(atomRdot[2]);
		//		}
		//	}
		//
		//	// Calculate dR		
		//	for(unsigned int j = 0; j < (ndofs); j++){
		//		dR[j] = R[j + (ndofs)] - R[j];
		//	}
		//
		//	// Calculate MSD
		//	SimTK::Real MSD = 0;
		//	if(dR.size() >= (ndofs)){
		//		MSD += magSq(dR);
		//		MSD /= (ndofs);
		//	}
		//
		//	// Calculate RRdot	
		//	SimTK::Real RRdot = 0;
		//	if(dR.size() >= (ndofs)){
		//		std::vector<SimTK::Real> tempDR = dR;
		//		std::vector<SimTK::Real> tempRdot = Rdot;
		//
		//		normalize(tempDR);
		//		normalize(tempRdot);
		//
		//		for(unsigned int j = 0; j < ndofs; j++){
		//			RRdot += tempDR[j] * tempRdot[j];
		//		}
		//	}
		//	std::cout << std::setprecision(10) << std::fixed;
		//	std::cout << "k= " << k << " RRdot= " << RRdot
		//		<< " MSD= " << MSD << std::endl;
		//
		//	// Pop current R and Rdot
		//	for(int j = ((ndofs) - 1); j >= 0; --j){
		//		R.pop_back();
		//		Rdot.pop_back();
		//	}
		//// Loop function END ////

	} // END main loop

	//// TODO: Final function BEGIN ////
	//// Pop initial R
	//for(int j = ((ndofs) - 1); j >= 0; --j){
	//	R.pop_back();
	//}
	//// Final function END ////

}

void HMCSampler::geomDihedral(){
    /* // INSTANT GEOMETRY
    SimTK::Vec3 a1pos, a2pos, a3pos, a4pos, a5pos;
    int a1, a2, a3, a4, a5;
    //DIHEDRAL 4 6 8 14 6 8 14 16
    //a1 = 16; a2 = 14; a3 = 0; a4 = 6; a5 = 8;
    a1 = 4; a2 = 6; a3 = 8; a4 = 14; a5 = 16;
    a1pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
    a2pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));
    a3pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)));
    a4pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a4)));
    a5pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a5)));
    int distA1, distA2, distA3, distA4;
    SimTK::Vec3 distA1pos, distA2pos;
    SimTK::Vec3 distA3pos, distA4pos;
    distA1 = 2; distA2 = 17; distA3 = 6; distA4 = 17;
//		    distA1pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA1)));
//		    distA2pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA2)));
//		    distA3pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA3)));
//		    distA4pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA4)));
//            std::cout << " dihedral elements"
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a1))
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a2))
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a3))
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a4))
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a5))
//              << std::endl;
//            std::cout << " poss: " << a1pos << ' ' << a2pos << ' ' << a3pos << ' ' << a4pos << ' ';
    std::cout << "geom "  << bDihedral(a1pos, a2pos, a3pos, a4pos) ;
    std::cout << " "  << bDihedral(a2pos, a3pos, a4pos, a5pos) ;
    //std::cout << " "  << 10 * (distA1pos - distA2pos).norm() ; // Ang
    //std::cout << " "  << 10 * (distA3pos - distA4pos).norm() ; // Ang
    std::cout << std::endl;
    // */
    // INSTANT GEOMETRY END
}

/** Store new configuration and energy terms**/
void HMCSampler::calcNewConfigurationAndEnergies(SimTK::State& someState)
{
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

    // Get new potential energy
    pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
    // TODO: replace with the following after checking is the same thing
    //pe_n = compoundSystem->calcPotentialEnergy(someState);

    logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);

    // Calculate total energy
    if(useFixman){
        etot_n = pe_n + ke_n + fix_n - (0.5 * RT * logSineSqrGamma2_n);
        etot_proposed = pe_o + ke_proposed + fix_o - (0.5 * RT * logSineSqrGamma2_o);
    }else{
        etot_n = pe_n + ke_n;
        etot_proposed = pe_o + ke_proposed;
    }
}

/** Set energies and configuration to new state **/
void HMCSampler::setSetConfigurationAndEnergiesToNew(SimTK::State& someState)
{
    setSetTVector(someState);
    pe_set = pe_n;
    fix_set = fix_n;
    logSineSqrGamma2_set = logSineSqrGamma2_n;
    ke_lastAccepted = ke_n;
    etot_set = pe_set + fix_set + ke_proposed + logSineSqrGamma2_set;
}

void HMCSampler::setSetConfigurationAndEnergiesToOld(SimTK::State& someState)
{
	assignConfFromSetTVector(someState);
	proposeExceptionCaught = false;
	acceptedStepsBuffer.push_back(0);
	acceptedStepsBuffer.pop_front();
}

SimTK::Real HMCSampler::MHAcceptProbability(SimTK::Real argEtot_proposed, SimTK::Real argEtot_n) const {
	if(argEtot_n < argEtot_proposed) {
		return 1;
	} else {
		// std::cout << "\tdiff=" << argEtot_n - argEtot_proposed << ", argEtot_n=" << argEtot_n << ", argEtot_proposed=" << argEtot_proposed << ", beta=" << beta << std::endl;
		return exp(-(argEtot_n - argEtot_proposed) * this->beta);
	}
}

/** Acception rejection step **/
bool HMCSampler::accRejStep(SimTK::State& someState) {

	bool validated;

	// Decide and get a new sample
	if ( getThermostat() == ThermostatName::ANDERSEN ) {
		// MD with Andersen thermostat
		this->acc = true;
		std::cout << "\tsample accepted (always with andersen thermostat)\n";
		update(someState);
	} else {
		// we do not consider this sample accepted unless it passes all checks
		this->acc = false;

		validated = validateProposal(); // TODO should be returned by propose()
		if (validated) {
			
			// Apply Metropolis-Hastings correction
			if(acceptSample()) {
				// sample is accepted
				this->acc = true;
				std::cout << "\tsample accepted\n";
				update(someState);
			} else {
				// sample is rejected
				std::cout << "\tsample rejected\n";
				setSetConfigurationAndEnergiesToOld(someState);
			}
			++nofSamples;
		}
		else {
			// std::cout << "\tsample not validated, reverting to previous configuration\n";
			setSetConfigurationAndEnergiesToOld(someState);
		}
	}

	return validated; // TODO should be returned by propose()
}

/** Checks if the proposal is valid **/
bool HMCSampler::validateProposal() const {
	// TODO should we check anything else?
	// TODO do we need to check the old values?

	if(proposeExceptionCaught) {
		std::cout << "\t[WARNING] invalid sample: propose exception caught!\n";
		return false;
	}

	// CHECK_IF_NAN(pe_o);
	CHECK_IF_NAN(pe_n);
	
	CHECK_IF_NAN(ke_proposed);
	CHECK_IF_NAN(ke_n);

	// CHECK_IF_NAN(fix_o);
	CHECK_IF_NAN(fix_n);

	// CHECK_IF_NAN(logSineSqrGamma2_o);
	CHECK_IF_NAN(logSineSqrGamma2_n);

	CHECK_IF_NAN(timestep);
	CHECK_IF_NAN(exp(-(etot_n - etot_proposed) / RT));

	CHECK_IF_NAN(etot_n);
	CHECK_IF_NAN(etot_proposed);

	return true;
}

/** Chooses whether to accept a sample or not based on a probability **/
bool HMCSampler::acceptSample() {
	const SimTK::Real rand_no = uniformRealDistribution(randomEngine);
	const SimTK::Real prob = MHAcceptProbability(etot_proposed, etot_n);

	// std::cout << "\trand_no=" << rand_no << ", prob=" << prob << ", beta=" << beta << std::endl;
	return rand_no < prob;
}

/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. **/
void HMCSampler::propose(SimTK::State& someState)
{

	// Store old configuration
	storeOldConfigurationAndPotentialEnergies(someState);

	// Initialize velocities according to the Maxwell-Boltzmann distribution
	initializeVelocities(someState);

	// Store the proposed energies 
	calcProposedKineticAndTotalEnergy(someState);

	// Adapt timestep
	if(shouldAdaptTimestep){
		adaptTimestep(someState);
	}
	
	// Adapt timestep
	bool shouldAdaptWorldBlocks = false;
	if(shouldAdaptWorldBlocks){
		adaptWorldBlocks(someState);
	}
	
	// Apply the L operator 
	integrateTrajectory(someState);
	//integrateTrajectoryOneStepAtATime(someState);

	calcNewConfigurationAndEnergies(someState);

	PrintDetailedEnergyInfo(someState);

}

/** Main function that contains all the 3 steps of HMC.
Implements the acception-rejection step and sets the state of the
compound to the appropriate conformation wether it accepted or not. **/
void HMCSampler::update(SimTK::State& someState)
{
	setSetConfigurationAndEnergiesToNew(someState); 
    ++acceptedSteps;
    acceptedStepsBuffer.push_back(1);
    acceptedStepsBuffer.pop_front();
}

/** Push Cartesian coordinates into R vector stored in Sampler.
Return the size of R **/
std::size_t HMCSampler::pushCoordinatesInR(SimTK::State& someState)
{
	for(const auto& topology : topologies){
		for(const auto& AtomList : topology.bAtomList){
			const auto aIx = AtomList.getCompoundAtomIndex();
			const auto& atomR = topology.calcAtomLocationInGroundFrame(someState, aIx);
			R.insert(R.end(), { atomR[0], atomR[1], atomR[2] });
		}
	}

	if(ndofs < std::numeric_limits<decltype(ndofs)>::max() / 2) {
		if(R.size() >= 2 * ndofs) {

			for(size_t j = 0; j < ndofs; j++){
				dR[j] = R[j + ndofs] - R[j];
			}

			// // Transfer upper to lower half and cleanup the upper half of R
			// for(int j = ((ndofs) - 1); j >= 0; --j){
			// 	//std::cout << "R size= " << R.size() << " j= " << j << " j+ndofs= " << j+ndofs << std::endl;
			// 	R[j] = R[j + ndofs];
			// 	R.pop_back();
			// }

			// If stuff breaks, look up for the original code. This also compacts memory
			// See https://stackoverflow.com/questions/7351899/remove-first-n-elements-from-a-stdvector for more details
			std::vector<decltype(R)::value_type>(R.begin() + ndofs, R.end()).swap(R);
		}
	} else {
		std::cout << "integer overflow at " << __LINE__ << " in " << __FILE__ << std::endl;
		throw exception();
	}

	return R.size();
}

/** Push velocities into Rdot vector stored in Sampler.
Return the size of Rdot **/
std::size_t HMCSampler::pushVelocitiesInRdot(SimTK::State& someState)
{
	for(const auto& topology : topologies){
		for(const auto& AtomList : topology.bAtomList){
			const auto aIx = AtomList.getCompoundAtomIndex();
			const auto& atomRdot = topology.calcAtomVelocityInGroundFrame(someState, aIx);
			Rdot.insert(Rdot.end(), { atomRdot[0], atomRdot[1], atomRdot[2] });
		}
	}

	// Calculate dRdot
	if(ndofs < std::numeric_limits<decltype(ndofs)>::max() / 2) {	
		if(Rdot.size() >= static_cast<size_t>(2 * ndofs)) {

			for(size_t j = 0; j < ndofs; j++){
				dRdot[j] = Rdot[j + ndofs] - Rdot[j];
			}

			// // Transfer upper to lower half and cleanup the upper half of Rdot
			// for(int j = ((ndofs) - 1); j >= 0; --j){
			// 	std::cout << "Rdotdot size= " << Rdot.size() << " j= " << j << " j+ndofs= " << j+ndofs << std::endl;
			// 	Rdot[j] = Rdot[j + ndofs];
			// 	Rdot.pop_back();
			// }

			// If stuff breaks, look up for the original code. This also compacts memory
			// See https://stackoverflow.com/questions/7351899/remove-first-n-elements-from-a-stdvector for more details
			std::vector<decltype(Rdot)::value_type>(Rdot.begin() + ndofs, Rdot.end()).swap(Rdot);
		}
	} else {
		std::cout << "integer overflow at " << __LINE__ << " in " << __FILE__ << std::endl;
		throw exception();
	}

	return Rdot.size();
}

bool HMCSampler::sample_iteration(SimTK::State& someState)
{
    std::cout << std::setprecision(10) << std::fixed;

	propose(someState); //TODO propose until validated
	if(accRejStep(someState)){

		if(this->acc) {
	
			// Add generalized coordinates to a buffer
			auto Q = someState.getQ(); // g++17 complains if this is auto& or const auto&
			QsBuffer.insert(QsBuffer.end(), Q.begin(), Q.end());
			QsBuffer.erase(QsBuffer.begin(), QsBuffer.begin() + Q.size());
	
			pushCoordinatesInR(someState);
			pushVelocitiesInRdot(someState);
	
			// Calculate MSD and RRdot to adapt the integration length
			std::cout << std::setprecision(10) << std::fixed;
			std::cout << "\tMSD= " << calculateMSD() << ", RRdot= " << calculateRRdot() << std::endl;
		}
		
		return this->acc;
	}else{
		return false;
	}
}


/** Calculate Mean Square Displacement based on stored R vectors **/
SimTK::Real HMCSampler::calculateMSD()
{
	SimTK::Real MSD = 0;
	if(dR.size() >= ndofs){
		MSD += magSq(dR);
		MSD /= static_cast<SimTK::Real>(ndofs); // TODO is this ok?
	}
	return MSD;
}

/** Calculate RRdot based on stored R and Rdot vectors **/
SimTK::Real HMCSampler::calculateRRdot()
{
	SimTK::Real RRdot = 0;
	if(dR.size() >= ndofs){
		std::vector<SimTK::Real> tempDR = dR;
		std::vector<SimTK::Real> tempRdot = Rdot;
		normalize(tempDR);
		normalize(tempRdot);

		for(std::size_t j = 0; j < ndofs; j++){
			RRdot += tempDR[j] * tempRdot[j];
		}
	}
	return RRdot;
}



int HMCSampler::getMDStepsPerSample() const {
    return MDStepsPerSample;
}

void HMCSampler::setMDStepsPerSample(int mdStepsPerSample) {
    MDStepsPerSample = mdStepsPerSample;
}

/** Print detailed energy information **/
void HMCSampler::PrintDetailedEnergyInfo(SimTK::State& someState)
{
    std::cout << std::setprecision(5) << std::fixed;
    std::cout
		<< "\tpe_o " << pe_o << ", pe_n " << pe_n << ", pe_nB " << /*" turned off for time being"*/ getPEFromEvaluator(someState)
        << "\n\tke_prop " << ke_proposed << ", ke_n " << ke_n
        << "\n\tfix_o " << fix_o << ", fix_n " << fix_n
        << "\n\tlogSineSqrGamma2_o " << logSineSqrGamma2_o << ", logSineSqrGamma2_n " << logSineSqrGamma2_n
        //<< " detmbat_n " << detmbat_n //<< " detmbat_o " << detmbat_o << " "
        << "\n\tts " << timestep  << ", exp(bdE) " << exp(-(etot_n - etot_proposed) / RT)
        << "\n\tetot_n " << etot_n  << ", etot_proposed " << etot_proposed
        << std::endl;
}


/** Modifies Q randomly
 **/
void HMCSampler::perturbQ(SimTK::State& someState)
{
    // Perturb Q
    //SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    //SimTK::Real rand_no = uniformRealDistribution_mpi_pi(randomEngine);
    int nq = someState.getNQ();
    //SimTK::Vector V(nq);
    for (int i=7; i < nq; ++i){
        //V[i] = uniformRealDistribution_mpi_pi(randomEngine);
        someState.updQ()[i] = uniformRealDistribution_mpi_pi(randomEngine);
    }
    std::cout << "perturbQ " << someState.getQ() << std::endl;
    system->realize(someState, SimTK::Stage::Position);

    // Get needed energies
    pe_o  = getOldPE();
    if(useFixman){
        fix_o = getOldFixman();
    }
    if(useFixman){
        fix_n = calcFixman(someState);
    }else{
        fix_n = 0.0;
    }

    pe_n = getPEFromEvaluator(someState); // OPENMM
    //std::cout << "Multibody PE " << getPEFromEvaluator(someState) << std::endl; // OPENMM
    //pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // ELIZA FULL

    setSetTVector(someState);
    setSetPE(pe_n);
    setSetFixman(fix_n);
    ++acceptedSteps;
    assignConfFromSetTVector(someState);

    std::cout << someState.getNU() << ' ' << 1 << ' '
              //<< getSetPE() + getREP() << ' ' << getLastAcceptedKE()
              << getSetPE() << ' ' << 0
              << ' ' << getSetFixman() << ' ' << fix_o << ' ' << fix_n << ' ';

    // Keep track of how many MC trials have been done
    ++nofSamples;
}


/** Load the map of mobods to joint types **/
/*
void HMCSampler::loadMbx2mobility(SimTK::State& someState)
{
	// Lop through topologies
	for(int i = 0; i < topologies.size(); i++){
		// Loop through atoms
		for (SimTK::Compound::AtomIndex aIx(0); aIx < topologies[i]->getNumAtoms(); ++aIx){
	
			// Get body, parentBody
			SimTK::MobilizedBodyIndex mbx = (topologies[i])->getAtomMobilizedBodyIndex(aIx);
			const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
			const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
			SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
			if(parentMbx != 0){
				// Get the neighbor atom in the parent mobilized body
				SimTK::Compound::AtomIndex chemParentAIx = (topologies[i])->getChemicalParent(matter, aIx);
			
				//std::cout << "mbx= " << mbx << " parentMbx= " << parentMbx
				//	<< " aIx= " << aIx << " chemParentAIx= " << chemParentAIx
				//	<< std::endl;
				// Get Top to parent frame
				SimTK::Compound::AtomIndex parentRootAIx = (topologies[i])->mbx2aIx[parentMbx];
			
				// Get mobility (joint type)
				bSpecificAtom *atom = (topologies[i])->updAtomByAtomIx(aIx);
				SimTK::BondMobility::Mobility mobility;
				bBond bond = (topologies[i])->getBond(
					(topologies[i])->getNumber(aIx), (topologies[i])->getNumber(chemParentAIx));
				mobility = bond.getBondMobility();
				mbx2mobility.insert(std::pair<SimTK::MobilizedBodyIndex, SimTK::BondMobility::Mobility>
                        		(mbx, mobility));
			
				(mobility == SimTK::BondMobility::Mobility::Torsion);
			
			} // Parent is not Ground
		
		} // END loop through atoms
	
	} // END loop through topologies
        for (SimTK::MobilizedBodyIndex mbx(2); mbx < matter->getNumBodies(); ++mbx){
                const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
                SimTK::QIndex qIx = mobod.getFirstQIndex(someState);
                int mobodNQ = mobod.getNumQ(someState);
                int mobodNU = mobod.getNumU(someState);
                //const SimTK::Transform T = mobod.getMobilizerTransform(someState);
                //const SimTK::MassProperties mp = mobod.getBodyMassProperties(someState);
                //const SimTK::UnitInertia unitInertia = mobod.getBodyUnitInertiaAboutBodyOrigin(someState);
                std::cout << "mbx= " << mbx << " mobility= " 
                //      << std::endl
                ;
                // Loop through body's Qs
                int internQIx = -1;
                for(int qi = qIx; qi < (mobodNQ + qIx); qi++){
                        internQIx++;
                        std::cout << " " << qi ;
                }
                std::cout << std::endl;
                // Loop through body's Us
                for(int ui = 0; ui < mobodNU; ui++){
                        //SimTK::SpatialVec H_FMCol = mobod.getH_FMCol(someState, SimTK::MobilizerUIndex(ui));
                        //std::cout << "H_FMCol= " << H_FMCol << std::endl;
                }
        }
}
*/
