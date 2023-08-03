/**@file
Implementation of MonteCarloSampler class. **/

#include "MonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor
MonteCarloSampler::MonteCarloSampler(World &argWorld,
		SimTK::CompoundSystem &argCompoundSystem,
		SimTK::SimbodyMatterSubsystem &argMatter,
		std::vector<Topology> &argTopologies, 
		SimTK::DuMMForceFieldSubsystem &argDumm,
		SimTK::GeneralForceSubsystem &argForces,
		SimTK::TimeStepper &argTimeStepper) :
        Sampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper)
{
	TVector = std::vector<SimTK::Transform>(matter->getNumBodies());
	SetTVector = std::vector<SimTK::Transform>(matter->getNumBodies());
	proposeExceptionCaught = false;

	for(int i = 0; i < acceptedStepsBufferSize; i++){ 
		acceptedStepsBuffer.push_back(0);
	}

}

// Destructor
MonteCarloSampler::~MonteCarloSampler()
{
}

// Seed the random number generator. Set simulation temperature,
// variables that store the configuration
// and variables that store the energies, both needed for the
// acception-rejection step. Also realize velocities and initialize
// the timestepper.
void MonteCarloSampler::initialize(SimTK::State& someState, SimTK::Real argTemperature, bool argUseFixman) 
{
	// Set the simulation temperature
	setTemperature(argTemperature); // Needed for Fixman


	// Store the configuration
	system->realize(someState, SimTK::Stage::Position);
	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		//const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState); // unused variable
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
	this->useFixman = argUseFixman;

	if(useFixman){
		std::cout << "Monte Carlo sampler: using Fixman potential." << std::endl;
		setOldFixman(calcFixman(someState));
		setSetFixman(getOldFixman());

		setOldLogSineSqrGamma2( ((Topology *)rootTopology)->calcLogSineSqrGamma2(someState));
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
	}else{
		setOldFixman(0.0);
		setSetFixman(getOldFixman());

		setOldLogSineSqrGamma2(0.0);
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
	}

	proposeExceptionCaught = false;
}

// Same as initialize 
void MonteCarloSampler::reinitialize(SimTK::State& someState, SimTK::Real argTemperature) 
{
	// Set the simulation temperature
	setTemperature(argTemperature); // Needed for Fixman

	// Store the configuration
	system->realize(someState, SimTK::Stage::Position);
	int i = 0;
	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		// const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState); // unused variable
		SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
		i++;
	}

	// Store potential energies
	//setOldPE(getPEFromEvaluator(someState));
	setSetPE(getOldPE());

	// Store Fixman potential
	if(useFixman){
		setOldFixman(calcFixman(someState));
		setSetFixman(getOldFixman());

		setOldLogSineSqrGamma2( ((Topology *)rootTopology)->calcLogSineSqrGamma2(someState));
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
	}else{
		setOldFixman(0.0);
		setSetFixman(getOldFixman());

		setOldLogSineSqrGamma2(0.0);
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
	}

	proposeExceptionCaught = false;
}

void MonteCarloSampler::useFixmanPotential()
{
    useFixman = true;
}


// Return true if use Fixman potential
bool MonteCarloSampler::isUsingFixmanPotential() const
{
    return useFixman;
}

// Is the sampler always accepting the proposed moves
//bool MonteCarloSampler::getAlwaysAccept(void) const
//{
//   return alwaysAccept;
//}

// Is the sampler always accepting the proposed moves
//void MonteCarloSampler::setAlwaysAccept(bool argAlwaysAccept)
//{
//    alwaysAccept = argAlwaysAccept;
//}


// Compute Fixman potential (should have been calcDetMInv ??)
SimTK::Real MonteCarloSampler::calcFixman(SimTK::State& someState){
    int nu = someState.getNU();
    SimTK::Vector V(nu);

    //for (int i=0; i < nu; ++i){
    //   V[i] = i;
    //}
    system->realize(someState, SimTK::Stage::Position);
    matter->realizeArticulatedBodyInertias(someState); // Move in calcDetM ?

    // Get M
    //SimTK::Matrix M(nu, nu);
    //matter->calcM(someState, M);

    // Get detM
    SimTK::Vector DetV(nu);
    SimTK::Real D0 = 1.0;

    // TODO: remove the request for Dynamics stage cache in SImbody files
    //std::cout << "MonteCarloSampler::calcFixman Stage: "<< matter->getStage(someState) << std::endl;
    matter->calcDetM(someState, V, DetV, &D0);

    //std::cout << "FixmanTorque: " << "MonteCarloSampler::calcFixman logdetM: " << std::setprecision(10) << std::log(D0) << std::setprecision(2) << std::endl;
    //std::cout << "MonteCarloSampler::calcFixman RT: " << RT << std::endl;
    // ---- Verify with Eigen ----------
    // Eigen M determinant
    //Eigen::MatrixXd EiM(nu, nu);
    //for(int i=0; i<nu; i++){
    //    for(int j=0; j<nu; j++){
    //        EiM(i, j) = M(i, j);
    //    }
    //}
    //SimTK::Real EiDetM = EiM.determinant();
    //std::cout << "EiDetM= " << EiDetM << std::endl;
    assert(RT > SimTK::TinyReal);
    //SimTK::Real result = 0.5 * RT * (((Topology *)rootTopology)->calcLogDetMBAT(someState) - std::log(D0));
    SimTK::Real result = 0.5 * RT * ( std::log(D0) - ((Topology *)rootTopology)->calcLogDetMBATInternal(someState) );
    //SimTK::Real result = 0.5 * RT * std::log(D0);
    //std::cout << std::setprecision(5) << std::fixed << "MonteCarloSampler::calcFixman 05kT logdetM: " << 0.5 * RT * std::log(D0) << std::endl;
    //std::cout << std::setprecision(5) << std::fixed << "MonteCarloSampler::calcFixman -05kT logdetMBATi: " << -0.5 * RT * ((Topology *)rootTopology)->calcLogDetMBATInternal(someState) << std::endl;

    if(SimTK::isInf(result)){
        result = 0.0;
    }
    
    return result;
}

// Compute Fixman potential numerically
SimTK::Real MonteCarloSampler::calcNumFixman(SimTK::State&){
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

// Compute mass matrix determinant numerically
// Not to be confused with the Fixman potential
SimTK::Real MonteCarloSampler::calcNumDetM(SimTK::State&){
    // args were SimTK::State& someState
    assert(!"Not implemented.");
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
    return EiM.determinant();
    */
    return 0;
}

// Get the set potential energy
SimTK::Real MonteCarloSampler::getSetPE() const
{
    return this->pe_set;
}

// Get the stored potential energy
SimTK::Real MonteCarloSampler::getOldPE() const
{
    return this->pe_o;
}

// Set set potential energy
void MonteCarloSampler::setSetPE(SimTK::Real argPE)
{
    this->pe_set = argPE;
}

// Set stored potential energy
void MonteCarloSampler::setOldPE(SimTK::Real argPE)
{
    this->pe_o = argPE;
}

// Set set Fixman potential
void MonteCarloSampler::setSetFixman(SimTK::Real argFixman)
{
    this->fix_set = argFixman;
}

// Set old Fixman potential
void MonteCarloSampler::setOldFixman(SimTK::Real argFixman)
{
    this->fix_o = argFixman;
}

// Set old Fixman potential
void MonteCarloSampler::setProposedFixman(SimTK::Real argFixman)
{
    this->fix_n = argFixman;
}

// Get set Fixman potential
SimTK::Real MonteCarloSampler::getSetFixman() const
{
    return this->fix_set;
}

// Get Fixman potential
SimTK::Real MonteCarloSampler::getOldFixman() const
{
    return this->fix_o;
}

// Get Fixman potential
SimTK::Real MonteCarloSampler::getProposedFixman() const
{
    return this->fix_n;
}




// Set/get External MBAT contribution potential
void MonteCarloSampler::setSetLogSineSqrGamma2(SimTK::Real argX){
    this->logSineSqrGamma2_set = argX;
}

SimTK::Real MonteCarloSampler::getSetLogSineSqrGamma2() const
{
    return this->logSineSqrGamma2_set;
}

// Set/get External MBAT contribution potential
void MonteCarloSampler::setOldLogSineSqrGamma2(SimTK::Real argX){
    this->logSineSqrGamma2_o = argX;
}

SimTK::Real MonteCarloSampler::getOldLogSineSqrGamma2() const
{
    return this->logSineSqrGamma2_o;
}

void MonteCarloSampler::setProposedLogSineSqrGamma2(SimTK::Real argFixman)
{
    this->logSineSqrGamma2_n = argFixman;
}

SimTK::Real MonteCarloSampler::getProposedLogSineSqrGamma2() const
{
    return this->logSineSqrGamma2_n;
}



// Set/get Residual Embedded Potential
void MonteCarloSampler::setREP(SimTK::Real inp)
{
    this->residualEmbeddedPotential = inp;
}

SimTK::Real MonteCarloSampler::getREP() const
{
    return this->residualEmbeddedPotential;
}

// Stores the configuration into an internal vector of transforms TVector
void MonteCarloSampler::setTVector(const SimTK::State& someState)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    TVector[i] = mobod.getMobilizerTransform(someState);
    i++;
  }
}

// Stores the configuration into an internal vector of transforms TVector
void MonteCarloSampler::setTVector(SimTK::Transform *inpTVector)
{
    // TODO pointer parameter is bad
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    TVector[i] = inpTVector[i];
    i++;
  }
}

// Stores the set configuration into an internal vector of transforms TVector
void MonteCarloSampler::setSetTVector(const SimTK::State& someState)
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
SimTK::Transform * MonteCarloSampler::getTVector()
{
    return &TVector[0];
}

// Stores the configuration into an internal vector of transforms TVector
// Get the stored configuration
SimTK::Transform * MonteCarloSampler::getSetTVector()
{
    return &SetTVector[0];
}

// Restores configuration from the internal set vector of transforms TVector
void MonteCarloSampler::assignConfFromSetTVector(SimTK::State& someState)
{
  //someState.invalidateAll(SimTK::Stage::Instance);
  //matter->invalidateArticulatedBodyInertias(someState);
  //system->realize(someState, SimTK::Stage::Position);
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    system->realize(someState, SimTK::Stage::Position);
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    mobod.setQToFitTransform(someState, SetTVector[i]);

    system->realize(someState, SimTK::Stage::Position);
    //matter->realizeArticulatedBodyInertias(someState);

    //std::cout <<  "Sampler after setQ State Cache Info: " << std::endl;
    //PrintSimbodyStateCache(someState);

    i++;
  }
  system->realize(someState, SimTK::Stage::Position);
  //matter->realizeArticulatedBodyInertias(someState);
}

// Restores configuration from the internal vector of transforms TVector
void MonteCarloSampler::assignConfFromTVector(SimTK::State& someState)
{
  //someState.invalidateAll(SimTK::Stage::Instance);
  //matter->invalidateArticulatedBodyInertias(someState);
  //system->realize(someState, SimTK::Stage::Position);
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    system->realize(someState, SimTK::Stage::Position);
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    mobod.setQToFitTransform(someState, TVector[i]);

    system->realize(someState, SimTK::Stage::Position);
    //matter->realizeArticulatedBodyInertias(someState);
    i++;
  }
  system->realize(someState, SimTK::Stage::Position);
  //matter->realizeArticulatedBodyInertias(someState);
}

// Assign random conformation
// Assign random conformation
// In torsional dynamics the first body has 7 Q variables for 6 dofs - one
// quaternion (q) and 3 Cartesian coordinates (x). updQ will return: 
// [qw, qx, qy, qz, x1, x2, x3]
bool MonteCarloSampler::proposeEquilibrium(SimTK::State& someState)
{
    for (int i = 1; i < matter->getNumBodies(); ++i){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(i));
        const SimTK::Real rand_no = uniformRealDistribution_0_2pi(randomEngine);

        for(int j=0; j<mobod.getNumQ(someState); j++){
            mobod.setOneQ(someState, j, rand_no);
            //someState.updQ()[i] = rand_no;
        }
    }

    system->realize(someState, SimTK::Stage::Position); // NECESSARY

    // TODO
    assert(!"What should we return here?");
    return true;

}

// The update step in Monte Carlo methods consists in:
// Acception - rejection step
void MonteCarloSampler::update(SimTK::State& someState){
    // SimTK_BOLTZMANN_CONSTANT_MD is long double, but it shouldn't be a problem. Should we add a higher precision constant?
    // mind that this is not the only occurence in our code base
    RT = getTemperature() * static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);

    // Get old energy
    pe_o = getOldPE();

    // Assign random configuration
    proposeEquilibrium(someState);

    // Send configuration to evaluator  
    sendConfToEvaluator(); // OPENMM

    // Get current potential energy from evaluator
    pe_n = getPEFromEvaluator(someState); // OPENMM

    // Apply Metropolis criterion
    assert(!std::isnan(pe_n));
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    if ((proposeExceptionCaught == false) || (pe_n < pe_o) || (rand_no < exp(-(pe_n - pe_o)/RT))){
        // Accept
        setTVector(someState);
        setOldPE(pe_n);
        ++acceptedSteps;
        this->acc = true;
    }else{ // Reject
        assignConfFromTVector(someState);
        this->acc = false;
    }
}

// Get the potential energy from an external source as far as the sampler
// is concerned - OPENMM has to be inserted here
SimTK::Real MonteCarloSampler::getPEFromEvaluator(SimTK::State& someState){
    if ( getThermostat() == ThermostatName::ANDERSEN ){
        return forces->getMultibodySystem().calcPotentialEnergy(someState);
        //return dumm->CalcFullPotEnergyIncludingRigidBodies(someState);
    }else{
        return forces->getMultibodySystem().calcPotentialEnergy(someState);
        //return dumm->CalcFullPotEnergyIncludingRigidBodies(someState);
    }
}

/*
// Get the desired simulation temperature. Not to be confused with 
// the instant temperature
SimTK::Real MonteCarloSampler::getTemperature(){
    return this->temperature;
}

// Set the desired simulation temperature. Not to be confused with 
// the instant temperature
void MonteCarloSampler::setTemperature(SimTK::Real argTemperature)
{
    std::cout << " MonteCarloSampler::setTemperature " << argTemperature << std::endl;
    this->temperature = argTemperature;
    if(this->temperature < 0){
        std::cerr << "Temperature set to " << this->temperature << std::endl;
    }
    RT = this->temperature * SimTK_BOLTZMANN_CONSTANT_MD;
}
*/

// Set a thermostat
void MonteCarloSampler::setThermostat(ThermostatName argThermostat){
    this->thermostat = argThermostat;
}
// Set a thermostat
void MonteCarloSampler::setThermostat(std::string thermoName){
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
void MonteCarloSampler::setThermostat(const char *argThermostat){
    setThermostat(std::string(argThermostat));
}

// Get the name of the thermostat
ThermostatName MonteCarloSampler::getThermostat() const{
    return thermostat;
}




// Send configuration to an external evaluator

void MonteCarloSampler::sendConfToEvaluator(){
    assert(!"Not implemented");
}

// Get the number of accpted conformations
int MonteCarloSampler::getAcceptedSteps() const
{
    return acceptedSteps;
}


/** Modifies Q randomly
 **/
void MonteCarloSampler::perturbQ(SimTK::State& someState)
{
	// Perturb Q
	//SimTK::Real rand_no = uniformRealDistribution(randomEngine);
	SimTK::Real rand_no = uniformRealDistribution_mpi_pi(randomEngine);
	int nq = someState.getNQ();
	//SimTK::Vector V(nq);

	//
	SimTK::Real q0_refVal, q1_refVal;
	if(rand_no < 0){
		q0_refVal = 2;
		q1_refVal = -1.7;
	}else{
		q0_refVal = 0;
		q1_refVal = 0;
	}
	// q1_refVal = -2.3
	SimTK::Real q0_actVal = someState.getQ()[0];
	SimTK::Real q1_actVal = someState.getQ()[1];

	// Draw from vonMises distribution
	SimTK::Real q0_vonMises = generateVonMisesSample(q0_actVal, 5);
	if(std::isnan(q0_vonMises)){q0_vonMises = q0_actVal;}
	std::cout << "vonMises " << q0_vonMises << "\n";

	// Add Gaussian noise
	//std::normal_distribution<double> Noiser(q0_refVal, 1);
	//double noise = Noiser(RandomCache.RandomEngine);

	std::cout << "HMCSampler::perturbQ " << someState.getQ() << " with " << q0_vonMises << std::endl;
	//someState.updQ()[which] = uniformRealDistribution_mpi_pi(randomEngine);
	someState.updQ()[0] = q0_vonMises;
	//someState.updQ()[1] = q1_vonMises;
	system->realize(someState, SimTK::Stage::Position);
	std::cout << "HMCSampler::perturbQ " << someState.getQ() << std::endl;

}


