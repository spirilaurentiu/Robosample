/* -------------------------------------------------------------------------- *
 *                       Robosample: Gmolmodel                                *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Copyright information.                                                     *
 * Authors: Laurentiu Spiridon                                                *
 * Contributors: Eliza C. Martin, Victor Ungureanu, Teodor A. Sulea           *
 *                                                                            *
 * Licensed under ; you may                                                   *
 * not use this file except                                                   *
 * You may obtain a copy of the License at                                    *
 * -------------------------------------------------------------------------- */


/** @file
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
//==============================================================================
//                           CONSTRUCTOR
//==============================================================================
// Description.
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
	// Ensure Sampler prerequisites
	assert(argCompoundSystem != nullptr);
	assert(argMatter != nullptr);
	assert(argDumm != nullptr);
	assert(argForces != nullptr);
	assert(argTimeStepper != nullptr);

	this->system = &argMatter->getSystem();


	//this->rootTopology = argResidue;

	if( !(topologies.size() > 0) ){
		std::cerr << "HMCSampler: No topologies found. Exiting...";
		throw std::exception();
		std::exit(1);
	}

	this->rootTopology = &topologies[0];

	// Set total number of atoms and dofs
	natoms = 0;
	for (const auto& topology: topologies){
		natoms += topology.getNumAtoms();
	}

	// This is just an assumption; We don't know if topology is realized
	// TODO: BUG: ndofs may not be equal with NU
	int ThreeFrom3D = 3;
	ndofs = natoms * ThreeFrom3D;
	// Sampler end

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
	//this->alwaysAccept = false;

	this->prevPrevAcceptance = SimTK::NaN;
	this->prevAcceptance = SimTK::NaN;
	this->acceptance = SimTK::NaN;


	this->prevPrevTimestep = SimTK::NaN; // ps
	this->prevTimestep = SimTK::NaN; // ps
	this->timestep = 0;

	this->temperature = 0.0;
	this->boostT = 0.0;
	this->boostRT = 0.0;
	this->sqrtBoostRT = 0.0;
	this->boostBeta = 0.0;

	MDStepsPerSample = 0;
	proposeExceptionCaught = false;
	shouldAdaptTimestep = false;

	dR.resize(ndofs, 0);
	dRdot.resize(ndofs, 0);

	// Put ridiculous values so it won't work by accident
	UScaleFactors.resize(ndofs, 99999999);
	UScaleFactorsNorm = std::sqrt(ndofs);
	InvUScaleFactors.resize(ndofs, 99999999);
	InvUScaleFactorsNorm = std::sqrt(ndofs);

	NormedUScaleFactors.resize(ndofs);
	DOFUScaleFactors.resize(ndofs);

	//NMARotation.resize(ndofs, std::vector<double>(ndofs, 0));
	//for(int i = 0; i < ndofs; i++){
	//	NMARotation[i][i] = 1.0;
	//}

	sampleGenerator = 0;
	integratorName = IntegratorName::EMPTY;

	// Non-equilibrium options
	DistortOpt = 0;
	FlowOpt = 0;
	WorkOpt = 0;

	// None-equilibrium params
	QScaleFactor = 1.0;

	MDStepsPerSampleStd = 0.5;

	NMAAltSign = 1.0;

	// Set total mass of the system to non-realistic value
	totalMass = 0;

	// Work
	work = 0.0;

}

/** Destructor **/
HMCSampler::~HMCSampler()
{
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
	//setOldPE(getPEFromEvaluator(someState));
	setOldPE(forces->getMultibodySystem().calcPotentialEnergy(someState));
	//setOldPE(dumm->CalcFullPotEnergyIncludingRigidBodies(someState));
	setSetPE(getOldPE());

	// Store Fixman potential
	if(useFixman){
		std::cout << "Hamiltonian Monte Carlo sampler: using Fixman potential.\n";
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

	for (int j=0; j < nu; ++j){
	UScaleFactors[j] = 1;
	InvUScaleFactors[j] = 1;
	}

	loadUScaleFactors(someState);

	// Buffer to hold Q means
	QsMeans.resize(nq), 0;

	// Total mass of the system
	this->totalMass = matter->calcSystemMass(someState);

	// Work
	work = 0.0;

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

		setOldLogSineSqrGamma2( ((Topology *)rootTopology)->calcLogSineSqrGamma2(someState));
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

	// Reset ndofs which was set to natoms*3 in constructors
	ndofs = nu;

	// kT
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

	for (int j=0; j < nu; ++j){
        UScaleFactors[j] = 1;
        InvUScaleFactors[j] = 1;
	}
	loadUScaleFactors(someState);

	// Total mass of the system
	this->totalMass = matter->calcSystemMass(someState);
	
	// Work
	work = 0.0;

}


SimTK::Real HMCSampler::getMDStepsPerSampleStd() const
{
	return MDStepsPerSampleStd;
}

void HMCSampler::setMDStepsPerSampleStd(SimTK::Real mdstd){
	MDStepsPerSampleStd = mdstd;
}

// Set the method of integration
void HMCSampler::setSampleGenerator(std::string& generatorNameNameArg)
{
	if (generatorNameNameArg == "EMPTY"){
		this->sampleGenerator = 0;
		setAlwaysAccept(true);
	}
	else if(generatorNameNameArg == "MC"){
		this->sampleGenerator = 1;
		setAlwaysAccept(false);
	}else{
		std::cerr << "Unknown sampling method.\n"; throw std::exception(); std::exit(1);
	}
}

void HMCSampler::setIntegratorName(IntegratorName integratorNameArg)
{
	this->integratorName = integratorNameArg;
}

void HMCSampler::setIntegratorName(std::string integratorNameArg)
{

	//this->integratorName = IntegratorNameS[integratorNameArg];

 	if(integratorNameArg == "OMMVV"){
		this->integratorName = IntegratorName::OMMVV;
	}else if (integratorNameArg == "VV"){
		integratorName = IntegratorName::VERLET;
	}else{
		integratorName = IntegratorName::EMPTY;
	}

}

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

		// OLD RESTORE TODO BOOST
		//RandomCache.sqrtRT = std::sqrt(RT);
		// NEW TRY TODO BOOST
		RandomCache.sqrtRT = std::sqrt(this->boostRT);

		RandomCache.nu = nu;

		// we don't get to use multithreading here
		RandomCache.FillWithGaussian();
	} else {
		// wait for random number generation to finish (should be done by this stage)
		RandomCache.task.wait();
	}

	// V[i] *= UScaleFactors[i] - note that V is already populated with random numbers
	//std::transform(RandomCache.V.begin(), RandomCache.V.end(), // apply an operation on this
	//	UScaleFactors.begin(), // and this
	//	RandomCache.V.begin(), // and store here
	//	std::multiplies<SimTK::Real>()); // this is the operation

	// Scale by square root of the inverse mass matrix
	matter->multiplyBySqrtMInv(someState, RandomCache.V, RandomCache.SqrtMInvV);

	// Set stddev according to temperature
	RandomCache.SqrtMInvV *= sqrtBoostRT;

	// Raise the temperature
	someState.updU() = RandomCache.SqrtMInvV;

	// Realize velocity
	system->realize(someState, SimTK::Stage::Velocity);

	// ask for a number of random numbers and check if we are done the next
	//  time we hit this function
	RandomCache.task = std::async(std::launch::async, RandomCache.FillWithGaussian);
}


/** Initialize velocities scaled by NMA factors.
 Coresponds to R operator in LAHMC **/
void HMCSampler::initializeNMAVelocities(SimTK::State& someState){

	// Get the number of dofs
	const int nu = someState.getNU();

	// Get sqrt of kT
	double sqrtRT = std::sqrt(RT);

	// U vector
	SimTK::Vector Us(nu, 1);

	// Vector multiplied by the sqrt of the inverse mass matrix
	SimTK::Vector sqrtMInvUs(nu, 1);

	if(DistortOpt == 1){ // For pure demonstrative purposes

		if((nofSamples % 2) == 0){ NMAAltSign = 1;}else{NMAAltSign = -1;}
		//if(((nofSamples) % 6) == 0){ NMAAltSign *= -1;}
		//std::cout << "NMAAltSign " << NMAAltSign << std::endl;

		Us *= NMAAltSign;

		// Multiply by the square root of the inverse mass matrix
		matter->multiplyBySqrtMInv(someState, Us, sqrtMInvUs);

		// Scale by NMA scale factors
		std::transform(Us.begin(), Us.end(), // apply an operation on this
			UScaleFactors.begin(), // and this
			Us.begin(), // and store here
			std::multiplies<SimTK::Real>()); // this is the operation

		// Restore vector length
		for(int i = 0; i < nu; i++){
			Us[i] = (Us[i] * (std::sqrt(nu) / UScaleFactorsNorm));
		}

		// Multiply by the square root of kT
		sqrtMInvUs *= sqrtRT;

	}else if(DistortOpt == 2){ // Use UscaleFactors with random sign

		// Get random -1 or 1
		SimTK::Real randSign;
		SimTK::Real randUni_m1_1 = uniformRealDistribution_m1_1(randomEngine);
		randSign = (randUni_m1_1 > 0) ? 1 : -1 ;
		Us *= randSign;

		//SimTK::Real totalMass = matter->calcSystemMass(someState);// DELETE

		///* Multiply by H sqrt(Mcartesian) ????

		for (SimTK::MobilizedBodyIndex mbx(2); mbx < matter->getNumBodies(); ++mbx){
			const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
			SimTK::Real M = mobod.getBodyMass(someState);
			SimTK::Real sqrtMobodMInv = 1 / std::sqrt(M);

			for(SimTK::UIndex uIx = mobod.getFirstUIndex(someState);
				uIx < mobod.getFirstUIndex(someState) + mobod.getNumU(someState);
				uIx++ ){
				sqrtMInvUs[int(uIx)] = Us[int(uIx)] * sqrtMobodMInv;
			}
		}

		std::cout << "sqrtMInvUs 1 "; for(int i = 0; i < nu; i++){ std::cout << sqrtMInvUs[i] << " " ;} std::cout << "\n";

		for(int i = 0; i < nu; i++){
			sqrtMInvUs[i] = (sqrtMInvUs[i] * 2.64);
		}

		// Scale by NMA scale factors
		std::transform(Us.begin(), Us.end(), // apply an operation on this
			UScaleFactors.begin(), // and this
			Us.begin(), // and store here
			std::multiplies<SimTK::Real>()); // this is the operation

		// Restore vector length
		for(int i = 0; i < nu; i++){
			Us[i] = (Us[i] * (std::sqrt(nu) / UScaleFactorsNorm));
		}

		// Multiply by the square root of kT
		sqrtMInvUs *= sqrtRT;

		//*/
	}else if(DistortOpt == 3){ // Set a random sign

		// Check if we can use our cache
		// Draw X from a normal distribution N(0, 1)
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

		Us = RandomCache.V;

		// Multiply by the square root of the inverse mass matrix
		matter->multiplyBySqrtMInv(someState, Us, sqrtMInvUs);

		// Scale by NMA scale factors
		std::transform(Us.begin(), Us.end(), // apply an operation on this
			UScaleFactors.begin(), // and this
			Us.begin(), // and store here
			std::multiplies<SimTK::Real>()); // this is the operation

		// Restore vector length
		for(int i = 0; i < nu; i++){
			Us[i] = (Us[i] * (std::sqrt(nu) / UScaleFactorsNorm));
		}

		// Multiply by the square root of kT
		sqrtMInvUs *= sqrtRT;

	}else if(DistortOpt == 4){ // Set a random sign

		// Direction vector UScaleFactors
		//SimTK::Vector NormedUScaleFactors(nu, 1);
		//SimTK::Vector DOFUScaleFactors(nu, 1);
		//for(int i = 0; i < nu; i++){
		//	NormedUScaleFactors[i] = UScaleFactors[i] / UScaleFactorsNorm;
		//}
		//for(int i = 0; i < nu; i++){
		//	DOFUScaleFactors[i] = NormedUScaleFactors[i] * std::sqrt(nu);
		//}

		//std::cout << "DOFUScaleFactors 0 "; for(int i = 0; i < nu; i++){ std::cout << DOFUScaleFactors[i] << " " ;} std::cout << "\n";

		// Check if we can use our cache
		// Draw X from a normal distribution N(0, 1)
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

		Us = RandomCache.V;

		std::cout << "Us 1 "; for(int i = 0; i < nu; i++){ std::cout << Us[i] << " " ;} std::cout << "\n";

		// Add noise
		std::transform(DOFUScaleFactors.begin(), DOFUScaleFactors.end(), // apply an operation on this
			Us.begin(), // and this
			DOFUScaleFactors.begin(), // and store here
			std::plus<SimTK::Real>()); // this is the operation

		std::cout << "DOFUScaleFactors 1 "; for(int i = 0; i < nu; i++){ std::cout << DOFUScaleFactors[i] << " " ;} std::cout << "\n";

		matter->multiplyBySqrtMInv(someState, DOFUScaleFactors, sqrtMInvUs);

		std::cout << "sqrtMInvUs 2 "; for(int i = 0; i < nu; i++){ std::cout << sqrtMInvUs[i] << " " ;} std::cout << "\n";

		// Scale by sqrt of kT
		sqrtMInvUs *= sqrtRT;

		std::cout << "sqrtMInvUs 3 "; for(int i = 0; i < nu; i++){ std::cout << sqrtMInvUs[i] << " " ;} std::cout << "\n";

		// Scale by NMA scale factors
		//std::transform(sqrtMInvUs.begin(), sqrtMInvUs.end(), // apply an operation on this
		//	NormedUScaleFactors.begin(), // and this
		//	sqrtMInvUs.begin(), // and store here
		//	std::plus<SimTK::Real>()); // this is the operation


		// Restore vector length
		//for(int i = 0; i < nu; i++){
		//	sqrtMInvUs[i] = (sqrtMInvUs[i] * (std::sqrt(nu) / UScaleFactorsNorm));
		//}



	}else if(DistortOpt == 5){
	    std::cout << "NMAOPTION 5 IS STILL UNDER DEVELOPMENT !!\n" <<std::flush;

		//proj(U, V, W);

		// Generate unit Von Mises Fisher around [1, 0, 0...]
		vector<double> X;
		vector<double> U;
		X.resize(ndofs, 0.0);
		U.resize(ndofs, 0.0);
		double concentration = 100;

		vonMisesFisher(X, concentration);
		std::cout << "X 0 "; for(int i = 0; i < nu; i++){ std::cout << X[i] << " " ;} std::cout << "\n";

		// Rotate to appropiate location
		bMulVecByMatrix(X, NMARotation, U);
		std::cout << "U 1 "; for(int i = 0; i < nu; i++){ std::cout << U[i] << " " ;} std::cout << "\n";
		//bPrintMat(NMARotation);

		// Get random -1 or 1 and generate bidirectional von Mises Fisher
		SimTK::Real randSign;
		SimTK::Real randUni_m1_1 = uniformRealDistribution_m1_1(randomEngine);
		randSign = (randUni_m1_1 > 0) ? 1 : -1 ;
		Us *= randSign;

		// Put c++ vMF sample in Simbody Vector
		for(int j = 0; j < ndofs; j++){Us[j] *= U[j];}
		std::cout << "Us 2 "; for(int i = 0; i < nu; i++){ std::cout << Us[i] << " " ;} std::cout << "\n";

		// Prolong vector to match temperature
		SimTK::Real velocityMacroLength = (std::sqrt(nu) * sqrtRT);
		for(int j = 0; j < ndofs; j++){Us[j] *= velocityMacroLength;}
		std::cout << "Us 3 "; for(int i = 0; i < nu; i++){ std::cout << Us[i] << " " ;} std::cout << "\n";

		// Add Gaussian noise to the length
		std::normal_distribution<double> lengthNoiser(1, 0.1);
		double lengthNoise = lengthNoiser(RandomCache.RandomEngine);
		std::cout << "lengthNoise " << lengthNoise << "\n";
		for(int j = 0; j < ndofs; j++){Us[j] *= lengthNoise;}
		std::cout << "Us 4 "; for(int i = 0; i < nu; i++){ std::cout << Us[i] << " " ;} std::cout << "\n";

		// Multiply by the square root of the mass matrix
		matter->multiplyBySqrtMInv(someState, Us, sqrtMInvUs);
		std::cout << "sqrtMInvUs 5 "; for(int i = 0; i < nu; i++){ std::cout << sqrtMInvUs[i] << " " ;} std::cout << "\n";

	}else if(DistortOpt == 6){

		//std::cout << "DOFUs 0 "; for(int i = 0; i < nu; i++){ std::cout << DOFUScaleFactors[i] << " " ;} std::cout << "\n";

		// Sample from N(0, 1)
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

		Us = RandomCache.V;
		//std::cout << "Us 0 "; for(int i = 0; i < nu; i++){ std::cout << Us[i] << " " ;} std::cout << "\n";

		// Get random -1 or 1
		SimTK::Real randSign;
		SimTK::Real randUni_m1_1 = uniformRealDistribution_m1_1(randomEngine);
		randSign = (randUni_m1_1 > 0) ? 1 : -1 ;

		// Shift average
		for(int j = 0; j < ndofs; j++){Us[j] += ((randSign) * DOFUScaleFactors[j]);}
		//std::cout << "Us 1 "; for(int i = 0; i < nu; i++){ std::cout << Us[i] << " " ;} std::cout << "\n";

		// Multiply by the square root of the mass matrix
		matter->multiplyBySqrtMInv(someState, Us, sqrtMInvUs);
		//std::cout << "sqrtMInvUs 1 "; for(int i = 0; i < nu; i++){ std::cout << sqrtMInvUs[i] << " " ;} std::cout << "\n";

		// Scale by sqrt of guidance kT
		//SimTK::Real sqrtBoostRT = std::sqrt(boostRT); // moved into private vars
		sqrtMInvUs *= sqrtBoostRT;
		//sqrtMInvUs *= sqrtRT;
		//std::cout << "sqrtRT " << sqrtRT
		//	<< " sqrtBoostRT " << sqrtBoostRT
		//	<< std::endl;
		//std::cout << "sqrtMInvUs 2 "; for(int i = 0; i < nu; i++){ std::cout << sqrtMInvUs[i] << " " ;} std::cout << "\n";

		//////////////////////
		// Regulate temperature of the guidance Hamiltonian
		sqrtMInvUs /= 1.4142135623730950488; // sqrt(2)
		//std::cout << "sqrtMInvUs 3 "; for(int i = 0; i < nu; i++){ std::cout << sqrtMInvUs[i] << " " ;} std::cout << "\n";
		//////////////////////

		// Probability terms
		SimTK::Real boostBetaPlus = this->boostBeta * 1.4142135623730950488;
		SimTK::Real localBoostFactor = boostBetaPlus / this->beta;
		//std::cout << "b+sqrt(2)/b " << localBoostFactor << std::endl;

		// Kinetic term
		SimTK::Vector X(sqrtMInvUs);
		SimTK::Vector XM(nu, 1);
		SimTK::Real XMX = 0.0;
		SimTK::Real bXMX = 0.0;

		matter->multiplyByM(someState, X, XM);
		for(int j = 0; j < ndofs; j++){XMX += (XM[j] * X[j]);}
		bXMX = localBoostFactor * XMX;

		// Bias term
		//SimTK::Vector DOFUScaleFactorsM(nu, 1);
		//SimTK::Real muMmu = 0.0;
		//SimTK::Real bmuMmu = 0.0;
		//matter->multiplyByM(someState, DOFUScaleFactors, DOFUScaleFactorsM);
		//for(int j = 0; j < ndofs; j++){muMmu += (DOFUScaleFactorsM[j] * DOFUScaleFactors[j]);}
		//bmuMmu = localBoostFactor * muMmu;

		// Correlation term
		SimTK::Real muMX = 0.0;
		SimTK::Real bmuMX = 0.0;
		SimTK::Real _bmuMX = 0.0;
		SimTK::Real bcorr = 0.0;
		for(int j = 0; j < ndofs; j++){muMX += (XM[j] * DOFUScaleFactors[j]);}
		bmuMX = boostBetaPlus * muMX;
		_bmuMX = -1.0 * bmuMX;

		if(bcorr < 16){ // safe limit to use exponentials
			bcorr = LSE2(bmuMX, _bmuMX);
		}else{
			bcorr = (bmuMX >= 0 ? bmuMX : _bmuMX); // approximation to LSE
		}
		bcorr /= this->beta;

		std::cout << "Keval terms bXMX muMmu bmuMX LSE(bmuMX) = "
			<< 0.5*bXMX << " " // << 0.5*bmuMmu << " "
			<< bmuMX << " " << bcorr << std::endl;

		ke_prop_nma6 = // (-1.0 * RT * std::log(0.5))
			+ (0.5*bXMX) // + (0.5* bmuMmu)
			- bcorr;
		// */

		/*// Kinetic energy proposed for the evaluation Hamiltonian
		ke_prop_nma6 = 0;
		SimTK::Vector v_miu(nu, 1);
		SimTK::Vector v_miuM(nu, 1);
		SimTK::Real leftExp = 0.0;
		SimTK::Real rightExp = 0.0;


		// Left exponential
		v_miu = (sqrtMInvUs - DOFUScaleFactors);
		std::cout << "v_miu 0 "; for(int i = 0; i < nu; i++){ std::cout << v_miu[i] << " " ;} std::cout << "\n";

		matter->multiplyByM(someState, v_miu, v_miuM);
		std::cout << "v_miuM 1 "; for(int i = 0; i < nu; i++){ std::cout << v_miuM[i] << " " ;} std::cout << "\n";

		for(int j = 0; j < ndofs; j++){leftExp += (v_miuM[j] * v_miu[j]);}
		std::cout << "leftExp 0 " << leftExp << "\n";

		//leftExp *= ((-0.5) * 1.4142135623730950488 * this->beta);
		leftExp *= ((-0.5) * 1.4142135623730950488 * this->boostBeta);
		std::cout << "leftExp 1 " << leftExp << "\n";

		leftExp = std::exp(leftExp);
		std::cout << "leftExp 2 " << leftExp << "\n";

		// Right exponential
		v_miu = (sqrtMInvUs + DOFUScaleFactors);
		std::cout << "v_miu 0 "; for(int i = 0; i < nu; i++){ std::cout << v_miu[i] << " " ;} std::cout << "\n";

		matter->multiplyByM(someState, v_miu, v_miuM);
		std::cout << "v_miuM 1 "; for(int i = 0; i < nu; i++){ std::cout << v_miuM[i] << " " ;} std::cout << "\n";

		for(int j = 0; j < ndofs; j++){rightExp += (v_miuM[j] * v_miu[j]);}
		std::cout << "rightExp 0 " << rightExp << "\n";

		//rightExp *= ((-0.5) * 1.4142135623730950488 * this->beta);
		std::cout << "beta boostBeta " << this->beta << " " << this->boostBeta << "\n" ;
		rightExp *= ((-0.5) * 1.4142135623730950488 * this->boostBeta);
		std::cout << "rightExp 1 " << rightExp << "\n";

		rightExp = std::exp(rightExp);
		std::cout << "rightExp 2 " << rightExp << "\n";

		ke_prop_nma6 = -1.0 * RT * std::log(0.5 * (leftExp + rightExp));
		// */

		std::cout << "ke_prop_nma6 " << ke_prop_nma6 << "\n";



	}

	// Set state velocities
	someState.updU() = sqrtMInvUs;

	// Realize velocity
	system->realize(someState, SimTK::Stage::Velocity);

	// ask for a number of random numbers and check if we are done the next
	//  time we hit this function
	RandomCache.task = std::async(std::launch::async, RandomCache.FillWithGaussian);
}

/*
* Integrate trajectory: Apply the L operator for non-equlibrium processes
*/
void HMCSampler::integrateTrajectory(SimTK::State& someState){
	//std::cout << "HMCSampler::integrateTrajectory: ts mds " << timestep << " " << MDStepsPerSample << std::endl;
	try {
		// Instant geometry
		//geomDihedral(someState);

		// Call Simbody TimeStepper to advance time
		this->world->ts->stepTo(someState.getTime() + (timestep*MDStepsPerSample));
		//this->timeStepper->stepTo(someState.getTime() + (timestep*MDStepsPerSample));
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

// Trajectory length has an average of MDStepsPerSample and a given std
void HMCSampler::integrateVariableTrajectory(SimTK::State& someState){
	try {

	float timeJump = 0;
	if(MDStepsPerSampleStd != 0){
		std::cout << "Got MDSTEPS_STD " << MDStepsPerSampleStd << '\n';

		// TODO moe to internal variables set
			std::normal_distribution<float> trajLenNoiser(1.0, MDStepsPerSampleStd);
			float trajLenNoise = trajLenNoiser(RandomCache.RandomEngine);
		timeJump = (timestep*MDStepsPerSample) * trajLenNoise;
	}else{
		std::cout << "Didn't get MDSTEPS_STD" << '\n';
		timeJump = (timestep*MDStepsPerSample);
	}

		this->timeStepper->stepTo(someState.getTime() + timeJump);

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

	try {
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

		int nu = someState.getNU();
		for(int k = 0; k < MDStepsPerSample; k++){
			this->timeStepper->stepTo(someState.getTime() + (timestep));
			system->realize(someState, SimTK::Stage::Position);

			const auto Q = someState.getQ(); // g++17 complains if this is auto& or const auto&
			const auto U = someState.getU(); // g++17 complains if this is auto& or const auto&
			SimTK::Vector MU(nu, 1.0);
			matter->multiplyByM(someState, U, MU);

			std::cout << "Q" << Q  << std::endl;
			std::cout << "P" << MU << std::endl;

			SimTK::Real KE = matter->calcKineticEnergy(someState);
			SimTK::Real PE = forces->getMultibodySystem().calcPotentialEnergy(someState);
			std::cout << "KE PE " << KE << " " << PE << std::endl;

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

bool b = false;

// ELIZA OPENMM FULLY FLEXIBLE INTEGRATION CODE

// ELIZA: Insert code here
void HMCSampler::OMM_setTemperature(double HMCBoostTemperature){
	assert(!"Not implemented");
}

// ELIZA: Insert code here
double HMCSampler::OMM_calcKineticEnergy(void){
	return dumm->OMM_calcKineticEnergy();
	// assert(!"Not implemented");
}

// ELIZA: Insert code here
double HMCSampler::OMM_calcPotentialEnergy(void){
	return dumm->OMM_calcPotentialEnergy();
}

void HMCSampler::OMM_integrateTrajectory(SimTK::State& someState){
	assert(!"Not implemented");

	try {
		// ELIZA: Insert code here
		dumm->OMM_integrateTrajectory(this->MDStepsPerSample );

	}catch(const std::exception&){
		// Send general message
		std::cout << "\t[WARNING] propose exception caught!\n";
		proposeExceptionCaught = true;

		// Restore configuration
		for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
			const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
			mobod.setQToFitTransform(someState, SetTVector[mbx - 1]);
		}
		system->realize(someState, SimTK::Stage::Position);
	}

}

// ELIZA: Check the code below
void HMCSampler::OMM_calcProposedKineticAndTotalEnergy(void){
	assert(!"Not implemented");

	this->ke_proposed = OMM_calcKineticEnergy();

	// Leave this here 
	this->etot_proposed = getOldPE() 
		+ getProposedKE() 
		+ getOldFixman() 
		+ getOldLogSineSqrGamma2();
}


// ELIZA: Check the code below
void HMCSampler::OMM_calcNewConfigurationAndEnergies(void){
	assert(!"Not implemented");

	// Get new Fixman potential
	if(useFixman){
		std::cerr 
		<< "Attempting Fixman potential calculation with OpenMM integrators.";
		throw std::exception();
		std::exit(1);
		//fix_n = calcFixman(someState);
		//logSineSqrGamma2_n = 
		//	((Topology *)rootTopology)->calcLogSineSqrGamma2(someState);
	}else{
		fix_n = 0.0;
		logSineSqrGamma2_n = 0.0;
	}

	// Get new kinetic energy
	ke_n = OMM_calcKineticEnergy();

	// Get new potential energy
	pe_n = OMM_calcPotentialEnergy();

	// Calculate total energy
	//if(useFixman){
	//	etot_n = pe_n + ke_n + fix_n - (0.5 * RT * logSineSqrGamma2_n);
	//	etot_proposed = pe_o + ke_proposed + fix_o - (0.5 * RT * logSineSqrGamma2_o);
	//}else{
		etot_n = pe_n + ke_n;
		etot_proposed = pe_o + ke_proposed;
	//}

}


/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. **/
bool HMCSampler::propose(SimTK::State& someState)
{

	// Adapt timestep
	if(shouldAdaptTimestep){
		adaptTimestep(someState);
	}

	// Store old configuration
	storeOldConfigurationAndPotentialEnergies(someState);

	if(integratorName == IntegratorName::OMMVV){

		// // ELIZA: Check the code below
		OMM_setTemperature(this->boostT);
		OMM_calcProposedKineticAndTotalEnergy();
		OMM_integrateTrajectory(someState);
		OMM_calcNewConfigurationAndEnergies();

	}else{

		// Initialize velocities according to the Maxwell-Boltzmann distribution
		initializeVelocities(someState);

		// Store the proposed energies
		calcProposedKineticAndTotalEnergy(someState);

		// Adapt timestep
		bool shouldAdaptWorldBlocks = false;
		if(shouldAdaptWorldBlocks){
			adaptWorldBlocks(someState);
		}

		// Apply the L operator
		integrateTrajectory(someState);
		//integrateTrajectoryOneStepAtATime(someState);

		calcNewConfigurationAndEnergies(someState);
	}

	return validateProposal();

}

/* 
* Compute mathematical, rather than robotic Jacobian.
* It translates generalized velocities u into Cartesian velocities
* for each atom on all topologies.
* Computational cost is high
*/
void HMCSampler::calcMathJacobian(const SimTK::State& someState,
	SimTK::Matrix& mathJ)
{
	// Mathematical Jacobian is 3N x nu dimensional
	unsigned int nu = someState.getNU();
	mathJ.resize(natoms * 3, nu);
	mathJ = 0;

	// Row number in math jacobian
	int i = 0;

	// Go through topologies
	for(auto& topology : topologies){
		// Go through atoms
		for(const auto& AtomList : topology.bAtomList){

			// Get atom indeces in Compound and Simbody
			const auto aIx = AtomList.getCompoundAtomIndex();
			const auto mbx = topology.getAtomMobilizedBodyIndexThroughDumm(
				aIx, *dumm);

				// Get atom station on mobod
				SimTK::Vec3 station = 
				topology.getAtomLocationInMobilizedBodyFrameThroughDumm(
				aIx, *dumm);
				
				// Get station Jacobian: 42*nt + 54*nb + 33*n flops
				SimTK::Matrix stationJ;
				matter->calcStationJacobian(someState, mbx, station, stationJ);
				
				// Add new three rows corresponding to this atom
				for(unsigned int j = 0; j < nu; j++){
					mathJ(i + 0, j) = stationJ(0, j);
					mathJ(i + 1, j) = stationJ(1, j);
					mathJ(i + 2, j) = stationJ(2, j);
				}

				i += 3; // increment row number - 3 rows per atom
		}
	}

}

/*
* Get the diagonal 3Nx3N matrix containing the atoms masses
*/
void HMCSampler::getCartesianMassMatrix(const SimTK::State& somestate,
	SimTK::Matrix& M)
{
	// Resize the matrix
	M.resize(natoms * 3, natoms *  3);
	M = 0;
	
	// Row number
	int i = 0;

	// Go through topologies
	for(auto& topology : topologies){
		// Go through atoms
		for(const auto& AtomList : topology.bAtomList){
			
			// Get atom indeces in Compound and DuMM
			const auto aIx = AtomList.getCompoundAtomIndex();
			const auto dAIx = topology.getDuMMAtomIndex(aIx);

			// Put atom mass on the diagonal
			SimTK::Real mass = dumm->getAtomMass(dAIx);
			M(i + 0, i + 0) = mass;
			M(i + 1, i + 1) = mass;
			M(i + 2, i + 2) = mass;

			// Increment row number - 3 rows per atom
			i += 3;
		}
	}


}

/*
* Get Joint type by examining hinge matrix H_FM
*/
int HMCSampler::getJointTypeFromH(const SimTK::State& someState,
	const SimTK::MobilizedBody& mobod)
{
	// Get the mobilized body
	//const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

	// Get number of degrees of freedom
	unsigned int nu = mobod.getNumU(someState);

	if(nu == 1){ // only one degree of freedom
		SpatialVec H_FMCol = mobod.getH_FMCol(someState,
			SimTK::MobilizerUIndex(0));
		std::cout << "H_FMCol " << H_FMCol << std::endl;

	}else{
		// Go through H_FM columns
		for (SimTK::MobilizerUIndex ux(0); ux < nu; ++ux){
			SpatialVec H_FMCol = mobod.getH_FMCol(someState, ux);
			std::cout << "H_FMCol " << H_FMCol << std::endl;
		}
	}


}

void HMCSampler::setDistortOption(const int& distortOptArg)
{
	this->DistortOpt = distortOptArg;
}

void HMCSampler::setQScaleFactor(const SimTK::Real& s)
{
	this->QScaleFactor = s;
}

/*
 * Shift all the generalized coordinates
 */
void HMCSampler::shiftQ(SimTK::State& someState)
{
	// Test
	//std::cout << "Got " << QScaleFactor << " scale factor\n";
	std::cout << "unshifted Q = " << someState.getQ() << std::endl;

	// 1. Update means of values before altering them
	world->updateTransformsMeans(someState);

	//world->PrintX_PFMeans();
	//world->PrintX_BMMeans();

	// 2. Get differences between current transforms and their means
	std::vector<SimTK::Real> X_PFdiffs;
	std::vector<SimTK::Real> X_BMdiffs;
	X_PFdiffs.resize(world->acosX_PF00_means.size(), 0);
	X_BMdiffs.resize(world->normX_BMp_means.size(), 0);

	for(unsigned int k = 0; k < X_PFdiffs.size(); k++){
		X_PFdiffs[k] = world->acosX_PF00[k] - world->acosX_PF00_means[k];
	}
	for(unsigned int k = 0; k < X_BMdiffs.size(); k++){
		X_BMdiffs[k] = world->normX_BMp[k] - world->normX_BMp_means[k];
	}

	// 3. Scale the differences
	int k = -1;
	for(auto& diff : X_PFdiffs){
		diff *= (QScaleFactor - 1);
	}
	for(auto& diff : X_BMdiffs){
		diff *= (QScaleFactor - 1);
	}

	// Print the differences	
/* 	for(unsigned int k = 0; k < X_PFdiffs.size(); k++){
		std::cout << "Excel X_PF " << k << " " << world->acosX_PF00[k] << " " << world->acosX_PF00_means[k] << " " 
			<< X_PFdiffs[k] << std::endl;
		std::cout << "Excel X_BMDiff " << k << " " << world->normX_BMp[k]  << " " << world->normX_BMp_means[k] << " "
			<< X_BMdiffs[k] << std::endl;
	}
	std::cout << "Excel END\n"; */

	// Ground and first mobod don't have internal coordinates
	int offset = 2; 
	for (SimTK::MobilizedBodyIndex mbx(offset);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobilized body
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		
		// we only allocated  X_PFs for non-Ground bodies
		mobod.setOneQ(someState, 0, -1.0 * X_PFdiffs[int(mbx) - 1]);
		mobod.setOneQ(someState, 1, -1.0 * X_BMdiffs[int(mbx) - 1]);
	}

	// Save changes by advancing to Position Stage
	system->realize(someState, SimTK::Stage::Dynamics);

	// Test
	std::cout << "shifted Q = " << someState.getQ() << std::endl;
	//world->getTransformsStatistics(someState);
}

/*
 * Test ground for SOA
 */
void HMCSampler::testSOA(SimTK::State& someState)
{	
	// Some SOA
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){
		
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

		const SpatialVec HColi = mobod.getHCol(someState, SimTK::MobilizerUIndex(0));
		//const SpatialVec H_FMCol = mobod.getH_FMCol(someState, SimTK::MobilizerUIndex(0));
		std::cout << "HColi\n" << HColi << std::endl; // H_PB_G (V_PB_G=H*u)
		//std::cout << "H_FMCol\n" << H_FMCol << std::endl;

		const SpatialInertia Mi_G = mobod.getBodySpatialInertiaInGround(someState); // about origin
		std::cout << "Mi_G\n" << Mi_G.toSpatialMat() << std::endl;

	}	

	// Get System Jacobian
	SimTK::Matrix J_G;
	matter->calcSystemJacobian(someState, J_G);
	//std::cout << "J_G\n" << J_G << std::endl;

	// Get mathematical jacobian
	SimTK::Matrix mathJ;
	calcMathJacobian(someState, mathJ);
	//std::cout << "mathJ\n" << mathJ << std::endl;

	// Get Cartesian mass matrix
	SimTK::Matrix cartM;
	getCartesianMassMatrix(someState, cartM);
	//std::cout << "cartM\n" << cartM << std::endl;

	// Our mass Matrix = mathJ x cartM x mathJ
	SimTK::Matrix myMassMatrix;
	myMassMatrix = mathJ.transpose() * cartM * mathJ;
	//std::cout << "myMassMatrix\n" << myMassMatrix << std::endl;	

 	// Get system mass matrix
	SimTK::Matrix M;
	matter->calcM(someState, M);
	std::cout << "M\n" << M << std::endl << std::flush;


	// Compare metric tensor
	SimTK::Matrix JtJ;
	JtJ = J_G.transpose() * J_G;
	//std::cout << "JtJ\n" << JtJ << std::endl;	

	SimTK::Matrix mathJtJ;
	mathJtJ = mathJ.transpose() * mathJ;
	std::cout << "mathJtJ\n" << mathJtJ << std::endl << std::flush;
	std::cout << "ndofs " << ndofs << std::endl << std::flush;

/* 	
	// linker problems ...
	SimTK::Eigen mathJtJEigen(mathJtJ);
	SimTK::Vector mathJtJEigenVals;
	mathJtJEigen.getAllEigenValues(mathJtJEigenVals);
	std::cout << "mathJtJEigenVals\n" << mathJtJEigenVals << std::endl; */

	// Get mathematical Jacobian square determinant
	SimTK::SymMat<12> smMathJtJ; // BUG: this should be a constant
	for(unsigned int i = 0; i < ndofs; i++){
		for (unsigned int j =0; j < ndofs; j++){
			smMathJtJ(i, j) = mathJtJ(i, j);
			smMathJtJ(j, i) = mathJtJ(i, j);
		}
	}
	SimTK::Real detMathJtJ = SimTK::det(smMathJtJ);
	std::cout << "mathJtJ determinant\n" << detMathJtJ << std::endl;

	// Get mass matrix determinant
	SimTK::SymMat<12> smM; // BUG: this should be a constant
	for(unsigned int i = 0; i < ndofs; i++){
		for (unsigned int j =0; j < ndofs; j++){
			smM(i, j) = M(i, j);
			smM(j, i) = M(i, j);
		}
	}
	SimTK::Real detM = SimTK::det(smM) ;
	std::cout << "M determinant\n" << detM << std::endl;

	SimTK::Real extractM = 0;
	for(unsigned int i = 0; i < cartM.nrow(); i++){
		extractM += std::log(cartM(i, i));
	}

	std::cout << "extractM " << extractM << std::endl;

}



bool HMCSampler::proposeNEHMC(SimTK::State& someState)
{

	// Store old configuration
	storeOldConfigurationAndPotentialEnergies(someState);

	// Apply Lambda protocol: T times
	int T = 1;
	for(int i = 0; i < T; i++){

		// WORK: perform work (alpha)
		if(DistortOpt == -1){
			shiftQ(someState);
		}
		calcNewConfigurationAndEnergies(someState);
		work += (pe_n - pe_o);
		//std::cout << "energies after shiftQ pe_o pe_n ke_prop "
		//	<< pe_o << " " << pe_n << " " << ke_proposed << '\n';

		/*
		// EQUILIBRIUM SIMULATION
		// Initialize velocities according to the Maxwell-Boltzmann distrib
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

		// Propagate with K (no heat if integrator is deterministic)
		//integrateTrajectoryOneStepAtATime(someState);
		integrateTrajectory(someState);
		*/

		calcNewConfigurationAndEnergies(someState);
	}

	return validateProposal();


}

bool HMCSampler::proposeNMA(SimTK::State& someState)
{

	// Store old configuration
	storeOldConfigurationAndPotentialEnergies(someState);

	// Initialize velocities according to the Maxwell-Boltzmann distribution
	initializeNMAVelocities(someState);

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
	integrateVariableTrajectory(someState);

	calcNewConfigurationAndEnergies(someState);

	return validateProposal();


}

/// Acception rejection step
bool HMCSampler::accRejStep(SimTK::State& someState) {
	// Decide and get a new sample
	//if ( getThermostat() == ThermostatName::ANDERSEN ) {
	if ( alwaysAccept == true ){
		// Simple molecular dynamics
		this->acc = true;
		std::cout << "\tsample accepted (simple molecular dynamics)\n";
		update(someState);
	} else {
		// we do not consider this sample accepted unless it passes all checks
		this->acc = false;

		// Apply Metropolis-Hastings correction
		if(acceptSample()) {
			// sample is accepted
			this->acc = true;
			std::cout << "\tsample accepted (metropolis-hastings)\n";
			update(someState);
		} else {
			// sample is rejected
			this->acc = false;
			std::cout << "\tsample rejected (metropolis-hastings)\n";
			setSetConfigurationAndEnergiesToOld(someState);
		}
	}

	return this->acc;
}

// The main function that generates a sample
bool HMCSampler::sample_iteration(SimTK::State& someState)
{
	std::cout << std::setprecision(10) << std::fixed;

	bool validated = false;

	// Propose
	std::cout << "DISTORT_OPTION " << DistortOpt << std::endl;
	for (int i = 0; i < 10; i++){
		if(DistortOpt > 0){
			validated = proposeNMA(someState);
		}else if(DistortOpt < 0){
			validated = proposeNEHMC(someState);
		}else{
			validated = propose(someState);
		}

		if (validated){
			break;
		}
	}

	if(validated){

		PrintDetailedEnergyInfo(someState);

		accRejStep(someState);

		if(this->acc) { // Only in this case ???
			// Add generalized coordinates to a buffer
			auto Q = someState.getQ(); // g++17 complains if this is auto& or const auto&
			//std::cout << "Q = " << Q << std::endl;
			QsBuffer.insert(QsBuffer.end(), Q.begin(), Q.end());
			QsBuffer.erase(QsBuffer.begin(), QsBuffer.begin() + Q.size());

			pushCoordinatesInR(someState);
			pushVelocitiesInRdot(someState);

			// Calculate MSD and RRdot to adapt the integration length
			std::cout << std::setprecision(10) << std::fixed;
			std::cout << "\tMSD= " << calculateMSD() << ", RRdot= " << calculateRRdot() << std::endl;
		}

		++nofSamples;

		return this->acc;
	}else{
		std::cout << "Warning: HMCSampler: Proposal was not validated after 10 tries." << std::endl;
		setSetConfigurationAndEnergiesToOld(someState);

		++nofSamples;

		return false;
	}

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
	if(argThermostat == ThermostatName::ANDERSEN){
		std::cout << "Adding Andersen thermostat.\n";
	}
	this->thermostat = argThermostat;
}
// Set a thermostat
void HMCSampler::setThermostat(std::string thermoName){
	thermoName.resize(thermoName.size());
	std::transform(thermoName.begin(), thermoName.end(), thermoName.begin(), ::tolower);

	std::cout << "Adding " << thermoName << " thermostat.\n";

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
SimTK::Real HMCSampler::getPEFromEvaluator(const SimTK::State& someState) const{
		return forces->getMultibodySystem().calcPotentialEnergy(someState);
	// Eliza's potential energy's calculations including rigid bodies
	// internal energy
	//return dumm->CalcFullPotEnergyIncludingRigidBodies(someState);

}

SimTK::Real HMCSampler::calcFixman(SimTK::State& someState){
	//std::cout << "DEBUG HMCSampler::calcFixman ";

	int nu = someState.getNU();
	SimTK::Vector V(nu);

	system->realize(someState, SimTK::Stage::Position);
	matter->realizeArticulatedBodyInertias(someState); // Move in calcDetM ?

	// Get detM
	SimTK::Vector DetV(nu);
	SimTK::Real D0 = 0.0;

	// TODO: remove the request for Dynamics stage cache in SImbody files
	matter->calcDetM(someState, V, DetV, &D0);
	//std::cout << "HMCSampler::calcFixman D0= "<< D0 << std::endl;

	assert(RT > SimTK::TinyReal);
	double detMBAT = ((Topology *)rootTopology)->calcLogDetMBATInternal(someState);
	//SimTK::Real result = 0.5 * RT * ( std::log(D0) - detMBAT ); // original
	SimTK::Real result = 0.5 * RT * ( D0 - detMBAT ); // log space already
	//std::cout << "detM detMBAT fixP " << D0 << " " << detMBAT << " " << result << std::endl;

	if(SimTK::isInf(result)){
		std::cout << "Fixman potential is infinite!\n";
		//result = 0.0;
	}

	// Get M
	//SimTK::Matrix M(nu, nu); // TODO: DELETE
	//matter->calcM(someState, M); // TODO: DELETE
    	//PrintBigMat(M, nu, nu, 12, "M");  // TODO: DELETE

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
  system->realize(someState, SimTK::Stage::Position);
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
	//system->realize(someState, SimTK::Stage::Position);
	const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
	mobod.setQToFitTransform(someState, SetTVector[i]);

	//system->realize(someState, SimTK::Stage::Position);


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
	InvUScaleFactors.resize(nu, 1);

	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		// const int mnu = mobod.getNumU(someState);
		const float scaleFactor = world->getMobodUScaleFactor(mbx);

		//std::cout << "loadUScaleFactors mbx scaleFactor uIxes " << int(mbx) << ' ' << scaleFactor;
		for(SimTK::UIndex uIx = mobod.getFirstUIndex(someState);
		uIx < mobod.getFirstUIndex(someState) + mobod.getNumU(someState);
		uIx++ ){
			//std::cout << ' ' << int(uIx) ;
			UScaleFactors[int(uIx)] = scaleFactor;
			InvUScaleFactors[int(uIx)] = 1.0 / scaleFactor;
		}
		//std::cout << '\n';
	}

	// Compute UScalefactors norm and its inverse
	// We don't want to normalize it though. Allow the user the liberty
	// to enforce its own velocities
	UScaleFactorsNorm = 0.0;
	InvUScaleFactorsNorm = 0.0;
	for (int i = 0; i < UScaleFactors.size(); ++i) {
		UScaleFactorsNorm += UScaleFactors[i] * UScaleFactors[i];
		InvUScaleFactorsNorm += InvUScaleFactors[i] * InvUScaleFactors[i];
	}
	UScaleFactorsNorm = std::sqrt(UScaleFactorsNorm);
	InvUScaleFactorsNorm = std::sqrt(InvUScaleFactorsNorm);

	// Direction vector UScaleFactors
	NormedUScaleFactors.resize(nu);
	DOFUScaleFactors.resize(nu);
	for(int i = 0; i < nu; i++){
		NormedUScaleFactors[i] = UScaleFactors[i] / UScaleFactorsNorm;
	}
	for(int i = 0; i < nu; i++){
		DOFUScaleFactors[i] = NormedUScaleFactors[i] * std::sqrt(nu);
	}


/* THIS IS A HUGE MEMORY PROBLEM
	// Set the NMA rotation matrix to identity
	NMARotation.resize(nu, std::vector<double>(nu, 0));

	for(int i = 0; i < nu; i++){
		NMARotation[i][i] = 1.0;
	}

	// Fill first vector of the NMA rotation matrix
	//std::cout << "NMARotation 1st column filled with \n";
	for(int i = 0; i < nu; i++){
		//std::cout << (UScaleFactors[i] / UScaleFactorsNorm) << " ";
		NMARotation[i][0] = (UScaleFactors[i] / UScaleFactorsNorm);
	}
	//std::cout << '\n';

	bMatrix tempNMARotation;
	tempNMARotation.resize(nu, std::vector<double>(nu, 0));
	bCopyMat(NMARotation, tempNMARotation);

	//gram_schmidt(tempNMARotation, NMARotation);
*/

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
	//std::cout << "HMC: boost temperature: " << this->boostT << std::endl;

	this->boostRT = this->boostT * SimTK_BOLTZMANN_CONSTANT_MD;
	//std::cout << "HMC: boostRT: " << this->boostRT << std::endl;

	this->sqrtBoostRT = std::sqrt(this->boostRT);
	//std::cout << "HMC: sqrtBoostRT: " << this->sqrtBoostRT << std::endl;

	this->boostBeta = 1.0 / boostRT;
	//std::cout << "HMC: boostBeta: " << this->boostBeta << std::endl;

	this->boostFactor = std::sqrt(this->boostT / this->temperature);
	this->unboostFactor = 1 / boostFactor;
	//std::cout << "HMC: boost velocity scale factor: " << this->boostFactor << std::endl;
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

/** Store the proposed energies **/
void HMCSampler::calcProposedKineticAndTotalEnergy(SimTK::State& someState){

	// Store proposed kinetic energy
	// setProposedKE(matter->calcKineticEnergy(someState));

	// OLD RESTORE TODO BOOST
	// this->ke_proposed = matter->calcKineticEnergy(someState);
	// NEW TRY TODO BOOST
	this->ke_proposed = this->unboostFactor * matter->calcKineticEnergy(someState);

	// Store proposed total energy
	this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();
}

// Stochastic optimization of the timestep using gradient descent
void HMCSampler::adaptTimestep(SimTK::State&)
{

	static int nofTimestepAdapt = 0;

	if( (nofSamples % acceptedStepsBufferSize) == (acceptedStepsBufferSize-1) ){

		nofTimestepAdapt += 1;

		std::cout << "Adapt BEGIN " << nofTimestepAdapt << ": ";

		//SimTK::Real idealAcceptance = 0.651;
		SimTK::Real idealAcceptance = 0.8;
		SimTK::Real newTimestep = SimTK::NaN;

		// Compute acceptance in the buffer
		int sum = std::accumulate(acceptedStepsBuffer.begin(),
			acceptedStepsBuffer.end(), 0);
		SimTK::Real newAcceptance = float(sum) / float(acceptedStepsBufferSize);

		// Advance timestep ruler
		prevPrevAcceptance = prevAcceptance;
		prevAcceptance = acceptance;
		acceptance = newAcceptance;

		std::cout << " ppAcc= " << prevPrevAcceptance
			<< " pAcc= " << prevAcceptance
			<< " acc= " << acceptance
			<< " ppTs= " << prevPrevTimestep
			<< " pTs= " << prevTimestep
			<< " ts= " << timestep
			<< std::endl;

		// Passed first two initial evaluations
		if( !SimTK::isNaN(prevPrevAcceptance) ){
			// Calculate gradients
			SimTK::Real a_n, a_n_1, a_n_2, t_n, t_n_1, t_n_2, f_n, f_n_1, f_n_2;
			a_n = acceptance; a_n_1 = prevAcceptance; a_n_2 = prevPrevAcceptance;
			t_n = timestep; t_n_1 = prevTimestep; t_n_2 = prevPrevTimestep;
			f_n   = a_n     - idealAcceptance;
			f_n_1 = a_n_1   - idealAcceptance;
			f_n_2 = a_n_2   - idealAcceptance;

			// f = a - aref and F = integral(f)
			// minimize F <=> f = 0 <=> a = aref
			SimTK::Real F_n =   (0.5 * (a_n * a_n))     - (a_n   * idealAcceptance);
			SimTK::Real F_n_1 = (0.5 * (a_n_1 * a_n_1)) - (a_n_1 * idealAcceptance);
			SimTK::Real F_n_2 = (0.5 * (a_n_2 * a_n_2)) - (a_n_2 * idealAcceptance);

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

			// Generate a new timestep adaptedTimestep
			SimTK::Real adaptedTimestep = t_n - (gamma * gradF_n);

			std::cout << "Adapt data:"
				<< " F_n= " << F_n  << " F_n_1= " << F_n_1  << " F_n_2= " << F_n_2
				<< " dF_n= " << dF_n << " dF_n_1= " << dF_n_1
				<< " dt_n= " << dt_n << " dt_n_1= " << dt_n_1
				<< " gradF_n= " << gradF_n << " gradF_n_1= " << gradF_n_1
				<< " num " << num << " " << denom << " gamma= " << gamma
				<< std::endl;

			// Acceptance is fine
			if(	   (std::abs(f_n) < 0.1)
				&& (std::abs(f_n_1) < 0.1)
				&& (std::abs(f_n_2) < 0.1)
			){
				std::cout << "Acceptance is nearly ideal.\n";

				newTimestep = timestep;

			// Not enough data to compute gradient
			//}else if((dF_n == 0.0) || (dF_n_1 == 0.0)){
			}else if(SimTK::isNaN(adaptedTimestep)){
				std::cout << "Adding noise to newTimestep to become: ";

				SimTK::Real r = uniformRealDistribution(randomEngine);
				if(acceptance < idealAcceptance){
					newTimestep = timestep - (0.3*r * timestep);
				}else{
					newTimestep = timestep + (0.3*r * timestep);
				}

				std::cout << newTimestep << "\n";

			} else { // Take the new timestep
				std::cout << "Setting newTimestep to adaptedTimestep.\n";
				newTimestep = adaptedTimestep;
			}

		} else { // Alter the intial timesteps to get a valid dF_n next time

			//if( SimTK::isNaN(prevPrevTimestep) ){
			if( !SimTK::isNaN(prevAcceptance) ){
				std::cout << "Third step\n";
				newTimestep = (timestep + prevTimestep) / 2;
			}else{
				std::cout << "Second step\n";
				newTimestep = 0.0009; // smaller than H vibration
			}
		}

		// Advance timestep ruler
		prevPrevTimestep = prevTimestep;
		prevTimestep = timestep;
		timestep = newTimestep;

		if(timestep < 0.0001){
			timestep = 0.0009;
		}

		// Set integrator timestep to newTimestep
		timeStepper->updIntegrator().setFixedStepSize(timestep);

		std::cout << "Adapt END: "  << " ppTs= " << prevPrevTimestep
			<< " pTs= " << prevTimestep << " ts= " << timestep << std::endl;

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


void HMCSampler::geomDihedral(SimTK::State& someState){
	 // INSTANT GEOMETRY
	SimTK::Vec3 a1pos, a2pos, a3pos, a4pos, a5pos;
	int a1, a2, a3, a4, a5;
	//DIHEDRAL 4 6 8 14 6 8 14 16
	//a1 = 16; a2 = 14; a3 = 0; a4 = 6; a5 = 8;
	a1 = 4; a2 = 6; a3 = 8; a4 = 14; a5 = 16;
	a1pos = ((Topology *)rootTopology)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
	a2pos = ((Topology *)rootTopology)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));
	a3pos = ((Topology *)rootTopology)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)));
	a4pos = ((Topology *)rootTopology)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a4)));
	//a5pos = ((Topology *)rootTopology)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a5)));
	int distA1, distA2, distA3, distA4;
	SimTK::Vec3 distA1pos, distA2pos;
	SimTK::Vec3 distA3pos, distA4pos;
	distA1 = 2; distA2 = 17; distA3 = 6; distA4 = 17;
		    distA1pos = ((Topology *)rootTopology)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA1)));
		    distA2pos = ((Topology *)rootTopology)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA2)));
		    distA3pos = ((Topology *)rootTopology)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA3)));
		    distA4pos = ((Topology *)rootTopology)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA4)));
//            std::cout << " dihedral elements"
//              << " " << rootTopology->getAtomElement(SimTK::Compound::AtomIndex(a1))
//              << " " << rootTopology->getAtomElement(SimTK::Compound::AtomIndex(a2))
//              << " " << rootTopology->getAtomElement(SimTK::Compound::AtomIndex(a3))
//              << " " << rootTopology->getAtomElement(SimTK::Compound::AtomIndex(a4))
//              << " " << rootTopology->getAtomElement(SimTK::Compound::AtomIndex(a5))
//              << std::endl;
//            std::cout << " poss: " << a1pos << ' ' << a2pos << ' ' << a3pos << ' ' << a4pos << ' ';
	std::cout << "geom "  << bDihedral(a1pos, a2pos, a3pos, a4pos) ;
//	std::cout << " "  << bDihedral(a2pos, a3pos, a4pos, a5pos) ;
	//std::cout << " "  << 10 * (distA1pos - distA2pos).norm() ; // Ang
	//std::cout << " "  << 10 * (distA3pos - distA4pos).norm() ; // Ang
	std::cout << std::endl;
	// */
	// INSTANT GEOMETRY END
}

/** Store new configuration and energy terms **/
void HMCSampler::calcNewConfigurationAndEnergies(SimTK::State& someState)
{
	// Get new Fixman potential
	if(useFixman){
		fix_n = calcFixman(someState);
		logSineSqrGamma2_n = ((Topology *)rootTopology)->calcLogSineSqrGamma2(someState);
	}else{
		fix_n = 0.0;
		logSineSqrGamma2_n = 0.0;
	}

	// Get new kinetic energy
	system->realize(someState, SimTK::Stage::Velocity);

	// OLD RESTORE TODO BOOST
	//ke_n = matter->calcKineticEnergy(someState);
	// NEW TRY TODO BOOST
	ke_n = this->unboostFactor * matter->calcKineticEnergy(someState);

	// Get new potential energy
	pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
	//pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState);
	// TODO: replace with the following after checking is the same thing
	//pe_n = compoundSystem->calcPotentialEnergy(someState);

	logSineSqrGamma2_n = ((Topology *)rootTopology)->calcLogSineSqrGamma2(someState);

	// Calculate total energy
	if(useFixman){
		etot_n = pe_n + ke_n + fix_n - (0.5 * RT * logSineSqrGamma2_n);
		etot_proposed = pe_o + ke_proposed + fix_o - (0.5 * RT * logSineSqrGamma2_o);
	}else{
		etot_n = pe_n + ke_n;
		etot_proposed = pe_o + ke_proposed;
	}

		/*// Kinetic energy new for the evaluation Hamiltonian
		const int nu = someState.getNU();
		ke_n_nma6 = 0;
		SimTK::Vector v_miu(nu, 1);
		SimTK::Vector v_miuM(nu, 1);

		SimTK::Real leftExp = 0.0;
		SimTK::Real rightExp = 0.0;

		// Left exponential
		v_miu = (someState.getU() - DOFUScaleFactors);
		std::cout << "v_miu 0 "; for(int i = 0; i < nu; i++){ std::cout << v_miu[i] << " " ;} std::cout << "\n";

		matter->multiplyByM(someState, v_miu, v_miuM);
		std::cout << "v_miuM 1 "; for(int i = 0; i < nu; i++){ std::cout << v_miuM[i] << " " ;} std::cout << "\n";

		for(int j = 0; j < ndofs; j++){leftExp += (v_miuM[j] * v_miu[j]);}
		std::cout << "leftExp 0 " << leftExp << "\n";

		//leftExp *= ((-0.5) * 1.4142135623730950488 * this->beta);
		leftExp *= ((-0.5) * 1.4142135623730950488 * this->boostBeta);
		std::cout << "leftExp 1 " << leftExp << "\n";

		leftExp = std::exp(leftExp);
		std::cout << "leftExp 2 " << leftExp << "\n";

		// Right exponential
		v_miu = (someState.getU() + DOFUScaleFactors);
		std::cout << "v_miu 0 "; for(int i = 0; i < nu; i++){ std::cout << v_miu[i] << " " ;} std::cout << "\n";

		matter->multiplyByM(someState, v_miu, v_miuM);
		std::cout << "v_miuM 1 "; for(int i = 0; i < nu; i++){ std::cout << v_miuM[i] << " " ;} std::cout << "\n";

		for(int j = 0; j < ndofs; j++){rightExp += (v_miuM[j] * v_miu[j]);}
		std::cout << "rightExp 0 " << rightExp << "\n";

		//rightExp *= ((-0.5) * 1.4142135623730950488 * this->beta);
		std::cout << "beta boostBeta " << this->beta << " " << this->boostBeta << "\n" ;
		rightExp *= ((-0.5) * 1.4142135623730950488 * this->boostBeta);
		std::cout << "rightExp 1 " << rightExp << "\n";

		rightExp = std::exp(rightExp);
		std::cout << "rightExp 2 " << rightExp << "\n";

		ke_n_nma6 = -1.0 * RT * std::log(0.5 * (leftExp + rightExp));
		// */

		// Probability terms
		const int nu = someState.getNU();
		SimTK::Real boostBetaPlus = this->boostBeta * 1.4142135623730950488;
		SimTK::Real localBoostFactor = boostBetaPlus / this->beta;

		// Kinetic term
		SimTK::Vector X(someState.getU());
		SimTK::Vector XM(nu, 1);
		SimTK::Real XMX = 0.0;
		SimTK::Real bXMX = 0.0;

		matter->multiplyByM(someState, X, XM);
		for(int j = 0; j < ndofs; j++){XMX += (XM[j] * X[j]);}
		bXMX = localBoostFactor * XMX;

		// Bias term
		//SimTK::Vector DOFUScaleFactorsM(nu, 1);
		//SimTK::Real muMmu = 0.0;
		//SimTK::Real bmuMmu = 0.0;
		//matter->multiplyByM(someState, DOFUScaleFactors, DOFUScaleFactorsM);
		//for(int j = 0; j < ndofs; j++){muMmu += (DOFUScaleFactorsM[j] * DOFUScaleFactors[j]);}
		//bmuMmu = localBoostFactor * muMmu;

		// Correlation term
		SimTK::Real muMX = 0.0;
		SimTK::Real bmuMX = 0.0;
		SimTK::Real _bmuMX = 0.0;
		SimTK::Real bcorr = 0.0;
		for(int j = 0; j < ndofs; j++){muMX += (XM[j] * DOFUScaleFactors[j]);}
		bmuMX = boostBetaPlus * muMX;
		_bmuMX = -1.0 * bmuMX;

		if(bcorr < 16){ // safe limit to use exponentials
			bcorr = LSE2(bmuMX, _bmuMX);
		}else{
			bcorr = (bmuMX >= 0 ? bmuMX : _bmuMX); // approximation to LSE
		}
		bcorr /= this->beta;

		std::cout << "Keval terms end bXMX muMmu bmuMX LSE(bmuMX) = "
			<< 0.5*bXMX << " " // << 0.5*bmuMmu << " "
			<< bmuMX << " " << bcorr << std::endl;

		ke_n_nma6 = // (-1.0 * RT * std::log(0.5))
			+ (0.5*bXMX) // + (0.5* bmuMmu)
			- bcorr;

		std::cout << "ke_n_nma6 " << ke_n_nma6 << "\n";

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
	SimTK::Real prob = 0.0;
	if(DistortOpt == 0){
		prob = MHAcceptProbability(etot_proposed, etot_n + work);
	}else if(DistortOpt > 0){

		if(useFixman){
			prob = MHAcceptProbability(pe_o + ke_prop_nma6 + fix_o,
									   pe_n + ke_n_nma6    + fix_n + work);
			//prob = MHAcceptProbability(pe_o + fix_o, pe_n + fix_n);
		}else{
			prob = MHAcceptProbability(pe_o + ke_prop_nma6, pe_n + ke_n_nma6 + work);
			//prob = MHAcceptProbability(pe_o, pe_n);
		}
	}else{
		prob = MHAcceptProbability(etot_proposed, etot_n + work);
	}

	// std::cout << "\trand_no=" << rand_no << ", prob=" << prob << ", beta=" << beta << std::endl;
	return rand_no < prob;
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
void HMCSampler::PrintDetailedEnergyInfo(const SimTK::State& someState) const
{
	std::cout << std::setprecision(5) << std::fixed;
	std::cout
		<< "\tpe_o " << pe_o << ", pe_n " << pe_n << ", pe_nB " << /*" turned off for time being"*/ getPEFromEvaluator(someState)
		<< "\n\tke_prop " << ke_proposed << ", ke_n " << ke_n
		<< "\n\tfix_o " << fix_o << ", fix_n " << fix_n
		<< "\n\tlogSineSqrGamma2_o " << logSineSqrGamma2_o << ", logSineSqrGamma2_n " << logSineSqrGamma2_n
		//<< " detmbat_n " << detmbat_n //<< " detmbat_o " << detmbat_o << " "
		<< "\n\tts " << timestep  //<< ", exp(bdE) " << exp(-(etot_n - etot_proposed) / RT)
		<< "\n\t, etot_n " << etot_n  << ", etot_proposed " << etot_proposed
		<< ", work " << work
		<< std::endl;
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
