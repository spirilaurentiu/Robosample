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

	// For RANDOM_WALK Docking Simulations
	this->sphereRadius=0.0;

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
	bendStretchJacobianDetLog = 0.0;

}

/** Destructor **/
HMCSampler::~HMCSampler()
{
}

/** ===============================
 * RANDOM NUMBERS
    =============================== */


/** Generate a random number from a uniform distribution
 * with limits L and R
*/
SimTK::Real HMCSampler::uniformRealDistributionRandTrunc(
	SimTK::Real L, SimTK::Real R)
{
	SimTK::Real r = uniformRealDistribution(randomEngine);
	
	return r * (R - L) + L;
}

/** Get the PDF of a random number from a uniform distribution
 * with limits L and R
*/
SimTK::Real HMCSampler::uniformRealDistributionPDFTrunc(
	SimTK::Real X, SimTK::Real L, SimTK::Real R)
{
	SimTK::Real pdf = SimTK::NaN;
	
	if ((X >= L) && (X <= R)){
		pdf = 1.0 / (R - L);
	}

	return pdf;
}

/** Get the CDF of a random number from a uniform distribution
 * with limits L and R
*/
SimTK::Real HMCSampler::uniformRealDistributionCDFTrunc(
	SimTK::Real X, SimTK::Real L, SimTK::Real R)
{
	SimTK::Real cdf = SimTK::NaN;
	
	if (X < L){
		cdf = 0.0;
	}
	else if(X > R){
		cdf = 1.0;
	}else{
		cdf = (X - L) / (R - L);
	}

	return cdf;
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

	// Get no of degrees of freedom
	const int nu = someState.getNU();

	// Potential energy and configuration need stage Position
	system->realize(someState, SimTK::Stage::Position);

	// Set old potential energy
	this->pe_o = this->pe_set = forces->getMultibodySystem().calcPotentialEnergy(
		someState);

	// Store the configuration
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx)
	{
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		SetTVector[mbx - 1] = 
			//TVector[mbx - 1] = 
			mobod.getMobilizerTransform(someState);
	}

	// Initialize QsBuffer with zeros
	int nq = matter->getNQ(someState);
	int totSize = QsBufferSize * nq;
	for(int i = 0; i < totSize; i++){
		//QsBuffer.push_back(SimTK::Vector(nq, SimTK::Real(0)));
		QsBuffer.push_back(SimTK::Real(0));
	}

/* 	// Store potential energies
	//setOldPE(getPEFromEvaluator(someState));
	setOldPE(forces->getMultibodySystem().calcPotentialEnergy(someState));
	//setOldPE(dumm->CalcFullPotEnergyIncludingRigidBodies(someState));
	setSetPE(getOldPE()); */

	// Store Fixman potential
	if(useFixman){
		std::cout << "Hamiltonian Monte Carlo sampler: using Fixman potential.\n";
		this->fix_o = calcFixman(someState);
		this->fix_set = this->fix_o;

		setOldLogSineSqrGamma2(
			((Topology *)rootTopology)->calcLogSineSqrGamma2(someState));
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
	}else{
		this->fix_o = 0.0;
		this->fix_set = this->fix_o;

		setOldLogSineSqrGamma2(0.0);
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
	}

	// Initialize velocities to temperature
/* 	double sqrtRT = std::sqrt(RT);
	SimTK::Vector V(nu);
	SimTK::Vector SqrtMInvV(nu);
	for (int j=0; j < nu; ++j){
		V[j] = gaurand(randomEngine);
	}
	matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);

	SqrtMInvV *= sqrtRT; // Set stddev according to temperature
	someState.updU() = SqrtMInvV;
	system->realize(someState, SimTK::Stage::Velocity); */

	initializeVelocities(someState);

/* 	// Store kinetic energies
	setOldKE(matter->calcKineticEnergy(someState));
	setSetKE(getProposedKE());

	// Store total energies
	this->etot_proposed = getOldPE() + getProposedKE() 
		+ getOldFixman() + getOldLogSineSqrGamma2();
	this->etot_set = this->etot_proposed; */

	// Set the generalized velocities scale factors
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
	bendStretchJacobianDetLog = 0.0;

}

/** Same as initialize **/
void HMCSampler::reinitialize(SimTK::State& someState)
{
	// After an event handler has made a discontinuous change to the
	// Integrator's "advanced state", this method must be called to
	// reinitialize the Integrator.
	//(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

	// Potential energy and configuration need stage Position
	system->realize(someState, SimTK::Stage::Position);

	// Set old potential energy
	this->pe_o = forces->getMultibodySystem().calcPotentialEnergy(
		someState);

	// Store the configuration
	int i = 0;
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		SetTVector[i] = mobod.getMobilizerTransform(someState); //= TVector[i]
		i++;
	}

	// Store potential energies
	setSetPE(getOldPE());

	// Store Fixman potential
	if(useFixman){
		this->fix_o = calcFixman(someState);
		this->fix_set = this->fix_o;

		setOldLogSineSqrGamma2( 
			((Topology *)rootTopology)->calcLogSineSqrGamma2(someState));
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
	}else{
		this->fix_o = 0.0;
		this->fix_set = this->fix_o;

		setOldLogSineSqrGamma2(0.0);
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
	}

	// Initialize velocities to temperature
	int nu = someState.getNU();

	// Reset ndofs
	ndofs = nu;

	// Initialize velocities to temperature
	/*// kT
	double sqrtRT = std::sqrt(RT);
	SimTK::Vector V(nu);
	SimTK::Vector SqrtMInvV(nu);
	for (int j=0; j < nu; ++j){
		V[j] = gaurand(randomEngine);
	}
	matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
	SqrtMInvV *= sqrtRT; // Set stddev according to temperature
	someState.updU() = SqrtMInvV;
	system->realize(someState, SimTK::Stage::Velocity); */

	//initializeVelocities(someState);

/* 	// Store kinetic energies
	setProposedKE(matter->calcKineticEnergy(someState));
	setLastAcceptedKE(getProposedKE());

	// Store total energies
	this->etot_proposed = getOldPE() + getProposedKE()
		+ getOldFixman() + getOldLogSineSqrGamma2();
	this->etot_set = this->etot_proposed; */

	for (int j=0; j < nu; ++j){
        UScaleFactors[j] = 1;
        InvUScaleFactors[j] = 1;
	}

	// Set the generalized velocities scale factors
	loadUScaleFactors(someState);

	// Total mass of the system
	this->totalMass = matter->calcSystemMass(someState);
	
	// Work
	bendStretchJacobianDetLog = 0.0;

}


SimTK::Real HMCSampler::getMDStepsPerSampleStd() const
{
	return MDStepsPerSampleStd;
}

void HMCSampler::setMDStepsPerSampleStd(SimTK::Real mdstd){
	MDStepsPerSampleStd = mdstd;
}

// Set the method of integration
void HMCSampler::setSampleGenerator(const std::string& generatorNameArg)
{
	if (generatorNameArg == "EMPTY"){
		this->sampleGenerator = 0;
		setAlwaysAccept(true);
	}
	else if(generatorNameArg == "MC"){
		this->sampleGenerator = 1;
		setAlwaysAccept(false);
	}else{
		std::cerr << "Unknown sampling method.\n";
		throw std::exception(); std::exit(1);
	}
}

void HMCSampler::setIntegratorName(IntegratorName integratorNameArg)
{
	this->integratorName = integratorNameArg;
}

void HMCSampler::setIntegratorName(const std::string integratorNameArg)
{

	//this->integratorName = IntegratorNameS[integratorNameArg];

 	if(integratorNameArg == "OMMVV"){
		this->integratorName = IntegratorName::OMMVV;

	}else if (integratorNameArg == "VV"){
		integratorName = IntegratorName::VERLET;

	}else if (integratorNameArg == "BOUND_WALK"){
		integratorName = IntegratorName::BOUND_WALK;
	
	}else if (integratorNameArg == "BOUND_HMC"){
		integratorName = IntegratorName::BOUND_HMC;
	}else if(integratorNameArg == "STATIONS_TASK"){
		integratorName = IntegratorName::STATIONS_TASK;
	}else{
		integratorName = IntegratorName::EMPTY;

	}

}

void HMCSampler::setVelocitiesToZero(SimTK::State& someState){
	
	// Set velocities to 0
	someState.updU() = 0.0;
}

/** Initialize velocities according to the Maxwell-Boltzmann
 * distribution.  Coresponds to R operator in LAHMC.
 * It takes the boost factor into account 
 * */
void HMCSampler::initializeVelocities(SimTK::State& someState){

	if(this->boostRT == 0){

		setVelocitiesToZero(someState);

	}else{
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
		matter->multiplyBySqrtMInv(someState,
			RandomCache.V,
			RandomCache.SqrtMInvV);

		// Set stddev according to temperature
		RandomCache.SqrtMInvV *= sqrtBoostRT;

		// Raise the temperature
		someState.updU() = RandomCache.SqrtMInvV;

	}

	// Realize velocity
	system->realize(someState, SimTK::Stage::Velocity);

	// Store kinetic energies
	this->ke_o = (this->unboostKEFactor) * matter->calcKineticEnergy(someState);
	this->ke_set = this->ke_o;

	// Update total energies
	this->etot_o = getOldPE() + getOldKE()
		+ getOldFixman() + getOldLogSineSqrGamma2();
	this->etot_set = this->etot_o;

	// ask for a number of random numbers and check if we are done the next
	//  time we hit this function
	RandomCache.task = std::async(std::launch::async,
		RandomCache.FillWithGaussian);
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

		generateVonMisesFisherSample(X, concentration);
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
	// assert(!"Not implemented");

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

	this->ke_o = OMM_calcKineticEnergy();

	// Leave this here 
	this->etot_o = getOldPE() 
		+ getOldKE() 
		+ getOldFixman() 
		+ getOldLogSineSqrGamma2();
}


// ELIZA: Check the code below
void HMCSampler::OMM_calcNewConfigurationAndEnergies(void){
	// assert(!"Not implemented");

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
		etot_o = pe_o + ke_o;
	//}

}

/* 
* Compute mathematical, rather than robotic Jacobian.
* It translates generalized velocities u into Cartesian velocities
* for each atom on all topologies.
* Computational cost is high
*/
SimTK::Matrix&
HMCSampler::calcMathJacobian(const SimTK::State& someState,
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

	return mathJ;

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
SimTK::Real HMCSampler::getNewPE(void) const
{
	return this->pe_n;
}

// Set set potential energy
void HMCSampler::setNewPE(SimTK::Real argPE)
{
	this->pe_n = argPE;
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

SimTK::Real HMCSampler::getDistortJacobianDetLog(void) const
{
	// Accumulate all transformation Jacobians in here
	SimTK::Real retValue = 0.0;

	// Bend-Stretch joints
	retValue += this->bendStretchJacobianDetLog;


	return retValue;
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


/** Stores the accepted kinetic energy. This should be set right after a
move is accepted. It's a component of the total energy stored. **/
void HMCSampler::setSetKE(SimTK::Real inpKE)
{
	this->ke_set = inpKE;
}

/** Sets the proposed kinetic energy before the proposal. This should be
set right after the velocities are initialized. **/
void HMCSampler::setOldKE(SimTK::Real inpKE)
{
	this->ke_o = inpKE;
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

	this->boostKEFactor = (this->boostT / this->temperature);
	this->unboostKEFactor = 1 / boostKEFactor;

	this->boostUFactor = std::sqrt(this->boostT / this->temperature);
	this->unboostUFactor = 1 / boostUFactor;
	//std::cout << "HMC: boost velocity scale factor: " << this->boostUFactor << std::endl;
}

void HMCSampler::setBoostMDSteps(int argMDSteps)
{
	this->boostMDSteps = argMDSteps;
	std::cout << "HMC: boost MD steps: " << this->boostMDSteps << std::endl;

}

/** 
 * Store configuration as a set of Transforms and 
 * potential energies
 */
void
HMCSampler::storeOldConfigurationAndPotentialEnergies(
	SimTK::State& someState)
{ 
	// Ensure stage Position is realized
	system->realize(someState, SimTK::Stage::Position);

	// Set old potential energy terms to the last set ones
/* 	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){

		TVector[mbx - 1] = SetTVector[mbx - 1];
	} */

	// Set old potential energy terms to the last set ones
	pe_o = pe_set;
	fix_o = fix_set;
	logSineSqrGamma2_o = logSineSqrGamma2_set;

}

/** Store the proposed energies **/
void HMCSampler::calcProposedKineticAndTotalEnergy(SimTK::State& someState){

	// Store proposed kinetic energy
	// setProposedKE(matter->calcKineticEnergy(someState));

	// OLD RESTORE TODO BOOST
	// this->ke_proposed = matter->calcKineticEnergy(someState);
	// NEW TRY TODO BOOST
	this->ke_o = (this->unboostKEFactor) * matter->calcKineticEnergy(someState);

	// Store proposed total energy
	this->etot_o = getOldPE() + getOldKE() 
		+ getOldFixman() + getOldLogSineSqrGamma2();
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
	ke_n = (this->unboostKEFactor) * matter->calcKineticEnergy(someState);

	// Get new potential energy
	pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
	//pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState);
	// TODO: replace with the following after checking is the same thing
	//pe_n = compoundSystem->calcPotentialEnergy(someState);

	// Calculate total energy
	if(useFixman){
		etot_n = pe_n + ke_n + fix_n - (0.5 * RT * logSineSqrGamma2_n);
		etot_o = pe_o + ke_o + fix_o - (0.5 * RT * logSineSqrGamma2_o);
	}else{
		etot_n = pe_n + ke_n;
		etot_o = pe_o + ke_o;
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
		std::cout << "v_miu 0 "; for(int i = 0; i < nu; i++){ 
			std::cout << v_miu[i] << " " ;} std::cout << "\n";

		matter->multiplyByM(someState, v_miu, v_miuM);
		std::cout << "v_miuM 1 "; for(int i = 0; i < nu; i++){
			std::cout << v_miuM[i] << " " ;} std::cout << "\n";

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
void HMCSampler::setSetConfigurationAndEnergiesToNew(
	SimTK::State& someState)
{
	// Store new configuration
	setSetTVector(someState);

	// Store new energies
	pe_set = pe_n;
	fix_set = fix_n;
	logSineSqrGamma2_set = logSineSqrGamma2_n;
	ke_set = ke_n;

	// Set final total energies
	etot_set = pe_set + fix_set + ke_set + logSineSqrGamma2_set;
}

void HMCSampler::setSetConfigurationToOld(
	SimTK::State& someState)
{
	assignConfFromSetTVector(someState);
	proposeExceptionCaught = false;
}

/** Restore configuration and set energies to old */
void HMCSampler::setSetConfigurationAndEnergiesToOld(
	SimTK::State& someState)
{
	// Restore old configuration
	assignConfFromSetTVector(someState);
	proposeExceptionCaught = false;

	// Set final energies to the precalculated old ones
	pe_set = pe_o;
	fix_set = fix_o;
	logSineSqrGamma2_set = logSineSqrGamma2_o;
	ke_set = ke_o;

	// Set the final total energy
	etot_set = pe_set + fix_set + ke_o + logSineSqrGamma2_set;

}

const int HMCSampler::getDistortOpt(void)
{
	return this->DistortOpt;
}

void HMCSampler::setDistortOption(const int& distortOptArg)
{
	this->DistortOpt = distortOptArg;
}

const SimTK::Real& HMCSampler::getBendStretchStdevScaleFactor(void)
{
	return this->QScaleFactor;
}

void HMCSampler::setBendStretchStdevScaleFactor(const SimTK::Real& s)
{
	this->QScaleFactor = s;
}

/* // TODO revise param1 and param2
SimTK::Real& HMCSampler::convoluteVariable(SimTK::Real& var,
		std::string distrib, SimTK::Real param1, SimTK::Real param2)
{
	// Bernoulli trial between the var and its inverse
	if(distrib == "BernoulliInverse"){
		SimTK::Real randomNumber_Unif;
		int randomSign;

		randomNumber_Unif = uniformRealDistribution(randomEngine);
		randomSign = int(std::floor(randomNumber_Unif * 2.0) - 1.0);
		if(randomSign < 0){
			var = 1.0 / var;
		}
	}
	// Bernoulli trial between the var and its reciprocal
	else if(distrib == "BernoulliReciprocal"){
		SimTK::Real randomNumber_Unif;
		int randomSign;

		randomNumber_Unif = uniformRealDistribution(randomEngine);
		randomSign = int(std::floor(randomNumber_Unif * 2.0) - 1.0);
		if(randomSign < 0){
			var = -1.0 * var;
		}
	}
	// Run Gaussian distribution
	else if(distrib == "normal"){

		var = var + (gaurand(randomEngine) * param1);
	}
	// Run truncated Gaussian distribution
	else if(distrib == "truncNormal"){

		SimTK::Real mean = var;
		var = var + (gaurand(randomEngine) * param1);

		if(var <= 0){
			var = mean;
		}
		if(var >= (2*mean)){
			var = mean ;
		}
	}
	// Run bimodal Gaussian distribution
	else if(distrib == "bimodalNormal"){

		var = convoluteVariable(var, "BernoulliInverse");
		var = var + gaurand(randomEngine);

	}
	// Gamma distribution
	else if(distrib == "gamma"){
		var = gammarand(randomEngine);
	}

	return var;

}

SimTK::Real HMCSampler::calcDeformationPotential(
		SimTK::Real& var,
		std::string distrib,
		SimTK::Real param1, SimTK::Real param2
)
{
	SimTK::Real retVal = 0.0;

	// Bernoulli trial between the var and its inverse
	if(distrib == "BernoulliInverse"){
		return 0;
	}
	// Bernoulli trial between the var and its reciprocal
	else if(distrib == "BernoulliReciprocal"){
		return 0;
	}
	// Run Gaussian distribution
	else if(distrib == "normal"){
		return 0;
	}
	// Run truncated Gaussian distribution
	else if(distrib == "truncNormal"){
		SimTK::Real var2 = var*var;

		SimTK::Real firstTerm = 0.0;
		SimTK::Real secondTerm = 0.0;

		firstTerm = (1.0 / (var2)) - var2;
		secondTerm = 2.0 * param1 * (var - (1.0/var));

		SimTK::Real numerator = firstTerm + secondTerm;
		SimTK::Real denominator = param2*param2;

		retVal = (1.0 / this->beta) * 0.5 * (numerator / denominator);

	}
	// Run bimodal Gaussian distribution
	else if(distrib == "bimodalNormal"){
		return 0;
	}
	// Gamma distribution
	else if(distrib == "gamma"){
		return 0;
	}

	return retVal;
	
} */


/*
 * Shift all the generalized coordinates to scale bonds and angles
 * through BendStretch joint
 */
void HMCSampler::setQToScaleBendStretch(SimTK::State& someState,
	std::vector<SimTK::Real>& scaleFactors)
{
	// Scaling factor is set by Context only in the begining
	std::cout << "shiftQ Got " << this->QScaleFactor << " scale factor ";
	std::cout << "and turned it into " << this->QScaleFactor << "\n";

	// Set the return scalingFactors to QScaleFactor
	scaleFactors.resize(world->acosX_PF00.size() + world->normX_BMp.size(),
		1.0);

	// Ground and first mobod don't have internal coordinates
	int mbxOffset = 2;
	int sfIxOffset = world->normX_BMp.size();
	SimTK::Real sf = 1.0;

	for (SimTK::MobilizedBodyIndex mbx(mbxOffset);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobilized body
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);	

		// Set Q for angle
		sf = this->QScaleFactor - 1.0;
		scaleFactors[ sfIxOffset + (int(mbx) - 1) ] = this->QScaleFactor;

		mobod.setOneQ(someState, 0,
			-1.0 * (sf * world->acosX_PF00[int(mbx) - 1])  );

		// Set Q for bond
		sf = this->QScaleFactor - 1.0;
		scaleFactors[ (int(mbx) - 1) ] = this->QScaleFactor;
		mobod.setOneQ(someState, 1,
			sf * world->normX_BMp[int(mbx) - 1]  );

	}

	// Save changes by advancing to Position Stage
	system->realize(someState, SimTK::Stage::Position);

	// Test
	std::cout << "shifted Q = " << someState.getQ() << std::endl;
}


/**
 * Calculate bond length and angle deviations from their means
*/ 
/* void HMCSampler::calcBendStretchDeviations(
	SimTK::State& someState,
	std::vector<SimTK::Real>& X_PFdiffs,
	std::vector<SimTK::Real>& X_BMdiffs
)
{

	// Make sure it has 
	X_PFdiffs.resize(world->acosX_PF00_means.size(), 0.0);
	X_BMdiffs.resize(world->normX_BMp_means.size(), 0.0);

	// 
	for(unsigned int k = 0; k < X_PFdiffs.size(); k++){
		X_PFdiffs[k] = world->acosX_PF00[k] - world->acosX_PF00_means[k];
	}
	for(unsigned int k = 0; k < X_BMdiffs.size(); k++){
		X_BMdiffs[k] = world->normX_BMp[k] - world->normX_BMp_means[k];
	}

} */

/*
 * Shift all the generalized coordinates to scale bonds and angles
 * standard deviations through BendStretch joint
 */
void HMCSampler::setQToShiftBendStretchStdev(SimTK::State& someState,
std::vector<SimTK::Real>& scaleFactors)
{
	/* 	world->PrintX_PFs();
	world->PrintX_BMs(); */
		// Print the scale factor
	std::cout << "shiftQ Got " << this->QScaleFactor << " scale factor "
		<< std::endl;

	// Prepare scaling factors for return
	scaleFactors.resize(world->acosX_PF00.size() + world->normX_BMp.size(),
		1.0);
		
	// 1. Get the deviations from their means
	std::vector<SimTK::Real> X_PFdiffs;
	std::vector<SimTK::Real> X_BMdiffs;
	world->calcBendStretchDeviations(someState, X_PFdiffs, X_BMdiffs);

	// Scale the differences with QScale. -1 is only here because Q is always 0
	int k = -1;
	for(auto& diff : X_PFdiffs){
		/* diff += QScaleFactor; */
		//std::cout << "diff= " << diff << std::endl;
	}
	for(auto& diff : X_BMdiffs){
		/* diff += QScaleFactor; */
		//std::cout << "diff= " << diff << std::endl;
	}

	// Ground and first mobod don't have internal coordinates
	int mbxOffset = 2;
	int sfIxOffset = world->normX_BMp.size();

	for (SimTK::MobilizedBodyIndex mbx(mbxOffset);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobilized body
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		
		// We only allocated  X_PFs for non-Ground bodies
		//RESTORE RESTORE RESTORE Check X_FM assignment 0

		// Get the scaleFactors too
		if(std::abs(world->normX_BMp[(int(mbx) - 1)]) > 0.00000001){

			mobod.setOneQ(someState, 1, this->QScaleFactor);

			scaleFactors[ (int(mbx) - 1) ] = 
			(world->normX_BMp[int(mbx) - 1] + this->QScaleFactor) /
				world->normX_BMp[int(mbx) - 1];
		}		

		if(std::abs(world->acosX_PF00[int(mbx) - 1]) > 0.00000001){

			mobod.setOneQ(someState, 0, -1.0 * this->QScaleFactor);
			
			scaleFactors[ sfIxOffset + (int(mbx) - 1) ] = 
			(world->acosX_PF00[int(mbx) - 1] + (-1.0 * this->QScaleFactor)) /
				world->acosX_PF00[int(mbx) - 1];
		}	

	}

	// Save changes by advancing to Position Stage
	system->realize(someState, SimTK::Stage::Position);

	// Test
	std::cout << "shifted Q = " << someState.getQ() << std::endl;
}

/**
 * Shift all the generalized coordinates to scale bonds and angles
 * standard deviations through BendStretch joint
*/
void HMCSampler::setQToScaleBendStretchStdev(SimTK::State& someState,
std::vector<SimTK::Real>& scaleFactors)
{
	// Print the scale factor
	std::cout << "shiftQ Got " << this->QScaleFactor << " scale factor "
		<< std::endl;

	//world->traceBendStretch(someState);
	//world->PrintAcosX_PFs();
	//world->PrintNormX_BMs();
	//world->PrintAcosX_PFMeans();
	//world->PrintNormX_BMMeans();

	// Resize scaling factors
	scaleFactors.resize(world->acosX_PF00.size() + world->normX_BMp.size(),
		1.0);
		
	// 1. Get the deviations from their means
	std::vector<SimTK::Real> Q_of_X_PFdiffs;
	std::vector<SimTK::Real> Q_of_X_BMdiffs;
	world->calcBendStretchDeviations(someState, Q_of_X_PFdiffs, Q_of_X_BMdiffs);

	// Scale differences with QScale-1. (Solve s(X-miu) + miu = Q + X)
	int k = -1;
	for(auto& diff : Q_of_X_PFdiffs){
		//std::cout << "Q(X_PFdiff)= " << diff << std::endl;
		
		diff *= (QScaleFactor - 1.0);
		//diff *= -1.0;
		
		if(world->visual){
			world->paraMolecularDecorator->updFCommVars(diff);
		}
	}
	for(auto& diff : Q_of_X_BMdiffs){
		//std::cout << "Q(X_BMdiff)= " << diff << std::endl;
		
		diff *= (QScaleFactor - 1.0);
		//diff *= -1.0;
		
		if(world->visual){
			world->paraMolecularDecorator->updBCommVars(diff);
		}
	}

	// Now set Qs to the scaled differences
	int mbxOffset = 2; // no internal coordinates for Ground and first mobod 
	int sfIxOffset = world->normX_BMp.size(); 

	for (SimTK::MobilizedBodyIndex mbx(mbxOffset);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobilized body
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		
		// Previous mbx
		int sfIx = (int(mbx) - 1);

		// Set the stretch 
		if(std::abs(world->normX_BMp[ sfIx ]) > 0.00000001){

			// Set Q to scaled diff
			if(mobod.getNumQ(someState) == 2){
				mobod.setOneQ(someState, 1, (+1.0) * Q_of_X_BMdiffs[sfIx]);
				//std::cout << "setQToScaleBendStretchStdev mobod.setOneQ 0\n";
			}else{
				mobod.setOneQ(someState, 0, (+1.0) * Q_of_X_BMdiffs[sfIx]);
			}

			// Record the scaleFactors too: (X + Q) / X
			scaleFactors[ sfIx ] = 
			(world->normX_BMp[ sfIx ] + ((+1.0) * Q_of_X_BMdiffs[ sfIx ])) /
				world->normX_BMp[ sfIx ];
		}

		// Set the bend rotation
		if(std::abs(world->acosX_PF00[ sfIx ]) > 0.00000001){

			// Set Q to scaled diff
			if(mobod.getNumQ(someState) == 2){
				mobod.setOneQ(someState, 0, (-1.0) * Q_of_X_PFdiffs[ sfIx ]);
			
				// Record the scaleFactors too: (X + Q) / X
				scaleFactors[ sfIxOffset + sfIx ] = 
				(world->acosX_PF00[ sfIx ] + ((+1.0) * Q_of_X_PFdiffs[ sfIx ])) /
					world->acosX_PF00[ sfIx ];
			}else{
				scaleFactors[ sfIxOffset + sfIx ] = 1.0;
			}

		}		
	}

	//std::cout << "scaleFactors: "; PrintCppVector(scaleFactors);

	// Save changes by advancing to Position Stage
	system->realize(someState, SimTK::Stage::Position);
	if(world->visual){
		std::this_thread::sleep_for(std::chrono::milliseconds(4000));
	}

	// Test
	std::cout << "shifted Q = " << someState.getQ() << std::endl;

}


/*
 * Shift all the generalized coordinates to scale bonds and angles
 * standard deviations through BendStretch joint
 */
void HMCSampler::setQToScaleBendStretchStdev_Old(SimTK::State& someState,
std::vector<SimTK::Real>& scaleFactors)
{
	// Scaling factor is set by Context only in the begining
	//this->QScaleFactor = 1.0;

	std::cout << "shiftQ Got " << this->QScaleFactor << " scale factor ";

	//convoluteVariable(this->QScaleFactor, "BernoulliInverse");
	//convoluteVariable(this->QScaleFactor, "truncNormal",
	//	0.1);
	
	std::cout << "and turned it into " << this->QScaleFactor << "\n";

	// Return the scaling factors
	scaleFactors.resize(world->acosX_PF00.size() + world->normX_BMp.size(),
		1.0);

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

	// 3. Scale the differences with QScale. -1 is only here because Q is always 0
	int k = -1;
	for(auto& diff : X_PFdiffs){
		diff *= QScaleFactor - 1.0;
		//std::cout << "diff= " << diff << std::endl;
	}
	for(auto& diff : X_BMdiffs){
		diff *= QScaleFactor - 1.0;
		//std::cout << "diff= " << diff << std::endl;
	}

	// Print the differences	
	/* for(unsigned int k = 0; k < X_PFdiffs.size(); k++){
		std::cout << "Excel X_PF " << k << " "
			<< world->acosX_PF00[k] << " "
			<< world->acosX_PF00_means[k] << " " 
			<< X_PFdiffs[k] << std::endl;
	}
	for(unsigned int k = 0; k < X_PFdiffs.size(); k++){
		std::cout << "Excel X_BM " << k << " "
			<< world->normX_BMp[k]  << " "
			<< world->normX_BMp_means[k] << " "
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
		mobod.setOneQ(someState, 1, X_BMdiffs[int(mbx) - 1]);

		// Get the scaleFactors too
		if(std::abs(world->normX_BMp[(int(mbx) - 1)]) > 0.00000001){
			scaleFactors[ (int(mbx) - 1) ] = 
			(world->normX_BMp[int(mbx) - 1] + X_BMdiffs[int(mbx) - 1]) /
				world->normX_BMp[int(mbx) - 1];

		/* std::cout << "i bm diff c "
			<< int(mbx) - 1 << " " 
			<< world->normX_BMp[int(mbx) - 1] << " "
			<< X_BMdiffs[int(mbx) - 1] << " "
			<< scaleFactors[ (int(mbx) - 1) ] << " "
			<< std::endl; */
		}		

		if(std::abs(world->acosX_PF00[int(mbx) - 1]) > 0.00000001){
			scaleFactors[ world->normX_BMp.size() + (int(mbx) - 1) ] = 
			(world->acosX_PF00[int(mbx) - 1] + (-1.0 * X_PFdiffs[int(mbx) - 1])) /
				world->acosX_PF00[int(mbx) - 1];
		}
		
		/* std::cout << "i pf diff c "
			<< world->normX_BMp.size() + (int(mbx) - 1) << " " 
			<< world->acosX_PF00[int(mbx) - 1] << " "
			<< X_PFdiffs[int(mbx) - 1] << " "
			<< scaleFactors[ (int(mbx) - 1) ] << " "
			<< std::endl; */

		// Set Q to a uniform distribution
		/*randomNumber_Unif = uniformRealDistribution(randomEngine);
		randomNumber_Unif = (randomNumber_Unif * 2.0) - 1.0;
		mobod.setOneQ(someState, 0, randomNumber_Unif * 0.0005);

		randomNumber_Unif = uniformRealDistribution(randomEngine);
		randomNumber_Unif = (randomNumber_Unif * 2.0) - 1.0;
		mobod.setOneQ(someState, 1, randomNumber_Unif * 0.0005); */

	}

	// Save changes by advancing to Position Stage
	system->realize(someState, SimTK::Stage::Position);

	// Test
	std::cout << "shifted Q = " << someState.getQ() << std::endl;

/* 	for(unsigned int k = 0; k < world->normX_BMp.size(); k++){
		if(world->normX_BMp[k] != 0.0){
			scaleFactors[k] = this->QScaleFactor;
			scaleFactors[k] *= 
				1.0 - (world->normX_BMp_means[k] / world->normX_BMp[k]);

			scaleFactors[k] +=
				(world->normX_BMp_means[k] / world->normX_BMp[k]);
		}
	}

	for(unsigned int k = 0; k < world->acosX_PF00.size(); k++){
		scaleFactors[world->normX_BMp.size() + k] = this->QScaleFactor;
		scaleFactors[world->normX_BMp.size() + k] *=
			1.0 - (world->acosX_PF00_means[k] / world->acosX_PF00[k]);

		scaleFactors[world->normX_BMp.size() + k] +=
			(world->acosX_PF00_means[k] / world->acosX_PF00[k]);
	} */

/* 	int nu = someState.getNU();

	Matrix mathJ;
	calcMathJacobian(someState, mathJ);
	PrintBigMat(mathJ, mathJ.nrow(), mathJ.ncol(), 2, "mathJacobian");

	SimTK::Matrix mathJtJ;
	mathJtJ = mathJ.transpose() * mathJ;

	PrintBigMat(mathJtJ, mathJtJ.nrow(), mathJtJ.ncol(), 2, "mathJacobianSquared");


	std::vector<SimTK::Real> tempM(mathJtJ.nrow() * mathJtJ.ncol());
	for(int i=0; i<mathJtJ.nrow(); i++){
		for(int j=0; j<mathJtJ.ncol(); j++){
			tempM[i * nu + j] = mathJtJ(i, j);
		}
	}

	SimTK::Real detMathJ = cstyle_det(&tempM[0], nu);
	std::cout << std::setprecision(10) << std::fixed;
	std::cout << "mathDeterminant " << detMathJ << std::endl;

	//SimTK::Lapack::getrf

	// Get System Jacobian
	SimTK::Matrix J_G;
	matter->calcSystemJacobian(someState, J_G);
	//PrintBigMat(J_G, J_G.nrow(), J_G.ncol(), 2, "systemJacobian");

	// Get Cartesian mass matrix
	SimTK::Matrix cartM;
	getCartesianMassMatrix(someState, cartM);
	PrintBigMat(cartM, cartM.nrow(), cartM.ncol(), 2, "CartesianMassMatrix");

	// Compare metric tensor
	SimTK::Matrix JtJ;
	JtJ = J_G.transpose() * J_G;
	PrintBigMat(JtJ, JtJ.nrow(), JtJ.ncol(), 2, "sysJacobianSquared");

	std::vector<SimTK::Real> EiM(nu * nu);
	for(int i=0; i<nu; i++){
		for(int j=0; j<nu; j++){
			EiM[i * nu + j] = JtJ(i, j);
		}
	}
	SimTK::Real detJ = cstyle_det(&EiM[0], nu);
	std::cout << "sysDeterminant " << detJ << std::endl; */

/* 	// Get System Jacobian
	SimTK::Matrix J_G;
	matter->calcSystemJacobian(someState, J_G);
	std::cout << "J_G\n" << J_G << std::endl;

	SimTK::Array_<SimTK::SpatialInertia, SimTK::MobilizedBodyIndex> R;
	const SimTK::ArticulatedInertia A;


	matter->calcCompositeBodyInertias(someState, R);
	
	std::cout << "Mass properties " << std::endl;
	int i = -1;
	for (SimTK::MobilizedBodyIndex mbx(1);
	mbx < matter->getNumBodies();
	++mbx){
		i += 1;
		PrintSpatialMat(R[mbx].toSpatialMat(),
			3, "Composite Body Inertia");

		const SimTK::ArticulatedInertia 
			A(matter->getArticulatedBodyInertia(someState, mbx));

		PrintSpatialMat(A.toSpatialMat(),
			3, "Articulated Body Inertia");

	} */

/* 	// Get mathematical Jacobian square determinant
	SimTK::Matrix mathJtJ;
	mathJtJ = mathJ.transpose() * mathJ;
	std::cout << "mathJtJ\n" << mathJtJ << std::endl << std::flush; */
	
	/* // linker problems ...
	SimTK::Eigen mathJtJEigen(mathJtJ);
	SimTK::Vector mathJtJEigenVals;
	mathJtJEigen.getAllEigenValues(mathJtJEigenVals);
	std::cout << "mathJtJEigenVals\n" << mathJtJEigenVals << std::endl; */

/* 	// Get mathematical Jacobian square determinant
	SimTK::SymMat<28> smMathJtJ; // BUG: this should be a constant
	for(unsigned int i = 0; i < ndofs; i++){
		for (unsigned int j = 0; j < ndofs; j++){
			smMathJtJ(i, j) = mathJtJ(i, j);
			smMathJtJ(j, i) = mathJtJ(i, j);
		}
	}
	SimTK::Real detMathJtJ = SimTK::det(smMathJtJ);
	std::cout << "mathJtJ determinant\n" << detMathJtJ << std::endl; */

}

/**
 * Get the log of the Jacobian of a bond-angle strtch
*/
SimTK::Real 
HMCSampler::calcBendStretchJacobianDetLog(SimTK::State& someState,
	std::vector<SimTK::Real> scaleFactors)
{

	int x_pf_k = 0;

	// Get log of the Cartesian->BAT Jacobian
	SimTK::Real logJacBAT_0 = 0.0; // log space
	//SimTK::Real logJacBAT_0 = 1.0; // normal space
	for(unsigned int k = 0; k < world->normX_BMp.size(); k++){

			// scaleFactors index is shifted
			x_pf_k = k + world->normX_BMp.size();

			// Checks
			/* std::cout << " k acosPF normBM sacosPF snormBM  " << k << " "
				<< world->acosX_PF00[k] << " " << world->normX_BMp[k] << " "
					<< world->acosX_PF00[k] * scaleFactors[x_pf_k] << " "
					<< world->normX_BMp[k] * scaleFactors[k] << " "
				<< std::endl << std::flush; */

			// Get bond term
			SimTK::Real d2 = world->normX_BMp[k];
			if(d2 != 0){
					d2 = 2.0 * std::log(d2); // log space
					//d2 = d2 * d2; // normal space
			}

			// Get the angle term
			SimTK::Real sineTheta = std::abs(world->acosX_PF00[k]);

			if(sineTheta != 0){
					sineTheta = std::log(std::sin(sineTheta)); // log space
					//sineTheta = std::sin(sineTheta); // normal space
			}

			// Accumulate result
			logJacBAT_0 += (d2 + sineTheta); // log space
			//logJacBAT_0 *= (d2 * sineTheta); // normal space
	}

	// Get log of the scaling Jacobian
	// Although the method seems stupid for now, it has generality
	// and allows insertion of new code
	SimTK::Real logJacScale = 0.0; // log space
	//SimTK::Real logJacScale = 1.0; // normal space
	for(unsigned int k = 0; k < scaleFactors.size(); k++){

			// Accumulate result
			if(scaleFactors[k] != 0){
					logJacScale += std::log(std::abs(scaleFactors[k])); // log space
					//logJacScale *= std::abs(scaleFactors[k]); // normal space
			}
	}

	// Get log of the BAT->Cartesian Jacobian after scaling
	SimTK::Real logJacBAT_tau = 0.0; // log space
	//SimTK::Real logJacBAT_tau = 1.0; // normal space
	for(unsigned int k = 0; k < world->normX_BMp.size(); k++){

			// scaleFactors index is shifted
			x_pf_k = k + world->normX_BMp.size();

			// Checks
			/* std::cout << " k acosPF normBM sacosPF snormBM  " << k << " "
				<< world->acosX_PF00[k] << " " << world->normX_BMp[k] << " "
					<< world->acosX_PF00[k] * scaleFactors[x_pf_k] << " "
					<< world->normX_BMp[k] * scaleFactors[k] << " "
				<< std::endl << std::flush; */

			// Get bond term
			SimTK::Real d2 = (world->normX_BMp[k] * scaleFactors[k]);

			if(d2 != 0){
					d2 = 2.0 * std::log(d2); // log space
					//d2 = d2 * d2; // normal space
			}

			// Get the angle term
			SimTK::Real sineTheta = std::abs(world->acosX_PF00[k] * scaleFactors[x_pf_k]);

			if(std::abs(sineTheta) > 0.0000001){
					sineTheta = std::log(std::sin(sineTheta)); // log space
					//sineTheta = std::sin(sineTheta); // normal space
  			}

			// Accumulate result
			logJacBAT_tau += (d2 + sineTheta); // log space
			//logJacBAT_tau *= (d2 * sineTheta); // normal space
	}

        // Final result
        //SimTK::Real logBendStretchJac = (-1.0 * logJacBAT_0) + logJacScale + logJacBAT_tau;
        //SimTK::Real logBendStretchJac = logJacScale; // 7000
        SimTK::Real logBendStretchJac = logJacBAT_0 + logJacScale + (-1.0 * logJacBAT_tau); // 8000
        //SimTK::Real logBendStretchJac = (-1.0 * logJacBAT_0) + (1.0 * logJacBAT_tau); // 9000

        std::cout << "logJacBAT " << logJacBAT_0
                << " logJacScale " << logJacScale
                << " logJacBATInv " << logJacBAT_tau
                << " logBendStretchJac " << logBendStretchJac
                << std::endl;

	return logBendStretchJac;

}

/**
 * It implements a non-equilibrium proposal move, 
 * then propagates the trajectory.
*/
bool HMCSampler::proposeNEHMC(SimTK::State& someState)
{

	/*// Store old configuration
	storeOldConfigurationAndPotentialEnergies(someState);*/

	// Resize the scale factors vector to be handed further
	std::vector<SimTK::Real> scaleFactors;
	scaleFactors.resize(world->acosX_PF00.size() + world->normX_BMp.size(),
		1.0);

	// Scale bonds and angles
	if(this->nofSamples >= 0){ // dont't take burn-in
		//setQToScaleBendStretch(someState, scaleFactors);
		setQToScaleBendStretchStdev(someState, scaleFactors);
	}

	// Get the log of the Jacobian of the change
	this->bendStretchJacobianDetLog =
		calcBendStretchJacobianDetLog(someState, scaleFactors);

	// Adapt timestep
	if(shouldAdaptTimestep){
		adaptTimestep(someState);
	}

	// Initialize velocities from Maxwell-Boltzmann distribution
	initializeVelocities(someState);

	// Store the proposed energies
	calcProposedKineticAndTotalEnergy(someState);

	// Adapt timestep
	bool shouldAdaptWorldBlocks = false;
	if(shouldAdaptWorldBlocks){
		adaptWorldBlocks(someState);
	}

	// Apply the L operator
	if(this->integratorName == IntegratorName::VERLET){
		integrateTrajectory(someState);
		//integrateTrajectoryOneStepAtATime(someState);
	}else{
		std::cout << "ProposeNEHMC: NON-VERLET integrator\n";
		system->realize(someState, SimTK::Stage::Dynamics);
	}

	calcNewConfigurationAndEnergies(someState);

	return validateProposal();

}

// Set Sphere Radius when doing RANDOM_WALK
void HMCSampler::setSphereRadius(float argSphereRadius)
{
	//std::cout << "HMCSampler::setSphereRadius:  Radius Set: " << sphereRadius << std::endl;	
	sphereRadius = argSphereRadius;
	std::cout << "HMCSampler::setSphereRadius:  Radius Set: " << sphereRadius << std::endl;
}

/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory. **/
bool HMCSampler::proposeEquilibrium(SimTK::State& someState)
{

/* 	// Store old configuration
	storeOldConfigurationAndPotentialEnergies(someState); */

	if(integratorName == IntegratorName::OMMVV){

		// // ELIZA: Check the code below
		OMM_setTemperature(this->boostT);
		OMM_calcProposedKineticAndTotalEnergy();
		OMM_integrateTrajectory(someState);
		OMM_calcNewConfigurationAndEnergies();

	}else if (integratorName == IntegratorName::VERLET){

		// Adapt timestep
		if(shouldAdaptTimestep){
			adaptTimestep(someState);
		}

		// Initialize velocities from Maxwell-Boltzmann distribution
		initializeVelocities(someState);

		// Store the proposed energies
		calcProposedKineticAndTotalEnergy(someState);

		// Adapt timestep
		bool shouldAdaptWorldBlocks = false;
		if(shouldAdaptWorldBlocks){
			adaptWorldBlocks(someState);
		}

		// Apply the L operator
		if(this->integratorName == IntegratorName::VERLET){
			integrateTrajectory(someState);
			//integrateTrajectoryOneStepAtATime(someState);
		}else{
			std::cout << "Propose: EMPTY integrator\n";
			system->realize(someState, SimTK::Stage::Dynamics);
		}

		calcNewConfigurationAndEnergies(someState);
	
	}else if (integratorName == IntegratorName::BOUND_WALK){

		// Hard-coded, need to be able to propagate from input
		// but for the moment this will have to do.
		
		// These are integers, since you'd never *really* need a long list
		// of bounders
		//const vector<vector<int>> bounders{{0},{2}};
		// These are strings, since you could  want all solvent
		// molecules to be bound, so the "s" option is implemented
		//const vector<vector<string>> boundeds{{"2"},{"s"}};
		// These are strings, since you could want the 
		//const vector<vector<string>> bounder_subset{{"74","88"},{"a"}};
		//const vector<float> bound_radii;




		// Set velocities to zero
		setVelocitiesToZero(someState);
		system->realize(someState, SimTK::Stage::Velocity);

		// Store the proposed energies
		calcProposedKineticAndTotalEnergy(someState);

		std::cout << "Propose: BOUND_WALK integrator" << std::endl;
		if(topologies.size() < 2){
			std::cout << "BOUND integrators should only be used over many molecules\n";
		}

		// Get the binding site center

		SimTK::Vec3 geometricCenter = world->getGeometricCenterOfSelection(someState);

		// We *assume* the last molecule is the ligand.
		const int nOfBodies = matter->getNumBodies();

		// Print the geometric center (for debugging purposes)
		std::cout << "HMCSampler Binding Site Center: \t" << geometricCenter << std::endl;
		std::cout << "HMCSampler Binding Site Sphere Radius: \t" << sphereRadius << " nm\n";	

		// Ligand
		const SimTK::MobilizedBody& mobod_L = matter->getMobilizedBody(
			SimTK::MobilizedBodyIndex(2) ); 
		const Vec3 COM_L = mobod_L.getBodyMassCenterStation(someState);
		const Vec3 COM_G = mobod_L.findMassCenterLocationInGround(someState);
		const Transform& X_FM = mobod_L.getMobilizerTransform(someState);
		const Transform& X_PF = mobod_L.getInboardFrame(someState);
		const Transform& X_GL = mobod_L.getBodyTransform(someState);

		// Water (Before rotation)
		const SimTK::MobilizedBody& mobod_W = matter->getMobilizedBody(
				SimTK::MobilizedBodyIndex(3) ); 
		const Vec3 COM_W = mobod_W.getBodyMassCenterStation(someState);
		const Transform& X_GW = mobod_W.getBodyTransform(someState);
		const Vec3 COM_LW = mobod_W.findMassCenterLocationInAnotherBody(someState, mobod_L); // Expressed in W (??)
		const Vec3 COM_GW = mobod_W.findMassCenterLocationInGround(someState);
		const Transform& X_PF_Wat = mobod_W.getInboardFrame(someState);


		std::cout << "COM_LW: " << COM_LW << " (" << COM_LW.norm() << ")\n" ;

		// It's best to first compute the rotation, then
		// compute the translation to fit the "New" X_FM
		
		// Generate random rotation quaternion
		double x,y,z, u,v,w, s;
		do {
			x = uniformRealDistribution_m1_1(randomEngine);
			y = uniformRealDistribution_m1_1(randomEngine);
			z = x*x + y*y;
		}while (z > 1);
		do {
			u = uniformRealDistribution_m1_1(randomEngine);
			v = uniformRealDistribution_m1_1(randomEngine);
			w = u*u + v*v;
		} while (w > 1);

		s = sqrt((1-z) / w);

		SimTK::Quaternion randQuat(x, y, s*u, s*v);
		SimTK::Quaternion currQuat(someState.updQ()[0], someState.updQ()[1],
			someState.updQ()[2], someState.updQ()[3]);
		SimTK::Quaternion resultQuat = multiplyQuaternions(randQuat, currQuat);

		mobod_L.setQToFitRotation(someState, SimTK::Rotation(resultQuat));
		//mobod_W.setQToFitRotation(someState, SimTK::Rotation(resultQuat));
		system->realize(someState, SimTK::Stage::Dynamics);

		// Translate
		SimTK::Vec3 New_LW = SimTK::Rotation(resultQuat) * COM_LW ; //Still expressed in L
		std::cout << "New_LW: " << New_LW << " (" << New_LW.norm() << ")\n" ;
		New_LW = X_GL.shiftBaseStationToFrame(New_LW); // Expressed in Ground
		New_LW = X_GW.shiftBaseStationToFrame(New_LW); // Expressed in Water
		std::cout << "New_LW (expressed in W): " << New_LW << " (" << New_LW.norm() << ")\n\n" ;
		//mobod_W.setQToFitTranslation(someState, New_LW);
		//system->realize(someState, SimTK::Stage::Dynamics);

		// Determine the distance between center of mass of ligand and geometricCenter
		SimTK::Vec3 ligandToSite = geometricCenter - COM_G;
		std::cout << "ligandToSite: " << ligandToSite << " ligandToSite Norm: " << ligandToSite.norm() << std::endl;

		// If ligand is too far reposition on a sphere centered on the last body
		// center of mass (also add a margin of error, so you don't stay too much on the surface
		// of the sphere)
		//if(ligandToSite.norm() > (sphereRadius)){
		if (0){

			SimTK::Vec3 BR={0,0,0};

			// Sample a random vector centered in 0 and expressed in G
			SimTK::Real theta = uniformRealDistribution_0_2pi(randomEngine);
			SimTK::Real phi = std::acos(2.0 * uniformRealDistribution(randomEngine) - 1.0);
			SimTK::Vec3 randVec={0,0,0};


			randVec[0] = sphereRadius * std::cos(theta) * std::sin(phi);
			randVec[1] = sphereRadius * std::sin(theta) * std::sin(phi);
			randVec[2] = sphereRadius * std::cos(phi);
			randVec *= uniformRealDistribution(randomEngine);

			// Move it in BindingSiteCenter (BS)
			SimTK::Vec3 GR;
			GR = randVec + (geometricCenter - X_PF.p());

			BR = mobod_L.expressGroundVectorInBodyFrame(someState, GR);
			
			// Account for COM_L
			BR = BR - COM_L;

			// Express BR in F
			BR = X_FM.R()*BR;
			mobod_L.setQToFitTranslation(someState, BR);

		// If not, just generate a random translation
		}else{
			someState.updQ()[4] += uniformRealDistribution_m1_1(randomEngine) * (1 * this->timestep);
			someState.updQ()[5] += uniformRealDistribution_m1_1(randomEngine) * (1 * this->timestep);
			someState.updQ()[6] += uniformRealDistribution_m1_1(randomEngine) * (1 * this->timestep);
		}

		//system->realize(someState, SimTK::Stage::Dynamics);


		/* // Water
		for (int waterIx = 3; waterIx < nOfBodies; waterIx++){ 
			const SimTK::MobilizedBody& mobod_W = matter->getMobilizedBody(
				SimTK::MobilizedBodyIndex(waterIx) ); 
			const Vec3 COM_W = mobod_W.getBodyMassCenterStation(someState);
			const Vec3 COM_GW = mobod_W.findMassCenterLocationInGround(someState);
			const Transform& X_FM_Water = mobod_W.getMobilizerTransform(someState);
			const Transform& X_PF_Water = mobod_W.getInboardFrame(someState);

			// Get distance between water and ligand
					
			const Vec3 WL = mobod_W.findStationLocationInAnotherBody(someState, COM_W, mobod_L);
			std::cout << "COM_L: " << COM_G << " COM_GW: " << COM_GW << std::endl;
			std::cout << "WL: " << WL << " WL_NORM: " << WL.norm() << std::endl;

			// If water is too far from ligand reposition on a sphere centered on the 
			// ligand center of mass
			float waterSphere = 1;
			if(WL.norm() > waterSphere){
			//if (1){

				std::cout << "Water (" << waterIx-1 << ") too far from ligand (" << WL.norm() << "), repositioning...\n";
				SimTK::Vec3 BR={0,0,0};

				// Sample a random vector centered in 0 and expressed in G
				SimTK::Real theta = uniformRealDistribution_0_2pi(randomEngine);
				SimTK::Real phi = std::acos(2.0 * uniformRealDistribution(randomEngine) - 1.0);
				SimTK::Vec3 randVec={0,0,0};

				randVec[0] = waterSphere * std::cos(theta) * std::sin(phi);
				randVec[1] = waterSphere * std::sin(theta) * std::sin(phi);
				randVec[2] = waterSphere * std::cos(phi);
				//randVec *= uniformRealDistribution(randomEngine);

				// Move it in BindingSiteCenter (BS)
				SimTK::Vec3 GR;
				//GR = randVec + (geometricCenter - X_PF.p());
				GR = randVec + (COM_G_PostTP - X_PF_Water.p());

				BR = mobod_W.expressGroundVectorInBodyFrame(someState, GR);
				
				// Account for COM_L
				BR = BR - COM_W;

				// Express BR in F
				BR = X_FM_Water.R()*BR;
				mobod_W.setQToFitTranslation(someState, BR);
				system->realize(someState, SimTK::Stage::Position);
				std::cout << "Water (" << waterIx-1
				          <<") repositioned in " << mobod_W.expressVectorInGroundFrame(someState, BR) 
						  << std::endl << std::endl;
				system->realize(someState, SimTK::Stage::Dynamics);

			}
		}
 */
		system->realize(someState, SimTK::Stage::Dynamics);
		calcNewConfigurationAndEnergies(someState);

	}else if (integratorName == IntegratorName::BOUND_HMC){

		// Set velocities to zero
		setVelocitiesToZero(someState);
		system->realize(someState, SimTK::Stage::Velocity);

		// Store the proposed energies
		calcProposedKineticAndTotalEnergy(someState);

		std::cout << "Propose: BOUND_HMC integrator" << std::endl;
		if(topologies.size() < 2){
			std::cout << "BOUND integrators should only be used over many molecules\n";
		}

		// Get the binding site center

		SimTK::Vec3 geometricCenter = world->getGeometricCenterOfSelection(someState);

		// We *assume* the last molecule is the ligand.
		const int nOfBodies = matter->getNumBodies();

		// Print the geometric center (for debugging purposes)
		std::cout << "HMCSampler Binding Site Center: \t" << geometricCenter << std::endl;
		std::cout << "HMCSampler Binding Site Sphere Radius: \t" << sphereRadius << " nm\n";	

		// Ligand
		const SimTK::MobilizedBody& mobod_L = matter->getMobilizedBody(
			SimTK::MobilizedBodyIndex(2) ); 
		const Vec3 COM_L = mobod_L.getBodyMassCenterStation(someState);
		const Vec3 COM_G = mobod_L.findMassCenterLocationInGround(someState);
		const Transform& X_FM = mobod_L.getMobilizerTransform(someState);
		const Transform& X_PF = mobod_L.getInboardFrame(someState);

		// Unlike "RANDOM_WALK", this integrator does not need to do 
		// random rotation, since we integrate the trajectory, which
		// includes rotation.

		// Determine the distance between center of mass of ligand and geometricCenter
		SimTK::Vec3 ligandToSite = geometricCenter - COM_G;
		std::cout << "ligandToSite(kick): " << ligandToSite << " ligandToSite Norm: " << ligandToSite.norm() << std::endl;

		// If ligand is too far reposition on a sphere centered on the last body
		// center of mass
		if(ligandToSite.norm() > sphereRadius){
		//if (1){

			std::cout << "Ligand too far from center (" << ligandToSite.norm() << "), repositioning...\n";
			SimTK::Vec3 BR={0,0,0};

			// Sample a random vector centered in 0 and expressed in G
			SimTK::Real theta = uniformRealDistribution_0_2pi(randomEngine);
			SimTK::Real phi = std::acos(2.0 * uniformRealDistribution(randomEngine) - 1.0);
			SimTK::Vec3 randVec={0,0,0};

			randVec[0] = sphereRadius * std::cos(theta) * std::sin(phi);
			randVec[1] = sphereRadius * std::sin(theta) * std::sin(phi);
			randVec[2] = sphereRadius * std::cos(phi);
			randVec *= uniformRealDistribution(randomEngine);

			// Move it in BindingSiteCenter (BS)
			SimTK::Vec3 GR;
			GR = randVec + (geometricCenter - X_PF.p());

			BR = mobod_L.expressGroundVectorInBodyFrame(someState, GR);
			
			// Account for COM_L
			BR = BR - COM_L;

			// Express BR in F
			BR = X_FM.R()*BR;
			mobod_L.setQToFitTranslation(someState, BR);

			system->realize(someState, SimTK::Stage::Dynamics);
		}

		// Water Part
		const SimTK::MobilizedBody& mobod_L_PostTP = matter->getMobilizedBody(
			SimTK::MobilizedBodyIndex(2) ); 
		const Vec3 COM_L_PostTP = mobod_L_PostTP.getBodyMassCenterStation(someState);
		const Vec3 COM_G_PostTP = mobod_L_PostTP.findMassCenterLocationInGround(someState);
		const Transform& X_FM_PostTP = mobod_L_PostTP.getMobilizerTransform(someState);
		const Transform& X_PF_PostTP = mobod_L_PostTP.getInboardFrame(someState);

		// Water
		for (int waterIx = 3; waterIx < nOfBodies; waterIx++){ 
			const SimTK::MobilizedBody& mobod_W = matter->getMobilizedBody(
				SimTK::MobilizedBodyIndex(waterIx) ); 
			const Vec3 COM_W = mobod_W.getBodyMassCenterStation(someState);
			const Vec3 COM_GW = mobod_W.findMassCenterLocationInGround(someState);
			const Transform& X_FM_Water = mobod_W.getMobilizerTransform(someState);
			const Transform& X_PF_Water = mobod_W.getInboardFrame(someState);

			// Get distance between water and ligand
					
			const Vec3 WL = mobod_W.findStationLocationInAnotherBody(someState, COM_W, mobod_L_PostTP);
			std::cout << "COM_L: " << COM_G << " COM_GW: " << COM_GW << std::endl;
			std::cout << "WL: " << WL << " WL_NORM: " << WL.norm() << std::endl;

			// If water is too far from ligand reposition on a sphere centered on the 
			// ligand center of mass
			float waterSphere = 1;
			if(WL.norm() > waterSphere){
			//if (1){

				std::cout << "Water (" << waterIx << ") too far from ligand (" << WL.norm() << "), repositioning...\n";
				SimTK::Vec3 BR={0,0,0};

				// Sample a random vector centered in 0 and expressed in G
				SimTK::Real theta = uniformRealDistribution_0_2pi(randomEngine);
				SimTK::Real phi = std::acos(2.0 * uniformRealDistribution(randomEngine) - 1.0);
				SimTK::Vec3 randVec={0,0,0};

				randVec[0] = waterSphere * std::cos(theta) * std::sin(phi);
				randVec[1] = waterSphere * std::sin(theta) * std::sin(phi);
				randVec[2] = waterSphere * std::cos(phi);
				//randVec *= uniformRealDistribution(randomEngine);

				// Move it in BindingSiteCenter (BS)
				SimTK::Vec3 GR;
				//GR = randVec + (geometricCenter - X_PF.p());
				GR = randVec + (COM_G_PostTP - X_PF_Water.p());

				BR = mobod_W.expressGroundVectorInBodyFrame(someState, GR);
				
				// Account for COM_L
				BR = BR - COM_W;

				// Express BR in F
				BR = X_FM_Water.R()*BR;
				mobod_W.setQToFitTranslation(someState, BR);
				system->realize(someState, SimTK::Stage::Position);
				std::cout << "Water (" << waterIx-1
				          <<") repositioned in " << mobod_W.expressVectorInGroundFrame(someState, BR) 
						  << std::endl << std::endl;

			}
		}
		system->realize(someState, SimTK::Stage::Dynamics);

		// Else, if not repositioned, integrate trajectory.
		if(ligandToSite.norm() <= sphereRadius) {
			initializeVelocities(someState);
			calcProposedKineticAndTotalEnergy(someState);

			integrateTrajectory(someState);
			system->realize(someState, SimTK::Stage::Dynamics);
		}

		calcNewConfigurationAndEnergies(someState);

	/* }else if (integratorName == IntegratorName::WATER_CAGE){

		// Set velocities to zero
		setVelocitiesToZero(someState);
		system->realize(someState, SimTK::Stage::Velocity);

		// Store the proposed energies
		calcProposedKineticAndTotalEnergy(someState);

		std::cout << "Propose: WATER_CAGE integrator" << std::endl;
		if(topologies.size() < 2){
			std::cout << "RANDOM_KICK integrators should only be used over many molecules\n";
		}

		// Get the binding site center

		// We *assume* the last molecule is the ligand.
		const int nOfBodies = matter->getNumBodies();

		SimTK::Vec3 geometricCenter = world->getGeometricCenterOfSelection(someState);

		// Print the geometric center (for debugging purposes)
		std::cout << "HMCSampler Binding Site Center: \t" << geometricCenter << std::endl;
		std::cout << "HMCSampler Binding Site Sphere Radius: \t" << sphereRadius << " nm\n"; 	

		// Ligand
		const SimTK::MobilizedBody& mobod_L = matter->getMobilizedBody(
			SimTK::MobilizedBodyIndex(2) ); 
		const Vec3 COM_L = mobod_L.getBodyMassCenterStation(someState);
		const Vec3 COM_G = mobod_L.findMassCenterLocationInGround(someState);
		const Transform& X_FM = mobod_L.getMobilizerTransform(someState);
		const Transform& X_PF = mobod_L.getInboardFrame(someState);
		
		// Unlike "RANDOM_WALK", this integrator does not need to do 
		// random rotation, since we integrate the trajectory, which
		// includes rotation.

		// Determine the distance between center of mass of ligand and geometricCenter
		SimTK::Vec3 ligandToSite = geometricCenter - COM_G;
		std::cout << "ligandToSite(kick): " << ligandToSite << " ligandToSite Norm: " << ligandToSite.norm() << std::endl;

		//if(ligandToSite.norm() > sphereRadius){
		if (1){

			std::cout << "Ligand too far from center (" << ligandToSite.norm() << "), repositioning...\n";
			SimTK::Vec3 BR={0,0,0};

			// Sample a random vector centered in 0 and expressed in G
			SimTK::Real theta = uniformRealDistribution_0_2pi(randomEngine);
			SimTK::Real phi = std::acos(2.0 * uniformRealDistribution(randomEngine) - 1.0);
			SimTK::Vec3 randVec={0,0,0};

			randVec[0] = sphereRadius * std::cos(theta) * std::sin(phi);
			randVec[1] = sphereRadius * std::sin(theta) * std::sin(phi);
			randVec[2] = sphereRadius * std::cos(phi);
			randVec *= uniformRealDistribution(randomEngine);

			// Move it in BindingSiteCenter (BS)
			SimTK::Vec3 GR;
			GR = randVec + (geometricCenter - X_PF.p());

			BR = mobod_L.expressGroundVectorInBodyFrame(someState, GR);
			
			// Account for COM_L
			BR = BR - COM_L;

			// Express BR in F
			BR = X_FM.R()*BR;
			mobod_L.setQToFitTranslation(someState, BR);

			system->realize(someState, SimTK::Stage::Dynamics);

			
		}

		const SimTK::MobilizedBody& mobod_L_PostTP = matter->getMobilizedBody(
			SimTK::MobilizedBodyIndex(2) ); 
		const Vec3 COM_L_PostTP = mobod_L_PostTP.getBodyMassCenterStation(someState);
		const Vec3 COM_G_PostTP = mobod_L_PostTP.findMassCenterLocationInGround(someState);
		const Transform& X_FM_PostTP = mobod_L_PostTP.getMobilizerTransform(someState);
		const Transform& X_PF_PostTP = mobod_L_PostTP.getInboardFrame(someState);

		// Water
		for (int waterIx = 1; waterIx < nOfBodies-2; waterIx++){ 
			const SimTK::MobilizedBody& mobod_W = matter->getMobilizedBody(
				SimTK::MobilizedBodyIndex(2+waterIx) ); 
			const Vec3 COM_W = mobod_W.getBodyMassCenterStation(someState);
			const Vec3 COM_GW = mobod_W.findMassCenterLocationInGround(someState);
			const Transform& X_FM_Water = mobod_W.getMobilizerTransform(someState);
			const Transform& X_PF_Water = mobod_W.getInboardFrame(someState);

			// Get distance between water and ligand
					
			const Vec3 WL = mobod_W.findStationLocationInAnotherBody(someState, COM_W, mobod_L);
			std::cout << "COM_L: " << COM_G << " COM_GW: " << COM_GW << std::endl;
			std::cout << "WL: " << WL << " WL_NORM: " << WL.norm() << std::endl;

			// If water is too far from ligand reposition on a sphere centered on the 
			// ligand center of mass
			float waterSphere = 2;
			if(WL.norm() > waterSphere){
			//if (1){

				std::cout << "Water (" << waterIx-1 << ") too far from ligand (" << WL.norm() << "), repositioning...\n";
				SimTK::Vec3 BR={0,0,0};

				// Sample a random vector centered in 0 and expressed in G
				SimTK::Real theta = uniformRealDistribution_0_2pi(randomEngine);
				SimTK::Real phi = std::acos(2.0 * uniformRealDistribution(randomEngine) - 1.0);
				SimTK::Vec3 randVec={0,0,0};

				randVec[0] = waterSphere * std::cos(theta) * std::sin(phi);
				randVec[1] = waterSphere * std::sin(theta) * std::sin(phi);
				randVec[2] = waterSphere * std::cos(phi);
				randVec *= uniformRealDistribution(randomEngine);

				// Move it in BindingSiteCenter (BS)
				SimTK::Vec3 GR;
				//GR = randVec + (geometricCenter - X_PF.p());
				GR = randVec + (COM_G_PostTP - X_PF_Water.p());

				BR = mobod_W.expressGroundVectorInBodyFrame(someState, GR);
				
				// Account for COM_L
				BR = BR - COM_W;

				// Express BR in F
				BR = X_FM_Water.R()*BR;
				mobod_W.setQToFitTranslation(someState, BR);
				system->realize(someState, SimTK::Stage::Position);
				std::cout << "Water (" << waterIx-1
				          <<") repositioned in " << mobod_W.expressVectorInGroundFrame(someState, BR) 
						  << std::endl << std::endl;

			}
		}
		system->realize(someState, SimTK::Stage::Dynamics);

		// Else, if not repositioned, integrate trajectory.
		//if(1) {
		else {
		initializeVelocities(someState);
		calcProposedKineticAndTotalEnergy(someState);

		integrateTrajectory(someState);
		} 

		system->realize(someState, SimTK::Stage::Dynamics);
		calcNewConfigurationAndEnergies(someState);
 */
	}else if(integratorName == IntegratorName::STATIONS_TASK){

		std::cout << "STATIONS_TASK\n";

		system->realize(someState, SimTK::Stage::Position);

		// Update target space
		world->updateTaskSpace(someState);
		SimTK::Array_<SimTK::Vec3> deltaStationP = world->getTaskSpaceDeltaStationP();
		std::cout << "deltaStationP " << deltaStationP << std::endl;

		SimTK::Matrix_<SimTK::Vec3> JS;
		world->calcStationJacobian(someState, JS);

		// Initialize velocities from Maxwell-Boltzmann distribution
		initializeVelocities(someState);
		//setVelocitiesToZero(someState);

		SimTK::Vector f;
		matter->multiplyByStationJacobianTranspose(someState, world->onBodyB[0], world->stationPInGuest[0], deltaStationP[0], f);

		SimTK::Vector givenU = SimTK::Vector(someState.getU().size());
		someState.setU(10*f);
		system->realize(someState, SimTK::Stage::Velocity);

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

	}else if (integratorName == IntegratorName::EMPTY){

/* 		// Alter Q in some way
		if(this->sampleGenerator == 1){
			
		} */

	}else{
		std::cout << "Warning UNKNOWN INTEGRATOR TREATED AS EMPTY\n";
	}

	return validateProposal();

}

bool HMCSampler::proposeNMA(SimTK::State& someState)
{

/* 	// Store old configuration
	storeOldConfigurationAndPotentialEnergies(someState); */

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

/**
 * Generate a trial move in the chain
*/
bool HMCSampler::generateProposal(SimTK::State& someState)
{
	// Return value
	bool validated = false;

	// TODO: delete
	std::cout << "DISTORT_OPTION " << this->DistortOpt << std::endl;

	// Propose
	for (int i = 0; i < 10; i++){

		// Normal Modes based trial
		if(DistortOpt > 0){
			validated = proposeNMA(someState);

		// Non-equilibrium trial - TFEP
		}else if(DistortOpt < 0){
			validated = proposeNEHMC(someState);

		// Equilibrium trial
		}else{
			validated = proposeEquilibrium(someState);
		}

		if (validated){
			break;
		}
	}

	// The sample was not validated
	if ( !validated ){ 
		// Warn user
		std::cout 
			<< "Warning: Proposal not validated after 10 tries.\n";

		// Reset
		this->acc = false;
		setSetConfigurationAndEnergiesToOld(someState);
	}

	return validated;

}

/** 
 * This is actually the Boltzmann factor not the probability 
 * because we don't have the partition function
*/
SimTK::Real HMCSampler::MHAcceptProbability(
	SimTK::Real argEtot_proposed, SimTK::Real argEtot_n) const 
{
	if(argEtot_n < argEtot_proposed) {
		return 1;
	} else {
		/* std::cout << "\tdiff=" << argEtot_n - argEtot_proposed 
			<< ", argEtot_n=" << argEtot_n 
			<< ", argEtot_proposed=" << argEtot_proposed
			<< ", beta=" << beta << std::endl; */
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

	CHECK_IF_NAN(ke_o);
	CHECK_IF_NAN(ke_n);

	// CHECK_IF_NAN(fix_o);
	CHECK_IF_NAN(fix_n);

	// CHECK_IF_NAN(logSineSqrGamma2_o);
	CHECK_IF_NAN(logSineSqrGamma2_n);

	CHECK_IF_NAN(timestep);
	CHECK_IF_NAN(exp(-(etot_n - etot_o) / RT));

	CHECK_IF_NAN(etot_n);
	CHECK_IF_NAN(etot_o);

	return true;
}

/** Chooses whether to accept a sample or not based on a probability **/
bool HMCSampler::acceptSample() {

	// The decision tree sets the value of internal variable acc
	if(this->alwaysAccept == true){ // Empty sampler
		this->acc = true;
	}else{ // Markov-Chain Monte Carlo
		const SimTK::Real rand_no = uniformRealDistribution(randomEngine);
		SimTK::Real prob = 0.0;
		if(DistortOpt == 0){
			prob = MHAcceptProbability(etot_o, etot_n);
		}else if(DistortOpt > 0){

			if(useFixman){
				prob = 
				MHAcceptProbability(pe_o + ke_prop_nma6 + fix_o,
									pe_n + ke_n_nma6 + fix_n);

			}else{
				prob =
				MHAcceptProbability(pe_o + ke_prop_nma6,
									pe_n + ke_n_nma6);

			}
		}else if(DistortOpt < 0){
			prob =
			MHAcceptProbability(pe_o + fix_o,
								pe_n + fix_n - this->bendStretchJacobianDetLog);
		}

		// std::cout << "\trand_no=" << rand_no << ", prob=" << prob 
		//<< ", beta=" << beta << std::endl;
		this->acc = (rand_no < prob);
	}

	return this->acc;
	
}

/**
 * Acception rejection step - not used
*/
bool HMCSampler::accRejStep(SimTK::State& someState) {

	// Empty sample generator
	if ( alwaysAccept == true ){
		// Simple molecular dynamics
		this->acc = true;
		std::cout << "\tsample accepted (simple molecular dynamics)\n";
		update(someState);
	} else {
		// we do not consider this sample accepted unless it passes all checks
		this->acc = false;

		// Apply Metropolis-Hastings criterion
		if(acceptSample()) {
			// sample is accepted
			this->acc = true;
			std::cout << "\tsample accepted (metropolis-hastings)\n";
			update(someState);
		} else {
			// sample is rejected
			this->acc = false;
			std::cout << "\tsample rejected (metropolis-hastings)\n";
			setSetConfigurationToOld(someState);
		}
	}

	return this->acc;
}

/**
 * The main function that generates a sample
*/
bool HMCSampler::sample_iteration(SimTK::State& someState)
{
/* std::cout << "Transforms before sample_iteration\n";
world->PrintAcosX_PFs();
world->PrintNormX_BMs();
world->PrintAcosX_PFMeans();
world->PrintNormX_BMMeans(); */

	// Set the number of decimals to be printed
	std::cout << std::setprecision(10) << std::fixed;

	// Store old configuration
	storeOldConfigurationAndPotentialEnergies(someState);

	// Generate a trial move in the stochastic chain
	//world->traceBendStretch(someState);

	// Set X_PF and X_BM means to whatever iti is now
	///////////////////////////////////////////////////
	if(nofSamples == 0){
		world->setTransformsMeansToCurrent(someState);
	}
	///////////////////////////////////////////////////

	// Calculate X_PFs and X_BMs
	world->getTransformsStatistics(someState);

	bool validated = true;
	validated = generateProposal(someState);
	
	// Apply the acceptance criterion
	if(validated){

		// Print all the energy terms first
		PrintDetailedEnergyInfo(someState);

		// Decide
		if(acceptSample()){
			// Print info
			std::cout << "\tsample accepted";
			if(this->alwaysAccept == true){
				std::cout << " (simple molecular dynamics)\n";
			}else{
				std::cout << " (metropolis-hastings)\n";
			}

			// UPDATE
			update(someState);

			// Deal with adaptive data
			storeAdaptiveData(someState);
			PrintAdaptiveData();

		}else{
			//Print info
			std::cout << "\tsample rejected (metropolis-hastings)\n";

			// RESET
			restore(someState);
		}

		// Calculate X_PFs and X_BMs
//world->traceBendStretch(someState);

		world->getTransformsStatistics(someState);

		// Update means of values before altering them
		if((this->nofSamples > 3000) && (this->nofSamples <= 6000)){
			//world->updateTransformsMeans(someState);
		}

	}

		// Increase the sample counter and return
		++nofSamples;
		return this->acc;

}

/**
*  Add generalized coordinates to a buffer
*/
void HMCSampler::updateQBuffer(const SimTK::State& someState)
{
	
	// g++17 complains if this is auto& or const auto&
	auto Q = someState.getQ(); 
	//std::cout << "Q = " << Q << std::endl;
	QsBuffer.insert(QsBuffer.end(), Q.begin(), Q.end());
	QsBuffer.erase(QsBuffer.begin(), QsBuffer.begin() + Q.size());
}

void HMCSampler::storeAdaptiveData(SimTK::State& someState)
{
	updateQBuffer(someState);
	pushCoordinatesInR(someState);
	pushVelocitiesInRdot(someState);
}

void HMCSampler::PrintAdaptiveData(void)
{
	// Calculate MSD and RRdot to adapt the integration length
	std::cout << std::setprecision(10) << std::fixed;
	std::cout << "\tMSD= " << calculateMSD() 
		<< ", RRdot= " << calculateRRdot() << std::endl;

}


/** Main function that contains all the 3 steps of HMC.
Implements the acception-rejection step and sets the state of the
compound to the appropriate conformation wether it accepted or not. **/
void HMCSampler::update(SimTK::State& someState)
{
	// Store final configuration and energy
	setSetConfigurationAndEnergiesToNew(someState);
	
	// Acceptance rate buffer
	++acceptedSteps;
	acceptedStepsBuffer.push_back(1);
	acceptedStepsBuffer.pop_front();
}

void HMCSampler::restore(SimTK::State& someState)
{
	// Restore configuration and energies to old
	setSetConfigurationAndEnergiesToOld(someState);

	// Update acceptance rate buffer
	acceptedStepsBuffer.push_back(0);
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
		<< "\tpe_o " << pe_o << ", pe_n " << pe_n << ", pe_nB "
		<< /*" turned off for time being"*/ getPEFromEvaluator(someState)
		<< "\n\tke_prop " << ke_o << ", ke_n " << ke_n
		<< "\n\tfix_o " << fix_o << ", fix_n " << fix_n
		<< "\n\tlogSineSqrGamma2_o " << logSineSqrGamma2_o
		<< ", logSineSqrGamma2_n " << logSineSqrGamma2_n
		//<< " detmbat_n " << detmbat_n //<< " detmbat_o " << detmbat_o << " "
		<< "\n\tts " << timestep  //<< ", exp(bdE) " << exp(-(etot_n - etot_proposed) / RT)
		<< "\n\t, etot_n " << etot_n  << ", etot_proposed " << etot_o
		<< ", JDetLog " << bendStretchJacobianDetLog
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
