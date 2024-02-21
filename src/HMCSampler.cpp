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

#define NAN_WARNING(n) \
	std::cout << "\t[WARNING]: " << #n << " is nan!\n"; \
	return false;

bool NAN_TO_INF(SimTK::Real& someNumber)
{
	if(std::isnan(someNumber)){
		//someNumber = SimTK::Infinity; // TODO: Victor check this please
		someNumber = std::numeric_limits<double>::max();
		return true;
	}else{
		return false;
	}
}


//** Constructor **/
//==============================================================================
//                           CONSTRUCTOR
//==============================================================================
// Description.
HMCSampler::HMCSampler(World &argWorld,
		SimTK::CompoundSystem &argCompoundSystem,
		SimTK::SimbodyMatterSubsystem &argMatter,
		std::vector<Topology> &argTopologies, 
		SimTK::DuMMForceFieldSubsystem &argDumm,
		SimTK::GeneralForceSubsystem &argForces,
		SimTK::TimeStepper &argTimeStepper) :
		Sampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper)
		//, MonteCarloSampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper)
{
	this->system = &argMatter.getSystem(); // TODO do we need this? was already defined in sampler.hpp

	//this->rootTopology = argResidue;

	// TODO do not throw in constructors
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

	// Acc
	this->acc = false;

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

/** Set simulation temperature,
velocities to desired temperature, variables that store the configuration
and variables that store the energies, both needed for the
acception-rejection step. Also realize velocities and initialize
the timestepper. **/
bool HMCSampler::initialize(SimTK::State& someState)
{
	//system->realize(someState, SimTK::Stage::Model);

	// After an event handler has made a discontinuous change to the
	// Integrator's "advanced state", 
	timeStepper->initialize(compoundSystem->getDefaultState());
	//timeStepper->initialize(someState);


	// Calculate Simbody configuration
	system->realize(someState, SimTK::Stage::Position);

	// Initialize QsBuffer with zeros
	if(this->nofSamples == 0){
		int totSize = QsBufferSize * matter->getNQ(someState);
		for(int i = 0; i < totSize; i++){
			QsBuffer.push_back(SimTK::Real(0));
		}
	}

	// Buffer to hold Q means
	if(this->nofSamples == 0){
		QsMeans.resize(matter->getNQ(someState), 0);
	}

	// Total mass of the system
	this->totalMass = matter->calcSystemMass(someState);

	return true;
}


/**
 * Set simulation temperature,
 * velocities to desired temperature, variables that store the configuration
 * and variables that store the energies, both needed for the
 * acception-rejection step. Also realize velocities and initialize
 * the timestepper.
 */
bool HMCSampler::reinitialize(SimTK::State& someState)
{
	// Set a validation flag
	bool validated = true;

	// After an event handler has made a discontinuous change to the
	// Integrator's "advanced state", this method must be called to
	// reinitialize the Integrator.
	if(this->nofSamples == 0){
		//timeStepper->initialize(compoundSystem->getDefaultState());
		timeStepper->initialize(someState);

		if(integratorName == IntegratorName::OMMVV){
		}
	}

	// Get no of degrees of freedom
	const int nu = someState.getNU();

	// Convenient variable to use for distortion detection
	pe_init = pe_set;


	// HMC: &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   X_O   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Calculate and set old configurationthis->natoms

	// Calculate Simbody configuration
	system->realize(someState, SimTK::Stage::Position);

	// Copy coordinates to OpenMM
	if(this->integratorName == IntegratorName::OMMVV){
		Simbody_To_OMM_setAtomsLocationsCartesian(someState);
	}

	// Initialize QsBuffer with zeros
	if(this->nofSamples == 0){
		int totSize = QsBufferSize * matter->getNQ(someState);
		for(int i = 0; i < totSize; i++){
			QsBuffer.push_back(SimTK::Real(0));
		}
	}

	// Store Simbody configuration
	storeSimbodyConfiguration_XFMs(someState);

	// Store OpenMM configuration
	if(this->integratorName == IntegratorName::OMMVV){
		omm_locations.resize(matter->getNumBodies());
		omm_locations_old.resize(matter->getNumBodies());

		OMM_storeOMMConfiguration_X(dumm->OMM_getPositions());
	}

	// Print Simbody
	//world->PrintFullTransformationGeometry(someState,
	//	true, true, true, true, true, true);


	// HMC: &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   PE_O   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Calculate and set old potential energy. This should call OpenMMPlugin to 
	// calculate energy and forces
	setOldPE(
		forces->getMultibodySystem().calcPotentialEnergy(someState)
		//dumm->CalcFullPotEnergyIncludingRigidBodies(someState) // NO OPENMM
	);


	// HMC: &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  FIX_O &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Calc old Fixman potential
	if(useFixman){
		
		// Internal Fixman potential
		setOldFixman(calcFixman(someState));

		// External Fixman potential
		setOldLogSineSqrGamma2( 
			((Topology *)rootTopology)->calcLogSineSqrGamma2(someState));

	}else{

		// Internal Fixman potential
		setOldFixman(0.0);

		// External Fixman potential
		setOldLogSineSqrGamma2(0.0);
	}

	// HMC: &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   PE_SET   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Store potential energies
	setSetPE(getOldPE());

	// HMC: &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  FIX_SET  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	// Set final Fixman potential to old
	if(useFixman){
		
		// Internal Fixman potential
		setSetFixman(getOldFixman());

		// External Fixman potential
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());

	}else{

		// Internal Fixman potential
		setSetFixman(getOldFixman());

		// External Fixman potential
		setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
	}

	// Initialize velocities to temperature
	// if(this->nofSamples == 0){
	// 	perturbVelocities(someState, VelocitiesPerturbMethod::TO_ZERO);
	// }

	// Reset ndofs
	ndofs = nu;

	// TODO replaced code beloe
	std::fill(UScaleFactors.begin(), UScaleFactors.end(), 1);
	std::fill(InvUScaleFactors.begin(), InvUScaleFactors.end(), 1);

	// // Initialize velocities to temperature
	// for (int j=0; j < nu; ++j){
	// 	UScaleFactors[j] = 1;
	// 	InvUScaleFactors[j] = 1;
	// }

	// Set the generalized velocities scale factors
	loadUScaleFactors(someState);

	// Buffer to hold Q means
	if(this->nofSamples == 0){
		QsMeans.resize(matter->getNQ(someState), 0);
	}

	// Total mass of the system
	this->totalMass = matter->calcSystemMass(someState);
	
	// Transformation Jacobian
	bendStretchJacobianDetLog = 0.0;

	// Set DuMM temperature : TODO: should propagate to OpenMM
	if(this->integratorName == IntegratorName::OMMVV){
		OMM_setDuMMTemperature(boostT);
	}
	

	return validated;

}

SimTK::Real HMCSampler::getMDStepsPerSampleStd() const
{
	return MDStepsPerSampleStd;
}

void HMCSampler::setMDStepsPerSampleStd(SimTK::Real mdstd){
	MDStepsPerSampleStd = mdstd;
}

// Set the method of integration
void HMCSampler::setSampleGenerator(SampleGenerator sampleGeneratorArg)
{
	if (SampleGenerator::EMPTY == sampleGeneratorArg) {
		setAlwaysAccept(true);
	}
	else if (SampleGenerator::MC == sampleGeneratorArg) {
		setAlwaysAccept(false);
	}
}

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

/** Store old and set kinetic and total energies */
void HMCSampler::storeOldAndSetKineticAndTotalEnergies(SimTK::State& someState)
{
	// Store kinetic energies
	if(this->integratorName == IntegratorName::OMMVV){
		this->ke_o = OMM_calcKineticEnergy();
	}else{
		this->ke_o = matter->calcKineticEnergy(someState);		
	}
	
	this->ke_o *= (this->unboostKEFactor);

	this->ke_set = this->ke_o;

	// Update total energies
	this->etot_o = getOldPE()
				 + getOldKE()
				 + getOldFixman()
				 + getOldLogSineSqrGamma2();

	this->etot_set = this->etot_o;

}

void HMCSampler::perturbPositions(SimTK::State& someState,
	PositionsPerturbMethod PPM)
{
	if(PPM == PositionsPerturbMethod::BENDSTRETCH){
	
		//std::cout << "Transforms before sample_iteration\n";
		//world->PrintAcosX_PFs();
		//world->PrintNormX_BMs();
		//world->PrintAcosX_PFMeans();
		//world->PrintNormX_BMMeans();
		//
		//world->traceBendStretch(someState);
		//
		// Set X_PF and X_BM means to whatever iti is now
		///////////////////////////////////////////////////
		// if(nofSamples == 0){
		//	world->setTransformsMeansToCurrent(someState);
		//}
		///////////////////////////////////////////////////

		// Resize the scale factors vector to be handed further
		//std::vector<SimTK::Real> scaleFactors;
		//scaleFactors.resize(world->acosX_PF00.size() + world->normX_BMp.size(),
		//	1.0);

		// Scale bonds and angles
		if(this->nofSamples >= 0){ // dont't take burn-in
			
			// Just for the Visualizer
			if(world->visual){
				//this->world->ts->stepTo(
				//	someState.getTime() + (0.00001));				
				std::cout << "Sleeping... before scaling " << std::flush;
				std::this_thread::sleep_for(std::chrono::milliseconds(6000));
				std::cout << "done.\n" << std::flush;
			}

			// Calculate X_PFs and X_BMs
			//world->getTransformsStatistics(someState);
			//world->updateTransformsMeans(someState);
			//
			//setQToScaleBendStretch(someState, scaleFactors);
			//setQToScaleBendStretchStdev(someState, scaleFactors);
			// DEBUG PRINT
			// World& myWorld = *world;
			// if(myWorld.ownWorldIndex == 1){
			// 	myWorld.PrintFullTransformationGeometry(someState,
			// 		false, false, false, true, false, false);
			// }

			calcSubZMatrixBATDeviations(someState);
			//PrintSubZMatrixBATDeviations(someState);

			SimTK::Real sJac =
				scaleSubZMatrixBATDeviations(someState, SimTK::Real(500.0/300.0));
			setDistortJacobianDetLog(sJac);
			//PrintSubZMatrixBATDeviations(someState);

			// Just for the Visualizer
			if(world->visual){
				//this->world->ts->stepTo(
				//	someState.getTime() + (0.00001));				
				std::cout << "Sleeping... after scaling " << std::flush;
				std::this_thread::sleep_for(std::chrono::milliseconds(6000));
				std::cout << "done.\n" << std::flush;
			}

		}

	}else{
		// Do nothing
	}
}


/**
 * Set velocities to 0
*/
void HMCSampler::setVelocitiesToZero(SimTK::State& someState)
{	
	// Set velocities to 0
	if(this->integratorName == IntegratorName::OMMVV){
		uint32_t seed = randomEngine() >> 32;
		dumm->setOpenMMvelocities(0, seed);
	}else{
		someState.updU() = 0;
	}

}

/** Set velocities according to the Maxwell-Boltzmann
distribution.  **/
void HMCSampler::setVelocitiesToGaussian(SimTK::State& someState)
{
	if (this->integratorName == IntegratorName::OMMVV){
		uint32_t seed = randomEngine() >> 32;
		dumm->setOpenMMvelocities(this->boostT, seed);
	} else {
		// Check if we can use our cache
		const int nu = someState.getNU();
		if (nu != RandomCache.getNU()) {
			// Rebuild the cache
			// We also get here if the cache is not initialized
			RandomCache.initialize(nu);
			sqrtMInvV.resize(nu);
		} else {
			// Wait for random number generation to finish (should be done by this stage)
			RandomCache.wait();
		}

		// Scale by square root of the inverse mass matrix
		matter->multiplyBySqrtMInv(someState, RandomCache.getV(), sqrtMInvV);

		// Set stddev according to temperature
		sqrtMInvV *= sqrtBoostRT;

		// Raise the temperature
		someState.updU() = sqrtMInvV;

		// Ask for a number of random numbers and check if we are done the next
		// time we hit this function
		RandomCache.generateGaussianVelocities();
	}
}


/** Initialize velocities according to the Maxwell-Boltzmann
 * distribution.  Coresponds to R operator in LAHMC.
 * It takes the boost factor into account 
 * */
void HMCSampler::perturbVelocities(SimTK::State& someState,
	VelocitiesPerturbMethod VPM){



	if((VPM == VelocitiesPerturbMethod::TO_ZERO) || (this->boostRT == 0)){

		// Set velocities to 0
		setVelocitiesToZero(someState);
	
	}else if(this->DistortOpt > 0){
		setVelocitiesToNMA(someState);
	}else{
		// Set velocities to Gaussian
		setVelocitiesToGaussian(someState);
	}

	// Realize velocity
	system->realize(someState, SimTK::Stage::Velocity);

	// Store kinetic and total energies
	storeOldAndSetKineticAndTotalEnergies(someState);

}

/** Initialize velocities scaled by NMA factors.
 Coresponds to R operator in LAHMC **/
void HMCSampler::setVelocitiesToNMA(SimTK::State& someState)
{
	// Get the number of dofs
	const int nu = someState.getNU();

	// Get sqrt of kT
	double sqrtRT = std::sqrt(RT);

	// U vector
	SimTK::Vector Us(nu, 1);

	// Vector multiplied by the sqrt of the inverse mass matrix
	// TODO move to class member to store memory
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
		if (nu != RandomCache.getNU()) {
			// Rebuild the cache
			// We also get here if the cache is not initialized
			RandomCache.initialize(nu);
		} else {
			// Wait for random number generation to finish (should be done by this stage)
			RandomCache.wait();
		}

		Us = RandomCache.getV();

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
		if (nu != RandomCache.getNU()) {
			// Rebuild the cache
			// We also get here if the cache is not initialized
			RandomCache.initialize(nu);
		} else {
			// Wait for random number generation to finish (should be done by this stage)
			RandomCache.wait();
		}

		Us = RandomCache.getV();

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
		double lengthNoise = lengthNoiser(randomEngine);
		std::cout << "lengthNoise " << lengthNoise << "\n";
		for(int j = 0; j < ndofs; j++){Us[j] *= lengthNoise;}
		std::cout << "Us 4 "; for(int i = 0; i < nu; i++){ std::cout << Us[i] << " " ;} std::cout << "\n";

		// Multiply by the square root of the mass matrix
		matter->multiplyBySqrtMInv(someState, Us, sqrtMInvUs);
		std::cout << "sqrtMInvUs 5 "; for(int i = 0; i < nu; i++){ std::cout << sqrtMInvUs[i] << " " ;} std::cout << "\n";

	}else if(DistortOpt == 6){

		//std::cout << "DOFUs 0 "; for(int i = 0; i < nu; i++){ std::cout << DOFUScaleFactors[i] << " " ;} std::cout << "\n";

		// Sample from N(0, 1)
		if (nu != RandomCache.getNU()) {
			// Rebuild the cache
			// We also get here if the cache is not initialized
			RandomCache.initialize(nu);
		} else {
			// Wait for random number generation to finish (should be done by this stage)
			RandomCache.wait();
		}

		Us = RandomCache.getV();
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
		// bool b = false;
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

		/* std::cout << "Keval terms bXMX muMmu bmuMX LSE(bmuMX) = "
			<< 0.5*bXMX << " " // << 0.5*bmuMmu << " "
			<< bmuMX << " " << bcorr << std::endl; */

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

		//std::cout << "ke_prop_nma6 " << ke_prop_nma6 << "\n";



	}

	// Set state velocities
	someState.updU() = sqrtMInvUs;

	// Realize velocity
	system->realize(someState, SimTK::Stage::Velocity);

	RandomCache.generateGaussianVelocities();
}

/*
* Integrate trajectory
*/
void HMCSampler::integrateTrajectory(SimTK::State& someState){

	// Adapt timestep
	if(shouldAdaptTimestep){
		adaptTimestep(someState);
	}

	if(this->integratorName == IntegratorName::VERLET){
		try {

			// Call Simbody TimeStepper to advance time
			this->world->ts->stepTo(
				someState.getTime() + (timestep*MDStepsPerSample));

			system->realize(someState, SimTK::Stage::Position);

		}catch(const std::exception&){

			std::cout << "\t[WARNING] propose exception caught!\n";
			proposeExceptionCaught = true;

			assignConfFromSetTVector(someState);

		}

	}else if(this->integratorName == IntegratorName::BOUND_WALK){
		try {

			// Call Simbody TimeStepper to advance time
			integrateTrajectory_Bounded(someState);

		}catch(const std::exception&){

			std::cout << "\t[WARNING] propose exception caught!\n";
			proposeExceptionCaught = true;

			assignConfFromSetTVector(someState);
			
		}

	}else if(this->integratorName == IntegratorName::BOUND_HMC){
		try {

			// Call Simbody TimeStepper to advance time
			integrateTrajectory_BoundHMC(someState);

		}catch(const std::exception&){

			std::cout << "\t[WARNING] propose exception caught!\n";
			proposeExceptionCaught = true;

			assignConfFromSetTVector(someState);
			
		}

	}else if(this->integratorName == IntegratorName::STATIONS_TASK){
		try {

			// Call Simbody TimeStepper to advance time
			integrateTrajectory_TaskSpace(someState);

		}catch(const std::exception&){

			std::cout << "\t[WARNING] propose exception caught!\n";
			proposeExceptionCaught = true;

			assignConfFromSetTVector(someState);
			
		}

	}else if(this->integratorName == IntegratorName::OMMVV){
		try {
			
			// This code works for updating simbody bodies
			// each body should be an atom
			assert(matter->getNumBodies() == dumm->getNumAtoms() + 1);

			// Actual openmm integration
			dumm->OMM_integrateTrajectory(this->MDStepsPerSample);

			// Somewhere, the topology gets ruined
			system->realizeTopology();

			// Transfer back to Simbody (TODO: might be redundant)
			OMM_To_Simbody_setAtomsLocations(someState); // COMPLETE

		}catch(const std::exception&){

			// Send general message
			std::cout << "\t[WARNING] propose exception caught!\n";
			proposeExceptionCaught = true;

			OMM_restoreConfiguration(someState);

			// Transfer back to Simbody (TODO: might be redundant)
			OMM_To_Simbody_setAtomsLocations(someState); // COMPLETE

		}

	}else if(this->integratorName == IntegratorName::EMPTY){
		try {

			// Advance to Position Stage
			system->realize(someState, SimTK::Stage::Dynamics);

		}catch(const std::exception&){

			std::cout << "\t[WARNING] propose exception caught!\n";
			proposeExceptionCaught = true;

			assignConfFromSetTVector(someState);
			
		}

	}else{ // Anything else is consodered Empty
		try {

			// Advance to Position Stage
			system->realize(someState, SimTK::Stage::Dynamics);

		}catch(const std::exception&){

			std::cout << "\t[WARNING] propose exception caught!\n";
			proposeExceptionCaught = true;
			
			assignConfFromSetTVector(someState);

		}

	}



}


void HMCSampler::OMM_integrateTrajectory(SimTK::State& someState)
{
	assert(!"Not implemented");
}

// Trajectory length has an average of MDStepsPerSample and a given std
void HMCSampler::integrateVariableTrajectory(SimTK::State& someState){
	try {

		float timeJump = 0;
		if(MDStepsPerSampleStd != 0){
			std::cout << "Got MDSTEPS_STD " << MDStepsPerSampleStd << '\n';

			// TODO moe to internal variables set
				std::normal_distribution<float> trajLenNoiser(1.0, MDStepsPerSampleStd);
				float trajLenNoise = trajLenNoiser(randomEngine);
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
			//SimTK::Real PE = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // DOESN'T WORK WITH OPENMM
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

/*
* Integrate trajectory: with boundary conditions using random 
* spherical distribution
	Hard-coded, need to be able to propagate from input
but for the moment this will have to do.

These are integers, since you'd never *really* need a long list
of bounders
const vector<vector<int>> bounders{{0},{2}};
These are strings, since you could  want all solvent
molecules to be bound, so the "s" option is implemented
const vector<vector<string>> boundeds{{"2"},{"s"}};
These are strings, since you could want the 
const vector<vector<string>> bounder_subset{{"74","88"},{"a"}};
const vector<float> bound_radii; */
void HMCSampler::integrateTrajectory_Bounded(SimTK::State& someState){

	assert(!"Not implemented");

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
		
		SimTK::Quaternion randQuat = generateRandomQuaternion();

		SimTK::Quaternion currQuat(someState.updQ()[0], someState.updQ()[1],
			someState.updQ()[2], someState.updQ()[3]);

		SimTK::Quaternion resultQuat = multiplyQuaternions(randQuat, currQuat);

		mobod_L.setQToFitRotation(someState, SimTK::Rotation(resultQuat));
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

		system->realize(someState, SimTK::Stage::Dynamics);
}

/*
* Integrate trajectory: with boundary conditions using random 
* spherical distribution and HMC
*/
void HMCSampler::integrateTrajectory_BoundHMC(SimTK::State& someState){
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
			perturbVelocities(someState);
			calcProposedKineticAndTotalEnergyOld(someState);

			integrateTrajectory(someState);
			system->realize(someState, SimTK::Stage::Dynamics);
		}
}


/** Integrate trajectory using task space forces */
void HMCSampler::integrateTrajectory_TaskSpace(SimTK::State& someState){

		/* std::cout << "STATIONS_TASK\n";
		system->realize(someState, SimTK::Stage::Position);

		// Update target space
		world->updateTaskSpace(someState);
		SimTK::Array_<SimTK::Vec3> deltaStationP = world->getTaskSpaceDeltaStationP();
		std::cout << "deltaStationP " << deltaStationP << std::endl;

		SimTK::Matrix_<SimTK::Vec3> JS;
		world->calcStationJacobian(someState, JS);

		SimTK::Vector taskForce;
		
		matter->multiplyByStationJacobianTranspose(someState,
			world->onBodyB[0], world->taskStationPInGuest[0],
			deltaStationP[0], taskForce);
		
		std::cout << "state numU " << someState.getNU()
			<< " taskForce " << -1.0 * taskForce
			<< std::endl;

		someState.setU(-1 * taskForce); */

		/* // Topologies and target atoms
		int hostTopology = 0;
		int guestTopology = 1;
		std::vector<int> bAtomIxs_host = {4}; // atoms on host topology
		std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology

		// Get stations in guest
		int topi = -1;
		for(auto& topology : (topologies)){ // topologies
			topi++;
			if(topi == guestTopology){
				int tz = -1;
				for (int bAtomIx : bAtomIxs_guest) { // atoms
					tz++;

					// Get body
					SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).compoundAtomIndex;
					SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(aIx);
					SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

					// Get mobod X_GB, X_PF, X_FM and X_BM
					const Transform& X_GB = mobod.getBodyTransform(someState);
					const Transform& X_PF = mobod.getInboardFrame(someState);
					const Transform& X_FM = mobod.getMobilizerTransform(someState);
					const Transform& X_BM = mobod.getOutboardFrame(someState);
					
					const Transform& X_GF = X_GB * X_BM * (~X_FM);
					SimTK::Vec3 G_S = deltaStationP[0];
					SimTK::Vec3 F_S = (~X_GF) * G_S;

					mobod.setOneU(someState, 0, F_S[0]);
					mobod.setOneU(someState, 1, F_S[1]);
					mobod.setOneU(someState, 2, F_S[2]);

				}
			}
		} */

		system->realize(someState, SimTK::Stage::Velocity);

}


// ELIZA OPENMM FULLY FLEXIBLE INTEGRATION CODE

// ELIZA: Insert code here
void HMCSampler::OMM_setDuMMTemperature(double HMCBoostTemperature){
	dumm->setDuMMTemperature(HMCBoostTemperature);
}

// ELIZA: Insert code here
double HMCSampler::OMM_calcKineticEnergy(void){
	return dumm->OMM_calcKineticEnergy();
}

// ELIZA: Insert code here
double HMCSampler::OMM_calcPotentialEnergy(void){
	return dumm->OMM_calcPotentialEnergy();
}

void HMCSampler::OMM_storeOMMConfiguration_X(const std::vector<OpenMM::Vec3>& positions)
{
		omm_locations_old[0] = SimTK::Vec3(0, 0, 0);

		for (int i = 0; i < positions.size(); i++) {
			omm_locations_old[i + 1] = SimTK::Vec3(
				positions[i][0],
				positions[i][1],
				positions[i][2]);
		}

}


void HMCSampler::OMM_restoreConfiguration(SimTK::State& someState)
{

	// New code
	// Restore OpenMM positions
	std::vector<SimTK::Vec3>::const_iterator it_begin = 
		omm_locations_old.begin() + size_t(1);
	std::vector<SimTK::Vec3>::const_iterator it_end = 
		omm_locations_old.end();

	const std::vector<SimTK::Vec3> omm_locations_old_1(it_begin, it_end);

	dumm->OMM_setOpenMMPositions(omm_locations_old_1);

	// Reset Simbody (may not be necessary)
	OMM_To_Simbody_setAtomsLocations(someState);

	// Old code
	/* // Restore configuration
	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		mobod.setQToFitTransform(someState, SetTVector[mbx - 1]);
	}
	system->realize(someState, SimTK::Stage::Position);

	omm_locations[0] = SimTK::Vec3(0, 0, 0);

	std::vector<OpenMM::Vec3>& positions = dumm->OMM_updatePositions();

	for (int i = 0; i < positions.size(); i++) {
		positions[i + 1] = OpenMM::Vec3(
			omm_locations[i][0],
			omm_locations[i][1],
			omm_locations[i][2]);
	} */


}


/**
 * Order should be DuMM::NonbondedAtomIx 
 * Update OpenMM position from a Simbody Cartesian world
 * */
void HMCSampler::Simbody_To_OMM_setAtomsLocationsCartesian(SimTK::State& someState,
	bool throughDumm)
{

	/* std::cout << "HMCSampler::Simbody_To_OMM "
		<< "someState.getTime() " << someState.getTime()
		<< std::endl << std::flush; */

	std::vector<SimTK::Vec3> startingPos;
	startingPos.resize(this->natoms);

	system->realize(someState, SimTK::Stage::Position);

	if(throughDumm){

		const Vector_<Vec3>& DuMMIncludedAtomStationsInG = 
			dumm->getIncludedAtomPositionsInG(someState);

		for(int atomCnt = 0; atomCnt < this->natoms; atomCnt++){
			startingPos[atomCnt] = DuMMIncludedAtomStationsInG[atomCnt];
		}

		/* //for(int atomCnt = 0; atomCnt < this->natoms; atomCnt++){
		for(int atomCnt = 0; atomCnt < std::min(10, static_cast<int>(this->natoms)); atomCnt++){
			std::cout << "atom " << atomCnt << " "
				<< startingPos[atomCnt][0] << " "
				<< startingPos[atomCnt][1] << " "
				<< startingPos[atomCnt][2] << " "
				<< std::endl;
		}

		// world->PrintFullTransformationGeometry(someState,
		//	false, false, false, true, true, true);
		std::cout << "numBodies " <<  matter->getNumBodies() << std::endl << std::flush;

		//for(SimTK::MobilizedBodyIndex mbx(1); mbx < (matter->getNumBodies()); mbx++){
		for(SimTK::MobilizedBodyIndex mbx(1);
		mbx < std::min(10, matter->getNumBodies());
		mbx++){

			SimTK::MobilizedBody mobod = matter->getMobilizedBody(mbx);

			//SimTK::Vec3 mobodOrig = mobod.expressVectorInGroundFrame(someState, SimTK::Vec3(0, 0, 0));
			SimTK::Vec3 mobodOrig = mobod.getBodyOriginLocation(someState);
			std::cout << "mobodOrig " << int(mbx) << " "
				<< mobodOrig[0] << " "
				<< mobodOrig[1] << " "
				<< mobodOrig[2] << " "
				<< std::endl;		
		} */

	}else{

		for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
			SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
			startingPos.push_back(mobod.getBodyOriginLocation(someState));
		}
	}

	// Apply
	dumm->OMM_setOpenMMPositions(startingPos);

}

// Transfer coordinates from openmm to simbody
void HMCSampler::OMM_To_Simbody_setAtomsLocations(SimTK::State& someState)
{

		/* std::cout << "HMCSampler::OMM_To_Simbody_setAtomsLocations BEFORE " << std::endl;
		world->PrintFullTransformationGeometry(someState); */

		//omm_locations.resize(matter->getNumBodies());

		omm_locations[0] = SimTK::Vec3(0, 0, 0);

		const std::vector<OpenMM::Vec3>& positions = dumm->OMM_getPositions();

		for (int i = 0; i < positions.size(); i++) {
			omm_locations[i + 1] = SimTK::Vec3(
				positions[i][0],
				positions[i][1],
				positions[i][2]);
		}

		// Invalidate all statges
		matter->invalidateSubsystemTopologyCache();

		for (int i = 0; i < dumm->getNumAtoms(); i++) {

			SimTK:DuMM::AtomIndex aix(i);
			SimTK::MobilizedBodyIndex mbx = dumm->getAtomBody(aix);
			SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
			
			const auto location = omm_locations[i + 1];
			const auto parent = omm_locations[mobod.getParentMobilizedBody().getMobilizedBodyIndex()];

			mobod.updateDefaultFrames(
				Transform(Rotation(), location),
				Transform(Rotation(), parent));

			mobod.setQToFitTransform(someState, Transform(Rotation()));

		}

		system->realizeTopology();

		compoundSystem->realize(someState, SimTK::Stage::Position);

		/* std::cout << "HMCSampler::OMM_To_Simbody_setAtomsLocations AFTER " << std::endl;
		world->PrintFullTransformationGeometry(someState); */


}


void HMCSampler::OMM_PrintLocations(void)
{
	const auto positions = dumm->OMM_getPositions();

	std::cout << "OMM locations" << std::endl;

	for (int i = 0; i < positions.size(); i++) {
		std::cout
			<< positions[i][0] << " "
			<< positions[i][1] << " "
			<< positions[i][2] << " "
			<< std::endl;
	}	
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

// Set old Fixman potential
void HMCSampler::setOldFixman(SimTK::Real argFixman)
{
	this->fix_o = argFixman;
}

// Get Fixman potential
SimTK::Real HMCSampler::getNewFixman(void) const
{
	return this->fix_n;
}


// Set old Fixman potential
void HMCSampler::setNewFixman(SimTK::Real argFixman)
{
	this->fix_n = argFixman;
}

// Get set Fixman potential
SimTK::Real HMCSampler::getSetFixman(void) const
{
	return this->fix_set;
}

// Set set Fixman potential
void HMCSampler::setSetFixman(SimTK::Real argFixman)
{
	this->fix_set = argFixman;
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

void HMCSampler::setDistortJacobianDetLog(SimTK::Real argJ)
{
	this->bendStretchJacobianDetLog = argJ;
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
	//return dumm->CalcFullPotEnergyIncludingRigidBodies(someState);// DOESN'T WORK WITH OPENMM

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
void HMCSampler::storeSimbodyConfiguration_XFMs(const SimTK::State& someState)
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
	const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
	mobod.setQToFitTransform(someState, SetTVector[i]);
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

	// TODO
    assert(!"What should we return here?");
    return -1;
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

	if(this->integratorName == IntegratorName::OMMVV){
		OMM_setDuMMTemperature(this->boostT);
	}

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
HMCSampler::storeOldPotentialEnergies(
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
void HMCSampler::calcProposedKineticAndTotalEnergyOld(SimTK::State& someState){

	// Get proposed kinetic energy
	if(integratorName == IntegratorName::OMMVV){
		this->ke_o = OMM_calcKineticEnergy();

	}else{
		this->ke_o = matter->calcKineticEnergy(someState);
	}

	// Unboost (from guidance to evaluation Hamiltonian)
	this->ke_o *= (this->unboostKEFactor);

	// Store proposed total energy
	this->etot_o = getOldPE()
				 + getOldKE() 
				 + getOldFixman()
				 + getOldLogSineSqrGamma2();
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

/** Cartesian Fixman potential */
SimTK::Real HMCSampler::CartesianFixmanPotential(void){
	std::cerr 
		<< "Attempting Fixman potential calculation with OpenMM integrators.";
	throw std::exception();
	std::exit(1);
}

/** Store new configuration and energy terms **/
void HMCSampler::calcNewEnergies(SimTK::State& someState)
{

	// Get new potential energy
	if(this->integratorName == IntegratorName::OMMVV){
		pe_n = OMM_calcPotentialEnergy();
	}else{
		pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
		//pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState);// DOESN'T WORK WITH OPENMM
		// TODO: replace with the following after checking is the same thing
		//pe_n = compoundSystem->calcPotentialEnergy(someState);
	}

	// Get new Fixman potential
	if(useFixman){
		if(this->integratorName == IntegratorName::OMMVV){
			fix_n = CartesianFixmanPotential();
		}else{
			fix_n = calcFixman(someState);
			logSineSqrGamma2_n = ((Topology *)rootTopology)->calcLogSineSqrGamma2(someState);
		}
	}else{
		fix_n = 0.0;
		logSineSqrGamma2_n = 0.0;
	}

	// Get new kinetic energy
	if(this->integratorName == IntegratorName::OMMVV){
		ke_n = OMM_calcKineticEnergy();
	}else{
		ke_n = matter->calcKineticEnergy(someState);
		system->realize(someState, SimTK::Stage::Velocity);
	}

	ke_n *= (this->unboostKEFactor); // TODO: Check

	// Get new total energy
	if(useFixman){
		etot_n = pe_n + ke_n + fix_n - (0.5 * RT * logSineSqrGamma2_n);
		etot_o = pe_o + ke_o + fix_o - (0.5 * RT * logSineSqrGamma2_o);
	}else{
		etot_n = pe_n + ke_n;
		etot_o = pe_o + ke_o;
	}

}


/** Set energies and configuration to new state **/
void HMCSampler::setSetConfigurationAndEnergiesToNew(
	SimTK::State& someState)
{		

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
SimTK::Real HMCSampler::setQToScaleBendStretch(SimTK::State& someState,
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


	// Set the Jacobian and return
	unsigned int startFromBody = 2;
	bendStretchJacobianDetLog =
		calcBendStretchJacobianDetLog(someState, scaleFactors, startFromBody);
	return bendStretchJacobianDetLog;

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
SimTK::Real HMCSampler::setQToShiftBendStretchStdev(SimTK::State& someState,
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


	// Set the Jacobian and return
	unsigned int startFromBody = 2;
	bendStretchJacobianDetLog =
		calcBendStretchJacobianDetLog(someState, scaleFactors, startFromBody);
	return bendStretchJacobianDetLog;

}

/**
 * Shift all the generalized coordinates to scale bonds and angles
 * standard deviations through BendStretch joint
*/
SimTK::Real
HMCSampler::setQToScaleBendStretchStdev(
	SimTK::State& someState,
	std::vector<SimTK::Real>& scaleFactors)
{

	// Print the scale factor
	std::cout << "shiftQ Got " << this->QScaleFactor << " scale factor "
		<< std::endl;

	if (this->QScaleFactor == 1){
		return 1;
	}
	
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

	// Iterate mobods
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
	} // every mobod

	//std::cout << "scaleFactors: "; PrintCppVector(scaleFactors);

	// Save changes by advancing to Position Stage
	system->realize(someState, SimTK::Stage::Topology);
	system->realize(someState, SimTK::Stage::Position);

	// Test
	if(0){
		matter->realizeArticulatedBodyInertias(someState);
		SimTK::Vector v(someState.getNU());
		SimTK::Vector MInvV(someState.getNU());
		SimTK::Real detM = 0.0;
		matter->calcDetM(someState, v, MInvV, &detM);
		std::cout << "logDetM " << detM << std::endl;
		/* //std::cout << "shifted Q = " << someState.getQ() << std::endl;
		// Get bonds and angles values
		for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
			const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
			const Transform& X_PF = mobod.getInboardFrame(someState);
			const Transform& X_FM = mobod.getMobilizerTransform(someState);
			const Transform& X_BM = mobod.getOutboardFrame(someState);

			std::cout << "HMCSampler check world " << world->ownWorldIndex << " " 
				<< "bond " << int(mbx) - 1 << " " << X_BM.p().norm() << " "
				<< std::endl;
		} */
	}

	// Set the Jacobian and return
	unsigned int startFromBody = 2;
	bendStretchJacobianDetLog =
		calcBendStretchJacobianDetLog(someState, scaleFactors, startFromBody);
	return bendStretchJacobianDetLog;

}

/**
 * Get the log of the Jacobian of a bond-angle strtch
*/
SimTK::Real 
HMCSampler::calcBendStretchJacobianDetLog(SimTK::State& someState,
	std::vector<SimTK::Real> scaleFactors,
	unsigned int startFromBody)
{

	int x_pf_k = 0;

	// Get log of the Cartesian->BAT Jacobian
	SimTK::Real logJacBAT_0 = 0.0; // log space
	startFromBody -= 1; // align with k (exclude Ground)
	for(unsigned int k = startFromBody; k < world->normX_BMp.size(); k++){

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
			}
			//std::cout << "k d^2 " << k << " " << d2 << std::endl;

			// Get the angle term
			SimTK::Real sineTheta = std::abs(world->acosX_PF00[k]);

			if(sineTheta != 0){
					sineTheta = std::log(std::sin(sineTheta)); // log space
			}
			//std::cout << "k sinTh " << k << " " << sineTheta << std::endl;


			// Accumulate result
			logJacBAT_0 += (d2 + sineTheta); // log space
	}

	// Get log of the scaling Jacobian
	// Although the method seems stupid for now, it has generality
	// and allows insertion of new code
	SimTK::Real logJacScale = 0.0; // log space
	/* for(unsigned int k = startFromBody; k < scaleFactors.size(); k++){

			// Accumulate result
			if(scaleFactors[k] != 0){
					logJacScale -= std::log(std::abs(scaleFactors[k])); // log space
			}

			//std::cout << "sf " << k << " " << scaleFactors[k] << std::endl; 
	} */

	int internNdofs = this->ndofs - 
		matter->getMobilizedBody(SimTK::MobilizedBodyIndex(1)).getNumU(someState);

	logJacScale =  internNdofs * std::log( ( this->QScaleFactor ) );

	// Get log of the BAT->Cartesian Jacobian after scaling
	SimTK::Real logJacBAT_tau = 0.0; // log space
	for(unsigned int k = startFromBody; k < world->normX_BMp.size(); k++){

			// scaleFactors index is shifted
			x_pf_k = k + world->normX_BMp.size();

			// Get bond term
			SimTK::Real d2 = (world->normX_BMp[k] * scaleFactors[k]);

			if(d2 != 0){
					d2 = 2.0 * std::log(d2); // log space
			}
			//std::cout << "k d^2 " << k << " " << d2 << std::endl;

			// Get the angle term
			SimTK::Real sineTheta = std::abs(world->acosX_PF00[k] * scaleFactors[x_pf_k]);

			if(std::abs(sineTheta) > 0.0000001){
					sineTheta = std::log(std::sin(sineTheta)); // log space
  			}
			//std::cout << "k sinTh " << k << " " << sineTheta << std::endl;

			// Accumulate result
			logJacBAT_tau += (d2 + sineTheta); // log space
	}

	// Final result
	//SimTK::Real logBendStretchJac = (-1.0 * logJacBAT_0) + logJacScale + logJacBAT_tau;
	//SimTK::Real logBendStretchJac = logJacScale;
	//SimTK::Real logBendStretchJac = -1.0 * logJacScale;
	//SimTK::Real logBendStretchJac = (-1.0 * logJacBAT_0) + (1.0 * logJacBAT_tau);
	//SimTK::Real logBendStretchJac = logJacBAT_0 - logJacScale + (-1.0 * logJacBAT_tau);
	//SimTK::Real logBendStretchJac = (-1.0 * logJacBAT_0) + logJacScale + logJacBAT_tau;
	SimTK::Real logBendStretchJac = logJacBAT_0 + logJacScale + (-1.0 * logJacBAT_tau);

	std::cout << "logJacBAT " << logJacBAT_0
			<< " logJacScale " << logJacScale
			<< " logJacBATInv " << logJacBAT_tau
			<< " logBendStretchJac " << logBendStretchJac
            << std::endl;


	//logBendStretchJac = std::log(this->QScaleFactor); 
	//std::cout << "LNJ_HARDCODED_s_lnJ " << this->QScaleFactor << " " <<  logBendStretchJac << std::endl;

	return logBendStretchJac;

}

// Set Sphere Radius when doing RANDOM_WALK
void HMCSampler::setSphereRadius(float argSphereRadius)
{
	sphereRadius = argSphereRadius;
	std::cout 
		<< "HMCSampler::setSphereRadius:  Radius Set: " 
		<< sphereRadius << std::endl;
}

/** Returns the 'how' argument of initializeVelocities */
VelocitiesPerturbMethod HMCSampler::velocitiesPerturbMethod(void)
{
	// How do we initialize velocities

	if(integratorName == IntegratorName::OMMVV){
		return VelocitiesPerturbMethod::TO_T;

	}else if (integratorName == IntegratorName::VERLET){
		return VelocitiesPerturbMethod::TO_T;

	}else if (integratorName == IntegratorName::BOUND_WALK){
		return VelocitiesPerturbMethod::TO_ZERO;

	}else if (integratorName == IntegratorName::BOUND_HMC){
		return VelocitiesPerturbMethod::TO_ZERO;

	}else if(integratorName == IntegratorName::STATIONS_TASK){
		return VelocitiesPerturbMethod::TO_ZERO;

	}else if (integratorName == IntegratorName::EMPTY){
		return VelocitiesPerturbMethod::TO_ZERO;

	}else{
		return VelocitiesPerturbMethod::TO_ZERO;
	}

	return VelocitiesPerturbMethod::TO_ZERO;

}

/** Returns the 'how' argument of initializeVelocities */
PositionsPerturbMethod HMCSampler::positionsPerturbMethod(void)
{
	
	PositionsPerturbMethod how = PositionsPerturbMethod::EMPTY;

	if(this->DistortOpt < 0){
		how = PositionsPerturbMethod::BENDSTRETCH;

	}else{
		how = PositionsPerturbMethod::EMPTY;

	}

	return how;
}

/**
 * It implements the proposal move in the Hamiltonian Monte Carlo
 * algorithm. It essentially propagates the trajectory.
 **/
bool HMCSampler::propose(SimTK::State& someState)
{

	bool proposalValidation = true;

	for (int i = 0; i < 10; i++){

		// Adapt Gibbs blocks (Transformer)
		bool shouldAdaptWorldBlocks = false;
		if(shouldAdaptWorldBlocks){
			adaptWorldBlocks(someState);
		}

		// Perturb positions 
		perturbPositions(someState, positionsPerturbMethod());

		// Perturb velocities
		perturbVelocities(someState, velocitiesPerturbMethod());

		// Perturb forces
		/* perturbForces(someState,
			forcePerturbationMethod(this->integratorName)); */

		// Store the proposed energies
		calcProposedKineticAndTotalEnergyOld(someState);

		// Integrate trajectory
		integrateTrajectory(someState);

		// Get all new energies after integration
		calcNewEnergies(someState);

		// TODO: Any validation should be inserted here

		// Check if any energies are nan
		// or any exception during perturbations or integrations
		proposalValidation = checkExceptionsAndEnergiesForNAN()
			&& proposalValidation;

		// Check if molecule was distorted during
		// perturbations or integrations
		proposalValidation = checkDistortionBasedOnE(pe_n - pe_o)
			&& proposalValidation;

		if(proposalValidation){
			return true;
		}

	}

	return false;

}

/** Checks if the proposal is valid **/
bool HMCSampler::checkExceptionsAndEnergiesForNAN() {
	// TODO should we check anything else?
	// TODO do we need to check the old values?

	if(proposeExceptionCaught) {
		std::cout << "\t[WARNING] invalid sample: propose exception caught!\n";
		return false;
	}

	if (NAN_TO_INF(pe_n)){
		NAN_WARNING(pe_n);
		return false;
	}

	if (NAN_TO_INF(ke_o)){
		NAN_WARNING(ke_o);
		return false;
	}
	if (NAN_TO_INF(ke_n)){
		NAN_WARNING(ke_n);
		return false;
	}

	if (NAN_TO_INF(fix_n)){
		NAN_WARNING(fix_n);
		return false;
	}

	if (NAN_TO_INF(logSineSqrGamma2_n)){
		NAN_WARNING(logSineSqrGamma2_n);
		return false;
	}

	if (NAN_TO_INF(timestep)){
		NAN_WARNING(timestep);
		return false;
	}

	CHECK_IF_NAN(exp(-(etot_n - etot_o) / RT)); // TODO:

	if (NAN_TO_INF(etot_o)){
		NAN_WARNING(etot_o);
		return false;
	}

	if (NAN_TO_INF(etot_n)){
		NAN_WARNING(etot_n);
		return false;
	}

	return true;
}

/** 
 * This is actually the Boltzmann factor not the probability 
 * because we don't have the partition function
*/
SimTK::Real HMCSampler::MetropolisHastings(
	SimTK::Real argEtot_proposed,
	SimTK::Real argEtot_n,
	SimTK::Real lnJ) const 
{

	SimTK::Real dE = argEtot_n - argEtot_proposed;
	SimTK::Real beta_dE = this->beta * dE;
	SimTK::Real beta_dE_mJ = beta_dE - lnJ;

	/* std::cout << "\tdiff=" << argEtot_n - argEtot_proposed 
		<< ", argEtot_n=" << argEtot_n 
		<< ", argEtot_proposed=" << argEtot_proposed
		<< ", beta=" << beta << std::endl; */

	if(beta_dE_mJ < 0) {
		return 1;
	} else {
		return exp(-1.0 * beta_dE_mJ );
	}

}

/** 
 * Chooses whether to accept a sample or not based on a probability 
 **/
bool HMCSampler::acceptSample() {

	// Return value
	bool acceptSampleReturn = false;

	// Local vars
	SimTK::Real hereE_o  = 0.0;
	SimTK::Real hereE_n  = 0.0;
	SimTK::Real here_lnj = 0.0;

	// The decision tree sets the value of internal variable acc
	if((this->alwaysAccept == true) || 
		((this->nofSamples < equilNofRounds) &&
			(world->ownWorldIndex == 0))) // DIRTY FIX
	{ // Empty sampler

		acceptSampleReturn = true;

	}else{ // Markov-Chain Monte Carlo

		const SimTK::Real rand_no = uniformRealDistribution(randomEngine);

		SimTK::Real prob = 0.0;

		if(DistortOpt == 0){
			hereE_o = etot_o;
			hereE_n = etot_n;
			here_lnj = 0.0;

		}else if(DistortOpt > 0){

			if(useFixman){
				hereE_o = pe_o + ke_prop_nma6 + fix_o;
				hereE_n = pe_n + ke_n_nma6    + fix_n;
				here_lnj = 0.0;

			}else{
				hereE_o = pe_o + ke_prop_nma6;
				hereE_n = pe_n + ke_n_nma6;
				here_lnj = 0.0;				

			}

		}else if(DistortOpt < 0){

			hereE_o = pe_o + fix_o;
			hereE_n = pe_n + fix_n;
			here_lnj = getDistortJacobianDetLog();
			// here_lnj = std::log(std::abs(this->QScaleFactor))

		}

		prob = MetropolisHastings(hereE_o, hereE_n, here_lnj);

		acceptSampleReturn = (rand_no < prob);

	}

	return acceptSampleReturn;
	
}

/**
 * Checks is there are any sudden jumps in potential energy which usually
 * indicate a distortion in the system
*/
bool HMCSampler::checkDistortionBasedOnE(SimTK::Real deltaPE)
{

	// Set an energy limit for the potential energy difference in kT
	SimTK::Real energyLimit = beta * 100 * ndofs;
	
	// Apply
	if((beta * deltaPE) > energyLimit){
		std::cout << "[WARNING] Reduced PE GT " << energyLimit << " "
			<< deltaPE << "." << std::endl;
		return false;
	}else{
		return true;
	}

}

/**
 * Print everything after proposal generation
*/
void HMCSampler::Print(const SimTK::State& someState,
	bool isTheSampleValid, bool isTheSampleAccepted)
{

	// Set precision
	std::cout << std::setprecision(10) << std::fixed;

	// Do we distort? // TODO remove
	std::cout << "DISTORT_OPTION " << this->DistortOpt << std::endl;

	// Print all the energy terms first
	PrintDetailedEnergyInfo(someState);

	//Print acc
	if(isTheSampleAccepted){
		std::cout << "\tacc";
	}else{
		std::cout << "\trej";
	}
	
	// Print simulation type
	if(this->alwaysAccept == true){
		std::cout << "\t (MD)";
	}else{
		std::cout << "\t (MH)";
	}

	// Print validity
	if(isTheSampleValid){
		std::cout << " \n";
	}else{
		std::cout << " invalid\n";
	}
	
}

/**
 * The main function that generates a sample
*/
bool HMCSampler::sample_iteration(SimTK::State& someState)
{

	// Flag if anything goes wrong during simulation
	bool validated = true;

	// Store old configuration
	storeOldPotentialEnergies(someState);

	// PROPOSE
	// Generate a trial move in the stochastic chain
	validated = propose(someState) && validated;

	// --- invalid --- //
	if ( !validated ){

				// Set status
				setAcc(false);

				// RESTORE
				restore(someState);

				// Deal with adaptive data
				storeAdaptiveData(someState); // PrintAdaptiveData();
				
				// Print
				Print(someState, validated, getAcc());

	// --- valid --- //
	}else{

				// --- accept --- //
				if(acceptSample()){

					// Set status
					setAcc(true);

					// Print
					Print(someState, validated, getAcc());

					// UPDATE
					update(someState);

					// Deal with adaptive data
					storeAdaptiveData(someState); // PrintAdaptiveData();

				// --- reject --- //
				}else{

					// Set status
					setAcc(false);

					// Print
					Print(someState, validated, getAcc());

					// RESTORE
					restore(someState);

					// Deal with adaptive data
					storeAdaptiveData(someState); // PrintAdaptiveData();
				}

	}

	// Increase the sample counter and return
	++nofSamples;

	// Return accepted
	return getAcc();

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
	// TODO THIS IS BAD, FIX IT
	//updateQBuffer(someState);
	// pushCoordinatesInR(someState);
	// pushVelocitiesInRdot(someState);
}

void HMCSampler::PrintAdaptiveData(void)
{
	// Calculate MSD and RRdot to adapt the integration length
	std::cout << std::setprecision(10) << std::fixed;
	std::cout << "\tMSD= " << calculateMSD() 
		<< ", RRdot= " << calculateRRdot() << std::endl;

}


///////////////////////////////////////////////////////
// RESTORE
///////////////////////////////////////////////////////

void HMCSampler::restoreConfiguration(
	SimTK::State& someState)
{

	if(integratorName == IntegratorName::OMMVV){
		OMM_restoreConfiguration(someState);
		OMM_To_Simbody_setAtomsLocations(someState); // COMPLETE
		
	}else{
		assignConfFromSetTVector(someState);
	}

	proposeExceptionCaught = false;
}

/** Restore energies */
void HMCSampler::restoreEnergies(void){

	// Set final energies to the precalculated old ones
	pe_set = pe_o;
	ke_set = ke_o;
	fix_set = fix_o;
	logSineSqrGamma2_set = logSineSqrGamma2_o;

	// Set the final total energy
	etot_set = pe_set + fix_set + ke_o + logSineSqrGamma2_set;
}

/** Restore */
void HMCSampler::restore(SimTK::State& someState)
{

	// Restore old configuration
	restoreConfiguration(someState);

	// Restore energies
	restoreEnergies();

	// Update acceptance rate buffer
	acceptedStepsBuffer.push_back(0);
	acceptedStepsBuffer.pop_front();
}

///////////////////////////////////////////////////////
// UPDATE
///////////////////////////////////////////////////////

/** Update energies */
void HMCSampler::updateEnergies(void)
{
	// Store new energies
	pe_set = pe_n;
	ke_set = ke_n;
	fix_set = fix_n;
	logSineSqrGamma2_set = logSineSqrGamma2_n;

	// Set final total energies
	etot_set = pe_set + fix_set + ke_set + logSineSqrGamma2_set;
}

/** Update **/
void HMCSampler::update(SimTK::State& someState)
{

	if(this->integratorName == IntegratorName::OMMVV){
		// Update Simbody too
		OMM_To_Simbody_setAtomsLocations(someState);
	}
	
	// Store final configuration and energy
	// Store new configuration
	storeSimbodyConfiguration_XFMs(someState);

	// Set the final energies to the new ones
	updateEnergies();
	
	// Acceptance rate buffer
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

const bool& Sampler::getAcc(void) const
{
	return this->acc;
}

bool& Sampler::updAcc(void)
{
	return this->acc;
}

void Sampler::setAcc(bool accArg)
{
	this->acc = accArg;
}

void HMCSampler::setOMMmass(SimTK::DuMM::NonbondAtomIndex nax, SimTK::Real mass) {
	// set if using openmm integrator
	if(integratorName == IntegratorName::OMMVV) {
		dumm->setOpenMMparticleMass(nax, mass);
	}
}

void HMCSampler::setNonequilibriumParameters(int distort, int work, int flow) {
	DistortOpt = distort;
	WorkOpt = work;
	FlowOpt = flow;
}

int HMCSampler::getDistortOption() const {
	return DistortOpt;
}

void HMCSampler::setGuidanceHamiltonian(SimTK::Real boostTemperature, int boostMDSteps) {
	setBoostTemperature(boostTemperature); // used for OpenMM and other minor stuff
	setBoostMDSteps(boostMDSteps); // not used
}



// ===========================================================================
// ===========================================================================
// ZMatrix BAT
// ===========================================================================
// ===========================================================================
/*!
 * <!--	Getter for the variableBATs map -->
*/
const std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>>&
HMCSampler::getSubZMatrixBATs() const {
	return subZMatrixBATs;
}

/*!
 * <!--	Getter for the variableBATs map -->
*/
std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>>&
HMCSampler::updSubZMatrixBATs() {
	return subZMatrixBATs;
}

/*!
 * <!--	Print BAT -->
*/
void HMCSampler::PrintSubZMatrixBAT() {
 
	scout("HMCSampler::PrintSubZMatrixBAT ") << eol;

    // Iterate over each key-value pair in the map
    for (const auto& pair : subZMatrixBATs) {

        // Print the MobilizedBodyIndex
        std::cout << "MobilizedBodyIndex: " << pair.first <<" ";

        // Print the vector of BAT values
        std::cout << "BAT Values: ";
        for (const auto& value : pair.second) {
            std::cout << value << " ";
        }

        std::cout << std::endl;
    }
}

/*!
 * <!--	 -->
*/
void
HMCSampler::calcSubZMatrixBATDeviations(
	SimTK::State& someState)
{

	SimTK::Real N = nofSamples + 1;

	// Calculate BAT means
	if(nofSamples == 0){
		
		// Iterate bats
		size_t bati = 0;
    	for (const auto& pair : subZMatrixBATs) {

			subZMatrixBATMeans.insert(pair);
			subZMatrixBATDiffs.insert(std::make_pair(pair.first, std::vector<SimTK::Real> {0, 0, 0}));

			bati++;
		}

	}else{

		// Iterate bats
		size_t bati = 0;
    	for (const auto& pair : subZMatrixBATs) {

			std::vector<SimTK::Real>& BAT      = subZMatrixBATs.at(pair.first);
			std::vector<SimTK::Real>& BATmeans = subZMatrixBATMeans.at(pair.first);
			std::vector<SimTK::Real>& BATdiffs = subZMatrixBATDiffs.at(pair.first);

			// Get BAT means
			SimTK::Real N_1_over_N = (N - 1.0) / N;
			SimTK::Real Ninv = 1.0 / N;
			
			BATmeans[0] = (N_1_over_N * BATmeans[0]) + (Ninv * BAT[0]); 
			BATmeans[1] = (N_1_over_N * BATmeans[1]) + (Ninv * BAT[1]); 
			BATmeans[2] = (N_1_over_N * BATmeans[2]) + (Ninv * BAT[2]); 

			// Get BAT deviations
			BATdiffs[0] = BAT[0] - BATmeans[0];
			BATdiffs[1] = BAT[1] - BATmeans[1];
			BATdiffs[2] = BAT[2] - BATmeans[2];

			bati++;
		}

	}

}

/*!
 * <!--	Print BATs -->
*/
void
HMCSampler::PrintSubZMatrixBATDeviations(
	SimTK::State& someState)
{

	scout("HMCSampler::PrintBATDeviations: ") << eol;

	// Iterate bats
	size_t bati = 0;
	for (const auto& pair : subZMatrixBATs) {

		std::vector<SimTK::Real>& BAT      = subZMatrixBATs.at(pair.first);
		std::vector<SimTK::Real>& BATmeans = subZMatrixBATMeans.at(pair.first);
		std::vector<SimTK::Real>& BATdiffs = subZMatrixBATDiffs.at(pair.first);

		SimTK::MobilizedBodyIndex mbx = pair.first;
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);


		// Print
		scout("PrintBATDeviations mbx qix BATs BATmeans BATdiffs ")
			<< pair.first <<" " << mobod.getFirstQIndex(someState) <<" | "
			<< std::setprecision(12)
			<< BAT[0] <<" " << BAT[1] <<" " << BAT[2] <<" "
			<< BATmeans[0] <<" " << BATmeans[1] <<" " << BATmeans[2] <<" "
			<< BATdiffs[0] <<" " << BATdiffs[1] <<" " << BATdiffs[2] <<" "
			<< eol;

		bati++;


		// 
		for(int qCnt = mobod.getFirstQIndex(someState);
		qCnt < mobod.getFirstQIndex(someState) + mobod.getNumQ(someState);
		qCnt++){
			scout("PrintBATDeviations stateq ") << qCnt <<" " << someState.getQ()[qCnt] << eol;
			scout("PrintBATDeviations mododq ") << mobod.getOneQ(someState, 0) << eol;
		}

	} // every variableBAT


}

/*!
 * <!--	Scale BATs -->
*/
SimTK::Real
HMCSampler::scaleSubZMatrixBATDeviations(
	SimTK::State& someState,
	SimTK::Real scalingFactor,
	std::vector<int> BATOrder,
	std::vector<SimTK::Real> BATSign
	)
{

	SimTK::Real scaleJacobian = 0.0;

	scout("scaleBATDeviations state before ") << std::setprecision(12) << someState.getQ() << eol;

	// Iterate bats
	size_t bati = 0;
	for (const auto& pair : subZMatrixBATs) {

		std::vector<SimTK::Real>& BAT      = subZMatrixBATs.at(pair.first);
		std::vector<SimTK::Real>& BATmeans = subZMatrixBATMeans.at(pair.first);
		std::vector<SimTK::Real>& BATdiffs = subZMatrixBATDiffs.at(pair.first);

		SimTK::MobilizedBodyIndex mbx = pair.first;
		SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

		// scale
		int mobodQCnt = 0;

		for(int qCnt = mobod.getFirstQIndex(someState);
		qCnt < mobod.getFirstQIndex(someState) + mobod.getNumQ(someState);
		qCnt++){

			int rearrMobodQCnt = BATOrder[mobodQCnt];

			if( !(std::isnan(BATdiffs[rearrMobodQCnt])) ){

				SimTK::Real& qEntry = someState.updQ()[qCnt];

				//qEntry += (BATdiffs[rearrMobodQCnt] * (scalingFactor - 1.0));



				// if(nofSamples % 2 == 0){
				//     qEntry += 0.01;
				// }else{
				//     qEntry -= 0.01;
				// }
				qEntry += 0.01;


				qEntry *= BATSign[rearrMobodQCnt];


				scaleJacobian += std::log( (BAT[rearrMobodQCnt] + BATdiffs[rearrMobodQCnt]) / (BAT[rearrMobodQCnt]) );

				scout("scaleBATDeviations ")
					<< "mbx " << mbx << " qCnt "	<< qCnt <<" " 
					<< "BAT[mbx]["<< rearrMobodQCnt << "] " << BAT[rearrMobodQCnt] <<" "
					<< "BATmeans[mbx]["<< rearrMobodQCnt << "] " << BATmeans[rearrMobodQCnt] <<" "
					<< "BATdiffs[mbx]["<< rearrMobodQCnt << "] " << BATdiffs[rearrMobodQCnt] <<" "
					<< "scalingFactor " << (BATdiffs[rearrMobodQCnt] * (scalingFactor - 1.0)) <<" "
					<< " qEntry " << qEntry 
					<< eol;

			}
			
			mobodQCnt++;
		} // every mobod q

		bati++;
	}

	system->realize(someState, SimTK::Stage::Dynamics);
	scout("scaleBATDeviations state after") << std::setprecision(12) << someState.getQ() << eol;

	//return someState;
	return scaleJacobian;

}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ZMatrix BAT
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


