#ifndef __HAMMONTECARLOSAMPLER_HPP__
#define __HAMMONTECARLOSAMPLER_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include "Sampler.hpp"
#include <thread>

// Just to remove the long syntax requirement
#ifndef pHMC
//#define pHMC(pSampler) dynamic_cast<HMCSampler *>(pSampler)
#define pHMC(pSampler) pSampler
#endif

// Other classes that we need
class Topology;
class IState;

struct RANDOM_CACHE {
	std::normal_distribution<> Gaussian;
	std::mt19937 RandomEngine; // mt19937_64 is faster in our case

	std::function<SimTK::Real()> GenerateGaussian = [this]() mutable {
		return this->Gaussian(this->RandomEngine);
	};

	std::function<void()> FillWithGaussian = [this]() mutable {
		std::generate(V.begin(), V.end(), GenerateGaussian);
	};

	std::future<void> task;

	SimTK::Vector V, SqrtMInvV;
	SimTK::Real sqrtRT;
	int nu = -1;

	RANDOM_CACHE() {
		// Fill all state 19k bits of Mersenne Twister
		std::vector<uint32_t> RandomData(624);
		std::random_device Source;
		std::generate(RandomData.begin(), RandomData.end(), std::ref(Source));
		std::seed_seq SeedSeq(RandomData.begin(), RandomData.end());
		RandomEngine = std::mt19937(SeedSeq);

		// We want a Gaussian distribution
		Gaussian = std::normal_distribution<>(0.0, 1.0);
	}
};

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
	const char *dirname, const char *prefix, int midlength,
	const char *sufix, double aTime);

/** A Generalized Coordiantes Hamiltonian Monte Carlo sampler as described in
J Chem Theory Comput. 2017 Oct 10;13(10):4649-4659. In short it consists
of the following steps:
   1. Initialize velocities from a random normal distribution with a 
	  covariance of kT sqrt(M) where M is the mass matrix tensor.
   2. Propagate the trial trajectory using a symplectic integrator provided
	  by Simbody
   3. An acception-rejection step which includes the Fixman potential 
	  if needed.
Step 1 and 2 are implemented in the fuction propose. Step 3 is implemented
in the function update.
**/
class HMCSampler : virtual public Sampler
{
friend class Context;
public:

	/** Constructor **/
	HMCSampler(World *argWorld,
		SimTK::CompoundSystem *argCompoundSystem,
		SimTK::SimbodyMatterSubsystem *argMatter,
		std::vector<Topology> &argTopologies, 
		SimTK::DuMMForceFieldSubsystem *argDumm,
		SimTK::GeneralForceSubsystem *argForces,
		SimTK::TimeStepper *argTimeStepper);

	/** Destructor **/
	virtual ~HMCSampler();

	// BEGIN MCSAMPLER
	// Get/Set a thermostat (even for MCMC)
	void setThermostat(ThermostatName);
	void setThermostat(std::string);
	void setThermostat(const char *);
	virtual ThermostatName getThermostat(void) const;
 
	void setIntegratorName(IntegratorName);
	const IntegratorName getIntegratorName(void){return integratorName;}
	void setIntegratorName(const std::string integratorName);

 	/* 
	* Compute mathematical, rather than robotic Jacobian.
	* It translates generalized velocities u into Cartesian velocities
	* for each atom on all topologies
	*/ 
	SimTK::Matrix& 
	calcMathJacobian(const SimTK::State& someState,
		SimTK::Matrix& mathJ);

	/*
	* Get the diagonal 3Nx3N matrix containing the atoms masses
	*/
	void getCartesianMassMatrix(const SimTK::State& somestate,
		SimTK::Matrix& M);

	// Return true if use Fixman potential
	void useFixmanPotential(void);
	bool isUsingFixmanPotential(void) const;

	// Compute Fixman potential
	SimTK::Real calcFixman(SimTK::State& someState);

	// Compute Fixman potential numerically
	SimTK::Real calcNumFixman(SimTK::State& someState);

	// Set/get Fixman potential
	void setOldFixman(SimTK::Real);
	SimTK::Real getOldFixman(void) const;

	// Set/get Fixman potential
	void setSetFixman(SimTK::Real);
	SimTK::Real getSetFixman(void) const;

	// Set/get External MBAT contribution potential
	void setOldLogSineSqrGamma2(SimTK::Real);
	SimTK::Real getOldLogSineSqrGamma2(void) const;

	// Set/get External MBAT contribution potential
	void setSetLogSineSqrGamma2(SimTK::Real);
	SimTK::Real getSetLogSineSqrGamma2(void) const;

	// 
	void setProposedLogSineSqrGamma2(SimTK::Real argFixman);

	// Evaluate the potential energy at current state
	SimTK::Real getPEFromEvaluator(const SimTK::State& someState) const;

	// Get/set current potential energy
	SimTK::Real getOldPE(void) const;
	void setOldPE(SimTK::Real argPE);

	// Get/set set potential energy
	SimTK::Real getNewPE(void) const;
	void setNewPE(SimTK::Real argPE); 

	// Get/set set potential energy
	SimTK::Real getSetPE(void) const;
	void setSetPE(SimTK::Real argPE); 

	// Set/get residual embedded potential energy: potential
	// stored inside rigid bodies
	void setREP(SimTK::Real); 
	SimTK::Real getREP(void) const; 

	// Set/get Fixman potential
	void setSetTVector(const SimTK::State& advanced); 
	SimTK::Transform * getSetTVector(void);
	void assignConfFromSetTVector(SimTK::State& advanced);

	// Store/restore the configuration from the internal transforms vector
	// TVector
	void setTVector(const SimTK::State& advanced);
	void setTVector(SimTK::Transform *);
	SimTK::Transform * getTVector(void);

	/** Calculate O(n2) the square root of the mass matrix inverse
	denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
	This is lower triangular **/
	void calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv) const;

	/** Calculate O(n2) the square root of the mass matrix inverse
	denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular **/
	void calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv) const;

	/** Calculate sqrt(M) using Eigen. For debug purposes. **/
	void calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper) const;

	/** Helper function for initialize velocities. Put generalized velocities
	scale factors into a fixed size array to avoid searching for them into a 
	map every time the velocities are intialized **/
	void loadUScaleFactors(SimTK::State& someState);

	/** Seed the random number GenerateGaussian. Set simulation temperature, 
	velocities to desired temperature, variables that store the configuration
	and variables that store the energies, both needed for the 
	acception-rejection step. Also realize velocities and initialize 
	the timestepper. **/
	virtual void initialize(SimTK::State& advanced);

	/** Same as initialize **/
	virtual void reinitialize(SimTK::State& advanced) ;

	/** Get the TimeStepper that manages the integrator **/
	const SimTK::TimeStepper * getTimeStepper();
	SimTK::TimeStepper * updTimeStepper();
	void setTimeStepper(SimTK::TimeStepper * someTimeStepper);

	/** Get/Set the timestep for integration **/
	virtual SimTK::Real getTimestep() const;
	virtual void setTimestep(SimTK::Real);

	/** Get/Set boost temperature **/
	SimTK::Real getBoostTemperature();
	void setBoostTemperature(SimTK::Real);
	void setBoostMDSteps(int);

	/** Store configuration and PE, Fixman potential and logsin gamma squared **/
	virtual void 
		storeOldConfigurationAndPotentialEnergies(
			SimTK::State& someState);

	// Set the method of integration
	void setSampleGenerator(const std::string& samplerNameArg);

	/** Initialize velocities according to the Maxwell-Boltzmann
	distribution.  Coresponds to R operator in LAHMC **/
	void setVelocitiesToZero(SimTK::State& someState);
	virtual void initializeVelocities(SimTK::State& someState);
	void initializeNMAVelocities(SimTK::State& someState);

	/** Store the proposed energies **/
	virtual void calcProposedKineticAndTotalEnergy(SimTK::State& someState);

	/** Apply the L operator **/
	virtual void integrateTrajectory(SimTK::State& someState);
	void integrateVariableTrajectory(SimTK::State& someState);

	/** Integrate trajectory one step at a time to compute quantities instantly **/
	virtual void integrateTrajectoryOneStepAtATime(SimTK::State& someState);

	/** Use stochastic optimization to adapt timestep **/
	virtual void adaptTimestep(SimTK::State& someState);

	/** Adapt Gibbs blocks definitions **/
	void adaptWorldBlocks(SimTK::State& someState);

	/** Store new configuration and energy terms**/
	virtual void calcNewConfigurationAndEnergies(SimTK::State& someState);

	/**  Restore old configuration and energies**/
	virtual void setSetConfigurationToOld(SimTK::State& someState);

	/** Update new configuration and energiees **/
	virtual void setSetConfigurationAndEnergiesToNew(
		SimTK::State& someState);

	/** Restore configuration and set energies to old */
	void setSetConfigurationAndEnergiesToOld(
		SimTK::State& someState);

	/** Metropolis-Hastings acceptance probability **/
	SimTK::Real MHAcceptProbability(SimTK::Real argEtot_proposed,
		SimTK::Real argEtot_n) const;

	/** Accetion rejection step **/
	virtual bool accRejStep(SimTK::State& someState);

	/** Checks if the proposal is valid **/
	bool validateProposal() const;

	/** Chooses whether to accept a sample or not based on a probability **/
	bool acceptSample();

	/*
	* Get Joint type by examining hinge matrix H_FM
	*/
	int getJointTypeFromH(const SimTK::State& someState,
		const SimTK::MobilizedBody& mobod);

	// 
	void setDistortOption(const int& distortOptArg);

	void setQScaleFactor(const SimTK::Real& s);

	/** Shift all the generalized coordinates
	 **/
	void shiftQ(SimTK::State& someState);

	// ELIZA OPENMM FULLY FLEXIBLE INTEGRATION CODE
	void OMM_setTemperature(double HMCBoostTemperature);

	double OMM_calcPotentialEnergy(void);

	double OMM_calcKineticEnergy(void);

	void OMM_calcProposedKineticAndTotalEnergy(void);

	void OMM_integrateTrajectory(SimTK::State&);

	void OMM_calcNewConfigurationAndEnergies(void);

	/** It implements the proposal move in the Hamiltonian Monte Carlo
	algorithm. It essentially propagates the trajectory after it stores
	the configuration and energies. Returns true if the proposal is 
	validatedTODO: break in two functions:initializeVelocities and 
	propagate/integrate **/
	bool proposeEquilibrium(SimTK::State& someState);
	bool proposeNEHMC(SimTK::State& someState);
	bool proposeNMA(SimTK::State& someState);

	/**
	 * Generate a trial move in the chain
	*/
	bool generateProposal(SimTK::State& someState);

	/** Calls setSetConfigurationAndEnergiesToNew **/
	void update(SimTK::State& someState);

	/** Calls setSetConfigurationAndEnergiesToOld **/
	void restore(SimTK::State& someState);

	virtual bool sample_iteration(SimTK::State& someState);

	/**
	*  Add generalized coordinates to a buffer
	*/
	void updateQBuffer(const SimTK::State& someState);

	/** Push Cartesian coordinates into R vector stored in Sampler.
	Return the size of R **/
	std::size_t pushCoordinatesInR(SimTK::State& someState);

	/** Push Cartesian velocities into Rdot vector stored in Sampler.
	Return the size of Rdot **/
	std::size_t pushVelocitiesInRdot(SimTK::State& someState);

	/** Push generalized coordinates into R vector stored in Sampler.
	Return the size of R **/
	std::size_t pushCoordinatesInQ(SimTK::State& someState);

	/** Push generalizedvelocities into Rdot vector stored in Sampler.
	Return the size of Rdot **/
	std::size_t pushVelocitiesInQdot(SimTK::State& someState);

	/** Push generalizedvelocities into Rdot vector stored in Sampler.
	Return the size of Rdot **/
	std::size_t pushVelocitiesInU(SimTK::State& someState);

	void storeAdaptiveData(SimTK::State& someState);

	/** 
	 * Print adaptive data
	*/
	void PrintAdaptiveData(void);

	/** Get the proposed kinetic energy. This is set right  after velocities
	are initialized. **/
	SimTK::Real getOldKE() { return this->ke_o; }
	
	/** Get the stored kinetic energy. This is set rightafter a move is
	accepted. It's a component of the total energy stored. **/
	SimTK::Real getLastAcceptedKE() { return this->ke_set; }
	
	/** Sets the proposed kinetic energy before the proposal. This should be
	set right after the velocities are initialized. **/
	void setOldKE(SimTK::Real);

	/** Stores the accepted kinetic energy. This should be set right after a 
	move is accepted. It's a component of the total energy stored. **/
	void setSetKE(SimTK::Real);

	int getMDStepsPerSample() const;

	void setMDStepsPerSample(int mdStepsPerSample);

	SimTK::Real getMDStepsPerSampleStd() const;
	void setMDStepsPerSampleStd(SimTK::Real mdstd = 0);

	/** Print detailed energy information **/
	void PrintDetailedEnergyInfo(const SimTK::State& someState) const;

	/** Calculate Mean Square Displacement based on stored R vectors **/
	SimTK::Real calculateMSD();

	/** Calculate RRdot based on stored R and Rdot vectors **/
	SimTK::Real calculateRRdot();

	/** **/
	void geomDihedral(SimTK::State& someState);

	/** Load the map of mobods to joint types **/
	//void loadMbx2mobility(SimTK::State& someState);

	/*
 	* Test ground mostly for SOA
 	*/
	void testSOA(SimTK::State& someState);

protected:

	RANDOM_CACHE RandomCache;

	// BEGIN MCSampler
	std::vector<SimTK::Transform> SetTVector; // Transform matrices
	std::vector<SimTK::Transform> TVector; // Transform matrices
	SimTK::Real pe_set = 0.0,
	    pe_o = 0.0,
	    pe_n = 0.0;

	SimTK::Real fix_set = 0.0,
	    fix_o = 0.0,
	    fix_n = 0.0;
	SimTK::Real detmbat_set = 0.0,
	    detmbat_o = 0.0,
	    detmbat_n = 0.0;
	SimTK::Real residualEmbeddedPotential = 0.0; // inside rigid bodies if weren't rigid

	SimTK::Real logSineSqrGamma2_o = 0.0, logSineSqrGamma2_n = 0.0, logSineSqrGamma2_set = 0.0;

	bool useFixman = false;
	//bool alwaysAccept = false;

	int acceptedSteps = 0;
	int acceptedStepsBufferSize = 500;
	std::deque<int> acceptedStepsBuffer;

	int QsBufferSize = 300;
	//std::list<SimTK::Vector> QsBuffer;
	std::deque<SimTK::Real> QsBuffer;

	// Buffer to hold Q means
	std::vector<SimTK::Real> QsMeans;

	SimTK::Real acceptance;
	SimTK::Real prevAcceptance;
	SimTK::Real prevPrevAcceptance;

	//
	int DistortOpt = 0;
	int FlowOpt = 0;
	int WorkOpt = 0;

	bool proposeExceptionCaught;
	// END MCSampler

	std::vector<SimTK::Real> R;
	std::vector<SimTK::Real> Rdot;

	std::vector<SimTK::Real> dR;
	std::vector<SimTK::Real> dRdot;

	// Integration
	IntegratorName integratorName = IntegratorName::EMPTY;

	// Sampling
	int sampleGenerator = 0;

	std::vector<SimTK::Real> UScaleFactors;
	SimTK::Real UScaleFactorsNorm = 0.0;
	std::vector<SimTK::Real> InvUScaleFactors;
	SimTK::Real InvUScaleFactorsNorm = 0.0;

	SimTK::Vector NormedUScaleFactors;
	SimTK::Vector DOFUScaleFactors;

	bMatrix NMARotation;

	SimTK::Real timestep;
	SimTK::Real prevTimestep;
	SimTK::Real prevPrevTimestep;

	SimTK::Real ke_set; // last accepted kinetic energy
	SimTK::Real ke_o; // proposed kinetic energy
	SimTK::Real ke_n; // new kinetic energy
	SimTK::Real etot_set; // stored total energy
	SimTK::Real etot_o; // last accepted total energ (same with stored)
	SimTK::Real etot_n;

	SimTK::Real ke_prop_nma6;
	SimTK::Real ke_n_nma6;

	SimTK::Real work;

	SimTK::Real boostT;
	SimTK::Real boostRT;
	SimTK::Real sqrtBoostRT; // vel init
	SimTK::Real boostBeta;
	SimTK::Real boostFactor;
	SimTK::Real unboostFactor;
	int boostMDSteps;

	int MDStepsPerSample;
	SimTK::Real MDStepsPerSampleStd;

	bool shouldAdaptTimestep;

	SimTK::Real NMAAltSign;

	SimTK::Real QScaleFactor;

};

#endif // __HAMMONTECARLOSAMPLER_HPP__
