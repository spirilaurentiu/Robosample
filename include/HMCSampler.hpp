#ifndef __HAMMONTECARLOSAMPLER_HPP__
#define __HAMMONTECARLOSAMPLER_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

/** @file
This defines the HMCSampler class, which (TODO)

Topology 0:                                                      :       
          :                                                      :
Position 0:                                                      :
          :          ┌─────────────────────────┐                 :
                     │                         │                 : X_o
                     │      REINITIALIZE       │                 : pe_o, pe_set
                     │                         │                 : fix_o, fix_set
                     └────────────┬────────────┘                 : 
                                  │                              :
                                  │                              :
                            ┌─────▼─────┐                        :
                            │  PROPOSE  │                        :     
                            └─────┬─────┘                        :
                                  │                              :
Dynamics 0                        │                              :
                                  │                              :
                             ┌────▼──────┐       NO              :
                             │ VALIDATE  │────────────────       :
                             └────┬──────┘               |       :
                                  │                      |       :
                                  │  YES                 |       :
                                  │                      |       :
                              ┌───▼─────┐                |       :
                              │ ACCEPT  │                |       :
                              └────┬────┘                |       :
                                   │                     |       :
                        YES        │       NO            |       :
             ┌────-────────────────┴─────────────────────┬       :
             |                                           │       :
        ┌────▼────-                                 ┌────▼─────┐ :          
        │ UPDATE  │                                 │ RESTORE  │ :
        └─────────┘                                 └──────────┘ :
         Dynamics 0                                  Position 0  :



The class should theoretically track the following vars:

pe ke    fix lnsingamma lnJ TVec fX_PF gX_BM

One iteration must include:

1. reinitialize
	1.1. calc initial vars      (pe_o         fix_o lnsingamma_o lnJ_o TVec_o fX_PF_o gX_BM_o)
2. proposal
    2.1. initialize velocities  (pe_o ke_prop fix_o lnsingamma_o lnJ_o TVec_o fX_PF_o gX_BM_o)
	2.2. integrate trajectories (pe_n ke_n    fix_n lnsingamma_n lnJ_n TVec_n fX_PF_n gX_BM_n)
3. acc-rej step                 (pe_s ke_s    fix_s lnsingamma_s lnJ_s TVec_s fX_PF_s gX_BM_s) = 
    3.1. accept                 (pe_n ke_n    fix_n lnsingamma_n lnJ_n TVec_n fX_PF_n gX_BM_n) OR
	3.2. reject                 (pe_o ke_prop fix_o lnsingamma_o lnJ_o TVec_o fX_PF_o gX_BM_o)

**/

//#include "Context.hpp"
#include "Sampler.hpp"
#include "TaskSpace.hpp"
#include <thread>

// Just to remove the long syntax requirement
#ifndef pHMC
//#define pHMC(pSampler) dynamic_cast<HMCSampler *>(pSampler)
#define pHMC(pSampler) pSampler
#endif

// Other classes that we need
class Topology;
class IState;
class Context;

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
	HMCSampler(World &argWorld,
		SimTK::CompoundSystem &argCompoundSystem,
		SimTK::SimbodyMatterSubsystem &argMatter,
		std::vector<Topology> &argTopologies, 
		SimTK::DuMMForceFieldSubsystem &argDumm,
		SimTK::GeneralForceSubsystem &argForces,
		SimTK::TimeStepper &argTimeStepper)
		//:Qmeans(nullptr), Qdiffs(nullptr), Qstds(nullptr)
;

	/** Destructor **/
	virtual ~HMCSampler();

	/** ===============================
	 * RANDOM NUMBERS
		=============================== */

	// Uniform distribution number generator
	SimTK::Real uniformRealDistributionRandTrunc(
		SimTK::Real L, SimTK::Real R);

	// Uniform distribution PDF
	SimTK::Real uniformRealDistributionPDFTrunc(
		SimTK::Real X, SimTK::Real L, SimTK::Real R);

	// Uniform distribution CDF
	SimTK::Real uniformRealDistributionCDFTrunc(
		SimTK::Real X, SimTK::Real L, SimTK::Real R);
		


	// BEGIN MCSAMPLER
	// Get/Set a thermostat (even for MCMC)
	void setThermostat(ThermostatName);
	void setThermostat(std::string);
	void setThermostat(const char *);
	virtual ThermostatName getThermostat(void) const;
 
	void setIntegratorName(IntegratorName);
	const IntegratorName getIntegratorName(void){return integratorName;}

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
	void setNewFixman(SimTK::Real);
	SimTK::Real getNewFixman(void) const;

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

	// Get/set Jacobians
	SimTK::Real getDistortJacobianDetLog(void) const;
	void setDistortJacobianDetLog(SimTK::Real argJ);

	// Set/get residual embedded potential energy: potential
	// stored inside rigid bodies
	void setREP(SimTK::Real); 
	SimTK::Real getREP(void) const; 

	// Set/get Fixman potential
	void storeSimbodyConfiguration_XFMs(const SimTK::State& advanced); 
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

	/** Get the TimeStepper that manages the integrator **/
	const SimTK::TimeStepper * getTimeStepper();
	SimTK::TimeStepper * updTimeStepper();
	void setTimeStepper(SimTK::TimeStepper * someTimeStepper);

	/** Get/Set the timestep for integration **/
	virtual SimTK::Real getTimestep() const;
	virtual void setTimestep(SimTK::Real ts, bool adaptive);

	/** Get/Set boost temperature **/
	SimTK::Real getBoostTemperature();
	void setBoostTemperature(SimTK::Real);
	void setBoostMDSteps(int);

	/** Store configuration and PE, Fixman potential and logsin gamma squared **/
	virtual void 
		storeOldPotentialEnergies(
			SimTK::State& someState);

	void storeOldAndSetKineticAndTotalEnergies(SimTK::State& someState);

	// Set the method of integration
	void setAcceptRejectMode(AcceptRejectMode sampleGeneratorArg);

	void perturbPositions(SimTK::State& someState, PositionsPerturbMethod);

	/** Set velocities to zero.  **/
	void setVelocitiesToZero(SimTK::State& someState);

	/** Set velocities according to the Maxwell-Boltzmann
	distribution.  **/
	void setVelocitiesToGaussian(SimTK::State& someState);

	virtual void perturbVelocities(SimTK::State& someState,
		VelocitiesPerturbMethod VPM = VelocitiesPerturbMethod::TO_T);

	void setVelocitiesToNMA(SimTK::State& someState);

	void perturbForces(SimTK::State& someState,
		ForcesPerturbMethod FPM);

	/** Store the proposed energies **/
	virtual void calcProposedKineticAndTotalEnergyOld(SimTK::State& someState);

	/** Apply the L operator **/
	virtual void integrateTrajectory(SimTK::State& someState);
	void integrateVariableTrajectory(SimTK::State& someState);

	/** Integrate trajectory one step at a time to compute quantities instantly **/
	virtual void integrateTrajectoryOneStepAtATime(SimTK::State& someState);

	/** BOUND WALK */
	void integrateTrajectory_Bounded(SimTK::State& someState);

	/** BOUND HMC */
	void integrateTrajectory_BoundHMC(SimTK::State& someState);

	/** Integrate trajectory using task space forces */
	void integrateTrajectory_TaskSpace(SimTK::State& someState);

	/** Use stochastic optimization to adapt timestep **/
	virtual void adaptTimestep(SimTK::State& someState);

	/** Adapt Gibbs blocks definitions **/
	void adaptWorldBlocks(SimTK::State& someState);

	/** Cartesian Fixman potential */
	SimTK::Real CartesianFixmanPotential(void);

	/** Store new configuration and energy terms**/
	virtual void calcNewEnergies(SimTK::State& someState);

	/**  Restore old configuration and energies**/
	virtual void restoreConfiguration(SimTK::State& someState);

	/** Update new configuration and energiees **/
	virtual void setSetConfigurationAndEnergiesToNew(
		SimTK::State& someState);

	/** Metropolis-Hastings acceptance probability **/
	SimTK::Real MetropolisHastings(
		SimTK::Real argEtot_o,
		SimTK::Real argEtot_n,
		SimTK::Real lnJ) const;

	/** Metropolis-Hastings acceptance probability **/
	SimTK::Real MetropolisHastings(
		SimTK::Real argEtot_proposed,
		SimTK::Real argEtot_n,
		SimTK::Real transProb,
		SimTK::Real invTransProb,
		SimTK::Real lnJ) const;

	/** Accetion rejection step **/
	//virtual bool accRejStep(SimTK::State& someState);

	/** Checks if the proposal is valid **/
	bool checkExceptionsAndEnergiesForNAN();

	/** Chooses whether to accept a sample or not based on a probability **/
	bool acceptSample();

	/*
	* Get Joint type by examining hinge matrix H_FM
	*/
	int getJointTypeFromH(const SimTK::State& someState,
		const SimTK::MobilizedBody& mobod);

	/** Set simulation temperature, 
	velocities to desired temperature, variables that store the configuration
	and variables that store the energies, both needed for the 
	acception-rejection step. Also realize velocities and initialize 
	the timestepper. **/
	//virtual void initialize(SimTK::State& advanced);

	/** Same as initialize **/
	virtual bool initialize(SimTK::State& advanced) ;
	virtual bool reinitialize(SimTK::State& advanced, std::stringstream& samplerOutStream) ;

	void PrintInitialParams(void);
	void getMsg_Header(std::stringstream& ss);
	void getMsg_InitialParams(std::stringstream& ss);
	void getMsg_EnergyDetails(
		std::stringstream& ss,
		const SimTK::State& someState,
		bool isTheSampleValid,
		bool isTheSampleAccepted
	);

	// ELIZA OPENMM FULLY FLEXIBLE INTEGRATION CODE
	void OMM_setDuMMTemperature(double HMCBoostTemperature);

	// Update OpenMM position from a Simbody Cartesian world
	void Simbody_To_OMM_setAtomsLocationsCartesian(SimTK::State& someState,
		bool throughDumm = true);


	double OMM_calcPotentialEnergy(void);

	double OMM_calcKineticEnergy(void);

	void OMM_storeOMMConfiguration_X(const std::vector<OpenMM::Vec3>& positions);

	void OMM_To_Simbody_setAtomsLocations(SimTK::State& someState);

	void OMM_PrintLocations(void);


	void OMM_integrateTrajectory(SimTK::State&);

	//void OMM_calcNewEnergies(void);

	void OMM_restoreConfiguration(SimTK::State& someState);

	void setSphereRadius(float argSphereRadius);

	///////////////////////////////////////////////////////
	// PROPOSE
	///////////////////////////////////////////////////////

	/** Returns the 'how' argument of perturbPositions */
	PositionsPerturbMethod positionsPerturbMethod(void);

	VelocitiesPerturbMethod velocitiesPerturbMethod(void);

	/** Returns the 'how' argument of perturbVelocities */
	ForcesPerturbMethod forcesPerturbMethod(void);

	// Docking funcitons
	void teleport(SimTK::State& someState);

	// Perturb Q, QDot or QDotDot
	void perturb_Q_QDot_QDotDot(SimTK::State& someState);

	/** It implements the proposal move in the Hamiltonian Monte Carlo
	algorithm. It essentially propagates the trajectory after it stores
	the configuration and energies. Returns true if the proposal is 
	validatedTODO: break in two functions:initializeVelocities and 
	propagate/integrate **/
	bool propose(SimTK::State& someState);

	///////////////////////////////////////////////////////
	// UPDATE
	///////////////////////////////////////////////////////

	/** Update energies */
	void updateEnergies(void);

	/** Calls setSetConfigurationAndEnergiesToNew **/
	void update(SimTK::State& someState);

	///////////////////////////////////////////////////////
	// RESTORE
	///////////////////////////////////////////////////////

	/** Restore energies */
	void restoreEnergies(void);
	/** Calls setSetConfigurationAndEnergiesToOld **/
	void restore(SimTK::State& someState);

	/**
	 * Checks is there are any sudden jumps in potential energy which usually
	 * indicate a distortion in the system
	*/
	bool checkDistortionBasedOnE(SimTK::Real deltaPE);

	virtual bool sample_iteration(SimTK::State& someState,
		std::stringstream& samplerOutStream);

	/**
	 * Print everything after proposal generation
	*/
	void Print(const SimTK::State& someState,
		bool isTheSampleValid, bool isTheSampleAccepted);

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
	SimTK::Real getNewKE() { return this->ke_n; }

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

	void setGuidanceHamiltonian(SimTK::Real boostTemperature, int boostMDSteps);
	
	//////////////////////////////////
	//////       Scaling        //////
	//////////////////////////////////

	void setNonequilibriumParameters(int distort, int work, int flow);
	int getDistortOption() const;

	// Are we performing work by modifying Q
	const int getDistortOpt();
	
	void setDistortOption(const int& distortOptArg);

	//------------------------------------------------------------------------------
    /** @name Scaling Q Directly
	*/

    /**@{**/

	/**@}**/

	//------------------------------------------------------------------------------
    /** @name Scaling Q BendStretch
	*/

    /**@{**/

	const SimTK::Real& getBendStretchStdevScaleFactor(void);
	void setBendStretchStdevScaleFactor(const SimTK::Real& s);

	/** Shift all the generalized coordinates and
	 * return the scale factors of angles and bonds
	 **/
	SimTK::Real setQToScaleBendStretch(SimTK::State& someState,
		std::vector<SimTK::Real>& scaleFactors);

	SimTK::Real setQToShiftBendStretchStdev(SimTK::State& someState,
		std::vector<SimTK::Real>& scaleFactors);

	/** Shift all the generalized coordinates and
	 * return the scale factors of angles and bonds
	 **/
	SimTK::Real setQToScaleBendStretchStdev(SimTK::State& someState,
		std::vector<SimTK::Real>& scaleFactors);

	void setQToScaleBendStretchStdev_Old(SimTK::State& someState,
		std::vector<SimTK::Real>& scaleFactors);

	/**
	 * Get the log of the Jacobian of a bond-angle stretch
	*/
	SimTK::Real calcBendStretchJacobianDetLog(SimTK::State& someState,
		std::vector<SimTK::Real> scaleFactors,
		unsigned int startFromBody = 0);

	/**@}**/
	// WORK Q PERTURB BEND STRETCH --------------------------------------------

	//------------------------------------------------------------------------------
    /** @name Z Matrix and BAT functions
	*/

    /**@{**/

	void PrintSubZMatrixBAT();

	void setSubZMatrixBATStats(
		std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATmeans,
		std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATdiffs,
		std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATstds,
		std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATstds_Alien
	);

	void
	PrintSubZMatrixBATAndRelated(
		SimTK::State& someState
	);

	SimTK::Real
	scaleSubZMatrixBATDeviations(
		SimTK::State& someState,
		SimTK::Real scalingFactor,
		bool BernoulliTrial = true,
		bool varianceBasedScalingFactor = true,		
		std::vector<int> BATOrder = {1, 0, 2}, 					// bendstretch {1, 0, 2};  spherical {2, 1, 0}
		std::vector<SimTK::Real> BATSign = {1, -1, 1}			// bendstretch {1, -1, 1}; spherical {1, 1, 1}
	);

	void
	updateSubZMatrixBAT(
		SimTK::State& someState,
		std::vector<int> BATOrder = {1, 0, 2},
		std::vector<SimTK::Real> BATSign = {1, -1, 1}
	);

	SimTK::Real calcBATJacobianDetLog(
		SimTK::State& someState,
		SimTK::BondMobility::Mobility bondMobility,
		std::vector<int> BATOrder = {1, 0, 2},
		std::vector<SimTK::Real> BATSign = {1, -1, 1}
	);

	// Calculate sub determinant of MBAT 
	SimTK::Real calcSubMBATDetLog(
		SimTK::State& someState);

	// Put in protected
	std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>&> subZMatrixBATs_ref;

	std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>> subZMatrixBATMeans;
	std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>> subZMatrixBATDiffs;
	std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>> subZMatrixBATVars;
	std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>> subZMatrixBATVars_Alien;

    // Updater getter for the map
    const std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>&>& getSubZMatrixBATsRef() const {
        return subZMatrixBATs_ref;
    }
	    std::map<SimTK::MobilizedBodyIndex, std::vector<SimTK::Real>&>& updSubZMatrixBATsRef() {
        return subZMatrixBATs_ref;
    }

	/**@}**/
	// BAT --------------------------------------------------------------------

	void set_dBMps(std::vector<SimTK::Real>& QArg){ prev_dBMps = &QArg; }

	void setPreviousQs(std::vector<SimTK::Real>& QArg){ previousQs = &QArg; }
	void setQmeans(std::vector<SimTK::Real>& QArg){ Qmeans = &QArg; }
	void setQdiffs(std::vector<SimTK::Real>& QArg){ Qdiffs = &QArg; }
	void setQvars(std::vector<SimTK::Real>& QArg){ Qvars = &QArg; }

	double calcMobodsMBAT(SimTK::State& someState);
	double studyBATScale(SimTK::State& someState);


protected:


	// Buffers to hold Q statistics
	std::vector<SimTK::Real>* prev_dBMps = nullptr;	
	std::vector<SimTK::Real>* previousQs = nullptr;
	std::vector<SimTK::Real>* Qmeans = nullptr;
	std::vector<SimTK::Real>* Qdiffs = nullptr;
	std::vector<SimTK::Real>* Qvars = nullptr;


	int equilNofRounds = 0;

	//bool NAN_TO_INF(SimTK::Real& someNumber);

	std::vector<SimTK::Vec3> omm_locations;
	std::vector<SimTK::Vec3> omm_locations_old;

	// BEGIN MCSampler
	std::vector<SimTK::Transform> SetTVector; // Transform matrices
	std::vector<SimTK::Transform> TVector; // Transform matrices
	SimTK::Real pe_init = 0.0,
		pe_set = 0.0,
	    pe_o = 0.0,
	    pe_n = 0.0;

	SimTK::Real fix_set = 0.0,
	    fix_o = 0.0,
	    fix_n = 0.0;
	SimTK::Real detmbat_set = 0.0,
	    detmbat_o = 0.0,
	    detmbat_n = 0.0;
	SimTK::Real residualEmbeddedPotential = 0.0; // inside rigid bodies if weren't rigid

	SimTK::Real logSineSqrGamma2_o = 0.0,
		logSineSqrGamma2_n = 0.0,
		logSineSqrGamma2_set = 0.0;

	bool useFixman = false;
	//bool alwaysAccept = false;

	int acceptedStepsBufferSize = 50;
	std::deque<int> acceptedStepsBuffer;
	SimTK::Real learningRate = 10e-6;
	SimTK::Real idealAcceptance = 0.651;
	SimTK::Real MDStepsPerSampleStd = 0.5;
	SimTK::Real timestep = std::numeric_limits<SimTK::Real>::max(),
		prevTimestep = std::numeric_limits<SimTK::Real>::max();
	int MDStepsPerSample = std::numeric_limits<int>::max(),
		prevMDStepsPerSample = std::numeric_limits<int>::max();
	bool shouldAdaptTimestep = false;

	int QsBufferSize = 300;
	//std::list<SimTK::Vector> QsBuffer;
	std::deque<SimTK::Real> QsBuffer;

	SimTK::Real acceptance = SimTK::NaN,
		prevAcceptance = SimTK::NaN;

	// Non-equilibrium options
	int DistortOpt = 0;
	int FlowOpt = 0;
	int WorkOpt = 0;

	bool proposeExceptionCaught = false;
	// END MCSampler

	std::vector<SimTK::Real> R;
	std::vector<SimTK::Real> Rdot;

	std::vector<SimTK::Real> dR;
	std::vector<SimTK::Real> dRdot;

	// Integration
	IntegratorName integratorName = IntegratorName::None;

	// For RANDOM_WALK Docking Simulations
	SimTK::Vec3 geometricCenter;
	float sphereRadius = SimTK::NaN;

	// Sampling
	AcceptRejectMode sampleGenerator = AcceptRejectMode::AlwaysAccept;

	std::vector<SimTK::Real> UScaleFactors;
	SimTK::Real UScaleFactorsNorm = 0.0;
	std::vector<SimTK::Real> InvUScaleFactors;
	SimTK::Real InvUScaleFactorsNorm = 0.0;

	SimTK::Vector NormedUScaleFactors;
	SimTK::Vector DOFUScaleFactors;

	bMatrix NMARotation;

	SimTK::Real ke_set; // last accepted kinetic energy
	SimTK::Real ke_o; // proposed kinetic energy
	SimTK::Real ke_n; // new kinetic energy
	SimTK::Real etot_set; // stored total energy
	SimTK::Real etot_o; // last accepted total energ (same with stored)
	SimTK::Real etot_n;

	SimTK::Real ke_prop_nma6;
	SimTK::Real ke_n_nma6;

	SimTK::Real bendStretchJacobianDetLog = 0.0;

	SimTK::Real boostT = SimTK::NaN,
		boostRT = SimTK::NaN,
		sqrtBoostRT = SimTK::NaN, // vel init
		boostBeta = SimTK::NaN;

	SimTK::Real boostKEFactor;
	SimTK::Real unboostKEFactor;

	SimTK::Real boostUFactor;
	SimTK::Real unboostUFactor;
	int boostMDSteps;

	SimTK::Real NMAAltSign = 1.0;
	SimTK::Real QScaleFactor = 1.0;

	// TODO explain what this does
	SimTK::Vector sqrtMInvV;

	// DELETE
	SimTK::Real debug_rand_no = 0.0;

};

#endif // __HAMMONTECARLOSAMPLER_HPP__
