#ifndef __HAMMONTECARLOSAMPLER_HPP__
#define __HAMMONTECARLOSAMPLER_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include "MonteCarloSampler.hpp"

// Just to remove the long syntax requirement
#ifndef pHMC
//#define pHMC(pSampler) dynamic_cast<HMCSampler *>(pSampler)
#define pHMC(pSampler) pSampler
#endif

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
	const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime);

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
class HMCSampler : virtual public MonteCarloSampler
{
friend class Context;
public:

	/** Constructor **/
	HMCSampler(SimTK::CompoundSystem *argCompoundSystem
			,SimTK::SimbodyMatterSubsystem *argMatter

			//,SimTK::Compound *argResidue
			,std::vector<Topology *>& argTopologies

			,SimTK::DuMMForceFieldSubsystem *argDumm
			,SimTK::GeneralForceSubsystem *forces
			,SimTK::TimeStepper *argTimeStepper
			);

	/** Destructor **/
	virtual ~HMCSampler();

	/** Calculate O(n2) the square root of the mass matrix inverse
	denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
	This is lower triangular **/
	void calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv);

	/** Calculate O(n2) the square root of the mass matrix inverse
	denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular **/
	void calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv);

	/** Calculate sqrt(M) using Eigen. For debug purposes. **/
	void calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper);

	/** Seed the random number generator. Set simulation temperature, 
	velocities to desired temperature, variables that store the configuration
	and variables that store the energies, both needed for the 
	acception-rejection step. Also realize velocities and initialize 
	the timestepper. **/
	virtual void initialize(SimTK::State& advanced);

	/** Same as initialize **/
	virtual void reinitialize(SimTK::State& advanced) ;

	/** Get the TimeStepper that manages the integrator **/
	const SimTK::TimeStepper * getTimeStepper(void);
	SimTK::TimeStepper * updTimeStepper(void);
	void setTimeStepper(SimTK::TimeStepper * someTimeStepper);

	/** Get/Set the timestep for integration **/
	virtual float getTimestep(void);
	virtual void setTimestep(float);

	/** Get/Set boost temperature **/
	SimTK::Real getBoostTemperature(void);
	void setBoostTemperature(SimTK::Real);
	void setBoostMDSteps(int);

	/** Store configuration and PE, Fixman potential and logsin gamma squared **/
	virtual void storeOldConfigurationAndPotentialEnergies(SimTK::State& someState);

	/** Initialize velocities according to the Maxwell-Boltzmann
	distribution.  Coresponds to R operator in LAHMC **/
	virtual void initializeVelocities(SimTK::State& someState);

	/** Store the proposed energies **/
	virtual void calcProposedKineticAndTotalEnergy(SimTK::State& someState);

	/** Apply the L operator **/
	virtual void integrateTrajectory(SimTK::State& someState);

	/** Integrate trajectory one step at a time to compute quantities instantly **/
	virtual void integrateTrajectoryOneStepAtATime(SimTK::State& someState);

	/** Use stochastic optimization to adapt timestep **/
	virtual void adaptTimestep(SimTK::State& someState);

	/** Adapt Gibbs blocks definitions **/
	void adaptWorldBlocks(SimTK::State& someState);

	/** Store new configuration and energy terms**/
	virtual void calcNewConfigurationAndEnergies(SimTK::State& someState);

	/**  Restore old configuration and energies**/
	virtual void setSetConfigurationAndEnergiesToOld(SimTK::State& someState);

	/** Update new configuration and energiees **/
	virtual void setSetConfigurationAndEnergiesToNew(SimTK::State& someState);

	/** Metropolis-Hastings acceptance probability **/
	virtual SimTK::Real MHAcceptProbability(SimTK::State& someState, 
	SimTK::Real argEtot_proposed, SimTK::Real argEtot_n);

	/** Accetion rejection step **/
	virtual bool accRejStep(SimTK::State& someState);

	/** It implements the proposal move in the Hamiltonian Monte Carlo
	algorithm. It essentially propagates the trajectory after it stores
	the configuration and energies. TODO: break in two functions:
	initializeVelocities and propagate/integrate **/
	void propose(SimTK::State& someState);

	/** Main function that contains all the 3 steps of HMC.
	Implements the acception-rejection step and sets the state of the 
	compound to the appropriate conformation wether it accepted or not. **/
	void update(SimTK::State& someState);

	virtual bool sample_iteration(SimTK::State& someState);

	/** Push Cartesian coordinates into R vector stored in Sampler.
	Return the size of R **/
	int pushCoordinatesInR(SimTK::State& someState);

	/** Push Cartesian velocities into Rdot vector stored in Sampler.
	Return the size of Rdot **/
	int pushVelocitiesInRdot(SimTK::State& someState);

	/** Push generalized coordinates into R vector stored in Sampler.
	Return the size of R **/
	int pushCoordinatesInQ(SimTK::State& someState);

	/** Push generalizedvelocities into Rdot vector stored in Sampler.
	Return the size of Rdot **/
	int pushVelocitiesInQdot(SimTK::State& someState);

	/** Push generalizedvelocities into Rdot vector stored in Sampler.
	Return the size of Rdot **/
	int pushVelocitiesInU(SimTK::State& someState);

	/** Modifies Q randomly
	 **/
	void perturbQ(SimTK::State& someState);

	/** Get the proposed kinetic energy. This is set right  after velocities
	are initialized. **/
	SimTK::Real getProposedKE(void) { return this->ke_proposed; }
	
	/** Get the stored kinetic energy. This is set rightafter a move is
	accepted. It's a component of the total energy stored. **/
	SimTK::Real getLastAcceptedKE(void) { return this->ke_lastAccepted; }
	
	/** Sets the proposed kinetic energy before the proposal. This should be
	set right after the velocities are initialized. **/
	void setProposedKE(SimTK::Real);

	/** Stores the accepted kinetic energy. This should be set right after a 
	move is accepted. It's a component of the total energy stored. **/
	void setLastAcceptedKE(SimTK::Real);

	int getMDStepsPerSample() const;

	void setMDStepsPerSample(int mdStepsPerSample);

	/** Print detailed energy information **/
	void PrintDetailedEnergyInfo(SimTK::State& someState);

	/** Calculate Mean Square Displacement based on stored R vectors **/
	SimTK::Real calculateMSD(void);

	/** Calculate RRdot based on stored R and Rdot vectors **/
	SimTK::Real calculateRRdot(void);

	/** **/
	void geomDihedral(void);

	/** Load the map of mobods to joint types **/
        //void loadMbx2mobility(SimTK::State& someState);

protected:

	float timestep;
	float prevTimestep;
	float prevPrevTimestep;
	bool shouldAdaptTimestep;

	SimTK::Real ke_lastAccepted; // last accepted kinetic energy
	SimTK::Real ke_proposed; // proposed kinetic energy
	SimTK::Real ke_n; // new kinetic energy
	SimTK::Real etot_set; // stored total energy
	SimTK::Real etot_proposed; // last accepted total energ (same with stored)
	SimTK::Real etot_n;

	SimTK::Real boostT;
	SimTK::Real boostFactor;
	SimTK::Real unboostFactor;
	int boostMDSteps;

	int MDStepsPerSample;

	std::vector<SimTK::Real> R;

	std::vector<SimTK::Real> Rdot;

	std::vector<SimTK::Real> dR;
	std::vector<SimTK::Real> dRdot;

};

#endif // __HAMMONTECARLOSAMPLER_HPP__

