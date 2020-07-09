#ifndef __LAHMCSAMPLER_HPP__
#define __LAHMCSAMPLER_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include "HMCSampler.hpp"

// Just to remove the long syntax requirement
#ifndef pHMC
//#define pHMC(pSampler) dynamic_cast<LAHMCSampler *>(pSampler)
#define pHMC(pSampler) pSampler
#endif

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
class LAHMCSampler : virtual public HMCSampler
{
friend class Context;
public:

    /** Constructor **/
    LAHMCSampler(SimTK::CompoundSystem *argCompoundSystem
                                 ,SimTK::SimbodyMatterSubsystem *argMatter
                                 ,SimTK::Compound *argResidue
                                 ,SimTK::DuMMForceFieldSubsystem *argDumm
                                 ,SimTK::GeneralForceSubsystem *forces
                                 ,SimTK::TimeStepper *argTimeStepper
                                 );

    /** Destructor **/
    virtual ~LAHMCSampler();

    /** Calculate O(n2) the square root of the mass matrix inverse
    denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
    This is lower triangular **/
    //void calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv);

    /** Calculate O(n2) the square root of the mass matrix inverse
    denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular **/
    //void calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv);

    /** Calculate sqrt(M) using Eigen. For debug purposes. **/
    //void calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper);

    /** Seed the random number generator. Set simulation temperature, 
    velocities to desired temperature, variables that store the configuration
    and variables that store the energies, both needed for the 
    acception-rejection step. Also realize velocities and initialize 
    the timestepper. **/
    virtual void initialize(SimTK::State& advanced);

    /** Same as initialize **/
    virtual void reinitialize(SimTK::State& advanced) ;

    /** Get the TimeStepper that manages the integrator **/
    //const SimTK::TimeStepper * getTimeStepper(void);
    //SimTK::TimeStepper * updTimeStepper(void);
    //void setTimeStepper(SimTK::TimeStepper * someTimeStepper);

    /** Get/Set the timestep for integration **/
    //virtual float getTimestep(void);
    //virtual void setTimestep(float);

    /** Get/Set boost temperature **/
    //SimTK::Real getBoostTemperature(void);
    //void setBoostTemperature(SimTK::Real);
    //void setBoostMDSteps(int);

    /** It implements the proposal move in the Hamiltonian Monte Carlo
    algorithm. It essentially propagates the trajectory after it stores
    the configuration and energies. TODO: break in two functions:
    initializeVelocities and propagate/integrate **/
    void propose(SimTK::State& someState);

    /** Main function that contains all the 3 steps of HMC.
    Implements the acception-rejection step and sets the state of the 
    compound to the appropriate conformation wether it accepted or not. **/
    bool update(SimTK::State& someState, SimTK::Real newBeta);

    /** Modifies Q randomly
     **/
    //void perturbQ(SimTK::State& someState);

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

protected:

    //float timestep;
    //bool adaptTimestep;
    //float timestepIncr = 0.0001;
    //SimTK::Real ke_lastAccepted; // last accepted kinetic energy
    //SimTK::Real ke_proposed; // proposed kinetic energy
    //SimTK::Real ke_n; // new kinetic energy
    //SimTK::Real etot_set; // stored total energy
    //SimTK::Real etot_proposed; // last accepted total energ (same with stored)
    //SimTK::Real etot_n;

    //SimTK::Real boostT;
    //SimTK::Real boostFactor;
    //SimTK::Real unboostFactor;
    //int boostMDSteps;

    //int MDStepsPerSample;

};

#endif // __LAHMCSAMPLER_HPP__

