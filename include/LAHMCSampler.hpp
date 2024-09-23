#ifndef __LAHMCSAMPLER_HPP__
#define __LAHMCSAMPLER_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include "HMCSampler.hpp"

// Just to remove the long syntax requirement
#ifndef pLAHMC
//#define pLAHMC(pSampler) dynamic_cast<LAHMCSampler *>(pSampler)
#define pLAHMC(pSampler) pSampler
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
    LAHMCSampler(World &argWorld,
		SimTK::CompoundSystem &argCompoundSystem,
		SimTK::SimbodyMatterSubsystem &argMatter,
		std::vector<Topology> &argTopologies, 
		SimTK::DuMMForceFieldSubsystem &argDumm,
		SimTK::GeneralForceSubsystem &argForces,
		SimTK::TimeStepper &argTimeStepper,
      unsigned int Kext);

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
    virtual bool initialize(SimTK::State& advanced);

    /** Same as initialize **/
    virtual bool reinitialize(SimTK::State& advanced) ;

    /** Get the TimeStepper that manages the integrator **/
    //const SimTK::TimeStepper * getTimeStepper();
    //SimTK::TimeStepper * updTimeStepper();
    //void setTimeStepper(SimTK::TimeStepper * someTimeStepper);

    /** Get/Set the timestep for integration **/
    //virtual SimTK::Real getTimestep() const;
    //virtual void setTimestep(SimTK::Real, bool adaptive);

    /** Get/Set boost temperature **/
    //SimTK::Real getBoostTemperature();
    //void setBoostTemperature(SimTK::Real);
    //void setBoostMDSteps(int);

    /** Store configuration and PE, Fixman potential and logsin gamma squared **/
    void storeOldPotentialEnergies(SimTK::State& someState);

    /** Initialize velocities according to the Maxwell-Boltzmann   
    distribution.  Coresponds to R operator in LAHMC **/
    void initializeVelocities(SimTK::State& someState);

    /** Store the proposed energies **/
    void calcProposedKineticAndTotalEnergy(SimTK::State& someState);

    /** Apply the L operator **/
    bool integrateTrajectory(SimTK::State& someState);

    /** Store new configuration and energy terms**/
    void calcNewConfigurationAndEnergies(SimTK::State& someState, int k);

    /*** Set C matrices entries to 0 ***/
    void resetCMatrices();

    /*** Return a submatrix of M with lesser rows and cols from the end ***/
    SimTK::Matrix extractFromTop(SimTK::Matrix M, int rowCut, int colCut);

    /*** Copies different sizes matrices entries and avoids Simbody resize ***/
    void injectFromTop(const std::vector<SimTK::Real>& src, std::vector<SimTK::Real>& dest);
    void injectFromTop(const SimTK::Matrix& src, SimTK::Matrix& dest);

    /*** Return the antidiagonal transpose of a matrix ***/
    SimTK::Matrix reverseMatrix(SimTK::Matrix M);

    /*** Set C matrices entry ***/
    void setCEntry(int i, int j, SimTK::Real entry);

    /*** Set C matrices entry ***/
    void setCtauEntry(int i, int j, SimTK::Real entry);

    /*** Convert C indeces to Ctau indeces ***/
    void C_to_Ctau_Indeces(int C_i, int C_j, int &Ctau_i, int &Ctau_j, int currSize);

    /*** Convert Ctau indeces to C indeces ***/
    void Ctau_to_C_Indeces(int Ctau_i, int Ctau_j, int &C_i, int &C_j, int currSize);

    /*** Transition probability from state o to state n ***/
     SimTK::Real MHAcceptProbability(SimTK::State& someState, SimTK::Real E_o, SimTK::Real E_n);

    /*** Compute cumulative transition probabilities ***/
    SimTK::Matrix& leap_prob_recurse_hard(SimTK::State& someState, std::vector<SimTK::Real> etot_ns_loc, SimTK::Matrix& CC_loc);

    void setSetConfigurationAndEnergiesToNew(SimTK::State& someState);

    /** It implements the proposal move in the Hamiltonian Monte Carlo
    algorithm. It essentially propagates the trajectory after it stores
    the configuration and energies. TODO: break in two functions:
    initializeVelocities and propagate/integrate **/
    bool propose(SimTK::State& someState);

    /** Main function that contains all the 3 steps of HMC.
    Implements the acception-rejection step and sets the state of the 
    compound to the appropriate conformation wether it accepted or not. **/
    void update(SimTK::State& someState);

    /** Accetion rejection step **/
    virtual bool accRejStep(SimTK::State& someState);

    bool sample_iteration(SimTK::State& someState);

    /** Modifies Q randomly
     **/
    //void perturbQ(SimTK::State& someState);

    /** Get the proposed kinetic energy. This is set right  after velocities
    are initialized. **/
    SimTK::Real getProposedKE() { return this->ke_o; }
      
    /** Get the stored kinetic energy. This is set rightafter a move is
    accepted. It's a component of the total energy stored. **/
    SimTK::Real getLastAcceptedKE() { return this->ke_set; }
    
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

    int K; // # times to apply the L operator

    // Energies of the new states
    std::vector<SimTK::Real> pe_ns;
    std::vector<SimTK::Real> fix_ns;
    std::vector<SimTK::Real> logSineSqrGamma2_ns;
    std::vector<SimTK::Real> ke_ns;
    std::vector<SimTK::Real> etot_ns;

    // Anti-diagonal transpose Ctau as denoted by 2007 Vasily Golyshev, 
    // Jan Stienstra. The same as Python C[:0:, :0:]
    SimTK::Matrix C, Ctau; // Cumulative probability matrices
    SimTK::Matrix P; // Ant-diagonal identity matrix

    SimTK::Matrix CC;


};

#endif // __LAHMCSAMPLER_HPP__

