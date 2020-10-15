#ifndef __MONTECARLOSAMPLER_HPP__
#define __MONTECARLOSAMPLER_HPP__

#include "Robo.hpp"
#include "Sampler.hpp"

//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
//#include <boost/random/normal_distribution.hpp>


// Just to remove the long syntax requirement
#ifndef pMC
//#define pMC(pSampler) dynamic_cast<MonteCarloSampler *>(pSampler)
#define pMC(pSampler) pSampler
#endif

// Other classes that we need
class Topology;
class IState;

class MonteCarloSampler : virtual public Sampler
{
    friend class Context;
public:

    // Constructor
    MonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                      SimTK::SimbodyMatterSubsystem *argMatter,

                      //SimTK::Compound *argResidue,
                      std::vector<Topology *>& argTopologies,

                      SimTK::DuMMForceFieldSubsystem *argDumm,
                      SimTK::GeneralForceSubsystem *forces,
                      SimTK::TimeStepper *argTimeStepper) ;

    // Destructor
    virtual ~MonteCarloSampler();

    // Simulation temperature related
    //SimTK::Real getTemperature(void);
    //void setTemperature(SimTK::Real);

    // Set a thermostat (even for MCMC)
    void setThermostat(ThermostatName);
    void setThermostat(std::string);
    void setThermostat(const char *);

    // Get the name of the thermostat
    virtual ThermostatName getThermostat(void) const;
    /** Seed the random number generator. Set simulation temperature,
    variables that store the configuration
    and variables that store the energies, both needed for the
    acception-rejection step. Also realize velocities and initialize
    the timestepper. **/
    virtual void initialize(SimTK::State& advanced, SimTK::Real argTemperature, bool argUseFixman = true) ;

    /** Same as initialize **/
    virtual void reinitialize(SimTK::State& advanced, SimTK::Real argTemperature) ;

    // Store/restore the configuration from the set nternal transforms vector
    // TVector
    void setSetTVector(const SimTK::State& advanced);
    SimTK::Transform * getSetTVector(void);
    void assignConfFromSetTVector(SimTK::State& advanced);

    // Store/restore the configuration from the internal transforms vector
    // TVector
    void setTVector(const SimTK::State& advanced);
    void setTVector(SimTK::Transform *);
    SimTK::Transform * getTVector(void);
    void assignConfFromTVector(SimTK::State& advanced);

    // Assign a random conformation
    void propose(SimTK::State& advanced);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
    void update(SimTK::State&);

    // Get/set set potential energy
    SimTK::Real getSetPE(void) const;
    void setSetPE(SimTK::Real argPE);

    // Get/set current potential energy
    SimTK::Real getOldPE(void) const;
    void setOldPE(SimTK::Real argPE);

    // Set/get Fixman potential
    void setSetFixman(SimTK::Real);
    SimTK::Real getSetFixman(void) const;

    // Set/get Fixman potential
    void setREP(SimTK::Real);
    SimTK::Real getREP(void) const;

    // Set/get Fixman potential
    void setOldFixman(SimTK::Real);
    SimTK::Real getOldFixman(void) const;



    // Set/get External MBAT contribution potential
    void setSetLogSineSqrGamma2(SimTK::Real);
    SimTK::Real getSetLogSineSqrGamma2(void) const;

    // Set/get External MBAT contribution potential
    void setOldLogSineSqrGamma2(SimTK::Real);
    SimTK::Real getOldLogSineSqrGamma2(void) const;

    //
    void setProposedLogSineSqrGamma2(SimTK::Real argFixman);
    SimTK::Real getProposedLogSineSqrGamma2(void) const;



    // Set/get Fixman potential
    void setProposedFixman(SimTK::Real);
    SimTK::Real getProposedFixman(void) const;

    // Evaluate the potential energy at current state
    SimTK::Real getPEFromEvaluator(SimTK::State& someState); 

    // Return true if use Fixman potential
    void useFixmanPotential(void);
    bool isUsingFixmanPotential(void) const;

    // Compute Fixman potential
    SimTK::Real calcFixman(SimTK::State& someState);

    // Compute Fixman potential numerically
    SimTK::Real calcNumFixman(SimTK::State& someState);

    // Compute mass matrix determinant numerically
    // Not to be confused with the Fixman potential
    SimTK::Real calcNumDetM(SimTK::State& someState);

    // Send configuration to an external evaluator
    void sendConfToEvaluator(void);
    
    // Is the sampler always accepting the proposed moves
    bool getAlwaysAccept(void) const;

    // Is the sampler always accepting the proposed moves
    void setAlwaysAccept(bool);

    // Get the number of accpted conformations
    int getAcceptedSteps(void) const;

protected:
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
    bool alwaysAccept = false;

    int acceptedSteps = 0;
    int acceptedStepsBufferSize = 30;
    std::list<int> acceptedStepsBuffer;

    float acceptance;
    float prevAcceptance;
    float prevPrevAcceptance;

    bool proposeExceptionCaught;
 
};

#endif // __MONTECARLOSAMPLER_HPP__

