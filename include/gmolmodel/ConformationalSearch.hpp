#ifndef __MONTECARLOSAMPLER_HPP__
#define __MONTECARLOSAMPLER_HPP__

#include "Robo.hpp"
#include "Sampler.hpp"

//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
//#include <boost/random/normal_distribution.hpp>

class Topology;
class IState;

class ConformationalSearch : public Sampler
{
public:

    // Constructor
    ConformationalSearch(SimTK::CompoundSystem *argCompoundSystem,
                      SimTK::SimbodyMatterSubsystem *argMatter,
                      SimTK::Compound *argResidue,
                      SimTK::DuMMForceFieldSubsystem *argDumm,
                      SimTK::GeneralForceSubsystem *forces,
                      SimTK::TimeStepper *argTimeStepper) ;

    // Destructor
    virtual ~ConformationalSearch();

    // Simulation temperature related
    SimTK::Real getTemperature(void);
    void setTemperature(SimTK::Real);

    // Set a thermostat (even for MCMC)
    void setThermostat(ThermostatName);
    void setThermostat(std::string);
    void setThermostat(const char *);

    // Get the name of the thermostat
    virtual ThermostatName getThermostat(void);
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

    // Get/set set potential energy
    SimTK::Real getSetPE(void);
    void setSetPE(SimTK::Real argPE);

    // Get/set current potential energy
    SimTK::Real getOldPE(void);
    void setOldPE(SimTK::Real argPE);

    // Set/get Fixman potential
    void setSetFixman(SimTK::Real);
    SimTK::Real getSetFixman(void);

    // Set/get Fixman potential
    void setREP(SimTK::Real);
    SimTK::Real getREP(void);

    // Set/get Fixman potential
    void setOldFixman(SimTK::Real);
    SimTK::Real getOldFixman(void);

    // Evaluate the potential energy at current state
    SimTK::Real getPEFromEvaluator(SimTK::State& someState); 

    // Return true if use Fixman potential
    void useFixmanPotential(void);
    bool isUsingFixmanPotential(void);

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
    bool getAlwaysAccept(void);

    // Is the sampler always accepting the proposed moves
    void setAlwaysAccept(bool);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
    void update(SimTK::State&);

    // Get the number of accpted conformations
    int getAcceptedSteps(void);

protected:
    SimTK::Transform *SetTVector; // Transform matrices
    SimTK::Transform *TVector; // Transform matrices
    SimTK::Real pe_set, pe_o;
    SimTK::Real temperature;
    SimTK::Real RT;

    bool useFixman;
    bool alwaysAccept;
    SimTK::Real fix_set, fix_o, fix_n;
    SimTK::Real residualEmbeddedPotential; // inside rigid bodies if weren't rigid

    int acceptedSteps;
 
    // Random number generators - not sure if I need two
    // Needs testing

    boost::random::mt19937 randomEngine = boost::random::mt19937();

    boost::random::uniform_real_distribution<double> uniformRealDistribution_0_2pi =
        boost::random::uniform_real_distribution<double>(SimTK::Zero, 2*SimTK::Pi);

    boost::random::uniform_real_distribution<double> uniformRealDistribution_mpi_pi =
        boost::random::uniform_real_distribution<double>((-1)*SimTK::Pi, SimTK::Pi);

    boost::random::uniform_real_distribution<double> uniformRealDistribution =
        boost::random::uniform_real_distribution<double>(SimTK::Zero, SimTK::One);

    boost::random::uniform_real_distribution<double> uniformRealDistribution_m1_1 =
        boost::random::uniform_real_distribution<double>((-1)*SimTK::One, SimTK::One);

    boost::normal_distribution<> gaurand = boost::normal_distribution<>(0.0, 1.0);
    //boost::math::normal gaurand;

};

#endif // __MONTECARLOSAMPLER_HPP__

