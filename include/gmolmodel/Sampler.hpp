#ifndef __SAMPLER_HPP__
#define __SAMPLER_HPP__

#include "Robo.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

#ifndef PRINT_BUFFER_SIZE
#define PRINT_BUFFER_SIZE 4096
#endif 

class Topology;

class Sampler
{
public:
    // Constructor
    Sampler(SimTK::CompoundSystem *argCompoundSystem,
            SimTK::SimbodyMatterSubsystem *argMatter,
            SimTK::Compound *argResidue,
            SimTK::DuMMForceFieldSubsystem *argDumm,
            SimTK::GeneralForceSubsystem *forces,
            SimTK::TimeStepper *argTimeStepper);


    // Destructor
    virtual ~Sampler();

    // Compute mass matrix determinant (O(n))
    // TODO Move
    SimTK::Real calcMassDeterminant(const SimTK::State&);
    SimTK::Real calcMassDeterminant(SimTK::State&);

    // Set / reset variables needed at the beginning of a simulation
    void initialize(SimTK::State& someState);
    void reinitialize(SimTK::State& someState);

    // Getter / setter for macroscopic temperature and RT
    // virtual void setTemperature(SimTK::Real) = 0; // RE
    SimTK::Real getTemperature() const;
    void setTemperature(SimTK::Real temperature);

    SimTK::Real getRT() const;

    /** Returns the number of samples extracted so far. **/
    int getNofSamples(void);

    // Get set the seed
    unsigned long long int getSeed(void);
    void setSeed(unsigned long long int);

    /** Generate a random number. **/
    SimTK::Real generateRandomNumber(GmolRandDistributionType);

    /** Propose a move **/
    virtual void propose(SimTK::State& someState) = 0;
    //virtual eval() = 0;
    virtual void update(SimTK::State& someState) = 0;

    // For debugging purposes
    void PrintSimbodyStateCache(SimTK::State& someState);

public:
    const SimTK::System *system;
    SimTK::CompoundSystem *compoundSystem;
    SimTK::SimbodyMatterSubsystem *matter;
    SimTK::Compound *residue;
    SimTK::DuMMForceFieldSubsystem *dumm;
    SimTK::GeneralForceSubsystem *forces;
    SimTK::TimeStepper *timeStepper;

    // Thermodynamics
    ThermostatName thermostat;
    SimTK::Real temperature;
    SimTK::Real RT;

    // Sampling
    int nofSamples;
    unsigned long long int seed;

    // Random number generators - not sure if I need two
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

 };


#endif // __SAMPLER_HPP__

