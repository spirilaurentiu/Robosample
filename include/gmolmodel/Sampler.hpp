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
class World;

// TODO: implement getCompoundSystem
class Sampler
{
public:
	// Constructor
	Sampler(World *argWorld,
		SimTK::CompoundSystem *argCompoundSystem,
		SimTK::SimbodyMatterSubsystem *argMatter,
		//SimTK::Compound *argResidue,
		std::vector<Topology> &argTopologies,
		SimTK::DuMMForceFieldSubsystem *argDumm,
		SimTK::GeneralForceSubsystem *argForces,
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
	void setTemperature(SimTK::Real temperature); // Also sets RT and beta

	SimTK::Real getRT() const;
	void setBeta(SimTK::Real argBeta);
	SimTK::Real getBeta() const;

	/** Load the map of mobods to joint types **/
	void loadMbx2mobility(SimTK::State& someState);

	/** Returns the number of samples extracted so far. **/
	int getNofSamples(void);

	// Get set the seed
	uint32_t getSeed(void) const;
	void setSeed(uint32_t);

	/** Generate a random number. **/
	SimTK::Real generateRandomNumber(GmolRandDistributionType);

	/** Propose a move **/
	virtual bool propose(SimTK::State& someState) = 0;
	//virtual eval() = 0;
	virtual void update(SimTK::State& someState) = 0;

	// For debugging purposes
	void PrintSimbodyStateCache(SimTK::State& someState);

public:
	World *world;
	const SimTK::System *system;
	SimTK::CompoundSystem *compoundSystem;
	SimTK::SimbodyMatterSubsystem *matter;

	//SimTK::Compound *rootTopology;
	Topology *rootTopology;

	std::vector<Topology>& topologies;
	std::size_t natoms;
	std::size_t ndofs;

	/** Joint types **/
	std::map< SimTK::MobilizedBodyIndex, SimTK::BondMobility::Mobility> mbx2mobility;
	std::map< SimTK::QIndex, JointType> qIndex2jointType;

	SimTK::DuMMForceFieldSubsystem *dumm;
	SimTK::GeneralForceSubsystem *forces;
	SimTK::TimeStepper *timeStepper;

	// Thermodynamics
	ThermostatName thermostat;
	SimTK::Real temperature;
	SimTK::Real RT;
	SimTK::Real beta;

	// Sampling
	int nofSamples;
	uint32_t seed;
	bool acc;

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

