#ifndef __SAMPLER_HPP__
#define __SAMPLER_HPP__

#include "Robo.hpp"

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

	// Is the sampler always accepting the proposed moves
	virtual bool getAlwaysAccept(void) const;

	// Is the sampler always accepting the proposed moves
	virtual void setAlwaysAccept(bool);

	// Getter / setter for macroscopic temperature and RT
	// virtual void setTemperature(SimTK::Real) = 0; // RE
	SimTK::Real getTemperature() const;
	void setTemperature(SimTK::Real temperature); // Also sets RT and beta

	SimTK::Real getRT() const;
	void setBeta(SimTK::Real argBeta);
	SimTK::Real getBeta() const;

	// Just for checking
	void checkAtomStationsThroughDumm(void);

	/** Load the map of mobods to joint types **/
	//void loadMbx2mobility(SimTK::State& someState); // SAFE
	void loadMbx2mobility(int whichWorld); // DANGER

	/** Returns the number of samples extracted so far. **/
	int getNofSamples();

	// Get set the seed
	int64_t getSeed() const;
	void setSeed(int64_t);

	// Draws one sample from vonMises distribution with concentration k
	// The algorithm is taken from 1979 Best, page 155
	SimTK::Real vonMises(SimTK::Real miu, SimTK::Real k);

	// Draws X from von Mises-Fisher distribution with concentration 
	// parameter k TODO: reference to the algorithm
	// X vector has the dimensions of the ndofs
	std::vector<double>& vonMisesFisher(std::vector<double>& X, double k);

	// Draws from chi distribution
	double chi(void);

	/** Generate a random number. **/
	SimTK::Real generateRandomNumber(GmolRandDistributionType);

	virtual void shiftQ ( SimTK::State& someState, SimTK::Real scalingFactor, int numIgnoredQs ) = 0;

	/** Propose a move **/
	virtual bool propose(SimTK::State& someState) = 0;
	//virtual eval() = 0;
	virtual void update(SimTK::State& someState) = 0;

	// For debugging purposes
	void PrintSimbodyStateCache(SimTK::State& someState);

public: 
	// Classes we need to access
	World *world;
	const SimTK::System *system;
	SimTK::CompoundSystem *compoundSystem;
	SimTK::SimbodyMatterSubsystem *matter;

	//SimTK::Compound *rootTopology;
	Topology *rootTopology;

	std::vector<Topology>& topologies;
	std::size_t natoms;
	std::size_t ndofs;

	// Total mass of the system
	SimTK::Real totalMass;

	/** Joint types **/
	std::map< SimTK::MobilizedBodyIndex, SimTK::BondMobility::Mobility> mbx2mobility;
	std::map< SimTK::QIndex, JointType> qIndex2jointType;

	SimTK::DuMMForceFieldSubsystem *dumm;
	SimTK::GeneralForceSubsystem *forces;
	SimTK::TimeStepper *timeStepper;

	// Thermodynamics
	bool alwaysAccept = false;
	ThermostatName thermostat;
	SimTK::Real temperature;
	SimTK::Real RT;
	SimTK::Real beta;

	// Sampling
	int nofSamples;
	int64_t seed;
	bool acc;

	// Random number generators - not sure if I need two
	using RANDOM_ENGINE = std::mt19937_64; // mt19937_64
	RANDOM_ENGINE randomEngine;

	// Use this to initialize randomEngine
	using RANDOM_ENGINE_INIT = std::minstd_rand; // ranlux48
	RANDOM_ENGINE_INIT randomEngineInit;

	using RANDOM_ENGINE_INIT_RESULT_TYPE = uint32_t; 

	std::uniform_real_distribution<SimTK::Real> uniformRealDistribution_0_2pi =
		    std::uniform_real_distribution<SimTK::Real>(SimTK::Zero, 2*SimTK::Pi);

	std::uniform_real_distribution<SimTK::Real> uniformRealDistribution_mpi_pi =
		    std::uniform_real_distribution<SimTK::Real>((-1)*SimTK::Pi, SimTK::Pi);

	std::uniform_real_distribution<SimTK::Real> uniformRealDistribution =
		    std::uniform_real_distribution<SimTK::Real>(SimTK::Zero, SimTK::One);

	std::uniform_real_distribution<SimTK::Real> uniformRealDistribution_m1_1 =
		    std::uniform_real_distribution<SimTK::Real>((-1)*SimTK::One, SimTK::One);

	// Gaussian random number distribution
	std::normal_distribution<> gaurand = std::normal_distribution<>(0.0, 1.0);

	std::gamma_distribution<double> gammarand = std::gamma_distribution<double>(1, 2);

 };


#endif // __SAMPLER_HPP__

