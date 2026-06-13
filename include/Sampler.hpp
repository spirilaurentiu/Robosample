#pragma once

#include "Random.hpp"
#include "bgeneral.hpp"

class Topology;
class World;

// TODO: implement getCompoundSystem
class Sampler {
    public:
    // Constructor
    Sampler(World& argWorld,
            SimTK::CompoundSystem& argCompoundSystem,
            SimTK::SimbodyMatterSubsystem& argMatter,
            Span<Topology> argTopologies,
            SimTK::DuMMForceFieldSubsystem& argDumm,
            SimTK::GeneralForceSubsystem& argForces,
            SimTK::TimeStepper& argTimeStepper);

    // Destructor
    virtual ~Sampler();

    // Compute mass matrix determinant (O(n))
    // TODO Move
    auto calcMassDeterminant(const SimTK::State& state) -> SimTK::Real;
    auto calcMassDeterminant(SimTK::State& state) -> SimTK::Real;

    // Set / reset variables needed at the beginning of a simulation
    void initialize(SimTK::State& someState);
    void reinitialize(SimTK::State& someState);

    // Is the sampler always accepting the proposed moves
    [[nodiscard]] virtual auto getAlwaysAccept() const -> bool;

    // Is the sampler always accepting the proposed moves
    virtual void setAlwaysAccept(bool);

    // Getter / setter for macroscopic temperature and RT
    // virtual void setTemperature(SimTK::Real) = 0; // RE
    [[nodiscard]] auto getTemperature() const -> SimTK::Real;
    void setTemperature(SimTK::Real temperature); // Also sets RT and beta

    [[nodiscard]] auto getRT() const -> SimTK::Real;
    void setBeta(SimTK::Real argBeta);
    [[nodiscard]] auto getBeta() const -> SimTK::Real;

    // Just for checking
    void checkAtomStationsThroughDumm();

    /** Load the map of mobods to joint types **/
    // void loadMbx2mobility(SimTK::State& someState); // SAFE
    void loadMbx2mobility(int whichWorld);

    /** Returns the number of samples extracted so far. **/
    [[nodiscard]] auto getNofSamples() -> int;

    // Get set the seed
    void setSeed(Random32& seeder);

    // Draws one sample from vonMises distribution with concentration k
    // The algorithm is taken from 1979 Best, page 155
    [[nodiscard]] auto generateVonMisesSample(SimTK::Real miu, SimTK::Real k) -> SimTK::Real;

    // Draws X from von Mises-Fisher distribution with concentration
    // parameter k TODO: reference to the algorithm
    // X vector has the dimensions of the ndofs
    void generateVonMisesFisherSample(std::vector<SimTK::Real>& X, SimTK::Real k);

    // Draws from chi distribution
    [[nodiscard]] auto generateChiSample() -> SimTK::Real;

    /** Generate a random number. **/
    [[nodiscard]] auto generateRandomNumber(GmolRandDistributionType) -> SimTK::Real;

    /** Generate random rotation quaternion */
    [[nodiscard]] auto generateRandomQuaternion() -> SimTK::Quaternion;

    /**
     * Take a variable and transforming according to some distribution
     * //TODO revise param1 and param2
     */
    auto convoluteVariable(SimTK::Real& var,
                           std::string distrib = "alternateInverse",
                           SimTK::Real param1 = 1.0,
                           SimTK::Real param2 = 0.0,
                           SimTK::Real param3 = 0.0) -> SimTK::Real;
    auto convoluteVariable(std::vector<SimTK::Real>& vvar,
                           std::string distrib = "alternateInverse",
                           SimTK::Real param1 = 1.0,
                           SimTK::Real param2 = 0.0) -> SimTK::Real;
    auto calcDeformationPotential(SimTK::Real& var,
                                  std::string distrib = "alternateInverse",
                                  SimTK::Real param1 = 1.0,
                                  SimTK::Real param2 = 0.0) -> SimTK::Real;
    virtual auto setQToScaleBendStretchStdev(SimTK::State& someState, std::vector<SimTK::Real>& scaleFactors)
        -> SimTK::Real = 0;

    // virtual void setIntegratorName(IntegratorName) = 0;
    // virtual void setIntegratorName(std::string integratorName) = 0;

    [[nodiscard]] virtual const bool& getAcc() const;
    virtual auto updAcc() -> bool&;
    virtual void setAcc(bool);

    // For debugging purposes
    void PrintSimbodyStateCache(SimTK::State& someState);

    // Classes we need to access
    std::reference_wrapper<World> world;
    std::reference_wrapper<const SimTK::System> system;
    std::reference_wrapper<SimTK::CompoundSystem> compoundSystem;
    std::reference_wrapper<SimTK::SimbodyMatterSubsystem> matter;
    std::reference_wrapper<SimTK::DuMMForceFieldSubsystem> dumm;
    std::reference_wrapper<SimTK::GeneralForceSubsystem> forces;
    std::reference_wrapper<SimTK::TimeStepper> timeStepper;

    // SimTK::Compound *rootTopology;
    Topology* rootTopology;

    Span<Topology> topologies;
    int natoms = 0;
    int numDegreesOfFreedom = 0;
    int acceptedSteps = 0;

    // Total mass of the system
    SimTK::Real totalMass = 0;

    /** Joint types **/
    std::map<SimTK::MobilizedBodyIndex, SimTK::BondMobility::Mobility> mbx2mobility;
    std::map<SimTK::QIndex, JointType> qIndex2jointType;

    // Thermodynamics
    bool alwaysAccept = false;
    ThermostatName thermostat = ThermostatName::None;
    SimTK::Real temperature = SimTK::NaN;
    SimTK::Real RT = SimTK::NaN;
    SimTK::Real beta = SimTK::NaN;

    // Sampling
    int nofSamples = 0;
    uint32_t seed = 0;
    bool acc = false;

    int numSamples_period = 0;
    int numAccepted_period = 0;

    // Random number generators
    Random64 randomEngine;
    RANDOM_CACHE RandomCache;

    std::uniform_real_distribution<SimTK::Real> uniformRealDistribution_0_2pi =
        std::uniform_real_distribution<SimTK::Real>(SimTK::Zero, 2 * SimTK::Pi);

    std::uniform_real_distribution<SimTK::Real> uniformRealDistribution_mpi_pi =
        std::uniform_real_distribution<SimTK::Real>((-1) * SimTK::Pi, SimTK::Pi);

    std::uniform_real_distribution<SimTK::Real> uniformRealDistribution =
        std::uniform_real_distribution<SimTK::Real>(SimTK::Zero, SimTK::One);

    std::uniform_real_distribution<SimTK::Real> uniformRealDistribution_m1_1 =
        std::uniform_real_distribution<SimTK::Real>((-1) * SimTK::One, SimTK::One);

    // Gaussian random number distribution
    std::normal_distribution<> gaurand = std::normal_distribution<>(0.0, 1.0);

    std::gamma_distribution<SimTK::Real> gammarand = std::gamma_distribution<SimTK::Real>(1, 2);

    std::lognormal_distribution<SimTK::Real> lognormal = std::lognormal_distribution<SimTK::Real>(0.0, 1.0);
};
