#pragma once

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

// #include "Context.hpp"
#include <cstdint>
#include <thread>

#include "EnergySnapshot.hpp"
#include "NUTS.hpp"
#include "OpenMM.hpp"
#include "Sampler.hpp"
#include "State.h"
#include "Vec3.h"
#include "bgeneral.hpp"
// #include "TaskSpace.hpp"

// Just to remove the long syntax requirement
#ifndef pHMC
// #define pHMC(pSampler) dynamic_cast<HMCSampler *>(pSampler)
#    define pHMC(pSampler) pSampler
#endif

// Other classes that we need
class Topology;
class Context;

void writePdb(SimTK::Compound& c,
              SimTK::State& advanced,
              const char* dirname,
              const char* prefix,
              int midlength,
              const char* sufix,
              SimTK::Real aTime);

/** A Generalized Coordinates Hamiltonian Monte Carlo sampler as described in
J Chem Theory Comput. 2017 Oct 10;13(10):4649-4659. In short it consists
of the following steps:
   1. Initialize velocities from a random normal distribution with a
      covariance of kT sqrt(M) where M is the mass matrix tensor.
   2. Propagate the trial trajectory using a symplectic integrator provided
      by Simbody
   3. An acception-rejection step which includes the Fixman potential
      if needed.
Step 1 and 2 are implemented in the function propose. Step 3 is implemented
in the function update.
**/
class HMCSampler : virtual public Sampler {
    friend class Context;

    public:
    [[nodiscard]] auto getPreviousEnergy() const -> const EnergySnapshot& {
        return previousEnergy;
    }

    [[nodiscard]] auto getCurrentEnergy() const -> const EnergySnapshot& {
        return currentEnergy;
    }

    [[nodiscard]] auto getNumSamples() const -> int {
        return numSamples;
    }

    [[nodiscard]] auto getNumAcceptedSamples() const -> int {
        return numAcceptedSamples;
    }

    [[nodiscard]] auto getAcceptance() const -> SimTK::Real {
        return numSamples > 0 ? static_cast<SimTK::Real>(numAcceptedSamples) / numSamples : SimTK::NaN;
    }

    void resetAcceptance() {
        numSamples = 0;
        numAcceptedSamples = 0;
    }

    [[nodiscard]] auto getTotalEnrgies() const -> const std::vector<SimTK::Real>& {
        return totalEnergiesBuffer;
    }

    [[nodiscard]] auto getNumDegreesOfFreedom() const -> int {
        return numDegreesOfFreedom;
    }

    std::vector<SimTK::Real> UCache, UDotCache, TorqueCache;

    /** Constructor **/
    HMCSampler(World& argWorld,
               SimTK::CompoundSystem& argCompoundSystem,
               SimTK::SimbodyMatterSubsystem& argMatter,
               Span<Topology> argTopologies,
               SimTK::DuMMForceFieldSubsystem& argDumm,
               SimTK::GeneralForceSubsystem& argForces,
               SimTK::TimeStepper& argTimeStepper);
    //: Qmeans(nullptr), Qdiffs(nullptr), Qstds(nullptr)

    /** Destructor **/
    virtual ~HMCSampler();

    /** Same as initialize **/
    virtual void initialize();
    virtual void reinitialize(SimTK::State& state, std::stringstream& samplerOutStream, bool verbose);

    /** ===============================
     * RANDOM NUMBERS
        =============================== */

    // Uniform distribution number generator
    SimTK::Real uniformRealDistributionRandTrunc(SimTK::Real L, SimTK::Real R);

    // Uniform distribution PDF
    SimTK::Real uniformRealDistributionPDFTrunc(SimTK::Real X, SimTK::Real L, SimTK::Real R);

    // Uniform distribution CDF
    SimTK::Real uniformRealDistributionCDFTrunc(SimTK::Real X, SimTK::Real L, SimTK::Real R);


    // BEGIN MCSAMPLER
    // Get/Set a thermostat (even for MCMC)
    void setThermostat(ThermostatName);
    void setThermostat(std::string);
    void setThermostat(const char*);
    [[nodiscard]] virtual auto getThermostat() const -> ThermostatName;

    void setIntegratorType(IntegratorType type);

    [[nodiscard]] auto getIntegratorType() const -> IntegratorType {
        return integratorType;
    }

    [[nodiscard]] auto getIntegratorTypeAsString() const -> std::string {
        switch (integratorType) {
            case IntegratorType::Empty:
                return "Empty";
            case IntegratorType::Verlet:
                return "Simbody Velocity Verlet";
            case IntegratorType::Euler:
                return "Euler";
            case IntegratorType::Euler2:
                return "Euler2";
            case IntegratorType::CPodes:
                return "CPodes";
            case IntegratorType::RungeKutta:
                return "RungeKutta";
            case IntegratorType::RungeKutta2:
                return "RungeKutta2";
            case IntegratorType::RungeKutta3:
                return "RungeKutta3";
            case IntegratorType::RungeKuttaFeldberg:
                return "RungeKuttaFeldberg";
            case IntegratorType::BendStretch:
                return "BendStretch";
            case IntegratorType::OpenMMVelocityVerlet:
                return "OpenMM Velocity Verlet";
            case IntegratorType::BoundWalk:
                return "BoundWalk";
            case IntegratorType::BoundHMC:
                return "BoundHMC";
            case IntegratorType::StationsTask:
                return "StationsTask";
            default:
                return "UnknownIntegratorType";
        }
    }

    void setUseNUTS(bool useNUTS) {
        this->useNUTS = useNUTS;
    }

    /*
     * Compute mathematical, rather than robotic Jacobian.
     * It translates generalized velocities u into Cartesian velocities
     * for each atom on all topologies
     */
    auto calcMathJacobian(const SimTK::State& someState, SimTK::Matrix& mathJ) -> SimTK::Matrix&;

    void PrintUDot(const SimTK::State& someState);
    auto GetUDot(const SimTK::State& someState) -> const SimTK::Vector&;

    /*
     * Get the diagonal 3Nx3N matrix containing the atoms masses
     */
    void getCartesianMassMatrix(const SimTK::State& somestate, SimTK::Matrix& M);

    // Return true if use Fixman potential
    void useFixmanPotential() {
        useFixman = true;
    }
    [[nodiscard]] auto isUsingFixmanPotential() const -> bool {
        return useFixman;
    }

    // Compute Fixman potential
    auto calcFixman(const SimTK::State& someState) -> SimTK::Real;

    // Get/set Jacobians
    [[nodiscard]] auto getDistortJacobianDetLog() const -> SimTK::Real;
    void setDistortJacobianDetLog(SimTK::Real argJ);

    // Set/get residual embedded potential energy: potential
    // stored inside rigid bodies
    void setREP(SimTK::Real);
    [[nodiscard]] auto getREP() const -> SimTK::Real;

    /** Calculate O(n2) the square root of the mass matrix inverse
    denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
    This is lower triangular **/
    void calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv) const;

    /** Calculate O(n2) the square root of the mass matrix inverse
    denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular **/
    void calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv) const;

    /** Helper function for initialize velocities. Put generalized velocities
    scale factors into a fixed size array to avoid searching for them into a
    map every time the velocities are initialized **/
    void loadUScaleFactors(const SimTK::State& someState);

    /** Get/Set the timestep for integration **/
    [[nodiscard]] virtual auto getTimestep() const -> SimTK::Real;
    virtual void setTimestep(SimTK::Real ts, bool adaptive);

    /** Get/Set boost temperature **/
    auto getBoostTemperature() -> SimTK::Real;
    void setBoostTemperature(SimTK::Real);
    void setBoostMDSteps(int);

    // Set the method of integration
    void setAcceptRejectMode(AcceptRejectMode acceptRejectMode);

#pragma region REBAS_TEST

    enum REBAS_MoleculeName_Ix {
        ETHANE,
        ALA1,
        TRPCH
    };
    std::vector<std::string> REBAS_MoleculeNames = {"ETHANE", "ALA1", "TRPCH"};

    std::vector<SimTK::Real>&
    dihedralSegmenter(int nofIntervals, SimTK::Real segHalfDiff, std::vector<SimTK::Real>& segLims);
    int findSegmentIndex(SimTK::Real value, const std::vector<SimTK::Real>& segLims);

    bool REBAS_Scale_Mbx(REBAS_MoleculeName_Ix molName, SimTK::MobilizedBodyIndex mbx);
    void perturbPositions(SimTK::State& someState, PositionsPerturbMethod);

#pragma endregion REBAS_TEST


    /** Set velocities to zero.  **/
    void setVelocitiesToZero(SimTK::State& someState);

    /** Set velocities according to the Maxwell-Boltzmann
    distribution.  **/
    void setVelocitiesToGaussian(SimTK::State& someState);

    virtual void perturbVelocities(SimTK::State& someState,
                                   VelocitiesPerturbMethod VPM = VelocitiesPerturbMethod::ToTemperature);

    void setVelocitiesToNMA(SimTK::State& someState);

    void perturbForces(SimTK::State& someState, ForcesPerturbMethod FPM);

    /** Apply the L operator **/
    void integrateVariableTrajectory(SimTK::State& someState);

    auto isUTurnOpenMM(const PhasePointOpenMM& minus, const PhasePointOpenMM& plus) -> bool;
    auto buildTreeWithOpenMM(std::vector<OpenMM::Vec3>& positions,
                             std::vector<OpenMM::Vec3>& velocities,
                             int depth,
                             NUTSDirection direction,
                             SimTK::Real logU,
                             SimTK::Real H0) -> NUTSNodeOpenMM;
    auto integrateNUTSWithOpenMM(SimTK::State& state, int maxDepth) -> NUTSResult;


    void buildTreeWithSimbody(SimTK::State& state,
                              int depth,
                              NUTSDirection direction,
                              SimTK::Real logU,
                              SimTK::Real H0,
                              NUTSNodeRef& nodeOut,
                              NUTSTrajectoryLog& log,
                              int& fwdLeafCount,
                              int& bwdLeafCount);
    auto integrateNUTSWithSimbody(const SimTK::State& someState, int maxDepth) -> NUTSResult;

    /** Integrate trajectory one step at a time to compute quantities instantly **/
    virtual void integrateTrajectoryOneStepAtATime(SimTK::State& someState);

    /** BOUND WALK */
    void integrateTrajectory_Bounded(SimTK::State& someState);

    /** BOUND HMC */
    SimTK::Vec3 sampleRandomVectorOnSphere(SimTK::Real radius);
    void integrateTrajectory_BoundHMC(SimTK::State& someState);

    /** Integrate trajectory using task space forces */
    void integrateTrajectory_TaskSpace(SimTK::State& someState);

    /** Store new configuration and energy terms**/
    virtual void calcNewEnergies(SimTK::State& someState);

    /** Update new configuration and energiees **/
    virtual void setSetConfigurationAndEnergiesToNew(SimTK::State& someState);

    /** Accetion rejection step **/
    // virtual bool accRejStep(SimTK::State& someState);

    /** Chooses whether to accept a sample or not based on a probability **/
    auto acceptSample(const EnergySnapshot& proposedEnergy, bool shouldPrint) -> bool;

    /*
     * Get Joint type by examining hinge matrix H_FM
     */
    int getJointTypeFromH(const SimTK::State& someState, const SimTK::MobilizedBody& mobod);

    /** Set simulation temperature,
    velocities to desired temperature, variables that store the configuration
    and variables that store the energies, both needed for the
    accept-rejection step. Also realize velocities and initialize
    the timestepper. **/
    // virtual void initialize(SimTK::State& advanced);

    void PrintInitialParams();

    void rebuildSimbodyTopologyFromOpenMMPositions(SimTK::State& someState);

    void setSphereRadius(float argSphereRadius);

    ///////////////////////////////////////////////////////
    // PROPOSE
    ///////////////////////////////////////////////////////

    /** Returns the 'how' argument of perturbPositions */
    PositionsPerturbMethod positionsPerturbMethod();

    VelocitiesPerturbMethod velocitiesPerturbMethod();

    /** Returns the 'how' argument of perturbVelocities */
    ForcesPerturbMethod forcesPerturbMethod();

    // Docking functions
    SimTK::Real getComComDistance(SimTK::State& someState,
                                  SimTK::MobilizedBodyIndex mbx1,
                                  SimTK::MobilizedBodyIndex mbx2);
    SimTK::Transform getRandomSphericalTransform(SimTK::Real radius);
    SimTK::Transform getRandomFM(SimTK::State& someState, SimTK::Real minDist, SimTK::Real maxDist);
    void teleport(SimTK::State& someState);

    // Perturb Q, QDot or QDotDot
    void perturb_Q_QDot_QDotDot(SimTK::State& someState);

    void printDrilling(SimTK::State& someState);

    auto sampleIteration(SimTK::State& state,
                         std::vector<SimTK::Compound::AtomTargetLocations>& proposedAtomTargetLocations,
                         bool shouldPrint) -> bool;

    [[nodiscard]] auto integrateWithOpenMM(SimTK::State& state) -> NUTSResult;

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
    void PrintAdaptiveData();

    int getMDStepsPerSample() const;

    void setMDStepsPerSample(int mdStepsPerSample);

    /** Calculate Mean Square Displacement based on stored R vectors **/
    SimTK::Real calculateMSD();

    /** Calculate RRdot based on stored R and Rdot vectors **/
    SimTK::Real calculateRRdot();

    /** Load the map of mobods to joint types **/
    // void loadMbx2mobility(SimTK::State& someState);

    /*
     * Test ground mostly for SOA
     */
    void testSOA(SimTK::State& someState);

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

    const SimTK::Real& getBendStretchStdevScaleFactor();
    void setBendStretchStdevScaleFactor(const SimTK::Real& s);

    /** Shift all the generalized coordinates and
     * return the scale factors of angles and bonds
     **/
    SimTK::Real setQToScaleBendStretch(SimTK::State& someState, std::vector<SimTK::Real>& scaleFactors);

    SimTK::Real setQToShiftBendStretchStdev(SimTK::State& someState, std::vector<SimTK::Real>& scaleFactors);

    /** Shift all the generalized coordinates and
     * return the scale factors of angles and bonds
     **/
    SimTK::Real setQToScaleBendStretchStdev(SimTK::State& someState, std::vector<SimTK::Real>& scaleFactors);

    void setQToScaleBendStretchStdev_Old(SimTK::State& someState, std::vector<SimTK::Real>& scaleFactors);

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

    void
    setSubZMatrixBATStats(std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATmeans,
                          std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATdiffs,
                          std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATstds,
                          std::map<SimTK::Compound::AtomIndex, std::vector<SimTK::Real>&> inBATstds_Alien);

    void PrintSubZMatrixBATAndRelated(SimTK::State& someState);

    SimTK::Real scaleSubZMatrixBATDeviations(
        SimTK::State& someState,
        SimTK::Real scalingFactor,
        bool BernoulliTrial = true,
        bool varianceBasedScalingFactor = true,
        std::vector<int> BATOrder = {1, 0, 2},        // bendstretch {1, 0, 2};  spherical {2, 1, 0}
        std::vector<SimTK::Real> BATSign = {1, -1, 1} // bendstretch {1, -1, 1}; spherical {1, 1, 1}
    );

    void updateSubZMatrixBAT(SimTK::State& someState,
                             std::vector<int> BATOrder = {1, 0, 2},
                             std::vector<SimTK::Real> BATSign = {1, -1, 1});

    SimTK::Real calcBATJacobianDetLog(SimTK::State& someState,
                                      SimTK::BondMobility::Mobility bondMobility,
                                      std::vector<int> BATOrder = {1, 0, 2},
                                      std::vector<SimTK::Real> BATSign = {1, -1, 1});

    // Calculate sub determinant of MBAT
    SimTK::Real calcSubMBATDetLog(SimTK::State& someState);

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

    void set_dBMps(std::vector<SimTK::Real>& QArg) {
        prev_dBMps = &QArg;
    }
    const std::vector<SimTK::Real>& get_dBMps() const {
        return *prev_dBMps;
    }
    void set_dPFrs(std::vector<SimTK::Real>& QArg) {
        prev_dPFrs = &QArg;
    }
    const std::vector<SimTK::Real>& get_dPFrs() const {
        return *prev_dPFrs;
    }

    void set_BMps_means(std::vector<SimTK::Real>& QArg) {
        prev_BMps_means = &QArg;
    }
    const std::vector<SimTK::Real>& get_BMps_means() const {
        return *prev_BMps_means;
    }
    void set_PFrs_means(std::vector<SimTK::Real>& QArg) {
        prev_PFrs_means = &QArg;
    }
    const std::vector<SimTK::Real>& get_PFrs_means() const {
        return *prev_PFrs_means;
    }

    void setPreviousQs(std::vector<SimTK::Real>& QArg) {
        previousQs = &QArg;
    }

    void setQmeans(std::vector<SimTK::Real>& QArg) {
        Qmeans = &QArg;
    }
    void setQdiffs(std::vector<SimTK::Real>& QArg) {
        Qdiffs = &QArg;
    }
    void setQvars(std::vector<SimTK::Real>& QArg) {
        Qvars = &QArg;
    }

    // Doesn't take masses into account
    SimTK::Real calcMobodsMBAT(const SimTK::State& someState);
    SimTK::Real calcMobodsBATJacobianDetLog_NEW(const SimTK::State& someState);
    SimTK::Real studyBATScale(const SimTK::State& someState);

#pragma region REBAS_TEST
    void setReplica(int thisReplica) {
        this->replicaIx = thisReplica;
    }
    void setThermodynamicState(int thisThermoStateIx) {
        this->thermoStateIx = thisThermoStateIx;
    }
#pragma endregion REBAS_TEST

    void setCartesianRandomSteps(int minSteps, int maxSteps) {
        cartesianRandomSteps = std::uniform_int_distribution<>(minSteps, maxSteps);
    }

    protected:
#pragma region REBAS_TEST
    int replicaIx;
    int thermoStateIx;
#pragma endregion REBAS_TEST

    // Buffers to hold Q statistics
    std::vector<SimTK::Real>* prev_BMps_means = nullptr;
    std::vector<SimTK::Real>* prev_PFrs_means = nullptr;
    std::vector<SimTK::Real>* prev_dBMps = nullptr;
    std::vector<SimTK::Real>* prev_dPFrs = nullptr;

    std::vector<SimTK::Real>* previousQs = nullptr;
    std::vector<SimTK::Real>* Qmeans = nullptr;
    std::vector<SimTK::Real>* Qdiffs = nullptr;
    std::vector<SimTK::Real>* Qvars = nullptr;

    // BEGIN MCSampler
    SimTK::Real detmbat_set = 0.0, detmbat_o = 0.0, detmbat_n = 0.0;
    SimTK::Real residualEmbeddedPotential = 0.0; // inside rigid bodies if weren't rigid

    bool useFixman = false;

    SimTK::Real MDStepsPerSampleStd = 0.5;
    SimTK::Real timestep = SimTK::NaN;
    int MDStepsPerSample = SimTK::NaN;

    int numSamples = 0;
    int numAcceptedSamples = 0;

    int QsBufferSize = 300;
    // std::list<SimTK::Vector> QsBuffer;
    std::deque<SimTK::Real> QsBuffer;

    // Non-equilibrium options
    int DistortOpt = 0;
    int FlowOpt = 0;
    int WorkOpt = 0;

    // END MCSampler

    std::vector<SimTK::Real> R;
    std::vector<SimTK::Real> Rdot;

    std::vector<SimTK::Real> dR;
    std::vector<SimTK::Real> dRdot;

    // Integration
    IntegratorType integratorType = IntegratorType::Empty;

    // For RANDOM_WALK Docking Simulations
    SimTK::Vec3 geometricCenter;
    SimTK::Real sphereRadius = SimTK::NaN;

    // Sampling
    int sampleGenerator = 0;

    std::vector<SimTK::Real> UScaleFactors;
    SimTK::Real UScaleFactorsNorm = 0.0;
    std::vector<SimTK::Real> InvUScaleFactors;
    SimTK::Real InvUScaleFactorsNorm = 0.0;

    SimTK::Vector NormedUScaleFactors;
    SimTK::Vector DOFUScaleFactors;

    std::vector<std::vector<SimTK::Real>> NMARotation;

    SimTK::Real ke_prop_nma6;
    SimTK::Real ke_n_nma6;

    // Transform Jacobian
    SimTK::Real bendStretchJacobianDetLog = 0.0;

    SimTK::Real boostT = SimTK::NaN, boostRT = SimTK::NaN,
                sqrtBoostRT = SimTK::NaN, // vel init
        boostBeta = SimTK::NaN;

    SimTK::Real boostKEFactor;
    SimTK::Real unboostKEFactor;

    SimTK::Real boostUFactor;
    SimTK::Real unboostUFactor;
    int boostMDSteps;

    SimTK::Real NMAAltSign = 1.0;
    SimTK::Real QScaleFactor = 1.0;

    EnergySnapshot previousEnergy, currentEnergy;

    // DELETE
    SimTK::Real debug_rand_no = 0.0;

    std::uniform_int_distribution<> cartesianRandomSteps{250, 2500};

    bool useNUTS = false;
    SimTK::Vector sqrtMInvV;
    SimTK::Vector dummyJointForces;
    SimTK::Vector dummyAccelerations;

    SimTK::Vector nutsOpenMMDeltaQ;
    SimTK::Vector nutsSimbodyDeltaQ;
    std::uniform_real_distribution<SimTK::Real> uniformReal01{0.0, 1.0};
    std::exponential_distribution<SimTK::Real> expDist{1.0};

    std::vector<SimTK::Real> totalEnergiesBuffer;

    NUTSWorkspaceSimbody nutsWorkspaceSimbody;
    int nutsMaxDepthSimbody = 16;
};
