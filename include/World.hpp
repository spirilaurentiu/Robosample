#pragma once

/* -------------------------------------------------------------------------- *
 *		                       Robosampling                           *
 * -------------------------------------------------------------------------- *
 * This is part of Robosample		                                      *
 */

#include <unistd.h>

#include <array>
#include <atomic>
#include <cmath>
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "Compound.h"
#include "CompoundSystem.h"
#include "Constraint.h"
#include "DuMMForceFieldSubsystem.h"
#include "FixmanTorque.hpp"
#include "TopologyElements.hpp"
#include "Transform.h"
#include "common.h"

#ifndef BaseSampler
#    define BaseSampler HMCSampler
#endif

#include "HMCSampler.hpp"
#include "Topology.hpp"

struct TopoAtom {
    int topoIx;
    SimTK::Compound::AtomIndex cAIx;
};

struct SpatialForceSample {
    int outboardPrmtopIndex;
    SimTK::Vec3 force;
    SimTK::Vec3 torque;
    SimTK::Real u;
    SimTK::Real uDot;
};

// ----------------------------------------------------------------------------
//  A rigid unit == a maximal set of atoms connected by non-ring Rigid tree
//  bonds. This reproduces CompoundSystem::buildUpRigidBody exactly: that
//  function recurses only over tree neighbours joined by Rigid bonds, and
//  ring-closing bonds are never tree edges, so they never merge bodies.
// ----------------------------------------------------------------------------
struct WorldRigidUnit {
    std::vector<int> atomCAIxs; // members (compound atom indices)
    int rootCAIx = -1;          // inboard atom of the unit
    int parentUnit = -1;        // -1 == attached to Ground
    int jointParentCAIx = -1;   // chemical parent of rootCAIx (in parent unit)
    SimTK::BondMobility::Mobility mobility = SimTK::BondMobility::Mobility::Rigid; // inboard joint mobility
    SimTK::MobilizedBodyIndex mbx; // filled when the body is created
};

// ---------------------------------------------------------------------------
//  Model-time cache. Built once (modelTopologies); reused every transfer.
// ---------------------------------------------------------------------------
struct FrameGraph {
    int totalAtoms = 0;
    std::vector<int> topoOffset; // size = numTopologies + 1; gAIx = topoOffset[t] + cAIx

    // Bucket G: full-geometry atoms (gp>=0 && rc>=0). Struct-of-arrays.
    std::vector<int> g_self, g_parent, g_gparent, g_refChild;

    // Bucket F: fallback atoms (root-child with no grandparent, OR leaf).
    std::vector<int> f_self, f_parent;

    // Bucket R: roots (one per molecule).
    std::vector<int> r_self;

    [[nodiscard]] auto numFull() const -> int {
        return (int)g_self.size();
    }
    [[nodiscard]] auto numFallback() const -> int {
        return (int)f_self.size();
    }
    [[nodiscard]] auto numRoot() const -> int {
        return (int)r_self.size();
    }
};


// Describes bonds involving a root atom
struct RootAtomBond {
    std::size_t topologyIndex = 0;

    SimTK::Compound::AtomIndex childCAIx;
    SimTK::Compound::AtomIndex parentCAIx;

    SimTK::MobilizedBodyIndex childMBIx;
    SimTK::MobilizedBodyIndex parentMBIx;
};

// Describes atom bonds and angles between two linked rigid bodies
struct RigidBodyAtomBond {
    std::size_t topologyIndex = 0;
    SimTK::BondMobility::Mobility mobility = SimTK::BondMobility::Default;

    std::size_t childAtomGlobalIndex = -1;
    std::size_t parentAtomGlobalIndex = -1;
    std::size_t grandParentAtomGlobalIndex = -1;

    SimTK::Compound::AtomIndex childCAIx;
    SimTK::Compound::AtomIndex parentCAIx;
    SimTK::Compound::AtomIndex grandParentCAIx;

    // Compound atom index of the root atom in the parent rigid body
    SimTK::Compound::AtomIndex parentMobodRootCAIx;

    SimTK::MobilizedBodyIndex childMBIx;
    SimTK::MobilizedBodyIndex parentMBIx;
    SimTK::MobilizedBodyIndex grandParentMBIx;
};

struct RigidBond {
    std::size_t topologyIndex = 0;
    SimTK::Compound::AtomIndex childCAIx;
    SimTK::Compound::AtomIndex parentCAIx;
    bool ringClosing = false;
};

struct RigidAngle {
    std::size_t topologyIndex = 0;
    SimTK::Compound::AtomIndex cAIx1, cAIx2, cAIx3;
    bool ringClosing = false;
};

struct RigidTorsion {
    std::size_t topologyIndex = 0;
    SimTK::Compound::AtomIndex cAIx1, cAIx2, cAIx3, cAIx4;
    bool ringClosing = false;
};

struct BondStretchKey {
    SimTK::DuMM::AtomClassIndex atomClassIndex1;
    SimTK::DuMM::AtomClassIndex atomClassIndex2;

    BondStretchKey(SimTK::DuMM::AtomClassIndex aCIx1, SimTK::DuMM::AtomClassIndex aCIx2) {
        atomClassIndex1 = std::min(aCIx1, aCIx2);
        atomClassIndex2 = std::max(aCIx1, aCIx2);
    }

    auto operator<(const BondStretchKey& other) const -> bool {
        if (atomClassIndex1 != other.atomClassIndex1) {
            return atomClassIndex1 < other.atomClassIndex1;
        }
        return atomClassIndex2 < other.atomClassIndex2;
    }
};

struct BondStretchValue {
    std::array<int, 2> globalAtomIndices;
    SimTK::Real stiffness;
    SimTK::Real length;

    auto operator==(const BondStretchValue& other) const -> bool {
        static constexpr SimTK::Real epsilon = 1e-9;
        return std::abs(stiffness - other.stiffness) < epsilon && std::abs(length - other.length) < epsilon;
    }

    auto operator!=(const BondStretchValue& other) const -> bool {
        return !(*this == other);
    }
};

struct BondBendKey {
    SimTK::DuMM::AtomClassIndex atomClassIndex1;
    SimTK::DuMM::AtomClassIndex atomClassIndex2; // center atom, stays fixed
    SimTK::DuMM::AtomClassIndex atomClassIndex3;

    BondBendKey(SimTK::DuMM::AtomClassIndex a1,
                SimTK::DuMM::AtomClassIndex a2,
                SimTK::DuMM::AtomClassIndex a3)
        : atomClassIndex2(a2) {
        atomClassIndex1 = std::min(a1, a3);
        atomClassIndex3 = std::max(a1, a3);
    }

    bool operator<(const BondBendKey& other) const {
        if (atomClassIndex1 != other.atomClassIndex1) {
            return atomClassIndex1 < other.atomClassIndex1;
        }
        if (atomClassIndex2 != other.atomClassIndex2) {
            return atomClassIndex2 < other.atomClassIndex2;
        }
        return atomClassIndex3 < other.atomClassIndex3;
    }
};

struct BondBendValue {
    std::array<int, 3> globalAtomIndices;
    SimTK::Real stiffness;
    SimTK::Real angleDeg;

    bool operator==(const BondBendValue& other) const {
        static constexpr SimTK::Real epsilon = 1e-9;
        return std::abs(stiffness - other.stiffness) < epsilon
               && std::abs(angleDeg - other.angleDeg) < epsilon;
    }

    bool operator!=(const BondBendValue& other) const {
        return !(*this == other);
    }
};

struct PeriodicTorsionKey {
    SimTK::DuMM::AtomClassIndex a1, a2, a3, a4;

    PeriodicTorsionKey(SimTK::DuMM::AtomClassIndex i,
                       SimTK::DuMM::AtomClassIndex j,
                       SimTK::DuMM::AtomClassIndex k,
                       SimTK::DuMM::AtomClassIndex l,
                       bool canonicalize) {
        // We don't canonicalize for improper torsions
        if (!canonicalize) {
            a1 = i;
            a2 = j;
            a3 = k;
            a4 = l;
            return;
        }

        // canonicalize (i,j,k,l) == (l,k,j,i)
        if (std::tie(i, j, k, l) <= std::tie(l, k, j, i)) {
            a1 = i;
            a2 = j;
            a3 = k;
            a4 = l;
        } else {
            a1 = l;
            a2 = k;
            a3 = j;
            a4 = i;
        }
    }

    bool operator<(const PeriodicTorsionKey& o) const {
        return std::tie(a1, a2, a3, a4) < std::tie(o.a1, o.a2, o.a3, o.a4);
    }
};

struct PeriodicTorsionValue {
    // store all 5 AMBER terms exactly as passed to Molmodel
    std::array<int, 4> globalAtomIndices;
    std::array<int, 5> periodicity;
    std::array<SimTK::Real, 5> amplitude;
    std::array<SimTK::Real, 5> phase;
    int numTerms;

    bool operator==(const PeriodicTorsionValue& o) const {
        if (numTerms != o.numTerms) {
            return false;
        }

        static constexpr SimTK::Real eps = 1e-9;
        for (int i = 0; i < 5; i++) {
            if (periodicity[i] != o.periodicity[i]) {
                return false;
            }
            if (std::abs(amplitude[i] - o.amplitude[i]) > eps) {
                return false;
            }
            if (std::abs(phase[i] - o.phase[i]) > eps) {
                return false;
            }
        }
        return true;
    }
    bool operator!=(const PeriodicTorsionValue& o) const {
        return !(*this == o);
    }
};

inline std::string getAtomDescription(const RoboAtom& atom) {
    return atom.identity.uniqueAtomName + " (atom class: " + atom.identity.atomClassName
           + ", atom class index " + std::to_string(atom.identity.atomClassIndex)
           + ", charged type: " + atom.identity.chargedAtomTypeName + ", charged type index "
           + std::to_string(atom.identity.chargedAtomTypeIndex) + ")";
}

class HarmonicImproperTorsionForce : public SimTK::DuMM::CustomBondTorsion {
    public:
    HarmonicImproperTorsionForce(SimTK::Real forceConstantInKJPerMol, SimTK::Real equilibriumAngleInRadians)
        : k(forceConstantInKJPerMol)
        , psi0(equilibriumAngleInRadians) {
    }

    SimTK::Real calcEnergy(SimTK::Real torsionInRadians) const override {
        const SimTK::Real dpsi = torsionInRadians - psi0;
        return k * dpsi * dpsi;
    }

    SimTK::Real calcTorque(SimTK::Real torsionInRadians) const override {
        const SimTK::Real dpsi = torsionInRadians - psi0;
        return -2.0 * k * dpsi;
    }

    private:
    SimTK::Real k;    // kψ in kJ/mol/rad^2
    SimTK::Real psi0; // equilibrium angle (rad)
};

//==============================================================================
//                   CLASS World
//==============================================================================
/**
 *  Contains a Symbody system and additional data that define a regimen
 **/
class World {
    public:
    void setAtomTargetLocationsToState(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets);
    void updateInboardAndOutboardFramesFromTopologies();

    explicit World(int worldIndex,
                   Span<Topology> topo,
                   bool testing,
                   const ZMatrix& _zMatrix,
                   bool wantSpatialForceHistory);

    [[nodiscard]] auto getAtomTargetLocationsCache() const
        -> const std::vector<SimTK::Compound::AtomTargetLocations>& {
        return atomTargetLocationsCache;
    }

    void setMobodLocks(const std::vector<std::vector<SimTK::MobilizedBodyIndex>>& mobodLocks) {
        this->mobodLocks = mobodLocks;
    }

    void generateDummParams(const std::vector<RoboAtom>& atoms,
                            const std::vector<RoboBond>& bonds,
                            const std::vector<RoboAngle>& angles,
                            const std::vector<RoboPeriodicTorsion>& properPeriodicTorsions,
                            const std::vector<RoboHarmonicImproperTorsion>& harmonicImproperTorsions);

    void modelTopologies(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets);

    SimTK::Real getRecommendedTimesteps();

    /** Add contact constraints to specific bodies **/
    void addRodConstraint(SimTK::State& someState);

    /** Add contact constraints to specific bodies **/
    const SimTK::State& addSpeedConstraint(int prmtopIndex);

    /** Calc station Jacobian */
    void calcStationJacobian(const SimTK::State& someState, SimTK::Matrix_<SimTK::Vec3>& JS) const;

    /** Get U scale factor for the mobilized body **/
    SimTK::Real getMobodUScaleFactor(SimTK::MobilizedBodyIndex&) const;
    //...............

    // Get geometric center of a subset of atoms
    // TEODOR
    SimTK::Vec3 getGeometricCenterOfSelection(const SimTK::State& state);

    /** Update Gmolmodel Atom Cartesian coordinates according to
    Molmodel Compound which in turn relizes Position and uses matter
     to calculate locations. **/
    void updateAtomListsFromSimbody(const SimTK::State& state);

    //.......................
    // --- Thermodynamics ---
    //.......................
    /** Get the World (macro) temperature **/
    SimTK::Real getTemperature();

    /** Set the World (macro) temperature **/
    void setTemperature(SimTK::Real);
    //...............

    /** Set the World (macro) temperature **/
    void setBoostTemperature(SimTK::Real);
    //...............

    // --- Simulation ---
    /** Get/Set seed for reproducibility. **/
    void setSeed(uint32_t argSeed);

    /** Use the Fixman torque as an additional force subsystem.
    Careful not have different temperatures for World and Fixman Torque. **/
    void addFixmanTorque();

    [[nodiscard]] auto getMultibodySystem() const -> const SimTK::CompoundSystem& {
        return *multibodySystem;
    }
    [[nodiscard]] auto updMultibodySystem() const -> SimTK::CompoundSystem& {
        return *multibodySystem;
    }

    [[nodiscard]] auto getMatterSubsystem() const -> const SimTK::SimbodyMatterSubsystem& {
        return *matter;
    }
    [[nodiscard]] auto updMatterSubsystem() const -> SimTK::SimbodyMatterSubsystem& {
        return *matter;
    }

    [[nodiscard]] auto getForces() const -> const SimTK::GeneralForceSubsystem& {
        return *forces;
    }
    [[nodiscard]] auto updForces() const -> SimTK::GeneralForceSubsystem& {
        return *forces;
    }

    const SimTK::DuMMForceFieldSubsystem& getForceField() const {
        return *forceField;
    }
    SimTK::DuMMForceFieldSubsystem& updForceField() {
        return *forceField;
    }

    const SimTK::VerletIntegrator& getIntegrator() const {
        return *integrator;
    }
    SimTK::VerletIntegrator& updIntegrator() {
        return *integrator;
    }

    const SimTK::TimeStepper& getTimeStepper() const {
        return *timeStepper;
    }
    SimTK::TimeStepper& updTimeStepper() {
        return *timeStepper;
    }

    /** Return true if the Fixman torque flag is set **/
    bool isUsingFixmanTorque() const;
    //...............

    // Calculate Fixman potential
    [[nodiscard]] auto getFixmanPotential() const -> SimTK::Real {
        return getSampler(0)->getCurrentEnergy().fixman;
    }

    /** Generate a number of samples **/
    auto generateSamples(int howMany,
                         std::stringstream& worldOutStream,
                         const std::string& header,
                         bool shouldPrint) -> bool;

    //...................
    // --- Statistics ---
    //...................
    /** How many samples did we have so far **/
    std::size_t getNofSamples() const;

    /** Sampler manipulation functions **/
    std::size_t getNofSamplers() const;

    /** Add a sampler to the World **/
    auto addSampler(SamplerName samplerName,
                    IntegratorType integratorType,
                    ThermostatName thermostatName,
                    bool useFixmanPotential,
                    bool useNUTS) -> bool;

    // TODO Use Sampler polymorphism
    /** Get a sampler based on its position in the samplers vector **/
    [[nodiscard]] auto getSampler(std::size_t which) const -> const BaseSampler* {
        return samplers[which].get();
    }

    /** Get a writable sampler based on its position in the samplers vector **/
    [[nodiscard]] auto updSampler(std::size_t which) const -> BaseSampler* {
        return samplers[which].get();
    }

    /** Get writble pointer to FixmanTorque implementation **/
    FixmanTorque* updFixmanTorque();

    /** Get pointer to FixmanTorque implementation **/
    FixmanTorque* getFixmanTorque() const;

    // Get the (potential) energy transfer
    // If any of the Q, U or tau is actively modifyied by the sampler
    // the Jacobian of that transformation will be included too
    SimTK::Real getWorkOrHeat();

    // Get the (potential) energy transfer in the form of work
    // If any of the Q, U or tau is actively modifyied by the sampler
    // the Jacobian of that transformation will be included too
    SimTK::Real getWork() const;

    // Set initial values of X_PF or X_BM
    void setTransformsMeansToIni();

    // Set initial values of X_PF or X_BM
    void setTransformsMeansToCurrent(SimTK::State& someState);

    // Set initial values of X_PF or X_BM
    void setTransformsMeans(const std::vector<SimTK::Real>& givenX_PF,
                            const std::vector<SimTK::Real>& givenX_BM);

    // Get X_PF and X_BM related values
    void getTransformsStatistics(SimTK::State& someState);

    /** Update transforms means given a previous value */
    void updateTransformsMeans(SimTK::State& someState);

    /** Get X_PF means */
    std::vector<SimTK::Real>& getX_PFMeans();

    /** Get X_BM means */
    std::vector<SimTK::Real>& getX_BMMeans();

    /**
     * Calculate bond length and angle deviations from their means
     */
    void calcBendStretchDeviations(SimTK::State& someState,
                                   std::vector<SimTK::Real>& X_PFdiffs,
                                   std::vector<SimTK::Real>& X_BMdiffs);

    // Print bond lengths and angle bends
    void traceBendStretch(SimTK::State& someState);

    // Print X_PF
    void PrintAcosX_PFs();

    // Print X_PF
    void PrintNormX_BMs();

    // Print X_PF means
    void PrintAcosX_PFMeans();

    // Print X_PF means
    void PrintNormX_BMMeans();

    //...............

    // REORIENT

    SimTK::Transform& getReorientTransformInAnotherBody(const SimTK::State& someState,
                                                        const SimTK::MobilizedBody& inBodyA,
                                                        const SimTK::MobilizedBody& ofBodyB,
                                                        const SimTK::Transform& reorientAB,
                                                        SimTK::Transform& X_FMprim);

    // RANDOM_WALK related functions; we don't need getter, since we only
    // use these values inside the scope of World.
    void setTopologyIXs(std::vector<int> topologyIXs);
    void setAmberAtomIXs(std::vector<std::vector<int>> AmberAtomIXs);

    const SimTK::Vector& getBMps();
    const SimTK::Vector& getPFrs();

    // Get Qs
    [[nodiscard]] auto getNQsFromAdvancedState() const -> int {
        return matter->getNQ(worldState);
    }

    [[nodiscard]] auto getNUsFromAdvancedState() const -> int {
        return matter->getNU(worldState);
    }

    [[nodiscard]] auto getAdvancedQs() const -> const SimTK::Vector& {
        return matter->getQ(worldState);
    }

    const SimTK::Vector& getAdvancedUs();

    void calcSimbodyBAT(std::vector<std::vector<int>>& ZMatrix,
                        std::vector<SimTK::Real>& BONDLengths,
                        std::vector<SimTK::Real>& ANGLEBends,
                        std::vector<SimTK::Real>& TORSIONAngles);

    void setFlexibilites(const std::vector<std::vector<SimTK::MobilizedBodyIndex>>& flexibilities_UNCHAINED) {
        this->flexibilities_UNCHAINED = flexibilities_UNCHAINED;
    }

    // This is non-copyable and non-movable, so we use a unique_ptr because World needs to be
    // copyable/movable.
    std::unique_ptr<SimTK::CompoundSystem> multibodySystem;

    // Subsystem->SimbodyMatterSubsystem
    // This is non-copyable and non-movable, so we use a unique_ptr because World needs to be
    // copyable/movable.
    std::unique_ptr<SimTK::SimbodyMatterSubsystem> matter;

    // Subsystem->ForceSubsystem->GeneralForceSubsystem
    // This is non-copyable and non-movable, so we use a unique_ptr because World needs to be
    // copyable/movable.
    std::unique_ptr<SimTK::GeneralForceSubsystem> forces;

    // Subsystem->ForceSubsystem->DuMMForceFieldSubsystem
    // This is non-copyable and non-movable, so we use a unique_ptr because World needs to be
    // copyable/movable.
    std::unique_ptr<SimTK::DuMMForceFieldSubsystem> forceField;

    // --- Simulation ---
    std::unique_ptr<SimTK::VerletIntegrator> integrator;
    std::unique_ptr<SimTK::TimeStepper> timeStepper;
    // TODO they belong to Sampler, not World

    std::vector<std::unique_ptr<BaseSampler>> samplers;

    std::vector<std::vector<SimTK::MobilizedBodyIndex>> flexibilities_UNCHAINED;

    SimTK::Vector BMps;
    SimTK::Vector PFrs;

    // std::vector<std::unique_ptr<SimTK::ConformationalController>> controller;
    // std::vector<std::unique_ptr<SimTK::Force::Custom>> controlForce;

    /** Nof molecules **/
    std::size_t numMolecules = 0;
    std::size_t numAtoms = 0;

    /** Molecules (topologies<-Compounds) objects **/
    Span<Topology> topologies;
    std::vector<std::string> roots;
    std::vector<std::string> rootMobilitiesStr;

    /** Joint types **/
    // std::map< SimTK::MobilizedBodyIndex, SimTK::BondMobility::Mobility> mbx2mobility;

    //
    std::vector<std::vector<int>> zMatrixTable;
    std::vector<std::vector<SimTK::Real>> zMatrixBAT;

    // --- Thermodynamics ---
    SimTK::Real temperature = SimTK::NaN;

    // --- Statistics ---
    std::vector<SimTK::Real> acosX_PF00;
    std::vector<SimTK::Real> normX_BMp;
    std::vector<SimTK::Real> acosX_PF00_means;
    std::vector<SimTK::Real> normX_BMp_means;

    // --- Mixing data ---
    int ownWorldIndex;

    /** Get writble pointer to Fixman Torque and other forces**/
    // std::unique_ptr<FixmanTorque> FixmanTorqueImpl;
    // std::unique_ptr<SimTK::Force::Custom> FixmanTorqueForce;

    std::unique_ptr<SimTK::Force::Custom> FixmanTorqueForce;
    // std::unique_ptr<FixmanTorque> FixmanTorqueImpl;
    std::unique_ptr<SimTK::Force::Custom> FixmanTorqueExtForce;
    // std::unique_ptr<FixmanTorqueExt> FixmanTorqueExtImpl;

    FixmanTorque* FixmanTorqueImpl = nullptr;
    // SimTK::Force::Custom* FixmanTorqueForce = nullptr;
    FixmanTorqueExt* FixmanTorqueExtImpl = nullptr;
    // SimTK::Force::Custom* FixmanTorqueExtForce = nullptr;

    // Task Space
    SimTK::Array_<SimTK::MobilizedBodyIndex> onBodyB;
    SimTK::Array_<SimTK::Vec3> taskStationPInGuest;
    SimTK::Array_<SimTK::Vec3> taskStationPInHost;
    SimTK::Array_<SimTK::Vec3> taskDeltaStationP;

    // Constraints
    std::vector<std::pair<SimTK::MobilizedBodyIndex, SimTK::MobilizedBodyIndex>> rodBodies;
    SimTK::Array_<SimTK::Vec3> conStationPInGuest;
    SimTK::Array_<SimTK::Vec3> conStationPInHost;
    SimTK::Array_<SimTK::Vec3> conDeltaStationP;
    SimTK::Array_<SimTK::Constraint::Rod> rodConstraints;

    // X axis to Z axis switch
    const SimTK::Transform X_to_Z = SimTK::Rotation(-90 * SimTK::Deg2Rad, SimTK::YAxis);
    const SimTK::Transform Z_to_X = ~X_to_Z;

    // Y axis to Z axis switch
    const SimTK::Transform Y_to_Z = SimTK::Transform(SimTK::Rotation(-90 * SimTK::Deg2Rad, SimTK::XAxis));
    const SimTK::Transform Z_to_Y = ~Y_to_Z;

    // X axis to X axis switch
    const SimTK::Transform Y_to_X = SimTK::Rotation(-90 * SimTK::Deg2Rad, SimTK::ZAxis);
    const SimTK::Transform X_to_Y = ~Y_to_X;

    //
    void setSamplesPerRound(int samples);
    [[nodiscard]] auto getSamplesPerRound() const -> int;

    void setDistortOption(int distort);
    [[nodiscard]] auto getDistortOption() const -> int;

    [[nodiscard]] auto getOwnIndex() const -> int {
        return ownWorldIndex;
    }

    [[nodiscard]] auto getRootAtomBonds() const -> const std::vector<RootAtomBond>& {
        return rootAtomBonds;
    }

    [[nodiscard]] auto getRigidBodyAtomBonds() const -> const std::vector<RigidBodyAtomBond>& {
        return rigidBodyAtomBonds;
    }

    [[nodiscard]] auto getMobodRootAtomIndex(SimTK::MobilizedBodyIndex mbIndex) const -> const TopoAtom& {
        return mbxRootCAIx[mbIndex];
    }

    void calcSpatialForces();
    void writeSpatialForces(const std::string& filename) const;
    [[nodiscard]] auto getWantSpatialForceHistory() const -> bool {
        return wantSpatialForceHistory;
    }

    // BAT --------------------------------------------------------------------

    /**
     * @brief Drill
     * @param
     */
    void printDrilling();

    // std::vector<std::vector<int>>& zMatrixTable;
    // std::vector<std::vector<SimTK::Real>>& zMatrixBAT;

    // setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value);
    // void calcZMatrixBAT(SimTK::State& someState);

    void checkCoordinateTransfer(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets) const;

    [[nodiscard]] auto getMbx(int topoIx, int cAIx) const -> SimTK::MobilizedBodyIndex {
        return topoAtomToMbx[topoIx][cAIx];
    }

    [[nodiscard]] auto getDAIx(int topoIx, int cAIx) const -> SimTK::DuMM::AtomIndex {
        return topoAtomToDAIX[topoIx][cAIx];
    }

    private:
    bool testing = false;

    std::vector<SimTK::Compound::AtomTargetLocations> atomTargetLocationsCache;
    std::vector<SimTK::Compound::AtomTargetLocations> atomTargetLocationsCacheOld;
    std::vector<std::pair<bool, SimTK::Real>> acceptanceRMSD;

    std::vector<RootAtomBond> rootAtomBonds;
    std::vector<RigidBodyAtomBond> rigidBodyAtomBonds;

    std::set<SimTK::MobilizedBodyIndex> interestingMobodIndices;
    std::map<SimTK::MobilizedBodyIndex, int> mbx2PrmtopInboardIndex;
    std::map<int, int> prmtopInboardIndex2PrmtopOutboardIndex;
    std::unordered_map<int, std::vector<SpatialForceSample>> spatialForceHistory;
    bool wantSpatialForceHistory = false;

    // Maps a generalized velocity scale factor for every mobod
    std::map<SimTK::MobilizedBodyIndex, SimTK::Real> mbx2uScale;

    // Binding Site Data: Topologies, AtomIx
    std::vector<int> topologyIXs;
    std::vector<std::vector<int>> amberAtomIXs;

    // Track the Stage of the system
    SimTK::Stage currStage;

    bool useFixmanTorque = false;
    int samplesPerRound = 0;

    Random32 randomEngine;

    std::vector<std::vector<SimTK::MobilizedBodyIndex>> mobodLocks;

    // Default return value for non-existing topology atom, pair
    std::pair<int, SimTK::Compound::AtomIndex> errorTopoAtomPair{
        -1,
        SimTK::Compound::AtomIndex(SimTK::InvalidIndex)};

    std::reference_wrapper<const ZMatrix> zMatrix;

    SimTK::State worldState;

    /**
     * Maps a specific atom within a topology to its assigned global Mobilized Body.
     *
     * Lookup: [topology_index][compound_atom_index] -> mobilized_body_index
     *
     * Use this to determine which rigid body an atom belongs to. Mobilized Body indices (mbx) are contiguous
     * [0, N) across the entire system. They are NOT local to the topology. Example: topo[0] might contain mbx
     * 0-10, and topo[1] starts at mbx 11.
     */
    std::vector<std::vector<SimTK::MobilizedBodyIndex>> topoAtomToMbx;

    /**
     * Maps a specific atom within a topology to its assigned global DuMM Atom Index.
     *
     * Lookup: [topology_index][compound_atom_index] -> dumm_atom_index
     *
     * Use this to determine which DuMM Atom an atom belongs to, which in turn determines its force field
     * parameters. DuMM Atom indices (aIx) are contiguous [0, N) across the entire system. They are NOT local
     * to the topology. Example: topo[0] might contain aIx 0-10, and topo[1] starts at aIx 11.
     */
    std::vector<std::vector<SimTK::DuMM::AtomIndex>> topoAtomToDAIX;

    std::vector<std::vector<bool>> topoAtomIsRigidBodyRoot;

    /**
     * Identifies the "Root" atom for a given Mobilized Body.
     *
     * Lookup: [mobilized_body_index] -> compound_atom_index
     *
     * The root atom defines the origin of the mobilized body's local frame.
     * Since cAIx is part of the static Topology definition, this value is constant across all worlds, even
     * though the mbx itself is world-specific.
     */
    std::vector<TopoAtom> mbxRootCAIx;


    auto decomposeRigidUnits(const Topology& topology) const -> std::vector<WorldRigidUnit>;
    void buildFrameGraph(const Span<Topology>& topologies);
    auto modelOneCompound(int topoIx, SimTK::RootMobility rootMobility) -> std::vector<WorldRigidUnit>;

    FrameGraph frameGraph;                   // built once in modelTopologies
    std::vector<SimTK::Vec3> targetFlat;     // [gAIx] reused each transfer
    std::vector<SimTK::Transform> frameFlat; // [gAIx] reused each transfer
    std::vector<SimTK::Transform> xpcBcFlat; // [gAIx] B == X_parentBC_childBC

    // flat -> (topo,cAIx) accessor; inline, no bounds work in release.
    const SimTK::Transform& F(int topoIx, SimTK::Compound::AtomIndex cAIx) const {
        return frameFlat[frameGraph.topoOffset[topoIx] + int(cAIx)];
    }
    const SimTK::Transform& Xpc(int topoIx, SimTK::Compound::AtomIndex cAIx) const {
        return xpcBcFlat[frameGraph.topoOffset[topoIx] + int(cAIx)];
    }


    std::vector<std::vector<SimTK::Transform>> atomFrameCache;
    std::vector<std::vector<SimTK::Transform>> topoAtomBodyFrame;
    // std::vector<std::vector<SimTK::Vec3>> atomStations;
    std::vector<SimTK::Real> clustersMass;
    std::vector<SimTK::Vec3> clustersCOM;
    std::vector<SimTK::Inertia> clustersInertia;
};
