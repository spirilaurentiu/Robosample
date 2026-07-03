#pragma once

#include <cstddef>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "OpenMM.h"
#include "TopologyElements.hpp"

// Thin Robosample-owned wrapper around a single OpenMM System/Context (a
// process-wide singleton). One OpenMM system for the whole run.
class OpenMMContext {
    public:
    static auto get() -> OpenMMContext& {
        static OpenMMContext omm;
        return omm;
    }

    auto initialize(const SystemTopology& systemTopology) -> bool;

    void shutdown() {
        integrator.reset();
        context.reset();
        system.reset();
        forceGroupLabels.clear();
        hasVirtualSites = false;
        initialized = false;
    }

    void setVelocitiesToTemperature(double temperature, int seed) {
        ensureInitialized();
        context->setVelocitiesToTemperature(temperature, seed);
    }

    // ---- NCMC alchemy (per-molecule intermolecular decoupling) -------------
    // enableAlchemy stores the decoupled atom-index SET (Region A, docs/specs/
    // ncmc-explicit-solvent/30-region-and-protocol-policy.md Sec.2, DECIDED:
    // arbitrary set, not only contiguous) and a flag; the correction force is
    // built in initialize() (which has the SystemTopology). PME/periodic is
    // rejected there. setAlchemicalLambda drives the global "lambda_inter".
    void enableAlchemy(const std::vector<int>& atomIndices);
    // Convenience: contiguous [atomBegin,atomEnd) Region A.
    void enableAlchemy(int atomBegin, int atomEnd);
    void setAlchemicalLambda(double lambdaInter);

    [[nodiscard]] auto getPotentialEnergy() const -> double {
        ensureInitialized();
        return potentialEnergy;
    }
    [[nodiscard]] auto getKineticEnergy() const -> double {
        ensureInitialized();
        return kineticEnergy;
    }

    // OpenMM's own degrees-of-freedom count for the whole system, matching the
    // convention OpenMM's StateDataReporter uses: 3 per particle with nonzero
    // mass (massless virtual sites carry no independent dof), minus the number
    // of SHAKE/SETTLE distance constraints, minus 3 if a CMMotionRemover force
    // (center-of-mass motion removal) is present. Used as the Cartesian world's
    // n_dof, since model_.nu is NOT the physical dof there (Cartesian worlds
    // collapse every atom into a single nu==1 body).
    [[nodiscard]] auto getNumDegreesOfFreedom() const -> int {
        ensureInitialized();
        int dof = 0;
        const int numParticles = system->getNumParticles();
        for (int i = 0; i < numParticles; ++i) {
            if (system->getParticleMass(i) > 0.0) {
                dof += 3;
            }
        }
        dof -= system->getNumConstraints();
        const int numForces = system->getNumForces();
        for (int i = 0; i < numForces; ++i) {
            if (dynamic_cast<const OpenMM::CMMotionRemover*>(&system->getForce(i)) != nullptr) {
                dof -= 3;
                break;
            }
        }
        return dof;
    }

    struct ForceGroupEnergy {
        int group = 0;
        std::string name;
        double energy = 0.0; // kJ/mol
    };

    void setSeparateForceGroups(bool enabled) {
        separateForceGroups = enabled;
    }
    [[nodiscard]] auto getSeparateForceGroups() const -> bool {
        return separateForceGroups;
    }

    // Multiple-timestep (r-RESPA) control. When enabled, initialize() builds an
    // MTSIntegrator instead of the single-rate VerletIntegrator: slow forces
    // (Nonbonded, GBSA, NBFIX) go to force group 0 and are evaluated once per
    // outer step; fast bonded forces go to group 1 and are evaluated
    // `innerSubsteps` times. Only affects the Cartesian world's on-device MD
    // (integrateTrajectory); the torsional worlds always get the full force sum.
    // MUST be called before initialize(). innerSubsteps < 2 disables MTS.
    void setMTS(bool enabled, int innerSubsteps) {
        useMTS = enabled && (innerSubsteps >= 2);
        mtsInnerSubsteps = innerSubsteps;
    }
    [[nodiscard]] auto getUseMTS() const -> bool {
        return useMTS;
    }

    void setPositions(const std::vector<OpenMM::Vec3>& positions) {
        ensureInitialized();
        context->setPositions(positions);
        if (hasVirtualSites) {
            context->computeVirtualSites();
        }
    }
    void getPositions(std::vector<OpenMM::Vec3>& out) const {
        ensureInitialized();
        out = context->getState(OpenMM::State::Positions, enforcePeriodicBox).getPositions();
    }

    // Whether OpenMM wraps coordinates into the primary box when we PULL state
    // (positions/forces/energy getState calls). For Robosample this MUST stay
    // false under explicit solvent: the robot engine consumes per-atom positions
    // and rebuilds each molecule's internal frames from contiguous geometry, so a
    // wrapped molecule that straddles a box face would yield a ~box-length "bond"
    // and corrupt the frame build. Energies/forces are unaffected by this flag
    // (OpenMM always applies the minimum image internally); wrapping is purely a
    // representation of the returned positions. Any periodic-image bookkeeping for
    // visualization is done as a WHOLE-MOLECULE rigid translation elsewhere, never
    // here. Default false.
    void setEnforcePeriodicBox(bool v) {
        enforcePeriodicBox = v;
    }
    [[nodiscard]] auto getEnforcePeriodicBox() const -> bool {
        return enforcePeriodicBox;
    }

    // True for the methods that require a periodic box (and thus exclude GBSA).
    [[nodiscard]] static auto isPeriodic(NonbondedMethod m) -> bool {
        return m == NonbondedMethod::CutoffPeriodic || m == NonbondedMethod::Ewald
               || m == NonbondedMethod::PME;
    }

    auto computePotentialEnergy(const std::vector<OpenMM::Vec3>& positions) -> double;
    auto computePotentialEnergyByGroup(const std::vector<OpenMM::Vec3>& positions)
        -> std::pair<double, std::vector<ForceGroupEnergy>>;
    void evaluateForcesFromPositionsCache(const std::vector<OpenMM::Vec3>& positions,
                                          std::vector<OpenMM::Vec3>& outForces) const;
    auto integrateTrajectory(int steps, double timeStepInPicoseconds) -> bool;
    auto computePeriodicBoxVectors_Context(double a_length,
                                           double b_length,
                                           double c_length,
                                           double alpha,
                                           double beta,
                                           double gamma)
        -> std::tuple<OpenMM::Vec3, OpenMM::Vec3, OpenMM::Vec3>;

    private:
    OpenMMContext() = default;

    void ensureInitialized() const {
        if (!initialized) {
            throw std::runtime_error("OpenMM subsystem not initialized. Call OpenMMContext::initialize() "
                                     "before using any other functions.");
        }
    }

    [[nodiscard]] auto createNonbondedForce(const SystemTopology& systemTopology) -> OpenMM::NonbondedForce*;
    [[nodiscard]] auto createGBSAOBCForce(const SystemTopology& systemTopology) -> OpenMM::GBSAOBCForce*;
    [[nodiscard]] auto createCustomNonbondedForce(const SystemTopology& systemTopology)
        -> OpenMM::CustomNonbondedForce*;
    // Alchemy correction force: total [begin,end) x rest pair energy becomes
    // lambda_inter * standard (LJ + Coulomb), via a (lambda_inter-1)*standard term
    // over an interaction group. Every excluded pair is intramolecular (never an
    // A x rest pair), so the exclusions below are energy-neutral; they exist only
    // to satisfy the CPU platform's shared-neighbor-list rule.
    [[nodiscard]] auto createAlchemyCorrectionForce(const SystemTopology& systemTopology)
        -> OpenMM::CustomNonbondedForce*;
    // Mirror the main NonbondedForce's exception pairs (all 1-2/1-3 exclusions and
    // 1-4 scaled pairs) as CustomNonbondedForce exclusions. The CPU platform shares
    // ONE neighbor list across every exclusion-using force and rejects the Context
    // ("All Forces must have identical exclusions") unless the lists match exactly;
    // CUDA/OpenCL route interaction-group custom forces around the shared list, so
    // this is a CPU-correctness requirement and a no-op on the energy elsewhere.
    static void addStandardExclusions(OpenMM::CustomNonbondedForce* force,
                                      const SystemTopology& systemTopology);
    [[nodiscard]] auto createHarmonicBondForce(const SystemTopology& systemTopology)
        -> OpenMM::HarmonicBondForce*;
    [[nodiscard]] auto createHarmonicAngleForce(const SystemTopology& systemTopology)
        -> OpenMM::HarmonicAngleForce*;
    [[nodiscard]] auto createPeriodicTorsionForce(const SystemTopology& systemTopology)
        -> OpenMM::PeriodicTorsionForce*;
    [[nodiscard]] auto createImproperHarmonicTorsionForce(const SystemTopology& systemTopology)
        -> OpenMM::CustomTorsionForce*;
    [[nodiscard]] auto createCMAPTorsionForce(const SystemTopology& systemTopology)
        -> OpenMM::CMAPTorsionForce*;
    [[nodiscard]] auto createUreyBradleyForce(const SystemTopology& systemTopology)
        -> OpenMM::HarmonicBondForce*;

    // PME/explicit-solvent NCMC. Reciprocal space cannot be localized to an
    // A x rest pair list, so electrostatics are scaled on the MAIN NonbondedForce
    // via a lambda_inter charge offset (PME-exact), and A's LJ is rebuilt as
    // soft-core A x rest + hard intra-A custom forces. Mutates `main`; returns the
    // two custom forces to add. Driven by the SAME lambda_inter global parameter,
    // so setAlchemicalLambda / ncmcMove need no changes.
    [[nodiscard]] auto createAlchemyDecouplingForces(const SystemTopology& systemTopology,
                                                     OpenMM::NonbondedForce* main)
        -> std::pair<OpenMM::CustomNonbondedForce*, OpenMM::CustomNonbondedForce*>;

    std::unique_ptr<OpenMM::Context> context;
    std::unique_ptr<OpenMM::System> system;
    std::unique_ptr<OpenMM::Integrator> integrator;

    std::size_t numAtoms = 0;
    double potentialEnergy = 0;
    double kineticEnergy = 0;

    bool separateForceGroups = false;
    std::vector<std::pair<int, std::string>> forceGroupLabels;

    // NCMC alchemy state. alchemyAtoms is Region A: an ascending, deduplicated
    // atom-index set (may be non-contiguous; enableAlchemy(int,int) builds the
    // contiguous case). alchemyAtomSet mirrors it as a std::set<int> for O(log n)
    // membership tests in the force builders.
    bool alchemyEnabled = false;
    std::vector<int> alchemyAtoms;
    std::set<int> alchemyAtomSet;
    OpenMM::CustomNonbondedForce* alchemyForce = nullptr; // owned by `system`

    // r-RESPA multiple-timestep state.
    bool useMTS = false;
    int mtsInnerSubsteps = 4;
    static constexpr int kMtsSlowGroup = 0; // Nonbonded, GBSA, NBFIX
    static constexpr int kMtsFastGroup = 1; // bonds, angles, torsions, CMAP, UB

    bool enforcePeriodicBox = false;
    bool hasVirtualSites = false;
    bool initialized = false;
};