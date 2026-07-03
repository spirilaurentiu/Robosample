#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "DCDWriter.hpp"
#include "OpenMMContext.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"

// Robosample top-level context. Constructed from Python as Context(base_name,
// seed); the Python subclass (context.py) adds the dihedral classifier and
// load_amber, which fills `systemTopology` in place.
//
// Owns the worlds (the Gibbs flexibility regimens), the replica temperature
// ladder, and the REX driver. World order defines the Gibbs sweep order, and
// per-atom Ground coordinates (nm) are the only currency passed between worlds.
class Context {
    public:
    Context(std::string baseName, std::uint32_t seed);

    // Bound to Python as `system_topology`. Filled in place by load_amber.
    SystemTopology systemTopology;

    // ---- modeling --------------------------------------------------------
    // Root mobility is a PER-WORLD property. systemTopology.rootMobilities is
    // only the build-time DEFAULT each world seeds from (set by load_amber);
    // there is deliberately no Context-level mutator. To change an individual
    // molecule's root attachment, call World::setRootMobility(...) on the world
    // returned by add*World (see World.hpp), which rebuilds that world's model
    // in isolation -- the same mechanism addDockingWorld already uses.

    // Add a Cartesian world (pure OpenMM MD on device) or a robotic/torsional
    // world (internal-coordinate HMC). Both return a reference to the new world
    // so Python can chain .add_sampler(...). The model is built immediately
    // from the current systemTopology + root mobilities.
    //
    // wantReactionReporter (docs/specs/reaction-force-monitoring.md Sec.2): opt
    // -in, off by default. Flags the new world the per-body applied-force
    // reporter (World::setReactionReporter). addCartesianWorld always THROWS
    // when true (Sec.3 integrator guard -- a Cartesian world's articulated
    // body indexing is not meaningful). addRoboticWorld derives the
    // interesting-body set from the just-built model: every non-Weld (flexed)
    // body plus its parent, excluding only Ground itself -- a Free-rooted
    // body (e.g. a receptor's root, directly attached to Ground) IS included
    // (Sec.2.1/Sec.4).
    World& addCartesianWorld(bool wantReactionReporter = false);
    World& addRoboticWorld(const Selection& sel, bool wantReactionReporter = false);

    // Add a DOCKING world. `ligandMoleculeIndices` lists which molecules are
    // ligands; their roots become Free (6 external DOF) and ALL their bonds stay
    // Rigid (rigid-body docking). Every other molecule is Welded to Ground and
    // rigid. Root mobility is thus a WORLD property here -- it is built from this
    // argument, NOT from systemTopology.rootMobilities, so the same system can
    // host a docking world and ordinary worlds at once. The binding-site centre
    // is the centroid of all non-ligand (receptor) atoms; pass a sphere radius
    // via add_sampler(...). Returns the world for chaining .add_sampler(...).
    World& addDockingWorld(const std::vector<int>& ligandMoleculeIndices);


    // Per-bond mobility selection. `bonds` empty/None => every eligible
    // (non-ring, non-terminal) bond gets `mobility`; otherwise only the listed
    // (i,j) bonds do. Ring-closing bonds always remain Rigid.
    Selection buildFlexibilities(const std::optional<std::vector<std::pair<int, int>>>& bonds,
                                 JointType mobility,
                                 bool flag);

    // ---- run -------------------------------------------------------------
    // Build the OpenMM system, set the replica temperature ladder, and seed
    // every replica with the reference coordinates. Empty list => single 300 K.
    void initialize(const std::vector<double>& temperatures);

    // Replica-exchange driver: `equil` then `prod` rounds; each round runs a
    // Gibbs sweep over all worlds for every replica, attempts adjacent-replica
    // swaps, and writes outputs every `writeFreq` production rounds.
    void runREX(int equilRounds, int prodRounds, int writeFreq, bool verbose);

    // Enable multiple-timestep (r-RESPA) integration for the Cartesian world's
    // on-device OpenMM MD: slow forces (Nonbonded, GBSA, ...) are evaluated once
    // per outer step, fast bonded forces `innerSubsteps` times. No effect on the
    // torsional worlds (their forces are always the full sum). MUST be called
    // before initialize(). Forwards to OpenMMContext.
    void setMTS(bool enabled, int innerSubsteps);


    // ---- energy ingestion / validation (unchanged) ----------------------
    auto initializeOpenMM() -> bool;
    [[nodiscard]] auto calcOpenMMPotentialEnergy() -> double;
    [[nodiscard]] auto computePotentialEnergyByGroup()
        -> std::pair<double, std::vector<OpenMMContext::ForceGroupEnergy>>;

    [[nodiscard]] auto getBaseName() const -> const std::string& {
        return baseName;
    }
    [[nodiscard]] auto getSeed() const -> std::uint32_t {
        return seed;
    }
    [[nodiscard]] int numWorlds() const {
        return static_cast<int>(worlds_.size());
    }

    private:
    double openmmPotential(const std::vector<robo::Vec3>& coords);
    void writeOutputs(int replica, int round, bool verbose);
    // Append a reporter world's captured per-body force rows (docs/specs/
    // reaction-force-monitoring.md Sec.4) to that replica's per-replica CSV,
    // tagged with the DCD frame index they pair with. No-op if rows is
    // empty. Always the fixed 10-column schema (Sec.1.2/4).
    void writeReactionRows(int replica, int frame, const std::vector<ReactionSample>& rows);

    // Refuse to start (or warn, if ROBO_ALLOW_BAD_START is set) when the input
    // geometry is non-finite or sterically clashing -- the docking world cannot
    // repair a frozen-receptor clash, so this would otherwise loop forever.
    void checkStartupGeometry();

    std::string baseName;
    std::uint32_t seed = 0;

    std::vector<std::unique_ptr<World>> worlds_;         // heap so World& stays valid
    std::vector<double> temperatures_;                   // one per replica
    std::vector<std::vector<robo::Vec3>> replicaCoords_; // per replica, nm, OpenMM order
    int writeCounter_ = 0;
    std::vector<dcd::Writer> dcdWriters_; // one trajectory per replica (baseName.<r>.dcd)
    std::vector<double> dcdScratch_;      // interleaved xyz, nm->Angstrom, reused
    std::mt19937_64 rexRng_;
    std::uniform_real_distribution<double> rexUniform_{0.0, 1.0};
};