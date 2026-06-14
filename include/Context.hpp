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
    // Override a molecule's root attachment to Ground (default Free, set by
    // load_amber). Writes systemTopology.rootMobilities, which the worlds read.
    void setRootMobility(int moleculeIndex, RootMobility mobility);

    // Add a Cartesian world (pure OpenMM MD on device) or a robotic/torsional
    // world (internal-coordinate HMC). Both return a reference to the new world
    // so Python can chain .add_sampler(...). The model is built immediately
    // from the current systemTopology + root mobilities.
    World& addCartesianWorld();
    World& addRoboticWorld(const Selection& sel);

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
                                 BondMobility mobility,
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