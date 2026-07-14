#include "Context.hpp"

#include <algorithm>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "DCDWriter.hpp"
#include "OpenMM.h"

Context::Context(std::string baseName, std::uint32_t seed)
    : baseName(std::move(baseName))
    , seed(seed)
    , rexRng_(seed ^ 0xD1B54A32D192ED03ULL) {
}

// Context::setRootMobility was removed: root mobility is now a per-world
// property. Use World::setRootMobility / World::setRootMobilities on the world
// returned by add*World instead (see World.cpp).

// ---------------------------------------------------------------------------
//  Worlds
// ---------------------------------------------------------------------------
World& Context::addCartesianWorld(bool wantReactionReporter) {
    const int idx = static_cast<int>(worlds_.size());
    worlds_.push_back(std::make_unique<World>(idx, /*cartesian*/ true, seed));
    World& world = *worlds_.back();
    world.buildModel(systemTopology, Selection{}, systemTopology.rootMobilities);
    if (wantReactionReporter) {
        // Always throws (Sec.3 integrator guard): a Cartesian world's
        // internal-coordinate articulated reaction is not meaningful.
        world.enableReactionReporter();
    }
    return world;
}

World& Context::addRoboticWorld(const Selection& sel, bool wantReactionReporter) {
    const int idx = static_cast<int>(worlds_.size());
    worlds_.push_back(std::make_unique<World>(idx, /*cartesian*/ false, seed));
    World& world = *worlds_.back();
    world.buildModel(systemTopology, sel, systemTopology.rootMobilities);
    if (wantReactionReporter) {
        world.enableReactionReporter();
    }
    return world;
}

World& Context::addDockingWorld(const std::vector<int>& ligandMoleculeIndices) {
    const int numMol = systemTopology.numMolecules;
    std::set<int> ligandSet(ligandMoleculeIndices.begin(), ligandMoleculeIndices.end());

    // Per-WORLD root mobilities: ligands Free, everything else Welded. This is
    // local to the docking world and does NOT touch systemTopology.rootMobilities.
    std::vector<JointType> rootMob(numMol, JointType::Rigid);
    for (int m : ligandMoleculeIndices) {
        if (m < 0 || m >= numMol) {
            throw std::out_of_range("addDockingWorld: ligand molecule index " + std::to_string(m)
                                    + " out of range (have " + std::to_string(numMol) + " molecules)");
        }
        rootMob[m] = JointType::Free;
    }

    // Rigid-body docking: every bond stays Rigid (each molecule = one body).
    Selection sel;
    sel.bondMobility.assign(systemTopology.numBonds, JointType::Rigid);

    const int idx = static_cast<int>(worlds_.size());
    worlds_.push_back(std::make_unique<World>(idx, /*cartesian*/ false, seed));
    worlds_.back()->buildModel(systemTopology, sel, rootMob);

    // One atom group PER ligand molecule (each gets its own auto-sized sphere and
    // is repositioned independently); the receptor atoms (all non-ligands) define
    // the sphere centre, as GLOBAL atom indices.
    std::vector<std::vector<int>> ligandGroups;
    std::vector<int> siteAtoms;
    for (int mol = 0; mol < numMol; ++mol) {
        const int begin = systemTopology.atomsBegin[mol];
        const int end = systemTopology.atomsEnd[mol];
        if (ligandSet.count(mol)) {
            std::vector<int> grp;
            for (int a = begin; a < end; ++a) {
                grp.push_back(a);
            }
            ligandGroups.push_back(std::move(grp));
        } else {
            for (int a = begin; a < end; ++a) {
                siteAtoms.push_back(a);
            }
        }
    }
    worlds_.back()->configureDocking(std::move(ligandGroups), std::move(siteAtoms));
    return *worlds_.back();
}

Selection Context::buildFlexibilities(const std::optional<std::vector<std::pair<int, int>>>& bonds,
                                      JointType mobility,
                                      bool /*flag*/) {
    Selection sel;
    sel.bondMobility.assign(systemTopology.numBonds, JointType::Rigid);

    std::set<std::pair<int, int>> want;
    if (bonds.has_value()) {
        for (const auto& b : *bonds) {
            want.insert({std::min(b.first, b.second), std::max(b.first, b.second)});
        }
    }

    for (int k = 0; k < systemTopology.numBonds; ++k) {
        if (systemTopology.bondsRingClosing[k]) {
            continue;
        }
        const int i = systemTopology.bondsI[k];
        const int j = systemTopology.bondsJ[k];
        const bool rotatable =
            systemTopology.atomsNumBondsInvolved[i] >= 2 && systemTopology.atomsNumBondsInvolved[j] >= 2;
        if (!rotatable) {
            continue;
        }
        if (bonds.has_value() && want.count({std::min(i, j), std::max(i, j)}) == 0) {
            continue;
        }
        sel.bondMobility[k] = mobility;
    }
    return sel;
}

void Context::setMTS(bool enabled, int innerSubsteps) {
    OpenMMContext::get().setMTS(enabled, innerSubsteps);
}

void Context::setSeparateForceGroups(bool enabled) {
    OpenMMContext::get().setSeparateForceGroups(enabled);
}

void Context::setEnforcePeriodicBox(bool enabled) {
    OpenMMContext::get().setEnforcePeriodicBox(enabled);
}

// ---------------------------------------------------------------------------
//  Initialize + run
// ---------------------------------------------------------------------------
void Context::initialize(const std::vector<double>& temperatures) {
    temperatures_ = temperatures.empty() ? std::vector<double>{300.0} : temperatures;

    if (!initializeOpenMM()) {
        throw std::runtime_error("Context::initialize: OpenMM initialization failed");
    }

    // ------------------------------------------------------------------------
    //  STARTUP HEALTH CHECK. A docking run inherits the input geometry as the
    //  carried-forward reference pose; if that pose is already broken (NaN, or a
    //  steric clash), the docking world can NEVER repair it (the receptor is
    //  welded -- rigid -- so its internal clashes are frozen), and every round
    //  degenerates to "kick, reject, restore, log the same PE". Catch it here,
    //  loudly, BEFORE any sampling, instead of letting it silently waste a run.
    //
    //  Two independent signals, because PE magnitude alone is ambiguous (it is
    //  EXTENSIVE -- a clean explicit-solvent box is legitimately ~-1e6 kJ/mol):
    //    1. Hard steric clashes: non-bonded atom pairs closer than kClashNm.
    //       This is INTENSIVE and solvent/size-independent -- the most reliable
    //       "not minimized" signal. A single such pair contributes ~1e6-1e8
    //       kJ/mol of r^-12 repulsion.
    //    2. Non-finite or extremely positive total PE.
    checkStartupGeometry();

    std::vector<robo::Vec3> ref(systemTopology.numAtoms);
    for (int a = 0; a < systemTopology.numAtoms; ++a) {
        ref[a] = robo::Vec3(systemTopology.atomsX[a], systemTopology.atomsY[a], systemTopology.atomsZ[a]);
    }
    replicaCoords_.assign(temperatures_.size(), ref);
    writeCounter_ = 0;

    // Truncate the per-run TEXT outputs so a rerun with the same base name starts
    // fresh. They are opened in append mode as frames are produced (writeReactionRows /
    // the energy + moves CSVs) with a "write the header only when the file is empty"
    // guard, so WITHOUT truncating here a second run appends its rows onto the first
    // run's file -- accumulating stale/duplicate frames (and inflating any downstream
    // analysis). The DCD writers already overwrite on initialize(); this brings the
    // CSVs in line. Opening an ofstream in trunc mode (and letting it close) empties
    // the file, or creates it empty.
    std::ofstream(baseName + ".moves.csv", std::ios::trunc);
    for (std::size_t r = 0; r < temperatures_.size(); ++r) {
        std::ofstream(baseName + "." + std::to_string(r) + ".csv", std::ios::trunc);
        std::ofstream(baseName + "." + std::to_string(r) + ".reactions.csv", std::ios::trunc);
    }

    dcdWriters_.clear();
    dcdWriters_.reserve(temperatures_.size());
    const bool periodicDcd =
        OpenMMContext::isPeriodic(systemTopology.nonbondedMethod) && systemTopology.boxVectors.size() == 9;
    for (std::size_t r = 0; r < temperatures_.size(); ++r) {
        dcdWriters_.emplace_back();
        dcdWriters_.back().initialize(baseName + "." + std::to_string(r) + ".dcd",
                                      systemTopology.numAtoms,
                                      /*withBox*/ periodicDcd);
    }
}

double Context::openmmPotential(const std::vector<robo::Vec3>& coords) {
    std::vector<OpenMM::Vec3> pos;
    pos.reserve(coords.size());
    for (const auto& c : coords) {
        pos.emplace_back(c[0], c[1], c[2]);
    }
    return OpenMMContext::get().computePotentialEnergy(pos);
}

// ---------------------------------------------------------------------------
//  OpenMM energy ingestion / validation (unchanged)
// ---------------------------------------------------------------------------
auto Context::initializeOpenMM() -> bool {
    return OpenMMContext::get().initialize(systemTopology);
}

auto Context::calcOpenMMPotentialEnergy() -> double {
    std::vector<OpenMM::Vec3> positions;
    positions.reserve(static_cast<std::size_t>(systemTopology.numAtoms));
    for (int i = 0; i < systemTopology.numAtoms; ++i) {
        positions.emplace_back(systemTopology.atomsX[i], systemTopology.atomsY[i], systemTopology.atomsZ[i]);
    }
    return OpenMMContext::get().computePotentialEnergy(positions);
}

auto Context::computePotentialEnergyByGroup()
    -> std::pair<double, std::vector<OpenMMContext::ForceGroupEnergy>> {
    std::vector<OpenMM::Vec3> positions;
    positions.reserve(static_cast<std::size_t>(systemTopology.numAtoms));
    for (int i = 0; i < systemTopology.numAtoms; ++i) {
        positions.emplace_back(systemTopology.atomsX[i], systemTopology.atomsY[i], systemTopology.atomsZ[i]);
    }
    return OpenMMContext::get().computePotentialEnergyByGroup(positions);
}