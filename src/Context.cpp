#include "Context.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "DCDWriter.hpp"
#include "OpenMM.h"

namespace {
constexpr double kBoltzmann_kJ = 0.0083144626; // kJ/mol/K
}

Context::Context(std::string baseName, std::uint32_t seed)
    : baseName(std::move(baseName))
    , seed(seed)
    , rexRng_(seed ^ 0xD1B54A32D192ED03ULL) {
}

void Context::setRootMobility(int moleculeIndex, RootMobility mobility) {
    if (moleculeIndex < 0 || moleculeIndex >= static_cast<int>(systemTopology.rootMobilities.size())) {
        throw std::out_of_range("setRootMobility: molecule index " + std::to_string(moleculeIndex)
                                + " out of range (have "
                                + std::to_string(systemTopology.rootMobilities.size()) + " molecules)");
    }
    systemTopology.rootMobilities[moleculeIndex] = mobility;
}

// ---------------------------------------------------------------------------
//  Worlds
// ---------------------------------------------------------------------------
World& Context::addCartesianWorld() {
    const int idx = static_cast<int>(worlds_.size());
    worlds_.push_back(std::make_unique<World>(idx, /*cartesian*/ true, seed));
    worlds_.back()->buildModel(systemTopology, Selection{}, systemTopology.rootMobilities);
    return *worlds_.back();
}

World& Context::addRoboticWorld(const Selection& sel) {
    const int idx = static_cast<int>(worlds_.size());
    worlds_.push_back(std::make_unique<World>(idx, /*cartesian*/ false, seed));
    worlds_.back()->buildModel(systemTopology, sel, systemTopology.rootMobilities);
    return *worlds_.back();
}

World& Context::addDockingWorld(const std::vector<int>& ligandMoleculeIndices) {
    const int numMol = systemTopology.numMolecules;
    std::set<int> ligandSet(ligandMoleculeIndices.begin(), ligandMoleculeIndices.end());

    // Per-WORLD root mobilities: ligands Free, everything else Welded. This is
    // local to the docking world and does NOT touch systemTopology.rootMobilities.
    std::vector<RootMobility> rootMob(numMol, RootMobility::Weld);
    for (int m : ligandMoleculeIndices) {
        if (m < 0 || m >= numMol) {
            throw std::out_of_range("addDockingWorld: ligand molecule index " + std::to_string(m)
                                    + " out of range (have " + std::to_string(numMol) + " molecules)");
        }
        rootMob[m] = RootMobility::Free;
    }

    // Rigid-body docking: every bond stays Rigid (each molecule = one body).
    Selection sel;
    sel.bondMobility.assign(systemTopology.numBonds, BondMobility::Rigid);

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
                                      BondMobility mobility,
                                      bool /*flag*/) {
    Selection sel;
    sel.bondMobility.assign(systemTopology.numBonds, BondMobility::Rigid);

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

// ---------------------------------------------------------------------------
//  Initialize + run
// ---------------------------------------------------------------------------
void Context::initialize(const std::vector<double>& temperatures) {
    temperatures_ = temperatures.empty() ? std::vector<double>{300.0} : temperatures;

    if (!initializeOpenMM()) {
        throw std::runtime_error("Context::initialize: OpenMM initialization failed");
    }

    std::vector<robo::Vec3> ref(systemTopology.numAtoms);
    for (int a = 0; a < systemTopology.numAtoms; ++a) {
        ref[a] = robo::Vec3(systemTopology.atomsX[a], systemTopology.atomsY[a], systemTopology.atomsZ[a]);
    }
    replicaCoords_.assign(temperatures_.size(), ref);
    writeCounter_ = 0;

    dcdWriters_.clear();
    dcdWriters_.reserve(temperatures_.size());
    for (std::size_t r = 0; r < temperatures_.size(); ++r) {
        dcdWriters_.emplace_back();
        dcdWriters_.back().initialize(baseName + "." + std::to_string(r) + ".dcd",
                                      systemTopology.numAtoms,
                                      /*withBox*/ false);
    }
}

void Context::runREX(int equilRounds, int prodRounds, int writeFreq, bool verbose) {
    const int R = static_cast<int>(temperatures_.size());
    const int total = equilRounds + prodRounds;

    // Per-world telemetry CSV (one row per world per replica per logged round).
    const std::string movesPath = baseName + ".moves.csv";
    {
        std::ofstream hdr(movesPath, std::ios::app);
        if (hdr.tellp() == 0) {
            hdr << "round,replica,T,world,type,kick,accepted,PE,KE,fixman,H\n";
        }
    }

    for (int round = 0; round < total; ++round) {
        const bool production = round >= equilRounds;
        const bool logThisRound = (writeFreq > 0) && (round % writeFreq == 0);

        // Gibbs sweep: each replica runs every world in order.
        for (int r = 0; r < R; ++r) {
            for (auto& w : worlds_) {
                w->setTemperature(temperatures_[r]);
                w->setAtomsLocationsInGround(replicaCoords_[r]);
                const bool accepted = w->generateSample();
                const robo::Vec3* p = w->getAtomsLocationsInGround();
                std::copy(p, p + systemTopology.numAtoms, replicaCoords_[r].begin());

                if (verbose) {
                    std::printf("[rex] round=%d replica=%d T=%.1f world=%d(%s) kick=%s acc=%s "
                                "PE=%.4f KE=%.4f Fix=%.4f H=%.4f\n",
                                round,
                                r,
                                temperatures_[r],
                                w->index(),
                                w->typeName(),
                                w->lastKickApplied() ? "Y" : "N",
                                accepted ? "ACC" : "rej",
                                w->lastPE(),
                                w->lastKE(),
                                w->lastFixman(),
                                w->lastTotalEnergy());
                }
                if (logThisRound) {
                    std::ofstream out(movesPath, std::ios::app);
                    if (out) {
                        out << round << ',' << r << ',' << temperatures_[r] << ',' << w->index() << ','
                            << w->typeName() << ',' << (w->lastKickApplied() ? 1 : 0) << ','
                            << (accepted ? 1 : 0) << ',' << w->lastPE() << ',' << w->lastKE() << ','
                            << w->lastFixman() << ',' << w->lastTotalEnergy() << '\n';
                    }
                }
            }
        }

        if (R > 1) {
            const int startPair = (round % 2 == 0) ? 0 : 1;
            for (int r = startPair; r + 1 < R; r += 2) {
                const double Ea = openmmPotential(replicaCoords_[r]);
                const double Eb = openmmPotential(replicaCoords_[r + 1]);
                const double betaA = 1.0 / (kBoltzmann_kJ * temperatures_[r]);
                const double betaB = 1.0 / (kBoltzmann_kJ * temperatures_[r + 1]);
                const double delta = (betaA - betaB) * (Ea - Eb);
                if (delta >= 0 || rexUniform_(rexRng_) < std::exp(delta)) {
                    std::swap(replicaCoords_[r], replicaCoords_[r + 1]);
                }
            }
        }

        if (production && writeFreq > 0 && ((round - equilRounds) % writeFreq == 0)) {
            for (int r = 0; r < R; ++r) {
                writeOutputs(r, round, verbose);
            }
            ++writeCounter_;
        }
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

void Context::writeOutputs(int replica, int round, bool verbose) {
    const double pe = openmmPotential(replicaCoords_[replica]);
    const std::string path = baseName + "." + std::to_string(replica) + ".csv";
    std::ofstream out(path, std::ios::app);
    if (out) {
        out << round << "," << replica << "," << temperatures_[replica] << "," << pe << "\n";
    }
    if (verbose) {
        std::printf("[rex] round=%d replica=%d T=%.1f PE=%.4f kJ/mol\n",
                    round,
                    replica,
                    temperatures_[replica],
                    pe);
    }

    if (replica >= 0 && replica < static_cast<int>(dcdWriters_.size())) {
        const auto& coords = replicaCoords_[replica];
        const int n = systemTopology.numAtoms;
        const auto& perm = systemTopology.atomsPrmtopIndex;
        const bool havePerm = (static_cast<int>(perm.size()) == n);
        dcdScratch_.resize(static_cast<std::size_t>(3 * n));
        for (int a = 0; a < n; ++a) {
            const int p = havePerm ? perm[a] : a;
            dcdScratch_[3 * p + 0] = coords[a][0] * 10.0;
            dcdScratch_[3 * p + 1] = coords[a][1] * 10.0;
            dcdScratch_[3 * p + 2] = coords[a][2] * 10.0;
        }
        dcdWriters_[replica].append(dcdScratch_);
    }
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