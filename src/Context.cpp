#include "Context.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
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

    dcdWriters_.clear();
    dcdWriters_.reserve(temperatures_.size());
    for (std::size_t r = 0; r < temperatures_.size(); ++r) {
        dcdWriters_.emplace_back();
        dcdWriters_.back().initialize(baseName + "." + std::to_string(r) + ".dcd",
                                      systemTopology.numAtoms,
                                      /*withBox*/ false);
    }
}

// ----------------------------------------------------------------------------
//  checkStartupGeometry -- refuse to start (or warn) on an unminimized system.
//  Set the env var ROBO_ALLOW_BAD_START=1 to downgrade the hard error to a
//  warning (e.g. if you intend to run a Cartesian relaxation world first).
// ----------------------------------------------------------------------------
void Context::checkStartupGeometry() {
    const int n = systemTopology.numAtoms;
    if (n <= 0) {
        return;
    }

    // Bonded/1-2 pairs are excluded from the clash scan (a bond length ~0.1 nm
    // is not a clash). We only exclude direct bonds here -- cheap and sufficient
    // to remove the obvious false positives.
    std::set<std::pair<int, int>> bonded;
    for (int k = 0; k < systemTopology.numBonds; ++k) {
        const int i = systemTopology.bondsI[k], j = systemTopology.bondsJ[k];
        bonded.insert({std::min(i, j), std::max(i, j)});
    }

    // Hard-clash distance. Two non-bonded heavy/H atoms closer than this are in
    // r^-12 overlap. 0.08 nm (0.8 A) is well inside any real contact (vdW
    // contacts are >= ~0.2 nm) yet above bonded H distances we already excluded.
    constexpr double kClashNm = 0.08;

    // O(N^2) is fine: this runs once. (For very large systems a grid would help.)
    int nClash = 0;
    double minNonbondedNm = 1e30;
    std::pair<int, int> worst{-1, -1};
    const auto& X = systemTopology.atomsX;
    const auto& Y = systemTopology.atomsY;
    const auto& Z = systemTopology.atomsZ;
    bool anyNaN = false;
    for (int a = 0; a < n; ++a) {
        if (!std::isfinite(X[a]) || !std::isfinite(Y[a]) || !std::isfinite(Z[a])) {
            anyNaN = true;
        }
    }
    for (int a = 0; a < n && !anyNaN; ++a) {
        for (int b = a + 1; b < n; ++b) {
            if (bonded.count({a, b})) {
                continue;
            }
            const double dx = X[a] - X[b], dy = Y[a] - Y[b], dz = Z[a] - Z[b];
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (d < minNonbondedNm) {
                minNonbondedNm = d;
                worst = {a, b};
            }
            if (d < kClashNm) {
                ++nClash;
            }
        }
    }

    const double pe = calcOpenMMPotentialEnergy();
    const bool peBad = !std::isfinite(pe) || pe > 1.0e4; // +1e4 kJ/mol is already pathological

    std::fprintf(stderr,
                 "[init] startup geometry check: nAtoms=%d  initial PE=%.4g kJ/mol  "
                 "min non-bonded distance=%.4f nm  hard clashes(<%.2f nm)=%d\n",
                 n,
                 pe,
                 (minNonbondedNm > 1e29 ? 0.0 : minNonbondedNm),
                 kClashNm,
                 nClash);

    if (!anyNaN && nClash == 0 && !peBad) {
        return; // healthy start
    }

    std::string msg = "Context::initialize: the input structure is NOT usable as-is.\n";
    if (anyNaN) {
        msg += "  * coordinates contain NaN/Inf.\n";
    }
    if (nClash > 0 && worst.first >= 0) {
        const auto nm = [&](int i) {
            return (i < (int)systemTopology.atomsUniqueName.size()) ? systemTopology.atomsUniqueName[i]
                                                                    : std::to_string(i);
        };
        char buf[256];
        std::snprintf(buf,
                      sizeof(buf),
                      "  * %d steric clash(es): non-bonded atoms closer than %.2f nm "
                      "(worst: %s -- %s at %.4f nm).\n",
                      nClash,
                      kClashNm,
                      nm(worst.first).c_str(),
                      nm(worst.second).c_str(),
                      minNonbondedNm);
        msg += buf;
    }
    if (peBad) {
        char buf[160];
        std::snprintf(buf,
                      sizeof(buf),
                      "  * initial potential energy is %.4g kJ/mol (clashing/unminimized).\n",
                      pe);
        msg += buf;
    }
    msg += "  The docking world welds the receptor RIGID, so it cannot relax these\n"
           "  clashes -- the run would loop forever rejecting kicks. Energy-MINIMIZE\n"
           "  the structure first (e.g. tleap/sander/OpenMM LocalEnergyMinimizer, or a\n"
           "  Cartesian relaxation world before the docking world). To proceed anyway,\n"
           "  set ROBO_ALLOW_BAD_START=1.\n";

    const char* allow = std::getenv("ROBO_ALLOW_BAD_START");
    if (allow != nullptr && allow[0] != '0' && allow[0] != '\0') {
        std::fprintf(stderr, "[init] WARNING (continuing, ROBO_ALLOW_BAD_START set):\n%s", msg.c_str());
        return;
    }
    throw std::runtime_error(msg);
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

    // -----------------------------------------------------------------------
    //  PRE-RUN INITIAL KICK: for every docking world that has
    //  maxInitialKickTries > 0, find a clash-free starting pose for each
    //  replica BEFORE round 0.  The found position is written back into
    //  replicaCoords_ so the first generateSample() receives it directly.
    //  Temperature is set to the replica's temperature for the sphere-size
    //  geometry (no thermal dependence, but keeps the API consistent).
    for (auto& w : worlds_) {
        if (!w->isDocking()) {
            continue;
        }
        for (int r = 0; r < R; ++r) {
            w->setTemperature(temperatures_[r]);
            w->setAtomsLocationsInGround(replicaCoords_[r]);
            w->findGoodStartingPose(); // no-op if maxInitialKickTries == 0
            // Write the (possibly repositioned) ligand coords back into
            // replicaCoords_ so every subsequent world in the sweep sees them.
            const robo::Vec3* p = w->getAtomsLocationsInGround();
            std::copy(p, p + systemTopology.numAtoms, replicaCoords_[r].begin());
        }
    }

    for (int round = 0; round < total; ++round) {
        const bool production = round >= equilRounds;
        const bool logThisRound = (writeFreq > 0) && (round % writeFreq == 0);

        // Flip every world into/out of equilibration mode. During equil all
        // worlds use AlwaysAccept; during production each world uses its own
        // configured acceptRejectMode. The timestep and mdSteps are unchanged.
        for (auto& w : worlds_) {
            w->setEquilPhase(!production);
        }
        if (round == 0 && equilRounds > 0) {
            std::fprintf(stderr,
                         "[rex] === EQUILIBRATION phase: rounds 0..%d "
                         "(AlwaysAccept on all worlds) ===\n",
                         equilRounds - 1);
        }
        if (round == equilRounds && equilRounds > 0) {
            std::fprintf(stderr,
                         "[rex] === PRODUCTION phase: rounds %d..%d "
                         "(each world uses its configured acceptRejectMode) ===\n",
                         equilRounds,
                         total - 1);
        }

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