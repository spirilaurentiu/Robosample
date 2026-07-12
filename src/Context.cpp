#include "Context.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "DCDWriter.hpp"
#include "OpenMM.h"
#include "PeriodicBox.hpp"

namespace {
constexpr double kBoltzmann_kJ = 0.0083144626; // kJ/mol/K

// Convert the three reduced lattice vectors (nm, lower-triangular, row-major
// a,b,c) into a CHARMM/DCD Box: side lengths in Angstrom and angles in degrees.
dcd::Box boxFromReducedVectors(const std::vector<double>& bv) {
    dcd::Box box;
    if (bv.size() != 9) {
        return box; // default 1 A cube
    }
    const double ax = bv[0];
    const double bx = bv[3], by = bv[4];
    const double cx = bv[6], cy = bv[7], cz = bv[8];
    const double a = ax;
    const double b = std::sqrt(bx * bx + by * by);
    const double c = std::sqrt(cx * cx + cy * cy + cz * cz);
    constexpr double kRad2Deg = 57.29577951308232;
    const double cosGamma = (a > 0 && b > 0) ? (ax * bx) / (a * b) : 0.0;
    const double cosBeta = (a > 0 && c > 0) ? (ax * cx) / (a * c) : 0.0;
    const double cosAlpha = (b > 0 && c > 0) ? (bx * cx + by * cy) / (b * c) : 0.0;
    constexpr double kNm2Ang = 10.0;
    box.sideA = a * kNm2Ang;
    box.sideB = b * kNm2Ang;
    box.sideC = c * kNm2Ang;
    box.angleAlpha = std::acos(std::max(-1.0, std::min(1.0, cosAlpha))) * kRad2Deg;
    box.angleBeta = std::acos(std::max(-1.0, std::min(1.0, cosBeta))) * kRad2Deg;
    box.angleGamma = std::acos(std::max(-1.0, std::min(1.0, cosGamma))) * kRad2Deg;
    return box;
}

} // namespace

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

    // Excluded pairs: not just 1-2 bonds, but the FULL intramolecular
    // non-interacting set -- the OpenMM exclusion list (1-2 and 1-3) and the 1-4
    // (scaling14) pairs. These are bonded geometry, never clashes. Using only
    // bonds is wrong for 4-point water: the extra point sits ~0.078 nm from each
    // hydrogen (a 1-3 pair that OpenMM excludes), so a bonds-only filter reports
    // two phantom "clashes" per water (2*Nwater of them). The clash scan is meant
    // to catch INTERMOLECULAR overlaps from a bad/unminimized placement, not the
    // internal geometry of a virtual-site water. Keyed as min*N+max in a hash set.
    std::unordered_set<long long> excluded;
    auto pairKey = [n](int i, int j) -> long long {
        if (i > j) {
            std::swap(i, j);
        }
        return static_cast<long long>(i) * n + j;
    };
    for (int k = 0; k < systemTopology.numBonds; ++k) {
        excluded.insert(pairKey(systemTopology.bondsI[k], systemTopology.bondsJ[k]));
    }
    for (int k = 0; k < systemTopology.numExclusions; ++k) {
        excluded.insert(pairKey(systemTopology.exclusionI[k], systemTopology.exclusionJ[k]));
    }
    for (int k = 0; k < systemTopology.numScaling14; ++k) {
        excluded.insert(pairKey(systemTopology.scaling14I[k], systemTopology.scaling14L[k]));
    }

    // Virtual sites / massless particles (e.g. a 4-point water's extra point) have
    // no Lennard-Jones term -- no steric presence -- so they cannot clash and are
    // skipped entirely. (Belt-and-suspenders with the exclusion set above.)
    std::vector<bool> isVirtual(static_cast<std::size_t>(n), false);
    for (int a = 0; a < n; ++a) {
        if (systemTopology.atomsMass[a] == 0.0) {
            isVirtual[static_cast<std::size_t>(a)] = true;
        }
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

    // Under explicit solvent the box is periodic, so a "distance" must be the
    // MINIMUM-IMAGE distance: two atoms on opposite faces are actually neighbours.
    // box_vectors are reduced (lower-triangular) a=(ax,0,0) b=(bx,by,0) c=(cx,cy,cz);
    // wrap the displacement by subtracting whole lattice vectors in c,b,a order
    // (the same order OpenMM reduces them). Without this the scan reports phantom
    // clashes (or misses real cross-boundary ones) on a solvated box.
    const bool periodic =
        OpenMMContext::isPeriodic(systemTopology.nonbondedMethod) && systemTopology.boxVectors.size() == 9;
    const auto& bv = systemTopology.boxVectors;
    auto minImage = [&](double& dx, double& dy, double& dz) {
        if (!periodic) {
            return;
        }
        // single source of truth: robo::pbc::minimumImage (PeriodicBox.hpp).
        robo::pbc::minimumImage(dx, dy, dz, bv.data());
    };

    bool anyNaN = false;
    for (int a = 0; a < n; ++a) {
        if (!std::isfinite(X[a]) || !std::isfinite(Y[a]) || !std::isfinite(Z[a])) {
            anyNaN = true;
        }
    }
    for (int a = 0; a < n && !anyNaN; ++a) {
        if (isVirtual[static_cast<std::size_t>(a)]) {
            continue;
        }
        for (int b = a + 1; b < n; ++b) {
            if (isVirtual[static_cast<std::size_t>(b)]) {
                continue;
            }
            if (excluded.count(pairKey(a, b)) != 0) {
                continue;
            }
            double dx = X[a] - X[b], dy = Y[a] - Y[b], dz = Z[a] - Z[b];
            minImage(dx, dy, dz);
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

                // Reaction-force reporter (docs/specs/reaction-force-monitoring.md):
                // one per-body force snapshot per reporter world, per replica,
                // on every DCD-write round -- the same cadence writeOutputs
                // just used.
                for (auto& w : worlds_) {
                    if (!w->isReactionReporter()) {
                        continue;
                    }
                    // Re-sync to replica r's FINAL coordinates (post adjacent-
                    // replica swap above, if any) -- the SAME conformation
                    // writeOutputs just wrote to the DCD frame -- so the CSV row
                    // and the DCD frame are guaranteed to agree (Sec.1 item 3)
                    // even when a swap relabelled replicaCoords_ this round.
                    // Read-only: the next thing that touches this world is
                    // always a fresh setAtomsLocationsInGround at the top of its
                    // next generateSample() turn, so this cannot perturb the
                    // sampler (Sec.3/6).
                    w->setAtomsLocationsInGround(replicaCoords_[r]);
                    w->captureReactionSnapshot();
                    writeReactionRows(r, writeCounter_, w->reactionSamples());
                }
            }
            ++writeCounter_;
        }
    }
}

// ---------------------------------------------------------------------------
//  Label-swap replica exchange (docs/specs/replica-exchange-nonequilibrium-
//  work.md). Stage 1: RUN_TYPE::REMC only.
// ---------------------------------------------------------------------------
void Context::setupReplicaExchange(RUN_TYPE runType) {
    runType_ = runType;

    const int R = static_cast<int>(temperatures_.size());
    const int W = static_cast<int>(worlds_.size());
    if (R == 0) {
        throw std::runtime_error("Context::RunREX: no replicas -- call initialize(temperatures) first");
    }

    replicas_.assign(R, Replica{});
    thermodynamicStates_.assign(R, ThermodynamicState{});
    for (int i = 0; i < R; ++i) {
        // Seeded from the same reference coordinates initialize() already put
        // in replicaCoords_ (B0: both drivers start from an identical state).
        replicas_[i].atomsLocations = replicaCoords_[i];

        ThermodynamicState& st = thermodynamicStates_[i];
        st.temperature = temperatures_[i];
        st.worldIndexes.resize(W);
        st.timeSteps.resize(W);
        st.mdSteps.resize(W);
        st.acceptRejectModes.resize(W);
        for (int w = 0; w < W; ++w) {
            st.worldIndexes[w] = w; // B0: every state schedules the SAME full ordered world list.
            st.timeSteps[w] = worlds_[w]->getTimeStep();
            st.mdSteps[w] = worlds_[w]->getMdSteps();
            st.acceptRejectModes[w] = worlds_[w]->getAcceptRejectMode();
        }
    }

    // PRE-RUN INITIAL KICK (S1 fix, mirrors runREX Context.cpp:409-429
    // verbatim): for every docking world with maxInitialKickTries>0, find a
    // clash-free starting pose for each replica BEFORE round 0. Without this
    // RunREX starts every replica from the unrepaired seed on a docking
    // config (e.g. FFAR1), diverging from the coordinate-swap oracle and
    // risking a clashing first generateSample(); it is a no-op for
    // non-docking systems (findGoodStartingPose no-ops when
    // maxInitialKickTries == 0), which is why the ala-dipeptide
    // INVARIANT-EQUIV run didn't exercise it.
    for (auto& w : worlds_) {
        if (!w->isDocking()) {
            continue;
        }
        for (int i = 0; i < R; ++i) {
            w->setTemperature(temperatures_[i]);
            w->setAtomsLocationsInGround(replicas_[i].atomsLocations);
            w->findGoodStartingPose(); // no-op if maxInitialKickTries == 0
            const robo::Vec3* p = w->getAtomsLocationsInGround();
            std::copy(p, p + systemTopology.numAtoms, replicas_[i].atomsLocations.begin());
        }
    }

    // Committed potentials, evaluated AFTER the pose-repair loop above so
    // INV-6 (stored potential == energy of stored coordinates) holds from
    // round 0 even when a docking world just moved the seed coordinates.
    for (int i = 0; i < R; ++i) {
        replicas_[i].potential = openmmPotential(replicas_[i].atomsLocations);
        replicas_[i].referencePotential = replicas_[i].potential;
    }

    // B0 NOTE: "the model is W shared worlds, R = T replicas/states". In
    // Stage 1 this is a STRUCTURAL invariant, not a runtime check:
    // replicas_ and thermodynamicStates_ are both `.assign(R, ...)` from the
    // SAME local `R` above, so replicas_.size() == thermodynamicStates_.size()
    // holds by construction and no input combination here can violate it --
    // an `if` guard here would be permanently dead code (reviewer S2). This
    // stops being vacuous, and SHALL regain a real runtime check, the moment
    // a future change (Stage 2's ThermodynamicState list, e.g. an explicit
    // addThermodynamicState() API per the original design) sources the state
    // count from a SECOND, independent input.

    // Identity maps (B0).
    replica2ThermoIxs_.resize(R);
    thermo2ReplicaIxs_.resize(R);
    for (int i = 0; i < R; ++i) {
        replica2ThermoIxs_[i] = i;
        thermo2ReplicaIxs_[i] = i;
    }

    nofAttemptedSwapsMatrix_.assign(R, std::vector<std::int64_t>(R, 0));
    nofAcceptedSwapsMatrix_.assign(R, std::vector<std::int64_t>(R, 0));
    exchangeRound_ = 0;
    exchangePairList_.clear();
}

void Context::swapThermodynamicStates(int thermoC, int thermoH) {
    const int X = thermo2ReplicaIxs_[thermoC];
    const int Y = thermo2ReplicaIxs_[thermoH];
    std::swap(replica2ThermoIxs_[X], replica2ThermoIxs_[Y]);
    std::swap(thermo2ReplicaIxs_[thermoC], thermo2ReplicaIxs_[thermoH]);
}

bool Context::attemptREXSwap(int thermoC, int thermoH) {
    // Fail loud on misuse (e.g. called from Python before RunREX/
    // setupReplicaExchange has built the objects, or with an out-of-range
    // state index) rather than reading/writing past the end of an empty
    // vector.
    if (thermodynamicStates_.empty()) {
        throw std::logic_error(
            "Context::attemptREXSwap: no thermodynamic states -- call RunREX (run_rex_label_swap) "
            "at least once first");
    }
    nofAttemptedSwapsMatrix_.at(thermoC).at(thermoH) += 1;
    nofAttemptedSwapsMatrix_.at(thermoH).at(thermoC) += 1;

    const int X = thermo2ReplicaIxs_.at(thermoC);
    const int Y = thermo2ReplicaIxs_.at(thermoH);
    const double betaC = 1.0 / (kBoltzmann_kJ * thermodynamicStates_.at(thermoC).temperature);
    const double betaH = 1.0 / (kBoltzmann_kJ * thermodynamicStates_.at(thermoH).temperature);
    const double refU_Xset = replicas_.at(X).referencePotential;
    const double refU_Yset = replicas_.at(Y).referencePotential;

    // correctionTerm (B6 step 4, D2/INV-9): the Hastings ratio of the
    // scale-factor proposal densities. == 1 (log == 0) PROVIDED the shared,
    // frozen, state-independent anchor precondition (INV-9) holds -- which is
    // exactly what Context::batAnchorStats_ (ONE global instance, not
    // per-state) and driveReplica's "one frozen snapshot per round, passed to
    // BOTH partners" discipline guarantee. Kept as an explicit named zero
    // (not simply omitted) so a future stochastic scale-factor randomiser
    // (perturbScalingFactor, deliberately NOT ported, D2) has an obvious
    // place to add its own log-density-ratio term instead of silently
    // reusing this one.
    const double logCorrectionTerm = 0.0;

    double logPAccept = 0.0;
    switch (runType_) {
        case RUN_TYPE::REMC:
            // ETerm_equal = -[(beta_H - beta_C)(refU_Xset - refU_Yset)] (B6 step 3).
            logPAccept = -((betaH - betaC) * (refU_Xset - refU_Yset));
            break;
        case RUN_TYPE::RENEMC: {
            // ETerm_nonequil (B6 step 3): the SAME parallel-tempering form,
            // but on the DRIVEN-ENDPOINT reference potentials -- no Jacobian
            // (INV-10: RENEMC's velocity/NMA drive is volume-preserving, so
            // there is none to carry). NOTE (Stage 2c TODO): the round-loop
            // that actually POPULATES referenceWORK_potential for RENEMC (the
            // velocity/NMA driven segment) is not wired yet (RunREX throws
            // for RUN_TYPE::RENEMC) -- this branch is exercised directly by
            // tests/TestRexAcceptanceAlgebra.cpp, which sets
            // referenceWORK_potential by hand.
            const double refU_Xtau = replicas_.at(X).referenceWORK_potential;
            const double refU_Ytau = replicas_.at(Y).referenceWORK_potential;
            logPAccept = -((betaH - betaC) * (refU_Xtau - refU_Ytau)) + logCorrectionTerm;
            break;
        }
        case RUN_TYPE::RENE:
        case RUN_TYPE::REBASONTOP: {
            // WTerm (B6 step 3, Derivation sketch; D7-simplified single-step
            // work since mdSteps==0 makes x^tau == x'):
            //   Work_partner = beta_target*U(x_partner^tau)
            //                - beta_source*U(x_partner^0) - lnJac_partner
            //   WTerm = -(Work_X + Work_Y)
            // X currently occupies thermoC (source beta_C, driving TOWARD
            // thermoH's temperature, so its endpoint is scored at beta_H);
            // Y occupies thermoH (source beta_H, driving toward thermoC,
            // scored at beta_C). This is the Ballard-Jarzynski / Nilmeier
            // deterministic-map acceptance (Derivation sketch; V8 proves
            // detailed balance for this exact form).
            const double refU_Xtau = replicas_.at(X).referenceWORK_potential;
            const double refU_Ytau = replicas_.at(Y).referenceWORK_potential;
            const double lnJacX = replicas_.at(X).WORK_Jacobian;
            const double lnJacY = replicas_.at(Y).WORK_Jacobian;
            const double workX = (betaH * refU_Xtau) - (betaC * refU_Xset) - lnJacX;
            const double workY = (betaC * refU_Ytau) - (betaH * refU_Yset) - lnJacY;
            logPAccept = -(workX + workY) + logCorrectionTerm;
            break;
        }
        default: // RUN_TYPE::Default
            // mixReplicas/runDrivenRound never reach here for Default (gated
            // out, B7); a direct call is a caller error, not a silent no-op.
            throw std::logic_error("Context::attemptREXSwap: called with RUN_TYPE::Default");
    }

    // NaN/inf fail-loud (reviewer N2): a non-finite acceptance exponent --
    // e.g. driveReplica forced WORK_Jacobian to -infinity after a Stage 2a
    // domain-invalid drive (r1<=0 or theta1 outside (0,pi)), or an OpenMM PE
    // blew up on a driven endpoint -- is an EXPLICIT automatic reject, not a
    // silent comparison. (IEEE754 already makes `NaN >= 0` and `u < exp(NaN)`
    // both false, and `-inf` already rejects naturally too, but a stray
    // "+inf, wrong sign" case from an unanticipated bug would NOT reject
    // naturally -- this guard forces the correct, conservative outcome
    // regardless of sign and makes the event visible.)
    bool forcedReject = false;
    if (std::isnan(logPAccept) || std::isinf(logPAccept)) {
        std::fprintf(stderr,
                     "[rexlabel] WARNING: non-finite acceptance exponent (logPAccept=%g) for "
                     "thermoC=%d thermoH=%d, runType=%d -- automatic reject (reviewer N2)\n",
                     logPAccept,
                     thermoC,
                     thermoH,
                     static_cast<int>(runType_));
        forcedReject = true;
    }

    const bool accept =
        !forcedReject && ((logPAccept >= 0.0) || (rexUniform_(rexRng_) < std::exp(logPAccept)));
    if (accept) {
        nofAcceptedSwapsMatrix_[thermoC][thermoH] += 1;
        nofAcceptedSwapsMatrix_[thermoH][thermoC] += 1;
        if (runType_ == RUN_TYPE::RENEMC || runType_ == RUN_TYPE::RENE || runType_ == RUN_TYPE::REBASONTOP) {
            // F4 atomic commit (INV-4): promote BOTH replicas' WORK_* trial
            // to committed together (coords + potential + referencePotential
            // + FixmanPotential), before the label swap. REMC/Default never
            // reach here (their accept is label-swap only, B6 step 7 -- this
            // IS the F4 fix: the original unconditionally ran the WORK commit
            // even for REMC, reverting coordinates from an unpopulated WORK
            // buffer).
            replicas_.at(X).commitWorkAsFinal();
            replicas_.at(Y).commitWorkAsFinal();
        }
        swapThermodynamicStates(thermoC, thermoH); // label swap only (INV-3)
    }
    return accept;
}

void Context::prepareExchangePairs(int round, int oddity) {
    exchangePairList_.clear();
    const int K = static_cast<int>(thermodynamicStates_.size());
    const int startIdx = (round + oddity) % 2;
    for (int thIx = startIdx; thIx + 1 < K; thIx += 2) {
        exchangePairList_.emplace_back(thIx, thIx + 1);
    }
}

void Context::mixAllReplicas(int nAttempts) {
    const int T = static_cast<int>(thermodynamicStates_.size());
    if (T <= 1) {
        return;
    }
    std::uniform_int_distribution<int> pick(0, T - 1);
    for (int a = 0; a < nAttempts; ++a) {
        int i = pick(rexRng_);
        int j = pick(rexRng_);
        while (j == i) {
            j = pick(rexRng_);
        }
        attemptREXSwap(i, j);
    }
}

void Context::mixReplicas(int mixi) {
    if (swapEvery_ <= 0 || (mixi % swapEvery_) != 0) {
        return;
    }
    const int T = static_cast<int>(thermodynamicStates_.size());
    if (runType_ == RUN_TYPE::Default || T <= 1) {
        return;
    }
    if (mixingScheme_ == ReplicaMixingScheme::Neighboring) {
        // Parity from the dedicated exchangeRound_ counter, NOT mixi (B7
        // revision-2 fix): keeps both pairing parities reachable regardless
        // of swapEvery_.
        prepareExchangePairs(exchangeRound_, /*oddity=*/0);
        for (const auto& pr : exchangePairList_) {
            attemptREXSwap(pr.first, pr.second);
        }
        ++exchangeRound_; // once per EXECUTED mix
    } else {
        mixAllReplicas(nSwapAttempts_);
    }
}

void Context::printSwapMatrix() const {
    const int T = static_cast<int>(thermodynamicStates_.size());
    std::fprintf(stderr, "[rexlabel] accepted/attempted swap matrix (by thermodynamic-state index):\n");
    for (int i = 0; i < T; ++i) {
        std::fprintf(stderr, "[rexlabel]  ");
        for (int j = 0; j < T; ++j) {
            std::fprintf(stderr,
                        "%4lld/%-4lld ",
                        static_cast<long long>(nofAcceptedSwapsMatrix_[i][j]),
                        static_cast<long long>(nofAttemptedSwapsMatrix_[i][j]));
        }
        std::fprintf(stderr, "\n");
    }
}

// ---------------------------------------------------------------------------
//  Stage 2b: INV-7/INV-10 preconditions, the driven (RENE/REBASONTOP) round,
//  and the REBASONTOP interleave (D4). NOT compiled or run (coordinator
//  directive, 2026-07-12) -- reviewed on paper only.
// ---------------------------------------------------------------------------
void Context::checkInv7AndInv10Guards(RUN_TYPE runType) const {
    if (runType != RUN_TYPE::RENE && runType != RUN_TYPE::RENEMC && runType != RUN_TYPE::REBASONTOP) {
        return; // INV-7 (as scoped by the Stage 2b directive) and INV-10 are driven-only preconditions
    }

    // INV-7/V9 (D3 biconditional): the swap acceptance excludes Fixman ONLY
    // correctly if Fixman IS enabled in every (non-Cartesian) sampler. A
    // Cartesian world has a flat metric (no Fixman term is meaningful there;
    // World::add_sampler auto-forces useFixman=false for it), so Cartesian
    // worlds are exempt.
    for (const auto& w : worlds_) {
        if (w->isCartesian()) {
            continue;
        }
        if (!w->getUseFixman()) {
            throw std::logic_error(
                "Context::checkInv7AndInv10Guards: INV-7/V9 violated -- world "
                + std::to_string(w->index())
                + " has Fixman disabled but the run type drives (RENE/RENEMC/REBASONTOP). The "
                  "swap acceptance excludes Fixman ONLY correctly when Fixman is enabled in "
                  "every sampler (D3 biconditional); a Fixman-off torsional world here would "
                  "target the wrong joint distribution.");
        }
    }

    // INV-10 (drive/run-type pairing): RENE/REBASONTOP SHALL drive with a
    // volume-changing BAT-scaling world (distortOption == ScaleBendStretch,
    // carrying lnJac); RENEMC SHALL drive with a volume-preserving velocity/
    // NMA world (distortOption == NMA, omitting it). Any world configured
    // with the WRONG drive for the active run type biases acceptance
    // (B9/INV-1) -- reject the pairing outright.
    bool anyDriven = false;
    for (const auto& w : worlds_) {
        const auto opt = w->getDistortOption();
        if (!opt.has_value()) {
            continue;
        }
        anyDriven = true;
        if ((runType == RUN_TYPE::RENE || runType == RUN_TYPE::REBASONTOP)
            && *opt != DistortOption::ScaleBendStretch) {
            throw std::logic_error(
                "Context::checkInv7AndInv10Guards: INV-10 violated -- world "
                + std::to_string(w->index())
                + " has a non-ScaleBendStretch distortOption under RUN_TYPE::RENE/REBASONTOP.");
        }
        if (runType == RUN_TYPE::RENEMC && *opt != DistortOption::NMA) {
            throw std::logic_error(
                "Context::checkInv7AndInv10Guards: INV-10 violated -- world "
                + std::to_string(w->index())
                + " has a non-NMA distortOption under RUN_TYPE::RENEMC.");
        }
    }
    if (!anyDriven) {
        throw std::logic_error(
            "Context::checkInv7AndInv10Guards: run type drives (RENE/RENEMC/REBASONTOP) but no "
            "world has a matching distortOption configured -- the drive would be silently inert.");
    }
}

void Context::driveReplica(int replicaIx,
                           int thermoIx,
                           double targetTemperature,
                           const robo::BatAnchorStats::Snapshot& anchor) {
    Replica& rep = replicas_.at(replicaIx);
    const ThermodynamicState& st = thermodynamicStates_.at(thermoIx);
    const double s = std::sqrt(targetTemperature / st.temperature); // B4 Q-scale-factor

    rep.WORK = 0.0;
    rep.WORK_Jacobian = 0.0;
    rep.WORK_atomsLocations = rep.atomsLocations; // start from the committed (equilibrium) endpoint x^0

    bool anyDriven = false;
    for (std::size_t pos = 0; pos < st.worldIndexes.size(); ++pos) {
        const int worldIx = st.worldIndexes[pos];
        World& w = *worlds_[worldIx];
        const auto distortOpt = w.getDistortOption();
        if (!distortOpt.has_value() || *distortOpt != DistortOption::ScaleBendStretch) {
            continue; // not a BAT-scaling driven world position (checkInv7AndInv10Guards already
                      // confirmed no OTHER-typed driven world exists for RENE/REBASONTOP)
        }
        anyDriven = true;
        const double uPrev = openmmPotential(rep.WORK_atomsLocations);
        w.setAtomsLocationsInGround(rep.WORK_atomsLocations);
        try {
            w.applyBatScalingDrive(s, anchor.meanR, anchor.meanTheta);
        } catch (const std::domain_error& e) {
            // S2 (reviewer, 2026-07-12): Stage 2a's r1<=0 / theta1 outside
            // (0,pi) domain guard fired for THIS (s, anchor) pair -- an
            // invalid scaled configuration. MUST NOT propagate and abort the
            // whole REX run: map it to an automatic reject of THIS swap
            // (reviewer N2 fail-loud, applied at the swap level) and let the
            // run continue. Force the reject by driving WORK_Jacobian to
            // -infinity: attemptREXSwap's Work_partner term then reads
            //   beta_target*U(x^tau) - beta_source*U(x^0) - (-inf) = +inf,
            // so WTerm = -(Work_X+Work_Y) = -inf and the swap always rejects
            // regardless of the (now-irrelevant) potential values -- setting
            // WORK_potential to a "plausible" finite number would NOT
            // reliably force rejection (the beta_target/beta_source weights
            // differ, so a finite-but-equal U(x^tau)==U(x^0) does not
            // generally give a negative logPAccept). The trial endpoint is
            // left at the last valid geometry (x^0 for this drive) rather
            // than the invalid target.
            std::fprintf(stderr,
                         "[rexlabel] WARNING: applyBatScalingDrive domain error for replica=%d "
                         "world=%d s=%g -- forcing automatic reject of this swap (%s)\n",
                         replicaIx,
                         worldIx,
                         s,
                         e.what());
            rep.WORK_Jacobian = -std::numeric_limits<double>::infinity();
            rep.WORK_potential = uPrev;
            rep.referenceWORK_potential = uPrev;
            return;
        }
        const robo::Vec3* p = w.getAtomsLocationsInGround();
        rep.WORK_atomsLocations.assign(p, p + systemTopology.numAtoms);
        const double uCurr = openmmPotential(rep.WORK_atomsLocations);
        rep.WORK += (uCurr - uPrev); // B5 per-driven-world work (Fixman excluded, D3)
        rep.WORK_Jacobian += w.getDistortJacobianDetLog(); // B12 live += accumulation
    }

    // D7: x^tau == x' (no post-scale MD) -- the driven endpoint's potential
    // is read directly off the final scaled geometry. If no driven world was
    // actually visited (e.g. an odd-T unpaired replica reaching here, or a
    // state schedule with none), fall back to the committed values so
    // referenceWORK_potential/WORK_Jacobian stay internally consistent (a
    // null drive == REMC's ETerm_equal limit, V4).
    rep.WORK_potential = anyDriven ? openmmPotential(rep.WORK_atomsLocations) : rep.potential;
    rep.referenceWORK_potential = rep.WORK_potential;
    if (!anyDriven) {
        rep.WORK_Jacobian = 0.0;
    }
}

void Context::runInterleavedRemcSubround() {
    const RUN_TYPE saved = runType_;
    runType_ = RUN_TYPE::REMC; // borrow attemptREXSwap's ETerm_equal branch (D4)
    for (int sub = 0; sub < rebasontopSubrounds_; ++sub) {
        prepareExchangePairs(exchangeRound_, /*oddity=*/sub % 2);
        for (const auto& pr : exchangePairList_) {
            attemptREXSwap(pr.first, pr.second);
        }
        ++exchangeRound_;
    }
    runType_ = saved;
}

void Context::runDrivenRound(int round, bool verbose) {
    const int R = static_cast<int>(replicas_.size());

    // (1) Equilibrium segment (B8): every replica's non-driven worlds, in
    // schedule order, exactly like the REMC sweep -- SKIP any world whose
    // distortOption is set (that is this state's nonequilibrium segment, run
    // separately below once pairing is known).
    for (int k = 0; k < R; ++k) {
        const int r = thermo2ReplicaIxs_[k];
        const ThermodynamicState& st = thermodynamicStates_[k];
        for (std::size_t pos = 0; pos < st.worldIndexes.size(); ++pos) {
            const int worldIx = st.worldIndexes[pos];
            World& w = *worlds_[worldIx];
            if (w.getDistortOption().has_value()) {
                continue; // driven world -- nonequilibrium segment, below
            }
            w.setTemperature(st.temperature);
            w.setTimeStep(st.timeSteps[pos]);
            w.setMdSteps(st.mdSteps[pos]);
            w.setAcceptRejectMode(st.acceptRejectModes[pos]);
            w.setAtomsLocationsInGround(replicas_[r].atomsLocations);
            const bool accepted = w.generateSample();
            const robo::Vec3* p = w.getAtomsLocationsInGround();
            std::copy(p, p + systemTopology.numAtoms, replicas_[r].atomsLocations.begin());

            // INV-9 accumulation: feed the shared/global anchor from THIS
            // equilibrium sample (B3/item 4). Safe/cheap on every world
            // (BatAnchorStats::accumulate no-ops on bodies with no scaled
            // DOF, e.g. every body in a Cartesian or purely-torsional world).
            accumulateBatAnchorStats(w);

            if (verbose) {
                std::printf("[rexlabel] round=%d state=%d replica=%d T=%.1f world=%d(%s) acc=%s "
                            "PE=%.4f KE=%.4f Fix=%.4f H=%.4f [equil]\n",
                            round,
                            k,
                            r,
                            st.temperature,
                            w.index(),
                            w.typeName(),
                            accepted ? "ACC" : "rej",
                            w.lastPE(),
                            w.lastKE(),
                            w.lastFixman(),
                            w.lastTotalEnergy());
            }
        }
        // INV-6: refresh the committed potential from the equilibrium
        // endpoint x^0 (D3: no Fixman/driven term in the physical potential).
        replicas_[r].potential = openmmPotential(replicas_[r].atomsLocations);
        replicas_[r].referencePotential = replicas_[r].potential;
    }

    // (2) Pairing FIRST (B4/B7): the Q-scale-factor needs the partner's
    // temperature before any drive can run.
    prepareExchangePairs(exchangeRound_, /*oddity=*/0);
    ++exchangeRound_; // once per EXECUTED (driven) round, mirrors mixReplicas' B7 rule

    // (3) ONE frozen anchor snapshot for this round (INV-9): every drive this
    // round -- both partners of every pair -- reads the SAME snapshot, taken
    // AFTER the equilibrium segment above has fed it.
    const robo::BatAnchorStats::Snapshot anchor = batAnchorSnapshot();

    // (4) Drive each paired replica toward its partner's temperature. A
    // domain-invalid drive (S2) is handled INSIDE driveReplica (forces that
    // replica's WORK_Jacobian to -infinity, an automatic reject at step (5)
    // below) -- it never throws out of this loop.
    for (const auto& pr : exchangePairList_) {
        const int thermoC = pr.first;
        const int thermoH = pr.second;
        const double tC = thermodynamicStates_[thermoC].temperature;
        const double tH = thermodynamicStates_[thermoH].temperature;
        const int x = thermo2ReplicaIxs_[thermoC];
        const int y = thermo2ReplicaIxs_[thermoH];
        driveReplica(x, thermoC, /*targetTemperature=*/tH, anchor);
        driveReplica(y, thermoH, /*targetTemperature=*/tC, anchor);
    }

    // (5) Attempt every pair's swap (main WTerm/ETerm_nonequil acceptance).
    for (const auto& pr : exchangePairList_) {
        attemptREXSwap(pr.first, pr.second);
    }

    // (6) REBASONTOP interleave (D4): periodic REMC sub-rounds "on top".
    if (runType_ == RUN_TYPE::REBASONTOP && interleaveRemcEvery_ > 0 && (round % interleaveRemcEvery_) == 0) {
        runInterleavedRemcSubround();
    }
}

void Context::RunREX(RUN_TYPE runType, int equilRounds, int prodRounds, int writeFreq, bool verbose) {
    if (runType == RUN_TYPE::RENEMC) {
        throw std::logic_error(
            "Context::RunREX: RUN_TYPE::RENEMC's driven round-loop (the velocity/NMA drive segment) "
            "is not wired yet -- Stage 2c (docs/specs/replica-exchange-nonequilibrium-work.md). Its "
            "acceptance formula (ETerm_nonequil) IS implemented in attemptREXSwap and is directly "
            "testable there; RUN_TYPE::Default/REMC/RENE/REBASONTOP are fully wired.");
    }

    setupReplicaExchange(runType);

    const bool driven = (runType == RUN_TYPE::RENE || runType == RUN_TYPE::REBASONTOP);
    if (driven) {
        checkInv7AndInv10Guards(runType); // INV-7/V9, INV-10 -- fail loud before any round runs
    }

    const int R = static_cast<int>(replicas_.size());
    const int total = equilRounds + prodRounds;

    for (int round = 0; round < total; ++round) {
        const bool production = round >= equilRounds;

        for (auto& w : worlds_) {
            w->setEquilPhase(!production);
        }

        if (driven) {
            // RENE/REBASONTOP (Stage 2b): the driven round (B6/B7/B8) --
            // equilibrium segment, pairing, drive, WTerm swap attempt, and
            // REBASONTOP's interleave -- is entirely encapsulated in
            // runDrivenRound (see its doc comment in Context.hpp).
            runDrivenRound(round, verbose);
        } else {
            // REMC/Default (Stage 1, INVARIANT-EQUIV-tested, UNCHANGED):
            // Gibbs sweep, iterate by THERMODYNAMIC STATE (not by replica
            // identity) and look up which replica currently occupies it
            // (thermo2ReplicaIxs_). A Gibbs sweep's per-replica draw order is
            // a free choice (each replica's draw is conditionally
            // independent given the others, B0) -- iterating by state,
            // rather than by replica, makes RunREX visit the shared worlds
            // in the SAME temporal order as the legacy coordinate-swap
            // runREX (which iterates by fixed temperature slot). Since a
            // shared world's RNG stream is consumed in call order regardless
            // of which coordinates are loaded, this ordering choice is what
            // makes label-swap and coordinate-swap REMC produce IDENTICAL
            // per-round trajectories (INVARIANT-EQUIV), not merely matching
            // statistics -- both are legitimate Gibbs sweep orders, but this
            // one is directly comparable to the retained oracle. The
            // schedule's per-world timestep/mdSteps/acceptRejectMode are
            // reset onto the shared worlds every position (Consequences
            // "Feasibility gap") -- in Stage 1 these are identical across
            // states (B0 NOTE), but resetting them unconditionally keeps
            // this loop correct.
            for (int k = 0; k < R; ++k) {
                const int r = thermo2ReplicaIxs_[k];
                const ThermodynamicState& st = thermodynamicStates_[k];
                for (std::size_t pos = 0; pos < st.worldIndexes.size(); ++pos) {
                    const int worldIx = st.worldIndexes[pos];
                    World& w = *worlds_[worldIx];
                    w.setTemperature(st.temperature);
                    w.setTimeStep(st.timeSteps[pos]);
                    w.setMdSteps(st.mdSteps[pos]);
                    w.setAcceptRejectMode(st.acceptRejectModes[pos]);
                    w.setAtomsLocationsInGround(replicas_[r].atomsLocations);
                    const bool accepted = w.generateSample();
                    const robo::Vec3* p = w.getAtomsLocationsInGround();
                    std::copy(p, p + systemTopology.numAtoms, replicas_[r].atomsLocations.begin());

                    if (verbose) {
                        std::printf("[rexlabel] round=%d state=%d replica=%d T=%.1f world=%d(%s) acc=%s "
                                    "PE=%.4f KE=%.4f Fix=%.4f H=%.4f\n",
                                    round,
                                    k,
                                    r,
                                    st.temperature,
                                    w.index(),
                                    w.typeName(),
                                    accepted ? "ACC" : "rej",
                                    w.lastPE(),
                                    w.lastKE(),
                                    w.lastFixman(),
                                    w.lastTotalEnergy());
                    }
                }
            }

            // INV-6: a replica's stored potential SHALL equal the energy of
            // its stored coordinates at swap time -- refresh before
            // attemptREXSwap reads referencePotential. REMC has no Fixman/
            // driven term (D3), so potential == referencePotential.
            for (int r = 0; r < R; ++r) {
                const double pe = openmmPotential(replicas_[r].atomsLocations);
                replicas_[r].potential = pe;
                replicas_[r].referencePotential = pe;
            }

            mixReplicas(round);
        }

        if (production && writeFreq > 0 && ((round - equilRounds) % writeFreq == 0)) {
            for (int k = 0; k < R; ++k) {
                const int r = thermo2ReplicaIxs_[k];
                writeOutputsCore(k, round, verbose, replicas_[r].atomsLocations, thermodynamicStates_[k].temperature);
            }
            ++writeCounter_;
        }
    }

    printSwapMatrix();
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
    writeOutputsCore(replica, round, verbose, replicaCoords_[replica], temperatures_[replica]);
}

// Pure extraction from the former writeOutputs body -- `idx` used to be
// `replica`/`replicaCoords_[replica]`/`temperatures_[replica]` verbatim.
// RunREX (the label-swap driver) reuses this on Replica-owned coordinates,
// indexed by THERMODYNAMIC STATE rather than replica-object identity (see
// the RunREX doc comment in Context.hpp).
void Context::writeOutputsCore(int idx, int round, bool verbose, const std::vector<robo::Vec3>& coords, double T) {
    const double pe = openmmPotential(coords);
    const std::string path = baseName + "." + std::to_string(idx) + ".csv";
    std::ofstream out(path, std::ios::app);
    if (out) {
        out << round << "," << idx << "," << T << "," << pe << "\n";
    }
    if (verbose) {
        std::printf("[rex] round=%d replica=%d T=%.1f PE=%.4f kJ/mol\n", round, idx, T, pe);
    }

    if (idx >= 0 && idx < static_cast<int>(dcdWriters_.size())) {
        const int n = systemTopology.numAtoms;
        const auto& perm = systemTopology.atomsPrmtopIndex;
        const bool havePerm = (static_cast<int>(perm.size()) == n);
        dcdScratch_.resize(static_cast<std::size_t>(3 * n));

        // Whole-molecule periodic imaging, applied ONLY to the output copy.
        // replicaCoords_ stays unwrapped (contiguous per molecule) so the next
        // setAtomsLocationsInGround/frame rebuild is unaffected. Each molecule is
        // shifted by integer lattice vectors so its center of mass lands in the
        // primary cell, then translated rigidly -- bonds never straddle a face.
        const auto& bv = systemTopology.boxVectors;
        const int numMol = systemTopology.numMolecules;
        const bool haveRanges = OpenMMContext::isPeriodic(systemTopology.nonbondedMethod) && bv.size() == 9
                                && static_cast<int>(systemTopology.atomsBegin.size()) == numMol
                                && static_cast<int>(systemTopology.atomsEnd.size()) == numMol;

        auto scatter = [&](int a, double sx, double sy, double sz) {
            const int p = havePerm ? perm[a] : a;
            dcdScratch_[3 * p + 0] = (coords[a][0] + sx) * 10.0;
            dcdScratch_[3 * p + 1] = (coords[a][1] + sy) * 10.0;
            dcdScratch_[3 * p + 2] = (coords[a][2] + sz) * 10.0;
        };

        if (!haveRanges) {
            // Non-periodic (or missing ranges): write coordinates verbatim.
            for (int a = 0; a < n; ++a) {
                scatter(a, 0.0, 0.0, 0.0);
            }
        } else {
            // Reduced lower-triangular box: a=(bv0,0,0) b=(bv3,bv4,0) c=(bv6,bv7,bv8).
            // Wrap in c -> b -> a order (same convention as the clash-scan minImage).
            for (int m = 0; m < numMol; ++m) {
                const int beg = systemTopology.atomsBegin[m];
                const int end = systemTopology.atomsEnd[m];
                if (beg >= end) {
                    continue;
                }

                // Mass-weighted COM in engine order (molecule is intact here).
                double cx = 0.0, cy = 0.0, cz = 0.0, mtot = 0.0;
                const int nMass = static_cast<int>(systemTopology.atomsMass.size());
                for (int a = beg; a < end; ++a) {
                    const double mass = (a < nMass) ? systemTopology.atomsMass[a] : 1.0;
                    cx += mass * coords[a][0];
                    cy += mass * coords[a][1];
                    cz += mass * coords[a][2];
                    mtot += mass;
                }
                if (mtot > 0.0) {
                    cx /= mtot;
                    cy /= mtot;
                    cz /= mtot;
                } else { // all-massless (e.g. pure virtual sites): use first atom
                    cx = coords[beg][0];
                    cy = coords[beg][1];
                    cz = coords[beg][2];
                }

                // Accumulate the rigid shift that brings the COM into [0, L).
                double sx = 0.0, sy = 0.0, sz = 0.0;
                const double nc = (bv[8] != 0.0) ? std::floor(cz / bv[8]) : 0.0;
                cx -= nc * bv[6];
                cy -= nc * bv[7];
                sx -= nc * bv[6];
                sy -= nc * bv[7];
                sz -= nc * bv[8];
                const double nb = (bv[4] != 0.0) ? std::floor(cy / bv[4]) : 0.0;
                cx -= nb * bv[3];
                sx -= nb * bv[3];
                sy -= nb * bv[4];
                const double na = (bv[0] != 0.0) ? std::floor(cx / bv[0]) : 0.0;
                sx -= na * bv[0];

                for (int a = beg; a < end; ++a) {
                    scatter(a, sx, sy, sz);
                }
            }
        }

        dcdWriters_[idx].append(dcdScratch_, boxFromReducedVectors(systemTopology.boxVectors));
    }
}

// docs/specs/reaction-force-monitoring.md Sec.4: one CSV per replica (the
// filename carries the replica, matching the per-replica .dcd), comma-
// delimited, appended as frames are produced. No-op on an empty row set (a
// non-write round never reaches here, and a reporter world with zero
// interesting bodies would otherwise write an empty file with only a header).
// Always the same 10 columns regardless of which term(s) the reporter world
// summed into `force`/`torque` (docs/specs/reaction-force-monitoring.md
// Sec.1.2) -- that choice is made once, at World::enableReactionReporter,
// and is invisible to the CSV format.
void Context::writeReactionRows(int replica, int frame, const std::vector<ReactionSample>& rows) {
    if (rows.empty()) {
        return;
    }
    const std::string path = baseName + "." + std::to_string(replica) + ".reactions.csv";
    std::ofstream out(path, std::ios::app);
    if (!out) {
        return;
    }
    if (out.tellp() == 0) {
        out << "# frame,replica,body_idx,atom_idx,fx,fy,fz,tx,ty,tz\n";
    }
    for (const auto& s : rows) {
        out << frame << ',' << replica << ',' << s.bodyIdx << ',' << s.atomIdx << ',' << s.force[0] << ','
            << s.force[1] << ',' << s.force[2] << ',' << s.torque[0] << ',' << s.torque[1] << ',' << s.torque[2]
            << '\n';
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