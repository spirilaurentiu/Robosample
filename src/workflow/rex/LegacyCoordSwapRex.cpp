#include "Context.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "GibbsSweep.hpp"
#include "ReplicaExchange.hpp"
#include "RexInternal.hpp"

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
                const bool accepted = robo::gibbs::stepWorldInGround(*w, replicaCoords_[r], systemTopology.numAtoms);

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
