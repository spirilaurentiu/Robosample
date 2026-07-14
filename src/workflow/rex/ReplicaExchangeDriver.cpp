#include "Context.hpp"

#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <utility>
#include <vector>

#include "GibbsSweep.hpp"
#include "ReplicaExchange.hpp"

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
                    const bool accepted =
                        robo::gibbs::stepWorldInGround(w, replicas_[r].atomsLocations, systemTopology.numAtoms);

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
