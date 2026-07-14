#include "Context.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

#include "GibbsSweep.hpp"
#include "ReplicaExchange.hpp"

// ---------------------------------------------------------------------------
//  Stage 2b: INV-7/INV-10 preconditions, the driven (RENE/REBASONTOP) round,
//  and the REBASONTOP interleave (D4). NOT compiled or run (coordinator
//  directive, 2026-07-12) -- reviewed on paper only.
// ---------------------------------------------------------------------------
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
            const bool accepted =
                robo::gibbs::stepWorldInGround(w, replicas_[r].atomsLocations, systemTopology.numAtoms);

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
