// ============================================================================
//  HmcMove.cpp - World's HMC move core.
//
//  Relocated verbatim from World.cpp (SPLIT-W9, pure code motion): the last
//  and most central World concern -- the momentum draw + starting-Hamiltonian
//  assembly (reinitialize), the total-energy assembly at a proposed state
//  (currentTotalEnergy), the acceptance test (metropolis), and the move
//  dispatcher (generateSample, which routes NCMC/Cartesian/torsional/docking
//  and runs the constrained-Verlet trajectory).
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "RobotIntegrator.hpp" // RobotEngine::stepTo/checkReversibility
#include "robo/world/detail/nma_debug.hpp"

using robo::Real;

// ----------------------------------------------------------------------------
//  Sampling
// ----------------------------------------------------------------------------
bool World::generateSample() {
    if (sampler_.moveType == MoveType::NcmcSwitch) {
        return ncmcMove();
    }
    if (cartesian_) {
        lastKickApplied_ = false;
        savedPosG_.assign(state_.atomPosG(), state_.atomPosG() + model_.numAtoms);
        bridge_.setAtomPositionsInGround(state_);
        const double peOld = bridge_.calcPotentialEnergy();
        bridge_.setVelocitiesToTemperature(temperature_, static_cast<int>(rng_()));
        bridge_.integrateTrajectoryOnDevice(state_, sampler_.mdSteps, sampler_.timeStep);
        const double peNew = bridge_.calcPotentialEnergy();
        state_.energy.pe = peNew;
        state_.energy.ke = 0.0; // device kinetic energy not pulled back here
        state_.energy.fixman = 0.0;
        state_.energy.logSineSqrGamma2 = 0.0;
        state_.energy.total = peNew;
        if (metropolis(peOld, peNew)) {
            lastAccepted_ = true;
            return true;
        }
        std::copy(savedPosG_.begin(), savedPosG_.end(), state_.atomPosG());
        state_.energy.pe = peOld;
        state_.energy.total = peOld;
        lastAccepted_ = false;
        return false;
    }

    // Internal-coordinate Generalized-Coordinate HMC over the (here: external)
    // DOF. The kick (above) is part of THIS move's proposal, not a separate
    // accept/reject: we save the pre-kick q + energy and reference Hold to the
    // pre-kick potential, so an overlapTorsiong placement is penalised by dH AND
    // caught by the energy validity check below -- and the whole move is rejected
    // back to the clean pre-kick state.
    savedQ_.assign(state_.q(), state_.q() + model_.nq); // pre-move (pre-kick) q
    lastKickApplied_ = false;
    double dockPotPre = 0.0; // pre-kick potential part of H
    double pePre = 0.0, fixPre = 0.0, lssPre = 0.0;
    if (docking_) {
        // Save the pre-kick CARTESIAN pose. The kick calls setAtomsLocationsInGround,
        // which redefines the body reference frames from the (kicked) positions and
        // zeroes q -- so restoring saved q would reconstruct the KICKED pose, not
        // this one. The pose, in Ground coordinates, is the reliable thing to keep.
        savedPosG_.assign(state_.atomPosG(), state_.atomPosG() + model_.numAtoms);

        RobotEngine::realizePosition(model_, state_);
        RobotEngine::realizeArticulatedBodyInertias(model_, state_);
        bridge_.evaluate(state_);
        pePre = bridge_.calcPotentialEnergy();
        if (sampler_.useFixman) {
            fixPre = calcFixman();
        }
        if (sampler_.useOrientationJacobian) {
            lssPre = calcLogSineSqrGamma2();
        }
        dockPotPre = pePre + fixPre - (0.5 * RT_ * lssPre);

        // Force a kick whenever the carried-forward state is unusable as a starting
        // pose. Three triggers, OR'd together:
        //   (a) non-finite pePre                -- nothing can integrate from NaN/Inf;
        //   (b) |pePre| > sampler_.maxStartPE   -- a clash; tied to the SAME absolute
        //       ceiling the acceptance gate uses below, so no admitted pose escapes
        //       the rescue (closes the old dead band where a +9800 clash sat under a
        //       1e4 rescue ceiling yet passed the relative validity gate);
        //   (c) dockingStuckCount_ >= maxStuckRounds -- a finite, sub-ceiling pose
        //       that is nonetheless LOCALLY NON-INTEGRABLE (every reseeded trajectory
        //       diverges). Neither (a) nor (b) catches it, so without this counter the
        //       move loops forever ("PE frozen, q restored"). This makes the trap
        //       escapable independent of energy magnitude.
        const bool stuck = (dockingStuckCount_ >= sampler_.maxStuckRounds);
        // A clash is a HIGH POSITIVE potential (steric overlap, r^-12). A bound
        // pose is strongly NEGATIVE (e.g. -2400 kJ/mol) and is exactly what we
        // want to keep -- it must NOT be flagged "bad". The old test
        // |pePre| > 1e3 force-kicked every well-bound pose, scrambling it every
        // round ("From -2400 -> wrong conf"). Gate on the positive ceiling only,
        // and use the configured maxStartPE (not a hard-coded 1e3) so a strained
        // but immovable receptor offset does not by itself trip the rescue.
        const bool preIsBad = !std::isfinite(pePre) || (pePre > sampler_.maxStartPE);
        std::fprintf(stderr,
                     "[dock] pre-kick: PE_old=%.2f  Fix_old=%.2f  H_old=%.2f kJ/mol  "
                     "triggers: always_kick=%s  preIsBad=%s(maxStartPE=%.0f)  "
                     "stuck=%s(%d/%d)  -> will_kick=%s\n",
                     pePre,
                     fixPre,
                     dockPotPre,
                     sampler_.alwaysKick ? "Y" : "N",
                     preIsBad ? "Y" : "N",
                     sampler_.maxStartPE,
                     stuck ? "Y" : "N",
                     dockingStuckCount_,
                     sampler_.maxStuckRounds,
                     (sampler_.alwaysKick || preIsBad || stuck) ? "Y" : "N");
        lastKickApplied_ = repositionLigands(sampler_.alwaysKick || preIsBad || stuck);
        if (stuck) {
            dockingStuckCount_ = 0; // fresh window after the forced shake
        }
    }

    reinitialize(); // seed u (KE), record Hold_ at the (post-kick) config
    if (docking_) {
        // Re-reference Hold to the PRE-kick potential. KE is config-independent
        // (= 1/2 RT |g|^2), so the only kick-dependent term is the potential.
        Hold_ = dockPotPre + state_.energy.ke;
        state_.energy.total = Hold_;
    }

    // Real pre-trajectory energy components, for an HONEST acceptance log and a
    // correct clash gate. reinitialize() already drew the momenta and folded the
    // resulting kinetic energy into Hold_ -- so KE_old is NOT zero; the metropolis
    // test compares the full Hold_ (PE + KE + Fixman + J) against Hnew. For a
    // torsional world these are the freshly-seeded values; for docking PE/Fix are
    // referenced to the pre-kick pose (KE is configuration-independent).
    const double peOld = docking_ ? pePre : state_.energy.pe;
    const double keOld = state_.energy.ke;
    const double fixOld = docking_ ? fixPre : state_.energy.fixman;

    const Real h = sampler_.timeStep;

    // PRE-STEP CLASH SCREEN (docking). reinitialize() has just evaluated the
    // post-kick forces/PE. If the kick drove the guest into a hard overlap, the
    // *change* in potential is enormous (or already non-finite). Handing such a
    // pose to the Verlet integrator is fatal: a single step against an ~Inf LJ
    // force turns finite q into NaN -- the "Particle coordinate is NaN" crash
    // seen when the proposal overlaps the receptor. So reject BEFORE stepTorsiong.
    // Gate on the move-INDUCED change (pePost - pePre), not the absolute total:
    // the rigid receptor may carry a large constant internal energy (e.g. an
    // unminimized protein at ~1e7 kJ/mol) that the ligand world neither created
    // nor can remove, and which cancels in the difference.
    bool stepsOk = true;
    if (docking_) {
        const double pePost = state_.energy.pe; // set by reinitialize()
        const double dPEpost = pePost - pePre;
        const bool screenedOut = !std::isfinite(pePost) || (dPEpost > sampler_.maxStartPE);
        std::fprintf(stderr,
                     "[dock] post-kick proposal (pre-MD):\n"
                     "[dock]   PE_old  = %12.2f kJ/mol\n"
                     "[dock]   PE_new  = %12.2f kJ/mol\n"
                     "[dock]   dPE     = %+12.2f kJ/mol  (new-old, clash ceiling=%.0f)\n"
                     "[dock]   pre-step screen: %s\n",
                     pePre,
                     pePost,
                     dPEpost,
                     sampler_.maxStartPE,
                     screenedOut ? "REJECT (skip MD, dPE>ceiling or non-finite)" : "pass -> run MD");
        if (screenedOut) {
            stepsOk = false;
        }
    }

    // Periodic reversibility probe (THEORY 5.7). Non-destructive: it integrates
    // mdSteps forward + back at this world's dt from the freshly seeded state and
    // restores the state, so the real proposal below is unaffected. It is a
    // smoke test for the CURRENT geometry only -- the always-on guard is the
    // per-step corrector throw inside verletStep. Disabled when interval == 0.
    const long revRound = generateSampleCalls_++;
    if (sampler_.reversibilityCheckInterval > 0 && stepsOk && sampler_.mdSteps > 0
        && (revRound % sampler_.reversibilityCheckInterval == 0)) {
        const Real revResid =
            RobotEngine::checkReversibility(model_, state_, bridge_, constraints_, sampler_.mdSteps, h);
        const Real revTol = Real(1e-6); // relative round-trip residual; ~1e-12 for a clean step
        const bool revBad = !std::isfinite(revResid) || revResid > revTol;
        std::fprintf(stderr,
                     "[rev] world %d round %ld: round-trip residual = %.3e over %d steps at dt=%.6g ps%s\n",
                     index_,
                     revRound,
                     (double)revResid,
                     sampler_.mdSteps,
                     (double)sampler_.timeStep,
                     revBad ? "  <-- WARNING: integrator not reversible at this dt/geometry; reduce timestep"
                            : "  (ok)");
    }

    for (int i = 0; stepsOk && i < sampler_.mdSteps; ++i) {
        // stepTo returns false if a non-finite force/coordinate appeared mid-step;
        // bail immediately so the broken pose is rejected, never carried forward.
        stepsOk = RobotEngine::stepTo(model_, state_, bridge_, constraints_, state_.time + h);
    }

    bool finite = stepsOk;
    if (finite) {
        const Real* qchk = state_.q();
        for (int i = 0; i < model_.nq; ++i) {
            if (!std::isfinite(qchk[i])) {
                finite = false;
                break;
            }
        }
    }

    bool accepted = false;
    if (finite) {
        const double Hnew = currentTotalEnergy(); // sets state_.energy (pe, ke, ...)
        const double peNew = state_.energy.pe;
        const double keNew = state_.energy.ke;
        const double fixNew = state_.energy.fixman;
        const double dPE = peNew - peOld;
        const double dKE = keNew - keOld;
        const double dFix = fixNew - fixOld;
        const double dH = Hnew - Hold_;
        const bool valid = std::isfinite(Hnew) && std::isfinite(peNew) && (dPE <= sampler_.maxStartPE);
        const bool mhPass = metropolis(Hold_, Hnew);
        const char* tag = docking_ ? "dock" : "hmc";
        std::fprintf(stderr,
                     "[%s] post-MD decision:\n"
                     "[%s]   PE_old  = %12.2f   PE_new  = %12.2f   dPE  = %+12.2f kJ/mol\n"
                     "[%s]   KE_old  = %12.2f   KE_new  = %12.2f   dKE  = %+12.2f kJ/mol\n"
                     "[%s]   Fix_old = %12.2f   Fix_new = %12.2f   dFix = %+12.2f kJ/mol\n"
                     "[%s]   H_old   = %12.2f   H_new   = %12.2f   dH   = %+12.2f kJ/mol\n"
                     "[%s]   valid(dPE<=%.0f)=%s  metropolis=%s  -> %s\n",
                     tag,
                     tag,
                     peOld,
                     peNew,
                     dPE,
                     tag,
                     keOld,
                     keNew,
                     dKE,
                     tag,
                     fixOld,
                     fixNew,
                     dFix,
                     tag,
                     Hold_,
                     Hnew,
                     dH,
                     tag,
                     sampler_.maxStartPE,
                     valid ? "Y" : "N",
                     mhPass ? "Y" : "N",
                     (valid && mhPass) ? "ACCEPT" : "reject");
        if (valid && mhPass) {
            RobotEngine::fillAtomPositionsFromBodies(model_, state_);
            accepted = true;
        }
    } else if (!stepsOk) {
        std::fprintf(stderr, "[hmc] post-MD: non-finite force mid-step or pre-screen -> reject\n");
    } else {
        std::fprintf(stderr, "[hmc] q NON-FINITE -> restoring\n");
    }

    if (accepted) {
        if (docking_) {
            dockingStuckCount_ = 0; // moved successfully -> not stuck
        }
        lastAccepted_ = true;
        return true;
    }

    // Rejected (non-finite, clash, or MH): restore the pre-kick state so nothing
    // broken is carried to the next round or reported in the log.
    if (docking_) {
        // Re-establish frames + q + positions from the saved pre-kick Cartesian
        // pose (same entry point the kick used). Restoring q alone is WRONG here:
        // the kick redefined the reference frames, so old q no longer maps to the
        // old pose. This is the fix for the "stuck in a clash forever" failure.
        setAtomsLocationsInGround(savedPosG_);
        state_.energy.pe = pePre;
        state_.energy.ke = 0.0;
        state_.energy.fixman = fixPre;
        state_.energy.logSineSqrGamma2 = lssPre;
        state_.energy.total = dockPotPre;
        // Count this stuck round. When it crosses maxStuckRounds the next entry
        // forces an unconditional kick (above), so no pose can trap the run.
        ++dockingStuckCount_;
    } else {
        std::copy(savedQ_.begin(), savedQ_.end(), state_.q());
        RobotEngine::realizePosition(model_, state_);
        RobotEngine::fillAtomPositionsFromBodies(model_, state_);
    }
    lastAccepted_ = false;
    return false;
}

void World::reinitialize() {
    RobotEngine::realizePosition(model_, state_);
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);

    // Seed u = sqrt(RT) * sqrt(M^-1) * g, g ~ N(0, I). multiplyBySqrtMInv never
    // forms M or M^-1: it sweeps the per-body DI (D^-1) blocks (Jain O(n)).
    const int nu = model_.nu;
    std::vector<Real> g(nu), seeded(nu);
    for (int i = 0; i < nu; ++i) {
        g[i] = gaussian_(rng_);
    }

    // simtk NMA Route B (DistortOption::NMA): draw the momentum from a symmetric
    // two-component Gaussian MIXTURE biased along the NMA direction instead of the
    // isotropic Gaussian. Us = z + s*mu, s = +/-1 uniform, mu = alpha * uhat where
    // uhat = uScaleFactors_/||uScaleFactors_|| is a UNIT direction and alpha =
    // nmaBiasScale is the directed push in thermal sigmas (so ||mu||^2 = alpha^2).
    // The map u = sqrt(RT) * M^-1/2 * Us below is unchanged; only the seed differs.
    // The bias steers proposals along soft directions; detailed balance is restored
    // in the acceptance by ke_mix = ke - RT*ln cosh(w.mu) (see nmaKineticCorrection).
    // With alpha=0 the mixture collapses to the plain draw. nullopt (None) => skip.
    if (sampler_.distortOption == DistortOption::NMA) {
        if (static_cast<int>(uScaleFactors_.size()) != nu) {
            uScaleFactors_.assign(nu, Real{1});
        }
        // Bias mu = alpha * uhat, where uhat = uScaleFactors_ / ||uScaleFactors_|| is a
        // UNIT direction and alpha = nmaBiasScale is the directed push in thermal-sigma
        // units. Hence ||mu||^2 = alpha^2 (NOT nu): the bias injects only ~1/2 RT alpha^2
        // of directed energy, so alpha controls boldness vs acceptance directly.
        Real norm2 = 0;
        for (int i = 0; i < nu; ++i) {
            norm2 += uScaleFactors_[i] * uScaleFactors_[i];
        }
        const Real alpha = sampler_.nmaBiasScale;
        const Real unitScale = (norm2 > 0) ? (alpha / std::sqrt(norm2)) : Real{0};
        nmaBias_.assign(nu, Real{0});
        for (int i = 0; i < nu; ++i) {
            nmaBias_[i] = uScaleFactors_[i] * unitScale;
        }

        // Symmetric mixture sign s = +/-1.
        const Real s = (uniform_(rng_) < 0.5) ? Real{-1} : Real{1};

        const bool trace = nmaDebugEnabled();
        Real zN2 = 0, muDotZ = 0;
        if (trace) {
            for (int i = 0; i < nu; ++i) {
                zN2 += g[i] * g[i];
                muDotZ += nmaBias_[i] * g[i];
            }
        }

        // Shift the white noise by the signed bias: Us = z + s*mu (in place in g).
        Real muDotUs = 0, UsN2 = 0;
        for (int i = 0; i < nu; ++i) {
            g[i] += s * nmaBias_[i];
            muDotUs += nmaBias_[i] * g[i];
            UsN2 += g[i] * g[i];
        }

        if (trace) {
            const int k = std::min(nu, 8);
            std::cout << "[nma] world " << index_ << ": Route B mixture draw, nu=" << nu << ", sign s=" << s
                      << ", alpha=" << alpha << "\n";
            std::cout << "[nma]   mu=alpha*uhat (first " << k << "): ";
            for (int i = 0; i < k; ++i) {
                std::cout << nmaBias_[i] << ' ';
            }
            std::cout << (k < nu ? "...\n" : "\n");
            std::cout << "[nma]   ||mu||^2=" << (alpha * alpha)
                      << " (==alpha^2; directed energy ~1/2 RT alpha^2=" << (0.5 * RT_ * alpha * alpha)
                      << "), ||z||^2=" << zN2 << ", mu.z=" << muDotZ << "\n";
            std::cout << "[nma]   Us=z+s*mu: ||Us||^2=" << UsN2 << ", mu.Us=" << muDotUs
                      << "  (mu.Us is the START w.mu; multiplyBySqrtM must reproduce it)\n"
                      << std::flush;
        }
    }

    RobotEngine::multiplyBySqrtMInv(model_, state_, g.data(), seeded.data());
    const Real scale = std::sqrt(RT_);
    Real* u = state_.u();
    for (int i = 0; i < nu; ++i) {
        u[i] = scale * seeded[i];
    }

    // Draw the Cartesian solvent velocities from the same Maxwell-Boltzmann
    // marginal (independent of the generalized draw -- flat diagonal metric).
    drawSolventVelocities();

    if (!constraints_.empty()) {
        RobotEngine::realizeVelocity(model_, state_);
        constraints_.enforceVelocityConstraints(model_, state_);
    }

    bridge_.evaluate(state_);
    RobotEngine::realizeVelocity(model_, state_);
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);
    RobotEngine::calcUDot(model_, state_);
    RobotEngine::calcQDot(model_, state_, state_.qdot());
    RobotEngine::calcQDotDot(model_, state_);

    const double pe = bridge_.calcPotentialEnergy();
    const double ke = RobotEngine::calcKineticEnergy(model_, state_);
    const double nmaCorr = nmaKineticCorrection(); // RT*ln cosh(w.mu); 0 unless Route B
    if (sampler_.distortOption == DistortOption::NMA && nmaDebugEnabled()) {
        // Physical KE = 1/2 u^T M u (equipartition target nu/2*RT), and the Route B
        // kinetic ke_mix = ke - nmaCorr that actually enters the acceptance H.
        std::cout << "[nma]   START: ke=1/2 u^T M u=" << ke << " (target nu/2*RT=" << (0.5 * nu * RT_)
                  << "), ke_mix=ke-corr=" << (ke - nmaCorr) << "\n"
                  << std::flush;
    }
    double fixman = 0.0;
    double logSineSqr = 0.0;
    if (sampler_.useFixman) {
        fixman = calcFixman();
    }
    if (sampler_.useOrientationJacobian) {
        logSineSqr = calcLogSineSqrGamma2();
    }
    const double keSolvent = calcSolventKE(); // 0 when no Cartesian solvent
    state_.energy.pe = pe;
    state_.energy.ke = ke;
    state_.energy.keSolvent = keSolvent;
    state_.energy.fixman = fixman;
    state_.energy.logSineSqrGamma2 = logSineSqr;
    // ke_mix = ke - nmaCorr replaces the kinetic term for the NMA Route B mixture
    // draw (nmaCorr == 0 for ordinary HMC, so Hold_ is unchanged off Route B).
    // keSolvent is the flat-space solvent kinetic energy (0 off solvent-relaxing
    // NCMC), so Hold_ stays identical to the welded path when there is no solvent.
    Hold_ = pe + ke + keSolvent + fixman - (0.5 * RT_ * logSineSqr) - nmaCorr;
    state_.energy.total = Hold_;
}

double World::currentTotalEnergy() {
    bridge_.evaluate(state_); // positions -> OpenMM -> forces (+ PE available)
    const double pe = bridge_.calcPotentialEnergy();
    RobotEngine::realizeVelocity(model_, state_);
    const double ke = RobotEngine::calcKineticEnergy(model_, state_);
    const double nmaCorr = nmaKineticCorrection(); // RT*ln cosh(w.mu) at the END; 0 unless Route B
    if (sampler_.distortOption == DistortOption::NMA && nmaDebugEnabled()) {
        std::cout << "[nma]   END:   ke=1/2 u^T M u=" << ke << ", ke_mix=ke-corr=" << (ke - nmaCorr)
                  << " (corr=RT*ln cosh(w.mu)=" << nmaCorr << ")\n"
                  << std::flush;
    }
    double fixman = 0.0;
    double logSineSqr = 0.0;
    if (sampler_.useFixman) {
        fixman = calcFixman(); // realizes ABI internally; position already current
    }
    if (sampler_.useOrientationJacobian) {
        logSineSqr = calcLogSineSqrGamma2();
    }
    const double keSolvent = calcSolventKE(); // 0 when no Cartesian solvent
    state_.energy.pe = pe;
    state_.energy.ke = ke;
    state_.energy.keSolvent = keSolvent;
    state_.energy.fixman = fixman;
    state_.energy.logSineSqrGamma2 = logSineSqr;
    // ke_mix = ke - nmaCorr (Route B mixture); nmaCorr == 0 off Route B. keSolvent
    // is the flat-space solvent KE (0 off solvent-relaxing NCMC).
    state_.energy.total = pe + ke + keSolvent + fixman - (0.5 * RT_ * logSineSqr) - nmaCorr;
    return state_.energy.total;
}

bool World::metropolis(double Hold, double Hnew) {
    // During burn-in (equilPhase_) every move is accepted regardless of the
    // world's configured acceptRejectMode. This lets the system relax from the
    // starting geometry (which may be far from equilibrium) without the
    // Metropolis gate blocking large-dH moves. The sampler configuration is
    // otherwise unchanged -- same timestep, same mdSteps, same kick logic --
    // so switching to production is a single flag flip with no other state change.
    if (equilPhase_ || sampler_.acceptRejectMode == AcceptRejectMode::AlwaysAccept) {
        return true;
    }
    const double dH = Hnew - Hold;
    if (dH <= 0) {
        return true;
    }
    return uniform_(rng_) < std::exp(-beta_ * dH);
}
