#include "Context.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "ReplicaExchange.hpp"
#include "RexInternal.hpp"

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
