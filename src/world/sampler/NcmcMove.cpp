// ============================================================================
//  NcmcMove.cpp - World's per-molecule NCMC move (docs/specs/
//  ncmc-explicit-solvent/).
//
//  Relocated verbatim from World.cpp (SPLIT-W4, pure code motion): a
//  lambda:1->0->1 alchemical decouple-move-recouple proposal, its
//  palindromic protocol schedule accessor, the lambda=0 trough teleport, the
//  Construction-II inner GHMC kernel, and the NCMC configuration entry
//  points. This is MoveType::NcmcSwitch, dispatched from generateSample (W9)
//  but otherwise a self-contained move regime.
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <vector>

#include "NCMCProtocol.hpp" // robo::ncmc::protocolLambda -- ONE schedule, prod + tests
#include "RobotIntegrator.hpp" // templated stepTo (ncmcInnerGhmcStep/ncmcMove instantiate it)
#include "engine_helpers.hpp"

using robo::Real;
using robo::Rotation;
using robo::Vec3;

namespace {

// Opt-in per-substep trace of Construction II's inner GHMC kernel
// (ncmcInnerGhmcStep). Set the env var ROBO_NCMC_DEBUG=1 to print, for every
// fixed-lambda substep, Hbefore/Hafter/dH/converged/accepted -- diagnosable
// evidence for whether the inner kernel is accepting genuine dt-dependent
// dynamics or freezing (dH pinned to a large/NaN, dt-independent value; see
// docs/specs/ncmc-explicit-solvent/). Off by default (up to ncmcSteps lines
// per move would otherwise flood stderr).
bool ncmcDebugEnabled() {
    static const bool on = [] {
        const char* e = std::getenv("ROBO_NCMC_DEBUG");
        return e != nullptr && e[0] != '0' && e[0] != '\0';
    }();
    return on;
}

} // namespace

void World::configureNcmc(std::vector<int> atomIndices, int ncmcSteps, double holdFraction) {
    // Region A is an arbitrary atom-index SET (docs/specs/ncmc-explicit-solvent/
    // 30-region-and-protocol-policy.md Sec.2). Sort + dedupe so ncmcTeleportRoot's
    // "first atom" read and the OpenMM-side index-set iteration are well defined.
    std::sort(atomIndices.begin(), atomIndices.end());
    atomIndices.erase(std::unique(atomIndices.begin(), atomIndices.end()), atomIndices.end());
    sampler_.ncmcAtomIndices = atomIndices;
    sampler_.ncmcSteps = ncmcSteps;
    sampler_.ncmcHoldFraction = holdFraction;
    sampler_.moveType = MoveType::NcmcSwitch; // must be set AFTER add_sampler
    bridge_.enableAlchemy(sampler_.ncmcAtomIndices);
}

void World::configureNcmc(int atomBegin, int atomEnd, int ncmcSteps, double holdFraction) {
    std::vector<int> atomIndices;
    if (atomEnd > atomBegin) {
        atomIndices.resize(static_cast<std::size_t>(atomEnd - atomBegin));
        std::iota(atomIndices.begin(), atomIndices.end(), atomBegin);
    }
    configureNcmc(std::move(atomIndices), ncmcSteps, holdFraction);
}

auto World::protocolLambda(int step) const -> double {
    // Delegate to the pure, unit-tested schedule (include/NCMCProtocol.hpp) so
    // production and the tests exercise ONE palindromic 1 -> 0 -> 1 schedule and
    // can never drift. Palindromic + endpoints pinned to lambda = 1 is what makes
    // the composed NCMC map F-reversible (see ncmcMove's acceptance comment).
    return robo::ncmc::protocolLambda(step, sampler_.ncmcSteps, sampler_.ncmcHoldFraction);
}

int World::ncmcTeleportRoot() const {
    // The tree root (parent == Ground) of the body carrying the FIRST (smallest-
    // index) atom of Region A, iff it is a Free joint (so a rigid reposition is
    // one quaternion+translation). ncmcAtomIndices is kept sorted ascending
    // (configureNcmc), so .front() is that atom for both the contiguous and the
    // general index-set case.
    if (sampler_.ncmcAtomIndices.empty()) {
        return -1;
    }
    const int firstAtom = sampler_.ncmcAtomIndices.front();
    if (firstAtom < 0 || firstAtom >= model_.numAtoms) {
        return -1;
    }
    int b = model_.atomBody[firstAtom];
    while (b > 0 && model_.bodyParent[b] != 0) {
        b = model_.bodyParent[b];
    }
    if (b > 0 && model_.bodyJoint[b] == JointType::Free) {
        return b;
    }
    return -1;
}

void World::ncmcApplyTroughTeleport(int rootBody) {
    // Rigid reposition of the region root at the λ=0 ghost trough. Mirrors the
    // engine-validated tests/TeleportMove.hpp teleportFreeRoot: set the root
    // orientation (Haar) and translation (uniform in the docking sphere), and
    // CO-ROTATE the root's angular AND linear speeds by ΔR = Rnew R_oldᵀ so KE is
    // exactly preserved (M_ang(q) is orientation-dependent). Operates on q/u
    // directly -- NO setAtomsLocationsInGround -- so the live momenta stay valid.
    robo::Real* q = state_.q();
    robo::Real* u = state_.u();
    const int qOff = model_.bodyQIndex[rootBody];
    const int uOff = model_.bodyUIndex[rootBody];

    const Rotation Rold = EngineHelpers::quatToRotation(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
    const Rotation Rnew = sampleUniformRotation();
    const Rotation dR(Rnew * Rold.transpose());

    const Vec3 site = atomSetCentroid(siteAtoms_);
    const double radius = atomSetRadius(siteAtoms_, site);
    const Vec3 target = site + sampleUniformInSphere(radius);

    Real w, x, y, z;
    EngineHelpers::rotationToQuaternion(Rnew, w, x, y, z);
    q[qOff + 0] = w;
    q[qOff + 1] = x;
    q[qOff + 2] = y;
    q[qOff + 3] = z;
    q[qOff + 4] = target[0];
    q[qOff + 5] = target[1];
    q[qOff + 6] = target[2];

    const Vec3 wR = dR * Vec3(u[uOff + 0], u[uOff + 1], u[uOff + 2]);
    const Vec3 vR = dR * Vec3(u[uOff + 3], u[uOff + 4], u[uOff + 5]);
    u[uOff + 0] = wR[0];
    u[uOff + 1] = wR[1];
    u[uOff + 2] = wR[2];
    u[uOff + 3] = vR[0];
    u[uOff + 4] = vR[1];
    u[uOff + 5] = vR[2];

    // re-realize geometry and re-seed the derivative chain for the next Verlet step
    RobotEngine::realizePosition(model_, state_);
    RobotEngine::fillAtomPositionsFromBodies(model_, state_);
    bridge_.evaluate(state_);
    RobotEngine::realizeVelocity(model_, state_);
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);
    RobotEngine::calcUDot(model_, state_);
    RobotEngine::calcQDot(model_, state_, state_.qdot());
    RobotEngine::calcQDotDot(model_, state_);
}

bool World::ncmcInnerGhmcStep(robo::Real h, int substepIndex, bool* acceptedOut) {
    // Construction II inner kernel (docs/specs/ncmc-explicit-solvent/
    // 10-acceptance-construction.md Sec.3; 20-inner-integrator.md Sec.2 F3).
    //
    // One fixed-lambda Verlet substep, proposed and then Metropolis accept/
    // reject'd against the FULL H_lambda = V_lambda + K(u;q) + K_s(v_s) + U_F(q)
    // - 1/2 RT ln sin^2(gamma2) - nmaCorr -- reusing currentTotalEnergy()'s own
    // assembly VERBATIM (F1/F2: a shortcut that dropped U_F/pitch here would make
    // the inner kernel non-pi_lambda-invariant and silently sample the WRONG
    // density; INV0 in tests/TestNcmcExplicitSolvent.cpp is the discriminating
    // oracle). On reject: restore (q,u,x_s,v_s) to their pre-substep values and
    // NEGATE (u,v_s) -- the standard GHMC reject-flip that keeps this kernel
    // F-reversible and hence exactly pi_lambda-invariant.
    //
    // F3: a velocity corrector that does not converge at this dt is ALSO an
    // automatic reject (never silently taken as Construction I's propagator
    // does) -- under Construction II a non-converged corrector need not be
    // F-reversible, and taking it as the GHMC proposal would break
    // pi_lambda-invariance with no visible symptom in the outer acceptance
    // (10-...:Sec.5 NOTE F3, 20-...:Sec.2 second bullet).
    if (acceptedOut) {
        *acceptedOut = false; // pessimistic default; set true only on a genuine accept below
    }
    const std::vector<int>& solv = state_.cartSolventAtoms();
    robo::Vec3* posG = state_.atomPosG();
    robo::Vec3* velG = state_.atomVelG();
    std::vector<Real> q0(state_.q(), state_.q() + model_.nq);
    std::vector<Real> u0(state_.u(), state_.u() + model_.nu);
    std::vector<robo::Vec3> xs0(solv.size()), vs0(solv.size());
    for (std::size_t j = 0; j < solv.size(); ++j) {
        xs0[j] = posG[solv[j]];
        vs0[j] = velG[solv[j]];
    }

    const double Hbefore = currentTotalEnergy(); // H_lambda at the CURRENT bridge lambda;
    // also refreshes bodyForceG/mobilityForce (bridge_.evaluate) and V_GB/qdot
    // (realizeVelocity) at (q0, the CURRENT lambda).

    // BUG FIX (state/energy-assembly inconsistency, freeze root cause): the
    // PERTURB substep in ncmcMove (lambda change) refreshes FORCES via
    // bridge_.evaluate but never recomputes udot/qdotdot -- it doesn't need to
    // for Construction I, whose stepTo call seeds that chain once at the top
    // of the move and never revisits it mid-substep. Construction II calls
    // THIS function once per substep, immediately after the SAME perturb
    // block may have just changed lambda, so without a local reseed here
    // `state_.udot()`/`qdotdot()` at verletStep's entry (the `a0` term of its
    // velocity-Verlet trapezoid) are evaluated at the STALE, PRE-perturb
    // lambda's forces -- a genuine assembly inconsistency baked directly into
    // the proposal itself (dt-INDEPENDENT: it does not shrink at small h the
    // way ordinary shadow work does, since it is a wrong-physics `a0`, not a
    // discretization-error `a0`). That corrupts the GHMC proposal and, given
    // alchemical decoupling can change intermolecular forces by orders of
    // magnitude at low lambda, is large enough to make every dH huge
    // regardless of dt -- exactly the observed 100x-dt-invariant freeze.
    // realizeArticulatedBodyInertias is called again UNCONDITIONALLY (cheap,
    // idempotent) rather than relying on calcFixman's internal call above, so
    // this reseed is correct even with useFixman=false. Scoped to
    // ncmcInnerGhmcStep (Construction-II-only code, never called by
    // Construction I), so Construction I's shared perturb block and stepTo
    // call are untouched -- bit-for-bit preserved.
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);
    RobotEngine::calcUDot(model_, state_);
    RobotEngine::calcQDot(model_, state_, state_.qdot());
    RobotEngine::calcQDotDot(model_, state_);

    bool converged = true;
    const bool stepOk =
        RobotEngine::stepTo(model_, state_, bridge_, constraints_, state_.time + h, &converged);
    if (!stepOk) {
        // Non-finite force/velocity: verletStep already restored its own
        // pre-step state internally. Mirror stepTo's contract exactly -- the
        // caller (ncmcMove) treats false as a hard abort of the whole move,
        // orthogonal to the F3 corrector-convergence guard below.
        if (ncmcDebugEnabled()) {
            std::fprintf(
                stderr,
                "[ncmc-inner] world %d substep %d: h=%.6g stepTo FAILED (non-finite) -> abort move\n",
                index_,
                substepIndex,
                (double)h);
        }
        return false;
    }

    bool accept = false;
    double Hafter = std::numeric_limits<double>::quiet_NaN();
    double dH = std::numeric_limits<double>::quiet_NaN();
    if (converged) {
        Hafter = currentTotalEnergy(); // same lambda; H_lambda at the proposed state
        dH = Hafter - Hbefore;
        accept = std::isfinite(dH) && (dH <= 0.0 || uniform_(rng_) < std::exp(-beta_ * dH));
    }
    // else: corrector did not converge at this dt -- F3 forces `accept = false`
    // rather than silently taking a proposal that need not be F-reversible.

    if (ncmcDebugEnabled()) {
        std::fprintf(stderr,
                     "[ncmc-inner] world %d substep %d: h=%.6g Hbefore=%.6f Hafter=%.6f dH=%+.6f "
                     "converged=%s accept=%s\n",
                     index_,
                     substepIndex,
                     (double)h,
                     Hbefore,
                     Hafter,
                     dH,
                     converged ? "true" : "false",
                     accept ? "true" : "false");
    }

    if (!accept) {
        std::copy(q0.begin(), q0.end(), state_.q());
        std::copy(u0.begin(), u0.end(), state_.u());
        for (std::size_t j = 0; j < solv.size(); ++j) {
            posG[solv[j]] = xs0[j];
            velG[solv[j]] = vs0[j];
        }
        // Reject-flip: negate the persistent momenta (GHMC).
        Real* u = state_.u();
        for (int i = 0; i < model_.nu; ++i) {
            u[i] = -u[i];
        }
        for (std::size_t j = 0; j < solv.size(); ++j) {
            robo::Vec3& v = velG[solv[j]];
            v = robo::Vec3(-v[0], -v[1], -v[2]);
        }
        // Re-seed the derivative chain at the restored (q0,-u0,x_s0,-v_s0) state,
        // mirroring ncmcApplyTroughTeleport's post-mutation reseed.
        RobotEngine::realizePosition(model_, state_);
        RobotEngine::fillAtomPositionsFromBodies(model_, state_);
        bridge_.evaluate(state_);
        RobotEngine::realizeVelocity(model_, state_);
        RobotEngine::realizeArticulatedBodyInertias(model_, state_);
        RobotEngine::calcUDot(model_, state_);
        RobotEngine::calcQDot(model_, state_, state_.qdot());
        RobotEngine::calcQDotDot(model_, state_);
    }
    if (acceptedOut) {
        *acceptedOut = accept;
    }
    return true;
}

bool World::ncmcMove() {
    // Per-molecule NCMC (Nilmeier, Crooks, Minh & Chodera 2011): a lambda:1->0->1
    // alchemical decouple-move-recouple proposal. During the switch the
    // INTERMOLECULAR nonbonded between this world's molecule and every other
    // molecule is softened (its intramolecular physics is untouched), so the
    // molecule strides free of the cage of contacting molecules near lambda=0 and
    // is recoupled as its mobile DOF relax to fit the (frozen) environment. The
    // environment itself relaxes in the Cartesian world of the Gibbs scan (DOF
    // coverage, THEORY 13.5); a co-mobilized local shell would let it relax inside
    // the move too (documented efficiency follow-up).
    savedQ_.assign(state_.q(), state_.q() + model_.nq);
    // Save the Cartesian solvent positions for a clean rollback on reject: the
    // body->atom fill skips these atoms, so restoring q alone would leave the
    // solvent at its (rejected) end pose. Empty/no-op off solvent-relaxing NCMC.
    {
        const std::vector<int>& solv = state_.cartSolventAtoms();
        const robo::Vec3* posG = state_.atomPosG();
        savedSolventPosG_.resize(solv.size());
        for (std::size_t j = 0; j < solv.size(); ++j) {
            savedSolventPosG_[j] = posG[solv[j]];
        }
    }
    bridge_.setAlchemicalLambda(1.0);
    reinitialize(); // draws p ~ N(0, RT M(q0)); sets Hold_ = V1 + K + F + J
    const double Hstart = Hold_;
    // Diagnostic-only (docs/specs/ncmc-explicit-solvent/40-reproducer-and-oracles.md
    // Sec.2/Sec.6): the Fixman/pitch state-function values at the START, so the
    // reproducer can compute dU_F = fixman_end - fixman_start and dJ from the log
    // without re-deriving them. Zero unless useFixman/useOrientationJacobian.
    const double fixmanStart = state_.energy.fixman;
    const double logSineSqrStart = state_.energy.logSineSqrGamma2;

    // No explicit momentum flip. Momenta are resampled from Maxwell-Boltzmann at
    // the start of every block (reinitialize, above), itself a Gibbs move on the
    // velocity marginal. Under that resampling the NCMC momentum-reversal is
    // unnecessary -- Nilmeier et al. 2011 explicitly sanction reinitializing
    // velocities after each NCMC step -- since the carried-over sign is overwritten
    // next block and never read. (Mirrors the existing MdHmc move, which also
    // resamples and accepts on dH without an explicit flip.)

    const robo::Real h = sampler_.timeStep;
    // DIAGNOSTIC ONLY (NOT used in acceptance): protocol work w = sum over the
    // PERTURBATION substeps of dV at fixed q (Nilmeier et al. 2011, Eq. 16).
    // NOTE: w is NOT equal to Hend - Hstart at finite dt. Two distinct effects open
    // the gap, and NEITHER biases acceptance (acceptance is on Hend - Hstart):
    //   (i)  integrator heat Q != 0 -- the fixed-lambda Verlet steps are only
    //        near-symplectic, so each pumps a little shadow energy that is in
    //        Hend - Hstart but never in w (w only sees the fixed-q perturbations);
    //   (ii) the Fixman + orientation-Jacobian state-function drift F(qend)-F(q0)
    //        and J(qend)-J(q0), which ARE in Hend - Hstart (currentTotalEnergy) but
    //        are excluded from w by construction.
    // So logging gap = w - (Hend - Hstart) is a free consistency check only: a
    // large gap flags a too-large dt / non-converged corrector, not a sampling bug.
    double work = 0.0;
    bool ok = true;
    double Vprev = bridge_.calcPotentialEnergy(); // V at lambda=1, current q
    double lamPrev = 1.0;
    const double peStart = state_.energy.pe; // for the per-move dPE/dKE breakdown
    const double keStart = state_.energy.ke;
    const double keSolvStart = state_.energy.keSolvent; // solvent Cartesian KE at start

    // λ=0 trough teleport (default off; tests/TestNcmcTeleport). Centered in the
    // λ=0 hold so the 1->0->1 protocol stays its own reverse; applied only for an
    // ACYCLIC Free-root region with a defined target region (docking site), where
    // there is no constraint-manifold branch ambiguity (CHMC Thm 3) to guard.
    const int teleHold = static_cast<int>(sampler_.ncmcHoldFraction * sampler_.ncmcSteps);
    const int teleRamp = std::max((sampler_.ncmcSteps - teleHold) / 2, 1);
    const int teleStep = teleRamp + teleHold / 2; // center of the hold
    const int teleRoot = ncmcTeleportRoot();
    const bool doTele =
        sampler_.ncmcTeleport && teleHold > 0 && teleRoot > 0 && !siteAtoms_.empty() && constraints_.empty();
    if (sampler_.ncmcTeleport && !doTele) {
        std::fprintf(stderr,
                     "[ncmc] teleport requested but inactive (needs hold>0, a Free-root region, a "
                     "docking site, and an acyclic system) -- using the plain uncaged stride\n");
    }

    // Per-move inner-GHMC acceptance count (Construction II only; stays 0/0 for
    // Construction I). Makes inner acceptance visible in the [ncmc] summary
    // line without inferring it from a frozen PE (the symptom that made the
    // Construction-II freeze a black box).
    int innerAccepted = 0;
    int innerAttempted = 0;

    for (int s = 0; ok && s < sampler_.ncmcSteps; ++s) {
        // (i) PERTURB: change lambda at FIXED q; accumulate work = V(lam_new) - V(lam_old).
        const double lam = protocolLambda(s);
        if (lam != lamPrev) {
            bridge_.setAlchemicalLambda(lam);
            bridge_.evaluate(state_); // recompute at same q, new lambda
            const double Vnew = bridge_.calcPotentialEnergy();
            if (!std::isfinite(Vnew)) {
                ok = false;
                break;
            }
            work += Vnew - Vprev;
            Vprev = Vnew;
            lamPrev = lam;
        }
        // (i.5) TELEPORT at the λ=0 trough (free: ghost, KE preserved by co-rotation).
        if (doTele && s == teleStep) {
            ncmcApplyTroughTeleport(teleRoot);
            Vprev = bridge_.calcPotentialEnergy(); // λ=0 => unchanged; refresh for the heat bookkeeping
        }
        // (ii) PROPAGATE one Verlet step at fixed lambda.
        if (sampler_.useMetropolizedInner) {
            // Construction II (10-acceptance-construction.md Sec.3): the substep
            // is an inner GHMC kernel, Metropolized against the FULL H_lambda, so
            // its shadow work is absorbed into inner rejections and never reaches
            // the outer acceptance. ncmcInnerGhmcStep returns false ONLY on a
            // genuinely unrecoverable non-finite condition (mirrors stepTo).
            bool innerAccept = false;
            ok = ncmcInnerGhmcStep(h, s, &innerAccept);
            ++innerAttempted;
            innerAccepted += innerAccept ? 1 : 0;
        } else {
            // Construction I (endpoint-DeltaH, unchanged): deterministic,
            // unadjusted, reversible Verlet.
            ok = RobotEngine::stepTo(model_, state_, bridge_, constraints_, state_.time + h);
        }
        if (ok) {
            Vprev = bridge_.calcPotentialEnergy(); // fixed-lambda V drift = heat, not work
        }
    }

    bool finite = ok;
    if (finite) {
        const robo::Real* q = state_.q();
        for (int i = 0; i < model_.nq; ++i) {
            if (!std::isfinite(q[i])) {
                finite = false;
                break;
            }
        }
    }

    bool accepted = false;
    if (finite) {
        bridge_.setAlchemicalLambda(1.0);
        const double Hend = currentTotalEnergy(); // V1(qend) + K + F + J at lambda=1
        const double fixmanEnd = state_.energy.fixman;
        const double logSineSqrEnd = state_.energy.logSineSqrGamma2;

        // ACCEPTANCE. Two exact constructions (docs/specs/ncmc-explicit-solvent/
        // 10-acceptance-construction.md); SHALL NOT mix them (Sec.2 CLAIM C1):
        //
        // Construction I (endpoint-DeltaH, sampler_.useMetropolizedInner == false,
        // unchanged): this NCMC trajectory is a valid HMC proposal -- a
        // DETERMINISTIC, volume-preserving map T on (q,p) that is momentum-flip
        // reversible, F T F == T^-1. T is the composition of the per-substep
        // fixed-lambda Verlet steps (the fixed-q lambda perturbations do not move
        // the state); each step is F-reversible, so the composition is reversible
        // BECAUSE the lambda schedule is a PALINDROME pinned to lambda = 1 at both
        // endpoints (protocolLambda / NCMCProtocol.hpp). We therefore accept on the
        // FULL Hamiltonian difference at the lambda=1 endpoints (Nilmeier et al.
        // 2011, Eq. 20; bistable-dimer Eq. 28), NOT on the work. GATE: lambda == 1
        // throughout => work == 0 and Hend - Hstart is the plain Verlet dH =>
        // reduces EXACTLY to the torsional-HMC metropolis test.
        //
        // Construction II (Metropolized-dynamics NCMC, useMetropolizedInner ==
        // true): every fixed-lambda substep was ALREADY Metropolized against the
        // full H_lambda inside the loop above (ncmcInnerGhmcStep), so shadow work
        // never reaches this acceptance -- accepting again on Hend - Hstart here
        // would double-count it and also re-admit the very bath shadow work the
        // construction exists to remove (10-...:Sec.3). The OUTER move instead
        // accepts on the protocol work ALONE: a = min(1, exp(-beta*W)) (Sec.3).
        // metropolis(0.0, work) computes exactly that (dH = work - 0 = work).
        const double dH = Hend - Hstart;
        // Per-move diagnostic breakdown (Phase 0): the gap = work - dH flags
        // integrator energy pumping (propagator: dt / mass-scale) plus the Fixman/
        // Jacobian drift, while the recoupling potential change dPE isolates the
        // irreducible reorganization (insertion) cost. dKE should stay ~0. This
        // breakdown is diagnostic under BOTH constructions (gap is never used in
        // Construction II's acceptance either -- only `work` is).
        const double dPE = state_.energy.pe - peStart;
        const double dKE = state_.energy.ke - keStart;
        const double dKEsolv = state_.energy.keSolvent - keSolvStart;
        std::fprintf(stderr,
                     "[ncmc] world %d: construction=%s Hstart=%.2f Hend=%.2f dH=%+.2f kJ/mol  "
                     "dPE=%+.2f dKE=%+.2f dKEsolv=%+.2f work(diag)=%.2f gap=%+.2f "
                     "fixman_start=%.4f fixman_end=%.4f logSineSqr_start=%.4f logSineSqr_end=%.4f "
                     "inner_accepted=%d/%d (steps=%d, teleport=%s)\n",
                     index_,
                     sampler_.useMetropolizedInner ? "II" : "I",
                     Hstart,
                     Hend,
                     dH,
                     dPE,
                     dKE,
                     dKEsolv,
                     work,
                     work - dH,
                     fixmanStart,
                     fixmanEnd,
                     logSineSqrStart,
                     logSineSqrEnd,
                     innerAccepted,
                     innerAttempted,
                     sampler_.ncmcSteps,
                     doTele ? "on" : "off");
        const bool moveAccept = sampler_.useMetropolizedInner
                                    ? (std::isfinite(work) && metropolis(0.0, work))
                                    : (std::isfinite(dH) && metropolis(Hstart, Hend));
        if (moveAccept) {
            RobotEngine::fillAtomPositionsFromBodies(model_, state_);
            accepted = true;
        }
    } else {
        std::fprintf(stderr, "[ncmc] non-finite during protocol -> reject\n");
    }

    if (!accepted) {
        bridge_.setAlchemicalLambda(1.0);
        std::copy(savedQ_.begin(), savedQ_.end(), state_.q());
        RobotEngine::realizePosition(model_, state_);
        RobotEngine::fillAtomPositionsFromBodies(model_, state_);
        // Restore the Cartesian solvent pose (the fill above skips these atoms).
        const std::vector<int>& solv = state_.cartSolventAtoms();
        robo::Vec3* posG = state_.atomPosG();
        for (std::size_t j = 0; j < solv.size(); ++j) {
            posG[solv[j]] = savedSolventPosG_[j];
        }
    }
    lastAccepted_ = accepted;
    return accepted;
}
