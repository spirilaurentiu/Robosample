#pragma once
// ============================================================================
//  SamplingHarness.hpp -- the two sampling loops every ensemble/NCMC test
//  file re-inlined (TESTS.md section 4): `for (i<N) drv.move(); observe...` for
//  the HMC ensemble suites, and the NCMC PERTURB/PROPAGATE work-chain loop
//  (`NcmcRun`/`runNcmcLoop`) for the NCMC work-identity suites. This header
//  lifts both so a converted test differs from its sibling only in the
//  assertion, per TESTS.md section 5.5.
//
//  The harness owns no engine state and no RNG: it drives the caller's
//  HmcDriver/bridge, which already carry the seeded streams (HmcDriver.hpp:76).
//  It never reseeds, so a converted call site draws the exact same sample
//  stream as its pre-harness inline loop.
// ============================================================================

#include "../HmcDriver.hpp"
#include "../RobotBuilders.hpp" // RobotModel/RobotState/robot_math via its own includes
#include "Constraints.hpp"
#include "NCMCProtocol.hpp"
#include "RobotEngine.hpp"
#include "RobotIntegrator.hpp" // stepTo

namespace rtest {

// ---------------------------------------------------------------------------
//  runHmcChain -- the shared HMC ensemble loop. Runs `nMoves` HmcDriver::move()
//  calls IN ORDER (no reseed, no skipped/extra moves), invoking
//  `observe(i, accepted)` every `stride` moves (i % stride == 0; stride=1
//  observes every move, matching the ensemble suites that never thinned).
//  `observe` pushes whatever scalar/vector observable the caller's own
//  accumulator/histogram(s) need -- runHmcChain does not know their shape, so
//  it takes the caller's accumulator only via the observer closure. The
//  returned Marginals mirrors HmcDriver's own accept/attempt counters
//  (HmcDriver.hpp:84-90) at the end of the chain; the harness does not
//  duplicate that bookkeeping.
// ---------------------------------------------------------------------------
struct Marginals {
    long attempted = 0;
    long accepted = 0;
};

template <class Bridge, class Observer>
Marginals runHmcChain(HmcDriver<Bridge>& drv, long nMoves, int stride, Observer&& observe) {
    for (long i = 0; i < nMoves; ++i) {
        const bool acc = drv.move();
        if (i % stride == 0) {
            observe(i, acc);
        }
    }
    Marginals out;
    out.attempted = drv.attempted();
    out.accepted = drv.accepted();
    return out;
}

// ---------------------------------------------------------------------------
//  The NCMC work-chain loop (moved verbatim from TestNCMCWork.cpp:236-336,
//  templatized over Bridge -- the one shape change TEST-004 makes, so
//  TestNcmcExplicitSolvent.cpp's HarmonicBridge<IntraLambdaInterPolicy>
//  construction path can share it; TestNCMCWork.cpp's own instantiation is
//  unchanged).
// ---------------------------------------------------------------------------

// Seed the derivative chain so the first verletStep has valid
// qdot0/udot0/qdd0 (identical preamble to TestIntegrator.cpp's
// seedDerivatives).
template <class Bridge>
void seedDerivatives(const RobotModel& m, RobotState& s, Bridge& bridge) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::fillAtomPositionsFromBodies(m, s);
    bridge.evaluate(s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
    RobotEngine::calcUDot(m, s);
    RobotEngine::calcQDot(m, s, s.qdot());
    RobotEngine::calcQDotDot(m, s);
}

// Generalized-coordinate kinetic energy 1/2 u^T M(q) u (the same
// calcKineticEnergy World uses in its Hamiltonian).
inline robo::Real kineticEnergy(const RobotModel& m, RobotState& s) {
    RobotEngine::realizeVelocity(m, s);
    return RobotEngine::calcKineticEnergy(m, s);
}

// Result of running the ncmcMove loop over a model potential: the telescope
// identity this preserves is Hend - Hstart == work + heat (TestNCMCWork.cpp's
// AccumulatorIsTheFixedQTelescope / WorkPlusHeatEqualsDeltaH), up to the
// integrator's shadow work captured in `heat`.
struct NcmcRun {
    double work = 0;   // the protocol-work accumulator (diagnostic in production)
    double Hstart = 0; // PE + KE at lambda=1, start
    double Hend = 0;   // PE + KE at lambda=1, end
    double heat = 0;   // sum of fixed-lambda propagate dH (the integrator's shadow work)
    bool ok = true;
};

// Reproduce ncmcMove's PERTURB/PROPAGATE loop EXACTLY (minus the Metropolis
// gate and momentum resampling, which are orthogonal to the work identity),
// over a lambda-aware analytic Bridge. H is PE+KE here (Fixman is constant in
// q-free perturb). nLambda controls how many distinct lambda values the
// schedule visits.
template <class Bridge>
NcmcRun runNcmcLoop(RobotModel& m,
                    RobotState& s,
                    Bridge& bridge,
                    robo::ConstraintSet& cs,
                    int ncmcSteps,
                    double holdFraction,
                    double h) {
    NcmcRun r;
    bridge.setLambda(1.0);
    seedDerivatives(m, s, bridge);
    r.Hstart = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);

    double work = 0.0, heat = 0.0;
    double Vprev = bridge.calcPotentialEnergy(s);
    double lamPrev = 1.0;

    for (int step = 0; r.ok && step < ncmcSteps; ++step) {
        // (i) PERTURB at FIXED q: change lambda, accumulate work = V_new - V_old.
        const double lam = robo::ncmc::protocolLambda(step, ncmcSteps, holdFraction);
        if (lam != lamPrev) {
            bridge.setLambda(static_cast<robo::Real>(lam));
            const double Vnew = bridge.calcPotentialEnergy(s); // SAME q, new lambda
            work += Vnew - Vprev;
            Vprev = Vnew;
            lamPrev = lam;
        }
        // (ii) PROPAGATE one Verlet step at fixed lambda; the fixed-lambda dH is heat.
        const double Hbefore = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);
        r.ok = RobotEngine::stepTo(m, s, bridge, cs, s.time + h);
        if (!r.ok) {
            break;
        }
        const double Hafter = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);
        heat += Hafter - Hbefore;
        Vprev = bridge.calcPotentialEnergy(s);
    }

    // Hend is evaluated at the protocol's ACTUAL final lambda: the work
    // accumulator only saw the lambda changes that fired, so the energy
    // telescope H_end - H_start == work + heat holds exactly only when Hend
    // uses the same final lambda the loop left the bridge in. With the
    // PALINDROMIC schedule the last index s = N-1 is pinned to lambda = 1
    // (both endpoints are), so this final lambda is 1 -- the same lambda
    // production's ncmcMove accepts on.
    r.Hend = bridge.calcPotentialEnergy(s) + kineticEnergy(m, s);
    r.work = work;
    r.heat = heat;
    return r;
}

} // namespace rtest
