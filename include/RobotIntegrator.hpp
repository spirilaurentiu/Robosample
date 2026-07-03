#pragma once
// ============================================================================
//  RobotIntegrator.hpp -- the fixed-step velocity-Verlet driver, as FUNCTION
//  TEMPLATES on the force-bridge type.
//
//  WHY A TEMPLATE. verletStep/stepTo/checkReversibility only ever touch the
//  bridge through `bridge.evaluate(s)`. Templating on the bridge type lets the
//  unit tests drive the REAL integrator with a test-only, OpenMM-free bridge
//  (tests/AnalyticForceBridge.hpp) -- no faking, no std::function indirection,
//  zero runtime cost in production (the only instantiation there is the real
//  ForceBridge). This is the minimal change that unlocks every integrator /
//  energy-conservation / reversibility / ensemble test without OpenMM.
//
//  These definitions were moved verbatim out of src/RobotEngine.cpp (only
//  `ForceBridge& bridge` -> `Bridge& bridge`); advanceQuatExp was a file-local
//  static there and moves here with them (its sole caller is verletStep). Note
//  RobotEngine / RobotModel / JointType are in the GLOBAL namespace (only the
//  math types live in `robo`), so these definitions and the using-block below
//  mirror the file-scope environment RobotEngine.cpp provided.
//
//  LEAF HEADER: include this from the TUs that actually call the integrator
//  (src/World.cpp and the integrator tests). It is intentionally NOT included by
//  RobotEngine.hpp, to avoid the RobotEngine.hpp <-> Constraints.hpp cycle
//  (Constraints.hpp includes RobotEngine.hpp; this header needs the COMPLETE
//  ConstraintSet, so it includes Constraints.hpp itself).
//
//  ROBO_CHECK: the optional per-step NaN-scan macro is a development aid defined
//  inside src/RobotEngine.cpp. The always-on correctness guards (forcesFinite /
//  velocitiesSane local lambdas, and the corrector-convergence throw) live in
//  the bodies below and are unaffected. When this header is included by a TU
//  that has NOT defined ROBO_CHECK (every caller), it is a no-op, which keeps all
//  instantiations of a given Bridge identical across TUs (no ODR surprise).
// ============================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "Constraints.hpp"
#include "RobotEngine.hpp"
#include "robot_math.hpp"

#ifndef ROBO_CHECK
#    define ROBO_CHECK(where) ((void)0)
#    define ROBO_INTEGRATOR_LOCAL_ROBO_CHECK
#endif

// The integrator bodies were moved from RobotEngine.cpp's global scope; mirror
// the math-type using-block it relied on (RobotEngine/RobotModel/JointType are
// already global).
using robo::ArticulatedInertia;
using robo::Mat33;
using robo::PhiMatrix;
using robo::Quaternion;
using robo::Real;
using robo::Rotation;
using robo::SpatialInertia;
using robo::SpatialVec;
using robo::SymMat33;
using robo::Transform;
using robo::Vec3;
using robo::Vec4;

// exact unit-quaternion advance (exponential map); sole caller is verletStep.
// -------- exact unit-quaternion advance (exponential map) ------------------
// Advance a unit quaternion under a constant angular velocity w_F expressed in
// the PARENT (F) frame, consistent with the engine's qdot = 1/2 (0,w_F) (x) q
// (left multiply; matches convertAngVelToQuaternionDot / N(q)). With theta =
// 1/2 |w| h: Dq = (cos theta, sin theta * w_hat) and q1 = Dq (x) q0. Both
// operands are unit, so q1 is unit BY CONSTRUCTION -- the |q|^2 -> Inf ->
// q/sqrt(Inf) = 0 overflow path of the linear "q += h*qdot" drift cannot occur.
// Reversible: negating w gives Dq^-1, so q0 = Dq^-1 (x) q1 (HMC needs this).
// As h -> 0 it reduces to q0 + h * (1/2 (0,w) (x) q0) = q0 + h * qdot0.
inline void advanceQuatExp(const Real* q0, const Vec3& wF, Real h, Real* q1) {
    const Real wn = std::sqrt(wF[0] * wF[0] + wF[1] * wF[1] + wF[2] * wF[2]);
    Real a, b, c, d; // Dq = (a, b, c, d)
    if (wn > Real(1e-12)) {
        const Real theta = Real(0.5) * wn * h;
        const Real ssc = std::sin(theta) / wn; // sin(theta) / |w|
        a = std::cos(theta);
        b = ssc * wF[0];
        c = ssc * wF[1];
        d = ssc * wF[2];
    } else {
        a = Real(1); // small-angle limit: Dq ~ (1, 1/2 h w)
        b = Real(0.5) * h * wF[0];
        c = Real(0.5) * h * wF[1];
        d = Real(0.5) * h * wF[2];
    }
    const Real w = q0[0], x = q0[1], y = q0[2], z = q0[3];
    // Hamilton product q1 = Dq (x) q0  (LEFT multiply), consistent with the
    // parent-frame map qdot = 1/2 (0,w) (x) q. (C-3 fix: the body previously
    // computed the right product q0 (x) Dq, which pairs with the body-frame map.)
    q1[0] = a * w - b * x - c * y - d * z;
    q1[1] = a * x + b * w + c * z - d * y;
    q1[2] = a * y + c * w + d * x - b * z;
    q1[3] = a * z + d * w + b * y - c * x;
}

template <class Bridge>
bool RobotEngine::verletStep(const RobotModel& m,
                             RobotState& s,
                             Bridge& bridge,
                             const robo::ConstraintSet& cset,
                             Real h,
                             bool* correctorConverged) {
    if (correctorConverged) {
        *correctorConverged = false; // pessimistic default; set true only on a converged step below
    }
    const int nq = m.nq, nu = m.nu;
    Real* q = s.q();
    Real* u = s.u();
    Real* qdot = s.qdot();
    Real* udot = s.udot();
    Real* qdd = s.qdotdot();

    std::vector<Real> q0(q, q + nq);
    std::vector<Real> u0(u, u + nu);
    std::vector<Real> qdot0(qdot, qdot + nq);
    std::vector<Real> udot0(udot, udot + nu);
    std::vector<Real> qdd0(qdd, qdd + nq);

    // ---- Cartesian-integrated solvent (solvent-relaxing NCMC) ----------------
    // A subset of atoms is advanced in FLAT Cartesian space by velocity-Verlet
    // driven by the same OpenMM force evaluation that drives the solute internal
    // step -- so the contact environment relaxes INSIDE the proposal instead of
    // being a welded wall (docs/specs/ncmc_solvent_relax.md). This shares the
    // existing integrator's structure exactly: position drift x1 = x0 + h v0 +
    // (h^2/2) a0, then the velocity half is the same implicit trapezoid the
    // generalized speeds use -- but the Cartesian force depends only on position
    // (no velocity-dependent terms), so that trapezoid is EXPLICIT and needs no
    // iteration: v1 = v0 + (h/2)(a0 + a1). The map is the standard symmetric,
    // time-reversible, volume-preserving Verlet, so the move's exact-dH HMC
    // acceptance argument is preserved (checkReversibility certifies it). The
    // block is skipped entirely when no atom is flagged (the welded engine),
    // making every existing path bit-identical.
    const std::vector<int>& solvAtoms = s.cartSolventAtoms();
    const std::vector<Real>& solvInvM = s.cartSolventInvMass();
    const int nSolv = static_cast<int>(solvAtoms.size());
    Vec3* posG = s.atomPosG();
    Vec3* velG = s.atomVelG();
    Vec3* frcG = s.atomForceG();
    std::vector<Vec3> xs0(nSolv), vs0(nSolv), fs0(nSolv);
    for (int j = 0; j < nSolv; ++j) {
        const int a = solvAtoms[j];
        xs0[j] = posG[a];  // x0
        vs0[j] = velG[a];  // v0
        fs0[j] = frcG[a];  // a0 = f0 * invMass (cached by the previous force eval)
    }

    // True iff every per-body spatial force is finite. OpenMM returns Inf forces
    // for a hard steric overlap (LJ r^-12); integrating against those is what
    // produced the "Particle coordinate is NaN" crash. Catch it here, at the
    // source, and signal the caller to reject this step's pose.
    auto forcesFinite = [&]() -> bool {
        const SpatialVec* bf = s.bodyForceG();
        for (int b = 1; b < m.numBodies; ++b) {
            if (!std::isfinite(bf[b][0][0]) || !std::isfinite(bf[b][0][1]) || !std::isfinite(bf[b][0][2])
                || !std::isfinite(bf[b][1][0]) || !std::isfinite(bf[b][1][1])
                || !std::isfinite(bf[b][1][2])) {
                return false;
            }
        }
        return true;
    };

    // forcesFinite() only catches OpenMM Inf. A CFL-unstable step (e.g. a stiff
    // explicit-solvent contact at too-large dt) amplifies u geometrically while
    // it is still finite, then overflows the quaternion. Cap ||u|| relative to
    // the freshly-seeded momentum norm so a runaway is rejected early, while a
    // merely hot trajectory (||u|| growing a few x) passes.
    Real uSeedNorm2 = 0;
    for (int i = 0; i < nu; ++i) {
        uSeedNorm2 += u0[i] * u0[i];
    }
    const Real uCap2 = (uSeedNorm2 + Real(1e-30)) * Real(1e6);
    auto velocitiesSane = [&]() -> bool {
        Real n2 = 0;
        for (int i = 0; i < nu; ++i) {
            if (!std::isfinite(u[i])) {
                return false;
            }
            n2 += u[i] * u[i];
        }
        return n2 <= uCap2;
    };

    auto restorePreStep = [&]() {
        std::copy(q0.begin(), q0.end(), q);
        std::copy(u0.begin(), u0.end(), u);
        std::copy(qdot0.begin(), qdot0.end(), qdot);
        std::copy(udot0.begin(), udot0.end(), udot);
        std::copy(qdd0.begin(), qdd0.end(), qdd);
        realizePosition(m, s);
        fillAtomPositionsFromBodies(m, s);
        // Cartesian solvent is owned by posG/velG, not by any body, so it is not
        // rolled back by the fill above -- restore it explicitly.
        for (int j = 0; j < nSolv; ++j) {
            const int a = solvAtoms[j];
            posG[a] = xs0[j];
            velG[a] = vs0[j];
            frcG[a] = fs0[j];
        }
    };

    // ---- position drift ----
    // Scalar DOFs (Torsion, Translation, and the Free TRANSLATION block): 2nd-order
    // Taylor q1 = q0 + h*qdot0 + (h^2/2)*qddot0 (local error O(h^3); global method
    // order 2 -- see THEORY 5.2).
    //
    // Quaternion DOFs (Free, Ball): EXACT exponential-map advance from the
    // midpoint angular velocity wHalf = u0 + (h/2)*udot0, applied as a left
    // Hamilton product onto q0. This DELIBERATELY diverges from Simbody's
    // linear-Taylor-plus-renormalize quaternion update, and is NOT a bug to be
    // "matched away": it keeps |q| = 1 by construction, is reversible, reduces to
    // the linear update as h -> 0 (so it changes only the proposal, not the target
    // -- THEORY 5.4), and -- crucially -- it advances the orientation purely from
    // the angular velocity, BYPASSING the quaternion second derivative qddot.
    // The linear-Taylor update instead leans on qddot (calcQDotDot); with the
    // port's reimplemented Free-joint kinematics that path injects kinetic energy
    // through the root body and (via the articulated recursion) the whole tree --
    // an observed ~+1300 kJ/mol/traj KE pump even at 1 fs. The exp-map is robust to
    // that latent N/qddot inconsistency. ROOT-CAUSE TODO: golden-test the Free-joint
    // quaternion kinematics (N, Ndot, qddot) against Simbody on a single free body
    // (THEORY 3.4); until that is closed, the exp-map is the correct propagator.
    for (int i = 0; i < nq; ++i) {
        q[i] = q0[i] + (h * qdot0[i]) + ((h * h / 2) * qdd0[i]);
    }

    for (int bodyIx = 1; bodyIx < m.numBodies; ++bodyIx) {
        if (!m.isQuaternionBody(bodyIx)) {
            continue;
        }
        const int qOff = m.bodyQIndex[bodyIx];
        const int uOff = m.bodyUIndex[bodyIx];
        // Half-step angular velocity in F. Ball/Free: the first 3 u ARE w_FM.
        // FreeLine: the 2 rotational speeds are (x,y) of w_FM expressed in M, so
        // w_FM = R_FM*(u0,u1,0) (R_FM from the step-start X_FM, still valid here).
        Vec3 wHalf;
        if (m.bodyJoint[bodyIx] == JointType::FreeLine) {
            const Rotation& R_FM = s.X_FM()[bodyIx].R();
            wHalf = R_FM
                    * Vec3(u0[uOff + 0] + (0.5 * h * udot0[uOff + 0]),
                           u0[uOff + 1] + (0.5 * h * udot0[uOff + 1]),
                           0);
        } else {
            wHalf = Vec3(u0[uOff + 0] + (0.5 * h * udot0[uOff + 0]),
                         u0[uOff + 1] + (0.5 * h * udot0[uOff + 1]),
                         u0[uOff + 2] + (0.5 * h * udot0[uOff + 2])); // w_FM in F, step start
        }
        std::array<Real, 4> qNew;
        advanceQuatExp(&q0[qOff], wHalf, h, qNew.data()); // overrides the 4 quaternion slots
        q[qOff + 0] = qNew[0];
        q[qOff + 1] = qNew[1];
        q[qOff + 2] = qNew[2];
        q[qOff + 3] = qNew[3];
        // translation block (Free q[qOff+4..6], FreeLine q[qOff+4..6]) keeps its
        // Taylor update from the loop above; Ball has no translation block.
    }
    normalizeQuaternions(m, s); // mops up ~1e-16 rounding after the exp-map; never rescues an overflow

    // Cartesian solvent position drift: x1 = x0 + h v0 + (h^2/2) a0, a0 = f0/m.
    // Written straight into posG; fillAtomPositionsFromBodies (below) skips these
    // atoms, so the drift survives and feeds the OpenMM force eval.
    {
        const Real h2half = Real(0.5) * h * h;
        for (int j = 0; j < nSolv; ++j) {
            const int a = solvAtoms[j];
            const Real im = solvInvM[j];
            posG[a] = Vec3(xs0[j][0] + h * vs0[j][0] + h2half * im * fs0[j][0],
                           xs0[j][1] + h * vs0[j][1] + h2half * im * fs0[j][1],
                           xs0[j][2] + h * vs0[j][2] + h2half * im * fs0[j][2]);
        }
    }

    auto refreshPos = [&]() {
        realizePosition(m, s);
        fillAtomPositionsFromBodies(m, s);
    };
    refreshPos();
    cset.enforcePositionConstraints(m, s, refreshPos); // localProjectQ

    // ---- velocity: implicit trapezoid + functional iteration ----
    for (int i = 0; i < nu; ++i) {
        u[i] = u0[i] + (h * udot0[i]); // u1_est
    }

    auto evalDerivs = [&]() -> bool {
        realizePosition(m, s);
        ROBO_CHECK("realizePosition");
        bridge.evaluate(s);
        ROBO_CHECK("bridge.evaluate");
        if (!forcesFinite()) {
            return false; // non-finite force -> abort before it corrupts udot/q
        }
        realizeVelocity(m, s);
        ROBO_CHECK("realizeVelocity");
        realizeArticulatedBodyInertias(m, s);
        calcUDot(m, s);
        ROBO_CHECK("calcUDot");
        for (int i = 0; i < nu; ++i) {
            if (!std::isfinite(udot[i])) {
                return false; // finite-but-diverging udot -> reject before it propagates
            }
        }
        calcQDot(m, s, qdot);
        calcQDotDot(m, s);
        return true;
    };
    if (!evalDerivs()) {
        restorePreStep();
        return false;
    }

    // Simbody: tol = min(1e-4, 0.1*accuracy). For fixed-step HMC (default accuracy
    // ~1e-3) this evaluates to 1e-4. Plain functional iteration, no under-relaxation,
    // max 10 sweeps -- matching VerletIntegrator::attemptDAEStep exactly.
    const Real tol = Real(1e-4);
    Real prevChange = std::numeric_limits<Real>::infinity();
    Real lastChange = std::numeric_limits<Real>::infinity(); // for the dt-too-large message
    int usedIters = 0;
    bool converged = false;

    for (int iter = 0; iter < 10; ++iter) {
        ++usedIters;
        Real num = 0;
        Real den = 0; // Simbody's relative 2-norm change

        for (int i = 0; i < nu; ++i) {
            const Real un = u0[i] + ((h / 2) * (udot0[i] + udot[i]));
            const Real d = un - u[i];
            num += d * d;
            den += u[i] * u[i];
            u[i] = un;
        }

        // Genuine non-finite / runaway. In Simbody this is the realize()/project()
        // exception path: caught, and in fixed-step mode the step still "succeeds",
        // propagating the bad state to the move-level energy validation, which
        // rejects the MOVE. We short-circuit to the same outcome by rejecting here.
        if (!velocitiesSane()) {
            restorePreStep();
            return false;
        }
        if (!evalDerivs()) {
            restorePreStep();
            return false;
        }

        const Real change = std::sqrt(num) / (std::sqrt(den) + Real(1e-30));
        lastChange = change;
        if (!std::isfinite(change) || change > Real(1e6)) {
            restorePreStep();
            return false;
        }

        if (change <= tol) {
            converged = true;
            break; // converged
        }

        // Functional iteration stopped contracting (after iter > 1, to skip the
        // crude forward-Euler seed's first non-monotone blip). We stop iterating
        // here; whether to take the step or reject is decided after the loop.
        if (iter > 1 && change > prevChange) {
            break;
        }

        prevChange = change;
    }

    // dt-too-large guard (deliberately STRICTER than Simbody's "take the step").
    // Reaching here without convergence means the implicit-trapezoid corrector
    // could not find its fixed point at this dt for the CURRENT configuration --
    // a FINITE, bounded solve that simply will not contract. This is distinct from
    // a steric clash (non-finite force / runaway u), which is handled above as a
    // move rejection (return false) because clashes are a normal, transient part of
    // sampling.
    //
    // THEORY (boilerplate "Propagation" / RATTLE): the fixed step is taken
    // UNCONDITIONALLY -- "a non-converged corrector does not shrink dt or reject
    // the step ... the best available velocity estimate is accepted", and "it is
    // the trajectory-level Metropolis test, not per-step control, that supplies
    // correctness." So on non-convergence we KEEP the advanced position q (set by
    // the position drift before this loop; the corrector only refines u) and the
    // last/best velocity iterate, and return success. A step that pumps energy
    // then shows up as a large dH and is rejected at the MOVE level (the caller
    // restores its saved q), which is the correct, theory-sanctioned behavior.
    //
    // The previous code instead called restorePreStep() here and still returned
    // true: that silently UNDID the step (q == q0) while reporting success, so the
    // move-level dH was always ~0 and every proposal was "accepted" as a no-op --
    // i.e. the trajectory froze in place. That both violated the theory above and
    // produced exactly that frozen-coordinate symptom; it is removed. The warning
    // is kept (fail-loud diagnostic): a non-converged corrector means dt is too
    // large for this geometry, so most such steps will be rejected by Metropolis
    // until the world's timestep is reduced.
    if (correctorConverged) {
        *correctorConverged = converged;
    }
    if (!converged) {
        std::fprintf(stderr,
                     "[verlet] world: velocity corrector did not converge at dt=%.6g ps "
                     "(relative change %.3e > tol %.3e after %d iterations). The timestep is "
                     "likely too large for the current configuration; the step is still taken "
                     "(THEORY: trajectory-level Metropolis supplies correctness) but expect "
                     "rejections until the timestep is reduced.\n",
                     (double)h,
                     (double)lastChange,
                     (double)tol,
                     usedIters);
    }

    cset.enforceVelocityConstraints(m, s); // localProjectU (RATTLE)
    realizeVelocity(m, s);                 // refresh V/KE at the projected u

    // Cartesian solvent velocity half: v1 = v0 + (h/2)(a0 + a1), a1 = f1/m. The
    // positions are frozen through the solute corrector, so f1 = frcG (evaluated
    // at the drifted x1 during evalDerivs) is the end-of-step force -- the update
    // is explicit and exact, no iteration. This completes the symmetric Verlet.
    for (int j = 0; j < nSolv; ++j) {
        const int a = solvAtoms[j];
        const Real hHalfInvM = Real(0.5) * h * solvInvM[j];
        velG[a] = Vec3(vs0[j][0] + hHalfInvM * (fs0[j][0] + frcG[a][0]),
                       vs0[j][1] + hHalfInvM * (fs0[j][1] + frcG[a][1]),
                       vs0[j][2] + hHalfInvM * (fs0[j][2] + frcG[a][2]));
    }

    s.time += h;
    return true;
}

template <class Bridge>
auto RobotEngine::stepTo(const RobotModel& model,
                         RobotState& state,
                         Bridge& bridge,
                         const robo::ConstraintSet& cset,
                         Real tEnd,
                         bool* correctorConverged) -> bool {
    const Real h = tEnd - state.time;
    if (h <= 0) {
        if (correctorConverged) {
            *correctorConverged = true; // no-op step: trivially "converged"
        }
        return true;
    }
    return verletStep(model, state, bridge, cset, h, correctorConverged);
}

template <class Bridge>
// ============================================================================
//  Reversibility diagnostic (HMC proposal sanity check)
// ============================================================================
// Integrate nSteps forward at fixed step h, flip every generalized speed
// (u -> -u), integrate nSteps "back", flip again. For a time-reversible map the
// state returns to its start to within ~machine epsilon amplified by the work
// done; a too-large h (the non-reversible, energy-pumTorsiong regime) returns a
// residual of order the trajectory size. Returns the RELATIVE round-trip
// residual ||(q,u)_returned - (q,u)_start|| / ||(q,u)_start||.
//
// Two properties of this probe matter for how it is used:
//   * It is NON-DESTRUCTIVE: the state s is restored to its entry value before
//     returning, so it can be called at startup or mid-run without perturbing
//     the chain.
//   * It certifies h ONLY for the configuration it is run from. The integrator
//     map depends on q through the configuration-dependent mass metric M(q) and
//     the local force stiffness, so the largest reversible h varies across
//     configuration space (THEORY 5.5). A startup call is therefore a smoke test
//     -- "is dt in the right ballpark for this geometry?" -- NOT a whole-run
//     guarantee. The ongoing, per-configuration guard is the corrector-
//     convergence throw inside verletStep, which tests the CURRENT geometry on
//     every step.
//
// If the corrector throws mid-probe (dt already non-convergent here), that is
// caught and reported as an infinite residual rather than aborting the probe.
//
// Quaternion DOF are compared with a double-cover-aware chordal distance
// min(|q - q0|, |q + q0|), since q and -q are the same rotation on S^3 (THEORY
// 3.4); scalar coordinates and the velocity vector use the plain Euclidean norm.
Real RobotEngine::checkReversibility(const RobotModel& m,
                                     RobotState& s,
                                     Bridge& bridge,
                                     const robo::ConstraintSet& cset,
                                     int nSteps,
                                     Real h) {
    const int nq = m.nq, nu = m.nu;

    // Snapshot the entry state (for comparison and for restoration).
    std::vector<Real> qStart(s.q(), s.q() + nq);
    std::vector<Real> uStart(s.u(), s.u() + nu);
    std::vector<Real> qdotStart(s.qdot(), s.qdot() + nq);
    std::vector<Real> udotStart(s.udot(), s.udot() + nu);
    std::vector<Real> qddStart(s.qdotdot(), s.qdotdot() + nq);
    const Real timeStart = s.time;

    // Cartesian solvent (solvent-relaxing NCMC): the joint map (q,u,x_s,v_s) must
    // be momentum-flip reversible too, so the probe carries the solvent DOF.
    const std::vector<int>& solvAtoms = s.cartSolventAtoms();
    const int nSolv = static_cast<int>(solvAtoms.size());
    Vec3* posG = s.atomPosG();
    Vec3* velG = s.atomVelG();
    std::vector<Vec3> xsStart(nSolv), vsStart(nSolv);
    for (int j = 0; j < nSolv; ++j) {
        xsStart[j] = posG[solvAtoms[j]];
        vsStart[j] = velG[solvAtoms[j]];
    }

    auto restoreStart = [&]() {
        std::copy(qStart.begin(), qStart.end(), s.q());
        std::copy(uStart.begin(), uStart.end(), s.u());
        std::copy(qdotStart.begin(), qdotStart.end(), s.qdot());
        std::copy(udotStart.begin(), udotStart.end(), s.udot());
        std::copy(qddStart.begin(), qddStart.end(), s.qdotdot());
        s.time = timeStart;
        realizePosition(m, s);
        realizeVelocity(m, s);
        fillAtomPositionsFromBodies(m, s);
        for (int j = 0; j < nSolv; ++j) {
            posG[solvAtoms[j]] = xsStart[j];
            velG[solvAtoms[j]] = vsStart[j];
        }
    };

    // Momentum-flip involution S: negate all generalized speeds and refresh the
    // velocity-dependent derived quantities against -u. Cartesian solvent speeds
    // are flipped in lockstep (the joint involution F).
    auto flipMomenta = [&]() {
        Real* u = s.u();
        for (int i = 0; i < nu; ++i) {
            u[i] = -u[i];
        }
        for (int j = 0; j < nSolv; ++j) {
            Vec3& v = velG[solvAtoms[j]];
            v = Vec3(-v[0], -v[1], -v[2]);
        }
        realizeVelocity(m, s);
        realizeArticulatedBodyInertias(m, s);
        calcUDot(m, s);
        calcQDot(m, s, s.qdot());
        calcQDotDot(m, s);
    };

    // Forward leg, flip, back leg, flip. Any non-finite rejection or corrector
    // throw means h is already not integrable/reversible here -> infinite residual.
    try {
        for (int i = 0; i < nSteps; ++i) {
            if (!verletStep(m, s, bridge, cset, h)) {
                restoreStart();
                return std::numeric_limits<Real>::infinity();
            }
        }
        flipMomenta();
        for (int i = 0; i < nSteps; ++i) {
            if (!verletStep(m, s, bridge, cset, h)) {
                restoreStart();
                return std::numeric_limits<Real>::infinity();
            }
        }
        flipMomenta();
    } catch (const std::exception&) {
        restoreStart();
        return std::numeric_limits<Real>::infinity();
    }

    // Round-trip residual. Quaternion blocks: double-cover chordal distance.
    const Real* q = s.q();
    const Real* u = s.u();
    Real resid2 = 0, scale2 = 0;

    std::vector<bool> isQuatSlot(nq, false);
    for (int b = 1; b < m.numBodies; ++b) {
        if (!m.isQuaternionBody(b)) {
            continue;
        }
        const int qOff = m.bodyQIndex[b];
        Real dPlus2 = 0, dMinus2 = 0;
        for (int k = 0; k < 4; ++k) {
            isQuatSlot[qOff + k] = true;
            const Real a = q[qOff + k], b0 = qStart[qOff + k];
            dPlus2 += (a - b0) * (a - b0);
            dMinus2 += (a + b0) * (a + b0);
            scale2 += b0 * b0;
        }
        resid2 += std::min(dPlus2, dMinus2); // q == -q on S^3
    }
    for (int i = 0; i < nq; ++i) {
        if (isQuatSlot[i]) {
            continue;
        }
        const Real d = q[i] - qStart[i];
        resid2 += d * d;
        scale2 += qStart[i] * qStart[i];
    }
    for (int i = 0; i < nu; ++i) {
        const Real d = u[i] - uStart[i];
        resid2 += d * d;
        scale2 += uStart[i] * uStart[i];
    }
    // Cartesian solvent positions and velocities (plain Euclidean -- flat space).
    for (int j = 0; j < nSolv; ++j) {
        const Vec3& x = posG[solvAtoms[j]];
        const Vec3& v = velG[solvAtoms[j]];
        for (int k = 0; k < 3; ++k) {
            const Real dx = x[k] - xsStart[j][k];
            const Real dv = v[k] - vsStart[j][k];
            resid2 += dx * dx + dv * dv;
            scale2 += xsStart[j][k] * xsStart[j][k] + vsStart[j][k] * vsStart[j][k];
        }
    }

    restoreStart(); // non-destructive probe
    return std::sqrt(resid2) / (std::sqrt(scale2) + Real(1e-30));
}

#ifdef ROBO_INTEGRATOR_LOCAL_ROBO_CHECK
#    undef ROBO_CHECK
#    undef ROBO_INTEGRATOR_LOCAL_ROBO_CHECK
#endif