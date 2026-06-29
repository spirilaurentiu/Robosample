// ============================================================================
//  TestBiasForces.cpp -- Phase 2: pin the VELOCITY-DEPENDENT bias terms of the
//  articulated-body dynamics in isolation -- gyroscopic, centrifugal /
//  centripetal, and Coriolis. These are computed in realizeVelocity
//  (s.gyro(), s.mobCoriolisA(), s.coriolisA()) and consumed by calcUDot.
//
//  WHY THIS FILE EXISTS. Until now these terms were touched only INDIRECTLY: the
//  one quantitative dynamics check (Dynamics.UDotEqualsMInvForceAtZeroVelocity,
//  TestMassMatrix.cpp) deliberately ZEROES velocity, so every bias term is 0 in
//  it and a regression there would be invisible. The realizeVelocity comments
//  document two historical energy-pumping bugs that live exactly here -- the
//  dropped Hdot_G u term ("EARLIER BUG") and the missing centripetal HDot_MB_F
//  term -- both of which conserve nothing yet pass a position-only OpenMM energy
//  check. This file exercises each term at NONZERO u, against an independent
//  re-derivation, plus a residual force balance and a system-level angular-
//  momentum invariant. Each test carries a must-FAIL-if-dropped guard so a
//  silent regression to the documented bugs turns the suite red.
//
//  No World / forcefield: generalized forces are set directly (calcUDot), and
//  the integrator tests use the OpenMM-free ZeroBridge below.
//
//  NOT PORTED (no operator in Robosample): reaction forces
//  (calcMobilizerReactionForces); gravity (none in Robosample).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "Constraints.hpp"
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotIntegrator.hpp" // templated verletStep (FreeBodyConservesAngularMomentum)
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::NearVec3;
using rtest::randomizeState;
using rtest::Rng;

namespace {

// ---- local tolerances (the spec's kFDbias / kEdrift; kTight/kLoose are shared) -
inline constexpr Real kFDbias = 1e-6; // central-FD bias residual (observed ~1e-9)
inline constexpr Real kEdrift = 1e-6; // free-body L drift over 2000 steps at small h

inline Real spatialNorm(const SpatialVec& v) {
    return std::sqrt(v[0].normSqr() + v[1].normSqr());
}

// makeDynForest-style mixed-joint builder (modeled on TestMassMatrix.cpp). Every
// body gets random-but-fixed joint frames; X_BM has a NONZERO translation, so
// r_MB = X_MB.p() != 0 and the centripetal HDot_MB_F term is active on every
// body (in particular the Torsion bodies -- the spec's "Torsion-with-nonzero-
// X_BM.p()" requirement). All ten joint types appear so CoriolisAccelMatches-
// ClassicForm exercises both the constant-H closed form (seven joints) and the
// q-dependent FD branch (BendStretch / SphericalCoords / FreeLine).
RobotModel makeBiasForest(Rng& rng) {
    auto F = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3(-0.6, 0.6)); // nonzero r_MB
        s.mass = rng.uniform(0.5, 2.5);
        s.com_B = rng.vec3(-0.25, 0.25);
        s.inertia_B = UnitInertia(rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7));
        return s;
    };
    std::vector<BodySpec> specs;
    specs.push_back(F(0, JointType::Free));            // 1  (robot 1 root)
    specs.push_back(F(1, JointType::Torsion));         // 2  Torsion, offset M frame, parent moving
    specs.push_back(F(2, JointType::Torsion));         // 3  Torsion (chain)
    specs.push_back(F(0, JointType::Ball));            // 4  (robot 2 root)
    specs.push_back(F(4, JointType::Cylinder));        // 5
    specs.push_back(F(0, JointType::Cartesian));       // 6  (robot 3 root)
    specs.push_back(F(6, JointType::BendStretch));     // 7  q-dependent H
    specs.push_back(F(6, JointType::SphericalCoords)); // 8  q-dependent H
    specs.push_back(F(0, JointType::FreeLine));        // 9  q-dependent H (robot 4 root)
    specs.push_back(F(1, JointType::Slider));          // 10
    return buildForest(specs);
}

void realizeForBias(const RobotModel& m, RobotState& s) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
}

// Independent re-derivation of the mobilizer bias a_mob[b] by finite difference,
// valid for EVERY joint type (constant-H and q-dependent alike):
//   a_tot[b] = d/dt V_GB[b] holding u fixed (udot==0), and
//   a_mob[b] = a_tot[b] - ~Phi[b] * a_tot[parent].
// We FD V_GB along q(t)=q0+t*qdot (quaternions renormalized, u restored each
// sample) and subtract the parent's shifted FD using the CENTER-point Phi.
// Returns a_mob_fd indexed by body. Restores s to its center configuration.
std::vector<SpatialVec> mobBiasByFD(const RobotModel& m, RobotState& s, Real h = 1e-6) {
    const int nb = m.numBodies;
    std::vector<Real> q0(s.q(), s.q() + m.nq);
    std::vector<Real> u(s.u(), s.u() + m.nu);
    std::vector<Real> qdot(m.nq);
    RobotEngine::calcQDot(m, s, qdot.data());
    std::vector<PhiMatrix> Phi(s.Phi(), s.Phi() + nb); // center-point Phi

    auto velAt = [&](Real t) {
        for (int i = 0; i < m.nq; ++i) {
            s.q()[i] = q0[i] + t * qdot[i];
        }
        RobotEngine::normalizeQuaternions(m, s);
        RobotEngine::realizePosition(m, s);
        std::copy(u.begin(), u.end(), s.u());
        RobotEngine::realizeVelocity(m, s);
        return std::vector<SpatialVec>(s.V_GB(), s.V_GB() + nb);
    };
    const std::vector<SpatialVec> Vp = velAt(+h);
    const std::vector<SpatialVec> Vm = velAt(-h);

    // restore center
    std::copy(q0.begin(), q0.end(), s.q());
    std::copy(u.begin(), u.end(), s.u());
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeVelocity(m, s);

    const Real inv2h = Real(1) / (Real(2) * h);
    std::vector<SpatialVec> dV(nb), aMob(nb, SpatialVec(Vec3(0), Vec3(0)));
    for (int b = 0; b < nb; ++b) {
        dV[b] = SpatialVec((Vp[b][0] - Vm[b][0]) * inv2h, (Vp[b][1] - Vm[b][1]) * inv2h);
    }
    for (int b = 1; b < nb; ++b) {
        const SpatialVec shifted = (~Phi[b]) * dV[m.bodyParent[b]];
        aMob[b] = SpatialVec(dV[b][0] - shifted[0], dV[b][1] - shifted[1]);
    }
    return aMob;
}

// OpenMM-free zero-force bridge matching the verletStep contract: clears every
// per-body spatial force and every mobility force, exactly like ForceBridge does
// before accumulating (here it accumulates nothing). With zero applied force the
// ONLY thing driving udot is the velocity-dependent bias -- which is the point.
struct ZeroBridge {
    const RobotModel& m;
    explicit ZeroBridge(const RobotModel& model)
        : m(model) {
    }
    void evaluate(RobotState& s) {
        for (int b = 0; b < m.numBodies; ++b) {
            s.bodyForceG()[b] = SpatialVec(Vec3(0), Vec3(0));
        }
        for (int i = 0; i < m.nu; ++i) {
            s.mobilityForce()[i] = Real(0);
        }
    }
};

// Spatial angular momentum of a single rigid body about its OWN COM, in Ground.
// Spatial momentum about the body origin Bo is h = Mk_G * V_GB = [L_Bo ; p];
// shifting the reference point Bo -> COM (offset c = Bo->COM in Ground):
//   L_com = L_Bo - c x p .
// For a torque-free free body this is the conserved quantity (it equals
// I_com * w in the inertial frame and is constant; w itself precesses).
Vec3 angularMomentumAboutCom(const RobotModel& m, RobotState& s, int b) {
    const SpatialVec hmom = s.Mk_G()[b] * s.V_GB()[b];
    const Vec3 p = hmom[1];                     // linear momentum m * v_com
    const Vec3 c = s.Mk_G()[b].getMassCenter(); // Bo -> COM in Ground
    return hmom[0] - (c % p);
}

} // namespace

// ---------------------------------------------------------------------------
//  B1: GyroscopicMatchesClosedForm. For every body of a random forest at random
//  u, s.gyro()[b] must equal the independent closed form
//      g = mass * [ w x (G w) ; w x (w x c) ]
//  built from Mk_G()[b].getUnitInertia()/getMassCenter()/getMass() and V_GB()[b].
//  FAIL guard: the body-frame variant mass * w x (I_body w) -- forgetting the
//  reexpress(~R_GB) into Ground -- must differ by > 1e-6 for at least one body
//  (otherwise the test could not catch that exact regression).
// ---------------------------------------------------------------------------
TEST(BiasForces, GyroscopicMatchesClosedForm) {
    Rng rng(0xB1A5);
    RobotModel m = makeBiasForest(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 20; ++rep) {
        randomizeState(m, s, rng);
        realizeForBias(m, s);

        const SpatialInertia* Mk = s.Mk_G();
        const SpatialVec* V = s.V_GB();
        const SpatialVec* gyro = s.gyro();

        Real maxWrong = 0; // body-frame (un-reexpressed) deviation, for the guard
        for (int b = 1; b < m.numBodies; ++b) {
            const Vec3 w = V[b][0];
            const Vec3 Gw = Mk[b].getUnitInertia() * w; // G w (Ground), per unit mass
            const Vec3 c = Mk[b].getMassCenter();       // COM offset in Ground
            const SpatialVec g = Mk[b].getMass() * SpatialVec(w % Gw, w % (w % c));

            EXPECT_TRUE(NearVec3(gyro[b][0], g[0], rtest::kTight)) << "ang body " << b << " rep " << rep;
            EXPECT_TRUE(NearVec3(gyro[b][1], g[1], rtest::kTight)) << "lin body " << b << " rep " << rep;

            // WRONG variant: unit inertia taken in the BODY frame (no reexpress).
            const Vec3 GwBody = m.bodyUnitInertia_B[b] * w;
            const Vec3 gWrongAng = Mk[b].getMass() * (w % GwBody);
            maxWrong = std::max(maxWrong, (gWrongAng - gyro[b][0]).norm());
        }
        EXPECT_GT(maxWrong, 1e-6) << "body-frame gyroscopic indistinguishable from Ground form, rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  B2: BiasVanishesAtZeroVelocity. With u = 0, every bias term is a product of
//  zero velocities, so gyro[b], mobCoriolisA[b], coriolisA[b] must all be
//  exactly Vec3(0) (to kTight) for every body. This is the boundary that the
//  pre-existing UDotEqualsMInvForceAtZeroVelocity relies on implicitly; pinning
//  it makes the assumption explicit.
// ---------------------------------------------------------------------------
TEST(BiasForces, BiasVanishesAtZeroVelocity) {
    Rng rng(0xB2A5);
    RobotModel m = makeBiasForest(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 10; ++rep) {
        randomizeState(m, s, rng);
        std::fill(s.u(), s.u() + m.nu, Real(0)); // zero velocity
        realizeForBias(m, s);

        for (int b = 1; b < m.numBodies; ++b) {
            EXPECT_TRUE(NearVec3(s.gyro()[b][0], Vec3(0), rtest::kTight)) << "gyro ang b" << b;
            EXPECT_TRUE(NearVec3(s.gyro()[b][1], Vec3(0), rtest::kTight)) << "gyro lin b" << b;
            EXPECT_TRUE(NearVec3(s.mobCoriolisA()[b][0], Vec3(0), rtest::kTight)) << "aMob ang b" << b;
            EXPECT_TRUE(NearVec3(s.mobCoriolisA()[b][1], Vec3(0), rtest::kTight)) << "aMob lin b" << b;
            EXPECT_TRUE(NearVec3(s.coriolisA()[b][0], Vec3(0), rtest::kTight)) << "aTot ang b" << b;
            EXPECT_TRUE(NearVec3(s.coriolisA()[b][1], Vec3(0), rtest::kTight)) << "aTot lin b" << b;
        }
    }
}

// ---------------------------------------------------------------------------
//  B3: BiasIsNonzeroWhenItMustBe -- the must-FAIL-if-dropped control. For bodies
//  whose mobilizer bias is genuinely nonzero (a Torsion child with an offset M
//  frame and a moving parent; an offset Ball/Free root) at random nonzero u,
//  ||mobCoriolisA[b]|| must exceed 1e-6. Dropping the centripetal_G / Hdot_G u
//  term (the documented "EARLIER BUG") collapses these to zero and turns this
//  red.
// ---------------------------------------------------------------------------
TEST(BiasForces, BiasIsNonzeroWhenItMustBe) {
    Rng rng(0xB3A5);
    // A focused forest: Free root (offset), Torsion child (offset, moving parent),
    // Ball root (offset). All three have a structurally nonzero a_mob.
    auto F = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3(-0.6, 0.6));
        s.mass = rng.uniform(0.8, 2.0);
        s.com_B = rng.vec3(-0.2, 0.2);
        s.inertia_B = UnitInertia(0.4, 0.5, 0.6);
        return s;
    };
    RobotModel m = buildForest({F(0, JointType::Free), F(1, JointType::Torsion), F(0, JointType::Ball)});
    RobotState s;
    s.allocateFull(m);
    const int freeRoot = 1, torsionChild = 2, ballRoot = 3;

    for (int rep = 0; rep < 20; ++rep) {
        randomizeState(m, s, rng); // u ~ gaussian, nonzero
        realizeForBias(m, s);
        EXPECT_GT(spatialNorm(s.mobCoriolisA()[freeRoot]), 1e-6) << "Free root, rep " << rep;
        EXPECT_GT(spatialNorm(s.mobCoriolisA()[torsionChild]), 1e-6) << "Torsion child, rep " << rep;
        EXPECT_GT(spatialNorm(s.mobCoriolisA()[ballRoot]), 1e-6) << "Ball root, rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  B4: CoriolisAccelMatchesClassicForm. Re-derive the mobilizer bias from the
//  classic centripetal + 2*Coriolis form using only parent body velocities and
//  static frames -- independent of the engine's extraAng/extraLin path:
//    a_mob.ang = w_GP x w_PB_G
//    a_mob.lin = w_GP x (v_GB - v_GP) + w_GP x v_PB_G + R_GF*(w_FM x (w_FM x r_MB_F))
//  This closed form is COMPLETE only when H_FM is constant in F (seven joints),
//  so it is asserted for those. For the three q-dependent joints (BendStretch /
//  SphericalCoords / FreeLine) the closed form misses the intrinsic dH_FM/dt u
//  term, so a_mob is instead checked against a finite difference of the body
//  velocity (== d/dt(H_PB_G u) after removing the parent's shifted acceleration).
// ---------------------------------------------------------------------------
TEST(BiasForces, CoriolisAccelMatchesClassicForm) {
    Rng rng(0xB4A5);
    RobotModel m = makeBiasForest(rng);
    RobotState s;
    s.allocateFull(m);

    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        realizeForBias(m, s);

        const SpatialVec* V = s.V_GB();
        const SpatialVec* V_FM = s.V_FM();
        const SpatialVec* V_PB_G = s.V_PB_G();
        const Transform* X_GB = s.X_GB();
        const Transform* X_FM = s.X_FM();
        const SpatialVec* aMob = s.mobCoriolisA();

        // FD reference for the q-dependent joints (also valid for the rest).
        const std::vector<SpatialVec> aMobFD = mobBiasByFD(m, s);

        for (int b = 1; b < m.numBodies; ++b) {
            const int p = m.bodyParent[b];
            if (RobotModel::jointHasConstantHFM(m.bodyJoint[b])) {
                const Vec3& w_GP = V[p][0];
                const Vec3& v_GP = V[p][1];
                const Vec3& v_GB = V[b][1];
                const Vec3& w_PB_G = V_PB_G[b][0];
                const Vec3& v_PB_G = V_PB_G[b][1];
                const Vec3& w_FM = V_FM[b][0];
                const Vec3 r_MB_F = X_FM[b].R() * (~m.X_BM[b]).p(); // Mo->Bo in F
                const Rotation R_GF = X_GB[p].R() * m.X_PF[b].R();
                const Vec3 centripetal_G = R_GF * (w_FM % (w_FM % r_MB_F));

                const Vec3 ang = w_GP % w_PB_G;
                const Vec3 lin = (w_GP % (v_GB - v_GP)) + (w_GP % v_PB_G) + centripetal_G;
                EXPECT_TRUE(NearVec3(aMob[b][0], ang, rtest::kTight))
                    << "classic ang body " << b << " rep " << rep;
                EXPECT_TRUE(NearVec3(aMob[b][1], lin, rtest::kTight))
                    << "classic lin body " << b << " rep " << rep;
            } else {
                // q-dependent: only the FD oracle captures the dH_FM/dt u term.
                EXPECT_TRUE(NearVec3(aMob[b][0], aMobFD[b][0], kFDbias))
                    << "FD ang body " << b << " joint " << static_cast<int>(m.bodyJoint[b]) << " rep " << rep;
                EXPECT_TRUE(NearVec3(aMob[b][1], aMobFD[b][1], kFDbias))
                    << "FD lin body " << b << " joint " << static_cast<int>(m.bodyJoint[b]) << " rep " << rep;
            }
        }
    }
}

// ---------------------------------------------------------------------------
//  B5: ResidualForceBalanceAtNonzeroVelocity. The missing companion to
//  UDotEqualsMInvForceAtZeroVelocity: drive calcUDot at NONZERO u with random
//  mobilityForce (tau) and zero bodyForceG, then verify the equation of motion
//      M udot = tau - c        (c = velocity-dependent gyroscopic+Coriolis bias)
//  in the cleaner equivalent form actually asserted here:
//      udot == multiplyByMInv(tau - c).
//  We obtain M^-1 (tau - c) without ever forming c explicitly: run calcUDot with
//  tau=0 to get udot_bias = M^-1 (-c), and use linearity of M^-1:
//      udot(tau) == multiplyByMInv(tau) + udot_bias.
//  FAIL guard: forcing bias to 0 (i.e. comparing udot(tau) to multiplyByMInv(tau)
//  alone) must NOT pass at nonzero u -- proving the bias actually enters udot.
// ---------------------------------------------------------------------------
TEST(Dynamics, ResidualForceBalanceAtNonzeroVelocity) {
    Rng rng(0xB5A5);
    RobotModel m = makeBiasForest(rng);
    RobotState s;
    s.allocateFull(m);

    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng); // nonzero u
        realizeForBias(m, s);
        std::fill(s.bodyForceG(), s.bodyForceG() + m.numBodies, SpatialVec(Vec3(0), Vec3(0)));

        std::vector<Real> tau(m.nu);
        for (int i = 0; i < m.nu; ++i) {
            tau[i] = rng.gaussian();
        }

        // udot at applied force tau (includes the bias).
        std::copy(tau.begin(), tau.end(), s.mobilityForce());
        RobotEngine::calcUDot(m, s);
        std::vector<Real> udotTau(s.udot(), s.udot() + m.nu);

        // udot at zero applied force = M^-1 (-c), the pure bias acceleration.
        std::fill(s.mobilityForce(), s.mobilityForce() + m.nu, Real(0));
        RobotEngine::calcUDot(m, s);
        std::vector<Real> udotBias(s.udot(), s.udot() + m.nu);

        // M^-1 tau via the explicit operator (no velocity terms).
        std::vector<Real> minvTau(m.nu);
        RobotEngine::multiplyByMInv(m, s, tau.data(), minvTau.data());

        // Equation of motion: udot(tau) == M^-1(tau) + M^-1(-c).
        Real maxBiasContribution = 0;
        for (int i = 0; i < m.nu; ++i) {
            EXPECT_NEAR(udotTau[i], minvTau[i] + udotBias[i], rtest::kLoose) << "i=" << i << " rep " << rep;
            maxBiasContribution = std::max(maxBiasContribution, std::abs(udotTau[i] - minvTau[i]));
        }
        // FAIL guard: bias-forced-to-0 (udotTau ?= minvTau) must NOT hold -- the
        // velocity-dependent term genuinely moves udot at nonzero u.
        EXPECT_GT(maxBiasContribution, 1e-6) << "bias did not enter udot, rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  B6: FreeBodyConservesAngularMomentum (system-level invariant). A single Free
//  body, NO applied force (so udot is driven entirely by the gyroscopic bias),
//  random initial spin, integrated 2000 steps with the production verletStep.
//  The spatial angular momentum about the COM in Ground is a constant of motion
//  for a torque-free body, so its relative drift must stay < kEdrift.
//  Sensitivity (FAIL) guard: a 10x-too-large h must drift past kEdrift.
// ---------------------------------------------------------------------------
TEST(BiasForces, FreeBodyConservesAngularMomentum) {
    auto relDrift = [](Real h, int nSteps) -> Real {
        // Deterministic single Free body with an ASYMMETRIC inertia (so w
        // precesses -- the gyroscopic term is genuinely exercised, not trivial).
        BodySpec spec;
        spec.parent = 0;
        spec.joint = JointType::Free;
        spec.X_PF = Transform();
        spec.X_BM = Transform();
        spec.mass = Real(2.0);
        spec.com_B = Vec3(0.07, -0.03, 0.05);
        spec.inertia_B = UnitInertia(Real(0.30), Real(0.55), Real(0.80));
        RobotModel m = buildForest({spec});

        RobotState s;
        s.allocateFull(m);
        Real* q = s.q();
        q[0] = 1; // identity quaternion
        q[1] = q[2] = q[3] = 0;
        q[4] = q[5] = q[6] = 0; // origin
        Real* u = s.u();
        u[0] = Real(1.3); // pure spin (w_FM), zero translational speed
        u[1] = Real(-0.8);
        u[2] = Real(0.5);
        u[3] = u[4] = u[5] = 0;

        ZeroBridge bridge(m);
        const ConstraintSet cset; // empty -> unconstrained

        // seed the derivative chain (mirrors TestIntegrator's seedDerivatives).
        RobotEngine::realizePosition(m, s);
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        bridge.evaluate(s);
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        RobotEngine::calcUDot(m, s);
        RobotEngine::calcQDot(m, s, s.qdot());
        RobotEngine::calcQDotDot(m, s);

        const Vec3 L0 = angularMomentumAboutCom(m, s, 1);
        for (int n = 0; n < nSteps; ++n) {
            const bool ok = RobotEngine::verletStep(m, s, bridge, cset, h);
            if (!ok) {
                return Real(1e9); // a rejected step is itself a failure of conservation
            }
        }
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);
        const Vec3 L1 = angularMomentumAboutCom(m, s, 1);
        return (L1 - L0).norm() / (L0.norm() + Real(1e-30));
    };

    const Real h = 5e-4;
    EXPECT_LT(relDrift(h, 2000), kEdrift) << "angular momentum drifted at the nominal step";
    // sensitivity: 10x the step must visibly break conservation.
    EXPECT_GT(relDrift(10 * h, 2000), kEdrift) << "10x step did not break conservation (insensitive)";
}