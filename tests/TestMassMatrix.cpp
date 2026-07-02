// ============================================================================
//  TestMassMatrix.cpp -- COMPOSITION of the articulated-body dynamics on a robot
//  forest: the matrix-free mass operators (M^-1, sqrt(M^-1), sqrt(M)), the
//  kinetic energy, the O(n) log-determinant, and forward dynamics (calcUDot).
//  These are the Robosample analogues of Simbody's matter-subsystem operator
//  tests. No World/forcefield: generalized forces are set directly.
//
//  The strong cross-checks here:
//    * sqrt(M) and sqrt(M^-1) are mutual inverses (Simbody guarantees this).
//    * M^-1 is symmetric (SPD), checked as a bilinear form.
//    * calcUDot at zero velocity equals M^-1 * tau == multiplyByMInv(tau): the
//      forward-dynamics solve and the explicit operator agree.
//    * KE(u = sqrt(M^-1) w) = 1/2 |w|^2 exactly (ties calcKineticEnergy to the
//      sqrt operator: if L L^T = M^-1 then u=Lw gives u^T M u = |w|^2).
//    * calcLogDetM equals -ln det(M^-1) built densely from multiplyByMInv.
// ============================================================================
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotLinearAlgebra.hpp" // jacobiSymEig / logDetSymPD (n <= 6 dense kernels)
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::randomizeState;
using rtest::Rng;

namespace {

// A forest with a healthy mix of dof so nu is large and M is non-trivial.
RobotModel makeDynForest(Rng& rng) {
    auto F = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        s.mass = rng.uniform(0.5, 2.5);
        s.com_B = rng.vec3(-0.25, 0.25);
        s.inertia_B = UnitInertia(rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7));
        return s;
    };
    std::vector<BodySpec> specs;
    specs.push_back(F(0, JointType::Free));      // 1
    specs.push_back(F(1, JointType::Torsion));   // 2
    specs.push_back(F(2, JointType::Torsion));   // 3 (chain)
    specs.push_back(F(0, JointType::Ball));      // 4 (second robot)
    specs.push_back(F(4, JointType::Cylinder));  // 5
    specs.push_back(F(0, JointType::Cartesian)); // 6 (third robot)
    return buildForest(specs);
}

void prep(const RobotModel& m, RobotState& s) {
    RobotEngine::realizePosition(m, s);
    RobotEngine::realizeVelocity(m, s);
    RobotEngine::realizeArticulatedBodyInertias(m, s);
}

Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real s = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        s += a[i] * b[i];
    }
    return s;
}

// log|det A| of a row-major n x n via partial-pivot Gaussian elimination.
Real logDet(std::vector<Real> A, int n) {
    Real ld = 0;
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r = col + 1; r < n; ++r) {
            if (std::abs(A[r * n + col]) > std::abs(A[piv * n + col])) {
                piv = r;
            }
        }
        if (piv != col) {
            for (int c = 0; c < n; ++c) {
                std::swap(A[piv * n + c], A[col * n + c]);
            }
        }
        const Real d = A[col * n + col];
        ld += std::log(std::abs(d));
        for (int r = col + 1; r < n; ++r) {
            const Real f = A[r * n + col] / d;
            for (int c = col; c < n; ++c) {
                A[r * n + c] -= f * A[col * n + c];
            }
        }
    }
    return ld;
}

// ---------------------------------------------------------------------------
//  Phase 4: an INDEPENDENT dense M, built from the system Jacobian -- the path
//  the existing suite never takes (M^-1 is only ever checked against udot and
//  against its own sqrt, never against a forward-formed M).
//
//  Column j of the system Jacobian J is V_GB over all bodies when u = e_j (one
//  realizeVelocity per basis vector; V_GB = J(q) u is linear in u). Then, in the
//  spatial-inertia metric,
//      M_ij = sum_b dot( V_GB^(i)[b], Mk_G[b] * V_GB^(j)[b] ).
//  dot(SpatialVec,SpatialVec) = w.tau + v.f is the existing pairing, and
//  Mk_G[b] * V_GB[b] is the spatial momentum. This uses ONLY forward kinematics
//  (realizeVelocity) and the body spatial inertias -- never the articulated-body
//  inverse recursion (multiplyByMInv / P / D / DI) -- so it is a genuinely
//  independent oracle for M.
//
//  The fixture is sized to nu == 6 so the dense SPD kernels jacobiSymEig /
//  logDetSymPD (fixed 6x6 buffers) apply directly.
// ---------------------------------------------------------------------------
RobotModel makeDenseForest(Rng& rng) {
    auto F = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        s.mass = rng.uniform(0.5, 2.5);
        s.com_B = rng.vec3(-0.25, 0.25);
        s.inertia_B = UnitInertia(rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7));
        return s;
    };
    // nu = 6: Ball (3) + three Torsions (1 each), a single coupled chain.
    return buildForest({F(0, JointType::Ball),
                        F(1, JointType::Torsion),
                        F(2, JointType::Torsion),
                        F(3, JointType::Torsion)});
}

// Build the dense nu x nu mass matrix from the system Jacobian (forward path).
// Requires realizePosition current; leaves u as it found it. Reads Mk_G from the
// state, so it reflects whatever spatial inertias the engine is currently using.
std::vector<Real> buildDenseM(const RobotModel& m, RobotState& s) {
    const int n = m.nu, B = m.numBodies;
    const SpatialInertia* Mk = s.Mk_G();

    // capture J's columns: V_GB when u = e_j (position held fixed throughout).
    std::vector<std::vector<SpatialVec>> col(static_cast<std::size_t>(n),
                                             std::vector<SpatialVec>(static_cast<std::size_t>(B)));
    std::vector<Real> uSave(s.u(), s.u() + n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            s.u()[i] = (i == j) ? Real(1) : Real(0);
        }
        RobotEngine::realizeVelocity(m, s);
        const SpatialVec* V = s.V_GB();
        for (int b = 0; b < B; ++b) {
            col[static_cast<std::size_t>(j)][static_cast<std::size_t>(b)] = V[b];
        }
    }
    std::copy(uSave.begin(), uSave.end(), s.u());

    std::vector<Real> M(static_cast<std::size_t>(n * n), 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Real acc = 0;
            for (int b = 1; b < B; ++b) { // body 0 == Ground: no inertia
                const SpatialVec MV = Mk[b] * col[static_cast<std::size_t>(j)][static_cast<std::size_t>(b)];
                acc += dot(col[static_cast<std::size_t>(i)][static_cast<std::size_t>(b)], MV);
            }
            M[static_cast<std::size_t>(i * n + j)] = acc;
        }
    }
    return M;
}

// smallest eigenvalue of a row-major symmetric n x n (n <= 6) via the in-house
// Jacobi solver -- used to assert M is positive-definite.
Real minEigSym(const std::vector<Real>& A, int n) {
    Real d[6], V[36];
    robo_linalg::jacobiSymEig(A.data(), n, d, V);
    Real lo = d[0];
    for (int k = 1; k < n; ++k) {
        lo = std::min(lo, d[k]);
    }
    return lo;
}

} // namespace

// ---------------------------------------------------------------------------
//  D1: sqrt(M) and sqrt(M^-1) are mutual inverses.
// ---------------------------------------------------------------------------
TEST(MassMatrix, SqrtOperatorsAreInverses) {
    Rng rng(0x511);
    RobotModel m = makeDynForest(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        prep(m, s);
        std::vector<Real> w(m.nu), u(m.nu), w2(m.nu);
        for (int i = 0; i < m.nu; ++i) {
            w[i] = rng.gaussian();
        }
        RobotEngine::multiplyBySqrtMInv(m, s, w.data(), u.data()); // u = L w
        RobotEngine::multiplyBySqrtM(m, s, u.data(), w2.data());   // w2 = L^-1 u
        for (int i = 0; i < m.nu; ++i) {
            EXPECT_NEAR(w2[i], w[i], rtest::kLoose) << "i=" << i << " rep " << rep;
        }
    }
}

// ---------------------------------------------------------------------------
//  D2: M^-1 is symmetric: x . M^-1 y == y . M^-1 x.
// ---------------------------------------------------------------------------
TEST(MassMatrix, MInvIsSymmetric) {
    Rng rng(0x522);
    RobotModel m = makeDynForest(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        prep(m, s);
        std::vector<Real> x(m.nu), y(m.nu), Mix(m.nu), Miy(m.nu);
        for (int i = 0; i < m.nu; ++i) {
            x[i] = rng.gaussian();
            y[i] = rng.gaussian();
        }
        RobotEngine::multiplyByMInv(m, s, x.data(), Mix.data());
        RobotEngine::multiplyByMInv(m, s, y.data(), Miy.data());
        EXPECT_NEAR(dot(y, Mix), dot(x, Miy), rtest::kLoose) << "rep " << rep;
        // SPD: x . M^-1 x > 0 for x != 0.
        EXPECT_GT(dot(x, Mix), 0.0) << "rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  D3: forward dynamics at zero velocity equals the explicit M^-1 operator.
//      M udot = tau (bias is zero at u=0)  =>  udot == multiplyByMInv(tau).
// ---------------------------------------------------------------------------
TEST(Dynamics, UDotEqualsMInvForceAtZeroVelocity) {
    Rng rng(0x533);
    RobotModel m = makeDynForest(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        std::fill(s.u(), s.u() + m.nu, Real(0)); // zero velocity -> no Coriolis/gyro bias

        std::vector<Real> tau(m.nu);
        for (int i = 0; i < m.nu; ++i) {
            tau[i] = rng.gaussian();
        }

        // forward dynamics path
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        std::fill(s.bodyForceG(), s.bodyForceG() + m.numBodies, SpatialVec(Vec3(0), Vec3(0)));
        std::copy(tau.begin(), tau.end(), s.mobilityForce());
        RobotEngine::calcUDot(m, s);
        std::vector<Real> udot(m.nu);
        std::copy(s.udot(), s.udot() + m.nu, udot.begin());

        // explicit operator path
        std::vector<Real> minvTau(m.nu);
        RobotEngine::multiplyByMInv(m, s, tau.data(), minvTau.data());

        for (int i = 0; i < m.nu; ++i) {
            EXPECT_NEAR(udot[i], minvTau[i], rtest::kLoose) << "i=" << i << " rep " << rep;
        }
    }
}

// ---------------------------------------------------------------------------
//  D4: KE(u = sqrt(M^-1) w) = 1/2 |w|^2 exactly. Plus basic KE algebra.
// ---------------------------------------------------------------------------
TEST(MassMatrix, KineticEnergyMatchesSqrtMInv) {
    Rng rng(0x544);
    RobotModel m = makeDynForest(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);

        std::vector<Real> w(m.nu), u(m.nu);
        Real wsq = 0;
        for (int i = 0; i < m.nu; ++i) {
            w[i] = rng.gaussian();
            wsq += w[i] * w[i];
        }
        RobotEngine::multiplyBySqrtMInv(m, s, w.data(), u.data());
        std::copy(u.begin(), u.end(), s.u());
        RobotEngine::realizeVelocity(m, s);
        const Real ke = RobotEngine::calcKineticEnergy(m, s);
        EXPECT_NEAR(ke, 0.5 * wsq, rtest::kLoose) << "rep " << rep;

        // KE(-u) == KE(u); KE(2u) == 4 KE(u).
        for (int i = 0; i < m.nu; ++i) {
            s.u()[i] = -u[i];
        }
        RobotEngine::realizeVelocity(m, s);
        EXPECT_NEAR(RobotEngine::calcKineticEnergy(m, s), ke, rtest::kLoose);
        for (int i = 0; i < m.nu; ++i) {
            s.u()[i] = 2 * u[i];
        }
        RobotEngine::realizeVelocity(m, s);
        EXPECT_NEAR(RobotEngine::calcKineticEnergy(m, s), 4 * ke, rtest::kLoose);
    }
}

// ---------------------------------------------------------------------------
//  D5: calcLogDetM == -ln det(M^-1), with M^-1 formed densely from the operator.
//      Cross-checks the O(n) articulated-body determinant against O(n^3) algebra.
// ---------------------------------------------------------------------------
TEST(MassMatrix, LogDetMatchesDenseInverse) {
    Rng rng(0x555);
    RobotModel m = makeDynForest(rng);
    RobotState s;
    s.allocateFull(m);
    const int n = m.nu;
    for (int rep = 0; rep < 10; ++rep) {
        randomizeState(m, s, rng);
        prep(m, s);

        // dense M^-1 by applying the operator to each basis vector.
        std::vector<Real> Minv(n * n, 0), e(n, 0), col(n, 0);
        for (int j = 0; j < n; ++j) {
            std::fill(e.begin(), e.end(), Real(0));
            e[j] = 1;
            RobotEngine::multiplyByMInv(m, s, e.data(), col.data());
            for (int i = 0; i < n; ++i) {
                Minv[i * n + j] = col[i];
            }
        }
        const Real lnDetMInv = logDet(Minv, n);
        const Real lnDetM = RobotEngine::calcLogDetM(m, s);
        EXPECT_NEAR(lnDetM, -lnDetMInv, 1e-7) << "rep " << rep;
    }
}

// ===========================================================================
//  Phase 4 -- the independent dense M (forward Jacobian path).
// ===========================================================================

// ---------------------------------------------------------------------------
//  P4.1: M (forward Jacobian assembly) is symmetric and positive-definite.
//        M_ij == M_ji to kLoose; every eigenvalue > 0 (jacobiSymEig).
// ---------------------------------------------------------------------------
TEST(MassMatrixDense, DenseMIsSymmetricSPD) {
    Rng rng(0x4D01);
    RobotModel m = makeDenseForest(rng);
    RobotState s;
    s.allocateFull(m);
    const int n = m.nu;
    ASSERT_LE(n, 6) << "dense SPD kernels assume nu <= 6";
    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s); // fixes X_GB, H, Mk_G for the Jacobian build
        const std::vector<Real> M = buildDenseM(m, s);

        Real maxAsym = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                maxAsym = std::max(maxAsym, std::abs(M[i * n + j] - M[j * n + i]));
            }
        }
        EXPECT_LT(maxAsym, rtest::kLoose) << "M not symmetric, rep " << rep;

        EXPECT_GT(minEigSym(M, n), 0.0) << "M not positive-definite, rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  P4.2: M * (multiplyByMInv(M e_j)) == e_j for every basis vector. This closes
//        the loop the current suite leaves open: M^-1 is finally checked against
//        an independently-formed M, not against udot or its own sqrt.
// ---------------------------------------------------------------------------
TEST(MassMatrixDense, MTimesMInvIsIdentity) {
    Rng rng(0x4D02);
    RobotModel m = makeDenseForest(rng);
    RobotState s;
    s.allocateFull(m);
    const int n = m.nu;
    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        const std::vector<Real> M = buildDenseM(m, s);     // forward path
        RobotEngine::realizeArticulatedBodyInertias(m, s); // inverse path inputs (P/G/DI)

        for (int j = 0; j < n; ++j) {
            std::vector<Real> e(static_cast<std::size_t>(n), 0), Mij(static_cast<std::size_t>(n), 0);
            e[static_cast<std::size_t>(j)] = 1;
            // M e_j
            for (int i = 0; i < n; ++i) {
                Real acc = 0;
                for (int k = 0; k < n; ++k) {
                    acc += M[i * n + k] * e[static_cast<std::size_t>(k)];
                }
                Mij[static_cast<std::size_t>(i)] = acc;
            }
            // M^-1 (M e_j) must be e_j
            std::vector<Real> back(static_cast<std::size_t>(n), 0);
            RobotEngine::multiplyByMInv(m, s, Mij.data(), back.data());
            for (int i = 0; i < n; ++i) {
                EXPECT_NEAR(back[static_cast<std::size_t>(i)], e[static_cast<std::size_t>(i)], rtest::kLoose)
                    << "M^-1 M e_" << j << " != e_" << j << " at " << i << ", rep " << rep;
            }
        }
    }
}

// ---------------------------------------------------------------------------
//  P4.3: calcKineticEnergy == 1/2 u^T M u at random u (the current KE test only
//        checks the sqrt-operator identity, never the quadratic form directly).
// ---------------------------------------------------------------------------
TEST(MassMatrixDense, KineticEnergyMatchesDenseForm) {
    Rng rng(0x4D03);
    RobotModel m = makeDenseForest(rng);
    RobotState s;
    s.allocateFull(m);
    const int n = m.nu;
    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng); // sets a random u as well
        RobotEngine::realizePosition(m, s);
        const std::vector<Real> M = buildDenseM(m, s); // restores u after probing columns

        // 1/2 u^T M u with the state's current u
        const Real* u = s.u();
        Real quad = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                quad += u[i] * M[i * n + j] * u[j];
            }
        }

        RobotEngine::realizeVelocity(m, s); // V_GB at the real u, for calcKineticEnergy
        const Real ke = RobotEngine::calcKineticEnergy(m, s);
        EXPECT_NEAR(ke, 0.5 * quad, rtest::kLoose) << "rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  P4.4: calcLogDetM == logDetSymPD(M), M built FORWARD. Strictly stronger than
//        LogDetMatchesDenseInverse, which builds M^-1 and inverts: here nothing
//        in the oracle touches the articulated-body inverse recursion.
// ---------------------------------------------------------------------------
TEST(MassMatrixDense, LogDetMMatchesDenseM) {
    Rng rng(0x4D04);
    RobotModel m = makeDenseForest(rng);
    RobotState s;
    s.allocateFull(m);
    const int n = m.nu;
    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        const std::vector<Real> M = buildDenseM(m, s);     // forward
        RobotEngine::realizeArticulatedBodyInertias(m, s); // calcLogDetM reads P
        const Real lnDetM = RobotEngine::calcLogDetM(m, s);
        const Real lnDense = robo_linalg::logDetSymPD(M.data(), n);
        EXPECT_NEAR(lnDetM, lnDense, 1e-7) << "rep " << rep;
    }
}

// ---------------------------------------------------------------------------
//  P4.5 (FAIL guard): the forward-M and inverse-M paths are now genuinely
//        independent, so a fault injected into ONE of them cannot be cancelled
//        by the other. Build the dense-M oracle from a CLEAN Mk_G, then corrupt
//        one Mk_G block in the state and rebuild the articulated inertias from
//        it. Only the engine's M^-1 / log-det now see the corruption; the oracle
//        still holds the clean M. Both MTimesMInvIsIdentity and LogDetMMatchesDenseM
//        must go red.
//
//        (Corrupting the MODEL mass and re-realizing EVERYTHING would feed the
//        same fault to both paths and they would still agree -- which is exactly
//        why this guard desynchronizes them instead.)
// ---------------------------------------------------------------------------
TEST(MassMatrixDense, CorruptingMkGBreaksBothDensePaths) {
    Rng rng(0x4D05);
    RobotModel m = makeDenseForest(rng);
    RobotState s;
    s.allocateFull(m);
    const int n = m.nu;
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    // oracle from the CLEAN spatial inertias
    const std::vector<Real> Mclean = buildDenseM(m, s);

    // sanity: clean paths agree (proves the corruption, not a pre-existing gap,
    // is what turns them red below).
    {
        RobotState sc;
        sc.allocateFull(m);
        std::copy(s.q(), s.q() + m.nq, sc.q());
        std::copy(s.u(), s.u() + m.nu, sc.u());
        RobotEngine::realizePosition(m, sc);
        RobotEngine::realizeArticulatedBodyInertias(m, sc);
        Real idErr = 0;
        for (int j = 0; j < n; ++j) {
            std::vector<Real> Mij(static_cast<std::size_t>(n), 0), back(static_cast<std::size_t>(n), 0);
            for (int i = 0; i < n; ++i) {
                Mij[static_cast<std::size_t>(i)] = Mclean[i * n + j];
            }
            RobotEngine::multiplyByMInv(m, sc, Mij.data(), back.data());
            for (int i = 0; i < n; ++i) {
                idErr = std::max(idErr, std::abs(back[static_cast<std::size_t>(i)] - (i == j ? 1.0 : 0.0)));
            }
        }
        ASSERT_LT(idErr, rtest::kLoose) << "clean paths already disagree; guard would be vacuous";
        ASSERT_NEAR(RobotEngine::calcLogDetM(m, sc), robo_linalg::logDetSymPD(Mclean.data(), n), 1e-7);
    }

    // Corrupt one Mk_G block in the STATE (engine view only), then rebuild the
    // articulated inertias so M^-1 / calcLogDetM are derived from the fault.
    SpatialInertia* Mk = s.Mk_G();
    const int badBody = 2;
    Mk[badBody] = SpatialInertia(Mk[badBody].getMass() * Real(1.7),
                                 Mk[badBody].getMassCenter(),
                                 Mk[badBody].getUnitInertia());
    RobotEngine::realizeArticulatedBodyInertias(m, s);

    // MTimesMInvIsIdentity must FAIL: M_clean^-1_engine M_clean e_j != e_j now.
    Real maxIdErr = 0;
    for (int j = 0; j < n; ++j) {
        std::vector<Real> Mij(static_cast<std::size_t>(n), 0), back(static_cast<std::size_t>(n), 0);
        for (int i = 0; i < n; ++i) {
            Mij[static_cast<std::size_t>(i)] = Mclean[i * n + j];
        }
        RobotEngine::multiplyByMInv(m, s, Mij.data(), back.data());
        for (int i = 0; i < n; ++i) {
            maxIdErr = std::max(maxIdErr, std::abs(back[static_cast<std::size_t>(i)] - (i == j ? 1.0 : 0.0)));
        }
    }
    EXPECT_GT(maxIdErr, 1e-3) << "corrupting Mk_G did NOT break M*M^-1==I -- paths are not independent";

    // LogDetMMatchesDenseM must FAIL: engine log-det now reflects the fault, the
    // dense oracle does not.
    const Real lnDetM = RobotEngine::calcLogDetM(m, s);
    const Real lnDense = robo_linalg::logDetSymPD(Mclean.data(), n);
    EXPECT_GT(std::abs(lnDetM - lnDense), 1e-3)
        << "corrupting Mk_G did NOT break the log-det match -- paths are not independent";
}

// ---------------------------------------------------------------------------
//  P5 (docs/specs/singular-dof-fixman.md S8 INVARIANT "path consistency", CC1/
//  C6): calcLogDetM's pseudo-determinant MUST exclude a structural phantom's
//  null hinge-inertia direction -- the entire point of the pseudoLogDet fix --
//  rather than silently reverting to the pre-fix "1e-300 clamp keeps it as a
//  finite residue" behavior. Nothing else in the suite pins this down: the
//  molecule oracle (TestRoboticsOracleMolecule.cpp) `continue`s past the
//  logDetM comparison whenever the tree is singular, and the
//  Cyclic1APQPhantomLogDetIsRunConstant LEMMA computes per-body ln(D_b)
//  directly, never through calcLogDetM.
//
//  Construction (spec S4 derivation, literal case): a leaf single-atom Torsion
//  body whose ONE atom sits exactly AT the body origin, on its own rotation
//  axis (c=0 => P_b = [[0,0],[0,0]] block-wise: UnitInertia(0,0,0) about the
//  origin has no off-axis mass to begin with) -- so D_b = ~H_b P_b H_b is
//  bit-exact 0.0, not merely tiny, no floating-point argument required. Paired
//  with a normal, well-conditioned Ball-jointed body so the tree is
//  non-trivial (a lone phantom would make the WHOLE tree singular, the
//  ill-posed case the molecule oracle explicitly guards against).
//
//  Oracle: calcLogDetM on the two-body model MUST equal calcLogDetM on a
//  SEPARATE one-body model containing only the well-conditioned body (built
//  with identical q) -- by definition, that reduced-model value IS
//  "sum over non-null bodies of ln det D_b", since the well-conditioned body
//  is the tree's only non-null body either way. This is a genuinely
//  independent comparison (different model topology, not comparing
//  calcLogDetM to itself under the same inputs) and it is DISCRIMINATING: if
//  pseudoLogDet regressed to the pre-fix floor-to-finite convention (simulated
//  below via the test-local, deliberately-unfixed robo_linalg::logDetSymPD,
//  tests/RobotLinearAlgebra.hpp S1.B), the two-body calcLogDetM would pick up
//  an extra ln(clamped D_phantom) ~= ln(1e-300) ~= -690.78 term that the
//  one-body reference never sees -- an ~O(690) swing, several orders past any
//  plausible tolerance, so a reversion cannot slip through as noise.
// ---------------------------------------------------------------------------
TEST(MassMatrix, LogDetMExcludesStructuralPhantomNullDirection) {
    // Body 1: well-conditioned Ball joint (3 dof), ordinary mass properties.
    BodySpec wellConditioned;
    wellConditioned.parent = 0;
    wellConditioned.joint = JointType::Ball;
    wellConditioned.X_PF = Transform(Rotation(Real(0.4), robo::XAxis), Vec3(0.2, -0.1, 0.05));
    wellConditioned.X_BM = Transform(Rotation(Real(-0.3), robo::YAxis), Vec3(0.05, 0.1, -0.05));
    wellConditioned.mass = Real(1.5);
    wellConditioned.com_B = Vec3(0.1, -0.05, 0.2);
    wellConditioned.inertia_B = UnitInertia(Real(0.4), Real(0.5), Real(0.6));

    // Body 2: structural phantom -- leaf single-atom Torsion, atom AT the body
    // origin (on the torsion axis by construction, spec S4). X_BM is identity
    // so the torsion axis (M frame z, H_FM's constant angular part) is exactly
    // the body's own z axis. X_PF carries a non-trivial ancestor rotation on
    // purpose (proves the D_b==0 result does not depend on R_GF, spec S4/
    // realizeArticulatedBodyInertias S3 comment).
    BodySpec phantom;
    phantom.parent = 0;
    phantom.joint = JointType::Torsion;
    phantom.X_PF = Transform(Rotation(Real(0.9), robo::XAxis), Vec3(0.3, -0.2, 0.1));
    phantom.X_BM = Transform();      // identity: torsion axis == body z axis
    phantom.mass = Real(12.0);
    phantom.com_B = Vec3(0, 0, 0);   // atom AT the body origin
    phantom.inertia_B = UnitInertia(Real(0), Real(0), Real(0)); // point mass at origin: exactly zero

    RobotModel mFull = buildForest({wellConditioned, phantom});
    rtest::attachAtoms(mFull, /*body=*/2, {Vec3(0, 0, 0)}, {Real(12.0)});
    RobotState sFull;
    sFull.allocateFull(mFull);
    std::fill(sFull.q(), sFull.q() + mFull.nq, Real(0));
    sFull.q()[0] = Real(1); // body 1's quaternion (qOff==0): identity (w,x,y,z)=(1,0,0,0)
    sFull.q()[4] = Real(0.37); // body 2's own torsion angle -- arbitrary, D_b is q-independent
    std::fill(sFull.u(), sFull.u() + mFull.nu, Real(0));
    prep(mFull, sFull);

    // Independent reduced model: the well-conditioned body ALONE, identical q.
    RobotModel mReduced = buildForest({wellConditioned});
    RobotState sReduced;
    sReduced.allocateFull(mReduced);
    std::fill(sReduced.q(), sReduced.q() + mReduced.nq, Real(0));
    sReduced.q()[0] = Real(1);
    std::fill(sReduced.u(), sReduced.u() + mReduced.nu, Real(0));
    prep(mReduced, sReduced);

    const Real lnDetFull = RobotEngine::calcLogDetM(mFull, sFull);
    const Real lnDetNonNullOnly = RobotEngine::calcLogDetM(mReduced, sReduced);

    ASSERT_TRUE(std::isfinite(lnDetFull)) << "calcLogDetM must stay finite in the presence of a structural phantom";

    // Sanity/precondition: the phantom's own D_b really is null (structurally,
    // not just "small at this q") -- otherwise this test would be vacuous.
    const int uOffPhantom = mFull.bodyUIndex[2];
    const SpatialVec HbPhantom = sFull.H()[uOffPhantom];
    const Real DbPhantom = dot(HbPhantom, sFull.P()[2] * HbPhantom);
    ASSERT_EQ(DbPhantom, Real(0)) << "atom-at-origin construction must give a bit-exact-null D_b (spec S4)";

    // THE assertion: calcLogDetM on the full (phantom-containing) model equals
    // the sum over non-null bodies only -- the phantom contributes 0, not
    // ln(residue).
    EXPECT_NEAR(lnDetFull, lnDetNonNullOnly, 1e-9)
        << "calcLogDetM must EXCLUDE the structural phantom's null direction (CC1/C6)";

    // Discrimination: what the PRE-FIX floor-to-finite convention would have
    // produced for the phantom's contribution (test-local reference kernel,
    // tests/RobotLinearAlgebra.hpp, deliberately NOT updated -- S1.B).
    const Real buggyPhantomTerm = robo_linalg::logDetSymPD(&DbPhantom, 1);
    const Real lnDetIfBuggy = lnDetNonNullOnly + buggyPhantomTerm;
    const Real swing = std::abs(lnDetFull - lnDetIfBuggy);
    std::cout << "[ MassMatrix ] LogDetMExcludesStructuralPhantomNullDirection: lnDetFull=" << lnDetFull
              << " lnDetNonNullOnly=" << lnDetNonNullOnly << " buggyPhantomTerm=" << buggyPhantomTerm
              << " lnDetIfBuggy=" << lnDetIfBuggy << " swing=" << swing << std::endl;
    EXPECT_GT(swing, 50.0) << "reverting calcLogDetM to the old 1e-300 floor must be clearly discriminated"
                               " (expected swing ~ -ln(1e-300) ~= 690.78 here)";
}
