// ============================================================================
//  TestQuaternion.cpp -- the quaternion surface of robo::Rotation / robo::Quat,
//  merged from the former TestQuaternionMaps.cpp and TestQuaternionRotation.cpp
//  (they exercised the same value types and shared fixtures; one file removes the
//  duplication while keeping the two distinct concerns as separate test suites).
//
//   * suite QuaternionMaps     -- the kinematic maps N(q), Ndot, NInv: the
//     parent-frame qdot = N(q) w map (Simbody calcUnnormalizedNForQuaternion),
//     tangency, round-trip, the C-4 parent-vs-body-frame regression guards.
//   * suite QuaternionRotation -- quaternion <-> rotation conversion:
//     fromQuaternion standard-map entries, exhaustive round-trips, the two
//     converters agreeing, all four Shepperd branches, normalization, and the
//     documented divergence that Rotation(Mat33) does NOT re-orthonormalize.
//
//  *** CONFLICT C-4 (parent-frame map) ***  The engine's generalized angular speed
//  w_FM is expressed in the parent frame F and Rotation::fromQuaternion builds the
//  standard R_FM, so the qdot map MUST be the PARENT-frame map
//      qdot_x = 1/2(qw*wx + qz*wy - qy*wz)   (left Hamilton product).
//  The body-frame map (off-diagonal signs flipped) pumped kinetic energy through
//  Free roots; these tests ENFORCE the parent-frame map and guard against regression.
// ============================================================================
#include <cmath>

#include "TestHelpers.hpp"
#include "engine_helpers.hpp"
#include "robot_math.hpp"

using namespace rtest;
using robo::crossMat;
using robo::Mat33;
using robo::Quat;
using robo::Rotation;
using robo::Vec3;
using robo::Vec4;

namespace {

// Build the 4x3 matrix N(q) by probing the engine map on the angular basis.
void buildN(const Vec4& q, Real N[4][3]) {
    for (int axis = 0; axis < 3; ++axis) {
        Vec3 e(0, 0, 0);
        e[axis] = 1;
        const Vec4 col = Rotation::convertAngVelToQuaternionDot(q, e);
        for (int r = 0; r < 4; ++r) {
            N[r][axis] = col[r];
        }
    }
}

// THEORY.md S3.4 parent-frame N(q), 4x3, rows {w,x,y,z}, cols {wx,wy,wz}:
//   N = 1/2 [ -qx -qy -qz ; qw qz -qy ; -qz qw qx ; qy -qx qw ]
void theoryParentN(const Vec4& q, Real N[4][3]) {
    const Real qw = q[0], qx = q[1], qy = q[2], qz = q[3];
    const Real m[4][3] = {{-qx, -qy, -qz}, {qw, qz, -qy}, {-qz, qw, qx}, {qy, -qx, qw}};
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 3; ++c) {
            N[r][c] = 0.5 * m[r][c];
        }
    }
}

// THEORY.md S3.4 NInv(q), 3x4, rows {wx,wy,wz}, cols {w,x,y,z}:
//   NInv = 2 [ -qx qw -qz qy ; -qy qz qw -qx ; -qz -qy qx qw ]
void theoryParentNInv(const Vec4& q, Real Ninv[3][4]) {
    const Real qw = q[0], qx = q[1], qy = q[2], qz = q[3];
    const Real m[3][4] = {{-qx, qw, -qz, qy}, {-qy, qz, qw, -qx}, {-qz, -qy, qx, qw}};
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 4; ++c) {
            Ninv[r][c] = 2.0 * m[r][c];
        }
    }
}

} // namespace

// ===========================================================================
//  suite QuaternionMaps -- the kinematic maps N, Ndot, NInv
// ===========================================================================

// --- Property 1: exact component signs == PARENT-FRAME (THEORY S3.4) ---------
TEST(QuaternionMaps, N_ComponentSigns_ParentFrame) {
    Rng rng(11);
    for (int t = 0; t < 200; ++t) {
        const Quat q = rng.unitQuat();
        const Real qw = q.elems[0], qx = q.elems[1], qy = q.elems[2], qz = q.elems[3];
        const Vec3 w = rng.vec3(-2, 2);
        const Vec4 got = Rotation::convertAngVelToQuaternionDot(Vec4(qw, qx, qy, qz), w);

        EXPECT_NEAR(got[0], 0.5 * (-qx * w[0] - qy * w[1] - qz * w[2]), kAlg);
        EXPECT_NEAR(got[1], 0.5 * (qw * w[0] + qz * w[1] - qy * w[2]), kAlg);  // +qz*wy -qy*wz
        EXPECT_NEAR(got[2], 0.5 * (-qz * w[0] + qw * w[1] + qx * w[2]), kAlg); // -qz*wx +qx*wz
        EXPECT_NEAR(got[3], 0.5 * (qy * w[0] - qx * w[1] + qw * w[2]), kAlg);  // +qy*wx -qx*wy
    }
}

// --- Property 1b: full 4x3 N matrix equals the theory parent-frame matrix ----
TEST(QuaternionMaps, N_MatrixMatchesTheoryParentFrame) {
    Rng rng(21);
    for (int t = 0; t < 100; ++t) {
        const Quat q = rng.unitQuat();
        const Vec4 q4(q.elems[0], q.elems[1], q.elems[2], q.elems[3]);
        Real Ncode[4][3], Nparent[4][3];
        buildN(q4, Ncode);
        theoryParentN(q4, Nparent);
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 3; ++c) {
                EXPECT_NEAR(Ncode[r][c], Nparent[r][c], kTight) << "N[" << r << "][" << c << "]";
            }
        }
    }
}

// --- Regression guard: N must NOT be the body-frame map ----------------------
TEST(QuaternionMaps, N_IsNotBodyFrame_RegressionGuard) {
    const Quat q = Rng(22).unitQuat();
    const Real qw = q.elems[0], qx = q.elems[1], qy = q.elems[2], qz = q.elems[3];
    const Vec3 w(0.7, -1.3, 0.4);
    const Vec4 got = Rotation::convertAngVelToQuaternionDot(Vec4(qw, qx, qy, qz), w);
    const Real bodyX = 0.5 * (qw * w[0] - qz * w[1] + qy * w[2]); // wrong (body-frame) x-component
    EXPECT_GT(std::abs(got[1] - bodyX), 1e-6)
        << "N matches the body-frame map -- C-4 has regressed (flip off-diagonal signs back)";
}

// --- Property 2: tangency q . qdot = 0 --------------------------------------
TEST(QuaternionMaps, N_Tangency) {
    Rng rng(7);
    for (int t = 0; t < 200; ++t) {
        const Quat q = rng.unitQuat();
        const Vec3 w = rng.vec3(-3, 3);
        const Vec4 qd =
            Rotation::convertAngVelToQuaternionDot(Vec4(q.elems[0], q.elems[1], q.elems[2], q.elems[3]), w);
        const Real dotp = q.elems[0] * qd[0] + q.elems[1] * qd[1] + q.elems[2] * qd[2] + q.elems[3] * qd[3];
        EXPECT_NEAR(dotp, 0.0, kAlg) << "qdot must be tangent to S^3";
    }
}

// --- Property 3: round-trip NInv(N w) = w (isometry) ------------------------
TEST(QuaternionMaps, N_RoundTripRecoversAngVel) {
    Rng rng(3);
    for (int t = 0; t < 200; ++t) {
        const Quat q = rng.unitQuat();
        const Vec3 w = rng.vec3(-3, 3);
        const Vec4 qd =
            Rotation::convertAngVelToQuaternionDot(Vec4(q.elems[0], q.elems[1], q.elems[2], q.elems[3]), w);
        Real N[4][3];
        buildN(Vec4(q.elems[0], q.elems[1], q.elems[2], q.elems[3]), N);
        Vec3 wRec(0, 0, 0);
        for (int a = 0; a < 3; ++a) {
            Real s = 0;
            for (int r = 0; r < 4; ++r) {
                s += N[r][a] * qd[r];
            }
            wRec[a] = 4 * s;
        }
        EXPECT_TRUE(NearVec3(wRec, w, kAlg));
    }
}

// --- Property 4: the convention is PARENT-frame (physical rotation rate) -----
TEST(QuaternionMaps, N_ConventionIsParentFrame) {
    Rng rng(5);
    int parentWins = 0, bodyWins = 0;
    const Real eps = 1e-7;
    for (int t = 0; t < 100; ++t) {
        const Quat q = rng.unitQuat();
        const Vec3 w = rng.vec3(-2, 2);
        const Vec4 q4(q.elems[0], q.elems[1], q.elems[2], q.elems[3]);
        const Vec4 qd = Rotation::convertAngVelToQuaternionDot(q4, w);

        Quat qStep(q4[0] + eps * qd[0], q4[1] + eps * qd[1], q4[2] + eps * qd[2], q4[3] + eps * qd[3]);
        qStep.normalize();
        const Mat33 R0 = Rotation::fromQuaternion(q);
        const Mat33 R1 = Rotation::fromQuaternion(qStep);

        Mat33 dR;
        for (int i = 0; i < 9; ++i) {
            dR.elems[static_cast<std::size_t>(i)] =
                (R1.elems[static_cast<std::size_t>(i)] - R0.elems[static_cast<std::size_t>(i)]) / eps;
        }

        const Mat33 parent = crossMat(w) * R0; // [w]x R
        const Mat33 body = R0 * crossMat(w);   // R [w]x
        Real dp = 0, db = 0;
        for (int i = 0; i < 9; ++i) {
            dp = std::max(
                dp,
                std::abs(dR.elems[static_cast<std::size_t>(i)] - parent.elems[static_cast<std::size_t>(i)]));
            db = std::max(
                db,
                std::abs(dR.elems[static_cast<std::size_t>(i)] - body.elems[static_cast<std::size_t>(i)]));
        }
        if (dp < db) {
            ++parentWins;
        } else {
            ++bodyWins;
        }
    }
    std::cerr << "[C-4 regression check] parent-frame=" << parentWins << " body-frame=" << bodyWins
              << " (of 100; expect parent=100)\n";
    EXPECT_EQ(parentWins, 100)
        << "orientation must advance as Rdot=[w]xR (parent frame); body-frame => C-4 regression";
    EXPECT_EQ(bodyWins, 0);
}

// --- Property 5: Ndot consistency via finite difference ---------------------
TEST(QuaternionMaps, Ndot_FiniteDifferenceConsistent) {
    Rng rng(9);
    for (int t = 0; t < 50; ++t) {
        const Quat q0 = rng.unitQuat();
        const Vec3 w = rng.vec3(-2, 2);
        const Vec3 wdot = rng.vec3(-2, 2);
        const Vec4 q4(q0.elems[0], q0.elems[1], q0.elems[2], q0.elems[3]);

        const Vec4 qddot = Rotation::convertAngVelDotToQuaternionDotDot(q4, w, wdot);

        auto qdotAt = [&](Real s) -> Vec4 {
            const Vec4 qd0 = Rotation::convertAngVelToQuaternionDot(q4, w);
            Quat qs(q4[0] + s * qd0[0], q4[1] + s * qd0[1], q4[2] + s * qd0[2], q4[3] + s * qd0[3]);
            qs.normalize();
            const Vec3 wAt(w[0] + s * wdot[0], w[1] + s * wdot[1], w[2] + s * wdot[2]);
            return Rotation::convertAngVelToQuaternionDot(
                Vec4(qs.elems[0], qs.elems[1], qs.elems[2], qs.elems[3]),
                wAt);
        };
        const Real h = 1e-6;
        const Vec4 fp = qdotAt(h), fm = qdotAt(-h);
        for (int r = 0; r < 4; ++r) {
            const Real fd = (fp[r] - fm[r]) / (2 * h);
            EXPECT_NEAR(qddot[r], fd, 1e-5) << "component " << r;
        }
    }
}

// --- Property 6: the two in-code implementations agree (single source) -------
TEST(QuaternionMaps, TwoImplementationsAgree) {
    Rng rng(13);
    for (int t = 0; t < 200; ++t) {
        const Quat q = rng.unitQuat();
        const Vec3 w = rng.vec3(-3, 3);
        const Quat a = Quat::angVelToQdot(q, w);
        const Vec4 b =
            Rotation::convertAngVelToQuaternionDot(Vec4(q.elems[0], q.elems[1], q.elems[2], q.elems[3]), w);
        for (int r = 0; r < 4; ++r) {
            EXPECT_NEAR(a.elems[static_cast<std::size_t>(r)], b[r], kTight) << "component " << r;
        }
    }
}

// --- Reference: theory parent-frame NInv . N == I_3 (and NInv = 4 N^T) -------
TEST(QuaternionMaps, NInv_TheoryParentFrameIsMutualInverse) {
    Rng rng(23);
    for (int t = 0; t < 100; ++t) {
        const Quat q = rng.unitQuat();
        const Vec4 q4(q.elems[0], q.elems[1], q.elems[2], q.elems[3]);
        Real N[4][3], Ninv[3][4];
        theoryParentN(q4, N);
        theoryParentNInv(q4, Ninv);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Real acc = 0;
                for (int k = 0; k < 4; ++k) {
                    acc += Ninv[i][k] * N[k][j];
                }
                EXPECT_NEAR(acc, (i == j) ? 1.0 : 0.0, kAlg) << "(" << i << "," << j << ")";
            }
        }
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < 4; ++k) {
                EXPECT_NEAR(Ninv[i][k], 4.0 * N[k][i], kTight);
            }
        }
    }
}

// ===========================================================================
//  suite QuaternionRotation -- quaternion <-> rotation conversion
// ===========================================================================

// --- fromQuaternion is the STANDARD R_FM map (exact entries) -----------------
TEST(QuaternionRotation, FromQuaternionStandardMap90Z) {
    const Real c = std::cos(M_PI / 4), s = std::sin(M_PI / 4);
    Quat q(c, 0, 0, s);
    const Mat33 R = Rotation::fromQuaternion(q);
    const Mat33 ref(0, -1, 0, 1, 0, 0, 0, 0, 1); // Rz(90): x->y, y->-x
    EXPECT_TRUE(NearMat33(R, ref, kAlg));
}

TEST(QuaternionRotation, FromQuaternion90AboutEachAxis) {
    const Real c = std::cos(M_PI / 4), s = std::sin(M_PI / 4);
    EXPECT_TRUE(NearMat33(Rotation::fromQuaternion(Quat(c, s, 0, 0)),
                          Mat33(1, 0, 0, 0, 0, -1, 0, 1, 0),
                          kAlg)); // Rx(90)
    EXPECT_TRUE(NearMat33(Rotation::fromQuaternion(Quat(c, 0, s, 0)),
                          Mat33(0, 0, 1, 0, 1, 0, -1, 0, 0),
                          kAlg)); // Ry(90)
    EXPECT_TRUE(NearMat33(Rotation::fromQuaternion(Quat(c, 0, 0, s)),
                          Mat33(0, -1, 0, 1, 0, 0, 0, 0, 1),
                          kAlg)); // Rz(90)
}

// --- exhaustive round-trips (quat -> R -> quat and R -> quat -> R) -----------
TEST(QuaternionRotation, RoundTripQuatRQuat) {
    Rng rng(1);
    for (int t = 0; t < 1000; ++t) {
        const Quat q = rng.unitQuat();
        const Mat33 R = Rotation::fromQuaternion(q);
        Real w, x, y, z;
        EngineHelpers::rotationToQuaternion(Rotation(R), w, x, y, z);
        EXPECT_LT(quatDist(q, Quat(w, x, y, z)), kAlg);
    }
}

TEST(QuaternionRotation, RoundTripRQuatR) {
    Rng rng(2);
    for (int t = 0; t < 1000; ++t) {
        const Mat33 R = rng.rotation();
        Real w, x, y, z;
        EngineHelpers::rotationToQuaternion(Rotation(R), w, x, y, z);
        const Mat33 R2 = Rotation::fromQuaternion(Quat(w, x, y, z));
        EXPECT_TRUE(NearMat33(R, R2, kAlg));
    }
}

// --- two converters agree: fromQuaternion vs quatToRotation (World.cpp) ------
TEST(QuaternionRotation, TwoConvertersAgree) {
    Rng rng(3);
    for (int t = 0; t < 500; ++t) {
        const Quat q = rng.unitQuat();
        const Mat33 a = Rotation::fromQuaternion(q);
        const Mat33 b = EngineHelpers::quatToRotation(q.elems[0], q.elems[1], q.elems[2], q.elems[3]);
        EXPECT_TRUE(NearMat33(a, b, kTight));
    }
}

// --- all four Shepperd branches (trace>0 and each of the 3 diagonal pivots) --
TEST(QuaternionRotation, ShepperdAllBranches) {
    const Real c = std::cos(M_PI / 2 * 0.99), s = std::sin(M_PI / 2 * 0.99); // ~178 deg
    const Quat cases[4] = {
        Quat(1, 0, 0, 0), // identity -> trace>0
        Quat(c, s, 0, 0), // ~180 about x -> m00 pivot
        Quat(c, 0, s, 0), // ~180 about y -> m11 pivot
        Quat(c, 0, 0, s), // ~180 about z -> m22 pivot
    };
    for (const Quat& q : cases) {
        const Mat33 R = Rotation::fromQuaternion(q);
        Real w, x, y, z;
        EngineHelpers::rotationToQuaternion(Rotation(R), w, x, y, z);
        EXPECT_LT(quatDist(q, Quat(w, x, y, z)), kAlg) << "Shepperd branch round-trip failed";
        EXPECT_TRUE(NearMat33(R, Rotation::fromQuaternion(Quat(w, x, y, z)), kAlg));
    }
}

// --- normalization semantics ------------------------------------------------
TEST(QuaternionRotation, NormalizeOnConstructionAndZeroGuard) {
    Quat q(Vec4(1, 2, 3, 4)); // ctor normalizes
    EXPECT_NEAR(q.norm(), 1.0, kTight);

    Quat zero(0, 0, 0, 0);
    zero.normalize(); // must map to identity, not divide by zero
    EXPECT_NEAR(zero.elems[0], 1.0, kTight);
    EXPECT_NEAR(zero.elems[1], 0.0, kTight);
    EXPECT_NEAR(zero.elems[2], 0.0, kTight);
    EXPECT_NEAR(zero.elems[3], 0.0, kTight);

    Quat tiny(1e-20, 0, 0, 0);
    tiny.normalize();
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(tiny.elems[static_cast<std::size_t>(i)]));
    }
}

// --- Rotation(Mat33) is a STRAIGHT COPY -- it does NOT re-orthonormalize. -----
// DIVERGENCE FROM SIMBODY: Simbody's Rotation(Mat33) builds the nearest orthonormal
// rotation; the current port copies verbatim. This pins that behavior so a future
// change is noticed; callers must feed already-orthonormal matrices.
TEST(QuaternionRotation, RotationFromMat33DoesNotReorthonormalize_DivergenceFromSimbody) {
    const Mat33 dirty(1.0, 0.02, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const Rotation R(dirty);
    for (int i = 0; i < 9; ++i) {
        EXPECT_EQ(R.elems[static_cast<std::size_t>(i)], dirty.elems[static_cast<std::size_t>(i)]);
    }
    EXPECT_FALSE(NearMat33(R * R.transpose(), Mat33::identity(), 1e-6))
        << "ctor unexpectedly re-orthonormalized; if a nearest-rotation ctor was "
           "added (matching Simbody), update this test and the frame-build assumptions";
}
