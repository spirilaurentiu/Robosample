// ============================================================================
//  RobotEngine.cpp - SimTK-free articulated-body dynamics.
//
//  This is a faithful port of the Simbody recursions (no MobilizedBody /
//  SimbodyMatterSubsystem / State / vtables). The SoA driver walks the body
//  tree (model.bodyParent[b] < b: outward = forward, inward = reverse) and
//  branches on model.bodyJoint[b] instead of virtual dispatch. The spatial
//  *value* algebra (SpatialVec, ArticulatedInertia, PhiMatrix, Transform,
//  SpatialInertia) is SimTK's own header-only, vtable-free code, so only the
//  recursion is reimplemented.
//
//  Source map (Simbody01/Simbody/src):
//    realizePosition / realizeVelocity ......... RigidBodyNodeSpec.h
//    calcBodyTransforms / H_PB_G ............... RigidBodyNodeSpec.{h,cpp}
//    Phi / gyroscopic / coriolis .............. RigidBodyNode.cpp
//    realizeArticulatedBodyInertias ........... RigidBodyNodeSpec.cpp
//    calcUDot pass1/2 ......................... RigidBodyNodeSpec.cpp
//    multiplyByMInv pass1/2 ................... RigidBodyNodeSpec.cpp
//    multiplyBySqrtMInvPassOutward ............ RigidBodyNodeSpec.cpp
//    per-joint H_FM / X_FM / qdot ............. RigidBodyNodeSpec_{Pin,Translation,Free}.h
//
//  VALIDATION GATE: every operator below must be diffed against the still-
//  present Simbody build (random q,u; all mobilizer types) to tolerance before
//  Simbody is removed. See validate_robotics notes.
// ============================================================================

#include "RobotEngine.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "Constraints.hpp"

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

namespace {

// log|det A| for a small SYMMETRIC POSITIVE-DEFINITE n x n matrix (row-major)
// via Cholesky: A = L L^T  =>  det A = prod(L_ii)^2  =>  log det = 2 sum log L_ii.
// D_b = ~H P H is SPD for any physical articulated body, so Cholesky is the
// right (cheapest, most stable) factorization here. Returns 0 for n == 0.
inline robo::Real logDetSymPD(const robo::Real* A, int n) {
    if (n <= 0) {
        return robo::Real(0);
    }
    robo::Real L[36]; // n <= 6
    for (int i = 0; i < n * n; ++i) {
        L[i] = robo::Real(0);
    }
    robo::Real logDet = robo::Real(0);
    for (int j = 0; j < n; ++j) {
        robo::Real sum = A[j * n + j];
        for (int k = 0; k < j; ++k) {
            sum -= L[j * n + k] * L[j * n + k];
        }
        // Guard against a non-finite / non-PD block (near-singular hinge): clamp
        // the pivot so the caller gets a large-but-finite Fixman rather than NaN,
        // which the Metropolis step can then reject cleanly.
        if (!(sum > robo::Real(1e-300))) {
            sum = robo::Real(1e-300);
        }
        const robo::Real Ljj = std::sqrt(sum);
        L[j * n + j] = Ljj;
        logDet += robo::Real(2) * std::log(Ljj);
        for (int i = j + 1; i < n; ++i) {
            robo::Real s = A[i * n + j];
            for (int k = 0; k < j; ++k) {
                s -= L[i * n + k] * L[j * n + k];
            }
            L[i * n + j] = s / Ljj;
        }
    }
    return logDet;
}


// ============================================================================
//  DEBUG LOGGING
//  ROBO_DEBUG   : master switch (0 disables every probe with zero overhead).
//  ROBO_VERBOSE : 0 = only the first-NaN dump + the D/DI dump at the failing
//                     body (recommended: pinpoints the origin, low spam);
//                 1 = + a one-line per-step growth summary;
//                 2 = + per-body X_GB/X_FM every step (huge; short runs only).
//  Override at compile time, e.g.  -DROBO_VERBOSE=1 .
// ============================================================================
#define ROBO_DEBUG 1
#define ROBO_VERBOSE 2

#ifndef ROBO_DEBUG
#    define ROBO_DEBUG 1
#endif
#ifndef ROBO_VERBOSE
#    define ROBO_VERBOSE 1
#endif

#if ROBO_DEBUG
namespace robodbg {

long evalCount = 0; // bumped once per full derivative evaluation
long stepCount = 0; // bumped once per Verlet step
bool firstNanDumped = false;

inline bool fin(Real x) {
    return std::isfinite(static_cast<double>(x));
}
inline bool finV(const Vec3& v) {
    return fin(v[0]) && fin(v[1]) && fin(v[2]);
}
inline bool finS(const SpatialVec& s) {
    return finV(s[0]) && finV(s[1]);
}

// Dump the q/u/position context of the body that first went non-finite.
void dumpBody(const RobotModel& m, RobotState& s, int b) {
    if (b < 1 || b >= m.numBodies) {
        return;
    }
    std::cout << "  body " << b << ": joint=" << static_cast<int>(m.bodyJoint[b])
              << " parent=" << m.bodyParent[b] << " qOff=" << m.bodyQIndex[b] << " uOff=" << m.bodyUIndex[b]
              << " dof=" << m.bodyNU[b] << "\n";
    const Real* q = s.q();
    const Real* u = s.u();
    const Real* ud = s.udot();
    std::cout << "    q =";
    for (int k = 0; k < m.bodyNQ[b]; ++k) {
        std::cout << ' ' << q[m.bodyQIndex[b] + k];
    }
    std::cout << "\n    u =";
    for (int k = 0; k < m.bodyNU[b]; ++k) {
        std::cout << ' ' << u[m.bodyUIndex[b] + k];
    }
    std::cout << "\n    udot =";
    for (int k = 0; k < m.bodyNU[b]; ++k) {
        std::cout << ' ' << ud[m.bodyUIndex[b] + k];
    }
    std::cout << "\n    X_GB.p = " << s.X_GB()[b].p() << "\n";
}

// Scan the core state for the first non-finite value after a given phase.
// On the very first detection in the whole run, dump rich context so the
// origin (which routine, which body) is unambiguous. Returns true if any
// non-finite value is present.
bool scan(const RobotModel& m, RobotState& s, const char* where) {
    int badBody = -1;
    const char* badArr = nullptr;
    int badIdx = -1;

    const int nq = m.nq, nu = m.nu, nb = m.numBodies, na = m.numAtoms;
    const Real* q = s.q();
    const Real* u = s.u();
    const Real* ud = s.udot();
    const Real* qd = s.qdot();
    const Real* qdd = s.qdotdot();
    const Transform* X = s.X_GB();
    const Vec3* pos = s.atomPosG();
    const SpatialVec* bf = s.bodyForceG();

    for (int i = 0; i < nq && !badArr; ++i) {
        if (!fin(q[i])) {
            badArr = "q";
            badIdx = i;
        }
    }
    for (int i = 0; i < nu && !badArr; ++i) {
        if (!fin(u[i])) {
            badArr = "u";
            badIdx = i;
        }
    }
    for (int i = 0; i < nu && !badArr; ++i) {
        if (!fin(ud[i])) {
            badArr = "udot";
            badIdx = i;
        }
    }
    for (int i = 0; i < nq && !badArr; ++i) {
        if (!fin(qd[i])) {
            badArr = "qdot";
            badIdx = i;
        }
    }
    for (int i = 0; i < nq && !badArr; ++i) {
        if (!fin(qdd[i])) {
            badArr = "qdotdot";
            badIdx = i;
        }
    }
    for (int b = 1; b < nb && !badArr; ++b) {
        if (!finV(X[b].p())) {
            badArr = "X_GB.p";
            badIdx = b;
            badBody = b;
        }
    }
    for (int a = 0; a < na && !badArr; ++a) {
        if (!finV(pos[a])) {
            badArr = "atomPosG";
            badIdx = a;
            badBody = m.atomBody[a];
        }
    }
    for (int b = 1; b < nb && !badArr; ++b) {
        if (!finS(bf[b])) {
            badArr = "bodyForceG";
            badIdx = b;
            badBody = b;
        }
    }

    if (!badArr) {
        return false;
    }

    if (firstNanDumped) {
        std::cout << "[robo] still non-finite @ " << where << " :: " << badArr << "[" << badIdx << "] (eval#"
                  << evalCount << " step#" << stepCount << ")\n"
                  << std::flush;
        return true;
    }
    firstNanDumped = true;
    std::cout << "\n================ FIRST NON-FINITE DETECTED ================\n"
              << "  phase  : " << where << "\n"
              << "  array  : " << badArr << "[" << badIdx << "]\n"
              << "  eval#  : " << evalCount << "   step# : " << stepCount << "\n";
    if (badBody < 0 && (std::string(badArr) == "u" || std::string(badArr) == "udot")) {
        // map a u-index back to its body for context
        for (int b = 1; b < nb; ++b) {
            if (badIdx >= m.bodyUIndex[b] && badIdx < m.bodyUIndex[b] + m.bodyNU[b]) {
                badBody = b;
                break;
            }
        }
    }
    if (badBody < 0
        && (std::string(badArr) == "q" || std::string(badArr) == "qdot"
            || std::string(badArr) == "qdotdot")) {
        for (int b = 1; b < nb; ++b) {
            if (badIdx >= m.bodyQIndex[b] && badIdx < m.bodyQIndex[b] + m.bodyNQ[b]) {
                badBody = b;
                break;
            }
        }
    }
    dumpBody(m, s, badBody);
    std::cout << "===========================================================\n" << std::flush;
    return true;
}

} // namespace robodbg
#    define ROBO_CHECK(where) (void)robodbg::scan(m, s, (where))
#else
#    define ROBO_CHECK(where) ((void)0)
#endif

// Spatial dot product:  ~a * b  ==  a.angular . b.angular + a.linear . b.linear.
inline Real spatialDot(const SpatialVec& a, const SpatialVec& b) {
    return dot(a[0], b[0]) + dot(a[1], b[1]);
}

// Dense inverse of an n x n block (n in {1,3,6}), row-major in/out.
// n == 1 is closed form (the Pin-joint hot path); n >= 2 goes through LAPACK
// LU (dgetrf/dgetri). Always returns a finite result (singular -> regularized).
static void jacobiSymEig(const Real* Ain, int n, Real* d, Real* V); // defined below

// Symmetric (Moore-Penrose) pseudo-inverse of the dof x dof articulated hinge
// matrix D = ~H P H, row-major in/out. D is symmetric and, for a well-posed
// hinge, SPD; this routine returns its inverse.
//
// CRITICAL ROBUSTNESS PROPERTY: a *degenerate* hinge -- a DOF with ~zero
// articulated inertia -- must be LOCKED, not regularized into a huge finite
// inverse. The canonical case is a single atom sitting exactly on its own
// torsion axis (D ~ 1e-33). It arises naturally once a ring-closing bond is
// removed and a ring atom (e.g. a proline/ring nitrogen) is left as a lone
// body whose root atom lies on the bond it rotates about; rotating a point
// mass on its own axis is unobservable, so the DOF has no dynamics.
//
// The previous code clamped such a D to +-1e-14 and returned 1/eps ~ 1e14.
// That inverse is enormous but FINITE, so it slips past the non-finite guard,
// multiplies the (tiny, noisy) force residual into a ~1e14 acceleration, and
// blows up on the next step -- exactly the observed failure. The correct
// action is a pseudo-inverse: invert the well-conditioned directions and set
// the inverse of any ~null direction to 0. DI=0 => G=0 => PPlus=P (full
// inertia passes through, i.e. the joint behaves as a weld) => udot=0 and the
// HMC seed velocity is 0 too. The DOF is frozen at its initial value, which is
// the physically exact treatment of a null coordinate.
void invertDense(const Real* A, int n, Real* Ainv) {
    const Real lockTol = Real(1e-12); // |eigenvalue| <= this is treated as null

    // Dominant case: 1-DOF Pin. Closed form, with a lock instead of a clamp.
    if (n == 1) {
        const Real d = A[0];
        Ainv[0] = (std::abs(d) > lockTol) ? (Real(1) / d) : Real(0);
        return;
    }

    // n in {3,6}: symmetric pseudo-inverse via the in-house Jacobi eigensolver
    // (self-contained; no LAPACK/Eigen link dependency). For a non-degenerate
    // hinge every eigenvalue is well above tol, so this equals the true inverse;
    // any null eigendirection contributes 0 (that direction is locked).
    Real d[6];
    Real V[36];
    jacobiSymEig(A, n, d, V);
    Real scale = 0;
    for (int k = 0; k < n; ++k) {
        scale = std::max(scale, std::abs(d[k]));
    }
    const Real tol = std::max(lockTol, scale * Real(1e-12));
    for (int i = 0; i < n * n; ++i) {
        Ainv[i] = Real(0);
    }
    for (int k = 0; k < n; ++k) {
        if (std::abs(d[k]) <= tol) {
            continue; // null direction -> pseudo-inverse contributes 0 (locked)
        }
        const Real inv = Real(1) / d[k];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Ainv[i * n + j] += inv * V[i * n + k] * V[j * n + k];
            }
        }
    }
}

// Symmetric eigensolver (cyclic Jacobi) for a small dense symmetric matrix
// (n <= 6 here). Row-major input A[n*n]; outputs eigenvalues d[n] and
// eigenvectors as COLUMNS of V (V[i*n+k] = component i of eigenvector k).
// Self-contained: no SimTK::Eigen / LAPACK dependency (that symbol is not in
// the link line and produced an undefined-symbol ImportError).
static void jacobiSymEig(const Real* Ain, int n, Real* d, Real* V) {
    Real a[36];
    for (int i = 0; i < n * n; ++i) {
        a[i] = Ain[i];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            V[i * n + j] = (i == j) ? Real(1) : Real(0);
        }
    }
    for (int sweep = 0; sweep < 100; ++sweep) {
        Real off = 0;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                off += a[p * n + q] * a[p * n + q];
            }
        }
        if (off <= Real(1e-30)) {
            break;
        }
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const Real apq = a[p * n + q];
                if (apq == Real(0)) {
                    continue;
                }
                const Real app = a[p * n + p];
                const Real aqq = a[q * n + q];
                const Real theta = (aqq - app) / (2 * apq);
                const Real t =
                    (theta >= 0 ? Real(1) : Real(-1)) / (std::abs(theta) + std::sqrt(theta * theta + 1));
                const Real cs = Real(1) / std::sqrt(t * t + 1);
                const Real sn = t * cs;
                // A <- J^T A J  (rotate columns then rows of the (p,q) plane).
                for (int i = 0; i < n; ++i) {
                    const Real aip = a[i * n + p];
                    const Real aiq = a[i * n + q];
                    a[i * n + p] = cs * aip - sn * aiq;
                    a[i * n + q] = sn * aip + cs * aiq;
                }
                for (int i = 0; i < n; ++i) {
                    const Real api = a[p * n + i];
                    const Real aqi = a[q * n + i];
                    a[p * n + i] = cs * api - sn * aqi;
                    a[q * n + i] = sn * api + cs * aqi;
                }
                // V <- V J  (accumulate eigenvectors as columns).
                for (int i = 0; i < n; ++i) {
                    const Real vip = V[i * n + p];
                    const Real viq = V[i * n + q];
                    V[i * n + p] = cs * vip - sn * viq;
                    V[i * n + q] = sn * vip + cs * viq;
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        d[i] = a[i * n + i];
    }
}

// Symmetric matrix square root of a dof x dof SPD block (for multiplyBySqrtMInv).
// dof==1 is a scalar sqrt; otherwise S = V diag(sqrt(lambda)) ~V, exactly as
// Simbody's multiplyBySqrtMInvPassOutward does (sqrtDI = U diag(sqrt) ~U).
void symSqrt(const Real* A, int n, Real* S) {
    if (n == 1) {
        S[0] = std::sqrt(std::max(Real(0), A[0]));
        return;
    }
    Real d[6];
    Real V[36];
    jacobiSymEig(A, n, d, V);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Real acc = 0;
            for (int k = 0; k < n; ++k) {
                acc += V[i * n + k] * std::sqrt(std::max(Real(0), d[k])) * V[j * n + k];
            }
            S[i * n + j] = acc;
        }
    }
}

// -------- per-joint cross-mobilizer transform X_FM(q) ----------------------
// Faithful to RigidBodyNodeSpec_{Pin,Translation,Free}.h / RigidBodyNode_Weld.
Transform jointX_FM(JointType jt, const Real* q, int qOff) {
    switch (jt) {
        case JointType::Weld:
            return Transform(); // identity
        case JointType::Pin: {
            Rotation R;
            R.setRotationFromAngleAboutZ(q[qOff]);
            return Transform(R, Vec3(0));
        }
        case JointType::Translation:
            return Transform(Rotation(), Vec3(q[qOff], q[qOff + 1], q[qOff + 2]));
        case JointType::Free: {
            // q = [q0 q1 q2 q3 x y z], quaternion (normalized) + translation.
            const Vec4 qv(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            Quaternion quat(qv); // normalizes
            Rotation R;
            R.setRotationFromQuaternion(quat);
            return Transform(R, Vec3(q[qOff + 4], q[qOff + 5], q[qOff + 6]));
        }
        default:
            return Transform();
    }
}

// -------- per-joint H_FM columns (in F), written into Hcol[0..dof-1] --------
void jointH_FM(JointType jt, SpatialVec* Hcol) {
    switch (jt) {
        case JointType::Weld:
            break;
        case JointType::Pin:
            Hcol[0] = SpatialVec(Vec3(0, 0, 1), Vec3(0));
            break;
        case JointType::Translation:
            Hcol[0] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            Hcol[1] = SpatialVec(Vec3(0), Vec3(0, 1, 0));
            Hcol[2] = SpatialVec(Vec3(0), Vec3(0, 0, 1));
            break;
        case JointType::Free:
            Hcol[0] = SpatialVec(Vec3(1, 0, 0), Vec3(0));
            Hcol[1] = SpatialVec(Vec3(0, 1, 0), Vec3(0));
            Hcol[2] = SpatialVec(Vec3(0, 0, 1), Vec3(0));
            Hcol[3] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            Hcol[4] = SpatialVec(Vec3(0), Vec3(0, 1, 0));
            Hcol[5] = SpatialVec(Vec3(0), Vec3(0, 0, 1));
            break;
        default:
            break;
    }
}

} // namespace

// ============================================================================
//  KINEMATICS
// ============================================================================
void RobotEngine::realizePosition(const RobotModel& m, RobotState& s) {
    Transform* X_GB = s.X_GB();
    Transform* X_FM = s.X_FM();
    Transform* X_PB = s.X_PB();
    PhiMatrix* Phi = s.Phi();
    SpatialInertia* Mk = s.Mk_G();
    Vec3* comG = s.comG();
    SpatialVec* HFM = s.H_FM();
    SpatialVec* H = s.H();
    const Real* q = s.q();

    X_GB[0] = Transform(); // Ground
    Mk[0] = SpatialInertia();

    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const JointType jt = m.bodyJoint[b];
        const int qOff = m.bodyQIndex[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];

        // X_FM, then X_PB = X_PF * X_FM * X_MB, X_GB = X_GP * X_PB.
        X_FM[b] = jointX_FM(jt, q, qOff);
        const Transform X_MB = ~m.X_BM[b];
        const Transform X_FB = X_FM[b] * X_MB;
        X_PB[b] = m.X_PF[b] * X_FB;
        X_GB[b] = X_GB[p] * X_PB[b];

        // H_FM (in F) and H_PB_G (in Ground): RigidBodyNodeSpec.cpp.
        jointH_FM(jt, &HFM[uOff]);
        const Rotation R_GF = X_GB[p].R() * m.X_PF[b].R();
        // r_MB is the vector from Mo to Bo expressed in M, i.e. X_MB.p() (Simbody
        // getX_MB().p()), NOT X_BM.p(). X_BM.p() is Mo's position in B (the
        // vector Bo->Mo expressed in B): wrong frame AND wrong sign. Using it
        // corrupts the LINEAR rows of the hinge matrix for every body whose M
        // frame is offset from Bo (all non-root torsion bodies, offset = bond
        // length), giving wrong cross-mobilizer velocities -> non-conserving
        // dynamics. The position path never touches H, so the OpenMM energy
        // check cannot catch this. X_MB == ~X_BM is already computed above.
        const Vec3 r_MB = X_MB.p();
        const Rotation& R_FM = X_FM[b].R();
        const Vec3 r_MB_F = R_FM * r_MB;
        for (int j = 0; j < dof; ++j) {
            const SpatialVec& h = HFM[uOff + j];
            // H_MB_F: top row 0, bottom row = -r_MB_F % h.angular.
            const SpatialVec hpb(h[0], h[1] - (r_MB_F % h[0]));
            H[uOff + j] = SpatialVec(R_GF * hpb[0], R_GF * hpb[1]);
        }

        // Joint-independent kinematics (Phi, COM_G, Mk_G).  RigidBodyNode.cpp.
        const Vec3 p_PB_G = X_GB[p].R() * X_PB[b].p();
        Phi[b] = PhiMatrix(p_PB_G);
        const Rotation& R_GB = X_GB[b].R();
        const Vec3& p_GB = X_GB[b].p();
        const robo::UnitInertia G_Bo_G = m.bodyUnitInertia_B[b].reexpress(~R_GB);
        const Vec3 p_BBc_G = R_GB * m.bodyCom_B[b];
        comG[b] = p_GB + p_BBc_G;
        Mk[b] = SpatialInertia(m.bodyMass[b], p_BBc_G, G_Bo_G);
    }

    // Per-atom Ground positions and stations (R_GB * station_B). The transfer
    // currency + the input to the OpenMM force bridge.
    Vec3* posG = s.atomPosG();
    Vec3* stG = s.atomStationG();
    for (int a = 0; a < m.numAtoms; ++a) {
        const int b = m.atomBody[a];
        const Vec3 st = X_GB[b].R() * m.atomStation_B[a];
        stG[a] = st;
        posG[a] = X_GB[b].p() + st;
    }
}

void RobotEngine::realizeVelocity(const RobotModel& m, RobotState& s) {
    calcQDot(m, s, s.qdot());

    const SpatialVec* H = s.H();
    const SpatialVec* HFM = s.H_FM();
    SpatialVec* V_FM = s.V_FM();
    SpatialVec* V_PB_G = s.V_PB_G();
    SpatialVec* V_GB = s.V_GB();
    SpatialVec* gyro = s.gyro();
    SpatialVec* a_tot = s.coriolisA();
    SpatialVec* a_mob = s.mobCoriolisA();
    const PhiMatrix* Phi = s.Phi();
    const SpatialInertia* Mk = s.Mk_G();
    const Transform* X_GB = s.X_GB(); // for R_GF (centripetal HDot_MB_F term)
    const Transform* X_FM = s.X_FM(); // for r_MB_F  (centripetal HDot_MB_F term)
    const Real* u = s.u();

    V_GB[0] = SpatialVec(Vec3(0), Vec3(0));
    a_tot[0] = SpatialVec(Vec3(0), Vec3(0));

    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];

        SpatialVec vfm(Vec3(0), Vec3(0));
        SpatialVec vpb(Vec3(0), Vec3(0));
        for (int j = 0; j < dof; ++j) {
            vfm += HFM[uOff + j] * u[uOff + j];
            vpb += H[uOff + j] * u[uOff + j];
        }
        V_FM[b] = vfm;
        V_PB_G[b] = vpb;

        // V_GB = ~Phi * V_GP + V_PB_G.
        const SpatialVec V_GBb = (~Phi[b]) * V_GB[p] + vpb;
        V_GB[b] = V_GBb;
        const Vec3& w_GB = V_GBb[0];
        const Vec3& v_GB = V_GBb[1];

        // Gyroscopic force b = mass * [ w x (G w) ; w x (w x c) ]  (Ground),
        // with G the unit inertia about Bo and c the COM offset (RigidBodyNode.cpp).
        const SpatialInertia& Mkb = Mk[b];
        const Vec3 Iw = Mkb.getUnitInertia() * w_GB;
        const Vec3 c = Mkb.getMassCenter();
        gyro[b] = Mkb.getMass() * SpatialVec(w_GB % Iw, w_GB % (w_GB % c));

        // Mobilizer coriolis (bias) acceleration  a_mob = d(~Phi)/dt V_GP + Hdot_G u.
        //
        // EARLIER BUG: this dropped Hdot_G u, arguing H_FM is constant in the M
        // frame. Wrong frame. H_PB_G is H_FM expressed in GROUND, and the
        // inboard frame F is fixed in the PARENT body, so R_GF (hence H_PB_G)
        // rotates with the parent angular velocity w_GP. Thus
        //   Hdot_G u = w_GP x (H_PB_G u) = w_GP x V_PB_G   (NONZERO).
        // Omitting it drops the angular Coriolis coupling entirely and half the
        // linear Coriolis, so the integrator stops conserving energy: velocities
        // grow every step from a thermal start until the trajectory explodes.
        // Full bias (equivalent to the classic centripetal + 2*Coriolis form
        //  w_GP x (w_GP x p_PB) + 2 w_GP x v_PB_G):
        //   a_mob = [ w_GP x w_PB_G ;
        //             w_GP x (v_GB - v_GP) + w_GP x v_PB_G ]
        const Vec3& w_GP = V_GB[p][0];
        const Vec3& v_GP = V_GB[p][1];
        const Vec3& w_PB_G = V_PB_G[b][0]; // cross-mobilizer angular vel in G
        const Vec3& v_PB_G = V_PB_G[b][1]; // cross-mobilizer linear  vel in G

        // MISSING-TERM FIX: the bias above is only complete when the mobilized
        // (M) frame origin coincides with the body (B) origin. Simbody's
        // calcParentToChildVelocityJacobianInGroundDot carries an extra HDot_MB_F
        // contribution whenever r_MB = X_MB.p() != 0 (true for every non-root
        // Pin/torsion body, where r_MB is the bond length). With the per-joint
        // H_FM constant in F (Pin/Free/Translation/Weld all have HDot_FM = 0),
        // that contribution reduces to the classic centripetal term
        //     R_GF * ( w_FM x (w_FM x r_MB_F) ),   r_MB_F = R_FM * r_MB,
        // and belongs in the LINEAR row of A_mob. Omitting it leaves the joint
        // bias acceleration inconsistent with H/Hdot and slowly pumps energy.
        // NOTE r_MB = X_MB.p() = (~X_BM).p(), the Mo->Bo vector in M (see the
        // matching fix in realizePosition); X_BM.p() would be wrong here too.
        const Vec3& w_FM = V_FM[b][0];                      // cross-mobilizer ang. vel in F
        const Vec3 r_MB_F = X_FM[b].R() * (~m.X_BM[b]).p(); // Mo->Bo, expressed in F
        const Rotation R_GF = X_GB[p].R() * m.X_PF[b].R();  // F orientation in Ground
        const Vec3 centripetal_G = R_GF * (w_FM % (w_FM % r_MB_F));

        const SpatialVec A_mob(w_GP % w_PB_G, w_GP % (v_GB - v_GP) + w_GP % v_PB_G + centripetal_G);
        a_mob[b] = A_mob;

        // Total coriolis accel a = ~Phi * a_parent + A_mob.
        a_tot[b] = (~Phi[b]) * a_tot[p] + A_mob;
    }
}

void RobotEngine::calcQDot(const RobotModel& m, const RobotState& s, Real* qdotOut) {
    const Real* u = s.u();
    const Transform* X_FM = s.X_FM();
    const Real* q = s.q();
    for (int b = 1; b < m.numBodies; ++b) {
        const int qOff = m.bodyQIndex[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        switch (m.bodyJoint[b]) {
            case JointType::Free: {
                // quaternion qdot = N(q)*w ; translation qdot = v.
                const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
                const Vec3 w_FM(u[uOff], u[uOff + 1], u[uOff + 2]);
                const Vec4 qd = Rotation::convertAngVelToQuaternionDot(quat, w_FM);
                qdotOut[qOff] = qd[0];
                qdotOut[qOff + 1] = qd[1];
                qdotOut[qOff + 2] = qd[2];
                qdotOut[qOff + 3] = qd[3];
                qdotOut[qOff + 4] = u[uOff + 3];
                qdotOut[qOff + 5] = u[uOff + 4];
                qdotOut[qOff + 6] = u[uOff + 5];
                break;
            }
            default: // Weld/Pin/Translation: qdot == u
                for (int j = 0; j < dof; ++j) {
                    qdotOut[qOff + j] = u[uOff + j];
                }
                break;
        }
    }
}

// ============================================================================
//  ARTICULATED-BODY INERTIAS (inward)   RigidBodyNodeSpec.cpp
// ============================================================================
void RobotEngine::realizeArticulatedBodyInertias(const RobotModel& m, RobotState& s) {
    ArticulatedInertia* P = s.P();
    ArticulatedInertia* PPlus = s.PPlus();
    const SpatialVec* H = s.H();
    SpatialVec* G = s.G();
    Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    const SpatialInertia* Mk = s.Mk_G();

    for (int b = m.numBodies - 1; b >= 1; --b) {
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];

        // P = Mk_G + sum_children Phi * PPlus_child * ~Phi.
        ArticulatedInertia Pb(Mk[b]);
        for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
            const int c = m.bodyChildren[ci];
            Pb += PPlus[c].shift(Phi[c].l());
        }
        P[b] = Pb;

        if (dof == 0) { // Weld: nothing felt through the (zero) mobility space
            PPlus[b] = Pb;
            continue;
        }

        // PH columns, D = ~H*P*H, DI = D^-1, G = P*H*DI.
        SpatialVec PH[6];
        for (int j = 0; j < dof; ++j) {
            PH[j] = Pb * H[uOff + j];
        }
        Real D[36];
        for (int i = 0; i < dof; ++i) {
            for (int j = 0; j < dof; ++j) {
                D[i * dof + j] = spatialDot(H[uOff + i], PH[j]);
            }
        }
        Real* DI = &DIpool[m.bodyUSqIndex[b]];
        invertDense(D, dof, DI);

#if ROBO_DEBUG
        {
            bool bad = false;
            for (int k = 0; k < dof * dof; ++k) {
                if (!robodbg::fin(D[k]) || !robodbg::fin(DI[k])) {
                    bad = true;
                    break;
                }
            }
            // Also flag a wildly large inverse (near-singular D), which is the
            // precursor to a blow-up even before it becomes a literal NaN.
            Real maxDI = 0;
            for (int k = 0; k < dof * dof; ++k) {
                maxDI = std::max(maxDI, std::abs(DI[k]));
            }
            if ((bad || maxDI > Real(1e8)) && !robodbg::firstNanDumped) {
                std::cout << "\n--- realizeABI: suspect hinge inverse at body " << b << " joint "
                          << static_cast<int>(m.bodyJoint[b]) << " dof " << dof << " (max|DI|=" << maxDI
                          << ", eval#" << robodbg::evalCount << " step#" << robodbg::stepCount << ")\n";
                std::cout << "  Mk: mass=" << Mk[b].getMass() << " com=" << Mk[b].getMassCenter() << "\n";
                std::cout << "  D =\n";
                for (int i = 0; i < dof; ++i) {
                    std::cout << "    ";
                    for (int j = 0; j < dof; ++j) {
                        std::cout << std::setw(13) << D[i * dof + j] << ' ';
                    }
                    std::cout << "\n";
                }
                std::cout << "  DI =\n";
                for (int i = 0; i < dof; ++i) {
                    std::cout << "    ";
                    for (int j = 0; j < dof; ++j) {
                        std::cout << std::setw(13) << DI[i * dof + j] << ' ';
                    }
                    std::cout << "\n";
                }
                std::cout << "  H cols (Ground):\n";
                for (int j = 0; j < dof; ++j) {
                    std::cout << "    H[" << j << "] = " << H[uOff + j] << "\n";
                }
                std::cout << std::flush;
            }
        }
#endif

        for (int j = 0; j < dof; ++j) {
            SpatialVec gj(Vec3(0), Vec3(0));
            for (int k = 0; k < dof; ++k) {
                gj += PH[k] * DI[k * dof + j];
            }
            G[uOff + j] = gj;
        }

        // PPlus = P - G*~PH, reconstructed in ArticulatedInertia block form
        // exactly as Simbody does (sum of outer products over columns).
        Mat33 massMoment(0);
        Mat33 mass(0);
        Mat33 inertia(0);
        for (int j = 0; j < dof; ++j) {
            const Vec3& Ga = G[uOff + j][0];
            const Vec3& Gl = G[uOff + j][1];
            const Vec3& PHa = PH[j][0];
            const Vec3& PHl = PH[j][1];
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    massMoment(r, c) += Ga[r] * PHl[c];
                    mass(r, c) += Gl[r] * PHl[c];
                    inertia(r, c) += Ga[r] * PHa[c];
                }
            }
        }
        const SymMat33 symMass(mass(0, 0),
                               (mass(1, 0) + mass(0, 1)) / 2,
                               mass(1, 1),
                               (mass(2, 0) + mass(0, 2)) / 2,
                               (mass(2, 1) + mass(1, 2)) / 2,
                               mass(2, 2));
        const SymMat33 symInertia(inertia(0, 0),
                                  (inertia(1, 0) + inertia(0, 1)) / 2,
                                  inertia(1, 1),
                                  (inertia(2, 0) + inertia(0, 2)) / 2,
                                  (inertia(2, 1) + inertia(1, 2)) / 2,
                                  inertia(2, 2));
        PPlus[b] = Pb - ArticulatedInertia(symMass, massMoment, symInertia);
    }

    // Articulated-body centrifugal force = P * a_mob + b  (seed for calcUDot
    // pass1). CRITICAL: Simbody (RigidBodyNode.cpp,
    // realizeArticulatedBodyVelocityCache) forms this from the MOBILIZER
    // coriolis acceleration (the per-joint incremental term A), NOT the TOTAL
    // coriolis acceleration (a = ~Phi*a_parent + A). Using the total here adds
    // a spurious, velocity^2-scaled centrifugal force on every non-root body;
    // it propagates inward through the pass-1 Phi*zPlus sum, corrupts udot,
    // and feeds back through Verlet as monotonic energy injection -> blow-up.
    const SpatialVec* a_mob = s.mobCoriolisA();
    const SpatialVec* gyro = s.gyro();
    SpatialVec* abcf = s.abCentrifugal();
    for (int b = 1; b < m.numBodies; ++b) {
        abcf[b] = P[b] * a_mob[b] + gyro[b];
    }
}

// ============================================================================
//  FORWARD DYNAMICS  (forces -> udot)   RigidBodyNodeSpec.cpp
// ============================================================================
void RobotEngine::calcUDot(const RobotModel& m, RobotState& s) {
    const SpatialVec* H = s.H();
    const SpatialVec* G = s.G();
    const ArticulatedInertia* P = s.P();
    const Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    const SpatialVec* abcf = s.abCentrifugal();
    const SpatialVec* a_mob = s.mobCoriolisA();
    const SpatialVec* bodyF = s.bodyForceG();
    const Real* jointF = s.mobilityForce();

    SpatialVec* Z = s.Z();
    SpatialVec* ZPlus = s.zPlus();
    Real* eps = s.eps();
    SpatialVec* A_GB = s.A_GB();
    Real* udot = s.udot();

    // Pass 1 inward: z = (P a + b) - F + sum Phi*zPlus_child ; eps = f - ~H z ;
    // zPlus = z + G eps.
    for (int b = m.numBodies - 1; b >= 1; --b) {
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        SpatialVec z = abcf[b] - bodyF[b];
        for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
            const int c = m.bodyChildren[ci];
            z += Phi[c] * ZPlus[c];
        }
        for (int j = 0; j < dof; ++j) {
            eps[uOff + j] = jointF[uOff + j] - spatialDot(H[uOff + j], z);
        }
        SpatialVec zp = z;
        for (int j = 0; j < dof; ++j) {
            zp += G[uOff + j] * eps[uOff + j];
        }
        Z[b] = z;
        ZPlus[b] = zp;
    }

    // Pass 2 outward: APlus = ~Phi A_parent; udot = DI eps - ~G APlus;
    // A_GB = APlus + H udot + a_mob.
    A_GB[0] = SpatialVec(Vec3(0), Vec3(0));
    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const SpatialVec APlus = (~Phi[b]) * A_GB[p];
        const Real* DI = &DIpool[m.bodyUSqIndex[b]];
        SpatialVec acc = APlus;
        for (int i = 0; i < dof; ++i) {
            Real ud = 0;
            for (int j = 0; j < dof; ++j) {
                ud += DI[i * dof + j] * eps[uOff + j];
            }
            ud -= spatialDot(G[uOff + i], APlus); // -~G*APlus, row i
            udot[uOff + i] = ud;
            acc += H[uOff + i] * ud;
        }
        A_GB[b] = acc + a_mob[b];
    }
}

// ============================================================================
//  M^-1 f   (two passes, no velocity terms)   RigidBodyNodeSpec.cpp
// ============================================================================
void RobotEngine::multiplyByMInv(const RobotModel& m, RobotState& s, const Real* in, Real* out) {
    const SpatialVec* H = s.H();
    const SpatialVec* G = s.G();
    const Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    SpatialVec* Z = s.Z();
    SpatialVec* ZPlus = s.zPlus();
    Real* eps = s.eps();
    SpatialVec* A_GB = s.A_GB();

    for (int b = m.numBodies - 1; b >= 1; --b) {
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        SpatialVec z(Vec3(0), Vec3(0));
        for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
            const int c = m.bodyChildren[ci];
            z += Phi[c] * ZPlus[c];
        }
        SpatialVec zp = z;
        for (int j = 0; j < dof; ++j) {
            eps[uOff + j] = in[uOff + j] - spatialDot(H[uOff + j], z);
        }
        for (int j = 0; j < dof; ++j) {
            zp += G[uOff + j] * eps[uOff + j];
        }
        Z[b] = z;
        ZPlus[b] = zp;
    }

    A_GB[0] = SpatialVec(Vec3(0), Vec3(0));
    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const SpatialVec APlus = (~Phi[b]) * A_GB[p];
        const Real* DI = &DIpool[m.bodyUSqIndex[b]];
        SpatialVec acc = APlus;
        for (int i = 0; i < dof; ++i) {
            Real v = 0;
            for (int j = 0; j < dof; ++j) {
                v += DI[i * dof + j] * eps[uOff + j];
            }
            v -= spatialDot(G[uOff + i], APlus);
            out[uOff + i] = v;
            acc += H[uOff + i] * v;
        }
        A_GB[b] = acc;
    }
}

// ============================================================================
//  sqrt(M^-1) f   (HMC velocity seeding)   RigidBodyNodeSpec.cpp
// ============================================================================
void RobotEngine::multiplyBySqrtMInv(const RobotModel& m, RobotState& s, const Real* in, Real* out) {
    const SpatialVec* H = s.H();
    const SpatialVec* G = s.G();
    const Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    SpatialVec* V = s.V_GB(); // reuse as scratch outward velocity

    V[0] = SpatialVec(Vec3(0), Vec3(0));
    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const SpatialVec VPlus = (~Phi[b]) * V[p];
        Real sqrtDI[36];
        symSqrt(&DIpool[m.bodyUSqIndex[b]], dof, sqrtDI);
        SpatialVec acc = VPlus;
        for (int i = 0; i < dof; ++i) {
            Real v = 0;
            for (int j = 0; j < dof; ++j) {
                v += sqrtDI[i * dof + j] * in[uOff + j];
            }
            v -= spatialDot(G[uOff + i], VPlus);
            out[uOff + i] = v;
            acc += H[uOff + i] * v;
        }
        V[b] = acc;
    }
}

robo::Real RobotEngine::calcLogDetM(const RobotModel& m, const RobotState& s) {
    using robo::Real;
    using robo::SpatialVec;
    const SpatialVec* H = s.H();
    const robo::ArticulatedInertia* P = s.P();

    Real logDet = Real(0);
    Real D[36]; // per-body dof x dof, dof <= 6; the global M is never formed
    for (int b = 1; b < m.numBodies; ++b) {
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        if (dof == 0) {
            continue; // Weld: contributes det = 1 -> log 0
        }
        // D_b = ~H_b * P_b * H_b  (exactly as in realizeArticulatedBodyInertias)
        for (int j = 0; j < dof; ++j) {
            const SpatialVec PHj = P[b] * H[uOff + j];
            for (int i = 0; i < dof; ++i) {
                D[i * dof + j] = spatialDot(H[uOff + i], PHj);
            }
        }
        logDet += logDetSymPD(D, dof); // ln det(D_b)
    }
    return logDet; // = ln|M_phi|
}


// ============================================================================
//  KINETIC ENERGY = 1/2 u^T M u  via the spatial sweep (Mk_G, V_GB).
// ============================================================================
Real RobotEngine::calcKineticEnergy(const RobotModel& m, const RobotState& s) {
    const SpatialInertia* Mk = s.Mk_G();
    const SpatialVec* V = s.V_GB();
    Real ke = 0;
    for (int b = 1; b < m.numBodies; ++b) {
        const SpatialVec MV = Mk[b] * V[b];
        ke += spatialDot(V[b], MV);
    }
    return Real(0.5) * ke;
}

// ============================================================================
//  TRANSFER: internal bodies -> per-atom Cartesian (accept path)
// ============================================================================
void RobotEngine::fillAtomPositionsFromBodies(const RobotModel& m, RobotState& s) {
    const Transform* X_GB = s.X_GB();
    Vec3* posG = s.atomPosG();
    for (int a = 0; a < m.numAtoms; ++a) {
        const int b = m.atomBody[a];
        posG[a] = X_GB[b].p() + X_GB[b].R() * m.atomStation_B[a];
    }
}

// ============================================================================
//  Quaternion renormalization (the no-constraint Verlet "projection").
// ============================================================================
void RobotEngine::normalizeQuaternions(const RobotModel& m, RobotState& s) {
    Real* q = s.q();
    for (int qs : m.quaternionQStart) {
        const Real n2 = q[qs] * q[qs] + q[qs + 1] * q[qs + 1] + q[qs + 2] * q[qs + 2] + q[qs + 3] * q[qs + 3];
        if (n2 <= Real(1e-24)) {
            // Degenerate (near-zero) quaternion: reset to identity instead of
            // dividing by ~0. Should not happen now that q is seeded to identity,
            // but this guarantees no 1/sqrt(0) -> NaN can ever recur.
            q[qs] = Real(1);
            q[qs + 1] = q[qs + 2] = q[qs + 3] = Real(0);
            continue;
        }
        const Real inv = Real(1) / std::sqrt(n2);
        q[qs] *= inv;
        q[qs + 1] *= inv;
        q[qs + 2] *= inv;
        q[qs + 3] *= inv;
    }
}

// ============================================================================
//  VERLET  (fixed step)   VerletIntegrator.cpp::attemptDAEStep
// ============================================================================
bool RobotEngine::verletStep(const RobotModel& m,
                             RobotState& s,
                             ForceBridge& bridge,
                             const robo::ConstraintSet& cset,
                             Real h) {
    const int nq = m.nq, nu = m.nu;
    Real* q = s.q();
    Real* u = s.u();
    Real* qdot = s.qdot();
    Real* udot = s.udot();
    Real* qdd = s.qdotdot();

    std::vector<Real> q0(q, q + nq), u0(u, u + nu), qdot0(qdot, qdot + nq), udot0(udot, udot + nu),
        qdd0(qdd, qdd + nq);

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
    auto restorePreStep = [&]() {
        std::copy(q0.begin(), q0.end(), q);
        std::copy(u0.begin(), u0.end(), u);
        std::copy(qdot0.begin(), qdot0.end(), qdot);
        std::copy(udot0.begin(), udot0.end(), udot);
        std::copy(qdd0.begin(), qdd0.end(), qdd);
        realizePosition(m, s);
        fillAtomPositionsFromBodies(m, s);
    };

    // ---- position: q1 = q0 + h*qdot0 + (h^2/2)*qddot0, normalize, then SHAKE ----
    for (int i = 0; i < nq; ++i) {
        q[i] = q0[i] + h * qdot0[i] + (h * h / 2) * qdd0[i];
    }
    normalizeQuaternions(m, s);

    auto refreshPos = [&]() {
        realizePosition(m, s);
        fillAtomPositionsFromBodies(m, s);
    };
    refreshPos();
    cset.enforcePositionConstraints(m, s, refreshPos); // localProjectQ

    // ---- velocity: implicit trapezoid + functional iteration ----
    for (int i = 0; i < nu; ++i) {
        u[i] = u0[i] + h * udot0[i]; // u1_est
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
        calcQDot(m, s, qdot);
        calcQDotDot(m, s);
        return true;
    };
    if (!evalDerivs()) {
        restorePreStep();
        return false;
    }

    const Real tol = Real(1e-4);
    Real prevChange = std::numeric_limits<Real>::infinity();
    for (int iter = 0; iter < 10; ++iter) {
        Real num = 0, den = 0; // Simbody's relative 2-norm change
        for (int i = 0; i < nu; ++i) {
            const Real un = u0[i] + (h / 2) * (udot0[i] + udot[i]);
            const Real d = un - u[i];
            num += d * d;
            den += u[i] * u[i];
            u[i] = un;
        }
        if (!evalDerivs()) {
            restorePreStep();
            return false;
        }
        const Real change = std::sqrt(num) / (std::sqrt(den) + Real(1e-30));
        if (change <= tol) {
            break; // converged
        }
        if (iter > 1 && change > prevChange) {
            break; // FIX: Simbody's break, only after i>1
        }
        prevChange = change;
    }

    cset.enforceVelocityConstraints(m, s); // localProjectU (RATTLE)
    realizeVelocity(m, s);                 // refresh V/KE at the projected u
    s.time += h;
    return true;
}

bool RobotEngine::stepTo(const RobotModel& m,
                         RobotState& s,
                         ForceBridge& bridge,
                         const robo::ConstraintSet& cset,
                         Real tEnd) {
    const Real h = tEnd - s.time;
    if (h <= 0) {
        return true;
    }
    return verletStep(m, s, bridge, cset, h);
}

// ---- calcQDotDot: qddot = N qddot-coupling. Pin/Translation: udot; Free: quat. ----
void RobotEngine::calcQDotDot(const RobotModel& m, RobotState& s) {
    const Real* u = s.u();
    const Real* udot = s.udot();
    const Real* q = s.q();
    Real* qdd = s.qdotdot();
    for (int b = 1; b < m.numBodies; ++b) {
        const int qOff = m.bodyQIndex[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        if (m.bodyJoint[b] == JointType::Free) {
            const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            const Vec3 w(u[uOff], u[uOff + 1], u[uOff + 2]);
            const Vec3 wd(udot[uOff], udot[uOff + 1], udot[uOff + 2]);
            const Vec4 qdd4 = Rotation::convertAngVelDotToQuaternionDotDot(quat, w, wd);
            for (int k = 0; k < 4; ++k) {
                qdd[qOff + k] = qdd4[k];
            }
            qdd[qOff + 4] = udot[uOff + 3];
            qdd[qOff + 5] = udot[uOff + 4];
            qdd[qOff + 6] = udot[uOff + 5];
        } else {
            for (int j = 0; j < dof; ++j) {
                qdd[qOff + j] = udot[uOff + j];
            }
        }
    }
}

// Per-joint kernels (X_FM, H_FM, qdot, qddot) are inlined directly into the
// realize* routines above; there are no separate helper symbols.