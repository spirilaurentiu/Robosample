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
//    per-joint H_FM / X_FM / qdot ............. RigidBodyNodeSpec_{Torsion,Translation,Free}.h
//
//  VALIDATION: Simbody has been REMOVED from the build (no find_package, no
//  SimTK/ Simbody01 tree, no live MobilizedBody code anywhere in src/ or tests/).
//  The historical "diff every operator against the still-present Simbody build"
//  gate therefore NO LONGER EXISTS and was never wired into the suite. The
//  oracle-of-record for these recursions is now:
//    * analytic closed forms (single/free rigid body, physical pendulum, KE),
//    * finite-difference cross-checks (H<->X_FM, V<->X_GB, A<->V, HDot<->H,
//      logDetM<->Cholesky of the forward-built dense M),
//    * operator round-trips (M*M^-1==I, sqrtMInv/sqrtM mutual inverse,
//      u^T M u == w^T w for u = multiplyBySqrtMInv(w)),
//    * statistical ensembles (equipartition, Fixman/Boltzmann marginal),
//    * frozen Simbody-golden CONSTANTS transcribed into TestSpatialAlgebra /
//      TestInertia (numbers, not a live diff), and
//    * a live OpenMM/ParmEd differential for POTENTIAL energy only.
//  Field/routine names below still track their Simbody originals (source map
//  above) purely as provenance, not as a runtime comparison.
// ============================================================================

#include "RobotEngine.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Constraints.hpp"
#include "JointKernels.hpp"
#include "robot_math.hpp"

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

// n x n symmetric eigensolver (cyclic Jacobi), defined below; forward-declared
// here so the null-space-lock helpers (this block) can use it before their
// point of definition in the file.
static void jacobiSymEig(const robo::Real* Ain, int n, robo::Real* d, robo::Real* V);

// ============================================================================
//  SHARED NULL-SPACE LOCK (docs/specs/singular-dof-fixman.md CC1/CC4/N1).
//
//  A direction of a dof x dof hinge-inertia block D_b = ~H_b P_b H_b is
//  "null" (a structural/gauge phantom, e.g. a leaf single-atom Torsion body
//  whose atom sits ON its own rotation axis) iff its magnitude is at or
//  below this RELATIVE lock, floored at an absolute constant so an
//  all-tiny block is still locked rather than scaled away to nothing. This
//  is the SINGLE source of truth: invertDense (dynamics pseudo-inverse),
//  pseudoLogDet (Fixman pseudo-determinant), and symSqrtInv (sqrt(M) for
//  NMA route B) all lock on exactly this test, so a direction removed by
//  the dynamics contributes factor 1 (0 to ln det) to the determinant and 0
//  to sqrt(D) -- never a floating-point residue treated as physically
//  present (the pre-fix inconsistency: invertDense zeroed the direction,
//  logDetSymPD's separate 1e-300 Cholesky-pivot clamp did not, so it added
//  ln(residue) instead of ln(1)==0).
// ============================================================================
constexpr robo::Real kNullLockAbs = robo::Real(1e-12);

inline robo::Real nullLockTol(robo::Real scale) {
    return std::max(kNullLockAbs, scale * kNullLockAbs);
}

// Eigen-decompose an n x n symmetric block and compute its null-space lock
// tolerance in one place, so every caller locks on IDENTICAL (d, V, tol) for
// the same input matrix. n == 1 is closed form (no eigensolve needed): the
// single "eigenvalue" is the scalar itself and the tolerance is the bare
// absolute constant (matches invertDense's pre-existing 1-dof fast path,
// unchanged by this refactor).
inline void eigDecompAndTol(const robo::Real* A, int n, robo::Real* d, robo::Real* V, robo::Real& tol) {
    if (n == 1) {
        d[0] = A[0];
        V[0] = robo::Real(1);
        tol = kNullLockAbs;
        return;
    }
    jacobiSymEig(A, n, d, V);
    robo::Real scale = 0;
    for (int k = 0; k < n; ++k) {
        scale = std::max(scale, std::abs(d[k]));
    }
    tol = nullLockTol(scale);
}

// ln|det D_b|, the PSEUDO-determinant: a null direction (per the shared lock
// above) contributes 0 to the sum (factor 1 to det), never ln(residue). This
// is what makes calcLogDetM agree with invertDense (CC1) -- a direction the
// dynamics pseudo-inverse removes must not still inflate/deflate ln|M_tree|
// by a floating-point accident. D_b is SPD for a well-posed hinge, so every
// non-null eigenvalue is expected positive; std::abs is defensive against
// roundoff-sign noise exactly at the lock boundary, not a silent sign flip
// of a real negative eigenvalue (that would be non-PD and is a separate bug,
// surfaced by NaN/Inf downstream, not masked here). Returns 0 for n == 0
// (Weld: det == 1).
inline robo::Real pseudoLogDet(const robo::Real* A, int n) {
    if (n <= 0) {
        return robo::Real(0);
    }
    robo::Real d[6];
    robo::Real V[36];
    robo::Real tol;
    eigDecompAndTol(A, n, d, V, tol);
    robo::Real logDet = robo::Real(0);
    for (int k = 0; k < n; ++k) {
        if (std::abs(d[k]) <= tol) {
            continue; // null direction: contributes 0 (pseudo-determinant), not ln(residue)
        }
        logDet += std::log(std::abs(d[k]));
    }
    return logDet;
}


// ============================================================================
//  DEBUG LOGGING
//  ROBO_DEBUG   : master switch (0 disables every probe with zero overhead).
//  ROBO_VERBOSE : 0 = only the first-NaN dump + the D/DI dump at the failing
//                     body (recommended: Torsionpoints the origin, low spam);
//                 1 = + a one-line per-step growth summary;
//                 2 = + per-body X_GB/X_FM every step (huge; short runs only).
//  Override at compile time, e.g.  -DROBO_VERBOSE=1 .
// ============================================================================
#define ROBO_DEBUG 1
#define ROBO_VERBOSE 2

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

// Symmetric (Moore-Penrose) pseudo-inverse of the dof x dof articulated hinge
// matrix D = ~H P H, row-major in/out. D is symmetric and, for a well-posed
// hinge, SPD; this routine returns its inverse. Returns the number of
// LOCKED (null) directions, 0 for a well-conditioned D -- the realize-ABI
// caller uses this to fire the CC1/CC4 fail-loud gate (docs/specs/
// singular-dof-fixman.md STEP 3) on any body whose lock was not expected.
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
int invertDense(const Real* A, int n, Real* Ainv) {
    // Dominant case: 1-DOF Torsion. Closed form, with a lock instead of a clamp.
    if (n == 1) {
        const Real d = A[0];
        const bool locked = std::abs(d) <= kNullLockAbs;
        Ainv[0] = locked ? Real(0) : (Real(1) / d);
        return locked ? 1 : 0;
    }

    // n in {3,6}: symmetric pseudo-inverse via the in-house Jacobi eigensolver
    // (self-contained; no LAPACK/Eigen link dependency). For a non-degenerate
    // hinge every eigenvalue is well above tol, so this equals the true inverse;
    // any null eigendirection contributes 0 (that direction is locked).
    Real d[6];
    Real V[36];
    Real tol;
    eigDecompAndTol(A, n, d, V, tol);
    for (int i = 0; i < n * n; ++i) {
        Ainv[i] = Real(0);
    }
    int numLocked = 0;
    for (int k = 0; k < n; ++k) {
        if (std::abs(d[k]) <= tol) {
            ++numLocked; // null direction -> pseudo-inverse contributes 0 (locked)
            continue;
        }
        const Real inv = Real(1) / d[k];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Ainv[i * n + j] += inv * V[i * n + k] * V[j * n + k];
            }
        }
    }
    return numLocked;
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

// Inverse symmetric square root of a dof x dof block A == DI (invertDense's
// pseudo-inverse of D): S = V diag(1/sqrt(lambda)) ~V. On every UNLOCKED
// direction this is the exact inverse of symSqrt(DI) (same eigenvectors,
// reciprocal sqrt eigenvalues), so symSqrtInv(DI)*symSqrt(DI) == I there --
// verified numerically to ~2e-13 over random SPD blocks (dof 1..6). Used by
// multiplyBySqrtM (the forward sqrt sweep, NMA route B).
//
// NULL-DIRECTION CONVENTION (S3, docs/specs/singular-dof-fixman.md): DI's
// eigenvalue is EXACTLY 0 on any direction invertDense locked (a structural
// phantom). The OLD code floored that 0 to 1e-300 and returned
// 1/sqrt(1e-300) ~ 3e149 -- an enormous but technically-finite regularizer
// that is NOT the pseudo-inverse convention every other null-space lock in
// this file uses, and is only harmless today because the caller's `in` on
// that direction is always bit-exact 0 (u_phantom is frozen by the SAME
// lock), so huge*0==0. That correctness depends on an invariant this
// function cannot see or enforce. Locking the direction to 0 HERE instead
// (same shared tol as invertDense/pseudoLogDet) makes symSqrtInv robust on
// its own terms -- multiplyBySqrtM's derivation only needs sqrt(D)==0 on a
// direction that mobility WILL NEVER be seeded on, and 0 is the correct
// "no information on this gauge direction" answer, not the largest float
// that still avoids inf.
void symSqrtInv(const Real* A, int n, Real* S) {
    if (n == 1) {
        const bool locked =
            A[0] <= kNullLockAbs; // DI >= 0 always; <=0 means locked (or non-PSD, caught downstream)
        S[0] = locked ? Real(0) : (Real(1) / std::sqrt(A[0]));
        return;
    }
    Real d[6];
    Real V[36];
    Real tol;
    eigDecompAndTol(A, n, d, V, tol);
    for (int i = 0; i < n * n; ++i) {
        S[i] = Real(0);
    }
    for (int k = 0; k < n; ++k) {
        if (std::abs(d[k]) <= tol) {
            continue; // null direction -> contributes 0 (locked; matches invertDense/pseudoLogDet)
        }
        const Real inv = Real(1) / std::sqrt(d[k]);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                S[i * n + j] += inv * V[i * n + k] * V[j * n + k];
            }
        }
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

    // Kinetic-metric preconditioning: fictitious per-body inertia scale applied
    // to Mk_G ONLY (draw + KE + Fixman all read Mk_G, so they stay consistent).
    // Empty => all 1.0 (physical). Hoisted out of the hot loop.
    const Real* massScale = m.bodyMassScale.empty() ? nullptr : m.bodyMassScale.data();

    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const JointType jt = m.bodyJoint[b];
        const int qOff = m.bodyQIndex[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];

        // X_FM, then X_PB = X_PF * X_FM * X_MB, X_GB = X_GP * X_PB.
        X_FM[b] = robo::jointX_FM(jt, q, qOff);
        const Transform X_MB = ~m.X_BM[b];
        const Transform X_FB = X_FM[b] * X_MB;
        X_PB[b] = m.X_PF[b] * X_FB;
        X_GB[b] = X_GB[p] * X_PB[b];

        // H_FM (in F) and H_PB_G (in Ground): RigidBodyNodeSpec.cpp.
        jointH_FM(jt, X_FM[b], &HFM[uOff]);
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
        // bodyMass is the only scaled argument: UnitInertia is mass-normalised, so
        // mass*scale propagates the factor through the whole spatial inertia (mass,
        // first moment, and inertia tensor) -> a uniformly heavier, SPD-valid copy
        // of the same rigid body. Bias-free (see RobotModel::bodyMassScale).
        const Real massB = massScale ? (m.bodyMass[b] * massScale[b]) : m.bodyMass[b];
        Mk[b] = SpatialInertia(massB, p_BBc_G, G_Bo_G);
    }

    // Per-atom Ground positions and stations (R_GB * station_B). The transfer
    // currency + the input to the OpenMM force bridge. Each atom is independent (no
    // reduction), so this vectorizes cleanly: __restrict removes the aliasing barrier
    // between the destination arrays and the (const) source arrays, and `omp simd`
    // (SIMD-only, -fopenmp-simd) lets the compiler pack the R*station FMAs across atoms.
    // Bitwise-identical to the scalar form.
    Vec3* __restrict posG = s.atomPosG();
    Vec3* __restrict stG = s.atomStationG();
    const int* __restrict atomBody = m.atomBody.data();
    const Vec3* __restrict station = m.atomStation_B.data();
#pragma omp simd
    for (int a = 0; a < m.numAtoms; ++a) {
        const int b = atomBody[a];
        const Vec3 st = X_GB[b].R() * station[a];
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
        const JointType jt = m.bodyJoint[b];

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
        // Torsion/torsion body, where r_MB is the bond length). With the per-joint
        // H_FM constant in F (Torsion/Free/Translation/Weld all have HDot_FM = 0),
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

        // For the seven joints with H_FM constant in F, dH_FM/dt = 0 and the
        // mobilizer bias is exactly the line above (this branch is skipped),
        // so those joints stay BIT-IDENTICAL to the prior engine. The three
        // q-dependent joints (BendStretch/SphericalCoords/FreeLine) add the
        // intrinsic dH_FM/dt*u term that the constant-H derivation dropped:
        //   extraAng_F = sum_j Hdot_ang(j) u_j
        //   extraLin_F = sum_j (Hdot_lin(j) - r_MB_F x Hdot_ang(j)) u_j
        // (The -rdot_MB_F x H_ang(j) part of the full d/dt(H_MB_F) is already
        //  captured by centripetal_G, so it is NOT re-added here.)
        Vec3 extraAng_F(0, 0, 0);
        Vec3 extraLin_F(0, 0, 0);

        if (!RobotModel::jointHasConstantHFM(jt)) {
            std::array<SpatialVec, 5> Hdot; // max dof = 5 (FreeLine)
            jointHDot_FM(jt, X_FM[b], V_FM[b], Hdot.data());
            for (int j = 0; j < dof; ++j) {
                const Real uj = u[uOff + j];
                extraAng_F += Hdot[j][0] * uj;
                extraLin_F += (Hdot[j][1] - (r_MB_F % Hdot[j][0])) * uj;
            }
        }

        const SpatialVec A_mob(w_GP % w_PB_G + R_GF * extraAng_F,
                               w_GP % (v_GB - v_GP) + w_GP % v_PB_G + centripetal_G + R_GF * extraLin_F);
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
            case JointType::Ball: {
                // 3 dof rotation: quaternion qdot = N(q) * w_FM. No translation.
                const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
                const Vec3 w_FM(u[uOff], u[uOff + 1], u[uOff + 2]);
                const Vec4 qd = Rotation::convertAngVelToQuaternionDot(quat, w_FM);
                qdotOut[qOff] = qd[0];
                qdotOut[qOff + 1] = qd[1];
                qdotOut[qOff + 2] = qd[2];
                qdotOut[qOff + 3] = qd[3];
                break;
            }
            case JointType::FreeLine: {
                // 2 rotational speeds are (x,y) of w_FM expressed in M, so
                // w_FM (in F) = R_FM * (u0, u1, 0); qdot_quat = N(q) * w_FM.
                // The 3 translational speeds (u2..4) are qdot of x,y,z directly.
                const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
                const Rotation& R_FM = X_FM[b].R();
                const Vec3 w_FM = R_FM * Vec3(u[uOff], u[uOff + 1], 0);
                const Vec4 qd = Rotation::convertAngVelToQuaternionDot(quat, w_FM);
                qdotOut[qOff] = qd[0];
                qdotOut[qOff + 1] = qd[1];
                qdotOut[qOff + 2] = qd[2];
                qdotOut[qOff + 3] = qd[3];
                qdotOut[qOff + 4] = u[uOff + 2];
                qdotOut[qOff + 5] = u[uOff + 3];
                qdotOut[qOff + 6] = u[uOff + 4];
                break;
            }
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
            default: // Weld/Torsion/Slider/Cylinder/Translation/BendStretch/SphericalCoords: qdot == u
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
// Position-only articulated-inertia factorization: P, PPlus, D, DI, G (the per-body
// Jacobi eigensolve / invertDense). Pure function of q (masses + geometry) -- NO velocity
// dependence -- so the verlet corrector hoists it to once/step (see the wrapper below and
// docs/specs/gpu-cartesian-kinematics/03-aba-parallelization.md Sec.0.5). The velocity-
// dependent centrifugal seed is split into seedArticulatedCentrifugal.
void RobotEngine::factorizeArticulatedInertias(const RobotModel& m, RobotState& s) {
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
        const int numLocked = invertDense(D, dof, DI);

        // STEP 3 fail-loud gate (docs/specs/singular-dof-fixman.md, Review
        // outcome): invertDense null-locked >=1 direction of this body's
        // hinge inertia. That is EXPECTED and safe only for a body whose
        // D_b is PROVABLY a run-constant (locked at this q => locked for
        // EVERY q, so freezing it is bias-free, C3) -- never merely
        // "small at this reference config" (B2's rejected false-positive
        // risk: an angle-flexible multi-dof joint, e.g. BendStretch/
        // SphericalCoords/Cartesian/FreeLine, can be genuinely collinear
        // only at isolated q and must NOT be treated as structural).
        //
        // A LEAF (no children) Torsion (1-dof) body qualifies unconditionally:
        // D_b = ~H_b P_b H_b with H_b (Ground) = R_GF * H_FM, H_FM CONSTANT
        // for Torsion (jointHasConstantHFM) and R_GF depending only on
        // ancestors (upstream of this body's own q); P_b == Mk_G[b] exactly
        // (leaf: no PPlus_child term). Writing Mk_G[b] = R_GB Mk_B R_GB^T
        // with R_GB = R_GF R_FM(q_own) R_MB, the R_GF factors cancel by
        // orthogonality in D_b = H_FM^T [R_FM(q_own) R_MB Mk_B R_MB^T
        // R_FM(q_own)^T] H_FM, which no longer contains R_GF at all -- D_b
        // is therefore independent of every ancestor config AND of the
        // body's own q (a revolute joint's own-axis inertia is invariant to
        // rotation about that same axis, the standard rigid-body fact).
        // Confirmed empirically to ~1e-14 by the Cyclic1APQPhantomLogDetIsRunConstant
        // LEMMA test (tests/TestRoboticsOracleMolecule.cpp) on 1APQ's three
        // leaf single-atom on-axis Torsion phantoms. STEP 4 (a build-time
        // weld removing these bodies outright, making this gate unreachable
        // for them) is DEFERRED -- see the coder checkpoint; this runtime
        // recognition is the provably-safe substitute while it is deferred.
        // Any OTHER locked body (non-leaf, non-Torsion, or multi-dof) is NOT
        // this shape and throws instead of silently locking (Rule 11).
        if (numLocked > 0
            && !(m.bodyChildrenBeg[b] == m.bodyChildrenEnd[b] && m.bodyJoint[b] == JointType::Torsion)) {
            Real eig[6];
            Real eigV[36];
            Real eigTol;
            eigDecompAndTol(D, dof, eig, eigV, eigTol);
            Real minEig = std::abs(eig[0]);
            for (int k = 1; k < dof; ++k) {
                minEig = std::min(minEig, std::abs(eig[k]));
            }
            std::string msg = "realizeArticulatedBodyInertias: invertDense null-locked "
                              + std::to_string(numLocked) + " direction(s) of D_b on body "
                              + std::to_string(b)
                              + " (JointType=" + std::to_string(static_cast<int>(m.bodyJoint[b]))
                              + ", dof=" + std::to_string(dof) + ", atoms=[";
            // bodyAtomsBeg/End are only populated for a fully-built (real-molecule)
            // RobotModel; a hand-built synthetic test model (tests/RobotBuilders.hpp
            // buildForest without attachAtoms) leaves them empty -- guard the index
            // so a genuine gate failure never masks itself behind an out-of-bounds
            // read while formatting the diagnostic.
            if (static_cast<std::size_t>(b) < m.bodyAtomsBeg.size()) {
                for (int a = m.bodyAtomsBeg[b]; a < m.bodyAtomsEnd[b]; ++a) {
                    if (a != m.bodyAtomsBeg[b]) {
                        msg += ",";
                    }
                    msg += std::to_string(m.bodyAtoms[a]);
                }
            } else {
                msg += "unavailable: synthetic model with no atom map";
            }
            msg += "], min-eig(D)=" + std::to_string(minEig)
                   + ") that is not a recognized structural phantom (leaf Torsion) -- "
                     "see docs/specs/singular-dof-fixman.md STEP 3";
            throw std::runtime_error(msg);
        }

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
}

// Velocity-dependent articulated centrifugal seed: abcf = P * a_mob + gyro (seed for
// calcUDot pass1). This is the ONLY part of the old realizeArticulatedBodyInertias that
// depends on u; splitting it out lets the corrector iterate just this + calcUDot while the
// (expensive, position-only) factorization above runs once/step. PRECONDITION: realizeVelocity
// current (a_mob/gyro) and factorizeArticulatedInertias current (P).
void RobotEngine::seedArticulatedCentrifugal(const RobotModel& m, RobotState& s) {
    // CRITICAL: Simbody (RigidBodyNode.cpp, realizeArticulatedBodyVelocityCache) forms this
    // from the MOBILIZER coriolis acceleration (the per-joint incremental term A), NOT the
    // TOTAL coriolis acceleration (a = ~Phi*a_parent + A). Using the total here adds a
    // spurious, velocity^2-scaled centrifugal force on every non-root body; it propagates
    // inward through the pass-1 Phi*zPlus sum, corrupts udot, and feeds back through Verlet
    // as monotonic energy injection -> blow-up.
    const ArticulatedInertia* P = s.P();
    const SpatialVec* a_mob = s.mobCoriolisA();
    const SpatialVec* gyro = s.gyro();
    SpatialVec* abcf = s.abCentrifugal();
    for (int b = 1; b < m.numBodies; ++b) {
        abcf[b] = P[b] * a_mob[b] + gyro[b];
    }
}

// Backward-compatible full pass = factorization + centrifugal seed, in the original order.
// Callers other than the verlet corrector (reinitialize, captureReactionSnapshot, ...) use
// this and are unchanged; the corrector (RobotIntegrator.hpp) calls the two halves
// separately, hoisting the factorization out of the u-only iteration.
void RobotEngine::realizeArticulatedBodyInertias(const RobotModel& m, RobotState& s) {
    factorizeArticulatedInertias(m, s);
    seedArticulatedCentrifugal(m, s);
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

// ============================================================================
//  sqrt(M) f   (exact inverse of multiplyBySqrtMInv; for NMA Route B acceptance)
// ============================================================================
// multiplyBySqrtMInv maps in -> out = M^(-1/2) in via the outward recurrence
//     VPlus_b = ~Phi_b V_p ;  out_b = sqrtDI_b in_b - ~G_b VPlus_b ;
//     V_b     = VPlus_b + H_b out_b .
// Inverting body-by-body (out_b is this sweep's INPUT, in_b its OUTPUT):
//     in_b = sqrtDI_b^{-1} ( out_b + ~G_b VPlus_b ) ,
// while V_b is rebuilt from the SAME (given) out_b, so V_b is identical to the
// forward sweep's. Substituting the forward out_b shows in_b is recovered
// exactly: sqrtDI^{-1}(sqrtDI in_b - ~G VPlus + ~G VPlus) = in_b. Hence this is
// the exact algebraic inverse for any H,G,DI,Phi -- no approximation. It uses a
// LOCAL scratch velocity (NOT s.V_GB()), so it is safe to call after a trajectory
// when V_GB holds the real end velocities needed by calcKineticEnergy.
void RobotEngine::multiplyBySqrtM(const RobotModel& m, RobotState& s, const Real* in, Real* out) {
    const SpatialVec* H = s.H();
    const SpatialVec* G = s.G();
    const Real* DIpool = s.DI();
    const PhiMatrix* Phi = s.Phi();
    std::vector<SpatialVec> V(m.numBodies, SpatialVec(Vec3(0), Vec3(0)));

    for (int b = 1; b < m.numBodies; ++b) {
        const int p = m.bodyParent[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const SpatialVec VPlus = (~Phi[b]) * V[p];
        Real sqrtDIinv[36];
        symSqrtInv(&DIpool[m.bodyUSqIndex[b]], dof, sqrtDIinv);
        // tmp_i = in_b[i] + ~G_b[i] . VPlus     (in[] plays the forward out_b)
        Real tmp[6];
        for (int i = 0; i < dof; ++i) {
            tmp[i] = in[uOff + i] + spatialDot(G[uOff + i], VPlus);
        }
        // out_b = sqrtDI^{-1} tmp
        for (int i = 0; i < dof; ++i) {
            Real v = 0;
            for (int j = 0; j < dof; ++j) {
                v += sqrtDIinv[i * dof + j] * tmp[j];
            }
            out[uOff + i] = v;
        }
        // V_b = VPlus + H_b in_b   (rebuilt from the given out_b == in[], matching forward)
        SpatialVec acc = VPlus;
        for (int i = 0; i < dof; ++i) {
            acc += H[uOff + i] * in[uOff + i];
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
        logDet += pseudoLogDet(D, dof); // pseudo-ln-det(D_b): null directions contribute 0 (CC1)
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
//  MOBILIZER REACTION FORCES   (port of calcMobilizerReactionForces)
//
//  The reaction transmitted to body b across its inboard mobilizer, in Ground,
//  is rigid Newton-Euler with the TRUE accelerations A_GB plus what the children
//  pass back inward:
//      reac_b@Bo = Mk_b A_GB_b + gyro_b - F_ext_b + sum_c Phi[c] reac_c@Bo
//  F_ext_b is the applied spatial BODY force on b only (bodyForceG from the
//  bridge). Matching Simbody (SimbodyMatterSubsystemRep::calcMobilizerReactionForces
//  / calcMobilizerReactionForcesUsingFreebodyMethod, SimbodyMatterSubsystemRep.cpp
//  :6061-6062,6067-6168): "any generalized forces applied at the mobilities end
//  up included in the reaction forces" -- i.e. applied mobility (generalized
//  joint) forces are NOT subtracted out here; the reported reaction is the one
//  actually transmitted across the joint given whatever generalized force was
//  applied. NOTE: this is behavior-neutral today because the force bridge
//  zeroes mobilityForce every step (include/ForceBridge.hpp:73-75), so mobF is
//  identically 0 and dropping the H*mobF term changes nothing numerically. It
//  becomes live once a Fixman/biasing generalized torque is introduced --  at
//  that point this convention (reaction includes actuation) is the intended one.
//  Phi[c] (offset = parent->child origin in Ground) shifts a child's force from
//  the child origin to b's origin -- identical to calcUDot's pass-1 transmission.
// ============================================================================
void RobotEngine::calcMobilizerReactionForces(const RobotModel& m,
                                              const RobotState& s,
                                              SpatialVec* reactionAtBoInG,
                                              SpatialVec* reactionAtMInG) {
    const SpatialInertia* Mk = s.Mk_G();
    const SpatialVec* A_GB = s.A_GB();
    const SpatialVec* gyro = s.gyro();
    const SpatialVec* bodyF = s.bodyForceG();
    const PhiMatrix* Phi = s.Phi();
    const Transform* X_GB = s.X_GB();

    // reaction at body origin, accumulated inward. Local scratch so the operator
    // is side-effect free on the cache (callers may not want a dedicated slot).
    std::vector<SpatialVec> reacBo(static_cast<std::size_t>(m.numBodies), SpatialVec(Vec3(0), Vec3(0)));

    for (int b = m.numBodies - 1; b >= 1; --b) {
        // applied spatial force on b: bridge body force only (see Simbody
        // convention note above -- generalized/mobility forces are NOT
        // subtracted here).
        const SpatialVec& fExt = bodyF[b];

        // rigid Newton-Euler residual at Bo, in Ground.
        SpatialVec reac = (Mk[b] * A_GB[b]) + gyro[b] - fExt;

        // add what the outboard children transmit inward (shift child Bo -> b Bo).
        for (int ci = m.bodyChildrenBeg[b]; ci < m.bodyChildrenEnd[b]; ++ci) {
            const int c = m.bodyChildren[ci];
            reac += Phi[c] * reacBo[static_cast<std::size_t>(c)];
        }

        reacBo[static_cast<std::size_t>(b)] = reac;
        if (reactionAtBoInG != nullptr) {
            reactionAtBoInG[b] = reac;
        }
        if (reactionAtMInG != nullptr) {
            // shift the spatial force from Bo to the outboard frame origin Mo:
            //   p_BoMo_G = R_GB * X_BM.p ;  [t;f]@Bo -> [t - p x f ; f]@Mo.
            const Vec3 p_BoMo_G = X_GB[b].R() * m.X_BM[b].p();
            reactionAtMInG[b] = SpatialVec(reac.angular - (p_BoMo_G % reac.linear), reac.linear);
        }
    }

    if (reactionAtBoInG != nullptr) {
        reactionAtBoInG[0] = SpatialVec(Vec3(0), Vec3(0)); // Ground: no inboard joint
    }
    if (reactionAtMInG != nullptr) {
        reactionAtMInG[0] = SpatialVec(Vec3(0), Vec3(0));
    }
}

auto RobotEngine::findMobilizerReactionOnBodyAtMInGround(const RobotModel& m, const RobotState& s, int body)
    -> SpatialVec {
    std::vector<SpatialVec> atM(static_cast<std::size_t>(m.numBodies), SpatialVec(Vec3(0), Vec3(0)));
    calcMobilizerReactionForces(m, s, nullptr, atM.data());
    if (body < 0 || body >= m.numBodies) {
        return SpatialVec(Vec3(0), Vec3(0));
    }
    return atM[static_cast<std::size_t>(body)];
}


// ============================================================================
//  TRANSFER: internal bodies -> per-atom Cartesian (accept path)
// ============================================================================
void RobotEngine::fillAtomPositionsFromBodies(const RobotModel& m, RobotState& s) {
    const Transform* X_GB = s.X_GB();
    Vec3* posG = s.atomPosG();
    // Atoms that are Cartesian-integrated inside the proposal (solvent-relaxing
    // NCMC) own their Ground positions directly in posG -- they are NOT placed by
    // any rigid body, so the body->atom fill must leave them untouched (otherwise
    // each step would snap them back onto the welded body and undo the relaxation).
    // mask == nullptr in the welded engine, so this is a no-op there.
    const char* cartMask = s.cartSolventMask();
    const int* __restrict atomBody = m.atomBody.data();
    const Vec3* __restrict station = m.atomStation_B.data();
    Vec3* __restrict posGr = posG;
    if (cartMask == nullptr) {
        // Welded engine (the common case): no Cartesian-solvent atoms to preserve, so
        // the whole loop vectorizes (see realizePosition for the same pattern).
#pragma omp simd
        for (int a = 0; a < m.numAtoms; ++a) {
            const int b = atomBody[a];
            posGr[a] = X_GB[b].p() + X_GB[b].R() * station[a];
        }
    } else {
        for (int a = 0; a < m.numAtoms; ++a) {
            if (cartMask[a]) {
                continue; // Cartesian-solvent atom owns posG directly; leave untouched
            }
            const int b = atomBody[a];
            posGr[a] = X_GB[b].p() + X_GB[b].R() * station[a];
        }
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


// ---- calcQDotDot: qddot = N qddot-coupling. Scalar joints: udot; quaternion ----
//      bodies (Ball/FreeLine/Free): quaternion second derivative. ----
void RobotEngine::calcQDotDot(const RobotModel& m, RobotState& s) {
    const Real* u = s.u();
    const Real* udot = s.udot();
    const Real* q = s.q();
    const Transform* X_FM = s.X_FM();
    Real* qdd = s.qdotdot();
    for (int b = 1; b < m.numBodies; ++b) {
        const int qOff = m.bodyQIndex[b];
        const int uOff = m.bodyUIndex[b];
        const int dof = m.bodyNU[b];
        const JointType jt = m.bodyJoint[b];
        if (jt == JointType::Ball) {
            const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            const Vec3 w(u[uOff], u[uOff + 1], u[uOff + 2]);
            const Vec3 wd(udot[uOff], udot[uOff + 1], udot[uOff + 2]);
            const Vec4 qdd4 = Rotation::convertAngVelDotToQuaternionDotDot(quat, w, wd);
            for (int k = 0; k < 4; ++k) {
                qdd[qOff + k] = qdd4[k];
            }
        } else if (jt == JointType::FreeLine) {
            // w_FM = R_FM*(u0,u1,0); wdot_FM = R_FM*(udot0,udot1,0) (the
            // R_FM_dot*(u0,u1,0) term is w_FM x w_FM = 0, see derivation).
            const Vec4 quat(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            const Rotation& R_FM = X_FM[b].R();
            const Vec3 w = R_FM * Vec3(u[uOff], u[uOff + 1], 0);
            const Vec3 wd = R_FM * Vec3(udot[uOff], udot[uOff + 1], 0);
            const Vec4 qdd4 = Rotation::convertAngVelDotToQuaternionDotDot(quat, w, wd);
            for (int k = 0; k < 4; ++k) {
                qdd[qOff + k] = qdd4[k];
            }
            qdd[qOff + 4] = udot[uOff + 2];
            qdd[qOff + 5] = udot[uOff + 3];
            qdd[qOff + 6] = udot[uOff + 4];
        } else if (jt == JointType::Free) {
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