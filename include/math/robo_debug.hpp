#pragma once

/**
 * @file robo_debug.hpp
 * @brief Development-only probes that scan realized @c RobotState for non-finite
 *        values and dump the offending body's coordinates.
 *
 * The probes are diagnostic: they read the model and state and write to
 * @c std::cout, but never mutate sampled coordinates, so including or excluding
 * this header does not change the sampled distribution. See @ref ROBO_DEBUG for
 * the switch semantics and the caveat that the switch is hard-defined on here.
 */

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>

#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "robot_math.hpp"

using robo::Real;
using robo::SpatialVec;
using robo::Transform;
using robo::Vec3;

/**
 * @def ROBO_DEBUG
 * @brief Master switch for the non-finite state probes.
 *
 * When nonzero, the @c robodbg helpers and the @ref ROBO_CHECK macro expand to
 * active scans; when zero, @ref ROBO_CHECK is a no-op and the probe functions
 * are not defined.
 * @warning This header hard-defines the switch to 1, so any translation unit
 *          that includes @c robo_debug.hpp compiles the probes in and cannot
 *          disable them from the command line (a @c -DROBO_DEBUG=0 would clash
 *          with this definition). The probes only read state and write
 *          diagnostics to @c std::cout; they never mutate sampled coordinates,
 *          so the sampled distribution is unaffected whether or not the header
 *          is included. "Compiled out" in production therefore means "the header
 *          is not included in that translation unit", not "the guard is off".
 */
#define ROBO_DEBUG 1
/**
 * @def ROBO_VERBOSE
 * @brief Intended verbosity level for the probes (0 = first-NaN dump only,
 *        1 = + per-step summary, 2 = + per-body transforms).
 * @note Hard-defined to 2 here and not read by any current probe; overriding it
 *       at compile time has no effect and clashes with this definition.
 */
#define ROBO_VERBOSE 2

#if ROBO_DEBUG
namespace robodbg {

/// Run counters the caller increments to tag diagnostics; process-global,
/// not thread-safe. @c evalCount is bumped once per full derivative evaluation,
/// @c stepCount once per Verlet step.
inline long evalCount = 0;
inline long stepCount = 0;
/// Latches true after the first non-finite value is reported, so the rich
/// context dump prints once per run and later detections print a one-liner.
inline bool firstNanDumped = false;

/**
 * @brief Finiteness predicate for a scalar.
 * @param[in] x Value to test.
 * @return @c true when @p x is neither NaN nor infinite.
 */
inline bool fin(Real x) {
    return std::isfinite(static_cast<double>(x));
}
/// @brief @c true when all three components of @p v are finite.
inline bool finV(const Vec3& v) {
    return fin(v[0]) && fin(v[1]) && fin(v[2]);
}
/// @brief @c true when both rows of the spatial vector @p s are finite.
inline bool finS(const SpatialVec& s) {
    return finV(s[0]) && finV(s[1]);
}

/**
 * @brief Print the q/u/udot/position context of one body to @c std::cout.
 * @param[in]  m Model providing this body's joint, parent, and coordinate offsets.
 * @param[in]  s State read for @p b's generalized coordinates and body origin.
 * @param[in]  b Body index to dump; out-of-range indices (@c b < 1 or
 *               @c b >= m.numBodies) print nothing.
 * @note Diagnostic output only; @p s is not modified.
 */
inline void dumpBody(const RobotModel& m, RobotState& s, int b) {
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

/**
 * @brief Scan the core realized state (q, u, udot, qdot, qdotdot, body origins,
 *        atom positions, body wrenches) for any non-finite value.
 * @param[in] m     Model supplying array sizes and the index->body mapping.
 * @param[in] s     State scanned; read only, never modified.
 * @param[in] where Phase label printed with the diagnostic (e.g. "calcUDot").
 * @return @c true when at least one scanned value is non-finite, @c false
 *         otherwise.
 * @post On the first non-finite detection of the run, prints a rich context
 *       block (phase, offending array and index, run counters, and the owning
 *       body via @ref dumpBody) and latches @ref firstNanDumped; on later
 *       detections prints a single-line notice. The only mutation is to the
 *       @ref firstNanDumped process-global flag and @c std::cout; sampled state
 *       is untouched.
 */
inline bool scan(const RobotModel& m, RobotState& s, const char* where) {
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
/**
 * @def ROBO_CHECK
 * @brief Scan the ambient state for non-finite values at a labelled phase.
 * @param where String literal naming the phase.
 *
 * Expands to a @ref robodbg::scan of the enclosing scope's @c m (a
 * @c RobotModel) and @c s (a @c RobotState) and discards the result, so it is a
 * statement usable anywhere those two names are in scope. Reduces to a no-op
 * when @ref ROBO_DEBUG is 0.
 */
#    define ROBO_CHECK(where) (void)robodbg::scan(m, s, (where))
#else
#    define ROBO_CHECK(where) ((void)0)
#endif
