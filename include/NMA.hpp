#pragma once

/**
 * @file NMA.hpp
 * @brief Route-B normal-mode analysis: the mass-weighted internal-coordinate
 *        Hessian and its softest non-rigid mode, used to bias HMC momenta.
 *
 * @ref robo::computeRouteBNMA builds the generalized stiffness
 * @f$ K_{ab} = \partial^2 V / \partial s_a \partial s_b @f$ by central
 * differencing the generalized force, mass-weights it into
 * @f$ \tilde H = N^\top K N @f$ with @f$ N = \sqrt{M^{-1}} @f$ (the same
 * articulated-body factor @c RobotEngine::multiplyBySqrtMInv applies, so
 * @f$ N^\top M N = I @f$), diagonalizes, and returns the softest non-trivial
 * mode as a unit-norm generalized-speed direction. It reads the model and state
 * and restores the state to the input configuration; it does not touch the
 * sampler, set momenta, or run acceptance.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "Constraints.hpp"   // ConstraintSet::{mapAtomForcesToGeneralizedForces, applyVelSpaceIncrementToQ}
#include "OpenMMContext.hpp" // force/energy singleton + OpenMM::Vec3
#include "RobotEngine.hpp"   // realizePosition, realizeArticulatedBodyInertias, multiplyBySqrtMInv
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "robot_math.hpp" // robo::Real, robo::Vec3

namespace robo {

/**
 * @brief Result of a Route-B normal-mode analysis over one configuration.
 *
 * A plain value bundle returned by @ref computeRouteBNMA. @c eigval holds the
 * mode frequencies @f$ \omega_k^2 @f$ in ascending order under the mass-weighted
 * u-metric; @c modeU[k] is the k-th mode in generalized-speed (u) space
 * (@f$ v_k = N y_k @f$); @c softMode indexes the chosen softest non-rigid mode;
 * @c nNearZero counts the rigid/external/noise modes below the zero tolerance
 * (about 6 per Free-rooted molecule); @c uScaleFactors is @c modeU[softMode]
 * renormalized to unit length, the direction the NMA momentum bias consumes.
 */
struct RouteBNMA {
    int nu = 0;
    int nNearZero = 0;                    // # eigenvalues treated as rigid/numerical-zero
    int softMode = -1;                    // index (in `eigval`) of the chosen soft mode
    std::vector<Real> eigval;             // omega^2, ASCENDING (mass-weighted u-metric)
    std::vector<std::vector<Real>> modeU; // modeU[k] = u-space mode v_k = N y_k (length nu)
    std::vector<Real> uScaleFactors;      // unit-norm direction of the soft mode -> NMA path
};

namespace detail {

/// @brief Copy the realized Ground-frame atom positions of @p s into @p posOut
///        as @c OpenMM::Vec3 (internal helper; @p s read only).
inline void buildOmmPositions(const RobotModel& m, const RobotState& s, std::vector<OpenMM::Vec3>& posOut) {
    const Vec3* p = s.atomPosG();
    posOut.resize(static_cast<std::size_t>(m.numAtoms));
    for (int a = 0; a < m.numAtoms; ++a) {
        posOut[a] = OpenMM::Vec3(p[a][0], p[a][1], p[a][2]);
    }
}

/**
 * @brief Generalized force @c g(q) at the current configuration (internal helper).
 * @pre @c RobotEngine::realizePosition(m, s) has run for the current @c q, so the
 *      atom positions the projection reads are valid.
 * @post @p gOut holds the @c nu-length generalized force. Virtual-site (mass==0)
 *       slots are zeroed before projection so their already-projected force is
 *       not double-counted (INV-2), matching the ForceBridge convention.
 */
inline void generalizedForce(const RobotModel& m,
                             const RobotState& s,
                             std::vector<OpenMM::Vec3>& posScratch,
                             std::vector<OpenMM::Vec3>& forceScratch,
                             std::vector<Vec3>& atomForce,
                             std::vector<Real>& gOut) {
    buildOmmPositions(m, s, posScratch);
    OpenMMContext::get().evaluateForcesFromPositionsCache(posScratch, forceScratch);

    atomForce.resize(static_cast<std::size_t>(m.numAtoms));
    for (int a = 0; a < m.numAtoms; ++a) {
        if (m.atomMass[a] == Real(0)) { // virtual site: parents already carry it
            atomForce[a] = Vec3(Real(0));
        } else {
            atomForce[a] = Vec3(forceScratch[a][0], forceScratch[a][1], forceScratch[a][2]);
        }
    }
    gOut.assign(static_cast<std::size_t>(m.nu), Real(0));
    ConstraintSet::mapAtomForcesToGeneralizedForces(m, s, atomForce.data(), gOut.data());
}

/**
 * @brief Symmetric eigensolver (cyclic Jacobi) with eigenvalues sorted ascending
 *        (internal helper).
 * @param[in]  A       Row-major symmetric @p n x @p n matrix, passed by value.
 * @param[in]  n       Dimension.
 * @param[out] evalOut Eigenvalues in ascending order (@p n entries).
 * @param[out] evecOut Eigenvectors as columns, row-major, reordered to match
 *                     @p evalOut (@c evecOut[r*n+k] is component @c r of the
 *                     @c k-th eigenvector).
 * @note Distinct from @ref robo::detail::jacobiSymEig (hinge_linalg), which
 *       leaves eigenvalues unsorted; this one sorts ascending because the mode
 *       selection depends on the ordering.
 */
inline void jacobiEigh(std::vector<Real> A, int n, std::vector<Real>& evalOut, std::vector<Real>& evecOut) {
    evecOut.assign(static_cast<std::size_t>(n) * n, Real(0));
    for (int i = 0; i < n; ++i) {
        evecOut[i * n + i] = Real(1);
    }

    const int maxSweeps = 100;
    for (int sweep = 0; sweep < maxSweeps; ++sweep) {
        Real off = 0;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                off += A[p * n + q] * A[p * n + q];
            }
        }
        if (off <= Real(1e-30)) {
            break;
        }

        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const Real apq = A[p * n + q];
                if (std::abs(apq) < Real(1e-300)) {
                    continue;
                }
                const Real app = A[p * n + p];
                const Real aqq = A[q * n + q];
                const Real phi = Real(0.5) * std::atan2(Real(2) * apq, aqq - app);
                const Real c = std::cos(phi), srot = std::sin(phi);
                for (int k = 0; k < n; ++k) {
                    const Real akp = A[k * n + p], akq = A[k * n + q];
                    A[k * n + p] = c * akp - srot * akq;
                    A[k * n + q] = srot * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {
                    const Real apk = A[p * n + k], aqk = A[q * n + k];
                    A[p * n + k] = c * apk - srot * aqk;
                    A[q * n + k] = srot * apk + c * aqk;
                }
                for (int k = 0; k < n; ++k) {
                    const Real vkp = evecOut[k * n + p], vkq = evecOut[k * n + q];
                    evecOut[k * n + p] = c * vkp - srot * vkq;
                    evecOut[k * n + q] = srot * vkp + c * vkq;
                }
            }
        }
    }

    std::vector<int> ord(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        ord[i] = i;
    }
    std::vector<Real> diag(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        diag[i] = A[i * n + i];
    }
    std::stable_sort(ord.begin(), ord.end(), [&](int x, int y) {
        return diag[x] < diag[y];
    });

    evalOut.assign(static_cast<std::size_t>(n), Real(0));
    std::vector<Real> evecSorted(static_cast<std::size_t>(n) * n, Real(0));
    for (int newc = 0; newc < n; ++newc) {
        const int oldc = ord[newc];
        evalOut[newc] = diag[oldc];
        for (int r = 0; r < n; ++r) {
            evecSorted[r * n + newc] = evecOut[r * n + oldc];
        }
    }
    evecOut.swap(evecSorted);
}

} // namespace detail


/**
 * @brief Compute the Route-B normal-mode analysis at the configuration currently
 *        held in @p state and return its softest non-rigid mode.
 * @param[in]     model   Immutable model (topology, masses, joint layout).
 * @param[in,out] state   State whose coordinates define the analysis point;
 *                        perturbed during finite differencing and restored to
 *                        the input configuration before return.
 * @param[in]     h       Finite-difference step along each generalized-speed
 *                        direction (radians for torsions).
 * @param[in]     zeroTol Eigenvalues with @c |omega^2| <= zeroTol are treated as
 *                        rigid/external/noise modes and skipped when selecting
 *                        the soft mode (about 6 near-zero per Free-rooted
 *                        molecule; 0 for an all-Rigid-rooted torsional world).
 * @return A @ref RouteBNMA with ascending @c eigval, the u-space modes, the
 *         chosen @c softMode, and the unit-norm @c uScaleFactors direction.
 * @pre OpenMM is initialized (post @c Context::initialize), @p model is a built
 *      internal-coordinate world (not Cartesian), and @p state holds the
 *      energy-minimized configuration @c q0.
 * @post @p state is left realized at the input @c q0; the analysis has no other
 *       observable effect on program state (it does not set the sampler's scale
 *       factors, draw momenta, or run acceptance). For @c nu == 0 the result is
 *       empty.
 * @note Cost is @c 2*nu OpenMM force evaluations plus @c nu sqrt(M^-1) sweeps.
 * @see RobotEngine::multiplyBySqrtMInv, VelocityDistortion
 */
inline RouteBNMA
computeRouteBNMA(const RobotModel& model, RobotState& state, Real h = Real(1e-5), Real zeroTol = Real(1e-6)) {
    const int nu = model.nu;
    RouteBNMA out;
    out.nu = nu;
    if (nu == 0) {
        return out;
    }

    // snapshot the minimized generalized coordinates
    std::vector<Real> q0(static_cast<std::size_t>(model.nq));
    {
        const Real* q = state.q();
        for (int i = 0; i < model.nq; ++i) {
            q0[i] = q[i];
        }
    }
    auto restoreQ0 = [&]() {
        Real* q = state.q();
        for (int i = 0; i < model.nq; ++i) {
            q[i] = q0[i];
        }
    };

    // scratch
    std::vector<OpenMM::Vec3> posScratch, forceScratch;
    std::vector<Vec3> atomForce;
    std::vector<Real> dvel(static_cast<std::size_t>(nu), Real(0));
    std::vector<Real> gPlus, gMinus;

    std::cout << "[RouteBNMA] computing mass-weighted internal-coordinate Hessian at q0..." << std::endl;

    // ---- 1. generalized stiffness K (column-by-column central difference) ----
    // K stored row-major nu x nu; column b is the response to perturbing s_b.
    std::vector<Real> K(static_cast<std::size_t>(nu) * nu, Real(0));
    for (int b = 0; b < nu; ++b) {
        restoreQ0();
        dvel[b] = +h;
        ConstraintSet::applyVelSpaceIncrementToQ(model, state, dvel.data());
        RobotEngine::realizePosition(model, state);
        detail::generalizedForce(model, state, posScratch, forceScratch, atomForce, gPlus);

        restoreQ0();
        dvel[b] = -h;
        ConstraintSet::applyVelSpaceIncrementToQ(model, state, dvel.data());
        RobotEngine::realizePosition(model, state);
        detail::generalizedForce(model, state, posScratch, forceScratch, atomForce, gMinus);

        dvel[b] = Real(0);
        const Real inv2h = Real(1) / (Real(2) * h);
        for (int a = 0; a < nu; ++a) {
            // g = -dV/ds  =>  K_ab = d2V/ds_a ds_b = -(gPlus - gMinus)/2h
            K[a * nu + b] = -(gPlus[a] - gMinus[a]) * inv2h;
        }
    }
    std::cout << "[RouteBNMA] ...done." << std::endl;

    // symmetrize (kills the antisymmetric part from numerical force noise)
    for (int i = 0; i < nu; ++i) {
        for (int j = i + 1; j < nu; ++j) {
            const Real avg = Real(0.5) * (K[i * nu + j] + K[j * nu + i]);
            K[i * nu + j] = avg;
            K[j * nu + i] = avg;
        }
    }

    std::cout << "[RouteBNMA] diagonalizing Htilde = N^T K N..." << std::endl;

    // ---- 2. mass-weighting factor N = sqrt(M^-1)  (N N^T = M^-1) -------------
    // Build N at the minimized configuration. N column j = sqrtMInv(e_j).
    restoreQ0();
    RobotEngine::realizePosition(model, state);
    RobotEngine::realizeArticulatedBodyInertias(model, state); // fills P, G, DI for sqrtMInv

    std::vector<Real> N(static_cast<std::size_t>(nu) * nu, Real(0));
    {
        std::vector<Real> ej(static_cast<std::size_t>(nu), Real(0));
        std::vector<Real> col(static_cast<std::size_t>(nu), Real(0));
        for (int j = 0; j < nu; ++j) {
            std::fill(ej.begin(), ej.end(), Real(0));
            ej[j] = Real(1);
            // NB: multiplyBySqrtMInv uses V_GB as scratch -- harmless here.
            RobotEngine::multiplyBySqrtMInv(model, state, ej.data(), col.data());
            for (int i = 0; i < nu; ++i) {
                N[i * nu + j] = col[i];
            }
        }
    }

    std::cout << "[RouteBNMA] ...done." << std::endl;

    // ---- 3. Htilde = N^T K N -------------------------------------------------
    std::vector<Real> KN(static_cast<std::size_t>(nu) * nu, Real(0)); // K @ N
    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j < nu; ++j) {
            Real acc = 0;
            for (int k = 0; k < nu; ++k) {
                acc += K[i * nu + k] * N[k * nu + j];
            }
            KN[i * nu + j] = acc;
        }
    }
    std::vector<Real> Htil(static_cast<std::size_t>(nu) * nu, Real(0)); // N^T @ KN
    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j < nu; ++j) {
            Real acc = 0;
            for (int k = 0; k < nu; ++k) {
                acc += N[k * nu + i] * KN[k * nu + j];
            }
            Htil[i * nu + j] = acc;
        }
    }
    for (int i = 0; i < nu; ++i) { // re-symmetrize against roundoff
        for (int j = i + 1; j < nu; ++j) {
            const Real avg = Real(0.5) * (Htil[i * nu + j] + Htil[j * nu + i]);
            Htil[i * nu + j] = avg;
            Htil[j * nu + i] = avg;
        }
    }

    std::cout << "[RouteBNMA] ...done." << std::endl;

    // ---- 4. diagonalize ------------------------------------------------------
    std::vector<Real> eval, evec;
    detail::jacobiEigh(Htil, nu, eval, evec); // eval ascending; evec col k = y_k
    out.eigval = eval;

    // u-space modes v_k = N y_k
    out.modeU.assign(static_cast<std::size_t>(nu), std::vector<Real>(static_cast<std::size_t>(nu), Real(0)));
    for (int k = 0; k < nu; ++k) {
        for (int i = 0; i < nu; ++i) {
            Real acc = 0;
            for (int j = 0; j < nu; ++j) {
                acc += N[i * nu + j] * evec[j * nu + k];
            }
            out.modeU[k][i] = acc;
        }
    }

    std::cout << "[RouteBNMA] ...done." << std::endl;

    // count near-zero (rigid/external/noise) modes and pick the softest above tol
    out.nNearZero = 0;
    for (int k = 0; k < nu; ++k) {
        if (std::abs(eval[k]) <= zeroTol) {
            ++out.nNearZero;
        }
    }
    out.softMode = -1;
    for (int k = 0; k < nu; ++k) {
        if (eval[k] > zeroTol) {
            out.softMode = k;
            break;
        }
    }
    if (out.softMode < 0) {
        out.softMode = nu - 1; // degenerate fallback
    }

    // ---- 5. hand-off vector for the NMA path (unit direction in u-space) -----
    out.uScaleFactors = out.modeU[static_cast<std::size_t>(out.softMode)];
    Real nrm = 0;
    for (Real x : out.uScaleFactors) {
        nrm += x * x;
    }
    nrm = std::sqrt(nrm);
    if (nrm > Real(0)) {
        for (Real& x : out.uScaleFactors) {
            x /= nrm;
        }
    }

    restoreQ0();
    RobotEngine::realizePosition(model, state); // leave caller's state at q0
    return out;
}

} // namespace robo