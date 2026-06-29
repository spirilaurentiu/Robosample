// ============================================================================
//  TestNMALinalg.cpp -- Phase 7: the OpenMM-free core of Route-B NMA.
//
//  computeRouteBNMA needs the OpenMM force singleton and is out of scope here,
//  but two pieces are pure and carry a LATENT INCONSISTENCY worth pinning: there
//  are now TWO cyclic-Jacobi symmetric eigensolvers in the codebase --
//      robo_linalg::jacobiSymEig   (tests/RobotLinearAlgebra.hpp; the engine's
//                                   hoisted hinge solver, tan-based rotation,
//                                   NO eigenvalue sort)
//      robo::detail::jacobiEigh    (include/NMA.hpp; atan2-based rotation,
//                                   eigenvalues sorted ASCENDING)
//  -- with different rotation formulas and different sort order. The NMA path
//  depends on eigval[0] being the softest mode, so the ascending sort is load-
//  bearing; and any future drift between the two conventions is exactly the bug
//  the duplication invites. We test jacobiEigh's contract directly and cross-
//  check the two solvers' spectra.
//
//  The third piece is the mass-weighting factor N = sqrt(M^-1): the header claims
//  N^T M N = I (it builds N column-by-column from multiplyBySqrtMInv(e_j), so
//  N N^T = M^-1 and hence N^T M N = I). We validate that against the Phase 4
//  dense M with ZERO OpenMM dependency.
//
//  OpenMM-FREE BUILD. NMA.hpp pulls OpenMMContext.hpp -> OpenMM.h for the force
//  evaluation, which these tests never touch. The test build supplies a minimal
//  OpenMM.h stub (declarations only; nothing linked), the same pattern the engine
//  TUs use with the tests/ForceBridge stub. detail::jacobiEigh reached this way
//  is the REAL production code, not a copy.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "NMA.hpp" // robo::detail::jacobiEigh (the OpenMM-free core)
#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotLinearAlgebra.hpp" // robo_linalg::jacobiSymEig (the other solver)
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::randomizeState;
using rtest::Rng;

namespace {

// nu = 6 coupled forest (Ball + 3 torsions): small enough for the n<=6 dense
// kernels and the same fixture the Phase 4 dense-M tests use.
RobotModel denseForest(Rng& rng) {
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
    return buildForest({F(0, JointType::Ball),
                        F(1, JointType::Torsion),
                        F(2, JointType::Torsion),
                        F(3, JointType::Torsion)});
}

// Dense nu x nu mass matrix from the system Jacobian (Phase 4 oracle): forward
// kinematics + body spatial inertias, no articulated-body inverse. Requires
// realizePosition current; restores u; uses realizeVelocity (writes V_GB).
std::vector<Real> buildDenseM(const RobotModel& m, RobotState& s) {
    const int n = m.nu, B = m.numBodies;
    const SpatialInertia* Mk = s.Mk_G();
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
            for (int b = 1; b < B; ++b) {
                const SpatialVec MV = Mk[b] * col[static_cast<std::size_t>(j)][static_cast<std::size_t>(b)];
                acc += dot(col[static_cast<std::size_t>(i)][static_cast<std::size_t>(b)], MV);
            }
            M[static_cast<std::size_t>(i * n + j)] = acc;
        }
    }
    return M;
}

} // namespace

// ---------------------------------------------------------------------------
//  1. JacobiEighReconstructsAndOrders: on random SPD A, detail::jacobiEigh gives
//     V diag(eval) V^T == A, V^T V == I (kLoose), and eigenvalues ASCENDING (the
//     NMA soft-mode path relies on eigval[0] being the softest).
// ---------------------------------------------------------------------------
TEST(NMALinalg, JacobiEighReconstructsAndOrders) {
    Rng rng(0x7001);
    for (int n : {1, 2, 3, 5, 6}) {
        for (int t = 0; t < 40; ++t) {
            Real flat[36];
            rng.spd(n, flat, 0.4, 4.0); // row-major SPD with eigenvalues in [0.4,4]
            std::vector<Real> A(flat, flat + n * n), eval, evec;
            robo::detail::jacobiEigh(A, n, eval, evec); // A passed by value (destroyed inside)

            // ascending eigenvalues
            for (int i = 1; i < n; ++i) {
                EXPECT_LE(eval[static_cast<std::size_t>(i - 1)],
                          eval[static_cast<std::size_t>(i)] + rtest::kLoose)
                    << "eigenvalues not ascending at i=" << i << " (n=" << n << ")";
            }

            // reconstruction V diag V^T == A and orthonormality V^T V == I
            Real maxRec = 0, maxOrth = 0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    Real rec = 0;
                    for (int k = 0; k < n; ++k) {
                        rec += evec[i * n + k] * eval[static_cast<std::size_t>(k)] * evec[j * n + k];
                    }
                    maxRec = std::max(maxRec, std::abs(rec - flat[i * n + j]));
                    Real orth = 0;
                    for (int k = 0; k < n; ++k) {
                        orth += evec[k * n + i] * evec[k * n + j];
                    }
                    maxOrth = std::max(maxOrth, std::abs(orth - (i == j ? Real(1) : Real(0))));
                }
            }
            EXPECT_LT(maxRec, rtest::kLoose) << "V diag V^T != A (n=" << n << ")";
            EXPECT_LT(maxOrth, rtest::kLoose) << "V^T V != I (n=" << n << ")";
        }
    }
}

// ---------------------------------------------------------------------------
//  2. TwoJacobiImplementationsAgree: the same SPD matrix through both solvers --
//     detail::jacobiEigh (sorted ascending) and robo_linalg::jacobiSymEig
//     (unsorted) -- gives the same SORTED eigenvalue spectrum to kLoose.
//     FAIL guard: if a future edit changes one solver's convention so the spectra
//     diverge, this goes red -- catching the drift the duplication invites.
// ---------------------------------------------------------------------------
TEST(NMALinalg, TwoJacobiImplementationsAgree) {
    Rng rng(0x7002);
    for (int n : {1, 2, 3, 5, 6}) {
        for (int t = 0; t < 40; ++t) {
            Real flat[36];
            rng.spd(n, flat, 0.3, 5.0);

            // detail::jacobiEigh -> already ascending
            std::vector<Real> A(flat, flat + n * n), evalJ, evecJ;
            robo::detail::jacobiEigh(A, n, evalJ, evecJ);

            // robo_linalg::jacobiSymEig -> unsorted; sort to compare
            Real d[6], V[36];
            robo_linalg::jacobiSymEig(flat, n, d, V);
            std::vector<Real> dl(d, d + n);
            std::sort(dl.begin(), dl.end());

            for (int i = 0; i < n; ++i) {
                EXPECT_NEAR(dl[static_cast<std::size_t>(i)],
                            evalJ[static_cast<std::size_t>(i)],
                            rtest::kLoose)
                    << "the two Jacobi solvers disagree on eigenvalue " << i << " (n=" << n << ")";
            }
        }
    }
}

// ---------------------------------------------------------------------------
//  3. MassWeightingIsOrthonormalizing: N = sqrt(M^-1) built column-by-column from
//     multiplyBySqrtMInv(e_j) satisfies the header's claim N^T M N = I, with M the
//     Phase 4 dense mass matrix. Validates the NMA mass-metric with no OpenMM.
// ---------------------------------------------------------------------------
TEST(NMALinalg, MassWeightingIsOrthonormalizing) {
    Rng rng(0x7003);
    RobotModel m = denseForest(rng);
    RobotState s;
    s.allocateFull(m);
    const int nu = m.nu;

    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);                // fills X_GB, H, Mk_G
        RobotEngine::realizeArticulatedBodyInertias(m, s); // P, G, DI for sqrtMInv

        // N column j = sqrt(M^-1) e_j (exactly the NMA construction).
        std::vector<Real> N(static_cast<std::size_t>(nu * nu), 0);
        {
            std::vector<Real> ej(static_cast<std::size_t>(nu), 0), col(static_cast<std::size_t>(nu), 0);
            for (int j = 0; j < nu; ++j) {
                std::fill(ej.begin(), ej.end(), Real(0));
                ej[static_cast<std::size_t>(j)] = 1;
                RobotEngine::multiplyBySqrtMInv(m, s, ej.data(), col.data());
                for (int i = 0; i < nu; ++i) {
                    N[static_cast<std::size_t>(i * nu + j)] = col[static_cast<std::size_t>(i)];
                }
            }
        }

        // dense M (its buildDenseM uses realizeVelocity as scratch -- fine, it does
        // not disturb Mk_G/X_GB which realizePosition above already fixed).
        const std::vector<Real> M = buildDenseM(m, s);

        // N^T M N
        std::vector<Real> MN(static_cast<std::size_t>(nu * nu), 0),
            NtMN(static_cast<std::size_t>(nu * nu), 0);
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nu; ++j) {
                Real acc = 0;
                for (int k = 0; k < nu; ++k) {
                    acc += M[i * nu + k] * N[k * nu + j];
                }
                MN[static_cast<std::size_t>(i * nu + j)] = acc;
            }
        }
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nu; ++j) {
                Real acc = 0;
                for (int k = 0; k < nu; ++k) {
                    acc += N[k * nu + i] * MN[k * nu + j];
                }
                NtMN[static_cast<std::size_t>(i * nu + j)] = acc;
            }
        }

        Real maxErr = 0;
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nu; ++j) {
                maxErr = std::max(maxErr, std::abs(NtMN[i * nu + j] - (i == j ? Real(1) : Real(0))));
            }
        }
        EXPECT_LT(maxErr, rtest::kLoose) << "N^T M N != I (rep " << rep << ", maxErr " << maxErr << ")";
    }
}