#pragma once
// ============================================================================
//  HmcDriver.hpp -- a header-only, OpenMM-free generalized-coordinate HMC move,
//  for the ENSEMBLE / statistical-correctness tests.
//
//  WHY THIS EXISTS. The C++ test suite never instantiates World: World hard-wires
//  the OpenMM ForceBridge and keeps its HMC internals (reinitialize, metropolis,
//  currentTotalEnergy) private. The dynamics tests instead drive RobotEngine +
//  the templated stepTo against the OpenMM-free AnalyticForceBridge (see
//  TestIntegrator.cpp / TestNcmcWork.cpp). HmcDriver packages that same recipe as
//  ONE reusable move, so the ensemble tests (orientation Haar-uniformity, Fixman
//  Boltzmann, equipartition, detailed balance, mass-scaling invariance) all share
//  a single vetted GC-HMC kernel.
//
//  WHAT ONE move() DOES, mirroring World::reinitialize + the MD leg + metropolis:
//    1. momentum refreshment:  u = sqrt(RT) * sqrt(M^-1) * g, g ~ N(0,I)
//       (RobotEngine::multiplyBySqrtMInv -- M is never formed), then RATTLE-project
//       u if there are loop closures, EXACTLY as reinitialize does.
//    2. Hold = H(q,u) with the SAME terms World uses:
//          H = PE + KE [+ Fixman] [+ orientation Jacobian]
//       each correction independently toggleable (useFixman / useOrientationJac),
//       which is the whole point: the orientation-Jacobian study (Group A) needs to
//       run the term ON and OFF without touching production code.
//    3. mdSteps fixed-step velocity-Verlet steps (RobotEngine::stepTo).
//    4. Metropolis accept on dH with beta = 1/RT; on reject restore q (the momenta
//       are discarded and redrawn next move, so only q must be rolled back).
//
//  The Fixman term here drops the q-independent ln|M_3N| Cartesian reference World
//  carries: it is a constant and cancels identically in dH and in any marginal, so
//  omitting it changes nothing measurable while keeping the driver self-contained.
//
//  The orientation-Jacobian term reads the quaternion straight out of s.q() for
//  every quaternion joint block (model.quaternionQStart). At engine level -- a
//  single body attached to Ground, no setAtomsLocationsInGround / recomputeGeometry
//  block handoff -- s.q() IS the body's true orientation, so the freeRootAbsRotation
//  atom-triplet reconstruction World needs (to dodge its per-block frame reset) is
//  unnecessary; the math under test is identical.
// ============================================================================

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "Constraints.hpp"
#include "RobotEngine.hpp"
#include "RobotIntegrator.hpp" // template definitions of stepTo / verletStep
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "engine_helpers.hpp"
#include "robot_math.hpp"

namespace rtest {

template <class Bridge>
class HmcDriver {
    public:
    using Real = robo::Real;

    HmcDriver(const RobotModel& model,
              RobotState& state,
              Bridge& bridge,
              robo::ConstraintSet& constraints,
              Real RT,
              Real timeStep,
              int mdSteps,
              std::uint64_t seed = 0x5DEECE66DULL)
        : m_(model)
        , s_(state)
        , bridge_(bridge)
        , cs_(constraints)
        , RT_(RT)
        , h_(timeStep)
        , mdSteps_(mdSteps)
        , gen_(seed) {
    }

    // term selectors (default: plain PE+KE HMC)
    bool useFixman = false;
    bool useOrientationJac = false;

    // diagnostics
    [[nodiscard]] long attempted() const {
        return attempted_;
    }
    [[nodiscard]] long accepted() const {
        return accepted_;
    }
    [[nodiscard]] double acceptanceRate() const {
        return attempted_ ? double(accepted_) / double(attempted_) : 0.0;
    }

    Real gaussian() {
        return std::normal_distribution<Real>(Real(0), Real(1))(gen_);
    }
    Real uniform01() {
        return std::uniform_real_distribution<Real>(Real(0), Real(1))(gen_);
    }

    // u = sqrt(RT) * sqrt(M^-1) * g, g ~ N(0,I); RATTLE-project if constrained.
    // Mirrors World::reinitialize's momentum draw exactly.
    void seedMomenta() {
        RobotEngine::realizePosition(m_, s_);
        RobotEngine::realizeArticulatedBodyInertias(m_, s_);
        const int nu = m_.nu;
        std::vector<Real> g(static_cast<std::size_t>(nu)), seeded(static_cast<std::size_t>(nu));
        for (int i = 0; i < nu; ++i) {
            g[i] = gaussian();
        }
        RobotEngine::multiplyBySqrtMInv(m_, s_, g.data(), seeded.data());
        const Real scale = std::sqrt(RT_);
        Real* u = s_.u();
        for (int i = 0; i < nu; ++i) {
            u[i] = scale * seeded[i];
        }
        if (!cs_.empty()) {
            RobotEngine::realizeVelocity(m_, s_);
            cs_.enforceVelocityConstraints(m_, s_);
        }
    }

    // 1/2 RT ( ln|M_tree| - ln det(G M^-1 G^T) ); the constant ln|M_3N| is dropped.
    Real fixmanTerm() {
        RobotEngine::realizeArticulatedBodyInertias(m_, s_);
        const Real lnDetM = RobotEngine::calcLogDetM(m_, s_);
        const Real lnDetZ = cs_.calcConstraintLogDet(m_, s_);
        return Real(0.5) * RT_ * (lnDetM - lnDetZ);
    }

    // -1/2 RT sum_blocks ln sin^2(gamma2). The Euler/spherical orientation Jacobian
    // (off by default in production; see World::calcLogSineSqrGamma2). gamma2 pitch
    // sine = 2(w*y - z*x) from each quaternion block q = [w x y z].
    Real orientationJacTerm() const {
        double acc = 0.0;
        for (int qs : m_.quaternionQStart) {
            const Real w = s_.q()[qs + 0];
            const Real x = s_.q()[qs + 1];
            const Real y = s_.q()[qs + 2];
            const Real z = s_.q()[qs + 3];
            const Real sinPitch = std::clamp(Real(2.0) * ((w * y) - (z * x)), Real(-1.0), Real(1.0));
            acc += EngineHelpers::safeLogSineSqr(std::asin(sinPitch));
        }
        return Real(-0.5) * RT_ * acc;
    }

    // H = PE + KE [+ Fixman] [+ orientation Jacobian]. Refreshes positions/forces/
    // velocity from the current (q,u) first, so it is safe to call at start and end.
    Real hamiltonian() {
        RobotEngine::realizePosition(m_, s_);
        RobotEngine::fillAtomPositionsFromBodies(m_, s_);
        bridge_.evaluate(s_);
        const Real pe = bridge_.calcPotentialEnergy(s_);
        RobotEngine::realizeVelocity(m_, s_);
        const Real ke = RobotEngine::calcKineticEnergy(m_, s_);
        Real H = pe + ke;
        if (useFixman) {
            H += fixmanTerm();
        }
        if (useOrientationJac) {
            H += orientationJacTerm();
        }
        return H;
    }

    // Fill the qdot0/udot0/qdotdot0 the first verletStep reads (identical preamble
    // to TestIntegrator/TestNcmcWork seedDerivatives).
    void seedDerivatives() {
        RobotEngine::realizePosition(m_, s_);
        RobotEngine::fillAtomPositionsFromBodies(m_, s_);
        bridge_.evaluate(s_);
        RobotEngine::realizeVelocity(m_, s_);
        RobotEngine::realizeArticulatedBodyInertias(m_, s_);
        RobotEngine::calcUDot(m_, s_);
        RobotEngine::calcQDot(m_, s_, s_.qdot());
        RobotEngine::calcQDotDot(m_, s_);
    }

    // One full GC-HMC move. Returns true if accepted. On reject q is restored.
    bool move() {
        ++attempted_;
        seedMomenta();
        std::vector<Real> q0(s_.q(), s_.q() + m_.nq);
        std::vector<Real> u0(s_.u(), s_.u() + m_.nu);
        seedDerivatives();
        const Real Hold = hamiltonian();

        bool ok = true;
        for (int k = 0; k < mdSteps_ && ok; ++k) {
            ok = RobotEngine::stepTo(m_, s_, bridge_, cs_, s_.time + h_);
        }
        if (!ok) {
            restore(q0, u0);
            return false;
        }

        const Real Hnew = hamiltonian();
        const Real dH = Hnew - Hold;
        const bool accept = (dH <= Real(0)) || (uniform01() < std::exp(-dH / RT_));
        if (accept) {
            ++accepted_;
            return true;
        }
        restore(q0, u0);
        return false;
    }

    private:
    void restore(const std::vector<Real>& q0, const std::vector<Real>& u0) {
        std::copy(q0.begin(), q0.end(), s_.q());
        std::copy(u0.begin(), u0.end(), s_.u());
        RobotEngine::realizePosition(m_, s_);
        RobotEngine::fillAtomPositionsFromBodies(m_, s_);
    }

    const RobotModel& m_;
    RobotState& s_;
    Bridge& bridge_;
    robo::ConstraintSet& cs_;
    Real RT_;
    Real h_;
    int mdSteps_;
    std::mt19937_64 gen_;
    long attempted_ = 0;
    long accepted_ = 0;
};

} // namespace rtest
