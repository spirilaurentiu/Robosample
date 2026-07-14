#pragma once
// ============================================================================
//  HarmonicBridge.hpp -- one OpenMM-free analytic-force-bridge template that
//  consolidates the six ad-hoc subclasses six different test files reinvented
//  (TESTS.md section 4): all are lambda/scaled variants of the SAME per-atom
//  isotropic harmonic anchor with the IDENTICAL per-atom -> per-body reduction
//  as AnalyticForceBridge.hpp/ForceBridge::getForcesFromOpenMM (clear
//  bodyForceG/mobilityForce every call; skip atomMass==0 virtual sites; moment
//  about the body origin, in Ground). `HarmonicBridge<Policy>::evaluate` runs
//  that ONE shared clear-then-apply sequence; each Policy supplies WHAT force
//  gets applied.
//
//  Policy contract (the per-atom-force policies -- everything except
//  PoisonForcePolicy): `Policy(const RobotModel&, const RobotState& s0,
//  <args...>)` captures anchors from the realized start state s0, exactly as
//  every subclass did; `Vec3 atomForce(const Vec3*, int) const` (0 for
//  atomMass==0); `Real energy(const RobotModel&, const RobotState&) const`;
//  `void apply(const RobotModel&, RobotState&)` runs applyPerAtomForce (below)
//  with *this, so the per-atom loop is written once and shared; `void
//  setLambda(Real)` on the lambda policies only -- a class-template member
//  function is only instantiated when a call site actually uses it, so a
//  policy without setLambda/energy never needs to define it.
//
//  PoisonForcePolicy is the one exception: InfForceBridge never ran a per-atom
//  loop -- it poisoned a single body force component directly. Routing it
//  through a per-atom force would ALSO poison that body's angular row (via r %
//  f) and any other atoms in the body, which the original never did. So
//  PoisonForcePolicy::apply writes the poison directly, matching bit-for-bit.
// ============================================================================

#include <limits>
#include <vector>

#include "../AnalyticForceBridge.hpp" // pulls RobotModel/RobotState/robot_math via its own includes

namespace rtest {

// Shared per-atom -> per-body reduction (clear() already ran in
// HarmonicBridge::evaluate). Reproduces AnalyticForceBridge.hpp:95-116 term-
// for-term: skip atomMass==0, moment about the body origin in Ground, and --
// when `cacheAtomForces` is true AND RobotState::wantsAtomForces() -- cache
// the raw per-atom force into atomForceG() (TwoRobotBridge's Cartesian-
// solvent cache contract, TestTwoRobotContact.cpp:250-260). Defaulting
// `cacheAtomForces` to false reproduces AnalyticForceBridge's behavior of
// never writing atomForceG(), regardless of wantsAtomForces().
template <class Policy>
void applyPerAtomForce(const RobotModel& model, RobotState& s, const Policy& policy, bool cacheAtomForces) {
    using Vec3 = robo::Vec3;
    using Real = robo::Real;
    const Vec3* posG = s.atomPosG();
    const robo::Transform* X_GB = s.X_GB();
    const bool cache = cacheAtomForces && s.wantsAtomForces();
    Vec3* atomForceOut = cache ? s.atomForceG() : nullptr;
    for (int a = 0; a < model.numAtoms; ++a) {
        if (model.atomMass[a] == Real(0)) {
            continue; // virtual site: skip (identical to the real bridge)
        }
        const int b = model.atomBody[a];
        const Vec3 f = policy.atomForce(posG, a);
        if (atomForceOut != nullptr) {
            atomForceOut[a] = f;
        }
        const Vec3 r = posG[a] - X_GB[b].p();
        s.bodyForceG()[b][1] += f; // linear
        s.bodyForceG()[b][0] += r % f; // angular (moment about origin)
    }
}

// NoForcePolicy -- k=0: clears forces, applies none. Replaces ZeroBridge
// (TestBiasForces.cpp) and the retired tests/ForceBridge.hpp no-op stub.
class NoForcePolicy {
    public:
    using Real = robo::Real;
    using Vec3 = robo::Vec3;

    NoForcePolicy(const RobotModel& /*model*/, const RobotState& /*s0*/) {}

    [[nodiscard]] Vec3 atomForce(const Vec3* /*pos*/, int /*a*/) const {
        return Vec3(0);
    }
    [[nodiscard]] Real energy(const RobotModel& /*model*/, const RobotState& /*s*/) const {
        return Real(0);
    }
    void apply(const RobotModel& /*model*/, RobotState& /*s*/) const {
        // Clear already ran in HarmonicBridge::evaluate; nothing further to apply.
    }
};

// SingleWellPolicy -- one isotropic harmonic anchor per atom, stiffness k.
// Replaces the plain AnalyticForceBridge-shaped uses and TwoRobotBridge's
// force term; `cacheAtomForces` (off by default) reproduces TwoRobotBridge's
// atomForceG() cache for the two-robot Cartesian-solvent suite.
class SingleWellPolicy {
    public:
    using Real = robo::Real;
    using Vec3 = robo::Vec3;

    SingleWellPolicy(const RobotModel& model, const RobotState& s0, Real k, bool cacheAtomForces = false)
        : model_(model)
        , k_(k)
        , cacheAtomForces_(cacheAtomForces)
        , anchor_(static_cast<std::size_t>(model.numAtoms)) {
        const Vec3* p = s0.atomPosG();
        for (int a = 0; a < model_.numAtoms; ++a) {
            anchor_[static_cast<std::size_t>(a)] = p[a];
        }
    }

    [[nodiscard]] Vec3 atomForce(const Vec3* pos, int a) const {
        if (model_.atomMass[a] == Real(0)) {
            return Vec3(0);
        }
        return (anchor_[static_cast<std::size_t>(a)] - pos[a]) * k_;
    }
    [[nodiscard]] Real energy(const RobotModel& /*model*/, const RobotState& s) const {
        const Vec3* pos = s.atomPosG();
        Real u = 0;
        for (int a = 0; a < model_.numAtoms; ++a) {
            if (model_.atomMass[a] == Real(0)) {
                continue;
            }
            const Vec3 d = pos[a] - anchor_[static_cast<std::size_t>(a)];
            u += robo::dot(d, d);
        }
        return Real(0.5) * k_ * u;
    }
    void apply(const RobotModel& model, RobotState& s) const {
        applyPerAtomForce(model, s, *this, cacheAtomForces_);
    }

    private:
    const RobotModel& model_;
    Real k_;
    bool cacheAtomForces_;
    std::vector<Vec3> anchor_;
};

// LambdaWellPolicy -- one anchor, setLambda scales the whole well: V = lambda
// * 1/2 k sum|r-anchor|^2. Replaces LambdaWellBridge (TestNcmcTeleport.cpp).
class LambdaWellPolicy {
    public:
    using Real = robo::Real;
    using Vec3 = robo::Vec3;

    LambdaWellPolicy(const RobotModel& model, const RobotState& s0, Real k)
        : model_(model)
        , k_(k)
        , anchor_(static_cast<std::size_t>(model.numAtoms)) {
        const Vec3* p = s0.atomPosG();
        for (int a = 0; a < model_.numAtoms; ++a) {
            anchor_[static_cast<std::size_t>(a)] = p[a];
        }
    }

    void setLambda(Real l) {
        lambda_ = l;
    }
    [[nodiscard]] Vec3 atomForce(const Vec3* pos, int a) const {
        if (model_.atomMass[a] == Real(0)) {
            return Vec3(0);
        }
        return (anchor_[static_cast<std::size_t>(a)] - pos[a]) * (k_ * lambda_);
    }
    [[nodiscard]] Real energy(const RobotModel& /*model*/, const RobotState& s) const {
        const Vec3* pos = s.atomPosG();
        Real u = 0;
        for (int a = 0; a < model_.numAtoms; ++a) {
            if (model_.atomMass[a] == Real(0)) {
                continue;
            }
            const Vec3 d = pos[a] - anchor_[static_cast<std::size_t>(a)];
            u += robo::dot(d, d);
        }
        return Real(0.5) * k_ * lambda_ * u;
    }
    void apply(const RobotModel& model, RobotState& s) const {
        applyPerAtomForce(model, s, *this, /*cacheAtomForces=*/false);
    }

    private:
    const RobotModel& model_;
    Real k_, lambda_ = 1.0;
    std::vector<Vec3> anchor_;
};

// IntraLambdaInterPolicy -- V = 1/2 kIntra sum|r-intraAnchor|^2 + lambda * 1/2
// kInter sum|r-interAnchor|^2, intraAnchor = start position, interAnchor =
// start position + (0.05,-0.03,0.04). Replaces LambdaAnalyticBridge
// (TestNCMCWork.cpp) and NcmcLambdaBridge (TestNcmcExplicitSolvent.cpp) -- the
// same two-anchor form under two names.
class IntraLambdaInterPolicy {
    public:
    using Real = robo::Real;
    using Vec3 = robo::Vec3;

    IntraLambdaInterPolicy(const RobotModel& model, const RobotState& s0, Real kIntra, Real kInter)
        : model_(model)
        , kIntra_(kIntra)
        , kInter_(kInter)
        , intraAnchor_(static_cast<std::size_t>(model.numAtoms))
        , interAnchor_(static_cast<std::size_t>(model.numAtoms)) {
        const Vec3* p = s0.atomPosG();
        for (int a = 0; a < model_.numAtoms; ++a) {
            intraAnchor_[static_cast<std::size_t>(a)] = p[a];
            interAnchor_[static_cast<std::size_t>(a)] = p[a] + Vec3(0.05, -0.03, 0.04);
        }
    }

    void setLambda(Real l) {
        lambda_ = l;
    }
    [[nodiscard]] Vec3 atomForce(const Vec3* pos, int a) const {
        if (model_.atomMass[a] == Real(0)) {
            return Vec3(0);
        }
        const Vec3 fi = (intraAnchor_[static_cast<std::size_t>(a)] - pos[a]) * kIntra_;
        const Vec3 fe = (interAnchor_[static_cast<std::size_t>(a)] - pos[a]) * (kInter_ * lambda_);
        return fi + fe;
    }
    [[nodiscard]] Real energy(const RobotModel& /*model*/, const RobotState& s) const {
        const Vec3* pos = s.atomPosG();
        Real ui = 0, ue = 0;
        for (int a = 0; a < model_.numAtoms; ++a) {
            if (model_.atomMass[a] == Real(0)) {
                continue;
            }
            const Vec3 di = pos[a] - intraAnchor_[static_cast<std::size_t>(a)];
            const Vec3 de = pos[a] - interAnchor_[static_cast<std::size_t>(a)];
            ui += robo::dot(di, di);
            ue += robo::dot(de, de);
        }
        return Real(0.5) * kIntra_ * ui + lambda_ * Real(0.5) * kInter_ * ue;
    }
    void apply(const RobotModel& model, RobotState& s) const {
        applyPerAtomForce(model, s, *this, /*cacheAtomForces=*/false);
    }

    private:
    const RobotModel& model_;
    Real kIntra_, kInter_, lambda_ = 1.0;
    std::vector<Vec3> intraAnchor_, interAnchor_;
};

// PoisonForcePolicy -- injects a non-finite body force to drive verletStep's
// reject-and-restore path without OpenMM. Replaces InfForceBridge
// (TestIntegrator.cpp). Does NOT go through applyPerAtomForce (see the
// file-level note): it writes body 1's linear[0] directly, exactly as the
// original did, leaving every other BF/mob entry at the clear()'d zero.
class PoisonForcePolicy {
    public:
    using Real = robo::Real;

    PoisonForcePolicy(const RobotModel& /*model*/, const RobotState& /*s0*/) {}

    void apply(const RobotModel& /*model*/, RobotState& s) const {
        s.bodyForceG()[1].linear[0] = std::numeric_limits<Real>::infinity();
    }
};

// HarmonicBridge<Policy> -- the ForceBridge-shaped contract (evaluate/
// calcPotentialEnergy/setLambda) every subclass reimplemented. evaluate() runs
// the ONE shared clear, then delegates to Policy::apply for what gets applied.
// calcPotentialEnergy/setLambda forward to the policy; per the file banner,
// they are only instantiated for policies whose call sites actually use them.
template <class Policy>
class HarmonicBridge {
    public:
    using Real = robo::Real;
    using Vec3 = robo::Vec3;

    template <class... PolicyArgs>
    HarmonicBridge(const RobotModel& model, const RobotState& s0, PolicyArgs&&... args)
        : model_(model)
        , policy_(model, s0, std::forward<PolicyArgs>(args)...) {}

    void setLambda(Real l) {
        policy_.setLambda(l);
    }
    [[nodiscard]] Real calcPotentialEnergy(const RobotState& s) const {
        return policy_.energy(model_, s);
    }
    void evaluate(RobotState& s) {
        robo::SpatialVec* BF = s.bodyForceG();
        for (int b = 0; b < model_.numBodies; ++b) {
            BF[b] = robo::SpatialVec(Vec3(0), Vec3(0));
        }
        Real* mob = s.mobilityForce();
        for (int i = 0; i < model_.nu; ++i) {
            mob[i] = Real(0);
        }
        policy_.apply(model_, s);
    }

    private:
    const RobotModel& model_;
    Policy policy_;
};

} // namespace rtest
