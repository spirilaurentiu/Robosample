#pragma once
// ============================================================================
//  AnalyticForceBridge.hpp -- a test-only, OpenMM-free stand-in for ForceBridge
//  with a KNOWN closed-form potential, so integrator / energy-conservation /
//  ensemble tests have an analytic oracle.
//
//  Potential: a per-atom isotropic harmonic anchor
//      U(s)  = 1/2 k * sum_{a : mass>0} | r_a - anchor_a |^2
//      F_a   = -k (r_a - anchor_a)                 (0 for mass==0 virtual sites)
//  where r_a = s.atomPosG()[a] and anchor_a is captured once from the start
//  configuration s0 (so U(start) == 0). Conservative, smooth, and analytic:
//  -dU/dr_a == F_a exactly, which is what makes energy conservation testable.
//
//  CONTRACT MATCH. evaluate(RobotState&) has the SAME call shape verletStep
//  invokes on the real bridge, and reduces per-atom Cartesian forces to per-body
//  spatial forces with the IDENTICAL rule as ForceBridge::getForcesFromOpenMM:
//    * clear bodyForceG and mobilityForce every call;
//    * skip atomMass==0 virtual sites;
//    * BF[b].linear  += f;  BF[b].angular += (r_a - X_GB[b].p()) % f
//      (moment about the BODY ORIGIN, in Ground).
//  Keeping the reduction identical means a future change to the real reduction
//  is mirrored here on purpose, not missed by accident. This class is NOT named
//  ForceBridge -- it is linked alongside it and passed to the templated
//  verletStep/stepTo/checkReversibility.
//
//  PRECONDITION: the constructor and evaluate() read s.atomPosG()/s.X_GB(), so
//  RobotEngine::realizePosition(model, s) (which fills them) must have run for
//  the state passed in. The integrator already calls realizePosition before
//  bridge.evaluate, so this holds on the hot path.
// ============================================================================

#include <vector>

#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "robot_math.hpp"

class AnalyticForceBridge {
    public:
    using Real = robo::Real;
    using Vec3 = robo::Vec3;

    // Capture anchors from the (already position-realized) start state s0.
    AnalyticForceBridge(const RobotModel& model, const RobotState& s0, Real k)
        : model_(model)
        , k_(k)
        , anchor_(static_cast<std::size_t>(model.numAtoms), Vec3(0)) {
        const Vec3* p = s0.atomPosG();
        for (int a = 0; a < model_.numAtoms; ++a) {
            anchor_[static_cast<std::size_t>(a)] = p[a];
        }
    }

    [[nodiscard]] Real stiffness() const {
        return k_;
    }
    [[nodiscard]] const Vec3& anchor(int a) const {
        return anchor_[static_cast<std::size_t>(a)];
    }

    // ---- force / energy as functions of an EXPLICIT position array -----------
    // (used by the gradient test; the RobotState overloads below just forward).

    // F_a = -k (r_a - anchor_a); 0 for virtual sites (mass==0).
    [[nodiscard]] Vec3 atomForceAt(const Vec3* positions, int a) const {
        if (model_.atomMass[a] == Real(0)) {
            return Vec3(0);
        }
        return (anchor_[static_cast<std::size_t>(a)] - positions[a]) * k_; // = -k (r - anchor)
    }

    // U = 1/2 k sum_{real a} |r_a - anchor_a|^2 over the given positions.
    [[nodiscard]] Real calcPotentialEnergyAt(const Vec3* positions, int nAtoms) const {
        Real u = 0;
        for (int a = 0; a < nAtoms; ++a) {
            if (model_.atomMass[a] == Real(0)) {
                continue;
            }
            const Vec3 d = positions[a] - anchor_[static_cast<std::size_t>(a)];
            u += robo::dot(d, d);
        }
        return Real(0.5) * k_ * u;
    }

    // ---- RobotState convenience overloads ------------------------------------
    [[nodiscard]] Vec3 atomForce(const RobotState& s, int a) const {
        return atomForceAt(s.atomPosG(), a);
    }
    [[nodiscard]] Real calcPotentialEnergy(const RobotState& s) const {
        return calcPotentialEnergyAt(s.atomPosG(), model_.numAtoms);
    }

    // ---- the bridge contract verletStep calls --------------------------------
    // Same signature + same reduction as ForceBridge::getForcesFromOpenMM.
    void evaluate(RobotState& s) {
        robo::SpatialVec* BF = s.bodyForceG();
        for (int b = 0; b < model_.numBodies; ++b) {
            BF[b] = robo::SpatialVec(Vec3(0), Vec3(0));
        }
        Real* mob = s.mobilityForce();
        for (int i = 0; i < model_.nu; ++i) {
            mob[i] = Real(0);
        }
        const Vec3* posG = s.atomPosG();
        const robo::Transform* X_GB = s.X_GB();
        for (int a = 0; a < model_.numAtoms; ++a) {
            if (model_.atomMass[a] == Real(0)) {
                continue; // virtual site: skip (identical to the real bridge)
            }
            const int b = model_.atomBody[a];
            const Vec3 f = atomForceAt(posG, a);
            const Vec3 r = posG[a] - X_GB[b].p(); // station in Ground, about body origin
            BF[b][1] += f;                        // linear (force)
            BF[b][0] += r % f;                    // angular (moment about origin); % == cross
        }
    }

    private:
    const RobotModel& model_;
    Real k_;
    std::vector<Vec3> anchor_; // Ground-frame anchor per atom (start positions)
};