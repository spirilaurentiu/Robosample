// ============================================================================
//  GeometryFitter.cpp - World::setAtomsLocationsInGround / recomputeGeometry.
//
//  Relocated verbatim from World.cpp (SPLIT-W6, pure code motion): the
//  geometry re-fitting concern. setAtomsLocationsInGround is the Gibbs-block
//  continuation entry point that copies in per-atom Ground coordinates,
//  rebuilds rigid-body frames, resets q, and retargets loop-closure
//  distances; recomputeGeometry is its worker, which rebuilds every body's
//  mass properties and the X_PF/X_BM joint transforms from the incoming
//  Cartesian geometry.
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

using robo::Real;
using robo::Rotation;
using robo::Transform;
using robo::Vec3;

namespace {

inline void atomFrameFromGeometry(const Vec3* tgt,
                                  int self,
                                  int parent,
                                  int gparent,
                                  int refChild,
                                  Transform& outFrame,
                                  Transform* outB) {
    const Vec3& pChild = tgt[self];
    const Vec3& pParent = tgt[parent];
    const Vec3& pGparent = tgt[gparent];
    const Vec3& pRefChild = tgt[refChild];

    const Transform G(Rotation(robo::UnitVec3(pChild - pParent),
                               robo::XAxis,
                               robo::UnitVec3(pGparent - pParent),
                               robo::YAxis),
                      pParent);
    const Real theta = robo::calcDihedralAngle(pGparent, pParent, pChild, pRefChild);
    const Real d = (pChild - pParent).norm();

    const Transform B = Transform(Rotation(theta, robo::XAxis)) * Transform(Vec3(d, 0, 0))
                        * Transform(Rotation(robo::Pi, robo::YAxis));
    const Transform C(Rotation(robo::Pi, robo::XAxis));

    outFrame = G * B * C;
    if (outB != nullptr) {
        *outB = B;
    }
}

inline void
computeAllAtomFrames(const RobotModel::FrameGraph& fg, const Vec3* tgt, Transform* frame, Transform* xpcBc) {
    for (int k = 0; k < (int)fg.r_self.size(); ++k) {
        const int s = fg.r_self[k];
        frame[s] = Transform(Rotation(), tgt[s]);
        xpcBc[s] = Transform();
    }
    for (int k = 0; k < (int)fg.g_self.size(); ++k) {
        atomFrameFromGeometry(tgt,
                              fg.g_self[k],
                              fg.g_parent[k],
                              fg.g_gparent[k],
                              fg.g_refChild[k],
                              frame[fg.g_self[k]],
                              &xpcBc[fg.g_self[k]]);
    }
    for (int k = 0; k < (int)fg.f_self.size(); ++k) {
        const int s = fg.f_self[k];
        const int p = fg.f_parent[k];
        const Vec3 dir = tgt[s] - tgt[p];
        Rotation R;
        R.setRotationFromOneAxis(robo::UnitVec3(dir), robo::XAxis);
        frame[s] = Transform(R, tgt[s]);
        const Real d = dir.norm();
        xpcBc[s] = Transform(Vec3(d, 0, 0)) * Transform(Rotation(robo::Pi, robo::YAxis));
    }
}

} // namespace

// ----------------------------------------------------------------------------
//  setAtomsLocationsInGround
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
//  setAtomsLocationsInGround  -- the GIBBS-BLOCK CONTINUATION entry point.
//
//  This is how one Gibbs block hands the full configuration to the next. The
//  incoming atomPosG is the previous block's POST-Metropolis state (accepted ->
//  new point; rejected -> the point we stayed at). Because Cartesian coordinates
//  fully encode every bond length, bond angle, and torsion, copying them in and
//  rebuilding this world's internal frames from them carries ALL degrees of
//  freedom forward unchanged -- exactly the Gibbs requirement that each block
//  condition on the current values of the coordinates it does not itself sample.
//
//  Important: the bonds/angles this world holds rigid are NOT hard-constrained to
//  idealized values. recomputeGeometry() rebuilds every rigid body's shape from
//  the ACTUAL incoming coordinates (dst), so a coordinate frozen here sits at
//  whatever value it currently has -- and if some other block (e.g. a future
//  bond-length world) moves it, this block picks up the new value on its next
//  visit. The freeze is per-block (Gibbs conditioning), never a permanent
//  constraint. The q=0 / identity-quaternion reset does NOT discard geometry: it
//  only declares "the current configuration is this block's reference", with the
//  geometry living in the frames just rebuilt from dst. q=0 reconstructs the
//  incoming pose exactly (exact for a tree), which is also why a rejected move
//  restores the received geometry bit-for-bit.
// ----------------------------------------------------------------------------
void World::setAtomsLocationsInGround(const std::vector<robo::Vec3>& atomPosG) {
    Vec3* dst = state_.atomPosG();
    const int n = std::min<int>(model_.numAtoms, (int)atomPosG.size());
    for (int a = 0; a < n; ++a) {
        dst[a] = atomPosG[a];
    }
    if (cartesian_) {
        return;
    }

    recomputeGeometry(dst); // rigid-body shapes from the carried-over geometry
    std::fill(state_.q(), state_.q() + model_.nq, Real(0));
    for (int qs : model_.quaternionQStart) {
        state_.q()[qs] = Real(1);
    }
    std::fill(state_.u(), state_.u() + model_.nu, Real(0));
    RobotEngine::realizePosition(model_, state_);

    // Loop-closure target = the CARRIED-OVER distance, not the force-field r0.
    // The closure distance is part of the state being continued, so it must track
    // whatever the previous block left (continuation), exactly as the bonds and
    // angles above do. Do NOT reset this to the equilibrium bond length r0 -- that
    // would stamp a fixed idealized geometry over the carried-over configuration
    // and break Gibbs continuation. (Across torsional blocks the closure is held
    // fixed, so this value simply propagates; a block that genuinely samples the
    // ring region would update it through the carried-over Cartesian.)
    for (auto& c : constraints_.distance) {
        c.restLength = (dst[c.atomA] - dst[c.atomB]).norm();
    }
}

// ----------------------------------------------------------------------------
//  recomputeGeometry
// ----------------------------------------------------------------------------
void World::recomputeGeometry(const robo::Vec3* targets) {
    computeAllAtomFrames(model_.frameGraph, targets, frameFlat_.data(), xpcFlat_.data());

    static const Transform X_to_Z(Rotation(-90 * robo::Deg2Rad, robo::YAxis));

    for (int b = 1; b < model_.numBodies; ++b) {
        const int root = model_.bodyRootAtom[b];
        const Transform& T_X_B = frameFlat_[root];
        const Transform B_X_T = ~T_X_B;

        Real mass = 0;
        Vec3 com(0);
        robo::Inertia inertia(0);
        for (int ci = model_.bodyAtomsBeg[b]; ci < model_.bodyAtomsEnd[b]; ++ci) {
            const int a = model_.bodyAtoms[ci];
            const Vec3 station = (B_X_T * frameFlat_[a]).p();
            model_.atomStation_B[a] = station;
            const Real m = model_.atomMass[a];
            mass += m;
            com += m * station;
            inertia += robo::Inertia(station, m);
        }
        model_.bodyMass[b] = mass;
        const robo::MassProperties mp(mass, mass > 0 ? Vec3(com / mass) : Vec3(0), inertia);
        model_.bodyCom_B[b] = mp.getMassCenter();
        model_.bodyUnitInertia_B[b] = mp.getUnitInertia();

        const int p = model_.bodyParent[b];
        if (p == 0) {
            model_.X_PF[b] = T_X_B;
            model_.X_BM[b] = Transform();
        } else {
            // Per-joint axis convention, ported from molmodel calc_XPF_XBM_new:
            //   group A  (bond axis -> joint Z): Torsion, Translation, Free,
            //            Cylinder, Ball, FreeLine  -> apply X_to_Z.
            //   group B  (bond axis stays on X):  Slider, BendStretch,
            //            SphericalCoords           -> no axis switch.
            // (Weld never reaches here as a non-root flexible joint.)
            const JointType jt = model_.bodyJoint[b];
            const bool axisToZ =
                (jt == JointType::Torsion || jt == JointType::Cartesian || jt == JointType::Free
                 || jt == JointType::Cylinder || jt == JointType::Ball || jt == JointType::FreeLine);
            const Transform Proot_X_root = (~frameFlat_[model_.bodyRootAtom[p]]) * T_X_B;
            model_.X_BM[b] = axisToZ ? Transform(xpcFlat_[root] * X_to_Z) : xpcFlat_[root];
            model_.X_PF[b] = Proot_X_root * model_.X_BM[b];
        }
    }

    // atomStation_B was just refit -> the fused CUDA path must re-upload the stations on its
    // next evaluate (they are otherwise treated as resident/unchanged across steps).
    bridge_.markStationsDirty();
}
