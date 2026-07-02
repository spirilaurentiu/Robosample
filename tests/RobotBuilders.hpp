// ============================================================================
//  RobotBuilders.hpp -- hand-assemble RobotModel "robots" for engine tests,
//  WITHOUT the World/Context/OpenMM/molecule pipeline. A robot here is exactly
//  what Simbody/Robosample mean structurally: a tree of rigid bodies rooted on
//  the shared Ground (body 0). A molecular *system* is just a forest of such
//  trees sharing Ground -- so these builders let one test exercise several
//  robots at once, which is the whole point of "each molecule is a robot".
//
//  We populate only the fields the RobotEngine sweeps read: the body tree +
//  children CSR, the q/u layout (from RobotModel::jointNQ/jointNU, the single
//  source of truth), the static joint frames X_PF/X_BM, and per-body mass
//  properties. No atoms, no frame graph, no z-matrix -- those belong to the
//  geometry/transfer layer, not to the multibody mechanics under test here.
// ============================================================================
#pragma once

#include <vector>

#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"
#include "robot_math.hpp"

namespace rtest {

using robo::Rotation;
using robo::Transform;
using robo::UnitInertia;
using robo::Vec3;

// One body to add to a robot: its parent (already-added) body index and joint.
struct BodySpec {
    int parent; // parent body index (0 == Ground)
    JointType joint;
    Transform X_PF; // inboard frame on the parent
    Transform X_BM; // outboard frame on this body
    robo::Real mass = robo::Real(1);
    Vec3 com_B = Vec3(0); // COM in body frame
    UnitInertia inertia_B = UnitInertia(robo::Real(0.4), robo::Real(0.5), robo::Real(0.6));
};

// Assemble a RobotModel from a list of body specs (index i -> body i+1; body 0
// is Ground). Parents must be < child (topological order), which is natural
// since a spec can only name an already-listed parent.
inline RobotModel buildForest(const std::vector<BodySpec>& specs) {
    RobotModel m;
    const int B = static_cast<int>(specs.size()) + 1; // + Ground
    m.numBodies = B;
    m.numAtoms = 0;

    m.bodyParent.assign(B, -1);
    m.bodyLevel.assign(B, 0);
    m.bodyJoint.assign(B, JointType::Rigid);
    m.bodyRootAtom.assign(B, -1);
    m.X_PF.assign(B, Transform());
    m.X_BM.assign(B, Transform());
    m.bodyMass.assign(B, robo::Real(0));
    m.bodyCom_B.assign(B, Vec3(0));
    m.bodyUnitInertia_B.assign(B, UnitInertia(robo::Real(0), robo::Real(0), robo::Real(0)));

    m.bodyQIndex.assign(B, 0);
    m.bodyNQ.assign(B, 0);
    m.bodyUIndex.assign(B, 0);
    m.bodyNU.assign(B, 0);
    m.bodyUSqIndex.assign(B, 0);

    int qc = 0, uc = 0, usq = 0;
    for (int i = 0; i < (int)specs.size(); ++i) {
        const int b = i + 1;
        const BodySpec& s = specs[i];
        m.bodyParent[b] = s.parent;
        m.bodyLevel[b] = m.bodyLevel[s.parent] + 1;
        m.bodyJoint[b] = s.joint;
        m.X_PF[b] = s.X_PF;
        m.X_BM[b] = s.X_BM;
        m.bodyMass[b] = s.mass;
        m.bodyCom_B[b] = s.com_B;
        m.bodyUnitInertia_B[b] = s.inertia_B;

        const int nq = RobotModel::jointNQ(s.joint);
        const int nu = RobotModel::jointNU(s.joint);
        m.bodyQIndex[b] = qc;
        m.bodyNQ[b] = nq;
        m.bodyUIndex[b] = uc;
        m.bodyNU[b] = nu;
        m.bodyUSqIndex[b] = usq;
        if (RobotModel::jointUsesQuaternion(s.joint)) {
            m.quaternionQStart.push_back(qc);
        }
        qc += nq;
        uc += nu;
        usq += nu * nu;
    }
    m.nq = qc;
    m.nu = uc;
    m.nuSq = usq;

    // children CSR (needed by the inward ABI / udot sweeps)
    std::vector<std::vector<int>> kids(B);
    for (int b = 1; b < B; ++b) {
        kids[m.bodyParent[b]].push_back(b);
    }
    m.bodyChildrenBeg.assign(B, 0);
    m.bodyChildrenEnd.assign(B, 0);
    for (int b = 0; b < B; ++b) {
        m.bodyChildrenBeg[b] = (int)m.bodyChildren.size();
        for (int c : kids[b]) {
            m.bodyChildren.push_back(c);
        }
        m.bodyChildrenEnd[b] = (int)m.bodyChildren.size();
    }
    return m;
}

// Convenience: a single body of joint type `jt` attached straight to Ground,
// with random-but-fixed joint frames (so the test exercises non-trivial X_PF/
// X_BM, not just identity).
inline RobotModel buildSingle(JointType jt, Rng& rng) {
    BodySpec s;
    s.parent = 0;
    s.joint = jt;
    s.X_PF = Transform(rng.rotation(), rng.vec3());
    s.X_BM = Transform(rng.rotation(), rng.vec3());
    return buildForest({s});
}

// A single chain Ground -> joints[0] -> joints[1] -> ... (joints[0] must be a
// legal root). Every body gets random-but-fixed frames and mass. This is the
// Robosample analogue of Simbody TestAngleConversions' "one of every mobilizer
// in a chain" stress configuration.
inline RobotModel buildChain(const std::vector<JointType>& joints, Rng& rng) {
    std::vector<BodySpec> specs;
    for (int i = 0; i < (int)joints.size(); ++i) {
        BodySpec s;
        s.parent = i; // body i+1's parent is body i (Ground for the first)
        s.joint = joints[i];
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        s.mass = rng.uniform(0.5, 2.5);
        s.com_B = rng.vec3(-0.25, 0.25);
        s.inertia_B = UnitInertia(rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7), rng.uniform(0.3, 0.7));
        specs.push_back(s);
    }
    return buildForest(specs);
}

// ----------------------------------------------------------------------------
//  attachAtoms -- give a hand-built robot atoms (buildForest leaves numAtoms=0).
//  `stations_B` are in the BODY frame (the convention realizePosition consumes:
//  posG = X_GB.p + R_GB*station_B); `masses` are daltons (0 == virtual site).
//  ORDERING CONTRACT: call BEFORE RobotState::allocateFull(model) -- the state
//  sizes its per-atom arrays from model.numAtoms at allocation time.
//
//  (TestAtomTransfer.cpp has a sibling random/mass-deriving variant; this is the
//  explicit-placement form TestBuilders.cpp needs to check the station formula.)
// ----------------------------------------------------------------------------
inline void attachAtoms(RobotModel& m,
                        int body,
                        const std::vector<Vec3>& stations_B,
                        const std::vector<robo::Real>& masses) {
    if (body < 0 || body >= m.numBodies) {
        return;
    }
    if (stations_B.size() != masses.size()) {
        return;
    }
    const int firstNew = m.numAtoms;
    for (std::size_t k = 0; k < stations_B.size(); ++k) {
        m.atomBody.push_back(body);
        m.atomMass.push_back(masses[k]);
        m.atomStation_B.push_back(stations_B[k]);
        m.atomDummIndex.push_back(-1);
        m.atomCompoundIndex.push_back(m.numAtoms);
        ++m.numAtoms;
    }
    if (m.bodyRootAtom[body] < 0 && !stations_B.empty()) {
        m.bodyRootAtom[body] = firstNew;
    }
    m.bodyAtoms.clear();
    m.bodyAtomsBeg.assign(static_cast<std::size_t>(m.numBodies), 0);
    m.bodyAtomsEnd.assign(static_cast<std::size_t>(m.numBodies), 0);
    for (int b = 0; b < m.numBodies; ++b) {
        m.bodyAtomsBeg[b] = static_cast<int>(m.bodyAtoms.size());
        for (int a = 0; a < m.numAtoms; ++a) {
            if (m.atomBody[a] == b) {
                m.bodyAtoms.push_back(a);
            }
        }
        m.bodyAtomsEnd[b] = static_cast<int>(m.bodyAtoms.size());
    }
}

// Fill q with a valid random configuration: scalars uniform, quaternion blocks a
// random unit quaternion (Shoemake). u is plain Gaussian.
inline void randomizeState(const RobotModel& m, RobotState& s, Rng& rng) {
    robo::Real* q = s.q();
    for (int i = 0; i < m.nq; ++i) {
        q[i] = rng.uniform(robo::Real(-1.5), robo::Real(1.5));
    }
    for (int qs : m.quaternionQStart) {
        const auto quat = rng.unitQuat();
        q[qs + 0] = quat.elems[0];
        q[qs + 1] = quat.elems[1];
        q[qs + 2] = quat.elems[2];
        q[qs + 3] = quat.elems[3];
    }
    robo::Real* u = s.u();
    for (int i = 0; i < m.nu; ++i) {
        u[i] = rng.gaussian(0, robo::Real(0.7));
    }
}

// ----------------------------------------------------------------------------
//  buildBentTorsionChain -- a Free 6-DOF root followed by nTorsions Torsion
//  joints in a zigzag chain, alternating the bend axis every link (extends
//  TestFixmanBoltzmann.cpp::twoTorsionChain's 90-degree-bend coupling idea to
//  many links) so consecutive torsion axes are never parallel and the
//  mass-metric tensor M(phi) is dense / strongly off-diagonal rather than
//  block-diagonal. n_dof = 6 + nTorsions (RobotModel::nu). Shared by
//  TestEquipartition.cpp (T0.1) and TestEnsembleValidation.cpp (T0.2) so both
//  tiers exercise the SAME fixture (ensemble-validation spec 10-tier0
//  T0.1/T0.2). Two atoms are attached per body so downstream Cartesian
//  reconstruction (fillAtomPositionsFromBodies) has something to transform,
//  even on call sites that never read PE.
// ----------------------------------------------------------------------------
inline RobotModel buildBentTorsionChain(int nTorsions, Rng& rng) {
    std::vector<BodySpec> specs;

    BodySpec root;
    root.parent = 0;
    root.joint = JointType::Free;
    root.mass = rng.uniform(robo::Real(1.2), robo::Real(2.0));
    root.com_B = rng.vec3(robo::Real(-0.05), robo::Real(0.05));
    root.inertia_B = UnitInertia(rng.uniform(robo::Real(0.4), robo::Real(0.6)),
                                 rng.uniform(robo::Real(0.4), robo::Real(0.6)),
                                 rng.uniform(robo::Real(0.4), robo::Real(0.6)));
    specs.push_back(root);

    for (int i = 0; i < nTorsions; ++i) {
        BodySpec b;
        b.parent = i; // body i+1's parent is body i (root == body 1)
        b.joint = JointType::Torsion;
        const robo::CoordinateAxis bendAxis = (i % 2 == 0) ? robo::XAxis : robo::YAxis;
        b.X_PF = Transform(Rotation(robo::Real(M_PI_2), bendAxis), Vec3(robo::Real(0.15), robo::Real(0.0), robo::Real(0.0)));
        b.X_BM = Transform();
        b.mass = rng.uniform(robo::Real(0.8), robo::Real(1.6));
        b.com_B = Vec3(rng.uniform(robo::Real(0.08), robo::Real(0.14)),
                       rng.uniform(robo::Real(-0.05), robo::Real(0.05)),
                       rng.uniform(robo::Real(-0.05), robo::Real(0.05)));
        b.inertia_B = UnitInertia(rng.uniform(robo::Real(0.3), robo::Real(0.5)),
                                  rng.uniform(robo::Real(0.3), robo::Real(0.5)),
                                  rng.uniform(robo::Real(0.3), robo::Real(0.5)));
        specs.push_back(b);
    }

    RobotModel m = buildForest(specs);
    for (int b = 1; b < m.numBodies; ++b) {
        attachAtoms(m, b, {Vec3(0, 0, 0), Vec3(robo::Real(0.05), robo::Real(0.02), robo::Real(-0.01))},
                    {robo::Real(12.0), robo::Real(1.0)});
    }
    return m;
}

} // namespace rtest