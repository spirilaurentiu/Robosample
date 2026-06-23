// ============================================================================
//  World.cpp - SimTK-free World.
//
//  CHANGES in this revision:
//   * boostMDSteps removed (was dead).
//   * Fixman potential + logSineSqrGamma2 wired into the acceptance Hamiltonian
//     of torsional worlds, summed over ALL bodies / ALL root bodies (the old
//     code applied logSineSqrGamma2 to molecule 0 only -- a multi-molecule bug).
//   * Rigid-kick docking move (RigidKick): a symmetric Cartesian proposal,
//     Metropolised on dU, optionally followed by MD-HMC relaxation.
//   * Root mobility is supplied per-world via buildModel's rootMobilities arg.
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>

using robo::Real;
using robo::Rotation;
using robo::Transform;
using robo::Vec3;

namespace {
constexpr double kBoltzmann_kJ = 0.0083144626; // kJ/mol/K

// Runtime toggle for per-kick docking diagnostics, read once. Set the env var
//   ROBO_DOCK_DEBUG=1
// to print, every kick, the sphere geometry and the sampled placement in GROUND
// (Cartesian) coordinates. All lengths below are nanometres (the engine's and
// OpenMM's native unit; AMBER Angstrom inputs are converted by ANG_TO_NM on load).
bool dockDebugEnabled() {
    static const bool on = [] {
        const char* e = std::getenv("ROBO_DOCK_DEBUG");
        return e != nullptr && e[0] != '0' && e[0] != '\0';
    }();
    return on;
}

struct DSU {
    std::vector<int> p;
    explicit DSU(int n)
        : p(n) {
        std::iota(p.begin(), p.end(), 0);
    }
    int find(int x) {
        return p[x] == x ? x : (p[x] = find(p[x]));
    }
    void join(int a, int b) {
        p[find(a)] = find(b);
    }
};

bool isFlexible(BondMobility m) {
    return m != BondMobility::Rigid;
}

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

// Rotation (row-major Mat33) -> unit quaternion (w,x,y,z). Shepperd's method.
inline void rotationToQuaternion(const Rotation& R, Real& w, Real& x, Real& y, Real& z) {
    const Real m00 = R(0, 0), m11 = R(1, 1), m22 = R(2, 2);
    const Real tr = m00 + m11 + m22;
    if (tr > 0) {
        Real s = std::sqrt(tr + 1.0) * 2.0;
        w = 0.25 * s;
        x = (R(2, 1) - R(1, 2)) / s;
        y = (R(0, 2) - R(2, 0)) / s;
        z = (R(1, 0) - R(0, 1)) / s;
    } else if (m00 > m11 && m00 > m22) {
        Real s = std::sqrt(1.0 + m00 - m11 - m22) * 2.0;
        w = (R(2, 1) - R(1, 2)) / s;
        x = 0.25 * s;
        y = (R(0, 1) + R(1, 0)) / s;
        z = (R(0, 2) + R(2, 0)) / s;
    } else if (m11 > m22) {
        Real s = std::sqrt(1.0 + m11 - m00 - m22) * 2.0;
        w = (R(0, 2) - R(2, 0)) / s;
        x = (R(0, 1) + R(1, 0)) / s;
        y = 0.25 * s;
        z = (R(1, 2) + R(2, 1)) / s;
    } else {
        Real s = std::sqrt(1.0 + m22 - m00 - m11) * 2.0;
        w = (R(1, 0) - R(0, 1)) / s;
        x = (R(0, 2) + R(2, 0)) / s;
        y = (R(1, 2) + R(2, 1)) / s;
        z = 0.25 * s;
    }
}

// log(sin^2(pitch)) with a floor that keeps it finite near the singularity.
inline Real safeLogSineSqr(Real pitch) {
    const Real s = std::sin(pitch);
    Real s2 = s * s;
    constexpr Real kFloor = 1e-12;
    if (s2 < kFloor) {
        s2 = kFloor;
    }
    return std::log(s2);
}

// Unit quaternion (w,x,y,z) -> rotation matrix (row-major).
inline Rotation quatToRotation(Real qw, Real qx, Real qy, Real qz) {
    const Real xx = qx * qx, yy = qy * qy, zz = qz * qz;
    const Real xy = qx * qy, xz = qx * qz, yz = qy * qz;
    const Real wx = qw * qx, wy = qw * qy, wz = qw * qz;
    return Rotation(robo::Mat33(1 - 2 * (yy + zz),
                                2 * (xy - wz),
                                2 * (xz + wy),
                                2 * (xy + wz),
                                1 - 2 * (xx + zz),
                                2 * (yz - wx),
                                2 * (xz - wy),
                                2 * (yz + wx),
                                1 - 2 * (xx + yy)));
}

} // namespace

World::World(int index, bool cartesian, std::uint32_t seed)
    : index_(index)
    , cartesian_(cartesian)
    , bridge_(model_)
    , rng_(static_cast<std::uint64_t>(seed) ^ (0x9e3779b97f4a7c15ULL * (index + 1))) {
}

World& World::add_sampler(double timeStep,
                          int mdSteps,
                          AcceptRejectMode mode,
                          bool useNuts,
                          double sphereFactor,
                          std::optional<bool> useFixman,
                          bool alwaysKick,
                          double clashThreshold,
                          int maxInitialKickTries) {
    sampler_.timeStep = timeStep;
    sampler_.mdSteps = mdSteps;
    sampler_.acceptRejectMode = mode;
    sampler_.useNuts = useNuts;
    sampler_.sphereFactor = (sphereFactor > 0.0) ? sphereFactor : 1.0;
    sampler_.alwaysKick = alwaysKick;
    sampler_.clashThreshold = (clashThreshold > 0.0) ? clashThreshold : 1.0e4;
    sampler_.maxInitialKickTries = (maxInitialKickTries > 0) ? maxInitialKickTries : 0;

    if (sampler_.timeStep == 0.0 || sampler_.mdSteps == 0) {
        // Pure proposal world: no internal dynamics. For docking this means a
        // rigid teleport only -- placement energy is the sole criterion.
        // timeStep=0 on a Cartesian world is the "OpenMM HMC reference" mode:
        // the docking world places the ligand, then the Cartesian world runs
        // its own integrator each round from that starting pose.
        sampler_.timeStep = 0.0;
        sampler_.mdSteps = 0;
        std::fprintf(stderr,
                     "[world %d] timeStep=0 / mdSteps=0: pure proposal mode "
                     "(no internal dynamics; %s)\n",
                     index_,
                     docking_ ? "rigid teleport only" : "OpenMM HMC reference mode");
    } else if (sampler_.timeStep > 0.005) {
        std::fprintf(stderr,
                     "[world %d] WARNING: timeStep=%.4g ps is large for all-atom MD; "
                     "the integrator may go NON-FINITE every step. "
                     "Typical stable values are 0.001-0.002 ps.\n",
                     index_,
                     sampler_.timeStep);
    }
    // AUTO default: Fixman + logSineSqr ON for non-Cartesian (torsional + docking)
    // worlds, OFF for Cartesian (flat space => constant metric => no correction).
    sampler_.useFixman = useFixman.has_value() ? (*useFixman && !cartesian_) : !cartesian_;
    sampler_.moveType = docking_ ? MoveType::RigidKick : MoveType::MdHmc;
    return *this;
}

void World::configureDocking(std::vector<std::vector<int>> ligandGroups, std::vector<int> siteAtoms) {
    docking_ = true;
    ligandGroups_ = std::move(ligandGroups);
    siteAtoms_ = std::move(siteAtoms);
    sampler_.moveType = MoveType::RigidKick;
}

void World::setTemperature(double T) {
    temperature_ = T;
    RT_ = kBoltzmann_kJ * T;
    beta_ = (RT_ > 0) ? 1.0 / RT_ : 0.0;
}

// ----------------------------------------------------------------------------
//  Kinetic-metric preconditioning (fictitious mass; sampling only).
//  Scales the spatial inertia used ONLY in the proposal (momentum draw, KE,
//  Fixman ln det M) -> raises the stable dt ~sqrt(scale) for the scaled body
//  with ZERO configurational bias (RobotModel::bodyMassScale documents why).
//  Call AFTER buildModel(). Re-sizes defensively if the model array is missing
//  (e.g. a Cartesian world that took the compact build path).
// ----------------------------------------------------------------------------
void World::setBodyMassScale(int body, double scale) {
    if ((int)model_.bodyMassScale.size() != model_.numBodies) {
        model_.bodyMassScale.assign(model_.numBodies, robo::Real(1));
    }
    if (body >= 1 && body < model_.numBodies) {
        model_.bodyMassScale[body] = robo::Real(scale);
    }
}

void World::setMassScaleByJoint(JointType jt, double scale) {
    if ((int)model_.bodyMassScale.size() != model_.numBodies) {
        model_.bodyMassScale.assign(model_.numBodies, robo::Real(1));
    }
    for (int b = 1; b < model_.numBodies; ++b) {
        if (model_.bodyJoint[b] == jt) {
            model_.bodyMassScale[b] = robo::Real(scale);
        }
    }
}

void World::setReversibilityCheck(int interval) {
    sampler_.reversibilityCheckInterval = (interval > 0) ? interval : 0;
    if (sampler_.reversibilityCheckInterval > 0) {
        // Typically set from add_*_world() BEFORE add_sampler(), so mdSteps and
        // timeStep are not known here yet; the per-round "[rev]" line logs the
        // actual step count and dt at run time.
        std::fprintf(stderr,
                     "[world %d] reversibility probe ENABLED: every %d round(s) "
                     "(non-destructive; residual logged per round; see THEORY 5.7).\n",
                     index_,
                     sampler_.reversibilityCheckInterval);
    }
}

// ----------------------------------------------------------------------------
//  buildModel  (unchanged structure; stores ln|M_3N| for the Fixman reference)
// ----------------------------------------------------------------------------
void World::buildModel(const SystemTopology& sys,
                       const Selection& sel,
                       const std::vector<RootMobility>& rootMobilities) {
    const int nAtoms = sys.numAtoms;
    model_.numAtoms = nAtoms;

    if (cartesian_) {
        model_.numBodies = 1;
        model_.numZRows = 0;
        model_.atomBody.assign(nAtoms, 0);
        state_.allocateCompact(model_);
        return;
    }

    DSU dsu(nAtoms);
    for (int k = 0; k < sys.numBonds; ++k) {
        if (sys.bondsRingClosing[k]) {
            continue;
        }
        const BondMobility mob =
            (k < (int)sel.bondMobility.size()) ? sel.bondMobility[k] : BondMobility::Rigid;
        if (!isFlexible(mob)) {
            dsu.join(sys.bondsI[k], sys.bondsJ[k]);
        }
    }

    std::vector<int> compToBody(nAtoms, -1);
    int nextBody = 1;
    std::vector<int> atomBody(nAtoms, 0);
    for (int a = 0; a < nAtoms; ++a) {
        const int c = dsu.find(a);
        if (compToBody[c] < 0) {
            compToBody[c] = nextBody++;
        }
        atomBody[a] = compToBody[c];
    }
    const int B = nextBody;
    model_.numBodies = B;
    model_.atomBody = atomBody;
    model_.atomMass.assign(nAtoms, 0);
    for (int a = 0; a < nAtoms; ++a) {
        model_.atomMass[a] = sys.atomsMass[a];
    }

    struct Edge {
        int bi, bj, ai, aj;
    };
    std::vector<Edge> jointEdges;
    for (int k = 0; k < sys.numBonds; ++k) {
        if (sys.bondsRingClosing[k]) {
            continue;
        }
        const BondMobility mob =
            (k < (int)sel.bondMobility.size()) ? sel.bondMobility[k] : BondMobility::Rigid;
        if (isFlexible(mob)) {
            const int ai = sys.bondsI[k], aj = sys.bondsJ[k];
            jointEdges.push_back({atomBody[ai], atomBody[aj], ai, aj});
        }
    }

    std::vector<std::vector<int>> adj(B);
    for (int e = 0; e < (int)jointEdges.size(); ++e) {
        adj[jointEdges[e].bi].push_back(e);
        adj[jointEdges[e].bj].push_back(e);
    }

    std::vector<int> rootBodyOfMol;
    std::vector<RootMobility> rootMobOfBody(B, RootMobility::Weld);
    std::vector<int> rootAtomOfBody(B, -1);
    for (int mol = 0; mol < sys.numMolecules; ++mol) {
        const int rootAtom = sys.atomsBegin[mol];
        const int rb = atomBody[rootAtom];
        rootBodyOfMol.push_back(rb);
        rootMobOfBody[rb] = (mol < (int)rootMobilities.size()) ? rootMobilities[mol] : RootMobility::Weld;
        rootAtomOfBody[rb] = rootAtom;
    }

    model_.bodyParent.assign(B, -1);
    model_.bodyLevel.assign(B, 0);
    model_.bodyJoint.assign(B, JointType::Weld);
    model_.bodyRootAtom.assign(B, -1);
    std::vector<bool> placed(B, false);
    placed[0] = true;
    std::queue<int> bfs;

    auto rootJointType = [](RootMobility rm) {
        switch (rm) {
            case RootMobility::Free:
                return JointType::Free;
            case RootMobility::Cartesian:
                return JointType::Translation;
            case RootMobility::Pin:
                return JointType::Pin;
            case RootMobility::Weld:
            default:
                return JointType::Weld;
        }
    };

    for (int rb : rootBodyOfMol) {
        if (placed[rb]) {
            continue;
        }
        model_.bodyParent[rb] = 0;
        model_.bodyLevel[rb] = 1;
        model_.bodyJoint[rb] = rootJointType(rootMobOfBody[rb]);
        model_.bodyRootAtom[rb] = rootAtomOfBody[rb];
        placed[rb] = true;
        bfs.push(rb);
        while (!bfs.empty()) {
            const int u = bfs.front();
            bfs.pop();
            for (int e : adj[u]) {
                const Edge& ed = jointEdges[e];
                const int v = (ed.bi == u) ? ed.bj : ed.bi;
                if (placed[v]) {
                    continue;
                }
                model_.bodyParent[v] = u;
                model_.bodyLevel[v] = model_.bodyLevel[u] + 1;
                model_.bodyJoint[v] = JointType::Pin;
                model_.bodyRootAtom[v] = (ed.bi == u) ? ed.aj : ed.ai;
                placed[v] = true;
                bfs.push(v);
            }
        }
    }

    // Relabel bodies into topological (parent<child) order.
    {
        std::vector<int> order;
        order.reserve(B > 1 ? B - 1 : 0);
        for (int b = 1; b < B; ++b) {
            order.push_back(b);
        }
        std::stable_sort(order.begin(), order.end(), [&](int x, int y) {
            if (model_.bodyLevel[x] != model_.bodyLevel[y]) {
                return model_.bodyLevel[x] < model_.bodyLevel[y];
            }
            return x < y;
        });
        std::vector<int> newId(B, 0);
        for (int i = 0; i < (int)order.size(); ++i) {
            newId[order[i]] = i + 1;
        }
        std::vector<int> np(B, -1), nlvl(B, 0), nra(B, -1);
        std::vector<JointType> njt(B, JointType::Weld);
        for (int b = 1; b < B; ++b) {
            const int nb = newId[b];
            const int op = model_.bodyParent[b];
            np[nb] = (op <= 0) ? 0 : newId[op];
            nlvl[nb] = model_.bodyLevel[b];
            njt[nb] = model_.bodyJoint[b];
            nra[nb] = model_.bodyRootAtom[b];
        }
        model_.bodyParent.swap(np);
        model_.bodyLevel.swap(nlvl);
        model_.bodyJoint.swap(njt);
        model_.bodyRootAtom.swap(nra);
        for (int a = 0; a < nAtoms; ++a) {
            atomBody[a] = newId[atomBody[a]];
        }
        model_.atomBody = atomBody;
    }

    std::vector<std::vector<int>> kids(B);
    for (int b = 1; b < B; ++b) {
        if (model_.bodyParent[b] >= 0) {
            kids[model_.bodyParent[b]].push_back(b);
        }
    }
    model_.bodyChildrenBeg.assign(B, 0);
    model_.bodyChildrenEnd.assign(B, 0);
    for (int b = 0; b < B; ++b) {
        model_.bodyChildrenBeg[b] = (int)model_.bodyChildren.size();
        for (int c : kids[b]) {
            model_.bodyChildren.push_back(c);
        }
        model_.bodyChildrenEnd[b] = (int)model_.bodyChildren.size();
    }

    auto nqOf = [](JointType jt) {
        switch (jt) {
            case JointType::Weld:
                return 0;
            case JointType::Pin:
                return 1;
            case JointType::Translation:
                return 3;
            case JointType::Free:
                return 7;
            default:
                return 0;
        }
    };
    auto nuOf = [](JointType jt) {
        switch (jt) {
            case JointType::Weld:
                return 0;
            case JointType::Pin:
                return 1;
            case JointType::Translation:
                return 3;
            case JointType::Free:
                return 6;
            default:
                return 0;
        }
    };
    model_.bodyQIndex.assign(B, 0);
    model_.bodyNQ.assign(B, 0);
    model_.bodyUIndex.assign(B, 0);
    model_.bodyNU.assign(B, 0);
    model_.bodyUSqIndex.assign(B, 0);
    int qc = 0, uc = 0, usq = 0;
    for (int b = 1; b < B; ++b) {
        const JointType jt = model_.bodyJoint[b];
        const int nq = nqOf(jt), nu = nuOf(jt);
        model_.bodyQIndex[b] = qc;
        model_.bodyNQ[b] = nq;
        model_.bodyUIndex[b] = uc;
        model_.bodyNU[b] = nu;
        model_.bodyUSqIndex[b] = usq;
        if (jt == JointType::Free) {
            model_.quaternionQStart.push_back(qc);
        }
        qc += nq;
        uc += nu;
        usq += nu * nu;
    }
    model_.nq = qc;
    model_.nu = uc;
    model_.nuSq = usq;
    model_.numZRows = 0;

    model_.bodyAtomsBeg.assign(B, 0);
    model_.bodyAtomsEnd.assign(B, 0);
    {
        std::vector<std::vector<int>> ba(B);
        for (int a = 0; a < nAtoms; ++a) {
            ba[atomBody[a]].push_back(a);
        }
        for (int b = 0; b < B; ++b) {
            model_.bodyAtomsBeg[b] = (int)model_.bodyAtoms.size();
            for (int a : ba[b]) {
                model_.bodyAtoms.push_back(a);
            }
            model_.bodyAtomsEnd[b] = (int)model_.bodyAtoms.size();
        }
    }

    {
        auto& fg = model_.frameGraph;
        fg = RobotModel::FrameGraph{};
        fg.topoOffset = {0, nAtoms};
        fg.totalAtoms = nAtoms;
        std::vector<int> par(nAtoms, -1), firstChild(nAtoms, -1);
        for (int k = 0; k < sys.numBonds; ++k) {
            if (sys.bondsRingClosing[k]) {
                continue;
            }
            const int p = sys.bondsI[k];
            const int ch = sys.bondsJ[k];
            par[ch] = p;
            if (firstChild[p] < 0 || ch < firstChild[p]) {
                firstChild[p] = ch;
            }
        }
        for (int a = 0; a < nAtoms; ++a) {
            if (par[a] < 0) {
                fg.r_self.push_back(a);
                continue;
            }
            const int gp = par[par[a]];
            const int rc = firstChild[a];
            if (gp >= 0 && rc >= 0) {
                fg.g_self.push_back(a);
                fg.g_parent.push_back(par[a]);
                fg.g_gparent.push_back(gp);
                fg.g_refChild.push_back(rc);
            } else {
                fg.f_self.push_back(a);
                fg.f_parent.push_back(par[a]);
            }
        }
    }

    frameFlat_.assign(nAtoms, Transform());
    xpcFlat_.assign(nAtoms, Transform());
    model_.X_PF.assign(B, Transform());
    model_.X_BM.assign(B, Transform());
    model_.bodyMass.assign(B, 0);
    model_.bodyCom_B.assign(B, Vec3(0));
    model_.bodyUnitInertia_B.assign(B, robo::UnitInertia(0, 0, 0));
    model_.bodyMassScale.assign(B, robo::Real(1)); // kinetic-metric preconditioning; 1.0 = physical
    model_.atomStation_B.assign(nAtoms, Vec3(0));

    state_.allocateFull(model_);

    std::vector<Vec3> ref(nAtoms);
    for (int a = 0; a < nAtoms; ++a) {
        ref[a] = Vec3(sys.atomsX[a], sys.atomsY[a], sys.atomsZ[a]);
    }
    std::copy(ref.begin(), ref.end(), state_.atomPosG());
    recomputeGeometry(state_.atomPosG());
    std::fill(state_.q(), state_.q() + model_.nq, Real(0));
    for (int qs : model_.quaternionQStart) {
        state_.q()[qs] = Real(1);
    }
    std::fill(state_.u(), state_.u() + model_.nu, Real(0));
    RobotEngine::realizePosition(model_, state_);
    RobotEngine::fillAtomPositionsFromBodies(model_, state_);

    // Constant Cartesian mass-matrix log-determinant: ln|M_3N| = 3 * sum ln m_a.
    lnDetMCartesian_ = 0.0;
    for (int a = 0; a < nAtoms; ++a) {
        if (model_.atomMass[a] > 0) {
            lnDetMCartesian_ += 3.0 * std::log(model_.atomMass[a]);
        }
    }

    constraints_.distance.clear();
    for (int k = 0; k < sys.numBonds; ++k) {
        if (!sys.bondsRingClosing[k]) {
            continue;
        }
        const int ai = sys.bondsI[k];
        const int aj = sys.bondsJ[k];
        if (model_.atomBody[ai] != model_.atomBody[aj]) {
            constraints_.distance.push_back(robo::DistanceConstraint{ai, aj, Real(0)});
        }
    }
}

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
            const Transform Proot_X_root = (~frameFlat_[model_.bodyRootAtom[p]]) * T_X_B;
            model_.X_BM[b] = xpcFlat_[root] * X_to_Z;
            model_.X_PF[b] = Proot_X_root * model_.X_BM[b];
        }
    }
}

// ----------------------------------------------------------------------------
//  Fixman / coordinate-Jacobian corrections (torsional worlds)
// ----------------------------------------------------------------------------
double World::calcFixman() {
    // Fixman compensating potential, Spiridon & Minh 2017 (JCTC 13:4649) Eq. 3:
    //   U_F = (1/2) RT ln( |M_{N_f}| / |M_3N| )
    // where |M_{N_f}| is the mass-metric determinant of the constrained system's
    // FLEXIBLE coordinates and |M_3N| (constant) is the Cartesian reference. U_F
    // makes the constrained-move marginal (their Eq. 2, rho ~ |M_{N_f}|^{1/2}
    // e^{-bU}) match the unconstrained Cartesian Boltzmann marginal.
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);

    // (a) Tree term: ln|M_tree| = sum_b ln det(D_b), the O(n) articulated-body
    //     determinant (Jain et al., refs 9 & 23). For an ACYCLIC molecule the
    //     flexible coordinates ARE the tree torsions, so |M_{N_f}| = |M_tree| and
    //     this term alone is exact -- the regime the paper validated.
    const double lnDetM = RobotEngine::calcLogDetM(model_, state_);

    // (b) Loop-closure term: for a CYCLIC molecule the ring is opened into the
    //     spanning tree and closed by a RATTLE distance constraint, which removes
    //     one flexible DOF per ring. The RATTLE momentum projection (G M^-1 p = 0)
    //     contributes det(G M^-1 G^T)^{-1/2} to the marginal, so the correct
    //     flexible determinant is |M_{N_f}| = |M_tree| / det(G M^-1 G^T). The Jain
    //     tree algorithm does not include this (it is a branched-molecule method);
    //     the paper never tested ring Boltzmann correctness (its macrocycle data,
    //     Sec. 3.5, measured efficiency only), so the term was simply missing for
    //     cyclic systems. calcConstraintLogDet returns 0 when there are no loop
    //     closures, making this a guaranteed no-op on every acyclic system.
    const double lnDetZ = constraints_.calcConstraintLogDet(model_, state_);

    // U_F = 1/2 RT ( ln|M_tree| - ln det(G M^-1 G^T) - ln|M_3N| ).
    // = 1/2 RT ( ln|M_{N_f}| - ln|M_3N| ), i.e. Eq. 3 with the cyclic |M_{N_f}|.
    return 0.5 * RT_ * (lnDetM - lnDetZ - lnDetMCartesian_);
}

// ---------------------------------------------------------------------------
//  Absolute orientation of a free-root body, built from a body-fixed atom
//  triplet in the CURRENT Ground frame.
//
//  WHY NOT X_GB[b].R(): the per-block frame handoff (setAtomsLocationsInGround
//  -> recomputeGeometry) rebuilds every root frame with IDENTITY rotation and
//  resets the joint coordinate q = 0 (computeAllAtomFrames r_self branch). So
//  X_GB[b].R() carries only the WITHIN-block deviation from that reset -- it is
//  exactly identity at the START of every move (the gimbal pole, sin(pitch)->0).
//  Reading it floors J(q) to its eps cap for every free root at H_old, every
//  round, independent of the real pose -- the +41352 vs +7352 artifact. The
//  external-rotation Jacobian needs the molecule's TRUE orientation in space
//  (Section 6: "body b's orientation"), which is frame-reset-invariant and read
//  here straight from the atom geometry.
//
//  Convention: x along (root -> reference atom), z along the triplet normal,
//  y = z x x. Any fixed convention is admissible: it shifts J by a per-body
//  constant that is IDENTICAL at H_old and H_new, so only the physical change in
//  orientation across the trajectory survives in dH. Returns false for a body
//  without 3 non-collinear real atoms (no 3D orientation DOF -> no J term).
static bool freeRootAbsRotation(const RobotModel& m, const RobotState& s, int b, Rotation& Rout) {
    const Vec3* P = s.atomPosG();
    const int a0 = m.bodyRootAtom[b];
    if (a0 < 0) {
        return false;
    }
    const Vec3 p0 = P[a0];
    // Pick the two real (mass>0) body atoms, distinct from the root, that span
    // the largest triangle -- the most numerically stable, orientation-defining
    // pair (a near-collinear pick would make the frame ill-conditioned).
    int bestA1 = -1, bestA2 = -1;
    Real bestArea = 0;
    for (int ci = m.bodyAtomsBeg[b]; ci < m.bodyAtomsEnd[b]; ++ci) {
        const int a1 = m.bodyAtoms[ci];
        if (a1 == a0 || m.atomMass[a1] <= Real(0)) {
            continue;
        }
        for (int cj = ci + 1; cj < m.bodyAtomsEnd[b]; ++cj) {
            const int a2 = m.bodyAtoms[cj];
            if (a2 == a0 || m.atomMass[a2] <= Real(0)) {
                continue;
            }
            const Real area = ((P[a1] - p0) % (P[a2] - p0)).norm();
            if (area > bestArea) {
                bestArea = area;
                bestA1 = a1;
                bestA2 = a2;
            }
        }
    }
    if (bestA1 < 0 || bestArea < Real(1e-10)) {
        return false;
    }
    const Vec3 ex = (P[bestA1] - p0) / (P[bestA1] - p0).norm();
    Vec3 ez = ex % (P[bestA2] - p0);
    ez = ez / ez.norm();
    const Vec3 ey = ez % ex;
    Rout = Rotation(robo::Mat33(ex[0], ey[0], ez[0], ex[1], ey[1], ez[1], ex[2], ey[2], ez[2]));
    return true;
}

double World::calcLogSineSqrGamma2() const {
    // Sum over EVERY free root body (parent == Ground, Free joint). gamma2_b is
    // the pitch of the body's TRUE orientation in space (freeRootAbsRotation),
    // NOT X_GB[b].R() -- which the per-block reset pins to identity (q = 0) at
    // the start of every move, flooring J for all roots (see helper).
    double acc = 0.0;
    for (int b = 1; b < model_.numBodies; ++b) {
        if (model_.bodyParent[b] != 0 || model_.bodyJoint[b] != JointType::Free) {
            continue;
        }
        Rotation Rabs;
        if (!freeRootAbsRotation(model_, state_, b, Rabs)) {
            continue; // < 3 non-collinear real atoms: no 3D orientation DOF
        }
        Real w, x, y, z;
        rotationToQuaternion(Rabs, w, x, y, z);
        const Real sinPitch = std::clamp(Real(2.0) * ((w * y) - (z * x)), Real(-1.0), Real(1.0));
        acc += safeLogSineSqr(std::asin(sinPitch));
    }
    return acc;
}

// ---------------------------------------------------------------------------
//  TEMPORARY DIAGNOSTIC -- remove once the J / logSineSqr term is confirmed.
//
//  Prints the per-free-root pitch sine that feeds the external-rotation
//  Jacobian J, now read from the body's TRUE orientation (freeRootAbsRotation),
//  the same quantity calcLogSineSqrGamma2 uses. Called once in reinitialize()
//  (START of move) and once in currentTotalEnergy() (END). After the fix BOTH
//  should show ordinary, non-floored values and nearly the same sum (the pose
//  barely moves in one block) -- i.e. floored=0 at "reinit", and J no longer
//  dominates dH.
static void dbgPrintFreeRootPitches(const char* where, const RobotModel& m, const RobotState& s) {
    int nFree = 0, nFloored = 0, shown = 0;
    double totLss = 0.0;
    std::fprintf(stderr, "[Jdbg %-6s] per-free-root sinPitch (first few shown):\n", where);
    for (int b = 1; b < m.numBodies; ++b) {
        if (m.bodyParent[b] != 0 || m.bodyJoint[b] != JointType::Free) {
            continue;
        }
        ++nFree;
        Rotation Rabs;
        if (!freeRootAbsRotation(m, s, b, Rabs)) {
            continue;
        }
        Real w, x, y, z;
        rotationToQuaternion(Rabs, w, x, y, z);
        const Real sinPitch = std::clamp(Real(2.0) * ((w * y) - (z * x)), Real(-1.0), Real(1.0));
        const Real s2 = sinPitch * sinPitch;
        const Real lss = safeLogSineSqr(std::asin(sinPitch));
        totLss += lss;
        if (s2 < Real(1e-12)) {
            ++nFloored;
        }
        if (shown < 5) {
            std::fprintf(stderr,
                         "[Jdbg %-6s]   body=%-5d quat=(%+.4f %+.4f %+.4f %+.4f) "
                         "sinPitch=%+.3e sin^2=%.3e ln(sin^2)=%+8.3f\n",
                         where,
                         b,
                         w,
                         x,
                         y,
                         z,
                         sinPitch,
                         s2,
                         lss);
            ++shown;
        }
    }
    std::fprintf(stderr,
                 "[Jdbg %-6s] freeRoots=%d  floored(sin^2<1e-12)=%d  sum ln(sin^2)=%.1f\n",
                 where,
                 nFree,
                 nFloored,
                 totLss);
}

// ----------------------------------------------------------------------------
//  Docking: auto-sized binding sphere + conditional reposition ("bound HMC")
//
//  The sphere is sized automatically, PER LIGAND, as
//      R_i = sphereFactor * ( R_receptor + R_ligand_i )
//  R_receptor = max distance from the receptor geometric centre to any receptor
//  atom; R_ligand_i = max distance from ligand i's geometric centre to any of its
//  atoms. The sphere radius is R_receptor + sphereFactor * R_ligand_i: the guest
//  may roam roughly sphereFactor ligand-radii beyond the receptor surface;
//  sphereFactor (from Python) scales only that ligand allowance.
//  Nothing is hard-coded; Robosample sizes it from the current geometry.
//
//  Move: every docking sample is one Generalized-Coordinate HMC move over the
//  guest's external DOF. The kick (reposition) is part of that move's PROPOSAL,
//  not a separate accept/reject: Hold is referenced to the pre-kick state, so a
//  placement that drives the ligand into the receptor explodes the potential and
//  the whole move is rejected on the single acceptance criterion (MH, or always
//  under AlwaysAccept), restoring the pre-kick pose. This is a heuristic move
//  (the teleport is not strictly reversible); for rigorous binding free energies
//  use a flat-bottom / funnel COM restraint added to U with reweighting.
// ----------------------------------------------------------------------------
robo::Vec3 World::atomSetCentroid(const std::vector<int>& atoms) const {
    Vec3 c(0);
    if (atoms.empty()) {
        return c;
    }
    const Vec3* P = state_.atomPosG();
    for (int a : atoms) {
        c += P[a];
    }
    return c / static_cast<Real>(atoms.size());
}

robo::Vec3 World::atomSetMassCenter(const std::vector<int>& atoms) const {
    const Vec3* P = state_.atomPosG();
    Vec3 c(0);
    Real m = 0;
    for (int a : atoms) {
        const Real ma = model_.atomMass[a];
        c += ma * P[a];
        m += ma;
    }
    return (m > 0) ? Vec3(c / m) : c;
}

double World::atomSetRadius(const std::vector<int>& atoms, const robo::Vec3& center) const {
    const Vec3* P = state_.atomPosG();
    Real rmax = 0;
    for (int a : atoms) {
        const Real d = (P[a] - center).norm();
        if (d > rmax) {
            rmax = d;
        }
    }
    return rmax;
}

double World::groupSphereRadius(int g) const {
    const double Rrec = atomSetRadius(siteAtoms_, atomSetCentroid(siteAtoms_));
    const double Rlig = atomSetRadius(ligandGroups_[g], atomSetCentroid(ligandGroups_[g]));
    return Rrec + sampler_.sphereFactor * Rlig;
}

robo::Vec3 World::sampleUniformInSphere(double radius) {
    const Real theta = 2.0 * robo::Pi * uniform_(rng_);
    const Real phi = std::acos(2.0 * uniform_(rng_) - 1.0);
    const Real r = radius * std::cbrt(uniform_(rng_)); // cube-root -> uniform in volume
    return Vec3(r * std::cos(theta) * std::sin(phi), r * std::sin(theta) * std::sin(phi), r * std::cos(phi));
}

robo::Rotation World::sampleUniformRotation() {
    // Shoemake: a uniform unit quaternion -> rotation matrix. Uniform on SO(3).
    const Real u1 = uniform_(rng_), u2 = uniform_(rng_), u3 = uniform_(rng_);
    const Real s1 = std::sqrt(1.0 - u1), s2 = std::sqrt(u1);
    const Real qx = s1 * std::sin(2.0 * robo::Pi * u2);
    const Real qy = s1 * std::cos(2.0 * robo::Pi * u2);
    const Real qz = s2 * std::sin(2.0 * robo::Pi * u3);
    const Real qw = s2 * std::cos(2.0 * robo::Pi * u3);
    return quatToRotation(qw, qx, qy, qz);
}

bool World::repositionLigands(bool forceAll) {
    const Vec3 site = atomSetCentroid(siteAtoms_);
    Vec3* P = state_.atomPosG();
    bool movedAny = false;

    // Receptor geometry (defines the sphere CENTRE and its radius contribution).
    // All lengths in nm (the engine's native unit; AMBER Å inputs are divided by
    // 10 at load time). Printed unconditionally to stderr so every kick is
    // traceable without rebuilding or setting an env var.
    const Vec3 siteMassCom = atomSetMassCenter(siteAtoms_);
    const double Rrec = atomSetRadius(siteAtoms_, site);
    std::fprintf(stderr,
                 "[dock] --- kick decision  always_kick=%s ---\n"
                 "[dock] receptor : nAtoms=%d  "
                 "centroid/sphere-center=(% .4f % .4f % .4f) nm  "
                 "massCOM=(% .4f % .4f % .4f) nm  "
                 "R_receptor=%.4f nm\n",
                 forceAll ? "always" : "containment",
                 (int)siteAtoms_.size(),
                 site[0],
                 site[1],
                 site[2],
                 siteMassCom[0],
                 siteMassCom[1],
                 siteMassCom[2],
                 Rrec);

    // The kick is a PURE proposal: it relocates the ligand and does nothing else.
    // There is NO acceptance here -- the whole move (kick + dynamics) is judged in
    // generateSample, where Hold is referenced to the pre-kick state and an
    // overlapping placement is rejected on the energy validity check.
    // CONTAINMENT (forceAll == false): fire only when the COM has left the sphere.
    // alwaysKick (forceAll == true): fire every round regardless.
    for (int g = 0; g < (int)ligandGroups_.size(); ++g) {
        const auto& grp = ligandGroups_[g];
        if (grp.empty()) {
            continue;
        }
        const Vec3 ligCentroid = atomSetCentroid(grp);
        const double Rlig = atomSetRadius(grp, ligCentroid);
        const double radius = groupSphereRadius(g); // Rrec + sphereFactor * Rlig
        const Vec3 com = atomSetMassCenter(grp);
        const double dist = (com - site).norm();
        const bool inside = dist <= radius;
        const bool willKick = forceAll || !inside;

        std::fprintf(stderr,
                     "[dock] ligand[%d]: nAtoms=%d  "
                     "R_ligand=%.4f nm  "
                     "COM=(% .4f % .4f % .4f) nm  "
                     "|COM-center|=%.4f nm\n"
                     "[dock]   sphere: center=(% .4f % .4f % .4f) nm  "
                     "radius=%.4f nm  "
                     "(R_rec=%.4f + sphereFactor=%.3f * R_lig=%.4f)\n"
                     "[dock]   inside=%s  should_kick=%s  reason=%s\n",
                     g,
                     (int)grp.size(),
                     Rlig,
                     com[0],
                     com[1],
                     com[2],
                     dist,
                     site[0],
                     site[1],
                     site[2],
                     radius,
                     Rrec,
                     sampler_.sphereFactor,
                     Rlig,
                     inside ? "Y" : "N",
                     willKick ? "Y" : "N",
                     forceAll ? "always_kick/rescue"
                              : (inside ? "containment:inside->skip" : "containment:escaped->kick"));

        if (!willKick) {
            continue;
        }

        // Perturb the full external q: uniform COM position in the sphere + uniform
        // reorientation (rigid: rotate the ligand about its COM, then translate).
        const Vec3 offset = sampleUniformInSphere(radius);
        const Vec3 target = site + offset;
        const Rotation R = sampleUniformRotation();
        for (int a : grp) {
            P[a] = target + (R * (P[a] - com));
        }
        movedAny = true;

        // Re-measure COM from the *updated* P[] to confirm the draw actually moved
        // the ligand (sanity: newCOM should equal target within floating-point noise).
        const Vec3 newCom = atomSetMassCenter(grp);
        std::fprintf(stderr,
                     "[dock]   KICK: offset=(% .4f % .4f % .4f) nm  "
                     "|offset|=%.4f nm  (max=radius=%.4f nm)\n"
                     "[dock]   new target=(% .4f % .4f % .4f) nm  "
                     "new COM=(% .4f % .4f % .4f) nm  [Ground/Cartesian, nm]\n",
                     offset[0],
                     offset[1],
                     offset[2],
                     offset.norm(),
                     radius,
                     target[0],
                     target[1],
                     target[2],
                     newCom[0],
                     newCom[1],
                     newCom[2]);
    }

    if (movedAny) {
        // Rebuild q/frames from the proposed Cartesian pose so the dynamics step
        // (and the pre-kick-referenced Hold) are consistent.
        std::vector<robo::Vec3> pos(P, P + model_.numAtoms);
        setAtomsLocationsInGround(pos);
    }
    return movedAny;
}

// ----------------------------------------------------------------------------
//  findGoodStartingPose
//
//  Called once before round 0 when sampler_.maxInitialKickTries > 0.
//  Keeps drawing random placements for every ligand until the immediate
//  post-kick PE change is below the clash ceiling (maxStartPE) for ALL
//  ligands simultaneously -- the same gate the normal per-round pre-step
//  screen uses.  When a clean pose is found it is committed to
//  state_.atomPosG() (and the replica coord array upstream via the
//  normal setAtomsLocationsInGround path), so round 0 starts from a
//  clash-free geometry rather than the raw input file position.
//
//  Returns the number of attempts used.  Throws std::runtime_error if the
//  budget is exhausted without finding a clean pose, so the user gets an
//  immediate, explicit failure rather than a run that silently wastes every
//  round on rejections.
// ----------------------------------------------------------------------------
int World::findGoodStartingPose() {
    if (!docking_ || sampler_.maxInitialKickTries <= 0) {
        return 0;
    }

    // Evaluate the current PE so we have a pePre baseline for dPE gating.
    // (We don't call the full reinitialize() here -- we just need the energy
    // to judge whether a candidate placement is clash-free.)
    RobotEngine::realizePosition(model_, state_);
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);
    bridge_.evaluate(state_);
    const double pePre = bridge_.calcPotentialEnergy();

    std::fprintf(stderr,
                 "[dock] findGoodStartingPose: pePre=%.2f kJ/mol  "
                 "clash_ceiling(maxStartPE)=%.0f  budget=%d tries\n",
                 pePre,
                 sampler_.maxStartPE,
                 sampler_.maxInitialKickTries);

    // Save the current positions so we can restore them if needed.
    std::vector<Vec3> saved(state_.atomPosG(), state_.atomPosG() + model_.numAtoms);

    for (int attempt = 1; attempt <= sampler_.maxInitialKickTries; ++attempt) {
        // Force-kick every ligand unconditionally (forceAll=true).
        repositionLigands(/*forceAll=*/true);

        // Evaluate the post-kick PE at the proposed Cartesian positions.
        // repositionLigands already called setAtomsLocationsInGround, which
        // rebuilt q/frames, so bridge_.evaluate sees the new geometry.
        bridge_.evaluate(state_);
        const double pePost = bridge_.calcPotentialEnergy();
        const double dPE = pePost - pePre;
        const bool clean = std::isfinite(pePost) && (dPE <= sampler_.maxStartPE);

        std::fprintf(stderr,
                     "[dock]   attempt %d/%d: pePost=%.2f  dPE=%+.2f kJ/mol  -> %s\n",
                     attempt,
                     sampler_.maxInitialKickTries,
                     pePost,
                     dPE,
                     clean ? "GOOD (accepted as start)" : "clash, retry");

        if (clean) {
            // Commit: leave state_.atomPosG() at this position.
            // Zero u so the first reinitialize() seeds from rest.
            std::fill(state_.u(), state_.u() + model_.nu, Real(0));
            std::fprintf(stderr,
                         "[dock] findGoodStartingPose: found clean pose in %d attempt(s). "
                         "pePost=%.2f kJ/mol  dPE=%+.2f kJ/mol\n",
                         attempt,
                         pePost,
                         dPE);
            return attempt;
        }

        // Restore the receptor+ligand positions before the next draw so that
        // the sphere-center geometry (centroid of siteAtoms_) is always correct.
        // Only the LIGAND atoms need to be reset -- receptor is welded and
        // setAtomsLocationsInGround doesn't move it -- but the simplest safe
        // approach is to restore everything and let repositionLigands pick a
        // fresh draw next iteration.
        setAtomsLocationsInGround(saved);
    }

    // Budget exhausted.
    char msg[256];
    std::snprintf(msg,
                  sizeof(msg),
                  "findGoodStartingPose: could not find a clash-free starting pose for "
                  "the ligand(s) after %d attempt(s) (maxStartPE=%.0f kJ/mol). "
                  "Check that ligand_molecule_indices is correct and that the receptor "
                  "is minimized. Increase max_initial_kick_tries if the binding site is "
                  "very occluded.",
                  sampler_.maxInitialKickTries,
                  sampler_.maxStartPE);
    throw std::runtime_error(msg);
}

// ----------------------------------------------------------------------------
//  Sampling
// ----------------------------------------------------------------------------
bool World::generateSample() {
    if (sampler_.moveType == MoveType::NcmcSwitch) {
        return ncmcMove();
    }
    if (cartesian_) {
        lastKickApplied_ = false;
        savedPosG_.assign(state_.atomPosG(), state_.atomPosG() + model_.numAtoms);
        bridge_.setAtomPositionsInGround(state_);
        const double peOld = bridge_.calcPotentialEnergy();
        bridge_.setVelocitiesToTemperature(temperature_, static_cast<int>(rng_()));
        bridge_.integrateTrajectoryOnDevice(state_, sampler_.mdSteps, sampler_.timeStep);
        const double peNew = bridge_.calcPotentialEnergy();
        state_.energy.pe = peNew;
        state_.energy.ke = 0.0; // device kinetic energy not pulled back here
        state_.energy.fixman = 0.0;
        state_.energy.logSineSqrGamma2 = 0.0;
        state_.energy.total = peNew;
        if (metropolis(peOld, peNew)) {
            lastAccepted_ = true;
            return true;
        }
        std::copy(savedPosG_.begin(), savedPosG_.end(), state_.atomPosG());
        state_.energy.pe = peOld;
        state_.energy.total = peOld;
        lastAccepted_ = false;
        return false;
    }

    // Internal-coordinate Generalized-Coordinate HMC over the (here: external)
    // DOF. The kick (above) is part of THIS move's proposal, not a separate
    // accept/reject: we save the pre-kick q + energy and reference Hold to the
    // pre-kick potential, so an overlapping placement is penalised by dH AND
    // caught by the energy validity check below -- and the whole move is rejected
    // back to the clean pre-kick state.
    savedQ_.assign(state_.q(), state_.q() + model_.nq); // pre-move (pre-kick) q
    lastKickApplied_ = false;
    double dockPotPre = 0.0; // pre-kick potential part of H
    double pePre = 0.0, fixPre = 0.0, lssPre = 0.0;
    if (docking_) {
        // Save the pre-kick CARTESIAN pose. The kick calls setAtomsLocationsInGround,
        // which redefines the body reference frames from the (kicked) positions and
        // zeroes q -- so restoring saved q would reconstruct the KICKED pose, not
        // this one. The pose, in Ground coordinates, is the reliable thing to keep.
        savedPosG_.assign(state_.atomPosG(), state_.atomPosG() + model_.numAtoms);

        RobotEngine::realizePosition(model_, state_);
        RobotEngine::realizeArticulatedBodyInertias(model_, state_);
        bridge_.evaluate(state_);
        pePre = bridge_.calcPotentialEnergy();
        if (sampler_.useFixman) {
            fixPre = calcFixman();
            lssPre = calcLogSineSqrGamma2();
        }
        dockPotPre = pePre + fixPre - (0.5 * RT_ * lssPre);

        // Force a kick whenever the carried-forward state is unusable as a starting
        // pose. Three triggers, OR'd together:
        //   (a) non-finite pePre                -- nothing can integrate from NaN/Inf;
        //   (b) |pePre| > sampler_.maxStartPE   -- a clash; tied to the SAME absolute
        //       ceiling the acceptance gate uses below, so no admitted pose escapes
        //       the rescue (closes the old dead band where a +9800 clash sat under a
        //       1e4 rescue ceiling yet passed the relative validity gate);
        //   (c) dockingStuckCount_ >= maxStuckRounds -- a finite, sub-ceiling pose
        //       that is nonetheless LOCALLY NON-INTEGRABLE (every reseeded trajectory
        //       diverges). Neither (a) nor (b) catches it, so without this counter the
        //       move loops forever ("PE frozen, q restored"). This makes the trap
        //       escapable independent of energy magnitude.
        const bool stuck = (dockingStuckCount_ >= sampler_.maxStuckRounds);
        // A clash is a HIGH POSITIVE potential (steric overlap, r^-12). A bound
        // pose is strongly NEGATIVE (e.g. -2400 kJ/mol) and is exactly what we
        // want to keep -- it must NOT be flagged "bad". The old test
        // |pePre| > 1e3 force-kicked every well-bound pose, scrambling it every
        // round ("From -2400 -> wrong conf"). Gate on the positive ceiling only,
        // and use the configured maxStartPE (not a hard-coded 1e3) so a strained
        // but immovable receptor offset does not by itself trip the rescue.
        const bool preIsBad = !std::isfinite(pePre) || (pePre > sampler_.maxStartPE);
        std::fprintf(stderr,
                     "[dock] pre-kick: PE_old=%.2f  Fix_old=%.2f  H_old=%.2f kJ/mol  "
                     "triggers: always_kick=%s  preIsBad=%s(maxStartPE=%.0f)  "
                     "stuck=%s(%d/%d)  -> will_kick=%s\n",
                     pePre,
                     fixPre,
                     dockPotPre,
                     sampler_.alwaysKick ? "Y" : "N",
                     preIsBad ? "Y" : "N",
                     sampler_.maxStartPE,
                     stuck ? "Y" : "N",
                     dockingStuckCount_,
                     sampler_.maxStuckRounds,
                     (sampler_.alwaysKick || preIsBad || stuck) ? "Y" : "N");
        lastKickApplied_ = repositionLigands(sampler_.alwaysKick || preIsBad || stuck);
        if (stuck) {
            dockingStuckCount_ = 0; // fresh window after the forced shake
        }
    }

    reinitialize(); // seed u (KE), record Hold_ at the (post-kick) config
    if (docking_) {
        // Re-reference Hold to the PRE-kick potential. KE is config-independent
        // (= 1/2 RT |g|^2), so the only kick-dependent term is the potential.
        Hold_ = dockPotPre + state_.energy.ke;
        state_.energy.total = Hold_;
    }

    // Real pre-trajectory energy components, for an HONEST acceptance log and a
    // correct clash gate. reinitialize() already drew the momenta and folded the
    // resulting kinetic energy into Hold_ -- so KE_old is NOT zero; the metropolis
    // test compares the full Hold_ (PE + KE + Fixman + J) against Hnew. For a
    // torsional world these are the freshly-seeded values; for docking PE/Fix are
    // referenced to the pre-kick pose (KE is configuration-independent).
    const double peOld = docking_ ? pePre : state_.energy.pe;
    const double keOld = state_.energy.ke;
    const double fixOld = docking_ ? fixPre : state_.energy.fixman;

    const Real h = sampler_.timeStep;

    // PRE-STEP CLASH SCREEN (docking). reinitialize() has just evaluated the
    // post-kick forces/PE. If the kick drove the guest into a hard overlap, the
    // *change* in potential is enormous (or already non-finite). Handing such a
    // pose to the Verlet integrator is fatal: a single step against an ~Inf LJ
    // force turns finite q into NaN -- the "Particle coordinate is NaN" crash
    // seen when the proposal overlaps the receptor. So reject BEFORE stepping.
    // Gate on the move-INDUCED change (pePost - pePre), not the absolute total:
    // the rigid receptor may carry a large constant internal energy (e.g. an
    // unminimized protein at ~1e7 kJ/mol) that the ligand world neither created
    // nor can remove, and which cancels in the difference.
    bool stepsOk = true;
    if (docking_) {
        const double pePost = state_.energy.pe; // set by reinitialize()
        const double dPEpost = pePost - pePre;
        const bool screenedOut = !std::isfinite(pePost) || (dPEpost > sampler_.maxStartPE);
        std::fprintf(stderr,
                     "[dock] post-kick proposal (pre-MD):\n"
                     "[dock]   PE_old  = %12.2f kJ/mol\n"
                     "[dock]   PE_new  = %12.2f kJ/mol\n"
                     "[dock]   dPE     = %+12.2f kJ/mol  (new-old, clash ceiling=%.0f)\n"
                     "[dock]   pre-step screen: %s\n",
                     pePre,
                     pePost,
                     dPEpost,
                     sampler_.maxStartPE,
                     screenedOut ? "REJECT (skip MD, dPE>ceiling or non-finite)" : "pass -> run MD");
        if (screenedOut) {
            stepsOk = false;
        }
    }

    // Periodic reversibility probe (THEORY 5.7). Non-destructive: it integrates
    // mdSteps forward + back at this world's dt from the freshly seeded state and
    // restores the state, so the real proposal below is unaffected. It is a
    // smoke test for the CURRENT geometry only -- the always-on guard is the
    // per-step corrector throw inside verletStep. Disabled when interval == 0.
    const long revRound = generateSampleCalls_++;
    if (sampler_.reversibilityCheckInterval > 0 && stepsOk && sampler_.mdSteps > 0
        && (revRound % sampler_.reversibilityCheckInterval == 0)) {
        const Real revResid =
            RobotEngine::checkReversibility(model_, state_, bridge_, constraints_, sampler_.mdSteps, h);
        const Real revTol = Real(1e-6); // relative round-trip residual; ~1e-12 for a clean step
        const bool revBad = !std::isfinite(revResid) || revResid > revTol;
        std::fprintf(stderr,
                     "[rev] world %d round %ld: round-trip residual = %.3e over %d steps at dt=%.6g ps%s\n",
                     index_,
                     revRound,
                     (double)revResid,
                     sampler_.mdSteps,
                     (double)sampler_.timeStep,
                     revBad ? "  <-- WARNING: integrator not reversible at this dt/geometry; reduce timestep"
                            : "  (ok)");
    }

    for (int i = 0; stepsOk && i < sampler_.mdSteps; ++i) {
        // stepTo returns false if a non-finite force/coordinate appeared mid-step;
        // bail immediately so the broken pose is rejected, never carried forward.
        stepsOk = RobotEngine::stepTo(model_, state_, bridge_, constraints_, state_.time + h);
    }

    bool finite = stepsOk;
    if (finite) {
        const Real* qchk = state_.q();
        for (int i = 0; i < model_.nq; ++i) {
            if (!std::isfinite(qchk[i])) {
                finite = false;
                break;
            }
        }
    }

    bool accepted = false;
    if (finite) {
        const double Hnew = currentTotalEnergy(); // sets state_.energy (pe, ke, ...)
        const double peNew = state_.energy.pe;
        const double keNew = state_.energy.ke;
        const double fixNew = state_.energy.fixman;
        const double dPE = peNew - peOld;
        const double dKE = keNew - keOld;
        const double dFix = fixNew - fixOld;
        const double dH = Hnew - Hold_;
        const bool valid = std::isfinite(Hnew) && std::isfinite(peNew) && (dPE <= sampler_.maxStartPE);
        const bool mhPass = metropolis(Hold_, Hnew);
        const char* tag = docking_ ? "dock" : "hmc";
        std::fprintf(stderr,
                     "[%s] post-MD decision:\n"
                     "[%s]   PE_old  = %12.2f   PE_new  = %12.2f   dPE  = %+12.2f kJ/mol\n"
                     "[%s]   KE_old  = %12.2f   KE_new  = %12.2f   dKE  = %+12.2f kJ/mol\n"
                     "[%s]   Fix_old = %12.2f   Fix_new = %12.2f   dFix = %+12.2f kJ/mol\n"
                     "[%s]   H_old   = %12.2f   H_new   = %12.2f   dH   = %+12.2f kJ/mol\n"
                     "[%s]   valid(dPE<=%.0f)=%s  metropolis=%s  -> %s\n",
                     tag,
                     tag,
                     peOld,
                     peNew,
                     dPE,
                     tag,
                     keOld,
                     keNew,
                     dKE,
                     tag,
                     fixOld,
                     fixNew,
                     dFix,
                     tag,
                     Hold_,
                     Hnew,
                     dH,
                     tag,
                     sampler_.maxStartPE,
                     valid ? "Y" : "N",
                     mhPass ? "Y" : "N",
                     (valid && mhPass) ? "ACCEPT" : "reject");
        if (valid && mhPass) {
            RobotEngine::fillAtomPositionsFromBodies(model_, state_);
            accepted = true;
        }
    } else if (!stepsOk) {
        std::fprintf(stderr, "[hmc] post-MD: non-finite force mid-step or pre-screen -> reject\n");
    } else {
        std::fprintf(stderr, "[hmc] q NON-FINITE -> restoring\n");
    }

    if (accepted) {
        if (docking_) {
            dockingStuckCount_ = 0; // moved successfully -> not stuck
        }
        lastAccepted_ = true;
        return true;
    }

    // Rejected (non-finite, clash, or MH): restore the pre-kick state so nothing
    // broken is carried to the next round or reported in the log.
    if (docking_) {
        // Re-establish frames + q + positions from the saved pre-kick Cartesian
        // pose (same entry point the kick used). Restoring q alone is WRONG here:
        // the kick redefined the reference frames, so old q no longer maps to the
        // old pose. This is the fix for the "stuck in a clash forever" failure.
        setAtomsLocationsInGround(savedPosG_);
        state_.energy.pe = pePre;
        state_.energy.ke = 0.0;
        state_.energy.fixman = fixPre;
        state_.energy.logSineSqrGamma2 = lssPre;
        state_.energy.total = dockPotPre;
        // Count this stuck round. When it crosses maxStuckRounds the next entry
        // forces an unconditional kick (above), so no pose can trap the run.
        ++dockingStuckCount_;
    } else {
        std::copy(savedQ_.begin(), savedQ_.end(), state_.q());
        RobotEngine::realizePosition(model_, state_);
        RobotEngine::fillAtomPositionsFromBodies(model_, state_);
    }
    lastAccepted_ = false;
    return false;
}

void World::reinitialize() {
    RobotEngine::realizePosition(model_, state_);
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);

    // Seed u = sqrt(RT) * sqrt(M^-1) * g, g ~ N(0, I). multiplyBySqrtMInv never
    // forms M or M^-1: it sweeps the per-body DI (D^-1) blocks (Jain O(n)).
    const int nu = model_.nu;
    std::vector<Real> g(nu), seeded(nu);
    for (int i = 0; i < nu; ++i) {
        g[i] = gaussian_(rng_);
    }
    RobotEngine::multiplyBySqrtMInv(model_, state_, g.data(), seeded.data());
    const Real scale = std::sqrt(RT_);
    Real* u = state_.u();
    for (int i = 0; i < nu; ++i) {
        u[i] = scale * seeded[i];
    }

    if (!constraints_.empty()) {
        RobotEngine::realizeVelocity(model_, state_);
        constraints_.enforceVelocityConstraints(model_, state_);
    }

    bridge_.evaluate(state_);
    RobotEngine::realizeVelocity(model_, state_);
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);
    RobotEngine::calcUDot(model_, state_);
    RobotEngine::calcQDot(model_, state_, state_.qdot());
    RobotEngine::calcQDotDot(model_, state_);

    const double pe = bridge_.calcPotentialEnergy();
    const double ke = RobotEngine::calcKineticEnergy(model_, state_);
    double fixman = 0.0;
    double logSineSqr = 0.0;
    if (sampler_.useFixman) {
        fixman = calcFixman();
        logSineSqr = calcLogSineSqrGamma2();
        // dbgPrintFreeRootPitches("reinit", model_, state_); // TEMP: J-term diagnostic
    }
    state_.energy.pe = pe;
    state_.energy.ke = ke;
    state_.energy.fixman = fixman;
    state_.energy.logSineSqrGamma2 = logSineSqr;
    Hold_ = pe + ke + fixman - (0.5 * RT_ * logSineSqr);
    state_.energy.total = Hold_;
}

double World::currentTotalEnergy() {
    bridge_.evaluate(state_); // positions -> OpenMM -> forces (+ PE available)
    const double pe = bridge_.calcPotentialEnergy();
    RobotEngine::realizeVelocity(model_, state_);
    const double ke = RobotEngine::calcKineticEnergy(model_, state_);
    double fixman = 0.0;
    double logSineSqr = 0.0;
    if (sampler_.useFixman) {
        fixman = calcFixman(); // realizes ABI internally; position already current
        logSineSqr = calcLogSineSqrGamma2();
        // dbgPrintFreeRootPitches("postMD", model_, state_); // TEMP: J-term diagnostic
    }
    state_.energy.pe = pe;
    state_.energy.ke = ke;
    state_.energy.fixman = fixman;
    state_.energy.logSineSqrGamma2 = logSineSqr;
    state_.energy.total = pe + ke + fixman - (0.5 * RT_ * logSineSqr);
    return state_.energy.total;
}

bool World::metropolis(double Hold, double Hnew) {
    // During burn-in (equilPhase_) every move is accepted regardless of the
    // world's configured acceptRejectMode. This lets the system relax from the
    // starting geometry (which may be far from equilibrium) without the
    // Metropolis gate blocking large-dH moves. The sampler configuration is
    // otherwise unchanged -- same timestep, same mdSteps, same kick logic --
    // so switching to production is a single flag flip with no other state change.
    if (equilPhase_ || sampler_.acceptRejectMode == AcceptRejectMode::AlwaysAccept) {
        return true;
    }
    const double dH = Hnew - Hold;
    if (dH <= 0) {
        return true;
    }
    return uniform_(rng_) < std::exp(-beta_ * dH);
}

void World::configureNcmc(int atomBegin, int atomEnd, int ncmcSteps, double holdFraction) {
    sampler_.ncmcAtomBegin = atomBegin;
    sampler_.ncmcAtomEnd = atomEnd;
    sampler_.ncmcSteps = ncmcSteps;
    sampler_.ncmcHoldFraction = holdFraction;
    sampler_.moveType = MoveType::NcmcSwitch; // must be set AFTER add_sampler
    bridge_.enableAlchemy(atomBegin, atomEnd);
}

double World::protocolLambda(int s) const {
    // s in [0, ncmcSteps): triangle 1 -> 0 -> 1 with an optional flat lambda=0 hold.
    const int N = sampler_.ncmcSteps;
    int hold = static_cast<int>(sampler_.ncmcHoldFraction * N);
    int ramp = (N - hold) / 2;
    if (ramp < 1) {
        ramp = 1;
    }
    if (s < ramp) {
        return 1.0 - static_cast<double>(s) / ramp; // 1 -> 0
    }
    if (s < ramp + hold) {
        return 0.0; // the uncaged stride
    }
    const int up = s - ramp - hold;
    const double l = static_cast<double>(up) / ramp; // 0 -> 1
    return (l > 1.0) ? 1.0 : l;
}

bool World::ncmcMove() {
    // Per-molecule NCMC (Nilmeier, Crooks, Minh & Chodera 2011): a lambda:1->0->1
    // alchemical decouple-move-recouple proposal. During the switch the
    // INTERMOLECULAR nonbonded between this world's molecule and every other
    // molecule is softened (its intramolecular physics is untouched), so the
    // molecule strides free of the cage of contacting molecules near lambda=0 and
    // is recoupled as its mobile DOF relax to fit the (frozen) environment. The
    // environment itself relaxes in the Cartesian world of the Gibbs scan (DOF
    // coverage, THEORY 13.5); a co-mobilized local shell would let it relax inside
    // the move too (documented efficiency follow-up).
    savedQ_.assign(state_.q(), state_.q() + model_.nq);
    bridge_.setAlchemicalLambda(1.0);
    reinitialize(); // draws p ~ N(0, RT M(q0)); sets Hold_ = V1 + K + F + J
    const double Hstart = Hold_;

    // No explicit momentum flip. Momenta are resampled from Maxwell-Boltzmann at
    // the start of every block (reinitialize, above), itself a Gibbs move on the
    // velocity marginal. Under that resampling the NCMC momentum-reversal is
    // unnecessary -- Nilmeier et al. 2011 explicitly sanction reinitializing
    // velocities after each NCMC step -- since the carried-over sign is overwritten
    // next block and never read. (Mirrors the existing MdHmc move, which also
    // resamples and accepts on dH without an explicit flip.)

    const robo::Real h = sampler_.timeStep;
    // DIAGNOSTIC ONLY (NOT used in acceptance): protocol work w = sum over the
    // PERTURBATION substeps of dV at fixed q (Nilmeier et al. 2011, Eq. 16). For a
    // deterministic, reversible, volume-preserving propagator (THEORY 5.5) the heat
    // satisfies dS = 0, so w == Hend - Hstart up to the small non-symplectic drift;
    // logging the gap w - (Hend - Hstart) is a free consistency check (a large gap
    // flags a non-converged corrector / dt too large for this configuration).
    double work = 0.0;
    bool ok = true;
    double Vprev = bridge_.calcPotentialEnergy(); // V at lambda=1, current q
    double lamPrev = 1.0;

    for (int s = 0; ok && s < sampler_.ncmcSteps; ++s) {
        // (i) PERTURB: change lambda at FIXED q; accumulate work = V(lam_new) - V(lam_old).
        const double lam = protocolLambda(s);
        if (lam != lamPrev) {
            bridge_.setAlchemicalLambda(lam);
            bridge_.evaluate(state_); // recompute at same q, new lambda
            const double Vnew = bridge_.calcPotentialEnergy();
            if (!std::isfinite(Vnew)) {
                ok = false;
                break;
            }
            work += Vnew - Vprev;
            Vprev = Vnew;
            lamPrev = lam;
        }
        // (ii) PROPAGATE one Verlet step at fixed lambda (deterministic, reversible).
        ok = RobotEngine::stepTo(model_, state_, bridge_, constraints_, state_.time + h);
        if (ok) {
            Vprev = bridge_.calcPotentialEnergy(); // fixed-lambda V drift = heat, not work
        }
    }

    bool finite = ok;
    if (finite) {
        const robo::Real* q = state_.q();
        for (int i = 0; i < model_.nq; ++i) {
            if (!std::isfinite(q[i])) {
                finite = false;
                break;
            }
        }
    }

    bool accepted = false;
    if (finite) {
        bridge_.setAlchemicalLambda(1.0);
        const double Hend = currentTotalEnergy(); // V1(qend) + K + F + J at lambda=1

        // Acceptance for DETERMINISTIC, reversible, volume-preserving propagation
        // (THEORY 5.5 => path-action dS = 0): accept on the FULL Hamiltonian
        // difference at the lambda=1 endpoints, NOT on the work (Nilmeier et al.
        // 2011, Eq. 20; bistable-dimer Eq. 28). For such integrators w == Hend -
        // Hstart, so adding work would DOUBLE-COUNT. The alchemical cost is already
        // inside Hend - Hstart: a slow protocol lets the mobile DOF relax as lambda
        // returns to 1, lowering Hend; a fast one does not. The thermodynamic
        // (alchemical) perturbation has unit coordinate Jacobian (alpha-ratio = 1)
        // and the 1->0->1 protocol is its own reverse (protocol ratio = 1), so no
        // extra factors enter dH. GATE: lambda == 1 throughout => work == 0 and
        // Hend - Hstart is the plain Verlet dH => reduces EXACTLY to the
        // torsional-HMC metropolis test, which is the existing metropolis() call.
        const double dH = Hend - Hstart;
        std::fprintf(stderr,
                     "[ncmc] world %d: Hstart=%.2f Hend=%.2f dH=%+.2f kJ/mol  "
                     "work(diag)=%.2f gap=%+.2f (steps=%d)\n",
                     index_,
                     Hstart,
                     Hend,
                     dH,
                     work,
                     work - dH,
                     sampler_.ncmcSteps);
        if (std::isfinite(dH) && metropolis(Hstart, Hend)) {
            RobotEngine::fillAtomPositionsFromBodies(model_, state_);
            accepted = true;
        }
    } else {
        std::fprintf(stderr, "[ncmc] non-finite during protocol -> reject\n");
    }

    if (!accepted) {
        bridge_.setAlchemicalLambda(1.0);
        std::copy(savedQ_.begin(), savedQ_.end(), state_.q());
        RobotEngine::realizePosition(model_, state_);
        RobotEngine::fillAtomPositionsFromBodies(model_, state_);
    }
    lastAccepted_ = accepted;
    return accepted;
}