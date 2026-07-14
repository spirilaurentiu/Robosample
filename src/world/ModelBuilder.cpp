// ============================================================================
//  ModelBuilder.cpp - World::buildModel / setRootMobility(ies).
//
//  Relocated verbatim from World.cpp (SPLIT-W1, pure code motion): turns the
//  shared SystemTopology + this world's per-bond mobilities + this world's
//  root mobilities into an immutable RobotModel (body forest by union-find +
//  BFS, topological relabel, q/u index tables, quaternion-slot registry,
//  z-matrix/BAT rows, body-atom ranges, frame graph, mass-property init,
//  loop-closure constraints), plus the two in-place rebuild entry points that
//  re-run it on a root-mobility change.
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>

using robo::Real;
using robo::Rotation;
using robo::Transform;
using robo::Vec3;

namespace {

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

bool isFlexible(JointType m) {
    return RobotModel::jointIsFlexible(m); // flexible == not Weld (== not Rigid)
}

} // namespace

// ----------------------------------------------------------------------------
//  buildModel  (unchanged structure; stores ln|M_3N| for the Fixman reference)
// ----------------------------------------------------------------------------
void World::buildModel(const SystemTopology& sys,
                       const Selection& sel,
                       const std::vector<JointType>& rootMobilities) {
    // Retain the inputs so setRootMobility(ies) can rebuild this world in place.
    // sys is Context-owned and outlives the World; sel / rootMobilities become
    // this world's own copies. Guard the self-rebuild aliasing case (a setter
    // passes &rootMobilities_ back in) so the copy stays well defined.
    sys_ = &sys;
    if (&sel != &sel_) {
        sel_ = sel;
    }
    if (&rootMobilities != &rootMobilities_) {
        rootMobilities_ = rootMobilities;
    }

    // These three are the only append-only (push_back, no prior clear/assign)
    // members below; clearing them here is what makes buildModel idempotent and
    // therefore safe to re-run on a root-mobility change.
    model_.bodyChildren.clear();
    model_.bodyAtoms.clear();
    model_.quaternionQStart.clear();

    const int nAtoms = sys.numAtoms;
    model_.numAtoms = nAtoms;

    if (cartesian_) {
        model_.numBodies = 1;
        model_.numZRows = 0;
        model_.bodyZRow.assign(1, -1);
        model_.atomBody.assign(nAtoms, 0);
        state_.allocateCompact(model_);
        return;
    }

    DSU dsu(nAtoms);
    for (int k = 0; k < sys.numBonds; ++k) {
        if (sys.bondsRingClosing[k]) {
            continue;
        }
        const JointType mob = (k < (int)sel.bondMobility.size()) ? sel.bondMobility[k] : JointType::Rigid;
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
        JointType jt; // the per-bond joint type (was previously discarded -> always Torsion)
    };
    std::vector<Edge> jointEdges;
    for (int k = 0; k < sys.numBonds; ++k) {
        if (sys.bondsRingClosing[k]) {
            continue;
        }
        const JointType mob = (k < (int)sel.bondMobility.size()) ? sel.bondMobility[k] : JointType::Rigid;
        if (isFlexible(mob)) {
            const int ai = sys.bondsI[k], aj = sys.bondsJ[k];
            jointEdges.push_back({atomBody[ai], atomBody[aj], ai, aj, mob});
        }
    }

    std::vector<std::vector<int>> adj(B);
    for (int e = 0; e < (int)jointEdges.size(); ++e) {
        adj[jointEdges[e].bi].push_back(e);
        adj[jointEdges[e].bj].push_back(e);
    }

    std::vector<int> rootBodyOfMol;
    std::vector<JointType> rootMobOfBody(B, JointType::Rigid);
    std::vector<int> rootAtomOfBody(B, -1);
    for (int mol = 0; mol < sys.numMolecules; ++mol) {
        const int rootAtom = sys.atomsBegin[mol];
        const int rb = atomBody[rootAtom];
        rootBodyOfMol.push_back(rb);
        const JointType rm = (mol < (int)rootMobilities.size()) ? rootMobilities[mol] : JointType::Rigid;
        // Validate: only the former-RootMobility subset is meaningful as a
        // molecule-root attachment to Ground (Slider/Cylinder/BendStretch/
        // SphericalCoords need a bond axis that does not exist at the Ground
        // hinge). One enum, one rule -- fail loud instead of degrading to Weld.
        if (!RobotModel::jointIsLegalRoot(rm)) {
            throw std::invalid_argument(
                "World::buildModel: molecule " + std::to_string(mol)
                + " has an illegal root JointType (value " + std::to_string(static_cast<int>(rm))
                + "); legal roots are Free, Translation(Cartesian), Weld(Rigid), FreeLine, Ball, Torsion");
        }
        rootMobOfBody[rb] = rm;
        rootAtomOfBody[rb] = rootAtom;
    }

    model_.bodyParent.assign(B, -1);
    model_.bodyLevel.assign(B, 0);
    model_.bodyJoint.assign(B, JointType::Rigid);
    model_.bodyRootAtom.assign(B, -1);
    std::vector<bool> placed(B, false);
    placed[0] = true;
    std::queue<int> bfs;

    // Root mobility IS a JointType now (validated above), so the former
    // RootMobility->JointType switch is gone -- it is used directly. Cartesian
    // is just the alias spelling of Translation; both arrive here as the same
    // value, so the old silent "FreeLine/Ball -> Weld" degradation cannot recur.

    for (int rb : rootBodyOfMol) {
        if (placed[rb]) {
            continue;
        }
        model_.bodyParent[rb] = 0;
        model_.bodyLevel[rb] = 1;
        model_.bodyJoint[rb] = rootMobOfBody[rb];
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
                model_.bodyJoint[v] = ed.jt; // the per-bond JointType (was hardcoded Torsion)
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
        std::vector<JointType> njt(B, JointType::Rigid);
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

    // q/u sizes come from the single source of truth (RobotModel::jointNQ/NU),
    // so every joint type is counted -- the former local lambdas silently
    // returned 0/0 for everything past Free, zeroing out Slider/Cylinder/
    // BendStretch/Ball/SphericalCoords/FreeLine.
    auto nqOf = [](JointType jt) {
        return RobotModel::jointNQ(jt);
    };
    auto nuOf = [](JointType jt) {
        return RobotModel::jointNU(jt);
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
        // Register the 4-wide quaternion slot of EVERY quaternion body so it is
        // renormalized each step. Previously only Free was pushed, so a Ball (or
        // FreeLine) body's orientation drifted off the unit sphere unchecked --
        // the latent bug isQuaternionBody() already advertised. The quaternion is
        // always the FIRST 4 q of the body (Ball: q0..3; FreeLine/Free: q0..3,
        // translation after), so qc is the quaternion start.
        if (RobotModel::jointUsesQuaternion(jt)) {
            model_.quaternionQStart.push_back(qc);
        }
        qc += nq;
        uc += nu;
        usq += nu * nu;
    }
    model_.nq = qc;
    model_.nu = uc;
    model_.nuSq = usq;

    // ---- z-matrix (BAT) rows, one per FLEXIBLE body (RobotModel.hpp doc) ---
    // bodyParent/bodyRootAtom are already in their FINAL topological order
    // (the relabel above), so this is a pure lookup, no further BFS needed.
    model_.zI.clear();
    model_.zJ.clear();
    model_.zK.clear();
    model_.zL.clear();
    model_.bodyZRow.assign(B, -1);
    for (int b = 1; b < B; ++b) {
        if (!isFlexible(model_.bodyJoint[b])) {
            continue; // Rigid/Weld bodies carry no BAT coordinate
        }
        const int p1 = model_.bodyParent[b];
        const int zI = model_.bodyRootAtom[b];
        const int zJ = (p1 > 0) ? model_.bodyRootAtom[p1] : -1;
        const int p2 = (p1 > 0) ? model_.bodyParent[p1] : -1;
        const int zK = (p2 > 0) ? model_.bodyRootAtom[p2] : -1;
        const int p3 = (p2 > 0) ? model_.bodyParent[p2] : -1;
        const int zL = (p3 > 0) ? model_.bodyRootAtom[p3] : -1;
        model_.bodyZRow[b] = static_cast<int>(model_.zI.size());
        model_.zI.push_back(zI);
        model_.zJ.push_back(zJ);
        model_.zK.push_back(zK);
        model_.zL.push_back(zL);
    }
    model_.numZRows = static_cast<int>(model_.zI.size());

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
//  setRootMobility / setRootMobilities  -- per-world root attachment override
// ----------------------------------------------------------------------------
//  Root mobility is a property of THIS world only. Both mutate this world's own
//  rootMobilities_ copy and rebuild the model from the retained sys_/sel_; the
//  shared SystemTopology is never touched, so the same molecule can have a
//  different root in another world (e.g. Free in a solvation shell, Welded in
//  the bulk). buildModel() is idempotent (it clears its append-only tables), so
//  the rebuild is a full, clean replacement of model_ + state_.
//
//  Must be called after buildModel() (which captures sys_) and BEFORE
//  add_sampler / mass scaling, because the rebuild resets per-body sampler
//  state (bodyMassScale, NMA factors) to their defaults.
void World::setRootMobility(int moleculeIndex, JointType mobility) {
    if (sys_ == nullptr) {
        throw std::logic_error("World::setRootMobility called before buildModel");
    }
    if (moleculeIndex < 0 || moleculeIndex >= (int)rootMobilities_.size()) {
        throw std::out_of_range("World::setRootMobility: molecule index " + std::to_string(moleculeIndex)
                                + " out of range (have " + std::to_string(rootMobilities_.size())
                                + " molecules)");
    }
    rootMobilities_[moleculeIndex] = mobility;
    buildModel(*sys_, sel_, rootMobilities_);
}

void World::setRootMobilities(const std::vector<JointType>& rootMobilities) {
    if (sys_ == nullptr) {
        throw std::logic_error("World::setRootMobilities called before buildModel");
    }
    if ((int)rootMobilities.size() != (int)rootMobilities_.size()) {
        throw std::invalid_argument("World::setRootMobilities: expected "
                                    + std::to_string(rootMobilities_.size()) + " entries, got "
                                    + std::to_string(rootMobilities.size()));
    }
    rootMobilities_ = rootMobilities;
    buildModel(*sys_, sel_, rootMobilities_);
}
