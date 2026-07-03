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
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>

#include "NCMCProtocol.hpp" // robo::ncmc::protocolLambda -- ONE schedule, prod + tests
#include "NMA.hpp"
#include "RobotIntegrator.hpp" // templated verletStep/stepTo/checkReversibility

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

// Opt-in trace of the NMA velocity-distortion draw (DistortOption::NMA). Set the
// env var ROBO_NMA_DEBUG=1 to see, per Gibbs block, exactly how the momentum draw
// is re-pointed. Off by default so the per-block reinitialize() path stays quiet.
bool nmaDebugEnabled() {
    return true;


    static const bool on = [] {
        const char* e = std::getenv("ROBO_NMA_DEBUG");
        return e != nullptr && e[0] != '0' && e[0] != '\0';
    }();
    return on;
}

// Opt-in per-substep trace of Construction II's inner GHMC kernel
// (ncmcInnerGhmcStep). Set the env var ROBO_NCMC_DEBUG=1 to print, for every
// fixed-lambda substep, Hbefore/Hafter/dH/converged/accepted -- diagnosable
// evidence for whether the inner kernel is accepting genuine dt-dependent
// dynamics or freezing (dH pinned to a large/NaN, dt-independent value; see
// docs/specs/ncmc-explicit-solvent/). Off by default (up to ncmcSteps lines
// per move would otherwise flood stderr).
bool ncmcDebugEnabled() {
    static const bool on = [] {
        const char* e = std::getenv("ROBO_NCMC_DEBUG");
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

bool isFlexible(JointType m) {
    return RobotModel::jointIsFlexible(m); // flexible == not Weld (== not Rigid)
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
                          int maxInitialKickTries,
                          std::optional<DistortOption> distortOption,
                          double nmaBiasScale) {
    sampler_.timeStep = timeStep;
    sampler_.mdSteps = mdSteps;
    sampler_.acceptRejectMode = mode;
    sampler_.useNuts = useNuts;
    sampler_.sphereFactor = (sphereFactor > 0.0) ? sphereFactor : 1.0;
    sampler_.alwaysKick = alwaysKick;
    sampler_.clashThreshold = (clashThreshold > 0.0) ? clashThreshold : 1.0e4;
    sampler_.maxInitialKickTries = (maxInitialKickTries > 0) ? maxInitialKickTries : 0;
    sampler_.distortOption = distortOption; // nullopt => no velocity distortion
    sampler_.nmaBiasScale = nmaBiasScale;   // NMA Route B bias magnitude (thermal sigmas)

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
    // AUTO default: Fixman ON for non-Cartesian (torsional + docking) worlds, OFF
    // for Cartesian (flat space => constant metric => no correction). The
    // orientational Jacobian (useOrientationJacobian) stays at its struct default
    // (OFF) -- it is the Euler-angle factor and is wrong for quaternion roots; see
    // SamplerConfig and World::calcLogSineSqrGamma2.
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

void World::setReactionReporter(std::vector<int> interestingBodies) {
    if (cartesian_) {
        throw std::runtime_error(
            "World::setReactionReporter: not supported on a Cartesian (OpenMMVelocityVerlet-equivalent) "
            "world -- it has no articulated body indexing for the interesting-body set to select from "
            "(docs/specs/reaction-force-monitoring.md Sec.3 integrator guard). Flag a torsional/"
            "robotic (add_robotic_world/add_torsional_world) world instead.");
    }
    reactionInterestingBodies_ = std::move(interestingBodies);
    reactionReporter_ = true;
}

void World::enableReactionReporter(bool reportFreeBodies, bool includeOpenmm, bool includeReaction) {
    if (cartesian_) {
        throw std::runtime_error(
            "World::enableReactionReporter: not supported on a Cartesian (OpenMMVelocityVerlet-equivalent) "
            "world -- it has no articulated body indexing for the interesting-body set to select from "
            "(docs/specs/reaction-force-monitoring.md Sec.3 integrator guard). Flag a torsional/"
            "robotic (add_robotic_world/add_torsional_world) world instead.");
    }
    // Interesting-body selection (docs/specs/reaction-force-monitoring.md
    // Sec.2.1/Sec.4): every flexed (non-Weld) body PLUS its parent -- the
    // "both sides of the flex point" domain-boundary picture the old writer
    // used (its childMBIx/parentMBIx pair). bodyForceG is well-defined for
    // ANY non-Ground body (it needs only that body's own atoms, not an
    // inboard atom), so the only exclusion is Ground itself (body 0, no
    // atom to serve as the CSV's atom_idx representative -- Terms: "Ground,
    // body 0, is never included"). In particular a body whose OWN inboard
    // joint attaches directly to Ground (e.g. a Free-rooted receptor's root
    // body) IS included here -- excluding it would silently drop that root
    // body's own applied-force reading (e.g. TM1 of a 7-body TM bundle).
    std::set<int> interesting;
    auto reportable = [](int body) { return body != 0; };
    for (int b = 1; b < model_.numBodies; ++b) {
        if (model_.bodyNU[b] == 0) {
            continue; // Weld: not a flexed joint
        }
        if (reportable(b)) {
            interesting.insert(b);
        }
        const int p = model_.bodyParent[b];
        if (reportable(p)) {
            interesting.insert(p);
        }
    }

    // reportFreeBodies == false: drop FREE-FLOATING RIGID bodies -- bodies
    // that are BOTH Ground-rooted (bodyParent[b] == 0) AND CHILDLESS (no
    // other body's parent is b). Such a body is a lone rigid molecule sitting
    // on a Free root (e.g. a lipid given a Free root but no internal DOFs):
    // it entered `interesting` above only because a Free root always carries
    // bodyNU == 6 != 0, tripping the "flexed" gate even though nothing about
    // it is actually flexed. A flexed molecule's OWN root (e.g. a receptor's
    // root TM body) is Ground-rooted too, but is NOT childless -- its
    // internally-jointed bodies branch from it -- so this predicate leaves it
    // untouched.
    if (!reportFreeBodies) {
        std::vector<int> childCount(static_cast<std::size_t>(model_.numBodies), 0);
        for (int b = 1; b < model_.numBodies; ++b) {
            ++childCount[static_cast<std::size_t>(model_.bodyParent[b])];
        }
        for (auto it = interesting.begin(); it != interesting.end();) {
            const int b = *it;
            if (model_.bodyParent[b] == 0 && childCount[static_cast<std::size_t>(b)] == 0) {
                it = interesting.erase(it);
            } else {
                ++it;
            }
        }
    }

    reactionIncludeOpenmm_ = includeOpenmm;
    reactionIncludeReaction_ = includeReaction;
    setReactionReporter(std::vector<int>(interesting.begin(), interesting.end()));
}

void World::captureReactionSnapshot() {
    reactionSamples_.clear();
    if (!reactionReporter_) {
        return;
    }
    if (!reactionIncludeOpenmm_ && !reactionIncludeReaction_) {
        throw std::runtime_error(
            "World::captureReactionSnapshot: both includeOpenmm and includeReaction are false -- "
            "nothing to report (docs/specs/reaction-force-monitoring.md Sec.1.2). Pass at least one "
            "true to World::enableReactionReporter.");
    }

    // 1. Realize position at the CURRENT (accepted) q. Idempotent/cheap
    // (mirrors currentConstraintLogDet()'s own pattern): guarantees X_GB is
    // current for THIS q regardless of which generateSample() branch
    // (accept, or one of the several reject-restore paths) produced it.
    // bodyForceG's reduction (ForceBridge::getForcesFromOpenMM) needs
    // atomPosG/X_GB only, not articulated-body inertias/velocities.
    RobotEngine::realizePosition(model_, state_);

    // 2. Refresh the field at q: bridge_.evaluate writes atomPosG -> OpenMM,
    // reads back per-atom forces, and reduces them to the per-body spatial
    // force `bodyForceG` (net force + torque about each body's origin, in
    // Ground -- docs/specs/reaction-force-monitoring.md Sec.1/3,
    // ForceBridge::getForcesFromOpenMM). Velocity-independent: this alone
    // cannot leak a trajectory velocity into the recorded quantity, and it
    // only overwrites arena scratch (bodyForceG/atomForceG/mobilityForce)
    // that the next round's evaluate() recomputes from scratch anyway -- no
    // save/restore dance is needed for THIS quantity (unlike the OPTIONAL
    // u=0 reaction term below, which does touch u and restores it). ALWAYS
    // run, regardless of includeOpenmm: `bodyForceG` is also the F_ext the
    // reaction term (calcMobilizerReactionForces) needs as a precondition.
    bridge_.evaluate(state_);

    // 3. OPTIONAL: the static (u=0) mobilizer reaction term (docs/specs/
    // reaction-force-monitoring.md Sec.1.2). Off by default
    // (reactionIncludeReaction_ == false): skipped entirely -- no forward-
    // dynamics step, no u touched, exactly the OpenMM-only path this
    // reporter shipped with.
    std::vector<robo::SpatialVec> reactionAtBo;
    if (reactionIncludeReaction_) {
        // Static (u=0) forward dynamics ("Why static", docs/specs/
        // reaction-force-monitoring.md): save u, zero it, re-realize velocity
        // so every u-derived cache (V_GB, gyro, Coriolis/centrifugal bias) is
        // EXACTLY the u=0 value, then solve forward dynamics so A_GB/zPlus
        // are consistent with u=0 and the field already evaluated in step 2.
        RobotEngine::realizeArticulatedBodyInertias(model_, state_);
        const int nu = model_.nu;
        std::vector<Real> uSave(state_.u(), state_.u() + nu);
        std::fill(state_.u(), state_.u() + nu, Real(0));
        RobotEngine::realizeVelocity(model_, state_);
        RobotEngine::calcUDot(model_, state_);

        // The static reaction: one O(n) inward sweep, at-Bo only.
        reactionAtBo.assign(static_cast<std::size_t>(model_.numBodies),
                            robo::SpatialVec(robo::Vec3(0), robo::Vec3(0)));
        RobotEngine::calcMobilizerReactionForces(model_, state_, reactionAtBo.data(), nullptr);

        // Restore u (and re-sync every u-derived cache to it) -- READ-ONLY
        // w.r.t. the sampler: the next round's generateSample() sees exactly
        // the q/u it would have seen had this snapshot never run.
        std::copy(uSave.begin(), uSave.end(), state_.u());
        RobotEngine::realizeVelocity(model_, state_);
        RobotEngine::calcUDot(model_, state_);
    }

    // 4. Buffer rows: for each interesting body, the SUM of whichever term(s)
    // are enabled, plus the VMD/DCD atom index (prmtop order) of its FIRST
    // atom as a representative atom for placement in vmd/arrows.tcl. Falls
    // back to the engine atom index itself if this world was built without a
    // prmtop permutation (sys_->atomsPrmtopIndex empty/unset).
    const bool havePerm =
        sys_ != nullptr && static_cast<int>(sys_->atomsPrmtopIndex.size()) == model_.numAtoms;
    auto vmdAtomOf = [&](int body) -> int {
        const int engineAtom = model_.bodyAtoms[static_cast<std::size_t>(model_.bodyAtomsBeg[body])];
        return havePerm ? sys_->atomsPrmtopIndex[static_cast<std::size_t>(engineAtom)] : engineAtom;
    };
    const robo::SpatialVec* bodyForceG = state_.bodyForceG();
    reactionSamples_.reserve(reactionInterestingBodies_.size());
    for (int b : reactionInterestingBodies_) {
        robo::Vec3 force(0);
        robo::Vec3 torque(0);
        if (reactionIncludeOpenmm_) {
            force += bodyForceG[static_cast<std::size_t>(b)].linear;
            torque += bodyForceG[static_cast<std::size_t>(b)].angular;
        }
        if (reactionIncludeReaction_) {
            force += reactionAtBo[static_cast<std::size_t>(b)].linear;
            torque += reactionAtBo[static_cast<std::size_t>(b)].angular;
        }
        reactionSamples_.push_back(ReactionSample{b, vmdAtomOf(b), force, torque});
    }
}

double World::setNMASoftModeFromHessian(const std::vector<double>& atomPosGFlat, double h, double zeroTol) {
    if (cartesian_) {
        return 0.0; // Route B is an internal-coordinate move; Cartesian world has nu==1 body
    }
    // Position the world at the minimized geometry. This rebuilds the rigid-body
    // frames (recomputeGeometry) and sets q<-0, u<-0, then realizePosition. The
    // minimum is therefore at q0==0 in these frames -- exactly where we want H.
    const int nA = model_.numAtoms;
    std::vector<robo::Vec3> pos(static_cast<std::size_t>(nA));
    for (int a = 0; a < nA; ++a) {
        pos[a] = robo::Vec3(atomPosGFlat[3 * a], atomPosGFlat[3 * a + 1], atomPosGFlat[3 * a + 2]);
    }
    setAtomsLocationsInGround(pos);

    robo::RouteBNMA nma = robo::computeRouteBNMA(model_, state_, h, zeroTol);

    // This is the hand-off: a length-nu vector here makes reinitialize() use it
    // (its all-ones fallback only fires when the size != nu).
    uScaleFactors_ = nma.uScaleFactors;

    return (nma.softMode >= 0) ? nma.eigval[nma.softMode] : 0.0;
}

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
    // External-rotation Jacobian J(q) = -(1/2) RT sum_b ln sin^2(gamma2_b), the
    // sin(theta) polar volume factor of an EULER/SPHERICAL parameterization of
    // orientation. This is gated by sampler_.useOrientationJacobian, DEFAULT OFF,
    // because Robosample's free/ball roots are parameterized by UNIT QUATERNIONS,
    // not Euler angles, and for a quaternion root the term does NOT belong:
    //
    //   The flat (uniform) measure on the unit 3-sphere S^3 -- i.e. drawing the
    //   quaternion uniformly subject to ||q|| = 1 -- pushes forward to EXACTLY
    //   the Haar (rotation-invariant) measure on SO(3). "Flat in quaternion
    //   space" and "Haar-uniform on rotations" are the same distribution; there
    //   is no Jacobian between them. The quaternion exp-map integrator
    //   (RobotIntegrator advanceQuatExp) is the geodesic flow of the
    //   left-invariant kinetic metric, and for a single free body with constant
    //   body inertia det M(q) is orientation-independent, so GC-HMC with U = 0
    //   already samples orientation Haar-uniformly with NO correction term. The
    //   sin^2(gamma2) factor is the correct Jacobian only if one sampled in Euler
    //   angles; applied on top of the quaternion parameterization it biases the
    //   marginal toward the poles. See tests/TestEnsembleOrientation.cpp, which
    //   demonstrates: OFF -> Haar-uniform; ON -> measurably biased.
    //
    // The term is retained behind the flag for any future Euler-parameterized
    // joint where it WOULD be correct. gamma2_b is the pitch of the body's TRUE
    // orientation in space (freeRootAbsRotation), NOT X_GB[b].R() -- which the
    // per-block reset Torsions to identity (q = 0) at the start of every move,
    // flooring J for all roots (see helper).
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
    // overlapTorsiong placement is rejected on the energy validity check.
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
    // pre-kick potential, so an overlapTorsiong placement is penalised by dH AND
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
        }
        if (sampler_.useOrientationJacobian) {
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
    // seen when the proposal overlaps the receptor. So reject BEFORE stepTorsiong.
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

void World::setCartesianSolvent(const std::vector<int>& atomIndices) {
    // Keep only real (massive) atoms: massless virtual sites/EPs have no
    // independent Cartesian DOF (OpenMM reconstructs them from their parents), so
    // they must NOT be Verlet-integrated here. The mask also drives the skip in
    // fillAtomPositionsFromBodies, so it must list exactly the integrated atoms.
    std::vector<int> kept;
    kept.reserve(atomIndices.size());
    for (int a : atomIndices) {
        if (a >= 0 && a < model_.numAtoms && model_.atomMass[a] > robo::Real(0)) {
            kept.push_back(a);
        }
    }

    // PRECONDITION P1 (docs/specs/two-robot-contact/10-mixed-integrator-
    // correctness.md Sec.2): buildModel assigns EVERY atom to some articulated
    // body, so a Cartesian-integrated atom is safe ONLY inside a 0-DOF
    // Weld/Rigid body whose atom membership the mask covers EXACTLY. If a
    // flagged atom instead sat in a body with bodyNU > 0 (e.g. a Free-rooted
    // second robot), that body's KE would be double-counted (the articulated
    // `ke` from calcKineticEnergy PLUS `keSolvent` from calcSolventKE) and its
    // momentum double-drawn (multiplyBySqrtMInv PLUS drawSolventVelocities) --
    // a silent, non-crashing Boltzmann corruption, not a crash. Likewise a body
    // only PARTIALLY covered by the mask breaks the "one rigid body, one
    // motion" invariant the zero cross-Jacobians of M1 depend on. Fail loud
    // instead of sampling the wrong density.
    {
        std::vector<char> flagged(static_cast<std::size_t>(model_.numAtoms), 0);
        for (int a : kept) {
            flagged[static_cast<std::size_t>(a)] = 1;
        }
        std::vector<char> bodyChecked(static_cast<std::size_t>(model_.numBodies), 0);
        for (int a : kept) {
            const int b = model_.atomBody[a];
            if (bodyChecked[static_cast<std::size_t>(b)]) {
                continue;
            }
            bodyChecked[static_cast<std::size_t>(b)] = 1;
            if (model_.bodyNU[b] != 0) {
                throw std::runtime_error(
                    "World::setCartesianSolvent: flagged atom " + std::to_string(a) + " belongs to body "
                    + std::to_string(b) + " with " + std::to_string(model_.bodyNU[b])
                    + " DOF; Cartesian-integrated atoms MUST sit in a 0-DOF Weld/Rigid body "
                      "(PRECONDITION P1, docs/specs/two-robot-contact/10-mixed-integrator-correctness.md "
                      "Sec.2). A DOF>0 body double-counts KE and double-draws momentum for its atoms.");
            }
            for (int ci = model_.bodyAtomsBeg[b]; ci < model_.bodyAtomsEnd[b]; ++ci) {
                const int bodyAtom = model_.bodyAtoms[ci];
                if (model_.atomMass[bodyAtom] <= robo::Real(0)) {
                    continue; // massless virtual site: never flagged, never required
                }
                if (!flagged[static_cast<std::size_t>(bodyAtom)]) {
                    throw std::runtime_error(
                        "World::setCartesianSolvent: body " + std::to_string(b)
                        + " is only PARTIALLY covered by the cartSolvent mask (atom "
                        + std::to_string(bodyAtom)
                        + " of this 0-DOF body is not flagged); the mask must cover a flagged body's atoms "
                          "EXACTLY (PRECONDITION P1, docs/specs/two-robot-contact/"
                          "10-mixed-integrator-correctness.md Sec.2).");
                }
            }
        }
    }

    state_.setCartSolvent(kept, model_.atomMass.data());
    std::fprintf(stderr,
                 "[ncmc] world %d: Cartesian-integrated solvent atoms = %d (of %d requested); the "
                 "contact environment now relaxes inside the proposal.\n",
                 index_,
                 static_cast<int>(kept.size()),
                 static_cast<int>(atomIndices.size()));
}

void World::drawSolventVelocities() {
    // Maxwell-Boltzmann v_s ~ N(0, RT/m_s) per Cartesian component. Mirrors the
    // generalized-momentum draw in reinitialize() (a Gibbs update of the velocity
    // marginal), so the move needs no explicit momentum flip. No-op when empty.
    const std::vector<int>& solv = state_.cartSolventAtoms();
    const std::vector<robo::Real>& invM = state_.cartSolventInvMass();
    robo::Vec3* velG = state_.atomVelG();
    for (std::size_t j = 0; j < solv.size(); ++j) {
        const robo::Real sigma = std::sqrt(RT_ * invM[j]); // sqrt(RT/m)
        robo::Vec3& v = velG[solv[j]];
        v = robo::Vec3(sigma * gaussian_(rng_), sigma * gaussian_(rng_), sigma * gaussian_(rng_));
    }
}

double World::calcSolventKE() const {
    // 1/2 sum_s m_s |v_s|^2 over the Cartesian-integrated atoms (flat metric).
    const std::vector<int>& solv = state_.cartSolventAtoms();
    const std::vector<robo::Real>& invM = state_.cartSolventInvMass();
    const robo::Vec3* velG = state_.atomVelG();
    double ke = 0.0;
    for (std::size_t j = 0; j < solv.size(); ++j) {
        const robo::Real m = (invM[j] > robo::Real(0)) ? (robo::Real(1) / invM[j]) : robo::Real(0);
        const robo::Vec3& v = velG[solv[j]];
        ke += m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    return 0.5 * ke;
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

    // simtk NMA Route B (DistortOption::NMA): draw the momentum from a symmetric
    // two-component Gaussian MIXTURE biased along the NMA direction instead of the
    // isotropic Gaussian. Us = z + s*mu, s = +/-1 uniform, mu = alpha * uhat where
    // uhat = uScaleFactors_/||uScaleFactors_|| is a UNIT direction and alpha =
    // nmaBiasScale is the directed push in thermal sigmas (so ||mu||^2 = alpha^2).
    // The map u = sqrt(RT) * M^-1/2 * Us below is unchanged; only the seed differs.
    // The bias steers proposals along soft directions; detailed balance is restored
    // in the acceptance by ke_mix = ke - RT*ln cosh(w.mu) (see nmaKineticCorrection).
    // With alpha=0 the mixture collapses to the plain draw. nullopt (None) => skip.
    if (sampler_.distortOption == DistortOption::NMA) {
        if (static_cast<int>(uScaleFactors_.size()) != nu) {
            uScaleFactors_.assign(nu, Real{1});
        }
        // Bias mu = alpha * uhat, where uhat = uScaleFactors_ / ||uScaleFactors_|| is a
        // UNIT direction and alpha = nmaBiasScale is the directed push in thermal-sigma
        // units. Hence ||mu||^2 = alpha^2 (NOT nu): the bias injects only ~1/2 RT alpha^2
        // of directed energy, so alpha controls boldness vs acceptance directly.
        Real norm2 = 0;
        for (int i = 0; i < nu; ++i) {
            norm2 += uScaleFactors_[i] * uScaleFactors_[i];
        }
        const Real alpha = sampler_.nmaBiasScale;
        const Real unitScale = (norm2 > 0) ? (alpha / std::sqrt(norm2)) : Real{0};
        nmaBias_.assign(nu, Real{0});
        for (int i = 0; i < nu; ++i) {
            nmaBias_[i] = uScaleFactors_[i] * unitScale;
        }

        // Symmetric mixture sign s = +/-1.
        const Real s = (uniform_(rng_) < 0.5) ? Real{-1} : Real{1};

        const bool trace = nmaDebugEnabled();
        Real zN2 = 0, muDotZ = 0;
        if (trace) {
            for (int i = 0; i < nu; ++i) {
                zN2 += g[i] * g[i];
                muDotZ += nmaBias_[i] * g[i];
            }
        }

        // Shift the white noise by the signed bias: Us = z + s*mu (in place in g).
        Real muDotUs = 0, UsN2 = 0;
        for (int i = 0; i < nu; ++i) {
            g[i] += s * nmaBias_[i];
            muDotUs += nmaBias_[i] * g[i];
            UsN2 += g[i] * g[i];
        }

        if (trace) {
            const int k = std::min(nu, 8);
            std::cout << "[nma] world " << index_ << ": Route B mixture draw, nu=" << nu << ", sign s=" << s
                      << ", alpha=" << alpha << "\n";
            std::cout << "[nma]   mu=alpha*uhat (first " << k << "): ";
            for (int i = 0; i < k; ++i) {
                std::cout << nmaBias_[i] << ' ';
            }
            std::cout << (k < nu ? "...\n" : "\n");
            std::cout << "[nma]   ||mu||^2=" << (alpha * alpha)
                      << " (==alpha^2; directed energy ~1/2 RT alpha^2=" << (0.5 * RT_ * alpha * alpha)
                      << "), ||z||^2=" << zN2 << ", mu.z=" << muDotZ << "\n";
            std::cout << "[nma]   Us=z+s*mu: ||Us||^2=" << UsN2 << ", mu.Us=" << muDotUs
                      << "  (mu.Us is the START w.mu; multiplyBySqrtM must reproduce it)\n"
                      << std::flush;
        }
    }

    RobotEngine::multiplyBySqrtMInv(model_, state_, g.data(), seeded.data());
    const Real scale = std::sqrt(RT_);
    Real* u = state_.u();
    for (int i = 0; i < nu; ++i) {
        u[i] = scale * seeded[i];
    }

    // Draw the Cartesian solvent velocities from the same Maxwell-Boltzmann
    // marginal (independent of the generalized draw -- flat diagonal metric).
    drawSolventVelocities();

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
    const double nmaCorr = nmaKineticCorrection(); // RT*ln cosh(w.mu); 0 unless Route B
    if (sampler_.distortOption == DistortOption::NMA && nmaDebugEnabled()) {
        // Physical KE = 1/2 u^T M u (equipartition target nu/2*RT), and the Route B
        // kinetic ke_mix = ke - nmaCorr that actually enters the acceptance H.
        std::cout << "[nma]   START: ke=1/2 u^T M u=" << ke << " (target nu/2*RT=" << (0.5 * nu * RT_)
                  << "), ke_mix=ke-corr=" << (ke - nmaCorr) << "\n"
                  << std::flush;
    }
    double fixman = 0.0;
    double logSineSqr = 0.0;
    if (sampler_.useFixman) {
        fixman = calcFixman();
    }
    if (sampler_.useOrientationJacobian) {
        logSineSqr = calcLogSineSqrGamma2();
    }
    const double keSolvent = calcSolventKE(); // 0 when no Cartesian solvent
    state_.energy.pe = pe;
    state_.energy.ke = ke;
    state_.energy.keSolvent = keSolvent;
    state_.energy.fixman = fixman;
    state_.energy.logSineSqrGamma2 = logSineSqr;
    // ke_mix = ke - nmaCorr replaces the kinetic term for the NMA Route B mixture
    // draw (nmaCorr == 0 for ordinary HMC, so Hold_ is unchanged off Route B).
    // keSolvent is the flat-space solvent kinetic energy (0 off solvent-relaxing
    // NCMC), so Hold_ stays identical to the welded path when there is no solvent.
    Hold_ = pe + ke + keSolvent + fixman - (0.5 * RT_ * logSineSqr) - nmaCorr;
    state_.energy.total = Hold_;
}

double World::nmaKineticCorrection() {
    // RT*ln cosh(w.mu), w = M^(1/2) u / sqrt(RT), mu = nmaBias_. ke_mix = ke - this
    // equals -RT*ln g(u|q) for the symmetric biased-mixture momentum draw (up to the
    // same q-independent constant the standard kinetic term drops), so subtracting it
    // from both Hold_ and Hnew makes the Route B move detailed-balanced. Returns 0
    // (strict no-op) for ordinary HMC, so it never perturbs the standard path.
    if (sampler_.distortOption != DistortOption::NMA) {
        return 0.0;
    }
    const int nu = model_.nu;
    if (static_cast<int>(nmaBias_.size()) != nu) {
        return 0.0; // bias not yet built (no reinitialize this block) => no correction
    }
    // multiplyBySqrtM needs the articulated-body inertias at the current config; it
    // uses local scratch, so it does NOT disturb V_GB (the velocities calcKineticEnergy
    // just read). Position is already realized at every call site.
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);
    std::vector<Real> sqrtMu(nu);
    RobotEngine::multiplyBySqrtM(model_, state_, state_.u(), sqrtMu.data()); // M^(1/2) u
    const Real inv = Real(1) / std::sqrt(RT_);
    Real wDotMu = 0.0;
    for (int i = 0; i < nu; ++i) {
        wDotMu += (sqrtMu[i] * inv) * nmaBias_[i];
    }
    // Numerically stable ln cosh(x) = |x| + log1p(exp(-2|x|)) - ln 2.
    const Real ax = std::abs(wDotMu);
    const Real lnCosh = ax + std::log1p(std::exp(-2.0 * ax)) - std::log(2.0);
    if (nmaDebugEnabled()) {
        // w.mu via multiplyBySqrtM; at the START this must match the mu.Us printed
        // by reinitialize (a live check that sqrt(M) inverts sqrt(M^-1) on the draw).
        std::cout << "[nma]   w.mu=" << wDotMu << " (via M^1/2), RT*ln cosh(w.mu)=" << (RT_ * lnCosh) << "\n"
                  << std::flush;
    }
    return RT_ * lnCosh;
}

double World::currentTotalEnergy() {
    bridge_.evaluate(state_); // positions -> OpenMM -> forces (+ PE available)
    const double pe = bridge_.calcPotentialEnergy();
    RobotEngine::realizeVelocity(model_, state_);
    const double ke = RobotEngine::calcKineticEnergy(model_, state_);
    const double nmaCorr = nmaKineticCorrection(); // RT*ln cosh(w.mu) at the END; 0 unless Route B
    if (sampler_.distortOption == DistortOption::NMA && nmaDebugEnabled()) {
        std::cout << "[nma]   END:   ke=1/2 u^T M u=" << ke << ", ke_mix=ke-corr=" << (ke - nmaCorr)
                  << " (corr=RT*ln cosh(w.mu)=" << nmaCorr << ")\n"
                  << std::flush;
    }
    double fixman = 0.0;
    double logSineSqr = 0.0;
    if (sampler_.useFixman) {
        fixman = calcFixman(); // realizes ABI internally; position already current
    }
    if (sampler_.useOrientationJacobian) {
        logSineSqr = calcLogSineSqrGamma2();
    }
    const double keSolvent = calcSolventKE(); // 0 when no Cartesian solvent
    state_.energy.pe = pe;
    state_.energy.ke = ke;
    state_.energy.keSolvent = keSolvent;
    state_.energy.fixman = fixman;
    state_.energy.logSineSqrGamma2 = logSineSqr;
    // ke_mix = ke - nmaCorr (Route B mixture); nmaCorr == 0 off Route B. keSolvent
    // is the flat-space solvent KE (0 off solvent-relaxing NCMC).
    state_.energy.total = pe + ke + keSolvent + fixman - (0.5 * RT_ * logSineSqr) - nmaCorr;
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

void World::configureNcmc(std::vector<int> atomIndices, int ncmcSteps, double holdFraction) {
    // Region A is an arbitrary atom-index SET (docs/specs/ncmc-explicit-solvent/
    // 30-region-and-protocol-policy.md Sec.2). Sort + dedupe so ncmcTeleportRoot's
    // "first atom" read and the OpenMM-side index-set iteration are well defined.
    std::sort(atomIndices.begin(), atomIndices.end());
    atomIndices.erase(std::unique(atomIndices.begin(), atomIndices.end()), atomIndices.end());
    sampler_.ncmcAtomIndices = atomIndices;
    sampler_.ncmcSteps = ncmcSteps;
    sampler_.ncmcHoldFraction = holdFraction;
    sampler_.moveType = MoveType::NcmcSwitch; // must be set AFTER add_sampler
    bridge_.enableAlchemy(sampler_.ncmcAtomIndices);
}

void World::configureNcmc(int atomBegin, int atomEnd, int ncmcSteps, double holdFraction) {
    std::vector<int> atomIndices;
    if (atomEnd > atomBegin) {
        atomIndices.resize(static_cast<std::size_t>(atomEnd - atomBegin));
        std::iota(atomIndices.begin(), atomIndices.end(), atomBegin);
    }
    configureNcmc(std::move(atomIndices), ncmcSteps, holdFraction);
}

auto World::protocolLambda(int step) const -> double {
    // Delegate to the pure, unit-tested schedule (include/NCMCProtocol.hpp) so
    // production and the tests exercise ONE palindromic 1 -> 0 -> 1 schedule and
    // can never drift. Palindromic + endpoints pinned to lambda = 1 is what makes
    // the composed NCMC map F-reversible (see ncmcMove's acceptance comment).
    return robo::ncmc::protocolLambda(step, sampler_.ncmcSteps, sampler_.ncmcHoldFraction);
}

int World::ncmcTeleportRoot() const {
    // The tree root (parent == Ground) of the body carrying the FIRST (smallest-
    // index) atom of Region A, iff it is a Free joint (so a rigid reposition is
    // one quaternion+translation). ncmcAtomIndices is kept sorted ascending
    // (configureNcmc), so .front() is that atom for both the contiguous and the
    // general index-set case.
    if (sampler_.ncmcAtomIndices.empty()) {
        return -1;
    }
    const int firstAtom = sampler_.ncmcAtomIndices.front();
    if (firstAtom < 0 || firstAtom >= model_.numAtoms) {
        return -1;
    }
    int b = model_.atomBody[firstAtom];
    while (b > 0 && model_.bodyParent[b] != 0) {
        b = model_.bodyParent[b];
    }
    if (b > 0 && model_.bodyJoint[b] == JointType::Free) {
        return b;
    }
    return -1;
}

void World::ncmcApplyTroughTeleport(int rootBody) {
    // Rigid reposition of the region root at the λ=0 ghost trough. Mirrors the
    // engine-validated tests/TeleportMove.hpp teleportFreeRoot: set the root
    // orientation (Haar) and translation (uniform in the docking sphere), and
    // CO-ROTATE the root's angular AND linear speeds by ΔR = Rnew R_oldᵀ so KE is
    // exactly preserved (M_ang(q) is orientation-dependent). Operates on q/u
    // directly -- NO setAtomsLocationsInGround -- so the live momenta stay valid.
    robo::Real* q = state_.q();
    robo::Real* u = state_.u();
    const int qOff = model_.bodyQIndex[rootBody];
    const int uOff = model_.bodyUIndex[rootBody];

    const Rotation Rold = quatToRotation(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
    const Rotation Rnew = sampleUniformRotation();
    const Rotation dR(Rnew * Rold.transpose());

    const Vec3 site = atomSetCentroid(siteAtoms_);
    const double radius = atomSetRadius(siteAtoms_, site);
    const Vec3 target = site + sampleUniformInSphere(radius);

    Real w, x, y, z;
    rotationToQuaternion(Rnew, w, x, y, z);
    q[qOff + 0] = w;
    q[qOff + 1] = x;
    q[qOff + 2] = y;
    q[qOff + 3] = z;
    q[qOff + 4] = target[0];
    q[qOff + 5] = target[1];
    q[qOff + 6] = target[2];

    const Vec3 wR = dR * Vec3(u[uOff + 0], u[uOff + 1], u[uOff + 2]);
    const Vec3 vR = dR * Vec3(u[uOff + 3], u[uOff + 4], u[uOff + 5]);
    u[uOff + 0] = wR[0];
    u[uOff + 1] = wR[1];
    u[uOff + 2] = wR[2];
    u[uOff + 3] = vR[0];
    u[uOff + 4] = vR[1];
    u[uOff + 5] = vR[2];

    // re-realize geometry and re-seed the derivative chain for the next Verlet step
    RobotEngine::realizePosition(model_, state_);
    RobotEngine::fillAtomPositionsFromBodies(model_, state_);
    bridge_.evaluate(state_);
    RobotEngine::realizeVelocity(model_, state_);
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);
    RobotEngine::calcUDot(model_, state_);
    RobotEngine::calcQDot(model_, state_, state_.qdot());
    RobotEngine::calcQDotDot(model_, state_);
}

bool World::ncmcInnerGhmcStep(robo::Real h, int substepIndex, bool* acceptedOut) {
    // Construction II inner kernel (docs/specs/ncmc-explicit-solvent/
    // 10-acceptance-construction.md Sec.3; 20-inner-integrator.md Sec.2 F3).
    //
    // One fixed-lambda Verlet substep, proposed and then Metropolis accept/
    // reject'd against the FULL H_lambda = V_lambda + K(u;q) + K_s(v_s) + U_F(q)
    // - 1/2 RT ln sin^2(gamma2) - nmaCorr -- reusing currentTotalEnergy()'s own
    // assembly VERBATIM (F1/F2: a shortcut that dropped U_F/pitch here would make
    // the inner kernel non-pi_lambda-invariant and silently sample the WRONG
    // density; INV0 in tests/TestNcmcExplicitSolvent.cpp is the discriminating
    // oracle). On reject: restore (q,u,x_s,v_s) to their pre-substep values and
    // NEGATE (u,v_s) -- the standard GHMC reject-flip that keeps this kernel
    // F-reversible and hence exactly pi_lambda-invariant.
    //
    // F3: a velocity corrector that does not converge at this dt is ALSO an
    // automatic reject (never silently taken as Construction I's propagator
    // does) -- under Construction II a non-converged corrector need not be
    // F-reversible, and taking it as the GHMC proposal would break
    // pi_lambda-invariance with no visible symptom in the outer acceptance
    // (10-...:Sec.5 NOTE F3, 20-...:Sec.2 second bullet).
    if (acceptedOut) {
        *acceptedOut = false; // pessimistic default; set true only on a genuine accept below
    }
    const std::vector<int>& solv = state_.cartSolventAtoms();
    robo::Vec3* posG = state_.atomPosG();
    robo::Vec3* velG = state_.atomVelG();
    std::vector<Real> q0(state_.q(), state_.q() + model_.nq);
    std::vector<Real> u0(state_.u(), state_.u() + model_.nu);
    std::vector<robo::Vec3> xs0(solv.size()), vs0(solv.size());
    for (std::size_t j = 0; j < solv.size(); ++j) {
        xs0[j] = posG[solv[j]];
        vs0[j] = velG[solv[j]];
    }

    const double Hbefore = currentTotalEnergy(); // H_lambda at the CURRENT bridge lambda;
    // also refreshes bodyForceG/mobilityForce (bridge_.evaluate) and V_GB/qdot
    // (realizeVelocity) at (q0, the CURRENT lambda).

    // BUG FIX (state/energy-assembly inconsistency, freeze root cause): the
    // PERTURB substep in ncmcMove (lambda change) refreshes FORCES via
    // bridge_.evaluate but never recomputes udot/qdotdot -- it doesn't need to
    // for Construction I, whose stepTo call seeds that chain once at the top
    // of the move and never revisits it mid-substep. Construction II calls
    // THIS function once per substep, immediately after the SAME perturb
    // block may have just changed lambda, so without a local reseed here
    // `state_.udot()`/`qdotdot()` at verletStep's entry (the `a0` term of its
    // velocity-Verlet trapezoid) are evaluated at the STALE, PRE-perturb
    // lambda's forces -- a genuine assembly inconsistency baked directly into
    // the proposal itself (dt-INDEPENDENT: it does not shrink at small h the
    // way ordinary shadow work does, since it is a wrong-physics `a0`, not a
    // discretization-error `a0`). That corrupts the GHMC proposal and, given
    // alchemical decoupling can change intermolecular forces by orders of
    // magnitude at low lambda, is large enough to make every dH huge
    // regardless of dt -- exactly the observed 100x-dt-invariant freeze.
    // realizeArticulatedBodyInertias is called again UNCONDITIONALLY (cheap,
    // idempotent) rather than relying on calcFixman's internal call above, so
    // this reseed is correct even with useFixman=false. Scoped to
    // ncmcInnerGhmcStep (Construction-II-only code, never called by
    // Construction I), so Construction I's shared perturb block and stepTo
    // call are untouched -- bit-for-bit preserved.
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);
    RobotEngine::calcUDot(model_, state_);
    RobotEngine::calcQDot(model_, state_, state_.qdot());
    RobotEngine::calcQDotDot(model_, state_);

    bool converged = true;
    const bool stepOk =
        RobotEngine::stepTo(model_, state_, bridge_, constraints_, state_.time + h, &converged);
    if (!stepOk) {
        // Non-finite force/velocity: verletStep already restored its own
        // pre-step state internally. Mirror stepTo's contract exactly -- the
        // caller (ncmcMove) treats false as a hard abort of the whole move,
        // orthogonal to the F3 corrector-convergence guard below.
        if (ncmcDebugEnabled()) {
            std::fprintf(stderr,
                         "[ncmc-inner] world %d substep %d: h=%.6g stepTo FAILED (non-finite) -> abort move\n",
                         index_,
                         substepIndex,
                         (double)h);
        }
        return false;
    }

    bool accept = false;
    double Hafter = std::numeric_limits<double>::quiet_NaN();
    double dH = std::numeric_limits<double>::quiet_NaN();
    if (converged) {
        Hafter = currentTotalEnergy(); // same lambda; H_lambda at the proposed state
        dH = Hafter - Hbefore;
        accept = std::isfinite(dH) && (dH <= 0.0 || uniform_(rng_) < std::exp(-beta_ * dH));
    }
    // else: corrector did not converge at this dt -- F3 forces `accept = false`
    // rather than silently taking a proposal that need not be F-reversible.

    if (ncmcDebugEnabled()) {
        std::fprintf(stderr,
                     "[ncmc-inner] world %d substep %d: h=%.6g Hbefore=%.6f Hafter=%.6f dH=%+.6f "
                     "converged=%s accept=%s\n",
                     index_,
                     substepIndex,
                     (double)h,
                     Hbefore,
                     Hafter,
                     dH,
                     converged ? "true" : "false",
                     accept ? "true" : "false");
    }

    if (!accept) {
        std::copy(q0.begin(), q0.end(), state_.q());
        std::copy(u0.begin(), u0.end(), state_.u());
        for (std::size_t j = 0; j < solv.size(); ++j) {
            posG[solv[j]] = xs0[j];
            velG[solv[j]] = vs0[j];
        }
        // Reject-flip: negate the persistent momenta (GHMC).
        Real* u = state_.u();
        for (int i = 0; i < model_.nu; ++i) {
            u[i] = -u[i];
        }
        for (std::size_t j = 0; j < solv.size(); ++j) {
            robo::Vec3& v = velG[solv[j]];
            v = robo::Vec3(-v[0], -v[1], -v[2]);
        }
        // Re-seed the derivative chain at the restored (q0,-u0,x_s0,-v_s0) state,
        // mirroring ncmcApplyTroughTeleport's post-mutation reseed.
        RobotEngine::realizePosition(model_, state_);
        RobotEngine::fillAtomPositionsFromBodies(model_, state_);
        bridge_.evaluate(state_);
        RobotEngine::realizeVelocity(model_, state_);
        RobotEngine::realizeArticulatedBodyInertias(model_, state_);
        RobotEngine::calcUDot(model_, state_);
        RobotEngine::calcQDot(model_, state_, state_.qdot());
        RobotEngine::calcQDotDot(model_, state_);
    }
    if (acceptedOut) {
        *acceptedOut = accept;
    }
    return true;
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
    // Save the Cartesian solvent positions for a clean rollback on reject: the
    // body->atom fill skips these atoms, so restoring q alone would leave the
    // solvent at its (rejected) end pose. Empty/no-op off solvent-relaxing NCMC.
    {
        const std::vector<int>& solv = state_.cartSolventAtoms();
        const robo::Vec3* posG = state_.atomPosG();
        savedSolventPosG_.resize(solv.size());
        for (std::size_t j = 0; j < solv.size(); ++j) {
            savedSolventPosG_[j] = posG[solv[j]];
        }
    }
    bridge_.setAlchemicalLambda(1.0);
    reinitialize(); // draws p ~ N(0, RT M(q0)); sets Hold_ = V1 + K + F + J
    const double Hstart = Hold_;
    // Diagnostic-only (docs/specs/ncmc-explicit-solvent/40-reproducer-and-oracles.md
    // Sec.2/Sec.6): the Fixman/pitch state-function values at the START, so the
    // reproducer can compute dU_F = fixman_end - fixman_start and dJ from the log
    // without re-deriving them. Zero unless useFixman/useOrientationJacobian.
    const double fixmanStart = state_.energy.fixman;
    const double logSineSqrStart = state_.energy.logSineSqrGamma2;

    // No explicit momentum flip. Momenta are resampled from Maxwell-Boltzmann at
    // the start of every block (reinitialize, above), itself a Gibbs move on the
    // velocity marginal. Under that resampling the NCMC momentum-reversal is
    // unnecessary -- Nilmeier et al. 2011 explicitly sanction reinitializing
    // velocities after each NCMC step -- since the carried-over sign is overwritten
    // next block and never read. (Mirrors the existing MdHmc move, which also
    // resamples and accepts on dH without an explicit flip.)

    const robo::Real h = sampler_.timeStep;
    // DIAGNOSTIC ONLY (NOT used in acceptance): protocol work w = sum over the
    // PERTURBATION substeps of dV at fixed q (Nilmeier et al. 2011, Eq. 16).
    // NOTE: w is NOT equal to Hend - Hstart at finite dt. Two distinct effects open
    // the gap, and NEITHER biases acceptance (acceptance is on Hend - Hstart):
    //   (i)  integrator heat Q != 0 -- the fixed-lambda Verlet steps are only
    //        near-symplectic, so each pumps a little shadow energy that is in
    //        Hend - Hstart but never in w (w only sees the fixed-q perturbations);
    //   (ii) the Fixman + orientation-Jacobian state-function drift F(qend)-F(q0)
    //        and J(qend)-J(q0), which ARE in Hend - Hstart (currentTotalEnergy) but
    //        are excluded from w by construction.
    // So logging gap = w - (Hend - Hstart) is a free consistency check only: a
    // large gap flags a too-large dt / non-converged corrector, not a sampling bug.
    double work = 0.0;
    bool ok = true;
    double Vprev = bridge_.calcPotentialEnergy(); // V at lambda=1, current q
    double lamPrev = 1.0;
    const double peStart = state_.energy.pe; // for the per-move dPE/dKE breakdown
    const double keStart = state_.energy.ke;
    const double keSolvStart = state_.energy.keSolvent; // solvent Cartesian KE at start

    // λ=0 trough teleport (default off; tests/TestNcmcTeleport). Centered in the
    // λ=0 hold so the 1->0->1 protocol stays its own reverse; applied only for an
    // ACYCLIC Free-root region with a defined target region (docking site), where
    // there is no constraint-manifold branch ambiguity (CHMC Thm 3) to guard.
    const int teleHold = static_cast<int>(sampler_.ncmcHoldFraction * sampler_.ncmcSteps);
    const int teleRamp = std::max((sampler_.ncmcSteps - teleHold) / 2, 1);
    const int teleStep = teleRamp + teleHold / 2; // center of the hold
    const int teleRoot = ncmcTeleportRoot();
    const bool doTele = sampler_.ncmcTeleport && teleHold > 0 && teleRoot > 0 && !siteAtoms_.empty()
                        && constraints_.empty();
    if (sampler_.ncmcTeleport && !doTele) {
        std::fprintf(stderr,
                     "[ncmc] teleport requested but inactive (needs hold>0, a Free-root region, a "
                     "docking site, and an acyclic system) -- using the plain uncaged stride\n");
    }

    // Per-move inner-GHMC acceptance count (Construction II only; stays 0/0 for
    // Construction I). Makes inner acceptance visible in the [ncmc] summary
    // line without inferring it from a frozen PE (the symptom that made the
    // Construction-II freeze a black box).
    int innerAccepted = 0;
    int innerAttempted = 0;

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
        // (i.5) TELEPORT at the λ=0 trough (free: ghost, KE preserved by co-rotation).
        if (doTele && s == teleStep) {
            ncmcApplyTroughTeleport(teleRoot);
            Vprev = bridge_.calcPotentialEnergy(); // λ=0 => unchanged; refresh for the heat bookkeeping
        }
        // (ii) PROPAGATE one Verlet step at fixed lambda.
        if (sampler_.useMetropolizedInner) {
            // Construction II (10-acceptance-construction.md Sec.3): the substep
            // is an inner GHMC kernel, Metropolized against the FULL H_lambda, so
            // its shadow work is absorbed into inner rejections and never reaches
            // the outer acceptance. ncmcInnerGhmcStep returns false ONLY on a
            // genuinely unrecoverable non-finite condition (mirrors stepTo).
            bool innerAccept = false;
            ok = ncmcInnerGhmcStep(h, s, &innerAccept);
            ++innerAttempted;
            innerAccepted += innerAccept ? 1 : 0;
        } else {
            // Construction I (endpoint-DeltaH, unchanged): deterministic,
            // unadjusted, reversible Verlet.
            ok = RobotEngine::stepTo(model_, state_, bridge_, constraints_, state_.time + h);
        }
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
        const double fixmanEnd = state_.energy.fixman;
        const double logSineSqrEnd = state_.energy.logSineSqrGamma2;

        // ACCEPTANCE. Two exact constructions (docs/specs/ncmc-explicit-solvent/
        // 10-acceptance-construction.md); SHALL NOT mix them (Sec.2 CLAIM C1):
        //
        // Construction I (endpoint-DeltaH, sampler_.useMetropolizedInner == false,
        // unchanged): this NCMC trajectory is a valid HMC proposal -- a
        // DETERMINISTIC, volume-preserving map T on (q,p) that is momentum-flip
        // reversible, F T F == T^-1. T is the composition of the per-substep
        // fixed-lambda Verlet steps (the fixed-q lambda perturbations do not move
        // the state); each step is F-reversible, so the composition is reversible
        // BECAUSE the lambda schedule is a PALINDROME pinned to lambda = 1 at both
        // endpoints (protocolLambda / NCMCProtocol.hpp). We therefore accept on the
        // FULL Hamiltonian difference at the lambda=1 endpoints (Nilmeier et al.
        // 2011, Eq. 20; bistable-dimer Eq. 28), NOT on the work. GATE: lambda == 1
        // throughout => work == 0 and Hend - Hstart is the plain Verlet dH =>
        // reduces EXACTLY to the torsional-HMC metropolis test.
        //
        // Construction II (Metropolized-dynamics NCMC, useMetropolizedInner ==
        // true): every fixed-lambda substep was ALREADY Metropolized against the
        // full H_lambda inside the loop above (ncmcInnerGhmcStep), so shadow work
        // never reaches this acceptance -- accepting again on Hend - Hstart here
        // would double-count it and also re-admit the very bath shadow work the
        // construction exists to remove (10-...:Sec.3). The OUTER move instead
        // accepts on the protocol work ALONE: a = min(1, exp(-beta*W)) (Sec.3).
        // metropolis(0.0, work) computes exactly that (dH = work - 0 = work).
        const double dH = Hend - Hstart;
        // Per-move diagnostic breakdown (Phase 0): the gap = work - dH flags
        // integrator energy pumping (propagator: dt / mass-scale) plus the Fixman/
        // Jacobian drift, while the recoupling potential change dPE isolates the
        // irreducible reorganization (insertion) cost. dKE should stay ~0. This
        // breakdown is diagnostic under BOTH constructions (gap is never used in
        // Construction II's acceptance either -- only `work` is).
        const double dPE = state_.energy.pe - peStart;
        const double dKE = state_.energy.ke - keStart;
        const double dKEsolv = state_.energy.keSolvent - keSolvStart;
        std::fprintf(stderr,
                     "[ncmc] world %d: construction=%s Hstart=%.2f Hend=%.2f dH=%+.2f kJ/mol  "
                     "dPE=%+.2f dKE=%+.2f dKEsolv=%+.2f work(diag)=%.2f gap=%+.2f "
                     "fixman_start=%.4f fixman_end=%.4f logSineSqr_start=%.4f logSineSqr_end=%.4f "
                     "inner_accepted=%d/%d (steps=%d, teleport=%s)\n",
                     index_,
                     sampler_.useMetropolizedInner ? "II" : "I",
                     Hstart,
                     Hend,
                     dH,
                     dPE,
                     dKE,
                     dKEsolv,
                     work,
                     work - dH,
                     fixmanStart,
                     fixmanEnd,
                     logSineSqrStart,
                     logSineSqrEnd,
                     innerAccepted,
                     innerAttempted,
                     sampler_.ncmcSteps,
                     doTele ? "on" : "off");
        const bool moveAccept = sampler_.useMetropolizedInner
                                     ? (std::isfinite(work) && metropolis(0.0, work))
                                     : (std::isfinite(dH) && metropolis(Hstart, Hend));
        if (moveAccept) {
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
        // Restore the Cartesian solvent pose (the fill above skips these atoms).
        const std::vector<int>& solv = state_.cartSolventAtoms();
        robo::Vec3* posG = state_.atomPosG();
        for (std::size_t j = 0; j < solv.size(); ++j) {
            posG[solv[j]] = savedSolventPosG_[j];
        }
    }
    lastAccepted_ = accepted;
    return accepted;
}