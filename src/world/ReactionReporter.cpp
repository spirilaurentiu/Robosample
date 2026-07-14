// ============================================================================
//  ReactionReporter.cpp - World::setReactionReporter / enableReactionReporter
//  / captureReactionSnapshot.
//
//  Relocated verbatim from World.cpp (SPLIT-W3, pure code motion): an opt-in,
//  orthogonal telemetry concern that snapshots, per selected body, the sum of
//  the OpenMM net applied force (bodyForceG) and/or the static (u=0)
//  mobilizer reaction, into a per-frame CSV row buffer. A world that never
//  enables it allocates nothing and does zero extra work.
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <vector>

using robo::Real;

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
    auto reportable = [](int body) {
        return body != 0;
    };
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
