// ============================================================================
//  DockingMove.cpp - World's docking concern: the auto-sized binding-sphere
//  geometry, the rigid-kick proposal (uniform position in the sphere +
//  uniform SO(3) reorientation of each ligand), and the clash-free
//  initial-placement search.
//
//  Relocated verbatim from World.cpp (SPLIT-W2, pure code motion). The move
//  itself (accept/reject) lives in HMC (W9); this TU is the proposal
//  machinery plus the docking configuration entry point.
//
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
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "engine_helpers.hpp"

using robo::Real;
using robo::Rotation;
using robo::Vec3;

namespace {

// Runtime toggle for per-kick docking diagnostics, read once. Set the env var
//   ROBO_DOCK_DEBUG=1
// to print, every kick, the sphere geometry and the sampled placement in GROUND
// (Cartesian) coordinates. All lengths below are nanometres (the engine's and
// OpenMM's native unit; AMBER Angstrom inputs are converted by ANG_TO_NM on load).
//
// NOTE (SPLIT-W2/CX-8): this helper appears unused elsewhere in the docking
// TU (no call site found by grep). Moved verbatim with the docking concept
// regardless -- dead code observation flagged for the reviewer, not deleted
// in this code-motion ticket.
bool dockDebugEnabled() {
    static const bool on = [] {
        const char* e = std::getenv("ROBO_DOCK_DEBUG");
        return e != nullptr && e[0] != '0' && e[0] != '\0';
    }();
    return on;
}

} // namespace

void World::configureDocking(std::vector<std::vector<int>> ligandGroups, std::vector<int> siteAtoms) {
    docking_ = true;
    ligandGroups_ = std::move(ligandGroups);
    siteAtoms_ = std::move(siteAtoms);
    sampler_.moveType = MoveType::RigidKick;
}

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
    return EngineHelpers::quatToRotation(qw, qx, qy, qz);
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
