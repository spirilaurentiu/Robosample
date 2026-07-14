// ============================================================================
//  CartesianSolvent.cpp - World::setCartesianSolvent / drawSolventVelocities
//  / calcSolventKE.
//
//  Relocated verbatim from World.cpp (SPLIT-W7, pure code motion): the
//  Cartesian-solvent concern (docs/specs/ncmc_solvent_relax.md,
//  docs/specs/two-robot-contact/) -- marking a set of atoms to be advanced in
//  flat Cartesian space by velocity-Verlet inside the proposal, plus the
//  Maxwell-Boltzmann velocity draw and the flat-metric kinetic energy for
//  those atoms. The atoms stay welded as 0-DOF rigid bodies (no
//  Fixman/Jacobian contribution); only their per-atom (x,v) move. An empty
//  set is the welded engine, bit-for-bit.
// ============================================================================

#include "World.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

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
