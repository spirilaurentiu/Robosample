#pragma once

// ============================================================================
//  RobotState - the single aligned slab holding every per-step quantity, SoA.
//  Replaces SimTK::State + the SBTree*Cache caches (Simbody01/Simbody/src/
//  SimbodyTreeState.h). Field names track those caches as provenance only;
//  Simbody has been removed from the build, so there is no longer a live
//  1:1 diff against it (see the VALIDATION note in RobotEngine.cpp for the
//  oracle-of-record that replaced it).
//
//  H_FM / H / G are stored as SpatialVec COLUMNS: column j of a body lives at
//  index (bodyUIndex + j). One SpatialVec == 2 Vec3, so an array of `nu`
//  SpatialVec is exactly the 2*nu Vec3 hinge storage Simbody uses.
//
//  DI (the dof x dof articulated hinge inverse) lives in a flat Real pool of
//  size model.nuSq, addressed by model.bodyUSqIndex[b] (row-major dof x dof).
// ============================================================================

#include <ostream>
#include <vector>

#include "MemoryArena.hpp"
#include "RobotModel.hpp"
#include "robot_math.hpp"

struct EnergySnapshot {
    robo::Real pe = 0;
    robo::Real ke = 0; // solute generalized KE  1/2 u^T M(q) u
    // Cartesian kinetic energy 1/2 sum_s m_s |v_s|^2 of the solvent (or contact
    // shell) atoms that are integrated in flat Cartesian space by OpenMM forces
    // INSIDE the proposal (solvent-relaxing NCMC, docs/specs/ncmc_solvent_relax.md).
    // Zero whenever no atom is Cartesian-integrated, so H reduces exactly to the
    // welded torsional-HMC Hamiltonian (T0).
    robo::Real keSolvent = 0;
    robo::Real fixman = 0;
    robo::Real logSineSqrGamma2 = 0;
    robo::Real total = 0;
};

class RobotState {
    public:
    RobotState() = default;

    auto allocateFull(const RobotModel& model) -> void {
        // Reset the arena first so this is safe to call more than once (a World
        // rebuild, e.g. after World::setRootMobility, re-runs buildModel which
        // re-allocates). The arena is one-shot -- reserve() throws after commit()
        // -- so without this a second allocate would abort. Move-assigning a
        // fresh arena frees the previous slab. No-op on the first (default) call.
        arena_ = MemoryArena{};

        nq_ = model.nq;
        nu_ = model.nu;
        nAtoms_ = model.numAtoms;
        nBodies_ = model.numBodies;
        nZRows_ = model.numZRows;
        nuSq_ = model.nuSq;

        hQ_ = arena_.reserve<robo::Real>(nq_);
        hU_ = arena_.reserve<robo::Real>(nu_);
        hQDot_ = arena_.reserve<robo::Real>(nq_);
        hUDot_ = arena_.reserve<robo::Real>(nu_);
        hQDotDot_ = arena_.reserve<robo::Real>(nq_);

        hAtomPosG_ = arena_.reserve<robo::Vec3>(nAtoms_);
        hAtomStationG_ = arena_.reserve<robo::Vec3>(nAtoms_);
        // Per-atom Cartesian velocity / force, used ONLY for atoms that are
        // Cartesian-integrated (solvent-relaxing NCMC). Allocated unconditionally
        // in the full layout (negligible vs the body-level slabs); untouched when
        // no atom is flagged Cartesian (see cartSolventAtoms_ / wantsAtomForces()).
        hAtomVelG_ = arena_.reserve<robo::Vec3>(nAtoms_);
        hAtomForceG_ = arena_.reserve<robo::Vec3>(nAtoms_);

        hX_GB_ = arena_.reserve<robo::Transform>(nBodies_);
        hX_FM_ = arena_.reserve<robo::Transform>(nBodies_);
        hX_PB_ = arena_.reserve<robo::Transform>(nBodies_);
        hPhi_ = arena_.reserve<robo::PhiMatrix>(nBodies_);
        hMk_G_ = arena_.reserve<robo::SpatialInertia>(nBodies_);
        hComG_ = arena_.reserve<robo::Vec3>(nBodies_);
        hH_FM_ = arena_.reserve<robo::SpatialVec>(nu_);
        hH_ = arena_.reserve<robo::SpatialVec>(nu_);
        hV_FM_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hV_PB_G_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hV_GB_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hA_GB_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hGyro_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hCoriolisA_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hMobCoriolisA_ = arena_.reserve<robo::SpatialVec>(nBodies_);

        hP_ = arena_.reserve<robo::ArticulatedInertia>(nBodies_);
        hPPlus_ = arena_.reserve<robo::ArticulatedInertia>(nBodies_);
        hG_ = arena_.reserve<robo::SpatialVec>(nu_);
        hDI_ = arena_.reserve<robo::Real>(nuSq_);
        hABCentrifugal_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hZ_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hZPlus_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hEps_ = arena_.reserve<robo::Real>(nu_);

        hBodyForceG_ = arena_.reserve<robo::SpatialVec>(nBodies_);
        hMobilityForce_ = arena_.reserve<robo::Real>(nu_);
        hBAT_ = arena_.reserve<robo::Vec3>(nZRows_);

        arena_.commit();
        full_ = true;
    }

    auto allocateCompact(const RobotModel& model) -> void {
        arena_ = MemoryArena{}; // re-callable: see allocateFull
        nAtoms_ = model.numAtoms;
        nZRows_ = model.numZRows;
        hAtomPosG_ = arena_.reserve<robo::Vec3>(nAtoms_);
        hBAT_ = arena_.reserve<robo::Vec3>(nZRows_);
        arena_.commit();
        full_ = false;
    }

    [[nodiscard]] auto q() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hQ_);
    }
    [[nodiscard]] auto u() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hU_);
    }
    [[nodiscard]] auto qdot() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hQDot_);
    }
    [[nodiscard]] auto udot() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hUDot_);
    }
    [[nodiscard]] auto qdotdot() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hQDotDot_);
    }

    [[nodiscard]] auto atomPosG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hAtomPosG_);
    }
    [[nodiscard]] auto atomStationG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hAtomStationG_);
    }
    // Ground-frame Cartesian velocity / force per atom (full layout only). Valid
    // for atoms listed in cartSolventAtoms(); other slots are unused scratch.
    [[nodiscard]] auto atomVelG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hAtomVelG_);
    }
    [[nodiscard]] auto atomForceG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hAtomForceG_);
    }

    [[nodiscard]] auto X_GB() const -> robo::Transform* {
        return arena_.ptr<robo::Transform>(hX_GB_);
    }
    [[nodiscard]] auto X_FM() const -> robo::Transform* {
        return arena_.ptr<robo::Transform>(hX_FM_);
    }
    [[nodiscard]] auto X_PB() const -> robo::Transform* {
        return arena_.ptr<robo::Transform>(hX_PB_);
    }
    [[nodiscard]] auto Phi() const -> robo::PhiMatrix* {
        return arena_.ptr<robo::PhiMatrix>(hPhi_);
    }
    [[nodiscard]] auto Mk_G() const -> robo::SpatialInertia* {
        return arena_.ptr<robo::SpatialInertia>(hMk_G_);
    }
    [[nodiscard]] auto comG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hComG_);
    }
    [[nodiscard]] auto H_FM() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hH_FM_);
    }
    [[nodiscard]] auto H() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hH_);
    }
    [[nodiscard]] auto V_FM() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hV_FM_);
    }
    [[nodiscard]] auto V_PB_G() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hV_PB_G_);
    }
    [[nodiscard]] auto V_GB() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hV_GB_);
    }
    [[nodiscard]] auto A_GB() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hA_GB_);
    }
    [[nodiscard]] auto gyro() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hGyro_);
    }
    [[nodiscard]] auto coriolisA() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hCoriolisA_);
    }
    [[nodiscard]] auto mobCoriolisA() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hMobCoriolisA_);
    }

    [[nodiscard]] auto P() const -> robo::ArticulatedInertia* {
        return arena_.ptr<robo::ArticulatedInertia>(hP_);
    }
    [[nodiscard]] auto PPlus() const -> robo::ArticulatedInertia* {
        return arena_.ptr<robo::ArticulatedInertia>(hPPlus_);
    }
    [[nodiscard]] auto G() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hG_);
    }
    [[nodiscard]] auto DI() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hDI_);
    }
    [[nodiscard]] auto abCentrifugal() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hABCentrifugal_);
    }
    [[nodiscard]] auto Z() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hZ_);
    }
    [[nodiscard]] auto zPlus() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hZPlus_);
    }
    [[nodiscard]] auto eps() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hEps_);
    }

    [[nodiscard]] auto bodyForceG() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hBodyForceG_);
    }
    [[nodiscard]] auto mobilityForce() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hMobilityForce_);
    }
    [[nodiscard]] auto BAT() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hBAT_);
    }

    [[nodiscard]] auto nq() const -> int {
        return nq_;
    }
    [[nodiscard]] auto nu() const -> int {
        return nu_;
    }
    [[nodiscard]] auto numAtoms() const -> int {
        return nAtoms_;
    }
    [[nodiscard]] auto numBodies() const -> int {
        return nBodies_;
    }
    [[nodiscard]] auto isFull() const -> bool {
        return full_;
    }

    EnergySnapshot energy;
    robo::Real time = 0;

    // ---- Cartesian-integrated solvent set (per-world, runtime) ---------------
    // The atoms whose flat-space Cartesian (x,v) are advanced by velocity-Verlet
    // inside the proposal (solvent-relaxing NCMC). Empty => the world behaves
    // exactly as the welded torsional engine (every integrator/energy path below
    // is guarded on emptiness), which is what keeps the existing tests bit-exact.
    //   cartSolventAtoms_   : the atom indices (drives the O(n_solv) Verlet loops)
    //   cartSolventInvMass_ : 1/m for each, parallel to cartSolventAtoms_
    //   cartSolventMask_    : per-atom 0/1 for O(1) skip in fillAtomPositionsFromBodies
    void setCartSolvent(const std::vector<int>& atoms, const robo::Real* atomMass) {
        cartSolventAtoms_ = atoms;
        cartSolventInvMass_.resize(atoms.size());
        cartSolventMask_.assign(static_cast<std::size_t>(nAtoms_), 0);
        for (std::size_t j = 0; j < atoms.size(); ++j) {
            const int a = atoms[j];
            cartSolventInvMass_[j] = (atomMass[a] > robo::Real(0)) ? (robo::Real(1) / atomMass[a]) : robo::Real(0);
            if (a >= 0 && a < nAtoms_) {
                cartSolventMask_[static_cast<std::size_t>(a)] = 1;
            }
        }
    }
    [[nodiscard]] auto cartSolventAtoms() const -> const std::vector<int>& {
        return cartSolventAtoms_;
    }
    [[nodiscard]] auto cartSolventInvMass() const -> const std::vector<robo::Real>& {
        return cartSolventInvMass_;
    }
    // nullptr when no atom is Cartesian-integrated (so callers do no work).
    [[nodiscard]] auto cartSolventMask() const -> const char* {
        return cartSolventMask_.empty() ? nullptr : cartSolventMask_.data();
    }
    // True iff per-atom OpenMM forces must be cached into atomForceG() (i.e. some
    // atom is Cartesian-integrated). Lets ForceBridge skip the copy otherwise.
    [[nodiscard]] auto wantsAtomForces() const -> bool {
        return !cartSolventAtoms_.empty();
    }

    auto copyTransferPayloadTo(RobotState& dst) const -> void {
        const robo::Vec3* srcPos = atomPosG();
        robo::Vec3* dstPos = dst.atomPosG();
        for (int i = 0; i < nAtoms_; ++i) {
            dstPos[i] = srcPos[i];
        }
        dst.energy = energy;
    }

    private:
    MemoryArena arena_;
    bool full_ = false;
    int nq_ = 0, nu_ = 0, nAtoms_ = 0, nBodies_ = 0, nZRows_ = 0, nuSq_ = 0;

    std::vector<int> cartSolventAtoms_;
    std::vector<robo::Real> cartSolventInvMass_;
    std::vector<char> cartSolventMask_;

    using Hnd = MemoryArena::Handle;
    Hnd hQ_{}, hU_{}, hQDot_{}, hUDot_{}, hQDotDot_{};
    Hnd hAtomPosG_{}, hAtomStationG_{}, hAtomVelG_{}, hAtomForceG_{};
    Hnd hX_GB_{}, hX_FM_{}, hX_PB_{}, hPhi_{}, hMk_G_{}, hComG_{}, hH_FM_{}, hH_{};
    Hnd hV_FM_{}, hV_PB_G_{}, hV_GB_{}, hA_GB_{}, hGyro_{}, hCoriolisA_{}, hMobCoriolisA_{};
    Hnd hP_{}, hPPlus_{}, hG_{}, hDI_{}, hABCentrifugal_{}, hZ_{}, hZPlus_{}, hEps_{};
    Hnd hBodyForceG_{}, hMobilityForce_{}, hBAT_{};
};

namespace robo {

inline auto operator<<(std::ostream& os, const Vec3& v) -> std::ostream& {
    return os << '[' << v[0] << ' ' << v[1] << ' ' << v[2] << ']';
}

inline auto operator<<(std::ostream& os, const SpatialVec& sv) -> std::ostream& {
    // SpatialVec == 2 Vec3: [0] angular, [1] linear
    return os << "{ang=" << sv[0] << " lin=" << sv[1] << '}';
}

} // namespace robo