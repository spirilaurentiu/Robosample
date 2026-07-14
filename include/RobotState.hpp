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

/**
 * @brief Additive energy decomposition of one state. Energies are in kJ/mol;
 *        @c logSineSqrGamma2 is dimensionless. Populated by the sampler's energy
 *        evaluation (HmcMove/NcmcMove), not by any realization stage.
 */
struct EnergySnapshot {
    /** @brief Potential energy (OpenMM force field). */
    robo::Real pe = 0;
    /** @brief Solute generalized kinetic energy 1/2 u^T M(q) u. */
    robo::Real ke = 0;
    /**
     * @brief Cartesian kinetic energy 1/2 sum_s m_s |v_s|^2 of the atoms
     *        integrated in flat Cartesian space inside the proposal
     *        (solvent-relaxing NCMC).
     * @note NCMC-solvent field (OQ-6). Zero whenever no atom is
     *       Cartesian-integrated, so the Hamiltonian reduces exactly to the
     *       welded torsional-HMC form.
     */
    robo::Real keSolvent = 0;
    /** @brief Fixman potential contribution (0 when the Fixman term is unused). */
    robo::Real fixman = 0;
    /** @brief log sin^2(gamma) angular-coordinate Jacobian term (dimensionless). */
    robo::Real logSineSqrGamma2 = 0;
    /** @brief Sum consumed by the Metropolis acceptance test for this state. */
    robo::Real total = 0;
};

/**
 * @brief Per-step mutable solver cache for one articulated-body world: every
 *        generalized-coordinate, kinematic, and dynamic quantity the RobotEngine
 *        reads and writes during a step, laid out struct-of-arrays over a single
 *        MemoryArena slab and exposed through single-letter accessors.
 *
 * @par Ownership and lifetime
 * Owned by @c World by value. Every pointer accessor returns a @b borrowed raw
 * pointer into the arena slab; it is valid only while this RobotState is alive
 * and only until the next allocateFull()/allocateCompact() call, which replaces
 * the arena, frees the previous slab, and invalidates every outstanding pointer.
 * The accessors are @c const yet return @b mutable pointers: the solver writes
 * the caches in place through them, so @c const here means the handle table is
 * unchanged, not that the pointed-to data is read-only.
 *
 * @par INV-4 stage-validity contract
 * A cache holds meaningful data only after the realization stage that writes it
 * has run for the @e current q()/u(); reading it earlier returns the previous
 * step's values or (on the first step) uninitialized slab memory, never an
 * error. The stages run in a fixed order, each depending on all earlier ones:
 *   - @b Stage-P (position): RobotEngine::realizePosition. Requires a current q().
 *   - @b Stage-V (velocity): RobotEngine::realizeVelocity. Requires Stage-P and u().
 *   - @b Stage-A (articulated-body inertias): factorizeArticulatedInertias
 *        (position-only: P, PPlus, G, DI) then seedArticulatedCentrifugal
 *        (velocity-coupled: abCentrifugal); realizeArticulatedBodyInertias runs
 *        both. Requires Stage-P for the factorization and Stage-V for the seed.
 *   - @b Stage-U (accelerations): RobotEngine::calcUDot. Requires Stage-A plus
 *        the force inputs bodyForceG() and mobilityForce().
 * Each accessor below states the earliest stage at which its cache is valid, or
 * marks itself a caller-owned input written before a stage rather than by one.
 * The contract is enforced by call-order convention, not by the type system
 * (ARCHITECTURE OQ-1); it is recovered from the write sites in
 * RobotEngine_kinematics.cpp and RobotEngine_dynamics.cpp, not from the names.
 *
 * @see MemoryArena for slab/handle mechanics and reset semantics.
 */
class RobotState {
    public:
    RobotState() = default;

    /**
     * @brief Allocate the full per-step cache layout (all backing arrays) sized
     *        from @p model.
     * @param[in] model immutable topology; borrowed for the call only. Supplies
     *        the nq/nu/numAtoms/numBodies/numZRows/nuSq array extents.
     * @post Every accessor returns a fresh arena pointer into uninitialized slab
     *       memory; each cache becomes meaningful only when its realization stage
     *       runs. isFull() == true.
     * @note Re-callable: resets the arena first, freeing the previous slab and
     *       invalidating all previously handed-out pointers. Used on a World
     *       rebuild (e.g. after World::setRootMobility re-runs buildModel).
     */
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

    /**
     * @brief Allocate the minimal transfer-only layout: atomPosG() and BAT() only.
     * @param[in] model borrowed; supplies numAtoms and numZRows.
     * @post isFull() == false; only atomPosG() and BAT() are backed by the arena.
     *       Accessing any full-layout cache is undefined.
     * @note Re-callable; resets the arena as allocateFull() does. Used for states
     *       that only carry the per-atom transfer payload between worlds.
     */
    auto allocateCompact(const RobotModel& model) -> void {
        arena_ = MemoryArena{}; // re-callable: see allocateFull
        nAtoms_ = model.numAtoms;
        nZRows_ = model.numZRows;
        hAtomPosG_ = arena_.reserve<robo::Vec3>(nAtoms_);
        hBAT_ = arena_.reserve<robo::Vec3>(nZRows_);
        arena_.commit();
        full_ = false;
    }

    /**
     * @brief Generalized coordinates (position DOF), @c nq entries.
     * @note Caller-owned input to Stage-P: the sampler/integrator writes q, then
     *       realizePosition consumes it. Quaternion blocks are kept
     *       unit-normalized (RobotEngine::normalizeQuaternions).
     */
    [[nodiscard]] auto q() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hQ_);
    }
    /**
     * @brief Generalized speeds, @c nu entries.
     * @note Caller-owned input to Stage-V: written by the integrator/sampler,
     *       read by realizeVelocity.
     */
    [[nodiscard]] auto u() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hU_);
    }
    /**
     * @brief Time derivatives of the generalized coordinates, @c nq entries.
     * @pre Stage-V (RobotEngine::calcQDot, run first inside realizeVelocity);
     *      maps u through the per-joint N-matrix.
     */
    [[nodiscard]] auto qdot() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hQDot_);
    }
    /**
     * @brief Generalized accelerations, @c nu entries.
     * @pre Stage-U (RobotEngine::calcUDot) for the current step's value; the
     *      integrator also reads the previous step's udot as a drift seed.
     */
    [[nodiscard]] auto udot() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hUDot_);
    }
    /**
     * @brief Second time derivatives of the generalized coordinates, @c nq entries.
     * @pre RobotEngine::calcQDotDot, which runs after Stage-U (it reads udot()).
     */
    [[nodiscard]] auto qdotdot() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hQDotDot_);
    }

    /**
     * @brief Per-atom Cartesian position in the Ground frame (nm), @c numAtoms
     *        entries.
     * @pre Stage-P for atoms placed by a rigid body. Slots for
     *      Cartesian-integrated NCMC-solvent atoms (cartSolventAtoms()) are owned
     *      directly by the integrator and are @b not overwritten by Stage-P.
     * @note This array is the inter-world transfer currency (INV-3) and the
     *       per-atom input to the OpenMM force bridge.
     */
    [[nodiscard]] auto atomPosG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hAtomPosG_);
    }
    /**
     * @brief Per-atom station vector R_GB * station_B in Ground (nm), @c numAtoms
     *        entries: the atom's offset from its body origin rotated into Ground.
     * @pre Stage-P.
     */
    [[nodiscard]] auto atomStationG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hAtomStationG_);
    }
    /**
     * @brief Per-atom Cartesian velocity in Ground (nm/ps), @c numAtoms entries;
     *        full layout only.
     * @note NCMC-solvent field (OQ-6). Meaningful only for atoms in
     *       cartSolventAtoms(); other slots are unused scratch. Owned and advanced
     *       by the integrator's Cartesian-solvent Verlet, not by any realization
     *       stage.
     */
    [[nodiscard]] auto atomVelG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hAtomVelG_);
    }
    /**
     * @brief Per-atom Cartesian force in Ground (kJ/mol/nm), @c numAtoms entries;
     *        full layout only.
     * @note NCMC-solvent field (OQ-6). Written by ForceBridge only when
     *       wantsAtomForces(); consumed by the Cartesian-solvent Verlet. Unused
     *       scratch for non-solvent slots.
     */
    [[nodiscard]] auto atomForceG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hAtomForceG_);
    }

    /**
     * @brief Body-to-Ground transforms, @c numBodies entries; index 0 is Ground
     *        (identity).
     * @pre Stage-P.
     */
    [[nodiscard]] auto X_GB() const -> robo::Transform* {
        return arena_.ptr<robo::Transform>(hX_GB_);
    }
    /**
     * @brief Cross-mobilizer transforms (inboard F frame to outboard M frame),
     *        @c numBodies entries.
     * @pre Stage-P.
     */
    [[nodiscard]] auto X_FM() const -> robo::Transform* {
        return arena_.ptr<robo::Transform>(hX_FM_);
    }
    /**
     * @brief Parent-body-to-child-body transforms, @c numBodies entries.
     * @pre Stage-P.
     */
    [[nodiscard]] auto X_PB() const -> robo::Transform* {
        return arena_.ptr<robo::Transform>(hX_PB_);
    }
    /**
     * @brief Rigid shift operators from each body to its parent (PhiMatrix),
     *        @c numBodies entries.
     * @pre Stage-P.
     */
    [[nodiscard]] auto Phi() const -> robo::PhiMatrix* {
        return arena_.ptr<robo::PhiMatrix>(hPhi_);
    }
    /**
     * @brief Per-body spatial inertia in Ground about the body origin,
     *        @c numBodies entries.
     * @pre Stage-P.
     * @note Carries any per-body kinetic-metric mass scaling
     *       (RobotModel::bodyMassScale); the draw, KE, and Fixman paths all read
     *       Mk_G so they stay mutually consistent.
     */
    [[nodiscard]] auto Mk_G() const -> robo::SpatialInertia* {
        return arena_.ptr<robo::SpatialInertia>(hMk_G_);
    }
    /**
     * @brief Per-body center-of-mass position in Ground (nm), @c numBodies entries.
     * @pre Stage-P.
     */
    [[nodiscard]] auto comG() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hComG_);
    }
    /**
     * @brief Joint (hinge) matrix columns in the mobilizer M frame, stored as
     *        SpatialVec columns indexed by RobotModel::bodyUIndex[b]+j; @c nu
     *        columns total.
     * @pre Stage-P.
     */
    [[nodiscard]] auto H_FM() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hH_FM_);
    }
    /**
     * @brief Joint matrix columns expressed in Ground (H_PB_G), stored as
     *        SpatialVec columns indexed by bodyUIndex[b]+j; @c nu columns total.
     * @pre Stage-P.
     */
    [[nodiscard]] auto H() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hH_);
    }
    /**
     * @brief Cross-mobilizer spatial velocity in the F frame, @c numBodies entries.
     * @pre Stage-V.
     */
    [[nodiscard]] auto V_FM() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hV_FM_);
    }
    /**
     * @brief Cross-mobilizer (parent-to-body) spatial velocity in Ground,
     *        @c numBodies entries.
     * @pre Stage-V.
     */
    [[nodiscard]] auto V_PB_G() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hV_PB_G_);
    }
    /**
     * @brief Body spatial velocity in Ground (angular, linear), @c numBodies
     *        entries; index 0 (Ground) is zero.
     * @pre Stage-V.
     */
    [[nodiscard]] auto V_GB() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hV_GB_);
    }
    /**
     * @brief Body spatial acceleration in Ground, @c numBodies entries; index 0
     *        (Ground) is zero.
     * @pre Stage-U (second, outward pass of calcUDot).
     */
    [[nodiscard]] auto A_GB() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hA_GB_);
    }
    /**
     * @brief Per-body gyroscopic (velocity-quadratic) spatial force in Ground,
     *        @c numBodies entries.
     * @pre Stage-V.
     */
    [[nodiscard]] auto gyro() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hGyro_);
    }
    /**
     * @brief Total Coriolis/centrifugal bias acceleration per body in Ground
     *        (~Phi*a_parent + mobilizer term), @c numBodies entries.
     * @pre Stage-V.
     */
    [[nodiscard]] auto coriolisA() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hCoriolisA_);
    }
    /**
     * @brief Mobilizer-incremental Coriolis bias acceleration per body in Ground,
     *        @c numBodies entries: the per-joint term A, not the propagated total.
     * @pre Stage-V.
     */
    [[nodiscard]] auto mobCoriolisA() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hMobCoriolisA_);
    }

    /**
     * @brief Articulated-body inertia per body, @c numBodies entries.
     * @pre Stage-A (factorizeArticulatedInertias); a pure function of q.
     */
    [[nodiscard]] auto P() const -> robo::ArticulatedInertia* {
        return arena_.ptr<robo::ArticulatedInertia>(hP_);
    }
    /**
     * @brief Articulated-body inertia with the mobility subspace removed (P+),
     *        @c numBodies entries.
     * @pre Stage-A (factorizeArticulatedInertias).
     */
    [[nodiscard]] auto PPlus() const -> robo::ArticulatedInertia* {
        return arena_.ptr<robo::ArticulatedInertia>(hPPlus_);
    }
    /**
     * @brief Articulated gain columns P*H*DI, stored as SpatialVec columns indexed
     *        by RobotModel::bodyUIndex[b]+j; @c nu columns total.
     * @pre Stage-A (factorizeArticulatedInertias).
     */
    [[nodiscard]] auto G() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hG_);
    }
    /**
     * @brief Inverse hinge inertia D^-1, a flat Real pool of @c nuSq entries;
     *        body b's dof x dof block starts at RobotModel::bodyUSqIndex[b],
     *        row-major.
     * @pre Stage-A (factorizeArticulatedInertias).
     */
    [[nodiscard]] auto DI() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hDI_);
    }
    /**
     * @brief Articulated centrifugal seed P*a_mob + gyro for the calcUDot inward
     *        pass, @c numBodies entries.
     * @pre Stage-A centrifugal seed (seedArticulatedCentrifugal), which requires
     *      Stage-V (a_mob, gyro) in addition to the factorization.
     */
    [[nodiscard]] auto abCentrifugal() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hABCentrifugal_);
    }
    /**
     * @brief Articulated residual spatial force per body (calcUDot inward pass),
     *        @c numBodies entries.
     * @pre Stage-U.
     */
    [[nodiscard]] auto Z() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hZ_);
    }
    /**
     * @brief Articulated residual force after the mobility projection (calcUDot
     *        inward pass), @c numBodies entries.
     * @pre Stage-U.
     */
    [[nodiscard]] auto zPlus() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hZPlus_);
    }
    /**
     * @brief Per-DOF hinge residual f - H^T z (calcUDot inward pass), @c nu entries.
     * @pre Stage-U.
     */
    [[nodiscard]] auto eps() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hEps_);
    }

    /**
     * @brief Applied per-body spatial wrench (angular, linear) in Ground about the
     *        body origin (kJ/mol, kJ/mol/nm), @c numBodies entries.
     * @note Caller-owned input to Stage-U: ForceBridge zeroes then accumulates the
     *       reduced OpenMM forces here (INV-1) before calcUDot reads it. Must be
     *       current for the step's q(). The integrator also scans it for
     *       non-finite forces (steric-clash rejection).
     */
    [[nodiscard]] auto bodyForceG() const -> robo::SpatialVec* {
        return arena_.ptr<robo::SpatialVec>(hBodyForceG_);
    }
    /**
     * @brief Applied generalized (joint) force per DOF, @c nu entries.
     * @note Caller-owned input to Stage-U: ForceBridge clears it every step (this
     *       model applies no direct joint torque, so it is identically zero unless
     *       a Fixman/bias term adds into it after the force eval); calcUDot reads
     *       it as f in eps = f - H^T z.
     */
    [[nodiscard]] auto mobilityForce() const -> robo::Real* {
        return arena_.ptr<robo::Real>(hMobilityForce_);
    }
    /**
     * @brief Per-Z-row bond/angle/torsion internal coordinates, @c numZRows entries.
     * @warning No caller in the current source tree reads or writes through this
     *          accessor; the backing slab is allocated (in both layouts) but
     *          otherwise unused. See findings (dead accessor).
     */
    [[nodiscard]] auto BAT() const -> robo::Vec3* {
        return arena_.ptr<robo::Vec3>(hBAT_);
    }

    /** @brief Number of generalized coordinates (length of q()/qdot()/qdotdot()). */
    [[nodiscard]] auto nq() const -> int {
        return nq_;
    }
    /** @brief Number of generalized speeds / DOF (length of u()/udot()/eps()). */
    [[nodiscard]] auto nu() const -> int {
        return nu_;
    }
    /** @brief Number of atoms (length of the per-atom arrays). */
    [[nodiscard]] auto numAtoms() const -> int {
        return nAtoms_;
    }
    /** @brief Number of bodies, including Ground at index 0. */
    [[nodiscard]] auto numBodies() const -> int {
        return nBodies_;
    }
    /**
     * @brief True after allocateFull(); false after allocateCompact()
     *        (transfer-only layout).
     */
    [[nodiscard]] auto isFull() const -> bool {
        return full_;
    }

    /**
     * @brief Latest energy decomposition for this state.
     * @note Written by the sampler's energy evaluation, not by a realization stage.
     */
    EnergySnapshot energy;
    /** @brief Simulation time accumulated by the integrator (ps). */
    robo::Real time = 0;

    // ---- Cartesian-integrated solvent set (per-world, runtime) ---------------
    // The atoms whose flat-space Cartesian (x,v) are advanced by velocity-Verlet
    // inside the proposal (solvent-relaxing NCMC). Empty => the world behaves
    // exactly as the welded torsional engine (every integrator/energy path below
    // is guarded on emptiness), which is what keeps the existing tests bit-exact.
    //   cartSolventAtoms_   : the atom indices (drives the O(n_solv) Verlet loops)
    //   cartSolventInvMass_ : 1/m for each, parallel to cartSolventAtoms_
    //   cartSolventMask_    : per-atom 0/1 for O(1) skip in fillAtomPositionsFromBodies
    /**
     * @brief Designate the atom subset advanced in flat Cartesian space by the
     *        in-proposal solvent-relaxing NCMC Verlet.
     * @param[in] atoms atom indices to Cartesian-integrate; copied into the state.
     * @param[in] atomMass borrowed per-atom mass array (Daltons), indexable at
     *        every entry of @p atoms; read once to precompute inverse masses.
     *        Zero-mass atoms receive inverse mass 0.
     * @post cartSolventAtoms(), cartSolventInvMass(), and cartSolventMask() reflect
     *       @p atoms; wantsAtomForces() becomes true iff @p atoms is non-empty.
     * @note NCMC-solvent field (OQ-6): this runtime set is concern-bleed into the
     *       core cache. An empty set leaves every downstream integrator/energy path
     *       bit-identical to the welded torsional engine.
     */
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
    /**
     * @brief The Cartesian-integrated atom indices set by setCartSolvent().
     * @note NCMC-solvent field (OQ-6). Empty in the welded engine.
     */
    [[nodiscard]] auto cartSolventAtoms() const -> const std::vector<int>& {
        return cartSolventAtoms_;
    }
    /**
     * @brief Inverse masses (1/Da) parallel to cartSolventAtoms().
     * @note NCMC-solvent field (OQ-6).
     */
    [[nodiscard]] auto cartSolventInvMass() const -> const std::vector<robo::Real>& {
        return cartSolventInvMass_;
    }
    /**
     * @brief Per-atom 0/1 membership mask over all @c numAtoms atoms, for an O(1)
     *        skip in fillAtomPositionsFromBodies.
     * @return borrowed pointer to @c numAtoms bytes, or @c nullptr when no atom is
     *         Cartesian-integrated (callers then do no work).
     * @note NCMC-solvent field (OQ-6).
     */
    [[nodiscard]] auto cartSolventMask() const -> const char* {
        return cartSolventMask_.empty() ? nullptr : cartSolventMask_.data();
    }
    /**
     * @brief True iff some atom is Cartesian-integrated, so ForceBridge must cache
     *        per-atom forces into atomForceG(); false lets it skip the copy.
     * @note NCMC-solvent field (OQ-6).
     */
    [[nodiscard]] auto wantsAtomForces() const -> bool {
        return !cartSolventAtoms_.empty();
    }

    /**
     * @brief Copy the inter-world transfer payload into @p dst: the per-atom Ground
     *        positions and the energy snapshot.
     * @param[out] dst destination state; must already be allocated with the same
     *        numAtoms() as this state. Only atomPosG() and energy are written; all
     *        other caches in @p dst are left untouched.
     * @pre Stage-P current here, so atomPosG() holds valid positions.
     */
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

/** @brief Stream a Vec3 as @c [x y z] (diagnostics only). */
inline auto operator<<(std::ostream& os, const Vec3& v) -> std::ostream& {
    return os << '[' << v[0] << ' ' << v[1] << ' ' << v[2] << ']';
}

/** @brief Stream a SpatialVec as @c {ang=[...] lin=[...]}; element 0 is angular,
 *  element 1 linear (diagnostics only). */
inline auto operator<<(std::ostream& os, const SpatialVec& sv) -> std::ostream& {
    // SpatialVec == 2 Vec3: [0] angular, [1] linear
    return os << "{ang=" << sv[0] << " lin=" << sv[1] << '}';
}

} // namespace robo