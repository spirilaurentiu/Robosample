#pragma once

/**
 * @brief The sole adapter between `robo::` types and OpenMM: the concrete
 *        OpenMM-host (and fused-CUDA) `Bridge` the templated integrator drives.
 *
 * @par Direction contract
 *      Robosample -> OpenMM: per-atom Cartesian positions in Ground (nm).
 *      OpenMM -> Robosample: per-atom Cartesian forces (kJ/mol/nm), which this
 *      bridge reduces to per-body spatial wrenches itself (via `ForceReducer`,
 *      INV-1/INV-2). Atom order == OpenMM particle order == global/BFS order, so
 *      the position transfer is an index-free copy and the force->body scatter
 *      uses only `model.atomBody[a]`.
 *
 * @par The Bridge concept (template requirement)
 *      `RobotIntegrator::verletStep<Bridge>` is templated on the bridge type and
 *      requires exactly one call: `bridge.evaluate(RobotState& s)`. Evidenced by
 *      all three instantiations - this class, `tests/AnalyticForceBridge.hpp`,
 *      and the OpenMM-free `tests/ForceBridge.hpp` stub - a conforming `Bridge`
 *      SHALL, in `evaluate(s)`: read `s.atomPosG()` and `s.X_GB()` (position
 *      already realized by the integrator), overwrite `s.bodyForceG()` with the
 *      per-body wrench under INV-1 (`angular = moment about the body origin,
 *      linear = net force`, Ground) and INV-2 (skip massless slots), zero
 *      `s.mobilityForce()` (no direct joint forces in this model), and leave the
 *      potential energy retrievable via `calcPotentialEnergy(s)`.
 *
 * @note Ownership: constructed by `World` by value; an observer of `RobotModel`
 *       (holds a reference) and an adapter to the `OpenMMContext` singleton. Owns
 *       only its staging buffers. Header-only: every method is defined in-class
 *       (implicitly inline), so there is one definition across all TUs.
 * @note The header still includes `OpenMMContext.hpp` concretely (the forwarders
 *       and `posCache_` use OpenMM types); the ARCHITECTURE 6.1
 *       forward-declaration resolution is not applied here - recorded in findings.
 */

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "OpenMMContext.hpp" // the existing OpenMM singleton (reused, not wrapped)
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "bridge/ForceReducer.hpp" // the single host force->wrench reduction (INV-1/INV-2)
#include "robot_math.hpp"

class ForceBridge {
    public:
    /// @brief Binds the bridge to @p model (observed by reference; must outlive
    ///        the bridge). Allocates nothing until first evaluate.
    explicit ForceBridge(const RobotModel& model)
        : model_(model) {
    }

    /**
     * @brief Caches @p s 's Ground-frame positions in OpenMM particle order.
     *
     * The actual `context->setPositions` happens inside `OpenMMContext` on the
     * next energy/force call, which consumes this cache. Marks the last
     * evaluation as non-fused so `calcPotentialEnergy` recomputes from the cache.
     *
     * @param[in] s  Position-realized state; `s.atomPosG()` is read.
     */
    void setAtomPositionsInGround(const RobotState& s) {
        const robo::Vec3* p = s.atomPosG();
        posCache_.resize(static_cast<std::size_t>(model_.numAtoms));
        for (int a = 0; a < model_.numAtoms; ++a) {
            posCache_[a] = OpenMM::Vec3(p[a][0], p[a][1], p[a][2]);
        }
        lastEvalWasFused_ = false; // host positions now current; PE comes from posCache_
    }

    /**
     * @brief Potential energy for the current positions [kJ/mol].
     * @return On the fused CUDA path, the energy computed on device during the
     *         last `evaluate()`; otherwise OpenMM recomputes it from the cached
     *         positions.
     * @pre A prior `setAtomPositionsInGround`/`evaluate` established the positions.
     */
    [[nodiscard]] robo::Real calcPotentialEnergy() const {
        if (lastEvalWasFused_) {
            return lastPE_;
        }
        return OpenMMContext::get().computePotentialEnergy(posCache_);
    }

    /**
     * @brief OpenMM -> Robosample: evaluates per-atom forces for the cached
     *        positions and reduces them to per-body spatial wrenches on @p s.
     *
     * The host force path. Overwrites `s.bodyForceG()` with `(moment about the
     * body origin, net force)` in Ground via `reduceAtomForcesToBodies`
     * (INV-1/INV-2), and zeros `s.mobilityForce()` (this model applies forces
     * only as Cartesian atom forces; a stale mobility force would feed `calcUDot`
     * garbage). When `s.wantsAtomForces()` (solvent-relaxing NCMC), the raw
     * per-atom Cartesian force is additionally cached into `s.atomForceG()` under
     * the same INV-2 skip.
     *
     * @param[in,out] s  Force targets `bodyForceG`/`mobilityForce` (and optionally
     *                   `atomForceG`) are written; positions come from the cache.
     * @pre `setAtomPositionsInGround(s)` ran with the current positions.
     */
    void getForcesFromOpenMM(RobotState& s) const {
        std::vector<OpenMM::Vec3> forces;
        OpenMMContext::get().evaluateForcesFromPositionsCache(posCache_, forces);

        robo::SpatialVec* BF = s.bodyForceG();
        for (int b = 0; b < model_.numBodies; ++b) {
            BF[b] = robo::SpatialVec(robo::Vec3(0), robo::Vec3(0));
        }
        // Generalized (joint) forces. calcUDot reads mobilityForce as the
        // applied per-DOF force (eps = f - H^T z). This model applies forces
        // ONLY as Cartesian atom forces (reduced into bodyForceG above); there
        // are no direct joint torques (Fixman unused), so the generalized force
        // is identically zero. It lives in the arena and is never written
        // elsewhere, so it MUST be cleared here every step -- otherwise
        // calcUDot reads uninitialized slab memory -> garbage udot -> NaN.
        // (Any future Fixman/biasing term should add into this cleared vector
        //  after evaluate() and before calcUDot.)
        robo::Real* mob = s.mobilityForce();
        for (int i = 0; i < model_.nu; ++i) {
            mob[i] = robo::Real(0);
        }
        // When some atoms are Cartesian-integrated inside the proposal
        // (solvent-relaxing NCMC), the velocity-Verlet for those atoms needs the
        // raw per-atom Cartesian force. Cache it here (the only place the per-atom
        // force vector exists) so the integrator can read it back. Skipped
        // entirely otherwise -- the welded engine never allocates/reads it.
        const bool cacheAtomForces = s.wantsAtomForces();
        robo::Vec3* atomForce = cacheAtomForces ? s.atomForceG() : nullptr;
        const robo::Vec3* posG = s.atomPosG();
        const robo::Transform* X_GB = s.X_GB();

        // Convert OpenMM's per-atom forces to robo::Vec3 (reduceAtomForcesToBodies
        // is OpenMM-free) and cache them for the solvent Verlet where requested,
        // under the SAME mass==0 skip the reduction below applies (a virtual site
        // never needs its own Cartesian-Verlet force cached).
        std::vector<robo::Vec3> forceG(static_cast<std::size_t>(model_.numAtoms));
        for (int a = 0; a < model_.numAtoms; ++a) {
            const robo::Vec3 f(forces[a][0], forces[a][1], forces[a][2]);
            forceG[static_cast<std::size_t>(a)] = f;
            if (atomForce && model_.atomMass[a] != robo::Real(0)) {
                atomForce[a] = f; // per-atom Cartesian force for the solvent Verlet
            }
        }

        // OpenMM -> Robosample body-wrench reduction (INV-1/INV-2), hoisted so the
        // CUDA reduceForces kernel is pinned against this one definition
        // (tests/TestForceReducer.cpp). Virtual sites (mass == 0, e.g. the
        // OPC/TIP4P M-site) carry no independent DOF: OpenMM's getState(Forces)
        // has ALREADY projected the force computed at the site onto its real
        // parent atoms (verify: sum over REAL atoms reproduces the system net
        // force / -dPE/dq to machine precision, while the raw force still left in
        // the site's own slot makes the global sum nonzero). Reducing that
        // leftover slot too would DOUBLE-COUNT the site force -- here it would
        // re-apply ~10^3 kJ/mol/nm per water at the M-site station, a
        // non-conservative kick that pumps the body's kinetic energy every step.
        // reduceAtomForcesToBodies skips it; the parents already carry the
        // contribution.
        reduceAtomForcesToBodies(forceG.data(), posG, model_.atomMass.data(), model_.atomBody.data(), X_GB,
                                 model_.numAtoms, model_.numBodies, BF);
    }

    /**
     * @brief The `Bridge` concept entry point: one full force evaluation
     *        (positions then forces) that leaves the per-body wrench on @p s.
     *
     * Selects the fused CUDA path (positions and force reduction on device, only
     * `O(numBodies)` host<->device traffic per step) when
     * `cudaKinematicsAvailable()` and `!s.wantsAtomForces()`; otherwise the exact
     * host path (`setAtomPositionsInGround` then `getForcesFromOpenMM`). Both
     * paths honor the same INV-1/INV-2 contract and zero `mobilityForce`.
     *
     * @param[in,out] s  Position-realized state; `bodyForceG` and `mobilityForce`
     *                   are written. Potential energy is then via
     *                   `calcPotentialEnergy()`.
     * @note The Cartesian MD world does not reach here; it uses
     *       `integrateTrajectoryOnDevice`.
     */
    void evaluate(RobotState& s) {
        if (OpenMMContext::get().cudaKinematicsAvailable() && !s.wantsAtomForces()) {
            evaluateFused(s);
        } else {
            setAtomPositionsInGround(s);
            getForcesFromOpenMM(s);
        }
    }

    /// @brief Marks the atom stations dirty so the fused path re-uploads them on
    ///        its next `evaluate`. Call after any refit of
    ///        `RobotModel::atomStation_B` (a coordinate transfer). No-op for the
    ///        host path.
    void markStationsDirty() {
        stationsDirty_ = true;
    }

    /// @brief Draws Maxwell-Boltzmann velocities on the OpenMM Context.
    /// @param[in] temperature  Kelvin. @param[in] seed  RNG seed.
    void setVelocitiesToTemperature(double temperature, int seed) const {
        OpenMMContext::get().setVelocitiesToTemperature(temperature, seed);
    }

    /// @brief Forwards Region A (arbitrary atom-index set) to the alchemy factory.
    ///        @see AlchemyForceFactory for lambda semantics. Call before the
    ///        OpenMM system is initialized.
    void enableAlchemy(const std::vector<int>& atomIndices) const {
        OpenMMContext::get().enableAlchemy(atomIndices);
    }
    /// @brief Convenience overload: Region A = contiguous `[atomBegin, atomEnd)`.
    void enableAlchemy(int atomBegin, int atomEnd) const {
        OpenMMContext::get().enableAlchemy(atomBegin, atomEnd);
    }
    /// @brief Sets the alchemical coupling `lambda_inter` (1 = coupled/unmodified
    ///        field, 0 = A<->rest decoupled). @see AlchemyForceFactory.
    void setAlchemicalLambda(double lambdaInter) const {
        OpenMMContext::get().setAlchemicalLambda(lambdaInter);
    }

    /**
     * @brief Runs on-device Cartesian MD for @p steps and writes the resulting
     *        positions back into @p s (the Cartesian-world path).
     *
     * Pushes `s.atomPosG()` to the live Context, integrates, then pulls the
     * post-integration coordinates back into `s.atomPosG()`.
     *
     * @param[in,out] s         Positions read then overwritten with the result.
     * @param[in]     steps     Integration steps.
     * @param[in]     timestep  Step size [ps].
     * @return `true` on success; `false` iff the integrator threw (positions are
     *         still read back from the post-step state).
     */
    bool integrateTrajectoryOnDevice(RobotState& s, int steps, robo::Real timestep) {
        setAtomPositionsInGround(s);                  // posCache_ <- s.atomPosG()
        OpenMMContext::get().setPositions(posCache_); // push to the live context
        const bool ok =
            OpenMMContext::get().integrateTrajectory(static_cast<int>(steps), static_cast<double>(timestep));
        OpenMMContext::get().getPositions(posCache_); // pull post-integration coords
        robo::Vec3* p = s.atomPosG();
        for (int a = 0; a < model_.numAtoms; ++a) {
            p[a] = robo::Vec3(posCache_[a][0], posCache_[a][1], posCache_[a][2]);
        }
        return ok;
    }

    private:
    // Allocate the reused staging buffers once and build the constant virtual-site mask.
    // stationFlat_ is only sized here; its CONTENTS are repacked every step (atomStation_B
    // is refit per coordinate transfer, so it changes between rounds).
    void ensureFusedStaging() {
        if (fusedStagingBuilt_) {
            return;
        }
        stationFlat_.assign(static_cast<std::size_t>(3 * model_.numAtoms), 0.0);
        isVirtualI_.resize(static_cast<std::size_t>(model_.numAtoms));
        for (int a = 0; a < model_.numAtoms; ++a) {
            isVirtualI_[static_cast<std::size_t>(a)] =
                (model_.atomMass[static_cast<std::size_t>(a)] == robo::Real(0)) ? 1 : 0;
        }
        xgbFlat_.assign(static_cast<std::size_t>(12 * model_.numBodies), 0.0);
        bodyForceFlat_.assign(static_cast<std::size_t>(6 * model_.numBodies), 0.0);
        fusedStagingBuilt_ = true;
    }

    // Fused CUDA evaluate: push X_GB*station into OpenMM's device posq, evaluate forces +
    // energy on device, and reduce forces to per-body spatial forces on device. Per robotics
    // STEP the only host<->device traffic is X_GB up (12*numBodies) and bodyForceG down
    // (6*numBodies); the atom stations (3*numAtoms) go up only when they change -- i.e. once
    // per round, on a coordinate transfer (World::recomputeGeometry -> markStationsDirty).
    void evaluateFused(RobotState& s) {
        ensureFusedStaging();
        OpenMMContext& omm = OpenMMContext::get();

        // atomStation_B is refit only on a coordinate transfer (once per round), so repack
        // the flat stations ONLY when they changed. Otherwise this per-atom pass is skipped.
        if (stationsDirty_) {
            for (int a = 0; a < model_.numAtoms; ++a) {
                const robo::Vec3& st = model_.atomStation_B[static_cast<std::size_t>(a)];
                stationFlat_[static_cast<std::size_t>(3 * a)] = st[0];
                stationFlat_[static_cast<std::size_t>(3 * a + 1)] = st[1];
                stationFlat_[static_cast<std::size_t>(3 * a + 2)] = st[2];
            }
        }
        omm.ensureKinematicsConstants(static_cast<const void*>(this),
                                      model_.numAtoms,
                                      model_.numBodies,
                                      stationFlat_.data(),
                                      model_.atomBody.data(),
                                      isVirtualI_.data(),
                                      model_.bodyAtomsBeg.data(),
                                      model_.bodyAtomsEnd.data(),
                                      model_.bodyAtoms.data(),
                                      /*stationsChanged=*/stationsDirty_);
        stationsDirty_ = false;

        // Pack X_GB as 12 doubles/body: 9 row-major rotation (robo::Mat33::elems) + 3
        // translation. Matches the kernel's row-major R*station + p.
        const robo::Transform* X_GB = s.X_GB();
        for (int b = 0; b < model_.numBodies; ++b) {
            const std::array<robo::Real, 9>& e = X_GB[b].R().elems;
            double* d = xgbFlat_.data() + static_cast<std::size_t>(12 * b);
            for (int k = 0; k < 9; ++k) {
                d[k] = e[static_cast<std::size_t>(k)];
            }
            const robo::Vec3& p = X_GB[b].p();
            d[9] = p[0];
            d[10] = p[1];
            d[11] = p[2];
        }

        omm.pushBodyTransforms(xgbFlat_.data());          // K1: posq <- X_GB*station
        lastPE_ = omm.computeForcesAndEnergyOnDevice();    // forces+energy on device
        lastEvalWasFused_ = true;                          // calcPotentialEnergy uses lastPE_
        omm.reduceForcesToBodies(bodyForceFlat_.data());   // K2: forces -> per-body

        robo::SpatialVec* BF = s.bodyForceG();
        for (int b = 0; b < model_.numBodies; ++b) {
            const double* f = bodyForceFlat_.data() + static_cast<std::size_t>(6 * b);
            BF[b] = robo::SpatialVec(robo::Vec3(f[0], f[1], f[2]),   // angular (moment about Bo)
                                     robo::Vec3(f[3], f[4], f[5]));  // linear (net force)
        }
        // Generalized (joint) forces are applied only as Cartesian atom forces, so the
        // per-DOF force is identically zero -- clear it (calcUDot reads it). Mirrors the
        // host getForcesFromOpenMM; any future Fixman/bias term adds here before calcUDot.
        robo::Real* mob = s.mobilityForce();
        for (int i = 0; i < model_.nu; ++i) {
            mob[i] = robo::Real(0);
        }
    }

    const RobotModel& model_;
    std::vector<OpenMM::Vec3> posCache_; // OpenMM particle order, nm

    // Fused CUDA robot-kinematics path state (see evaluate/evaluateFused). All inert
    // when the fused path is never taken (host path leaves lastEvalWasFused_ false).
    bool lastEvalWasFused_ = false;
    double lastPE_ = 0.0;
    bool fusedStagingBuilt_ = false;
    bool stationsDirty_ = true;          // stations refit per coordinate transfer; re-upload then
    std::vector<double> stationFlat_;    // 3*numAtoms, body-frame stations (per round)
    std::vector<int> isVirtualI_;        // numAtoms, 1 => skip (virtual site)
    std::vector<double> xgbFlat_;        // 12*numBodies, repacked + uploaded each step
    std::vector<double> bodyForceFlat_;  // 6*numBodies, downloaded each step
};