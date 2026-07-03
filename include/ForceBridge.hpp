#pragma once

// ============================================================================
//  ForceBridge -- Robosample's DIRECT contract with OpenMM. No DuMM, no Simbody
//  subsystem, no NonBondedMappings gather.
//
//  Contract:
//    Robosample --> OpenMM : per-atom Cartesian positions in Ground (nm)
//    OpenMM     --> Robosample : per-atom Cartesian forces (kJ/mol/nm), which
//                   Robosample reduces to per-body spatial forces itself.
//  Atom order == OpenMM particle order == global/BFS order, so the transfer is
//  an index-free copy; the force->body scatter uses only model.atomBody[a].
//
//  HEADER-ONLY: every method is defined in-class (implicitly inline), so there
//  is exactly one definition across all translation units and no separate
//  ForceBridge.cpp / CMake entry is needed. (The previous declared-only methods
//  produced the `undefined symbol: ForceBridge::setAtomPositionsInGround` link
//  failure at import time.)
// ============================================================================

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "OpenMMContext.hpp" // the existing OpenMM singleton (reused, not wrapped)
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "robot_math.hpp"

class ForceBridge {
    public:
    explicit ForceBridge(const RobotModel& model)
        : model_(model) {
    }

    // Robosample -> OpenMM: cache the Ground-frame positions in OpenMM particle
    // order. The actual context->setPositions happens inside OpenMMContext on
    // the next energy/force call (both of those take the position cache).
    void setAtomPositionsInGround(const RobotState& s) {
        const robo::Vec3* p = s.atomPosG();
        posCache_.resize(static_cast<std::size_t>(model_.numAtoms));
        for (int a = 0; a < model_.numAtoms; ++a) {
            posCache_[a] = OpenMM::Vec3(p[a][0], p[a][1], p[a][2]);
        }
        lastEvalWasFused_ = false; // host positions now current; PE comes from posCache_
    }

    // Potential energy from OpenMM for the current positions [kJ/mol]. On the fused
    // CUDA path the energy was already computed on device during evaluate() (positions
    // live in posq, not posCache_), so return that; otherwise recompute from posCache_.
    [[nodiscard]] robo::Real calcPotentialEnergy() const {
        if (lastEvalWasFused_) {
            return lastPE_;
        }
        return OpenMMContext::get().computePotentialEnergy(posCache_);
    }

    // OpenMM -> Robosample: read per-atom Cartesian forces and assemble per-body
    // spatial forces about each body's origin, expressed in Ground:
    //   bodyForceG[b].linear  += f_a
    //   bodyForceG[b].angular += (r_a - origin_b) x f_a
    // This is the exact reduction the articulated-body solver consumes; it
    // replaces the whole DuMM NonBondedMappings machinery.
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
        for (int a = 0; a < model_.numAtoms; ++a) {
            // Virtual sites (mass == 0, e.g. the OPC/TIP4P M-site) carry no
            // independent DOF: OpenMM's getState(Forces) has ALREADY projected
            // the force computed at the site onto its real parent atoms (verify:
            // sum over REAL atoms reproduces the system net force / -dPE/dq to
            // machine precision, while the raw force still left in the site's own
            // slot makes the global sum nonzero). Reducing that leftover slot too
            // would DOUBLE-COUNT the site force -- here it would re-apply ~10^3
            // kJ/mol/nm per water at the M-site station, a non-conservative kick
            // that pumps the body's kinetic energy every step. Skip it; the
            // parents already carry the contribution.
            if (model_.atomMass[a] == robo::Real(0)) {
                continue;
            }
            const int b = model_.atomBody[a];
            const robo::Vec3 f(forces[a][0], forces[a][1], forces[a][2]);
            if (atomForce) {
                atomForce[a] = f; // per-atom Cartesian force for the solvent Verlet
            }
            const robo::Vec3 r = posG[a] - X_GB[b].p(); // station in Ground, about body origin
            BF[b][1] += f;                              // linear (force)
            BF[b][0] += r % f;                          // angular (moment about origin); SimTK % == cross
        }
    }

    // One full force evaluation: positions then forces (Dynamics-stage realize).
    // Fused CUDA path (opt-in, robot worlds): compute posG = X_GB*station straight
    // into OpenMM's device posq, evaluate forces+energy on device, and reduce forces
    // to per-body spatial forces on device -- so per step only O(numBodies) transforms
    // go up and O(numBodies) forces come down (no per-atom host round trip). Falls back
    // to the exact host path when the toggle is off, not a CUDA build, or when per-atom
    // forces are needed on the host (wantsAtomForces, e.g. NCMC solvent Verlet). The
    // Cartesian MD world never reaches here (it uses integrateTrajectoryOnDevice).
    void evaluate(RobotState& s) {
        if (OpenMMContext::get().cudaKinematicsAvailable() && !s.wantsAtomForces()) {
            evaluateFused(s);
        } else {
            setAtomPositionsInGround(s);
            getForcesFromOpenMM(s);
        }
    }

    // Mark the atom stations (body-frame positions) as changed so the fused path re-uploads
    // them on its next evaluate. Call after any refit of RobotModel::atomStation_B (a
    // coordinate transfer -- World::recomputeGeometry). No-op for the host path.
    void markStationsDirty() {
        stationsDirty_ = true;
    }

    void setVelocitiesToTemperature(double temperature, int seed) const {
        OpenMMContext::get().setVelocitiesToTemperature(temperature, seed);
    }

    // NCMC passthroughs (forward only ints; no SystemTopology coupling here).
    // Region A as an arbitrary atom-index set (docs/specs/ncmc-explicit-solvent/
    // 30-region-and-protocol-policy.md Sec.2).
    void enableAlchemy(const std::vector<int>& atomIndices) const {
        OpenMMContext::get().enableAlchemy(atomIndices);
    }
    // Convenience: contiguous [atomBegin,atomEnd) Region A.
    void enableAlchemy(int atomBegin, int atomEnd) const {
        OpenMMContext::get().enableAlchemy(atomBegin, atomEnd);
    }
    void setAlchemicalLambda(double lambdaInter) const {
        OpenMMContext::get().setAlchemicalLambda(lambdaInter);
    }

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