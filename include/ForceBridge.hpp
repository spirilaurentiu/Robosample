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
    }

    // Potential energy from OpenMM for the cached positions [kJ/mol].
    [[nodiscard]] robo::Real calcPotentialEnergy() const {
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
    void evaluate(RobotState& s) {
        setAtomPositionsInGround(s);
        getForcesFromOpenMM(s);
    }

    void setVelocitiesToTemperature(double temperature, int seed) const {
        OpenMMContext::get().setVelocitiesToTemperature(temperature, seed);
    }

    // NCMC passthroughs (forward only ints; no SystemTopology coupling here).
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
    const RobotModel& model_;
    std::vector<OpenMM::Vec3> posCache_; // OpenMM particle order, nm
};