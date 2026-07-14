#include "bridge/ForceReducer.hpp"

void reduceAtomForcesToBodies(const robo::Vec3* atomForceG,
                              const robo::Vec3* atomPosG,
                              const robo::Real* atomMass,
                              const int* atomBody,
                              const robo::Transform* X_GB,
                              int numAtoms,
                              int numBodies,
                              robo::SpatialVec* bodyForceG) {
    (void) numBodies; // bodyForceG is caller-sized and caller-zeroed
    for (int a = 0; a < numAtoms; ++a) {
        // Virtual sites (mass == 0, e.g. the OPC/TIP4P M-site) carry no
        // independent DOF: OpenMM's getState(Forces) has ALREADY projected
        // the force computed at the site onto its real parent atoms (verify:
        // sum over REAL atoms reproduces the system net force / -dPE/dq to
        // machine precision, while the raw force still left in the site's own
        // slot makes the global sum nonzero). Reducing that leftover slot too
        // would DOUBLE-COUNT the M-site force -- here it would re-apply ~10^3
        // kJ/mol/nm per water at the M-site station, a non-conservative kick
        // that pumps the body's kinetic energy every step. Skip it; the
        // parents already carry the contribution.
        if (atomMass[a] == robo::Real(0)) {
            continue;
        }
        const int b = atomBody[a];
        const robo::Vec3& f = atomForceG[a];
        const robo::Vec3 r = atomPosG[a] - X_GB[b].p(); // station in Ground, about body origin
        bodyForceG[b][1] += f;                          // linear (force)
        bodyForceG[b][0] += r % f;                      // angular (moment about origin); SimTK % == cross
    }
}
