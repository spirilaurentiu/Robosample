#pragma once

#include "RobotModel.hpp"
#include "robot_math.hpp"

namespace robo {

// -------- per-joint cross-mobilizer transform X_FM(q) ----------------------
// Faithful to RigidBodyNodeSpec_{Torsion,Slider,Cylinder,Translation,Ball,
// BendStretch,SphericalCoords,FreeLine,Free}.h / RigidBodyNode_Weld.
auto jointX_FM(JointType jt, const Real* q, int qOff) -> Transform;

// -------- per-joint H_FM columns (in F), written into Hcol[0..dof-1] --------
// X_FM is needed only by the three q-dependent joints (BendStretch,
// SphericalCoords, FreeLine); the other seven ignore it (their H_FM is constant
// in F). Mirrors each node's calcAcrossJointVelocityJacobian.
void jointH_FM(JointType jt, const Transform& X_FM, SpatialVec* Hcol);

// -------- per-joint dH_FM/dt columns (in F), written into Hdot[0..dof-1] -----
// Zero for the seven joints with H_FM constant in F (jointHasConstantHFM); the
// engine never calls this for them. For BendStretch / SphericalCoords / FreeLine
// these are q,u-dependent (the M axes rotate), ported verbatim from each node's
// calcAcrossJointVelocityJacobianDot. V_FM is the cross-mobilizer spatial
// velocity in F: V_FM[0] = w_FM, V_FM[1] = v_FM.
//
// VALIDATION NOTE: these three feed the mobilizer-bias (Coriolis) acceleration,
// the same hot path whose HDot/Coriolis terms have historically hidden silent
// energy-pumTorsiong bugs (see realizeVelocity). They are faithful to Simbody on
// paper but have NOT been golden-tested against a Simbody single-body reference
// in this port; do that (energy conservation on a free BendStretch / Spherical /
// FreeLine body) before trusting these joints in production.
void jointHDot_FM(JointType jt, const Transform& X_FM, const SpatialVec& V_FM, SpatialVec* Hdot);

} // namespace robo