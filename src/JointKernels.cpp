#include "JointKernels.hpp"

namespace robo {

// -------- per-joint cross-mobilizer transform X_FM(q) ----------------------
// Faithful to RigidBodyNodeSpec_{Torsion,Slider,Cylinder,Translation,Ball,
// BendStretch,SphericalCoords,FreeLine,Free}.h / RigidBodyNode_Weld.
auto jointX_FM(JointType jt, const Real* q, int qOff) -> Transform {
    switch (jt) {
        case JointType::Rigid:
            return {}; // identity
        case JointType::Torsion: {
            Rotation R;
            R.setRotationFromAngleAboutZ(q[qOff]);
            return Transform(R, Vec3(0));
        }
        case JointType::Slider:
            // 1 dof: translation along F's x axis (the bond axis stays on x).
            return Transform(Rotation(), Vec3(q[qOff], 0, 0));
        case JointType::Cylinder: {
            // 2 dof: rotation about z (q0) + translation along z (q1).
            Rotation R;
            R.setRotationFromAngleAboutZ(q[qOff]);
            return Transform(R, Vec3(0, 0, q[qOff + 1]));
        }
        case JointType::BendStretch: {
            // 2 dof: rotation about z (q0) + translation along M's x (q1),
            // re-expressed in F. RigidBodyNodeSpec_Derived.cpp::RBNodeBendStretch.
            Rotation R;
            R.setRotationFromAngleAboutZ(q[qOff]);
            return Transform(R, R * Vec3(q[qOff + 1], 0, 0));
        }
        case JointType::Cartesian:
            return Transform(Rotation(), Vec3(q[qOff], q[qOff + 1], q[qOff + 2]));
        case JointType::Ball: {
            // 3 dof rotation; q = [q0 q1 q2 q3] unit quaternion, no translation.
            const Vec4 qv(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            Quaternion quat(qv); // normalizes
            Rotation R;
            R.setRotationFromQuaternion(quat);
            return Transform(R, Vec3(0));
        }
        case JointType::SphericalCoords: {
            // 3 dof BAT: body-fixed Z-Y sequence (azimuth q0 about Fz, zenith q1
            // about the new My) + radius q2 along the RADIAL axis. We implement
            // the molmodel default (BondMobility::Spherical): radial axis = Mz,
            // zero offsets, no negation. RigidBodyNodeSpec_SphericalCoords.h.
            // (OrthoSpherical / radial-X would set axisT = Mx; not represented by
            // the single SphericalCoords value -- see RobotModel note.)
            Rotation Rz;
            Rz.setRotationFromAngleAboutAxis(q[qOff], robo::ZAxis);
            Rotation Ry;
            Ry.setRotationFromAngleAboutAxis(q[qOff + 1], robo::YAxis);
            const Rotation R = Rz * Ry; // body-fixed 3-2
            const Vec3 Mz_F = R * Vec3(0, 0, 1);
            return Transform(R, Mz_F * q[qOff + 2]);
        }
        case JointType::FreeLine: {
            // 5 dof: orientation as unit quaternion q[0..3] (R_FM), translation
            // q[4..6]. The missing 6th dofs (Torsion about the body's own line, M's
            // z) lives only in the velocity Jacobian, not here.
            const Vec4 qv(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            Quaternion quat(qv);
            Rotation R;
            R.setRotationFromQuaternion(quat);
            return Transform(R, Vec3(q[qOff + 4], q[qOff + 5], q[qOff + 6]));
        }
        case JointType::Free: {
            // q = [q0 q1 q2 q3 x y z], quaternion (normalized) + translation.
            const Vec4 qv(q[qOff], q[qOff + 1], q[qOff + 2], q[qOff + 3]);
            Quaternion quat(qv); // normalizes
            Rotation R;
            R.setRotationFromQuaternion(quat);
            return Transform(R, Vec3(q[qOff + 4], q[qOff + 5], q[qOff + 6]));
        }
        default:
            return Transform();
    }
}

// -------- per-joint H_FM columns (in F), written into Hcol[0..dof-1] --------
// X_FM is needed only by the three q-dependent joints (BendStretch,
// SphericalCoords, FreeLine); the other seven ignore it (their H_FM is constant
// in F). Mirrors each node's calcAcrossJointVelocityJacobian.
void jointH_FM(JointType jt, const Transform& X_FM, SpatialVec* Hcol) {
    switch (jt) {
        case JointType::Rigid:
            break;
        case JointType::Torsion:
            Hcol[0] = SpatialVec(Vec3(0, 0, 1), Vec3(0));
            break;
        case JointType::Slider:
            Hcol[0] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            break;
        case JointType::Cylinder:
            Hcol[0] = SpatialVec(Vec3(0, 0, 1), Vec3(0)); // rotation about z
            Hcol[1] = SpatialVec(Vec3(0), Vec3(0, 0, 1)); // translation along z
            break;
        case JointType::Cartesian:
            Hcol[0] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            Hcol[1] = SpatialVec(Vec3(0), Vec3(0, 1, 0));
            Hcol[2] = SpatialVec(Vec3(0), Vec3(0, 0, 1));
            break;
        case JointType::Ball:
            Hcol[0] = SpatialVec(Vec3(1, 0, 0), Vec3(0));
            Hcol[1] = SpatialVec(Vec3(0, 1, 0), Vec3(0));
            Hcol[2] = SpatialVec(Vec3(0, 0, 1), Vec3(0));
            break;
        case JointType::BendStretch: {
            // u0: rotation about Fz; u1: translation along M's x (in F).
            const Vec3 p_FM = X_FM.p();
            const Vec3 Mx_F = X_FM.R() * Vec3(1, 0, 0);
            Hcol[0] = SpatialVec(Vec3(0, 0, 1), Vec3(0, 0, 1) % p_FM);
            Hcol[1] = SpatialVec(Vec3(0), Mx_F);
            break;
        }
        case JointType::SphericalCoords: {
            // Z-radial default: u0 about Fz, u1 about My, u2 along Mz.
            const Rotation& R = X_FM.R();
            const Vec3 p_FM = X_FM.p();
            const Vec3 sFz(0, 0, 1);
            const Vec3 sMy = R * Vec3(0, 1, 0);
            const Vec3 sMt = R * Vec3(0, 0, 1); // radial axis = Mz
            Hcol[0] = SpatialVec(sFz, sFz % p_FM);
            Hcol[1] = SpatialVec(sMy, sMy % p_FM);
            Hcol[2] = SpatialVec(Vec3(0), sMt);
            break;
        }
        case JointType::FreeLine: {
            // u0,u1: (x,y) angular velocity of M in F, expressed in M (so the
            // columns are M's x/y axes in F); u2..4: translation in F.
            const Rotation& R = X_FM.R();
            const Vec3 Mx_F = R * Vec3(1, 0, 0);
            const Vec3 My_F = R * Vec3(0, 1, 0);
            Hcol[0] = SpatialVec(Mx_F, Vec3(0));
            Hcol[1] = SpatialVec(My_F, Vec3(0));
            Hcol[2] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            Hcol[3] = SpatialVec(Vec3(0), Vec3(0, 1, 0));
            Hcol[4] = SpatialVec(Vec3(0), Vec3(0, 0, 1));
            break;
        }
        case JointType::Free:
            Hcol[0] = SpatialVec(Vec3(1, 0, 0), Vec3(0));
            Hcol[1] = SpatialVec(Vec3(0, 1, 0), Vec3(0));
            Hcol[2] = SpatialVec(Vec3(0, 0, 1), Vec3(0));
            Hcol[3] = SpatialVec(Vec3(0), Vec3(1, 0, 0));
            Hcol[4] = SpatialVec(Vec3(0), Vec3(0, 1, 0));
            Hcol[5] = SpatialVec(Vec3(0), Vec3(0, 0, 1));
            break;
        default:
            break;
    }
}

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
void jointHDot_FM(JointType jt, const Transform& X_FM, const SpatialVec& V_FM, SpatialVec* Hdot) {
    const Rotation& R = X_FM.R();
    const Vec3 p_FM = X_FM.p();
    const Vec3& w_FM = V_FM[0];
    const Vec3& v_FM = V_FM[1];
    switch (jt) {
        case JointType::BendStretch:
            Hdot[0] = SpatialVec(Vec3(0), Vec3(0, 0, 1) % v_FM);
            Hdot[1] = SpatialVec(Vec3(0), w_FM % (R * Vec3(1, 0, 0)));
            break;
        case JointType::SphericalCoords: {
            const Vec3 sFz(0, 0, 1);
            const Vec3 sMy = R * Vec3(0, 1, 0);
            const Vec3 sMt = R * Vec3(0, 0, 1);
            const Vec3 dsMy = w_FM % sMy;
            const Vec3 dsMt = w_FM % sMt;
            Hdot[0] = SpatialVec(Vec3(0), sFz % v_FM);
            Hdot[1] = SpatialVec(dsMy, (dsMy % p_FM) + (sMy % v_FM));
            Hdot[2] = SpatialVec(Vec3(0), dsMt);
            break;
        }
        case JointType::FreeLine:
            Hdot[0] = SpatialVec(w_FM % (R * Vec3(1, 0, 0)), Vec3(0));
            Hdot[1] = SpatialVec(w_FM % (R * Vec3(0, 1, 0)), Vec3(0));
            Hdot[2] = SpatialVec(Vec3(0), Vec3(0));
            Hdot[3] = SpatialVec(Vec3(0), Vec3(0));
            Hdot[4] = SpatialVec(Vec3(0), Vec3(0));
            break;
        default:
            break;
    }
}

} // namespace robo