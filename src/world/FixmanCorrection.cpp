// ============================================================================
//  FixmanCorrection.cpp - World::calcFixman / World::calcLogSineSqrGamma2.
//
//  Relocated verbatim from World.cpp (SPLIT-W5, pure code motion): the
//  Fixman / coordinate-Jacobian corrections for torsional worlds -- the
//  compensating potential U_F = 1/2 RT (ln|M_tree| - ln det(G M^-1 G^T) -
//  ln|M_3N|) (Spiridon & Minh 2017) and the external-rotation Jacobian
//  -1/2 RT Sigma_b ln sin^2(gamma2_b) (default OFF; wrong for the engine's
//  unit-quaternion roots). Pure functions of (RobotModel, RobotState,
//  ConstraintSet) plus the constant lnDetMCartesian_.
// ============================================================================

#include "World.hpp"

#include <algorithm>
#include <cmath>

#include "engine_helpers.hpp"

using robo::Real;
using robo::Rotation;
using robo::Vec3;

namespace {

// ---------------------------------------------------------------------------
//  Absolute orientation of a free-root body, built from a body-fixed atom
//  triplet in the CURRENT Ground frame.
//
//  WHY NOT X_GB[b].R(): the per-block frame handoff (setAtomsLocationsInGround
//  -> recomputeGeometry) rebuilds every root frame with IDENTITY rotation and
//  resets the joint coordinate q = 0 (computeAllAtomFrames r_self branch). So
//  X_GB[b].R() carries only the WITHIN-block deviation from that reset -- it is
//  exactly identity at the START of every move (the gimbal pole, sin(pitch)->0).
//  Reading it floors J(q) to its eps cap for every free root at H_old, every
//  round, independent of the real pose -- the +41352 vs +7352 artifact. The
//  external-rotation Jacobian needs the molecule's TRUE orientation in space
//  (Section 6: "body b's orientation"), which is frame-reset-invariant and read
//  here straight from the atom geometry.
//
//  Convention: x along (root -> reference atom), z along the triplet normal,
//  y = z x x. Any fixed convention is admissible: it shifts J by a per-body
//  constant that is IDENTICAL at H_old and H_new, so only the physical change in
//  orientation across the trajectory survives in dH. Returns false for a body
//  without 3 non-collinear real atoms (no 3D orientation DOF -> no J term).
bool freeRootAbsRotation(const RobotModel& m, const RobotState& s, int b, Rotation& Rout) {
    const Vec3* P = s.atomPosG();
    const int a0 = m.bodyRootAtom[b];
    if (a0 < 0) {
        return false;
    }
    const Vec3 p0 = P[a0];
    // Pick the two real (mass>0) body atoms, distinct from the root, that span
    // the largest triangle -- the most numerically stable, orientation-defining
    // pair (a near-collinear pick would make the frame ill-conditioned).
    int bestA1 = -1, bestA2 = -1;
    Real bestArea = 0;
    for (int ci = m.bodyAtomsBeg[b]; ci < m.bodyAtomsEnd[b]; ++ci) {
        const int a1 = m.bodyAtoms[ci];
        if (a1 == a0 || m.atomMass[a1] <= Real(0)) {
            continue;
        }
        for (int cj = ci + 1; cj < m.bodyAtomsEnd[b]; ++cj) {
            const int a2 = m.bodyAtoms[cj];
            if (a2 == a0 || m.atomMass[a2] <= Real(0)) {
                continue;
            }
            const Real area = ((P[a1] - p0) % (P[a2] - p0)).norm();
            if (area > bestArea) {
                bestArea = area;
                bestA1 = a1;
                bestA2 = a2;
            }
        }
    }
    if (bestA1 < 0 || bestArea < Real(1e-10)) {
        return false;
    }
    const Vec3 ex = (P[bestA1] - p0) / (P[bestA1] - p0).norm();
    Vec3 ez = ex % (P[bestA2] - p0);
    ez = ez / ez.norm();
    const Vec3 ey = ez % ex;
    Rout = Rotation(robo::Mat33(ex[0], ey[0], ez[0], ex[1], ey[1], ez[1], ex[2], ey[2], ez[2]));
    return true;
}

} // namespace

// ----------------------------------------------------------------------------
//  Fixman / coordinate-Jacobian corrections (torsional worlds)
// ----------------------------------------------------------------------------
double World::calcFixman() {
    // Fixman compensating potential, Spiridon & Minh 2017 (JCTC 13:4649) Eq. 3:
    //   U_F = (1/2) RT ln( |M_{N_f}| / |M_3N| )
    // where |M_{N_f}| is the mass-metric determinant of the constrained system's
    // FLEXIBLE coordinates and |M_3N| (constant) is the Cartesian reference. U_F
    // makes the constrained-move marginal (their Eq. 2, rho ~ |M_{N_f}|^{1/2}
    // e^{-bU}) match the unconstrained Cartesian Boltzmann marginal.
    RobotEngine::realizeArticulatedBodyInertias(model_, state_);

    // (a) Tree term: ln|M_tree| = sum_b ln det(D_b), the O(n) articulated-body
    //     determinant (Jain et al., refs 9 & 23). For an ACYCLIC molecule the
    //     flexible coordinates ARE the tree torsions, so |M_{N_f}| = |M_tree| and
    //     this term alone is exact -- the regime the paper validated.
    const double lnDetM = RobotEngine::calcLogDetM(model_, state_);

    // (b) Loop-closure term: for a CYCLIC molecule the ring is opened into the
    //     spanning tree and closed by a RATTLE distance constraint, which removes
    //     one flexible DOF per ring. The RATTLE momentum projection (G M^-1 p = 0)
    //     contributes det(G M^-1 G^T)^{-1/2} to the marginal, so the correct
    //     flexible determinant is |M_{N_f}| = |M_tree| / det(G M^-1 G^T). The Jain
    //     tree algorithm does not include this (it is a branched-molecule method);
    //     the paper never tested ring Boltzmann correctness (its macrocycle data,
    //     Sec. 3.5, measured efficiency only), so the term was simply missing for
    //     cyclic systems. calcConstraintLogDet returns 0 when there are no loop
    //     closures, making this a guaranteed no-op on every acyclic system.
    const double lnDetZ = constraints_.calcConstraintLogDet(model_, state_);

    // U_F = 1/2 RT ( ln|M_tree| - ln det(G M^-1 G^T) - ln|M_3N| ).
    // = 1/2 RT ( ln|M_{N_f}| - ln|M_3N| ), i.e. Eq. 3 with the cyclic |M_{N_f}|.
    return 0.5 * RT_ * (lnDetM - lnDetZ - lnDetMCartesian_);
}

double World::calcLogSineSqrGamma2() const {
    // External-rotation Jacobian J(q) = -(1/2) RT sum_b ln sin^2(gamma2_b), the
    // sin(theta) polar volume factor of an EULER/SPHERICAL parameterization of
    // orientation. This is gated by sampler_.useOrientationJacobian, DEFAULT OFF,
    // because Robosample's free/ball roots are parameterized by UNIT QUATERNIONS,
    // not Euler angles, and for a quaternion root the term does NOT belong:
    //
    //   The flat (uniform) measure on the unit 3-sphere S^3 -- i.e. drawing the
    //   quaternion uniformly subject to ||q|| = 1 -- pushes forward to EXACTLY
    //   the Haar (rotation-invariant) measure on SO(3). "Flat in quaternion
    //   space" and "Haar-uniform on rotations" are the same distribution; there
    //   is no Jacobian between them. The quaternion exp-map integrator
    //   (RobotIntegrator advanceQuatExp) is the geodesic flow of the
    //   left-invariant kinetic metric, and for a single free body with constant
    //   body inertia det M(q) is orientation-independent, so GC-HMC with U = 0
    //   already samples orientation Haar-uniformly with NO correction term. The
    //   sin^2(gamma2) factor is the correct Jacobian only if one sampled in Euler
    //   angles; applied on top of the quaternion parameterization it biases the
    //   marginal toward the poles. See tests/TestEnsembleOrientation.cpp, which
    //   demonstrates: OFF -> Haar-uniform; ON -> measurably biased.
    //
    // The term is retained behind the flag for any future Euler-parameterized
    // joint where it WOULD be correct. gamma2_b is the pitch of the body's TRUE
    // orientation in space (freeRootAbsRotation), NOT X_GB[b].R() -- which the
    // per-block reset Torsions to identity (q = 0) at the start of every move,
    // flooring J for all roots (see helper).
    double acc = 0.0;
    for (int b = 1; b < model_.numBodies; ++b) {
        if (model_.bodyParent[b] != 0 || model_.bodyJoint[b] != JointType::Free) {
            continue;
        }
        Rotation Rabs;
        if (!freeRootAbsRotation(model_, state_, b, Rabs)) {
            continue; // < 3 non-collinear real atoms: no 3D orientation DOF
        }
        Real w, x, y, z;
        EngineHelpers::rotationToQuaternion(Rabs, w, x, y, z);
        const Real sinPitch = std::clamp(Real(2.0) * ((w * y) - (z * x)), Real(-1.0), Real(1.0));
        acc += EngineHelpers::safeLogSineSqr(std::asin(sinPitch));
    }
    return acc;
}
