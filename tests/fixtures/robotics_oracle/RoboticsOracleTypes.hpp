#pragma once

// ============================================================================
//  RoboticsOracleTypes.hpp -- shared, framework-free POD schema for the
//  robotics-oracle differential fixture (docs/specs/robotics-oracle-
//  differential.md, S4.1/S9).
//
//  This header has NO dependency on either engine (no SimTK/*, no robo::) so
//  it compiles unmodified on both sides of the differential:
//    - the clone-side generator (Robosample/tools/gen_robotics_oracle.cpp),
//      which fills an OracleCase from a live Simbody run and prints it as a
//      generated C++ source (S9);
//    - the disasm-side test (tests/TestRoboticsOracle.cpp), which replays the
//      SAME model spec + inputs through RobotEngine and diffs the S7 outputs
//      against the baked reference.
//
//  S4.1 correspondence: the model spec (X_PF/X_BM/mass/com/inertia) and the
//  scripted (q,u,bodyForceG,mobilityForce) inputs are stored ONCE here and
//  consumed identically by both sides -- neither side draws its own RNG for
//  the model, which is what makes the two builds comparable at all.
//
//  Field layout conventions (binding, must match S7.2 exactly):
//    * Spatial vectors are [angular; linear] (torque;force for body forces).
//    * Rotations are row-major 3x3 (R[row*3+col]).
//    * ArticulatedInertia blocks follow robo::ArticulatedInertia /
//      SimTK::ArticulatedInertia_: J (angAng, SymMat33 packed xx,xy,yy,xz,yz,zz),
//      F (angLin, full 3x3 row-major), M (linLin, SymMat33 packed).
//    * DI, Mdense are row-major nu x nu (nu <= kMaxDof).
//    * G is nu columns, each a spatial vector (ang, lin) -- G[j] is the gain
//      for mobility j (matches the port's SpatialVec G[uIndex+j] and
//      Simbody's Mat<2,dof,Vec3> column j).
// ============================================================================

namespace robotics_oracle {

inline constexpr int kSchemaVersion = 1;
inline constexpr int kMaxDof = 6; // Free is the largest Scope-A joint (dof=6)
inline constexpr int kMaxQ = 7;   // Free/FreeLine (quaternion-inflated)
inline constexpr int kMaxStates = 3;

// A single scripted (q, u, force) state and every S7 output computed from it.
struct OracleState {
    const char* label = "";

    // ---- inputs (S7.1) ----
    double q[kMaxQ] = {}; // quaternion blocks stored PRE-NORMALIZED (S6 stage 2)
    int nq = 0;
    double u[kMaxDof] = {};
    int nu = 0;
    double bodyForceTorque[3] = {}; // bodyForceG angular (torque) part, Bo, Ground
    double bodyForceForce[3] = {};  // bodyForceG linear (force) part, Bo, Ground
    double mobilityForce[kMaxDof] = {};

    // ---- stage 1: kinematics (Ground, invariant; rtol 1e-10) ----
    double X_GB_R[9] = {};
    double X_GB_p[3] = {};
    // frame-equality gate (S4.2): port X_FM must equal Simbody
    // getMobilizerTransform(state) to machine precision before ANY
    // frame-dependent (stage 3/4) comparison is trusted.
    double X_FM_R[9] = {};
    double X_FM_p[3] = {};

    // ---- stage 2: velocity (Ground, invariant; rtol 1e-9) ----
    double V_GB_ang[3] = {};
    double V_GB_lin[3] = {};
    double qdot[kMaxQ] = {};

    // ---- stage 3: articulated body inertia (frame-dependent; rtol 1e-8;
    //      gated on S6.1 min-eig(D) lock check before any element diff) ----
    double P_J[6] = {}; // SymMat33 packed: xx,xy,yy,xz,yz,zz
    double P_F[9] = {}; // Mat33 row-major
    double P_M[6] = {}; // SymMat33 packed
    double PPlus_J[6] = {};
    double PPlus_F[9] = {};
    double PPlus_M[6] = {};
    double DI[kMaxDof * kMaxDof] = {}; // row-major nu x nu
    double G_ang[kMaxDof][3] = {};
    double G_lin[kMaxDof][3] = {};
    double minEigD = 0; // S6.1 lock gate (computed from Simbody's own D)

    // ---- stage 4: acceleration pass (mixed; udot/A_GB invariant; rtol 1e-8) ----
    double Z_ang[3] = {};
    double Z_lin[3] = {};
    double ZPlus_ang[3] = {};
    double ZPlus_lin[3] = {};
    double eps[kMaxDof] = {};
    double udot[kMaxDof] = {};
    double A_GB_ang[3] = {};
    double A_GB_lin[3] = {};

    // ---- stage 5: dense M / logDetM / reactions (invariant; rtol 1e-8,
    //      1e-6 eig) ----
    double Mdense[kMaxDof * kMaxDof] = {}; // row-major nu x nu
    double logDetM = 0;
    double reactionBoAng[3] = {};
    double reactionBoLin[3] = {};
    double reactionMoAng[3] = {};
    double reactionMoLin[3] = {};
};

// The static model spec (S4.1) + its scripted states.
struct OracleCase {
    const char* name = "";

    double X_PF_R[9] = {};
    double X_PF_p[3] = {};
    double X_BM_R[9] = {};
    double X_BM_p[3] = {};
    double mass = 0;
    double com_B[3] = {};
    double unitInertia_B[6] = {}; // xx,xy,yy,xz,yz,zz

    int nq = 0;
    int nu = 0;
    int numStates = 0;
    OracleState states[kMaxStates];
};

// ============================================================================
//  Phase 1b (docs/specs/robotics-oracle-differential.md §8/§9): multi-body
//  structural cases (mixed chain, forest, wide-star hub, zero-DOF Rigid mid-
//  chain, duplicate molecules, applied-force discriminators) and the two
//  §8.1 stress cases (depth, conditioning), stored per the §9 storage rule:
//  structural cases get full per-body arrays; stress cases get ONLY
//  aggregate/extremal invariants (min-eig(D), logDetM, KE, ||udot||, and one
//  "report body"'s kinematics -- the deepest body for the depth chain, the
//  most-extreme body for the conditioning case).
// ============================================================================

inline constexpr int kMaxStructBodies = 12; // widest structural case: 8-child star hub (9 bodies)
inline constexpr int kMaxStructTotal = 24;  // widest total nq/nu: the mixed chain (nu=15, nq=17)
inline constexpr int kMaxStressBodies = 64;         // depth-stress needs >=50
inline constexpr int kMaxStressTotal = 7 * kMaxStressBodies; // generous q/u bound (7 = quaternion joint's nq)

// One body's model spec (S4.1 correspondence): parent index into the SAME
// per-case body list (0 == Ground, else 1-based), joint type (int, matching
// robo::JointType's underlying values 0..9 exactly -- RobotModel.hpp's enum
// is stable and this header cannot include it, being engine-agnostic), the
// two joint frames, and mass properties -- identical fields to OracleCase's
// single-body spec, just per-body.
struct BodyModelSpec {
    int parent = 0;
    int joint = 0; // robo::JointType underlying value
    double X_PF_R[9] = {};
    double X_PF_p[3] = {};
    double X_BM_R[9] = {};
    double X_BM_p[3] = {};
    double mass = 0;
    double com_B[3] = {};
    double unitInertia_B[6] = {};
};

// Per-body S7 outputs (same fields as OracleState's stage 1-4 block, minus
// the per-state inputs which live once per OracleMultiState, not per body).
struct OracleBodyOutput {
    double X_GB_R[9] = {};
    double X_GB_p[3] = {};
    double X_FM_R[9] = {};
    double X_FM_p[3] = {};

    double V_GB_ang[3] = {};
    double V_GB_lin[3] = {};
    double qdot[kMaxQ] = {}; // §8.2 #3: NOT compared raw for FreeLine

    double P_J[6] = {};
    double P_F[9] = {};
    double P_M[6] = {};
    double PPlus_J[6] = {};
    double PPlus_F[9] = {};
    double PPlus_M[6] = {};
    double DI[kMaxDof * kMaxDof] = {};
    double G_ang[kMaxDof][3] = {};
    double G_lin[kMaxDof][3] = {};
    double minEigD = 0;

    double Z_ang[3] = {};
    double Z_lin[3] = {};
    double ZPlus_ang[3] = {};
    double ZPlus_lin[3] = {};
    double eps[kMaxDof] = {};
    double udot[kMaxDof] = {};
    double A_GB_ang[3] = {};
    double A_GB_lin[3] = {};

    double reactionBoAng[3] = {};
    double reactionBoLin[3] = {};
    double reactionMoAng[3] = {};
    double reactionMoLin[3] = {};
};

// One scripted (q,u,force) state over the WHOLE multi-body case.
struct OracleMultiState {
    const char* label = "";

    double q[kMaxStructTotal] = {}; // system-wide, body order
    int nq = 0;
    double u[kMaxStructTotal] = {};
    int nu = 0;
    // bodyForce*[i] is body (i+1) (0-based index into the case's body list);
    // ground*[..] is the force applied to Ground (body 0) -- §8.2 #7's
    // "force on Ground must not leak into udot" discriminator.
    double bodyForceTorque[kMaxStructBodies][3] = {};
    double bodyForceForce[kMaxStructBodies][3] = {};
    double groundForceTorque[3] = {};
    double groundForceForce[3] = {};
    double mobilityForce[kMaxStructTotal] = {};

    OracleBodyOutput body[kMaxStructBodies];

    double Mdense[kMaxStructTotal * kMaxStructTotal] = {}; // row-major nu x nu, system-wide
    double logDetM = 0;
};

struct OracleMultiCase {
    const char* name = "";
    BodyModelSpec bodies[kMaxStructBodies];
    int numBodies = 0;
    int numStates = 0;
    OracleMultiState states[kMaxStates];
};

// ---- Stress cases (§8.1/§9): aggregate/extremal invariants only for the
// OUTPUTS (per §9's storage rule); the scripted INPUTS (q,u,bodyForceG) must
// still be stored in full -- both engines need them to reproduce the state.
struct OracleAggregateState {
    const char* label = "";

    double q[kMaxStressTotal] = {};
    int nq = 0;
    double u[kMaxStressTotal] = {};
    int nu = 0;
    double bodyForceTorque[kMaxStressBodies][3] = {};
    double bodyForceForce[kMaxStressBodies][3] = {};

    double logDetM = 0;
    double totalKE = 0;
    double normUdot = 0; // ||udot||_2 over the whole system
    double minEigD = 0;  // worst (smallest) min-eig(D) over all bodies

    // The "report body" (1-based index into the case's body list): the
    // deepest body for the depth-stress chain, the most mass-extreme body
    // for the conditioning-stress case.
    int reportBody = 0;
    double reportX_GB_R[9] = {};
    double reportX_GB_p[3] = {};
    double reportV_GB_ang[3] = {};
    double reportV_GB_lin[3] = {};
    double reportA_GB_ang[3] = {};
    double reportA_GB_lin[3] = {};
};

struct OracleAggregateCase {
    const char* name = "";
    BodyModelSpec bodies[kMaxStressBodies];
    int numBodies = 0;
    int numStates = 0;
    OracleAggregateState states[kMaxStates];
};

// ============================================================================
//  §8.3 randomized fuzz batch (docs/specs/robotics-oracle-differential.md
//  §8.3/§9): batch 1 draws N in [32,64] seeded random (q,u) states over each
//  of the existing Scope-A topologies; batch 2 draws a small seeded batch of
//  random topologies (depth in [2,60], random branching/joint/mass ratio).
//  Both store the SAME aggregate/extremal invariants as OracleAggregateState
//  (logDetM, KE, ||udot||, worst min-eig(D)) PLUS the frame-invariant
//  end-products udot and A_GB PER BODY (the physically meaningful
//  cross-check §8.3 asks for) -- still never full per-body DI/PPlus/G.
//
//  This is a NEW, additive struct rather than a reuse of OracleAggregateCase:
//  kMaxStates=3 is sized for the hand-authored rest/random/force battery and
//  cannot hold N in [32,64] fuzz states without either bumping a constant
//  every existing (already-green) fixture is silently exposed to, or adding
//  fields no existing case populates. A dedicated struct is schema-inert for
//  every other case kind (Rule 3).
// ============================================================================
inline constexpr int kMaxFuzzStates = 64; // §8.3: N in [32,64]

struct OracleFuzzState {
    const char* label = "";

    double q[kMaxStressTotal] = {};
    int nq = 0;
    double u[kMaxStressTotal] = {};
    int nu = 0;
    // The fuzz battery draws (q,u) only (applied-force discriminators are
    // already covered by the hand-authored ForceDiscriminators case) -- these
    // are kept (always zero) purely for on-disk schema uniformity with the
    // other case kinds, never populated with a nonzero value.
    double bodyForceTorque[kMaxStressBodies][3] = {};
    double bodyForceForce[kMaxStressBodies][3] = {};

    double logDetM = 0;
    double totalKE = 0;
    double normUdot = 0; // ||udot||_2 over the whole system
    double minEigD = 0;  // worst (smallest) min-eig(D) over all bodies, THIS engine

    // Frame-invariant end-products, PER BODY (§8.3 task 1): `udot` is the
    // system-wide vector (same body-order/qIndex-uIndex convention as `q`/
    // `u`, sliced per body by both sides via the shared BodyModelSpec/
    // RobotModel DOF table, §5 Scope-A correspondence); `A_GB` is a genuine
    // per-body spatial acceleration (defined even for a 0-dof body).
    double udot[kMaxStressTotal] = {};
    double A_GB_ang[kMaxStressBodies][3] = {};
    double A_GB_lin[kMaxStressBodies][3] = {};
};

struct OracleFuzzCase {
    const char* name = "";
    BodyModelSpec bodies[kMaxStressBodies];
    int numBodies = 0;
    int numStates = 0;
    OracleFuzzState states[kMaxFuzzStates];
    // §8.3/§9: "record the seed... in the fixture metadata" / "count and
    // record how many were resampled in the manifest". -1 is the "not
    // present" sentinel the loader uses for a non-fuzz-adjacent read.
    long long rngSeed = -1;
    int resampleCount = -1;
};

} // namespace robotics_oracle
