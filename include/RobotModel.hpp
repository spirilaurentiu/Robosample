#pragma once

// ============================================================================
//  RobotModel — the immutable, precomputed, per-World graph. Built ONCE in
//  World::build() from the SoA SystemTopology + this world's per-bond
//  mobilities. Never rebuilt on a coordinate transfer (this is what kills the
//  per-transfer realizeTopology() cost in the current code).
//
//  This is a faithful SoA port of what World currently precomputes across
//  decomposeRigidUnits / buildFrameGraph / modelOneCompound / modelTopologies,
//  with the joint frames and mass properties stored as CONSTANTS rather than
//  Simbody topology-stage defaults.
//
//  SoA, no vtables. Joint identity is an enum per body, branched in the engine
//  (replacing Simbody's RigidBodyNode virtual dispatch).
//
//  Tables here are built once and never resized; std::vector is fine. They can
//  later be frozen into a MemoryArena alongside RobotState if profiling wants
//  the sweep inputs cache-packed — the engine reads them through pointers only.
// ============================================================================

#include <cstdint>
#include <vector>

#include "robot_math.hpp"

// Joint type per body. Mirrors the mobility->MobilizedBody mapping in
// World::modelOneCompound. NO SimTK enum on the hot path.
enum class JointType : std::uint8_t {
    Weld = 0,    // 0 dof
    Pin,         // 1 dof  (Torsion bond mobility)
    Slider,      // 1 dof
    Cylinder,    // 2 dof
    BendStretch, // 2 dof
    Translation, // 3 dof  (Cartesian)
    Ball,        // 3 dof, q = 4 (quaternion) by default
    SphericalCoords,
    FreeLine, // 5 dof
    Free,     // 6 dof, q = 7 (quaternion + translation)
};

// How a molecule's root body attaches to Ground (mirrors SimTK::RootMobility).
// RootMobility is the SAME enum the topology + Python bindings already use
// (SimTK::RootMobility, defined in molmodel CompoundSystem.h with
// Free=0,Cartesian=1,Weld=2,FreeLine=3,Ball=4,Pin=5). We forward-declare it
// (opaque scoped enum with fixed underlying type == a complete type) and alias
// it so the engine shares ONE enum with the bound topology -- no parallel
// values, and no Molmodel header pulled into this foundational file. The full
// definition arrives via TopologyElements.hpp in the .cpp that need enumerators.
enum class RootMobility : std::uint8_t {
    Free,
    Cartesian,
    Weld,
    FreeLine,
    Ball,
    Pin
};

// ----------------------------------------------------------------------------
//  INDEX / ORDER CONVENTIONS (must match python/robosample/molecule_prototype.py
//  and z_matrix.py exactly -- see MIGRATION_PLAN_v2_corrections.md Part B 5/6):
//   * Atoms arrive in BFS/compound order from the root (compound 0 == root).
//     The array position IS the atom index; Molmodel cAIx is retired.
//   * Bonds are oriented PARENT->CHILD (root-ward first):  bondsI=parent,
//     bondsJ=child. Ring-closing bonds are flagged (transmitted explicitly)
//     and are NEVER tree edges (skipped in decomposeRigidUnits/buildFrameGraph).
//   * Angles oriented i->k (j center); proper torsions oriented i->l; all
//     root-ward-first, matching Python.
//   * Z-MATRIX uses the standard BAT convention and is the OPPOSITE direction:
//     zI = atom PLACED; zJ/zK/zL = back-references TOWARD root
//     (parent / grandparent / great-grandparent). BAT = dist(zI,zJ),
//     angle(zI,zJ,zK), dihedral(zI,zJ,zK,zL). -1 sentinels in leading rows.
//   * Frame-graph ancestry names: self / parent / grandParent (root-ward).
//     (refChild only exists if the legacy A1 orientation switch is kept.)
// ----------------------------------------------------------------------------
struct RobotModel {
    // ---- system sizes ------------------------------------------------------
    int numBodies = 0; // includes Ground == body 0
    int numAtoms = 0;
    int nq = 0; // total generalized coordinates (quaternion-inflated)
    int nu = 0; // total generalized speeds (== DOFs)

    // ---- body tree, TOPOLOGICAL ORDER (parent index < child index) ---------
    // Sweeps: inward = reverse iteration, outward = forward iteration.
    std::vector<int> bodyParent;      // [numBodies]   parent body (Ground's parent = -1)
    std::vector<int> bodyLevel;       // [numBodies]   depth from Ground (for future level-parallel sweeps)
    std::vector<JointType> bodyJoint; // [numBodies]
    std::vector<int> bodyRootAtom;    // [numBodies]   the inboard ("root") atom of this body
    std::vector<int> bodyChildrenBeg; // [numBodies]   CSR into bodyChildren
    std::vector<int> bodyChildrenEnd; // [numBodies]
    std::vector<int> bodyChildren;    // flattened children, CSR

    // ---- generalized-coordinate offsets ------------------------------------
    std::vector<int> bodyQIndex;   // [numBodies] first q of this body
    std::vector<int> bodyNQ;       // [numBodies] number of q (4 or 7 if quaternion)
    std::vector<int> bodyUIndex;   // [numBodies] first u of this body
    std::vector<int> bodyNU;       // [numBodies] number of u == dof
    std::vector<int> bodyUSqIndex; // [numBodies] offset into the dof*dof DI pool
    int nuSq = 0;                  // sum of dof^2 over bodies (DI pool size)
    // q-slots that are quaternions and must be renormalized each Verlet step
    // (port of the no-constraint q-projection). Stored as (qStart) of each quat.
    std::vector<int> quaternionQStart; // q index where a 4-wide quaternion begins

    // ---- static joint frames (CONSTANT once built) -------------------------
    // X_PF: inboard frame fixed on the PARENT body. X_BM: outboard frame fixed
    // on THIS body. These define the mobilizer; in the current code they were
    // Simbody topology defaults rewritten every transfer — here they are fixed.
    std::vector<robo::Transform> X_PF; // [numBodies]
    std::vector<robo::Transform> X_BM; // [numBodies]

    // ---- mass properties (RECOMPUTED every transfer; about body origin)
    std::vector<robo::Real> bodyMass;  // [numBodies]
    std::vector<robo::Vec3> bodyCom_B; // [numBodies] center of mass in body frame
    std::vector<robo::UnitInertia>
        bodyUnitInertia_B; // [numBodies] unit inertia about body origin, body frame

    // ---- atom <-> body maps ------------------------------------------------
    std::vector<int> atomBody;             // [numAtoms]  body index of each atom
    std::vector<robo::Real> atomMass;      // [numAtoms]  per-atom mass (Daltons), constant
    std::vector<robo::Vec3> atomStation_B; // [numAtoms]  atom in body frame (per transfer)
    std::vector<int> atomDummIndex;        // [numAtoms]  DuMM atom index (for OpenMM mapping)
    std::vector<int> atomCompoundIndex;    // [numAtoms]  cAIx (analysis / IO compatibility)
    std::vector<int> bodyAtomsBeg;         // [numBodies] CSR into a body-sorted atom list
    std::vector<int> bodyAtomsEnd;         // [numBodies]
    std::vector<int> bodyAtoms;            // flattened atoms per body, CSR

    // ---- frame graph (your existing Cartesian->frame kernel, kept verbatim) -
    struct FrameGraph {
        int totalAtoms = 0;
        std::vector<int> topoOffset;                              // [numTopologies+1]
        std::vector<int> g_self, g_parent, g_gparent, g_refChild; // full-geometry bucket
        std::vector<int> f_self, f_parent;                        // fallback bucket
        std::vector<int> r_self;                                  // roots
    } frameGraph;

    // ---- OpenMM force->body reduction (NO NonBondedMappings) ----------------
    // DuMM's NonBondedMappings (dummAtomIndex/includedAtomIndex/bodyIndex/
    // bodyStart) is DROPPED. Atom index == OpenMM particle index == SoA
    // position, so positions hand over with no gather. The force->body
    // reduction needs ONLY atomBody[a] (above): scatter into the small
    // bodyForce array (numBodies SpatialVecs, hot in L1). The body-sorted
    // atom CSR (bodyAtoms/bodyAtomsBeg/End, above) is OPTIONAL and only used
    // if a very large single robot wants a strict sequential reduction.

    // ---- z-matrix rows (BAT) ----------------------------------------------
    std::vector<int> zI, zJ, zK, zL; // [numZRows]; -1 sentinels in leading rows
    int numZRows = 0;

    // ---- convenience -------------------------------------------------------
    [[nodiscard]] bool isQuaternionBody(int b) const {
        return bodyJoint[b] == JointType::Ball || bodyJoint[b] == JointType::Free;
    }
};