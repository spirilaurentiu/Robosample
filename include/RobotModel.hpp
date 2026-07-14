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

/**
 * @brief Joint (mobilizer) identity for one body: the single per-body tag the
 *        engine branches on in place of virtual dispatch.
 * @note The DOF count, q width, quaternion use, constant-H_FM property, and
 *       root-legality of each value are given by the static predicates on
 *       RobotModel; those predicates are the single source of truth for
 *       joint-type facts.
 */
enum class JointType : std::uint8_t {
    Rigid = 0,       // 0 dof
    Torsion,         // 1 dof  rotation about the bond axis (the canonical dihedral)
    Slider,          // 1 dof  translation along the bond axis
    Cylinder,        // 2 dof  rotation + translation about/along the bond axis
    BendStretch,     // 2 dof  rotation perpendicular to the bond + translation along it
    Cartesian,       // 3 dof  3 translations
    Ball,            // 3 dof  rotation, q = 4 (quaternion) by default
    SphericalCoords, // 3 dof  BAT (azimuth, zenith, radius)
    FreeLine,        // 5 dof  2 rotations (no spin about own line) + 3 translations, q = 7
    Free,            // 6 dof  q = 7 (quaternion + translation)
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
/**
 * @brief Immutable, precomputed articulated-robot topology for one World: the
 *        body forest, generalized-coordinate index tables, static joint frames,
 *        mass properties, atom<->body maps, loop-closure (z-matrix) records, and
 *        the joint-fact predicates the solver queries.
 *
 * @par Ownership and immutability boundary
 * Owned by @c World by value and built once by ModelBuilder (World::buildModel),
 * the sole writer. After the build it is @b const for the life of the World and
 * every solver function takes it as @c const @c RobotModel&; it is rebuilt only
 * on a root-mobility change. Two field classes are exceptions to "set once at
 * build": bodyMass / bodyCom_B / bodyUnitInertia_B / atomStation_B are
 * recomputed on every coordinate transfer (marked in their group comments), and
 * bodyMassScale is a sampling-time run-constant. The topology and index tables
 * (sizes, tree, q/u offsets, joint frames, z-matrix) are fixed after the build.
 *
 * @par Index and ordering conventions
 * Atoms and bodies are in topological/BFS order (parent index < child index,
 * Ground == body 0); the array position is the index. The q/u offset tables
 * (bodyQIndex/bodyNQ/bodyUIndex/bodyNU/bodyUSqIndex) define the
 * generalized-coordinate layout the mass metric (INV-5) and constraints (INV-6)
 * read; @c nq exceeds @c nu only for quaternion bodies. See the field-group
 * comments for each array's index space and unit; ModelBuilder is the producer.
 *
 * @note SoA, no vtables: joint behavior is dispatched on bodyJoint[b] via the
 *       static predicates below, not through polymorphism.
 */
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

    // ---- kinetic-metric preconditioning (fictitious mass; SAMPLING only) ----
    // Per-body scale applied to the body spatial inertia Mk_G in the KINETIC
    // metric ONLY: the momentum draw (sqrt(M^-1)), the kinetic energy
    // (1/2 u^T M u), and the Fixman log-det (ln det M). Default 1.0 == physical.
    // MUST be a run-constant (never a function of q). Does NOT enter V.
    //
    // Why this is bias-free: Mk_G is the single source for P/D (hence both the
    // sqrt(M^-1) draw and ln det M) and for calcKineticEnergy, so scaling it by a
    // constant s_b shifts ln det M by the constant sum_b dof_b * ln s_b, which
    // cancels in dH; and the marginalised det(M)^{1/2} is cancelled exactly by
    // the Fixman exp(-beta F) for the SAME (scaled) M. Net effect: the sampled
    // configurational distribution is unchanged; only the proposal dynamics are
    // rescaled, raising the stable dt ~sqrt(s_b) for the scaled body. Used to
    // tame stiff explicit-solvent libration. See World::setMassScaleByJoint /
    // setBodyMassScale. May be left empty (treated as all-1.0).
    std::vector<robo::Real> bodyMassScale; // [numBodies]  default 1.0

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
    // One row per FLEXIBLE (non-Rigid) body: zI = this body's root ("placed")
    // atom; zJ/zK/zL = parent/grandparent/great-grandparent BODY's root atom
    // (root-ward, standard Z-matrix convention -- see the file-header index
    // convention block above). -1 sentinels where the ancestor chain reaches
    // Ground before 3/4 atoms are collected (root-adjacent bodies). Populated
    // by World::buildModel (docs/specs/replica-exchange-nonequilibrium-work.md
    // D6/D5 -- the BAT-scaling Jacobian reads r=dist(zI,zJ), theta=
    // angle(zI,zJ,zK) from these rows). A scaled joint (Slider/Cylinder/
    // BendStretch/SphericalCoords) can never be a molecule root
    // (jointIsLegalRoot excludes them), so those bodies' zJ is always valid.
    std::vector<int> zI, zJ, zK, zL; // [numZRows]; -1 sentinels in leading rows
    int numZRows = 0;
    std::vector<int> bodyZRow; // [numBodies] index into zI/zJ/zK/zL, or -1 (Rigid bodies)

    // ---- convenience -------------------------------------------------------
    /**
     * @brief Whether body @p b stores its orientation as a 4-wide unit quaternion
     *        (and so must be renormalized each Verlet step).
     * @param[in] b body index in [0, numBodies).
     */
    [[nodiscard]] bool isQuaternionBody(int b) const {
        return jointUsesQuaternion(bodyJoint[b]);
    }

    // ---- joint-type facts: the SINGLE source of truth -----------------------
    // Every place that needs per-joint sizes/flags (the builder's dof counting,
    // the engine's q/u layout, the quaternion renormalizer) reads these, so a
    // new joint is defined in exactly one spot. nq differs from nu only for the
    // quaternion bodies (orientation stored as a 4-wide unit quaternion).
    /**
     * @brief Number of generalized speeds (DOF) a joint of type @p jt contributes.
     * @return 0 (Rigid) up to 6 (Free); the value used to size bodyNU and lay out u.
     */
    [[nodiscard]] static constexpr int jointNU(JointType jt) {
        switch (jt) {
            case JointType::Rigid:
                return 0;
            case JointType::Torsion:
            case JointType::Slider:
                return 1;
            case JointType::Cylinder:
            case JointType::BendStretch:
                return 2;
            case JointType::Cartesian:
            case JointType::Ball:
            case JointType::SphericalCoords:
                return 3;
            case JointType::FreeLine:
                return 5;
            case JointType::Free:
                return 6;
        }
        return 0;
    }
    /**
     * @brief Number of generalized coordinates (q width) for a joint of type @p jt.
     * @return jointNU(jt) for non-quaternion joints; the quaternion-inflated width
     *         for Ball (4), FreeLine (7), and Free (7).
     */
    [[nodiscard]] static constexpr int jointNQ(JointType jt) {
        // Quaternion bodies inflate the orientation block from 3 (rotational u)
        // to 4 (unit quaternion q): Ball (4), FreeLine (4 + 3 trans = 7), Free
        // (4 + 3 trans = 7). All others have nq == nu.
        switch (jt) {
            case JointType::Ball:
                return 4;
            case JointType::FreeLine:
            case JointType::Free:
                return 7;
            default:
                return jointNU(jt);
        }
    }
    /**
     * @brief Whether type @p jt stores orientation as a unit quaternion (Ball,
     *        FreeLine, Free), making @c nq > @c nu for that body.
     */
    [[nodiscard]] static constexpr bool jointUsesQuaternion(JointType jt) {
        return jt == JointType::Ball || jt == JointType::FreeLine || jt == JointType::Free;
    }
    /**
     * @brief Whether type @p jt contributes any DOF; false only for Rigid, which
     *        welds its two atoms into one rigid unit.
     */
    [[nodiscard]] static constexpr bool jointIsFlexible(JointType jt) {
        return jt != JointType::Rigid; // Weld/Rigid welds its two atoms into one rigid unit
    }
    // H_FM is CONSTANT in the F (inboard) frame for these, so HDot_FM == 0 and
    // the engine's mobilizer-bias acceleration uses the cheap centripetal path.
    // The other three (BendStretch, SphericalCoords, FreeLine) have q-dependent
    // H_FM and go through the general HDot_FM*u term.
    /**
     * @brief Whether the joint matrix H_FM is constant in the inboard F frame, so
     *        HDot_FM == 0.
     * @return true for the seven constant-H_FM joints (engine uses the cheap
     *         centripetal mobilizer-bias path); false for BendStretch,
     *         SphericalCoords, FreeLine (q-dependent H_FM, general HDot_FM*u term).
     * @note Consumed in RobotEngine::realizeVelocity to select the bias path.
     */
    [[nodiscard]] static constexpr bool jointHasConstantHFM(JointType jt) {
        switch (jt) {
            case JointType::Rigid:
            case JointType::Torsion:
            case JointType::Slider:
            case JointType::Cylinder:
            case JointType::Cartesian:
            case JointType::Ball:
            case JointType::Free:
                return true;
            default: // BendStretch, SphericalCoords, FreeLine
                return false;
        }
    }
    // Only these joints are meaningful as a molecule-root attachment to Ground
    // (the former RootMobility set). Slider/Cylinder/BendStretch/SphericalCoords
    // are defined against a BOND axis that does not exist at the Ground hinge, so
    // they are rejected there (validated in World::buildModel). This is a rule on
    // ONE enum, not a second enum.
    /**
     * @brief Whether type @p jt may attach a molecule root to Ground.
     * @return true for Free, Cartesian, Rigid, FreeLine, Ball, Torsion; false for
     *         the bond-axis joints (Slider/Cylinder/BendStretch/SphericalCoords),
     *         which have no bond axis at the Ground hinge.
     * @note Enforced in World::buildModel; a scaled joint can never be a root.
     */
    [[nodiscard]] static constexpr bool jointIsLegalRoot(JointType jt) {
        switch (jt) {
            case JointType::Free:
            case JointType::Cartesian: // == Cartesian
            case JointType::Rigid:     // == Rigid
            case JointType::FreeLine:
            case JointType::Ball:
            case JointType::Torsion:
                return true;
            default:
                return false;
        }
    }
};