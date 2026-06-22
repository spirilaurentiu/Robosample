#pragma once

// ============================================================================
//  World - reimplemented, SimTK-free. Each World owns ONE robot (immutable
//  RobotModel + live RobotState), a RobotEngine (stateless ops), and a
//  ForceBridge to OpenMM.
//
//  Regimes (selected by moveType / cartesian_):
//    Cartesian world  : plain OpenMM MD on the device (integrateTrajectoryOnDevice).
//    Torsional/robotic: internal-coordinate Generalized-Coordinate HMC. Sampling
//                       = seed u via sqrt(M^-1), Verlet stepTo, Metropolis on the
//                       FULL Hamiltonian H = PE + KE + Fixman - 1/2 RT logSineSqr.
//    Docking world    : a torsional world whose ligand molecule roots are Free
//                       and everything else is Welded/Rigid. Its move is a rigid
//                       KICK of the ligand (a symmetric Cartesian proposal,
//                       Metropolised on dU), optionally followed by the standard
//                       MD-HMC relaxation sub-move. Root mobility here is a
//                       *World* property, supplied by Context::addDockingWorld.
// ============================================================================

#include <cstdint>
#include <optional>
#include <random>
#include <vector>

#include "Constraints.hpp"
#include "ForceBridge.hpp"
#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TopologyElements.hpp"
#include "robot_math.hpp"

enum class AcceptRejectMode : std::uint8_t {
    AlwaysAccept = 0,
    MetropolisHastings
};

// How a world proposes a sample.
//   MdHmc     : draw momenta, run constrained Verlet, Metropolis on full H.
//   RigidKick : rigid-body kick of the ligand (docking). A symmetric Cartesian
//               proposal Metropolised on dU; if mdSteps > 0 it is followed by an
//               MdHmc relaxation sub-move on the same world.
enum class MoveType : std::uint8_t {
    MdHmc = 0,
    RigidKick
};

enum class BondMobility : std::uint8_t {
    Rigid = 0,
    Torsion,
    Free,
    Ball,
    Pin,
    Slider,
    Cylinder,
    BendStretch
};

// Per-bond mobility selection produced by Context::build_flexibilities.
struct Selection {
    std::vector<BondMobility> bondMobility; // [numBonds]; default Rigid
};

struct SamplerConfig {
    double timeStep = 0.001;
    int mdSteps = 0;
    AcceptRejectMode acceptRejectMode = AcceptRejectMode::AlwaysAccept;
    bool useFixman = false; // include the Fixman potential + logSineSqr in H
    MoveType moveType = MoveType::MdHmc;
    bool useNuts = false;

    // Docking: the binding sphere is sized automatically, per ligand, as
    //   R_i = R_receptor + sphereFactor * R_ligand_i
    // (R_receptor / R_ligand_i = max extent of each set from its centroid, from
    // current geometry). sphereFactor scales only the ligand allowance: the guest
    // may roam roughly sphereFactor ligand-radii beyond the receptor surface.
    double sphereFactor = 1.0;

    // alwaysKick == false (default): CONTAINMENT -- the kick fires for a ligand
    // only when its COM has left the sphere; otherwise the move is pure dynamics.
    // alwaysKick == true: perturb every ligand every round regardless. Either way
    // the perturbation relocates the guest (uniform position in the sphere +
    // uniform reorientation = a full perturbation of its external q) and the
    // dynamics step then reseeds qdot/qdotdot.
    bool alwaysKick = false;

    // A proposed conformation is VALID only if its potential energy is finite and
    // |peNew| <= clashThreshold * max(1, |pePre|). In other words, reject if the
    // proposed PE is more than clashThreshold times larger (in absolute value) than
    // the pre-kick PE -- an order-of-magnitude relative check. This handles both
    // free-space (~0 kJ/mol) and bound (~-2000 kJ/mol) regimes without a hard-coded
    // absolute ceiling. Default 10 = one order of magnitude. Enforced in EVERY mode
    // (AlwaysAccept means skip the Metropolis probability, not allow broken geometry).
    double clashThreshold = 10.0; // dimensionless factor

    // ABSOLUTE physical ceiling on a pose's potential energy [kJ/mol], shared by
    // the acceptance gate and the forced-kick rescue. The relative clashThreshold
    // alone is dangerous in the bound regime: at pePre ~ -2400 it admits poses up
    // to ~24000, so a clash at e.g. +9800 passes validity AND sits below the old
    // hard-coded 1e4 rescue ceiling -- a strained, non-integrable pose that becomes
    // an absorbing state (PE frozen, q restored every round). Tying both gates to
    // ONE absolute ceiling removes that dead band: nothing the acceptance gate lets
    // in can exceed what the rescue treats as bad. A pose this far above the bound
    // minimum is physically a clash, never a resting state.
    double maxStartPE = 5.0e3; // kJ/mol; entry ceiling == rescue ceiling

    // INITIAL KICK: before round 0, keep re-drawing the ligand position until
    // the proposal passes the clash gate (dPE <= maxStartPE). This guarantees
    // the simulation starts from a clash-free pose without burning real rounds.
    // maxInitialKickTries caps the retry loop; if no clean pose is found within
    // that budget the run throws rather than starting from a clashing geometry.
    // Set to 0 (default) to disable -- the first round handles placement as usual.
    int maxInitialKickTries = 0; // 0 = disabled; >0 = retry budget

    // ESCAPE HATCH (docking): max consecutive rejected docking moves tolerated from
    // the SAME carried-forward pose before a kick is FORCED regardless of energy
    // magnitude. A pose that is finite and below maxStartPE can still be locally
    // non-integrable (every reseeded trajectory diverges), so neither the energy
    // gate nor the rescue ceiling would fire -- and the move loops forever. This
    // counter guarantees no pose is a permanent trap: after this many stuck rounds
    // the ligand is relocated unconditionally, then the counter resets.
    int maxStuckRounds = 25;
};

class World {
    public:
    World(int index, bool cartesian, std::uint32_t seed);

    // Build the immutable RobotModel once from the shared topology + this
    // world's per-bond mobilities + this world's root mobilities. NOTE: the
    // root-mobility vector is now a PER-WORLD argument (Context passes a global
    // default for ordinary worlds and a docking-specific one for docking
    // worlds), which is what makes "root mobility is a world property" true.
    void buildModel(const SystemTopology& sys,
                    const Selection& sel,
                    const std::vector<RootMobility>& rootMobilities);

    // Python: world.add_sampler(timeStep, mdSteps, acceptRejectMode, use_nuts,
    //                           sphere_factor=1.0, use_fixman=None)
    // sphere_factor scales the AUTO radius R_i = R_receptor + sphere_factor*
    // R_ligand_i per ligand; larger => guest may roam further past the receptor
    // surface. use_fixman == nullopt -> AUTO: ON for non-Cartesian worlds
    // (torsional + docking), forced OFF for Cartesian. boostMDSteps and the
    // kick_stride / local_kick_* knobs were removed (radius is now automatic and
    // the kick is conditional on leaving the sphere).
    World& add_sampler(double timeStep,
                       int mdSteps,
                       AcceptRejectMode mode,
                       bool useNuts,
                       double sphereFactor,
                       std::optional<bool> useFixman,
                       bool alwaysKick,
                       double clashThreshold,
                       int maxInitialKickTries);

    // Mark this world as a docking world. ligandGroups[i] is the global atom
    // index list of ligand molecule i (each gets its own auto-sized sphere and is
    // repositioned independently); siteAtoms is the receptor atom list whose
    // centroid is the sphere centre. Called by Context::addDockingWorld.
    void configureDocking(std::vector<std::vector<int>> ligandGroups, std::vector<int> siteAtoms);

    void setAtomsLocationsInGround(const std::vector<robo::Vec3>& atomPosG);
    [[nodiscard]] const robo::Vec3* getAtomsLocationsInGround() const {
        return state_.atomPosG();
    }

    // One Gibbs sample on this world. Returns true if accepted.
    bool generateSample();

    // Docking initialisation: keep drawing random ligand placements until the
    // post-kick PE change is below the clash ceiling, or the retry budget
    // (sampler_.maxInitialKickTries) is exhausted.  Returns the number of
    // attempts needed.  Throws if no clean pose was found within the budget.
    // Called by Context::runREX before round 0 when maxInitialKickTries > 0.
    int findGoodStartingPose();

    [[nodiscard]] const RobotModel& model() const {
        return model_;
    }
    [[nodiscard]] RobotState& state() {
        return state_;
    }
    [[nodiscard]] int index() const {
        return index_;
    }
    [[nodiscard]] bool isCartesian() const {
        return cartesian_;
    }
    [[nodiscard]] bool isDocking() const {
        return docking_;
    }

    // --- last-move telemetry (for the [rex] log) ---
    [[nodiscard]] bool lastAccepted() const {
        return lastAccepted_;
    }
    // Set by Context::runREX each round. During equilibration every world
    // behaves as AlwaysAccept regardless of its sampler_.acceptRejectMode,
    // so the burn-in warms up from the starting geometry without the
    // Metropolis gate blocking large-dH moves from an unrelaxed structure.
    void setEquilPhase(bool equil) {
        equilPhase_ = equil;
    }

    [[nodiscard]] bool lastKickApplied() const {
        return lastKickApplied_;
    }
    [[nodiscard]] double lastPE() const {
        return state_.energy.pe;
    }
    [[nodiscard]] double lastKE() const {
        return state_.energy.ke;
    }
    [[nodiscard]] double lastFixman() const {
        return state_.energy.fixman;
    }
    [[nodiscard]] double lastTotalEnergy() const {
        return state_.energy.total;
    }
    [[nodiscard]] const char* typeName() const {
        return cartesian_ ? "cartesian" : (docking_ ? "docking" : "torsional");
    }

    void setTemperature(double T);

    // Kinetic-metric preconditioning (fictitious mass; sampling only). Inflates
    // the spatial inertia used ONLY in the proposal (momentum draw, KE, Fixman
    // ln det M) for selected bodies, raising their stable dt ~sqrt(scale) with
    // zero configurational bias (see RobotModel::bodyMassScale). Call after
    // buildModel(). Typical use: setMassScaleByJoint(JointType::Free, 16.0) on a
    // solvent world to tame water libration.
    void setBodyMassScale(int body, double scale);
    void setMassScaleByJoint(JointType jt, double scale);

    private:
    // --- internal (torsional) HMC pieces ---
    void reinitialize(); // seed velocities, record initial H (incl. Fixman)
    bool metropolis(double Hold, double Hnew);
    double currentTotalEnergy(); // PE (OpenMM) + KE (engine) [+ Fixman - 1/2 RT logSineSqr]
    void recomputeGeometry(const robo::Vec3* targets);

    // --- Fixman / coordinate-Jacobian corrections (torsional worlds only) ---
    // calcFixman   = 1/2 RT ( ln|M_phi| - ln|M_3N| ); ln|M_phi| via the O(n) ABA
    //                D-determinant (matrix-free, ALL bodies). |M_3N| is constant
    //                and cancels in dH, kept for an absolute value.
    // calcLogSineSqrGamma2 = sum over EVERY root body (not just molecule 0) of
    //                log(sin^2(pitch)) from the root-atom quaternion.
    double calcFixman();
    double calcLogSineSqrGamma2() const;

    // --- docking ---
    // Reposition ligands (pure proposal; relocate + reorient, no acceptance here).
    // forceAll=false: containment -- move a ligand only if its COM has left its
    // sphere. true (always_kick): perturb every ligand every round. The move is
    // judged in generateSample (Hold at pre-kick + energy validity check).
    // Returns true if any ligand was perturbed.
    bool repositionLigands(bool forceAll);
    double groupSphereRadius(int g) const;                             // R_receptor + sphereFactor*R_ligand_g
    robo::Vec3 sampleUniformInSphere(double radius);                   // uniform in a ball
    robo::Rotation sampleUniformRotation();                            // uniform on SO(3)
    robo::Vec3 atomSetCentroid(const std::vector<int>& atoms) const;   // unweighted centre
    robo::Vec3 atomSetMassCenter(const std::vector<int>& atoms) const; // mass-weighted COM
    double atomSetRadius(const std::vector<int>& atoms, const robo::Vec3& center) const;

    int index_;
    bool cartesian_;
    bool docking_ = false;
    RobotModel model_;
    RobotState state_;
    ForceBridge bridge_;
    SamplerConfig sampler_;
    robo::ConstraintSet constraints_;

    // Docking bookkeeping (global atom indices, OpenMM order). One atom group per
    // ligand molecule; the receptor atoms define the sphere centre.
    std::vector<std::vector<int>> ligandGroups_;
    std::vector<int> siteAtoms_;
    bool lastAccepted_ = false;
    bool lastKickApplied_ = false;
    bool equilPhase_ = false; // true during burn-in: AlwaysAccept overrides MH

    // Consecutive rejected docking moves from the current carried-forward pose.
    // Reset on any acceptance; when it reaches sampler_.maxStuckRounds a kick is
    // forced unconditionally so a locally non-integrable pose cannot trap the run.
    int dockingStuckCount_ = 0;

    double temperature_ = 300.0;
    double RT_ = 0;
    double beta_ = 0;

    // Constant Cartesian mass-matrix log-determinant ln|M_3N| = 3 * sum ln m_a.
    double lnDetMCartesian_ = 0.0;

    std::mt19937_64 rng_;
    std::normal_distribution<double> gaussian_{0.0, 1.0};
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};

    std::vector<robo::Vec3> savedPosG_;
    std::vector<double> savedQ_;
    double Hold_ = 0;

    std::vector<robo::Transform> frameFlat_;
    std::vector<robo::Transform> xpcFlat_;
};