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
    RigidKick,
    NcmcSwitch // NEW: lambda:1->0->1 alchemical decouple-move-recouple (per-molecule)
};

// Velocity-distortion option for the HMC momentum draw, ported from the simtk
// Robosample "NMA scaling" (HMCSampler::setVelocitiesToNMA). Before the
// sqrt(M^-1) map, each generalized-speed component of the white-noise draw is
// scaled by its per-DOF NMA factor and the overall vector length is restored by
// sqrt(nu)/||uScale|| so the kinetic temperature is preserved. The source read
// the per-body factors from World::getMobodUScaleFactor but never populated them,
// so the effective factors are unity and the distortion reduces to the plain
// Gaussian draw unless non-unit factors are supplied. A nullopt option (Python
// None, the default) means the distortion is not applied at all.
enum class DistortOption : std::uint8_t {
    NMA = 0
};

// Per-bond mobility selection produced by Context::build_flexibilities. The
// per-bond joint type is now a JointType (the former BondMobility enum is gone;
// Rigid/Torsion remain usable as JointType aliases). Default Weld (== Rigid)
// welds the two bonded atoms into the same rigid unit.
struct Selection {
    std::vector<JointType> bondMobility; // [numBonds]; default Weld (== Rigid)
};

// Per-frame, per-BODY spatial force row (docs/specs/reaction-force-
// monitoring.md Sec.1/3/4): `force`/`torque`, about the body origin Bo, in
// Ground, is the SUM of whichever term(s) World::enableReactionReporter
// enabled -- the OpenMM net applied force (`bodyForceG`, read directly off
// `ForceBridge::evaluate`) and/or the static (u=0) mobilizer reaction
// (`RobotEngine::calcMobilizerReactionForces`). Both terms are spatial
// forces about the SAME point (Bo, in Ground), so summing them is a valid
// spatial-force addition, not an apples-to-oranges combination; which
// term(s) are included is a per-reporter-world CHOICE, not encoded in the
// row itself (docs/specs/reaction-force-monitoring.md Sec.1.2). RAW (native
// units, no normalization -- normalization, if wanted, is a render-only
// transform applied downstream of this CSV, in vmd/arrows.tcl, never here).
// `bodyIdx` is the RobotModel body index; `atomIdx` is the VMD/DCD atom
// index (prmtop order, SystemTopology::atomsPrmtopIndex) of a
// REPRESENTATIVE atom of that body, for placement in vmd/arrows.tcl.
struct ReactionSample {
    int bodyIdx;
    int atomIdx;
    robo::Vec3 force;
    robo::Vec3 torque;
};

struct SamplerConfig {
    double timeStep = 0.001;
    int mdSteps = 0;
    AcceptRejectMode acceptRejectMode = AcceptRejectMode::AlwaysAccept;
    bool useFixman = false; // include the Fixman potential (1/2 RT ln det M) in H
    // External-rotation Jacobian -(1/2)RT ln sin^2(gamma2) for free/ball roots.
    // DEFAULT OFF: it is the Euler-angle (sin theta) volume factor and is WRONG
    // for Robosample's UNIT-QUATERNION roots, where the flat S^3 measure already
    // equals Haar on SO(3) (no correction needed). See World::calcLogSineSqrGamma2
    // and tests/TestEnsembleOrientation.cpp. Kept as a flag for any future
    // Euler-parameterized joint where the term would be correct.
    bool useOrientationJacobian = false;
    MoveType moveType = MoveType::MdHmc;
    bool useNuts = false;

    // Velocity-distortion applied at momentum-draw time (World::reinitialize).
    // nullopt (Python None, the default) => plain Gaussian draw, no distortion.
    // DistortOption::NMA => apply the per-DOF NMA velocity scaling (see the
    // DistortOption enum). Ported from the simtk HMCSampler DistortOpt path.
    std::optional<DistortOption> distortOption = std::nullopt;
    // NMA Route B bias magnitude alpha (DistortOption::NMA only). The momentum
    // mixture is Us = z +/- alpha * uhat, uhat a UNIT direction, so alpha is the
    // directed push in thermal-sigma units along that direction (||bias||^2 =
    // alpha^2, NOT nu). Larger alpha => bolder proposals, lower acceptance. See
    // add_sampler docstring for guidance; 1.0 is one thermal sigma.
    double nmaBiasScale = 1.0;
    // NCMC λ=0 trough teleport (tests/TestNcmcTeleport): at the decoupled (ghost)
    // trough, rigidly reposition the region root with KE-preserving velocity
    // co-rotation, giving the implicit-docking long-range teleport with explicit-
    // solvent correctness. DEFAULT OFF (the existing uncaged-stride behaviour).
    bool ncmcTeleport = false;
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

    // Reversibility diagnostic cadence. 0 = OFF (default; preserves zero overhead).
    // N > 0 = run RobotEngine::checkReversibility every N generateSample() calls
    // (the first call is round 0, so N>0 also performs the startup check). The probe
    // integrates mdSteps forward + back at the world's timeStep from the freshly
    // seeded state, logs the relative round-trip residual, and warns if it is large
    // or non-finite. It is a NON-DESTRUCTIVE smoke test that certifies dt only for
    // the current configuration (the safe dt is configuration dependent); the
    // always-on guard remains the per-step corrector throw in verletStep. See
    // THEORY 5.7. Set from Python via the reversibility_check_every=N argument to
    // context.add_*_world() (which forwards to World::setReversibilityCheck).
    int reversibilityCheckInterval = 0;

    // ---- NCMC (MoveType::NcmcSwitch) -----------------------------------
    // Total protocol substeps; lambda ramps 1->0 then 0->1, one Verlet step per
    // substep (the torsional stride happens near the lambda~0 midpoint where the
    // mobile molecule is uncaged from all others). 0 => not an NCMC world.
    int ncmcSteps = 0;
    // Fraction of ncmcSteps held at lambda=0 between the down- and up-ramps.
    double ncmcHoldFraction = 0.0;
    // Atom-index SET (global/OpenMM order, ASCENDING + deduplicated) of Region A --
    // the substructure decoupled from every other atom (docs/specs/
    // ncmc-explicit-solvent/30-region-and-protocol-policy.md Sec.2). Set by
    // World::configureNcmc; a contiguous [begin,end) range is the common case (the
    // convenience overload builds this vector from it) but the set need not be
    // contiguous.
    std::vector<int> ncmcAtomIndices;
    // Construction II (Metropolized-dynamics NCMC, docs/specs/
    // ncmc-explicit-solvent/10-acceptance-construction.md): when true, each
    // fixed-lambda propagate substep is itself Metropolis accept/reject'd against
    // the full H_lambda (an inner GHMC kernel), and the OUTER move accepts on the
    // protocol work W alone (min(1,exp(-beta*W))) instead of the endpoint
    // Hend-Hstart. Default false preserves Construction I (endpoint-DeltaH, the
    // original behavior) bit-for-bit.
    bool useMetropolizedInner = false;
};

class World {
    public:
    World(int index, bool cartesian, std::uint32_t seed);

    // Build the immutable RobotModel once from the shared topology + this
    // world's per-bond mobilities + this world's root mobilities. NOTE: the
    // root-mobility vector is now a PER-WORLD argument (Context passes a global
    // default for ordinary worlds and a docking-specific one for docking
    // worlds), which is what makes "root mobility is a world property" true.
    void
    buildModel(const SystemTopology& sys, const Selection& sel, const std::vector<JointType>& rootMobilities);

    // Override this world's root attachment to Ground. Root mobility is a
    // PER-WORLD property: these mutate only THIS world's copy of the
    // root-mobility vector and rebuild this world's RobotModel in place; they
    // never touch the shared SystemTopology, so different worlds can give the
    // same molecule different roots (e.g. a solvation shell where near waters
    // are Free and the bulk is Welded). buildModel(...) MUST have been called
    // once (by Context::add*World) first. Call these BEFORE add_sampler / mass
    // scaling, since rebuilding resets the per-body sampler state.
    //
    //   setRootMobility   : change ONE molecule (rebuilds once).
    //   setRootMobilities : replace the WHOLE vector and rebuild ONCE -- the
    //                       O(N)-not-O(N^2) path for setting many molecules
    //                       (e.g. thousands of solvent molecules at once).
    void setRootMobility(int moleculeIndex, JointType mobility);
    void setRootMobilities(const std::vector<JointType>& rootMobilities);

    // Solvent-relaxing NCMC (docs/specs/ncmc_solvent_relax.md). Mark a set of
    // atoms (global/OpenMM order) to be advanced in FLAT Cartesian space by
    // velocity-Verlet driven by OpenMM forces INSIDE the proposal, so the contact
    // environment relaxes during the move instead of being a welded wall. The
    // atoms stay welded as rigid bodies (zero generalized DOF, no Fixman/Jacobian
    // contribution); only their per-atom Cartesian (x,v) move. Empty set (default)
    // == the welded engine, bit-for-bit. Call AFTER add_sampler (the rebuild in
    // setRootMobilities clears per-body state but not this runtime set).
    void setCartesianSolvent(const std::vector<int>& atomIndices);

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
                       int maxInitialKickTries,
                       std::optional<DistortOption> distortOption = std::nullopt,
                       double nmaBiasScale = 1.0);

    // Mark this world as a docking world. ligandGroups[i] is the global atom
    // index list of ligand molecule i (each gets its own auto-sized sphere and is
    // repositioned independently); siteAtoms is the receptor atom list whose
    // centroid is the sphere centre. Called by Context::addDockingWorld.
    void configureDocking(std::vector<std::vector<int>> ligandGroups, std::vector<int> siteAtoms);

    // Mark this world as a per-molecule NCMC world: during the lambda:1->0->1
    // switch the intermolecular nonbonded between Region A (atomIndices) and
    // every other atom is alchemically softened (intramolecular physics
    // untouched). atomIndices is an arbitrary atom-index set (need not be
    // contiguous, docs/specs/ncmc-explicit-solvent/30-region-and-protocol-
    // policy.md Sec.2); it is sorted + deduplicated internally. Call AFTER
    // add_sampler (it sets moveType last).
    void configureNcmc(std::vector<int> atomIndices, int ncmcSteps, double holdFraction);
    // Convenience overload: contiguous [atomBegin,atomEnd) Region A (the common
    // case -- a single contiguous solute block). Forwards to the index-set form.
    void configureNcmc(int atomBegin, int atomEnd, int ncmcSteps, double holdFraction);

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

    // Number of velocity/momentum coordinates actually drawn for this world
    // (ensemble-validation foundations spec, Sec. 2). Internal/torsional world:
    // nu - n_C, where n_C is the SAME loop-closure constraint count consumed by
    // calcConstraintLogDet (constraints_.numConstraints()) -- the removed
    // holonomic velocity constraints. Cartesian world: model_.nu is NOT the
    // physical dof (Cartesian worlds are modelled with nu==1 body, see
    // buildModel), so this instead reads OpenMM's own dof count
    // (3*N_real - constraints - 3*[CMMotionRemover present]) off the shared
    // OpenMMContext system.
    [[nodiscard]] int nDof() const {
        if (cartesian_) {
            return OpenMMContext::get().getNumDegreesOfFreedom();
        }
        return model_.nu - constraints_.numConstraints();
    }

    // Number of loop-closure DistanceConstraints (ConstraintSet::distance) this
    // world's molecule graph produced -- 0 on every acyclic molecule, >=1 once a
    // ring-closing bond was cut into a RATTLE constraint (World.cpp ~:679). Exposed
    // for the Tier-2 cyclic-path structural guard (docs/specs/
    // fixman-idealized-chains-validation.md T2.0): no other Python-reachable
    // signal distinguishes "the loop path never fired" from "it fired and is 0".
    [[nodiscard]] int numLoopConstraints() const {
        return constraints_.numConstraints();
    }

    // Loop-closure Fixman term ln det(G M^-1 G^T) at the CURRENT geometry
    // (Constraints.hpp::calcConstraintLogDet) -- the same quantity calcFixman()
    // folds into the acceptance Hamiltonian, exposed standalone (no sampling
    // round required) so a structural test can assert it is non-zero and
    // conformation-dependent on a cyclic molecule. Realizes position + articulated
    // -body inertias at the state's CURRENT q first (idempotent, cheap: this is
    // exactly calcFixman()'s own precondition). Returns 0 verbatim on an acyclic
    // world (ConstraintSet::calcConstraintLogDet's own no-op contract).
    [[nodiscard]] double currentConstraintLogDet() {
        RobotEngine::realizePosition(model_, state_);
        RobotEngine::realizeArticulatedBodyInertias(model_, state_);
        return constraints_.calcConstraintLogDet(model_, state_);
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

    // Enable/disable the periodic reversibility probe (THEORY 5.7). interval<=0
    // disables it; interval>0 runs it every `interval` generateSample() calls,
    // starting at round 0. Bound to Python as world.set_reversibility_check(...).
    void setReversibilityCheck(int interval);

    // Enable the NCMC λ=0 trough teleport (default off). See SamplerConfig::ncmcTeleport.
    void setNcmcTeleport(bool on) {
        sampler_.ncmcTeleport = on;
    }

    // Select Construction II (Metropolized-dynamics NCMC) for the outer
    // acceptance. Default false == Construction I (endpoint-DeltaH, unchanged).
    // See SamplerConfig::useMetropolizedInner and docs/specs/
    // ncmc-explicit-solvent/10-acceptance-construction.md.
    void setNcmcConstructionII(bool on) {
        sampler_.useMetropolizedInner = on;
    }

    // ---- Reaction-force reporter (docs/specs/reaction-force-monitoring.md) --
    // Mark this world the "reporter" for per-body force monitoring.
    // `interestingBodies` (RobotModel body indices, Ground (body 0) excluded --
    // Sec.2.1/§4) is the set whose spatial force `captureReactionSnapshot()`
    // records. Off by default -- a world that never calls this allocates
    // nothing and does zero extra work. Throws if this world is Cartesian
    // (OpenMMVelocityVerlet-equivalent): the internal-coordinate articulated
    // body indexing this reporter relies on is not meaningful there (Sec.3
    // integrator guard). Does NOT touch reactionIncludeOpenmm_/
    // reactionIncludeReaction_ (they stay whatever they already were,
    // defaults true/false) -- those flags are set exclusively by
    // enableReactionReporter's includeOpenmm/includeReaction arguments.
    void setReactionReporter(std::vector<int> interestingBodies);
    // Self-contained convenience: derive the interesting-body set from THIS
    // world's CURRENT model_ (every flexed/non-Weld body plus its parent,
    // excluding only Ground itself -- bodyForceG needs no inboard atom, so a
    // body whose own inboard joint attaches directly to Ground, e.g. a
    // Free-rooted receptor's root body, IS included -- same rule as
    // enableReactionReporter's own doc comment) and call setReactionReporter
    // with it. Safe to call at ANY point after buildModel/setRootMobility(ies) has
    // produced the FINAL model this world will sample with -- a rebuild
    // changes body indices, so calling this BEFORE the LAST rebuild derives a
    // stale set (e.g. add_ncmc_world's set_root_mobilities rebuild happens
    // AFTER add_robotic_world, docs/specs/reaction-force-monitoring.md
    // Sec.2.0).
    //
    // `reportFreeBodies` (default true, non-breaking) controls whether
    // FREE-FLOATING RIGID bodies survive the derivation: a body that is BOTH
    // Ground-rooted (bodyParent[b] == 0) AND CHILDLESS (no other body's
    // parent is b) is a lone rigid molecule sitting on a Free root -- e.g. a
    // lipid given a Free root but no internal DOFs. When false, such bodies
    // are dropped from the derived set AFTER collection, even though they
    // are "flexed" (bodyNU != 0 from the Free root) and would otherwise pass
    // the loop above. A flexed molecule's OWN root (e.g. a receptor's root
    // TM body) is Ground-rooted too but is NOT childless -- its internally-
    // jointed bodies branch from it -- so it is unaffected and always kept.
    //
    // `includeOpenmm` (default true) / `includeReaction` (default false)
    // select which spatial-force TERM(S) captureReactionSnapshot() SUMS into
    // each body's single (force, torque) row (docs/specs/
    // reaction-force-monitoring.md Sec.1.2): includeOpenmm adds the OpenMM
    // net applied force `bodyForceG`; includeReaction adds the static (u=0)
    // mobilizer reaction. Both terms are spatial forces about the SAME point
    // (Bo, in Ground), so summing them is a valid addition. DEFAULT
    // (true, false) reproduces the OpenMM-only output this reporter shipped
    // with, byte-for-byte -- no forward-dynamics step, no u save/zero/
    // restore. At least one of the two MUST be true (enforced in
    // captureReactionSnapshot, not here, since the flags can still change via
    // a later enableReactionReporter call before any snapshot is taken).
    // Bound to Python as world.enable_reaction_reporter(report_free_bodies=True,
    // include_openmm=True, include_reaction=False).
    void enableReactionReporter(bool reportFreeBodies = true,
                                bool includeOpenmm = true,
                                bool includeReaction = false);
    [[nodiscard]] bool isReactionReporter() const {
        return reactionReporter_;
    }
    [[nodiscard]] const std::vector<int>& reactionInterestingBodies() const {
        return reactionInterestingBodies_;
    }

    // Snapshot, per interesting body, the SUM of whichever term(s)
    // enableReactionReporter enabled -- the OpenMM net applied force
    // (`bodyForceG`, velocity-free) and/or the static (u=0) mobilizer
    // reaction -- on this world's CURRENT q (docs/specs/
    // reaction-force-monitoring.md Sec.3) into reactionSamples(). Intended
    // call point: AFTER generateSample() has resolved this round's
    // accept/reject, so the q read here is the same accepted conformation
    // that gets written to the paired DCD frame (Sec.1 item 3). No-op
    // (leaves reactionSamples() empty) unless isReactionReporter(). Throws
    // if BOTH includeOpenmm and includeReaction are false (nothing to
    // report). READ-ONLY w.r.t. the sampler: `bodyForceG`'s refresh only
    // touches ForceBridge/arena scratch (recomputed fresh at the top of the
    // next round anyway); the OPTIONAL u=0 reaction term DOES zero `u` for
    // its forward-dynamics step, but saves and restores it (and every
    // u-derived cache) before returning -- so either way this cannot perturb
    // the next round's generateSample() call (the no-perturbation property,
    // Sec.3/6).
    void captureReactionSnapshot();
    // The force rows captured by the last captureReactionSnapshot() call
    // (empty until that has been called at least once, or on a non-reporter
    // world). Read-only -- Context flushes this to the per-replica
    // reactions.csv at the DCD cadence (Sec.4).
    [[nodiscard]] const std::vector<ReactionSample>& reactionSamples() const {
        return reactionSamples_;
    }

    // Compute the Route B mass-weighted internal-coordinate Hessian at the given
    // (minimized) Ground-frame coordinates and store the softest non-trivial mode
    // as the per-DOF NMA direction consumed by reinitialize(). atomPosGFlat is
    // [x0,y0,z0, x1,...] in nm, global/OpenMM atom order. Returns omega^2 of the
    // chosen mode (sanity: should be > 0 at a real minimum). No-op on Cartesian worlds.
    double setNMASoftModeFromHessian(const std::vector<double>& atomPosGFlat,
                                     double h = 1e-5,
                                     double zeroTol = 1e-6);

    private:
    // --- internal (torsional) HMC pieces ---
    void reinitialize(); // seed velocities, record initial H (incl. Fixman)
    bool metropolis(double Hold, double Hnew);
    double currentTotalEnergy(); // PE (OpenMM) + KE (engine) [+ Fixman - 1/2 RT logSineSqr]
    // Cartesian solvent helpers (no-ops when no atom is Cartesian-integrated).
    void drawSolventVelocities(); // v_s ~ N(0, RT/m_s) Maxwell-Boltzmann
    double calcSolventKE() const; // 1/2 sum_s m_s |v_s|^2  -> state_.energy.keSolvent
    // NMA Route B kinetic correction RT*ln cosh(w.mu), w = M^(1/2) u / sqrt(RT).
    // Subtracted from the kinetic term in both Hold_ and Hnew so the acceptance
    // uses ke_mix = -RT*ln g(u|q) for the biased-mixture momentum draw. Returns 0
    // (strict no-op) unless DistortOption::NMA is active and nmaBias_ is set.
    double nmaKineticCorrection();
    bool ncmcMove();                    // the lambda-protocol HMC move (NcmcSwitch)
    double protocolLambda(int s) const;     // lambda schedule, s in [0, ncmcSteps)
    int ncmcTeleportRoot() const;           // Free-joint tree root of the NCMC region, or -1
    void ncmcApplyTroughTeleport(int rootBody); // rigid λ=0 teleport + KE co-rotation
    // Construction II inner kernel (docs/specs/ncmc-explicit-solvent/
    // 10-acceptance-construction.md Sec.3, 20-inner-integrator.md Sec.2 F3): one
    // fixed-lambda Verlet substep, Metropolized against the FULL H_lambda (the
    // SAME assembly currentTotalEnergy() uses -- F1/F2), reject-flip (u,v_s) on
    // reject. A non-converged velocity corrector is ALSO an automatic reject
    // (F3). Returns false ONLY on a genuinely unrecoverable non-finite
    // force/velocity (mirrors RobotEngine::stepTo's own contract); the caller
    // treats false exactly as it already treats a raw stepTo failure (abort the
    // whole NCMC move). A Metropolis reject or non-convergence is absorbed
    // internally and reported as true (the state is left valid, momentum-flipped).
    // `substepIndex` is diagnostic only (per-substep trace, ROBO_NCMC_DEBUG=1).
    // `acceptedOut`, if non-null, is set to whether the inner GHMC accepted this
    // substep (false on an F3 non-convergence reject too), so the caller can
    // accumulate a per-move inner-acceptance count.
    bool ncmcInnerGhmcStep(robo::Real h, int substepIndex, bool* acceptedOut = nullptr);
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
    long generateSampleCalls_ = 0; // for the reversibility-check cadence (SamplerConfig)
    bool equilPhase_ = false;      // true during burn-in: AlwaysAccept overrides MH

    // Reaction-force reporter state (docs/specs/reaction-force-monitoring.md).
    // reactionReporter_/reactionInterestingBodies_ are set once by
    // setReactionReporter(); reactionIncludeOpenmm_/reactionIncludeReaction_
    // are set by enableReactionReporter's includeOpenmm/includeReaction
    // arguments (not by setReactionReporter -- see its doc comment); they
    // select which term(s) captureReactionSnapshot() sums into each row
    // (Sec.1.2); reactionSamples_ is captureReactionSnapshot()'s output
    // buffer.
    bool reactionReporter_ = false;
    bool reactionIncludeOpenmm_ = true;
    bool reactionIncludeReaction_ = false;
    std::vector<int> reactionInterestingBodies_;
    std::vector<ReactionSample> reactionSamples_;

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

    // Per-DOF NMA velocity-scale factors for DistortOption::NMA. Sized to nu and
    // lazily initialised to 1.0 in reinitialize(); the simtk source read these
    // per mobilized body (getMobodUScaleFactor) but never populated them, so the
    // effective default is unity (NMA distortion == plain Gaussian draw).
    std::vector<robo::Real> uScaleFactors_;

    // NMA Route B bias mu = uScaleFactors_ * sqrt(nu)/||uScaleFactors_|| (length-
    // restored; unit factors => mu = ones). Set per block in reinitialize(); the
    // momentum draw is the mixture Us = z +/- mu, and the acceptance correction
    // RT*ln cosh(w.mu) uses this same mu.
    std::vector<robo::Real> nmaBias_;

    std::vector<robo::Vec3> savedPosG_;
    // Cartesian solvent positions saved at the start of an NCMC move and restored
    // on reject (the body->atom fill does not touch Cartesian-integrated atoms).
    std::vector<robo::Vec3> savedSolventPosG_;
    std::vector<double> savedQ_;
    double Hold_ = 0;

    std::vector<robo::Transform> frameFlat_;
    std::vector<robo::Transform> xpcFlat_;

    // Build inputs retained so setRootMobility(ies) can rebuild this world's
    // model in isolation. sys_ points at the Context-owned SystemTopology
    // (outlives every World); sel_ and rootMobilities_ are this world's own
    // copies. Captured on the first buildModel() call.
    const SystemTopology* sys_ = nullptr;
    Selection sel_;
    std::vector<JointType> rootMobilities_;
};