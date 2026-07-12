#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "BatScaling.hpp"
#include "DCDWriter.hpp"
#include "OpenMMContext.hpp"
#include "ReplicaExchange.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"

// Robosample top-level context. Constructed from Python as Context(base_name,
// seed); the Python subclass (context.py) adds the dihedral classifier and
// load_amber, which fills `systemTopology` in place.
//
// Owns the worlds (the Gibbs flexibility regimens), the replica temperature
// ladder, and the REX driver. World order defines the Gibbs sweep order, and
// per-atom Ground coordinates (nm) are the only currency passed between worlds.
class Context {
    public:
    Context(std::string baseName, std::uint32_t seed);

    // Bound to Python as `system_topology`. Filled in place by load_amber.
    SystemTopology systemTopology;

    // ---- modeling --------------------------------------------------------
    // Root mobility is a PER-WORLD property. systemTopology.rootMobilities is
    // only the build-time DEFAULT each world seeds from (set by load_amber);
    // there is deliberately no Context-level mutator. To change an individual
    // molecule's root attachment, call World::setRootMobility(...) on the world
    // returned by add*World (see World.hpp), which rebuilds that world's model
    // in isolation -- the same mechanism addDockingWorld already uses.

    // Add a Cartesian world (pure OpenMM MD on device) or a robotic/torsional
    // world (internal-coordinate HMC). Both return a reference to the new world
    // so Python can chain .add_sampler(...). The model is built immediately
    // from the current systemTopology + root mobilities.
    //
    // wantReactionReporter (docs/specs/reaction-force-monitoring.md Sec.2): opt
    // -in, off by default. Flags the new world the per-body applied-force
    // reporter (World::setReactionReporter). addCartesianWorld always THROWS
    // when true (Sec.3 integrator guard -- a Cartesian world's articulated
    // body indexing is not meaningful). addRoboticWorld derives the
    // interesting-body set from the just-built model: every non-Weld (flexed)
    // body plus its parent, excluding only Ground itself -- a Free-rooted
    // body (e.g. a receptor's root, directly attached to Ground) IS included
    // (Sec.2.1/Sec.4).
    World& addCartesianWorld(bool wantReactionReporter = false);
    World& addRoboticWorld(const Selection& sel, bool wantReactionReporter = false);

    // Add a DOCKING world. `ligandMoleculeIndices` lists which molecules are
    // ligands; their roots become Free (6 external DOF) and ALL their bonds stay
    // Rigid (rigid-body docking). Every other molecule is Welded to Ground and
    // rigid. Root mobility is thus a WORLD property here -- it is built from this
    // argument, NOT from systemTopology.rootMobilities, so the same system can
    // host a docking world and ordinary worlds at once. The binding-site centre
    // is the centroid of all non-ligand (receptor) atoms; pass a sphere radius
    // via add_sampler(...). Returns the world for chaining .add_sampler(...).
    World& addDockingWorld(const std::vector<int>& ligandMoleculeIndices);


    // Per-bond mobility selection. `bonds` empty/None => every eligible
    // (non-ring, non-terminal) bond gets `mobility`; otherwise only the listed
    // (i,j) bonds do. Ring-closing bonds always remain Rigid.
    Selection buildFlexibilities(const std::optional<std::vector<std::pair<int, int>>>& bonds,
                                 JointType mobility,
                                 bool flag);

    // ---- run -------------------------------------------------------------
    // Build the OpenMM system, set the replica temperature ladder, and seed
    // every replica with the reference coordinates. Empty list => single 300 K.
    void initialize(const std::vector<double>& temperatures);

    // Replica-exchange driver: `equil` then `prod` rounds; each round runs a
    // Gibbs sweep over all worlds for every replica, attempts adjacent-replica
    // swaps, and writes outputs every `writeFreq` production rounds.
    //
    // COORDINATE-swap REMC (INV-3 non-compliant by construction: it swaps
    // `replicaCoords_`, not labels). Retained verbatim as the
    // INVARIANT-EQUIV oracle for `RunREX` below (docs/specs/
    // replica-exchange-nonequilibrium-work.md, "Port-target baseline"). Do
    // NOT extend this method -- new run types and the label-swap object
    // model live in `RunREX`.
    void runREX(int equilRounds, int prodRounds, int writeFreq, bool verbose);

    // ---- label-swap replica exchange (docs/specs/replica-exchange-
    // nonequilibrium-work.md) ----------------------------------------------
    // Stage 1: RUN_TYPE::REMC (parallel tempering, B6). RUN_TYPE::Default
    // runs independent replicas (no exchange attempts). Stage 2b: RUN_TYPE::
    // RENE/REBASONTOP (driven BAT-scaling exchange, B5/B6/B8, WORK_*
    // accumulation, F4 atomic commit) are fully wired. RUN_TYPE::RENEMC's
    // ACCEPTANCE formula (ETerm_nonequil, B6) is wired in attemptREXSwap, but
    // its own driven ROUND-LOOP (the velocity/NMA drive segment) is a
    // Stage 2c TODO -- RunREX THROWS std::logic_error for RENEMC (the
    // acceptance math is directly testable via attemptREXSwap in isolation,
    // see tests/TestRexAcceptanceAlgebra.cpp).
    //
    // NONE of the RENE/REBASONTOP driven-round code below (this method, the
    // private runDrivenRound/driveReplica/runInterleavedRemcSubround/
    // checkInv7AndInv10Guards helpers, and attemptREXSwap's RENE/REBASONTOP/
    // RENEMC branches) has been compiled or run (coordinator directive,
    // 2026-07-12 "drop compiling and running entirely") -- treat every claim
    // below as a reviewed-on-paper design, not a build-confirmed one.
    //
    // Builds R = temperatures_.size() Replica/ThermodynamicState objects (B0,
    // asserting R == T) from the current worlds' schedule and the reference
    // coordinates seeded by `initialize()`, then runs `equil + prod` Gibbs
    // sweeps. Each round: every replica is propagated through its
    // ThermodynamicState's world schedule (temperature/timeStep/mdSteps/
    // acceptRejectMode reset onto the shared worlds via the runtime setters,
    // World::setTimeStep/setMdSteps/setAcceptRejectMode -- Consequences
    // "Feasibility gap"); `mixReplicas` (REMC/Default) or `runDrivenRound`
    // (RENE/REBASONTOP) then attempts a round of exchanges (B7), swapping
    // LABELS (`swapThermodynamicStates`, INV-3), never coordinates. Output
    // CSV/DCD files are indexed by THERMODYNAMIC STATE (not by replica-object
    // identity), matching the coordinate-swap `runREX`'s convention that a
    // fixed output slot is a fixed temperature -- this is what makes the two
    // drivers' per-state statistics directly comparable (INVARIANT-EQUIV, a
    // Stage-1 claim unaffected by the Stage-2b additions below).
    //
    // INV-7/V9 and INV-10 preconditions (checkInv7AndInv10Guards) are
    // asserted for RENE/REBASONTOP right after setup, before the round loop.
    void RunREX(RUN_TYPE runType, int equilRounds, int prodRounds, int writeFreq, bool verbose);

    // Attempt one swap between thermodynamic states thermoC/thermoH (B6).
    // REMC: accept iff ETerm_equal = -(beta_H - beta_C)(refU_X - refU_Y) >= 0
    // or U(0,1) < exp(ETerm_equal); on accept, swap LABELS only
    // (`swapThermodynamicStates`, INV-3).
    // RENEMC: ETerm_nonequil, the SAME PT form on the DRIVEN-ENDPOINT
    // reference potentials (referenceWORK_potential), no Jacobian (INV-10:
    // volume-preserving drive).
    // RENE/REBASONTOP: WTerm = -(Work_X + Work_Y), Work_partner =
    // beta_target*U(x_partner^tau) - beta_source*U(x_partner^0) -
    // lnJac_partner (B6/D7; correctionTerm == 1 under INV-9, D2). On accept,
    // F4 atomically commits BOTH replicas' WORK_* trial to committed
    // (Replica::commitWorkAsFinal) BEFORE the label swap; REMC/Default never
    // call it (their accept is label-swap only, INV-3/B6 -- this is the F4
    // fix: the original unconditionally ran the WORK commit even for REMC,
    // reverting coordinates from an unpopulated WORK buffer).
    // A non-finite acceptance exponent (NaN or +-inf, e.g. from a Stage 2a
    // domain-invalid drive endpoint or a blown-up OpenMM PE) is an EXPLICIT
    // automatic reject (reviewer N2 fail-loud), logged to stderr, not a
    // silent NaN comparison.
    // THROWS for RUN_TYPE::Default (a direct call is a caller error).
    // Returns true iff accepted. Exposed (not private) so a reproducer can
    // drive it directly (V1/INVARIANT-EQUIV do so indirectly via RunREX; V3/
    // V4/V8 in tests/TestRexAcceptanceAlgebra.cpp call it directly, per the
    // spec's acceptance-algebra oracles -- this is how RENEMC's acceptance
    // formula is tested despite its round-loop being unwired, Stage 2c).
    bool attemptREXSwap(int thermoC, int thermoH);

    // INV-7/V9 (Fixman-in-sampler <=> Fixman-out-of-acceptance, D3) and
    // INV-10 (drive/run-type pairing) preconditions for a DRIVEN run type.
    // Takes `runType` explicitly (not the `runType_` member) and reads only
    // `worlds_` (no replicas_/thermodynamicStates_ needed), so a reproducer
    // can call this directly on a hand-built Context (worlds added, no
    // initialize()/RunREX needed) to test the guard in isolation --
    // tests/TestRexAcceptanceAlgebra.cpp's INV-10 case does exactly this.
    // THROWS std::logic_error on violation; no-op for REMC/Default (INV-7 as
    // stated is a driven-only precondition per the Stage 2b directive; REMC's
    // own Fixman-off risk is a separate, not-yet-enforced concern, see the
    // Context.cpp comment). Public so it is independently testable.
    void checkInv7AndInv10Guards(RUN_TYPE runType) const;

    // Mixing configuration (Interface I3). swapEvery gates exchange
    // frequency (attempt a mix only when round % swapEvery == 0);
    // nSwapAttempts is the draw count for ReplicaMixingScheme::All;
    // swapFixman is an OFF-by-default DIAGNOSTIC (D3): Fixman is computed by
    // the sampler but SHALL NOT enter `attemptREXSwap`'s acceptance exponent
    // (INV-7) -- this flag exists only for the port-target interface parity
    // required by I3 and has no effect on Stage 1's REMC acceptance.
    void setReplicaMixingScheme(ReplicaMixingScheme scheme) {
        mixingScheme_ = scheme;
    }
    void setSwapEvery(int n) {
        swapEvery_ = (n > 0) ? n : 1;
    }
    void setNSwapAttempts(int n) {
        nSwapAttempts_ = (n > 0) ? n : 1;
    }
    void setSwapFixman(bool enabled) {
        swapFixman_ = enabled;
    }

    // REBASONTOP interleave (D4): "RENE work-swaps plus periodic REMC
    // neighbour-swap sub-rounds on top". Every `interleaveRemcEvery_` driven
    // rounds, `rebasontopSubrounds_` REMC-style (ETerm_equal, committed
    // potentials only, no WORK_* touched) neighbour-swap sub-rounds run in
    // addition to (not instead of) the main WTerm swap attempt that round.
    // Defaults (10, 6) are this port's choice -- the original's own count
    // ("six neighbour-swap sub-rounds", B10) is preserved for the subround
    // count; the original never specified an every-N cadence for entering
    // the sub-loop (the whole block was dead code, B10), so 10 is a new,
    // documented default, not a ported constant.
    void setInterleaveRemcEvery(int n) {
        interleaveRemcEvery_ = (n > 0) ? n : 1;
    }
    void setRebasontopSubrounds(int n) {
        rebasontopSubrounds_ = (n > 0) ? n : 1;
    }

    // Symmetric attempted/accepted swap counts, indexed by thermodynamic
    // state (T x T; nofAttemptedSwapsMatrix_[i][j] == [j][i]). Populated by
    // `attemptREXSwap`; V1's "non-zero attempted-swap count and a plausible
    // acceptance ratio" oracle reads these.
    [[nodiscard]] const std::vector<std::vector<std::int64_t>>& attemptedSwapsMatrix() const {
        return nofAttemptedSwapsMatrix_;
    }
    [[nodiscard]] const std::vector<std::vector<std::int64_t>>& acceptedSwapsMatrix() const {
        return nofAcceptedSwapsMatrix_;
    }

    // ---- BAT-scaling shared/global anchor (INV-9, D2 revision 2) ---------
    // Context (NOT ThermodynamicState) owns the single running-mean anchor
    // the RENE/REBASONTOP drive scales the deviation around (B4/B3), so it
    // is state-INDEPENDENT by construction -- the precondition the paired
    // scaling map needs to be an exact involution (INV-9). Stage 2a provides
    // the accumulation + frozen-snapshot API; wiring it into the equilibrium
    // round loop (calling accumulateBatAnchorStats after every equilibrium
    // world visit) and freezing exactly one snapshot per exchange round are
    // Stage 2b concerns.
    //
    // Reads `world`'s CURRENT committed geometry (getAtomsLocationsInGround)
    // -- callers SHALL call this only after an EQUILIBRIUM move (distortOption
    // == nullopt), never on a driven (ScaleBendStretch) world's output (that
    // would feed the anchor from nonequilibrium samples, biasing it).
    void accumulateBatAnchorStats(const World& world) {
        const robo::Vec3* p = world.getAtomsLocationsInGround();
        const std::vector<robo::Vec3> pos(p, p + world.model().numAtoms);
        batAnchorStats_.accumulate(world.model(), pos);
    }
    // Frozen snapshot (INV-9): take ONE and pass it to every drive in a
    // round -- both partners of a swap pair SHALL read the SAME Snapshot.
    [[nodiscard]] robo::BatAnchorStats::Snapshot batAnchorSnapshot() const {
        return batAnchorStats_.snapshot();
    }
    void resetBatAnchorStats() {
        batAnchorStats_.reset();
    }

    // Enable multiple-timestep (r-RESPA) integration for the Cartesian world's
    // on-device OpenMM MD: slow forces (Nonbonded, GBSA, ...) are evaluated once
    // per outer step, fast bonded forces `innerSubsteps` times. No effect on the
    // torsional worlds (their forces are always the full sum). MUST be called
    // before initialize(). Forwards to OpenMMContext.
    void setMTS(bool enabled, int innerSubsteps);


    // ---- energy ingestion / validation (unchanged) ----------------------
    auto initializeOpenMM() -> bool;
    [[nodiscard]] auto calcOpenMMPotentialEnergy() -> double;
    [[nodiscard]] auto computePotentialEnergyByGroup()
        -> std::pair<double, std::vector<OpenMMContext::ForceGroupEnergy>>;

    [[nodiscard]] auto getBaseName() const -> const std::string& {
        return baseName;
    }
    [[nodiscard]] auto getSeed() const -> std::uint32_t {
        return seed;
    }
    [[nodiscard]] int numWorlds() const {
        return static_cast<int>(worlds_.size());
    }

    private:
    double openmmPotential(const std::vector<robo::Vec3>& coords);
    void writeOutputs(int replica, int round, bool verbose);
    // Shared CSV+DCD writer core, factored out of writeOutputs so RunREX's
    // state-indexed output can reuse the same periodic-imaging/DCD-scatter
    // logic on `coords` (a Replica's committed coordinates) instead of
    // `replicaCoords_`. `idx` is both the output-file slot (baseName.<idx>.csv
    // / dcdWriters_[idx]) and the row's "replica" column -- for `runREX` that
    // is the replica-object index (== its fixed temperature slot); for
    // `RunREX` it is the THERMODYNAMIC-STATE index (see RunREX doc comment).
    // Pure extraction: writeOutputs(replica, round, verbose) behaves exactly
    // as before.
    void writeOutputsCore(int idx, int round, bool verbose, const std::vector<robo::Vec3>& coords, double T);
    // Append a reporter world's captured per-body force rows (docs/specs/
    // reaction-force-monitoring.md Sec.4) to that replica's per-replica CSV,
    // tagged with the DCD frame index they pair with. No-op if rows is
    // empty. Always the fixed 10-column schema (Sec.1.2/4).
    void writeReactionRows(int replica, int frame, const std::vector<ReactionSample>& rows);

    // Refuse to start (or warn, if ROBO_ALLOW_BAD_START is set) when the input
    // geometry is non-finite or sterically clashing -- the docking world cannot
    // repair a frozen-receptor clash, so this would otherwise loop forever.
    void checkStartupGeometry();

    // ---- label-swap replica exchange (RunREX) private helpers ------------
    // (Re)build replicas_/thermodynamicStates_/the inverse maps/the swap
    // matrices from temperatures_ + the current per-world schedule + the
    // reference coordinates seeded by initialize() (B0, identity maps).
    // Asserts R == T (B0 NOTE) -- throws std::logic_error otherwise.
    void setupReplicaExchange(RUN_TYPE runType);
    // Label swap only (INV-3): swap(replica2ThermoIxs_[X], [Y]),
    // swap(thermo2ReplicaIxs_[thermoC], [thermoH]) where X/Y are the replicas
    // currently occupying thermoC/thermoH. No Z-matrix pointer to swap in
    // Stage 1 (B1: the Z-matrix table is shared/global already; per-replica
    // BAT values are a Stage 2 concept).
    void swapThermodynamicStates(int thermoC, int thermoH);
    // Neighbouring-pairs schedule (B7): startIdx = (round + oddity) % 2, then
    // (startIdx, startIdx+1), (startIdx+2, startIdx+3), ... over the T
    // thermodynamic states. Fills exchangePairList_.
    void prepareExchangePairs(int round, int oddity);
    // ReplicaMixingScheme::All (B7): draw nAttempts distinct random
    // thermo-state pairs and attemptREXSwap each.
    void mixAllReplicas(int nAttempts);
    // One round's exchange attempt (B7), gated by swapEvery_ and runType_.
    // Parity comes from the dedicated exchangeRound_ counter (incremented
    // once per EXECUTED mix), never from the raw round/mixi index (B7
    // revision-2 fix, reviewer Should-fix R5) -- see prepareExchangePairs.
    void mixReplicas(int mixi);
    // Diagnostic: print the attempted/accepted swap matrices to stderr
    // (I4: "a PrintNofAcceptedSwapsMatrix-style acceptance matrix becomes
    // available"). Called once at the end of RunREX.
    void printSwapMatrix() const;

    // ---- Stage 2b: driven (RENE/REBASONTOP) round -------------------------
    // One driven round (B6/B7/B8): (1) every replica's EQUILIBRIUM worlds run
    // (a Gibbs sweep exactly like mixReplicas' REMC path, but SKIPPING any
    // world whose distortOption is set -- B8's equilibrium/nonequilibrium
    // partition, read per-world via World::getDistortOption rather than a
    // precomputed N1_wCnt split index, since the current engine's schedule is
    // homogeneous across states, B0), refreshing potential/referencePotential
    // and accumulating the BAT anchor (INV-9) from each equilibrium visit;
    // (2) the neighbour pairing is prepared (B7) BEFORE any drive runs, since
    // the Q-scale-factor s = sqrt(T_target/T_source) (B4) needs the PARTNER's
    // temperature; (3) ONE frozen anchor snapshot (INV-9) is taken and passed
    // to every drive this round; (4) each paired replica is driven toward its
    // partner's temperature (driveReplica); (5) attemptREXSwap runs for every
    // pair; (6) REBASONTOP's interleaved REMC sub-rounds (D4) run every
    // interleaveRemcEvery_ rounds.
    void runDrivenRound(int round, bool verbose);
    // Drive replica `replicaIx` (currently at thermodynamic state `thermoIx`)
    // toward `targetTemperature` (B4: s = sqrt(targetTemperature /
    // thermodynamicStates_[thermoIx].temperature)). Resets WORK/WORK_Jacobian
    // (INV-5), then for every world position in that state's schedule with
    // distortOption == ScaleBendStretch (RunREX's INV-10 guard already
    // confirmed no OTHER-typed driven world exists for RENE/REBASONTOP),
    // applies World::applyBatScalingDrive with the SHARED frozen `anchor`,
    // accumulating WORK (B5, live += accumulation, Fixman excluded per D3)
    // and WORK_Jacobian (B12) across every driven world visited. A
    // std::domain_error from applyBatScalingDrive (Stage 2a's r1<=0/theta1
    // outside (0,pi) guard) is caught here and converted to a forced-reject
    // WORK_Jacobian = -infinity (so attemptREXSwap's Work_partner term goes
    // to +infinity and the swap always rejects) rather than propagating an
    // exception out of the round loop (reviewer N2 fail-loud, applied at the
    // swap level, not a crash).
    void driveReplica(int replicaIx,
                       int thermoIx,
                       double targetTemperature,
                       const robo::BatAnchorStats::Snapshot& anchor);
    // D4: REBASONTOP's interleaved REMC neighbour-swap sub-rounds. Runs
    // rebasontopSubrounds_ full (alternating-parity) REMC-style
    // (ETerm_equal, committed potentials only) sweeps by temporarily
    // borrowing attemptREXSwap's REMC branch (runType_ flipped and restored
    // around the loop) -- these touch only potential/referencePotential,
    // never WORK_*, so they compose safely with the driven main swap.
    void runInterleavedRemcSubround();

    std::string baseName;
    std::uint32_t seed = 0;

    std::vector<std::unique_ptr<World>> worlds_;         // heap so World& stays valid
    std::vector<double> temperatures_;                   // one per replica
    std::vector<std::vector<robo::Vec3>> replicaCoords_; // per replica, nm, OpenMM order
    int writeCounter_ = 0;
    std::vector<dcd::Writer> dcdWriters_; // one trajectory per replica (baseName.<r>.dcd)
    std::vector<double> dcdScratch_;      // interleaved xyz, nm->Angstrom, reused
    std::mt19937_64 rexRng_;
    std::uniform_real_distribution<double> rexUniform_{0.0, 1.0};

    // ---- label-swap replica exchange state (RunREX) -----------------------
    RUN_TYPE runType_ = RUN_TYPE::Default;
    std::vector<Replica> replicas_;                     // R persistent configurations (B1)
    std::vector<ThermodynamicState> thermodynamicStates_; // T temperatures + schedules (B1)
    // Two mutually-inverse maps (B0/B6): replica2ThermoIxs_[replicaIx] = the
    // thermodynamic state currently simulating that replica's coordinates;
    // thermo2ReplicaIxs_[thermoIx] = the replica currently occupying that
    // state. Both initialised to identity by setupReplicaExchange (B0).
    std::vector<int> replica2ThermoIxs_;
    std::vector<int> thermo2ReplicaIxs_;

    ReplicaMixingScheme mixingScheme_ = ReplicaMixingScheme::Neighboring;
    int swapEvery_ = 1;
    int nSwapAttempts_ = 1;
    bool swapFixman_ = false; // OFF-by-default diagnostic (D3); never enters acceptance.
    int exchangeRound_ = 0;   // B7 revision-2: parity source, incremented once per executed mix.
    std::vector<std::pair<int, int>> exchangePairList_;
    std::vector<std::vector<std::int64_t>> nofAttemptedSwapsMatrix_; // T x T, symmetric
    std::vector<std::vector<std::int64_t>> nofAcceptedSwapsMatrix_;  // T x T, symmetric

    // Shared/global BAT-scaling anchor (INV-9) -- ONE instance, not per-state.
    robo::BatAnchorStats batAnchorStats_;

    // D4 REBASONTOP interleave configuration (setInterleaveRemcEvery/
    // setRebasontopSubrounds).
    int interleaveRemcEvery_ = 10;
    int rebasontopSubrounds_ = 6;
};