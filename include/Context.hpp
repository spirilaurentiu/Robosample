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

/**
 * @brief Robosample top-level run orchestrator: owns the Worlds, the replica
 *        temperature ladder, and the replica-exchange driver, and delegates a
 *        run to it.
 *
 * The engine is organized around replica exchange; a single-replica run is the
 * degenerate `R = 1` case of the default REMC driver. World order defines the
 * Gibbs sweep order, and per-atom Ground-frame coordinates (nm) are the only
 * currency passed between Worlds (INV-3).
 *
 * Lifecycle (the required call order): construct -> `add*World(...)` (one or
 * more) and `add_sampler` on each -> optional pre-`initialize` configuration
 * (`setMTS`, `setSeparateForceGroups`, `setEnforcePeriodicBox`, mixing setters)
 * -> `initialize(temperatures)` -> a run entry (`RunREX` or the legacy
 * `runREX`). Adding a World after `initialize` is unsupported; a run entry
 * before `initialize` throws.
 *
 * @note Constructed from Python as `Context(base_name, seed)`; the Python
 *       subclass (context.py) supplies the dihedral classifier and `load_amber`,
 *       which fills @ref systemTopology in place before any World is added.
 * @see docs/architecture, MODULES.md (workflow layer).
 */
class Context {
    public:
    /**
     * @brief Construct an empty context with a run base name and master seed.
     * @param[in] baseName Output-file stem; per-replica CSV/DCD paths are
     *                     `baseName.<idx>.{csv,dcd,reactions.csv}` and
     *                     `baseName.moves.csv`. Taken by value (moved in).
     * @param[in] seed     Master RNG seed. The REX driver's RNG is derived from
     *                     it; each World is seeded from it at `add*World`.
     * @post @ref systemTopology is default-constructed and empty; no Worlds
     *       exist. The caller (load_amber) fills @ref systemTopology next.
     */
    Context(std::string baseName, std::uint32_t seed);

    /// The molecular system definition (topology, coordinates, box, bonded
    /// sets). Bound to Python as `system_topology` and filled in place by
    /// `load_amber` before any World is added. A public member because it is
    /// the Python<->engine payload; every `add*World` reads it to build a model.
    SystemTopology systemTopology;

    // ---- modeling --------------------------------------------------------
    // Root mobility is a PER-WORLD property. systemTopology.rootMobilities is
    // only the build-time DEFAULT each world seeds from (set by load_amber);
    // there is deliberately no Context-level mutator. To change an individual
    // molecule's root attachment, call World::setRootMobility(...) on the world
    // returned by add*World (see World.hpp), which rebuilds that world's model
    // in isolation -- the same mechanism addDockingWorld already uses.

    /**
     * @brief Append a Cartesian World (pure on-device OpenMM MD) and build its
     *        model immediately from the current @ref systemTopology.
     * @param[in] wantReactionReporter Opt-in per-body applied-force reporter.
     *            @b Always throws when `true`: a Cartesian World has no
     *            meaningful articulated-body indexing for the reporter.
     * @return Reference to the newly added World, heap-owned by this Context so
     *         the reference stays valid across sweeps. Chain `.add_sampler(...)`
     *         on it.
     * @pre Call before `initialize`; @ref systemTopology already filled.
     * @post One World appended; World order (and thus Gibbs sweep order) is
     *       append order.
     */
    World& addCartesianWorld(bool wantReactionReporter = false);
    /**
     * @brief Append a robotic/torsional World (internal-coordinate HMC) with the
     *        given per-bond flexibility selection, building its model at once.
     * @param[in] sel Per-bond mobility (see @ref buildFlexibilities). Borrowed;
     *            copied into the built model.
     * @param[in] wantReactionReporter Opt-in per-body applied-force reporter.
     *            When `true`, flags every flexed (non-Weld) body plus its parent,
     *            excluding only Ground; a Free-rooted body (attached directly to
     *            Ground) is included.
     * @return Reference to the newly added World (heap-owned, stable across
     *         sweeps). Chain `.add_sampler(...)`.
     * @pre Call before `initialize`.
     */
    World& addRoboticWorld(const Selection& sel, bool wantReactionReporter = false);

    /**
     * @brief Append a rigid-body docking World: listed molecules become mobile
     *        6-DOF ligands, every other molecule is welded rigid to Ground.
     *
     * Root mobility here is a World property built from @p ligandMoleculeIndices,
     * not from `systemTopology.rootMobilities`, so a docking World and ordinary
     * Worlds can coexist over one system. Each ligand becomes a Free root with
     * all its bonds Rigid; the binding-site centre is the centroid of all
     * non-ligand (receptor) atoms; pass the sphere radius via `add_sampler(...)`.
     *
     * @param[in] ligandMoleculeIndices Molecule indices to treat as ligands.
     * @return Reference to the new World (heap-owned, stable). Chain
     *         `.add_sampler(...)`.
     * @throws std::out_of_range if any index is outside `[0, numMolecules)`.
     * @pre Call before `initialize`.
     */
    World& addDockingWorld(const std::vector<int>& ligandMoleculeIndices);


    /**
     * @brief Build a per-bond mobility selection for `addRoboticWorld`.
     * @param[in] bonds    If empty/None, every eligible bond gets @p mobility;
     *            otherwise only the listed `(i,j)` bonds do. A bond is eligible
     *            when it is non-ring-closing and both atoms have >= 2 bonds
     *            (rotatable). Order within a pair is irrelevant.
     * @param[in] mobility Joint type assigned to selected bonds; all others stay
     *            Rigid.
     * @param[in] flag     Currently unused.
     * @return A Selection with per-bond mobility; ring-closing bonds always
     *         remain Rigid regardless of @p bonds.
     */
    Selection buildFlexibilities(const std::optional<std::vector<std::pair<int, int>>>& bonds,
                                 JointType mobility,
                                 bool flag);

    // ---- run -------------------------------------------------------------
    /**
     * @brief Bring up OpenMM, run the startup geometry check, set the replica
     *        temperature ladder, and seed every replica with the reference
     *        coordinates.
     *
     * The transition from configuration to runnable state. Must be called after
     * every `add*World` and before any run entry.
     *
     * @param[in] temperatures Per-replica temperature ladder (K); an empty list
     *            means a single replica at 300 K. Its size fixes the replica
     *            count `R` for the run.
     * @throws std::runtime_error if OpenMM initialization fails, or (via
     *         `checkStartupGeometry`) if the input geometry is non-finite or
     *         sterically clashing and `ROBO_ALLOW_BAD_START` is not set.
     * @post `R = temperatures.size()` replicas seeded from `systemTopology`'s
     *       coordinates; per-replica CSV outputs truncated and DCD writers
     *       (re)created; `writeCounter` reset. Adding a World after this is
     *       unsupported.
     */
    void initialize(const std::vector<double>& temperatures);

    /**
     * @brief Legacy coordinate-swap replica-exchange driver, retained frozen as
     *        the differential oracle for the label-swap `RunREX`.
     *
     * Runs `equilRounds` then `prodRounds` rounds; each round is a Gibbs sweep
     * over all Worlds for every replica, followed by alternating-parity
     * adjacent-replica exchange attempts, then outputs every `writeFreq`
     * production rounds. On an accepted exchange it swaps the two replicas'
     * coordinate buffers directly (it moves configurations, not labels), so it
     * is INV-3 non-compliant by construction. The exchange accept test is the
     * REMC criterion `(beta_a - beta_b)(E_a - E_b) >= 0` or `u < exp(delta)`.
     *
     * Its role is the contract: `runREX` and `RunREX` produce equivalent
     * sampling (INV-8), and `runREX` @b is the equivalent-sampling reference
     * `RunREX` is validated against (Python `test_rex_label_swap_equivalence`).
     *
     * @param[in] equilRounds Equilibration rounds (all Worlds AlwaysAccept).
     * @param[in] prodRounds  Production rounds (each World's configured accept
     *                        mode).
     * @param[in] writeFreq   Output cadence in rounds; `<= 0` disables periodic
     *                        output. Also gates the per-move telemetry CSV.
     * @param[in] verbose     Emit per-World stdout telemetry each move.
     * @pre `initialize` has run.
     * @warning Frozen oracle: do @b not extend this method. New run types and
     *          the label-swap object model live in `RunREX`; extending `runREX`
     *          would break its role as the differential baseline. Driver scripts
     *          historically call it at a single replica; that is context, not a
     *          recommendation.
     */
    void runREX(int equilRounds, int prodRounds, int writeFreq, bool verbose);

    /**
     * @brief Label-swap replica-exchange driver: the default center of the
     *        engine. Runs `R` replicas on the temperature ladder and exchanges
     *        their thermodynamic-state labels, leaving configurations in place.
     *
     * Builds `R = temperatures.size()` Replica/ThermodynamicState objects
     * (asserting `R == T`) from the current Worlds' schedule and the coordinates
     * `initialize` seeded, then runs `equilRounds + prodRounds` rounds. Each
     * round propagates every replica through its state's World schedule (a Gibbs
     * sweep, iterated by thermodynamic-state index so the shared-World RNG is
     * consumed in the same order as the legacy `runREX`), refreshes committed
     * potentials (INV-6), then attempts a round of neighbour swaps via
     * `mixReplicas` -> `attemptREXSwap`, which swaps @b labels only (INV-3),
     * never coordinates. Output CSV/DCD are indexed by thermodynamic state (a
     * fixed output slot is a fixed temperature), matching `runREX` so the two
     * drivers' per-state statistics are directly comparable.
     *
     * Contract (INV-8): `RunREX` and `runREX` produce equivalent sampling; the
     * two index maps stay mutually inverse and the swap matrix keeps its
     * alternating-parity structure. The label-swap/coordinate-swap equivalence
     * is validated (4-replica REMC on CUDA; Python equivalence test).
     *
     * @param[in] runType     Acceptance rule / round structure. Runtime-verified:
     *            `Default`, `REMC`. `RENE`/`REBASONTOP` dispatch to the
     *            uncompiled driven round (OQ-5). `RENEMC` throws
     *            (`std::logic_error`): its acceptance formula is wired and
     *            unit-tested but its driven round-loop is unimplemented.
     * @param[in] equilRounds  Equilibration rounds.
     * @param[in] prodRounds   Production rounds.
     * @param[in] writeFreq    Output cadence in rounds; `<= 0` disables output.
     * @param[in] verbose      Emit per-move stdout telemetry.
     * @throws std::runtime_error if `initialize` has not seeded replicas.
     * @throws std::logic_error for `RUN_TYPE::RENEMC`, or (driven types) if the
     *         INV-7/INV-10 guards fail before the round loop.
     * @note This is @b not established as a proven stationary-distribution
     *       guarantee: the end-to-end detailed-balance / stationary-distribution
     *       oracle is a pending separate track. Today's evidence is the
     *       acceptance algebra (`attemptREXSwap`, tested) plus the label-swap
     *       equivalence to `runREX` (INV-8). Cite it as equivalence-to-`runREX`,
     *       not as a proven per-replica Boltzmann guarantee.
     */
    void RunREX(RUN_TYPE runType, int equilRounds, int prodRounds, int writeFreq, bool verbose);

    /**
     * @brief Attempt one exchange between thermodynamic states @p thermoC and
     *        @p thermoH under the active run type; on accept, swap labels only.
     *
     * The theory-derived acceptance algebra. Increments the attempted-swap
     * matrix, computes a log-acceptance exponent per run type, then accepts iff
     * the exponent is non-negative or `u < exp(exponent)`. On accept it
     * increments the accepted-swap matrix and swaps the two states' labels
     * (INV-3); it never moves coordinates.
     *
     * Per-run-type exponent (`beta = 1/(k_B T)`; `X`,`Y` are the replicas at
     * `thermoC`,`thermoH`):
     * - REMC: `ETerm_equal = -(beta_H - beta_C)(refU_X - refU_Y)`.
     * - RENEMC: `ETerm_nonequil`, the same form on the driven-endpoint reference
     *   potentials, no Jacobian (INV-10: volume-preserving drive).
     * - RENE / REBASONTOP: `WTerm = -(Work_X + Work_Y)` where
     *   `Work_p = beta_target*U(x_p^tau) - beta_source*U(x_p^0) - lnJac_p`
     *   (Ballard-Jarzynski / Nilmeier deterministic-map acceptance). On accept,
     *   both replicas' trial state is committed atomically
     *   (`Replica::commitWorkAsFinal`) @b before the label swap; REMC/Default
     *   never commit.
     *
     * @param[in] thermoC Cold-side thermodynamic-state index.
     * @param[in] thermoH Hot-side thermodynamic-state index.
     * @return `true` iff the swap was accepted.
     * @throws std::logic_error if no thermodynamic states exist (called before a
     *         run built them) or if the active run type is `RUN_TYPE::Default`.
     * @warning A non-finite acceptance exponent (NaN or +-inf, e.g. a
     *          domain-invalid driven endpoint whose sentinel Jacobian is -inf,
     *          or a blown-up OpenMM PE) forces an explicit automatic reject
     *          logged to stderr, never a silent NaN comparison (INV-8 guard).
     * @note Public (not private) so a reproducer can drive it directly; the
     *       REMC exponent (INV-8 detailed balance), Jacobian sign-flip
     *       blindness of the symmetric paired swap, RENEMC `ETerm_nonequil`, and
     *       the domain-error sentinel are pinned by
     *       tests/TestRexAcceptanceAlgebra.cpp. The RENEMC/RENE/REBASONTOP
     *       branches are otherwise uncompiled at run scope (OQ-5).
     */
    bool attemptREXSwap(int thermoC, int thermoH);

    /**
     * @brief Enforce the driven-run preconditions INV-7 (Fixman enabled in every
     *        non-Cartesian sampler) and INV-10 (each driven World's distortion
     *        matches the run type), or throw.
     *
     * A no-op for `REMC`/`Default`; for `RENE`/`RENEMC`/`REBASONTOP` it requires
     * that (a) every non-Cartesian World has Fixman enabled, because the swap
     * acceptance excludes Fixman correctly only then (D3 biconditional); (b)
     * every driven World's `distortOption` is the type the run demands
     * (`ScaleBendStretch` for RENE/REBASONTOP, `NMA` for RENEMC); and (c) at
     * least one driven World exists, so the drive is not silently inert.
     *
     * @param[in] runType Run type to check against, passed explicitly (not the
     *            `runType_` member) so it is callable in isolation.
     * @throws std::logic_error on any INV-7 or INV-10 violation, or if a driven
     *         run type has no matching driven World.
     * @note Reads only the Worlds (no replicas/states needed), so a reproducer
     *       can call it on a hand-built Context without `initialize`/`RunREX`;
     *       tests/TestRexAcceptanceAlgebra.cpp's INV-10 case does. The guards
     *       themselves are unexercised at run scope (OQ-5: driven runs never run).
     */
    void checkInv7AndInv10Guards(RUN_TYPE runType) const;

    /**
     * @brief Configure replica-exchange mixing, applied by the next `RunREX`.
     *
     * - `setReplicaMixingScheme`: `Neighboring` (default, alternating-parity
     *   adjacent pairs) or `All` (random pairs).
     * - `setSwapEvery(n)`: attempt a mix only when `round % n == 0` (clamped to
     *   `>= 1`).
     * - `setNSwapAttempts(n)`: draw count for `ReplicaMixingScheme::All`
     *   (clamped to `>= 1`).
     * - `setSwapFixman(enabled)`: off-by-default diagnostic only. Fixman is
     *   computed by the sampler but never enters `attemptREXSwap`'s acceptance
     *   exponent (INV-7); this flag exists for interface parity and does not
     *   affect REMC acceptance.
     */
    void setReplicaMixingScheme(ReplicaMixingScheme scheme) {
        mixingScheme_ = scheme;
    }
    /// Attempt a mix only every @p n th round (clamped `>= 1`). @see setReplicaMixingScheme
    void setSwapEvery(int n) {
        swapEvery_ = (n > 0) ? n : 1;
    }
    /// Draw count for `ReplicaMixingScheme::All` (clamped `>= 1`). @see setReplicaMixingScheme
    void setNSwapAttempts(int n) {
        nSwapAttempts_ = (n > 0) ? n : 1;
    }
    /// Off-by-default diagnostic; never enters acceptance (INV-7). @see setReplicaMixingScheme
    void setSwapFixman(bool enabled) {
        swapFixman_ = enabled;
    }

    /**
     * @brief Configure REBASONTOP's interleaved REMC sub-rounds.
     *
     * Every `interleaveRemcEvery` driven rounds, `rebasontopSubrounds`
     * REMC-style neighbour-swap sub-rounds (ETerm_equal, committed potentials
     * only, no trial state touched) run in addition to that round's main WTerm
     * swap. Both clamp to `>= 1`. Defaults are (10, 6).
     *
     * @note Only meaningful under `RUN_TYPE::REBASONTOP`, which is uncompiled at
     *       run scope (OQ-5).
     */
    void setInterleaveRemcEvery(int n) {
        interleaveRemcEvery_ = (n > 0) ? n : 1;
    }
    /// Number of interleaved REMC sub-rounds per interleave (clamped `>= 1`).
    /// @see setInterleaveRemcEvery
    void setRebasontopSubrounds(int n) {
        rebasontopSubrounds_ = (n > 0) ? n : 1;
    }

    /**
     * @brief Symmetric per-state-pair swap tallies accumulated across a `RunREX`.
     *
     * Both matrices are `T x T` and symmetric (`m[i][j] == m[j][i]`), indexed by
     * thermodynamic-state index, populated by `attemptREXSwap`. The attempted
     * count and the accepted/attempted ratio are the acceptance-rate diagnostic.
     * @return Borrowed references valid until the next `RunREX`/`initialize`.
     */
    [[nodiscard]] const std::vector<std::vector<std::int64_t>>& attemptedSwapsMatrix() const {
        return nofAttemptedSwapsMatrix_;
    }
    /// @copydoc attemptedSwapsMatrix
    [[nodiscard]] const std::vector<std::vector<std::int64_t>>& acceptedSwapsMatrix() const {
        return nofAcceptedSwapsMatrix_;
    }

    // ---- BAT-scaling shared/global anchor (driven run types) -------------
    /**
     * @brief Fold one World's current committed geometry into the single
     *        running-mean BAT anchor the driven drive scales deviations around.
     *
     * Context owns one anchor, shared across all thermodynamic states, so it is
     * state-independent by construction -- the precondition the paired scaling
     * map needs to be an exact involution (INV-9).
     *
     * @param[in] world World whose current Ground geometry is accumulated.
     *            Borrowed, read-only.
     * @pre Call only after an @b equilibrium move (`distortOption == nullopt`),
     *      never on a driven World's output, which would bias the anchor with
     *      nonequilibrium samples.
     * @note Part of the uncompiled driven path (OQ-5).
     */
    void accumulateBatAnchorStats(const World& world) {
        const robo::Vec3* p = world.getAtomsLocationsInGround();
        const std::vector<robo::Vec3> pos(p, p + world.model().numAtoms);
        batAnchorStats_.accumulate(world.model(), pos);
    }
    /**
     * @brief Freeze the current anchor into an immutable snapshot for one
     *        exchange round.
     * @return A value snapshot of the running means. Take exactly one per round
     *         and pass the same snapshot to both partners of every swap pair
     *         (INV-9). @note Uncompiled driven path (OQ-5).
     */
    [[nodiscard]] robo::BatAnchorStats::Snapshot batAnchorSnapshot() const {
        return batAnchorStats_.snapshot();
    }
    /// Clear the accumulated BAT-anchor running means. @note Uncompiled (OQ-5).
    void resetBatAnchorStats() {
        batAnchorStats_.reset();
    }

    /**
     * @brief Enable multiple-timestep (r-RESPA) integration for the Cartesian
     *        World's on-device OpenMM MD.
     * @param[in] enabled       Turn MTS on/off.
     * @param[in] innerSubsteps Fast bonded-force evaluations per outer step; slow
     *            forces (Nonbonded, GBSA, ...) are evaluated once per outer step.
     * @pre Call before `initialize`. No effect on torsional Worlds (they always
     *      use the full force sum). Forwards to the OpenMM singleton.
     */
    void setMTS(bool enabled, int innerSubsteps);

    /**
     * @brief Enable/disable separate OpenMM force groups per Force.
     * @param[in] enabled Whether each Force gets its own group (needed for
     *            per-group energy decomposition).
     * @note Thin forwarder onto the OpenMM singleton, present so callers hold
     *       only a `Context` and never reach the singleton directly.
     */
    void setSeparateForceGroups(bool enabled);

    /**
     * @brief Set whether OpenMM wraps coordinates into the primary box when
     *        state is pulled back to the engine.
     * @param[in] enabled Must stay `false` (the default) under explicit solvent
     *            so the robot engine receives whole molecules; energies/forces
     *            are unaffected (minimum image is always applied internally).
     * @note Thin forwarder onto the OpenMM singleton (same rationale as
     *       `setSeparateForceGroups`).
     */
    void setEnforcePeriodicBox(bool enabled);


    // ---- energy ingestion / validation ----------------------------------
    /// Bring up the OpenMM singleton from @ref systemTopology.
    /// @return `true` on success. Called by `initialize`.
    auto initializeOpenMM() -> bool;
    /// @return OpenMM potential energy (kJ/mol) of @ref systemTopology's current
    /// reference coordinates.
    [[nodiscard]] auto calcOpenMMPotentialEnergy() -> double;
    /// @return Total potential energy and its per-force-group decomposition for
    /// the reference coordinates. @pre `setSeparateForceGroups(true)` for a
    /// meaningful split.
    [[nodiscard]] auto computePotentialEnergyByGroup()
        -> std::pair<double, std::vector<OpenMMContext::ForceGroupEnergy>>;

    /// @return The run output-file base name.
    [[nodiscard]] auto getBaseName() const -> const std::string& {
        return baseName;
    }
    /// @return The master RNG seed.
    [[nodiscard]] auto getSeed() const -> std::uint32_t {
        return seed;
    }
    /// @return The number of Worlds added so far.
    [[nodiscard]] int numWorlds() const {
        return static_cast<int>(worlds_.size());
    }

    private:
    /// @return OpenMM potential energy (kJ/mol) of an arbitrary coordinate set.
    /// @param[in] coords Per-atom Ground coordinates (nm, engine order), borrowed.
    double openmmPotential(const std::vector<robo::Vec3>& coords);
    /// Write one replica's outputs for @ref runREX. Thin wrapper over
    /// @ref writeOutputsCore using `replicaCoords_[replica]` and
    /// `temperatures_[replica]`.
    void writeOutputs(int replica, int round, bool verbose);
    /**
     * @brief Emit one frame for output slot @p idx: append an energy row to the
     *        per-slot CSV and one imaged frame to the per-slot DCD.
     *
     * The single output core shared by both drivers. Recomputes the OpenMM PE of
     * @p coords, appends `round,idx,T,PE` to `baseName.<idx>.csv`, and (if a DCD
     * writer exists for @p idx) scatters @p coords into the DCD scratch buffer in
     * prmtop atom order, converts nm->Angstrom, and appends a frame.
     *
     * Whole-molecule periodic imaging is applied @b here and only to the DCD
     * copy: each molecule is shifted by integer lattice vectors so its mass
     * -weighted COM lands in the primary cell, then translated rigidly (bonds
     * never straddle a face). The imaged coordinates are output-only; the
     * sampled `coords`/`replicaCoords_` stay unwrapped and never re-enter
     * sampling. Non-periodic systems (or missing molecule ranges) write
     * coordinates verbatim.
     *
     * @param[in] idx     Output-file slot and the row's "replica" column. For
     *            `runREX` the replica-object index (== its fixed temperature
     *            slot); for `RunREX` the thermodynamic-state index -- either way
     *            a fixed slot is a fixed temperature.
     * @param[in] round   Round number written into the CSV row.
     * @param[in] verbose Also print the energy line to stdout.
     * @param[in] coords  Coordinates to image and write (nm, engine order),
     *            borrowed and not modified.
     * @param[in] T       Temperature written into the CSV row.
     */
    void writeOutputsCore(int idx, int round, bool verbose, const std::vector<robo::Vec3>& coords, double T);
    /**
     * @brief Append a reporter World's captured per-body reaction-force rows to
     *        that replica's `baseName.<replica>.reactions.csv`.
     *
     * Writes the header only when the file is empty, then one row per sample.
     * The schema is always the fixed 10 columns
     * `frame,replica,body_idx,atom_idx,fx,fy,fz,tx,ty,tz` regardless of which
     * force term(s) the reporter World summed.
     *
     * @param[in] replica Replica index; the CSV slot and the row's replica column.
     * @param[in] frame   DCD frame index these rows pair with (for CSV/DCD
     *            alignment).
     * @param[in] rows    Captured samples, borrowed. Empty => no-op (no empty
     *            header-only file is produced).
     */
    void writeReactionRows(int replica, int frame, const std::vector<ReactionSample>& rows);

    /**
     * @brief Startup geometry sanity scan: refuse to start (or warn) when the
     *        input structure is unusable.
     *
     * A precondition check run once inside `initialize`, before any sampling. It
     * never mutates sampled state. It flags three conditions: non-finite
     * coordinates; hard steric clashes (a non-excluded, non-virtual atom pair
     * closer than a fixed threshold, under minimum-image distance when the box
     * is periodic); and a non-finite or pathologically positive initial
     * potential energy. Bonded/excluded pairs (1-2/1-3/1-4) and massless virtual
     * sites are excluded from the clash scan. It always prints a one-line
     * summary to stderr.
     *
     * @throws std::runtime_error describing the failure (NaN, clashes, and/or
     *         bad PE) when any condition trips and `ROBO_ALLOW_BAD_START` is
     *         unset or `"0"`. When that env var is set to anything else, the
     *         same message is logged as a warning and the run continues.
     * @note A healthy structure (finite, no hard clashes, sane PE) returns
     *       silently after the summary line. Rationale: a docking World welds
     *       the receptor rigid and cannot relax a frozen clash, so an unchecked
     *       bad start would loop forever rejecting kicks.
     */
    void checkStartupGeometry();

    // ---- label-swap replica exchange (RunREX) private helpers ------------
    /**
     * @brief (Re)build the replica/state objects, the identity index maps, and
     *        the zeroed swap matrices from `temperatures_` and the current World
     *        schedule.
     *
     * Seeds each replica from the coordinates `initialize` staged, evaluates each
     * committed potential (INV-6), and applies the pre-run docking pose-repair
     * kick (a no-op for non-docking systems). Sets both inverse index maps to
     * identity.
     * @param[in] runType Stored as the active run type for the run.
     * @throws std::runtime_error if no replicas exist (`initialize` not called).
     */
    void setupReplicaExchange(RUN_TYPE runType);
    /**
     * @brief Exchange the labels of thermodynamic states @p thermoC and
     *        @p thermoH, keeping both inverse maps consistent (INV-3).
     * @post `replica2ThermoIxs_` and `thermo2ReplicaIxs_` remain mutually
     *       inverse; no coordinates move.
     */
    void swapThermodynamicStates(int thermoC, int thermoH);
    /**
     * @brief Fill `exchangePairList_` with the neighbouring-pairs schedule for a
     *        round: `(startIdx,startIdx+1),(startIdx+2,startIdx+3),...` where
     *        `startIdx = (round + oddity) % 2`.
     * @param[in] round  Parity source (the driver passes `exchangeRound_`).
     * @param[in] oddity Parity offset (0, or `sub%2` for interleaved sub-rounds).
     */
    void prepareExchangePairs(int round, int oddity);
    /// Draw @p nAttempts distinct random thermodynamic-state pairs and attempt a
    /// swap on each (`ReplicaMixingScheme::All`). No-op when `T <= 1`.
    void mixAllReplicas(int nAttempts);
    /**
     * @brief Run one round's exchange attempts, gated by `swapEvery_` and the
     *        run type.
     *
     * No-op when the mix is not due, the run type is `Default`, or `T <= 1`.
     * Otherwise attempts the `Neighboring` schedule (parity from the dedicated
     * `exchangeRound_` counter, incremented once per executed mix, so both
     * parities stay reachable for any `swapEvery_`) or `mixAllReplicas`.
     * @param[in] mixi Current round index, tested against `swapEvery_`.
     */
    void mixReplicas(int mixi);
    /// Print the accepted/attempted swap matrices (by thermodynamic-state index)
    /// to stderr. Called once at the end of `RunREX`.
    void printSwapMatrix() const;

    // ---- driven (RENE/REBASONTOP) round -----------------------------------
    /**
     * @brief Run one driven exchange round: equilibrium sweep, pairing, BAT
     *        drive, work-based swap, then REBASONTOP's optional interleave.
     *
     * Intended structure: (1) sweep every replica's non-driven Worlds (skipping
     * any World with a `distortOption`), refreshing committed potentials (INV-6)
     * and feeding the BAT anchor from each equilibrium visit; (2) prepare the
     * neighbour pairing before any drive, since the Q-scale factor
     * `s = sqrt(T_target/T_source)` needs the partner's temperature; (3) take one
     * frozen anchor snapshot and share it across every drive this round (INV-9);
     * (4) drive each paired replica toward its partner's temperature; (5) attempt
     * every pair's swap; (6) for REBASONTOP, run interleaved REMC sub-rounds
     * every `interleaveRemcEvery_` rounds.
     *
     * @param[in] round   Round index (parity and interleave cadence).
     * @param[in] verbose Emit per-move stdout telemetry.
     * @note Assumed: this method and the whole driven path are documented from
     *       source only; they are compiled-but-not-runtime-exercised (OQ-5) and
     *       have no test. The structure above is a reviewed-on-paper design, not
     *       build-confirmed behavior. Recorded as an OPEN-QUESTION in findings.
     */
    void runDrivenRound(int round, bool verbose);
    /**
     * @brief Drive one replica from its current state toward @p targetTemperature
     *        via BAT scaling, accumulating nonequilibrium work and Jacobian into
     *        its trial block.
     *
     * Intended behavior: computes `s = sqrt(targetTemperature / T_source)`,
     * resets the trial work fields (INV-5), reseeds the trial coordinates from
     * the committed endpoint, then for each `ScaleBendStretch` World in the
     * schedule applies the scaling drive with the shared frozen @p anchor,
     * accumulating `WORK` (Fixman excluded) and `WORK_Jacobian`. A
     * `std::domain_error` from the drive (an invalid scaled geometry) is caught
     * and converted to a forced-reject sentinel `WORK_Jacobian = -infinity` so
     * `attemptREXSwap` rejects the swap, rather than propagating out of the round.
     *
     * @param[in] replicaIx        Replica to drive.
     * @param[in] thermoIx         Its current thermodynamic state (source temp).
     * @param[in] targetTemperature Partner's temperature (drive target).
     * @param[in] anchor           Shared frozen BAT anchor for this round (INV-9),
     *            borrowed.
     * @note Assumed: uncompiled/untested driven path (OQ-5). The consuming side
     *       of the domain-error sentinel is unit-checked in
     *       tests/TestRexAcceptanceAlgebra.cpp, but this try/catch itself is not
     *       exercised end-to-end. Recorded in findings.
     */
    void driveReplica(int replicaIx,
                       int thermoIx,
                       double targetTemperature,
                       const robo::BatAnchorStats::Snapshot& anchor);
    /**
     * @brief Run REBASONTOP's interleaved REMC neighbour-swap sub-rounds.
     *
     * Intended behavior: runs `rebasontopSubrounds_` alternating-parity
     * REMC-style sweeps by temporarily borrowing `attemptREXSwap`'s REMC branch
     * (the run type is flipped and restored around the loop); these touch only
     * committed potentials, never the trial work fields, so they compose with the
     * driven main swap.
     * @note Assumed: uncompiled/untested driven path (OQ-5). Recorded in findings.
     */
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
    /// Two mutually-inverse permutation maps realizing the label swap.
    /// `replica2ThermoIxs_[replicaIx]` is the thermodynamic state currently
    /// simulating that replica's coordinates; `thermo2ReplicaIxs_[thermoIx]` is
    /// the replica currently occupying that state.
    /// @invariant They are inverses at all times:
    ///            `thermo2ReplicaIxs_[replica2ThermoIxs_[r]] == r`. Identity at
    ///            setup; every `swapThermodynamicStates` preserves the property.
    std::vector<int> replica2ThermoIxs_;
    std::vector<int> thermo2ReplicaIxs_; ///< @see replica2ThermoIxs_

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