#pragma once

/**
 * @file ReplicaExchange.hpp
 * @brief Replica-exchange value types: the run-type/topology enums and the
 *        `Replica` / `ThermodynamicState` records the label-swap driver owns.
 *
 * These are plain data carriers with no behavior beyond `commitWorkAsFinal`;
 * all sampling logic lives in Context (ReplicaExchangeDriver.cpp,
 * SwapAcceptance.cpp, DrivenRexDriver.cpp). The label-swap design keeps a
 * replica's configuration in place and swaps only its thermodynamic-state
 * label, so the two inverse index maps that realize the swap live on Context,
 * not on these types.
 *
 * @note Runtime-verified surface: only `RUN_TYPE::REMC` and `RUN_TYPE::Default`
 *       are exercised by a run. The RENE/REBASONTOP/RENEMC fields
 *       (`WORK_*`, `referenceWORK_potential`) and their acceptance paths are
 *       compiled-but-not-runtime-exercised (ARCHITECTURE OQ-5); their
 *       contracts below are stated from source, not from a run.
 * @see docs/specs/replica-exchange-nonequilibrium-work.md
 */

#include <cstdint>
#include <vector>

#include "World.hpp"     // AcceptRejectMode
#include "robot_math.hpp" // robo::Vec3

/**
 * @brief Selects the outer Markov chain over thermodynamic-state permutations:
 *        which acceptance rule (if any) a replica-exchange round applies.
 *
 * Passed to `Context::RunREX`. Chooses the acceptance branch in
 * `Context::attemptREXSwap` and the round structure in the driver.
 *
 * @note Runtime-exercised: `Default`, `REMC`. `RENE`/`REBASONTOP` are wired but
 *       uncompiled/untested (OQ-5). `RENEMC`'s acceptance formula is wired and
 *       unit-tested in isolation, but its driven round-loop is unimplemented:
 *       `Context::RunREX(RENEMC, ...)` throws `std::logic_error`.
 */
enum class RUN_TYPE : std::uint8_t {
    Default = 0, ///< Independent replicas; no exchange attempts.
    REMC,        ///< Label-swap parallel tempering, `ETerm_equal` acceptance (B6).
    RENEMC,      ///< `ETerm_nonequil` acceptance wired and tested; round-loop drive unimplemented (throws).
    RENE,        ///< Driven BAT-scaling exchange, work-based `WTerm` acceptance (uncompiled).
    REBASONTOP   ///< RENE work-swaps plus interleaved REMC sub-rounds (uncompiled).
};

/**
 * @brief Exchange topology: how `mixReplicas` chooses which state pairs to
 *        attempt each round.
 *
 * @note `Neighboring` (the default) attempts adjacent thermodynamic-state pairs
 *       with round-alternating parity; `All` draws random distinct pairs
 *       (`nSwapAttempts` draws).
 */
enum class ReplicaMixingScheme : std::uint8_t {
    All = 0,
    Neighboring = 1
};

/**
 * @brief One replica: a persistent molecular configuration, the energies
 *        measured on it, and (for driven run types) the nonequilibrium trial
 *        state held apart until accept/reject.
 *
 * Owned by Context in a `vector<Replica>` sized `R = temperatures.size()`. A
 * replica's `atomsLocations` is the only per-replica state that survives a
 * Gibbs sweep; Worlds hold none (INV-3). The committed block and the trial
 * block are kept separate so a driven swap can commit atomically on accept.
 *
 * @invariant INV-6: `potential` equals the OpenMM potential energy of
 *            `atomsLocations`. The driver refreshes it before every swap
 *            attempt.
 *
 * @note The `WORK_*` / `referenceWORK_*` trial fields are populated only during
 *       a driven round (RENE/REBASONTOP; RENEMC when its round-loop lands) and
 *       are uncompiled/untested (OQ-5). REMC/Default never touch them.
 */
class Replica {
    public:
    // ---- committed (equilibrium) state -----------------------------------
    /// Committed Cartesian coordinates (nm, engine/OpenMM atom order). The
    /// configuration whose label a swap moves.
    std::vector<robo::Vec3> atomsLocations;

    /// OpenMM potential energy of @ref atomsLocations (INV-6). @ref potential
    /// and @ref referencePotential coincide for every run type: Fixman is
    /// excluded from both unconditionally (D3). The pair exists for interface
    /// parity with the original two-potential design, not because they diverge.
    double potential = 0.0;
    double referencePotential = 0.0; ///< @see potential
    /// Committed Fixman potential, diagnostic only: computed by the sampler but
    /// never added into @ref potential / @ref referencePotential or any
    /// acceptance exponent (D3/INV-7). Tracked so the atomic commit has all
    /// four committed quantities to promote together.
    double FixmanPotential = 0.0;

    // ---- nonequilibrium trial state (driven run types only) --------------
    std::vector<robo::Vec3> WORK_atomsLocations; ///< Driven endpoint x^tau (== x' when mdSteps==0).
    double WORK_potential = 0.0;                 ///< U(x^tau), unreduced physical PE.
    double referenceWORK_potential = 0.0;        ///< == @ref WORK_potential (Fixman excluded, D3).
    double WORK_FixmanPotential = 0.0;           ///< Trial Fixman (diagnostic only, D3).
    /// Accumulated nonequilibrium work over the driven range: sum of per-driven
    /// -world (U_curr - U_prev), Fixman excluded (D3/B5). Reset once per driven
    /// range, accumulated once per driven World visited.
    double WORK = 0.0;
    /// Accumulated log-Jacobian over the driven range: sum of per-driven-world
    /// `getDistortJacobianDetLog()` (D6). A forced-reject sentinel of
    /// `-infinity` set by `driveReplica` on a domain-invalid drive.
    double WORK_Jacobian = 0.0;

    /**
     * @brief Promote the trial block to committed atomically (all four
     *        quantities: coordinates, potential, referencePotential,
     *        FixmanPotential).
     *
     * @pre Called only on an accepted driven swap, after the `WORK_*` trial has
     *      been populated by `Context::driveReplica`. Never called on reject
     *      (the committed state already stands) and never for REMC/Default
     *      (their accept is a label swap only).
     * @post `atomsLocations`, `potential`, `referencePotential`,
     *       `FixmanPotential` equal their `WORK_*` counterparts. No partial
     *       state is ever observable.
     * @note Part of the uncompiled driven path (OQ-5).
     */
    void commitWorkAsFinal() {
        atomsLocations = WORK_atomsLocations;
        potential = WORK_potential;
        referencePotential = referenceWORK_potential;
        FixmanPotential = WORK_FixmanPotential;
    }
};

/**
 * @brief One thermodynamic state: a target temperature plus the per-World
 *        simulation schedule a replica runs while it occupies this state.
 *
 * Defines the target Boltzmann distribution (via @ref temperature) for whichever
 * replica the driver's inverse index maps currently place here. Owned by
 * Context in a `vector<ThermodynamicState>` sized `T`, canonically `T == R`
 * (Context asserts this at setup, never assumes it).
 *
 * @note The schedule fields let the driver reset each shared World's
 *       temperature / timestep / MD-step-count / accept-mode every round from
 *       the state currently occupying it. Populated from the Worlds' own
 *       configuration at setup and identical across states for the runtime-
 *       exercised REMC/Default path; they diverge only on the uncompiled driven
 *       path (OQ-5).
 */
class ThermodynamicState {
    public:
    double temperature = 300.0; ///< Target temperature (K); fixes this state's Boltzmann distribution.

    /// World visitation order for a Gibbs sweep at this state; entry `pos` names
    /// the World index run at schedule position `pos`. Held explicit (not
    /// implicit 0..W-1) so a future permutation/subset schedule needs no
    /// interface change; currently `[0..W-1]` for every state.
    std::vector<int> worldIndexes;

    /// Per-schedule-position runtime overrides, one entry per @ref worldIndexes
    /// position, applied to the shared World before it runs at this state.
    std::vector<double> timeSteps;
    std::vector<int> mdSteps;                        ///< @see timeSteps
    std::vector<AcceptRejectMode> acceptRejectModes; ///< @see timeSteps
};
