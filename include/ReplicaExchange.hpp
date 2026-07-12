#pragma once

// ============================================================================
//  Replica-exchange object model (docs/specs/replica-exchange-nonequilibrium-
//  work.md). Stage 1: REMC (parallel-tempering) label-swap. Stage 2b/2c:
//  RENE/REBASONTOP driven exchange (B5/B6, WORK_* fields + F4 atomic commit)
//  and RENEMC's acceptance algebra (B6/INV-10) are implemented; RENEMC's OWN
//  velocity/NMA driven-round-loop wiring is a Stage 2c TODO (see
//  Context::RunREX) -- it reuses the ALREADY-EXISTING DistortOption::NMA
//  momentum-draw distortion (World::reinitialize), not a new mechanism, but
//  the round-loop partitioning + INV-10 pairing enforcement for it is not
//  wired end-to-end yet. NONE of the Stage 2b/2c code in this header/its
//  Context.cpp counterpart has been compiled or run (coordinator directive,
//  "drop compiling and running entirely") -- treat it as reviewed-on-paper
//  only until a build confirms it.
// ============================================================================

#include <cstdint>
#include <vector>

#include "World.hpp"     // AcceptRejectMode
#include "robot_math.hpp" // robo::Vec3

// Exchange acceptance rule selecting the outer Markov chain over
// thermodynamic-state permutations (spec B9). RUN_TYPE::Default runs
// independent replicas with no exchange attempts.
enum class RUN_TYPE : std::uint8_t {
    Default = 0,
    REMC,       // Stage 1: label-swap parallel tempering (B6 ETerm_equal).
    RENEMC,     // Stage 2b: acceptance (ETerm_nonequil) wired; round-loop drive is Stage 2c TODO.
    RENE,       // Stage 2b: driven (BAT-scaling, work-based, WTerm) -- fully wired.
    REBASONTOP  // Stage 2b: RENE work-swaps + interleaved REMC sub-rounds (D4) -- fully wired.
};

// Exchange topology (spec B7): Neighboring pairs adjacent thermodynamic
// states with alternating parity; All draws random pairs.
enum class ReplicaMixingScheme : std::uint8_t {
    All = 0,
    Neighboring = 1
};

// A persistent molecular configuration: committed Cartesian coordinates plus
// the energies measured on them (spec B1 ownership table, INV-6), PLUS the
// nonequilibrium TRIAL state (B2) held apart from the committed state until
// accept/reject. There are R of them (spec B0); worlds hold no persistent
// per-replica state between rounds (INV-3) -- a replica's coordinates are the
// only thing that survives across a Gibbs sweep.
class Replica {
    public:
    // ---- committed (equilibrium) state -----------------------------------
    // Committed Cartesian coordinates, nm, engine/OpenMM atom order.
    std::vector<robo::Vec3> atomsLocations;

    // U(atomsLocations), the OpenMM potential energy of the committed
    // configuration (INV-6: "a replica's stored potential SHALL equal the
    // energy of its stored coordinates"). REMC/Default never populate a
    // Fixman term here (D3: Fixman never enters the swap acceptance), so
    // `potential`/`referencePotential` coincide for those run types; RENE/
    // REBASONTOP still keep them equal for the SAME reason (D3's exclusion is
    // unconditional, not REMC-specific) -- the split exists purely for
    // interface parity with the ORIGINAL's two-potential design (B1), not
    // because this port ever diverges them.
    double potential = 0.0;
    double referencePotential = 0.0;
    // Committed Fixman potential (diagnostic only, D3/F6: computed by the
    // sampler for the Boltzmann-marginal biconditional INV-7, but NEVER added
    // into `potential`/`referencePotential` or any acceptance exponent).
    // Tracked so F4's atomic commit (INV-4) has a real quantity to promote
    // atomically alongside coordinates + potential + referencePotential.
    double FixmanPotential = 0.0;

    // ---- nonequilibrium trial state (B2, B5) -----------------------------
    // Populated only during a DRIVEN round (RENE/REBASONTOP; RENEMC once its
    // round-loop lands, Stage 2c). Reset to 0 (WORK/WORK_Jacobian) and
    // reseeded from the committed endpoint (WORK_atomsLocations) at the start
    // of each driven range (INV-5) by Context::driveReplica.
    std::vector<robo::Vec3> WORK_atomsLocations; // x^tau, the driven endpoint (== x' under D7)
    double WORK_potential = 0.0;                 // U(x^tau), unreduced physical PE
    double referenceWORK_potential = 0.0;        // == WORK_potential (D3 Fixman exclusion, mirrors referencePotential)
    double WORK_FixmanPotential = 0.0;           // trial Fixman (diagnostic only, D3)
    // Accumulated per-driven-range work/Jacobian (B5, B12 "live += accumulation"):
    // reset once per driven range, accumulated once per driven world visited
    // in that range (a range is normally a single D7 mdSteps=0 scaling world,
    // but B12 requires correctness for a multi-driven-world schedule too).
    double WORK = 0.0;          // sum of per-driven-world (U_curr - U_prev) (B5; Fixman excluded, D3)
    double WORK_Jacobian = 0.0; // sum of per-driven-world getDistortJacobianDetLog() (B5/D6)

    // F4 atomic commit (INV-4): on an ACCEPTED driven swap, promote the WORK_*
    // trial to committed, ALL FOUR quantities together -- coordinates,
    // potential, referencePotential, FixmanPotential -- so no partial/
    // inconsistent state is ever visible (the F4 Critical bug this replaces:
    // the original committed coords + potential only, silently leaving
    // referencePotential/FixmanPotential stale). Never called on reject
    // (nothing to promote; the committed state already stands, INV-3/B6 step
    // 7 "On reject: nothing").
    void commitWorkAsFinal() {
        atomsLocations = WORK_atomsLocations;
        potential = WORK_potential;
        referencePotential = referenceWORK_potential;
        FixmanPotential = WORK_FixmanPotential;
    }
};

// A temperature and a per-world simulation schedule (spec B1, B0). There are
// T of them, canonically R == T (B0) -- Context asserts this, never assumes
// it.
class ThermodynamicState {
    public:
    double temperature = 300.0;

    // Same world-index ordering for every state (B0 NOTE): [0..W-1]. Kept
    // explicit (rather than implicit 0..W-1) so a future permutation/subset
    // schedule (B12) does not require an interface change.
    std::vector<int> worldIndexes;

    // Per-world schedule, one entry per position in worldIndexes (spec I1/I3
    // "World runtime setters" Consequences bullet: these let the exchange
    // driver reset each world's timestep/MD-step-count/accept-mode every
    // round from the state currently occupying it). Stage 1 copies these
    // from the worlds' own configuration at setup time and they are
    // identical across states (B0: "per-world distort/flow/integrator
    // vectors are identical across states"); Stage 2's driven segment (D7:
    // mdSteps == 0 for distortOption < 0 worlds) is where they diverge.
    std::vector<double> timeSteps;
    std::vector<int> mdSteps;
    std::vector<AcceptRejectMode> acceptRejectModes;
};
