#pragma once

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include "OpenMM.h"
#include "TopologyElements.hpp"

/**
 * @brief NCMC per-molecule intermolecular decoupling: the alchemy force builders,
 *        the enable control surface, and the Region-A state (SPLIT-O4).
 *
 * Two construction paths, both driven by the single `lambda_inter` OpenMM global
 * parameter (set through `OpenMMContext::setAlchemicalLambda`):
 *   - vacuum/implicit (@ref createAlchemyCorrectionForce): exact linear A x rest
 *     scaling via a `(lambda_inter - 1) * standard` interaction group;
 *   - explicit solvent, PME/Ewald/CutoffPeriodic (@ref createAlchemyDecouplingForces):
 *     PME-exact charge scaling on the main `NonbondedForce` plus a soft-core
 *     A x rest force and a hard intra-A force.
 *
 * @par Lambda contract
 *      `lambda_inter == 1` reproduces the unmodified force field exactly (both
 *      paths); this endpoint fidelity is what NCMC/HMC acceptance relies on - the
 *      switching protocol between the endpoints is guidance only, since
 *      acceptance is evaluated at full Hamiltonian at `lambda == 1`.
 *      `lambda_inter == 0` fully removes the A<->rest coupling. The vacuum path
 *      is exactly linear in lambda; the explicit-solvent path additionally
 *      annihilates intra-A electrostatics for `lambda < 1` (a documented
 *      departure from pure decoupling, still exact and unbiased at `lambda == 1`).
 *
 * `OpenMMContext` owns one instance and forwards `enableAlchemy` to it;
 * `setAlchemicalLambda` stays on `OpenMMContext` (it drives the live Context,
 * which this factory never touches) and only queries @ref enabled here. This
 * factory owns the Region-A state so both builders can read it directly.
 */
class AlchemyForceFactory {
    public:
    /**
     * @brief Enables alchemy and records Region A as an ascending, deduplicated
     *        atom-index set (an arbitrary, possibly non-contiguous set).
     * @param[in] atomIndices  Region A atom indices; sorted and deduplicated
     *                         internally.
     */
    void enableAlchemy(const std::vector<int>& atomIndices);
    /// @brief Convenience overload: Region A = the contiguous range
    ///        `[atomBegin, atomEnd)`. @see enableAlchemy(const std::vector<int>&)
    void enableAlchemy(int atomBegin, int atomEnd);

    /// @brief Whether alchemy has been enabled (a `createAlchemy*` builder should run).
    [[nodiscard]] auto enabled() const -> bool {
        return alchemyEnabled_;
    }

    /**
     * @brief Builds the vacuum/implicit correction force so the total A x rest
     *        pair energy becomes `lambda_inter * standard` (LJ + Coulomb).
     *
     * Implemented as a `(lambda_inter - 1) * standard` term over an A x rest
     * interaction group. Every excluded pair is intramolecular (never A x rest),
     * so the mirrored exclusions are energy-neutral and exist only to satisfy the
     * CPU platform's shared-neighbor-list rule.
     *
     * @param[in] systemTopology  SoA topology. Borrowed.
     * @return Heap-owned force; ownership transfers to the `System` on `addForce`.
     */
    [[nodiscard]] auto createAlchemyCorrectionForce(const SystemTopology& systemTopology)
        -> OpenMM::CustomNonbondedForce*;
    /**
     * @brief Builds the explicit-solvent (PME/Ewald/CutoffPeriodic) decoupling
     *        forces and scales electrostatics in place on @p main.
     *
     * Reciprocal space cannot be localized to an A x rest pair list, so A's charge
     * is scaled on @p main via a `lambda_inter` charge offset (PME-exact) and A's
     * epsilon is zeroed there; A's LJ is rebuilt as a soft-core A x rest force and
     * a hard intra-A force. All three are driven by the same `lambda_inter`
     * parameter, so `setAlchemicalLambda` and `NcmcMove` need no change.
     *
     * @param[in,out] main            Main `NonbondedForce`; mutated (charge offset,
     *                                zeroed A epsilon). Not owned here.
     * @param[in]     systemTopology  SoA topology. Borrowed.
     * @return `{soft-core A x rest force, hard intra-A force}`, both heap-owned;
     *         ownership transfers to the `System` on `addForce`.
     */
    [[nodiscard]] auto createAlchemyDecouplingForces(const SystemTopology& systemTopology,
                                                     OpenMM::NonbondedForce* main)
        -> std::pair<OpenMM::CustomNonbondedForce*, OpenMM::CustomNonbondedForce*>;

    /**
     * @brief Records which force object the selected path built (or nullptr under
     *        explicit-solvent decoupling, where lambda is driven through the main
     *        force's charge offset).
     *
     * Called by `OpenMMContext::initialize()` right after it picks a path. The
     * stored handle is not read elsewhere; it mirrors the pre-split
     * `OpenMMContext::alchemyForce` write.
     *
     * @param[in] force  The built alchemy force, already owned by the `System`,
     *                   or nullptr. Borrowed as a handle only.
     */
    void setAlchemyForce(OpenMM::CustomNonbondedForce* force) {
        alchemyForce_ = force;
    }

    private:
    // NCMC alchemy state. alchemyAtoms_ is Region A: an ascending, deduplicated
    // atom-index set (may be non-contiguous; enableAlchemy(int,int) builds the
    // contiguous case). alchemyAtomSet_ mirrors it as a std::set<int> for O(log n)
    // membership tests in the two builders above.
    bool alchemyEnabled_ = false;
    std::vector<int> alchemyAtoms_;
    std::set<int> alchemyAtomSet_;
    OpenMM::CustomNonbondedForce* alchemyForce_ = nullptr; // owned by `system`
};
