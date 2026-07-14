#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "OpenMM.h"
#include "TopologyElements.hpp"
#include "bridge/AlchemyForceFactory.hpp"

/**
 * @file OpenMMSystemBuilder.hpp
 * @brief The OpenMM System/Integrator/Context construction procedure hoisted out
 *        of `OpenMMContext::initialize()` (SPLIT-O5).
 *
 * `OpenMMSystemBuilder::build` constructs, in a fixed order that OpenMM requires:
 * particles; the periodic box (before Context creation, so PME can build its
 * reciprocal grid); virtual sites (before Context creation) plus the
 * massless-real-particle safety net; GBSA usability gating; force-group
 * assignment; every Force (via `robo::forcefactory` and the caller's
 * `AlchemyForceFactory`); the integrator (MTS r-RESPA or plain Verlet); and
 * finally the platform and Context. Ownership of the built objects transfers to
 * the caller in the returned result; `OpenMMContext` then moves them into its
 * own members.
 */

/**
 * @brief Owning result of `OpenMMSystemBuilder::build`; the caller moves these
 *        into its own members.
 *
 * @note Partial-write-on-failure: @ref system, @ref integrator,
 *       @ref hasVirtualSites, @ref numAtoms, and @ref forceGroupLabels are
 *       populated even when @ref success is false, because they are fully built
 *       before the Context construction that may throw. @ref context and
 *       @ref openMMVersion are valid only when @ref success is true.
 */
struct OpenMMSystemBuildResult {
    std::unique_ptr<OpenMM::System> system;      ///< Built System (owned). Populated even on failure.
    std::unique_ptr<OpenMM::Integrator> integrator; ///< Built integrator (owned). Populated even on failure.
    std::unique_ptr<OpenMM::Context> context;    ///< Built Context (owned); null iff `!success`.
    std::vector<std::pair<int, std::string>> forceGroupLabels; ///< `(group, name)` per added force.
    bool hasVirtualSites = false;                ///< True iff any virtual site was declared.
    std::size_t numAtoms = 0;                    ///< Particle count added to the System.
    /// OpenMM version string, returned (not printed) so the caller logs it after
    /// the `ROBO_CUDA_KINEMATICS` env read. Empty iff `!success`.
    std::string openMMVersion;
    bool success = false;                        ///< False iff Context construction threw.
};

/// @brief Stateless builder namespace-struct for the OpenMM bring-up sequence.
struct OpenMMSystemBuilder {
    /**
     * @brief Builds the OpenMM System, Integrator, and Context from
     *        @p systemTopology in the required construction order.
     *
     * @param[in]     systemTopology     SoA topology; box vectors arrive already
     *                                   reduced (lower-triangular). Borrowed.
     * @param[in]     useMTS             Build an `MTSIntegrator` (else Verlet);
     *                                   also forces the two-tier force-group split.
     * @param[in]     mtsInnerSubsteps   Fast-group substeps per outer step when
     *                                   @p useMTS.
     * @param[in]     separateForceGroups  When set (and not @p useMTS), give each
     *                                   force its own group.
     * @param[in,out] alchemyFactory     Read for Region-A state by the two
     *                                   `createAlchemy*Force` builders; its
     *                                   `setAlchemyForce` records which force the
     *                                   vacuum/implicit path built (nullptr under
     *                                   explicit-solvent decoupling).
     * @return The built objects; check `.success`.
     *
     * @throws std::runtime_error  on periodic setup without 9 box vectors, a
     *         cutoff exceeding half the smallest box width, a massless
     *         non-virtual-site particle, NCMC alchemy combined with NBFIX, or
     *         more than 32 separate force groups. Context construction failure is
     *         reported via `.success == false`, not thrown.
     * @post Forces are placed per the force-group policy: MTS -> slow group 0 /
     *       fast group 1; else `separateForceGroups` -> one group per force; else
     *       all in group 0.
     */
    [[nodiscard]] static auto build(const SystemTopology& systemTopology,
                                    bool useMTS,
                                    int mtsInnerSubsteps,
                                    bool separateForceGroups,
                                    AlchemyForceFactory& alchemyFactory) -> OpenMMSystemBuildResult;
};
