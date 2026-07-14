#pragma once

/**
 * @file ForceFactory.hpp
 * @brief The standard (non-alchemy) OpenMM Force builders, hoisted out of
 *        `OpenMMContext` (SPLIT-O3).
 *
 * Each `create*Force` builder is a pure function of `SystemTopology` (plus, for
 * `createCustomNonbondedForce`, `OpenMMContext::isPeriodic`): it heap-allocates
 * an `OpenMM::*Force`, populates it from the SoA topology, and returns it. It
 * reads no `OpenMMContext` mutable state.
 *
 * @par Ownership
 *      Every builder returns a raw owning pointer. Ownership transfers to the
 *      OpenMM `System` on `System::addForce`; until then the caller owns it and
 *      must not leak it. The builders do not assign a force group - the caller
 *      (`OpenMMSystemBuilder::build`) sets each force's group per the
 *      force-group policy after construction.
 *
 * @note These builders define the potential (the energy/force field), not the
 *       per-body reduction; the force->wrench invariants INV-1/INV-2 live in
 *       `ForceReducer` / the CUDA `reduceForces` kernel and are untouched here.
 *       The alchemy/decoupling builders live in `AlchemyForceFactory` (SPLIT-O4)
 *       because they read Region-A state.
 */

#include "OpenMM.h"
#include "TopologyElements.hpp"

namespace robo::forcefactory {

/// @brief Builds the main `NonbondedForce` (LJ + electrostatics) with the
///        topology's nonbonded method, cutoff, exceptions, and 1-4 pairs. The
///        PME/explicit path returns the force alchemy later scales in place.
/// @param[in] systemTopology  SoA topology. Borrowed.
/// @return Heap-owned force; caller transfers ownership via `System::addForce`.
/// @throws std::invalid_argument on an unsupported nonbonded method.
[[nodiscard]] auto createNonbondedForce(const SystemTopology& systemTopology) -> OpenMM::NonbondedForce*;
/// @brief Builds the implicit-solvent `GBSAOBCForce` (OBC2 + ACE surface term).
/// @param[in] systemTopology  SoA topology. Borrowed.
/// @return Heap-owned force; caller transfers ownership via `System::addForce`.
/// @throws std::invalid_argument if the method is periodic (explicit solvent).
[[nodiscard]] auto createGBSAOBCForce(const SystemTopology& systemTopology) -> OpenMM::GBSAOBCForce*;
/// @brief Builds the NBFIX `CustomNonbondedForce` (tabulated A/B coefficients),
///        mirroring the main force's exclusions to avoid double-counting LJ.
/// @param[in] systemTopology  SoA topology; NBFIX table must be well-formed.
/// @return Heap-owned force; caller transfers ownership via `System::addForce`.
/// @throws std::runtime_error on a malformed NBFIX table (bad size or index).
[[nodiscard]] auto createCustomNonbondedForce(const SystemTopology& systemTopology)
    -> OpenMM::CustomNonbondedForce*;
/**
 * @brief Adds the main `NonbondedForce`'s exception pairs (all 1-2/1-3
 *        exclusions and 1-4 scaled pairs) as exclusions on @p force.
 *
 * The CPU platform shares one neighbor list across every exclusion-using force
 * and rejects the Context ("All Forces must have identical exclusions") unless
 * the lists match exactly; CUDA/OpenCL route interaction-group custom forces
 * around the shared list. This is therefore a CPU-correctness requirement and
 * energy-neutral on every platform.
 *
 * @param[in,out] force           Custom force to receive the exclusions.
 * @param[in]     systemTopology  Source of the exclusion/1-4 pair lists.
 */
void addStandardExclusions(OpenMM::CustomNonbondedForce* force, const SystemTopology& systemTopology);
/// @brief Builds the bonded `HarmonicBondForce` (stiffness doubled to OpenMM's
///        `k(r-r0)^2` convention). @return Heap-owned force; caller adds to System.
[[nodiscard]] auto createHarmonicBondForce(const SystemTopology& systemTopology) -> OpenMM::HarmonicBondForce*;
/// @brief Builds the `HarmonicAngleForce` (stiffness doubled). @return Heap-owned
///        force; caller transfers ownership via `System::addForce`.
[[nodiscard]] auto createHarmonicAngleForce(const SystemTopology& systemTopology) -> OpenMM::HarmonicAngleForce*;
/// @brief Builds the `PeriodicTorsionForce`. @return Heap-owned force; caller
///        transfers ownership via `System::addForce`.
[[nodiscard]] auto createPeriodicTorsionForce(const SystemTopology& systemTopology)
    -> OpenMM::PeriodicTorsionForce*;
/// @brief Builds improper harmonic torsions as a `CustomTorsionForce` (wrapped
///        `k*dtheta^2`). @return Heap-owned force; caller adds to System.
[[nodiscard]] auto createImproperHarmonicTorsionForce(const SystemTopology& systemTopology)
    -> OpenMM::CustomTorsionForce*;
/// @brief Builds the CHARMM `CMAPTorsionForce` from the topology's grids.
/// @return Heap-owned force; caller transfers ownership via `System::addForce`.
[[nodiscard]] auto createCMAPTorsionForce(const SystemTopology& systemTopology) -> OpenMM::CMAPTorsionForce*;
/// @brief Builds Urey-Bradley 1-3 terms as a `HarmonicBondForce` (stiffness
///        doubled). @return Heap-owned force; caller adds to System.
[[nodiscard]] auto createUreyBradleyForce(const SystemTopology& systemTopology) -> OpenMM::HarmonicBondForce*;

} // namespace robo::forcefactory
