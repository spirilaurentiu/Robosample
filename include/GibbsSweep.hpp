#pragma once

#include <vector>

#include "robot_math.hpp"

class World;

namespace robo::gibbs {

/**
 * @brief Advance one Gibbs block: load a configuration into a World, take one
 *        sampling move, and read the post-move configuration back out.
 *
 * The single sweep primitive shared by all replica-exchange drivers. One call
 * is exactly `setAtomsLocationsInGround(coords) -> generateSample() ->
 * getAtomsLocationsInGround() -> coords`. A driver builds a Gibbs sweep by
 * calling this once per World in schedule order on the same @p coords buffer.
 *
 * @param[in,out] w        The World for this Gibbs block. Borrowed; the move
 *                         reads @p coords into it, samples, and leaves the
 *                         World's transient sampler state overwritten. The
 *                         World retains no per-replica state that survives this
 *                         call (INV-3): the returned @p coords is the only
 *                         thing carried to the next block.
 * @param[in,out] coords   Per-atom Ground-frame Cartesian coordinates (nm,
 *                         engine/OpenMM atom order), the sole inter-world
 *                         currency (INV-3). On entry the configuration to
 *                         sample from; on return the post-move configuration
 *                         (unchanged on a rejected move). Must hold at least
 *                         @p numAtoms elements.
 * @param[in] numAtoms     Atom count to copy back out of the World. Must equal
 *                         the World's atom count and `coords.size()`.
 * @return `true` if `generateSample()` accepted the move, `false` if rejected.
 *         The accept flag is telemetry only; @p coords already reflects the
 *         accepted-or-restored configuration either way.
 */
bool stepWorldInGround(World& w, std::vector<robo::Vec3>& coords, int numAtoms);

} // namespace robo::gibbs
