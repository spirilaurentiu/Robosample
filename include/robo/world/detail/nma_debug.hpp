#pragma once

// ============================================================================
//  nma_debug.hpp -- shared opt-in trace toggle for the NMA Route B velocity
//  distortion (DistortOption::NMA).
//
//  Relocated verbatim from World.cpp's anonymous namespace (SPLIT-W8, pure
//  code motion): hoisted to a header because it is used by both
//  VelocityDistortion.cpp's nmaKineticCorrection and World.cpp's
//  reinitialize/currentTotalEnergy -- two separate translation units after
//  the World split, which an anonymous-namespace symbol cannot serve. Set
//  the env var ROBO_NMA_DEBUG=1 to see, per Gibbs block, exactly how the
//  momentum draw is re-pointed. Off by default so the per-block
//  reinitialize() path stays quiet.
// ============================================================================

#include <cstdlib>

/**
 * @brief Reports whether the NMA Route B velocity-distortion trace is enabled.
 *
 * Shared opt-in toggle for the NMA Route B momentum re-pointing
 * (DistortOption::NMA). When on, each Gibbs block logs exactly how the momentum
 * draw is re-pointed. Off by default so the per-block reinitialize() path stays
 * quiet.
 *
 * @return true iff the environment variable @c ROBO_NMA_DEBUG is set to a value
 *         that is neither empty nor "0". The environment is read once on the
 *         first call and cached for the process lifetime.
 */
inline bool nmaDebugEnabled() {
    static const bool on = [] {
        const char* e = std::getenv("ROBO_NMA_DEBUG");
        return e != nullptr && e[0] != '0' && e[0] != '\0';
    }();
    return on;
}
