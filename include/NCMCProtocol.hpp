#pragma once
/**
 * @file NCMCProtocol.hpp
 * @brief The pure lambda(step) coupling schedule for the per-molecule NCMC
 *        switch (Nilmeier, Crooks, Minh & Chodera 2011).
 *
 * @c World::protocolLambda delegates to @ref robo::ncmc::protocolLambda, so
 * production and tests exercise one definition. The coupling parameter follows
 * the convention @c lambda == 1 fully coupled, @c lambda == 0 fully decoupled;
 * the alchemy force this schedule drives SHALL use the same convention, so both
 * NCMC endpoints sit in the physical (fully-coupled) ensemble the work
 * accounting assumes.
 */

#include <algorithm>
#include <cmath>

namespace robo::ncmc {

/**
 * @brief Coupling parameter @c lambda for NCMC protocol substep @p step: a
 *        palindromic tent @c 1 -> 0 -> 1 with an optional flat @c lambda = 0 hold.
 * @param[in] step         Protocol substep index; conventionally in
 *                         @c [0, ncmcSteps).
 * @param[in] ncmcSteps    Total number of switching substeps.
 * @param[in] holdFraction Fraction of @p ncmcSteps spent in a flat
 *                         @c lambda = 0 hold centered on the trough
 *                         (the uncaged stride); 0 for a plain triangle.
 * @return @c lambda in @c [0, 1].
 * @post @c protocolLambda(0, N, .) == protocolLambda(N-1, N, .) == 1: both
 *       endpoints are fully coupled.
 * @post Palindromic: @c protocolLambda(s, N, .) == protocolLambda(N-1-s, N, .).
 *       The value decreases monotonically to the trough (@c lambda == 0 at the
 *       center) and increases monotonically back to 1; it is symmetric, not
 *       globally monotone. This palindrome is what makes accept-on-endpoint the
 *       exact Metropolis test for the reversible NCMC map.
 * @note Defined for any integer @p step: indices outside @c [0, N) clamp into
 *       @c [0, 1]. @c ncmcSteps <= 1 returns 1 (a single fully-coupled substep);
 *       a @p holdFraction wide enough to consume the whole protocol returns 0
 *       away from the pinned endpoints.
 */
[[nodiscard]] inline double protocolLambda(int step, int ncmcSteps, double holdFraction) {
    if (ncmcSteps <= 1) {
        return 1.0; // degenerate: a single, fully-coupled (lambda = 1) substep
    }
    const double center = (ncmcSteps - 1) / 2.0;          // palindrome axis (may be half-integer)
    const double holdHalf = (holdFraction * ncmcSteps) / 2.0; // half-width of the flat lambda=0 hold
    const double denom = center - holdHalf;               // ramp half-width in index units
    if (denom <= 0.0) {
        return 0.0; // degenerate: the hold spans the whole protocol
    }
    // Symmetric tent: 1 at the endpoints (|step-center| == center), 0 inside the
    // hold (|step-center| <= holdHalf), linear in between. Palindromic in step.
    const double lam = (std::abs(static_cast<double>(step) - center) - holdHalf) / denom;
    return std::clamp(lam, 0.0, 1.0);
}

} // namespace robo::ncmc