#pragma once
// ============================================================================
//  NcmcProtocol.hpp -- the pure lambda schedule for the per-molecule NCMC switch
//  (MoveType::NcmcSwitch). Extracted from World::protocolLambda so the schedule
//  is a free function of (step, ncmcSteps, holdFraction) with NO dependency on
//  World / OpenMM / RobotEngine, and can therefore be unit-tested directly.
//  World::protocolLambda delegates here, so production and tests exercise the
//  SAME code (no parallel re-implementation that could drift).
//
//  Schedule (Nilmeier, Crooks, Minh & Chodera 2011): a PALINDROMIC triangle
//  lambda: 1 -> 0 -> 1 across step in [0, ncmcSteps), with an OPTIONAL flat
//  lambda = 0 hold of round(holdFraction * ncmcSteps) steps centered at the
//  midpoint (the "uncaged stride" where the decoupled molecule moves freely).
//
//  REVERSIBILITY (why the schedule shape is load-bearing, not cosmetic). The
//  whole NCMC trajectory is one map T on (q,p): per substep we set lambda at
//  fixed q (a parameter change that does NOT move the state) then take ONE
//  velocity-Verlet step at that fixed lambda. Each fixed-lambda Verlet step is
//  itself momentum-flip reversible (F V_lambda F == V_lambda^-1), so the
//  composition T = V_{lambda_{N-1}} ... V_{lambda_0} satisfies F T F == T^-1
//  -- the condition that makes accept-on-endpoint-H an EXACT Metropolis test --
//  IFF the lambda sequence is a true palindrome: lambda_s == lambda_{N-1-s},
//  with BOTH endpoints pinned to lambda = 1 (s = 0 and s = N-1). The schedule
//  below is palindromic by construction (a symmetric tent in |step - center|),
//  so there is NO out-of-loop "final jump to lambda = 1" -- that off-by-one
//  jump is exactly what broke F-reversibility at O(1/ncmcSteps) before.
// ============================================================================

#include <algorithm>
#include <cmath>

namespace robo::ncmc {

// lambda for protocol substep `step` in [0, ncmcSteps). A palindromic tent:
// lambda(s) == lambda(N-1-s), lambda(0) == lambda(N-1) == 1, descending to a
// flat lambda = 0 hold centered at the trough. Defined for any integer step
// (steps outside [0, N) clamp into [0, 1]); ncmcSteps <= 1 is the degenerate
// single fully-coupled substep.
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