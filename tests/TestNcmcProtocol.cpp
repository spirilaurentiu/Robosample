// ============================================================================
//  TestNcmcProtocol.cpp -- the pure lambda schedule (robo::ncmc::protocolLambda),
//  split out of TestNCMCWork.cpp (TEST-006) as the FAST half: pure protocol-
//  schedule algebra, no sampling, no engine stepping. Runs on every push.
//
//  THEORY UNDER TEST (Fix 1): the whole NCMC trajectory is one map T on (q,p).
//  Accept-on-endpoint-H is an EXACT Metropolis test only if T is momentum-flip
//  reversible, F T F == T^-1. Since T composes per-substep fixed-lambda Verlet
//  steps (each itself F-reversible), the composition is reversible IFF the lambda
//  schedule is a true PALINDROME, lambda(s) == lambda(N-1-s), with BOTH endpoints
//  pinned to lambda = 1. These tests encode exactly that. They FAIL on the old
//  offset-by-one schedule (which pinned only s=0 to 1 and was symmetric about the
//  trough INDEX, not about the array center) and PASS on the palindromic one.
//
//  The work-chain / integration-level consequences of this schedule (protocol
//  work, Jacobian invariance, F-reversibility of the REAL integrator, Crooks
//  sign flip, Fixman block-diagonal invariance) are the SLOW half, split into
//  TestNcmcWorkChain.cpp.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

#include "NCMCProtocol.hpp"

using namespace robo;

// lambda(0) is exactly 1 (the protocol starts fully coupled at the physical
// lambda=1 state acceptance uses).
TEST(NcmcProtocol, StartsAtOne) {
    for (int n : {2, 5, 10, 33, 100, 1000}) {
        EXPECT_DOUBLE_EQ(ncmc::protocolLambda(0, n, 0.0), 1.0) << "n=" << n;
    }
}

// BOTH ENDPOINTS PINNED: lambda(0) == lambda(N-1) == 1. This is the schedule-level
// statement that the trajectory begins AND ends at the fully-coupled lambda=1
// physical state, with NO out-of-loop "final jump to 1". The OLD schedule failed
// this (its last grid point sat at 1 - 2/N), which is what broke F-reversibility
// at O(1/N). Checked across parities and with a hold.
TEST(NcmcProtocol, BothEndpointsPinnedToOne) {
    for (int n : {2, 3, 4, 5, 20, 21, 100, 101}) {
        for (double hold : {0.0, 0.1, 0.3}) {
            EXPECT_DOUBLE_EQ(ncmc::protocolLambda(0, n, hold), 1.0) << "start n=" << n << " hold=" << hold;
            EXPECT_DOUBLE_EQ(ncmc::protocolLambda(n - 1, n, hold), 1.0)
                << "end n=" << n << " hold=" << hold;
        }
    }
}

// PALINDROME: lambda(s) == lambda(N-1-s) for EVERY step and EVERY hold. This is
// the exact discrete condition for F T F == T^-1 (the protocol is its own
// time-reverse), so the protocol ratio is 1 and no extra factor enters dH. The
// OLD schedule was symmetric about the trough INDEX (ramp), not about (N-1)/2, so
// it violated this for even N -- the reversibility bug Fix 1 closes.
TEST(NcmcProtocol, IsPalindrome) {
    for (int n : {2, 3, 4, 7, 16, 17, 64, 101, 256}) {
        for (double hold : {0.0, 0.1, 0.4}) {
            for (int s = 0; s < n; ++s) {
                EXPECT_DOUBLE_EQ(ncmc::protocolLambda(s, n, hold), ncmc::protocolLambda(n - 1 - s, n, hold))
                    << "palindrome broken at n=" << n << " hold=" << hold << " s=" << s;
            }
        }
    }
}

// lambda stays in [0,1] for every step, including any step at or beyond ncmcSteps
// (defensive: ncmcMove indexes s in [0,ncmcSteps)).
TEST(NcmcProtocol, RangeIsUnitInterval) {
    for (int n : {2, 7, 16, 64, 257}) {
        for (double hold : {0.0, 0.1, 0.5, 0.9}) {
            for (int s = 0; s <= n + 3; ++s) {
                const double l = ncmc::protocolLambda(s, n, hold);
                EXPECT_GE(l, 0.0) << "n=" << n << " hold=" << hold << " s=" << s;
                EXPECT_LE(l, 1.0) << "n=" << n << " hold=" << hold << " s=" << s;
            }
        }
    }
}

// The trough reaches EXACTLY 0 (the fully-decoupled, uncaged state) at the
// palindrome center when N is odd (a single center index s=(N-1)/2). The trough
// is the global minimum of the schedule either way.
TEST(NcmcProtocol, TroughIsMinimumAndReachesZero) {
    for (int n : {5, 11, 51, 201}) { // odd -> exact center index exists
        const int center = (n - 1) / 2;
        EXPECT_DOUBLE_EQ(ncmc::protocolLambda(center, n, 0.0), 0.0) << "n=" << n;
        for (int s = 0; s < n; ++s) {
            EXPECT_GE(ncmc::protocolLambda(s, n, 0.0), ncmc::protocolLambda(center, n, 0.0))
                << "center is not the minimum, n=" << n << " s=" << s;
        }
    }
}

// Monotone non-increasing on the down-ramp (s up to the center), then
// non-decreasing on the up-ramp.
TEST(NcmcProtocol, MonotoneDownThenUp) {
    const int n = 101;
    const int center = (n - 1) / 2;
    for (int s = 1; s <= center; ++s) {
        EXPECT_LE(ncmc::protocolLambda(s, n, 0.0), ncmc::protocolLambda(s - 1, n, 0.0)) << "down s=" << s;
    }
    for (int s = center + 1; s < n; ++s) {
        EXPECT_GE(ncmc::protocolLambda(s, n, 0.0), ncmc::protocolLambda(s - 1, n, 0.0)) << "up s=" << s;
    }
}

// The flat lambda=0 hold is CENTERED on the palindrome axis: every pinned-at-0
// step is symmetric about (N-1)/2 (so the hold cannot itself break the palindrome,
// the uncaged stride straddles the trough). We assert the zero-set is palindromic
// and non-empty for a real hold.
TEST(NcmcProtocol, HoldZeroBlockIsCentered) {
    const int n = 100;
    const double hold = 0.4;
    int zeros = 0;
    for (int s = 0; s < n; ++s) {
        const bool isZero = (ncmc::protocolLambda(s, n, hold) == 0.0);
        if (isZero) {
            ++zeros;
            EXPECT_DOUBLE_EQ(ncmc::protocolLambda(n - 1 - s, n, hold), 0.0)
                << "zero-set not symmetric at s=" << s;
        }
    }
    EXPECT_GT(zeros, 0) << "a real hold produced no lambda=0 steps";
}

// Degenerate guards: N<=1 is a single fully-coupled substep, and holdFraction>=1
// (whole protocol held) still yields finite lambdas in [0,1].
TEST(NcmcProtocol, DegenerateBudgetsAreSafe) {
    EXPECT_DOUBLE_EQ(ncmc::protocolLambda(0, 1, 0.0), 1.0);
    for (int n : {0, 1, 2, 3}) {
        for (double hold : {0.0, 1.0, 2.0}) {
            for (int s = 0; s <= std::max(n, 1) + 1; ++s) {
                const double l = ncmc::protocolLambda(s, n, hold);
                EXPECT_TRUE(std::isfinite(l)) << "n=" << n << " hold=" << hold << " s=" << s;
                EXPECT_GE(l, 0.0);
                EXPECT_LE(l, 1.0);
            }
        }
    }
}
