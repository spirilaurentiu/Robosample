/**
 * TestSimTKCommon_Random.cpp
 *
 * Google Test suite for SimTK::Random::Uniform and SimTK::Random::Gaussian.
 *
 * What the original authors were testing:
 *
 * 1. Default construction parameters:
 *    - Uniform defaults to [0.0, 1.0); Gaussian defaults to mean=0, stddev=1.
 *
 * 2. Statistical correctness of generated samples:
 *    - Values produced by getValue() / fillArray() must fall within declared
 *      bounds on every call.
 *    - Bin counts across many samples must stay within 4 standard deviations
 *      of the analytically expected count (a χ²-style sanity check).
 *
 * 3. Seed determinism:
 *    - Re-seeding with the same value must reproduce the identical sequence.
 *    - fillArray() must produce the same sequence as repeated getValue() calls
 *      under the same seed.
 *
 * 4. Seed independence:
 *    - Two distinct seeds must produce different sequences.
 *    - Two independently constructed (unseeded) objects must start with
 *      different seeds, i.e. the library auto-seeds from a global source.
 *
 * 5. Range mutation:
 *    - Getters reflect values set via setMin/setMax (Uniform) and
 *      setMean/setStdDev (Gaussian).
 *    - Changing the range shifts the distribution accordingly.
 *
 * 6. Integer generation:
 *    - getIntValue() returns integers in [min, max) and their distribution
 *      is uniform.
 *
 * 7. Buffer safety:
 *    - fillArray() must not write beyond the requested element count.
 */

#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "SimTKcommon.h"

using SimTK::Random;
using SimTK::Real;

// ---------------------------------------------------------------------------
// Statistical helpers
// ---------------------------------------------------------------------------

namespace {

/**
 * Verify that each bin's observed count lies within 4 standard deviations of
 * its expected count.  Uses sqrt(expected) as the Poisson standard deviation,
 * which is the standard approach for large-sample histogram tests.
 */
auto verifyBinCounts(const std::vector<int>& expected, const std::vector<int>& found) -> void {
    ASSERT_EQ(expected.size(), found.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const double dev = std::sqrt(static_cast<double>(expected[i]));
        EXPECT_GE(found[i], (expected[i] - (4 * dev))) << "Bin " << i << " underflows";
        EXPECT_LE(found[i], (expected[i] + (4 * dev))) << "Bin " << i << " overflows";
    }
}

/**
 * Verify that every element of `values` lies in [min, max) and that the
 * 10-bin histogram is consistent with a uniform distribution over that range.
 */
auto verifyUniformDistribution(Real min, Real max, const std::vector<Real>& values) -> void {
    const int n = static_cast<int>(values.size());
    std::vector<int> expected(10, n / 10);
    std::vector<int> found(10, 0);

    for (const Real v : values) {
        ASSERT_GE(v, min) << "Value below declared minimum";
        ASSERT_LT(v, max) << "Value at or above declared maximum";
        const int index = static_cast<int>(((v - min) * 10) / (max - min));
        ++found[index];
    }
    verifyBinCounts(expected, found);
}

/**
 * Verify that every element of `values` lies in [min, max) and that the
 * per-integer histogram is consistent with a discrete uniform distribution.
 */
auto verifyUniformIntDistribution(int min, int max, const std::vector<int>& values) -> void {
    const int range = max - min;
    const int n = static_cast<int>(values.size());
    std::vector<int> expected(range, n / range);
    std::vector<int> found(range, 0);

    for (const int v : values) {
        ASSERT_GE(v, min) << "Integer value below declared minimum";
        ASSERT_LT(v, max) << "Integer value at or above declared maximum";
        ++found[v - min];
    }
    verifyBinCounts(expected, found);
}

/**
 * Verify that `values` follow a Gaussian distribution with the given
 * mean and standard deviation.  Uses a 6-bin scheme covering
 * (-∞, -2σ), [-2σ, -σ), [-σ, 0), [0, σ), [σ, 2σ), [2σ, +∞).
 */
auto verifyGaussianDistribution(Real mean, Real stddev, const std::vector<Real>& values) -> void {
    const int n = static_cast<int>(values.size());

    // Analytically derived probabilities for each bin under N(0,1).
    const int tail = static_cast<int>(0.0228 * n);
    const int outer = static_cast<int>((0.1587 * n) - tail);
    const int inner = static_cast<int>((0.5 * n) - outer);
    std::vector<int> expected = {tail, outer, inner, inner, outer, tail};
    std::vector<int> found(6, 0);

    for (const Real v : values) {
        const Real z = ((v - mean) / stddev);
        if (z < -2) {
            ++found[0];
        } else if (z < -1) {
            ++found[1];
        } else if (z < 0) {
            ++found[2];
        } else if (z < 1) {
            ++found[3];
        } else if (z < 2) {
            ++found[4];
        } else {
            ++found[5];
        }
    }
    verifyBinCounts(expected, found);
}

} // namespace

// ===========================================================================
// SimTKCommon_Random_Uniform
// ===========================================================================

// ---------------------------------------------------------------------------
// Default construction
// ---------------------------------------------------------------------------

/**
 * A default-constructed Uniform generator must report [0, 1) as its range.
 */
TEST(SimTKCommon_Random_Uniform, DefaultParametersAreZeroToOne) {
    const Random::Uniform rand;
    EXPECT_EQ(rand.getMin(), 0.0);
    EXPECT_EQ(rand.getMax(), 1.0);
}

// ---------------------------------------------------------------------------
// Statistical distribution over [0, 1)
// ---------------------------------------------------------------------------

/**
 * 2 000 samples drawn with a fixed seed must lie in [0, 1) and be
 * statistically uniform across 10 equal-width bins.
 */
TEST(SimTKCommon_Random_Uniform, SamplesAreUniformlyDistributedOverDefaultRange) {
    Random::Uniform rand;
    rand.setSeed(1);

    std::vector<Real> values(2000);
    for (Real& v : values) {
        v = rand.getValue();
    }
    verifyUniformDistribution(0.0, 1.0, values);
}

// ---------------------------------------------------------------------------
// Seed determinism – getValue
// ---------------------------------------------------------------------------

/**
 * Re-seeding with the same integer must reproduce the identical sequence from
 * getValue(), confirming that the PRNG is fully deterministic given its seed.
 */
TEST(SimTKCommon_Random_Uniform, ReseedingWithSameSeedReproducesGetValueSequence) {
    Random::Uniform rand;
    rand.setSeed(1);

    std::vector<Real> first(2000);
    for (Real& v : first) {
        v = rand.getValue();
    }

    rand.setSeed(1);
    for (int i = 0; i < 2000; ++i) {
        EXPECT_EQ(first[i], rand.getValue()) << "Mismatch at index " << i;
    }
}

// ---------------------------------------------------------------------------
// fillArray determinism
// ---------------------------------------------------------------------------

/**
 * fillArray() with the same seed must produce an element-wise identical array
 * to repeated getValue() calls under the same seed.
 */
TEST(SimTKCommon_Random_Uniform, FillArrayMatchesGetValueUnderSameSeed) {
    Random::Uniform rand;
    rand.setSeed(1);

    std::vector<Real> byGet(2000);
    for (Real& v : byGet) {
        v = rand.getValue();
    }

    std::vector<Real> byFill(2000);
    rand.setSeed(1);
    rand.fillArray(byFill.data(), static_cast<int>(byFill.size()));

    for (int i = 0; i < 2000; ++i) {
        EXPECT_EQ(byGet[i], byFill[i]) << "Mismatch at index " << i;
    }
}

// ---------------------------------------------------------------------------
// Seed independence
// ---------------------------------------------------------------------------

/**
 * Filling an array under seed 1 then seed 2 must yield at least one
 * differing element, confirming that distinct seeds produce distinct streams.
 */
TEST(SimTKCommon_Random_Uniform, DifferentSeedsProduceDifferentSequences) {
    Random::Uniform rand;

    std::vector<Real> seed1(2000);
    rand.setSeed(1);
    rand.fillArray(seed1.data(), static_cast<int>(seed1.size()));

    std::vector<Real> seed2(2000);
    rand.setSeed(2);
    rand.fillArray(seed2.data(), static_cast<int>(seed2.size()));

    // At least one element must differ; in practice all will differ.
    bool anyDifference = false;
    for (int i = 0; i < 2000; ++i) {
        if (seed1[i] != seed2[i]) {
            anyDifference = true;
            break;
        }
    }
    EXPECT_TRUE(anyDifference) << "Sequences under seed 1 and seed 2 are identical";
}

// ---------------------------------------------------------------------------
// Auto-seeding
// ---------------------------------------------------------------------------

/**
 * Two independently default-constructed Uniform objects must produce
 * different sequences, verifying that the library auto-seeds each instance
 * from a unique source (e.g. clock or counter).
 */
TEST(SimTKCommon_Random_Uniform, TwoDefaultInstancesProduceDifferentSequences) {
    Random::Uniform rand1;
    Random::Uniform rand2;

    std::vector<Real> v1(2000);
    rand1.fillArray(v1.data(), static_cast<int>(v1.size()));

    std::vector<Real> v2(2000);
    rand2.fillArray(v2.data(), static_cast<int>(v2.size()));

    bool anyDifference = false;
    for (int i = 0; i < 2000; ++i) {
        if (v1[i] != v2[i]) {
            anyDifference = true;
            break;
        }
    }
    EXPECT_TRUE(anyDifference) << "Two default-constructed instances produced identical sequences";
}

// ---------------------------------------------------------------------------
// Custom range – Real
// ---------------------------------------------------------------------------

/**
 * After calling setMin(5.0) and setMax(20.0), getters must reflect the new
 * range and 2 000 samples must be uniformly distributed over [5, 20).
 */
TEST(SimTKCommon_Random_Uniform, CustomRangeIsReflectedAndSamplesAreUniform) {
    Random::Uniform rand;
    rand.setMin(5.0);
    rand.setMax(20.0);

    EXPECT_EQ(rand.getMin(), 5.0);
    EXPECT_EQ(rand.getMax(), 20.0);

    std::vector<Real> values(2000);
    rand.fillArray(values.data(), static_cast<int>(values.size()));
    verifyUniformDistribution(5.0, 20.0, values);
}

// ---------------------------------------------------------------------------
// Integer generation
// ---------------------------------------------------------------------------

/**
 * getIntValue() must return integers in [min, max) following a uniform
 * discrete distribution when the generator range is [5.0, 20.0).
 */
TEST(SimTKCommon_Random_Uniform, IntValuesAreUniformlyDistributedInRange) {
    Random::Uniform rand;
    rand.setMin(5.0);
    rand.setMax(20.0);
    rand.setSeed(3);

    std::vector<int> values(2000);
    for (int& v : values) {
        v = rand.getIntValue();
    }
    verifyUniformIntDistribution(5, 20, values);
}

// ---------------------------------------------------------------------------
// Buffer safety
// ---------------------------------------------------------------------------

/**
 * fillArray(ptr, N) must not write to ptr[N]; a sentinel placed there before
 * the call must be unchanged afterwards.
 */
TEST(SimTKCommon_Random_Uniform, FillArrayDoesNotOverwriteBeyondRequestedCount) {
    Random::Uniform rand;
    rand.setSeed(1);

    std::array<Real, 2001> buffer{};
    buffer[2000] = 123.4;
    rand.fillArray(buffer.data(), 2000);

    EXPECT_EQ(buffer[2000], 123.4) << "fillArray() wrote past the end of the requested region";
}

// ===========================================================================
// SimTKCommon_Random_Gaussian
// ===========================================================================

// ---------------------------------------------------------------------------
// Default construction
// ---------------------------------------------------------------------------

/**
 * A default-constructed Gaussian generator must report mean=0, stddev=1.
 */
TEST(SimTKCommon_Random_Gaussian, DefaultParametersAreMeanZeroStdDevOne) {
    const Random::Gaussian rand;
    EXPECT_EQ(rand.getMean(), 0.0);
    EXPECT_EQ(rand.getStdDev(), 1.0);
}

// ---------------------------------------------------------------------------
// Statistical distribution N(0, 1)
// ---------------------------------------------------------------------------

/**
 * 2 000 samples from the default N(0, 1) generator must pass a 6-bin
 * Gaussian histogram test.
 */
TEST(SimTKCommon_Random_Gaussian, SamplesFollowStandardNormalDistribution) {
    Random::Gaussian rand;
    rand.setSeed(1);

    std::vector<Real> values(2000);
    for (Real& v : values) {
        v = rand.getValue();
    }
    verifyGaussianDistribution(0.0, 1.0, values);
}

// ---------------------------------------------------------------------------
// fillArray determinism
// ---------------------------------------------------------------------------

/**
 * fillArray() with the same seed must produce an element-wise identical array
 * to repeated getValue() calls under the same seed.
 */
TEST(SimTKCommon_Random_Gaussian, FillArrayMatchesGetValueUnderSameSeed) {
    Random::Gaussian rand;
    rand.setSeed(1);

    std::vector<Real> byGet(2000);
    for (Real& v : byGet) {
        v = rand.getValue();
    }

    std::vector<Real> byFill(2000);
    rand.setSeed(1);
    rand.fillArray(byFill.data(), static_cast<int>(byFill.size()));

    for (int i = 0; i < 2000; ++i) {
        EXPECT_EQ(byGet[i], byFill[i]) << "Mismatch at index " << i;
    }
}

// ---------------------------------------------------------------------------
// Custom parameters
// ---------------------------------------------------------------------------

/**
 * After setMean(10.0) and setStdDev(7.0), getters must reflect the new
 * values and 2 000 samples must follow N(10, 7).
 */
TEST(SimTKCommon_Random_Gaussian, CustomParametersAreReflectedAndSamplesAreCorrect) {
    Random::Gaussian rand;
    rand.setMean(10.0);
    rand.setStdDev(7.0);

    EXPECT_EQ(rand.getMean(), 10.0);
    EXPECT_EQ(rand.getStdDev(), 7.0);

    std::vector<Real> values(2000);
    rand.fillArray(values.data(), static_cast<int>(values.size()));
    verifyGaussianDistribution(10.0, 7.0, values);
}

// ---------------------------------------------------------------------------
// Buffer safety
// ---------------------------------------------------------------------------

/**
 * fillArray(ptr, N) must not write to ptr[N]; a sentinel placed there before
 * the call must be unchanged afterwards.
 */
TEST(SimTKCommon_Random_Gaussian, FillArrayDoesNotOverwriteBeyondRequestedCount) {
    Random::Gaussian rand;
    rand.setSeed(1);

    std::array<Real, 2001> buffer{};
    buffer[2000] = 123.4;
    rand.fillArray(buffer.data(), 2000);

    EXPECT_EQ(buffer[2000], 123.4) << "fillArray() wrote past the end of the requested region";
}