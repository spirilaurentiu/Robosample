#include <algorithm>
#include <gtest/gtest.h>
#include <vector>

#include "SimTKcommon.h"

using namespace SimTK;

namespace {

/**
 * @class SetFlagTask
 * @brief Task implementation to track visits to matrix indices.
 * * Verified behaviors:
 * 1. Correctness of 2D index mapping.
 * 2. Thread-safety of accumulation via localCount.
 */
class SetFlagTask : public Parallel2DExecutor::Task {
    public:
    explicit SetFlagTask(Array_<Array_<int>>& flags, int& count)
        : m_flags(flags)
        , m_count(count) {
    }

    auto execute(int i, int j) -> void override {
        m_flags[i][j]++;
        m_localCount++;
    }

    auto initialize() -> void override {
        m_localCount = 0;
    }

    auto finish() -> void override {
        m_count += m_localCount;
    }

    private:
    Array_<Array_<int>>& m_flags;
    int& m_count;
    static thread_local int m_localCount;
};

thread_local int SetFlagTask::m_localCount = 0;

/**
 * @brief Helper to reset the flag matrix state.
 */
auto clearFlags(Array_<Array_<int>>& flags, int size) -> void {
    flags.resize(size);
    for (int i = 0; i < size; ++i) {
        flags[i].resize(size);
        for (int j = 0; j < size; ++j) {
            flags[i][j] = 0;
        }
    }
}

} // namespace

/**
 * @details Validates that FullMatrix covers the entire i*j domain.
 */
TEST(SimTKCommon_Parallel2DExecutor_FullMatrix, ProcessesEveryElementInGrid) {
    const int numFlags = 100;
    const int iterations = 10;
    Array_<Array_<int>> flags;
    clearFlags(flags, numFlags);

    Parallel2DExecutor executor(numFlags);

    for (int iter = 0; iter < iterations; ++iter) {
        int count = 0;
        SetFlagTask task(flags, count);
        executor.execute(task, Parallel2DExecutor::FullMatrix);

        const int expectedTotal = (numFlags * numFlags);
        EXPECT_EQ(count, expectedTotal);

        for (int i = 0; i < numFlags; ++i) {
            for (int j = 0; j < numFlags; ++j) {
                EXPECT_EQ(flags[i][j], (iter + 1));
            }
        }
    }
}

/**
 * @details Validates strict lower triangle: j < i.
 */
TEST(SimTKCommon_Parallel2DExecutor_HalfMatrix, ProcessesOnlyStrictLowerTriangle) {
    const int numFlags = 100;
    Array_<Array_<int>> flags;
    clearFlags(flags, numFlags);

    Parallel2DExecutor executor(numFlags);
    int count = 0;
    SetFlagTask task(flags, count);

    executor.execute(task, Parallel2DExecutor::HalfMatrix);

    const int expectedCount = ((numFlags * (numFlags - 1)) / 2);
    EXPECT_EQ(count, expectedCount);

    for (int i = 0; i < numFlags; ++i) {
        for (int j = 0; j < numFlags; ++j) {
            const int expectedVal = (j < i) ? 1 : 0;
            EXPECT_EQ(flags[i][j], expectedVal);
        }
    }
}

/**
 * @details Validates lower triangle plus diagonal: j <= i.
 */
TEST(SimTKCommon_Parallel2DExecutor_HalfPlusDiagonal, ProcessesLowerTriangleAndDiagonal) {
    const int numFlags = 100;
    Array_<Array_<int>> flags;
    clearFlags(flags, numFlags);

    Parallel2DExecutor executor(numFlags);
    int count = 0;
    SetFlagTask task(flags, count);

    executor.execute(task, Parallel2DExecutor::HalfPlusDiagonal);

    const int expectedCount = ((numFlags * (numFlags + 1)) / 2);
    EXPECT_EQ(count, expectedCount);

    for (int i = 0; i < numFlags; ++i) {
        for (int j = 0; j < numFlags; ++j) {
            const int expectedVal = (j <= i) ? 1 : 0;
            EXPECT_EQ(flags[i][j], expectedVal);
        }
    }
}

/**
 * @details Ensures that serial execution (1 thread) results in the same coverage.
 */
TEST(SimTKCommon_Parallel2DExecutor_SingleThread, RespectsRangesWithOneWorker) {
    const int numFlags = 50;
    Array_<Array_<int>> flags;

    // Explicitly 1 thread
    Parallel2DExecutor executor(numFlags, 1);

    // Test Full Matrix Range
    clearFlags(flags, numFlags);
    int countFull = 0;
    SetFlagTask taskFull(flags, countFull);
    executor.execute(taskFull, Parallel2DExecutor::FullMatrix);
    EXPECT_EQ(countFull, (numFlags * numFlags));

    // Test Half Matrix Range
    clearFlags(flags, numFlags);
    int countHalf = 0;
    SetFlagTask taskHalf(flags, countHalf);
    executor.execute(taskHalf, Parallel2DExecutor::HalfMatrix);
    EXPECT_EQ(countHalf, ((numFlags * (numFlags - 1)) / 2));
}
