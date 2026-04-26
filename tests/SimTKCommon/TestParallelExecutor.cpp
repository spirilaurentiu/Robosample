#include <algorithm>
#include <gtest/gtest.h>
#include <numeric>
#include <vector>

#include "SimTKcommon.h"

using namespace SimTK;

namespace {

/**
 * @class SetFlagTask
 * @brief A helper task to verify the ParallelExecutor lifecycle.
 * * This task simulates work by incrementing flags and using thread-local storage
 * to verify that initialize(), execute(), and finish() are called in sequence.
 */
class SetFlagTask : public ParallelExecutor::Task {
    public:
    SetFlagTask(std::vector<int>& flags, int& count, bool expectedParallel)
        : m_flags(flags)
        , m_count(count)
        , m_expectedParallel(expectedParallel) {
    }

    auto execute(int index) -> void override {
        m_flags[static_cast<size_t>(index)]++;
        m_localCount++;
        // Verify we are in the expected thread context (worker vs master)
        EXPECT_EQ(ParallelExecutor::isWorkerThread(), m_expectedParallel);
    }

    auto initialize() -> void override {
        m_localCount = 0;
        EXPECT_EQ(ParallelExecutor::isWorkerThread(), m_expectedParallel);
    }

    auto finish() -> void override {
        m_count += m_localCount;
        EXPECT_EQ(ParallelExecutor::isWorkerThread(), m_expectedParallel);
    }

    private:
    std::vector<int>& m_flags;
    int& m_count;
    const bool m_expectedParallel;
    static thread_local int m_localCount;
};

thread_local int SetFlagTask::m_localCount = 0;

} // namespace

/**
 * @test SimTKCommon_ParallelExecutor_Execution.IteratesCorrectNumberOfTimes
 * @details Validates that the executor hits exactly the requested number of indices
 * and that the thread-local counts are correctly aggregated in the finish() step.
 */
TEST(SimTKCommon_ParallelExecutor_Execution, IteratesCorrectNumberOfTimes) {
    const int numFlags = 100;
    const int tasksToExecute = (numFlags - 10);
    const bool isParallel = (ParallelExecutor::getNumProcessors() > 1);

    std::vector<int> flags(numFlags, 0);
    ParallelExecutor executor;

    // Repeat the execution to catch potential race conditions or state leaks
    for (int i = 0; i < 100; ++i) {
        int totalCount = 0;
        std::fill(flags.begin(), flags.end(), 0);

        SetFlagTask task(flags, totalCount, isParallel);
        executor.execute(task, tasksToExecute);

        EXPECT_EQ(totalCount, tasksToExecute);
        for (int j = 0; j < numFlags; ++j) {
            if (j < tasksToExecute) {
                EXPECT_EQ(flags[static_cast<size_t>(j)], 1);
            } else {
                EXPECT_EQ(flags[static_cast<size_t>(j)], 0);
            }
        }
    }
}

/**
 * @test SimTKCommon_ParallelExecutor_WorkerThread.IdentifiesWorkerContextCorrectly
 * @details Ensures that ParallelExecutor::isWorkerThread() accurately distinguishes
 * between the master thread calling execute() and the threads in the pool.
 */
TEST(SimTKCommon_ParallelExecutor_WorkerThread, IdentifiesWorkerContextCorrectly) {
    ParallelExecutor executor;

    // Should be false on the main test thread
    EXPECT_FALSE(ParallelExecutor::isWorkerThread());

    const int numTasks = 10;
    std::vector<int> flags(numTasks, 0);
    int count = 0;
    const bool isParallel = (ParallelExecutor::getNumProcessors() > 1);

    SetFlagTask task(flags, count, isParallel);
    executor.execute(task, numTasks);

    // Should remain false after execution returns
    EXPECT_FALSE(ParallelExecutor::isWorkerThread());
}

/**
 * @test SimTKCommon_ParallelExecutor_SingleThreaded.BehavesCorrectlyWithOneThread
 * @details Validates that providing '1' to the constructor forces serial-like
 * execution while still following the Task lifecycle.
 */
TEST(SimTKCommon_ParallelExecutor_SingleThreaded, BehavesCorrectlyWithOneThread) {
    const int numFlags = 100;
    const int tasksToExecute = (numFlags - 10);
    std::vector<int> flags(numFlags, 0);
    int totalCount = 0;

    // Explicitly request a single-threaded executor
    ParallelExecutor executor(1);

    // In a single-threaded pool, SimTK typically executes on the caller thread
    SetFlagTask task(flags, totalCount, false);

    executor.execute(task, tasksToExecute);

    EXPECT_EQ(totalCount, tasksToExecute);
    for (int j = 0; j < tasksToExecute; ++j) {
        EXPECT_EQ(flags[static_cast<size_t>(j)], 1);
    }
}

/**
 * @test SimTKCommon_ParallelExecutor_Resize.ReportsCorrectMaxThreads
 * @details Verifies the constructor correctly sets and reports the thread pool size limit.
 */
TEST(SimTKCommon_ParallelExecutor_Resize, ReportsCorrectMaxThreads) {
    for (int x = 1; x < 100; ++x) {
        ParallelExecutor executor(x);
        EXPECT_EQ(executor.getMaxThreads(), x);
    }
}
