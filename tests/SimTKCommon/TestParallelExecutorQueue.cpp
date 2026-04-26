#include <algorithm>
#include <gtest/gtest.h>
#include <numeric>
#include <vector>

#include "SimTKcommon.h"

using namespace SimTK;

namespace {

/**
 * @class SetFlagTask
 * @brief Task implementation to mark a specific index in a flag array.
 *
 * Original authors' intent: Ensure each task is executed exactly once.
 * The assertion within execute() catches race conditions where a task might
 * be picked up by multiple workers.
 */
class SetFlagTask : public ParallelWorkQueue::Task {
    public:
    SetFlagTask(std::vector<int>& flags, int index)
        : m_flags(flags)
        , m_index(index) {
    }

    auto execute() -> void override {
        // If a task is executed twice, this expectation will fail.
        // We use static_cast to avoid signed/unsigned comparison warnings.
        const auto idx = static_cast<std::size_t>(m_index);

        EXPECT_EQ(m_flags[idx], 0) << "Task at index " << m_index << " executed more than once!";
        m_flags[idx] = 1;
    }

    private:
    std::vector<int>& m_flags;
    const int m_index;
};

} // namespace

/**
 * @test SimTKCommon_ParallelExecutorQueue_TaskExecution, ProcessesAllSubmittedTasks
 * @details This test populates a queue with 490 tasks (numFlags - 10) and verifies
 * that after flush() is called, exactly those indices are marked.
 */
TEST(SimTKCommon_ParallelExecutorQueue_TaskExecution, ProcessesAllSubmittedTasks) {
    const int numFlags = 500;
    const int tasksToSubmit = (numFlags - 10);

    // Using std::vector<int> for C++17 compatibility and clarity
    std::vector<int> flags(static_cast<std::size_t>(numFlags), 0);

    {
        // Create a queue with 10 worker threads
        ParallelWorkQueue queue(10);

        for (int i = 0; i < tasksToSubmit; ++i) {
            queue.addTask(new SetFlagTask(flags, i));
        }

        // flush() must block until all workers are finished
        queue.flush();
    }

    // Verify results
    for (int i = 0; i < numFlags; ++i) {
        const auto idx = static_cast<std::size_t>(i);
        if (i < tasksToSubmit) {
            EXPECT_EQ(flags[idx], 1) << "Task " << i << " was not executed.";
        } else {
            EXPECT_EQ(flags[idx], 0) << "Task " << i << " was executed but should not have been.";
        }
    }
}

/**
 * @test SimTKCommon_ParallelExecutorQueue_QueuePersistence, HandlesConsecutiveTaskBatches
 * @details Verifies that the queue can be reused for multiple batches of tasks
 * and successfully flushes each batch.
 */
TEST(SimTKCommon_ParallelExecutorQueue_QueuePersistence, HandlesConsecutiveTaskBatches) {
    const int numFlags = 100;
    std::vector<int> flags(static_cast<std::size_t>(numFlags), 0);
    ParallelWorkQueue queue(4);

    // Batch 1
    for (int i = 0; i < (numFlags / 2); ++i) {
        queue.addTask(new SetFlagTask(flags, i));
    }
    queue.flush();

    // Batch 2
    for (int i = (numFlags / 2); i < numFlags; ++i) {
        queue.addTask(new SetFlagTask(flags, i));
    }
    queue.flush();

    // All should be marked
    const auto allSet = std::all_of(flags.begin(), flags.end(), [](int f) {
        return f == 1;
    });
    EXPECT_TRUE(allSet);
}

/**
 * @test SimTKCommon_ParallelExecutorQueue_TaskSafety, PreventsDuplicateExecution
 * @details Validates the internal worker logic ensures no two threads execute the
 * same Task pointer.
 */
TEST(SimTKCommon_ParallelExecutorQueue_TaskSafety, PreventsDuplicateExecution) {
    const int numTasks = 1000;
    std::vector<int> flags(static_cast<std::size_t>(numTasks), 0);

    // High concurrency to increase chance of race condition
    ParallelWorkQueue queue(20);

    for (int i = 0; i < numTasks; ++i) {
        queue.addTask(new SetFlagTask(flags, i));
    }

    queue.flush();

    // If execute() expectation failed, GTest will report it.
    // We double-check the sum to ensure no tasks were missed.
    const int sum = std::accumulate(flags.begin(), flags.end(), 0);
    EXPECT_EQ(sum, numTasks);
}