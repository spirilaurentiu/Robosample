#include <gtest/gtest.h>
#include <iostream>
#include <thread>

#include "SimTKcommon.h"

using namespace SimTK;

/**
 * @class TestTask
 * @brief Logic to verify that the executing thread remains constant throughout
 * the lifecycle of a ParallelTask.
 * * Current ParallelExecutor design relies on TLS to store task data. If
 * an initialized task moves between threads, the TLS will point to invalid
 * data. Test that the task division logic does not move ParallelTasks
 * between threads.
 */
class TestTask : public Parallel2DExecutor::Task {
    public:
    TestTask() = default;

    auto execute(int /*row*/, int /*col*/) -> void override {
        // Verification: The thread executing the task must be the same one
        // that ran the initialize() method.
        EXPECT_EQ(m_threadId, std::this_thread::get_id())
            << "Task moved to a different thread during execution!";
    }

    auto initialize() -> void override {
        m_threadId = std::this_thread::get_id();
    }

    auto finish() -> void override {
        // Verification: The thread finishing the task must be the same one
        // that ran initialize().
        EXPECT_EQ(m_threadId, std::this_thread::get_id())
            << "Task moved to a different thread during finish()!";

        // Reset the TLS ID for this thread to avoid cross-pollination between test runs
        m_threadId = std::thread::id{};
    }

    private:
    // TLS to store the ID of the thread that initialized this specific task instance
    static thread_local std::thread::id m_threadId;
};

// Definition of the thread-local static member
thread_local std::thread::id TestTask::m_threadId;

/**
 * @test SimTKCommon_NoMovingTasks_StableThreadId.DefaultExecutorMaintainsAffinity
 * @expected_behavior The Parallel2DExecutor should ensure that for various matrix sizes,
 * every task stays pinned to its initializing thread when using the default internal executor.
 */
TEST(SimTKCommon_NoMovingTasks_StableThreadId, DefaultExecutorMaintainsAffinity) {
    const int numCpus = ParallelExecutor::getNumProcessors();
    const int maxRange = (numCpus * 2) + 1;

    for (int size = 1; size < maxRange; ++size) {
        Parallel2DExecutor executor(size);
        TestTask task;

        // Execute the task over the full matrix range
        EXPECT_NO_THROW(executor.execute(task, Parallel2DExecutor::FullMatrix));
    }
}

/**
 * @test SimTKCommon_NoMovingTasks_StableThreadId.FixedPoolMaintainsAffinity
 * @expected_behavior When using an external ParallelExecutor pool, the 2D wrapper
 * must not cause tasks to migrate between the pool's worker threads.
 */
TEST(SimTKCommon_NoMovingTasks_StableThreadId, FixedPoolMaintainsAffinity) {
    const int numCpus = ParallelExecutor::getNumProcessors();
    const int poolSize = (numCpus > 1) ? numCpus : 1;
    const int maxRange = (numCpus * 2) + 1;

    ParallelExecutor pool(poolSize);

    for (int size = 1; size < maxRange; ++size) {
        Parallel2DExecutor executor(size, pool);
        TestTask task;

        EXPECT_NO_THROW(executor.execute(task, Parallel2DExecutor::FullMatrix));
    }
}

/**
 * @test SimTKCommon_NoMovingTasks_StableThreadId.SingleThreadMaintainsAffinity
 * @expected_behavior In a single-threaded execution context, the lifecycle must
 * remain strictly on the calling/worker thread.
 */
TEST(SimTKCommon_NoMovingTasks_StableThreadId, SingleThreadMaintainsAffinity) {
    const int numCpus = ParallelExecutor::getNumProcessors();
    const int maxRange = (numCpus * 2) + 1;

    // Force a single-threaded pool
    ParallelExecutor singleThreadExecutor(1);

    for (int size = 1; size < maxRange; ++size) {
        Parallel2DExecutor executor(size, singleThreadExecutor);
        TestTask task;

        EXPECT_NO_THROW(executor.execute(task, Parallel2DExecutor::FullMatrix));
    }
}

auto main(int argc, char** argv) -> int {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
