// -----------------------------------------------------------------------------
// TestAtomic.cpp
//
// Original history
// ----------------
// This test originally exercised the SimTK::AtomicInteger class template.
// When SimTK::AtomicInteger was retired in favour of C++11 std::atomic the
// low-level operator tests became trivially true (the standard guarantees
// atomicity), so the test was repurposed.
//
// What is actually being tested here
// -----------------------------------
// These tests verify that SimTK::ParallelExecutor correctly schedules N
// independent tasks across its thread pool AND that the std::atomic post-
// increment (++) and compound-addition (+=) operators survive genuine data
// races without losing updates.  Two invariants are checked:
//
//   1. Post-increment (++):  5 000 tasks each post-increment a shared atomic
//      counter and use the *returned* pre-increment value to stamp exactly one
//      slot in a flag array.  Because the operation is atomic, every slot
//      0..4 999 must be stamped exactly once per outer iteration; slots
//      5 000..9 999 must remain untouched.
//
//   2. Compound-addition (+=):  5 000 tasks each add 2 atomically to a
//      shared counter.  The final value must be exactly 10 000; any value
//      below that would indicate a lost write caused by a non-atomic RMW.
//
// Both scenarios repeat for 100 outer iterations to increase the probability
// of triggering latent races on any scheduler.
// -----------------------------------------------------------------------------

#include <atomic>
#include <vector>

#include "gtest/gtest.h"

#include "SimTKcommon.h"

using namespace SimTK;

// =============================================================================
// Anonymous namespace — internal task helpers
// =============================================================================
namespace {

// -----------------------------------------------------------------------------
// SetFlagTask
//
// Each call to execute() atomically post-increments `index_` and uses the
// *pre-increment* value to bump the corresponding slot in `flags_`.
// If ++ were not atomic two threads could obtain the same index and write
// to the same slot, causing a count greater than expected and/or leaving
// other slots at zero.
// -----------------------------------------------------------------------------
class SetFlagTask final : public ParallelExecutor::Task {
    public:
    SetFlagTask(std::vector<int>& flags, std::atomic<int>& index)
        : flags_(flags)
        , index_(index) {
    }

    auto execute(int /*taskIndex*/) -> void override {
        // Atomic post-increment: fetch current value, increment, return old.
        // The old value is a unique slot guaranteed to no other thread.
        flags_[index_++]++;
    }

    private:
    std::vector<int>& flags_;
    std::atomic<int>& index_;
};

// -----------------------------------------------------------------------------
// IncrementTask
//
// Each call to execute() atomically adds `kStep` to `index_`.
// If += were not atomic concurrent read-modify-writes would overlap and some
// additions would be silently discarded, producing a total below
// kTaskCount * kStep.
// -----------------------------------------------------------------------------
class IncrementTask final : public ParallelExecutor::Task {
    public:
    explicit IncrementTask(std::atomic<int>& index)
        : index_(index) {
    }

    auto execute(int /*taskIndex*/) -> void override {
        index_ += kStep;
    }

    // Exposed so that callers can compute the expected total without
    // duplicating the magic constant.
    static constexpr int kStep = 2;

    private:
    std::atomic<int>& index_;
};

} // namespace

// =============================================================================
// Test: post-increment operator is atomic under parallel execution
//
// Naming: SimTKCommon_Atomic_PostIncrementOperator
//         UniqueSlotAssignedToEveryParallelTask
//
// Setup
// -----
//   flags  : array of 10 000 ints, all initialised to 0
//   index  : atomic counter starting at 0
//   tasks  : 5 000 — each task post-increments index and stamps flags[old]
//
// Expected invariants after each of 100 iterations (iteration i, 0-based)
// -----------------------------------------------------------------------
//   index          == kTaskCount          (counter advanced exactly once/task)
//   flags[j]       == i + 1   for j in [0, kTaskCount)  (each slot hit once)
//   flags[j]       == 0       for j in [kTaskCount, kFlagCount) (untouched)
// =============================================================================
TEST(SimTKCommon_Atomic_PostIncrementOperator, UniqueSlotAssignedToEveryParallelTask) {
    constexpr int kFlagCount = 10'000;
    constexpr int kTaskCount = 5'000;
    constexpr int kIterations = 100;

    std::atomic<int> index{0};
    std::vector<int> flags(kFlagCount, 0);

    ParallelExecutor executor;

    for (int i = 0; i < kIterations; ++i) {
        SetFlagTask task(flags, index);
        index = 0;
        executor.execute(task, kTaskCount);

        // The atomic post-increment must have advanced the counter to exactly
        // kTaskCount.  Any deviation means a task obtained a duplicate or
        // out-of-range index.
        ASSERT_EQ(index.load(), kTaskCount) << "Atomic index mismatch on iteration " << i;

        for (int j = 0; j < kFlagCount; ++j) {
            // Slots below kTaskCount accumulate one visit per iteration.
            // Slots at or above kTaskCount must never be touched.
            const int expected = (j < kTaskCount) ? (i + 1) : 0;
            EXPECT_EQ(flags[j], expected) << "flags[" << j << "] wrong on iteration " << i;
        }
    }
}

// =============================================================================
// Test: compound-addition operator is atomic under parallel execution
//
// Naming: SimTKCommon_Atomic_CompoundAdditionOperator
//         TotalEqualsTaskCountTimesStepAfterParallelExecution
//
// Setup
// -----
//   index  : atomic counter starting at 0
//   tasks  : 5 000 — each task adds IncrementTask::kStep (= 2) atomically
//
// Expected invariant after each of 100 iterations
// ------------------------------------------------
//   index == kTaskCount * IncrementTask::kStep   (no update must be lost)
// =============================================================================
TEST(SimTKCommon_Atomic_CompoundAdditionOperator, TotalEqualsTaskCountTimesStepAfterParallelExecution) {
    constexpr int kTaskCount = 5'000;
    constexpr int kIterations = 100;
    constexpr int kExpected = kTaskCount * IncrementTask::kStep;

    std::atomic<int> index{0};

    ParallelExecutor executor;

    for (int i = 0; i < kIterations; ++i) {
        IncrementTask task(index);
        index = 0;
        executor.execute(task, kTaskCount);

        // If += were not atomic, overlapping reads would cause some
        // additions to be computed on a stale value and then overwrite a
        // more recent result, so the total would fall below kExpected.
        ASSERT_EQ(index.load(), kExpected) << "+= total wrong on iteration " << i;
    }
}
