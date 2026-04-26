#include <atomic>
#include <chrono>
#include <gtest/gtest.h>
#include <iostream>
#include <thread>

#include "SimTKsimbody.h"

using namespace SimTK;

namespace {

std::atomic<int> g_parallelCalls{0};
std::atomic<int> g_nonParallelCalls{0};

class ParallelForceImpl : public Force::Custom::Implementation {
    public:
    [[nodiscard]] auto shouldBeParallelIfPossible() const -> bool override {
        return true;
    }

    void calcForce(const State&, Vector_<SpatialVec>&, Vector_<Vec3>&, Vector&) const override {
        g_parallelCalls.fetch_add(1, std::memory_order_relaxed);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    [[nodiscard]] auto calcPotentialEnergy(const State& state) const -> Real override {
        return 0.0;
    }
};

class NonParallelForceImpl : public Force::Custom::Implementation {
    public:
    bool shouldBeParallelIfPossible() const override {
        return false;
    }

    void calcForce(const State&, Vector_<SpatialVec>&, Vector_<Vec3>&, Vector&) const override {
        g_nonParallelCalls.fetch_add(1, std::memory_order_relaxed);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    [[nodiscard]] auto calcPotentialEnergy(const State& state) const -> Real override {
        return 0.0;
    }
};

void runParallelSystem() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    for (int i = 0; i < 50; ++i) {
        Force::Custom(forces, new ParallelForceImpl());
    }

    system.realizeTopology();
    State state = system.getDefaultState();
    system.realize(state, Stage::Dynamics);
}

void runNonParallelSystem() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    for (int i = 0; i < 50; ++i) {
        Force::Custom(forces, new NonParallelForceImpl());
    }

    system.realizeTopology();
    State state = system.getDefaultState();
    system.realize(state, Stage::Dynamics);
}

} // namespace

TEST(Simbody_ParallelForces_CorrectInvocationCount, AllForcesAreEvaluatedExactlyOncePerSystem) {
    g_parallelCalls = 0;
    g_nonParallelCalls = 0;

    runParallelSystem();
    runNonParallelSystem();

    EXPECT_EQ(g_parallelCalls.load(), 50);
    EXPECT_EQ(g_nonParallelCalls.load(), 50);
}

TEST(Simbody_ParallelForces_ParallelExecution_IsFasterThanSerial, ParallelVersionShouldReduceWallClockTime) {
    if (std::thread::hardware_concurrency() <= 1) {
        GTEST_SKIP() << "Insufficient hardware concurrency for parallel test.";
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    runParallelSystem();
    auto t1 = std::chrono::high_resolution_clock::now();

    auto parallelMs = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    t0 = std::chrono::high_resolution_clock::now();
    runNonParallelSystem();
    t1 = std::chrono::high_resolution_clock::now();

    auto serialMs = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    // Relaxed threshold to avoid flakiness due to scheduling noise.
    EXPECT_LT(parallelMs, serialMs);
}

TEST(Simbody_ParallelForces_SystemRealization_DoesNotCrash, SystemCompletesDynamicsRealizationSuccessfully) {
    EXPECT_NO_THROW({
        runParallelSystem();
        runNonParallelSystem();
    });
}