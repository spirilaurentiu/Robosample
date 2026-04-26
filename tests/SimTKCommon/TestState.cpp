// =============================================================================
// TestState.cpp
//
// Google Test conversion of the original SimTK CTest-based StateTest.
//
// What the authors tested (three logical domains):
//
//  1. Cache-entry validity & prerequisites  (testCacheValidity)
//     - Allocation metadata (stage, initial value) is recorded faithfully.
//     - A cache entry is inaccessible until its depends-on stage is realized.
//     - Unsound cache-to-cache prerequisite edges are rejected at allocation.
//     - Dependency lists of z, discrete variables, and upstream cache entries
//       are populated correctly after Instance stage is realized.
//     - Marking a cache entry realized makes it readable; subsequently
//       modifying any prerequisite invalidates it without changing system stage.
//     - Copy-construction reconstructs all dependency lists identically.
//     - Modifying a Model-stage discrete variable invalidates and *deallocates*
//       Instance-stage cache entries, removing them from prerequisite lists.
//     - Re-advancing through the guaranteed-valid stage makes a cache entry
//       automatically accessible without an explicit markCacheValueRealized().
//
//  2. Basic State operations  (testMisc)
//     - Time is exactly 0 after advancing to Topology stage.
//     - Q allocation returns contiguous, ascending indices.
//     - Event-trigger allocation returns distinct, non-negative indices.
//     - A discrete variable retains its value after an explicit update.
//
//  3. Inter-State consistency  (testConsistent)
//     - Two identically-structured empty states are consistent at Instance stage.
//     - isConsistent() throws below Instance stage.
//     - Differing Q / U / Z / QErr / UErr / UDotErr sizes break consistency;
//       matching sizes restore it.
//     - Per-stage event-trigger counts participate in the consistency check.
//     - A different number of subsystems makes two states inconsistent.
// =============================================================================

#include <gtest/gtest.h>
#include <string>

#include "SimTKcommon.h"

using namespace SimTK;

// -----------------------------------------------------------------------------
// Shared helper: advance the State one stage forward (from stage-1 to `stage`)
// for every subsystem, then advance the system itself.
// -----------------------------------------------------------------------------
auto advanceStage(State& state, Stage stage) -> void {
    for (SubsystemIndex sx(0); sx < state.getNumSubsystems(); ++sx) {
        state.advanceSubsystemToStage(sx, stage);
    }
    state.advanceSystemToStage(stage);
}

// =============================================================================
// Domain 1 – Cache-entry validity & prerequisites
// =============================================================================

// -----------------------------------------------------------------------------
// After allocating discrete variables and cache entries at Topology stage,
// the stored allocation stage, invalidation stage, and initial values must
// match what was requested.
// -----------------------------------------------------------------------------
TEST(SimTKCommon_State_CacheEntry, AllocationMetadataIsRecordedCorrectly) {
    const SubsystemIndex Sub0(0);
    const SubsystemIndex Sub1(1);

    State s;
    s.setNumSubsystems(2);

    EXPECT_EQ(s.getSystemStage(), Stage::Empty);

    // Allocate at Topology stage a Model-stage-invalidating state variable.
    const DiscreteVariableIndex dvx1TopoModel =
        s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    EXPECT_EQ(s.getDiscreteVarAllocationStage(Sub1, dvx1TopoModel), Stage::Topology);
    EXPECT_EQ(s.getDiscreteVarInvalidatesStage(Sub1, dvx1TopoModel), Stage::Model);
    EXPECT_EQ(Value<Real>::downcast(s.getDiscreteVariable(Sub1, dvx1TopoModel)), Real(2));

    // Allocate at Topology stage a cache entry that depends on Model stage
    // and is guaranteed to be valid at Time stage. In between (at Model or
    // Instance stage) it *may* be valid if explicitly marked so.
    const CacheEntryIndex cx0TopoModel =
        s.allocateCacheEntry(Sub0, Stage::Model, Stage::Time, new Value<int>(41));

    EXPECT_EQ(s.getCacheEntryAllocationStage(Sub0, cx0TopoModel), Stage::Topology);
}

// -----------------------------------------------------------------------------
// Accessing a cache entry before its depends-on stage has been realized must
// throw, even when its allocation stage (Topology) has already been reached.
// -----------------------------------------------------------------------------
TEST(SimTKCommon_State_CacheEntry, AccessBeforeDependsOnStageThrows) {
    const SubsystemIndex Sub0(0);
    const SubsystemIndex Sub1(1);

    State s;
    s.setNumSubsystems(2);

    // Provide the Model-invalidating variable so Topology can be realized.
    s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    // depends-on = Model, guaranteed-valid = Time.
    const CacheEntryIndex cx0TopoModel =
        s.allocateCacheEntry(Sub0, Stage::Model, Stage::Time, new Value<int>(41));

    advanceStage(s, Stage::Topology);

    // Topology is realized, but the depends-on stage (Model) has not been
    // reached yet — access must throw.
    EXPECT_THROW(s.getCacheEntry(Sub0, cx0TopoModel), std::exception);
}

// -----------------------------------------------------------------------------
// When a proposed prerequisite cache entry is invalidated at a finer stage
// than the dependent entry, the dependency is unsound and must be rejected
// at allocation time.
//
// Specifically: cx0TopoVelocity is invalidated when Velocity stage changes.
// A dependent entry that only gets invalidated at Position could appear valid
// even after Velocity changes — so the edge must be rejected.
// -----------------------------------------------------------------------------
TEST(SimTKCommon_State_CacheEntry, UnsoundCacheToCachePrerequisiteIsRejected) {
    const SubsystemIndex Sub0(0);
    const SubsystemIndex Sub1(1);

    State s;
    s.setNumSubsystems(2);

    s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    // Here is a cache entry allocated at Topology stage, with depends-on
    // Velocity and no good-by guarantee.
    const CacheEntryIndex cx0TopoVelocity =
        s.allocateCacheEntry(Sub0, Stage::Velocity, Stage::Infinity, new Value<char>('v'));

    advanceStage(s, Stage::Topology);
    advanceStage(s, Stage::Model);

    // This attempt to create a cache-to-cache dependency should fail because
    // the prerequisite gets invalidated when Velocity stage changes but
    // the "dependent" doesn't get invalidated unless Position stage does. Thus
    // a velocity could change, invalidating the prereq, but the downstream
    // cache entry still looks valid.
    EXPECT_THROW(s.allocateCacheEntryWithPrerequisites(Sub1,
                                                       Stage::Position,
                                                       Stage::Infinity,
                                                       false,
                                                       false,
                                                       false,
                                                       {},
                                                       {CacheEntryKey(Sub0, cx0TopoVelocity)},
                                                       new Value<int>(-1)),
                 std::exception);
}

// -----------------------------------------------------------------------------
// After Instance stage is realized, each cache entry with prerequisites must
// appear in the dependency list of every prerequisite it declared (z vector,
// discrete variable, and upstream cache entry).
// -----------------------------------------------------------------------------
TEST(SimTKCommon_State_CacheEntry, DependencyListsArePopulatedAtInstanceStage) {
    const SubsystemIndex Sub0(0);
    const SubsystemIndex Sub1(1);

    State s;
    s.setNumSubsystems(2);

    s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    // Allocate at Topology stage a cache entry that depends on Model stage
    // and is guaranteed to be valid at Time stage.
    const CacheEntryIndex cx0TopoModel =
        s.allocateCacheEntry(Sub0, Stage::Model, Stage::Time, new Value<int>(41));

    advanceStage(s, Stage::Topology);

    // Allocate at Model stage a Position-invalidating state variable.
    const DiscreteVariableIndex dvx0ModelPos =
        s.allocateDiscreteVariable(Sub0, Stage::Position, new Value<int>(31));

    advanceStage(s, Stage::Model);

    // Allocate a cache entry at Instance stage that has Time as depends-on and
    // also has a cross-subsystem dependency on discrete variable dvx0ModelPos.
    const CacheEntryIndex cx1InstanceTime =
        s.allocateCacheEntryWithPrerequisites(Sub1,
                                              Stage::Time,
                                              Stage::Velocity,
                                              false,
                                              false,
                                              true, // depends on z
                                              {DiscreteVarKey(Sub0, dvx0ModelPos)},
                                              {CacheEntryKey(Sub0, cx0TopoModel)},
                                              new Value<std::string>("hasPrereqs_Time"));

    // Same but with depends-on Instance.
    const CacheEntryIndex cx1InstInst =
        s.allocateCacheEntryWithPrerequisites(Sub1,
                                              Stage::Instance,
                                              Stage::Velocity,
                                              false,
                                              false,
                                              true, // depends on z
                                              {DiscreteVarKey(Sub0, dvx0ModelPos)},
                                              {CacheEntryKey(Sub0, cx0TopoModel)},
                                              new Value<std::string>("hasPrereqs_Instance"));

    advanceStage(s, Stage::Instance);

    const CacheEntryKey ckey1(Sub1, cx1InstanceTime);
    const CacheEntryKey ckey2(Sub1, cx1InstInst);

    // z dependency list
    EXPECT_EQ(s.getZDependents().size(), 2u);
    EXPECT_TRUE(s.getZDependents().contains(ckey1));
    EXPECT_TRUE(s.getZDependents().contains(ckey2));

    // Discrete variable dependency list
    const DiscreteVarInfo& dvinfo = s.getDiscreteVarInfo(DiscreteVarKey(Sub0, dvx0ModelPos));
    EXPECT_EQ(dvinfo.getDependents().size(), 2u);
    EXPECT_TRUE(dvinfo.getDependents().contains(ckey1));
    EXPECT_TRUE(dvinfo.getDependents().contains(ckey2));

    // Upstream cache entry dependency list
    const CacheEntryInfo& ceinfo = s.getCacheEntryInfo(CacheEntryKey(Sub0, cx0TopoModel));
    EXPECT_EQ(ceinfo.getDependents().size(), 2u);
    EXPECT_TRUE(ceinfo.getDependents().contains(ckey1));
    EXPECT_TRUE(ceinfo.getDependents().contains(ckey2));
}

// -----------------------------------------------------------------------------
// Modifying any prerequisite (z, discrete variable, or upstream cache entry)
// must make the dependent cache entry inaccessible again without altering
// the system stage. Re-marking restores access.
// -----------------------------------------------------------------------------
TEST(SimTKCommon_State_CacheEntry, PrerequisiteModificationInvalidatesWithoutStageChange) {
    const SubsystemIndex Sub0(0);
    const SubsystemIndex Sub1(1);

    State s;
    s.setNumSubsystems(2);

    s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    const CacheEntryIndex cx0TopoModel =
        s.allocateCacheEntry(Sub0, Stage::Model, Stage::Time, new Value<int>(41));

    advanceStage(s, Stage::Topology);

    const DiscreteVariableIndex dvx0ModelPos =
        s.allocateDiscreteVariable(Sub0, Stage::Position, new Value<int>(31));

    advanceStage(s, Stage::Model);

    const CacheEntryIndex cx1InstanceTime =
        s.allocateCacheEntryWithPrerequisites(Sub1,
                                              Stage::Time,
                                              Stage::Velocity,
                                              false,
                                              false,
                                              true, // depends on z
                                              {DiscreteVarKey(Sub0, dvx0ModelPos)},
                                              {CacheEntryKey(Sub0, cx0TopoModel)},
                                              new Value<std::string>("hasPrereqs_Time"));

    const CacheEntryIndex cx1InstInst =
        s.allocateCacheEntryWithPrerequisites(Sub1,
                                              Stage::Instance,
                                              Stage::Velocity,
                                              false,
                                              false,
                                              true, // depends on z
                                              {DiscreteVarKey(Sub0, dvx0ModelPos)},
                                              {CacheEntryKey(Sub0, cx0TopoModel)},
                                              new Value<std::string>("hasPrereqs_Instance"));

    advanceStage(s, Stage::Instance);

    // Although cx0TopoModel *could* be valid at this point,
    // no one has said so, so we expect it to throw.
    EXPECT_THROW(s.getCacheEntry(Sub0, cx0TopoModel), std::exception);

    // If we say it is valid, we should be able to obtain its value.
    s.markCacheValueRealized(Sub0, cx0TopoModel);
    EXPECT_EQ(Value<int>::downcast(s.getCacheEntry(Sub0, cx0TopoModel)), 41);

    advanceStage(s, Stage::Time);

    // cx1InstanceTime isn't automatically valid but can be now.
    EXPECT_THROW(s.getCacheEntry(Sub1, cx1InstanceTime), std::exception);
    EXPECT_THROW(s.getCacheEntry(Sub1, cx1InstInst), std::exception);

    s.markCacheValueRealized(Sub1, cx1InstanceTime);
    s.markCacheValueRealized(Sub1, cx1InstInst);

    EXPECT_EQ(Value<std::string>::downcast(s.getCacheEntry(Sub1, cx1InstanceTime)).get(), "hasPrereqs_Time");
    EXPECT_EQ(Value<std::string>::downcast(s.getCacheEntry(Sub1, cx1InstInst)).get(), "hasPrereqs_Instance");

    // That cache entry does not depend on q or u so changes to them should
    // have no effect on stage or cache entry validity.
    s.updQ() = 0.;
    s.updU() = 0.;
    EXPECT_EQ(s.getSystemStage(), Stage::Time);
    EXPECT_EQ(Value<std::string>::downcast(s.getCacheEntry(Sub1, cx1InstanceTime)).get(), "hasPrereqs_Time");

    // Changing prerequisites should make the cache entry inaccessible again,
    // although the stage should not change and the cache entry should not get
    // deallocated.
    s.updZ() = 0.; // z is a prerequisite
    EXPECT_EQ(s.getSystemStage(), Stage::Time);
    EXPECT_TRUE(s.hasCacheEntry(CacheEntryKey(Sub1, cx1InstanceTime)));
    EXPECT_THROW(s.getCacheEntry(Sub1, cx1InstanceTime), std::exception);

    // Re-marking restores access.
    s.markCacheValueRealized(Sub1, cx1InstanceTime);
    EXPECT_EQ(Value<std::string>::downcast(s.getCacheEntry(Sub1, cx1InstanceTime)).get(), "hasPrereqs_Time");

    // Modify discrete-variable prerequisite — both entries must be invalidated.
    Value<int>::updDowncast(s.updDiscreteVariable(Sub0, dvx0ModelPos)).upd() = 99;
    EXPECT_EQ(s.getSystemStage(), Stage::Time);
    EXPECT_TRUE(s.hasCacheEntry(CacheEntryKey(Sub1, cx1InstanceTime)));
    EXPECT_THROW(s.getCacheEntry(Sub1, cx1InstanceTime), std::exception);
    EXPECT_THROW(s.getCacheEntry(Sub1, cx1InstInst), std::exception);

    s.markCacheValueRealized(Sub1, cx1InstanceTime);
    s.markCacheValueRealized(Sub1, cx1InstInst);

    // Invalidate upstream cache entry prerequisite (modifying value alone is
    // not enough — it must be explicitly marked not-realized).
    s.markCacheValueNotRealized(Sub0, cx0TopoModel);
    EXPECT_THROW(s.getCacheEntry(Sub1, cx1InstanceTime), std::exception);
    EXPECT_EQ(s.getSystemStage(), Stage::Time);
}

// -----------------------------------------------------------------------------
// Copy-constructing a State must faithfully reconstruct all dependency lists
// for z, discrete variables, and upstream cache entries.
// -----------------------------------------------------------------------------
TEST(SimTKCommon_State_CacheEntry, CopyConstructionPreservesDependencyLists) {
    const SubsystemIndex Sub0(0);
    const SubsystemIndex Sub1(1);

    State s;
    s.setNumSubsystems(2);

    s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    const CacheEntryIndex cx0TopoModel =
        s.allocateCacheEntry(Sub0, Stage::Model, Stage::Time, new Value<int>(41));

    advanceStage(s, Stage::Topology);

    const DiscreteVariableIndex dvx0ModelPos =
        s.allocateDiscreteVariable(Sub0, Stage::Position, new Value<int>(31));

    advanceStage(s, Stage::Model);

    const CacheEntryIndex cx1InstanceTime =
        s.allocateCacheEntryWithPrerequisites(Sub1,
                                              Stage::Time,
                                              Stage::Velocity,
                                              false,
                                              false,
                                              true,
                                              {DiscreteVarKey(Sub0, dvx0ModelPos)},
                                              {CacheEntryKey(Sub0, cx0TopoModel)},
                                              new Value<std::string>("hasPrereqs_Time"));

    const CacheEntryIndex cx1InstInst =
        s.allocateCacheEntryWithPrerequisites(Sub1,
                                              Stage::Instance,
                                              Stage::Velocity,
                                              false,
                                              false,
                                              true,
                                              {DiscreteVarKey(Sub0, dvx0ModelPos)},
                                              {CacheEntryKey(Sub0, cx0TopoModel)},
                                              new Value<std::string>("hasPrereqs_Instance"));

    advanceStage(s, Stage::Instance);
    advanceStage(s, Stage::Time);

    s.markCacheValueRealized(Sub0, cx0TopoModel);
    s.markCacheValueRealized(Sub1, cx1InstanceTime);
    s.markCacheValueRealized(Sub1, cx1InstInst);

    const CacheEntryKey ckey1(Sub1, cx1InstanceTime);
    const CacheEntryKey ckey2(Sub1, cx1InstInst);

    // Test state copying. Should copy through at least Instance stage.
    State s2(s); // copy construction

    EXPECT_EQ(s2.getNumSubsystems(), s.getNumSubsystems());
    EXPECT_GE(s2.getSystemStage(), Stage::Instance);
    EXPECT_TRUE(s2.hasCacheEntry(ckey1));
    EXPECT_TRUE(s2.hasCacheEntry(ckey2));

    // Dependency lists should have been reconstructed in the copy.
    EXPECT_EQ(s2.getZDependents().size(), 2u);
    EXPECT_TRUE(s2.getZDependents().contains(ckey1));
    EXPECT_TRUE(s2.getZDependents().contains(ckey2));

    const DiscreteVarInfo& dvinfo2 = s2.getDiscreteVarInfo(DiscreteVarKey(Sub0, dvx0ModelPos));
    EXPECT_EQ(dvinfo2.getDependents().size(), 2u);
    EXPECT_TRUE(dvinfo2.getDependents().contains(ckey1));
    EXPECT_TRUE(dvinfo2.getDependents().contains(ckey2));

    const CacheEntryInfo& ceinfo2 = s2.getCacheEntryInfo(CacheEntryKey(Sub0, cx0TopoModel));
    EXPECT_EQ(ceinfo2.getDependents().size(), 2u);
    EXPECT_TRUE(ceinfo2.getDependents().contains(ckey1));
    EXPECT_TRUE(ceinfo2.getDependents().contains(ckey2));

    // The copied state must independently allow marking and reading entries.
    s2.markCacheValueRealized(Sub1, cx1InstInst);
    EXPECT_EQ(Value<std::string>::downcast(s2.getCacheEntry(Sub1, cx1InstInst)).get(), "hasPrereqs_Instance");
}

// -----------------------------------------------------------------------------
// Modifying a Model-stage discrete variable drops the system to Topology,
// deallocates all Instance-stage cache entries, and removes them from their
// prerequisite dependency lists.
// -----------------------------------------------------------------------------
TEST(SimTKCommon_State_CacheEntry, ModelStageVariableChangeDeallocatesInstanceCache) {
    const SubsystemIndex Sub0(0);
    const SubsystemIndex Sub1(1);

    State s;
    s.setNumSubsystems(2);

    // Allocate at Topology stage a Model stage-invalidating state variable.
    const DiscreteVariableIndex dvx1TopoModel =
        s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    const CacheEntryIndex cx0TopoModel =
        s.allocateCacheEntry(Sub0, Stage::Model, Stage::Time, new Value<int>(41));

    advanceStage(s, Stage::Topology);

    const DiscreteVariableIndex dvx0ModelPos =
        s.allocateDiscreteVariable(Sub0, Stage::Position, new Value<int>(31));

    advanceStage(s, Stage::Model);

    const CacheEntryIndex cx1InstanceTime =
        s.allocateCacheEntryWithPrerequisites(Sub1,
                                              Stage::Time,
                                              Stage::Velocity,
                                              false,
                                              false,
                                              true,
                                              {DiscreteVarKey(Sub0, dvx0ModelPos)},
                                              {CacheEntryKey(Sub0, cx0TopoModel)},
                                              new Value<std::string>("hasPrereqs_Time"));

    advanceStage(s, Stage::Instance);
    advanceStage(s, Stage::Time);

    const CacheEntryKey ckey1(Sub1, cx1InstanceTime);

    // Capture a reference to the upstream cache entry's info so we can later
    // verify that its dependent list was cleared.
    const CacheEntryInfo& ceinfo = s.getCacheEntryInfo(CacheEntryKey(Sub0, cx0TopoModel));

    // Now modify the Model-stage variable again.
    Value<Real>::updDowncast(s.updDiscreteVariable(Sub1, dvx1TopoModel)) = 9;

    EXPECT_EQ(s.getSystemStage(), Stage::Topology);
    EXPECT_THROW(s.getCacheEntry(Sub0, cx0TopoModel), std::exception);

    // Unallocating the cache entry should have removed it from its
    // prerequisite's dependency list.
    EXPECT_FALSE(s.hasCacheEntry(ckey1));
    EXPECT_TRUE(s.getZDependents().empty());
    EXPECT_TRUE(ceinfo.getDependents().empty());
}

// -----------------------------------------------------------------------------
// After a Model-stage variable change and re-realization, a cache entry's
// value must be recomputed and explicitly stored. When the guaranteed-valid
// stage (Time) is subsequently reached, the entry becomes automatically
// accessible without calling markCacheValueRealized().
// -----------------------------------------------------------------------------
TEST(SimTKCommon_State_CacheEntry, CacheValueReflectsRecomputationAndAutoValidAtGuaranteedStage) {
    const SubsystemIndex Sub0(0);
    const SubsystemIndex Sub1(1);

    State s;
    s.setNumSubsystems(2);

    // Allocate at Topology stage a Model stage-invalidating state variable.
    const DiscreteVariableIndex dvx1TopoModel =
        s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    // Cache entry: depends-on Model, guaranteed-valid Time.
    const CacheEntryIndex cx0TopoModel =
        s.allocateCacheEntry(Sub0, Stage::Model, Stage::Time, new Value<int>(0));

    advanceStage(s, Stage::Topology);
    s.allocateDiscreteVariable(Sub0, Stage::Position, new Value<int>(31));
    advanceStage(s, Stage::Model);

    EXPECT_FALSE(s.isCacheValueRealized(Sub0, cx0TopoModel));
    EXPECT_THROW(s.getCacheEntry(Sub0, cx0TopoModel), std::exception);

    // "calculate" the cache entry (2 * dvx1TopoModel = 2 * 9) and mark it valid.
    Value<Real>::updDowncast(s.updDiscreteVariable(Sub1, dvx1TopoModel)) = 9;
    advanceStage(s, Stage::Model); // re-realize after variable change

    Value<int>::updDowncast(s.updCacheEntry(Sub0, cx0TopoModel)) =
        static_cast<int>(2 * Value<Real>::downcast(s.getDiscreteVariable(Sub1, dvx1TopoModel)));
    s.markCacheValueRealized(Sub0, cx0TopoModel);

    EXPECT_EQ(Value<int>::downcast(s.getCacheEntry(Sub0, cx0TopoModel)), 18);

    // Now modify the Model-stage variable again, but realize through Time stage.
    // We should be able to access the cache entry without explicitly marking it
    // valid (because Time is the guaranteed-valid stage).
    Value<Real>::updDowncast(s.updDiscreteVariable(Sub1, dvx1TopoModel)) = -100;
    advanceStage(s, Stage::Model);
    advanceStage(s, Stage::Instance);
    advanceStage(s, Stage::Time);

    // The previously computed value (18) is still present because nobody
    // re-computed it; what changed is only that the entry is now auto-valid.
    EXPECT_EQ(Value<int>::downcast(s.getCacheEntry(Sub0, cx0TopoModel)), 18);
}

// =============================================================================
// Domain 2 – Basic State operations
// =============================================================================

// Advancing to Stage::Topology must set t = 0.
TEST(SimTKCommon_State_Basics, TimeIsZeroAfterTopologyStage) {
    State s;
    s.setNumSubsystems(1);
    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Topology);
    s.advanceSystemToStage(Stage::Topology);

    EXPECT_EQ(s.getTime(), 0.0);
}

// Two successive Q allocations of sizes 3 and 2 must return indices 0 and 3.
TEST(SimTKCommon_State_Basics, QAllocationReturnsContiguousAscendingIndices) {
    State s;
    s.setNumSubsystems(1);
    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Topology);
    s.advanceSystemToStage(Stage::Topology);

    const QIndex q1 = s.allocateQ(SubsystemIndex(0), Vector(3));
    const QIndex q2 = s.allocateQ(SubsystemIndex(0), Vector(2));

    EXPECT_EQ(static_cast<int>(q1), 0);
    EXPECT_EQ(static_cast<int>(q2), 3);
}

// Event triggers allocated at different stages must have distinct,
// non-negative stage-local indices.
TEST(SimTKCommon_State_Basics, EventTriggerAllocationReturnsDistinctNonNegativeIndices) {
    State s;
    s.setNumSubsystems(1);
    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Topology);
    s.advanceSystemToStage(Stage::Topology);

    const EventTriggerByStageIndex e1 = s.allocateEventTrigger(SubsystemIndex(0), Stage::Position, 3);
    const EventTriggerByStageIndex e2 = s.allocateEventTrigger(SubsystemIndex(0), Stage::Instance, 2);

    EXPECT_GE(static_cast<int>(e1), 0);
    EXPECT_GE(static_cast<int>(e2), 0);
    EXPECT_NE(static_cast<int>(e1), static_cast<int>(e2));
}

// Writing a new value through updDiscreteVariable must be visible via
// getDiscreteVariable immediately afterwards.
TEST(SimTKCommon_State_Basics, DiscreteVariableUpdateRetainsNewValue) {
    State s;
    s.setNumSubsystems(1);
    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Topology);
    s.advanceSystemToStage(Stage::Topology);

    const DiscreteVariableIndex dv =
        s.allocateDiscreteVariable(SubsystemIndex(0), Stage::Dynamics, new Value<int>(5));

    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Model);
    s.advanceSystemToStage(Stage::Model);

    Value<int>::updDowncast(s.updDiscreteVariable(SubsystemIndex(0), dv)) = 71;

    EXPECT_EQ(Value<int>::downcast(s.getDiscreteVariable(SubsystemIndex(0), dv)), 71);
}

// =============================================================================
// Domain 3 – Inter-State consistency
// =============================================================================

// Helper: advance both states to Instance stage, then assert they ARE
// consistent in both directions.
auto expectConsistent(State& sA, State& sB) -> void {
    while (sA.getSystemStage() < Stage::Instance) {
        advanceStage(sA, sA.getSystemStage().next());
    }
    while (sB.getSystemStage() < Stage::Instance) {
        advanceStage(sB, sB.getSystemStage().next());
    }
    EXPECT_TRUE(sA.isConsistent(sB));
    EXPECT_TRUE(sB.isConsistent(sA));
}

// Helper: advance both states to Instance stage, then assert they are NOT
// consistent in both directions.
auto expectNotConsistent(State& sA, State& sB) -> void {
    while (sA.getSystemStage() < Stage::Instance) {
        advanceStage(sA, sA.getSystemStage().next());
    }
    while (sB.getSystemStage() < Stage::Instance) {
        advanceStage(sB, sB.getSystemStage().next());
    }
    EXPECT_FALSE(sA.isConsistent(sB));
    EXPECT_FALSE(sB.isConsistent(sA));
}

// Two freshly created states with the same subsystem count must be consistent
// once both reach Instance stage. Calling isConsistent() earlier must throw.
TEST(SimTKCommon_State_Consistency, IdenticalEmptyStatesAreConsistentAtInstanceStage) {
    State sA;
    State sB;

    sA.setNumSubsystems(3);
    sB.setNumSubsystems(3);

    // Must realize to Instance to check consistency.
    EXPECT_THROW(sA.isConsistent(sB), SimTK::Exception::StageTooLow);

    expectConsistent(sA, sB);
}

// States with a different number of subsystems must not be consistent;
// equalizing the count must restore consistency.
TEST(SimTKCommon_State_Consistency, DifferentSubsystemCountIsNotConsistent) {
    State sA;
    State sB;

    sA.setNumSubsystems(3);
    sB.setNumSubsystems(3);

    expectConsistent(sA, sB);

    sA.setNumSubsystems(4);
    expectNotConsistent(sA, sB);

    sB.setNumSubsystems(4);
    expectConsistent(sA, sB);
}

// Allocating Q, U, or Z with differing sizes in one state must break
// consistency; matching sizes must restore it.
TEST(SimTKCommon_State_Consistency, DifferentContinuousVectorSizesAreNotConsistent) {
    State sA;
    State sB;

    const int numSubsys = 3;
    sA.setNumSubsystems(numSubsys);
    sB.setNumSubsystems(numSubsys);

    expectConsistent(sA, sB);

    for (int i = 0; i < numSubsys; ++i) {
        // Helper lambda: reallocate in sA only → expect inconsistency;
        // then in sB with the same size → expect consistency.
        auto checkAllocPair = [&](auto allocFn, int sizePerSubsys) {
            sA.invalidateAll(Stage::Model);
            (sA.*allocFn)(SubsystemIndex(i), Vector(sizePerSubsys));
            expectNotConsistent(sA, sB);

            sB.invalidateAll(Stage::Model);
            (sB.*allocFn)(SubsystemIndex(i), Vector(sizePerSubsys));
            expectConsistent(sA, sB);
        };

        checkAllocPair(&State::allocateQ, (5 + i));
        checkAllocPair(&State::allocateU, (3 + i));
        checkAllocPair(&State::allocateZ, (8 + i));
    }
}

// Allocating QErr, UErr, or UDotErr with differing counts in one state must
// break consistency; matching counts must restore it.
TEST(SimTKCommon_State_Consistency, DifferentConstraintVectorSizesAreNotConsistent) {
    State sA;
    State sB;

    const int numSubsys = 3;
    sA.setNumSubsystems(numSubsys);
    sB.setNumSubsystems(numSubsys);

    expectConsistent(sA, sB);

    for (int i = 0; i < numSubsys; ++i) {
        auto checkAllocPair = [&](auto allocFn, int count) {
            sA.invalidateAll(Stage::Model);
            (sA.*allocFn)(SubsystemIndex(i), count);
            expectNotConsistent(sA, sB);

            sB.invalidateAll(Stage::Model);
            (sB.*allocFn)(SubsystemIndex(i), count);
            expectConsistent(sA, sB);
        };

        checkAllocPair(&State::allocateQErr, (2 + i));
        checkAllocPair(&State::allocateUErr, (7 + i));
        checkAllocPair(&State::allocateUDotErr, (1 + i));
    }
}

// Per-stage event-trigger counts participate in the consistency check:
// mismatched allocations must break consistency, and matching ones restore it.
TEST(SimTKCommon_State_Consistency, DifferentEventTriggerCountsAreNotConsistent) {
    State sA;
    State sB;

    const int numSubsys = 3;
    sA.setNumSubsystems(numSubsys);
    sB.setNumSubsystems(numSubsys);

    expectConsistent(sA, sB);

    for (int i = 0; i < numSubsys; ++i) {
        for (Stage stage = Stage::LowestValid; stage <= Stage::HighestRuntime; ++stage) {
            sA.invalidateAll(Stage::Model);
            sA.allocateEventTrigger(SubsystemIndex(i), stage, (12 + i));
            expectNotConsistent(sA, sB);

            sB.invalidateAll(Stage::Model);
            sB.allocateEventTrigger(SubsystemIndex(i), stage, (12 + i));
            expectConsistent(sA, sB);
        }
    }
}