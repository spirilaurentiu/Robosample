#pragma once

#include <array>
#include <gtest/gtest.h>

#include "Simbody.h"
using namespace SimTK;

#define EXPECT_SIMTK_SIZE(expected, actual, n)                            \
    EXPECT_TRUE(SimTK::Test::numericallyEqual((expected), (actual), (n))) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\nSize:      " << (n) << "\n"

#define EXPECT_NEAR_CUSTOM_TOL_SIMTK(expected, actual, tol)                              \
    EXPECT_TRUE(SimTK::Test::numericallyEqual((expected), (actual), /*scale*/ 1, (tol))) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\nTolerance: " << (tol) << "\n"

#define EXPECT_NEAR_DEFAULT_TOL_SIMTK(expected, actual) \
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((expected), (actual), 1e-6)

#define EXPECT_NOT_NEAR_DEFAULT_TOL_SIMTK(expected, actual)                        \
    EXPECT_FALSE(SimTK::Test::numericallyEqual((expected), (actual), /*scale*/ 1)) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\n"


// Reaction force information: results are the joint reaction, at the F frame
// on the parent and M frame on the child, expressed in the parent or child
// frame, resp. Note that Gazebo's GetForceTorque() method uses the negation of
// the joint reaction, and Gazebo's results are ordered (force,torque) rather
// than (torque,force) as in a Simbody SpatialVec.
struct ReactionPair {
    SpatialVec reactionOnParentInParent;
    SpatialVec reactionOnChildInChild;
};

static inline auto getReactionPair(const State& state, const MobilizedBody& mobod) -> ReactionPair {
    const SpatialVec parent = mobod.findMobilizerReactionOnParentAtFInGround(state);
    const SpatialVec child = mobod.findMobilizerReactionOnBodyAtMInGround(state);

    const Rotation& R_GC = mobod.getBodyRotation(state);
    const Rotation& R_GP = mobod.getParentMobilizedBody().getBodyRotation(state);

    ReactionPair pair;
    pair.reactionOnChildInChild = ~R_GC * child;    // from Ground to Child
    pair.reactionOnParentInParent = ~R_GP * parent; // from Ground to Parent

    return pair;
}
