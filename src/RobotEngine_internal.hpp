#pragma once

#include "robot_math.hpp"

namespace robo::detail {

// Spatial dot product:  ~a * b  ==  a.angular . b.angular + a.linear . b.linear.
inline Real spatialDot(const SpatialVec& a, const SpatialVec& b) {
    return dot(a[0], b[0]) + dot(a[1], b[1]);
}

} // namespace robo::detail
