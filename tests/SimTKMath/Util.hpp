#pragma once

#include <gtest/gtest.h>
#include <vector>

#include "SimTKmath.h"

template <class T>
void assertEqual(T val1, T val2, const SimTK::Real TOL = 1e-4) {
    const T sz = std::max(std::abs(val1), std::abs(val2));
    const SimTK::Real tol = std::max(TOL, sz * TOL);
    EXPECT_NEAR(val1, val2, tol);
}

template <int N>
void assertEqual(const SimTK::Vec<N>& val1, const SimTK::Vec<N>& val2, const SimTK::Real TOL = 1e-4) {
    for (int i = 0; i < N; ++i) {
        const SimTK::Real sz = std::max(std::abs(val1[i]), std::abs(val2[i]));
        const SimTK::Real tol = std::max(TOL, sz * TOL);
        EXPECT_NEAR(val1[i], val2[i], tol);
    }
}

inline void assertEqual(const SimTK::UnitVec3& v1, const SimTK::UnitVec3& v2) {
    assertEqual(v1.asVec3(), v2.asVec3());
}

/**
* This function computes a standard central difference dy/dx.
* If extrap_endpoints is set to 1, then the derivative at the
* end points is estimated by linearly extrapolating the dy/dx
* values beside the end points
*
* @param x domain vector
* @param y range vector
& @param extrap_endpoints:
*   (false)   Endpoints of the returned vector will be zero,
*             because a central difference is undefined at
*             these endpoints
*    (true)   Endpoints are computed by linearly extrapolating
*             using a first difference from the neighboring 2
*             points
*
* @returns dy/dx computed using central differences
*/
inline auto getCentralDifference(const SimTK::Vector& x, const SimTK::Vector& y, bool extrap_endpoints)
    -> SimTK::Vector {
    SimTK::Vector dy(x.size());
    double dx1;
    double dx2;
    double dy1;
    double dy2;
    int size = x.size();
    for (int i = 1; i < x.size() - 1; i++) {
        dx1 = x(i) - x(i - 1);
        dx2 = x(i + 1) - x(i);
        dy1 = y(i) - y(i - 1);
        dy2 = y(i + 1) - y(i);
        dy(i) = (0.5 * dy1 / dx1) + (0.5 * dy2 / dx2);
    }

    if (extrap_endpoints) {
        dy1 = dy(2) - dy(1);
        dx1 = x(2) - x(1);
        dy(0) = dy(1) + ((dy1 / dx1) * (x(0) - x(1)));

        dy2 = dy(size - 2) - dy(size - 3);
        dx2 = x(size - 2) - x(size - 3);
        dy(size - 1) = dy(size - 2) + ((dy2 / dx2) * (x(size - 1) - x(size - 2)));
    }

    return dy;
}
