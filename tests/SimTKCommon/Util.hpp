#pragma once

#include <gtest/gtest.h>
#include <string>

#include "Simbody.h"

using namespace SimTK;

// Helper to perform element-wise numerical comparison with detailed error logging
template <typename T1, typename T2>
auto AssertSimTKEqual(const std::string& m_expr,
                      const std::string& expected_expr,
                      const T1& actual,
                      const T2& expected) -> ::testing::AssertionResult {
    if (SimTK::Test::numericallyEqual(actual, expected, 1)) {
        return ::testing::AssertionSuccess();
    }

    return ::testing::AssertionFailure() << "Value of: " << m_expr << "\n"
                                         << "  Actual: " << actual << "\n"
                                         << "Expected: " << expected_expr << "\n"
                                         << "  Which is: " << expected;
}
