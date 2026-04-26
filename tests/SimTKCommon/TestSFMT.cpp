#include <array>
#include <cstdint>
#include <gtest/gtest.h>
#include <memory>
#include <vector>

// Assuming these are the headers for the SFMT library used in Simbody
#include "../../Random/src/SFMT.h"

using namespace SimTK_SFMT;

class SimTKCommon_SFTM : public ::testing::Test {
    protected:
    // Using unique_ptr with a custom deleter for automatic cleanup
    struct SFMTDeleter {
        void operator()(SFMTData* p) const {
            deleteSFMTData(p);
        }
    };

    std::unique_ptr<SFMTData, SFMTDeleter> sfmt{createSFMTData()};

    // Expected values constant across tests
    static constexpr std::array<uint32_t, 5> expected = {3440181298U,
                                                         1564997079U,
                                                         1510669302U,
                                                         2930277156U,
                                                         1452439940U};
};

// Test Case 1: Sequential generation verification
TEST_F(SimTKCommon_SFTM, GeneratesCorrectSequentialValues) {
    init_gen_rand(1234, *sfmt);

    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(gen_rand32(*sfmt), expected[i]) << "Mismatch at index " << i;
    }
}

// Test Case 2: Bulk array filling verification
TEST_F(SimTKCommon_SFTM, MatchesExpectedValuesWhenFillingArray) {
    const int min_size = get_min_array_size32();
    const int length = std::max(5, min_size);

    std::vector<uint32_t> values(length);

    init_gen_rand(1234, *sfmt);
    fill_array32(values.data(), length, *sfmt);

    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(values[i], expected[i]) << "Array mismatch at index " << i;
    }
}
