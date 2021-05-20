#include "util/util.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {
using namespace ::testing;

TEST(MatReadWrite, Empty) {
    Matd empty;
    write_mat("test", empty);
    Matd res = read_mat("test");
    ASSERT_TRUE(res.empty());

    std::filesystem::remove("test");
}

TEST(MatReadWrite, Idempotent) {
    constexpr uint32_t size = 10;
    Matd mat(size, size);
    for (uint32_t i = 0; i < size; ++i) {
        for (uint32_t j = 0; j < size; ++j) {
            mat(i, j) = sqrt(i + j);
        }
    }
    write_mat("test", mat);

    Matd mat2 = read_mat("test");

    for (uint32_t i = 0; i < size; ++i) {
        for (uint32_t j = 0; j < size; ++j) {
            ASSERT_NEAR(mat(i, j), sqrt(i + j), 1e-5);
        }
    }
}

TEST(ReadPositions, SomeData) {
    std::vector<std::vector<uint32_t>> positions = read_positions("data/cosmic_small.vcf");
    std::vector<std::vector<uint32_t>> expected_pos
            = { { 65797, 66373, 66814, 69224, 69345, 69486, 69518, 69538, 69559, 69569, 69598,
                  69666 },
                {},
                {},
                {},
                {},
                {},
                {},
                {},
                {},
                { 7197537, 7197566, 7197575, 7197584, 7197593, 7197616 },
                {
                        2133595,
                        2133597,
                        2133605,
                        2133606,
                } };
    ASSERT_EQ(11, positions.size());
    ASSERT_EQ(positions, expected_pos);
}

} // namespace
