#include "util/util.hpp"

#include <gtest/gtest.h>

namespace {
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
} // namespace
