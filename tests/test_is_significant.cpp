#include "util/is_significant.hpp"

#include "util/util.hpp"

#include <gtest/gtest.h>

namespace {

bool is_significant_helper(const std::string &bases, double theta = 0.01) {
    std::array<uint16_t, 4> base_count = { 0, 0, 0, 0 };
    for (uint8_t c : bases) {
        base_count[CharToInt[c]]++;
    }
    return is_significant(base_count, theta);
}

TEST(LogFact, Zero) {
    ASSERT_EQ(0, log_fact(0));
}

TEST(LogFact, One) {
    ASSERT_EQ(0, log_fact(1));
}

TEST(LogFact, SomeValues) {
    double values[] = { 0.69314718056, 1.79175946923, 3.17805383035, 4.78749174278 };
    for (uint32_t i = 0; i < 4; ++i) {
        ASSERT_NEAR(values[i], log_fact(i + 2), 1e-10);
    }
}

TEST(LogFact, LargeValues) {
    double values[] = { 716.86, 722.01, 727.17, 732.33 };
    for (uint32_t i = 0; i < 4; ++i) {
        ASSERT_NEAR(values[i], log_fact(i + 172), 1e-2);
    }
}

TEST(Preprocess, Cov52OneDifferent) {
    std::string bases = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCC";
    ASSERT_FALSE(is_significant_helper(bases));
}

TEST(Preprocess, Cov52TenDifferent) {
    std::string bases = "CCACGTACGTACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCC";
    ASSERT_TRUE(is_significant_helper(bases));
}

TEST(Preprocess, Cov59TwoDifferent) {
    std::string bases = "tttttTTTTaTTTttTaTtTTTTtTTtTTTtTttTTtTtTtttTTttttTTttTTTTtt";
    ASSERT_FALSE(is_significant_helper(bases));
}

TEST(Preprocess, AtLimit) {
    std::string bases = "CcccccccCcccCCCCcaCcccCccACccccCCCcCcCCCC";
    ASSERT_TRUE(is_significant_helper(bases, 0.001));
}

/** Make sure we round the length of the string properly (for computing the proper threshold K),
 * i.e. length 65 in this case should generate index 5 */
TEST(Preprocess, LengthAtLimitRoundDown) {
    std::string bases = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    ASSERT_TRUE(is_significant_helper(bases, 0.001));
}

/** Make sure we round the length of the string properly (for computing the proper threshold K),
 * i.e. length 55 in this case should generate index 5 (because we emulate rounding to nearest even
 */
TEST(Preprocess, LengthAtLimitRoundUp) {
    std::string bases = "GGGGGGGGGGGGGTGGGGGGGGGGGAGGGGGGGGGGGGGGGTGGGGGGGGGGGGG";
    ASSERT_FALSE(is_significant_helper(bases, 0.001));
}

} // namespace
