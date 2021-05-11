#include "util/is_significant.hpp"

#include "util/util.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {
using namespace ::testing;

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

TEST(Filter, Empty) {
    auto [pos_data, coverage] = filter({}, {}, {}, "", {}, 1);
    ASSERT_TRUE(pos_data.empty());
}

TEST(Filter, OnePosSignificant) {
    std::vector<uint32_t> read_ids = { 0, 5, 9 };
    std::vector<uint16_t> cell_ids_and_bases = { 1 << 2 | 0, 2 << 2 | 1, 3 << 2 | 2 };
    PosData pd = { 1, read_ids, cell_ids_and_bases };
    auto [filtered, coverage] = filter({ { pd } }, { 0 }, { 0 }, "", 1e-3, 1);
    ASSERT_EQ(coverage, 3.0);
    std::vector<PosData> chromosome_data = { pd };
    ASSERT_THAT(filtered, ElementsAre(chromosome_data));
}

PosData assemble(uint32_t pos,
                 std::vector<uint32_t> read_ids,
                 std::vector<uint16_t> cell_ids,
                 std::vector<uint8_t> bases) {
    std::vector<uint16_t> cell_ids_and_bases(bases.size());
    for (uint32_t i = 0; i < bases.size(); ++i) {
        cell_ids_and_bases[i] = cell_ids[i] << 2 | bases[i];
    }
    return { pos, read_ids, cell_ids_and_bases };
}

TEST(Filter, OnePosNotSignificant) {
    std::vector<uint32_t> read_ids = { 0, 5, 9 };
    std::vector<uint8_t> bases = { 0, 0, 0 }; // all bases are the same
    std::vector<uint16_t> cell_ids = { 1, 3, 5 };
    PosData pd = assemble(1, read_ids, cell_ids, bases);

    auto [filtered, coverage] = filter({ { pd } }, { 0 }, { 0 }, "", 1e-3, 1);
    ASSERT_EQ(coverage, 0);
    ASSERT_EQ(1, filtered.size());
    ASSERT_TRUE(filtered[0].empty());
}

TEST(Filter, AllSignificant) {
    std::vector<std::vector<PosData>> pos_data;
    for (uint32_t chr = 0; chr < 23; ++chr) {
        std::vector<uint32_t> read_ids = { 0, 5, 9 };
        std::vector<uint8_t> bases = { 0, 1, 2 }; // all bases are different, positions are kept
        std::vector<uint16_t> cell_ids = { 1, 3, 5 };
        std::vector<PosData> chromosome_data;
        for (uint32_t i = 0; i < 100; ++i) {
            chromosome_data.push_back(assemble(i + 1, read_ids, cell_ids, bases));
        }
        pos_data.push_back(chromosome_data);
    }
    std::vector<uint16_t> id_to_group(100);
    std::iota(id_to_group.begin(), id_to_group.end(), 0);
    std::vector<uint32_t> id_to_pos(100);
    std::iota(id_to_pos.begin(), id_to_pos.end(), 0);

    std::vector<std::vector<PosData>> filtered;
    double coverage;
    std::tie(filtered, coverage) = filter(pos_data, id_to_group, id_to_pos, "", 1e-3, 2);
    ASSERT_EQ(coverage, 3);
    ASSERT_EQ(23, filtered.size());
    ASSERT_EQ(filtered, pos_data);
}

TEST(Filter, NoneSignificant) {
    std::vector<uint32_t> read_ids = { 0, 5, 9 };
    std::vector<uint8_t> bases = { 0, 0, 0 }; // all bases are the same, position is removed
    std::vector<uint16_t> cell_ids = { 1, 3, 5 };
    std::vector<PosData> chromosome_data;
    for (uint32_t i = 0; i < 100; ++i) {
        chromosome_data.push_back(assemble(i + 1, read_ids, cell_ids, bases));
    }
    std::vector<uint16_t> id_to_group(100);
    std::iota(id_to_group.begin(), id_to_group.end(), 0);
    std::vector<uint32_t> id_to_pos(100);
    std::iota(id_to_pos.begin(), id_to_pos.end(), 0);

    auto [filtered, coverage] = filter({ chromosome_data }, id_to_group, id_to_pos, "", 1e-3, 2);
    ASSERT_EQ(coverage, 0);
    ASSERT_EQ(1, filtered.size());
    ASSERT_TRUE(filtered[0].empty());
}

} // namespace
