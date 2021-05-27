#include "util/is_significant.hpp"

#include "util/util.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {
using namespace ::testing;

bool is_significant_helper(const std::string &bases, double theta = 0.01) {
    Filter filter(theta);
    std::array<uint16_t, 4> base_count = { 0, 0, 0, 0 };
    for (uint8_t c : bases) {
        base_count[CharToInt[c]]++;
    }
    return filter.is_significant(base_count);
}

TEST(LogFact, Zero) {
    Filter filter(0.01);
    ASSERT_EQ(0, filter.log_fact(0));
}

TEST(LogFact, One) {
    Filter filter(0.01);
    ASSERT_EQ(0, filter.log_fact(1));
}

TEST(LogFact, SomeValues) {
    Filter filter(0.01);
    double values[] = { 0.69314718056, 1.79175946923, 3.17805383035, 4.78749174278 };
    for (uint32_t i = 0; i < 4; ++i) {
        ASSERT_NEAR(values[i], filter.log_fact(i + 2), 1e-10);
    }
}

TEST(LogFact, LargeValues) {
    Filter filter(0.01);
    double values[] = { 716.86, 722.01, 727.17, 732.33 };
    for (uint32_t i = 0; i < 4; ++i) {
        ASSERT_NEAR(values[i], filter.log_fact(i + 172), 1e-2);
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
    // this is now considered insignificant, as 2 bases are not "convincing" enough
    ASSERT_FALSE(is_significant_helper(bases, 0.001));
}

//TEST(Preprocess, TwoSigmasAwayTrue) {
//    std::string bases
//            = "GGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
//              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
//    ASSERT_TRUE(is_significant_helper(bases, 0.01));
//}

TEST(Preprocess, TwoSigmasAwayFalse) {
    std::string bases
            = "GGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    ASSERT_FALSE(is_significant_helper(bases, 0.01));
}

TEST(Preprocess, Paradox) {
    std::string bases = "GGGA";
    ASSERT_FALSE(is_significant_helper(bases, 0.001));
}


// this test is removed, as we stopped emulating the python version after introducing
// is_two_sigmas_away()
/** Make sure we round the length of the string properly (for computing the proper threshold K),
 * i.e. length 55 in this case should generate index 5 (because we emulate rounding to nearest even
 */
// TEST(Preprocess, LengthAtLimitRoundUp) {
//    std::string bases = "GGGGGGGGGGGGGTGGGGGGGGGGGAGGGGGGGGGGGGGGGTGGGGGGGGGGGGG";
//    ASSERT_FALSE(is_significant_helper(bases, 0.001));
//}

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

TEST(Filter, Empty) {
    Filter filter(1e-3);
    auto [pos_data, coverage] = filter.filter({}, {}, {}, "", 1);
    ASSERT_TRUE(pos_data.empty());
}

TEST(Filter, OnePosSignificant) {
    // Filtering is pretty drastic in order to increase signal to noise ratio
    // Even with a seq error rate of 1e-3 we need at least 10 out of 100 bases to be different in
    // order to consider a position significant. This does mean that small clusters won't be seen
    uint32_t coverage = 100;
    uint32_t num_diff = 10; // bases different from the others

    std::vector<uint32_t> read_ids;
    std::vector<uint8_t> bases;
    std::vector<uint16_t> cell_ids;

    // building a position like: AAAAACCCCCCCCCCCC...
    for (uint32_t i = 0; i < coverage; ++i) {
        read_ids.push_back(2*i+1);
        if (i < num_diff) {
            bases.push_back(0);
        } else {
            bases.push_back(1);
        }
        cell_ids.push_back(i);
    }

    PosData pd = assemble(1, read_ids, cell_ids, bases);

    std::vector<uint16_t> id_to_group(coverage);
    std::iota(id_to_group.begin(), id_to_group.end(), 0);
    std::vector<uint32_t> id_to_pos(coverage);
    std::iota(id_to_pos.begin(), id_to_pos.end(), 0);

    Filter filter(1e-3);
    auto [filtered, avg_coverage] = filter.filter({ { pd } }, id_to_group, id_to_pos, "", 1);
    ASSERT_EQ(avg_coverage, 100.0);
    std::vector<PosData> chromosome_data = { pd };
    ASSERT_THAT(filtered, ElementsAre(chromosome_data));
}

TEST(Filter, OnePosNotSignificant) {
    std::vector<uint32_t> read_ids = { 0, 5, 9 };
    std::vector<uint8_t> bases = { 0, 0, 0 }; // all bases are the same
    std::vector<uint16_t> cell_ids = { 0, 1, 2 };
    PosData pd = assemble(1, read_ids, cell_ids, bases);

    Filter filter(1e-3);
    auto [filtered, coverage] = filter.filter({ { pd } }, { 0, 1, 2 }, { 0, 1, 2 }, "", 1);
    ASSERT_EQ(coverage, 0);
    ASSERT_EQ(1, filtered.size());
    ASSERT_TRUE(filtered[0].empty());
}

TEST(Filter, AllSignificant) {
    logger()->set_level(spdlog::level::trace);
    std::vector<std::vector<PosData>> pos_data;
    uint32_t coverage = 100;
    uint32_t num_diff = 10; // bases different from the others

    for (uint32_t chr = 0; chr < 23; ++chr) {
        std::vector<uint32_t> read_ids;
        std::vector<uint8_t> bases;
        std::vector<uint16_t> cell_ids;

        // building a position like: AAAAACCCCCCCCCCCC...
        for (uint32_t i = 0; i < coverage; ++i) {
            read_ids.push_back(2*i+1);
            if (i < num_diff) {
                bases.push_back(0);
            } else {
                bases.push_back(1);
            }
            cell_ids.push_back(i);
        }

        std::vector<PosData> chromosome_data;
        for (uint32_t i = 0; i < 100; ++i) {
            chromosome_data.push_back(assemble(i + 1, read_ids, cell_ids, bases));
        }
        pos_data.push_back(chromosome_data);
    }
    std::vector<uint16_t> id_to_group(coverage);
    std::iota(id_to_group.begin(), id_to_group.end(), 0);
    std::vector<uint32_t> id_to_pos(coverage);
    std::iota(id_to_pos.begin(), id_to_pos.end(), 0);

    std::vector<std::vector<PosData>> filtered;
    double avg_coverage;
    Filter filter(1e-3);
    std::tie(filtered, avg_coverage) = filter.filter(pos_data, id_to_group, id_to_pos, "", 2);
    ASSERT_EQ(avg_coverage, coverage);
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
    std::vector<uint16_t> id_to_group(10);
    std::iota(id_to_group.begin(), id_to_group.end(), 0);
    std::vector<uint32_t> id_to_pos(10);
    std::iota(id_to_pos.begin(), id_to_pos.end(), 0);

    Filter filter(1e-3);
    auto [filtered, coverage] = filter.filter({ chromosome_data }, id_to_group, id_to_pos, "", 2);
    ASSERT_EQ(coverage, 0);
    ASSERT_EQ(1, filtered.size());
    ASSERT_TRUE(filtered[0].empty());
}

} // namespace
