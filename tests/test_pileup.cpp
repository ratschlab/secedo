#include "pileup.hpp"
#include "util/pileup_reader.hpp"

#include <gtest/gtest.h>

namespace {

TEST(pileup, read) {
    uint32_t chromosome_id = 0;
    uint32_t max_coverage = 10;
    std::vector<PosData> data
            = pileup_bams({ "data/test1.bam", "data/test2.bam" }, "data/test_pileup", true,
                          chromosome_id, max_coverage, 1, 1);
    ASSERT_EQ(10, data.size());
    for (uint32_t i = 0; i < 5; ++i) {
        ASSERT_EQ(2, data[i].size());
        ASSERT_EQ(0, data[i].read_ids[0]);
        ASSERT_EQ(1, data[i].read_ids[1]);
        ASSERT_EQ(0, data[i].base(0));
        ASSERT_EQ(2, data[i].base(1));
        ASSERT_EQ(0, data[i].cell_id(0));
        ASSERT_EQ(1, data[i].cell_id(1));
    }
    for (uint32_t i = 5; i < 10; ++i) {
        ASSERT_EQ(2, data[i].size());
        ASSERT_EQ(0, data[i].read_ids[0]);
        ASSERT_EQ(1, data[i].read_ids[1]);
        ASSERT_EQ(1, data[i].base(0));
        ASSERT_EQ(3, data[i].base(1));
        ASSERT_EQ(0, data[i].cell_id(0));
        ASSERT_EQ(1, data[i].cell_id(1));
    }
}

TEST(pileup, read_file) {
    uint32_t chromosome_id = 0;
    uint32_t max_coverage = 10;
    pileup_bams({ "data/test1.bam", "data/test2.bam" }, "data/test_pileup", true, chromosome_id,
                max_coverage, 1, 1);
    auto [data, cell_ids, max_len] = read_pileup("data/test_pileup.bin", {0,1});
    ASSERT_EQ(424, max_len);
    ASSERT_EQ(10, data.size());
    for (uint32_t i = 0; i < 5; ++i) {
        ASSERT_EQ(2, data[i].size());
        ASSERT_EQ(0, data[i].read_ids[0]);
        ASSERT_EQ(1, data[i].read_ids[1]);
        ASSERT_EQ(0, data[i].base(0));
        ASSERT_EQ(2, data[i].base(1));
        ASSERT_EQ(0, data[i].cell_id(0));
        ASSERT_EQ(1, data[i].cell_id(1));
    }
    for (uint32_t i = 5; i < 10; ++i) {
        ASSERT_EQ(2, data[i].size());
        ASSERT_EQ(0, data[i].read_ids[0]);
        ASSERT_EQ(1, data[i].read_ids[1]);
        ASSERT_EQ(1, data[i].base(0));
        ASSERT_EQ(3, data[i].base(1));
        ASSERT_EQ(0, data[i].cell_id(0));
        ASSERT_EQ(1, data[i].cell_id(1));
    }
}


} // namespace
