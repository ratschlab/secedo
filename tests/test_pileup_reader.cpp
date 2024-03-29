#include "util/pileup_reader.hpp"

#include "util/util.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {
using namespace ::testing;

/**
 * Tests that the given pileup file contains the expected data.
 */
void check_binary(const std::string &fname,
                  uint32_t merge_count,
                  const std::string &group_fname,
                  const std::vector<PosData> &data,
                  uint16_t cell_count,
                  uint32_t max_len) {
    std::vector<PosData> data2;
    uint16_t cell_count2;
    uint32_t max_len2;
    std::tie(data2, cell_count2, max_len2)
            = read_pileup(fname, get_grouping(merge_count, group_fname));

    ASSERT_EQ(max_len2, max_len);
    ASSERT_EQ(cell_count2, cell_count);

    ASSERT_EQ(data2.size(), data.size());
    for (uint32_t i = 0; i < data.size(); ++i) {
        ASSERT_EQ(data2[i].position, data[i].position);
        ASSERT_EQ(data2[i].size(), data[i].size());
        for (uint32_t j = 0; j < data[i].read_ids.size(); ++j) {
            ASSERT_EQ(data2[i].base(j), data[i].base(j));
            ASSERT_EQ(data2[i].group_id(j), data[i].group_id(j));
        }
    }
    std::filesystem::remove_all(fname); // clean up
}

TEST(Reader, empty) {
    std::vector<PosData> data;
    uint16_t cell_count;
    uint32_t max_len;
    std::tie(data, cell_count, max_len) = read_pileup("data/empty.pileup", get_grouping());
    ASSERT_TRUE(data.empty());
    ASSERT_EQ(cell_count, 0);
}

TEST(Reader, one_row) {
    std::vector<PosData> data;
    uint16_t cell_count;
    uint32_t max_len;
    std::tie(data, cell_count, max_len) = read_pileup("data/one_row.pileup", get_grouping());
    std::vector<uint16_t> expected_cell_ids
            = { 95, 437, 458, 695, 887, 1011, 1216, 1223, 1522, 1612, 1795, 1924, 2163, 2163 };
    ASSERT_EQ(cell_count, 13);

    ASSERT_EQ(1, data.size());

    EXPECT_EQ(16050060, data[0].position);
    ASSERT_EQ(14, data[0].read_ids.size());
    ASSERT_EQ(14, data[0].group_ids_bases.size());

    std::vector<uint8_t> expected_bases
            = { 'C', 'A', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C' };

    for (uint32_t i = 0; i < data[0].size(); ++i) {
        ASSERT_EQ(data[0].read_ids[i], i);
        ASSERT_EQ(data[0].base(i), CharToInt[expected_bases[i]]);
        ASSERT_EQ(data[0].group_id(i), expected_cell_ids[i]);
    }

    check_binary("data/one_row.pileup.bin", 1, "", data, 2164, max_len);
}

TEST(Reader, three_rows) {
    std::vector<PosData> data;
    uint16_t cell_count;
    uint32_t max_len;
    std::tie(data, cell_count, max_len) = read_pileup("data/three_rows.pileup", get_grouping());
    std::vector<uint16_t> all_cell_ids = { 1, 2, 3, 4, 9 };
    ASSERT_EQ(cell_count, 5);

    ASSERT_EQ(3, data.size());

    std::vector<uint64_t> expected_positions = { 1, 2, 3 };
    std::vector<std::vector<uint16_t>> expected_cell_ids
            = { { 1, 3 }, { 1, 3, 2 }, { 1, 2, 3, 9, 4 } };
    std::vector<std::vector<uint8_t>> expected_bases
            = { { 'C', 'A' }, { 'A', 'C', 'T' }, { 'A', 'A', 'A', 'C', 'A' } };
    std::vector<std::vector<uint32_t>> expected_read_ids
            = { { 0, 1 }, { 0, 1, 2 }, { 0, 2, 1, 3, 4 } };

    for (uint32_t i : { 0, 1, 2 }) {
        ASSERT_EQ(expected_positions[i], data[i].position);

        for (uint32_t j = 0; j < data[i].size(); ++j) {
            ASSERT_EQ(data[i].read_ids[j], expected_read_ids[i][j]);
            ASSERT_EQ(data[i].base(j), CharToInt[expected_bases[i][j]]);
            ASSERT_EQ(data[i].group_id(j), expected_cell_ids[i][j]);
        }
    }

    check_binary("data/three_rows.pileup.bin", 1, "", data, 10, max_len);
}

TEST(Reader, filter_positions) {
    std::vector<PosData> data;
    uint16_t cell_count;
    uint32_t max_len;
    std::tie(data, cell_count, max_len) = read_pileup("data/three_rows.pileup", get_grouping(),
                                                      [](uint32_t) {}, 100U, { 0, 2, 9, 11 });
    std::vector<uint16_t> all_cell_ids = { 1, 2, 3 };
    ASSERT_EQ(cell_count, 3);

    ASSERT_EQ(1, data.size());

    std::vector<std::vector<uint16_t>> expected_cell_ids = { { 1, 3, 2 } };
    std::vector<std::vector<uint8_t>> expected_bases = { { 'A', 'C', 'T' } };
    std::vector<std::vector<uint32_t>> expected_read_ids = { { 0, 1, 2 } };

    ASSERT_EQ(2, data[0].position);

    for (uint32_t j = 0; j < data[0].size(); ++j) {
        ASSERT_EQ(data[0].read_ids[j], expected_read_ids[0][j]);
        ASSERT_EQ(data[0].base(j), CharToInt[expected_bases[0][j]]);
        ASSERT_EQ(data[0].group_id(j), expected_cell_ids[0][j]);
    }

    check_binary("data/three_rows.pileup.bin", 1, "", data, 4, max_len);
}

/**
 * Reads a pileup file containing data about 6 cells, numbered 0..5. Data from 2 consecutive cells
 * is merged, resulting in 3 effective cells.
 */
TEST(Reader, group_by_two) {
    std::vector<PosData> data;
    uint16_t cell_count;
    uint32_t max_len;
    std::tie(data, cell_count, max_len) = read_pileup("data/six_cells.pileup", get_grouping(2));
    std::vector<uint16_t> all_cell_ids = { 0, 1, 2, 3, 4, 5 };
    ASSERT_EQ(cell_count, 6);

    ASSERT_EQ(3, data.size());

    std::vector<uint64_t> expected_positions = { 1, 19, 103 };

    for (uint32_t i : { 0, 1, 2 }) {
        ASSERT_EQ(expected_positions[i], data[i].position);


        std::vector<std::vector<uint16_t>> expected_cell_ids
                = { { 0, 0, 1, 1, 2 }, { 0, 0, 1, 2 }, { 2, 2, 0, 1, 1, 0 } };
        std::vector<std::vector<uint8_t>> expected_bases = { { 'A', 'A', 'A', 'A', 'C' },
                                                             { 'A', 'C', 'G', 'T' },
                                                             { 'T', 'G', 'C', 'A', 'T', 'A' } };
        std::vector<std::vector<uint32_t>> expected_read_ids
                = { { 0, 1, 2, 3, 4 }, { 0, 1, 2, 3 }, { 0, 1, 2, 3, 4, 5 } };

        for (uint32_t j = 0; j < data[i].size(); ++j) {
            ASSERT_EQ(data[i].read_ids[j], expected_read_ids[i][j]);
            ASSERT_EQ(data[i].base(j), CharToInt[expected_bases[i][j]]);
            ASSERT_EQ(data[i].group_id(j), expected_cell_ids[i][j]);
        }
    }

    check_binary("data/six_cells.pileup.bin", 2, "", data, cell_count, max_len);
}

/**
 * Reads a pileup file containing data about 6 cells and a grouping file that merges the cells into
 * 2 groups of 3 cells each.
 */
TEST(Reader, group_using_file) {
    std::vector<PosData> data;
    uint16_t cell_count;
    uint32_t max_len;
    std::tie(data, cell_count, max_len)
            = read_pileup("data/six_cells.pileup", get_grouping(1, "data/six_cells.pileup.group"));
    std::vector<uint16_t> all_cell_ids = { 0, 1, 2, 3, 4, 5 };
    ASSERT_EQ(cell_count, 6);

    ASSERT_EQ(3, data.size());

    std::vector<uint64_t> expected_positions = { 1, 19, 103 };

    for (uint32_t i : { 0, 1, 2 }) {
        ASSERT_EQ(expected_positions[i], data[i].position);


        std::vector<std::vector<uint16_t>> expected_cell_ids
                = { { 0, 0, 0, 1, 1 }, { 0, 0, 1, 1 }, { 1, 1, 0, 1, 0, 0 } };
        std::vector<std::vector<uint8_t>> expected_bases = { { 'A', 'A', 'A', 'A', 'C' },
                                                             { 'A', 'C', 'G', 'T' },
                                                             { 'T', 'G', 'C', 'A', 'T', 'A' } };
        std::vector<std::vector<uint32_t>> expected_read_ids
                = { { 0, 1, 2, 3, 4 }, { 0, 1, 2, 3 }, { 0, 1, 2, 3, 4, 5 } };

        for (uint32_t j = 0; j < data[i].size(); ++j) {
            ASSERT_EQ(data[i].read_ids[j], expected_read_ids[i][j]);
            ASSERT_EQ(data[i].base(j), CharToInt[expected_bases[i][j]]);
            ASSERT_EQ(data[i].group_id(j), expected_cell_ids[i][j]);
        }
    }

    check_binary("data/six_cells.pileup.bin", 1, "data/six_cells.pileup.group", data, cell_count,
                 max_len);
}

} // namespace
