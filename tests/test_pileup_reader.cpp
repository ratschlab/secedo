#include "pileup_reader.hpp"

#include "util.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {
using namespace ::testing;

TEST(Reader, empty) {
    std::vector<PosData> data;
    std::unordered_set<uint32_t> read_ids;
    std::tie(data, read_ids) = read_pileup("data/empty.pileup");
    ASSERT_TRUE(data.empty());
    ASSERT_TRUE(read_ids.empty());
}

TEST(Reader, one_row) {
    std::vector<PosData> data;
    std::unordered_set<uint32_t> read_ids;
    std::tie(data, read_ids) = read_pileup("data/one_row.pileup");
    std::vector<uint16_t> expected_cell_ids
            = { 95, 437, 458, 695, 887, 1011, 1216, 1223, 1522, 1612, 1795, 1924, 2163, 2163 };
    std::vector<uint16_t> distinct_cell_ids
            = { 95, 437, 458, 695, 887, 1011, 1216, 1223, 1522, 1612, 1795, 1924, 2163 };
    ASSERT_THAT(read_ids, UnorderedElementsAreArray(distinct_cell_ids));

    ASSERT_EQ(1, data.size());

    EXPECT_EQ(16050060, data[0].position);
    ASSERT_EQ(14, data[0].cells_data.size());

    std::vector<uint8_t> expected_bases
            = { 'C', 'A', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C' };
    std::vector<std::string> expected_read_ids = {
        "A00228:202:HCCGKDMXX:2:1121:20048:18865", "A00228:202:HCCGKDMXX:2:2220:19904:9940",
        "A00228:202:HCCGKDMXX:2:1147:11198:23954", "A00228:202:HCCGKDMXX:1:1334:10791:27007",
        "A00228:202:HCCGKDMXX:2:1456:17065:20682", "A00228:202:HCCGKDMXX:2:2284:17833:6887",
        "A00228:202:HCCGKDMXX:2:2460:7771:12383",  "A00228:202:HCCGKDMXX:2:2336:29179:10535",
        "A00228:202:HCCGKDMXX:1:1273:30364:7451",  "A00228:202:HCCGKDMXX:2:1337:22797:15452",
        "A00228:202:HCCGKDMXX:2:1369:32425:18474", "A00228:202:HCCGKDMXX:2:2464:15926:5776",
        "A00228:202:HCCGKDMXX:2:2469:6451:23312",  "A00228:202:HCCGKDMXX:1:1204:9986:4539"
    };

    for (uint32_t i = 0; i < data[0].cells_data.size(); ++i) {
        ASSERT_EQ(data[0].cells_data[i].read_id, expected_read_ids[i]);
        ASSERT_EQ(data[0].cells_data[i].base, CharToInt[expected_bases[i]]);
        ASSERT_EQ(data[0].cells_data[i].cell_id, expected_cell_ids[i]);
    }
}

TEST(Reader, three_rows) {
    std::vector<PosData> data;
    std::unordered_set<uint32_t> cell_ids;
    std::tie(data, cell_ids) = read_pileup("data/three_rows.pileup");
    std::vector<uint16_t> all_cell_ids = { 1, 2, 3, 4, 9 };
    ASSERT_THAT(cell_ids, UnorderedElementsAreArray(all_cell_ids));

    ASSERT_EQ(3, data.size());

    std::vector<uint64_t> expected_positions = { 1, 2, 3 };

    for (uint32_t i : { 0, 1, 2 }) {
        ASSERT_EQ(expected_positions[i], data[i].position);


        std::vector<std::vector<uint16_t>> expected_cell_ids
                = { { 1, 3 }, { 1, 3, 2 }, { 1, 2, 3, 9, 4 } };
        std::vector<std::vector<uint8_t>> expected_bases
                = { { 'C', 'A' }, { 'A', 'C', 'T' }, { 'A', 'A', 'A', 'C', 'A' } };
        std::vector<std::vector<std::string>> expected_read_ids
                = { { "R1", "R2" }, { "R1", "R2", "R3" }, { "R1", "R3", "R2", "R4", "R5" } };

        for (uint32_t j = 0; j < data[i].cells_data.size(); ++j) {
            ASSERT_EQ(data[i].cells_data[j].read_id, expected_read_ids[i][j]);
            ASSERT_EQ(data[i].cells_data[j].base, CharToInt[expected_bases[i][j]]);
            ASSERT_EQ(data[i].cells_data[j].cell_id, expected_cell_ids[i][j]);
        }
    }
}

} // namespace
