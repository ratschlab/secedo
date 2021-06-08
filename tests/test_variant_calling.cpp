#include "variant_calling.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <fstream>
#include <vector>

namespace {

struct Vcf {
    uint8_t chromosome;
    uint32_t pos;
    char ref;
    char alt;
    std::string gt;
};

std::vector<Vcf> read_vcf_file(const std::string &path) {
    assert(std::filesystem::exists(path));
    std::ifstream f(path);
    std::string line;
    std::vector<Vcf> result;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::vector<std::string> cols = split(line, '\t');
        result.push_back({ (uint8_t)std::stoi(cols[0]), (uint32_t)std::stoul(cols[1]), cols[3][0],
                           cols[4][0], cols[9] });
    }
    return result;
}

class VariantCalling : public ::testing::Test {
  protected:
    virtual ~VariantCalling() {}

    virtual void SetUp() { std::filesystem::remove_all("data/" + name()); }
    virtual void TearDown() { std::filesystem::remove_all("data/" + name()); }

    std::string name() { return ::testing::UnitTest::GetInstance()->current_test_info()->name(); }
};


TEST(ReadFasta, Empty) {
    std::ifstream f("data/empty.pileup");
    std::vector<uint8_t> chr_data;
    get_next_chromosome(f, {}, &chr_data);
    ASSERT_EQ(0, chr_data.size());
}

TEST(ReadFasta, FemalGenome) {
    std::ifstream f("data/genome_diploid_female.fa");
    std::vector<uint8_t> chr_data;
    std::string expected_maternal = "AAAAAGGGGG";
    std::string expected_paternal = "CCCCCTTTTT";
    get_next_chromosome(f, {}, &chr_data);
    ASSERT_EQ(10, chr_data.size());
    for (uint32_t i = 0; i < chr_data.size(); ++i) {
        ASSERT_EQ(CharToInt[(int)expected_paternal[i]], chr_data[i] & 7);
        ASSERT_EQ(CharToInt[(int)expected_maternal[i]], chr_data[i] >> 3);
    }

    get_next_chromosome(f, {}, &chr_data);
    ASSERT_EQ(8, chr_data.size());

    expected_paternal = "AAAAGGGG";
    expected_maternal = "CCCCNNNN";
    for (uint32_t i = 0; i < chr_data.size(); ++i) {
        ASSERT_EQ(CharToInt[(int)expected_paternal[i]], chr_data[i] & 7);
        ASSERT_EQ(CharToInt[(int)expected_maternal[i]], chr_data[i] >> 3);
    }
}

TEST(ReadFasta, MaleGenome) {
    std::ifstream f("data/genome_diploid_male.fa");
    std::vector<uint8_t> chr_data;
    std::string expected_maternal = "AAAAAGGGGG";
    std::string expected_paternal = "CCCCCTTTTT";
    get_next_chromosome(f, {}, &chr_data);
    ASSERT_EQ(10, chr_data.size());
    for (uint32_t i = 0; i < chr_data.size(); ++i) {
        ASSERT_EQ(CharToInt[(int)expected_paternal[i]], chr_data[i] & 7);
        ASSERT_EQ(CharToInt[(int)expected_maternal[i]], chr_data[i] >> 3);
    }

    get_next_chromosome(f, {}, &chr_data);
    ASSERT_EQ(8, chr_data.size());

    std::string expected_X = "CCCCNNNN";
    for (uint32_t i = 0; i < chr_data.size(); ++i) {
        ASSERT_EQ(CharToInt[(int)expected_X[i]], chr_data[i] & 7);
        ASSERT_EQ(CharToInt[(int)expected_X[i]], chr_data[i] >> 3);
    }

    get_next_chromosome(f, {}, &chr_data);
    ASSERT_EQ(8, chr_data.size());

    std::string expected_Y = "AAAAGGGG";
    for (uint32_t i = 0; i < chr_data.size(); ++i) {
        ASSERT_EQ(CharToInt[(int)expected_Y[i]], chr_data[i] & 7);
        ASSERT_EQ(CharToInt[(int)expected_Y[i]], chr_data[i] >> 3);
    }
}

TEST(ReadFasta, MaleGenomeMapped) {
    std::ifstream f("data/genome_diploid_male.fa");
    std::vector<uint8_t> chr_data;
    std::string expected_maternal = "ANAAAAGGGGG";
    std::string expected_paternal = "NCCCCCTTTTT";
    auto map = read_map("data/genome_diploid_male.map");
    get_next_chromosome(f, map, &chr_data);
    ASSERT_EQ(11, chr_data.size());
    for (uint32_t i = 0; i < chr_data.size(); ++i) {
        ASSERT_EQ(CharToInt[(int)expected_paternal[i]], chr_data[i] & 7);
        ASSERT_EQ(CharToInt[(int)expected_maternal[i]], chr_data[i] >> 3);
    }

    get_next_chromosome(f, map, &chr_data);
    ASSERT_EQ(7, chr_data.size());

    std::string expected_X = "CCCNNNN";
    for (uint32_t i = 0; i < chr_data.size(); ++i) {
        ASSERT_EQ(CharToInt[(int)expected_X[i]], chr_data[i] & 7);
        ASSERT_EQ(CharToInt[(int)expected_X[i]], chr_data[i] >> 3);
    }

    get_next_chromosome(f, map, &chr_data);
    ASSERT_EQ(6, chr_data.size());

    std::string expected_Y = "AAAGGG";
    for (uint32_t i = 0; i < chr_data.size(); ++i) {
        ASSERT_EQ(CharToInt[(int)expected_Y[i]], chr_data[i] & 7);
        ASSERT_EQ(CharToInt[(int)expected_Y[i]], chr_data[i] >> 3);
    }
}

TEST_F(VariantCalling, Empty) {
    variant_calling({ {} }, {}, "data/genome_diploid_female.fa", "", 1e-3, 1e-3, "data/" + name());

    ASSERT_FALSE(std::filesystem::exists("data/" + name() + "/cluster_0.vcf"));
}

TEST_F(VariantCalling, EmptyPos) {
    variant_calling({ {} }, { 1, 1, 1, 2, 2, 2 }, "data/genome_diploid_female.fa", "", 1e-3, 1e-3,
                    "data/" + name());

    ASSERT_TRUE(read_vcf_file("data/" + name() + "/cluster_0.vcf").empty());
}

/**
 * On pos 0 of chromosome 0, the reference genome is heterozygous AC. The new genome will be
 * homozygous AA and the VCF should reflect this.
 */
TEST_F(VariantCalling, OnePosOneVariant) {
    const uint32_t num_cells = 10;
    std::vector<uint32_t> read_ids(num_cells);
    std::vector<uint16_t> cell_ids_bases(num_cells);
    std::iota(read_ids.begin(), read_ids.end(), 0);
    for (uint32_t i = 0; i < num_cells; ++i) {
        cell_ids_bases[i] = i << 2; // all cells have an 'A' in this position
    }
    PosData pd(1, read_ids, cell_ids_bases);
    const std::vector<std::vector<PosData>> pds = { { pd } };
    std::vector<uint16_t> clusters(num_cells, 1);

    variant_calling(pds, clusters, "data/genome_diploid_female.fa", "", 1e-3, 1e-3,
                    "data/" + name());

    std::vector<Vcf> vcfs = read_vcf_file("data/" + name() + "/cluster_0.vcf");

    ASSERT_EQ(1, vcfs.size());
    ASSERT_EQ(vcfs[0].ref, 'C');
    ASSERT_EQ(vcfs[0].alt, 'A');
    ASSERT_EQ(vcfs[0].gt, "1/1");
}

/**
 * On pos 0 of chromosome 0, the reference genome is heterozygous AC. The new genome will also be
 * AC, so the VCF should be empty
 */
TEST_F(VariantCalling, OnePosNoVariant) {
    const uint32_t num_cells = 10;
    std::vector<uint32_t> read_ids(num_cells);
    std::vector<uint16_t> cell_ids_bases(num_cells);
    std::iota(read_ids.begin(), read_ids.end(), 0);
    for (uint32_t i = 0; i < num_cells; ++i) {
        cell_ids_bases[i] = i % 2 ? i << 2 : i << 2 | 1; // cells have A or C at this position
    }
    PosData pd(1, read_ids, cell_ids_bases);
    const std::vector<std::vector<PosData>> pds = { { pd } };
    std::vector<uint16_t> clusters(num_cells, 1);

    variant_calling(pds, clusters, "data/genome_diploid_female.fa", "", 1e-3, 1e-3,
                    "data/" + name());

    std::vector<Vcf> vcfs = read_vcf_file("data/" + name() + "/cluster_0.vcf");

    ASSERT_EQ(0, vcfs.size());
}

/**
 * On pos 0 of chromosome 0, the reference genome is heterozygous AC. The new genome will be
 * GT, so the VCF should contain two rows, one for A-G and one for C-T.
 */
TEST_F(VariantCalling, OnePosTwoVariants) {
    const uint32_t num_cells = 10;
    std::vector<uint32_t> read_ids(num_cells);
    std::vector<uint16_t> cell_ids_bases(num_cells);
    std::iota(read_ids.begin(), read_ids.end(), 0);
    for (uint32_t i = 0; i < num_cells; ++i) {
        cell_ids_bases[i] = i % 2 ? i << 2 | 2 : i << 2 | 3; // cells have G or T at this position
    }
    PosData pd(1, read_ids, cell_ids_bases);
    const std::vector<std::vector<PosData>> pds = { { pd } };
    std::vector<uint16_t> clusters(num_cells, 1);

    variant_calling(pds, clusters, "data/genome_diploid_female.fa", "", 1e-3, 1e-3,
                    "data/" + name());

    std::vector<Vcf> vcfs = read_vcf_file("data/" + name() + "/cluster_0.vcf");

    ASSERT_EQ(2, vcfs.size());

    ASSERT_EQ(vcfs[0].ref, 'A');
    ASSERT_EQ(vcfs[0].alt, 'G');
    ASSERT_EQ(vcfs[0].gt, "1/1");

    ASSERT_EQ(vcfs[1].ref, 'C');
    ASSERT_EQ(vcfs[1].alt, 'T');
    ASSERT_EQ(vcfs[1].gt, "1/1");
}

/**
 * On pos 0 of chromosome 0, the reference genome is heterozygous AC. The new genome will be
 * AT, so the VCF should contain one row for C->T.
 * On pos 0 of chromosome 1, the reference genome is heterozygous CA. The new genome will be
 * AA, so the VCF should contain one row for C->A.
 */
TEST_F(VariantCalling, TwoPosTwoVariants) {
    const uint32_t num_cells = 10;
    std::vector<uint32_t> read_ids(num_cells);
    std::vector<uint16_t> cell_ids_bases1(num_cells);
    std::vector<uint16_t> cell_ids_bases2(num_cells);
    std::iota(read_ids.begin(), read_ids.end(), 0);
    for (uint32_t i = 0; i < num_cells; ++i) {
        cell_ids_bases1[i] = i % 2 ? i << 2 : i << 2 | 3; // cells have AT genotype
    }
    for (uint32_t i = 0; i < num_cells; ++i) {
        cell_ids_bases2[i] = i << 2; // cells have A at this position
    }
    PosData pd1(1, read_ids, cell_ids_bases1);
    PosData pd2(1, read_ids, cell_ids_bases2);
    const std::vector<std::vector<PosData>> pds = { { pd1 }, { pd2 } };
    std::vector<uint16_t> clusters(num_cells, 1);

    variant_calling(pds, clusters, "data/genome_diploid_female.fa", "", 1e-3, 1e-3,
                    "data/" + name());

    std::vector<Vcf> vcfs = read_vcf_file("data/" + name() + "/cluster_0.vcf");

    ASSERT_EQ(2, vcfs.size());

    ASSERT_EQ(vcfs[0].ref, 'C');
    ASSERT_EQ(vcfs[0].alt, 'T');
    ASSERT_EQ(vcfs[0].gt, "1/1");

    ASSERT_EQ(vcfs[1].ref, 'C');
    ASSERT_EQ(vcfs[1].alt, 'A');
    ASSERT_EQ(vcfs[1].gt, "1/1");
}

TEST(ReadMap, EmptyName) {
    ASSERT_TRUE(read_map("").empty());
}

TEST(ReadMap, EmptyFile) {
    std::ofstream f("data/empty.map");
    f.close();
    ASSERT_TRUE(read_map("data/empty.map").empty());
    std::filesystem::remove("data/empty.map");
}

TEST(ReadMap, File) {
    std::unordered_map<std::string, std::vector<ChrMap>> map = read_map("data/test.map");
    ASSERT_EQ(2, map.size());
    ASSERT_EQ(6, map["1_maternal"].size());
    std::vector<uint32_t> expected_pat_pos = { 63734, 66269, 82130, 91547, 285647, 289210, 708137 };
    std::vector<char> expected_pat_tr = { 'D', 'I', 'D', 'I', 'D', 'I', 'I' };
    std::vector<uint32_t> expected_pat_len = { 3, 1, 2, 1, 4, 1, 1 };
    ASSERT_EQ(expected_pat_pos.size(), map["1_paternal"].size());
    for (uint32_t i = 0; i < expected_pat_pos.size(); ++i) {
        ASSERT_EQ(expected_pat_pos[i], map["1_paternal"][i].start_pos);
        ASSERT_EQ(expected_pat_tr[i], map["1_paternal"][i].tr);
        ASSERT_EQ(expected_pat_len[i], map["1_paternal"][i].len);
    }
}

TEST(ApplyMap, EmptyMap) {
    std::vector<uint8_t> chromosome = { 0, 0, 0, 0 }; // AAAA
    std::vector<uint8_t> result;
    apply_map({}, chromosome, &result);
    ASSERT_EQ(chromosome, result);
}

TEST(ApplyMap, OneDeletionBeg) {
    std::vector<uint8_t> chromosome = { 0, 0, 0, 0 }; // AAAA
    std::vector<uint8_t> expected = { 5, 0, 0, 0, 0 }; // NAAAA
    std::vector<uint8_t> result;
    ChrMap map_el = { 1, 0, 1, 'D' };
    apply_map({ map_el }, chromosome, &result);
    ASSERT_EQ(expected, result);
}

TEST(ApplyMap, OneDeletionMid) {
    std::vector<uint8_t> chromosome = { 0, 0, 0, 0 }; // AAAA
    std::vector<uint8_t> expected = { 0, 0, 5, 0, 0 }; // AANAA
    std::vector<uint8_t> result;
    ChrMap map_el = { 1, 2, 1, 'D' };
    apply_map({ map_el }, chromosome, &result);
    ASSERT_EQ(expected, result);
}

TEST(ApplyMap, OneInsertionBeg) {
    std::vector<uint8_t> chromosome = { 0, 1, 2, 3 }; // ACGT
    std::vector<uint8_t> expected = { 1, 2, 3 }; // CGT
    std::vector<uint8_t> result;
    ChrMap map_el = { 1, 0, 1, 'I' };
    apply_map({ map_el }, chromosome, &result);
    ASSERT_EQ(expected, result);
}

TEST(ApplyMap, OneInsertionMid) {
    std::vector<uint8_t> chromosome = { 0, 1, 2, 3 }; // ACGT
    std::vector<uint8_t> expected = { 0, 1, 3 }; // ACT
    std::vector<uint8_t> result;
    ChrMap map_el = { 1, 2, 1, 'I' };
    apply_map({ map_el }, chromosome, &result);
    ASSERT_EQ(expected, result);
}

TEST(ApplyMap, InsertionBegDeletionMid) {
    std::vector<uint8_t> chromosome = { 0, 1, 2, 3 };
    std::vector<uint8_t> expected = { 5, 0, 1, 3 };
    std::vector<uint8_t> result;
    ChrMap map_el1 = { 1, 0, 1, 'D' };
    ChrMap map_el2 = { 1, 2, 1, 'I' };
    apply_map({ map_el1, map_el2 }, chromosome, &result);
    ASSERT_EQ(expected, result);
}

TEST(ApplyMap, InsertionMidDeletionBeg) {
    std::vector<uint8_t> chromosome = { 0, 1, 2, 3 };
    std::vector<uint8_t> expected = { 1, 5, 2, 3 };
    std::vector<uint8_t> result;
    ChrMap map_el1 = { 1, 0, 1, 'I' };
    ChrMap map_el2 = { 1, 2, 1, 'D' };
    apply_map({ map_el1, map_el2 }, chromosome, &result);
    ASSERT_EQ(expected, result);
}


} // namespace
