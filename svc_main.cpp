#include "logger.hpp"
#include "pileup_reader.hpp"
#include "preprocess.hpp"
#include "sequenced_data.hpp"
#include "similarity_matrix.hpp"
#include "spectral_clustering.hpp"
#include "util.hpp"
#include "variant_calling.hpp"

#include <gflags/gflags.h>
#include <progress_bar/progress_bar.hpp>

#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <vector>

DEFINE_double(seq_error_rate, 0.01, "Sequencing errors rate, denoted by theta");
DEFINE_double(mutation_rate,
              0.01,
              "epsilon, estimated frequency of mutated loci in the pre-processed data set");
// estimate of how many positions are actually homozygous germline, were only included because of
// sequencing (or alignment!) errors
DEFINE_double(
        homozygous_prob,
        0.15,
        "The probability that a locus is homozygous, (not filtered correctly in the first step");

DEFINE_string(i,
              "",
              "Input file or directory containing 'pileup' textual or binary format from"
              " an alignment, as written by preprocessing.py");

DEFINE_string(o, "./", "Directory where the output will be written.");

DEFINE_uint32(num_threads, 8, "Number of threads to use");

DEFINE_string(labels_file, "", "Input file containing labels");

DEFINE_string(chromosomes,
              "The chromosomes on which to run the algorithm",
              "1,2,3,4,5,6,7,8,9,10,11,12,13.14,15,16,17,18,19,20,21,22,X");

DEFINE_string(log_level,
              "trace",
              "The log verbosity: debug, trace, info, warn, error, critical, off");

DEFINE_uint32(merge_count,
              1,
              "Pool data from  merge_count consecutive cells as if they were a single cell, in"
              "  order to artificially increase coverage. Only work on synthetic data where we know"
              " consecutive cells are part of the same cluster!");

DEFINE_string(merge_file,
              "",
              "File containing cell grouping. Cells in the same group are treated as if they were "
              "a single cell. Useful for artificially increasing coverage for testing.");

DEFINE_string(clustering_type,
              "SPECTRAL6",
              "How to perform spectral clustering. One of FIEDLER, SPECTRAL2, SPECTRAL6, "
              "GMM_ASSIGN, GMM_PROB. See spectral_clustering.hpp for details.");

static bool ValidateClusteringType(const char *flagname, const std::string &value) {
    if (value != "FIEDLER" && value != "SPECTRAL2" && value != "SPECTRAL6" && value != "GMM_PROB"
        && value != "GMM_ASSIGN") {
        printf("Invalid value for --%s: %s.\nShould be one of FIEDLER, SPECTRAL2, SPECTRAL6, "
               "GMM_ASSIGN, GMM_PROB\n",
               flagname, value.c_str());
        return false;
    }
    return true;
}

DEFINE_validator(clustering_type, ValidateClusteringType);

static bool ValidateNormalization(const char *flagname, const std::string &value) {
    if (value != "ADD_MIN" && value != "EXPONENTIATE" && value != "SCALE_MAX_1") {
        printf("Invalid value for --%s: %s.\nShould be one of ADD_MIN, EXPONENTIATE, SCALE_MAX_1\n",
               flagname, value.c_str());
        return false;
    }
    return true;
}
DEFINE_string(normalization,
              "ADD_MIN",
              "How to normalize the similarity matrix. One of ADD_MIN, EXPONENTIATE, SCALE_MAX_1");
DEFINE_validator(normalization, &ValidateNormalization);

DEFINE_uint32(max_cell_count,
              10'000,
              "Maximum expected cell count, used to initialize the cell grouping");

DEFINE_uint32(max_coverage,
              100,
              "Positions with higher coverage are considered anomalies and discarded");

constexpr uint16_t NO_POS = std::numeric_limits<uint16_t>::max();

void divide(const std::vector<std::vector<PosData>> &pos_data,
            uint32_t max_read_length,
            const std::vector<uint16_t> &id_to_group,
            const std::vector<uint32_t> &id_to_pos,
            const std::vector<uint32_t> &pos_to_id,
            double mutation_rate,
            double heterozygous_rate,
            double seq_error_rate,
            const uint32_t num_threads,
            const std::string &out_dir,
            const std::string &normalization,
            const std::string marker) {
    std::cout << std::endl << std::endl;
    if (!marker.empty()) {
        logger()->info("Performing clustering of sub-cluster {} with {} elements", marker,
                       pos_to_id.size());
    }

    logger()->info("Computing similarity matrix...");
    Matd sim_mat = computeSimilarityMatrix(pos_data, pos_to_id.size(), max_read_length, id_to_pos,
                                           mutation_rate, heterozygous_rate, seq_error_rate,
                                           num_threads, FLAGS_o, marker, FLAGS_normalization);

    logger()->info("Performing spectral clustering...");
    std::vector<double> cluster;
    bool is_done = spectral_clustering(sim_mat, FLAGS_clustering_type, Termination::AIC, FLAGS_o,
                                       marker, &cluster);
    if (is_done) {
        return;
    }

    std::vector<uint16_t> id_to_cluster(id_to_group.size());
    for (uint16_t i = 0; i < id_to_group.size(); ++i) {
        uint16_t pos = id_to_pos[id_to_group[i]];
        id_to_cluster[i] = pos == NO_POS ? NO_POS : cluster[pos];
    }
    write_vec(std::filesystem::path(out_dir) / ("spectral_clustering" + marker), id_to_cluster);

    logger()->info("Performing clustering refinement via expectation maximization...");
    expectation_maximization(pos_data, id_to_pos, FLAGS_num_threads, FLAGS_seq_error_rate,
                             &cluster);

    for (uint16_t i = 0; i < id_to_group.size(); ++i) {
        uint32_t pos = id_to_pos[id_to_group[i]];
        id_to_cluster[i] = pos == NO_POS ? NO_POS : cluster[pos];
    }
    write_vec(std::filesystem::path(out_dir) / ("expectation_maximization" + marker),
              id_to_cluster);

    std::vector<std::vector<PosData>> pos_data_a;
    std::vector<std::vector<PosData>> pos_data_b;
    std::vector<uint32_t> pos_to_id_a;
    std::vector<uint32_t> pos_to_id_b;
    std::vector<uint32_t> id_to_pos_a(id_to_pos.size(), NO_POS);
    std::vector<uint32_t> id_to_pos_b(id_to_pos.size(), NO_POS);

    for (uint32_t cell_idx = 0; cell_idx < cluster.size(); ++cell_idx) {
        uint32_t cell_id = pos_to_id[cell_idx];
        if (cluster[cell_idx] < 0.05) {
            id_to_pos_a[cell_id] = pos_to_id_a.size();
            pos_to_id_a.push_back(cell_id);
        } else if (cluster[cell_idx] > 0.95) {
            id_to_pos_b[cell_id] = pos_to_id_b.size();
            pos_to_id_b.push_back(cell_id);
        }
    }

    if (pos_to_id_a.size() < 30 || pos_to_id_b.size() < 30) {
        logger()->trace("Cluster size is too small. Stopping.");
        return; // TODO: hack - figure it out
    }

    uint32_t total_coverage_a = 0;
    uint32_t total_coverage_b = 0;
    uint32_t total_positions_a = 0;
    uint32_t total_positions_b = 0;

    for (uint32_t chr_idx = 0; chr_idx < pos_data.size(); ++chr_idx) {
        std::vector<PosData> positions_a;
        std::vector<PosData> positions_b;
        for (uint32_t pos_idx = 0; pos_idx < pos_data[chr_idx].size(); ++pos_idx) {
            std::vector<CellData> cd_a;
            std::vector<CellData> cd_b;
            for (uint32_t cell_idx = 0; cell_idx < pos_data[chr_idx][pos_idx].cells_data.size();
                 ++cell_idx) {
                uint16_t cell_id = pos_data[chr_idx][pos_idx].cells_data[cell_idx].cell_id;
                if (cluster[id_to_pos[cell_id]] < 0.05) {
                    cd_a.push_back(pos_data[chr_idx][pos_idx].cells_data[cell_idx]);
                } else if (cluster[id_to_pos[cell_id]] > 0.95) {
                    cd_b.push_back(pos_data[chr_idx][pos_idx].cells_data[cell_idx]);
                }
            }
            uint32_t coverage;
            PosData pda = { pos_data[chr_idx][pos_idx].position, cd_a };
            if (is_significant(pda, seq_error_rate, &coverage)) {
                positions_a.push_back(pda);
                total_coverage_a += coverage;
            }
            PosData pdb = { pos_data[chr_idx][pos_idx].position, cd_b };
            if (is_significant(pdb, seq_error_rate, &coverage)) {
                positions_b.push_back(pdb);
                total_coverage_b += coverage;
            }
        }
        pos_data_a.push_back(positions_a);
        pos_data_b.push_back(positions_b);
        total_positions_a += positions_a.size();
        total_positions_b += positions_b.size();
    }

    double coverage_a = total_positions_a == 0
            ? 0
            : static_cast<double>(total_coverage_a) / total_positions_a;
    double coverage_b = total_positions_b == 0
            ? 0
            : static_cast<double>(total_coverage_b) / total_positions_b;
    logger()->trace("Avg coverage for cluster {}: {}. Total positions: {}", marker + 'A',
                    coverage_a, total_positions_a);
    logger()->trace("Avg coverage for cluster {}: {}. Total positions: {}", marker + 'B',
                    coverage_b, total_positions_b);
    // recursively try to further divide each of the new clusters
    if (coverage_a < 9 || coverage_b < 9) {
        logger()->trace("Coverage for at least one cluster is lower than 9. Stopping.");
        return;
    }

    divide(pos_data_a, max_read_length, id_to_group, id_to_pos_a, pos_to_id_a, mutation_rate,
           heterozygous_rate, seq_error_rate, num_threads, out_dir, normalization, marker + 'A');
    divide(pos_data_b, max_read_length, id_to_group, id_to_pos_b, pos_to_id_b, mutation_rate,
           heterozygous_rate, seq_error_rate, num_threads, out_dir, normalization, marker + 'B');
}

int main(int argc, char *argv[]) {
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    spdlog::set_level(spdlog::level::from_str(FLAGS_log_level));

    std::vector<std::filesystem::path> mpileup_files = { FLAGS_i };
    // if the input is a directory, get all pileup files in the directory
    if (std::filesystem::is_directory(FLAGS_i)) {
        mpileup_files = get_files(FLAGS_i, ".bin");
        if (mpileup_files.empty()) {
            logger()->info("No binary pileup files found. Looking for textual pileup files...");
            mpileup_files = get_files(FLAGS_i, ".pileup");
        }
        spdlog::info("Found {} .pileup files in '{}'", mpileup_files.size(), FLAGS_i);
    }

    if (mpileup_files.empty()) {
        spdlog::info("No pileup files found in {}. Bailing out.", FLAGS_i);
        std::exit(0);
    }

    std::vector<uint16_t> id_to_group
            = get_grouping(FLAGS_merge_count, FLAGS_merge_file, FLAGS_max_cell_count);

    uint64_t total_size = 0;
    for (const auto &f : mpileup_files) {
        total_size += std::filesystem::file_size(f);
    }
    ProgressBar read_progress(total_size, "Reading progress", std::cout);

    // read input files in parallel
    std::vector<std::vector<PosData>> pos_data(mpileup_files.size());
    std::vector<std::unordered_set<uint32_t>> cell_ids(mpileup_files.size());
    std::vector<uint32_t> max_read_lengths(mpileup_files.size());
#pragma omp parallel for num_threads(FLAGS_num_threads)
    for (uint32_t i = 0; i < pos_data.size(); ++i) {
        std::tie(pos_data[i], cell_ids[i], max_read_lengths[i])
                = read_pileup(mpileup_files[i], id_to_group,
                              [&read_progress](uint32_t progress) { read_progress += progress; }, FLAGS_max_coverage);
    }
    uint32_t max_read_length = *std::max_element(max_read_lengths.begin(), max_read_lengths.end());

    std::unordered_set<uint32_t> all_cell_ids;
    for (uint32_t i = 0; i < cell_ids.size(); ++i) {
        std::copy(cell_ids[i].begin(), cell_ids[i].end(),
                  std::inserter(all_cell_ids, all_cell_ids.end()));
    }
    uint32_t num_cells = *std::max_element(all_cell_ids.begin(), all_cell_ids.end()) + 1;

    if (FLAGS_merge_file.empty()) {
        // now that we know the actual number of cells, resize the mapping
        id_to_group.resize(num_cells);
    } else if (num_cells != id_to_group.size()) {
        logger()->error(
                "Invalid merge file {}. Merge files contains {} cell ids, data has {} cell ids",
                FLAGS_merge_file, id_to_group.size(), num_cells);
    }

    uint32_t num_groups = *std::max_element(id_to_group.begin(), id_to_group.end()) + 1;

    // this will maps the cell id to its actual position and vice-versa, as we divide cells into
    // smaller and smaller clusters
    std::vector<uint32_t> cell_id_map(num_groups);
    std::iota(cell_id_map.begin(), cell_id_map.end(), 0);

    divide(pos_data, max_read_length, id_to_group, cell_id_map, cell_id_map, FLAGS_mutation_rate,
           FLAGS_homozygous_prob, FLAGS_seq_error_rate, FLAGS_num_threads, FLAGS_o,
           FLAGS_normalization, "");


    logger()->info("Done.");
}
