#include "pileup_reader.hpp"
#include "sequenced_data.hpp"
#include "similarity_matrix.hpp"
#include "spectral_clustering.hpp"
#include "util.hpp"
#include "variant_calling.hpp"

#include <gflags/gflags.h>
#include <spdlog/spdlog.h>

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

DEFINE_double(seq_error_rate, 0.001, "Sequencing errors rate, denoted by theta");
DEFINE_double(mutation_rate,
              0,
              "epsilon, estimated frequency of mutated loci in the pre-processed data set");
// estimate of how many positions are actually homozygous germline, were only included because of
// sequencing (or alignment!) errors
DEFINE_double(
        hzygous_prob,
        0,
        "The probability that a loci is homozygous, (not filtered correctly in the first step");

DEFINE_string(mpileup_file,
              "",
              "Input file containing 'pileup' textual format from an alignment, as written by "
              "preprocessing.py");
DEFINE_string(cells_file,
              "",
              "File with identifiers of cells (numbers between 0 and (num_of_cells - 1); the "
              "matrices will be computed only for these cells; if absent, all cells are used.");
DEFINE_string(cell_group_file,
              "",
              "For each cell, it contains a group to which it belongs. Cells from one group will "
              "be treated as one cell. Useful e.g. when creating data with artificially higher "
              "coverage; if absent, each cell has its own group");

DEFINE_string(out_dir, "./", "Directory where the similarity matrices will be written to");
// TODO: do we need this or leave as constexpr?
DEFINE_uint32(max_read_size,
              1000,
              "Maximal considered insert size (for paired-end sequencing). In plain English, this "
              "is the maximal length of a read we consider. Reads that are longer will not be "
              "mapped correctly.");

DEFINE_uint32(num_threads, 8, "Number of threads to use");

DEFINE_string(labels_file, "", "Input file containing labels");

DEFINE_string(chromosomes,
              "The chromosomes on which to run the algorithm",
              "1,2,3,4,5,6,7,8,9,10,11,12,13.14,15,16,17,18,19,20,21,22,X");

DEFINE_string(log_level,
              "trace",
              "The log verbosity: debug, trace, info, warn, error, critical, off");

int main(int argc, char *argv[]) {
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    spdlog::set_level(spdlog::level::from_str(FLAGS_log_level));

    std::vector<std::filesystem::path> mpileup_files = { FLAGS_mpileup_file };
    // if the input is a directory, get all pileup files in the directory
    if (std::filesystem::is_directory(FLAGS_mpileup_file)) {
        mpileup_files = get_files(FLAGS_mpileup_file, ".pileup");
        spdlog::info("Found {} .pileup files in '{}'", mpileup_files.size(), FLAGS_mpileup_file);
    }

    spdlog::trace("Reading data...");
    std::vector<std::vector<PosData>> pos_data(mpileup_files.size());
    std::vector<std::unordered_set<uint32_t>> cell_ids(mpileup_files.size());
#pragma omp parallel for num_threads(FLAGS_num_threads)
    for (uint32_t i = 0; i < pos_data.size(); ++i) {
        std::tie(pos_data[i], cell_ids[i]) = read_pileup(mpileup_files[i]);
    }
    std::unordered_set<uint32_t> all_cell_ids;
    for (uint32_t i = 0; i < cell_ids.size(); ++i) {
        spdlog::trace("{} has {} cell ids.", mpileup_files[i].string(), cell_ids[i].size());
        std::copy(cell_ids[i].begin(), cell_ids[i].end(),
                  std::inserter(all_cell_ids, all_cell_ids.end()));
    }
    uint32_t num_cells = *std::max_element(all_cell_ids.begin(), all_cell_ids.end()) + 1;

    spdlog::trace("Computing similarity matrix...");
    std::vector<uint32_t> cell_id_map;
    if (!FLAGS_cells_file.empty()) {
        cell_id_map = int_split<uint32_t>(read_file(FLAGS_cells_file), ' ');
        // num_cells = len(cells_ids);
    } else {
        cell_id_map.resize(num_cells);
        std::iota(cell_id_map.begin(), cell_id_map.end(), 0);
    }

    std::vector<uint32_t> cell_group_map;
    if (!FLAGS_cell_group_file.empty()) {
        cell_group_map = int_split<uint32_t>(read_file(FLAGS_cells_file), ',');
    } else {
        cell_group_map.resize(num_cells);
        std::iota(cell_group_map.begin(), cell_group_map.end(), 0);
    }

    computeSimilarityMatrix(pos_data, num_cells, cell_id_map, cell_group_map, FLAGS_mutation_rate,
                            FLAGS_hzygous_prob, FLAGS_seq_error_rate, FLAGS_num_threads,
                            FLAGS_out_dir);
    spdlog::info("Done.");
}
