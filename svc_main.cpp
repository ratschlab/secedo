#include "sequenced_data.hpp"
#include "similarity_matrix.hpp"
#include "spectral_clustering.hpp"
#include "util/is_significant.hpp"
#include "util/logger.hpp"
#include "util/pileup_reader.hpp"
#include "util/util.hpp"
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
              "1,2,3,4,5,6,7,8,9,10,11,12,13.14,15,16,17,18,19,20,21,22,X,Y");

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

static bool ValidateTermination(const char *flagname, const std::string &value) {
    if (value != "AIC" && value != "BIC") {
        printf("Invalid value for --%s: %s.\nShould be one of AIC, BIC\n", flagname, value.c_str());
        return false;
    }
    return true;
}
DEFINE_string(termination,
              "BIC",
              "Which criteria to use for determining if a Simple or Multivariate Gaussian matches "
              "the data better (AIC/BIC)");
DEFINE_validator(termination, ValidateTermination);

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

DEFINE_bool(
        arma_kmeans,
        false,
        "Whether to use Armadillo's k-means implementation or our own (which is very simple, but "
        "it weights the Fiedler vector higher, reducing the effect of the curse of dimensionality");

DEFINE_string(pos_file,
              "",
              "When present, only consider positions in this file. The file must have 2 columns, "
              "first one is chromosome id, second is position.");

/**
 * Recursively divides cells into 2 sub-clusters until a termination criteria is met.
 * N - number of cells
 * G - number of cell groups (normally N=G, but in case we group cells to artificially increase
 * coverage we have G < N)
 * NC - number of cells in the current sub-cluster. At the first call of divide() NC=G.
 * @param pos_data
 * @param max_read_length length of the longest fragment (typically around 500)
 * @param id_to_group of size N maps cell ids to cell groups. Data from cells in the same group is
 * treated as if it came from one cell. Used to artificially increase coverage when testing
 * @param id_to_pos of size G maps a cell id to its position in the similarity matrix as we
 * subdivide into smaller and smaller clusters. At the beginning, this is the identity permutation.
 * If a cell with id 'cell_id' is not in the current cluster, then id_to_pos[cell_id]==NO_POS. The
 * position of a cell in the similarity matrix is given by id_to_pos[id_to_group[cell_id]].
 * @param pos_to_id of size NC the inverse of #id_to_pos, it maps each position 0...pos in the
 * current subgroup to the actual cell it corresponds to
 * @param mutation_rate epsilon, estimated frequency of mutated loci in the pre-processed data set
 * @param homozygous_rate  the probability that a locus is homozygous, (not filtered correctly in
 * the first step)
 * @param seq_error_rate error rate of the sequencer, e.g. 1e-3 if using Illumina reads with base
 * quality >30
 * @param num_threads number of threads to use for the computation
 * @param out_dir where to output the clustering results
 * @param normalization the type of normalization to use for the similiarity matrix (see the flag
 * with the same name)
 * @param marker marks the current sub-cluster; for example AB means we are in the second
 * sub-cluster (B) of the first cluster (A)
 */
void variant_call(const std::vector<std::vector<PosData>> &pds,
                  uint32_t max_read_length,
                  const std::vector<uint16_t> &id_to_group,
                  const std::vector<uint32_t> &id_to_pos,
                  const std::vector<uint32_t> &pos_to_id,
                  double mutation_rate,
                  double homozygous_rate,
                  double seq_error_rate,
                  const uint32_t num_threads,
                  const std::string &out_dir,
                  const std::string &normalization,
                  const std::string marker) {
    if (!marker.empty()) {
        logger()->info("\n\nPerforming clustering of sub-cluster {} with {} elements", marker,
                       pos_to_id.size());
    }
    logger()->info("Filtering significant positions...");
    Filter filter;
    auto [pos_data, coverage]
            = filter.filter(pds, id_to_group, id_to_pos, marker, seq_error_rate, num_threads);
    if (coverage < 9) {
        logger()->trace("Coverage of cluster {} is lower than 9. Stopping.", marker);
        return;
    }

    logger()->info("Computing similarity matrix...");
    uint32_t n_cells_subcluster = pos_to_id.size();
    uint32_t n_cells_total = id_to_group.size();
    Matd sim_mat = computeSimilarityMatrix(pos_data, n_cells_subcluster, max_read_length, id_to_pos,
                                           mutation_rate, homozygous_rate, seq_error_rate,
                                           num_threads, FLAGS_o, marker, FLAGS_normalization);

    logger()->info("Performing spectral clustering...");
    std::vector<double> cluster; // size n_cells
    Termination termination = parse_termination(FLAGS_termination);
    ClusteringType clustering_type = parse_clustering_type(FLAGS_clustering_type);
    bool is_done = spectral_clustering(sim_mat, clustering_type, termination, FLAGS_o, marker,
                                       FLAGS_arma_kmeans, &cluster);
    if (is_done) {
        return;
    }

    std::vector<uint16_t> id_to_cluster(n_cells_total);
    for (uint16_t cell_id = 0; cell_id < n_cells_total; ++cell_id) {
        uint16_t pos = id_to_pos[id_to_group[cell_id]];
        id_to_cluster[cell_id] = pos == NO_POS ? NO_POS : cluster[pos];
    }
    write_vec(std::filesystem::path(out_dir) / ("spectral_clustering" + marker), id_to_cluster);

    logger()->info("Performing clustering refinement via expectation maximization...");
    expectation_maximization(pos_data, id_to_pos, FLAGS_num_threads, FLAGS_seq_error_rate,
                             &cluster);

    for (uint16_t i = 0; i < n_cells_total; ++i) {
        uint32_t pos = id_to_pos[id_to_group[i]];
        id_to_cluster[i] = pos == NO_POS ? NO_POS : cluster[pos];
    }
    write_vec(std::filesystem::path(out_dir) / ("expectation_maximization" + marker),
              id_to_cluster);

    std::vector<uint32_t> pos_to_id_a;
    std::vector<uint32_t> pos_to_id_b;
    std::vector<uint32_t> id_to_pos_a(id_to_pos.size(), NO_POS);
    std::vector<uint32_t> id_to_pos_b(id_to_pos.size(), NO_POS);

    for (uint32_t cell_idx = 0; cell_idx < n_cells_subcluster; ++cell_idx) {
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

    variant_call(pds, max_read_length, id_to_group, id_to_pos_a, pos_to_id_a, mutation_rate,
                 homozygous_rate, seq_error_rate, num_threads, out_dir, normalization,
                 marker + 'A');
    variant_call(pds, max_read_length, id_to_group, id_to_pos_b, pos_to_id_b, mutation_rate,
                 homozygous_rate, seq_error_rate, num_threads, out_dir, normalization,
                 marker + 'B');
}

//============================================================================
int main(int argc, char *argv[]) {
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    spdlog::set_level(spdlog::level::from_str(FLAGS_log_level));

    std::vector<uint16_t> id_to_group
            = get_grouping(FLAGS_merge_count, FLAGS_merge_file, FLAGS_max_cell_count);

    std::vector<std::filesystem::path> input_files = { FLAGS_i };
    // if the input is a directory, get all pileup files in the directory
    if (std::filesystem::is_directory(FLAGS_i)) {
        input_files = get_files(FLAGS_i, ".bin");
        if (input_files.empty()) {
            logger()->info("No binary pileup files found. Looking for textual pileup files...");
            input_files = get_files(FLAGS_i, ".pileup");
        }
        std::sort(input_files.begin(), input_files.end());
        logger()->info("Found {} input files in '{}'", input_files.size(), FLAGS_i);
    }

    if (input_files.empty()) {
        logger()->info("No input files found in {}. Nothing to do.", FLAGS_i);
        std::exit(0);
    }

    std::vector<std::vector<uint32_t>> positions;
    if (!FLAGS_pos_file.empty()) {
        positions = read_positions(FLAGS_pos_file);
        if (positions.size() < input_files.size()) {
            //TODO: this won't work for the Y chromosome
            logger()->error(
                    "Number of chromosomes in {} ({}) does not match number of input files ({})",
                    FLAGS_pos_file, positions.size(), input_files.size());
            std::exit(1);
        }
    }

    uint64_t total_size = 0;
    for (const auto &f : input_files) {
        total_size += std::filesystem::file_size(f);
    }

    // read input files in parallel
    std::vector<std::vector<PosData>> pos_data(input_files.size());
    std::vector<std::unordered_set<uint32_t>> cell_ids(input_files.size());
    std::vector<uint32_t> max_read_lengths(input_files.size());

    ProgressBar read_progress(total_size, "Reading progress", std::cout);
    read_progress.SetFrequencyUpdate(total_size / 100);
#pragma omp parallel for num_threads(FLAGS_num_threads)
    for (uint32_t i = 0; i < pos_data.size(); ++i) {
        std::tie(pos_data[i], cell_ids[i], max_read_lengths[i]) = read_pileup(
                input_files[i], id_to_group,
                [&read_progress](uint32_t progress) { read_progress += progress; },
                FLAGS_max_coverage, positions[i]);
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

    // this maps the cell id to its actual position as we divide cells into smaller and smaller
    // clusters
    std::vector<uint32_t> cell_id_map(num_groups);
    std::iota(cell_id_map.begin(), cell_id_map.end(), 0);

    variant_call(pos_data, max_read_length, id_to_group, cell_id_map, cell_id_map,
                 FLAGS_mutation_rate, FLAGS_homozygous_prob, FLAGS_seq_error_rate,
                 FLAGS_num_threads, FLAGS_o, FLAGS_normalization, "");


    logger()->info("Done.");
}
