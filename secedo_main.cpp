#include "sequenced_data.hpp"
#include "spectral_clustering.hpp"
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
// The filtering step attempts to eliminate positions not consistent with the "all cells have the
// same  (heterozygous or homozygous) genotype". This value approximates the percentage of positions
// that snuck through although they are homozygous and have the same genotype (they were mis-labeled
// as 'interesting' due to sequencing errors)
DEFINE_double(homozygous_filtered_rate,
              0.5,
              "The faction of homozygous non-informative loci, i.e. the probability that cells "
              "at a filtered locus actually have identical homozygous genotype.");

// While the value above refers to the filtered loci, this value refers to all loci and it simply
// estimates the rate of heterozygous loci in the human genome
DEFINE_double(heterozygous_prob,
              1e-3,
              "The probability that a locus is heterozygous in the human genome");

DEFINE_string(i,
              "./",
              "Input file or directory containing 'pileup' textual or binary format from"
              " an alignment, as written by preprocessing.py");

DEFINE_string(o, "./", "Directory where the output will be written.");

DEFINE_uint32(num_threads, 8, "Number of threads to use");

DEFINE_string(labels_file, "", "Input file containing labels");

DEFINE_string(chromosomes,
              "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X",
              "The chromosomes on which to run the algorithm");

DEFINE_string(log_level,
              "trace",
              "The log verbosity: debug, trace, info, warn, error, critical, off");

DEFINE_uint32(merge_count,
              1,
              "Pool data from  merge_count consecutive cells as if they were a single cell, in"
              "  order to artificially increase coverage. Only works on synthetic data where we "
              "know consecutive cells are part of the same cluster!");

DEFINE_string(merge_file,
              "",
              "File containing cell grouping. Cells in the same group are treated as if they were "
              "a single cell. Useful for artificially increasing coverage for testing.");

DEFINE_string(clustering_type,
              "SPECTRAL6",
              "How to perform spectral clustering. One of FIEDLER, SPECTRAL2, SPECTRAL6. "
              "See spectral_clustering.hpp for details.");

DEFINE_bool(expectation_maximization,
            false,
            "If true, the spectral clustering results will be refined using an Expectation "
            "Maximization algorithm");

DEFINE_string(
        reference_genome,
        "",
        "The genome against which variants are called, if provided. The fasta file must have 2x24 "
        "sections, with the maternal chromosome being followed by the paternal chromosome");
DEFINE_string(map_file,
              "",
              "If not empty, maps positions in --reference_genome to positions in the haploid "
              "genome that --reference_genome is based on (e.g. to GRCh38)");

static bool ValidateClusteringType(const char *flagname, const std::string &value) {
    if (value != "FIEDLER" && value != "SPECTRAL2" && value != "SPECTRAL6") {
        printf("Invalid value for --%s: %s.\nShould be one of FIEDLER, SPECTRAL2, SPECTRAL6\n",
               flagname, value.c_str());
        return false;
    }
    return true;
}
DEFINE_validator(clustering_type, ValidateClusteringType);

DEFINE_uint32(min_cluster_size,
              100,
              "Stop clustering when the size of a cluster is below this value");

DEFINE_bool(compute_read_stats,
            false,
            "If set to true, the algorithm will compute the maximum fragment length and the "
            "longest fragment for each pileup file (expensive)");

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

DEFINE_string(clustering,
              "",
              "If provided read the clustering from this file and only perform variant calling");

//============================================================================
int main(int argc, char *argv[]) {
    if (argc == 1) {
        gflags::ShowUsageWithFlags(argv[0]);
        exit(1);
    }
    gflags::SetUsageMessage("secedo -i <input_dir> -o <output_dir> arguments");
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    spdlog::set_level(spdlog::level::from_str(FLAGS_log_level));

    std::vector<uint16_t> id_to_group
            = get_grouping(FLAGS_merge_count, FLAGS_merge_file, FLAGS_max_cell_count);

    std::vector<std::string> chromosomes_str = split(FLAGS_chromosomes, ',');
    std::vector<uint32_t> chromosome_ids;
    for (const auto &chr : chromosomes_str) {
        chromosome_ids.push_back(chromosome_to_id(chr));
    }

    std::vector<std::filesystem::path> input_files = { FLAGS_i };
    // if the input is a directory, get all pileup files in the directory
    if (std::filesystem::is_directory(FLAGS_i)) {
        logger()->info("Looking for binary pileup files in {}", std::filesystem::absolute(FLAGS_i));
        input_files = get_files(FLAGS_i, ".bin");
        if (input_files.empty()) {
            logger()->info("No binary pileup files found. Looking for textual pileup files...");
            input_files = get_files(FLAGS_i, ".pileup");
            logger()->info("Found {} textual pileup files files {}", input_files.size(),
                           join_vec(input_files));
        } else {
            logger()->info("Found {} binary pileup files: {}", input_files.size(),
                           join_vec(input_files));
        }
        std::sort(input_files.begin(), input_files.end());
    }

    if (input_files.empty()) {
        logger()->info("No input files found in {}. Nothing to do.", FLAGS_i);
        std::exit(0);
    }

    std::vector<std::vector<uint32_t>> positions;
    if (!FLAGS_pos_file.empty()) {
        positions = read_positions(FLAGS_pos_file);
        if (positions.size() < input_files.size()) {
            // TODO: this won't work for the Y chromosome
            logger()->error(
                    "Number of chromosomes in {} ({}) does not match number of input files ({})",
                    FLAGS_pos_file, positions.size(), input_files.size());
            std::exit(1);
        }
    } else {
        positions.resize(input_files.size());
    }

    uint64_t total_size = 0;
    std::unordered_set<uint32_t> available_chromosome_ids;
    for (const auto &f : input_files) {
        total_size += std::filesystem::file_size(f);
        available_chromosome_ids.insert(get_chromosome(f));
    }

    for (auto chr_id : chromosome_ids) {
        auto it = available_chromosome_ids.find(chr_id);
        if (it == available_chromosome_ids.end()) {
            logger()->error(
                    "Chromosome {} specified with --chromosomes={}, but no input file for it was "
                    "found",
                    id_to_chromosome(chr_id), FLAGS_chromosomes);
            std::exit(1);
        }
    }

    constexpr uint32_t chr_count = 24; // total number of chromosomes (including X and Y)

    // read input files in parallel
    std::vector<std::vector<PosData>> pos_data(chr_count);
    std::vector<uint16_t> num_cells_chr(chr_count);
    std::vector<uint32_t> max_read_lengths(chr_count);

    ProgressBar read_progress(total_size, "Reading progress", std::cout);
    read_progress.SetFrequencyUpdate(total_size / 100);

#pragma omp parallel for num_threads(FLAGS_num_threads)
    for (uint32_t i = 0; i < input_files.size(); ++i) {
        uint32_t chromosome_id = get_chromosome(input_files[i]);
        if (std::find(chromosome_ids.begin(), chromosome_ids.end(), chromosome_id)
            == chromosome_ids.end()) {
            logger()->trace("Skipping {} (not in --chromosomes)", input_files[i]);
            continue;
        }
        std::tie(pos_data[chromosome_id], num_cells_chr[chromosome_id],
                 max_read_lengths[chromosome_id])
                = read_pileup(
                        input_files[i], id_to_group,
                        [&read_progress](uint32_t progress) { read_progress += progress; },
                        FLAGS_max_coverage, positions[chromosome_id], FLAGS_compute_read_stats);
    }
    uint32_t max_read_length = *std::max_element(max_read_lengths.begin(), max_read_lengths.end());
    uint32_t num_cells = *std::max_element(num_cells_chr.begin(), num_cells_chr.end());

    if (FLAGS_merge_file.empty()) {
        if (id_to_group.size() < num_cells) {
            logger()->error(
                    "--max_cell_count is {}, but number of cells is {}. Please add "
                    "--max_cell_count={} to the command line",
                    FLAGS_max_cell_count, num_cells, FLAGS_max_cell_count);
            std::exit(1);
        }
        // now that we know the actual number of cells, resize the mapping
        id_to_group.resize(num_cells);
    } else if (num_cells != id_to_group.size()) {
        logger()->error(
                "Invalid merge file {}. Merge files contains {} cell ids, data has {} cell ids",
                FLAGS_merge_file, id_to_group.size(), num_cells);
        std::exit(1);
    }

    uint32_t num_groups = *std::max_element(id_to_group.begin(), id_to_group.end()) + 1;

    // this maps the cell id to its actual position as we divide cells into smaller and smaller
    // clusters
    std::vector<uint32_t> cell_id_map(num_groups);
    std::iota(cell_id_map.begin(), cell_id_map.end(), 0);

    std::vector<uint16_t> clusters(num_cells); // contains the final clustering
    uint16_t cluster_idx = 1; // 0 means "no cluster"
    if (FLAGS_clustering.empty()) {
        divide_cluster(pos_data, max_read_length, id_to_group, cell_id_map, cell_id_map,
                       FLAGS_mutation_rate, FLAGS_homozygous_filtered_rate, FLAGS_seq_error_rate,
                       FLAGS_num_threads, FLAGS_o, FLAGS_normalization, FLAGS_termination,
                       FLAGS_clustering_type, FLAGS_arma_kmeans, FLAGS_expectation_maximization,
                       FLAGS_min_cluster_size, "", &clusters, &cluster_idx);
    } else {
        logger()->info("Using provided clustering file {}", FLAGS_clustering);
        clusters = int_split<uint16_t>(read_file(FLAGS_clustering), ',');
        if (clusters.size() != num_cells) {
            logger()->error("Number of clusters ({}) doesn't match number of cells ({})",
                            clusters.size(), num_cells);
        }
    }

    if (!FLAGS_reference_genome.empty()) {
        logger()->info("Performing variant calling against {}", FLAGS_reference_genome);
        variant_calling(pos_data, clusters, FLAGS_reference_genome, FLAGS_map_file, 1e-3,
                        FLAGS_seq_error_rate, FLAGS_o);
    } else {
        logger()->info("Skipping variant calling, because no reference genome was provided");
    }

    logger()->info("Done.");
}
