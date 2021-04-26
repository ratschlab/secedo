#include "sequenced_data.hpp"
#include "util/bam_reader.hpp"
#include "util/logger.hpp"
#include "util/util.hpp"

#include <gflags/gflags.h>

#include <cstdint>
#include <unordered_set>
#include <vector>

DEFINE_double(seq_error_rate, 0.01, "Sequencing errors rate, denoted by theta");
DEFINE_string(i, "", "Input directory containing BAM files");

DEFINE_string(o,
              "./",
              "File prefix where the pileup output will be written. Written files will be of the "
              "form <out_prefix>_<chromosome_id>.pileup.bin");

DEFINE_uint32(num_threads, 8, "Number of threads to use");

DEFINE_string(chromosomes,
              "The chromosomes on which to run the algorithm",
              "1,2,3,4,5,6,7,8,9,10,11,12,13.14,15,16,17,18,19,20,21,22,X,Y");

DEFINE_string(log_level,
              "trace",
              "The log verbosity: debug, trace, info, warn, error, critical, off");

DEFINE_uint32(min_base_quality,
              13,
              "Minimum Phred quality score for the sequencer. A quality score of 20 corresponds to "
              "a sequencing error rate of 0.01");


DEFINE_uint32(max_coverage,
              100,
              "Positions with higher coverage are considered anomalies and discarded during "
              "preprocessing (also at clustering time)");

uint32_t chromosome_to_id(const std::string &chromosome) {
    char *p;
    uint32_t converted = strtol(chromosome.c_str(), &p, 10);
    if (*p) {
        if (chromosome != "X" && chromosome != "Y") {
            logger()->error("Invalid chromosome: {}. Must be 1..22, X, Y", chromosome);
            std::exit(1);
        }
        return (chromosome == "X") ? 22 : 23;
    } else {
        return converted - 1;
    }
}

//============================================================================
int main(int argc, char *argv[]) {
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    spdlog::set_level(spdlog::level::from_str(FLAGS_log_level));

    std::vector<std::filesystem::path> input_files = { FLAGS_i };
    if (!std::filesystem::is_directory(FLAGS_i)) {
        logger()->error("-i <input_dir> must be a directory, not a file.");
        std::exit(1);
    }
    input_files = get_files(FLAGS_i, ".bam");
    if (input_files.empty()) {
        logger()->info("No BAM files found in {}. Done.", FLAGS_i);
        std::exit(0);
    }
    std::sort(input_files.begin(), input_files.end());
    logger()->info("Found {} input files in '{}'", input_files.size(), FLAGS_i);

    if (std::filesystem::is_directory(FLAGS_o)) {
        logger()->error("-o <output_dir> must be a file prefix, not a directory");
        std::exit(1);
    }

    for (const auto &chromosome : split(FLAGS_chromosomes, ',')) {
        auto output_file = std::filesystem::path(FLAGS_o) / ("_" + chromosome + ".pileup");
        read_bam(input_files, output_file, true, chromosome_to_id(chromosome), FLAGS_max_coverage,
                 FLAGS_min_base_quality, FLAGS_num_threads, FLAGS_seq_error_rate);
        logger()->trace("Done processing chromosome {}", chromosome);
    }
}
