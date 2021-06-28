#include "pileup.hpp"
#include "sequenced_data.hpp"
#include "util/logger.hpp"
#include "util/util.hpp"

#include <gflags/gflags.h>

#include <cstdint>
#include <unordered_set>
#include <vector>

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

DEFINE_uint32(min_different,
              3,
              "Minimum number of bases that need to differ from the majority in order to add the "
              "locus to the pileup");

//============================================================================
int main(int argc, char *argv[]) {
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    spdlog::set_level(spdlog::level::from_str(FLAGS_log_level));

    std::vector<std::filesystem::path> input_files = { FLAGS_i };
    if (std::filesystem::is_directory(FLAGS_i)) {
        input_files = get_files(FLAGS_i, ".bam");
        if (input_files.empty()) {
            logger()->info("No BAM files found in {}. Done.", FLAGS_i);
            std::exit(0);
        }
        std::sort(input_files.begin(), input_files.end());
        logger()->info("Found {} input files in '{}'", input_files.size(), FLAGS_i);

        // write the mapping between the actual cell id and the cell idx used in the pileup
        auto f_map_name = std::filesystem::path(FLAGS_o + "_" + FLAGS_chromosomes + ".map");
        std::ofstream f_map(f_map_name);
        for (uint32_t i = 0; i < input_files.size(); ++i) {
            std::string fname = input_files[i].filename().replace_extension().string();
            f_map << fname.substr(0, fname.rfind("_")) << "\t" << i << std::endl;
        }
    }

    if (std::filesystem::is_directory(FLAGS_o)) {
        logger()->error("-o <output_dir> must be a file prefix, not a directory");
        std::exit(1);
    }

    for (const auto &chromosome : split(FLAGS_chromosomes, ',')) {
        auto output_file = std::filesystem::path(FLAGS_o + "_" + chromosome + ".pileup");
        pileup_bams(input_files, output_file, true, chromosome_to_id(chromosome),
                    FLAGS_max_coverage, FLAGS_min_base_quality, FLAGS_num_threads,
                    FLAGS_min_different);
        logger()->trace("Done processing chromosome {}", chromosome);
    }
}
