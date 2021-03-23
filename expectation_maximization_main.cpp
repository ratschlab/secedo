#include "expectation_maximization.hpp"

#include "pileup_reader.hpp"
#include "util.hpp"

#include <gflags/gflags.h>

#include <fstream>
#include <vector>

DEFINE_string(mpileup_file,
              "",
              "Input file containing 'pileup' textual format from an alignment, as written by "
              "preprocessing.py");
DEFINE_string(labels_file, "", "Input file containing labels");

DEFINE_double(seq_error_rate, 0.001, "Sequencing errors rate, denoted by theta");

DEFINE_uint32(num_threads, 8, "Number of threads to use");

DEFINE_string(chromosomes,
              "The chromosomes on which to run the algorithm",
              "1,2,3,4,5,6,7,8,9,10,11,12,13.14,15,16,17,18,19,20,21,22,X");

int main(int argc, char *argv[]) {
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    std::vector<double> prob_cluster_b = double_split(read_file(FLAGS_labels_file), ',');
    const std::vector<std::string> chromosomes = split(FLAGS_chromosomes, ',');

    std::vector<std::vector<PosData>> pos_data(chromosomes.size());
    std::unordered_set<uint32_t> unused(chromosomes.size());

    // read the data from all chromosomes in parallel
#pragma omp parallel for num_threads(FLAGS_num_threads)
    for (uint32_t idx = 0; idx < chromosomes.size(); ++idx) {
        std::string chromosome = chromosomes[idx];
        std::tie(pos_data[idx], unused) = read_pileup(FLAGS_mpileup_file + chromosome);
    }

    expectation_maximization(pos_data, FLAGS_seq_error_rate, &prob_cluster_b);
}
