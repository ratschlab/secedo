#include "expectation_maximization.hpp"
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

std::vector<PosData> read_pileup_files(const std::string fname) {
    std::vector<PosData> result;

    std::ifstream g(fname);
    std::string line;
    // process position by position
    while (std::getline(g, line)) {
        // parse the data, split by tabs
        std::vector<std::string> splitLine = split(line, '\t');

        // 4th column: bases/reads for the given position
        std::string bases = splitLine[3];
        // 5th column: cell ids the cell from which the reads are coming, comma separated
        std::vector<uint32_t> cell_ids = int_split(splitLine[4], ',');

        std::vector<CellData> cells_data;
        for (uint32_t j = 0; j < bases.size(); ++j) {
            uint16_t cell_index = cell_ids[j] / 4; // TODO: remove this once coverage works
            cells_data.push_back({ cell_index, CharToInt[(uint8_t)bases[j]] });
        }
        result.push_back({ cells_data });
    }
    return result;
}


int main(int argc, char *argv[]) {
    std::vector<double> prob_cluster_b = double_split(read_file(FLAGS_labels_file), ',');
    const std::vector<std::string> chromosomes = split(FLAGS_chromosomes, ',');

    std::vector<std::vector<PosData>> pos_data(chromosomes.size());

    // read the data from all chromosomes in parallel
#pragma omp parallel for num_threads(FLAGS_num_threads)
    for (uint32_t idx = 0; idx < chromosomes.size(); ++idx) {
        std::string chromosome = chromosomes[idx];
        pos_data[idx] = read_pileup_files(FLAGS_mpileup_file + chromosome);
    }

    expectation_maximization(pos_data, FLAGS_seq_error_rate, &prob_cluster_b);
}
