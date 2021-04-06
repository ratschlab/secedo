#pragma once

#include "mat.hpp"
#include "sequenced_data.hpp"

#include <cstdint>
#include <vector>

enum class Normalization {
    // add the minimum value to the matrix
    ADD_MIN,
    // compuate 1/(1+exp(mat_sim)) (equivalent to Equation 5.6 in the thesis)
    EXPONENTIATE,
    // scale the exponentiated version so that max. element (excluding diagonal) is 1
    SCALE_MAX_1,

};

/**
 * Compute mat_same (the matrix giving probabilities of cells i and j given they are in
 * the same cluster) and mat_diff (prob. of cells i and j given they are in different clusters)
 *
 * @param pos_data for each chromosome and each position, all the reads for that position (from all
 * cells)
 * @param max_read_length the maximum read length in the processed files (typically around 600
 * bases). This value is used to figure out when reads have ended and can be processed.
 * @param cell_id_to_cell_idx
 * @param mutation_rate estimated mutation rate
 * @param heterozygous_rate estimated probability that a loci is heterozygous
 * @param seq_error_rate estimated error rate in the sequencing technology
 * @param num_threads number of threads to use for parallelizing the computation
 * @param mpileup_file a file that contains all the reads from all the cells at a
 * @param out_dir directory where the output files (mat_same, mat_diff, combs_xs_xd) will be written
 * @param marker string that marks the subset that is currently being processed, e.g. "AAB"
 * @param normalization the type of normalization to perform on the similarity matrix
 *
 * @return a pair containing the
 */
Matd computeSimilarityMatrix(const std::vector<std::vector<PosData>> &pos_data,
                             uint32_t num_cells,
                             uint32_t max_read_length,
                             const std::vector<uint32_t> &cell_id_to_cell_idx,
                             double mutation_rate,
                             double heterozygous_rate,
                             double seq_error_rate,
                             const uint32_t num_threads,
                             const std::string &out_dir,
                             const std::string &marker,
                             const std::string &normalization);
