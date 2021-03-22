#pragma once

#include "sequenced_data.hpp"

#include <cstdint>
#include <vector>

/**
 * Compute mat_same (the matrix giving probabilities of cells i and j given they are in
 * the same cluster) and mat_diff (prob. of cells i and j given they are in different clusters)
 *
 * @param pos_data for each chromosome and each position, all the reads for that position (from all
 * cells)
 * @param num_cells number of cells that were sequenced
 * @param cell_ids
 * @param cell_groups
 * @param mutation_rate estimated mutation rate
 * @param heterozygous_rate estimated probability that a loci is heterozygous
 * @param seq_error_rate estimated error rate in the sequencing technology
 * @param num_threads number of threads to use for parallelizing the computation
 * @param mpileup_file a file that contains all the reads from all the cells at a
 * @param out_dir directory where the output files (mat_same, mat_diff, combs_xs_xd) will be written
 */
void computeSimilarityMatrix(const std::vector<std::vector<PosData>> &pos_data,
                             uint32_t num_cells,
                             const std::vector<uint32_t> &cell_ids,
                             const std::vector<uint32_t> &cell_groups,
                             double mutation_rate,
                             double heterozygous_rate,
                             double seq_error_rate,
                             const uint32_t num_threads,
                             const std::string &out_dir);
