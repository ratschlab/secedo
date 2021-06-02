#pragma once

#include "sequenced_data.hpp"
#include "util/mat.hpp"

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
 * @param max_fragment_length the maximum DNA fragment length in the processed files (typically
 * around 600 bases). This value is used to figure out when reads have ended and can be processed. A
 * note on terminology: Illumina sequencing involves chopping the DNA into *fragments*, which vary
 * in size from ~50 to ~600. Each fragment is then being read from both ends (that's called a paired
 * end read). The machine reads 100 base paires from each end (so 200 total). If the DNA fragment
 * happens to be shorter than 200, then some reads will overlap. If the DAN fragment happens to be
 * longer than 200, then the middle of the fragment won't be processed. This unprocessed part is
 * called the "insert" and its length the "insert length".
 * @param cell_id_to_cell_idx re-maps the original cell ids (as they appear in #pos_data) so that
 * they form a continuous sequence starting with 0. At the first clustering step, this mapping is
 * simply the identity. At subsequent clustering steps, it maps the cells in the cluster to
 * 0..cell_count.
 * @param mutation_rate estimated mutation rate
 * @param homozygous_rate estimated probability that a loci is homozygous (it wasn't filtered
 * correctly due to sequencing errors being mistaken for heterozygosity)
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
                             uint32_t max_fragment_length,
                             const std::vector<uint32_t> &cell_id_to_cell_idx,
                             double mutation_rate,
                             double homozygous_rate,
                             double seq_error_rate,
                             const uint32_t num_threads,
                             const std::string &marker,
                             const std::string &normalization);
