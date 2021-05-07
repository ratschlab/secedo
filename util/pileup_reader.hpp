#pragma once

#include "sequenced_data.hpp"

#include <functional>
#include <string>
#include <unordered_set>

/**
 * Reads a pileup file (typically containing single cell data for one chromosome) and returns all
 * the bases read at each position.
 * @param fname file name to read from
 * @param sequencing_error_rate the error rate of the sequencer (used to filter out non-informative
 * positions)
 * @param merge_count artificially increases the coverage by grouping data from merge_count
 * consecutive cells together (this assumes that consecutive cell ids are part of the same cluster)
 * @param progress callback that is invoked to report the number of bytes processed (e.g. for
 * reporting progress in the caller)
 * @return a tuple containing:
 *    1. a vector with the reads at each position
 *    2. all the cell ids
 *    3. the longest DNA fragment length (this will help us later decide when a fragment is fully
 * processed)
 */
std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t> read_pileup(
        const std::string fname,
        double sequencing_error_rate,
        const std::vector<uint16_t> &id_to_group,
        const std::function<void(uint64_t)> &progress = [](uint32_t) {},
        uint32_t max_coverage = 100);


/**
 * Creates the cell grouping based on the merge_count and merge_file. If merge_count is 1 and no
 * merge_file is specified, the identity permutation (no grouping) is returned.
 * @param merge_count merge_count artificially increases the coverage by grouping data from
 * merge_count consecutive cells together (this assumes that consecutive cell ids are part of the
 * same cluster)
 * @param merge_file file name containing comma separated values indicating the group of each cell
 * @param max_cell_count max number of cells (used to initialize the identity mapping)
 * @return a vector v, where v[i] contains the group cell number i belongs to
 */
std::vector<uint16_t> get_grouping(uint16_t merge_count = 1,
                                   const std::string &merge_file = "",
                                   uint16_t max_cell_count = 10'000);
