#pragma once

#include "sequenced_data.hpp"

#include <functional>
#include <string>
#include <unordered_set>

/**
 * Reads a pileup file (typically containing single cell data for one chromosome) and returns all
 * the bases read at each position.
 * @param fname file name to read from
 * @param merge_count artificially increases the coverage by grouping data from merge_count
 * consecutive cells together (this assumes that consecutive cell ids are part of the same cluster)
 * @param progress callback that is invoked to report the number of bytes processed (e.g. for
 * reporting progress in the caller)
 * @param max_coverage remove positions with coverage higher than max_coverage
 * @param positions if not empty, only consider positions listed here
 * @return a tuple containing:
 *    1. a vector with the reads at each position
 *    2. all the cell ids
 *    3. the longest DNA fragment length (this will help us later decide when a fragment is fully
 * processed)
 *
 * IMPLEMENTATION NOTE: We chose to read all pileup data into memory once, and then apply the
 * filtering step in memory for each clustering step. A more economical approach is to apply the
 * filtering (#is_significant()) when reading and only load into memory the positions that are
 * significant for the current cluster. This would reduce memory usage by a factor of 20-30x (for
 * the artificial dataset with 2500 cells the filtered data is 5GB vs 135GB) at the cost of having
 * to read and re-filter all the data again at every sub-clustering step (about 12 minutes for 2500
 * cells at coverage 0.04x and ETH's slow network disk).
 */
std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t> read_pileup(
        const std::string fname,
        const std::vector<uint16_t> &id_to_group,
        const std::function<void(uint64_t)> &progress = [](uint32_t) {},
        uint32_t max_coverage = 100,
        const std::vector<uint32_t> positions = {});


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
