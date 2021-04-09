#pragma once

#include "sequenced_data.hpp"

#include <string>
#include <unordered_set>

/**
 * Reads a pileup file (typically containing single cell data for one chromosome) and returns all
 * the bases read at each position.
 * @param fname file name to read from
 * @param merge_count artificially increases the coverage by grouping data from merge_count
 * consecutive cells together (this assumes that consecutive cell ids are part of the same cluster)
 * @return a tuple containing a vector with the reads at each position, all the cell ids and the
 * maximum read length
 */
std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t>
read_pileup(const std::string fname, const std::vector<uint16_t> &id_to_group);


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
