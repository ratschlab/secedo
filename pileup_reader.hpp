#pragma once

#include "sequenced_data.hpp"

#include <string>
#include <unordered_set>

/**
 * Reads a pileup file (typically containing single cell data for one chromosome) and returns all
 * the bases read at each position.
 * @param fname file name to read from
 * @param coverage_factor artificially increases the coverage by grouping data from coverage_factor
 * consecutive cells together (this assumes that consecutive cell ids are part of the same cluster)
 * @return a tuple containing a vector with the reads at each position, all the cell ids and the
 * maximum read length
 */
std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t>
read_pileup(const std::string fname, uint16_t coverage_factor = 1);
