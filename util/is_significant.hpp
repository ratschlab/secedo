#pragma once

#include "sequenced_data.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>


/**
 * Computes the log factorial using a table for small values or Stirling's formula for larger
 * values.
 */
double log_fact(uint32_t n);

/**
 * Decides if a given position is worth keeping, i.e. it will be useful in distinguishing cell
 * genotypes.
 * @param base_count counts of A,C,G, and T in the pooled data at a fixed position
 * @param theta sequencing error rate (e.g. ~0.01 on Illumina machines)
 * @return true if the position is kept
 */
bool is_significant(std::array<uint16_t, 4> &base_count, double theta);

bool is_significant(const PosData &pos_data, double theta, uint16_t *coverage);
