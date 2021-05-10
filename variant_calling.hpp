#pragma once

#include "expectation_maximization.hpp"
#include "util/util.hpp"


#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

constexpr uint8_t NO_GENOTYPE = std::numeric_limits<uint8_t>::max();

/**
 * Find the most likely genotype, if we observe bases counts
 * @parm nBases (array of length 4: number of As, Cs, etc.)
 * @param heteroPrior the prior on heterozygous genotype
 * @param theta sequencing error rate
 */
uint8_t mostLikelyGenotype(const std::array<uint32_t, 4> &nBases, double heteroPrior, double theta);

void variant_calling(const std::vector<std::vector<PosData>> &pos_data,
                     const std::vector<double> cell_labels,
                     double hetero_prior,
                     double theta,
                     const std::filesystem::path &out_dir);
