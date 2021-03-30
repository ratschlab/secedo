#pragma once

#include "sequenced_data.hpp"

#include <cstdint>
#include <string>
#include <vector>

/*
 * Expectation-Maximization applied to refine the result of spectral clustering, i.e. the
 * probability that each cell belongs to one of two clusters. Unlike a typical
 * Expectation-Maximization, our goal is not to estimate the parameters of the Mixture Model, i.e.
 * (μ_0, μ_1, Σ_0, Σ_1) in cause of a Gaussian mixture (also, this is not a Gaussian mixture), but
 * rather to refine the probability of each cell belonging to a cluster or the other - so the
 * computation is a bit simpler.
 * @param pos_data the list of bases (A,C,G,T) sequences for each chromosome at each position in the
 * genome; pos_data[i][j] contains all the reads for chromosome i at some position p_j (the position
 * itself is not stored, as it's not needed)
 * @param cell_id_to_cell_pos maps cell ids from #pos_data (ranging from 0 to number of sequenced
 * cells) to positions in prob_cluster_b
 * @param theta the sequencing error rate (typically 1e-3)
 * @param[in, out] the probability that cell i belongs to the second cluster; the initial value of
 * this vector is obtained via spectral clustering and is being now refined using expectation
 * maximization
 */
void expectation_maximization(const std::vector<std::vector<PosData>> &pos_data,
                              const std::vector<uint32_t> &cell_id_to_cell_pos,
                              double theta,
                              std::vector<double> *prob_cluster_b);
