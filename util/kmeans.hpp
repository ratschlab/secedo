#pragma once

#include <armadillo>

#include <cstdint>
#include <vector>

using namespace std;


/**
 * @param points the points to be clustered, one per row
 * @return the cluster for each point
 */
std::vector<uint32_t> kmeans(const arma::mat &points, uint32_t K, uint32_t max_iter);
