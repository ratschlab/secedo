#pragma once

#include <armadillo>

#include <cstdint>
#include <vector>

using namespace std;


class KMeans {
  private:
    double inertia_;
    std::vector<uint32_t> labels_;

  public:
    double inertia() { return inertia_; }
    const std::vector<uint32_t> &labels() { return labels_; }

    /**
     * @param points the points to be clustered, one per row
     * @return the cluster for each point
     */
    std::vector<uint32_t>
    run(const arma::mat &points, uint32_t K, uint32_t max_iter, uint32_t num_tries);
};
