#include "kmeans.hpp"

#include "logger.hpp"

#include <unordered_set>

double weighted_dist(const arma::rowvec &a, const ::arma::rowvec &b) {
    arma::rowvec diff = a - b;
    diff[1] *= 1.2;
    return arma::norm(diff);
}

double weighted_dist2(const arma::rowvec &a, const ::arma::rowvec &b) {
    arma::rowvec diff = a - b;
    diff[1] *= 1.2;
    return arma::dot(diff, diff);
}

uint32_t nearest_cluster(const std::vector<arma::rowvec> &centroids, const arma::rowvec &point) {
    double min_dist = weighted_dist(centroids[0], point);
    uint32_t result = 0;

    for (uint32_t i = 1; i < centroids.size(); i++) {
        double dist = weighted_dist(centroids[i], point);

        if (dist < min_dist) {
            min_dist = dist;
            result = i;
        }
    }

    return result;
}

/**
 * @param points the points to be clustered, one per row
 * @return the cluster for each point
 */
std::vector<uint32_t>
KMeans::run(const arma::mat &points, uint32_t K, uint32_t max_iter, uint32_t num_tries) {
    uint32_t n_points = points.n_rows;
    if (n_points == 0) {
        logger()->warn("No points to cluster");
        return {};
    }

    std::vector<uint32_t> idx_to_cluster(n_points);
    std::vector<uint32_t> best_idx_to_cluster(n_points);
    double min_inertia = std::numeric_limits<double>::max();
    std::vector<arma::rowvec> centroids(K, arma::rowvec(points.n_cols));

    std::default_random_engine generator;
    std::uniform_int_distribution<uint32_t> rand_idx(0, n_points);


    for (uint32_t run = 0; run < num_tries; ++run) {
        // Initializing Clusters
        std::unordered_set<uint32_t> used_pointIds;

        for (uint32_t i = 0; i < K; i++) {
            while (true) {
                uint32_t index = rand_idx(generator);

                if (used_pointIds.find(index) == used_pointIds.end()) {
                    used_pointIds.insert(index);
                    centroids[i] = points.row(i);
                    break;
                }
            }
        }

        for (uint32_t iter = 0; iter < max_iter; ++iter) {
            bool done = true;

            // Add all points to their nearest cluster
            for (uint32_t i = 0; i < n_points; i++) {
                uint32_t new_cluster = nearest_cluster(centroids, points.row(i));
                if (new_cluster != idx_to_cluster[i]) {
                    done = false;
                }
                idx_to_cluster[i] = new_cluster;
            }

            for (auto &c : centroids) {
                c.zeros();
            }
            std::vector<uint32_t> counts(K, 0);
            // Recalculating the center of each cluster
            for (uint32_t i = 0; i < n_points; i++) {
                centroids[idx_to_cluster[i]] += points.row(i);
                counts[idx_to_cluster[i]]++;
            }

            for (uint32_t i = 0; i < K; ++i) {
                if (counts[i] > 0) {
                    centroids[i] /= counts[i];
                }
            }

            if (done) {
                break;
            }
        }
        double inertia = 0;
        for (uint32_t i = 0; i < n_points; i++) {
            inertia += weighted_dist2(points.row(i), centroids[idx_to_cluster[i]]);
        }
        if (inertia < min_inertia) {
            best_idx_to_cluster = idx_to_cluster;
            min_inertia = inertia;
        }
    }
    labels_ = best_idx_to_cluster;
    inertia_ = min_inertia;
    return best_idx_to_cluster;
}
