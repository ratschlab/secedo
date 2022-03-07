#include "spectral_clustering.hpp"
#include "expectation_maximization.hpp"
#include "similarity_matrix.hpp"

#include "util/is_significant.hpp"
#include "util/kmeans.hpp"
#include "util/logger.hpp"
#include "util/util.hpp"

#include <armadillo>

#include <cstdint>
#include <vector>

// set to true in order to log the unnormalized similarity matrices
constexpr bool log_matrices = false;

ClusteringType parse_clustering_type(const std::string &clustering_type) {
    if (clustering_type == "FIEDLER")
        return ClusteringType::FIEDLER;
    else if (clustering_type == "SPECTRAL2")
        return ClusteringType::SPECTRAL2;
    else if (clustering_type == "SPECTRAL6")
        return ClusteringType::SPECTRAL6;
    else
        throw std::invalid_argument(clustering_type);
}

Termination parse_termination(const std::string &str_termination) {
    return str_termination == "AIC" ? Termination::AIC : Termination::BIC;
}

Matd laplacian(const Matd &a) {
    std::vector<double> diag(a.rows());
    for (uint32_t r = 0; r < a.rows(); ++r) {
        assert(a(r, r) == 0); // diagonal elements MUST be zero
        for (uint32_t c = 0; c < a.cols(); ++c) {
            assert(a(r, c) == a(c, r));
            diag[r] += a(r, c);
        }
    }
    std::transform(diag.begin(), diag.end(), diag.begin(),
                   [](double v) { return v == 0 ? 0 : 1 / std::sqrt(v); });
    Matd result = Matd::zeros(a.rows(), a.rows());
    for (uint32_t r = 0; r < a.rows(); ++r) {
        for (uint32_t c = 0; c <= r; ++c) {
            result(r, c) = (r == c ? 1 : 0) - diag[r] * diag[c] * a(r, c);
            result(c, r) = result(r, c);
        }
    }
    return result;
}

/**
 * Returns the number of parameters for the given Gaussian Mixture Model. Assuming we have K
 * Gaussians and the data has dimension N, the parameters in a GMM are:
 *   - the means ( K * N values)
 *   - the covariances (K * N * (N + 1)/2 values, since the covariance matrix is symmetric)
 *   - the priors for each Gaussian minus one (since the priors must sum up to 1, the last one is
 *  fully determined by the previous priors (K-1)
 *
 * @return the number of parameters
 */
uint32_t num_params(const arma::gmm_full &gmm) {
    uint32_t data_dim = gmm.n_dims(); // the dimensionality of the data
    uint32_t cov_params = gmm.n_gaus() * data_dim * (data_dim + 1) / 2;
    uint32_t mean_params = gmm.n_gaus() * data_dim;
    return cov_params + mean_params + gmm.n_gaus() - 1;
}

/**
 * Akaike information criterion for the given gaussian mixture model and data.
 * @param gmm the GMM used to explain #data
 * @param data the data to compute the AIC for
 * @return the AIC for the given GMM and data, the lower the better
 */
double aic(const arma::gmm_full &gmm, const arma::mat &data) {
    return 2 * num_params(gmm) - 2 * data.n_cols * gmm.avg_log_p(data);
}

/**
 * Compute the probabilities of each data point according to a specific Gaussian in the GMM.
 *
 * @param gmm the Gaussian Mixture Model used to compute the probabilities for #data.
 * @param data the data points to compute probabilties for; each column represents one data point,
 * so data.n_rows == gmm.n_dims()
 * @param idx which of the Gaussians in the GMM to use for computing the probabilities. 0 <= idx <=
 * gmm.n_gaus()
 * @return a vector of size data.n_cols, representing the probabilities for each sample in data
 * according to Gaussian number g in the model
 */
std::vector<double>
get_probabilities(const arma::gmm_full &gmm, const arma::mat &data, uint32_t idx) {
    arma::Row<double> likelihoods = arma::exp(gmm.log_p(data, idx)) * gmm.hefts[idx];
    arma::Row<double> sum = arma::exp(gmm.log_p(data, 0)) * gmm.hefts[0];
    for (uint32_t i = 1; i < gmm.n_gaus(); ++i) {
        sum += arma::exp(gmm.log_p(data, i)) * gmm.hefts[i];
    }
    return arma::conv_to<std::vector<double>>::from(likelihoods / sum);
}

std::vector<double>
get_assignments(const arma::gmm_full &gmm, const arma::mat &data, uint32_t idx) {
    return arma::conv_to<std::vector<double>>::from(gmm.assign(data, arma::prob_dist));
}

/**
 * Bayesian information criterion for the given gaussian mixture model and data.
 * @param gmm the GMM used to explain #data
 * @param data the data to compute the AIC for
 * @return the BIC for the given GMM and data, the lower the better
 */
double bic(const arma::gmm_full &gmm, const arma::mat &data) {
    return num_params(gmm) * std::log(data.n_cols) - 2 * data.n_cols * gmm.avg_log_p(data);
}

uint32_t spectral_clustering(const Matd &similarity,
                             const ClusteringType &clustering,
                             const Termination &termination,
                             const std::string &out_dir,
                             const std::string &marker,
                             bool use_arma_kmeans,
                             std::vector<double> *cluster) {
    cluster->resize(similarity.rows());

    // compute graph Laplacian, the eigenvalues and eigenvectors
    Matd L = laplacian(similarity);
    arma::mat lap(similarity.rows(), similarity.rows());
    for (uint32_t r = 0; r < similarity.rows(); ++r) {
        for (uint32_t c = 0; c < similarity.rows(); ++c) {
            lap(r, c) = L(r, c);
        }
    }

    assert(lap.is_symmetric());
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, lap);

    // save the first 20 eigenvalues
    std::ofstream f(out_dir + "sim_mat_eigenvalues" + marker + ".csv");
    f << eigenvalues.cols(0ull, std::min(20ull, eigenvalues.n_cols - 1ull));
    f.close();

    if (log_matrices) {
        // write the eigenvectors to a file (useful for visualization)
        if (eigenvectors.n_cols > 2) {
            f.open(out_dir + "sim_mat_eigenvectors" + marker + ".csv");
            f << eigenvectors << std::endl;
            f.close();
        }
    }

    // decide if it seems there are two or more different clusters, by fitting a Gaussian Mixture
    // Model to the first 5 eigenvectors. We fit one model with 1 component, one model with 2
    // components, etc. and then we compare which model fits the data best using AIC and BIC we
    // assume the coordinates for cells within one cluster are normally distributed
    // ->fit a Gaussian Mixture Model with 1 component and 2 components, compare using AIC, BIC

    if (eigenvectors.n_cols < 2) {
        std::cout << "Perfect decomposition, all cells in same cluster" << std::endl;
        return true; // one single cluster, done
    }

    // the first 5 non-trivial eigenvectors used for GMM
    arma::mat cell_coord_1_5 = eigenvectors.cols(1, std::min(5ull, eigenvectors.n_cols - 1)).t();

    constexpr uint32_t max_clusters = 4;

    // the first 2 eigenvectors (including trivial), used for the k-means classification
    arma::mat cell_coord_0_2
            = eigenvectors.cols(0, std::min(2U, static_cast<uint32_t>(eigenvectors.n_cols - 1)));

    std::vector<arma::gmm_full> gmms(max_clusters);
    std::vector<double> aics(max_clusters);
    std::vector<double> bics(max_clusters);
    std::vector<double> gmm_probs(max_clusters);
    std::vector<double> inertia(max_clusters);
    std::vector<double> gaps(max_clusters - 1);
    std::vector<bool> statuses(max_clusters);

    for (uint32_t i = 0; i < max_clusters; ++i) {
        statuses[i] = gmms[i].learn(cell_coord_1_5, i + 1 /* components */, arma::eucl_dist,
                                    arma::random_subset, 10, 5, 1e-10, false);
        aics[i] = aic(gmms[i], cell_coord_1_5);
        bics[i] = bic(gmms[i], cell_coord_1_5);
        gmm_probs[i] = gmms[i].avg_log_p(cell_coord_1_5);
        KMeans kmeans;
        kmeans.run(cell_coord_0_2, i + 1, 100, 10);
        inertia[i] = kmeans.inertia();
    }

    // the gaps are only used to decide if we attempt to cluster in 2, 3 or 4 groups
    for (uint32_t i = 1; i < max_clusters; ++i) {
        gaps[i - 1] = inertia[i - 1] - inertia[i];
    }
    uint32_t cluster_count = 2;
    for (uint32_t i = 1; i < max_clusters - 1; ++i) {
        if (gaps[i] > 0.75 * gaps[i - 1]) {
            cluster_count = i + 2;
        } else {
            break;
        }
    }

    logger()->trace(
            "Avg log-likelihood for GMM 1/2/3/4 {}/{}/{}/{}\taic 1/2/3/4 {}/{}/{}/{}\tbic 1/2/3/4 "
            "{}/{}/{}/{}",
            gmm_probs[0], gmm_probs[1], gmm_probs[2], gmm_probs[3], aics[0], aics[1], aics[2],
            aics[3], bics[0], bics[1], bics[2], bics[3]);

    logger()->trace("Attempting to cluster into {} clusters", cluster_count);

    switch (clustering) {
        case ClusteringType::FIEDLER: {
            // TODO: use the min-sparsity cut described in
            // https://people.eecs.berkeley.edu/~jrs/189s17/lec/22.pdf rather than the 0 cut
            arma::vec fiedler = eigenvectors.col(1);
            std::sort(fiedler.begin(), fiedler.end());

            // corner case: if the smallest value is zero, treat all zero values as negative, and
            // all positive values as port of the other cluster
            double threshold = fiedler[0] == 0 ? std::numeric_limits<double>::min() : 0;
            for (uint32_t i = 0; i < similarity.rows(); ++i) {
                cluster->at(i) = eigenvectors(i, 1) >= threshold;
            }
            break;
        }
        case ClusteringType::SPECTRAL2:
        case ClusteringType::SPECTRAL6: {
            // TODO: this is totally weird, but it works better if considering the zeroth column,
            // at lest in the tests. Maybe this comment from scikit.SpectralClustering is relevant:
            // "The first eigenvector is constant only for fully connected graphs and should be kept
            // for spectral clustering"
            uint32_t col_idx = clustering == ClusteringType::SPECTRAL2 ? 2 : 6;
            arma::mat ev = eigenvectors.cols(
                    0, std::min(col_idx, static_cast<uint32_t>(eigenvectors.n_cols - 1)));
            // normalize each row
            for (uint32_t i = 0; i < ev.n_rows; ++i) {
                double norm = arma::norm(ev.row(i));
                if (norm > 0) {
                    ev.row(i) = ev.row(i) / norm;
                }
            }
            f.open(out_dir + "sim_mat_eigenvectors_norm" + marker + ".csv");
            f << ev << std::endl;
            f.close();

            if (use_arma_kmeans) {
                ev = ev.t(); // k-means expects each column to be one sample
                arma::mat means;
                bool status = arma::kmeans(means, ev, cluster_count, arma::random_spread,
                                           100 /* iterations */, false);
                if (!status) {
                    logger()->error("K-means clustering failed.");
                } else {
                    for (uint32_t i = 0; i < similarity.rows(); ++i) {
                        double min_dist = arma::norm(means.col(0) - ev.col(i));
                        uint32_t idx = 0;
                        for (uint32_t c = 1; c < cluster_count; ++c) {
                            double dist = arma::norm(means.col(c) - ev.col(i));
                            if (dist < min_dist) {
                                min_dist = dist;
                                idx = c;
                            }
                        }
                        cluster->at(i) = idx;
                    }

                    f.open(out_dir + "centroids" + marker + ".csv");
                    f << means << std::endl;
                }
            } else {
                KMeans kmeans;
                std::vector<uint32_t> labels = kmeans.run(ev, cluster_count, 100, 10);
                for (uint32_t i = 0; i < similarity.rows(); ++i) {
                    cluster->at(i) = labels[i];
                }
            }
            f.close();
            break;
        }
    }

    // stop if either we couldn't fit the 2/3-component GMMs or if the 1-component GMM fits better
    bool is_done = (!statuses[1] && !statuses[2] && !statuses[3]) || termination == Termination::AIC
            ? (aics[0] < std::min({ aics[1], aics[2], aics[3] }))
            : (bics[0] < std::min({ bics[1], bics[2], bics[3] }));
    if (is_done) {
        logger()->trace("Simple Gaussian Model matches data better - stopping the clustering");
        return 1;
    }
    for (uint32_t i = 0; i < max_clusters; ++i) {
        uint32_t count = std::count_if(cluster->begin(), cluster->end(),
                                       [i](double v) { return std::abs(v - i) < 1e-3; });
        logger()->trace("Cluster {} has {} cells", i + 1, count);
    }
    return cluster_count;
}

uint32_t get_chromosome(const std::filesystem::path &filename) {
    std::string fname = filename.filename().replace_extension().replace_extension();
    std::vector<std::string> parts = split(fname, '_');
    if (parts.size() != 2) {
        logger()->error("Invalid pileup filename. Must be <bla>_chromosome.*");
        std::exit(1);
    }
    return chromosome_to_id(parts[1]);
}

void divide_cluster(const std::vector<std::vector<PosData>> &pds,
                    uint32_t max_read_length,
                    const std::vector<uint16_t> &id_to_group,
                    const std::vector<uint32_t> &id_to_pos,
                    const std::vector<uint32_t> &pos_to_id,
                    double mutation_rate,
                    double homozygous_rate,
                    double seq_error_rate,
                    const uint32_t num_threads,
                    const std::string &out_dir,
                    const std::string &normalization,
                    const std::string &termination_str,
                    const std::string &clustering_type_str,
                    bool use_arma_kmeans,
                    bool use_expectation_maximization,
                    uint32_t min_cluster_size,
                    const std::string marker,
                    std::vector<uint16_t> *clusters,
                    uint16_t *cluster_idx) {
    if (!marker.empty()) {
        logger()->info("\n\nPerforming clustering of sub-cluster {} with {} elements", marker,
                       pos_to_id.size());
    }
    logger()->info("Filtering significant positions...");
    Filter filter(seq_error_rate);
    auto [pos_data, coverage] = filter.filter(pds, id_to_pos, marker, num_threads);
    std::ofstream filtered(std::filesystem::path(out_dir) / ("significant_positions" + marker));
    for (uint32_t i = 0; i < pos_data.size(); ++i) {
        for (uint32_t j = 0; j < pos_data[i].size(); ++j) {
            filtered << id_to_chromosome(i) << '\t' << pos_data[i][j].position << std::endl;
        }
    }
    filtered.close();

    if (coverage < 9) {
        logger()->trace("Coverage of cluster {} is lower than 9. Stopping.", marker);
        return;
    }

    logger()->info("Computing similarity matrix...");
    uint32_t n_groups_subcluster = pos_to_id.size();
    uint32_t n_groups_total = id_to_group.size();
    Matd sim_mat = computeSimilarityMatrix(pos_data, n_groups_subcluster, max_read_length,
                                           id_to_pos, mutation_rate, homozygous_rate,
                                           seq_error_rate, num_threads, marker, normalization);

    logger()->info("Performing spectral clustering...");
    std::vector<double> cluster; // size n_groups_subcluster
    Termination termination = parse_termination(termination_str);
    ClusteringType clustering_type = parse_clustering_type(clustering_type_str);
    uint32_t num_clusters = spectral_clustering(sim_mat, clustering_type, termination, out_dir,
                                                marker, use_arma_kmeans, &cluster);
    if (num_clusters == 1) {
        return;
    }

    std::vector<uint16_t> id_to_cluster(n_groups_total);
    for (uint16_t cell_id = 0; cell_id < n_groups_total; ++cell_id) {
        uint16_t pos = id_to_pos[id_to_group[cell_id]];
        id_to_cluster[cell_id] = (pos == NO_POS ? NO_POS : cluster[pos]);
    }
    write_vec(std::filesystem::path(out_dir) / ("spectral_clustering" + marker), id_to_cluster);

    if (use_expectation_maximization && num_clusters == 2) {
        logger()->info("Performing clustering refinement via expectation maximization...");
        expectation_maximization(pos_data, id_to_pos, num_threads, seq_error_rate, &cluster);
        for (uint16_t i = 0; i < n_groups_total; ++i) {
            uint32_t pos = id_to_pos[id_to_group[i]];
            id_to_cluster[i] = pos == NO_POS ? NO_POS : cluster[pos];
        }
        write_vec(std::filesystem::path(out_dir) / ("expectation_maximization" + marker),
                  id_to_cluster);
    } else {
        logger()->info("Skipping clustering refinement via expectation maximization...");
    }

    std::vector<std::vector<uint32_t>> pos_to_id_new(num_clusters);
    std::vector<std::vector<uint32_t>> id_to_pos_new(
            num_clusters, std::vector<uint32_t>(id_to_pos.size(), NO_POS));

    // TODO: move this computation out of the function
    std::unordered_map<uint16_t, std::vector<uint16_t>> group_id_to_cell_ids;
    for (uint32_t i = 0; i < id_to_group.size(); ++i) {
        group_id_to_cell_ids[id_to_group[i]].push_back(i);
    }
    for (uint32_t cell_idx = 0; cell_idx < n_groups_subcluster; ++cell_idx) {
        uint32_t group_id = pos_to_id[cell_idx];
        bool assigned = false;
        for (uint32_t c = 0; !assigned && c < num_clusters; ++c) {
            if (std::abs(cluster[cell_idx] - c) < 0.05) {
                id_to_pos_new[c][group_id] = pos_to_id_new[c].size();
                pos_to_id_new[c].push_back(group_id);
                for (uint32_t cell_id : group_id_to_cell_ids[group_id]) {
                    clusters->at(cell_id) = *cluster_idx + c;
                }
                assigned = true;
            }
        }
        if (!assigned) { // set to 0 cluster id of cells that couldn't be assigned to a cluster
            for (uint32_t cell_id : group_id_to_cell_ids[group_id]) {
                clusters->at(cell_id) = 0;
            }
        }
    }
    *cluster_idx += num_clusters;
    write_vec(std::filesystem::path(out_dir) / "clustering", *clusters);

    for (uint32_t c = 0; c < num_clusters; ++c) {
        if (pos_to_id_new[c].size() < min_cluster_size) {
            logger()->trace("Cluster {} size is too small ({} vs {}). Stopping.", c,
                            pos_to_id_new[c].size(), min_cluster_size);
        } else if (n_groups_subcluster - pos_to_id_new[c].size() < min_cluster_size) {
            logger()->trace("Cluster {} size is too large relative to total ({} vs {}). Stopping.",
                            c, pos_to_id_new[c].size(), n_groups_subcluster);
        } else {
            divide_cluster(pds, max_read_length, id_to_group, id_to_pos_new[c], pos_to_id_new[c],
                           mutation_rate, homozygous_rate, seq_error_rate, num_threads, out_dir,
                           normalization, termination_str, clustering_type_str, use_arma_kmeans,
                           use_expectation_maximization, min_cluster_size,
                           marker + static_cast<char>('A' + c), clusters, cluster_idx);
        }
    }
}
