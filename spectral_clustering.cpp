#include "spectral_clustering.hpp"

#include "logger.hpp"
#include "util.hpp"

#include <armadillo>

#include <cstdint>
#include <vector>

Matd laplacian(const Matd &a) {
    std::vector<double> diag(a.rows());
    for (uint32_t r = 0; r < a.rows(); ++r) {
        assert(a(r, r) == 0); // diagonal elements MUST be zero
        for (uint32_t c = 0; c < a.cols(); ++c) {
            diag[r] += a(r, c);
        }
    }
    std::transform(diag.begin(), diag.end(), diag.begin(),
                   [](double v) { return 1 / std::sqrt(v); });
    Matd result = Matd::zeros(a.rows(), a.rows());
    for (uint32_t r = 0; r < a.rows(); ++r) {
        for (uint32_t c = 0; c < a.cols(); ++c) {
            result(r, c) = (r == c ? 1 : 0) - diag[r] * diag[c] * a(r, c);
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

bool spectral_clustering(const Matd &similarity,
                         const ClusteringType &clustering,
                         const Termination &termination,
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
    std::ofstream f("laplacian_eigenvalues");
    f << eigenvalues(0ull, std::min(20ull, eigenvalues.n_cols - 1));
    f.close();

    // print the second and third smallest eigenVectors
    if (eigenvectors.n_cols > 2) {
        f.open("laplacian_eigenvectors");
        f << eigenvectors.cols(0, 2) << std::endl;
        f.close();
    }

    // decide if it seems there are two different clusters, by fitting a Gaussian Mixture Model to
    // the first 5 eigenvectors. We fit one model with 1 component, and one model with 2 components,
    // and then we compare which model fits the data best using AIC and BIC we assume the
    // coordinates for cells within one cluster are normally distributed
    // ->fit a Gaussian Mixture Model with 1 component and 2 components, compare using AIC, BIC

    if (eigenvectors.n_cols < 2) {
        std::cout << "Perfect decomposition, all cells in same cluster" << std::endl;
        return true; // one single cluster, done
    }

    // the first 5 non-trivial eigenvectors
    arma::mat cell_coord = eigenvectors.cols(1, std::min(5ull, eigenvectors.n_cols - 1)).t();
    // std::cout << cell_coord << std::endl << std::endl;
    // GMM with 1 component
    arma::gmm_full gmm1;
    gmm1.learn(cell_coord, 1 /* components */, arma::eucl_dist, arma::random_subset, 10, 5, 1e-10,
               false);

    double aic1 = aic(gmm1, cell_coord);
    double bic1 = bic(gmm1, cell_coord);

    // GMM with 2 components
    arma::gmm_full gmm2;

    // Using Mahalanobis distance here, but worth trying Euclidean on real data
    bool status2 = gmm2.learn(cell_coord, 2 /* components */, arma::maha_dist, arma::random_subset,
                              10, 5, 1e-10, false);

    double aic2 = aic(gmm2, cell_coord);
    double bic2 = bic(gmm2, cell_coord);

    logger()->trace("Avg log-likelihood for GMM 1/2 {}/{}\taic 1/2 {}/{}\tbic 1/2 {}/{}",
                    gmm1.avg_log_p(cell_coord), gmm2.avg_log_p(cell_coord), aic1, aic2, bic1, bic2);

    // TODO: investigate using gmm with tied variances as in the Python version

    if (clustering == ClusteringType::GMM_ASSIGN) {
        *cluster = get_assignments(gmm2, cell_coord, 0);
    } else if (clustering == ClusteringType::GMM_PROB) {
        *cluster = get_probabilities(gmm2, cell_coord, 0);
    } else if (clustering == ClusteringType::SPECTRAL2) {
        arma::mat ev = eigenvectors.cols(0, std::min(1ull, eigenvectors.n_cols - 1));
        // normalize each row
        for (uint32_t i = 0; i < ev.n_rows; ++i) {
            double norm = arma::norm(ev.row(i));
            if (norm > 0) {
                ev.row(i) = ev.row(i) / norm;
            }
        }
        ev = ev.t(); // kmeans expects each column to be one sample
        arma::mat means;
        arma::kmeans(means, ev, 2, arma::random_spread, 10 /* iterations */, false);
        for (uint32_t i = 0; i < similarity.rows(); ++i) {
            (*cluster)[i]
                    = arma::norm(means.col(0) - ev.col(i)) > arma::norm(means.col(1) - ev.col(i));
        }
    } else if (clustering == ClusteringType::FIEDLER) {
        // TODO: use the min-sparsity cut described in
        // https://people.eecs.berkeley.edu/~jrs/189s17/lec/22.pdf rather than the 0 cut
        for (uint32_t i = 0; i < similarity.rows(); ++i) {
            cluster->at(i) = eigenvectors(i, 1) >= 0;
        }
    }

    // stop if either we couldn't fit the 2-component GMM or if the 1-component GMM fits better
    bool is_done = !status2 || termination == Termination::AIC ? aic1 < aic2 : bic1 < bic2;
    if (is_done) {
        logger()->trace("Simple Gaussian Model matches data better - stopping the clustering");
    } else {
        uint32_t count = std::count(cluster->begin(), cluster->end(), 0);
        logger()->trace("First cluster has {} cells, second cluster has {} cells", count,
                        cluster->size() - count);
    }
    return is_done;
}
