#include "expectation_maximization.hpp"

#include "util/logger.hpp"
#include "util/util.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

/**
 * Computes the log cluster center for a given position in the chromosome.
 * @param pos_data all the reads from all cells (that have a read) at a certain position
 * @param prob_cluster the probability of each cell belonging to a cluster
 * @param theta sequencing error rate
 * @return the log of cluster center, obtained as a weighted sum of base counts (C_j[k,A] on page 31
 * of the master thesis)
 */
std::array<double, 4>
cluster_center(const PosData &pos_data, const std::vector<double> &prob_cluster, double theta) {
    std::array<double, 4> center = { 0, 0, 0, 0 };
    // sum the bases over all cells belonging to the cluster(weighted by prob_cluster_a)
    for (const CellData &cd : pos_data.cells_data) {
        center[cd.base] += prob_cluster[cd.cell_id];
    }
    // normalize the composition at each position
    double s = center[0] + center[1] + center[2] + center[3];
    if (s == 0) { // the probability is theta for each of ACGT, so normalized is 1/4
        return { log(0.25), log(0.25), log(0.25), log(0.25) };
    }
    std::transform(center.begin(), center.end(), center.begin(),
                   [&](double v) { return v / s > theta ? v / s : theta; });
    // re-normalize the composition
    s = center[0] + center[1] + center[2] + center[3];
    std::transform(center.begin(), center.end(), center.begin(),
                   [&](double v) { return log(v / s); });

    return center;
}

/**
 * The M step of the Expectation Maximization algorithm, computes the log likelihood of the observed
 * data given the probability of a cell being in each cluster.
 *
 * In a typical Expectation Maximization, this step involves computing our best new guess for the
 * model parameters (μ_0, μ_1, Σ_0, Σ_1) of the Gaussians that form the mixture model, and then
 * compute the likelihoods given the new μ and Σ. Since our distributions are not Gaussian and also
 * we don't care about these parameters, we just compute the log of the PDF for each cell in
 * the cluster.
 * @param prob_cluster_b the probability of each cell being in cluster b
 * @param pos_data contains the bases sequenced at each chromosome and each position
 * @param theta the probability of a sequencing error
 * @return the log likelihood for each cell being in cluster a or b, i.e.
 * {ln P(reads from cell i | C_0), ln P(reads from cell i | C_1)}
 */
void maximization_step(const std::vector<double> &prob_cluster_b,
                       const std::vector<PosData> &pos_data,
                       const std::vector<uint32_t> &id_to_pos,
                       double theta,
                       std::vector<double> *log_likelihood_a,
                       std::vector<double> *log_likelihood_b) {
    uint32_t n_cells = prob_cluster_b.size();

    // probability of each cell being in clusters a or b
    std::vector<double> prob_cluster_a(n_cells);
    for (uint32_t i = 0; i < n_cells; ++i) {
        prob_cluster_a[i] = 1 - prob_cluster_b[i];
    }

    log_likelihood_a->resize(n_cells);
    log_likelihood_b->resize(n_cells);

    // process position by position
    for (const PosData &pd : pos_data) {
        // compute the logarithm of the cluster centers
        std::array<double, 4> center_a = cluster_center(pd, prob_cluster_a, theta);
        std::array<double, 4> center_b = cluster_center(pd, prob_cluster_b, theta);

        // compute the log likelihoods given the new cluster centers
        for (const CellData &cd : pd.cells_data) {
            log_likelihood_a->at(id_to_pos[cd.cell_id]) += center_a[cd.base];
            log_likelihood_b->at(id_to_pos[cd.cell_id]) += center_b[cd.base];
        }
    }
}

/**
 * The E (expectation) step of expectation maximization, see
 * https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm#E_step
 *
 * Computes the probability that each cell belongs to one of the two clusters (T_2,i in the wiki
 * article above).
 * @param log_likelihood_a the logarithm of the first PDF evaluated for each cell, i.e.
 * ln P(reads from cell i | C_0)
 * @param log_likelihood_b the logarithm of the second PDF evaluated for each cell, i.e.
 * ln P(reads from cell i | C_1)
 * @param[in, out] prob_cluster_b the probability that a cell belongs to the 2nd cluster (T_2,i in
 * the wiki article)
 * @return true if the values in prob_cluster_b didn't change much, so the algorithm converged
 */
bool expectation_step(const std::vector<double> &log_likelihood_a,
                      const std::vector<double> &log_likelihood_b,
                      std::vector<double> *prob_cluster_b) {
    constexpr double EPSILON = 1e-2;
    assert(prob_cluster_b->size() == log_likelihood_a.size()
           && log_likelihood_a.size() == log_likelihood_b.size());

    uint32_t n_cells = log_likelihood_a.size();
    assert(n_cells > 0);

    // compute the priors, i.e. the probability that a point belongs to the first or 2nd cluster
    double prior_b = sum(*prob_cluster_b) / n_cells;
    double prior_a = 1 - prior_b;

    bool done = true;
    for (uint32_t i = 0; i < n_cells; ++i) {
        // the likelihood of cluster b divided by the likelihood of cluster a
        double odds_b_a = exp(std::clamp(log_likelihood_b[i] - log_likelihood_a[i], -100., 100.));
        double prob = 1 - 1 / (1 + odds_b_a * prior_b / prior_a);
        done &= std::abs(prob - (*prob_cluster_b)[i]) < EPSILON;
        (*prob_cluster_b)[i] = prob;
    }
    return done;
}

void expectation_maximization(const std::vector<std::vector<PosData>> &pos_data,
                              const std::vector<uint32_t> &id_to_pos,
                              uint32_t, // num_threads,
                              double theta,
                              std::vector<double> *prob_cluster_b) {
    std::vector<double> log_likelihood_a(prob_cluster_b->size(), 0);
    std::vector<double> log_likelihood_b(prob_cluster_b->size(), 0);
    do {
        // perform the maximization step for each chromosome in parallel
        // #pragma omp parallel for num_threads(num_threads)
        for (uint32_t idx = 0; idx < pos_data.size(); ++idx) {
            std::vector<double> log_likelihood_a_chr;
            std::vector<double> log_likelihood_b_chr;
            maximization_step(*prob_cluster_b, pos_data[idx], id_to_pos, theta,
                              &log_likelihood_a_chr, &log_likelihood_b_chr);
#pragma omp critical
            {
                for (uint32_t i = 0; i < log_likelihood_a.size(); ++i) {
                    log_likelihood_a[i] += log_likelihood_a_chr[i];
                    log_likelihood_b[i] += log_likelihood_b_chr[i];
                }
            }
        }

        // perform the expectation step and stop if the probabilities didn't change much
    } while (!expectation_step(log_likelihood_a, log_likelihood_b, prob_cluster_b));

    uint32_t count_a = std::count_if(prob_cluster_b->begin(), prob_cluster_b->end(),
                                     [](double v) { return v < 0.05; });
    uint32_t count_b = std::count_if(prob_cluster_b->begin(), prob_cluster_b->end(),
                                     [](double v) { return v > 0.95; });
    logger()->trace(
            "After refinement, first cluster has {} cells, second cluster has {} cells, {} "
            "cells are undecided",
            count_a, count_b, prob_cluster_b->size() - count_a - count_b);
}
