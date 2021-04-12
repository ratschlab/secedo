#pragma once

#include "util.hpp"

/**
 * Whether to use the Akaike Information Criterion or the Bayesian Information criterion to decide
 * between the 1-component and the 2-component Gaussian Mixture Models in the termination criterion.
 */
enum class Termination { AIC, BIC };

/**
 * How to cluster the cells:
 *   FIEDLER - compute the Fiedler vector (2nd eigenvector) of the similarity matrix Laplacian, and
 * place cells corresponding to positive values in one cluster, the others in a different cluster
 *   SPECTRAL5 - perform k-Means clustering on the first 5 eigenvectors of the Laplacian (spectral
 * clustering)
 *   GMM - use the log likelihoods given by one of the 2 Gaussians in the 2-component Gaussian
 * Mixture Model used as a termination condition
 */
enum class ClusteringType {
    /**
     * Use only the Fiedler vector for clustering.
     */
    FIEDLER,
    /**
     * Spectral clustering with the first 2 eigenvectors
     */
    SPECTRAL2,
    /**
     * Spectral clustering with the first 6 eigenvectors
     */
    SPECTRAL6,
    /**
     * Use the probabilities given by the Gaussian Mixture Model
     */
    GMM_PROB,
    /**
     * Use the discrete assignments given by the Gaussian Mixture Model (choose the index of the
     * Gaussian whose center/mean is the closest to the data point)
     */
    GMM_ASSIGN,

};

/**
 * Computes the symmetric normalized Laplacian of matrix a
 * @param a input matrix
 * @return the symmetric normalized Laplacian
 */
Matd laplacian(const Matd &a);

/**
 * Clusters a group of points defined by the similarity matrix into 2 groups.
 * @param similarity square matrix with elements in [0,1]. similarity[i][j] measures the similarity
 * between point i and point j, where 1 means identical and 0 means not at all similar.
 * @param clustering the clustering variant to use
 * @param termination the criteria to use as a stopping criterion, e.g. the criterion to decide if
 * dividing into 2 clusters explains the data better than using a single cluster
 * @param cluster a vector of size similarity.size(), with values 0 or 1 assigning clusters to each
 * point
 * @return true if clustering into 2 groups explains the data better, false otherwise
 */
bool spectral_clustering(const Matd &similarity,
                         const std::string &clustering,
                         const Termination &termination,
                         std::vector<double> *cluster);
