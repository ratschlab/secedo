#pragma once

#include "sequenced_data.hpp"
#include "util/util.hpp"

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
};

ClusteringType parse_clustering_type(const std::string &clustering_type);

Termination parse_termination(const std::string &str_termination);

/**
 * Computes the symmetric normalized Laplacian of matrix a
 * @param a input matrix
 * @return the symmetric normalized Laplacian, I - D^0.5 * a * D^0.5, where D is the diagonal degree
 * matrix. The returned matrix is guaranteed to be symmetrical.
 */
Matd laplacian(const Matd &a);

/**
 * Clusters a group of points defined by the similarity matrix into 2 groups.
 * @param similarity square matrix with elements in [0,1]. similarity[i][j] measures the similarity
 * between point i and point j, where 1 means identical and 0 means not at all similar.
 * @param clustering the clustering variant to use
 * @param termination the criteria to use as a stopping criterion, e.g. the criterion to decide if
 * dividing into 2 clusters explains the data better than using a single cluster
 * @param use_arma_kmeans if true, use the armadillo library for clustering, otherwise use our own
 * primitive k-means clustering, which assigns a higher weight to the Fiedler vector
 * @param cluster a vector of size similarity.size(), with values 0 or 1 assigning clusters to each
 * point
 * @return the number of clusters that best explain the data (1 means that no further clustering is
 * possible)
 */
uint32_t spectral_clustering(const Matd &similarity,
                             const ClusteringType &clustering,
                             const Termination &termination,
                             const std::string &out_dir,
                             const std::string &marker,
                             bool use_arma_kmeans,
                             std::vector<double> *cluster);

/**
 * Recursively divides cells into 2 sub-clusters until a termination criteria is met.
 * N - number of cells
 * G - number of cell groups (normally N=G, but in case we group cells to artificially increase
 * coverage we have G < N)
 * NC - number of cells in the current sub-cluster. At the first call of divide() NC=G.
 * @param pos_data
 * @param max_read_length length of the longest fragment (typically around 500)
 * @param id_to_group of size N, maps cell ids to cell groups. Data from cells in the same group is
 * treated as if it came from one cell. Used to artificially increase coverage when testing
 * @param id_to_pos of size G, maps a cell group id to its position in the similarity matrix as we
 * subdivide into smaller and smaller clusters. At the beginning, this is the identity permutation.
 * If a cell with id 'cell_id' is not in the current cluster, then id_to_pos[cell_id]==NO_POS. The
 * position of a cell in the similarity matrix is given by id_to_pos[id_to_group[cell_id]].
 * @param pos_to_id of size NC the inverse of #id_to_pos, it maps each position 0...pos in the
 * current subgroup to the actual cell it corresponds to
 * @param mutation_rate epsilon, estimated frequency of mutated loci in the pre-processed data set
 * @param homozygous_rate  the probability that a locus is homozygous, (not filtered correctly in
 * the first step)
 * @param seq_error_rate error rate of the sequencer, e.g. 1e-3 if using Illumina reads with base
 * quality >30
 * @param num_threads number of threads to use for the computation
 * @param out_dir where to output the clustering results
 * @param normalization the type of normalization to use for the similiarity matrix (see the flag
 * with the same name)
 * @param marker marks the current sub-cluster; for example AB means we are in the second
 * sub-cluster (B) of the first cluster (A)
 * @param[in, out] clusters vector of size num_cells, will contain the final clustering assignment.
 * Positions marked as zero indicate that a cluster couldn't be assigned.
 * @param[in, out] cluster_idx keeps track of the current cluster index, in order to assign distinct
 * cluster values to each subcluster
 */
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
                    uint16_t *cluster_idx);

// TODO(ddanciu): this is brittle - just write the chromosome id into the binary pileup file
uint32_t get_chromosome(const std::filesystem::path &filename);
