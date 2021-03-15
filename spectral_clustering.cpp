#include "util.hpp"

#include <cstdint>
#include <vector>

Matd laplacian(const Matd &a) {
    std::vector<double> diag(a.size());
    for (uint32_t r = 0; r < a.size(); ++r) {
        assert(a[r][r] == 0); // diagonal elements MUST be zero
        for (uint32_t c = 0; c < a[r].size(); ++c) {
            diag[r] += a[r][c];
        }
    }
    std::transform(diag.begin(), diag.end(), diag.begin(),
                   [](double v) { return 1 / std::sqrt(v); });
    Matd result = newMat(a.size(), a[0].size(), 0.0);
    for (uint32_t r = 0; r < a.size(); ++r) {
        for (uint32_t c = 0; c < a[r].size(); ++c) {
            result[r][c] = (r == c ? 1 : 0) - diag[r] * diag[c] * a[r][c];
        }
    }
    return result;
}

/*
void compute_aic_bic(const Matd &sim_mat, uint32_t index) {
    // compute graph Laplacian, the eigenvalues and eigenvectors
    constexpr bool normed = true;
    Matd L = laplacian(sim_mat, normed);
    eigenvalues, eigenvectors = np.linalg.eig(L);
    idx = eigenvalues.argsort();
    eigenvalues = eigenvalues[idx;] eigenvectors = eigenvectors[:, idx];
    // save the first 20 eigenvalues
with open("eigenvalues_" + str(index), 'w') as f:
f.write(','.join([str(x) for x in eigenvalues[:20]]));
// print the second and third smallest eigenVectors
with open("eigenvectors_" + str(index), 'w') as rf:
for i in range(sim_mat.shape[0]) {
    rf.write(str(eigenvectors[i, 1]) + "," + str(eigenvectors[i, 2]) + "\n");
}

// decide if it seems there are two different clusters
// coordinates of cells... value of the first eigenvector
// we assume the coordinates for cells within one cluster are normally distributed
// ->fit a Gaussian Mixture Model with 1 component and 2 components, compare using AIC, BIC
cell_coord = eigenvectors[:, [ 1, 2, 3, 4, 5 ]];
// cell_coord = eigenvectors[ :, 1].reshape(-1, 1)
// GMM with 1 component
gmm_1 = GaussianMixture(n_components = 1, n_init = 1, random_state = 29);
gmm_1.fit(cell_coord);
aic_1 = gmm_1.aic(cell_coord);
bic_1 = gmm_1.bic(cell_coord);

// GMM with 2 components
gmm_2 = GaussianMixture(n_components = 2, n_init = 3, random_state = 29);
gmm_2.fit(cell_coord);
aic_2 = gmm_2.aic(cell_coord);
bic_2 = gmm_2.bic(cell_coord);

// GMM with 2 components and tied variance
gmm_2_tied = GaussianMixture(n_components = 2, n_init = 3, random_state = 29,
                             covariance_type = 'tied');
gmm_2_tied.fit(cell_coord);
aic_2_tied = gmm_2_tied.aic(cell_coord);
bic_2_tied = gmm_2_tied.bic(cell_coord);

// k - means on the first five eigenvectors(excluding the trivial one)
// X = eigenvectors[ :, [1, 2, 3, 4, 5]]
// print(X.shape)
// kmeans = KMeans(n_clusters = 2, random_state = 0).fit(X)
// np.savetxt("labels_specClus_5eigen.csv", kmeans.labels_, delimiter = ',')
// return all the AIC 's and BIC' s
return aic_1, aic_2, aic_2_tied, bic_1, bic_2, bic_2_tied;
}
*/
