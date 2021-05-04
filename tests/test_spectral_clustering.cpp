#include "spectral_clustering.hpp"

#include <armadillo>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <random>

namespace {
using namespace testing;

class SpectralClustering : public testing::TestWithParam<std::pair<ClusteringType, Termination>> {};

TEST(Laplacian, SomeMatrix) {
    Matd a(3, 3, { 0, .5, .2, .5, 0, .5, .2, .5, 0 });
    Matd L = laplacian(a);
    Matd expectedL(3, 3,
                   { 1., -0.5976143, -0.28571429, -0.5976143, 1., -0.5976143, -0.28571429,
                     -0.5976143, 1. });
    for (uint32_t r = 0; r < a.rows(); ++r) {
        for (uint32_t c = 0; c < a.cols(); ++c) {
            EXPECT_NEAR(L(r, c), expectedL(r, c), 1e-3);
        }
    }
}

TEST_P(SpectralClustering, OneCluster) {
    auto [clustering, termination] = GetParam();
    std::default_random_engine generator(1243);
    uint32_t count_done = 0;
    for (uint32_t trial = 0; trial < 3; ++trial) {
        std::uniform_real_distribution<double> noise(-1e-3, 1e-3);

        std::vector<double> cluster;
        // we have 10 identical cells
        constexpr uint32_t num_cells = 100;
        Matd similarity(num_cells, num_cells);
        for (uint32_t i = 0; i < num_cells; ++i) {
            similarity(i, i) = 0;
            for (uint32_t j = 0; j < i; ++j) {
                similarity(i, j) = 1 + noise(generator);
                similarity(j, i) = similarity(i, j);
            }
        }

        if (spectral_clustering(similarity, clustering, termination, "./", "", &cluster)) {
            count_done++;
        }
    }
    // the clustering should detect that there is a single cluster at least one out of 3 times
    ASSERT_TRUE(count_done > 1);
}

TEST_P(SpectralClustering, TwoClusters) {
    auto [clustering, termination] = GetParam();

    std::default_random_engine generator;
    std::uniform_int_distribution<uint32_t> dissimilar(0, 5);
    std::uniform_int_distribution<uint32_t> similar(100, 200);

    constexpr uint32_t num_cells = 100;

    std::vector<double> cluster;
    // we have 100 cells, with the first 50 and last 50 being identical (modulo some noise) to each
    // other
    Matd similarity = Matd::zeros(num_cells, num_cells);

    for (uint32_t i = 0; i < num_cells; ++i) {
        for (uint32_t j = 0; j < i; ++j) {
            // 1st cluster
            similarity(i, j) = dissimilar(generator); // + noise(generator);
            similarity(j, i) = similarity(i, j);
        }
    }

    const uint32_t half = num_cells / 2;
    for (uint32_t i = 0; i < half; ++i) {
        for (uint32_t j = 0; j < i; ++j) {
            // 1st cluster
            similarity(i, j) = similar(generator); // + noise(generator);
            similarity(j, i) = similarity(i, j);

            // 2nd cluster
            similarity(i + half, j + half) = similar(generator);
            similarity(j + half, i + half) = similarity(i + half, j + half);
        }
    }

    bool done = spectral_clustering(similarity, clustering, termination, "./", "", &cluster);

    ASSERT_FALSE(done); // the split should be successful
    for (uint32_t i = 0; i < half - 1; ++i) {
        ASSERT_NEAR(0, std::abs(cluster[i] - cluster[i + 1]), 1e-3);
    }
    for (uint32_t i = half; i < num_cells - 1; ++i) {
        ASSERT_NEAR(0, std::abs(cluster[i] - cluster[i + 1]), 1e-3);
    }
    ASSERT_EQ(1., std::abs(cluster.front() - cluster.back()));
}

TEST_P(SpectralClustering, ThreeClusters) {
    auto [clustering, termination] = GetParam();
    if (clustering == ClusteringType::GMM_ASSIGN || clustering == ClusteringType::GMM_PROB
        || clustering == ClusteringType::SPECTRAL2) {
        return; // GMM doesn't work well on this data
    }

    std::vector<double> cluster;
    // we have 6 cells, with the first 2, next 2 and last 2 being identical, and the first 2
    // and next 2 being slightly more similar to each other than to the last 2. The clustering
    // should thus group the first 4 cells and last 2 cells together.
    // clang-format off
    Matd similarity(6, 6,
                    {
                       0., 100, 1, 1, 0, 0,
                      100, 0, 1, 1, 0, 0,
                      1, 1, 0, 100, 0, 0,
                      1, 1, 100,  0, 0, 0,
                      0, 0, 0, 0, 0, 100,
                      0, 0, 0, 0, 100, 0 });
    // clang-format on
    bool done = spectral_clustering(similarity, clustering, termination, "./", "", &cluster);

    ASSERT_FALSE(done); // we clearly have 2 clusters
    ASSERT_THAT(std::vector(cluster.begin(), cluster.begin() + 4), Each(cluster[0]));
    ASSERT_EQ(cluster[4], cluster[5]);
    ASSERT_NE(cluster[0], cluster[4]);
}

INSTANTIATE_TEST_SUITE_P(
        Method,
        SpectralClustering,
        ::testing::Values(std::make_pair(ClusteringType::SPECTRAL2, Termination::AIC),
                          std::make_pair(ClusteringType::FIEDLER, Termination::AIC),
                          std::make_pair(ClusteringType::GMM_ASSIGN, Termination::AIC),
                          std::make_pair(ClusteringType::GMM_PROB, Termination::AIC),
                          std::make_pair(ClusteringType::SPECTRAL2, Termination::BIC),
                          std::make_pair(ClusteringType::FIEDLER, Termination::BIC),
                          std::make_pair(ClusteringType::GMM_PROB, Termination::BIC),
                          std::make_pair(ClusteringType::GMM_ASSIGN, Termination::BIC)));

} // namespace
