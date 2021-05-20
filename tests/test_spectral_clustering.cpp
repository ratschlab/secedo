#include "spectral_clustering.hpp"

#include <armadillo>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <random>

namespace {
using namespace testing;

class SpectralClustering
    : public testing::TestWithParam<std::tuple<ClusteringType, Termination, bool>> {};

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
    auto [clustering, termination, use_arma_kmeans] = GetParam();
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

        if (spectral_clustering(similarity, clustering, termination, "./", "", use_arma_kmeans,
                                &cluster)
            == 1) {
            count_done++;
        }
    }
    // the clustering should detect that there is a single cluster at least one out of 3 times
    ASSERT_TRUE(count_done > 1);
}

TEST_P(SpectralClustering, TwoClusters) {
    logger()->set_level(spdlog::level::trace);
    auto [clustering, termination, use_arma_kmeans] = GetParam();

    std::default_random_engine generator;
    std::uniform_int_distribution<uint32_t> dissimilar(0, 5);
    std::uniform_int_distribution<uint32_t> similar(100, 200);

    constexpr uint32_t num_cells = 100;

    std::vector<double> cluster;
    // we have 100 cells, with the first 50 and last 50 being identical (modulo some noise) to each
    // other
    Matd similarity = Matd::zeros(num_cells, num_cells);

    // set everything to dissimilar first
    for (uint32_t i = 0; i < num_cells; ++i) {
        for (uint32_t j = 0; j < i; ++j) {
            if (similar(generator) % 5) {
                similarity(i, j) = dissimilar(generator) + 20;
            } else {
                similarity(i, j) = dissimilar(generator);
            }
            similarity(j, i) = similarity(i, j);
        }
    }

    // set first half to similar
    const uint32_t half = num_cells / 2;
    for (uint32_t i = 0; i < half; ++i) {
        for (uint32_t j = 0; j < i; ++j) {
            if (similar(generator) % 2) {
                // 1st cluster
                similarity(i, j) = similar(generator);
                similarity(j, i) = similarity(i, j);
            } else {
                // 2nd cluster
                similarity(i + half, j + half) = similar(generator);
                similarity(j + half, i + half) = similarity(i + half, j + half);
            }
        }
    }

    uint32_t num_clusters = spectral_clustering(similarity, clustering, termination, "./", "",
                                                use_arma_kmeans, &cluster);

    ASSERT_EQ(num_clusters, 2);
    uint32_t mismatches = 0;
    for (uint32_t i = 0; i < half - 1; ++i) {
        if (std::abs(cluster[i] - cluster[i + 1]) > 1e-3) {
            mismatches++;
        }
    }
    for (uint32_t i = half; i < num_cells - 1; ++i) {
        if (std::abs(cluster[i] - cluster[i + 1]) > 1e-3) {
            mismatches++;
        }
    }
    ASSERT_LT(mismatches, 4);
    ASSERT_EQ(1., std::abs(cluster.front() - cluster.back()));
}

TEST_P(SpectralClustering, ThreeClusters) {
    logger()->set_level(spdlog::level::trace);
    auto [clustering, termination, use_arma_kmeans] = GetParam();

    if (clustering != ClusteringType::SPECTRAL2 && clustering != ClusteringType::SPECTRAL6) {
        return; // clustering in 3 doesn't really work well with the other methods
    }

    std::default_random_engine generator;
    std::uniform_int_distribution<uint32_t> dissimilar(0, 5);
    std::uniform_int_distribution<uint32_t> similar(100, 200);

    constexpr uint32_t num_cells = 99;

    std::vector<double> cluster;
    // we have 100 cells, with the first 50 and last 50 being identical (modulo some noise) to each
    // other
    Matd similarity = Matd::zeros(num_cells, num_cells);

    // set everything to dissimilar first
    for (uint32_t i = 0; i < num_cells; ++i) {
        for (uint32_t j = 0; j < i; ++j) {
            if (similar(generator) % 5) {
                similarity(i, j) = dissimilar(generator) + 20;
            } else {
                similarity(i, j) = dissimilar(generator);
            }
            similarity(j, i) = similarity(i, j);
        }
    }

    // set first half to similar
    const uint32_t third = num_cells / 3;
    for (uint32_t i = 0; i < third; ++i) {
        for (uint32_t j = 0; j < i; ++j) {
            if (similar(generator) % 2) {
                // 1st cluster
                similarity(i, j) = similar(generator);
                similarity(j, i) = similarity(i, j);
            } else {
                // 2nd cluster
                similarity(i + third, j + third) = similar(generator);
                similarity(j + third, i + third) = similarity(i + third, j + third);

                similarity(i + 2 * third, j + 2 * third) = similar(generator);
                similarity(j + 2 * third, i + 2 * third) = similarity(i + 2 * third, j + 2 * third);
            }
        }
    }

    uint32_t num_clusters = spectral_clustering(similarity, clustering, termination, "./", "",
                                                use_arma_kmeans, &cluster);
    ASSERT_TRUE(2 == num_clusters || 3 == num_clusters);

    // for 2 clusters check that the split is 33/66, for 3 clusters 33/33/33
    for (uint32_t i = 0; i < num_clusters; ++i) {
        uint32_t count = std::count(cluster.begin(), cluster.end(), i);
        // 0 or 33 or 66
        ASSERT_TRUE(labs(static_cast<int32_t>(count)) < 2
                    || labs(static_cast<int32_t>(count - 33U)) < 2
                    || labs(static_cast<int32_t>(count - 66U)) < 2);
    }
}

INSTANTIATE_TEST_SUITE_P(
        Method,
        SpectralClustering,
        ::testing::Values(std::make_tuple(ClusteringType::SPECTRAL2, Termination::AIC, false),
                          std::make_tuple(ClusteringType::FIEDLER, Termination::AIC, false),
                          std::make_tuple(ClusteringType::GMM_ASSIGN, Termination::AIC, false),
                          std::make_tuple(ClusteringType::GMM_PROB, Termination::AIC, false),
                          std::make_tuple(ClusteringType::SPECTRAL2, Termination::BIC, false),
                          std::make_tuple(ClusteringType::FIEDLER, Termination::BIC, false),
                          std::make_tuple(ClusteringType::GMM_PROB, Termination::BIC, false),
                          std::make_tuple(ClusteringType::GMM_ASSIGN, Termination::BIC, false),
                          std::make_tuple(ClusteringType::SPECTRAL2, Termination::AIC, true),
                          std::make_tuple(ClusteringType::FIEDLER, Termination::AIC, true),
                          std::make_tuple(ClusteringType::GMM_ASSIGN, Termination::AIC, true),
                          std::make_tuple(ClusteringType::GMM_PROB, Termination::AIC, true),
                          std::make_tuple(ClusteringType::SPECTRAL2, Termination::BIC, true),
                          std::make_tuple(ClusteringType::FIEDLER, Termination::BIC, true),
                          std::make_tuple(ClusteringType::GMM_PROB, Termination::BIC, true),
                          std::make_tuple(ClusteringType::GMM_ASSIGN, Termination::BIC, true)));

} // namespace
