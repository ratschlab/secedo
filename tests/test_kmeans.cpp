#include "util/kmeans.hpp"

#include "util/util.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {
using namespace ::testing;
TEST(kmeans, Empty) {
    ASSERT_TRUE(kmeans(arma::mat(), 2, 100).empty());
}

TEST(kmeans, OnePoint) {
    arma::mat points = { { 1, 2, 3 } };
    std::vector<uint32_t> clusters = kmeans(points, 1, 100);
    ASSERT_EQ(1, clusters.size());
    ASSERT_EQ(0, clusters[0]);
}

TEST(kmeans, TwoPointsTwoClusters) {
    arma::mat points = { { 1, 2, 3 }, { 4, 5, 6 } };
    std::vector<uint32_t> clusters = kmeans(points, 2, 100);
    ASSERT_EQ(2, clusters.size());
    ASSERT_EQ(0, clusters[0]);
    ASSERT_EQ(1, clusters[1]);
}

TEST(kmeans, TwoPointsOneCluster) {
    arma::mat points = { { 1, 2, 3 }, { 4, 5, 6 } };
    std::vector<uint32_t> clusters = kmeans(points, 1, 100);
    ASSERT_EQ(2, clusters.size());
    ASSERT_EQ(0, clusters[0]);
    ASSERT_EQ(0, clusters[1]);
}

TEST(kmeans, ThreePointsTwoClusters) {
    arma::mat points = { { 1, 2, 3 }, { 4, 5, 6 }, { 1.1, 2, 3 } };
    std::vector<uint32_t> clusters = kmeans(points, 2, 100);
    ASSERT_EQ(3, clusters.size());
    ASSERT_EQ(clusters[2], clusters[0]);
    ASSERT_NE(clusters[2], clusters[1]);
}

TEST(kmeans, TwoClustersRandom) {
    std::default_random_engine generator;
    std::uniform_real_distribution<uint32_t> noise(0, 1);

    arma::mat points(100, 3);
    for (uint32_t i = 0; i < 50; ++i) {
        points.row(i) = { 1. + noise(generator), 2. + noise(generator), 3. + noise(generator) };
    }
    for (uint32_t i = 50; i < 100; ++i) {
        points.row(i) = { 1. + noise(generator), 2. + noise(generator), 3. + noise(generator) };
    }
    std::vector<uint32_t> clusters = kmeans(points, 2, 100);
    for (uint32_t i = 1; i < 50; ++i) {
        ASSERT_EQ(clusters[0], clusters[i]);
    }
    for (uint32_t i = 51; i < 100; ++i) {
        ASSERT_EQ(clusters[50], clusters[i]);
    }
    ASSERT_NE(clusters[0], clusters[50]);
}

TEST(kmeans, TwoClustersSpectralClustering) {
    ifstream f("data/kmeans.csv");
    std::string line;
    arma::mat points(2500, 4);
    uint32_t l = 0;
    while (std::getline(f, line)) {
        std::vector<double> point = double_split(line, ' ');
        for(uint32_t c = 0; c < 4/*point.size()*/; ++c) {
            points(l,c) = point[c];
        }
        l++;
    }

    std::vector<uint32_t> clusters = kmeans(points, 2, 100);

    for (uint32_t i = 1; i < 500; ++i) {
        ASSERT_EQ(clusters[0], clusters[i]) << "Point " << i;
    }
    for (uint32_t i = 501; i < 2500; ++i) {
        ASSERT_EQ(clusters[500], clusters[i]) << "Point " << i;
    }
    ASSERT_NE(clusters[0], clusters[500]);
}


} // namespace
