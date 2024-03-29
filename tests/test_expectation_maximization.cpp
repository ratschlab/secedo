#include "expectation_maximization.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {
using namespace testing;

constexpr double theta = 1e-3;

/**
 * One cell and no reads - the EM will converge immediately and not change the input probabilities
 */
TEST(EM, OneCell) {
    const std::vector<std::vector<PosData>> pos_data = { {} };
    std::vector<double> prob_cluster_b = { 1.0 };
    expectation_maximization(pos_data, {}, 1, theta, &prob_cluster_b);
    ASSERT_EQ(1, prob_cluster_b.size());
    ASSERT_EQ(1.0, prob_cluster_b[0]);
}

/**
 * We have 2 cells, each with 3 identical reads - the probabilities of the cells belonging to either
 * cluster should converge to equal values. This is counter-intuitive: the cells won't belong to the
 * same cluster with probability 1.0. Instead, since the cells are identical, all we can say is that
 * they will have identical probabilities of belonging to either cluster.
 */
TEST(EM, TwoCellsSame) {
    std::vector<uint32_t> read_ids = { 1000, 1001 };
    std::vector<uint16_t> cell_ids_bases = { 0 << 2 | 1, 1 << 2 | 1 };
    PosData one_pos = { 1234, read_ids, cell_ids_bases };
    const std::vector<std::vector<PosData>> pos_data = { { one_pos, one_pos, one_pos } };
    std::vector<double> prob_cluster_b = { 0.3, 0.4 };
    expectation_maximization(pos_data, { 0, 1 }, 1, theta, &prob_cluster_b);
    ASSERT_NEAR(0, prob_cluster_b[1] - prob_cluster_b[0], 1e-3);
}

/**
 * We have 2 cells, each with 3 different reads. Even though the initial probabilities put both
 * cells with high likelihood in the same cluster, in the end the probabilities should converge
 * towards different clusters.
 */
TEST(EM, TwoCellsDifferent) {
    std::vector<uint32_t> read_ids = { 1000, 1001 };
    std::vector<uint16_t> cell_ids_bases = { 0 << 2 | 1, 1 << 2 | 2 };
    PosData one_pos = { 1234, read_ids, cell_ids_bases };
    const std::vector<std::vector<PosData>> pos_data = { { one_pos, one_pos, one_pos } };
    std::vector<double> prob_cluster_b = { 0.01, 0.02 };
    expectation_maximization(pos_data, { 0, 1 }, 1, theta, &prob_cluster_b);
    ASSERT_NEAR(std::abs(prob_cluster_b[0] - prob_cluster_b[1]), 1.0, 1e-3);
}

/**
 * We have 4 cells, each with 3 different reads. The first 2 cells and the last 2 cells have
 * identical reads. Even though the initial probabilities put cells 0 and 3 with high likelihood in
 * the same cluster, in the end the probabilities should converge towards the correct clusters.
 */
TEST(EM, FourCellsTwoGroups22) {
    std::vector<uint32_t> read_ids = { 1000, 1001, 1002, 1003 };
    std::vector<uint16_t> cell_ids_bases = { 0 << 2 | 1, 1 << 2 | 1, 2 << 2 | 2, 3 << 2 | 2   };
    PosData one_pos = { 1234, read_ids, cell_ids_bases };
    const std::vector<std::vector<PosData>> pos_data = { { one_pos, one_pos, one_pos } };
    std::vector<double> prob_cluster_b = { 0.9, 0.02, 0.03, 0.9 };
    expectation_maximization(pos_data, { 0, 1, 2, 3 }, 1, theta, &prob_cluster_b);
    ASSERT_THAT(prob_cluster_b,
                ElementsAre(DoubleNear(0, 1e-3), DoubleNear(0, 1e-3), DoubleNear(1, 1e-3),
                            DoubleNear(1, 1e-3)));
}

/**
 * We have 4 cells, each with 3 different reads. Cells 0, 2 and 3 have identical reads, while Cell 1
 * is completely different Even though the initial probabilities put cells 0 and 1 with high
 * likelihood in the same cluster, in the end the probabilities should converge towards the correct
 * clusters.
 */
TEST(EM, FourCellsTwoGroups31) {
    std::vector<uint32_t> read_ids = { 1000, 1001, 1002, 1003 };
    std::vector<uint16_t> cell_ids_bases = { 0 << 2 | 2, 1 << 2 | 1, 2 << 2 | 2, 3 << 2 | 2   };
    PosData one_pos = { 1234, read_ids, cell_ids_bases };
    const std::vector<std::vector<PosData>> pos_data = { { one_pos, one_pos, one_pos } };
    std::vector<double> prob_cluster_b = { 0.9, 0.9, 0.03, 0.1 };
    expectation_maximization(pos_data, { 0, 1, 2, 3 }, 1, theta, &prob_cluster_b);
    ASSERT_THAT(prob_cluster_b,
                ElementsAre(DoubleNear(0, 1e-3), DoubleNear(1, 1e-3), DoubleNear(0, 1e-3),
                            DoubleNear(0, 1e-3)));
}

} // namespace
