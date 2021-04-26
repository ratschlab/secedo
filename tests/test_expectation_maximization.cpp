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
    CellData cell_data_a = { 1000, 0, 1 }; // cell id '0' has a read of 'C'
    CellData cell_data_b = { 1001, 1, 1 }; // cell id '1' also has a read of 'C'
    PosData one_pos = { 1234, { cell_data_a, cell_data_b } };
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
    CellData cell_data_a = { 1000, 0, 1 }; // cell id '0' has a read of 'C'
    CellData cell_data_b = { 1001, 1, 2 }; // cell id '1' has a read of 'G'
    PosData one_pos = { 1234, { cell_data_a, cell_data_b } };
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
    CellData cell_data_a = { 1000, 0, 1 }; // cell id '0' has a read of 'C'
    CellData cell_data_b = { 1001, 1, 1 }; // cell id '1' also has a read of 'C'
    CellData cell_data_c = { 1002, 2, 2 }; // cell id '2' has a read of 'G'
    CellData cell_data_d = { 1003, 3, 2 }; // cell id '3' also has a read of 'G'
    PosData one_pos = { 1234, { cell_data_a, cell_data_b, cell_data_c, cell_data_d } };
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
    CellData cell_data_a = { 1000, 0, 2 }; // cell id '0' has a read of 'C'
    CellData cell_data_b = { 1001, 1, 1 }; // cell id '1' also has a read of 'C'
    CellData cell_data_c = { 1002, 2, 2 }; // cell id '2' has a read of 'G'
    CellData cell_data_d = { 1003, 3, 2 }; // cell id '3' also has a read of 'G'
    PosData one_pos = { 1234, { cell_data_a, cell_data_b, cell_data_c, cell_data_d } };
    const std::vector<std::vector<PosData>> pos_data = { { one_pos, one_pos, one_pos } };
    std::vector<double> prob_cluster_b = { 0.9, 0.9, 0.03, 0.1 };
    expectation_maximization(pos_data, { 0, 1, 2, 3 }, 1, theta, &prob_cluster_b);
    ASSERT_THAT(prob_cluster_b,
                ElementsAre(DoubleNear(0, 1e-3), DoubleNear(1, 1e-3), DoubleNear(0, 1e-3),
                            DoubleNear(0, 1e-3)));
}

} // namespace
