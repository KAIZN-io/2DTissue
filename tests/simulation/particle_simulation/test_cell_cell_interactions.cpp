// author: @Jan-Piotraschke
// date: 2023-07-05
// license: Apache License 2.0
// version: 0.1.0

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <particle_simulation/cell_cell_interactions.h>


class RepulsiveAdhesionTest : public ::testing::Test {
protected:
    const double k = 10;
    const double σ = 1.4166666666666667;
    const double r_adh = 1;
    const double k_adh = 0.75;
};

TEST_F(RepulsiveAdhesionTest, Test1) {
    Eigen::Vector2d dist_v(0.0203217, 0.010791);
    double dist = 0.953489;
    Eigen::Vector2d expected_result(-0.141406, -0.075088);
    Eigen::Vector2d result = repulsive_adhesion_motion(k, σ, dist, r_adh, k_adh, dist_v);

    ASSERT_NEAR(result[0], expected_result[0], 1e-5);
    ASSERT_NEAR(result[1], expected_result[1], 1e-5);
}

TEST_F(RepulsiveAdhesionTest, Test2) {
    Eigen::Vector2d dist_v(-0.0203217, -0.010791);
    double dist = 0.953489;
    Eigen::Vector2d expected_result(0.141406, 0.075088);
    Eigen::Vector2d result = repulsive_adhesion_motion(k, σ, dist, r_adh, k_adh, dist_v);

    ASSERT_NEAR(result[0], expected_result[0], 1e-5);
    ASSERT_NEAR(result[1], expected_result[1], 1e-5);
}