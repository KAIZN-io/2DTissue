// author: @Jan-Piotraschke
// date: 2023-07-18
// license: Apache License 2.0
// version: 0.2.0

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <iostream>

#include <Simulator.h>

// Define the test fixture class
class SimulatorTest : public ::testing::Test {
protected:
    SimulatorTest()
        : v0(0), k(10), σ(1.4166666666666667), μ(0), r_adh(1), k_adh(0.75), step_size(0),
          r_UV(Eigen::Matrix<double, Eigen::Dynamic, 2>()),
          r_dot(Eigen::Matrix<double, Eigen::Dynamic, 2>()),
          n(Eigen::VectorXd()),
          vertices_3D_active(std::vector<int>()),
          distance_matrix(Eigen::MatrixXd()),
          dist_length(Eigen::MatrixXd()),
          sim(r_UV, r_dot, n, vertices_3D_active, distance_matrix, dist_length, v0, k, σ, μ, r_adh, k_adh, step_size)
    {}

    double v0, k, σ, μ, r_adh, k_adh, step_size;
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV, r_dot;
    Eigen::VectorXd n;
    std::vector<int> vertices_3D_active;
    Eigen::MatrixXd distance_matrix, dist_length;
    Simulator sim;
};


// Utilize the test fixture in your tests
TEST_F(SimulatorTest, ThrowsWhenInputIsEmpty) {
    std::cout << "Testing mean_unit_circle_vector_angle_degrees" << std::endl;
    std::vector<double> empty;
    EXPECT_THROW(sim.mean_unit_circle_vector_angle_degrees(empty), std::invalid_argument);
}

TEST_F(SimulatorTest, CorrectlyCalculatesMeanAngle) {
    std::vector<double> angles {0, 180, 90};
    double expected_mean_angle = 90.0;
    double actual_mean_angle = sim.mean_unit_circle_vector_angle_degrees(angles);

    EXPECT_NEAR(expected_mean_angle, actual_mean_angle, 1e-5); // 1e-5 is the allowed error
}

TEST_F(SimulatorTest, CorrectlyHandlesNegativeAngles) {
    std::vector<double> angles {-45, -90, -135, -180, -225, -270, -315};
    double expected_mean_angle = 180.0;
    double actual_mean_angle = sim.mean_unit_circle_vector_angle_degrees(angles);

    EXPECT_NEAR(expected_mean_angle, actual_mean_angle, 1e-5); // 1e-5 is the allowed error
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
