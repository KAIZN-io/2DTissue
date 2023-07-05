// author: @Jan-Piotraschke
// date: 2023-07-05
// license: Apache License 2.0
// version: 0.1.0

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <particle_simulation/motion.h>


TEST(SymmetricMatrixTest, BasicTest) {
    Eigen::MatrixXd mat(3, 3);
    mat << 1, 2, 3,
           4, 5, 6,
           7, 8, 9;

    transform_into_symmetric_matrix(mat);

    // Expected symmetric matrix
    Eigen::MatrixXd expected_mat(3, 3);
    expected_mat << 1, 2, 3,
                    2, 5, 6,
                    3, 6, 9;

    ASSERT_TRUE(mat.isApprox(expected_mat));
}

TEST(SymmetricMatrixTest, ZeroTest) {
    Eigen::MatrixXd mat(3, 3);
    mat << 1, 0, 3,
           4, 5, 6,
           7, 8, 0;

    transform_into_symmetric_matrix(mat);

    // Expected symmetric matrix
    Eigen::MatrixXd expected_mat(3, 3);
    expected_mat << 1, 0, 3,
                    0, 5, 6,
                    3, 6, 0;

    ASSERT_TRUE(mat.isApprox(expected_mat));
}

TEST(SymmetricMatrixTest, AllZerosTest) {
    Eigen::MatrixXd mat(3, 3);
    mat << 0, 0, 0,
           0, 0, 0,
           0, 0, 0;

    transform_into_symmetric_matrix(mat);

    // Expected symmetric matrix
    Eigen::MatrixXd expected_mat(3, 3);
    expected_mat << 0, 0, 0,
                    0, 0, 0,
                    0, 0, 0;

    ASSERT_TRUE(mat.isApprox(expected_mat));
}


TEST(UnitCircleVectorTest, ThrowsWhenInputIsEmpty) {
    std::vector<double> empty;
    EXPECT_THROW(mean_unit_circle_vector_angle_degrees(empty), std::invalid_argument);
}

TEST(UnitCircleVectorTest, CorrectlyCalculatesMeanAngle) {
    std::vector<double> angles {0, 180, 90};
    double expected_mean_angle = 90.0;
    double actual_mean_angle = mean_unit_circle_vector_angle_degrees(angles);

    EXPECT_NEAR(expected_mean_angle, actual_mean_angle, 1e-5); // 1e-5 is the allowed error
}

TEST(UnitCircleVectorTest, CorrectlyHandlesNegativeAngles) {
    std::vector<double> angles {-45, -90, -135, -180, -225, -270, -315};
    double expected_mean_angle = 180.0;
    double actual_mean_angle = mean_unit_circle_vector_angle_degrees(angles);

    EXPECT_NEAR(expected_mean_angle, actual_mean_angle, 1e-5); // 1e-5 is the allowed error
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
