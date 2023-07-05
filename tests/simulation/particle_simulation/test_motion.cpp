// author: @Jan-Piotraschke
// date: 2023-07-05
// license: Apache License 2.0
// version: 0.1.0

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <particle_simulation/motion.h>


// Assuming that the function is declared in a namespace called 'your_namespace'
void transform_into_symmetric_matrix(Eigen::MatrixXd &A);

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

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
