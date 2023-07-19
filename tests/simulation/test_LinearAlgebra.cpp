// author: @Jan-Piotraschke
// date: 2023-07-19
// license: Apache License 2.0
// version: 0.1.0

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <LinearAlgebra.h>

std::unique_ptr<LinearAlgebra> linear_algebra_ptr = std::make_unique<LinearAlgebra>();

TEST(LinearAlgebraTest, BasicAssertions) {
    // Set up the test input
    Eigen::VectorXd avg_n(6);
    avg_n << 0, 45, 90, 180, 270, 360;

    // Call the function
    Eigen::Matrix<double, Eigen::Dynamic, 2> n_vec = linear_algebra_ptr->angles_to_unit_vectors(avg_n);

    // Check the results
    ASSERT_NEAR(n_vec(0, 0), 1, 1e-9);
    ASSERT_NEAR(n_vec(0, 1), 0, 1e-9);

    ASSERT_NEAR(n_vec(1, 0), sqrt(2) / 2, 1e-9);
    ASSERT_NEAR(n_vec(1, 1), sqrt(2) / 2, 1e-9);

    ASSERT_NEAR(n_vec(2, 0), 0, 1e-9);
    ASSERT_NEAR(n_vec(2, 1), 1, 1e-9);

    ASSERT_NEAR(n_vec(3, 0), -1, 1e-9);
    ASSERT_NEAR(n_vec(3, 1), 0, 1e-9);

    ASSERT_NEAR(n_vec(4, 0), 0, 1e-9);
    ASSERT_NEAR(n_vec(4, 1), -1, 1e-9);

    ASSERT_NEAR(n_vec(5, 0), 1, 1e-9);
    ASSERT_NEAR(n_vec(5, 1), 0, 1e-9);
}
