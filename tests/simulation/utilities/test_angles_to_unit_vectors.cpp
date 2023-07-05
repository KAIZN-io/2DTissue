// author: @Jan-Piotraschke
// date: 2023-07-05
// license: Apache License 2.0
// version: 0.1.0

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <utilities/angles_to_unit_vectors.h>


TEST(AngleToUnitVectorTest, BasicAssertions) {
    // Set up the test input
    Eigen::VectorXd avg_n(6);
    avg_n << 0, 45, 90, 180, 270, 360;

    // Call the function
    Eigen::Matrix<double, Eigen::Dynamic, 2> n_vec = angles_to_unit_vectors(avg_n);

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

// Run all the tests
int main(int argc, char **argv) {
    // The function ::testing::InitGoogleTest(&argc, argv) initializes the Google Test framework.
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
