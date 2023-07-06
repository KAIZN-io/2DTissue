// author: @Jan-Piotraschke
// date: 2023-07-05
// license: Apache License 2.0
// version: 0.1.0

#include <gtest/gtest.h>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <particle_simulation/motion.h>


/**
 * @brief Test the function transform_into_symmetric_matrix
*/
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


/**
 * @brief Test the function get_dist_vect
*/
TEST(GetDistVectTest, BasicTest) {
    // 3x2 matrix
    Eigen::Matrix<double, Eigen::Dynamic, 2> r(3, 2);
    r << 1, 2,
         3, 4,
         5, 6;

    std::vector<Eigen::MatrixXd> dist_vect = get_dist_vect(r);

    Eigen::MatrixXd expected_diff_x(3, 3);
    expected_diff_x <<  0, -2, -4,
                        2,  0, -2,
                        4,  2,  0;

    Eigen::MatrixXd expected_diff_y(3, 3);
    expected_diff_y <<  0, -2, -4,
                        2,  0, -2,
                        4,  2,  0;

    ASSERT_TRUE(dist_vect[0].isApprox(expected_diff_x));
    ASSERT_TRUE(dist_vect[1].isApprox(expected_diff_y));
}

TEST(GetDistVectTest, ZeroMatrixTest) {
    // 3x2 matrix
    Eigen::Matrix<double, Eigen::Dynamic, 2> r(3, 2);
    r << 0, 0,
         0, 0,
         0, 0;

    std::vector<Eigen::MatrixXd> dist_vect = get_dist_vect(r);

    Eigen::MatrixXd expected_diff(3, 3);
    expected_diff << 0, 0, 0,
                     0, 0, 0,
                     0, 0, 0;

    ASSERT_TRUE(dist_vect[0].isApprox(expected_diff));
    ASSERT_TRUE(dist_vect[1].isApprox(expected_diff));
}

TEST(GetDistVectTest, OneDimensionTest) {
    // 1x2 matrix
    Eigen::Matrix<double, Eigen::Dynamic, 2> r(1, 2);
    r << 1, 2;

    std::vector<Eigen::MatrixXd> dist_vect = get_dist_vect(r);

    Eigen::MatrixXd expected_diff(1, 1);
    expected_diff << 0;

    ASSERT_TRUE(dist_vect[0].isApprox(expected_diff));
    ASSERT_TRUE(dist_vect[1].isApprox(expected_diff));
}

TEST(GetDistVectTest, HandlesSquareMatrixCorrectly) {
    Eigen::Matrix<double, Eigen::Dynamic, 2> r(2, 2);
    r << 1.0, 2.0, 3.0, 4.0;

    std::vector<Eigen::MatrixXd> dist_vect = get_dist_vect(r);

    Eigen::MatrixXd expected_diff_x(2, 2);
    expected_diff_x << 0.0, -2.0, 2.0, 0.0;

    Eigen::MatrixXd expected_diff_y(2, 2);
    expected_diff_y << 0.0, -2.0, 2.0, 0.0;

    double tolerance = 1e-5;
    ASSERT_TRUE(dist_vect[0].isApprox(expected_diff_x, tolerance));
    ASSERT_TRUE(dist_vect[1].isApprox(expected_diff_y, tolerance));
}

TEST(GetDistVectTest, TwoDimensionTest) {
    // 2x2 matrix
    Eigen::Matrix<double, Eigen::Dynamic, 2> r(10, 2);
    r << 0.448453,  0.365021,
        0.252378,  0.477139,
        0.0309307,  0.327166,
        0.903785,  0.160117,
        0.257268,  0.436529,
        0.289008,   0.49713,
        0.844151,  0.209798,
        0.268783,  0.352196,
        0.968116,  0.673458,
        0.188375,  0.567491;

    std::vector<Eigen::MatrixXd> dist_vect = get_dist_vect(r);

    Eigen::MatrixXd expected_diff_x(10, 10);
    expected_diff_x << 0,    0.196075,    0.417522,   -0.455332,    0.191185,    0.159445,   -0.395698,     0.17967,   -0.519663,    0.260078,
                -0.196075,           0,    0.221447,   -0.651407, -0.00488968,  -0.0366297,   -0.591773,  -0.0164047,   -0.715738,   0.0640027,
                -0.417522,   -0.221447,           0,   -0.872854,   -0.226337,   -0.258077,    -0.81322,   -0.237852,   -0.937185,   -0.157445,
                0.455332,    0.651407,    0.872854,           0,    0.646517,    0.614777,   0.0596337,    0.635002,   -0.064331,    0.715409,
                -0.191185,  0.00488968,    0.226337,   -0.646517,           0,    -0.03174,   -0.586883,   -0.011515,   -0.710848,   0.0688923,
                -0.159445,   0.0366297,    0.258077,   -0.614777,     0.03174,           0,   -0.555143,    0.020225,   -0.679108,    0.100632,
                0.395698,    0.591773,     0.81322,  -0.0596337,    0.586883,    0.555143,           0,    0.575368,   -0.123965,    0.655776,
                -0.17967,   0.0164047,    0.237852,   -0.635002,    0.011515,   -0.020225,   -0.575368,           0,   -0.699333,   0.0804073,
                0.519663,    0.715738,    0.937185,    0.064331,    0.710848,    0.679108,    0.123965,    0.699333,           0,     0.77974,
                -0.260078,  -0.0640027,    0.157445,   -0.715409,  -0.0688923,   -0.100632,   -0.655776,  -0.0804073,    -0.77974,           0;

    Eigen::MatrixXd expected_diff_y(10, 10);
    expected_diff_y << 0,  -0.112118,  0.0378543,   0.204904, -0.0715087,  -0.132109,   0.155223,  0.0128243,  -0.308438,  -0.202471,
                0.112118,          0,   0.149973,   0.317022,  0.0406097,  -0.019991,   0.267341,   0.124943,  -0.196319, -0.0903523,
                -0.0378543,  -0.149973,          0,   0.167049,  -0.109363,  -0.169964,   0.117369,   -0.02503,  -0.346292,  -0.240325,
                -0.204904,  -0.317022,  -0.167049,          0,  -0.276412,  -0.337013, -0.0496807,  -0.192079,  -0.513341,  -0.407374,
                0.0715087, -0.0406097,   0.109363,   0.276412,          0, -0.0606007,   0.226732,   0.084333,  -0.236929,  -0.130962,
                0.132109,   0.019991,   0.169964,   0.337013,  0.0606007,          0,   0.287332,   0.144934,  -0.176328, -0.0703613,
                -0.155223,  -0.267341,  -0.117369,  0.0496807,  -0.226732,  -0.287332,          0,  -0.142399,  -0.463661,  -0.357694,
                -0.0128243,  -0.124943,    0.02503,   0.192079,  -0.084333,  -0.144934,   0.142399,          0,  -0.321262,  -0.215295,
                0.308438,   0.196319,   0.346292,   0.513341,   0.236929,   0.176328,   0.463661,   0.321262,          0,   0.105967,
                0.202471,  0.0903523,   0.240325,   0.407374,   0.130962,  0.0703613,   0.357694,   0.215295,  -0.105967,          0;

    double tolerance = 1e-5;
    ASSERT_TRUE(dist_vect[0].isApprox(expected_diff_x, tolerance));
    ASSERT_TRUE(dist_vect[1].isApprox(expected_diff_y, tolerance));
}


/**
 * @brief Test the function mean_unit_circle_vector_angle_degrees
*/
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
