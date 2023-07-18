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


// Test the repulsive_adhesion_motion function
TEST_F(SimulatorTest, RepulsiveAdhesionTest1) {
    Eigen::Vector2d dist_v(0.0203217, 0.010791);
    double dist = 0.953489;
    Eigen::Vector2d expected_result(-0.141406, -0.075088);
    Eigen::Vector2d result = sim.repulsive_adhesion_motion(k, σ, dist, r_adh, k_adh, dist_v);

    ASSERT_NEAR(result[0], expected_result[0], 1e-5);
    ASSERT_NEAR(result[1], expected_result[1], 1e-5);
}

TEST_F(SimulatorTest, RepulsiveAdhesionTest2) {
    Eigen::Vector2d dist_v(-0.0203217, -0.010791);
    double dist = 0.953489;
    Eigen::Vector2d expected_result(0.141406, 0.075088);
    Eigen::Vector2d result = sim.repulsive_adhesion_motion(k, σ, dist, r_adh, k_adh, dist_v);

    ASSERT_NEAR(result[0], expected_result[0], 1e-5);
    ASSERT_NEAR(result[1], expected_result[1], 1e-5);
}


// Test the mean_unit_circle_vector_angle_degrees function
TEST_F(SimulatorTest, ThrowsWhenInputIsEmpty) {
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


/**
 * @brief Test the function transform_into_symmetric_matrix
*/
TEST_F(SimulatorTest, SymmetricMatrixTestBasicTest) {
    Eigen::MatrixXd mat(3, 3);
    mat << 1, 2, 3,
           4, 5, 6,
           7, 8, 9;

    sim.transform_into_symmetric_matrix(mat);

    // Expected symmetric matrix
    Eigen::MatrixXd expected_mat(3, 3);
    expected_mat << 1, 2, 3,
                    2, 5, 6,
                    3, 6, 9;

    ASSERT_TRUE(mat.isApprox(expected_mat));
}

TEST_F(SimulatorTest, SymmetricMatrixTestZeroTest) {
    Eigen::MatrixXd mat(3, 3);
    mat << 1, 0, 3,
           4, 5, 6,
           7, 8, 0;

    sim.transform_into_symmetric_matrix(mat);

    // Expected symmetric matrix
    Eigen::MatrixXd expected_mat(3, 3);
    expected_mat << 1, 0, 3,
                    0, 5, 6,
                    3, 6, 0;

    ASSERT_TRUE(mat.isApprox(expected_mat));
}

TEST_F(SimulatorTest, SymmetricMatrixTestAllZerosTest) {
    Eigen::MatrixXd mat(3, 3);
    mat << 0, 0, 0,
           0, 0, 0,
           0, 0, 0;

    sim.transform_into_symmetric_matrix(mat);

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
TEST_F(SimulatorTest, GetDistVectTestBasicTest) {
    // 3x2 matrix
    Eigen::Matrix<double, Eigen::Dynamic, 2> r(3, 2);
    r << 1, 2,
         3, 4,
         5, 6;

    std::vector<Eigen::MatrixXd> dist_vect = sim.get_dist_vect(r);

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

TEST_F(SimulatorTest, GetDistVectTestZeroMatrixTest) {
    // 3x2 matrix
    Eigen::Matrix<double, Eigen::Dynamic, 2> r(3, 2);
    r << 0, 0,
         0, 0,
         0, 0;

    std::vector<Eigen::MatrixXd> dist_vect = sim.get_dist_vect(r);

    Eigen::MatrixXd expected_diff(3, 3);
    expected_diff << 0, 0, 0,
                     0, 0, 0,
                     0, 0, 0;

    ASSERT_TRUE(dist_vect[0].isApprox(expected_diff));
    ASSERT_TRUE(dist_vect[1].isApprox(expected_diff));
}

TEST_F(SimulatorTest, GetDistVectTestOneDimensionTest) {
    // 1x2 matrix
    Eigen::Matrix<double, Eigen::Dynamic, 2> r(1, 2);
    r << 1, 2;

    std::vector<Eigen::MatrixXd> dist_vect = sim.get_dist_vect(r);

    Eigen::MatrixXd expected_diff(1, 1);
    expected_diff << 0;

    ASSERT_TRUE(dist_vect[0].isApprox(expected_diff));
    ASSERT_TRUE(dist_vect[1].isApprox(expected_diff));
}

TEST_F(SimulatorTest, GetDistVectTestHandlesSquareMatrixCorrectly) {
    Eigen::Matrix<double, Eigen::Dynamic, 2> r(2, 2);
    r << 1.0, 2.0, 3.0, 4.0;

    std::vector<Eigen::MatrixXd> dist_vect = sim.get_dist_vect(r);

    Eigen::MatrixXd expected_diff_x(2, 2);
    expected_diff_x << 0.0, -2.0, 2.0, 0.0;

    Eigen::MatrixXd expected_diff_y(2, 2);
    expected_diff_y << 0.0, -2.0, 2.0, 0.0;

    double tolerance = 1e-5;
    ASSERT_TRUE(dist_vect[0].isApprox(expected_diff_x, tolerance));
    ASSERT_TRUE(dist_vect[1].isApprox(expected_diff_y, tolerance));
}

TEST_F(SimulatorTest, GetDistVectTestTenDimensionTest) {
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

    std::vector<Eigen::MatrixXd> dist_vect = sim.get_dist_vect(r);

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
 * @brief Test the function calculate_average_n_within_distance
*/
void CompareMatrices(const Eigen::MatrixXd& expected, const Eigen::MatrixXd& actual, double tolerance) {
    ASSERT_EQ(expected.rows(), actual.rows());
    ASSERT_EQ(expected.cols(), actual.cols());

    for(int i = 0; i < expected.rows(); ++i){
        for(int j = 0; j < expected.cols(); ++j){
            EXPECT_NEAR(expected(i, j), actual(i, j), tolerance);
        }
    }
}

TEST_F(SimulatorTest, AverageNWithinDistanceTest1){
    double σ = 1.4166666666666667;

    Eigen::MatrixXd diff_x(10, 10);
    diff_x << 0,  -0.041037,    0.22485,   0.550749,   0.614337,    0.27906,   0.312221,  -0.183676,   0.422344,   0.344173,
        0.041037,          0,   0.265887,   0.591786,   0.655374,   0.320097,   0.353258,  -0.142639,   0.463381,    0.38521,
        -0.22485,  -0.265887,          0,   0.325899,   0.389487,  0.0542097,  0.0873707,  -0.408526,   0.197494,   0.119323,
        -0.550749,  -0.591786,  -0.325899,          0,  0.0635878,   -0.27169,  -0.238529,  -0.734425,  -0.128405,  -0.206576,
        -0.614337,  -0.655374,  -0.389487, -0.0635878,          0,  -0.335277,  -0.302116,  -0.798013,  -0.191993,  -0.270164,
        -0.27906,  -0.320097, -0.0542097,    0.27169,   0.335277,          0,   0.033161,  -0.462735,   0.143284,  0.0651137,
        -0.312221,  -0.353258, -0.0873707,   0.238529,   0.302116,  -0.033161,          0,  -0.495896,   0.110123,  0.0319527,
        0.183676,   0.142639,   0.408526,   0.734425,   0.798013,   0.462735,   0.495896,          0,    0.60602,   0.527849,
        -0.422344,  -0.463381,  -0.197494,   0.128405,   0.191993,  -0.143284,  -0.110123,   -0.60602,          0, -0.0781707,
        -0.344173,   -0.38521,  -0.119323,   0.206576,   0.270164, -0.0651137, -0.0319527,  -0.527849,  0.0781707,          0;

    Eigen::MatrixXd diff_y(10, 10);
    diff_y << 0,  0.0483417,   0.289696,   -0.34716,   0.175601,  0.0146273,  -0.164702,  -0.198517,  -0.109245,  0.0868453,
        -0.0483417,          0,   0.241355,  -0.395502,   0.127259, -0.0337143,  -0.213044,  -0.246859,  -0.157587,  0.0385037,
        -0.289696,  -0.241355,          0,  -0.636857,  -0.114096,  -0.275069,  -0.454398,  -0.488213,  -0.398941,  -0.202851,
        0.34716,   0.395502,   0.636857,          0,   0.522761,   0.361788,   0.182458,   0.148643,   0.237915,   0.434006,
        -0.175601,  -0.127259,   0.114096,  -0.522761,          0,  -0.160973,  -0.340303,  -0.374118,  -0.284846, -0.0887553,
        -0.0146273,  0.0337143,   0.275069,  -0.361788,   0.160973,          0,  -0.179329,  -0.213144,  -0.123872,   0.072218,
        0.164702,   0.213044,   0.454398,  -0.182458,   0.340303,   0.179329,          0,  -0.033815,   0.055457,   0.251547,
        0.198517,   0.246859,   0.488213,  -0.148643,   0.374118,   0.213144,   0.033815,          0,   0.089272,   0.285362,
        0.109245,   0.157587,   0.398941,  -0.237915,   0.284846,   0.123872,  -0.055457,  -0.089272,          0,    0.19609,
        -0.0868453, -0.0385037,   0.202851,  -0.434006,  0.0887553,  -0.072218,  -0.251547,  -0.285362,   -0.19609,          0;

    std::vector<Eigen::MatrixXd> dist_vect;
    dist_vect.push_back(diff_x);
    dist_vect.push_back(diff_y);

    Eigen::MatrixXd dist_length(10, 10);
    dist_length << 0, 1.94061, 5.60903, 8.28046, 8.47736, 11.0131, 14.2291,  6.0693, 12.7292,  10.761,
            1.94061,       0, 4.51284, 6.31737, 7.23179, 11.6579,  12.934,  5.8962,  10.894, 10.6601,
            5.60903, 4.51284,       0, 5.55914, 2.81003, 9.00323, 10.8524, 10.1511, 7.80554, 6.77901,
            8.28046, 6.31737, 5.55914,       0, 5.75308, 14.5522,  6.5117, 7.98363, 5.19999, 11.8869,
            8.47736, 7.23179, 2.81003, 5.75308,       0, 9.30527, 9.01998, 12.5087, 5.71498, 6.17054,
            11.0131, 11.6579, 9.00323, 14.5522, 9.30527,       0, 9.45922, 14.2346, 12.4246, 3.59292,
            14.2291,  12.934, 10.8524,  6.5117, 9.01998, 9.45922,       0, 10.0873, 3.28317, 12.8653,
            6.0693,  5.8962, 10.1511, 7.98363, 12.5087, 14.2346, 10.0873,       0,  11.829, 16.5712,
            12.7292,  10.894, 7.80554, 5.19999, 5.71498, 12.4246, 3.28317,  11.829,       0, 11.2371,
            10.761, 10.6601, 6.77901, 11.8869, 6.17054, 3.59292, 12.8653, 16.5712, 11.2371,       0;

    Eigen::VectorXd n(10);
    n << 168, 154, 290, 83, 110, 46, 48, 144, 227, 48;

    sim.calculate_average_n_within_distance(dist_vect, dist_length, n, σ);

    Eigen::VectorXd expected_n(10);
    expected_n << 161, 161, 191.31, 83, 191.31, 46, 48, 144, 227, 48;

    CompareMatrices(expected_n, n, 2);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
