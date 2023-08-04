// author: @Jan-Piotraschke
// date: 2023-08-04
// license: Apache License 2.0
// version: 0.1.0

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

#include <Compass.h>

class CompassTest : public ::testing::Test {
protected:
    Eigen::Vector2d original_pole{0.1021, 0.0347795};
    Compass compass{original_pole};
};

TEST_F(CompassTest, TestCalculateNPole) {
    Eigen::Matrix<double, Eigen::Dynamic, 2> positions(3, 2);
    positions << 0.3021, 0.0347795,
                 0.3, 0.4,
                 0.1021, 0.3;
    Eigen::VectorXd orientations(3);
    orientations << 0.0, 90.0, 270.0;

    Eigen::VectorXd result = compass.calculate_n_pole(positions, orientations);

    double tolerance = 3;
    EXPECT_DOUBLE_EQ(result(0), 270.0);
    EXPECT_NEAR(result(1), 120, tolerance);
    EXPECT_DOUBLE_EQ(result(2), 270);
}

TEST_F(CompassTest, TestCalculateDelta) {
    Eigen::Matrix<double, Eigen::Dynamic, 2> start_points;
    start_points.resize(6, Eigen::NoChange);
    start_points << 0.5, 0.5,
                    0.2, 0.4,
                    0.2, 0.8,
                    0.1, 0.5,
                    0.3, 0.4,
                    0.1, 0.4;
    Eigen::Vector2d end_point(0.2, 0.5);

    Eigen::VectorXd result = compass.calculate_delta(start_points, end_point);
    EXPECT_DOUBLE_EQ(result(0), 270.0);
    EXPECT_DOUBLE_EQ(result(1), 0.0);
    EXPECT_DOUBLE_EQ(result(2), 180.0);
    EXPECT_DOUBLE_EQ(result(3), 90.0);
    EXPECT_DOUBLE_EQ(result(4), 315.0);
    EXPECT_DOUBLE_EQ(result(5), 45.0);
}

TEST_F(CompassTest, TestCalculateN) {
    Eigen::Matrix<double, Eigen::Dynamic, 2> positions(4, 2);
    positions << 0.3021, 0.0347795,
                 0.3, 0.4,
                 0.1021, 0.3,
                 0.83421, 0.9347795;
    Eigen::VectorXd n_start(4);
    n_start << 10.0, 90.0, 270.0, 133.1;

    Eigen::VectorXd n_pole = compass.calculate_n_pole(positions, n_start);
    Eigen::VectorXd result = compass.assign_n(positions, original_pole, n_pole);

    double tolerance = 3;
    EXPECT_NEAR(result(0), 10.0, tolerance);
    EXPECT_NEAR(result(1), 90.0, tolerance);
    EXPECT_NEAR(result(2), 270.0, tolerance);
    EXPECT_NEAR(result(3), 133.1, tolerance);
}
