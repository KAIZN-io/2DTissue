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
    Eigen::Matrix<double, Eigen::Dynamic, 2> positions(2, 2);
    positions << 0.5, 0.5,
                 0.7, 0.7;
    Eigen::VectorXd orientations(2);
    orientations << 0.0, 90.0;

    Eigen::VectorXd result = compass.calculate_n_pole(positions, orientations);
    EXPECT_DOUBLE_EQ(result(0), 45.0);
    EXPECT_DOUBLE_EQ(result(1), 45.0);
}

TEST_F(CompassTest, TestAssignNPoleOrientation) {
    Eigen::Matrix<double, Eigen::Dynamic, 2> newPositions(2, 2);
    newPositions << 0.938641, 0.627895,
                    0.419308, 0.601002;
    Eigen::VectorXd relativeAngles(2);
    relativeAngles << 45.0, 135.0;

    Eigen::Vector2d virtual_pole(0.8, 0.8);
    Eigen::VectorXd result = compass.assign_n_pole_orientation(newPositions, relativeAngles, virtual_pole);
    EXPECT_DOUBLE_EQ(result(0), 315.0);
    EXPECT_DOUBLE_EQ(result(1), 225.0);
}

TEST_F(CompassTest, TestVectorAngle) {
    Eigen::Vector2d position(0.15, 0.15);
    Eigen::Vector2d pole(0.5, 0.5);

    double result = compass.vectorAngle(position, pole);
    EXPECT_DOUBLE_EQ(result, 135.0);
}
