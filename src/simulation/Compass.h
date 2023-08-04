// Compass.h

#pragma once

#include <map>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <cstdint>
#include <gtest/gtest_prod.h>

class Compass {
public:
    Compass(const Eigen::Vector2d& original_pole)
        : original_pole_(original_pole) {}

    Eigen::VectorXd calculate_n_pole(const Eigen::Matrix<double, Eigen::Dynamic, 2>& positions, const Eigen::VectorXd& orientations) const {
        Eigen::VectorXd n_pole_(positions.rows());

        for(int i = 0; i < positions.rows(); i++) {
            Eigen::Vector2d position = positions.row(i);
            double angle = vectorAngle(position, original_pole_); // use original_pole
            n_pole_(i) = relativeAngle(orientations(i), angle);
        }

        return n_pole_;
    }

    Eigen::VectorXd calculate_n(const Eigen::Matrix<double, Eigen::Dynamic, 2>& positions, const Eigen::Vector2d& pole, const Eigen::VectorXd& n_pole_) const {
        Eigen::VectorXd n_values(positions.rows());

        // Calculate the angle between the vertical 0-degree line and the pole vector for each position
        Eigen::VectorXd deltas = calculate_delta(positions, pole);

        for(int i = 0; i < positions.rows(); i++) {
            n_values(i) = fmod(std::abs(deltas(i) - n_pole_(i)) + 360, 360);
        }

        return n_values;
    }

    Eigen::VectorXd assign_n_pole_orientation(const Eigen::Matrix<double, Eigen::Dynamic, 2>& newPositions, const Eigen::VectorXd& n_pole_, Eigen::Vector2d virtual_pole_) const {
        Eigen::VectorXd newOrientations(newPositions.rows());

        for(int i = 0; i < newPositions.rows(); i++) {
            Eigen::Vector2d newPosition = newPositions.row(i);
            double angle = vectorAngle(newPosition, virtual_pole_); // use virtual_pole
            newOrientations(i) = fmod(angle - n_pole_(i) + 360, 360);
        }

        return newOrientations;
    }

private:
    Eigen::Vector2d original_pole_;

    double vectorAngle(const Eigen::Vector2d& position, const Eigen::Vector2d& pole) const {
        double dx = position.x() - pole.x();
        double dy = position.y() - pole.y();
        double rad = atan2(dy, dx);
        double deg = rad * (180 / M_PI);
        return fmod(deg + 360, 360);
    }

    double relativeAngle(double orientation, double vectorAngle) const {
        double relative = vectorAngle - orientation;
        return fmod(relative + 360, 360);
    }

    Eigen::VectorXd calculate_delta(const Eigen::Matrix<double, Eigen::Dynamic, 2>& start_points, const Eigen::Vector2d& end_point) const {
        Eigen::VectorXd angles(start_points.rows());
        Eigen::Vector2d upVector(0, 1); // A vector pointing directly up

        for(int i = 0; i < start_points.rows(); i++) {
            Eigen::Vector2d position = start_points.row(i);
            Eigen::Vector2d vectorToPole = end_point - position;

            double dotProduct = upVector.dot(vectorToPole);
            double crossProduct = upVector.x() * vectorToPole.y() - upVector.y() * vectorToPole.x();

            double angle_rad = std::atan2(crossProduct, dotProduct);
            double angle_deg = angle_rad * 180 / M_PI;

            // Since we want the clockwise angle, subtract the angle from 360 if it's positive
            if (angle_deg > 0) {
                angle_deg = 360 - angle_deg;
            } else {
                angle_deg = -angle_deg;
            }

            angles(i) = angle_deg;
        }

        return angles;
    }

FRIEND_TEST(CompassTest, TestVectorAngle);
FRIEND_TEST(CompassTest, TestCalculateDelta);

};
