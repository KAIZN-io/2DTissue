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

    Eigen::VectorXi calculate_n_pole(const Eigen::Matrix<double, Eigen::Dynamic, 2>& positions, const Eigen::VectorXi& orientations) const {
        Eigen::VectorXi n_pole_(positions.rows());

        for(int i = 0; i < positions.rows(); i++) {
            Eigen::Vector2d position = positions.row(i);
            Eigen::Vector2d orientationVector(sin(orientations(i) * M_PI / 180), cos(orientations(i) * M_PI / 180));
            double angle = calculate_angle(position, original_pole_, orientationVector);
            n_pole_(i) = fmod(angle + 360, 360);
        }

        return n_pole_;
    }

    Eigen::VectorXi calculate_delta(const Eigen::Matrix<double, Eigen::Dynamic, 2>& start_points, const Eigen::Vector2d& end_point) const {
        Eigen::VectorXi angles(start_points.rows());
        Eigen::Vector2d upVector(0, 1); // A vector pointing directly up

        for(int i = 0; i < start_points.rows(); i++) {
            Eigen::Vector2d position = start_points.row(i);
            angles(i) = calculate_angle(position, end_point, upVector);
        }

        return angles;
    }

    Eigen::VectorXi assign_n(const Eigen::Matrix<double, Eigen::Dynamic, 2>& positions, const Eigen::Vector2d& pole, const Eigen::VectorXi& n_pole_) const {
        Eigen::VectorXi n_values(positions.rows());

        // Calculate the angle between the vertical 0-degree line and the pole vector for each position
        Eigen::VectorXi deltas = calculate_delta(positions, pole);

        for(int i = 0; i < positions.rows(); i++) {
            n_values(i) = fmod(deltas(i) - n_pole_(i) + 360, 360);
        }

        return n_values;
    }

private:
    Eigen::Vector2d original_pole_;

    double calculate_angle(const Eigen::Vector2d& position, const Eigen::Vector2d& pole, const Eigen::Vector2d& orientationVector) const {
        Eigen::Vector2d vectorToPole = pole - position;

        double dotProduct = orientationVector.dot(vectorToPole);
        double crossProduct = orientationVector.x() * vectorToPole.y() - orientationVector.y() * vectorToPole.x();

        double angle_rad = std::atan2(crossProduct, dotProduct);
        double angle_deg = angle_rad * 180 / M_PI;

        // Since we want the clockwise angle, subtract the angle from 360 if it's positive
        if (angle_deg > 0) {
            angle_deg = 360 - angle_deg;
        } else {
            angle_deg = -angle_deg;
        }

        return angle_deg;
    }
};
