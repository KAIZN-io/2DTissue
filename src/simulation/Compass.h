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

    double calculate_n(const Eigen::Vector2d& position, const Eigen::Vector2d& pole, double n_pole_) const {

        // Gete the angle between the vertical 0-degree line and the pole vector
        double delta = calculate_delta(position, pole);

        return fmod(std::abs(delta - n_pole_) + 360, 360);
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

    double calculate_delta(const Eigen::Vector2d& position, const Eigen::Vector2d& pole) const {
        Eigen::Vector2d vectorToPole = pole - position;
        Eigen::Vector2d upVector(0, 1); // A vector pointing directly up

        double dotProduct = vectorToPole.dot(upVector);
        double magnitudeA = vectorToPole.norm();
        double magnitudeB = upVector.norm();

        double cosineAngle = dotProduct / (magnitudeA * magnitudeB);

        // Clamp the value to the range [-1, 1] to handle potential numerical inaccuracies
        cosineAngle = std::max(-1.0, std::min(1.0, cosineAngle));

        double rad = acos(cosineAngle);
        double deg = rad * (180 / M_PI);

        return deg;
    }

FRIEND_TEST(CompassTest, TestVectorAngle);

};
