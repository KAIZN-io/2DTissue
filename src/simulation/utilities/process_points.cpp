// author: @Jan-Piotraschke
// date: 2023-06-26
// license: Apache License 2.0
// version: 0.1.0

#include "utilities/process_points.h"

// function to calculate y given x
double interpolateY(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double x) {
    return pointA[1] + ((x - pointA[0]) * (pointB[1] - pointA[1])) / (pointB[0] - pointA[0]);
}

// function to calculate x given y
double interpolateX(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double y) {
    return pointA[0] + ((y - pointA[1]) * (pointB[0] - pointA[0])) / (pointB[1] - pointA[1]);
}

// Function to calculate steepness switch
int calculateSteepnessSwitch(double steepness) {
    if (steepness > 0) return -1;
    else if (steepness < 0) return 1;
    return 0;
}

std::pair<Eigen::Vector2d, double> processPoints(const Eigen::Vector2d& pointA, const Eigen::Vector2d& point_outside, double n) {
    Eigen::Vector2d entry_angle(1, 1);
    Eigen::Vector2d new_point(2, 1);
    auto delta_x = point_outside[0] - pointA[0];
    auto delta_y = point_outside[1] - pointA[1];
    auto steepness = delta_y / delta_x;
    int steepness_switch = calculateSteepnessSwitch(steepness);

    if (point_outside[0] < 0 || point_outside[0] > 1 || point_outside[1] < 0 || point_outside[1] > 1) {
        // oben oder rechts
        if (delta_x >= 0 && delta_y >= 0){
            double x = 1;
            double y = interpolateY(pointA, point_outside, x);
            // For the case that the point is on the edge
            if (y == 1) {
                y -= 0.0001;
            }
            else if (y == 0) {
                y += 0.0001;
            }
            // rechte Grenze passiert
            if (y < 1 && y > 0){
                Eigen::Vector2d exit_point(1, y);
                Eigen::Vector2d entry_point(y, 1);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y, because we also switched the entry coordinates before
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(0) *= steepness_switch;  // has to be variable
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n -= 90;
            }
            // obere Grenze passiert
            else {
                double y_back = 1;
                double x_back = interpolateX(pointA, point_outside, y_back);
                if (x_back == 1) {
                    x_back -= 0.0001;
                }
                else if (x_back == 0) {
                    x_back += 0.0001;
                }

                Eigen::Vector2d exit_point(x_back, 1);
                Eigen::Vector2d entry_point(1, x_back);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);

                entry_angle.row(1) *= steepness_switch;

                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n += 90;
            }
        }
        // unten oder rechts
        else if (delta_x > 0 && delta_y < 0){
            double x = 1;
            double y = interpolateY(pointA, point_outside, x);
            if (y == 1) {
                y -= 0.0001;
            }
            else if (y == 0) {
                y += 0.0001;
            }
            // rechte Grenze passiert
            if (y < 1 && y > 0){

                Eigen::Vector2d exit_point(1, y);
                Eigen::Vector2d entry_point(y, 1);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);

                entry_angle.row(0) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array()  * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n -= 90;
            }
            // unten Grenze passiert
            else {
                double y_back_neg = 0;
                double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);
                if (x_back_neg == 1) {
                    x_back_neg -= 0.0001;
                }
                else if (x_back_neg == 0) {
                    x_back_neg += 0.0001;
                }
                Eigen::Vector2d exit_point(x_back_neg, 0);
                Eigen::Vector2d entry_point(0, x_back_neg);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(1) *= steepness_switch;

                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n += 90;
            }
        }
        // oben oder links
        else if (delta_x < 0 && delta_y > 0){
            double x = 0;
            double y = interpolateY(pointA, point_outside, x);
            if (y == 1) {
                y -= 0.0001;
            }
            else if (y == 0) {
                y += 0.0001;
            }
            // linke Grenze passiert
            if (y < 1 && y > 0){
                Eigen::Vector2d exit_point(0, y);
                Eigen::Vector2d entry_point(y, 0);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(0) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n -= 90;
            }
            // obere Grenze passiert
            else {
                double y_back = 1;
                double x_back = interpolateX(pointA, point_outside, y_back);
                if (x_back == 1) {
                    x_back -= 0.0001;
                }
                else if (x_back == 0) {
                    x_back += 0.0001;
                }
                Eigen::Vector2d exit_point(x_back, 1);
                Eigen::Vector2d entry_point(1, x_back);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(1) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n += 90;
            }
        }
        // unten oder links
        else {
            double x = 0;
            double y = interpolateY(pointA, point_outside, x);
            if (y == 1) {
                y -= 0.0001;
            }
            else if (y == 0) {
                y += 0.0001;
            }
            // linke Grenze passiert
            if (y < 1 && y > 0){
                Eigen::Vector2d exit_point(0, y);
                Eigen::Vector2d entry_point(y, 0);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();

                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);

                entry_angle.row(0) *= steepness_switch;

                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n -= 90;
            }
            // unten Grenze passiert
            else {
                double y_back_neg = 0;
                double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);
                if (x_back_neg == 1) {
                    x_back_neg -= 0.0001;
                }
                else if (x_back_neg == 0) {
                    x_back_neg += 0.0001;
                }
                Eigen::Vector2d exit_point(x_back_neg, 0);
                Eigen::Vector2d entry_point(0, x_back_neg);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(1) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n += 90;
            }
        }
    }
    else {
        new_point = point_outside;
    }
    return std::make_pair(new_point, n);
}
