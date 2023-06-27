// process_points.h

#pragma once

#include <Eigen/Dense>

std::tuple<Eigen::Vector2d, double, Eigen::Vector2d> processPoints(const Eigen::Vector2d& pointA, const Eigen::Vector2d& point_outside, double n);