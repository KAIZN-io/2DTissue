// process_points.h

#pragma once

#include <Eigen/Dense>

std::pair<Eigen::Vector2d, double> processPoints(const Eigen::Vector2d& pointA, const Eigen::Vector2d& point_outside, double n);