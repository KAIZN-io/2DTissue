// angles_to_unit_vectors.h

#pragma once

#include <Eigen/Dense>

Eigen::Matrix<double, Eigen::Dynamic, 2> angles_to_unit_vectors(const Eigen::VectorXd& avg_n);