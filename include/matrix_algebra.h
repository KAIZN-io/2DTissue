// matrix_algebra.h
#pragma once

#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd normalize_3D_matrix(const Eigen::MatrixXd &A);

Eigen::MatrixXd calculate_3D_cross_product(
    const Eigen::MatrixXd &A,
    const Eigen::MatrixXd &B
);