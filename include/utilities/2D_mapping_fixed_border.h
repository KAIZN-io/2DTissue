// 2D_mapping_fixed_border.h

#pragma once

#include <Eigen/Dense>

void opposite_seam_edges_square_border(Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_new);
void diagonal_seam_edges_square_border(
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_new,
    Eigen::MatrixXd& n_UV_new
);