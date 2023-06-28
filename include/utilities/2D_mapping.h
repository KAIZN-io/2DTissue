// 2D_mapping.h

#pragma once

#include <Eigen/Dense>

void opposite_seam_edges_square_border(Eigen::MatrixXd& r_UV_new);
void diagonal_seam_edges_square_border(
    Eigen::MatrixXd r_UV,
    Eigen::MatrixXd& r_UV_new,
    Eigen::MatrixXd& n_UV_new
);