// 2D_mapping.h

#pragma once

#include <Eigen/Dense>

void opposite_seam_edges(Eigen::MatrixXd& r_UV_new);
void diagonal_seam_edges(
    Eigen::MatrixXd r_UV,
    Eigen::MatrixXd& r_UV_new
);