// splay_state.h

#pragma once

#include <tuple>
#include <vector>
#include <Eigen/Dense>

using Matrix3Xi = Eigen::Matrix<int, Eigen::Dynamic, 3>;

std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 2>, std::vector<int>> get_splay_state_vertices(
    const Matrix3Xi& mesh_loaded_uv,
    const Eigen::MatrixXd& halfedges_uv,
    int modula_mode = 10
);

std::vector<int> get_3D_splay_vertices(
    Eigen::MatrixXd distance_matrix,
    int modula_mode
);
