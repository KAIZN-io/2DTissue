// barycentric_coord.h

#pragma once

#include <Eigen/Dense>

std::pair<Eigen::Vector3d, int> calculate_barycentric_3D_coord(
    const Eigen::MatrixXd r,
    const Eigen::MatrixXd halfedges_uv,
    const Eigen::MatrixXi faces_uv,
    Eigen::MatrixXd vertices_uv,
    Eigen::MatrixXd vertices_3D,
    std::vector<int64_t> h_v_mapping,
    int interator
);

Eigen::Vector3d calculate_barycentric_2D_coord(
    const Eigen::MatrixXd& start_3D_points,
    const Eigen::MatrixXi& faces_3D_static,
    const Eigen::MatrixXd& vertices_uv,
    const Eigen::MatrixXd& vertices_3D,
    std::vector<int64_t>& h_v_mapping,
    int iterator
);
