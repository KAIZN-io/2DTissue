// barycentric_coord.h

#pragma once

#include <Eigen/Dense>

Eigen::Vector3d calculate_barycentric_3D_coord(
    const Eigen::MatrixXd r,
    const Eigen::MatrixXd halfedges_uv,
    const Eigen::MatrixXi faces_uv,
    Eigen::MatrixXd vertices_uv,
    Eigen::MatrixXd vertices_3D,
    int interator
);