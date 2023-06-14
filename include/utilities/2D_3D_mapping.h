// 2D_3D_mapping.h
#pragma once

#include <vector>
#include <cstdint>
#include <Eigen/Dense>


std::pair<Eigen::MatrixXd, std::vector<int>> get_r3d(
    const Eigen::MatrixXd& r,
    const Eigen::MatrixXd& halfedges_uv,
    const Eigen::MatrixXi faces_uv,
    const Eigen::MatrixXd vertices_uv,
    const Eigen::MatrixXd vertices_3D,
    std::vector<int64_t> h_v_mapping
);

std::vector<int> find_vertice_rows_index(
    std::vector<int64_t> h_v_mapping_vector,
    std::vector<int> r3d_vertices
);

Eigen::MatrixXd get_coordinates(
    std::vector<int> indices,
    Eigen::MatrixXd coord
);