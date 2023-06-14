// 2D_3D_mapping.h
#pragma once

#include <vector>
#include <cstdint>
#include <Eigen/Dense>

std::vector<int64_t> get_first_uv_halfedge_from_3D_vertice_id(
    const std::vector<int64_t>& _vertice_3D_id,
    const std::vector<int64_t>& _halfedge_vertices_mapping
);

Eigen::MatrixXd get_r_from_halfedge_id(
    const std::vector<int64_t>& halfedge_id,
    const Eigen::MatrixXd& halfedges_uv
);

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