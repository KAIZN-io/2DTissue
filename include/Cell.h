// Cell.h
#pragma once

#include <vector>
#include <cstdint>
#include <Eigen/Dense>

#include <IO.h>

std::pair<Eigen::MatrixXd, std::vector<int>> get_r3d(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> r,
    const Eigen::MatrixXd halfedges_uv,
    const Eigen::MatrixXi faces_uv,
    const Eigen::MatrixXd vertices_uv,
    const Eigen::MatrixXd vertices_3D,
    std::vector<int64_t> h_v_mapping
);

Eigen::Matrix<double, Eigen::Dynamic, 2> get_r2d(
    const Eigen::MatrixXd r,
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

std::pair<Eigen::Vector3d, int> calculate_barycentric_3D_coord(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> r,
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

void init_particle_position(
    const Eigen::MatrixXi faces_uv,
    const Eigen::MatrixXd halfedges_uv,
    int num_part,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r,
    Eigen::VectorXd& n
);