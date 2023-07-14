// 2D_surface.h
#pragma once
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>

#include <utilities/mesh_descriptor.h>


void calculate_distances(
    _3D::Mesh mesh,
    _3D::vertex_descriptor start_node,
    std::vector<_3D::vertex_descriptor>& predecessor_pmap,
    std::vector<int>& distance
);

_3D::vertex_descriptor find_farthest_vertex(
    const _3D::Mesh mesh,
    _3D::vertex_descriptor start_node,
    const std::vector<int> distance
);

std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> create_uv_surface(
    std::string mesh_file_path,
    int32_t start_node_int
);

std::vector<_3D::edge_descriptor> set_UV_border_edges(
    const std::string mesh_file_path,
    _3D::vertex_descriptor start_node
);

std::string get_mesh_name(
   const std::string mesh_3D_path
);