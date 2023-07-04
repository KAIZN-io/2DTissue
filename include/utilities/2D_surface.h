// 2D_surface.h
#pragma once
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>

std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> create_uv_surface(
    std::string mesh_file_path,
    int32_t start_node_int
);

std::string get_mesh_name(
   const std::string mesh_3D_path
);