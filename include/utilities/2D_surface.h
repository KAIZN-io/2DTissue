// 2D_surface.h
#pragma once
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>

std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> create_uv_surface_intern(
    std::string mesh_3D = "Ellipsoid",
    int32_t start_node_int = 0
);