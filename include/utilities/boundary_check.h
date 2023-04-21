// boundary_check.h
#pragma once

#include <vector>
#include <Eigen/Dense>


// Find the indices of vertices that are inside the UV parametrization bounds
std::vector<int> find_inside_uv_vertices_id(const Eigen::MatrixXd& r);

std::vector<int> set_difference(int num_part, const std::vector<int>& inside_uv_ids);
