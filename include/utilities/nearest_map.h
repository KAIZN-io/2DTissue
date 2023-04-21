// nearest_map.h

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cstdint>

#include <utilities/sim_structs.h>


std::pair<Eigen::MatrixXd, std::vector<int64_t>> find_nearest_vertice_map(
    int target_vertex,
    const Eigen::MatrixXd distance_matrix,
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map
);