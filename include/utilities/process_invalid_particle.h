// process_invalid_particle.h
#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <map>
#include <vector>

#include "utilities/sim_structs.h"

void process_if_not_valid(
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
    std::vector<int> old_vertices_3D,
    std::vector<VertexData>& vertex_data,
    int num_part,
    Eigen::MatrixXd& distance_matrix_v,
    Eigen::VectorXd& n,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt
);