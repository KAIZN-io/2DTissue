// sim_structs.h
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cstdint>

struct VertexData {
    Eigen::MatrixXd old_particle_pos;
    Eigen::MatrixXd next_particle_pos;
    bool valid;
    int uv_mesh_id;
};

struct Mesh_UV_Struct {
    int start_vertice_id;
    Eigen::MatrixXd mesh;
    std::vector<int64_t> h_v_mapping;
    Eigen::MatrixXd vertices_UV;
    Eigen::MatrixXd vertices_3D;
    std::string mesh_file_path;
};