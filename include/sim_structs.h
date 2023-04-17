// sim_structs.h
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cstdint>

struct VertexData {
    int64_t old_id;
    int64_t next_id;
    bool valid;
    int uv_mesh_id;
};

struct Mesh_UV_Struct {
    int start_vertice_id;
    Eigen::MatrixXd mesh;
    std::vector<int64_t> h_v_mapping;
};