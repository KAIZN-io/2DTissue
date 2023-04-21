// update.h

#pragma once

#include <Eigen/Dense>
#include <vector>

#include <utilities/sim_structs.h>

std::vector<VertexData> update_vertex_data(
    const std::vector<int>& vertices_3D_active,
    const Eigen::VectorXd& vertice_3D_id,
    const std::vector<int>& inside_uv_ids
);

void update_if_valid(
    std::vector<VertexData>& vertex_data,
    const Eigen::MatrixXd& r_new,
    const Eigen::VectorXd& vertice_3D_id,
    int start_id
);