// update.h

#pragma once

#include <Eigen/Dense>
#include <vector>

#include <utilities/sim_structs.h>

std::vector<VertexData> update_vertex_data(
    const Eigen::MatrixXd& old_r_3D_coord,
    const Eigen::MatrixXd& new_r_3D_coord,
    const std::vector<int>& inside_uv_ids,
    int start_id
);

void update_if_valid(
    std::vector<VertexData>& vertex_data,
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_coord,
    const Eigen::MatrixXd& r_3D_coord,
    int start_id
);