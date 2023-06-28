// error_checking.h

#pragma once

#include <Eigen/Dense>
#include <vector>

#include <utilities/sim_structs.h>

void error_unvalid_vertices(
    std::vector<VertexData> vertex_data
);

void error_invalid_values(
    Eigen::MatrixXd r_new
);

void error_lost_particles(
    Eigen::MatrixXd r_UV_new,
    int num_part
);