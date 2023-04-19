// init_particle_position.h
#pragma once

#include <Eigen/Dense>

void init_particle_position(
    const Eigen::MatrixXi faces_uv,
    const Eigen::MatrixXd halfedges_uv,
    int num_part,
    Eigen::MatrixXd& r,
    Eigen::MatrixXd& n
);