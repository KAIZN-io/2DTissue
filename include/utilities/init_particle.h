// init_particle.h
#pragma once

#include <Eigen/Dense>

void init_particle_position(
    const Eigen::MatrixXi faces_uv,
    const Eigen::MatrixXd halfedges_uv,
    int num_part,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r,
    Eigen::MatrixXd& n
);