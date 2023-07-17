// dye_particle.h
#pragma once

#include <Eigen/Core>
#include <vector>

void count_particle_neighbors(
    std::vector<int>& particles_color,
    const Eigen::MatrixXd dist_length,
    double σ
);