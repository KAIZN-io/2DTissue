// forces.h
#pragma once

#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd calculate_forces_between_particles(
    const std::vector<Eigen::MatrixXd>& dist_vect,
    const Eigen::MatrixXd& dist_length,
    double k,
    double σ,
    double r_adh,
    double k_adh
);
