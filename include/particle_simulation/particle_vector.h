// particle_vector.h
#pragma once

#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd correct_n(
    const Eigen::MatrixXd& r_dot,
    const Eigen::MatrixXd n,
    double τ,
    double dt
);

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> calculate_particle_vectors(
    Eigen::MatrixXd &r_dot,
    Eigen::MatrixXd n,
    double dt
);