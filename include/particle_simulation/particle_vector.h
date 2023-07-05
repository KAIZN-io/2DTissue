// particle_vector.h
#pragma once

#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd correct_n(
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot,
    const Eigen::MatrixXd n,
    double τ,
    double dt
);

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> calculate_particle_vectors(
    Eigen::Matrix<double, Eigen::Dynamic, 2> &r_dot,
    Eigen::MatrixXd n,
    double dt
);