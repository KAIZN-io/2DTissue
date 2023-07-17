// particle_vector.h
#pragma once

#include <vector>
#include <Eigen/Dense>

Eigen::VectorXd correct_n(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> r_dot,
    const Eigen::VectorXd n,
    double τ,
    double dt
);

// std::pair<Eigen::VectorXd, Eigen::VectorXd> calculate_particle_vectors(
//     Eigen::Matrix<double, Eigen::Dynamic, 2> &r_dot,
//     Eigen::VectorXd n,
//     double dt
// );