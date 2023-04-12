// flight_of_the_particle.h
#pragma once

#include <tuple>
#include <vector>
#include <Eigen/Dense>

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> simulate_flight(
    Eigen::MatrixXd& r,
    Eigen::MatrixXd& n,
    std::vector<int>& vertices_3D_active,
    Eigen::MatrixXd distance_matrix_v,
    double v0,
    double k,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt
);