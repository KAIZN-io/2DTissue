// motion.h
#pragma once

#include <tuple>
#include <vector>
#include <Eigen/Dense>

std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 2>, Eigen::MatrixXd, Eigen::MatrixXd> simulate_flight(
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r,
    Eigen::VectorXd& n,
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