// motion.h
#pragma once

#include <tuple>
#include <vector>
#include <Eigen/Dense>

void transform_into_symmetric_matrix(Eigen::MatrixXd &A);

std::vector<Eigen::MatrixXd> get_dist_vect(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r);

double mean_unit_circle_vector_angle_degrees(std::vector<double> angles);

void calculate_average_n_within_distance(
    const std::vector<Eigen::MatrixXd> dist_vect,
    const Eigen::MatrixXd dist_length,
    Eigen::VectorXd& n,
    double σ
);

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> simulate_flight(
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