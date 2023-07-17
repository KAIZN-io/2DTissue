// Cell.h

#pragma once

#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <utilities/angles_to_unit_vectors.h>

class Cell {
public:
    Cell(
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot,
        Eigen::VectorXd& n,
        std::vector<int> vertices_3D_active,
        Eigen::MatrixXd distance_matrix_v,
        double v0,
        double k,
        double σ,
        double μ,
        double r_adh,
        double k_adh,
        double step_size
    );

    Eigen::MatrixXd simulate_flight();

private:
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot;
    Eigen::VectorXd& n;
    std::vector<int> vertices_3D_active;
    Eigen::MatrixXd distance_matrix_v;
    double v0, k, σ, μ, r_adh, k_adh, step_size;

    static void transform_into_symmetric_matrix(Eigen::MatrixXd &A);
    static std::vector<Eigen::MatrixXd> get_dist_vect(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r);
    static Eigen::MatrixXd get_distances_between_particles(
        Eigen::Matrix<double, Eigen::Dynamic, 2> r,
        Eigen::MatrixXd distance_matrix,
        std::vector<int> vertice_3D_id
    );
    static double mean_unit_circle_vector_angle_degrees(std::vector<double> angles);
    static void calculate_average_n_within_distance(
        const std::vector<Eigen::MatrixXd> dist_vect,
        const Eigen::MatrixXd dist_length,
        Eigen::VectorXd& n,
        double σ
    );
    Eigen::Vector2d repulsive_adhesion_motion(
        double k,
        double σ,
        double dist,
        double r_adh,
        double k_adh,
        const Eigen::Vector2d dist_v
    );
    Eigen::Matrix<double, Eigen::Dynamic, 2> calculate_forces_between_particles(
        const std::vector<Eigen::MatrixXd> dist_vect,
        const Eigen::MatrixXd dist_length,
        double k,
        double σ,
        double r_adh,
        double k_adh
    );
};