// Simulator.h

#pragma once

#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <gtest/gtest_prod.h>
#include <utilities/angles_to_unit_vectors.h>

class Simulator {
public:
    Simulator(
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot,
        Eigen::VectorXd& n,
        std::vector<int>& vertices_3D_active,
        Eigen::MatrixXd& distance_matrix,
        Eigen::MatrixXd& dist_length,
        double v0,
        double k,
        double σ,
        double μ,
        double r_adh,
        double k_adh,
        double step_size
    );

    void simulate_flight();

private:
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot;
    Eigen::VectorXd& n;
    std::vector<int>& vertices_3D_active;
    Eigen::MatrixXd& distance_matrix;
    Eigen::MatrixXd& dist_length;
    double v0, k, σ, μ, r_adh, k_adh, step_size;

    Eigen::Matrix<double, Eigen::Dynamic, 2> F_track;

    void resize_F_track();
    static void transform_into_symmetric_matrix(Eigen::MatrixXd &A);
    static std::vector<Eigen::MatrixXd> get_dist_vect(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r);
    static void get_distances_between_particles(
        Eigen::MatrixXd& dist_length,
        Eigen::Matrix<double, Eigen::Dynamic, 2> r,
        Eigen::MatrixXd distance_matrix,
        std::vector<int> vertice_3D_id
    );
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
    void calculate_forces_between_particles(const std::vector<Eigen::MatrixXd> dist_vect);
    static double mean_unit_circle_vector_angle_degrees(std::vector<double> angles);

FRIEND_TEST(SimulatorTest, ThrowsWhenInputIsEmpty);
FRIEND_TEST(SimulatorTest, CorrectlyCalculatesMeanAngle);
FRIEND_TEST(SimulatorTest, CorrectlyHandlesNegativeAngles);
FRIEND_TEST(SimulatorTest, SymmetricMatrixTestBasicTest);
FRIEND_TEST(SimulatorTest, SymmetricMatrixTestZeroTest);
FRIEND_TEST(SimulatorTest, SymmetricMatrixTestAllZerosTest);
FRIEND_TEST(SimulatorTest, GetDistVectTestBasicTest);
FRIEND_TEST(SimulatorTest, GetDistVectTestZeroMatrixTest);
FRIEND_TEST(SimulatorTest, GetDistVectTestOneDimensionTest);
FRIEND_TEST(SimulatorTest, GetDistVectTestHandlesSquareMatrixCorrectly);
FRIEND_TEST(SimulatorTest, GetDistVectTestTenDimensionTest);
};
