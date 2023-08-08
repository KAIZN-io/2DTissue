// Simulator.h

#pragma once

#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <gtest/gtest_prod.h>

#include "LinearAlgebra.h"

class Simulator {
public:
    Simulator(
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot,
        Eigen::VectorXi& n,
        std::vector<int>& vertices_3D_active,
        Eigen::MatrixXd& distance_matrix,
        Eigen::MatrixXd& dist_length,
        double v0,
        double k,
        double σ,
        double μ,
        double r_adh,
        double k_adh,
        double step_size,
        std::unique_ptr<LinearAlgebra> linear_algebra_ptr
    );
    static void get_distances_between_particles(
        Eigen::MatrixXd& dist_length,
        Eigen::MatrixXd distance_matrix,
        std::vector<int> vertice_3D_id
    );
    void simulate_flight();
    void opposite_seam_edges_square_border();
    void diagonal_seam_edges_square_border();

private:
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot;
    Eigen::VectorXi& n;
    std::vector<int>& vertices_3D_active;
    Eigen::MatrixXd& distance_matrix;
    Eigen::MatrixXd& dist_length;
    double v0, k, σ, μ, r_adh, k_adh, step_size;
    std::unique_ptr<LinearAlgebra> linear_algebra_ptr;

    Eigen::Matrix<double, Eigen::Dynamic, 2> F_track;

    void resize_F_track();
    static void transform_into_symmetric_matrix(Eigen::MatrixXd &A);
    static std::vector<Eigen::MatrixXd> get_dist_vect(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r);
    static void calculate_average_n_within_distance(
        const std::vector<Eigen::MatrixXd> dist_vect,
        const Eigen::MatrixXd dist_length,
        Eigen::VectorXi& n,
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
    double interpolateX(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double y);
    double interpolateY(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double x);
    int calculateSteepnessSwitch(double steepness);
    std::tuple<Eigen::Vector2d, double, Eigen::Vector2d> processPoints(const Eigen::Vector2d& pointA, const Eigen::Vector2d& point_outside, double n);
    void map_between_arbitrary_seam_edges();

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
FRIEND_TEST(SimulatorTest, AverageNWithinDistanceTest1);
FRIEND_TEST(SimulatorTest, RepulsiveAdhesionTest1);
FRIEND_TEST(SimulatorTest, RepulsiveAdhesionTest2);
};
