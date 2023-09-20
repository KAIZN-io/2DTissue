// Locomotion.h
#pragma once

#include <tuple>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <gtest/gtest_prod.h>

#include "LocomotionHelperInterface.h"
#include "LinearAlgebra.h"

class Locomotion {
public:
    Locomotion(
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
    void simulate_flight();
    static void get_distances_between_particles(
        Eigen::MatrixXd& dist_length,
        Eigen::MatrixXd distance_matrix,
        std::vector<int> vertice_3D_id
    );

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

FRIEND_TEST(LocomotionTest, ThrowsWhenInputIsEmpty);
FRIEND_TEST(LocomotionTest, CorrectlyCalculatesMeanAngle);
FRIEND_TEST(LocomotionTest, CorrectlyHandlesNegativeAngles);
FRIEND_TEST(LocomotionTest, SymmetricMatrixTestBasicTest);
FRIEND_TEST(LocomotionTest, SymmetricMatrixTestZeroTest);
FRIEND_TEST(LocomotionTest, SymmetricMatrixTestAllZerosTest);
FRIEND_TEST(LocomotionTest, GetDistVectTestBasicTest);
FRIEND_TEST(LocomotionTest, GetDistVectTestZeroMatrixTest);
FRIEND_TEST(LocomotionTest, GetDistVectTestOneDimensionTest);
FRIEND_TEST(LocomotionTest, GetDistVectTestHandlesSquareMatrixCorrectly);
FRIEND_TEST(LocomotionTest, GetDistVectTestTenDimensionTest);
FRIEND_TEST(LocomotionTest, AverageNWithinDistanceTest1);
FRIEND_TEST(LocomotionTest, RepulsiveAdhesionTest1);
FRIEND_TEST(LocomotionTest, RepulsiveAdhesionTest2);
};
