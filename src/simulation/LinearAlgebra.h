// LinearAlgebra.h
#pragma once

#include <vector>
#include <Eigen/Dense>


class LinearAlgebra
{
public:
    void calculate_order_parameter(
        Eigen::VectorXd& v_order,
        Eigen::Matrix<double, Eigen::Dynamic, 2> r,
        Eigen::Matrix<double, Eigen::Dynamic, 2> r_dot,
        int tt
    );

    Eigen::Matrix<double, Eigen::Dynamic, 2> angles_to_unit_vectors(const Eigen::VectorXi n);

    Eigen::MatrixXd normalize_3D_matrix(const Eigen::MatrixXd A);

    Eigen::MatrixXd calculate_3D_cross_product(
        const Eigen::MatrixXd A,
        const Eigen::MatrixXd B
    );
};