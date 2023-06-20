// author: @Jan-Piotraschke
// date: 2023-04-14
// license: Apache License 2.0
// version: 0.1.0

#include <Eigen/Dense>

#include <utilities/analytics.h>
#include <utilities/matrix_algebra.h>


void calculate_order_parameter(
    Eigen::VectorXd& v_order, 
    Eigen::MatrixXd r, 
    Eigen::MatrixXd r_dot, 
    int current_step
) {
    int num_part = r.rows();
    // Define a vector normal to position vector and velocity vector
    Eigen::MatrixXd v_tp = calculate_3D_cross_product(r, r_dot);

    // Normalize v_tp
    Eigen::MatrixXd v_norm = v_tp.rowwise().normalized();

    // Sum v_tp vectors and divide by number of particle to obtain order parameter of collective motion for spheroids
    v_order(current_step - 1) = (1.0 / num_part) * v_norm.colwise().sum().norm();
}
