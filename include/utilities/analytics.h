// analytics.h
#pragma once

#include "Eigen/Dense"

void calculate_order_parameter(
    Eigen::VectorXd& v_order, 
    Eigen::MatrixXd r, 
    Eigen::MatrixXd r_dot, 
    int tt
);