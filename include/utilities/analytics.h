// analytics.h
#pragma once

#include "Eigen/Dense"

void calculate_order_parameter(
    Eigen::VectorXd& v_order, 
    Eigen::Matrix<double, Eigen::Dynamic, 2> r, 
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_dot, 
    int tt
);