// author: @Jan-Piotraschke
// date: 2023-07-05
// license: Apache License 2.0
// version: 0.1.0

#include <cmath>
#include <Eigen/Dense>

#include <utilities/angles_to_unit_vectors.h>


/**
 * @brief Convert the angle degree to 2D unit vectors
 *
 * @info: Unittest implemented
*/
Eigen::Matrix<double, Eigen::Dynamic, 2> angles_to_unit_vectors(const Eigen::VectorXd& avg_n) {
    if (avg_n.cols() != 1) {
        throw std::invalid_argument("The input matrix must have exactly 1 column.");
    }

    // Initialize an Eigen::MatrixXd to store the 2D unit vectors
    Eigen::Matrix<double, Eigen::Dynamic, 2> n_vec(avg_n.rows(), 2);

    for (int i = 0; i < avg_n.rows(); ++i) {
        double angle_degrees = avg_n(i);
        double angle_radians = angle_degrees * M_PI / 180.0;

        // Convert the angle to a 2D unit vector
        Eigen::Vector2d vec(cos(angle_radians), sin(angle_radians));
        n_vec.row(i) = vec;
    }

    return n_vec;
}