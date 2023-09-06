/**
 * @file        LinearAlgebra.cpp
 * @brief       Mathematical operations on vectors and matrices
 *
 * @author      Jan-Piotraschke
 * @date        2023-Jul-19
 * @version     0.1.0
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        -
 */


#include <LinearAlgebra.h>


// ========================================
// Public Functions
// ========================================

/**
 * @brief Convert the angle degree to 2D unit vectors
 *
 * @info: Unittest implemented
*/
Eigen::Matrix<double, Eigen::Dynamic, 2> LinearAlgebra::angles_to_unit_vectors(const Eigen::VectorXi avg_n) {
    if (avg_n.cols() != 1) {
        throw std::invalid_argument("The input matrix must have exactly 1 column.");
    }

    // Initialize an Eigen::MatrixXd to store the 2D unit vectors
    Eigen::Matrix<double, Eigen::Dynamic, 2> n_vec(avg_n.rows(), 2);

    for (int i = 0; i < avg_n.rows(); ++i) {
        double angle_degrees = avg_n(i);
        double angle_radians = angle_degrees * DEG_TO_RAD;

        // Convert the angle to a 2D unit vector
        Eigen::Vector2d vec(cos(angle_radians), sin(angle_radians));
        n_vec.row(i) = vec;
    }

    return n_vec;
}

/*
normalize a 3D matrix A
*/
Eigen::MatrixXd LinearAlgebra::normalize_3D_matrix(const Eigen::MatrixXd A) {
    Eigen::VectorXd row_norms_t = A.rowwise().norm(); // Compute row-wise Euclidean norms
    Eigen::MatrixXd row_norms = row_norms_t.replicate(1, 3); // Repeat each norm for each column

    return row_norms; // Normalize each row
}


/**
 * Calculate the cross product of two 3D matrices A and B.
 *
 * @param A The first input matrix.
 * @param B The second input matrix.
 *
 * @return The cross product of A and B, computed for each row of the matrices.
 */
Eigen::MatrixXd LinearAlgebra::calculate_3D_cross_product(
    const Eigen::MatrixXd A,
    const Eigen::MatrixXd B
){
    // Ensure that A and B have the correct size
    assert(A.cols() == 3 && B.cols() == 3 && A.rows() == B.rows());

    // Preallocate output matrix with the same size and type as A
    const int num_rows = A.rows();
    Eigen::MatrixXd new_A = Eigen::MatrixXd::Zero(num_rows, 3);

    // Compute cross product for each row and directly assign the result to the output matrix
    for (int i = 0; i < num_rows; ++i) {
        // Get the i-th row of matrices A and B
        Eigen::Vector3d A_row = A.row(i);
        Eigen::Vector3d B_row = B.row(i);

        // Compute the cross product of the i-th rows of A and B
        new_A.row(i) = A_row.cross(B_row);
    }

    return new_A;
}


void LinearAlgebra::calculate_order_parameter(
    Eigen::VectorXd& v_order,
    Eigen::Matrix<double, Eigen::Dynamic, 2> r,
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_dot,
    int current_step
) {
    int num_part = r.rows();
    // Define a vector normal to position vector and velocity vector
    Eigen::MatrixXd v_tp = calculate_3D_cross_product(r, r_dot);

    // Normalize v_tp
    Eigen::MatrixXd v_norm = v_tp.rowwise().normalized();

    // Sum v_tp vectors and divide by number of particle to obtain order parameter of collective motion for spheroids
    v_order(current_step) = (1.0 / num_part) * v_norm.colwise().sum().norm();
}
